# -*- coding: utf-8 -*-
"""
Elliptic curves over the rational numbers

AUTHORS:

- William Stein (2005): first version

- William Stein (2006-02-26): fixed Lseries_extended which didn't work
  because of changes elsewhere in Sage.

- David Harvey (2006-09): Added padic_E2, padic_sigma, padic_height,
  padic_regulator methods.

- David Harvey (2007-02): reworked padic-height related code

- Christian Wuthrich (2007): added padic sha computation

- David Roe (2007-09): moved sha, l-series and p-adic functionality to
  separate files.

- John Cremona (2008-01)

- Tobias Nagel and Michael Mardaus (2008-07): added integral_points

- John Cremona (2008-07): further work on integral_points

- Christian Wuthrich (2010-01): moved Galois reps and modular
  parametrization in a separate file

- Simon Spicer (2013-03): Added code for modular degrees and congruence
  numbers of higher level

"""

##############################################################################
#       Copyright (C) 2005,2006,2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

import constructor
import BSD
from   ell_generic import EllipticCurve_generic, is_EllipticCurve
import ell_modular_symbols
from   ell_number_field import EllipticCurve_number_field
import ell_point
import ell_tate_curve
import ell_torsion
import heegner
from   gp_simon import simon_two_descent
from   lseries_ell import Lseries_ell
import mod5family
from   modular_parametrization import ModularParameterization
import padic_lseries
import padics

from sage.modular.modsym.modsym import ModularSymbols

import sage.modular.modform.constructor
import sage.modular.modform.element
import sage.libs.mwrank.all as mwrank
import sage.databases.cremona

import sage.rings.arith as arith
import sage.rings.all as rings
from sage.rings.all import (
    PowerSeriesRing, LaurentSeriesRing, O,
    infinity as oo,
    ZZ, QQ,
    Integer,
    IntegerRing, RealField,
    ComplexField, RationalField)

import sage.misc.misc as misc
from sage.misc.all import verbose

from sage.misc.functional import log
from sage.sets.set import Set

# Use some interval arithmetic to guarantee correctness.  We assume
# that alpha is computed to the precision of a float.
# IR = rings.RIF
#from sage.rings.interval import IntervalRing; IR = IntervalRing()

import sage.matrix.all as matrix
from   sage.libs.pari.all import pari, PariError
from sage.functions.other import gamma_inc
from math import sqrt
from sage.interfaces.all import gp
from sage.misc.cachefunc import cached_method
from copy import copy

factor = arith.factor
mul = misc.mul
next_prime = arith.next_prime

Q = RationalField()
C = ComplexField()
R = RealField()
Z = IntegerRing()
IR = rings.RealIntervalField(20)

_MAX_HEIGHT=21

# complex multiplication dictionary:
# CMJ is a dict of pairs (j,D) where j is a rational CM j-invariant
# and D is the corresponding quadratic discriminant

CMJ={ 0: -3, 54000: -12, -12288000: -27, 1728: -4, 287496: -16,
      -3375: -7, 16581375: -28, 8000: -8, -32768: -11, -884736: -19,
      -884736000: -43, -147197952000: -67, -262537412640768000: -163}



class EllipticCurve_rational_field(EllipticCurve_number_field):
    r"""
    Elliptic curve over the Rational Field.

    INPUT:

    - ``ainvs`` (list or string) -- either `[a_1,a_2,a_3,a_4,a_6]` or
      `[a_4,a_6]` (with `a_1=a_2=a_3=0`) or a valid label from the
      database.

    .. note::

       See constructor.py for more variants.

    EXAMPLES:

    Construction from Weierstrass coefficients (`a`-invariants), long form::

        sage: E = EllipticCurve([1,2,3,4,5]); E
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field

    Construction from Weierstrass coefficients (`a`-invariants),
    short form (sets `a_1=a_2=a_3=0`)::

        sage: EllipticCurve([4,5]).ainvs()
        (0, 0, 0, 4, 5)

    Construction from a database label::

        sage: EllipticCurve('389a1')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    """
    def __init__(self, ainvs, extra=None):
        r"""
        Constructor for the EllipticCurve_rational_field class.

        INPUT:

        - ``ainvs`` (list or string) -- either `[a_1,a_2,a_3,a_4,a_6]`
          or `[a_4,a_6]` (with `a_1=a_2=a_3=0`) or a valid label from
          the database.

        .. note::

           See constructor.py for more variants.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,4,5]); E
            Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field

        Constructor from `[a_4,a_6]` sets `a_1=a_2=a_3=0`::

            sage: EllipticCurve([4,5]).ainvs()
            (0, 0, 0, 4, 5)

        Constructor from a Cremona label::

            sage: EllipticCurve('389a1')
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

        Constructor from an LMFDB label::

            sage: EllipticCurve('462.f3')
            Elliptic Curve defined by y^2 + x*y = x^3 - 363*x + 1305 over Rational Field

        TESTS:

        When constructing a curve from the large database using a
        label, we must be careful that the copied generators have the
        right curve (see #10999: the following used not to work when
        the large database was installed)::

            sage: E=EllipticCurve('389a1')
            sage: [P.curve() is E for P in E.gens()]
            [True, True]

        """
        if extra != None:   # possibility of two arguments (the first would be the field)
            ainvs = extra
        self.__np = {}
        self.__gens = {}
        self.__rank = {}
        self.__regulator = {}
        self.__generalized_modular_degree = {}
        self.__generalized_congruence_number = {}
        self._isoclass = {}
        if isinstance(ainvs, str):
            label = ainvs
            X = sage.databases.cremona.CremonaDatabase()[label]
            EllipticCurve_number_field.__init__(self, Q, list(X.a_invariants()))
            for attr in ['rank', 'torsion_order', 'cremona_label', 'conductor',
                         'modular_degree', 'gens', 'regulator']:
                s = "_EllipticCurve_rational_field__"+attr
                if hasattr(X,s):
                    if attr == 'gens': # see #10999
                        gens_dict = getattr(X, s)
                        for boo in gens_dict.keys():
                            gens_dict[boo] = [self(P) for P in gens_dict[boo]]
                            setattr(self, s, gens_dict)
                    else:
                        setattr(self, s, getattr(X, s))
            if hasattr(X,'_lmfdb_label'):
                self._lmfdb_label = X._lmfdb_label
            return
        EllipticCurve_number_field.__init__(self, Q, ainvs)
        if self.base_ring() != Q:
            raise TypeError, "Base field (=%s) must be the Rational Field."%self.base_ring()

    def _set_rank(self, r):
        """
        Internal function to set the cached rank of this elliptic curve to
        r.

        .. warning::

           No checking is done! Not intended for use by users.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: E._set_rank(99)  # bogus value -- not checked
            sage: E.rank()         # returns bogus cached value
            99
            sage: E._EllipticCurve_rational_field__rank={} # undo the damage
            sage: E.rank()         # the correct rank
            1
        """
        self.__rank = {}
        self.__rank[True] = Integer(r)

    def _set_torsion_order(self, t):
        """
        Internal function to set the cached torsion order of this elliptic
        curve to t.

        .. warning::

           No checking is done! Not intended for use by users.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E._set_torsion_order(99)  # bogus value -- not checked
            sage: E.torsion_order()         # returns bogus cached value
            99
            sage: T = E.torsion_subgroup()  # causes actual torsion to be computed
            sage: E.torsion_order()         # the correct value
            1
        """
        self.__torsion_order = Integer(t)

    def _set_cremona_label(self, L):
        """
        Internal function to set the cached label of this elliptic curve to
        L.

        .. warning::

           No checking is done! Not intended for use by users.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E._set_cremona_label('bogus')
            sage: E.label()
            'bogus'
            sage: E.database_curve().label()
            '37a1'
            sage: E.label() # no change
            'bogus'
            sage: E._set_cremona_label(E.database_curve().label())
            sage: E.label() # now it is correct
            '37a1'
        """
        self.__cremona_label = L

    def _set_conductor(self, N):
        """
        Internal function to set the cached conductor of this elliptic
        curve to N.

        .. warning::

           No checking is done! Not intended for use by users.
           Setting to the wrong value will cause strange problems (see
           examples).

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E._set_conductor(99)      # bogus value -- not checked
            sage: E.conductor()             # returns bogus cached value
            99
        """
        self.__conductor_pari = Integer(N)

    def _set_modular_degree(self, deg):
        """
        Internal function to set the cached modular degree of this elliptic
        curve to deg.

        .. warning::

           No checking is done!

        EXAMPLES::

            sage: E=EllipticCurve('5077a1')
            sage: E.modular_degree()
            1984
            sage: E._set_modular_degree(123456789)
            sage: E.modular_degree()
            123456789
            sage: E._set_modular_degree(E.database_curve().modular_degree())
            sage: E.modular_degree()
            1984
        """
        self.__modular_degree = Integer(deg)

    def _set_gens(self, gens):
        """
        Internal function to set the cached generators of this elliptic
        curve to gens.

        .. warning::

           No checking is done!

        EXAMPLES::

            sage: E=EllipticCurve('5077a1')
            sage: E.rank()
            3
            sage: E.gens() # random
            [(-2 : 3 : 1), (-7/4 : 25/8 : 1), (1 : -1 : 1)]
            sage: E._set_gens([]) # bogus list
            sage: E.rank()        # unchanged
            3
            sage: E._set_gens([E(-2,3), E(-1,3), E(0,2)])
            sage: E.gens()
            [(-2 : 3 : 1), (-1 : 3 : 1), (0 : 2 : 1)]
        """
        self.__gens = {}
        self.__gens[True] = [self.point(x, check=True) for x in gens]
        self.__gens[True].sort()


    def is_p_integral(self, p):
        r"""
        Returns True if this elliptic curve has `p`-integral
        coefficients.

        INPUT:


        -  ``p`` - a prime integer


        EXAMPLES::

            sage: E=EllipticCurve(QQ,[1,1]); E
            Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
            sage: E.is_p_integral(2)
            True
            sage: E2=E.change_weierstrass_model(2,0,0,0); E2
            Elliptic Curve defined by y^2 = x^3 + 1/16*x + 1/64 over Rational Field
            sage: E2.is_p_integral(2)
            False
            sage: E2.is_p_integral(3)
            True
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        if self.is_integral():
            return True
        return bool(misc.mul([x.valuation(p) >= 0 for x in self.ainvs()]))

    def is_integral(self):
        """
        Returns True if this elliptic curve has integral coefficients (in
        Z)

        EXAMPLES::

            sage: E=EllipticCurve(QQ,[1,1]); E
            Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
            sage: E.is_integral()
            True
            sage: E2=E.change_weierstrass_model(2,0,0,0); E2
            Elliptic Curve defined by y^2 = x^3 + 1/16*x + 1/64 over Rational Field
            sage: E2.is_integral()
            False
        """
        try:
            return self.__is_integral
        except AttributeError:
            one = Integer(1)
            self.__is_integral = bool(misc.mul([x.denominator() == 1 for x in self.ainvs()]))
            return self.__is_integral


    def mwrank(self, options=''):
        r"""
        Run Cremona's mwrank program on this elliptic curve and return the
        result as a string.

        INPUT:


        -  ``options`` (string) -- run-time options passed when starting mwrank.
           The format is as follows (see below for examples of usage):

           - ``-v n``    (verbosity level)       sets verbosity to n (default=1)
           - ``-o``      (PARI/GP style output flag)  turns ON extra PARI/GP short output (default is OFF)
           - ``-p n``    (precision)       sets precision to `n` decimals (default=15)
           - ``-b n``    (quartic bound)   bound on quartic point search (default=10)
           - ``-x n``    (n_aux)           number of aux primes used for sieving (default=6)
           - ``-l``      (generator list flag)            turns ON listing of points (default ON unless v=0)
           - ``-s``      (selmer_only flag)     if set, computes Selmer rank only (default: not set)
           - ``-d``      (skip_2nd_descent flag)        if set, skips the second descent for curves with 2-torsion (default: not set)
           - ``-S n``    (sat_bd)          upper bound on saturation primes (default=100, -1 for automatic)

        OUTPUT:

        -  ``string`` - output of mwrank on this curve


        .. note::

           The output is a raw string and completely illegible using
           automatic display, so it is recommended to use print for
           legible output.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: E.mwrank() #random
            ...
            sage: print E.mwrank()
            Curve [0,0,1,-1,0] :        Basic pair: I=48, J=-432
            disc=255744
            ...
            Generator 1 is [0:-1:1]; height 0.05111...

            Regulator = 0.05111...

            The rank and full Mordell-Weil basis have been determined unconditionally.
            ...

        Options to mwrank can be passed::

            sage: E = EllipticCurve([0,0,0,877,0])

        Run mwrank with 'verbose' flag set to 0 but list generators if
        found

        ::

            sage: print E.mwrank('-v0 -l')
            Curve [0,0,0,877,0] :   0 <= rank <= 1
            Regulator = 1

        Run mwrank again, this time with a higher bound for point searching
        on homogeneous spaces::

            sage: print E.mwrank('-v0 -l -b11')
            Curve [0,0,0,877,0] :   Rank = 1
            Generator 1 is [29604565304828237474403861024284371796799791624792913256602210:-256256267988926809388776834045513089648669153204356603464786949:490078023219787588959802933995928925096061616470779979261000]; height 95.980371987964
            Regulator = 95.980371987964
        """
        if options == "":
            from sage.interfaces.all import mwrank
        else:
            from sage.interfaces.all import Mwrank
            mwrank = Mwrank(options=options)
        return mwrank(list(self.a_invariants()))

    def conductor(self, algorithm="pari"):
        """
        Returns the conductor of the elliptic curve.

        INPUT:


        -  ``algorithm`` - str, (default: "pari")

           -  ``"pari"`` - use the PARI C-library ellglobalred
              implementation of Tate's algorithm

           -  ``"mwrank"`` - use Cremona's mwrank implementation
              of Tate's algorithm; can be faster if the curve has integer
              coefficients (TODO: limited to small conductor until mwrank gets
              integer factorization)

           -  ``"gp"`` - use the GP interpreter.

           -  ``"generic"`` - use the general number field
              implementation

           -  ``"all"`` - use all four implementations, verify
              that the results are the same (or raise an error), and output the
              common value.


        EXAMPLE::

            sage: E = EllipticCurve([1, -1, 1, -29372, -1932937])
            sage: E.conductor(algorithm="pari")
            3006
            sage: E.conductor(algorithm="mwrank")
            3006
            sage: E.conductor(algorithm="gp")
            3006
            sage: E.conductor(algorithm="generic")
            3006
            sage: E.conductor(algorithm="all")
            3006

        .. note::

           The conductor computed using each algorithm is cached
           separately. Thus calling ``E.conductor('pari')``, then
           ``E.conductor('mwrank')`` and getting the same result
           checks that both systems compute the same answer.
        """

        if algorithm == "pari":
            try:
                return self.__conductor_pari
            except AttributeError:
                self.__conductor_pari = Integer(self.pari_mincurve().ellglobalred()[0])
            return self.__conductor_pari

        elif algorithm == "gp":
            try:
                return self.__conductor_gp
            except AttributeError:
                self.__conductor_gp = Integer(gp.eval('ellglobalred(ellinit(%s,0))[1]'%list(self.a_invariants())))
                return self.__conductor_gp

        elif algorithm == "mwrank":
            try:
                return self.__conductor_mwrank
            except AttributeError:
                if self.is_integral():
                    self.__conductor_mwrank = Integer(self.mwrank_curve().conductor())
                else:
                    self.__conductor_mwrank = Integer(self.minimal_model().mwrank_curve().conductor())
            return self.__conductor_mwrank

        elif algorithm == "generic":
            try:
                return self.__conductor_generic
            except AttributeError:
                self.__conductor_generic = sage.schemes.elliptic_curves.ell_number_field.EllipticCurve_number_field.conductor(self).gen()
                return self.__conductor_generic

        elif algorithm == "all":
            N1 = self.conductor("pari")
            N2 = self.conductor("mwrank")
            N3 = self.conductor("gp")
            N4 = self.conductor("generic")
            if N1 != N2 or N2 != N3 or N2 != N4:
                raise ArithmeticError, "PARI, mwrank, gp and Sage compute different conductors (%s,%s,%s,%3) for %s"%(
                    N1, N2, N3, N4, self)
            return N1
        else:
            raise RuntimeError, "algorithm '%s' is not known."%algorithm

    ####################################################################
    #  Access to PARI curves related to this curve.
    ####################################################################

    def pari_curve(self, prec = None, factor = 1):
        """
        Return the PARI curve corresponding to this elliptic curve.

        INPUT:


        -  ``prec`` - The precision of quantities calculated for the
           returned curve, in bits. If None, defaults to factor
           multiplied by the precision of the largest cached curve (or
           the default real precision if none yet computed).

        -  ``factor`` - The factor by which to increase the
           precision over the maximum previously computed precision. Only used
           if prec (which gives an explicit precision) is None.


        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: e = E.pari_curve()
            sage: type(e)
            <type 'sage.libs.pari.gen.gen'>
            sage: e.type()
            't_VEC'
            sage: e.ellan(10)
            [1, -2, -3, 2, -2, 6, -1, 0, 6, 4]

        ::

            sage: E = EllipticCurve(RationalField(), ['1/3', '2/3'])
            sage: e = E.pari_curve(prec = 100)
            sage: E._pari_curve.has_key(100)
            True
            sage: e.type()
            't_VEC'
            sage: e[:5]
            [0, 0, 0, 1/3, 2/3]

        This shows that the bug uncovered by trac #3954 is fixed::

            sage: E._pari_curve.has_key(100)
            True

        ::

            sage: E = EllipticCurve('37a1').pari_curve()
            sage: E[14].python().prec()
            64
            sage: [a.precision() for a in E]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4] # 32-bit
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3] # 64-bit

        This shows that the bug uncovered by trac #4715 is fixed::

            sage: Ep = EllipticCurve('903b3').pari_curve()
        """
        try:
            # if the PARI curve has already been computed to this
            # precision, returned the cached copy
            return self._pari_curve[prec]
        except AttributeError:
            # no PARI curves have been computed for this elliptic curve
            self._pari_curve = {}
            if prec is None:
                prec = rings.RealField().precision()
        except KeyError:
            # PARI curves are cached for this elliptic curve, but they
            # are not of the requested precision (or prec = None)
            if prec is None:
                L = self._pari_curve.keys()
                L.sort()
                if factor == 1:
                    return self._pari_curve[L[-1]]
                else:
                    prec = int(factor * L[-1])
        self._pari_curve[prec] = pari(list(self.a_invariants())).ellinit(precision=prec)
        return self._pari_curve[prec]

    def pari_mincurve(self, prec = None, factor = 1):
        """
        Return the PARI curve corresponding to a minimal model for this
        elliptic curve.

        INPUT:


        -  ``prec`` - The precision of quantities calculated for the
           returned curve, in bits. If None, defaults to factor
           multiplied by the precision of the largest cached curve (or
           the default real precision if none yet computed).

        -  ``factor`` - The factor by which to increase the
           precision over the maximum previously computed precision. Only used
           if prec (which gives an explicit precision) is None.


        EXAMPLES::

            sage: E = EllipticCurve(RationalField(), ['1/3', '2/3'])
            sage: e = E.pari_mincurve()
            sage: e[:5]
            [0, 0, 0, 27, 486]
            sage: E.conductor()
            47232
            sage: e.ellglobalred()
            [47232, [1, 0, 0, 0], 2]
        """
        try:
            # if the PARI curve has already been computed to this
            # precision, returned the cached copy
            return self._pari_mincurve[prec]
        except AttributeError:
            # no PARI curves have been computed for this elliptic curve
            self._pari_mincurve = {}
        except KeyError:
            # PARI curves are cached for this elliptic curve, but they
            # are not of the requested precision (or prec = None)
            if prec is None:
                L = self._pari_mincurve.keys()
                L.sort()
                if factor == 1:
                    return self._pari_mincurve[L[-1]]
                else:
                    prec = int(factor * L[-1])
        e = self.pari_curve(prec)
        mc, change = e.ellminimalmodel()
        self._pari_mincurve[prec] = mc
        # self.__min_transform = change
        return mc

    def database_curve(self):
        """
        Return the curve in the elliptic curve database isomorphic to this
        curve, if possible. Otherwise raise a RuntimeError exception.

        EXAMPLES::

            sage: E = EllipticCurve([0,1,2,3,4])
            sage: E.database_curve()
            Elliptic Curve defined by y^2  = x^3 + x^2 + 3*x + 5 over Rational Field

        .. note::

           The model of the curve in the database can be different
           from the Weierstrass model for this curve, e.g., database
           models are always minimal.
        """
        try:
            return self.__database_curve
        except AttributeError:
            misc.verbose("Looking up %s in the database."%self)
            D = sage.databases.cremona.CremonaDatabase()
            ainvs = list(self.minimal_model().ainvs())
            try:
                self.__database_curve = D.elliptic_curve_from_ainvs(ainvs)
            except RuntimeError:
                raise RuntimeError, "Elliptic curve %s not in the database."%self
            return self.__database_curve


    def Np(self, p):
        r"""
        The number of points on `E` modulo `p`.

        INPUT:

        - ``p`` (int) -- a prime, not necessarily of good reduction.


        OUTPUT:

        (int) The number ofpoints on the reduction of `E` modulo `p`
        (including the singular point when `p` is a prime of bad
        reduction).

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.Np(2)
            5
            sage: E.Np(3)
            5
            sage: E.conductor()
            11
            sage: E.Np(11)
            11

        This even works when the prime is large::

            sage: E = EllipticCurve('37a')
            sage: E.Np(next_prime(10^30))
            1000000000000001426441464441649
        """
        if self.conductor() % p == 0:
            return p + 1 - self.ap(p)
        return p+1 - self.ap(p)

    #def __pari_double_prec(self):
    #    EllipticCurve_number_field._EllipticCurve__pari_double_prec(self)
    #    try:
    #        del self._pari_mincurve
    #    except AttributeError:
    #        pass

    ####################################################################
    #  Access to mwrank
    ####################################################################
    def mwrank_curve(self, verbose=False):
        """
        Construct an mwrank_EllipticCurve from this elliptic curve

        The resulting mwrank_EllipticCurve has available methods from John
        Cremona's eclib library.

        EXAMPLES::

            sage: E=EllipticCurve('11a1')
            sage: EE=E.mwrank_curve()
            sage: EE
            y^2+ y = x^3 - x^2 - 10*x - 20
            sage: type(EE)
            <class 'sage.libs.mwrank.interface.mwrank_EllipticCurve'>
            sage: EE.isogeny_class()
            ([[0, -1, 1, -10, -20], [0, -1, 1, -7820, -263580], [0, -1, 1, 0, 0]],
            [[0, 5, 5], [5, 0, 0], [5, 0, 0]])
        """
        try:
            return self.__mwrank_curve
        except AttributeError:
            pass
        self.__mwrank_curve = mwrank.mwrank_EllipticCurve(
            list(self.ainvs()), verbose=verbose)
        return self.__mwrank_curve

    def two_descent(self, verbose=True,
                    selmer_only = False,
                    first_limit = 20,
                    second_limit = 8,
                    n_aux = -1,
                    second_descent = 1):
        """
        Compute 2-descent data for this curve.

        INPUT:


        -  ``verbose`` - (default: True) print what mwrank is
           doing. If False, **no output** is printed.

        -  ``selmer_only`` - (default: False) selmer_only
           switch

        -  ``first_limit`` - (default: 20) firstlim is bound
           on x+z second_limit- (default: 8) secondlim is bound on log max
           x,z , i.e. logarithmic

        -  ``n_aux`` - (default: -1) n_aux only relevant for
           general 2-descent when 2-torsion trivial; n_aux=-1 causes default
           to be used (depends on method)

        -  ``second_descent`` - (default: True)
           second_descent only relevant for descent via 2-isogeny


        OUTPUT:

        Returns ``True`` if the descent succeeded, i.e. if the lower bound and
        the upper bound for the rank are the same. In this case, generators and
        the rank are cached. A return value of ``False`` indicates that either
        rational points were not found, or that Sha[2] is nontrivial and mwrank
        was unable to determine this for sure.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.two_descent(verbose=False)
            True

        """
        misc.verbose("Calling mwrank C++ library.")
        C = self.mwrank_curve()
        C.two_descent(verbose, selmer_only,
                        first_limit, second_limit,
                        n_aux, second_descent)
        if C.certain():
            self.__gens[True] = [self.point(x, check=True) for x in C.gens()]
            self.__gens[True].sort()
            self.__rank[True] = len(self.__gens[True])
        return C.certain()

    ####################################################################
    #  Etc.
    ####################################################################

    def aplist(self, n, python_ints=False):
        r"""
        The Fourier coefficients `a_p` of the modular form
        attached to this elliptic curve, for all primes `p\leq n`.

        INPUT:


        -  ``n`` - integer

        -  ``python_ints`` - bool (default: False); if True
           return a list of Python ints instead of Sage integers.


        OUTPUT: list of integers

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: e.aplist(1)
            []
            sage: e.aplist(2)
            [-2]
            sage: e.aplist(10)
            [-2, -3, -2, -1]
            sage: v = e.aplist(13); v
            [-2, -3, -2, -1, -5, -2]
            sage: type(v[0])
            <type 'sage.rings.integer.Integer'>
            sage: type(e.aplist(13, python_ints=True)[0])
            <type 'int'>
        """
        e = self.pari_mincurve()
        v = e.ellaplist(n, python_ints=True)
        if python_ints:
            return v
        else:
            return [Integer(a) for a in v]



    def anlist(self, n, python_ints=False):
        """
        The Fourier coefficients up to and including `a_n` of the
        modular form attached to this elliptic curve. The i-th element of
        the return list is a[i].

        INPUT:


        -  ``n`` - integer

        -  ``python_ints`` - bool (default: False); if True
           return a list of Python ints instead of Sage integers.


        OUTPUT: list of integers

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.anlist(3)
            [0, 1, -2, -1]

        ::

            sage: E = EllipticCurve([0,1])
            sage: E.anlist(20)
            [0, 1, 0, 0, 0, 0, 0, -4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 8, 0]
        """
        n = int(n)
        e = self.pari_mincurve()
        if n >= 2147483648:
            raise RuntimeError, "anlist: n (=%s) must be < 2147483648."%n

        v = [0] + e.ellan(n, python_ints=True)
        if not python_ints:
            v = [Integer(x) for x in v]
        return v


        # There is some overheard associated with coercing the PARI
        # list back to Python, but it's not bad.  It's better to do it
        # this way instead of trying to eval the whole list, since the
        # int conversion is done very sensibly.  NOTE: This would fail
        # if a_n won't fit in a C int, i.e., is bigger than
        # 2147483648; however, we wouldn't realistically compute
        # anlist for n that large anyway.
        #
        # Some relevant timings:
        #
        # E <--> [0, 1, 1, -2, 0]   389A
        #  E = EllipticCurve([0, 1, 1, -2, 0]);   // Sage or MAGMA
        #  e = E.pari_mincurve()
        #  f = ellinit([0,1,1,-2,0]);
        #
        #  Computation                                              Time (1.6Ghz Pentium-4m laptop)
        #  time v:=TracesOfFrobenius(E,10000);  // MAGMA            0.120
        #  gettime;v=ellan(f,10000);gettime/1000                    0.046
        #  time v=e.ellan (10000)                                   0.04
        #  time v=E.anlist(10000)                                   0.07

        #  time v:=TracesOfFrobenius(E,100000);  // MAGMA           1.620
        #  gettime;v=ellan(f,100000);gettime/1000                   0.676
        #  time v=e.ellan (100000)                                  0.7
        #  time v=E.anlist(100000)                                  0.83

        #  time v:=TracesOfFrobenius(E,1000000);  // MAGMA          20.850
        #  gettime;v=ellan(f,1000000);gettime/1000                  9.238
        #  time v=e.ellan (1000000)                                 9.61
        #  time v=E.anlist(1000000)                                 10.95  (13.171 in cygwin vmware)

        #  time v:=TracesOfFrobenius(E,10000000);  //MAGMA          257.850
        #  gettime;v=ellan(f,10000000);gettime/1000      FAILS no matter how many allocatemem()'s!!
        #  time v=e.ellan (10000000)                                139.37
        #  time v=E.anlist(10000000)                                136.32
        #
        #  The last Sage comp retries with stack size 40MB,
        #  80MB, 160MB, and succeeds last time.  It's very interesting that this
        #  last computation is *not* possible in GP, but works in py_pari!
        #

    def q_expansion(self, prec):
        r"""
        Return the `q`-expansion to precision prec of the newform
        attached to this elliptic curve.

        INPUT:


        -  ``prec`` - an integer


        OUTPUT:

        a power series (in the variable 'q')

        .. note::

           If you want the output to be a modular form and not just a
           `q`-expansion, use :meth:`.modular_form`.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.q_expansion(20)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + 4*q^10 - 5*q^11 - 6*q^12 - 2*q^13 + 2*q^14 + 6*q^15 - 4*q^16 - 12*q^18 + O(q^20)
        """
        return PowerSeriesRing(Q, 'q')(self.anlist(prec), prec, check=True)

    def modular_form(self):
        r"""
        Return the cuspidal modular form associated to this elliptic
        curve.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: f = E.modular_form()
            sage: f
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)

        If you need to see more terms in the `q`-expansion::

            sage: f.q_expansion(20)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + 4*q^10 - 5*q^11 - 6*q^12 - 2*q^13 + 2*q^14 + 6*q^15 - 4*q^16 - 12*q^18 + O(q^20)

        .. note::

           If you just want the `q`-expansion, use
           :meth:`.q_expansion`.
        """
        try:
            return self.__modular_form
        except AttributeError:
            M = sage.modular.modform.constructor.ModularForms(self.conductor(),weight=2)
            f = sage.modular.modform.element.ModularFormElement_elliptic_curve(M, self)
            self.__modular_form = f
            return f

    def modular_symbol_space(self, sign=1, base_ring=Q, bound=None):
        r"""
        Return the space of cuspidal modular symbols associated to this
        elliptic curve, with given sign and base ring.

        INPUT:


        -  ``sign`` - 0, -1, or 1

        -  ``base_ring`` - a ring


        EXAMPLES::

            sage: f = EllipticCurve('37b')
            sage: f.modular_symbol_space()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(37) of weight 2 with sign 1 over Rational Field
            sage: f.modular_symbol_space(-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(37) of weight 2 with sign -1 over Rational Field
            sage: f.modular_symbol_space(0, bound=3)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field

        .. note::

           If you just want the `q`-expansion, use
           :meth:`.q_expansion`.
        """
        typ = (sign, base_ring)
        try:
            return self.__modular_symbol_space[typ]
        except AttributeError:
            self.__modular_symbol_space = {}
        except KeyError:
            pass
        M = ell_modular_symbols.modular_symbol_space(self, sign, base_ring, bound=bound)
        self.__modular_symbol_space[typ] = M
        return M

    def modular_symbol(self, sign=1, use_eclib = False, normalize = "L_ratio"):
        r"""
        Return the modular symbol associated to this elliptic curve,
        with given sign and base ring.  This is the map that sends `r/s`
        to a fixed multiple of the integral of `2 \pi i f(z) dz`
        from `\infty` to `r/s`, normalized so that all values of this map take
        values in `\QQ`.

        The normalization is such that for sign +1,
        the value at the cusp 0 is equal to the quotient of `L(E,1)`
        by the least positive period of `E` (unlike in ``L_ratio``
        of ``lseries()``, where the value is also divided by the
        number of connected components of `E(\RR)`). In particular the
        modular symbol depends on `E` and not only the isogeny class of `E`.

        INPUT:

        -  ``sign`` - 1 (default) or -1

        -  ``use_eclib`` - (default: False); if True the computation is
           done with John Cremona's implementation of modular
           symbols in ``eclib``

        -  ``normalize`` - (default: 'L_ratio'); either 'L_ratio',
           'period', or 'none';
           For 'L_ratio', the modular symbol tries to normalized correctly
           as explained above by comparing it to ``L_ratio`` for the
           curve and some small twists.
           The normalization 'period' is only available if
           ``use_eclib=False``. It uses the ``integral_period_map`` for modular
           symbols and is known to be equal to the above normalization
           up to the sign and a possible power of 2.
           For 'none', the modular symbol is almost certainly
           not correctly normalized, i.e. all values will be a
           fixed scalar multiple of what they should be.  But
           the initial computation of the modular symbol is
           much faster if ``use_eclib=False``, though evaluation of
           it after computing it won't be any faster.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: M=E.modular_symbol(); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: M(1/2)
            0
            sage: M(1/5)
            1

        ::

            sage: E=EllipticCurve('121b1')
            sage: M=E.modular_symbol()
            Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by -1, 2 or -2.
            sage: M(1/7)
            -1/2

        ::

            sage: E=EllipticCurve('11a1')
            sage: E.modular_symbol()(0)
            1/5
            sage: E=EllipticCurve('11a2')
            sage: E.modular_symbol()(0)
            1
            sage: E=EllipticCurve('11a3')
            sage: E.modular_symbol()(0)
            1/25

        ::

            sage: E=EllipticCurve('11a2')
            sage: E.modular_symbol(use_eclib=True, normalize='L_ratio')(0)
            1
            sage: E.modular_symbol(use_eclib=True, normalize='none')(0)
            2/5
            sage: E.modular_symbol(use_eclib=True, normalize='period')(0)
            Traceback (most recent call last):
            ...
            ValueError: no normalization 'period' known for modular symbols using John Cremona's eclib
            sage: E.modular_symbol(use_eclib=False, normalize='L_ratio')(0)
            1
            sage: E.modular_symbol(use_eclib=False, normalize='none')(0)
            1
            sage: E.modular_symbol(use_eclib=False, normalize='period')(0)
            1

        ::

            sage: E=EllipticCurve('11a3')
            sage: E.modular_symbol(use_eclib=True, normalize='L_ratio')(0)
            1/25
            sage: E.modular_symbol(use_eclib=True, normalize='none')(0)
            2/5
            sage: E.modular_symbol(use_eclib=True, normalize='period')(0)
            Traceback (most recent call last):
            ...
            ValueError: no normalization 'period' known for modular symbols using John Cremona's eclib
            sage: E.modular_symbol(use_eclib=False, normalize='L_ratio')(0)
            1/25
            sage: E.modular_symbol(use_eclib=False, normalize='none')(0)
            1
            sage: E.modular_symbol(use_eclib=False, normalize='period')(0)
            1/25

        """
        typ = (sign, normalize, use_eclib)
        try:
            return self.__modular_symbol[typ]
        except AttributeError:
            self.__modular_symbol = {}
        except KeyError:
            pass
        if use_eclib :
            M = ell_modular_symbols.ModularSymbolECLIB(self, sign, normalize=normalize)
        else :
            M = ell_modular_symbols.ModularSymbolSage(self, sign, normalize=normalize)
        self.__modular_symbol[typ] = M
        return M

    padic_lseries = padics.padic_lseries

    def newform(self):
        r"""
        Same as ``self.modular_form()``.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.newform()
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: E.newform() == E.modular_form()
            True
        """
        return self.modular_form()

    def q_eigenform(self, prec):
        r"""
        Synonym for ``self.q_expansion(prec)``.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.q_eigenform(10)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + O(q^10)
            sage: E.q_eigenform(10) == E.q_expansion(10)
            True
        """
        return self.q_expansion(prec)

    def analytic_rank(self, algorithm="pari", leading_coefficient=False):
        r"""
        Return an integer that is *probably* the analytic rank of this
        elliptic curve.  If leading_coefficient is ``True`` (only implemented
        for PARI), return a tuple `(rank, lead)` where `lead` is the value of
        the first non-zero derivative of the L-function of the elliptic
        curve.

        INPUT:

        - algorithm -

          - 'pari' (default) - use the PARI library function.

          - 'sympow' -use Watkins's program sympow

          - 'rubinstein' - use Rubinstein's L-function C++ program lcalc.

          - 'magma' - use MAGMA

          - 'all' - compute with all other free algorithms, check that
            the answers agree, and return the common answer.

        .. note::

           If the curve is loaded from the large Cremona database,
           then the modular degree is taken from the database.

        Of the three above, probably Rubinstein's is the most
        efficient (in some limited testing I've done).

        .. note::

           It is an open problem to *prove* that *any* particular
           elliptic curve has analytic rank `\geq 4`.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: E.analytic_rank(algorithm='pari')
            2
            sage: E.analytic_rank(algorithm='rubinstein')
              ***   Warning:...new stack size = ...
            2
            sage: E.analytic_rank(algorithm='sympow')
            2
            sage: E.analytic_rank(algorithm='magma')    # optional - magma
            2
            sage: E.analytic_rank(algorithm='all')
              ***   Warning:...new stack size = ...
            2

        With the optional parameter leading_coefficient set to ``True``, a
        tuple of both the analytic rank and the leading term of the
        L-series at `s = 1` is returned::

            sage: EllipticCurve([0,-1,1,-10,-20]).analytic_rank(leading_coefficient=True)
            (0, 0.25384186085591068...)
            sage: EllipticCurve([0,0,1,-1,0]).analytic_rank(leading_coefficient=True)
            (1, 0.30599977383405230...)
            sage: EllipticCurve([0,1,1,-2,0]).analytic_rank(leading_coefficient=True)
            (2, 1.518633000576853...)
            sage: EllipticCurve([0,0,1,-7,6]).analytic_rank(leading_coefficient=True)
            (3, 10.39109940071580...)
            sage: EllipticCurve([0,0,1,-7,36]).analytic_rank(leading_coefficient=True)
            (4, 196.170903794579...)

        TESTS:

        When the input is horrendous, some of the algorithms just bomb out with a RuntimeError::

            sage: EllipticCurve([1234567,89101112]).analytic_rank(algorithm='rubinstein')
            Traceback (most recent call last):
            ...
            RuntimeError: unable to compute analytic rank using rubinstein algorithm ('unable to convert x (= 6.19283e+19 and is too large) to an integer')
            sage: EllipticCurve([1234567,89101112]).analytic_rank(algorithm='sympow')
            Traceback (most recent call last):
            ...
            RuntimeError: failed to compute analytic rank
        """
        if algorithm == 'pari':
            rank_lead = self.pari_curve().ellanalyticrank()
            if leading_coefficient:
                return (rings.Integer(rank_lead[0]), rank_lead[1].python())
            else:
                return rings.Integer(self.pari_curve().ellanalyticrank()[0])
        elif algorithm == 'rubinstein':
            if leading_coefficient:
                raise NotImplementedError, "Cannot compute leading coefficient using rubinstein algorithm"
            try:
                from sage.lfunctions.lcalc import lcalc
                return lcalc.analytic_rank(L=self)
            except TypeError,msg:
                raise RuntimeError, "unable to compute analytic rank using rubinstein algorithm ('%s')"%msg
        elif algorithm == 'sympow':
            if leading_coefficient:
                raise NotImplementedError, "Cannot compute leading coefficient using sympow"
            from sage.lfunctions.sympow import sympow
            return sympow.analytic_rank(self)[0]
        elif algorithm == 'magma':
            if leading_coefficient:
                raise NotImplementedError, "Cannot compute leading coefficient using magma"
            from sage.interfaces.all import magma
            return rings.Integer(magma(self).AnalyticRank())
        elif algorithm == 'all':
            if leading_coefficient:
                S = set([self.analytic_rank('pari', True)])
            else:
                S = set([self.analytic_rank('pari'),
                    self.analytic_rank('rubinstein'), self.analytic_rank('sympow')])
            if len(S) != 1:
                raise RuntimeError, "Bug in analytic_rank; algorithms don't agree! (E=%s)"%E
            return list(S)[0]
        else:
            raise ValueError, "algorithm %s not defined"%algorithm


    def simon_two_descent(self, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Computes a lower bound for the rank of the Mordell-Weil group `E(Q)`,
        the rank of the 2-Selmer group, and a list of independent points on
        `E(Q)/2E(Q)`.

        INPUT:


        -  ``verbose`` - integer, 0,1,2,3; (default: 0), the
           verbosity level

        -  ``lim1`` - (default: 5) limite des points triviaux
           sur les quartiques

        -  ``lim3`` - (default: 50) limite des points sur les
           quartiques ELS

        -  ``limtriv`` - (default: 10) limite des points
           triviaux sur la courbe elliptique

        -  ``maxprob`` - (default: 20)

        -  ``limbigprime`` - (default: 30) to distinguish
           between small and large prime numbers. Use probabilistic tests for
           large primes. If 0, don't any probabilistic tests.


        OUTPUT:


        -  ``integer`` - lower bound on the rank of self

        -  ``integer`` - the 2-rank of the Selmer group

        -  ``list`` - list of independent points on the
           quotient `E(Q)/2E(Q)`.


        IMPLEMENTATION: Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/

        EXAMPLES: These computations use pseudo-random numbers, so we set
        the seed for reproducible testing.

        We compute the ranks of the curves of lowest known conductor up to
        rank `8`. Amazingly, each of these computations finishes
        almost instantly!

        ::

            sage: E = EllipticCurve('11a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (0, 0, [])
            sage: E = EllipticCurve('37a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (1, 1, [(0 : 0 : 1)])
            sage: E = EllipticCurve('389a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (2, 2, [(1 : 0 : 1), (-11/9 : -55/27 : 1)])
            sage: E = EllipticCurve('5077a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (3, 3, [(1 : -1 : 1), (2 : 0 : 1), (0 : 2 : 1)])

        In this example Simon's program does not find any points, though it
        does correctly compute the rank of the 2-Selmer group.

        ::

            sage: E = EllipticCurve([1, -1, 0, -751055859, -7922219731979])
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (1, 1, [])

        The rest of these entries were taken from Tom Womack's page
        http://tom.womack.net/maths/conductors.htm

        ::

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (4, 4, [(4 : 3 : 1), (5 : -2 : 1), (6 : -1 : 1), (8 : 7 : 1)])
            sage: E = EllipticCurve([0, 0, 1, -79, 342])
            sage: set_random_seed(0)
            sage: E.simon_two_descent()  # long time (9s on sage.math, 2011)
            (5, 5, [(5 : 8 : 1), (4 : 9 : 1), (3 : 11 : 1), (-1 : 20 : 1), (-6 : -25 : 1)])
            sage: E = EllipticCurve([1, 1, 0, -2582, 48720])
            sage: set_random_seed(0)
            sage: r, s, G = E.simon_two_descent(); r,s
            (6, 6)
            sage: E = EllipticCurve([0, 0, 0, -10012, 346900])
            sage: set_random_seed(0)
            sage: r, s, G = E.simon_two_descent(); r,s
            (7, 7)
            sage: E = EllipticCurve([0, 0, 1, -23737, 960366])
            sage: set_random_seed(0)
            sage: r, s, G = E.simon_two_descent(); r,s
            (8, 8)

        Example from trac 10832::

            sage: E = EllipticCurve([1,0,0,-6664,86543])
            sage: E.simon_two_descent()
            (2, 3, [(173 : 1943 : 1), (-73 : -394 : 1), (323/4 : 1891/8 : 1)])
            sage: E.rank()
            2
            sage: E.gens()
            [(-73 : -394 : 1), (323/4 : 1891/8 : 1)]

        Example from Trac #11372::

            sage: E = EllipticCurve([1, 1, 0, -23611790086, 1396491910863060])
            sage: E.simon_two_descent()
            (1, 2, [(88716 : -44358 : 1)])
            sage: E.rank()
            1
            sage: E.gens()
            [(4311692542083/48594841 : -13035144436525227/338754636611 : 1)]
        """
        verbose = int(verbose)
        t = simon_two_descent(self, verbose=verbose, lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime)
        if t=='fail':
            raise RuntimeError, 'Run-time error in simon_two_descent.'
        if verbose>0:
            print "simon_two_descent returns", t
        rank_low_bd = rings.Integer(t[0])
        two_selmer_rank = rings.Integer(t[1])
        gens_mod_two = [self(P) for P in t[2]]
        if rank_low_bd == two_selmer_rank - self.two_torsion_rank():
            if verbose>0:
                print "Rank determined successfully, saturating..."
            gens = [P for P in gens_mod_two if P.has_infinite_order()]
            gens = self.saturation(gens)[0]
            if len(gens) == rank_low_bd:
                self.__gens[True] = gens
                self.__gens[True].sort()
            self.__rank[True] = rank_low_bd
        return rank_low_bd, two_selmer_rank, gens_mod_two

    two_descent_simon = simon_two_descent

    def three_selmer_rank(self, algorithm='UseSUnits'):
        r"""
        Return the 3-selmer rank of this elliptic curve, computed using
        Magma.

        INPUT:


        -  ``algorithm`` - 'Heuristic' (which is usually much
           faster in large examples), 'FindCubeRoots', or 'UseSUnits'
           (default)


        OUTPUT: nonnegative integer

        EXAMPLES: A rank 0 curve::

            sage: EllipticCurve('11a').three_selmer_rank()       # optional - magma
            0

        A rank 0 curve with rational 3-isogeny but no 3-torsion

        ::

            sage: EllipticCurve('14a3').three_selmer_rank()      # optional - magma
            0

        A rank 0 curve with rational 3-torsion::

            sage: EllipticCurve('14a1').three_selmer_rank()      # optional - magma
            1

        A rank 1 curve with rational 3-isogeny::

            sage: EllipticCurve('91b').three_selmer_rank()       # optional - magma
            2

        A rank 0 curve with nontrivial 3-Sha. The Heuristic option makes
        this about twice as fast as without it.

        ::

            sage: EllipticCurve('681b').three_selmer_rank(algorithm='Heuristic')   # long time (10 seconds); optional - magma
            2
        """
        from sage.interfaces.all import magma
        E = magma(self)
        return Integer(E.ThreeSelmerGroup(MethodForFinalStep = magma('"%s"'%algorithm)).Ngens())

    def rank(self, use_database=False, verbose=False,
                   only_use_mwrank=True,
                   algorithm='mwrank_lib',
                   proof=None):
        """
        Return the rank of this elliptic curve, assuming no conjectures.

        If we fail to provably compute the rank, raises a RuntimeError
        exception.

        INPUT:


        -  ``use_database (bool)`` - (default: False), if
           True, try to look up the regulator in the Cremona database.

        -  ``verbose`` - (default: None), if specified changes
           the verbosity of mwrank computations. algorithm -

        -  ``- 'mwrank_shell'`` - call mwrank shell command

        -  ``- 'mwrank_lib'`` - call mwrank c library

        -  ``only_use_mwrank`` - (default: True) if False try
           using analytic rank methods first.

        -  ``proof`` - bool or None (default: None, see
           proof.elliptic_curve or sage.structure.proof). Note that results
           obtained from databases are considered proof = True


        OUTPUT:


        -  ``rank (int)`` - the rank of the elliptic curve.


        IMPLEMENTATION: Uses L-functions, mwrank, and databases.

        EXAMPLES::

            sage: EllipticCurve('11a').rank()
            0
            sage: EllipticCurve('37a').rank()
            1
            sage: EllipticCurve('389a').rank()
            2
            sage: EllipticCurve('5077a').rank()
            3
            sage: EllipticCurve([1, -1, 0, -79, 289]).rank()   # This will use the default proof behavior of True
            4
            sage: EllipticCurve([0, 0, 1, -79, 342]).rank(proof=False)
            5
            sage: EllipticCurve([0, 0, 1, -79, 342]).simon_two_descent()[0]  # long time (7s on sage.math, 2012)
            5

        Examples with denominators in defining equations::

            sage: E = EllipticCurve([0, 0, 0, 0, -675/4])
            sage: E.rank()
            0
            sage: E = EllipticCurve([0, 0, 1/2, 0, -1/5])
            sage: E.rank()
            1
            sage: E.minimal_model().rank()
            1

        A large example where mwrank doesn't determine the result with certainty::

            sage: EllipticCurve([1,0,0,0,37455]).rank(proof=False)
            0
            sage: EllipticCurve([1,0,0,0,37455]).rank(proof=True)
            Traceback (most recent call last):
            ...
            RuntimeError: Rank not provably correct.
        """
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)
        try:
            return self.__rank[proof]
        except KeyError:
            if proof is False and self.__rank.has_key(True):
                return self.__rank[True]
        if use_database:
            try:
                self.__rank[True] = self.database_curve().rank()
                return self.__rank[True]
            except (AttributeError, RuntimeError):
                pass
        if not only_use_mwrank:
            N = self.conductor()
            prec = int(4*float(sqrt(N))) + 10
            if self.root_number() == 1:
                L, err = self.lseries().at1(prec)
                if abs(L) > err + R(0.0001):  # definitely doesn't vanish
                    misc.verbose("rank 0 because L(E,1)=%s"%L)
                    self.__rank[proof] = 0
                    return self.__rank[proof]
            else:
                Lprime, err = self.lseries().deriv_at1(prec)
                if abs(Lprime) > err + R(0.0001):  # definitely doesn't vanish
                    misc.verbose("rank 1 because L'(E,1)=%s"%Lprime)
                    self.__rank[proof] = 1
                    return self.__rank[proof]

        if algorithm == 'mwrank_lib':
            misc.verbose("using mwrank lib")
            if self.is_integral(): E = self
            else: E = self.integral_model()
            C = E.mwrank_curve()
            C.set_verbose(verbose)
            r = C.rank()
            if C.certain():
                proof = True
            else:
                if proof:
                    print "Unable to compute the rank with certainty (lower bound=%s)."%C.rank()
                    print "This could be because Sha(E/Q)[2] is nontrivial."
                    print "Try calling something like two_descent(second_limit=13) on the"
                    print "curve then trying this command again.  You could also try rank"
                    print "with only_use_mwrank=False."
                    del E.__mwrank_curve
                    raise RuntimeError, 'Rank not provably correct.'
                else:
                    misc.verbose("Warning -- rank not proven correct", level=1)
            self.__rank[proof] = r
        elif algorithm == 'mwrank_shell':
            misc.verbose("using mwrank shell")
            X = self.mwrank()
            if 'determined unconditionally' not in X or 'only a lower bound of' in X:
                if proof:
                    X= "".join(X.split("\n")[-4:-2])
                    print X
                    raise RuntimeError, 'Rank not provably correct.'
                else:
                    misc.verbose("Warning -- rank not proven correct", level=1)

                s = "lower bound of"
                X = X[X.rfind(s)+len(s)+1:]
                r = Integer(X.split()[0])
            else:
                if proof is False:
                    proof = True #since we actually provably found the rank
                match = 'Rank ='
                i = X.find(match)
                if i == -1:
                    match = 'found points of rank'
                    i = X.find(match)
                    if i == -1:
                        raise RuntimeError, "%s\nbug -- tried to find 'Rank =' or 'found points of rank' in mwrank output but couldn't."%X
                j = i + X[i:].find('\n')
                r = Integer(X[i+len(match)+1:j])
            self.__rank[proof] = r

        return self.__rank[proof]

    def gens(self, verbose=False, rank1_search=10,
             algorithm='mwrank_lib',
             only_use_mwrank=True,
             proof = None,
             use_database = True,
             descent_second_limit=12,
             sat_bound = 1000):
        """
        Compute and return generators for the Mordell-Weil group E(Q)
        *modulo* torsion.

        .. warning::

           If the program fails to give a provably correct result, it
           prints a warning message, but does not raise an
           exception. Use the gens_certain command to find out if this
           warning message was printed.

        INPUT:


        -  ``verbose`` - (default: None), if specified changes
           the verbosity of mwrank computations.

        -  ``rank1_search`` - (default: 10), if the curve has
           analytic rank 1, try to find a generator by a direct search up to
           this logarithmic height. If this fails the usual mwrank procedure
           is called. algorithm -

        -  ``- 'mwrank_shell' (default)`` - call mwrank shell
           command

        -  ``- 'mwrank_lib'`` - call mwrank c library

        -  ``only_use_mwrank`` - bool (default True) if
           False, attempts to first use more naive, natively implemented
           methods.

        -  ``proof`` - bool or None (default None, see
           proof.elliptic_curve or sage.structure.proof).

        -  ``use_database`` - bool (default True) if True,
           attempts to find curve and gens in the (optional) database

        -  ``descent_second_limit`` - (default: 12)- used in 2-descent

        - ``sat_bound`` - (default: 1000) - bound on primes used in
           saturation.  If the computed bound on the index of the
           points found by two-descent in the Mordell-Weil group is
           greater than this, a warning message will be displayed.

        OUTPUT:


        -  ``generators`` - List of generators for the
           Mordell-Weil group modulo torsion.


        IMPLEMENTATION: Uses Cremona's mwrank C library.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: E.gens()                 # random output
            [(-1 : 1 : 1), (0 : 0 : 1)]

        A non-integral example::

            sage: E = EllipticCurve([-3/8,-2/3])
            sage: E.gens() # random (up to sign)
            [(10/9 : 29/54 : 1)]

        A non-minimal example::

            sage: E = EllipticCurve('389a1')
            sage: E1 = E.change_weierstrass_model([1/20,0,0,0]); E1
            Elliptic Curve defined by y^2 + 8000*y = x^3 + 400*x^2 - 320000*x over Rational Field
            sage: E1.gens() # random (if database not used)
            [(-400 : 8000 : 1), (0 : -8000 : 1)]
        """
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)

        # If the gens are already cached, return them:
        try:
            return list(self.__gens[proof])  # return copy so not changed
        except AttributeError:
            pass
        except KeyError:
            if proof is False and self.__gens.has_key(True):
                return self.__gens[True]

        # At this point, either self.__gens does not exist, or
        # self.__gens[False] exists but not self.__gens[True], and
        # proof is True

        # If the optional extended database is installed and an
        # isomorphic curve is in the database then its gens will be
        # known; if only the default database is installed, the rank
        # will be known but not the gens.

        if use_database:
            try:
                E = self.database_curve()
                iso = E.isomorphism_to(self)
                try:
                    self.__gens[True] = [iso(P) for P in E.__gens[True]]
                    return self.__gens[True]
                except (KeyError,AttributeError): # database curve does not have the gens
                    pass
            except (RuntimeError, KeyError):  # curve or gens not in database
                pass

        if self.conductor() > 10**7:
            only_use_mwrank = True

        if not only_use_mwrank:
            try:
                misc.verbose("Trying to compute rank.")
                r = self.rank(only_use_mwrank = False)
                misc.verbose("Got r = %s."%r)
                if r == 0:
                    misc.verbose("Rank = 0, so done.")
                    self.__gens[True] = []
                    return self.__gens[True]
                if r == 1 and rank1_search:
                    misc.verbose("Rank = 1, so using direct search.")
                    h = 6
                    while h <= rank1_search:
                        misc.verbose("Trying direct search up to height %s"%h)
                        G = self.point_search(h, verbose)
                        G = [P for P in G if P.order() == oo]
                        if len(G) > 0:
                            misc.verbose("Direct search succeeded.")
                            G, _, _ = self.saturation(G, verbose=verbose)
                            misc.verbose("Computed saturation.")
                            self.__gens[True] = G
                            return self.__gens[True]
                        h += 2
                    misc.verbose("Direct search FAILED.")
            except RuntimeError:
                pass
        # end if (not_use_mwrank)
        if algorithm == "mwrank_lib":
            misc.verbose("Calling mwrank C++ library.")
            if not self.is_integral():
                xterm = 1; yterm = 1
                ai = self.a_invariants()
                for a in ai:
                    if not a.is_integral():
                       for p, _ in a.denom().factor():
                          e  = min([(ai[i].valuation(p)/[1,2,3,4,6][i]) for i in range(5)]).floor()
                          ai = [ai[i]/p**(e*[1,2,3,4,6][i]) for i in range(5)]
                          xterm *= p**(2*e)
                          yterm *= p**(3*e)
                E = constructor.EllipticCurve(list(ai))
            else:
                E = self; xterm = 1; yterm = 1
            C = E.mwrank_curve(verbose)
            if not (verbose is None):
                C.set_verbose(verbose)
            C.two_descent(verbose=verbose, second_limit=descent_second_limit)
            C.saturate(bound=sat_bound)
            G = C.gens()
            if proof is True and C.certain() is False:
                del self.__mwrank_curve
                raise RuntimeError, "Unable to compute the rank, hence generators, with certainty (lower bound=%s, generators found=%s).  This could be because Sha(E/Q)[2] is nontrivial."%(C.rank(),G) + \
                      "\nTry increasing descent_second_limit then trying this command again."
            else:
                proof = C.certain()
            G = [[x*xterm,y*yterm,z] for x,y,z in G]
        else:
            # when gens() calls mwrank it passes the command-line
            # parameter "-p 100" which helps curves with large
            # coefficients and 2-torsion and is otherwise harmless.
            # This is pending a more intelligent handling of mwrank
            # options in gens() (which is nontrivial since gens() needs
            # to parse the output from mwrank and this is seriously
            # affected by what parameters the user passes!).
            # In fact it would be much better to avoid the mwrank console at
            # all for gens() and just use the library. This is in
            # progress (see trac #1949).
            X = self.mwrank('-p 100 -S '+str(sat_bound))
            misc.verbose("Calling mwrank shell.")
            if not 'The rank and full Mordell-Weil basis have been determined unconditionally' in X:
                msg = 'Generators not provably computed.'
                if proof:
                    raise RuntimeError, '%s\n%s'%(X,msg)
                else:
                    misc.verbose("Warning -- %s"%msg, level=1)
            elif proof is False:
                proof = True
            G = []
            i = X.find('Generator ')
            while i != -1:
                j = i + X[i:].find(';')
                k = i + X[i:].find('[')
                G.append(eval(X[k:j].replace(':',',')))
                X = X[j:]
                i = X.find('Generator ')
        ####
        self.__gens[proof] = [self.point(x, check=True) for x in G]
        self.__gens[proof].sort()
        self.__rank[proof] = len(self.__gens[proof])
        return self.__gens[proof]

    def gens_certain(self):
        """
        Return True if the generators have been proven correct.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.gens()                   # random (up to sign)
            [(0 : -1 : 1)]
            sage: E.gens_certain()
            True
        """
        return self.__gens.has_key(True)

    def ngens(self, proof = None):
        """
        Return the number of generators of this elliptic curve.

        .. note::

           See :meth:'.gens' for further documentation. The function
           :meth:`.ngens` calls :meth:`.gens` if not already done, but
           only with default parameters.  Better results may be
           obtained by calling ``mwrank()`` with carefully chosen
           parameters.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.ngens()
            1

        TO DO: This example should not cause a run-time error.

        ::

            sage: E=EllipticCurve([0,0,0,877,0])
            sage: # E.ngens()  ######## causes run-time error

        ::

            sage: print E.mwrank('-v0 -b12 -l')
            Curve [0,0,0,877,0] :   Rank = 1
            Generator 1 is [29604565304828237474403861024284371796799791624792913256602210:-256256267988926809388776834045513089648669153204356603464786949:490078023219787588959802933995928925096061616470779979261000]; height 95.980371987964
            Regulator = 95.980...
        """
        return len(self.gens(proof = proof))

    def regulator(self, use_database=True, proof=None, precision=None,
                        descent_second_limit=12, verbose=False):
        """
        Returns the regulator of this curve, which must be defined over Q.

        INPUT:


        -  ``use_database`` - bool (default: False), if True,
           try to look up the generators in the Cremona database.

        -  ``proof`` - bool or None (default: None, see
           proof.[tab] or sage.structure.proof). Note that results from
           databases are considered proof = True

        -  ``precision`` - int or None (default: None): the
           precision in bits of the result (default real precision if None)

        -  ``descent_second_limit`` - (default: 12)- used in 2-descent

        -  ``verbose`` - whether to print mwrank's verbose output

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.regulator()
            0.0511114082399688
            sage: EllipticCurve('11a').regulator()
            1.00000000000000
            sage: EllipticCurve('37a').regulator()
            0.0511114082399688
            sage: EllipticCurve('389a').regulator()
            0.152460177943144
            sage: EllipticCurve('5077a').regulator()
            0.41714355875838...
            sage: EllipticCurve([1, -1, 0, -79, 289]).regulator()
            1.50434488827528
            sage: EllipticCurve([0, 0, 1, -79, 342]).regulator(proof=False)  # long time (6s on sage.math, 2011)
            14.790527570131...
        """
        if precision is None:
            RR = rings.RealField()
            precision = RR.precision()
        else:
            RR = rings.RealField(precision)

        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)

        # We return a cached value if it exists and has sufficient precision:
        try:
            reg = self.__regulator[proof]
            if reg.parent().precision() >= precision:
                return RR(reg)
            else: # Found regulator value but precision is too low
                pass
        except KeyError:
            if proof is False and self.__regulator.has_key(True):
                reg = self.__regulator[True]
                if reg.parent().precision() >= precision:
                    return RR(reg)
                else: # Found regulator value but precision is too low
                    pass

        # Next we find the gens, taking them from the database if they
        # are there and use_database is True, else computing them:

        G = self.gens(proof=proof, use_database=use_database, descent_second_limit=descent_second_limit, verbose=verbose)

        # Finally compute the regulator of the generators found:

        self.__regulator[proof] = self.regulator_of_points(G,precision=precision)
        return self.__regulator[proof]

    def saturation(self, points, verbose=False, max_prime=0, odd_primes_only=False):
        """
        Given a list of rational points on E, compute the saturation in
        E(Q) of the subgroup they generate.

        INPUT:


        -  ``points (list)`` - list of points on E

        -  ``verbose (bool)`` - (default: False), if True, give
           verbose output

        -  ``max_prime (int)`` - (default: 0), saturation is
           performed for all primes up to max_prime. If max_prime==0,
           perform saturation at *all* primes, i.e., compute the true
           saturation.

        -  ``odd_primes_only (bool)`` - only do saturation at
           odd primes


        OUTPUT:


        -  ``saturation (list)`` - points that form a basis for
           the saturation

        -  ``index (int)`` - the index of the group generated
           by points in their saturation

        -  ``regulator (real with default precision)`` -
           regulator of saturated points.


        ALGORITHM: Uses Cremona's ``mwrank`` package. With ``max_prime=0``,
        we call ``mwrank`` with successively larger prime bounds until the full
        saturation is provably found. The results of saturation at the
        previous primes is stored in each case, so this should be
        reasonably fast.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: P=E(0,0)
            sage: Q=5*P; Q
            (1/4 : -5/8 : 1)
            sage: E.saturation([Q])
            ([(0 : 0 : 1)], 5, 0.0511114082399688)

        TESTS:

        See #10590.  This example would loop for ever at default precision::

            sage: E = EllipticCurve([1, 0, 1, -977842, -372252745])
            sage: P = E([-192128125858676194585718821667542660822323528626273/336995568430319276695106602174283479617040716649, 70208213492933395764907328787228427430477177498927549075405076353624188436/195630373799784831667835900062564586429333568841391304129067339731164107, 1])
            sage: P.height()
            113.302910926080
            sage: E.saturation([P])
            ([(-192128125858676194585718821667542660822323528626273/336995568430319276695106602174283479617040716649 : 70208213492933395764907328787228427430477177498927549075405076353624188436/195630373799784831667835900062564586429333568841391304129067339731164107 : 1)], 1, 113.302910926080)
            sage: E.saturation([2*P]) ## needs higher precision
            After 10 attempts at enlargement, giving up!
            ...
            ([(1755450733726721618440965414535034458701302721700399/970334851896750960577261378321772998240802013604 : -59636173615502879504846810677646864329901430096139563516090202443694810309127/955833935771565601591243078845907133814963790187832340692216425242529192 : 1)], 2, 113.302910926080)

        See #10840.  This causes eclib to crash since the curve is
        non-minimal at 2::

            sage: E = EllipticCurve([0,0,0,-13711473216,0])
            sage: P = E([-19992,16313472])
            sage: Q = E([-24108,-17791704])
            sage: R = E([-97104,-20391840])
            sage: S = E([-113288,-9969344])
            sage: E.saturation([P,Q,R,S])
            ([(-19992 : 16313472 : 1), (-24108 : -17791704 : 1), (-97104 : -20391840 : 1), (-113288 : -9969344 : 1)], 1, 172.792031341679)

        """
        if not isinstance(points, list):
            raise TypeError, "points (=%s) must be a list."%points
        if len(points) == 0:
            return [], None, R(1)

        v = []
        for P in points:
            if not isinstance(P, ell_point.EllipticCurvePoint_field):
                P = self(P)
            elif P.curve() != self:
                raise ArithmeticError, "point (=%s) must be %s."%(P,self)

        minimal = True
        if not self.is_minimal():
            minimal = False
            Emin = self.minimal_model()
            phi = self.isomorphism_to(Emin)
            points = [phi(P) for P in points]
        else:
            Emin = self

        for P in points:
            x, y = P.xy()
            d = x.denominator().lcm(y.denominator())
            v.append((x*d, y*d, d))

        c = Emin.mwrank_curve()
        mw = mwrank.mwrank_MordellWeil(c, verbose)
        mw.process(v)
        repeat_until_saturated = False
        if max_prime == 0:
            repeat_until_saturated = True
            max_prime = 9973
        from sage.libs.all import mwrank_get_precision, mwrank_set_precision
        prec0 = mwrank_get_precision()
        prec = 100
        if prec0<prec:
            mwrank_set_precision(prec)
        else:
            prec = prec0
        while True:
            ok, index, unsat = mw.saturate(max_prime=max_prime, odd_primes_only = odd_primes_only)
            reg = mw.regulator()
            if ok or not repeat_until_saturated: break
            max_prime = arith.next_prime(max_prime + 1000)
            prec += 50
            #print "Increasing precision to ",prec," and max_prime to ",max_prime
            mwrank_set_precision(prec)
        if prec!=prec0: mwrank_set_precision(prec0)
        #print "precision was originally ",prec0," and is now ",mwrank_get_precision()
        sat = mw.points()
        sat = [Emin(P) for P in sat]
        if not minimal:
            phi_inv = ~phi
            sat = [phi_inv(P) for P in sat]
        reg = self.regulator_of_points(sat)
        return sat, index, R(reg)


    def CPS_height_bound(self):
        r"""
        Return the Cremona-Prickett-Siksek height bound. This is a
        floating point number B such that if P is a rational point on
        the curve, then `h(P) \le \hat{h}(P) + B`, where `h(P)` is
        the naive logarithmic height of `P` and `\hat{h}(P)` is the
        canonical height.

        SEE ALSO: silverman_height_bound for a bound that also works for
        points over number fields.

        EXAMPLES::

            sage: E = EllipticCurve("11a")
            sage: E.CPS_height_bound()
            2.8774743273580445
            sage: E = EllipticCurve("5077a")
            sage: E.CPS_height_bound()
            0.0
            sage: E = EllipticCurve([1,2,3,4,1])
            sage: E.CPS_height_bound()
            Traceback (most recent call last):
            ...
            RuntimeError: curve must be minimal.
            sage: F = E.quadratic_twist(-19)
            sage: F
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 + 1376*x - 130 over Rational Field
            sage: F.CPS_height_bound()
            0.6555158376972852

        IMPLEMENTATION:
            Call the corresponding mwrank C++ library function.  Note that
            the formula in the [CPS] paper is given for number fields.  It's
            only the implementation in Sage that restricts to the rational
            field.
        """
        if not self.is_minimal():
            raise RuntimeError, "curve must be minimal."
        return self.mwrank_curve().CPS_height_bound()


    def silverman_height_bound(self, algorithm='default'):
        r"""
        Return the Silverman height bound.  This is a positive real
        (floating point) number B such that for all points `P` on the
        curve over any number field, `|h(P) - \hat{h}(P)| \leq B`,
        where `h(P)` is the naive logarithmic height of `P` and
        `\hat{h}(P)` is the canonical height.

        INPUT:

            - ``algorithm`` --

                 - 'default' (default) -- compute using a Python
                   implementation in Sage

                 - 'mwrank' -- use a C++ implementation in the mwrank
                   library

        NOTES:

           - The CPS_height_bound is often better (i.e. smaller) than
             the Silverman bound, but it only applies for points over
             the base field, whereas the Silverman bound works over
             all number fields.

           - The Silverman bound is also fairly straightforward to
             compute over number fields, but isn't implemented here.

           - Silverman's paper is 'The Difference Between the Weil
             Height and the Canonical Height on Elliptic Curves',
             Math. Comp., Volume 55, Number 192, pages 723-743.  We
             use a correction by Bremner with 0.973 replaced by 0.961,
             as explained in the source code to mwrank (htconst.cc).

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.silverman_height_bound()
            4.825400758180918
            sage: E.silverman_height_bound(algorithm='mwrank')
            4.825400758180918
            sage: E.CPS_height_bound()
            0.16397076103046915
        """
        if algorithm == 'default':
            Delta   = self.discriminant()
            j       = self.j_invariant()
            b2      = self.b2()
            twostar = 2 if b2 else 1
            from math import log
            def h(x):
                return log(max(abs(x.numerator()), abs(x.denominator())))
            def h_oo(x):
                return log(max(abs(x),1))
            mu    = h(Delta)/12 + h_oo(j)/12 + h_oo(b2/12)/2 + log(twostar)/2
            lower = 2*(-h(j)/24 - mu - 0.961)
            upper = 2*(mu + 1.07)
            return max(abs(lower), abs(upper))
        elif algorithm == 'mwrank':
            return self.mwrank_curve().silverman_bound()
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm



    def point_search(self, height_limit, verbose=False, rank_bound=None):
        """
        Search for points on a curve up to an input bound on the naive
        logarithmic height.

        INPUT:


        -  ``height_limit (float)`` - bound on naive height

        -  ``verbose (bool)`` - (default: False)

           If True, report on each point as found together with linear
           relations between the points found and the saturation process.

           If False, just return the result.

        -  ``rank_bound (bool)`` - (default: None)

           If provided, stop searching for points once we find this many
           independent nontorsion points.

        OUTPUT: points (list) - list of independent points which generate
        the subgroup of the Mordell-Weil group generated by the points
        found and then saturated.

        .. warning::

           height_limit is logarithmic, so increasing by 1 will cause
           the running time to increase by a factor of approximately
           4.5 (=exp(1.5)).

        IMPLEMENTATION: Uses Michael Stoll's ratpoints library.

        EXAMPLES::

            sage: E=EllipticCurve('389a1')
            sage: E.point_search(5, verbose=False)
            [(-1 : 1 : 1), (-3/4 : 7/8 : 1)]

        Increasing the height_limit takes longer, but finds no more
        points::

            sage: E.point_search(10, verbose=False)
            [(-1 : 1 : 1), (-3/4 : 7/8 : 1)]

        In fact this curve has rank 2 so no more than 2 points will ever be
        output, but we are not using this fact.

        ::

            sage: E.saturation(_)
            ([(-1 : 1 : 1), (-3/4 : 7/8 : 1)], 1, 0.152460177943144)

        What this shows is that if the rank is 2 then the points listed do
        generate the Mordell-Weil group (mod torsion). Finally,

        ::

            sage: E.rank()
            2

        If we only need one independent generator::

            sage: E.point_search(5, verbose=False, rank_bound=1)
            [(-2 : 0 : 1)]

        """
        from sage.libs.ratpoints import ratpoints
        from sage.functions.all import exp
        from sage.rings.arith import GCD
        H = exp(float(height_limit)) # max(|p|,|q|) <= H, if x = p/q coprime
        coeffs = [16*self.b6(), 8*self.b4(), self.b2(), 1]
        points = []
        a1 = self.a1()
        a3 = self.a3()
        new_H = H*2 # since we change the x-coord by 2 below
        for X,Y,Z in ratpoints(coeffs, new_H, verbose):
            if Z == 0: continue
            z = 2*Z
            x = X/2
            y = (Y/z - a1*x - a3*z)/2
            d = GCD((x,y,z))
            x = x/d
            if max(x.numerator().abs(), x.denominator().abs()) <= H:
                y = y/d
                z = z/d
                points.append(self((x,y,z)))
                if rank_bound is not None:
                    points = self.saturation(points, verbose=verbose)[0]
                    if len(points) == rank_bound:
                        break
        if rank_bound is None:
            points = self.saturation(points, verbose=verbose)[0]
        return points


    def selmer_rank(self):
        """
        The rank of the 2-Selmer group of the curve.

        EXAMPLE: The following is the curve 960D1, which has rank 0, but
        Sha of order 4.

        ::

            sage: E = EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.selmer_rank()
            3

        Here the Selmer rank is equal to the 2-torsion rank (=1) plus
        the 2-rank of Sha (=2), and the rank itself is zero::

            sage: E.rank()
            0

        In contrast, for the curve 571A, also with rank 0 and Sha of
        order 4, we get a worse bound::

            sage: E = EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.selmer_rank()
            2
            sage: E.rank_bound()
            2

        To establish that the rank is in fact 0 in this case, we would
        need to carry out a higher descent::

            sage: E.three_selmer_rank() # optional: magma
            0

        Or use the L-function to compute the analytic rank::

            sage: E.rank(only_use_mwrank=False)
            0

        """
        try:
            return self.__selmer_rank
        except AttributeError:
            C = self.mwrank_curve()
            self.__selmer_rank = C.selmer_rank()
            return self.__selmer_rank


    def rank_bound(self):
        """
        Upper bound on the rank of the curve, computed using
        2-descent.  In many cases, this is the actual rank of the
        curve.  If the curve has no 2-torsion it is the same as the
        2-selmer rank.

        EXAMPLE: The following is the curve 960D1, which has rank 0, but
        Sha of order 4.

        ::

            sage: E = EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.rank_bound()
            0

        It gives 0 instead of 2, because it knows Sha is nontrivial. In
        contrast, for the curve 571A, also with rank 0 and Sha of order 4,
        we get a worse bound::

            sage: E = EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.rank_bound()
            2
            sage: E.rank(only_use_mwrank=False)   # uses L-function
            0

        """
        try:
            return self.__rank_bound
        except AttributeError:
            C = self.mwrank_curve()
            self.__rank_bound = C.rank_bound()
            return self.__rank_bound


    def an(self, n):
        """
        The n-th Fourier coefficient of the modular form corresponding to
        this elliptic curve, where n is a positive integer.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: [E.an(n) for n in range(20) if n>0]
            [1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0]
        """
        return Integer(self.pari_mincurve().ellak(n))

    def ap(self, p):
        """
        The p-th Fourier coefficient of the modular form corresponding to
        this elliptic curve, where p is prime.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: [E.ap(p) for p in prime_range(50)]
            [-2, -3, -2, -1, -5, -2, 0, 0, 2, 6, -4, -1, -9, 2, -9]
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        return Integer(self.pari_mincurve().ellap(p))

    def quadratic_twist(self, D):
        """
       Return the quadratic twist of this elliptic curve by D.

       D must be a nonzero rational number.

       .. note::

          This function overrides the generic ``quadratic_twist()``
          function for elliptic curves, returning a minimal model.

       EXAMPLES::

           sage: E=EllipticCurve('37a1')
           sage: E2=E.quadratic_twist(2); E2
           Elliptic Curve defined by y^2  = x^3 - 4*x + 2 over Rational Field
           sage: E2.conductor()
           2368
           sage: E2.quadratic_twist(2) == E
           True
       """
        return EllipticCurve_number_field.quadratic_twist(self, D).minimal_model()

    def minimal_model(self):
        r"""
        Return the unique minimal Weierstrass equation for this elliptic
        curve. This is the model with minimal discriminant and
        `a_1,a_2,a_3 \in \{0,\pm 1\}`.

        EXAMPLES::

            sage: E=EllipticCurve([10,100,1000,10000,1000000])
            sage: E.minimal_model()
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + x + 1 over Rational Field
        """
        try:
            return self.__minimal_model
        except AttributeError:
            F = self.pari_mincurve()
            self.__minimal_model = EllipticCurve_rational_field([Q(F[i]) for i in range(5)])
            return self.__minimal_model

    def is_minimal(self):
        r"""
        Return True iff this elliptic curve is a reduced minimal model.

        The unique minimal Weierstrass equation for this elliptic curve.
        This is the model with minimal discriminant and
        `a_1,a_2,a_3 \in \{0,\pm 1\}`.

        TO DO: This is not very efficient since it just computes the
        minimal model and compares. A better implementation using the Kraus
        conditions would be preferable.

        EXAMPLES::

            sage: E=EllipticCurve([10,100,1000,10000,1000000])
            sage: E.is_minimal()
            False
            sage: E=E.minimal_model()
            sage: E.is_minimal()
            True
        """
        return self.ainvs() == self.minimal_model().ainvs()

    def is_p_minimal(self, p):
        """
        Tests if curve is p-minimal at a given prime p.

        INPUT: p - a primeOUTPUT: True - if curve is p-minimal


        -  ``False`` - if curve isn't p-minimal


        EXAMPLES::

            sage: E = EllipticCurve('441a2')
            sage: E.is_p_minimal(7)
            True

        ::

            sage: E = EllipticCurve([0,0,0,0,(2*5*11)**10])
            sage: [E.is_p_minimal(p) for p in prime_range(2,24)]
            [False, True, False, True, False, True, True, True, True]
        """
        if not p.is_prime():
            raise ValueError,"p must be prime"
        if not self.is_p_integral(p):
            return False
        if p > 3:
            return ((self.discriminant().valuation(p) < 12) or (self.c4().valuation(p) < 4))
        # else p = 2,3
        Emin = self.minimal_model()
        return self.discriminant().valuation(p) == Emin.discriminant().valuation(p)

#    Duplicate!
#
#    def is_integral(self):
#        for n in self.ainvs():
#            if n.denominator() != 1:
#                return False
#        return True

    def kodaira_type(self, p):
        """
        Local Kodaira type of the elliptic curve at `p`.

        INPUT:

        -   p, an integral prime


        OUTPUT:


        - the Kodaira type of this elliptic curve at p,
          as a KodairaSymbol.


        EXAMPLES::

            sage: E = EllipticCurve('124a')
            sage: E.kodaira_type(2)
            IV
        """
        return self.local_data(p).kodaira_symbol()

    kodaira_symbol = kodaira_type

    def kodaira_type_old(self, p):
        """
        Local Kodaira type of the elliptic curve at `p`.

        INPUT:


        -   p, an integral prime


        OUTPUT:

        - the kodaira type of this elliptic curve at p,
          as a KodairaSymbol.

        EXAMPLES::

            sage: E = EllipticCurve('124a')
            sage: E.kodaira_type_old(2)
            IV
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        try:
            self.__kodaira_type
        except AttributeError:
            self.__kodaira_type = {}
            self.__tamagawa_number = {}
        if not self.__kodaira_type.has_key(p):
            v = self.pari_mincurve().elllocalred(p)
            from kodaira_symbol import KodairaSymbol
            self.__kodaira_type[p] = KodairaSymbol(v[1])
            self.__tamagawa_number[p] = Integer(v[3])
        return self.__kodaira_type[p]

    def tamagawa_number(self, p):
        r"""
        The Tamagawa number of the elliptic curve at `p`.

        This is the order of the component group
        `E(\QQ_p)/E^0(\QQ_p)`.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: E.tamagawa_number(11)
            5
            sage: E = EllipticCurve('37b')
            sage: E.tamagawa_number(37)
            3
        """
        return self.local_data(p).tamagawa_number()

    def tamagawa_number_old(self, p):
        r"""
        The Tamagawa number of the elliptic curve at `p`.

        This is the order of the component group
        `E(\QQ_p)/E^0(\QQ_p)`.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: E.tamagawa_number_old(11)
            5
            sage: E = EllipticCurve('37b')
            sage: E.tamagawa_number_old(37)
            3
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        try:
            return self.__tamagawa_number[p]
        except (AttributeError, KeyError):
            self.kodaira_type_old(p)
            return self.__tamagawa_number[p]

    def tamagawa_exponent(self, p):
        """
        The Tamagawa index of the elliptic curve at `p`.

        This is the index of the component group
        `E(\QQ_p)/E^0(\QQ_p)`. It equals the
        Tamagawa number (as the component group is cyclic) except for types
        `I_m^*` (`m` even) when the group can be
        `C_2 \times C_2`.

        EXAMPLES::

            sage: E = EllipticCurve('816a1')
            sage: E.tamagawa_number(2)
            4
            sage: E.tamagawa_exponent(2)
            2
            sage: E.kodaira_symbol(2)
            I2*

        ::

            sage: E = EllipticCurve('200c4')
            sage: E.kodaira_symbol(5)
            I4*
            sage: E.tamagawa_number(5)
            4
            sage: E.tamagawa_exponent(5)
            2

        See #4715::

            sage: E=EllipticCurve('117a3')
            sage: E.tamagawa_exponent(13)
            4
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        cp = self.tamagawa_number(p)
        if not cp==4:
            return cp
        ks = self.kodaira_type(p)
        if ks._roman==1 and ks._n%2==0 and ks._starred:
            return 2
        return 4

    def tamagawa_product(self):
        """
        Returns the product of the Tamagawa numbers.

        EXAMPLES::

            sage: E = EllipticCurve('54a')
            sage: E.tamagawa_product ()
            3
        """
        try:
            return self.__tamagawa_product
        except AttributeError:
            self.__tamagawa_product = Integer(self.pari_mincurve().ellglobalred()[2].python())
            return self.__tamagawa_product

    def real_components(self):
        """
        Returns 1 if there is 1 real component and 2 if there are 2.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.real_components ()
            2
            sage: E = EllipticCurve('37b')
            sage: E.real_components ()
            2
            sage: E = EllipticCurve('11a')
            sage: E.real_components ()
            1
        """
        invs = self.short_weierstrass_model().ainvs()
        x = rings.polygen(self.base_ring())
        f = x**3 + invs[3]*x + invs[4]
        if f.discriminant() > 0:
            return 2
        else:
            return 1

    def has_good_reduction_outside_S(self,S=[]):
        r"""
        Tests if this elliptic curve has good reduction outside `S`.

        INPUT:

            -  S - list of primes (default: empty list).

        .. note::

            Primality of elements of S is not checked, and the output
            is undefined if S is not a list or contains non-primes.

            This only tests the given model, so should only be applied to
            minimal models.

        EXAMPLES::

            sage: EllipticCurve('11a1').has_good_reduction_outside_S([11])
            True
            sage: EllipticCurve('11a1').has_good_reduction_outside_S([2])
            False
            sage: EllipticCurve('2310a1').has_good_reduction_outside_S([2,3,5,7])
            False
            sage: EllipticCurve('2310a1').has_good_reduction_outside_S([2,3,5,7,11])
            True
        """
        return self.discriminant().is_S_unit(S)

    def period_lattice(self, embedding=None):
        r"""
        Returns the period lattice of the elliptic curve with respect to
        the differential `dx/(2y + a_1x + a_3)`.

        INPUT:

        -  ``embedding`` - ignored (for compatibility with the
           period_lattice function for elliptic_curve_number_field)

        OUTPUT:

        (period lattice) The PeriodLattice_ell object associated to
        this elliptic curve (with respect to the natural embedding of
        `\QQ` into `\RR`).

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        try:
            return self._period_lattice
        except AttributeError:
            from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            self._period_lattice = PeriodLattice_ell(self)
            return self._period_lattice

    def elliptic_exponential(self, z, embedding=None):
        r"""
        Computes the elliptic exponential of a complex number with respect to the elliptic curve.

        INPUT:

        - ``z`` (complex) -- a complex number

        -  ``embedding`` - ignored (for compatibility with the
           period_lattice function for elliptic_curve_number_field)

        OUTPUT:

        The image of `z` modulo `L` under the Weierstrass parametrization
        `\CC/L \to E(\CC)`.

        .. note::

           The precision is that of the input ``z``, or the default
           precision of 53 bits if ``z`` is exact.

        EXAMPLES::

            sage: E = EllipticCurve([1,1,1,-8,6])
            sage: P = E([1,-2])
            sage: z = P.elliptic_logarithm() # default precision is 100 here
            sage: E.elliptic_exponential(z)
            (1.0000000000000000000000000000 : -2.0000000000000000000000000000 : 1.0000000000000000000000000000)
            sage: z = E([1,-2]).elliptic_logarithm(precision=201)
            sage: E.elliptic_exponential(z)
            (1.00000000000000000000000000000000000000000000000000000000000 : -2.00000000000000000000000000000000000000000000000000000000000 : 1.00000000000000000000000000000000000000000000000000000000000)

        ::

            sage: E = EllipticCurve('389a')
            sage: Q = E([3,5])
            sage: E.elliptic_exponential(Q.elliptic_logarithm())
            (3.0000000000000000000000000000 : 5.0000000000000000000000000000 : 1.0000000000000000000000000000)
            sage: P = E([-1,1])
            sage: P.elliptic_logarithm()
            0.47934825019021931612953301006 + 0.98586885077582410221120384908*I
            sage: E.elliptic_exponential(P.elliptic_logarithm())
            (-1.0000000000000000000000000000 : 1.0000000000000000000000000000 : 1.0000000000000000000000000000)


        Some torsion examples::

            sage: w1,w2 = E.period_lattice().basis()
            sage: E.two_division_polynomial().roots(CC,multiplicities=False)
            [-2.0403022002854..., 0.13540924022175..., 0.90489296006371...]
            sage: [E.elliptic_exponential((a*w1+b*w2)/2)[0] for a,b in [(0,1),(1,1),(1,0)]]
            [-2.0403022002854..., 0.13540924022175..., 0.90489296006371...]

            sage: E.division_polynomial(3).roots(CC,multiplicities=False)
            [-2.88288879135...,
            1.39292799513...,
            0.078313731444316... - 0.492840991709...*I,
            0.078313731444316... + 0.492840991709...*I]
            sage: [E.elliptic_exponential((a*w1+b*w2)/3)[0] for a,b in [(0,1),(1,0),(1,1),(2,1)]]
            [-2.8828887913533..., 1.39292799513138,
            0.0783137314443... - 0.492840991709...*I,
            0.0783137314443... + 0.492840991709...*I]

        Observe that this is a group homomorphism (modulo rounding error)::

            sage: z = CC.random_element()
            sage: 2 * E.elliptic_exponential(z)
            (-1.52184235874404 - 0.0581413944316544*I : 0.948655866506124 - 0.0381469928565030*I : 1.00000000000000)
            sage: E.elliptic_exponential(2 * z)
            (-1.52184235874404 - 0.0581413944316562*I : 0.948655866506128 - 0.0381469928565034*I : 1.00000000000000)
        """
        return self.period_lattice().elliptic_exponential(z)

    def lseries(self):
        """
        Returns the L-series of this elliptic curve.

        Further documentation is available for the functions which apply to
        the L-series.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.lseries()
            Complex L-series of the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        try:
            return self.__lseries
        except AttributeError:
            from lseries_ell import Lseries_ell
            self.__lseries = Lseries_ell(self)
            return self.__lseries

    def Lambda(self, s, prec):
        r"""
        Returns the value of the Lambda-series of the elliptic curve E at
        s, where s can be any complex number.

        IMPLEMENTATION: Fairly *slow* computation using the definitions
        and implemented in Python.

        Uses prec terms of the power series.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: E.Lambda(1.4+0.5*I, 50)
            -0.354172680517... + 0.874518681720...*I
        """
        from sage.all import pi

        s = C(s)
        N = self.conductor()
        pi = R(pi)
        a = self.anlist(prec)
        eps = self.root_number()
        sqrtN = float(N.sqrt())
        def _F(n, t):
            return gamma_inc(t+1, 2*pi*n/sqrtN) * C(sqrtN/(2*pi*n))**(t+1)
        return sum([a[n]*(_F(n,s-1) + eps*_F(n,1-s)) for n in xrange(1,prec+1)])

    def is_local_integral_model(self,*p):
        r"""
        Tests if self is integral at the prime `p`, or at all the
        primes if `p` is a list or tuple of primes

        EXAMPLES::

            sage: E=EllipticCurve([1/2,1/5,1/5,1/5,1/5])
            sage: [E.is_local_integral_model(p) for p in (2,3,5)]
            [False, True, False]
            sage: E.is_local_integral_model(2,3,5)
            False
            sage: Eint2=E.local_integral_model(2)
            sage: Eint2.is_local_integral_model(2)
            True
        """
        if len(p)==1: p=p[0]
        if isinstance(p,(tuple,list)):
            return misc.forall(p, lambda x : self.is_local_integral_model(x))[0]
        assert p.is_prime(), "p must be prime in is_local_integral_model()"
        return misc.forall(self.ainvs(), lambda x : x.valuation(p) >= 0)[0]

    def local_integral_model(self,p):
        r"""
        Return a model of self which is integral at the prime `p`.

        EXAMPLES::

            sage: E=EllipticCurve([0, 0, 1/216, -7/1296, 1/7776])
            sage: E.local_integral_model(2)
            Elliptic Curve defined by y^2 + 1/27*y = x^3 - 7/81*x + 2/243 over Rational Field
            sage: E.local_integral_model(3)
            Elliptic Curve defined by y^2 + 1/8*y = x^3 - 7/16*x + 3/32 over Rational Field
            sage: E.local_integral_model(2).local_integral_model(3) == EllipticCurve('5077a1')
            True
        """
        assert p.is_prime(), "p must be prime in local_integral_model()"
        ai = self.a_invariants()
        e  = min([(ai[i].valuation(p)/[1,2,3,4,6][i]) for i in range(5)]).floor()
        return constructor.EllipticCurve([ai[i]/p**(e*[1,2,3,4,6][i]) for i in range(5)])

    def is_global_integral_model(self):
        r"""
        Return true iff self is integral at all primes.

        EXAMPLES::

            sage: E=EllipticCurve([1/2,1/5,1/5,1/5,1/5])
            sage: E.is_global_integral_model()
            False
            sage: Emin=E.global_integral_model()
            sage: Emin.is_global_integral_model()
            True
        """
        return self.is_integral()

    def global_integral_model(self):
        r"""
        Return a model of self which is integral at all primes.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1/216, -7/1296, 1/7776])
            sage: F = E.global_integral_model(); F
            Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
            sage: F == EllipticCurve('5077a1')
            True
        """
        ai = self.a_invariants()
        for a in ai:
            if not a.is_integral():
               for p, _ in a.denom().factor():
                  e  = min([(ai[i].valuation(p)/[1,2,3,4,6][i]) for i in range(5)]).floor()
                  ai = [ai[i]/p**(e*[1,2,3,4,6][i]) for i in range(5)]
        for z in ai:
            assert z.denominator() == 1, "bug in global_integral_model: %s" % ai
        return constructor.EllipticCurve(list(ai))

    integral_model = global_integral_model

    def integral_short_weierstrass_model(self):
        r"""
        Return a model of the form `y^2 = x^3 + ax + b` for this
        curve with `a,b\in\ZZ`.

        EXAMPLES::

            sage: E = EllipticCurve('17a1')
            sage: E.integral_short_weierstrass_model()
            Elliptic Curve defined by y^2  = x^3 - 11*x - 890 over Rational Field
        """
        F = self.minimal_model().short_weierstrass_model()
        _,_,_,A,B = F.ainvs()
        for p in [2,3]:
            e=min(A.valuation(p)/4,B.valuation(p)/6).floor()
            A /= Integer(p**(4*e))
            B /= Integer(p**(6*e))
        return constructor.EllipticCurve([A,B])

    # deprecated function replaced by integral_short_weierstrass_model, see trac 3974.
    def integral_weierstrass_model(self):
        r"""
        Return a model of the form `y^2 = x^3 + ax + b` for this
        curve with `a,b\in\ZZ`.

        Note that this function is deprecated, and that you should use
        integral_short_weierstrass_model instead as this will be
        disappearing in the near future.

        EXAMPLES::

            sage: E = EllipticCurve('17a1')
            sage: E.integral_weierstrass_model() #random
            doctest:1: DeprecationWarning: integral_weierstrass_model is deprecated, use integral_short_weierstrass_model instead!
            Elliptic Curve defined by y^2  = x^3 - 11*x - 890 over Rational Field
        """
        from sage.misc.superseded import deprecation
        deprecation(3974, "integral_weierstrass_model is deprecated, use integral_short_weierstrass_model instead!")
        return self.integral_short_weierstrass_model()


    def _generalized_congmod_numbers(self, M, invariant="both"):
        """
        Internal method to compute the generalized modular degree and congruence number
        at level `MN`, where `N` is the conductor of `E`.
        Values obtained are cached.

        This function is called by self.modular_degree() and self.congruence_number() when
        `M>1`. Since so much of the computation of the two values is shared, this method
        by default computes and caches both.

        INPUT:

        - ``M`` - Non-negative integer; this function is only ever called on M>1, although
          the algorithm works fine for the case `M==1`

        - ``invariant`` - String; default "both". Options are:

          - "both" - Both modular degree and congruence number at level `MN` are computed

          - "moddeg" - Only modular degree is computed

          - "congnum" - Only congruence number is computed

        OUTPUT:

        - A dictionary containing either the modular degree (a positive integer) at index "moddeg",
          or the congruence number (a positive integer) at index "congnum", or both.

        As far as we know there is no other implementation for this algorithm, so as yet
        there is nothing to check the below examples against.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: for M in range(2,8): print(M,E.modular_degree(M=M),E.congruence_number(M=M)) # long time (22s on 2009 MBP)
            (2, 5, 20)
            (3, 7, 28)
            (4, 50, 400)
            (5, 32, 128)
            (6, 1225, 19600)
            (7, 63, 252)
        """
        # Check invariant specification before we get going
        if invariant not in ["moddeg", "congnum", "both"]:
            raise ValueError("Invalid invariant specification")

        # Cuspidal space at level MN
        N = self.conductor()
        S = ModularSymbols(N*M,sign=1).cuspidal_subspace()

        # Cut out the subspace by hitting it with T_p for enough p
        A = S
        d = self.dimension()*arith.sigma(M,0)
        p = 2
        while A.dimension() > d:
            while N*M % p == 0:
                p = next_prime(p)
            Tp = A.hecke_operator(p)
            A = (Tp - self.ap(p)).kernel()
            p = next_prime(p)
        B = A.complement().cuspidal_submodule()

        L = {}
        if invariant in ["moddeg", "both"]:
            V = A.integral_structure()
            W = B.integral_structure()
            moddeg  = (V + W).index_in(S.integral_structure())
            L["moddeg"] = moddeg
            self.__generalized_modular_degree[M] = moddeg

        if invariant in ["congnum", "both"]:
            congnum = A.congruence_number(B)
            L["congnum"] = congnum
            self.__generalized_congruence_number[M] = congnum

        return L


    def modular_degree(self, algorithm='sympow', M=1):
        r"""
        Return the modular degree at level `MN` of this elliptic curve. The case
        `M==1` corresponds to the classical definition of modular degree.

        When `M>1`, the function returns the degree of the map from `X_0(MN) \to A`, where
        A is the abelian variety generated by embeddings of `E` into `J_0(MN)`.

        The result is cached. Subsequent calls, even with a different
        algorithm, just returned the cached result. The algorithm argument is ignored
        when `M>1`.

        INPUT:

        -  ``algorithm`` - string:

        -  ``'sympow'`` - (default) use Mark Watkin's (newer) C
           program sympow

        -  ``'magma'`` - requires that MAGMA be installed (also
           implemented by Mark Watkins)

        -  ``M`` - Non-negative integer; the modular degree at level `MN` is returned
                   (see above)

        .. note::

            On 64-bit computers ec does not work, so Sage uses sympow
            even if ec is selected on a 64-bit computer.

        The correctness of this function when called with algorithm "sympow"
        is subject to the following three hypothesis:


        -  Manin's conjecture: the Manin constant is 1

        -  Steven's conjecture: the `X_1(N)`-optimal quotient is
           the curve with minimal Faltings height. (This is proved in most
           cases.)

        -  The modular degree fits in a machine double, so it better be
           less than about 50-some bits. (If you use sympow this constraint
           does not apply.)


        Moreover for all algorithms, computing a certain value of an
        `L`-function 'uses a heuristic method that discerns when
        the real-number approximation to the modular degree is within
        epsilon [=0.01 for algorithm='sympow'] of the same integer for 3
        consecutive trials (which occur maybe every 25000 coefficients or
        so). Probably it could just round at some point. For rigour, you
        would need to bound the tail by assuming (essentially) that all the
        `a_n` are as large as possible, but in practice they
        exhibit significant (square root) cancellation. One difficulty is
        that it doesn't do the sum in 1-2-3-4 order; it uses
        1-2-4-8--3-6-12-24-9-18- (Euler product style) instead, and so you
        have to guess ahead of time at what point to curtail this
        expansion.' (Quote from an email of Mark Watkins.)

        .. note::

            If the curve is loaded from the large Cremona database,
            then the modular degree is taken from the database.

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E.modular_degree()
            1
            sage: E = EllipticCurve('5077a')
            sage: E.modular_degree()
            1984
            sage: factor(1984)
            2^6 * 31

        ::

            sage: EllipticCurve([0, 0, 1, -7, 6]).modular_degree()
            1984
            sage: EllipticCurve([0, 0, 1, -7, 6]).modular_degree(algorithm='sympow')
            1984
            sage: EllipticCurve([0, 0, 1, -7, 6]).modular_degree(algorithm='magma')  # optional - magma
            1984

        We compute the modular degree of the curve with rank 4 having
        smallest (known) conductor::

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: factor(E.conductor())  # conductor is 234446
            2 * 117223
            sage: factor(E.modular_degree())
            2^7 * 2617

        Higher level cases::

            sage: E = EllipticCurve('11a')
            sage: for M in range(1,11): print(E.modular_degree(M=M)) # long time (20s on 2009 MBP)
            1
            1
            3
            2
            7
            45
            12
            16
            54
            245
        """
        # Case 1: standard modular degree
        if M==1:
            try:
                return self.__modular_degree

            except AttributeError:
                if algorithm == 'sympow':
                    from sage.lfunctions.all import sympow
                    m = sympow.modular_degree(self)
                elif algorithm == 'magma':
                    from sage.interfaces.all import magma
                    m = rings.Integer(magma(self).ModularDegree())
                else:
                    raise ValueError, "unknown algorithm %s"%algorithm
                self.__modular_degree = m
                return m

        # Case 2: M > 1
        else:
            try:
                return self.__generalized_modular_degree[M]
            except KeyError:
                # self._generalized_congmod_numbers() also populates cache
                return self._generalized_congmod_numbers(M)["moddeg"]


    def modular_parametrization(self):
        r"""
        Returns the modular parametrization of this elliptic curve, which is
        a map from `X_0(N)` to self, where `N` is the conductor of self.

        EXAMPLES::

            sage: E = EllipticCurve('15a')
            sage: phi = E.modular_parametrization(); phi
            Modular parameterization from the upper half plane to Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
            sage: z = 0.1 + 0.2j
            sage: phi(z)
            (8.20822465478531 - 13.1562816054682*I : -8.79855099049364 + 69.4006129342200*I : 1.00000000000000)

        This map is actually a map on `X_0(N)`, so equivalent representatives
        in the upper half plane map to the same point::

            sage: phi((-7*z-1)/(15*z+2))
            (8.20822465478524 - 13.1562816054681*I : -8.79855099049... + 69.4006129342...*I : 1.00000000000000)

        We can also get a series expansion of this modular parameterization::

            sage: E=EllipticCurve('389a1')
            sage: X,Y=E.modular_parametrization().power_series()
            sage: X
            q^-2 + 2*q^-1 + 4 + 7*q + 13*q^2 + 18*q^3 + 31*q^4 + 49*q^5 + 74*q^6 + 111*q^7 + 173*q^8 + 251*q^9 + 379*q^10 + 560*q^11 + 824*q^12 + 1199*q^13 + 1773*q^14 + 2548*q^15 + 3722*q^16 + 5374*q^17 + O(q^18)
            sage: Y
            -q^-3 - 3*q^-2 - 8*q^-1 - 17 - 33*q - 61*q^2 - 110*q^3 - 186*q^4 - 320*q^5 - 528*q^6 - 861*q^7 - 1383*q^8 - 2218*q^9 - 3472*q^10 - 5451*q^11 - 8447*q^12 - 13020*q^13 - 19923*q^14 - 30403*q^15 - 46003*q^16 + O(q^17)

        The following should give 0, but only approximately::

            sage: q = X.parent().gen()
            sage: E.defining_polynomial()(X,Y,1) + O(q^11) == 0
            True
        """
        return ModularParameterization(self)

    def congruence_number(self, M=1):
        r"""
        The case `M==1` corresponds to the classical definition of congruence number:
        Let `X` be the subspace of `S_2(\Gamma_0(N))` spanned by the newform
        associated with this elliptic curve, and `Y` be orthogonal compliment
        of `X` under the Petersson inner product. Let `S_X` and `S_Y` be the
        intersections of `X` and `Y` with `S_2(\Gamma_0(N), \ZZ)`. The congruence
        number is defined to be `[S_X \oplus S_Y : S_2(\Gamma_0(N),\ZZ)]`.
        It measures congruences between `f` and elements of `S_2(\Gamma_0(N),\ZZ)`
        orthogonal to `f`.

        The congruence number for higher levels, when M>1, is defined as above, but
        instead considers `X` to be the subspace of `S_2(\Gamma_0(MN))` spanned by
        embeddings into `S_2(\Gamma_0(MN))` of the newform associated with this
        elliptic curve; this subspace has dimension `\sigma_0(M)`, i.e. the number
        of divisors of `M`. Let `Y` be the orthogonal complement in `S_2(\Gamma_0(MN))`
        of `X` under the Petersson inner product, and `S_X` and `S_Y` the intersections
        of `X` and `Y` with `S_2(\Gamma_0(MN), \ZZ)` respectively. Then the congruence
        number at level `MN` is `[S_X \oplus S_Y : S_2(\Gamma_0(MN),\ZZ)]`.

        INPUT:

        -  ``M`` - Non-negative integer; congruence number is computed at level `MN`,
                   where `N` is the conductor of self.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.congruence_number()
            2
            sage: E.congruence_number()
            2
            sage: E = EllipticCurve('54b')
            sage: E.congruence_number()
            6
            sage: E.modular_degree()
            2
            sage: E = EllipticCurve('242a1')
            sage: E.modular_degree()
            16
            sage: E.congruence_number()  # long time (4s on sage.math, 2011)
            176

        Higher level cases::

            sage: E = EllipticCurve('11a')
            sage: for M in range(1,11): print(E.congruence_number(M)) # long time (20s on 2009 MBP)
            1
            1
            3
            2
            7
            45
            12
            4
            18
            245

        It is a theorem of Ribet that the congruence number (at level `N`) is equal
        to the modular degree in the case of square free conductor. It is a conjecture
        of Agashe, Ribet, and Stein that `ord_p(c_f/m_f) \le ord_p(N)/2`.

        TESTS::

            sage: E = EllipticCurve('11a')
            sage: E.congruence_number()
            1
        """
        # Case 1: M==1
        if M==1:
            try:
                return self.__congruence_number
            except AttributeError:
                pass
            # Currently this is *much* faster to compute
            m = self.modular_degree()
            if self.conductor().is_squarefree():
                self.__congruence_number = m
            else:
                W = self.modular_symbol_space(sign=1)
                V = W.complement().cuspidal_subspace()
                self.__congruence_number = W.congruence_number(V)
                if not m.divides(self.__congruence_number):
                    # We should never get here
                    raise ValueError, "BUG in modular degree or congruence number computation of: %s" % self
            return self.__congruence_number

        # Case 2: M > 1
        else:
            try:
                return self.__generalized_congruence_number[M]
            except KeyError:
                # self._generalized_congmod_numbers() also populates cache
                return self._generalized_congmod_numbers(M)["congnum"]


    def cremona_label(self, space=False):
        """
        Return the Cremona label associated to (the minimal model) of this
        curve, if it is known. If not, raise a RuntimeError exception.

        EXAMPLES::

            sage: E=EllipticCurve('389a1')
            sage: E.cremona_label()
            '389a1'

        The default database only contains conductors up to 10000, so any
        curve with conductor greater than that will cause an error to be
        raised. The optional package 'database_cremona_ellcurve-20120606'
        contains more curves, with conductors up to 240000.

        ::

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: E.conductor()
            234446
            sage: E.cremona_label()  # optional - database_cremona_ellcurve
            '234446a1'
            sage: E = EllipticCurve((0, 0, 1, -79, 342))
            sage: E.conductor()
            19047851
            sage: E.cremona_label()
            Traceback (most recent call last):
            ...
            RuntimeError: Cremona label not known for Elliptic Curve defined by y^2 + y = x^3 - 79*x + 342 over Rational Field.
        """
        try:
            if not space:
                return self.__cremona_label.replace(' ','')
            return self.__cremona_label
        except AttributeError:
            try:
                X = self.database_curve()
            except RuntimeError:
                raise RuntimeError, "Cremona label not known for %s."%self
            self.__cremona_label = X.__cremona_label
            return self.cremona_label(space)

    label = cremona_label

    def reduction(self,p):
       """
       Return the reduction of the elliptic curve at a prime of good
       reduction.

       .. note::

          The actual reduction is done in ``self.change_ring(GF(p))``;
          the reduction is performed after changing to a model which
          is minimal at p.

       INPUT:

       -  ``p`` - a (positive) prime number


       OUTPUT: an elliptic curve over the finite field GF(p)

       EXAMPLES::

           sage: E = EllipticCurve('389a1')
           sage: E.reduction(2)
           Elliptic Curve defined by y^2 + y = x^3 + x^2 over Finite Field of size 2
           sage: E.reduction(3)
           Elliptic Curve defined by y^2 + y = x^3 + x^2 + x over Finite Field of size 3
           sage: E.reduction(5)
           Elliptic Curve defined by y^2 + y = x^3 + x^2 + 3*x over Finite Field of size 5
           sage: E.reduction(38)
           Traceback (most recent call last):
           ...
           AttributeError: p must be prime.
           sage: E.reduction(389)
           Traceback (most recent call last):
           ...
           AttributeError: The curve must have good reduction at p.
           sage: E=EllipticCurve([5^4,5^6])
           sage: E.reduction(5)
           Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5
       """
       p = rings.Integer(p)
       if not p.is_prime():
           raise AttributeError, "p must be prime."
       disc = self.discriminant()
       if not disc.valuation(p) == 0:
           local_data=self.local_data(p)
           if local_data.has_good_reduction():
               return local_data.minimal_model().change_ring(rings.GF(p))
           raise AttributeError, "The curve must have good reduction at p."
       return self.change_ring(rings.GF(p))

    def torsion_order(self):
        """
        Return the order of the torsion subgroup.

        EXAMPLES::

            sage: e = EllipticCurve('11a')
            sage: e.torsion_order()
            5
            sage: type(e.torsion_order())
            <type 'sage.rings.integer.Integer'>
            sage: e = EllipticCurve([1,2,3,4,5])
            sage: e.torsion_order()
            1
            sage: type(e.torsion_order())
            <type 'sage.rings.integer.Integer'>
        """
        try:
            return self.__torsion_order
        except AttributeError:
            self.__torsion_order = self.torsion_subgroup().order()
            return self.__torsion_order

    def _torsion_bound(self,number_of_places = 20):
        r"""
        Computes an upper bound on the order of the torsion group of the
        elliptic curve by counting points modulo several primes of good
        reduction. Note that the upper bound returned by this function is a
        multiple of the order of the torsion group.

        INPUT:


        -  ``number_of_places (default = 20)`` - the number
           of places that will be used to find the bound


        OUTPUT:


        -  ``integer`` - the upper bound


        EXAMPLES:
        """
        E = self
        bound = Integer(0)
        k = 0
        p = Integer(2)   # will run through odd primes
        while k < number_of_places :
            p = p.next_prime()
            # check if the formal group at the place is torsion-free
            # if so the torsion injects into the reduction
            while not E.is_local_integral_model(p) or not E.is_good(p): p = p.next_prime()
            bound = arith.gcd(bound,E.reduction(p).cardinality())
            if bound == 1:
                return bound
            k += 1
        return bound


    def torsion_subgroup(self, algorithm="pari"):
        """
        Returns the torsion subgroup of this elliptic curve.

        INPUT:


        -  ``algorithm`` - string:

        -  ``"pari"`` - (default) use the PARI library

        -  ``"doud"`` - use Doud's algorithm

        -  ``"lutz_nagell"`` - use the Lutz-Nagell theorem


        OUTPUT: The EllipticCurveTorsionSubgroup instance associated to
        this elliptic curve.

        .. note::

           To see the torsion points as a list, use :meth:`.torsion_points`.

        EXAMPLES::

            sage: EllipticCurve('11a').torsion_subgroup()
            Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: EllipticCurve('37b').torsion_subgroup()
            Torsion Subgroup isomorphic to Z/3 associated to the Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field

        ::

            sage: e = EllipticCurve([-1386747,368636886]);e
            Elliptic Curve defined by y^2  = x^3 - 1386747*x + 368636886 over Rational Field
            sage: G = e.torsion_subgroup(); G
            Torsion Subgroup isomorphic to Z/2 + Z/8 associated to the Elliptic
            Curve defined by y^2 = x^3 - 1386747*x + 368636886 over
            Rational Field
            sage: G.0
            (1227 : 22680 : 1)
            sage: G.1
            (282 : 0 : 1)
            sage: list(G)
            [(0 : 1 : 0), (1227 : 22680 : 1), (2307 : -97200 : 1), (8787 : 816480 : 1), (1011 : 0 : 1), (8787 : -816480 : 1), (2307 : 97200 : 1), (1227 : -22680 : 1), (282 : 0 : 1), (-933 : 29160 : 1), (-285 : -27216 : 1), (147 : 12960 : 1), (-1293 : 0 : 1), (147 : -12960 : 1), (-285 : 27216 : 1), (-933 : -29160 : 1)]
        """
        try:
            return self.__torsion_subgroup
        except AttributeError:
            self.__torsion_subgroup = ell_torsion.EllipticCurveTorsionSubgroup(self, algorithm)
            self.__torsion_order = self.__torsion_subgroup.order()
            return self.__torsion_subgroup

    def torsion_points(self, algorithm="pari"):
        """
        Returns the torsion points of this elliptic curve as a sorted
        list.

        INPUT:


        -  ``algorithm`` - string:

           -  "pari" - (default) use the PARI library

           -  "doud" - use Doud's algorithm

           -  "lutz_nagell" - use the Lutz-Nagell theorem


        OUTPUT: A list of all the torsion points on this elliptic curve.

        EXAMPLES::

            sage: EllipticCurve('11a').torsion_points()
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1), (16 : -61 : 1), (16 : 60 : 1)]
            sage: EllipticCurve('37b').torsion_points()
            [(0 : 1 : 0), (8 : -19 : 1), (8 : 18 : 1)]

        ::

            sage: E=EllipticCurve([-1386747,368636886])
            sage: T=E.torsion_subgroup(); T
            Torsion Subgroup isomorphic to Z/2 + Z/8 associated to the Elliptic Curve defined by y^2  = x^3 - 1386747*x + 368636886 over Rational Field
            sage: T == E.torsion_subgroup(algorithm="doud")
            True
            sage: T == E.torsion_subgroup(algorithm="lutz_nagell")
            True
            sage: E.torsion_points()
            [(-1293 : 0 : 1),
            (-933 : -29160 : 1),
            (-933 : 29160 : 1),
            (-285 : -27216 : 1),
            (-285 : 27216 : 1),
            (0 : 1 : 0),
            (147 : -12960 : 1),
            (147 : 12960 : 1),
            (282 : 0 : 1),
            (1011 : 0 : 1),
            (1227 : -22680 : 1),
            (1227 : 22680 : 1),
            (2307 : -97200 : 1),
            (2307 : 97200 : 1),
            (8787 : -816480 : 1),
            (8787 : 816480 : 1)]
        """
        return sorted(self.torsion_subgroup(algorithm).points())

    ## def newform_eval(self, z, prec):
##         """
##         The value of the newform attached to this elliptic curve at
##         the point z in the complex upper half plane, computed using
##         prec terms of the power series expansion.  Note that the power
##         series need not converge well near the real axis.
##         """
##         raise NotImplementedError

    @cached_method
    def root_number(self, p=None):
        """
        Returns the root number of this elliptic curve.

        This is 1 if the order of vanishing of the L-function L(E,s) at 1
        is even, and -1 if it is odd.

        INPUT::

             - `p` -- optional, default (None); if given, return the local
                      root number at `p`

        EXAMPLES::

            sage: EllipticCurve('11a1').root_number()
            1
            sage: EllipticCurve('37a1').root_number()
            -1
            sage: EllipticCurve('389a1').root_number()
            1
            sage: type(EllipticCurve('389a1').root_number())
            <type 'sage.rings.integer.Integer'>

            sage: E = EllipticCurve('100a1')
            sage: E.root_number(2)
            -1
            sage: E.root_number(5)
            1
            sage: E.root_number(7)
            1

        The root number is cached::

            sage: E.root_number(2) is E.root_number(2)
            True
            sage: E.root_number()
            1
        """
        e = self.pari_mincurve()
        if p is None:
            return Integer(e.ellrootno())
        else:
            return Integer(e.ellrootno(p))

    def has_cm(self):
        """
        Returns True iff this elliptic curve has Complex Multiplication.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.has_cm()
            False
            sage: E=EllipticCurve('32a1')
            sage: E.has_cm()
            True
            sage: E.j_invariant()
            1728
        """

        return CMJ.has_key(self.j_invariant())

    def cm_discriminant(self):
        """
        Returns the associated quadratic discriminant if this elliptic
        curve has Complex Multiplication.

        A ValueError is raised if the curve does not have CM (see the
        function has_cm()).

        EXAMPLES::

            sage: E=EllipticCurve('32a1')
            sage: E.cm_discriminant()
            -4
            sage: E=EllipticCurve('121b1')
            sage: E.cm_discriminant()
            -11
            sage: E=EllipticCurve('37a1')
            sage: E.cm_discriminant()
            Traceback (most recent call last):
            ...
            ValueError: Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field does not have CM
        """

        try:
            return CMJ[self.j_invariant()]
        except KeyError:
            raise ValueError, "%s does not have CM"%self


    def quadratic_twist(self, D):
        """
        Return the global minimal model of the quadratic twist of this
        curve by D.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E7=E.quadratic_twist(7); E7
            Elliptic Curve defined by y^2  = x^3 - 784*x + 5488 over Rational Field
            sage: E7.conductor()
            29008
            sage: E7.quadratic_twist(7) == E
            True
        """
        return EllipticCurve_number_field.quadratic_twist(self, D).minimal_model()

    def minimal_quadratic_twist(self):
        r"""
        Determines a quadratic twist with minimal conductor. Returns a
        global minimal model of the twist and the fundamental
        discriminant of the quadratic field over which they are
        isomorphic.

        .. note::

           If there is more than one curve with minimal conductor, the
           one returned is the one with smallest label (if in the
           database), or the one with minimal `a`-invariant list
           (otherwise).

        EXAMPLES::

            sage: E = EllipticCurve('121d1')
            sage: E.minimal_quadratic_twist()
            (Elliptic Curve defined by y^2 + y = x^3 - x^2 over Rational Field, -11)
            sage: Et, D = EllipticCurve('32a1').minimal_quadratic_twist()
            sage: D
            1

            sage: E = EllipticCurve('11a1')
            sage: Et, D = E.quadratic_twist(-24).minimal_quadratic_twist()
            sage: E == Et
            True
            sage: D
            -24

            sage: E = EllipticCurve([0,0,0,0,1000])
            sage: E.minimal_quadratic_twist()
            (Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field, 40)
            sage: E = EllipticCurve([0,0,0,1600,0])
            sage: E.minimal_quadratic_twist()
            (Elliptic Curve defined by y^2 = x^3 + 4*x over Rational Field, 5)


        """
        j = self.j_invariant()
        if j!=0 and j!=1728:
            # the constructor from j will give the minimal twist
            Et = constructor.EllipticCurve_from_j(j)
        else:
            if j==0:
                c = -2*self.c6()
                for p in c.support():
                    e = c.valuation(p)//3
                    c /= p**(3*e)
                E1 = constructor.EllipticCurve([0,0,0,0,c])
            elif j==1728:
                c = -3*self.c4()
                for p in c.support():
                    e = c.valuation(p)//2
                    c /= p**(2*e)
                E1 = constructor.EllipticCurve([0,0,0,c,0])
            tw = [-1,2,-2,3,-3,6,-6]
            Elist = [E1] + [E1.quadratic_twist(t) for t in tw]
            crv_cmp = lambda E,F: cmp(E.conductor(),F.conductor())
            Elist.sort(cmp=crv_cmp)
            Et = Elist[0]

        Et = Et.minimal_model()

        D = self.is_quadratic_twist(Et) # 1 or square-free
        if D % 4 != 1:
            D *= 4

        return Et, D


    ##########################################################
    # Isogeny class
    ##########################################################
    def isogeny_class(self, algorithm="sage", verbose=False, fill_matrix=True, return_maps=False, order=None, use_tuple=True):
        r"""
        Returns all curves over `\QQ` isogenous to this elliptic curve.

        INPUT:

        -  ``algorithm`` - string: one of the following:

           - "mwrank" - use the mwrank C++ library

           - "database" - use the Cremona database (only works if
             curve is isomorphic to a curve in the database)

           - "sage" (default) - use the native Sage implementation.

        - ``verbose`` -- bool (default: False): ignored unless
          ``algorithm`` == "mwrank", in which case it is passed to the
          isogeny class function in mwrank.

        - ``fill_matrix`` -- bool (default: True): See below.

        - ``return_maps`` -- bool (default: False): Ignored unless
          ``algorithm`` == "sage"; then if True, also returns a list of
          lists indexed by pairs of curves in the class such that the
          `(i,j)` entry is the isogeny from curve `i` to curve `j` if
          that isogeny has prime degree, else 0.

        - ``order`` -- None, string, or list of curves (default:
          None): If not None then the curves in the class are
          reordered after being computed.  Note that if the order is
          None then the resulting order will depend on the algorithm.

          - if ``order`` is "mwrank", "database", or "sage" then the
            reordering is so that the order of curves lines up with
            the order produced by that algorithm.

          - if ``order`` is "lmfdb" then the curves are sorted
            lexicographically by a-invariants.

          - if ``order`` is a list of curves, then the curves in the
            class are reordered to be isomorphic with the specified
            list of curves.

        - ``use_tuple`` -- bool (default: True).  Controls the output
          format: see below.  ``use_tuple==True`` is deprecated.

        OUTPUT:

        If ``use_tuple`` is False, returns a
        :class:`sage.schemes.elliptic_curves.isogeny_class.IsogenyClass_EC_Rational`
        instance.  This object models a list of minimal models (with
        containment, index, etc based on isomorphism classes).  It
        also has methods for computing the isogeny matrix and the list
        of isogenies between curves in this class.

        If ``use_tuple`` is True (deprecated), then returns a
        tuple (isogeny class, matrix of integers) or (isogeny class,
        matrix of integers, matrix of isogenies). The sorted list of
        all curves isogenous to self is returned. If ``algorithm`` is
        not "database", the isogeny matrix is also returned, otherwise
        None is returned as the second return value.  When fill_matrix
        is True (default) the `(i,j)` entry of the matrix is a
        positive integer giving the least degree of a cyclic isogeny
        from curve `i` to curve `j` (in particular, the diagonal
        entries are all `1`).  When fill_matrix is False, the
        non-prime entries are replaced by `0` (used in the
        ``isogeny_graph()`` function).

        .. note::

            The curves in the isogeny class are all standard minimal models.

        .. warning::

            With ``algorithm`` "mwrank", the result is *not*
            provably correct, in the sense that when the numbers are
            huge isogenies could be missed because of precision
            issues.  Using ``algorithm`` "sage" avoids this problem,
            though is slower.

        .. note::

           TO DO: when composition of isogenies is implemented, the
           optional returned matrix of isogenies should be completed
           to include all the isogenies, not just those of prime
           degree.

        EXAMPLES::

            sage: isocls = EllipticCurve('37b').isogeny_class('mwrank',order="lmfdb",use_tuple=False)
            sage: isocls
            Elliptic curve isogeny class 37b
            sage: isocls.curves
            (Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 3*x + 1 over Rational Field)
            sage: isocls.matrix()
            [1 3 9]
            [3 1 3]
            [9 3 1]

        ::

            sage: isocls = EllipticCurve('37b').isogeny_class('database', order="lmfdb", use_tuple=False); isocls.curves
            (Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 3*x + 1 over Rational Field)

        This is an example of a curve with a `37`-isogeny::

            sage: E = EllipticCurve([1,1,1,-8,6])
            sage: isocls = E.isogeny_class(use_tuple=False); isocls
            Isogeny class of Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 8*x + 6 over Rational Field
            sage: isocls.matrix()
            [ 1 37]
            [37  1]
            sage: print "\n".join([repr(E) for E in isocls.curves])
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 8*x + 6 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 208083*x - 36621194 over Rational Field

        This curve had numerous `2`-isogenies::

            sage: e=EllipticCurve([1,0,0,-39,90])
            sage: isocls = e.isogeny_class(use_tuple=False); isocls.matrix()
            [1 2 4 4 8 8]
            [2 1 2 2 4 4]
            [4 2 1 4 8 8]
            [4 2 4 1 2 2]
            [8 4 8 2 1 4]
            [8 4 8 2 4 1]

        See http://math.harvard.edu/~elkies/nature.html for more
        interesting examples of isogeny structures.

        ::

            sage: E = EllipticCurve(j = -262537412640768000)
            sage: isocls = E.isogeny_class(algorithm="sage", use_tuple=False); isocls.matrix()
            [  1 163]
            [163   1]
            sage: print "\n".join([repr(C) for C in isocls.curves])
            Elliptic Curve defined by y^2 + y = x^3 - 2174420*x + 1234136692 over Rational Field
            Elliptic Curve defined by y^2 + y = x^3 - 57772164980*x - 5344733777551611 over Rational Field

        For large examples, the "mwrank" algorithm may fail to find
        some isogenies since it works in fixed precision::

            sage: E1 = E.quadratic_twist(6584935282)
            sage: E1.isogeny_class(algorithm="mwrank", use_tuple=False).matrix()
            [1]

        Using algorithm="sage" gives another isogeny (even though
        results are cached: the cache depends on the algorithm)::

            sage: isocls = E1.isogeny_class(algorithm="sage", use_tuple=False); isocls.matrix()
            [  1 163]
            [163   1]
            sage: E1.conductor()
            18433092966712063653330496

        ::

            sage: E = EllipticCurve('14a1')
            sage: isocls = E.isogeny_class(algorithm="sage", use_tuple=False); isocls.matrix()
            [ 1  2  3  3  6  6]
            [ 2  1  6  6  3  3]
            [ 3  6  1  9  2 18]
            [ 3  6  9  1 18  2]
            [ 6  3  2 18  1  9]
            [ 6  3 18  2  9  1]
            sage: print "\n".join([repr(C) for C in isocls.curves])
            Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 36*x - 70 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 171*x - 874 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 11*x + 12 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 2731*x - 55146 over Rational Field
            sage: isocls2 = isocls.reorder('lmfdb'); isocls2.matrix()
            [ 1  2  3  9 18  6]
            [ 2  1  6 18  9  3]
            [ 3  6  1  3  6  2]
            [ 9 18  3  1  2  6]
            [18  9  6  2  1  3]
            [ 6  3  2  6  3  1]
            sage: print "\n".join([repr(C) for C in isocls2.curves])
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 2731*x - 55146 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 171*x - 874 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 36*x - 70 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 11*x + 12 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field

        ::

            sage: E = EllipticCurve('11a1')
            sage: isocls = E.isogeny_class(algorithm="sage", use_tuple=False); isocls.matrix()
            [ 1  5  5]
            [ 5  1 25]
            [ 5 25  1]
            sage: f = isocls.isogenies()[0][1]; f.kernel_polynomial()
            x^2 + x - 29/5
        """
        try:
            isoclass = self._isoclass[algorithm]
        except KeyError:
            from sage.schemes.elliptic_curves.isogeny_class import IsogenyClass_EC_Rational
            if hasattr(self, "_lmfdb_label") and self._lmfdb_label:
                label = self._lmfdb_label[:-1]
            elif hasattr(self, "_EllipticCurve_rational_field__cremona_label") and self.__cremona_label:
                label = self.__cremona_label[:-1]
            else:
                label = None
            isoclass = IsogenyClass_EC_Rational(self, algorithm, verbose, label)
            self._isoclass[algorithm] = isoclass
        if order:
            isoclass = isoclass.reorder(order)
        if use_tuple:
            from sage.misc.superseded import deprecation
            deprecation(12768, """For elliptic curves E over Q, isogeny_class(use_tuple=False) returns a class
that has methods producing the isogeny graph and list of lists of isogenies.
use_tuple=True (currently default) is deprecated.""")
            # After a year or so (in May 2013) we should switch the default to use_tuple=False as the default and deprecate
            # the keyword argument
            if return_maps:
                return isoclass, isoclass.matrix(fill_matrix), isoclass.isogenies()
            else:
                return isoclass, isoclass.matrix(fill_matrix)
        else:
            return isoclass

    def isogenies_prime_degree(self, l=None):
        r"""
        Returns a list of `\ell`-isogenies from self, where `\ell` is a
        prime.

        INPUT:

        - ``l`` -- either None or a prime or a list of primes.

        OUTPUT:

        (list) `\ell`-isogenies for the given `\ell` or if `\ell` is None, all
        `\ell`-isogenies.

        .. note::

           The codomains of the isogenies returned are standard
           minimal models.  This is because the functions
           :meth:`isogenies_prime_degree_genus_0()` and
           :meth:`isogenies_sporadic_Q()` are implemented that way for
           curves defined over `\QQ`.

        EXAMPLES::

            sage: E = EllipticCurve([45,32])
            sage: E.isogenies_prime_degree()
            []
            sage: E = EllipticCurve(j = -262537412640768000)
            sage: E.isogenies_prime_degree()
            [Isogeny of degree 163 from Elliptic Curve defined by y^2 + y = x^3 - 2174420*x + 1234136692 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - 57772164980*x - 5344733777551611 over Rational Field]
            sage: E1 = E.quadratic_twist(6584935282)
            sage: E1.isogenies_prime_degree()
            [Isogeny of degree 163 from Elliptic Curve defined by y^2 = x^3 - 94285835957031797981376080*x + 352385311612420041387338054224547830898 over Rational Field to Elliptic Curve defined by y^2 = x^3 - 2505080375542377840567181069520*x - 1526091631109553256978090116318797845018020806 over Rational Field]

            sage: E = EllipticCurve('14a1')
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - 36*x - 70 over Rational Field]
            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - x over Rational Field, Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 - 171*x - 874 over Rational Field]
            sage: E.isogenies_prime_degree(5)
            []
            sage: E.isogenies_prime_degree(11)
            []
            sage: E.isogenies_prime_degree(29)
            []
            sage: E.isogenies_prime_degree(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not prime.

        """
        from isogeny_small_degree import isogenies_prime_degree_genus_0, isogenies_sporadic_Q

        if l in [2, 3, 5, 7, 13]:
            return isogenies_prime_degree_genus_0(self, l)
        elif l != None and type(l) != list:
            try:
                if l.is_prime(proof=False):
                    return isogenies_sporadic_Q(self, l)
                else:
                    raise ValueError("%s is not prime."%l)
            except AttributeError:
                raise ValueError("%s is not prime."%l)
        if l == None:
            isogs = isogenies_prime_degree_genus_0(self)
            if isogs != []:
                return isogs
            else:
                return isogenies_sporadic_Q(self)
        if type(l) == list:
            isogs = []
            i = 0
            while i<len(l):
                isogenies = [f for f in self.isogenies_prime_degree(l[i]) if not f in isogs]
                isogs.extend(isogenies)
                i = i+1
            return isogs

    def is_isogenous(self, other, proof=True, maxp=200):
        """
        Returns whether or not self is isogenous to other.

        INPUT:

        - ``other`` -- another elliptic curve.

        - ``proof`` (default True) -- If ``False``, the function will
          return ``True`` whenever the two curves have the same
          conductor and are isogenous modulo `p` for `p` up to ``maxp``.
          If ``True``, this test is followed by a rigorous test (which
          may be more time-consuming).

        - ``maxp`` (int, default 200) -- The maximum prime `p` for
          which isogeny modulo `p` will be checked.

        OUTPUT:

        (bool) True if there is an isogeny from curve ``self`` to
        curve ``other``.

        METHOD:

        First the conductors are compared as well as the Traces of
        Frobenius for good primes up to ``maxp``.  If any of these
        tests fail, ``False`` is returned.  If they all pass and
        ``proof`` is ``False`` then ``True`` is returned, otherwise a
        complete set of curves isogenous to ``self`` is computed and
        ``other`` is checked for isomorphism with any of these,

        EXAMPLES::

            sage: E1 = EllipticCurve('14a1')
            sage: E6 = EllipticCurve('14a6')
            sage: E1.is_isogenous(E6)
            True
            sage: E1.is_isogenous(EllipticCurve('11a1'))
            False

        ::

            sage: EllipticCurve('37a1').is_isogenous(EllipticCurve('37b1'))
            False

        ::

            sage: E = EllipticCurve([2, 16])
            sage: EE = EllipticCurve([87, 45])
            sage: E.is_isogenous(EE)
            False
        """
        if not is_EllipticCurve(other):
            raise ValueError, "Second argument is not an Elliptic Curve."
        if not other.base_field() is QQ:
            raise ValueError, "If first argument is an elliptic curve over QQ then the second argument must be also."

        if self.is_isomorphic(other):
            return True

        E1 = self.minimal_model()
        E2 = other.minimal_model()
        D1 = E1.discriminant()
        D2 = E2.discriminant()

        if any([E1.change_ring(rings.GF(p)).cardinality() != E2.change_ring(rings.GF(p)).cardinality() for p in [p for p in rings.prime_range(2,maxp) if D1.valuation(p) == 0 and D2.valuation(p) == 0]]):
            return False

        if E1.conductor() != E2.conductor():
            return False

        if not proof:
            return True
        else:
            return  E2 in E1.isogeny_class(use_tuple=False).curves

    def isogeny_degree(self, other):
        """
        Returns the minimal degree of an isogeny between self and
        other.

        INPUT:

        - ``other`` -- another elliptic curve.

        OUTPUT:

        (int) The minimal degree of an isogeny from ``self`` to
        ``other``, or 0 if the curves are not isogenous.

        EXAMPLES::

            sage: E = EllipticCurve([-1056, 13552])
            sage: E2 = EllipticCurve([-127776, -18037712])
            sage: E.isogeny_degree(E2)
            11

        ::

            sage: E1 = EllipticCurve('14a1')
            sage: E2 = EllipticCurve('14a2')
            sage: E3 = EllipticCurve('14a3')
            sage: E4 = EllipticCurve('14a4')
            sage: E5 = EllipticCurve('14a5')
            sage: E6 = EllipticCurve('14a6')
            sage: E3.isogeny_degree(E1)
            3
            sage: E3.isogeny_degree(E2)
            6
            sage: E3.isogeny_degree(E3)
            1
            sage: E3.isogeny_degree(E4)
            9
            sage: E3.isogeny_degree(E5)
            2
            sage: E3.isogeny_degree(E6)
            18

        ::

            sage: E1 = EllipticCurve('30a1')
            sage: E2 = EllipticCurve('30a2')
            sage: E3 = EllipticCurve('30a3')
            sage: E4 = EllipticCurve('30a4')
            sage: E5 = EllipticCurve('30a5')
            sage: E6 = EllipticCurve('30a6')
            sage: E7 = EllipticCurve('30a7')
            sage: E8 = EllipticCurve('30a8')
            sage: E1.isogeny_degree(E1)
            1
            sage: E1.isogeny_degree(E2)
            2
            sage: E1.isogeny_degree(E3)
            3
            sage: E1.isogeny_degree(E4)
            4
            sage: E1.isogeny_degree(E5)
            4
            sage: E1.isogeny_degree(E6)
            6
            sage: E1.isogeny_degree(E7)
            12
            sage: E1.isogeny_degree(E8)
            12

        ::

            sage: E1 = EllipticCurve('15a1')
            sage: E2 = EllipticCurve('15a2')
            sage: E3 = EllipticCurve('15a3')
            sage: E4 = EllipticCurve('15a4')
            sage: E5 = EllipticCurve('15a5')
            sage: E6 = EllipticCurve('15a6')
            sage: E7 = EllipticCurve('15a7')
            sage: E8 = EllipticCurve('15a8')
            sage: E1.isogeny_degree(E1)
            1
            sage: E7.isogeny_degree(E2)
            8
            sage: E7.isogeny_degree(E3)
            2
            sage: E7.isogeny_degree(E4)
            8
            sage: E7.isogeny_degree(E5)
            16
            sage: E7.isogeny_degree(E6)
            16
            sage: E7.isogeny_degree(E8)
            4

        0 is returned when the curves are not isogenous::

            sage: A = EllipticCurve('37a1')
            sage: B = EllipticCurve('37b1')
            sage: A.isogeny_degree(B)
            0
            sage: A.is_isogenous(B)
            False
        """
        E1 = self.minimal_model()
        E2 = other.minimal_model()

        if not E1.is_isogenous(E2, proof=False):
            return Integer(0)

        isocls = E1.isogeny_class(use_tuple=False)
        try:
            return isocls.matrix(fill=True)[0,isocls.index(E2)]
        except ValueError:
            return Integer(0)

#
#     The following function can be implemented once composition of
#     isogenies has been implemented.
#
#     def contruct_isogeny(self, other):
#         """
#         Returns an isogeny from self to other if the two curves are in
#         the same isogeny class.
#         """


    def optimal_curve(self):
        """
        Given an elliptic curve that is in the installed Cremona
        database, return the optimal curve isogenous to it.

        EXAMPLES:

        The following curve is not optimal::

            sage: E = EllipticCurve('11a2'); E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field
            sage: E.optimal_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E.optimal_curve().cremona_label()
            '11a1'

        Note that 990h is the special case where the optimal curve
        isn't the first in the Cremona labeling::

            sage: E = EllipticCurve('990h4'); E
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 + 6112*x - 41533 over Rational Field
            sage: F = E.optimal_curve(); F
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 1568*x - 4669 over Rational Field
            sage: F.cremona_label()
            '990h3'
            sage: EllipticCurve('990a1').optimal_curve().cremona_label()   # a isn't h.
            '990a1'

        If the input curve is optimal, this function returns that
        curve (not just a copy of it or a curve isomorphic to it!)::

            sage: E = EllipticCurve('37a1')
            sage: E.optimal_curve() is E
            True

        Also, if this curve is optimal but not given by a minimal
        model, this curve will still be returned, so this function
        need not return a minimal model in general.

        ::

            sage: F = E.short_weierstrass_model(); F
            Elliptic Curve defined by y^2  = x^3 - 16*x + 16 over Rational Field
            sage: F.optimal_curve()
            Elliptic Curve defined by y^2  = x^3 - 16*x + 16 over Rational Field
        """
        label = self.cremona_label()
        N, isogeny, number = sage.databases.cremona.parse_cremona_label(label)
        if N == 990 and isogeny == 'h':
            optimal_label = '990h3'
        else:
            optimal_label = '%s%s1'%(N,isogeny)
        if optimal_label == label: return self
        return constructor.EllipticCurve(optimal_label)

    def isogeny_graph(self, order=None):
        r"""
        Returns a graph representing the isogeny class of this elliptic
        curve, where the vertices are isogenous curves over
        `\QQ` and the edges are prime degree isogenies

        .. note:

            The vertices are labeled 1 to n rather than 0 to n-1 to
            correspond to LMFDB and Cremona labels.

        EXAMPLES::

            sage: LL = []
            sage: for e in cremona_optimal_curves(range(1, 38)):
            ...    G = e.isogeny_graph()
            ...    already = False
            ...    for H in LL:
            ...        if G.is_isomorphic(H):
            ...            already = True
            ...            break
            ...    if not already:
            ...        LL.append(G)
            ...
            sage: graphs_list.show_graphs(LL)

        ::

            sage: E = EllipticCurve('195a')
            sage: G = E.isogeny_graph()
            sage: for v in G: print v, G.get_vertex(v)
            ...
            1 Elliptic Curve defined by y^2 + x*y  = x^3 - 110*x + 435 over Rational Field
            2 Elliptic Curve defined by y^2 + x*y  = x^3 - 115*x + 392 over Rational Field
            3 Elliptic Curve defined by y^2 + x*y  = x^3 + 210*x + 2277 over Rational Field
            4 Elliptic Curve defined by y^2 + x*y  = x^3 - 520*x - 4225 over Rational Field
            5 Elliptic Curve defined by y^2 + x*y  = x^3 + 605*x - 19750 over Rational Field
            6 Elliptic Curve defined by y^2 + x*y  = x^3 - 8125*x - 282568 over Rational Field
            7 Elliptic Curve defined by y^2 + x*y  = x^3 - 7930*x - 296725 over Rational Field
            8 Elliptic Curve defined by y^2 + x*y  = x^3 - 130000*x - 18051943 over Rational Field
            sage: G.plot(edge_labels=True)
        """
        return self.isogeny_class(algorithm='mwrank', order=order, use_tuple=False).graph()

    def manin_constant(self):
        r"""
        Return the Manin constant of this elliptic curve. If `\phi: X_0(N) \to E` is the modular
        parametrization of minimal degree, then the Manin constant `c` is defined to be the rational
        number `c` such that `\phi^*(\omega_E) = c\cdot \omega_f` where `\omega_E` is a Neron differential and `\omega_f = f(q) dq/q` is the differential on `X_0(N)` corresponding to the
        newform `f` attached to the isogeny class of `E`.

        It is known that the Manin constant is an integer. It is conjectured that in each class there is at least one, more precisely the so-called strong Weil curve or `X_0(N)`-optimal curve, that has Manin constant `1`.

        OUTPUT:

        an integer

        This function only works if the curve is in the installed
        Cremona database.  Sage includes by default a small databases;
        for the full database you have to install an optional package.

        EXAMPLES::

            sage: EllipticCurve('11a1').manin_constant()
            1
            sage: EllipticCurve('11a2').manin_constant()
            1
            sage: EllipticCurve('11a3').manin_constant()
            5

        Check that it works even if the curve is non-minimal::

            sage: EllipticCurve('11a3').change_weierstrass_model([1/35,0,0,0]).manin_constant()
            5

        Rather complicated examples (see #12080) ::

            sage: [ EllipticCurve('27a%s'%i).manin_constant() for i in [1,2,3,4]]
            [1, 1, 3, 3]
            sage: [ EllipticCurve('80b%s'%i).manin_constant() for i in [1,2,3,4]]
            [1, 2, 1, 2]

        """
        from sage.databases.cremona import CremonaDatabase

        if self.conductor() > CremonaDatabase().largest_conductor():
            raise NotImplementedError, "The Manin constant can only be evaluated for curves in Cremona's tables. If you have not done so, you may wish to install the optional large database."

        E = self.minimal_model()
        C = self.optimal_curve()
        m = C.isogeny_class(use_tuple=False).matrix()
        ma = max(max(x) for x in m)
        OmC = C.period_lattice().basis()
        OmE = E.period_lattice().basis()
        q_plus = QQ(gp.bestappr(OmE[0]/OmC[0],ma+1) )
        n_plus = q_plus.numerator()

        cinf_E = E.real_components()
        cinf_C = C.real_components()
        OmC_minus = OmC[1].imag()
        if cinf_C == 1:
            OmC_minus *= 2
        OmE_minus = OmE[1].imag()
        if cinf_E == 1:
            OmE_minus *= 2
        q_minus = QQ(gp.bestappr(OmE_minus/OmC_minus, ma+1))
        n_minus = q_minus.numerator()
        n = ZZ(n_minus * n_plus)

        if cinf_C == cinf_E:
            return n
        # otherwise they have different number of connected component and we have to adjust for this
        elif cinf_C > cinf_E:
            if ZZ(n_plus) % 2 == 0 and ZZ(n_minus) % 2 == 0:
                return n // 2
            else:
                return n
        else: #if cinf_C < cinf_E:
            if q_plus.denominator() % 2 == 0 and q_minus.denominator() % 2 == 0:
                return n
            else:
                return n*2

    def _shortest_paths(self):
        r"""
        Technical internal function that returns the list of isogenies
        curves and corresponding dictionary of shortest isogeny paths
        from self to each other curve in the isogeny class.

        .. warning::

           The result is *not* provably correct, in the
           sense that when the numbers are huge isogenies could be
           missed because of precision issues.

        OUTPUT:

        list, dict

        EXAMPLES::

            sage: EllipticCurve('11a1')._shortest_paths()
            ((Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field,
              Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field,
              Elliptic Curve defined by y^2 + y = x^3 - x^2 over Rational Field),
             {0: 0, 1: 5, 2: 5})
            sage: EllipticCurve('11a2')._shortest_paths()
            ((Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field,
              Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field,
              Elliptic Curve defined by y^2 + y = x^3 - x^2 over Rational Field),
             {0: 0, 1: 5, 2: 25})
        """
        from sage.graphs.graph import Graph
        isocls = self.isogeny_class(algorithm='mwrank',use_tuple=False)
        M = isocls.matrix(fill=True).change_ring(rings.RR)
        # see trac #4889 for nebulous M.list() --> M.entries() change...
        # Take logs here since shortest path minimizes the *sum* of the weights -- not the product.
        M = M.parent()([a.log() if a else 0 for a in M.list()])
        G = Graph(M, format='weighted_adjacency_matrix')
        G.set_vertices(dict([(v,isocls[v]) for v in G.vertices()]))
        v = G.shortest_path_lengths(0, by_weight=True, weight_sums=True)
        # Now exponentiate and round to get degrees of isogenies
        v = dict([(i, j.exp().round() if j else 0) for i,j in v.iteritems()])
        return isocls.curves, v

    def _multiple_of_degree_of_isogeny_to_optimal_curve(self):
        r"""
        Internal function returning an integer m such that the degree of
        the isogeny between this curve and the optimal curve in its
        isogeny class is a divisor of m.

        .. warning::

           The result is *not* provably correct, in the
           sense that when the numbers are huge isogenies could be
           missed because of precision issues.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: E._multiple_of_degree_of_isogeny_to_optimal_curve()
            5
            sage: E = EllipticCurve('11a2')
            sage: E._multiple_of_degree_of_isogeny_to_optimal_curve()
            25
            sage: E = EllipticCurve('11a3')
            sage: E._multiple_of_degree_of_isogeny_to_optimal_curve()
            25
        """
        _, v = self._shortest_paths()
        # Compute the degree of an isogeny from self to anything else
        # in the isogeny class of self.  Assuming the isogeny
        # enumeration is complete (which need not be the case a priori!), the LCM
        # of these numbers is a multiple of the degree of the isogeny
        # to the optimal curve.
        v = [deg for num, deg in v.iteritems() if deg]  # get just the degrees
        return arith.LCM(v)

    ##########################################################
    # Galois Representations
    ##########################################################

    def galois_representation(self):
        r"""
        The compatible family of the Galois representation
        attached to this elliptic curve.

        Given an elliptic curve `E` over `\QQ`
        and a rational prime number `p`, the `p^n`-torsion
        `E[p^n]` points of `E` is a representation of the
        absolute Galois group of `\QQ`. As `n` varies
        we obtain the Tate module `T_p E` which is a
        a representation of `G_K` on a free `\ZZ_p`-module
        of rank `2`. As `p` varies the representations
        are compatible.

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: rho.is_irreducible(7)
            True
            sage: rho.is_irreducible(5)
            False
            sage: rho.is_surjective(11)
            True
            sage: rho.non_surjective()
            [5]
            sage: rho = EllipticCurve('37a1').galois_representation()
            sage: rho.non_surjective()
            []
            sage: rho = EllipticCurve('27a1').galois_representation()
            sage: rho.is_irreducible(7)
            True
            sage: rho.non_surjective()   # cm-curve
            [0]

       """
        try:
            return self.__rho
        except AttributeError:
            from gal_reps import GaloisRepresentation
            self.__rho = GaloisRepresentation(self)
        return self.__rho

    # deprecated as it should be the is_reducible for a scheme (and hence return False always).
    def is_reducible(self, p):
        """
        Return True if the mod-p representation attached to E is
        reducible.

        Note that this function is deprecated, and that you should use
        galois_representation().is_reducible(p) instead as this will be
        disappearing in the near future.

        EXAMPLES::

            sage: EllipticCurve('20a1').is_reducible(3) #random
            doctest:1: DeprecationWarning: is_reducible is deprecated, use galois_representation().is_reducible(p) instead!
            True

        """
        from sage.misc.superseded import deprecation
        deprecation(8118, "is_reducible is deprecated, use galois_representation().is_reducible(p) instead!")
        return self.galois_representation().is_reducible(p)

    # deprecated as it should be the is_irreducible for a scheme (and hence return True always).
    def is_irreducible(self, p):
        """
        Return True if the mod p representation is irreducible.

        Note that this function is deprecated, and that you should use
        galois_representation().is_irreducible(p) instead as this will be
        disappearing in the near future.

        EXAMPLES::

            sage: EllipticCurve('20a1').is_irreducible(7) #random
            doctest:1: DeprecationWarning: is_irreducible is deprecated, use galois_representation().is_irreducible(p) instead!
            True

        """
        from sage.misc.superseded import deprecation
        deprecation(8118, "is_irreducible is deprecated, use galois_representation().is_irreducible(p) instead!")
        return self.galois_representation().is_irreducible(p)

    # deprecated
    def is_surjective(self, p, A=1000):
        r"""
        Returns true if the mod p representation is surjective

        Note that this function is deprecated, and that you should use
        galois_representation().is_surjective(p) instead as this will be
        disappearing in the near future.

        EXAMPLES::

            sage: EllipticCurve('20a1').is_surjective(7) #random
            doctest:1: DeprecationWarning: is_surjective is deprecated, use galois_representation().is_surjective(p) instead!
            True

        """
        from sage.misc.superseded import deprecation
        deprecation(8118, "is_surjective is deprecated, use galois_representation().is_surjective(p) instead!")
        return self.galois_representation().is_surjective(p,A)

    # deprecated
    def reducible_primes(self):
        r"""
        Returns a list of reducible primes.

        Note that this function is deprecated, and that you should use
        galois_representation().reducible_primes() instead as this will be
        disappearing in the near future.

        EXAMPLES::

            sage: EllipticCurve('20a1').reducible_primes() #random
            doctest:1: DeprecationWarning: reducible_primes is deprecated, use galois_representation().reducible_primes() instead!
            [2,3]

       """
        from sage.misc.superseded import deprecation
        deprecation(8118, "reducible_primes is deprecated, use galois_representation().reducible_primes() instead!")
        return self.galois_representation().reducible_primes()

    # deprecated
    def non_surjective(self, A=1000):
        r"""
        Returns a list of primes p for which the Galois representation mod p is not surjective.

        Note that this function is deprecated, and that you should use
        galois_representation().non_surjective() instead as this will be
        disappearing in the near future.

        EXAMPLES::

            sage: EllipticCurve('20a1').non_surjective() #random
            doctest:1: DeprecationWarning: non_surjective is deprecated, use galois_representation().non_surjective() instead!
            [2,3]

        """
        from sage.misc.superseded import deprecation
        deprecation(8118, "non_surjective is deprecated, use galois_representation().non_surjective() instead!")
        return self.galois_representation().non_surjective()

    def is_semistable(self):
        """
        Return True iff this elliptic curve is semi-stable at all primes.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.is_semistable()
            True
            sage: E=EllipticCurve('90a1')
            sage: E.is_semistable()
            False
        """
        return self.conductor().is_squarefree()

    def is_ordinary(self, p, ell=None):
        """
        Return True precisely when the mod-p representation attached to
        this elliptic curve is ordinary at ell.

        INPUT:

        -  ``p`` - a prime ell - a prime (default: p)

        OUTPUT: bool

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.is_ordinary(37)
            True
            sage: E=EllipticCurve('32a1')
            sage: E.is_ordinary(2)
            False
            sage: [p for p in prime_range(50) if E.is_ordinary(p)]
            [5, 13, 17, 29, 37, 41]

        """
        if ell is None:
            ell = p
        return self.ap(ell) % p != 0

    def is_good(self, p, check=True):
        """
        Return True if `p` is a prime of good reduction for
        `E`.

        INPUT:

        -  ``p`` - a prime

        OUTPUT: bool

        EXAMPLES::

            sage: e = EllipticCurve('11a')
            sage: e.is_good(-8)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
            sage: e.is_good(-8, check=False)
            True

        """
        if check:
            if not arith.is_prime(p):
                raise ValueError, "p must be prime"
        return self.conductor() % p != 0


    def is_supersingular(self, p, ell=None):
        """
        Return True precisely when p is a prime of good reduction and the
        mod-p representation attached to this elliptic curve is
        supersingular at ell.

        INPUT:

        -  ``p`` - a prime ell - a prime (default: p)

        OUTPUT: bool

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.is_supersingular(37)
            False
            sage: E=EllipticCurve('32a1')
            sage: E.is_supersingular(2)
            False
            sage: E.is_supersingular(7)
            True
            sage: [p for p in prime_range(50) if E.is_supersingular(p)]
            [3, 7, 11, 19, 23, 31, 43, 47]

        """
        if ell is None:
            ell = p
        return self.is_good(p) and not self.is_ordinary(p, ell)

    def supersingular_primes(self, B):
        """
        Return a list of all supersingular primes for this elliptic curve
        up to and possibly including B.

        EXAMPLES::

            sage: e = EllipticCurve('11a')
            sage: e.aplist(20)
            [-2, -1, 1, -2, 1, 4, -2, 0]
            sage: e.supersingular_primes(1000)
            [2, 19, 29, 199, 569, 809]

        ::

            sage: e = EllipticCurve('27a')
            sage: e.aplist(20)
            [0, 0, 0, -1, 0, 5, 0, -7]
            sage: e.supersingular_primes(97)
            [2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89]
            sage: e.ordinary_primes(97)
            [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97]
            sage: e.supersingular_primes(3)
            [2]
            sage: e.supersingular_primes(2)
            [2]
            sage: e.supersingular_primes(1)
            []
        """
        v = self.aplist(max(B, 3))
        P = rings.prime_range(max(B,3)+1)
        N = self.conductor()
        return [P[i] for i in [0,1] if P[i] <= B and v[i]%P[i]==0 and N%P[i] != 0] + \
                      [P[i] for i in range(2,len(v)) if v[i] == 0 and N%P[i] != 0]

    def ordinary_primes(self, B):
        """
        Return a list of all ordinary primes for this elliptic curve up to
        and possibly including B.

        EXAMPLES::

            sage: e = EllipticCurve('11a')
            sage: e.aplist(20)
            [-2, -1, 1, -2, 1, 4, -2, 0]
            sage: e.ordinary_primes(97)
            [3, 5, 7, 11, 13, 17, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
            sage: e = EllipticCurve('49a')
            sage: e.aplist(20)
            [1, 0, 0, 0, 4, 0, 0, 0]
            sage: e.supersingular_primes(97)
            [3, 5, 13, 17, 19, 31, 41, 47, 59, 61, 73, 83, 89, 97]
            sage: e.ordinary_primes(97)
            [2, 11, 23, 29, 37, 43, 53, 67, 71, 79]
            sage: e.ordinary_primes(3)
            [2]
            sage: e.ordinary_primes(2)
            [2]
            sage: e.ordinary_primes(1)
            []
        """
        v = self.aplist(max(B, 3) )
        P = rings.prime_range(max(B,3) +1)
        return [P[i] for i in [0,1] if P[i] <= B and v[i]%P[i]!=0] +\
               [P[i] for i in range(2,len(v)) if v[i] != 0]

    def eval_modular_form(self, points, prec):
        """
        Evaluate the modular form of this elliptic curve at points in CC

        INPUT:


        -  ``points`` - a list of points in the half-plane of
           convergence

        -  ``prec`` - precision


        OUTPUT: A list of values L(E,s) for s in points

        .. note::

           Better examples are welcome.

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: E.eval_modular_form([1.5+I,2.0+I,2.5+I],0.000001)
            [0, 0, 0]
        """
        if not isinstance(points, (list,xrange)):
            try:
                points = list(points)
            except TypeError:
                return self.eval_modular_form([points],prec)
        an = self.pari_mincurve().ellan(prec)
        s = 0
        c = pari('2 * Pi * I')
        ans = []
        for z in points:
            s = pari(0)
            r0 = (c*z).exp()
            r = r0
            for n in xrange(1,prec):
                s += an[n-1]*r
                r *= r0
            ans.append(s.python())
        return ans


    ########################################################################
    # The Tate-Shafarevich group
    ########################################################################

    def sha(self):
        """
        Return an object of class
        'sage.schemes.elliptic_curves.sha_tate.Sha' attached to this
        elliptic curve.

        This can be used in functions related to bounding the order of Sha
        (The Tate-Shafarevich group of the curve).

        EXAMPLES::

            sage: E=EllipticCurve('37a1')
            sage: S=E.sha()
            sage: S
            Tate-Shafarevich group for the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: S.bound_kolyvagin()
            ([2], 1)
        """
        try:
            return self.__sha
        except AttributeError:
            from sha_tate import Sha
            self.__sha = Sha(self)
            return self.__sha

    #################################################################################
    # Functions related to Heegner points#################################################################################
    heegner_point = heegner.ell_heegner_point
    kolyvagin_point = heegner.kolyvagin_point

    heegner_discriminants = heegner.ell_heegner_discriminants
    heegner_discriminants_list = heegner.ell_heegner_discriminants_list
    satisfies_heegner_hypothesis = heegner.satisfies_heegner_hypothesis

    heegner_point_height = heegner.heegner_point_height

    heegner_index = heegner.heegner_index
    _adjust_heegner_index = heegner._adjust_heegner_index
    heegner_index_bound = heegner.heegner_index_bound
    _heegner_index_in_EK = heegner._heegner_index_in_EK

    heegner_sha_an = heegner.heegner_sha_an

    _heegner_forms_list = heegner._heegner_forms_list
    _heegner_best_tau = heegner._heegner_best_tau

    #################################################################################
    # p-adic functions
    #################################################################################

    padic_regulator = padics.padic_regulator

    padic_height_pairing_matrix = padics.padic_height_pairing_matrix

    padic_height = padics.padic_height
    padic_height_via_multiply = padics.padic_height_via_multiply

    padic_sigma = padics.padic_sigma
    padic_sigma_truncated = padics.padic_sigma_truncated

    padic_E2 = padics.padic_E2

    matrix_of_frobenius = padics.matrix_of_frobenius

    # def weierstrass_p(self):
    #         # TODO: add allowing negative valuations for power series
    #         return 1/t**2 + a1/t + rings.frac(1,12)*(a1-8*a2) -a3*t \
    #                - (a4+a1*a3)*t**2  + O(t**3)


    def mod5family(self):
        """
        Return the family of all elliptic curves with the same mod-5
        representation as self.

        EXAMPLES::

            sage: E=EllipticCurve('32a1')
            sage: E.mod5family()
            Elliptic Curve defined by y^2  = x^3 + 4*x over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        E = self.short_weierstrass_model()
        a = E.a4()
        b = E.a6()
        return mod5family.mod5family(a,b)

    def tate_curve(self, p):
        r"""
        Creates the Tate Curve over the `p`-adics associated to
        this elliptic curves.

        This Tate curve a `p`-adic curve with split multiplicative
        reduction of the form `y^2+xy=x^3+s_4 x+s_6` which is
        isomorphic to the given curve over the algebraic closure of
        `\QQ_p`. Its points over `\QQ_p`
        are isomorphic to `\QQ_p^{\times}/q^{\ZZ}`
        for a certain parameter `q\in\ZZ_p`.

        INPUT:

        p - a prime where the curve has multiplicative reduction.

        EXAMPLES::

            sage: e = EllipticCurve('130a1')
            sage: e.tate_curve(2)
            2-adic Tate curve associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field

        The input curve must have multiplicative reduction at the prime.

        ::

            sage: e.tate_curve(3)
            Traceback (most recent call last):
            ...
            ValueError: The elliptic curve must have multiplicative reduction at 3

        We compute with `p=5`::

            sage: T = e.tate_curve(5); T
            5-adic Tate curve associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field

        We find the Tate parameter `q`::

            sage: T.parameter(prec=5)
            3*5^3 + 3*5^4 + 2*5^5 + 2*5^6 + 3*5^7 + O(5^8)

        We compute the `\mathcal{L}`-invariant of the curve::

            sage: T.L_invariant(prec=10)
            5^3 + 4*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 3*5^8 + 5^9 + O(5^10)
        """
        try:
            return self._tate_curve[p]
        except AttributeError:
            self._tate_curve = {}
        except KeyError:
            pass

        Eq = ell_tate_curve.TateCurve(self,p)
        self._tate_curve[p] = Eq
        return Eq

    def height(self, precision=None):
        """
        Returns the real height of this elliptic curve. This is used in
        integral_points()

        INPUT:


        -  ``precision`` - desired real precision of the result
           (default real precision if None)


        EXAMPLES::

            sage: E=EllipticCurve('5077a1')
            sage: E.height()
            17.4513334798896
            sage: E.height(100)
            17.451333479889612702508579399
            sage: E=EllipticCurve([0,0,0,0,1])
            sage: E.height()
            1.38629436111989
            sage: E=EllipticCurve([0,0,0,1,0])
            sage: E.height()
            7.45471994936400
        """
        if precision is None:
            precision = RealField().precision()
        R = RealField(precision)
        c4 = self.c4()
        c6 = self.c6()
        j = self.j_invariant()
        log_g2 = R((c4/12)).abs().log()
        log_g3 = R((c6/216)).abs().log()

        if j == 0:
            h_j = R(1)
        else:
            h_j = max(log(R(abs(j.numerator()))), log(R(abs(j.denominator()))))
        if (self.c4() != 0) and (self.c6() != 0):
            h_gs = max(1, log_g2, log_g3)
        elif c4 == 0:
            if c6 == 0:
                return max(1,h_j)
            h_gs = max(1, log_g3)
        else:
            h_gs = max(1, log_g2)
        return max(R(1),h_j, h_gs)

    def antilogarithm(self, z, max_denominator=None):
        r"""
        Returns the rational point (if any) associated to this complex
        number; the inverse of the elliptic logarithm function.

        INPUT:

        -  ``z`` -- a complex number representing an element of
           `\CC/L` where `L` is the period lattice of the elliptic curve

        - ``max_denominator`` (int or None) -- parameter controlling
          the attempted conversion of real numbers to rationals.  If
          None, ``simplest_rational()`` will be used; otherwise,
          ``nearby_rational()`` will be used with this value of
          ``max_denominator``.

        OUTPUT:

        - point on the curve: the rational point which is the
          image of `z` under the Weierstrass parametrization, if it
          exists and can be determined from `z` and the given value
          of max_denominator (if any); otherwise a ValueError exception
          is raised.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E(-1,1)
            sage: z = P.elliptic_logarithm()
            sage: E.antilogarithm(z)
            (-1 : 1 : 1)
            sage: Q = E(0,-1)
            sage: z = Q.elliptic_logarithm()
            sage: E.antilogarithm(z)
            Traceback (most recent call last):
            ...
            ValueError: approximated point not on the curve
            sage: E.antilogarithm(z, max_denominator=10)
            (0 : -1 : 1)

            sage: E = EllipticCurve('11a1')
            sage: w1,w2 = E.period_lattice().basis()
            sage: [E.antilogarithm(a*w1/5,1) for a in range(5)]
            [(0 : 1 : 0), (16 : -61 : 1), (5 : -6 : 1), (5 : 5 : 1), (16 : 60 : 1)]
        """
        if z.is_zero():
            return self(0)
        expZ = self.elliptic_exponential(z)
        xy = [t.real() for t in expZ[:2]]
        if max_denominator is None:
            xy = [t.simplest_rational() for t in xy]
        else:
            xy = [t.nearby_rational(max_denominator=max_denominator) for t in xy]
        try:
            return self(xy)
        except TypeError:
            raise ValueError, "approximated point not on the curve"

    def integral_x_coords_in_interval(self,xmin,xmax):
        r"""
        Returns the set of integers `x` with `xmin\le x\le xmax` which are
        `x`-coordinates of rational points on this curve.

        INPUT:

        - ``xmin``, ``xmax`` (integers) -- two integers.

        OUTPUT:

        (set) The set of integers `x` with `xmin\le x\le xmax` which
        are `x`-coordinates of rational points on the elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -7, 6])
            sage: xset = E.integral_x_coords_in_interval(-100,100)
            sage: xlist = list(xset); xlist.sort(); xlist
            [-3, -2, -1, 0, 1, 2, 3, 4, 8, 11, 14, 21, 37, 52, 93]

        TODO: re-implement this using the much faster point searching
        implemented in Stoll's ``ratpoints`` program.

        """
        xmin=Integer(xmin)
        xmax=Integer(xmax)
        ans = set([])
        x = xmin
        while x<=xmax:
            if self.is_x_coord(x):
                ans.add(x)
            x+=1
        return ans

    prove_BSD = BSD.prove_BSD

    def integral_points(self, mw_base='auto', both_signs=False, verbose=False):
        """
        Computes all integral points (up to sign) on this elliptic curve.

        INPUT:


        -  ``mw_base`` - list of EllipticCurvePoint generating
           the Mordell-Weil group of E (default: 'auto' - calls self.gens())

        -  ``both_signs`` - True/False (default False): if
           True the output contains both P and -P, otherwise only one of each
           pair.

        -  ``verbose`` - True/False (default False): if True,
           some details of the computation are output


        OUTPUT: A sorted list of all the integral points on E (up to sign
        unless both_signs is True)

        .. note::

           The complexity increases exponentially in the rank of curve
           E. The computation time (but not the output!) depends on
           the Mordell-Weil basis. If mw_base is given but is not a
           basis for the Mordell-Weil group (modulo torsion), integral
           points which are not in the subgroup generated by the given
           points will almost certainly not be listed.

        EXAMPLES: A curve of rank 3 with no torsion points

        ::

            sage: E=EllipticCurve([0,0,1,-7,6])
            sage: P1=E.point((2,0)); P2=E.point((-1,3)); P3=E.point((4,6))
            sage: a=E.integral_points([P1,P2,P3]); a
            [(-3 : 0 : 1), (-2 : 3 : 1), (-1 : 3 : 1), (0 : 2 : 1), (1 : 0 : 1), (2 : 0 : 1), (3 : 3 : 1), (4 : 6 : 1), (8 : 21 : 1), (11 : 35 : 1), (14 : 51 : 1), (21 : 95 : 1), (37 : 224 : 1), (52 : 374 : 1), (93 : 896 : 1), (342 : 6324 : 1), (406 : 8180 : 1), (816 : 23309 : 1)]

        ::

            sage: a = E.integral_points([P1,P2,P3], verbose=True)
            Using mw_basis  [(2 : 0 : 1), (3 : -4 : 1), (8 : -22 : 1)]
            e1,e2,e3:  -3.0124303725933... 1.0658205476962... 1.94660982489710
            Minimal eigenvalue of height pairing matrix:  0.63792081458500...
            x-coords of points on compact component with  -3 <=x<= 1
            [-3, -2, -1, 0, 1]
            x-coords of points on non-compact component with  2 <=x<= 6
            [2, 3, 4]
            starting search of remaining points using coefficient bound  5
            x-coords of extra integral points:
            [2, 3, 4, 8, 11, 14, 21, 37, 52, 93, 342, 406, 816]
            Total number of integral points: 18

        It is not necessary to specify mw_base; if it is not provided,
        then the Mordell-Weil basis must be computed, which may take
        much longer.

        ::

            sage: E=EllipticCurve([0,0,1,-7,6])
            sage: a=E.integral_points(both_signs=True); a
            [(-3 : -1 : 1), (-3 : 0 : 1), (-2 : -4 : 1), (-2 : 3 : 1), (-1 : -4 : 1), (-1 : 3 : 1), (0 : -3 : 1), (0 : 2 : 1), (1 : -1 : 1), (1 : 0 : 1), (2 : -1 : 1), (2 : 0 : 1), (3 : -4 : 1), (3 : 3 : 1), (4 : -7 : 1), (4 : 6 : 1), (8 : -22 : 1), (8 : 21 : 1), (11 : -36 : 1), (11 : 35 : 1), (14 : -52 : 1), (14 : 51 : 1), (21 : -96 : 1), (21 : 95 : 1), (37 : -225 : 1), (37 : 224 : 1), (52 : -375 : 1), (52 : 374 : 1), (93 : -897 : 1), (93 : 896 : 1), (342 : -6325 : 1), (342 : 6324 : 1), (406 : -8181 : 1), (406 : 8180 : 1), (816 : -23310 : 1), (816 : 23309 : 1)]

        An example with negative discriminant::

            sage: EllipticCurve('900d1').integral_points()
            [(-11 : 27 : 1), (-4 : 34 : 1), (4 : 18 : 1), (16 : 54 : 1)]

        Another example with rank 5 and no torsion points::

            sage: E=EllipticCurve([-879984,319138704])
            sage: P1=E.point((540,1188)); P2=E.point((576,1836))
            sage: P3=E.point((468,3132)); P4=E.point((612,3132))
            sage: P5=E.point((432,4428))
            sage: a=E.integral_points([P1,P2,P3,P4,P5]); len(a)  # long time (18s on sage.math, 2011)
            54

        TESTS:

        The bug reported on trac #4525 is now fixed::

            sage: EllipticCurve('91b1').integral_points()
            [(-1 : 3 : 1), (1 : 0 : 1), (3 : 4 : 1)]

        ::

            sage: [len(e.integral_points(both_signs=False)) for e in cremona_curves([11..100])]  # long time (15s on sage.math, 2011)
            [2, 0, 2, 3, 2, 1, 3, 0, 2, 4, 2, 4, 3, 0, 0, 1, 2, 1, 2, 0, 2, 1, 0, 1, 3, 3, 1, 1, 4, 2, 3, 2, 0, 0, 5, 3, 2, 2, 1, 1, 1, 0, 1, 3, 0, 1, 0, 1, 1, 3, 6, 1, 2, 2, 2, 0, 0, 2, 3, 1, 2, 2, 1, 1, 0, 3, 2, 1, 0, 1, 0, 1, 3, 3, 1, 1, 5, 1, 0, 1, 1, 0, 1, 2, 0, 2, 0, 1, 1, 3, 1, 2, 2, 4, 4, 2, 1, 0, 0, 5, 1, 0, 1, 2, 0, 2, 2, 0, 0, 0, 1, 0, 3, 1, 5, 1, 2, 4, 1, 0, 1, 0, 1, 0, 1, 0, 2, 2, 0, 0, 1, 0, 1, 1, 4, 1, 0, 1, 1, 0, 4, 2, 0, 1, 1, 2, 3, 1, 1, 1, 1, 6, 2, 1, 1, 0, 2, 0, 6, 2, 0, 4, 2, 2, 0, 0, 1, 2, 0, 2, 1, 0, 3, 1, 2, 1, 4, 6, 3, 2, 1, 0, 2, 2, 0, 0, 5, 4, 1, 0, 0, 1, 0, 2, 2, 0, 0, 2, 3, 1, 3, 1, 1, 0, 1, 0, 0, 1, 2, 2, 0, 2, 0, 0, 1, 2, 0, 0, 4, 1, 0, 1, 1, 0, 1, 2, 0, 1, 4, 3, 1, 2, 2, 1, 1, 1, 1, 6, 3, 3, 3, 3, 1, 1, 1, 1, 1, 0, 7, 3, 0, 1, 3, 2, 1, 0, 3, 2, 1, 0, 2, 2, 6, 0, 0, 6, 2, 2, 3, 3, 5, 5, 1, 0, 6, 1, 0, 3, 1, 1, 2, 3, 1, 2, 1, 1, 0, 1, 0, 1, 0, 5, 5, 2, 2, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1]

        The bug reported at #4897 is now fixed::

            sage: [P[0] for P in EllipticCurve([0,0,0,-468,2592]).integral_points()]
            [-24, -18, -14, -6, -3, 4, 6, 18, 21, 24, 36, 46, 102, 168, 186, 381, 1476, 2034, 67246]

        .. note::

           This function uses the algorithm given in [Co1].

        REFERENCES:

        - [Co1] Cohen H., Number Theory Vol I: Tools and Diophantine
          Equations GTM 239, Springer 2007

        AUTHORS:

        - Michael Mardaus (2008-07)

        - Tobias Nagel (2008-07)

        - John Cremona (2008-07)
        """
        #####################################################################
        # INPUT CHECK #######################################################
        if not self.is_integral():
            raise ValueError, "integral_points() can only be called on an integral model"

        if mw_base=='auto':
            mw_base = self.gens()
            r = len(mw_base)
        else:
            try:
                r = len(mw_base)
            except TypeError:
                raise TypeError, 'mw_base must be a list'
            if not all([P.curve() is self for P in mw_base]):
                raise ValueError, "points are not on the correct curve"

        tors_points = self.torsion_points()

        # END INPUT-CHECK####################################################
        #####################################################################

        #####################################################################
        # INTERNAL FUNCTIONS ################################################

        ############################## begin ################################
        def point_preprocessing(free,tor):
            r"""
            Transforms the mw_basis ``free`` into a `\ZZ`-basis for
            `E(\QQ)\cap E^0(`\RR)`. If there is a torsion point on the
            "egg" we add it to any of the gens on the egg; otherwise
            we replace the free generators with generators of a
            subgroup of index 2.
            """
            r = len(free)
            newfree = [Q for Q in free] # copy
            tor_egg = [T for T in tor if not T.is_on_identity_component()]
            free_id = [P.is_on_identity_component() for P in free]
            if any(tor_egg):
                T = tor_egg[0]
                for i in range(r):
                    if not free_id[i]:
                        newfree[i] += T
            else:
                if not all(free_id):
                    i0 = free_id.index(False)
                    P = free[i0]
                    for i in range(r):
                        if not free_id[i]:
                            if i==i0:
                                newfree[i] = 2*newfree[i]
                            else:
                                newfree[i] += P
            return newfree
        ##############################  end  ################################

        # END Internal functions #############################################
        ######################################################################

        if (r == 0):
            int_points = [P for P in tors_points if not P.is_zero()]
            int_points = [P for P in int_points if P[0].is_integral()]
            if not both_signs:
                xlist = set([P[0] for P in int_points])
                int_points = [self.lift_x(x) for x in xlist]
            int_points.sort()
            if verbose:
                print 'Total number of integral points:',len(int_points)
            return int_points

        if verbose:
            import sys  # so we can flush stdout for debugging

        g2 = self.c4()/12
        g3 = self.c6()/216
        disc = self.discriminant()
        j = self.j_invariant()
        b2 = self.b2()

        Qx = rings.PolynomialRing(RationalField(),'x')
        pol = Qx([-self.c6()/216,-self.c4()/12,0,4])
        if disc > 0: # two real component -> 3 roots in RR
            #on curve 897e4, only one root is found with default precision!
            RR = R
            prec = RR.precision()
            ei = pol.roots(RR,multiplicities=False)
            while len(ei)<3:
                prec*=2
                RR=RealField(prec)
                ei = pol.roots(RR,multiplicities=False)
            e1,e2,e3 = ei
            if r >= 1: #preprocessing of mw_base only necessary if rank > 0
                mw_base = point_preprocessing(mw_base, tors_points)
                  #at most one point in E^{egg}

        elif disc < 0: # one real component => 1 root in RR (=: e3),
                       # 2 roots in C (e1,e2)
            roots = pol.roots(C,multiplicities=False)
            e3 = pol.roots(R,multiplicities=False)[0]
            roots.remove(e3)
            e1,e2 = roots

        from sage.all import pi
        e = R(1).exp()
        pi = R(pi)

        M = self.height_pairing_matrix(mw_base)
        mw_base, U = self.lll_reduce(mw_base,M)
        M = U.transpose()*M*U

        if verbose:
            print "Using mw_basis ",mw_base
            print "e1,e2,e3: ",e1,e2,e3
            sys.stdout.flush()

        # Algorithm presented in [Co1]
        h_E = self.height()
        w1, w2 = self.period_lattice().basis()
        mu = R(disc).abs().log() / 6
        if j!=0:
            mu += max(R(1),R(j).abs().log()) / 6
        if b2!=0:
            mu += max(R(1),R(b2).abs().log())
            mu += log(R(2))
        else:
            mu += 1

        c1 = (mu + 2.14).exp()
        c2 = min(M.charpoly ().roots(multiplicities=False))
        if verbose:
            print "Minimal eigenvalue of height pairing matrix: ", c2
            sys.stdout.flush()

        c3 = (w1**2)*R(b2).abs()/48 + 8
        c5 = (c1*c3).sqrt()
        c7 = abs((3*pi)/((w1**2) * (w1/w2).imag()))

        mw_base_log = [] #contains \Phi(Q_i)
        mod_h_list = []  #contains h_m(Q_i)
        c9_help_list = []
        for i in range(0,r):
            mw_base_log.append(mw_base[i].elliptic_logarithm().abs())
            mod_h_list.append(max(mw_base[i].height(),h_E,c7*mw_base_log[i]**2))
            c9_help_list.append((mod_h_list[i]).sqrt()/mw_base_log[i])
        c8 = max(e*h_E,max(mod_h_list))
        c9 = e/c7.sqrt() * min(c9_help_list)
        n=r+1
        c10 = R(2 * 10**(8+7*n) * R((2/e)**(2 * n**2)) * (n+1)**(4 * n**2 + 10 * n) * log(c9)**(-2*n - 1) * misc.prod(mod_h_list))

        top = Z(128) #arbitrary first upper bound
        bottom = Z(0)
        log_c9=log(c9); log_c5=log(c5)
        log_r_top = log(R(r*(10**top)))
#        if verbose:
#            print "[bottom,top] = ",[bottom,top]

        while R(c10*(log_r_top+log_c9)*(log(log_r_top)+h_E+log_c9)**(n+1)) > R(c2/2 * (10**top)**2 - log_c5):
            #initial bound 'top' too small, upshift of search interval
            bottom = top
            top = 2*top
        while top >= bottom: #binary-search like search for fitting exponent (bound)
#            if verbose:
#                print "[bottom,top] = ",[bottom,top]
            bound = (bottom + (top - bottom)/2).floor()
            log_r_bound = log(R(r*(10**bound)))
            if R(c10*(log_r_bound+log_c9)*(log(log_r_bound)+h_E+log_c9)**(n+1)) > R(c2/2 * (10**bound)**2 - log_c5):
                bottom = bound + 1
            else:
                top = bound - 1

        H_q = R(10)**bound
        break_cond = 0 #at least one reduction step
        #reduction via LLL
        M = matrix.MatrixSpace(Z,n)
        while break_cond < 0.9: #as long as the improvement of the new bound in comparison to the old is greater than 10%
            c = R((H_q**n)*10)  #c has to be greater than H_q^n
            m = copy(M.identity_matrix())
            for i in range(r):
                m[i, r] = R(c*mw_base_log[i]).round()
            m[r,r] = max(Z(1),R(c*w1).round()) #ensures that m isn't singular

            #LLL - implemented in sage - operates on rows not on columns
            m_LLL = m.LLL()
            m_gram = m_LLL.gram_schmidt()[0]
            b1_norm = R(m_LLL.row(0).norm())

            #compute constant c1 ~ c1_LLL of Corollary 2.3.17 and hence d(L,0)^2 ~ d_L_0
            c1_LLL = -1
            for i in range(n):
                tmp = R(b1_norm/(m_gram.row(i).norm()))
                if tmp > c1_LLL:
                    c1_LLL = tmp

            if c1_LLL < 0:
                raise RuntimeError, 'Unexpected intermediate result. Please try another Mordell-Weil base'

            d_L_0 = R(b1_norm**2 / c1_LLL)

            #Reducing of upper bound
            Q = r * H_q**2
            T = (1 + (3/2*r*H_q))/2
            if d_L_0 < R(T**2+Q):
                d_L_0 = 10*(T**2*Q)
            low_bound = R(((d_L_0 - Q).sqrt() - T)/c)

            #new bound according to low_bound and upper bound
            #[c_5 exp((-c_2*H_q^2)/2)] provided by Corollary 8.7.3
            if low_bound != 0:
                H_q_new = R((log(low_bound/c5)/(-c2/2))).sqrt()
                H_q_new = H_q_new.ceil()
                if H_q_new == 1:
                    break_cond = 1 # stops reduction
                else:
                    break_cond = R(H_q_new/H_q)
                H_q = H_q_new
            else:
                break_cond = 1 # stops reduction, so using last H_q > 0
            #END LLL-Reduction loop

        b2_12 = b2/12
        if disc > 0:
            ##Points in egg have X(P) between e1 and e2 [X(P)=x(P)+b2/12]:
            x_int_points = self.integral_x_coords_in_interval((e1-b2_12).ceil(), (e2-b2_12).floor())
            if verbose:
                print 'x-coords of points on compact component with ',(e1-b2_12).ceil(),'<=x<=',(e2-b2_12).floor()
                L = list(x_int_points) # to have the order
                L.sort()               # deterministic for doctests!
                print L
                sys.stdout.flush()
        else:
            x_int_points = set()

        ##Points in noncompact component with X(P)< 2*max(|e1|,|e2|,|e3|) , espec. X(P)>=e3
        x0 = (e3-b2_12).ceil()
        x1 = (2*max(abs(e1),abs(e2),abs(e3)) - b2_12).ceil()
        x_int_points2 = self.integral_x_coords_in_interval(x0, x1)
        x_int_points = x_int_points.union(x_int_points2)
        if verbose:
            print 'x-coords of points on non-compact component with ',x0,'<=x<=',x1-1
            L = list(x_int_points2)
            L.sort()
            print L
            sys.stdout.flush()

        if verbose:
            print 'starting search of remaining points using coefficient bound ',H_q
            sys.stdout.flush()
        x_int_points3 = integral_points_with_bounded_mw_coeffs(self,mw_base,H_q)
        x_int_points = x_int_points.union(x_int_points3)
        if verbose:
            print 'x-coords of extra integral points:'
            L = list(x_int_points3)
            L.sort()
            print L
            sys.stdout.flush()

        if len(tors_points)>1:
            x_int_points_t = set()
            for x in x_int_points:
                P = self.lift_x(x)
                for T in tors_points:
                    Q = P+T
                    if not Q.is_zero() and Q[0].is_integral():
                        x_int_points_t = x_int_points_t.union([Q[0]])
            x_int_points = x_int_points.union(x_int_points_t)

        # Now we have all the x-coordinates of integral points, and we
        # construct the points, depending on the parameter both_signs:
        if both_signs:
            int_points = sum([self.lift_x(x,all=True) for x in x_int_points],[])
        else:
            int_points = [self.lift_x(x) for x in x_int_points]
        int_points.sort()
        if verbose:
            print 'Total number of integral points:',len(int_points)
        return int_points

    def S_integral_points(self, S, mw_base='auto', both_signs=False, verbose=False, proof=None):
        """
        Computes all S-integral points (up to sign) on this elliptic curve.

        INPUT:

        - ``S`` -  list of primes

        - ``mw_base`` - list of EllipticCurvePoint generating the
          Mordell-Weil group of E (default: 'auto' - calls
          :meth:`.gens`)

        - ``both_signs`` - True/False (default False): if True the
          output contains both P and -P, otherwise only one of each
          pair.

        - ``verbose`` - True/False (default False): if True, some
          details of the computation are output.

        - ``proof`` - True/False (default True): if True ALL
          S-integral points will be returned.  If False, the MW basis
          will be computed with the proof=False flag, and also the
          time-consuming final call to
          S_integral_x_coords_with_abs_bounded_by(abs_bound) is
          omitted.  Use this only if the computation takes too long,
          but be warned that then it cannot be guaranteed that all
          S-integral points will be found.

        OUTPUT:

        A sorted list of all the S-integral points on E (up to sign
        unless both_signs is True)

        .. note::

           The complexity increases exponentially in the rank of curve
           E and in the length of S.  The computation time (but not
           the output!) depends on the Mordell-Weil basis.  If mw_base
           is given but is not a basis for the Mordell-Weil group
           (modulo torsion), S-integral points which are not in the
           subgroup generated by the given points will almost
           certainly not be listed.

        EXAMPLES:

        A curve of rank 3 with no torsion points::

            sage: E=EllipticCurve([0,0,1,-7,6])
            sage: P1=E.point((2,0)); P2=E.point((-1,3)); P3=E.point((4,6))
            sage: a=E.S_integral_points(S=[2,3], mw_base=[P1,P2,P3], verbose=True);a
            max_S: 3 len_S: 3 len_tors: 1
            lambda 0.485997517468...
            k1,k2,k3,k4 6.68597129142710e234 1.31952866480763 3.31908110593519e9 2.42767548272846e17
            p= 2 : trying with p_prec =  30
            mw_base_p_log_val =  [2, 2, 1]
            min_psi =  2 + 2^3 + 2^6 + 2^7 + 2^8 + 2^9 + 2^11 + 2^12 + 2^13 + 2^16 + 2^17 + 2^19 + 2^20 + 2^21 + 2^23 + 2^24 + 2^28 + O(2^30)
            p= 3 : trying with p_prec =  30
            mw_base_p_log_val =  [1, 2, 1]
            min_psi =  3 + 3^2 + 2*3^3 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + 2*3^11 + 2*3^12 + 2*3^13 + 3^15 + 2*3^16 + 3^18 + 2*3^19 + 2*3^22 + 2*3^23 + 2*3^24 + 2*3^27 + 3^28 + 3^29 + O(3^30)
            mw_base [(1 : -1 : 1), (2 : 0 : 1), (0 : -3 : 1)]
            mw_base_log [0.667789378224099, 0.552642660712417, 0.818477222895703]
            mp [5, 7]
            mw_base_p_log [[2^2 + 2^3 + 2^6 + 2^7 + 2^8 + 2^9 + 2^14 + 2^15 + 2^18 + 2^19 + 2^24 + 2^29 + O(2^30), 2^2 + 2^3 + 2^5 + 2^6 + 2^9 + 2^11 + 2^12 + 2^14 + 2^15 + 2^16 + 2^18 + 2^20 + 2^22 + 2^23 + 2^26 + 2^27 + 2^29 + O(2^30), 2 + 2^3 + 2^6 + 2^7 + 2^8 + 2^9 + 2^11 + 2^12 + 2^13 + 2^16 + 2^17 + 2^19 + 2^20 + 2^21 + 2^23 + 2^24 + 2^28 + O(2^30)], [2*3^2 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^9 + 2*3^10 + 3^12 + 2*3^14 + 3^15 + 3^17 + 2*3^19 + 2*3^23 + 3^25 + 3^28 + O(3^30), 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 2*3^7 + 2*3^8 + 3^10 + 2*3^12 + 3^13 + 2*3^14 + 3^15 + 3^18 + 3^22 + 3^25 + 2*3^26 + 3^27 + 3^28 + O(3^30), 3 + 3^2 + 2*3^3 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + 2*3^11 + 2*3^12 + 2*3^13 + 3^15 + 2*3^16 + 3^18 + 2*3^19 + 2*3^22 + 2*3^23 + 2*3^24 + 2*3^27 + 3^28 + 3^29 + O(3^30)]]
            k5,k6,k7 0.321154513240... 1.55246328915... 0.161999172489...
            initial bound 2.6227097483365...e117
            bound_list [58, 58, 58]
            bound_list [8, 9, 9]
            bound_list [8, 7, 7]
            bound_list [8, 7, 7]
            starting search of points using coefficient bound  8
            x-coords of S-integral points via linear combination of mw_base and torsion:
            [-3, -26/9, -8159/2916, -2759/1024, -151/64, -1343/576, -2, -7/4, -1, -47/256, 0, 1/4, 4/9, 9/16, 58/81, 7/9, 6169/6561, 1, 17/16, 2, 33/16, 172/81, 9/4, 25/9, 3, 31/9, 4, 25/4, 1793/256, 8, 625/64, 11, 14, 21, 37, 52, 6142/81, 93, 4537/36, 342, 406, 816, 207331217/4096]
            starting search of extra S-integer points with absolute value bounded by 3.89321964979420
            x-coords of points with bounded absolute value
            [-3, -2, -1, 0, 1, 2]
            Total number of S-integral points: 43
            [(-3 : 0 : 1), (-26/9 : 28/27 : 1), (-8159/2916 : 233461/157464 : 1), (-2759/1024 : 60819/32768 : 1), (-151/64 : 1333/512 : 1), (-1343/576 : 36575/13824 : 1), (-2 : 3 : 1), (-7/4 : 25/8 : 1), (-1 : 3 : 1), (-47/256 : 9191/4096 : 1), (0 : 2 : 1), (1/4 : 13/8 : 1), (4/9 : 35/27 : 1), (9/16 : 69/64 : 1), (58/81 : 559/729 : 1), (7/9 : 17/27 : 1), (6169/6561 : 109871/531441 : 1), (1 : 0 : 1), (17/16 : -25/64 : 1), (2 : 0 : 1), (33/16 : 17/64 : 1), (172/81 : 350/729 : 1), (9/4 : 7/8 : 1), (25/9 : 64/27 : 1), (3 : 3 : 1), (31/9 : 116/27 : 1), (4 : 6 : 1), (25/4 : 111/8 : 1), (1793/256 : 68991/4096 : 1), (8 : 21 : 1), (625/64 : 14839/512 : 1), (11 : 35 : 1), (14 : 51 : 1), (21 : 95 : 1), (37 : 224 : 1), (52 : 374 : 1), (6142/81 : 480700/729 : 1), (93 : 896 : 1), (4537/36 : 305425/216 : 1), (342 : 6324 : 1), (406 : 8180 : 1), (816 : 23309 : 1), (207331217/4096 : 2985362173625/262144 : 1)]

        It is not necessary to specify mw_base; if it is not provided,
        then the Mordell-Weil basis must be computed, which may take
        much longer.

        ::

            sage: a = E.S_integral_points([2,3])
            sage: len(a)
            43

        An example with negative discriminant::

            sage: EllipticCurve('900d1').S_integral_points([17], both_signs=True)
            [(-11 : -27 : 1), (-11 : 27 : 1), (-4 : -34 : 1), (-4 : 34 : 1), (4 : -18 : 1), (4 : 18 : 1), (2636/289 : -98786/4913 : 1), (2636/289 : 98786/4913 : 1), (16 : -54 : 1), (16 : 54 : 1)]

        Output checked with Magma (corrected in 3 cases)::

            sage: [len(e.S_integral_points([2], both_signs=False)) for e in cremona_curves([11..100])] # long time (17s on sage.math, 2011)
            [2, 0, 2, 3, 3, 1, 3, 1, 3, 5, 3, 5, 4, 1, 1, 2, 2, 2, 3, 1, 2, 1, 0, 1, 3, 3, 1, 1, 5, 3, 4, 2, 1, 1, 5, 3, 2, 2, 1, 1, 1, 0, 1, 3, 0, 1, 0, 1, 1, 3, 7, 1, 3, 3, 3, 1, 1, 2, 3, 1, 2, 3, 1, 2, 1, 3, 3, 1, 1, 1, 0, 1, 3, 3, 1, 1, 7, 1, 0, 1, 1, 0, 1, 2, 0, 3, 1, 2, 1, 3, 1, 2, 2, 4, 5, 3, 2, 1, 1, 6, 1, 0, 1, 3, 1, 3, 3, 1, 1, 1, 1, 1, 3, 1, 5, 1, 2, 4, 1, 1, 1, 1, 1, 0, 1, 0, 2, 2, 0, 0, 1, 0, 1, 1, 6, 1, 0, 1, 1, 0, 4, 3, 1, 2, 1, 2, 3, 1, 1, 1, 1, 8, 3, 1, 2, 1, 2, 0, 8, 2, 0, 6, 2, 3, 1, 1, 1, 3, 1, 3, 2, 1, 3, 1, 2, 1, 6, 9, 3, 3, 1, 1, 2, 3, 1, 1, 5, 5, 1, 1, 0, 1, 1, 2, 3, 1, 1, 2, 3, 1, 3, 1, 1, 1, 1, 0, 0, 1, 3, 3, 1, 3, 1, 1, 2, 2, 0, 0, 6, 1, 0, 1, 1, 1, 1, 3, 1, 2, 6, 3, 1, 2, 2, 1, 1, 1, 1, 7, 5, 4, 3, 3, 1, 1, 1, 1, 1, 1, 8, 5, 1, 1, 3, 3, 1, 1, 3, 3, 1, 1, 2, 3, 6, 1, 1, 7, 3, 3, 4, 5, 9, 6, 1, 0, 7, 1, 1, 3, 1, 1, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 7, 8, 2, 3, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1]

        An example from [PZGH]::

            sage: E = EllipticCurve([0,0,0,-172,505])
            sage: E.rank(), len(E.S_integral_points([3,5,7]))  # long time (5s on sage.math, 2011)
            (4, 72)

        This is curve "7690e1" which failed until \#4805 was fixed::

            sage: EllipticCurve([1,1,1,-301,-1821]).S_integral_points([13,2])
            [(-13 : 16 : 1),
            (-9 : 20 : 1),
            (-7 : 4 : 1),
            (21 : 30 : 1),
            (23 : 52 : 1),
            (63 : 452 : 1),
            (71 : 548 : 1),
            (87 : 756 : 1),
            (2711 : 139828 : 1),
            (7323 : 623052 : 1),
            (17687 : 2343476 : 1)]

        REFERENCES:

        - [PZGH] Petho A., Zimmer H.G., Gebel J. and Herrmann E.,
          Computing all S-integral points on elliptic curves
          Math. Proc. Camb. Phil. Soc. (1999), 127, 383-402

        - Some parts of this implementation are partially based on the
          function integral_points()

        AUTHORS:

        - Tobias Nagel (2008-12)

        - Michael Mardaus (2008-12)

        - John Cremona (2008-12)
        """
        # INPUT CHECK #######################################################

        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)


        if not self.is_integral():
            raise ValueError, "S_integral_points() can only be called on an integral model"
        if not all([self.is_p_minimal(s) for s in S]):
            raise ValueError, "%s must be p-minimal for all primes in S"%self

        try:
            len_S = len(S)
            if len_S == 0:
                return self.integral_points(mw_base, both_signs, verbose)
            if not all([s.is_prime() for s in S]):
                raise ValueError, "All elements of S must be prime"
            S.sort()
        except TypeError:
            raise TypeError, 'S must be a list of primes'
        except AttributeError:#catches: <tuple>.sort(), <!ZZ>.is_prime()
            raise AttributeError,'S must be a list of primes'

        if mw_base=='auto':
            if verbose:
                print "Starting computation of MW basis"
            mw_base = self.gens(proof=proof)
            r = len(mw_base)
            if verbose:
                print "Finished computation of MW basis; rank is ",r
        else:
            try:
                r = len(mw_base)
            except TypeError:
                raise TypeError, 'mw_base must be a list'
            if not all([P.curve() is self for P in mw_base]):
                raise ValueError, "mw_base-points are not on the correct curve"

        #End Input-Check ######################################################

        #Internal functions ###################################################
        def reduction_at(p):
            r"""
            Reducing the bound `H_q` at the finite place p in S via LLL
            """
            indexp = S.index(p)
            pc = Z(p**(R(c.log()/log(p,e)).ceil()))
            m = copy(M.identity_matrix())
            for i in range(r):
                try:
                    m[i, r] = Z((beta[indexp][i])%pc)
                except ZeroDivisionError:  #If Inverse doesn't exist, change denominator (which is only approx)
                    val_nu = (beta[indexp][i]).numerator()
                    val_de = (beta[indexp][i]).denominator()
                    m[i, r] = Z((val_nu/(val_de+1))%pc)
            m[r,r] = max(Z(1), pc)

            #LLL - implemented in sage - operates on rows not on columns
            m_LLL = m.LLL()
            m_gram = m_LLL.gram_schmidt()[0]
            b1_norm = R(m_LLL.row(0).norm())

            c1_LLL = -1
            for i in range(n):
                tmp = R(b1_norm/(m_gram.row(i).norm()))
                if tmp > c1_LLL:
                    c1_LLL = tmp
            if c1_LLL < 0:
                raise RuntimeError, 'Unexpected intermediate result. Please try another Mordell-Weil base'
            d_L_0 = R(b1_norm**2 / c1_LLL)

            #Reducing of upper bound
            Q = r * H_q**2
            T = (1 + (3/2*r*H_q))/2
            if d_L_0 < R(T**2+Q):
                d_L_0 = 10*(T**2*Q)
            low_bound = R(((d_L_0 - Q).sqrt() - T)/c)

            ##new bound according to low_bound and upper bound
            ##[k5*k6 exp(-k7**H_q^2)]
            if low_bound != 0:
                H_q_infinity = R(((low_bound/(k6)).log()/(-k7)).sqrt())
                return (H_q_infinity.ceil())
            else:
                return (H_q)
    #<-------------------------------------------------------------------------
    #>-------------------------------------------------------------------------
        def S_integral_points_with_bounded_mw_coeffs():
            r"""
            Returns the set of S-integers x which are x-coordinates of
            points on the curve which are linear combinations of the
            generators (basis and torsion points) with coefficients
            bounded by `H_q`.  The bound `H_q` will be computed at
            runtime.
            (Modified version of integral_points_with_bounded_mw_coeffs() in
             integral_points() )

            TODO: Make this more efficient.  In the case ``S=[]`` we
            worked over the reals and only computed a combination
            exactly if the real coordinates were approximately
            integral.  We need a version of this which works for
            S-integral points, probably by finding a bound on the
            denominator.
            """
            from sage.groups.generic import multiples
            xs=set()
            N=H_q

            def test(P):
                """
                Record x-coord of a point if S-integral.
                """
                if not P.is_zero():
                    xP = P[0]
                    if xP.is_S_integral(S):
                        xs.add(xP)

            def test_with_T(R):
                """
                Record x-coords of a 'point+torsion' if S-integral.
                """
                for T in tors_points:
                    test(R+T)

         # For small rank and small H_q perform simple search
            if r==1 and N<=10:
                for P in multiples(mw_base[0],N+1):
                    test_with_T(P)
                return xs

         # explicit computation and testing linear combinations
         # ni loops through all tuples (n_1,...,n_r) with |n_i| <= N
         # stops when (0,0,...,0) is reached because after that, only inverse points of
         # previously tested points would be tested

            E0=E(0)
            ni = [-N for i in range(r)]
            mw_baseN = [-N*P for P in mw_base]
            Pi = [0 for j in range(r)]
            Pi[0] = mw_baseN[0]
            for i in range(1,r):
                Pi[i] = Pi[i-1] + mw_baseN[i]

            while True:
                if all([n==0 for n in ni]):
                    test_with_T(E0)
                    break

                # test the ni-combination which is Pi[r-1]
                test_with_T(Pi[r-1])

                # increment indices and stored points
                i0 = r-1
                while ni[i0]==N:
                    ni[i0] = -N
                    i0 -= 1
                ni[i0] += 1
                if all([n==0 for n in ni[0:i0+1]]):
                    Pi[i0] = E0
                else:
                    Pi[i0] += mw_base[i0]
                for i in range(i0+1,r):
                    Pi[i] = Pi[i-1] + mw_baseN[i]

            return xs
    #<-------------------------------------------------------------------------
    #>-------------------------------------------------------------------------
        def S_integral_x_coords_with_abs_bounded_by(abs_bound):
            r"""
            Extra search of points with `|x|< ` abs_bound, assuming
            that `x` is `S`-integral and `|x|\ge|x|_q` for all primes
            `q` in `S`. (Such points are not covered by the main part
            of the code).  We know

            .. math::

               x=\frac{\xi}{\p_1^{\alpha_1} \cdot \dots \cdot \p_s^{\alpha_s}},\ (gcd(\xi,\p_i)=1),\ p_i \in S

            so a bound of `\alpha_i` can be found in terms of
            abs_bound. Additionally each `\alpha` must be even, giving
            another restriction.  This gives a finite list of
            denominators to test, and for each, a range of numerators.
            All candidates for `x` resulting from this theory are then
            tested, and a list of the ones which are `x`-coordinates
            of (`S`-integral) points is returned.

            TODO: Make this more efficient.  If we had an efficient
            function for searching for integral points (for example,
            by wrapping Stoll's ratpoint program) then it should be
            better to scale the equation by the maximum denominator
            and search for integral points on the scaled model.

            """
            x_min = min(self.two_division_polynomial().roots(R,multiplicities=False))
            x_min_neg = bool(x_min<0)
            x_min_pos = not x_min_neg
            log_ab = R(abs_bound.log())
            alpha = [(log_ab/R(log(p,e))).floor() for p in S]
            if all([alpha_i <= 1 for alpha_i in alpha]): # so alpha_i must be 0 to satisfy that denominator is a square
                return set([x for x  in range(-abs_bound,abs_bound) if E.is_x_coord(x)])
            else:
                xs = []
                alpha_max_even = [y-y%2 for y in alpha]
                p_pow_alpha = []
                list_alpha = []
                for i in range(len_S-1):
                    list_alpha.append(range(0,alpha_max_even[i]+2,2))
                    p_pow_alpha.append([S[i]**list_alpha[i][j] for j in range(len(list_alpha[i]))])
                if verbose:
                    print list_alpha, p_pow_alpha
                # denom_maxpa is a list of pairs (d,q) where d runs
                # through possible denominators, and q=p^a is the
                # maximum prime power divisor of d:
                denom_maxpa = [(misc.prod(tmp),max(tmp)) for tmp in cartesian_product_iterator(p_pow_alpha)]
#               The maximum denominator is this (not used):
#                denom = [misc.prod([pp[-1] for pp in p_pow_alpha],1)]
                for de,maxpa in denom_maxpa:
                    n_max = (abs_bound*de).ceil()
                    n_min = maxpa*de
                    if x_min_pos:
                        pos_n_only = True
                        if x_min > maxpa:
                            n_min = (x_min*de).floor()
                    else:
                        pos_n_only = False
                        neg_n_max = (x_min.abs()*de).ceil()

#                   if verbose:
#                       print "testing denominator ",de
#                       print "numerator bounds = ",(n_min,n_max)

                    for n in misc.xsrange(n_min,n_max+1):
                        tmp = n/de  # to save time, do not check de is the exact denominator
                        if E.is_x_coord(tmp):
                            xs+=[tmp]
                        if not pos_n_only:
                            if n <= neg_n_max:
                                if E.is_x_coord(-tmp):
                                    xs+=[-tmp]

                return set(xs)
    #<-------------------------------------------------------------------------
        #End internal functions ###############################################
        from sage.misc.all import cartesian_product_iterator

        E = self
        tors_points = E.torsion_points()

        if (r==0):#only Torsionpoints to consider
            int_points = [P for P in tors_points if not P.is_zero()]
            int_points = [P for P in int_points if P[0].is_S_integral(S)]
            if not both_signs:
                xlist = set([P[0] for P in int_points])
                int_points = [E.lift_x(x) for x in xlist]
            int_points.sort()
            if verbose:
                print 'Total number of S-integral points:',len(int_points)
            return int_points

        if verbose:
            import sys  # so we can flush stdout for debugging

        e = R(1).exp()
        a1, a2, a3, a4, a6 = E.a_invariants()
        b2, b4, b6, b8 = E.b_invariants()
        c4, c6 = E.c_invariants()
        disc = E.discriminant()
        #internal function is doing only a comparison of E and E.short_weierstass_model() so the following is easier
        if a1 == a2 == a3 == 0:
            is_short = True
        else:
            is_short = False

        w1, w2 = E.period_lattice().basis()

        Qx = rings.PolynomialRing(RationalField(),'x')
        pol = Qx([-54*c6,-27*c4,0,1])
        if disc > 0: # two real component -> 3 roots in RR
            # it is possible that only one root is found with default precision! (see integral_points())
            RR = R
            prec = RR.precision()
            ei = pol.roots(RR,multiplicities=False)
            while len(ei)<3:
                prec*=2
                RR=RealField(prec)
                ei = pol.roots(RR,multiplicities=False)
            e1,e2,e3 = ei
        elif disc < 0: # one real component => 1 root in RR (=: e3),
                       # 2 roots in C (e1,e2)
            roots = pol.roots(C,multiplicities=False)
            e3 = pol.roots(R,multiplicities=False)[0]
            roots.remove(e3)
            e1,e2 = roots

        len_tors = len(tors_points)
        n = r + 1

        M = E.height_pairing_matrix(mw_base)
        mw_base, U = E.lll_reduce(mw_base,M)
        M = U.transpose()*M*U

        # NB "lambda" is a reserved word in Python!
        lamda = min(M.charpoly(algorithm="hessenberg").roots(multiplicities = False))
        max_S = max(S)
        len_S += 1 #Counting infinity (always "included" in S)
        if verbose:
            print 'max_S:',max_S,'len_S:',len_S,'len_tors:',len_tors
            print 'lambda',lamda
            sys.stdout.flush()

        if is_short:
            disc_0_abs = R((4*a4**3 + 27*a6**2).abs())
            k4 = R(10**4 * max(16*a4**2, 256*disc_0_abs.sqrt()**3))
            k3 = R(32/3 * disc_0_abs.sqrt() * (8 + 0.5*disc_0_abs.log())**4)
        else:
            disc_sh = R(E.short_weierstrass_model().discriminant()) #computes y^2=x^3 -27c4x -54c6
            k4 = R(20**4 * max(3**6 * c4**2, 16*(disc_sh.abs().sqrt())**3))
            k3 = R(32/3 * disc_sh.abs().sqrt() * (8 + 0.5*disc_sh.abs().log())**4)


        k2 = max(R(b2.abs()), R(b4.abs().sqrt()), R(b6.abs()**(1/3)), R(b8.abs()**(1/4))).log()
        k1 = R(7 * 10**(38*len_S+49)) * R(len_S**(20*len_S+15)) * max_S**24 * R(max(1,log(max_S, e))**(4*len_S - 2)) * k3 * k3.log()**2 * ((20*len_S - 19)*k3 + (e*k4).log()) + 2*R(2*b2.abs()+6).log()

        if verbose:
            print 'k1,k2,k3,k4',k1,k2,k3,k4
            sys.stdout.flush()
        #H_q -> [PZGH]:N_0 (due to consistency to integral_points())
        H_q = R(((k1/2+k2)/lamda).sqrt())

        #computation of logs
        mw_base_log = [(pts.elliptic_logarithm().abs())*(len_tors/w1) for pts in mw_base]
        mw_base_p_log = []
        beta = []
        mp=[]
        tmp = 0
        for p in S:
            Np = E.Np(p)
            cp = E.tamagawa_exponent(p)
            mp_temp = Z(len_tors).lcm(cp*Np)
            mp.append(mp_temp) #only necessary because of verbose below
            p_prec=30+E.discriminant().valuation(p)
            p_prec_ok=False
            while not p_prec_ok:
                if verbose:
                    print "p=",p,": trying with p_prec = ",p_prec
                try:
                    mw_base_p_log.append([mp_temp*(pts.padic_elliptic_logarithm(p,absprec=p_prec)) for pts in mw_base])
                    p_prec_ok=True
                except ValueError:
                    p_prec *= 2
            #reorder mw_base_p: last value has minimal valuation at p
            mw_base_p_log_val = [mw_base_p_log[tmp][i].valuation() for i in range(r)]
            if verbose:
                print "mw_base_p_log_val = ",mw_base_p_log_val
            min_index = mw_base_p_log_val.index(min(mw_base_p_log_val))
            min_psi = mw_base_p_log[tmp][min_index]
            if verbose:
                print "min_psi = ",min_psi
            mw_base_p_log[tmp].remove(min_psi)
            mw_base_p_log[tmp].append(min_psi)
            #beta needed for reduction at p later on
            try:
                beta.append([-mw_base_p_log[tmp][j]/min_psi for j in range(r)])
            except ValueError:
                # e.g. mw_base_p_log[tmp]==[0]:  can occur e.g. [?]'172c6, S=[2]
                beta.append([0] for j in range(r))
            tmp +=1

        if verbose:
            print 'mw_base',mw_base
            print 'mw_base_log', mw_base_log
            print 'mp', mp
            print 'mw_base_p_log',mw_base_p_log
            #print 'beta', beta
            sys.stdout.flush()

        #constants in reduction (not needed to be computed every reduction step)
        k5 = R((2*len_tors)/(3*w1))
        k6 = R((k2/len_S).exp())
        k7 = R(lamda/len_S)

        if verbose:
            print 'k5,k6,k7',k5,k6,k7
            sys.stdout.flush()

        break_cond = 0
        M = matrix.MatrixSpace(Z,n)
   #Reduction of initial bound
        if verbose:
            print 'initial bound',H_q
            sys.stdout.flush()

        while break_cond < 0.9:
         #reduction at infinity
            bound_list=[]
            c = R((H_q**n)*100)
            m = copy(M.identity_matrix())
            for i in range(r):
                m[i, r] = R(c*mw_base_log[i]).round()
            m[r,r] = max(Z(1), R(c*w1).round())
            #LLL - implemented in sage - operates on rows not on columns
            m_LLL = m.LLL()
            m_gram = m_LLL.gram_schmidt()[0]
            b1_norm = R(m_LLL.row(0).norm())

            #compute constant c1_LLL (cf. integral_points())
            c1_LLL = -1
            for i in range(n):
                tmp = R(b1_norm/(m_gram.row(i).norm()))
                if tmp > c1_LLL:
                    c1_LLL = tmp
            if c1_LLL < 0:
                raise RuntimeError, 'Unexpected intermediate result. Please try another Mordell-Weil base'
            d_L_0 = R(b1_norm**2 / c1_LLL)

            #Reducing of upper bound
            Q = r * H_q**2
            T = (1 + (3/2*r*H_q))/2
            if d_L_0 < R(T**2+Q):
                d_L_0 = 10*(T**2*Q)
            low_bound = R(((d_L_0 - Q).sqrt() - T)/c)

            ##new bound according to low_bound and upper bound
            ##[k5*k6 exp(-k7**H_q^2)]
            if low_bound != 0:
                H_q_infinity = R(((low_bound/(k5*k6)).log()/(-k7)).abs().sqrt())
                bound_list.append(H_q_infinity.ceil())
            else:
                bound_list.append(H_q)

         ##reduction for finite places in S
            for p in S:
                bound_list.append(reduction_at(p))

            if verbose:
                print 'bound_list',bound_list
                sys.stdout.flush()

            H_q_new = max(bound_list)
            if (H_q_new > H_q): #no improvement
                break_cond = 1 #stop reduction
            elif (H_q_new == 1): #best possible H_q
                H_q = H_q_new
                break_cond = 1 #stop
            else:
                break_cond = R(H_q_new/H_q)
                H_q = H_q_new
    #end of reductions

    #search of S-integral points
        #step1: via linear combination and H_q
        x_S_int_points = set()
        if verbose:
            print 'starting search of points using coefficient bound ',H_q
            sys.stdout.flush()
        x_S_int_points1 = S_integral_points_with_bounded_mw_coeffs()
        x_S_int_points = x_S_int_points.union(x_S_int_points1)
        if verbose:
            print 'x-coords of S-integral points via linear combination of mw_base and torsion:'
            L = list(x_S_int_points1)
            L.sort()
            print L
            sys.stdout.flush()

        #step 2: Extra search
        if e3 < 0:
            M = R( max((27*c4).abs().sqrt(), R((54*c6).abs()**(1/3)) / R(2**(1/3))-1) )
        else:
            M = R(0)
        e0 = max(e1+e2, 2*e3) + M
        abs_bound = R((max(0,e0)+6*b2.abs())/36)

        if proof:
            if verbose:
                print 'starting search of extra S-integer points with absolute value bounded by',abs_bound
                sys.stdout.flush()
            if abs_bound != 0:
                x_S_int_points2 = S_integral_x_coords_with_abs_bounded_by(abs_bound)
                x_S_int_points = x_S_int_points.union(x_S_int_points2)
                if verbose:
                    print 'x-coords of points with bounded absolute value'
                    L = list(x_S_int_points2)
                    L.sort()
                    print L
                    sys.stdout.flush()

        if len(tors_points)>1:
            x_S_int_points_t = set()
            for x in x_S_int_points:
                P = E.lift_x(x)
                for T in tors_points:
                    Q = P+T
                    if not Q.is_zero() and Q[0].is_S_integral(S):
                        x_S_int_points_t = x_S_int_points_t.union([Q[0]])
            x_S_int_points = x_S_int_points.union(x_S_int_points_t)

        # All x values collected, now considering "both_signs"
        if both_signs:
            S_int_points = sum([self.lift_x(x,all=True) for x in x_S_int_points],[])
        else:
            S_int_points = [self.lift_x(x) for x in x_S_int_points]
        S_int_points.sort()
        if verbose:
            print 'Total number of S-integral points:',len(S_int_points)
        return S_int_points


def cremona_curves(conductors):
    """
    Return iterator over all known curves (in database) with conductor
    in the list of conductors.

    EXAMPLES::

        sage: [(E.label(), E.rank()) for E in cremona_curves(srange(35,40))]
        [('35a1', 0),
        ('35a2', 0),
        ('35a3', 0),
        ('36a1', 0),
        ('36a2', 0),
        ('36a3', 0),
        ('36a4', 0),
        ('37a1', 1),
        ('37b1', 0),
        ('37b2', 0),
        ('37b3', 0),
        ('38a1', 0),
        ('38a2', 0),
        ('38a3', 0),
        ('38b1', 0),
        ('38b2', 0),
        ('39a1', 0),
        ('39a2', 0),
        ('39a3', 0),
        ('39a4', 0)]
    """
    if isinstance(conductors, (int,long, rings.RingElement)):
        conductors = [conductors]
    return sage.databases.cremona.CremonaDatabase().iter(conductors)

def cremona_optimal_curves(conductors):
    """
    Return iterator over all known optimal curves (in database) with
    conductor in the list of conductors.

    EXAMPLES::

        sage: [(E.label(), E.rank()) for E in cremona_optimal_curves(srange(35,40))]
        [('35a1', 0),
        ('36a1', 0),
        ('37a1', 1),
        ('37b1', 0),
        ('38a1', 0),
        ('38b1', 0),
        ('39a1', 0)]

    There is one case -- 990h3 -- when the optimal curve isn't labeled with a 1::

        sage: [e.cremona_label() for e in cremona_optimal_curves([990])]
        ['990a1', '990b1', '990c1', '990d1', '990e1', '990f1', '990g1', '990h3', '990i1', '990j1', '990k1', '990l1']

    """
    if isinstance(conductors, (int,long,rings.RingElement)):
        conductors = [conductors]
    return sage.databases.cremona.CremonaDatabase().iter_optimal(conductors)

def integral_points_with_bounded_mw_coeffs(E, mw_base, N):
    r"""
    Returns the set of integers `x` which are
    `x`-coordinates of points on the curve `E` which
    are linear combinations of the generators (basis and torsion
    points) with coefficients bounded by `N`.
    """
    from sage.groups.generic import multiples
    xs=set()
    tors_points = E.torsion_points()
    r = len(mw_base)

    def use(P):
        """
        Helper function to record x-coord of a point if integral.
        """
        if not P.is_zero():
            xP = P[0]
            if xP.is_integral():
                xs.add(xP)

    def use_t(R):
        """
        Helper function to record x-coords of a point +torsion if
        integral.
        """
        for T in tors_points:
            use(R+T)

    # We use a naive method when the number of possibilities is small:

    if r==1 and N<=10:
        for P in multiples(mw_base[0],N+1):
            use_t(P)
        return xs

    # Otherwise it is very very much faster to first compute
    # the linear combinations over RR, and only compute them as
    # rational points if they are approximately integral.

    # Note: making eps larger here will dramatically increase
    # the running time.  If evidence arises that integral
    # points are being missed, it would be better to increase
    # the real precision than to increase eps.

    def is_approx_integral(P):
        r"""
        Local function, returns True if the real point `P` is approximately integral.
        """
        eps = 0.0001
        return (abs(P[0]-P[0].round()))<eps and (abs(P[1]-P[1].round()))<eps

    RR = RealField(100) #(100)
    ER = E.change_ring(RR)
    ER0 = ER(0)

    # Note: doing [ER(P) for P in mw_base] sometimes fails.  The
    # following way is harder, since we have to make sure we don't use
    # -P instead of P, but is safer.

    Rgens = [ER.lift_x(P[0]) for P in mw_base]
    for i in range(r):
        if abs(Rgens[i][1]-mw_base[i][1])>abs((-Rgens[i])[1]-mw_base[i][1]):
            Rgens[i] = -Rgens[i]

    # the ni loop through all tuples (a1,a2,...,ar) with
    # |ai|<=N, but we stop immediately after using the tuple
    # (0,0,...,0).

    # Initialization:
    ni = [-N for i in range(r)]
    RgensN = [-N*P for P in Rgens]
    # RPi[i] = -N*(Rgens[0]+...+Rgens[i])
    RPi = [0 for j in range(r)]
    RPi[0] = RgensN[0]
    for i in range(1,r):
        RPi[i] = RPi[i-1] + RgensN[i]

    tors_points_R = map(ER, tors_points)
    while True:
        if all([n==0 for n in ni]):
             use_t(E(0))
             break

        # test the ni-combination which is RPi[r-1]
        RP = RPi[r-1]

        for T, TR in zip(tors_points, tors_points_R):
            if is_approx_integral(RP + TR):
                 P = sum([ni[i]*mw_base[i] for i in range(r)],T)
                 use(P)

        # increment indices and stored points
        i0 = r-1
        while ni[i0]==N:
            ni[i0] = -N
            i0 -= 1
        ni[i0] += 1
        # The next lines are to prevent rounding error: (-P)+P behaves
        # badly for real points!
        if all([n==0 for n in ni[0:i0+1]]):
            RPi[i0] = ER0
        else:
            RPi[i0] += Rgens[i0]
        for i in range(i0+1,r):
            RPi[i] = RPi[i-1] + RgensN[i]

    return xs

