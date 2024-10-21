# sage_setup: distribution = sagemath-schemes
r"""
Local data for elliptic curves over number fields

Let `E` be an elliptic curve over a number field `K` (including `\QQ`).
There are several local invariants at a finite place `v` that
can be computed via Tate's algorithm (see [Sil1994]_ IV.9.4 or [Tate1975]_).

These include the type of reduction (good, additive, multiplicative),
a minimal equation of `E` over `K_v`,
the Tamagawa number `c_v`, defined to be the index `[E(K_v):E^0(K_v)]`
of the points with good reduction among the local points, and the
exponent of the conductor `f_v`.

The functions in this file will typically be called by using ``local_data``.

EXAMPLES::

    sage: # needs sage.rings.number_field
    sage: x = polygen(ZZ, 'x')
    sage: K.<i> = NumberField(x^2 + 1)
    sage: E = EllipticCurve([(2+i)^2, (2+i)^7])
    sage: pp = K.fractional_ideal(2+i)
    sage: da = E.local_data(pp)
    sage: da.has_bad_reduction()
    True
    sage: da.has_multiplicative_reduction()
    False
    sage: da.kodaira_symbol()
    I0*
    sage: da.tamagawa_number()
    4
    sage: da.minimal_model()
    Elliptic Curve defined by y^2 = x^3 + (4*i+3)*x + (-29*i-278)
     over Number Field in i with defining polynomial x^2 + 1

An example to show how the Neron model can change as one extends the field::

    sage: E = EllipticCurve([0,-1])
    sage: E.local_data(2)
    Local data at Principal ideal (2) of Integer Ring:
      Reduction type: bad additive
      Local minimal model: Elliptic Curve defined by y^2 = x^3 - 1 over Rational Field
      Minimal discriminant valuation: 4
      Conductor exponent: 4
      Kodaira Symbol: II
      Tamagawa Number: 1

    sage: # needs sage.rings.number_field
    sage: EK = E.base_extend(K)
    sage: EK.local_data(1+i)
    Local data at Fractional ideal (i + 1):
      Reduction type: bad additive
      Local minimal model: Elliptic Curve defined by y^2 = x^3 + (-1)
                           over Number Field in i with defining polynomial x^2 + 1
      Minimal discriminant valuation: 8
      Conductor exponent: 2
      Kodaira Symbol: IV*
      Tamagawa Number: 3

Or how the minimal equation changes::

    sage: E = EllipticCurve([0,8])
    sage: E.is_minimal()
    True

    sage: # needs sage.rings.number_field
    sage: EK = E.base_extend(K)
    sage: da = EK.local_data(1+i)
    sage: da.minimal_model()
    Elliptic Curve defined by y^2 = x^3 + (-i)
     over Number Field in i with defining polynomial x^2 + 1

AUTHORS:

- John Cremona: First version 2008-09-21 (refactoring code from
  ``ell_number_field.py`` and ``ell_rational_field.py``)

- Chris Wuthrich: more documentation 2010-01
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.verbose import verbose

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal

from sage.rings.number_field.number_field_base import NumberField
from sage.rings.ideal import Ideal_generic

from .constructor import EllipticCurve
from .kodaira_symbol import KodairaSymbol


class EllipticCurveLocalData(SageObject):
    r"""
    The class for the local reduction data of an elliptic curve.

    Currently supported are elliptic curves defined over `\QQ`, and
    elliptic curves defined over a number field, at an arbitrary prime
    or prime ideal.

    INPUT:

    - ``E`` -- an elliptic curve defined over a number field, or `\QQ`

    - ``P`` -- a prime ideal of the field, or a prime integer if the field is `\QQ`

    - ``proof`` -- boolean; if ``True``, only use provably correct
      methods (default: controlled by global proof module).  Note
      that the proof module is number_field, not elliptic_curves,
      since the functions that actually need the flag are in
      number fields.

    - ``algorithm`` -- string (default: ``'pari'``); ignored unless the
      base field is `\QQ`.  If "pari", use the PARI C-library
      ``ellglobalred`` implementation of Tate's algorithm over
      `\QQ`. If "generic", use the general number field
      implementation.

    .. NOTE::

        This function is not normally called directly by users, who
        may access the data via methods of the EllipticCurve
        classes.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
        sage: E = EllipticCurve('14a1')
        sage: EllipticCurveLocalData(E,2)
        Local data at Principal ideal (2) of Integer Ring:
          Reduction type: bad non-split multiplicative
          Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6
                               over Rational Field
          Minimal discriminant valuation: 6
          Conductor exponent: 1
          Kodaira Symbol: I6
          Tamagawa Number: 2
    """

    def __init__(self, E, P, proof=None, algorithm='pari', globally=False):
        r"""
        Initialize the reduction data for the elliptic curve `E` at the prime `P`.

        INPUT:

        - ``E`` -- an elliptic curve defined over a number field, or `\QQ`

        - ``P`` -- a prime ideal of the field, or a prime integer if the field is `\QQ`

        - ``proof`` -- boolean; if ``True``, only use provably correct
          methods (default: controlled by global proof module).  Note
          that the proof module is number_field, not elliptic_curves,
          since the functions that actually need the flag are in
          number fields.

        - ``algorithm`` -- string (default: ``'pari'``); ignored unless the
          base field is `\QQ`.  If "pari", use the PARI C-library
          ``ellglobalred`` implementation of Tate's algorithm over
          `\QQ`. If "generic", use the general number field
          implementation.

        - ``globally`` -- boolean (default: ``False``); if ``True``, the algorithm
          uses the generators of principal ideals rather than an arbitrary
          uniformizer.

        .. NOTE::

            This function is not normally called directly by users, who
            may access the data via methods of the EllipticCurve
            classes.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve('14a1')
            sage: EllipticCurveLocalData(E, 2)
            Local data at Principal ideal (2) of Integer Ring:
              Reduction type: bad non-split multiplicative
              Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6
                                   over Rational Field
              Minimal discriminant valuation: 6
              Conductor exponent: 1
              Kodaira Symbol: I6
              Tamagawa Number: 2

        ::

            sage: EllipticCurveLocalData(E, 2, algorithm='generic')
            Local data at Principal ideal (2) of Integer Ring:
              Reduction type: bad non-split multiplicative
              Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6
                                   over Rational Field
              Minimal discriminant valuation: 6
              Conductor exponent: 1
              Kodaira Symbol: I6
              Tamagawa Number: 2

        ::

            sage: EllipticCurveLocalData(E, 2, algorithm='pari')
            Local data at Principal ideal (2) of Integer Ring:
              Reduction type: bad non-split multiplicative
              Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6
                                   over Rational Field
              Minimal discriminant valuation: 6
              Conductor exponent: 1
              Kodaira Symbol: I6
              Tamagawa Number: 2

        ::

            sage: EllipticCurveLocalData(E, 2, algorithm='unknown')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'pari', 'generic'

        ::

            sage: EllipticCurveLocalData(E, 3)
            Local data at Principal ideal (3) of Integer Ring:
              Reduction type: good
              Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6
                                   over Rational Field
              Minimal discriminant valuation: 0
              Conductor exponent: 0
              Kodaira Symbol: I0
              Tamagawa Number: 1

        ::

            sage: EllipticCurveLocalData(E, 7)
            Local data at Principal ideal (7) of Integer Ring:
              Reduction type: bad split multiplicative
              Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6
                                   over Rational Field
              Minimal discriminant valuation: 3
              Conductor exponent: 1
              Kodaira Symbol: I3
              Tamagawa Number: 3
        """
        self._curve = E
        K = E.base_field()
        p = check_prime(K,P) # error handling done in that function
        if algorithm != "pari" and algorithm != "generic":
            raise ValueError("algorithm must be one of 'pari', 'generic'")

        self._reduction_type = None
        if K is QQ:
            self._prime = ZZ.ideal(p)
        else:
            self._prime = p

        if algorithm == "pari" and K is QQ:
            Eint = E.integral_model()
            data = Eint.pari_curve().elllocalred(p)
            self._fp = data[0].sage()
            self._KS = KodairaSymbol(data[1].sage())
            self._cp = data[3].sage()
            # We use a global minimal model since we can:
            self._Emin_reduced = Eint.minimal_model()
            self._val_disc = self._Emin_reduced.discriminant().valuation(p)
            if self._fp > 0:
                self._reduction_type = Eint.ap(p) # = 0,-1 or +1
        else:
            self._Emin, _, self._val_disc, self._fp, self._KS, self._cp, self._split = self._tate(proof, globally)
            if self._fp > 0:
                if self._Emin.c4().valuation(p) > 0:
                    self._reduction_type = 0
                elif self._split:
                    self._reduction_type = +1
                else:
                    self._reduction_type = -1

    def __repr__(self):
        r"""
        Return the string representation of this reduction data.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve('14a1')
            sage: EllipticCurveLocalData(E,2).__repr__()
            'Local data at Principal ideal (2) of Integer Ring:\nReduction type: bad non-split multiplicative\nLocal minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field\nMinimal discriminant valuation: 6\nConductor exponent: 1\nKodaira Symbol: I6\nTamagawa Number: 2'
            sage: EllipticCurveLocalData(E,3).__repr__()
            'Local data at Principal ideal (3) of Integer Ring:\nReduction type: good\nLocal minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field\nMinimal discriminant valuation: 0\nConductor exponent: 0\nKodaira Symbol: I0\nTamagawa Number: 1'
            sage: EllipticCurveLocalData(E,7).__repr__()
            'Local data at Principal ideal (7) of Integer Ring:\nReduction type: bad split multiplicative\nLocal minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field\nMinimal discriminant valuation: 3\nConductor exponent: 1\nKodaira Symbol: I3\nTamagawa Number: 3'
        """
        red_type = "good"
        if self._reduction_type is not None:
            red_type = ["bad non-split multiplicative","bad additive","bad split multiplicative"][1+self._reduction_type]
        return "Local data at %s:\nReduction type: %s\nLocal minimal model: %s\nMinimal discriminant valuation: %s\nConductor exponent: %s\nKodaira Symbol: %s\nTamagawa Number: %s" % (self._prime,red_type,self.minimal_model(),self._val_disc,self._fp,self._KS,self._cp)

    def minimal_model(self, reduce=True):
        """
        Return the (local) minimal model from this local reduction data.

        INPUT:

        - ``reduce`` -- boolean (default: ``True``); if set to ``True`` and if
          the initial elliptic curve had globally integral
          coefficients, then the elliptic curve returned by Tate's
          algorithm will be "reduced" as specified in _reduce_model()
          for curves over number fields.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2 = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E, 2)
            sage: data.minimal_model()
            Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
            sage: data.minimal_model() == E.local_minimal_model(2)
            True

        To demonstrate the behaviour of the parameter ``reduce``::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: E = EllipticCurve(K, [0, 0, a, 0, 1])
            sage: E.local_data(K.ideal(a-1)).minimal_model()
            Elliptic Curve defined by y^2 + a*y = x^3 + 1
             over Number Field in a with defining polynomial x^3 + x + 1
            sage: E.local_data(K.ideal(a-1)).minimal_model(reduce=False)
            Elliptic Curve defined by y^2 + (a+2)*y = x^3 + 3*x^2 + 3*x + (-a+1)
             over Number Field in a with defining polynomial x^3 + x + 1

            sage: E = EllipticCurve([2, 1, 0, -2, -1])
            sage: E.local_data(ZZ.ideal(2), algorithm='generic').minimal_model(reduce=False)
            Elliptic Curve defined by y^2 + 2*x*y + 2*y = x^3 + x^2 - 4*x - 2 over Rational Field
            sage: E.local_data(ZZ.ideal(2), algorithm='pari').minimal_model(reduce=False)
            Traceback (most recent call last):
            ...
            ValueError: the argument reduce must not be False if algorithm=pari is used
            sage: E.local_data(ZZ.ideal(2), algorithm='generic').minimal_model()
            Elliptic Curve defined by y^2 = x^3 - x^2 - 3*x + 2 over Rational Field
            sage: E.local_data(ZZ.ideal(2), algorithm='pari').minimal_model()
            Elliptic Curve defined by y^2 = x^3 - x^2 - 3*x + 2 over Rational Field

        :issue:`14476`::

            sage: # needs sage.rings.number_field
            sage: t = QQ['t'].0
            sage: K.<g> = NumberField(t^4 - t^3-3*t^2 - t +1)
            sage: E = EllipticCurve([-2*g^3 + 10/3*g^2 + 3*g - 2/3,
            ....:                    -11/9*g^3 + 34/9*g^2 - 7/3*g + 4/9,
            ....:                    -11/9*g^3 + 34/9*g^2 - 7/3*g + 4/9, 0, 0])
            sage: vv = K.fractional_ideal(g^2 - g - 2)
            sage: E.local_data(vv).minimal_model()
            Elliptic Curve defined by
             y^2 + (-2*g^3+10/3*g^2+3*g-2/3)*x*y + (-11/9*g^3+34/9*g^2-7/3*g+4/9)*y
              = x^3 + (-11/9*g^3+34/9*g^2-7/3*g+4/9)*x^2
             over Number Field in g with defining polynomial t^4 - t^3 - 3*t^2 - t + 1
        """
        if reduce:
            try:
                return self._Emin_reduced
            except AttributeError:
                pass
            # trac 14476 we only reduce if the coefficients are globally integral
            if all(a.is_integral() for a in self._Emin.a_invariants()):
                self._Emin_reduced = self._Emin._reduce_model()
                return self._Emin_reduced
            else:
                return self._Emin
        else:
            try:
                return self._Emin
            except AttributeError:
                raise ValueError("the argument reduce must not be False if algorithm=pari is used")

    def prime(self):
        """
        Return the prime ideal associated with this local reduction data.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2 = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.prime()
            Principal ideal (2) of Integer Ring
        """
        return self._prime

    def conductor_valuation(self):
        """
        Return the valuation of the conductor from this local reduction data.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2 = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.conductor_valuation()
            2
        """
        return self._fp

    def discriminant_valuation(self):
        """
        Return the valuation of the minimal discriminant from this local reduction data.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2 = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.discriminant_valuation()
            4
        """
        return self._val_disc

    def kodaira_symbol(self):
        r"""
        Return the Kodaira symbol from this local reduction data.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2 = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.kodaira_symbol()
            IV
        """
        return self._KS

    def tamagawa_number(self):
        r"""
        Return the Tamagawa number from this local reduction data.

        This is the index `[E(K_v):E^0(K_v)]`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2 = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.tamagawa_number()
            3
        """
        return self._cp

    def tamagawa_exponent(self):
        r"""
        Return the Tamagawa index from this local reduction data.

        This is the exponent of `E(K_v)/E^0(K_v)`; in most cases it is
        the same as the Tamagawa index.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve('816a1')
            sage: data = EllipticCurveLocalData(E, 2)
            sage: data.kodaira_symbol()
            I2*
            sage: data.tamagawa_number()
            4
            sage: data.tamagawa_exponent()
            2

            sage: E = EllipticCurve('200c4')
            sage: data = EllipticCurveLocalData(E, 5)
            sage: data.kodaira_symbol()
            I4*
            sage: data.tamagawa_number()
            4
            sage: data.tamagawa_exponent()
            2
        """
        cp = self._cp
        if cp != 4:
            return cp
        ks = self._KS
        if ks._roman == 1 and ks._n % 2 == 0 and ks._starred:
            return ZZ(2)
        return ZZ(4)

    def bad_reduction_type(self):
        r"""
        Return the type of bad reduction of this reduction data.

        OUTPUT:

        integer or ``None``:

        - +1 for split multiplicative reduction
        - -1 for non-split multiplicative reduction
        - 0  for additive reduction
        - ``None`` for good reduction

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p,E.local_data(p).bad_reduction_type()) for p in prime_range(15)]
            [(2, -1), (3, None), (5, None), (7, 1), (11, None), (13, None)]

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p,E.local_data(p).bad_reduction_type()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), None), (Fractional ideal (2*a + 1), 0)]
        """
        return self._reduction_type

    def has_good_reduction(self):
        r"""
        Return ``True`` if there is good reduction.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_good_reduction()) for p in prime_range(15)]
            [(2, False), (3, True), (5, True), (7, False), (11, True), (13, True)]

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p,E.local_data(p).has_good_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), True),
             (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type is None

    def has_bad_reduction(self):
        r"""
        Return ``True`` if there is bad reduction.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_bad_reduction()) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p,E.local_data(p).has_bad_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), True)]
        """
        return self._reduction_type is not None

    def has_multiplicative_reduction(self):
        r"""
        Return ``True`` if there is multiplicative reduction.

        .. NOTE::

            See also ``has_split_multiplicative_reduction()`` and
            ``has_nonsplit_multiplicative_reduction()``.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.local_data(p).has_multiplicative_reduction()) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p,E.local_data(p).has_multiplicative_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type in (-1,+1)

    def has_split_multiplicative_reduction(self):
        r"""
        Return ``True`` if there is split multiplicative reduction.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.local_data(p).has_split_multiplicative_reduction())
            ....:  for p in prime_range(15)]
            [(2, False), (3, False), (5, False), (7, True), (11, False), (13, False)]

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p,E .local_data(p).has_split_multiplicative_reduction())
            ....:  for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type == +1

    def has_nonsplit_multiplicative_reduction(self):
        r"""
        Return ``True`` if there is non-split multiplicative reduction.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.local_data(p).has_nonsplit_multiplicative_reduction())
            ....:  for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, False), (11, False), (13, False)]

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p, E.local_data(p).has_nonsplit_multiplicative_reduction())
            ....:  for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type == -1

    def has_additive_reduction(self):
        r"""
        Return ``True`` if there is additive reduction.

        EXAMPLES::

            sage: E = EllipticCurve('27a1')
            sage: [(p, E.local_data(p).has_additive_reduction()) for p in prime_range(15)]
            [(2, False), (3, True), (5, False), (7, False), (11, False), (13, False)]

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0, 0, 0, 0, 2*a+1])
            sage: [(p, E.local_data(p).has_additive_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), True)]
        """
        return self._reduction_type == 0

    def _tate(self, proof=None, globally=False):
        r"""
        Tate's algorithm for an elliptic curve over a number field.

        Computes both local reduction data at a prime ideal and a
        local minimal model.

        The model is not required to be integral on input.  If `P` is
        principal, uses a generator as uniformizer, so it will not
        affect integrality or minimality at other primes.  If `P` is not
        principal, the minimal model returned will preserve
        integrality at other primes, but not minimality.

        The optional argument globally, when set to True, tells the algorithm to use the generator of the prime ideal if it is principal. Otherwise just any uniformizer will be used.

        .. NOTE::

            Called only by ``EllipticCurveLocalData.__init__()``.

        OUTPUT:

        (tuple) ``(Emin, p, val_disc, fp, KS, cp)`` where:

        - ``Emin`` -- :class:`EllipticCurve`; a model (integral and) minimal at P
        - ``p`` -- integer; the residue characteristic
        - ``val_disc`` -- integer; the valuation of the local minimal discriminant
        - ``fp`` -- integer; the valuation of the conductor
        - ``KS`` -- string; the Kodaira symbol
        - ``cp`` -- integer; the Tamagawa number

        EXAMPLES (this raised a type error in sage prior to 4.4.4, see :issue:`7930`) ::

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('99d1')
            sage: R.<X> = QQ[]
            sage: K.<t> = NumberField(X^3 + X^2 - 2*X - 1)
            sage: L.<s> = NumberField(X^3 + X^2 - 36*X - 4)
            sage: EK = E.base_extend(K)
            sage: toK = EK.torsion_order()
            sage: da = EK.local_data()  # indirect doctest
            sage: EL = E.base_extend(L)
            sage: da = EL.local_data()  # indirect doctest

        EXAMPLES:

        The following example shows that the bug at :issue:`9324` is fixed::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - x + 6)
            sage: E = EllipticCurve([0,0,0,-53160*a-43995,-5067640*a+19402006])
            sage: E.conductor() # indirect doctest
            Fractional ideal (18, 6*a)

        The following example shows that the bug at :issue:`9417` is fixed::

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(x^2 + 18*x + 1)
            sage: E = EllipticCurve(K, [0, -36, 0, 320, 0])
            sage: E.tamagawa_number(K.ideal(2))
            4

        This is to show that the bug :issue:`11630` is fixed. (The computation of the class group would produce a warning)::

            sage: # needs sage.rings.number_field
            sage: K.<t> = NumberField(x^7 - 2*x + 177)
            sage: E = EllipticCurve([0,1,0,t,t])
            sage: P = K.ideal(2, t^3 + t + 1)
            sage: E.local_data(P).kodaira_symbol()
            II
        """
        E = self._curve
        P = self._prime
        K = E.base_ring()
        OK = K.maximal_order()
        t = verbose("Running Tate's algorithm with P = %s" % P, level=1)
        F = OK.residue_field(P)
        p = F.characteristic()

        # In case P is not principal we mostly use a uniformiser which
        # is globally integral (with positive valuation at some other
        # primes); for this to work, it is essential that we can
        # reduce (mod P) elements of K which are not integral (but are
        # P-integral).  However, if the model is non-minimal and we
        # end up dividing a_i by pi^i then at that point we use a
        # uniformiser pi which has nonpositive valuation at all other
        # primes, so that we can divide by it without losing
        # integrality at other primes.

        if globally:
            principal_flag = P.is_principal()
        else:
            principal_flag = False

        if (K is QQ) or principal_flag:
            pi = P.gens_reduced()[0]
            verbose("P is principal, generator pi = %s" % pi, t, 1)
        else:
            pi = K.uniformizer(P, 'positive')
            verbose("uniformizer pi = %s" % pi, t, 1)
        pi2 = pi * pi
        pi3 = pi * pi2
        pi4 = pi * pi3
        pi_neg = None
        prime = pi if K is QQ else P

        pval = lambda x: x.valuation(prime)
        pdiv = lambda x: x.is_zero() or pval(x) > 0
        # Since ResidueField is cached in a way that
        # does not care much about embeddings of number
        # fields, it can happen that F.p.ring() is different
        # from K. This is a problem: If F.p.ring() has no
        # embedding but K has, then there is no coercion
        # from F.p.ring().maximal_order() to K. But it is
        # no problem to do an explicit conversion in that
        # case (Simon King, github issue #8800).

        from sage.categories.pushout import pushout, CoercionException
        try:
            if hasattr(F.p.ring(), 'maximal_order'): # it is not ZZ
                pushout(F.p.ring().maximal_order(), K)
            pinv = lambda x: F.lift(~F(x))
            proot = lambda x,e: F.lift(F(x).nth_root(e, extend=False, all=True)[0])
            preduce = lambda x: F.lift(F(x))
        except CoercionException: # the pushout does not exist, we need conversion
            pinv = lambda x: K(F.lift(~F(x)))
            proot = lambda x,e: K(F.lift(F(x).nth_root(e, extend=False, all=True)[0]))
            preduce = lambda x: K(F.lift(F(x)))

        def _pquadroots(a, b, c):
            r"""
            Local function returning ``True`` iff `ax^2 + bx + c` has roots
            modulo `P`.
            """
            (a, b, c) = (F(a), F(b), F(c))
            if a == 0:
                return (b != 0) or (c == 0)
            elif p == 2:
                return len(PolynomialRing(F, "x")([c,b,a]).roots()) > 0
            else:
                return (b**2 - 4*a*c).is_square()

        def _pcubicroots(b, c, d):
            r"""
            Local function returning the number of roots of `x^3 +
            b*x^2 + c*x + d` modulo `P`, counting multiplicities
            """

            return sum([rr[1] for rr in PolynomialRing(F, 'x')([F(d), F(c), F(b), F(1)]).roots()],0)

        if p == 2:
            halfmodp = OK(Integer(0))
        else:
            halfmodp = pinv(Integer(2))

        A = E.a_invariants()
        A = [0, A[0], A[1], A[2], A[3], 0, A[4]]
        indices = [1,2,3,4,6]
        if min([pval(a) for a in A if a != 0]) < 0:
            verbose("Non-integral model at P: valuations are %s; making integral" % ([pval(a) for a in A if a != 0]), t, 1)
            e = 0
            for i in range(7):
                if A[i] != 0:
                    e = max(e, (-pval(A[i])/i).ceil())
            pie = pi**e
            for i in range(7):
                if A[i] != 0:
                    A[i] *= pie**i
            verbose("P-integral model is %s, with valuations %s" % ([A[i] for i in indices], [pval(A[i]) for i in indices]), t, 1)

        split = None # only relevant for multiplicative reduction

        (a1, a2, a3, a4, a6) = (A[1], A[2], A[3], A[4], A[6])
        while True:
            C = EllipticCurve([a1, a2, a3, a4, a6])
            (b2, b4, b6, b8) = C.b_invariants()
            (c4, c6) = C.c_invariants()
            delta = C.discriminant()
            val_disc = pval(delta)

            if val_disc == 0:
                ## Good reduction already
                cp = 1
                fp = 0
                KS = KodairaSymbol("I0")
                break #return

            # Otherwise, we change coordinates so that p | a3, a4, a6
            if p == 2:
                if pdiv(b2):
                    r = proot(a4, 2)
                    t = proot(((r + a2)*r + a4)*r + a6, 2)
                else:
                    temp = pinv(a1)
                    r = temp * a3
                    t = temp * (a4 + r*r)
            elif p == 3:
                if pdiv(b2):
                    r = proot(-b6, 3)
                else:
                    r = -pinv(b2) * b4
                t = a1 * r + a3
            else:
                if pdiv(c4):
                    r = -pinv(12) * b2
                else:
                    r = -pinv(12*c4) * (c6 + b2 * c4)
                t = -halfmodp * (a1 * r + a3)
            r = preduce(r)
            t = preduce(t)
            verbose("Before first transform C = %s" % C)
            verbose("[a1,a2,a3,a4,a6] = %s" % ([a1, a2, a3, a4, a6]))
            C = C.rst_transform(r, 0, t)
            (a1, a2, a3, a4, a6) = C.a_invariants()
            (b2, b4, b6, b8) = C.b_invariants()
            if min([pval(a) for a in (a1, a2, a3, a4, a6) if a != 0]) < 0:
                raise RuntimeError("Non-integral model after first transform!")
            verbose("After first transform %s\n, [a1,a2,a3,a4,a6] = %s\n, valuations = %s" % ([r, 0, t], [a1, a2, a3, a4, a6], [pval(a1), pval(a2), pval(a3), pval(a4), pval(a6)]), t, 2)
            if pval(a3) == 0:
                raise RuntimeError("p does not divide a3 after first transform!")
            if pval(a4) == 0:
                raise RuntimeError("p does not divide a4 after first transform!")
            if pval(a6) == 0:
                raise RuntimeError("p does not divide a6 after first transform!")

            # Now we test for Types In, II, III, IV
            # NB the c invariants never change.

            if not pdiv(c4):
                # Multiplicative reduction: Type In (n = val_disc)
                split = False
                if _pquadroots(1, a1, -a2):
                    cp = val_disc
                    split = True
                elif Integer(2).divides(val_disc):
                    cp = 2
                else:
                    cp = 1
                KS = KodairaSymbol("I%s" % val_disc)
                fp = 1
                break #return

            # Additive reduction

            if pval(a6) < 2:
                ## Type II
                KS = KodairaSymbol("II")
                fp = val_disc
                cp = 1
                break #return
            if pval(b8) < 3:
                ## Type III
                KS = KodairaSymbol("III")
                fp = val_disc - 1
                cp = 2
                break #return
            if pval(b6) < 3:
                ## Type IV
                cp = 1
                a3t = preduce(a3/pi)
                a6t = preduce(a6/pi2)
                if _pquadroots(1, a3t, -a6t):
                    cp = 3
                KS = KodairaSymbol("IV")
                fp = val_disc - 2
                break #return

            # If our curve is none of these types, we change coords so that
            # p | a1, a2;  p^2 | a3, a4;  p^3 | a6
            if p == 2:
                s = proot(a2, 2)        # so s^2=a2 (mod pi)
                t = pi*proot(a6/pi2, 2) # so t^2=a6 (mod pi^3)
            elif p == 3:
                s = a1       # so a1'=2s+a1=3a1=0 (mod pi)
                t = a3       # so a3'=2t+a3=3a3=0 (mod pi^2)
            else:
                s = -a1*halfmodp   # so a1'=2s+a1=0 (mod pi)
                t = -a3*halfmodp   # so a3'=2t+a3=0 (mod pi^2)
            C = C.rst_transform(0, s, t)
            (a1, a2, a3, a4, a6) = C.a_invariants()

            verbose("After second transform %s\n[a1, a2, a3, a4, a6] = %s\nValuations: %s" % ([0, s, t], [a1,a2,a3,a4,a6],[pval(a1),pval(a2),pval(a3),pval(a4),pval(a6)]), t, 2)
            if pval(a1) == 0:
                raise RuntimeError("p does not divide a1 after second transform!")
            if pval(a2) == 0:
                raise RuntimeError("p does not divide a2 after second transform!")
            if pval(a3) < 2:
                raise RuntimeError("p^2 does not divide a3 after second transform!")
            if pval(a4) < 2:
                raise RuntimeError("p^2 does not divide a4 after second transform!")
            if pval(a6) < 3:
                raise RuntimeError("p^3 does not divide a6 after second transform!")
            if min(pval(a1), pval(a2), pval(a3), pval(a4), pval(a6)) < 0:
                raise RuntimeError("Non-integral model after second transform!")

            # Analyze roots of the cubic T^3 + bT^2 + cT + d = 0 mod P, where
            # b = a2/p, c = a4/p^2, d = a6/p^3
            b = preduce(a2/pi)
            c = preduce(a4/pi2)
            d = preduce(a6/pi3)
            bb = b*b
            cc = c*c
            bc = b*c
            w = 27*d*d - bb*cc + 4*b*bb*d - 18*bc*d + 4*c*cc
            x = 3*c - bb
            if pdiv(w):
                if pdiv(x):
                    sw = 3
                else:
                    sw = 2
            else:
                sw = 1
            verbose("Analyzing roots of cubic T^3 + %s*T^2 + %s*T + %s, case %s" % (b, c, d, sw), t, 1)
            if sw == 1:
                ## Three distinct roots - Type I*0
                verbose("Distinct roots", t, 1)
                KS = KodairaSymbol("I0*")
                cp = 1 + _pcubicroots(b, c, d)
                fp = val_disc - 4
                break #return
            elif sw == 2:
                ## One double root - Type I*m for some m
                verbose("One double root", t, 1)
                ## Change coords so that the double root is T = 0 mod p
                if p == 2:
                    r = proot(c, 2)
                elif p == 3:
                    r = c * pinv(b)
                else:
                    r = (bc - 9*d)*pinv(2*x)
                r = pi * preduce(r)
                C = C.rst_transform(r, 0, 0)
                (a1, a2, a3, a4, a6) = C.a_invariants()

                # The rest of this branch is just to compute cp, fp, KS.
                # We use pi to keep transforms integral.
                ix = 3
                iy = 3
                mx = pi2
                my = mx
                while True:
                    a2t = preduce(a2 / pi)
                    a3t = preduce(a3 / my)
                    a4t = preduce(a4 / (pi*mx))
                    a6t = preduce(a6 / (mx*my))
                    if pdiv(a3t*a3t + 4*a6t):
                        if p == 2:
                            t = my*proot(a6t, 2)
                        else:
                            t = my*preduce(-a3t*halfmodp)
                        C = C.rst_transform(0, 0, t)
                        (a1, a2, a3, a4, a6) = C.a_invariants()

                        my *= pi
                        iy += 1
                        a2t = preduce(a2 / pi)
                        a3t = preduce(a3/my)
                        a4t = preduce(a4/(pi*mx))
                        a6t = preduce(a6/(mx*my))
                        if pdiv(a4t*a4t - 4*a6t*a2t):
                            if p == 2:
                                r = mx*proot(a6t*pinv(a2t), 2)
                            else:
                                r = mx*preduce(-a4t*pinv(2*a2t))
                            C = C.rst_transform(r, 0, 0)
                            (a1, a2, a3, a4, a6) = C.a_invariants()

                            mx *= pi
                            ix += 1 # and stay in loop
                        else:
                            if _pquadroots(a2t, a4t, a6t):
                                cp = 4
                            else:
                                cp = 2
                            break # exit loop
                    else:
                        if _pquadroots(1, a3t, -a6t):
                            cp = 4
                        else:
                            cp = 2
                        break
                KS = KodairaSymbol("I%s*" % (ix+iy-5))
                fp = val_disc - ix - iy + 1
                break #return
            else: # sw == 3
                ## The cubic has a triple root
                verbose("Triple root", t, 1)
                ## First we change coordinates so that T = 0 mod p
                if p == 2:
                    r = b
                elif p == 3:
                    r = proot(-d, 3)
                else:
                    r = -b * pinv(3)
                r = pi*preduce(r)
                C = C.rst_transform(r, 0, 0)
                (a1, a2, a3, a4, a6) = C.a_invariants()

                verbose("After third transform %s\n[a1,a2,a3,a4,a6] = %s\nValuations: %s" % ([r,0,0],[a1,a2,a3,a4,a6],[pval(ai) for ai in [a1,a2,a3,a4,a6]]), t, 2)
                if min(pval(ai) for ai in [a1,a2,a3,a4,a6]) < 0:
                    raise RuntimeError("Non-integral model after third transform!")
                if pval(a2) < 2 or pval(a4) < 3 or pval(a6) < 4:
                    raise RuntimeError("Cubic after transform does not have a triple root at 0")
                a3t = preduce(a3/pi2)
                a6t = preduce(a6/pi4)
                # We test for Type IV*
                if not pdiv(a3t*a3t + 4*a6t):
                    cp = 3 if _pquadroots(1, a3t, -a6t) else 1
                    KS = KodairaSymbol("IV*")
                    fp = val_disc - 6
                    break #return
                # Now change coordinates so that p^3|a3, p^5|a6
                if p == 2:
                    t = -pi2*proot(a6t, 2)
                else:
                    t = pi2*preduce(-a3t*halfmodp)
                C = C.rst_transform(0, 0, t)
                (a1, a2, a3, a4, a6) = C.a_invariants()

                # We test for types III* and II*
                if pval(a4) < 4:
                    ## Type III*
                    KS = KodairaSymbol("III*")
                    fp = val_disc - 7
                    cp = 2
                    break #return
                if pval(a6) < 6:
                    ## Type II*
                    KS = KodairaSymbol("II*")
                    fp = val_disc - 8
                    cp = 1
                    break #return
                if pi_neg is None:
                    if principal_flag:
                        pi_neg = pi
                    else:
                        pi_neg = K.uniformizer(P, 'negative')
                    pi_neg2 = pi_neg*pi_neg
                    pi_neg3 = pi_neg*pi_neg2
                    pi_neg4 = pi_neg*pi_neg3
                    pi_neg6 = pi_neg4*pi_neg2
                a1 /= pi_neg
                a2 /= pi_neg2
                a3 /= pi_neg3
                a4 /= pi_neg4
                a6 /= pi_neg6
                verbose("Non-minimal equation, dividing out...\nNew model is %s" % ([a1, a2, a3, a4, a6]), t, 1)
        return (C, p, val_disc, fp, KS, cp, split)


def check_prime(K, P):
    r"""
    Function to check that `P` determines a prime of `K`, and return that ideal.

    INPUT:

    - ``K`` -- a number field (including `\QQ`)

    - ``P`` -- an element of ``K`` or a (fractional) ideal of ``K``

    OUTPUT: if ``K`` is `\QQ`: the prime integer equal to or which generates `P`

    - If ``K`` is not `\QQ`: the prime ideal equal to or generated by `P`.

    .. NOTE::

        If `P` is not a prime and does not generate a prime, a :exc:`TypeError`
        is raised.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_local_data import check_prime
        sage: check_prime(QQ, 3)
        3
        sage: check_prime(QQ, QQ(3))
        3
        sage: check_prime(QQ, ZZ.ideal(31))
        31

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 5)
        sage: check_prime(K, a)
        Fractional ideal (a)
        sage: check_prime(K, a + 1)
        Fractional ideal (a + 1)
        sage: [check_prime(K, P) for P in K.primes_above(31)]
        [Fractional ideal (5/2*a + 1/2), Fractional ideal (5/2*a - 1/2)]
        sage: L.<b> = NumberField(x^2 + 3)
        sage: check_prime(K, L.ideal(5))
        Traceback (most recent call last):
        ...
        TypeError: The ideal Fractional ideal (5) is not a prime ideal of
        Number Field in a with defining polynomial x^2 - 5
        sage: check_prime(K, L.ideal(b))
        Traceback (most recent call last):
        ...
        TypeError: No compatible natural embeddings found for
        Number Field in a with defining polynomial x^2 - 5 and
        Number Field in b with defining polynomial x^2 + 3
    """
    if K is QQ:
        if P in ZZ or isinstance(P, (Integer, int)):
            P = Integer(P)
            if P.is_prime():
                return P
            else:
                raise TypeError("The element %s is not prime" % (P,))
        elif P in QQ:
            raise TypeError("The element %s is not prime" % (P,))
        elif isinstance(P, Ideal_generic) and P.base_ring() is ZZ:
            if P.is_prime():
                return P.gen()
            else:
                raise TypeError("The ideal %s is not a prime ideal of %s" % (P, ZZ))
        else:
            raise TypeError("%s is neither an element of QQ or an ideal of %s" % (P, ZZ))

    if not isinstance(K, NumberField):
        raise TypeError("%s is not a number field" % (K,))

    if isinstance(P, NumberFieldFractionalIdeal) or P in K:
        # if P is an ideal, making sure it is a fractional ideal of K
        P = K.fractional_ideal(P)
        if P.is_prime():
            return P
        else:
            raise TypeError("The ideal %s is not a prime ideal of %s" % (P, K))

    raise TypeError("%s is not a valid prime of %s" % (P, K))
