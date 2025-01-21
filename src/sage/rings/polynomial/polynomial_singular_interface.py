"""
Polynomial Interfaces to Singular

AUTHORS:

- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-04-21)
- Robert Bradshaw: Re-factor to avoid multiple inheritance vs. Cython (2007-09)
- Syed Ahmad Lavasani: Added function field to _singular_init_ (2011-12-16);
  Added non-prime finite fields to _singular_init_ (2012-1-22)

TESTS::

    sage: R = PolynomialRing(GF(2**8,'a'), 10, 'x', order='invlex')                     # needs sage.rings.finite_rings
    sage: R == loads(dumps(R))                                                          # needs sage.rings.finite_rings
    True
    sage: P.<a,b> = PolynomialRing(GF(7), 2)
    sage: f = (a^3 + 2*b^2*a)^7; f
    a^21 + 2*a^7*b^14
"""
#################################################################
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
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
#
######################################################################

import sage.rings.fraction_field
import sage.rings.abc
import sage.rings.number_field as number_field

from sage.rings.rational_field import RationalField
from sage.rings.function_field.function_field_rational import RationalFunctionField
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field_base import NumberField

import sage.rings.finite_rings.finite_field_constructor


def _do_singular_init_(singular, base_ring, char, _vars, order):
    r"""
    Implementation of :meth:`PolynomialRing_singular_repr._singular_init_`.

    This code was extracted from :class:`PolynomialRing_singular_repr`
    to make it callable from other places, in particular
    :class:`MPolynomialRing_libsingular`.

    TESTS::

        sage: from sage.rings.polynomial.polynomial_singular_interface import _do_singular_init_
        sage: _do_singular_init_(singular, ZZ, 0, 'X', 'dp')                            # needs sage.libs.singular
        (polynomial ring, over a domain, global ordering
         // coefficients: ZZ
         // number of vars : 1
         //        block   1 : ordering dp
         //                  : names    X
         //        block   2 : ordering C,
         None)
    """
    make_ring = lambda s: singular.ring(s, _vars, order=order)

    if base_ring is ZZ:
        return make_ring("(ZZ)"), None

    if isinstance(base_ring, sage.rings.rational_field.RationalField):
        return make_ring("(QQ)"), None

    elif isinstance(base_ring, sage.rings.abc.RealField):
        # singular converts to bits from base_10 in mpr_complex.cc by:
        #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
        precision = base_ring.precision()
        digits = (2*precision + 4) // 7
        return make_ring(f"(real,{digits},0)"), None

    elif isinstance(base_ring, sage.rings.abc.ComplexField):
        # singular converts to bits from base_10 in mpr_complex.cc by:
        #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
        precision = base_ring.precision()
        digits = (2*precision + 4) // 7
        return make_ring(f"(complex,{digits},0,I)"), None

    elif isinstance(base_ring, sage.rings.abc.RealDoubleField):
        # singular converts to bits from base_10 in mpr_complex.cc by:
        #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
        return make_ring("(real,15,0)"), None

    elif isinstance(base_ring, sage.rings.abc.ComplexDoubleField):
        # singular converts to bits from base_10 in mpr_complex.cc by:
        #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
        return make_ring("(complex,15,0,I)"), None

    elif isinstance(base_ring, sage.rings.abc.IntegerModRing):
        char = base_ring.characteristic()
        if isinstance(base_ring, FiniteField) and char <= 2147483647:
            return make_ring(str(char)), None
        if char.is_power_of(2):
            return make_ring(f"(integer,2,{char.nbits()-1})"), None
        return make_ring(f"(integer,{char})"), None

    elif isinstance(base_ring, FiniteField):
        # not the prime field!
        gen = str(base_ring.gen())
        R = make_ring(f"({char},{gen})")

        minpoly = str(base_ring.modulus()).replace("x",gen).replace(" ","")
        if singular.eval('minpoly') != f"({minpoly})":
            singular.eval(f"minpoly={minpoly}")
            minpoly = singular.eval('minpoly')[1:-1]

        return R, minpoly

    elif isinstance(base_ring, NumberField) and base_ring.is_absolute():
        # not the rationals!
        gen = str(base_ring.gen())
        poly = base_ring.polynomial()
        poly_gen = str(poly.parent().gen())
        poly_str = str(poly).replace(poly_gen,gen)
        R = make_ring(f"({char},{gen})")

        minpoly = poly_str.replace(" ","")
        if singular.eval('minpoly') != f"({minpoly})":
            singular.eval(f"minpoly={minpoly}")
            minpoly = singular.eval('minpoly')[1:-1]

        return R, minpoly

    elif isinstance(base_ring, sage.rings.fraction_field.FractionField_generic):
        if base_ring.ngens() == 1:
            gens = str(base_ring.gen())
        else:
            gens = str(base_ring.gens())

        B = base_ring.base_ring()
        base_char = base_ring.characteristic()

        if B.is_prime_field() or B is ZZ:
            return make_ring(f"({base_char},{gens})"), None

        if isinstance(B, FiniteField) and B.characteristic() <= 2147483647:
            ext_gen = str(B.gen())
            _vars = '(' + ext_gen + ', ' + _vars[1:]

            R = make_ring(f"({base_char},{gens})")

            base_ring.__minpoly = (str(B.modulus()).replace("x",ext_gen)).replace(" ","")
            singular.eval('setring ' + R._name)

            from sage.misc.stopgap import stopgap
            stopgap("Denominators of fraction field elements are sometimes dropped without warning.", 17696)

            return singular(f"std(ideal({base_ring.__minpoly}))", type='qring'), None

    elif isinstance(base_ring, RationalFunctionField) \
            and base_ring.constant_field().is_prime_field():
        gen = str(base_ring.gen())
        return make_ring(f"({base_ring.characteristic()},{gen})"), None

    raise TypeError("no conversion to a Singular ring defined")


class PolynomialRing_singular_repr:
    """
    Implement methods to convert polynomial rings to Singular.

    This class is a base class for all univariate and multivariate
    polynomial rings which support conversion from and to Singular
    rings.
    """
    def _singular_(self, singular=None):
        r"""
        Return a Singular ring for this polynomial ring.

        Currently `\QQ`, `{\rm GF}(p), {\rm GF}(p^n)`, `\CC`, `\RR`, `\ZZ` and
        `\ZZ/n\ZZ` are supported.

        INPUT:

        - ``singular`` -- Singular instance

        OUTPUT: Singular ring matching this ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(CC)                                          # needs sage.rings.real_mpfr
            sage: singular(R)                                                           # needs sage.libs.singular sage.rings.real_mpfr
            polynomial ring, over a field, global ordering
            // coefficients: real[I](complex:15 digits, additional 0 digits)/(I^2+1)
            // number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R.<x,y> = PolynomialRing(RealField(100))                              # needs sage.rings.real_mpfr
            sage: singular(R)                                                           # needs sage.libs.singular sage.rings.real_mpfr
            polynomial ring, over a field, global ordering
            // coefficients: Float()
            // number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: w = polygen(ZZ, 'w')
            sage: R.<x> = PolynomialRing(NumberField(w^2 + 1, 's'))                     # needs sage.rings.number_field
            sage: singular(R)                                                           # needs sage.libs.singular sage.rings.number_field
            polynomial ring, over a field, global ordering
            //   coefficients: QQ[s]/(s^2+1)
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127), 'x', implementation='singular')           # needs sage.libs.singular
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/127
            //   number of vars : 1
            //        block   1 : ordering dp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ, 'x', implementation='singular')                # needs sage.libs.singular
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering dp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,'x')
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127), 'x')
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = Frac(ZZ['a,b'])['x,y']
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: QQ(a, b)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C


            sage: R = IntegerModRing(1024)['x,y']
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a ring (with zero-divisors), global ordering
            //   coefficients: ZZ/(2^10)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = IntegerModRing(15)['x,y']
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a ring (with zero-divisors), global ordering
            //   coefficients: ZZ/...(15)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = ZZ['x,y']
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a domain, global ordering
            //   coefficients: ZZ
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = ZZ['x']
            sage: singular(R)                                                           # needs sage.libs.singular
            polynomial ring, over a domain, global ordering
            // coefficients: ZZ
            // number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = FiniteField(25)
            sage: R = k['x']
            sage: K = R.fraction_field()
            sage: S = K['y']
            sage: singular(S)                                                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/5(x)
            //   number of vars : 2
            //        block   1 : ordering lp
            //                  : names    a y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=a2-a+2

        .. warning::

            - If the base ring is a finite extension field or a number field
              the ring will not only be returned but also be set as the current
              ring in Singular.
            - Singular represents precision of floating point numbers base 10
              while Sage represents floating point precision base 2.
        """
        if singular is None:
            from sage.interfaces.singular import singular
        try:
            R = self.__singular
            if R.parent() is not singular:
                raise ValueError
            R._check_valid()
            if self.base_ring() is ZZ or self.base_ring().is_prime_field():
                return R
            if isinstance(self.base_ring(), FiniteField) or \
                    (isinstance(self.base_ring(), NumberField) and self.base_ring().is_absolute()):
                R.set_ring()  # sorry for that, but needed for minpoly
                if singular.eval('minpoly') != f"({self.__minpoly})":
                    singular.eval(f"minpoly={self.__minpoly}")
                    self.__minpoly = singular.eval('minpoly')[1:-1]
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=None):
        """
        Return a newly created Singular ring matching this ring.

        EXAMPLES::

            sage: PolynomialRing(QQ,'u_ba')._singular_init_()                           # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    u_ba
            //        block   2 : ordering C
        """
        if not can_convert_to_singular(self):
            raise TypeError("no conversion of this ring to a Singular ring defined")

        if self.ngens() == 1:
            _vars = f'({self.gen()})'
            if "*" in _vars:  # 1.000...000*x
                _vars = _vars.split("*")[1]
            order = 'lp'
        else:
            _vars = str(self.gens())
            order = self.term_order().singular_str()

        if singular is None:
            from sage.interfaces.singular import singular

        self.__singular, self.__minpoly = _do_singular_init_(singular, self.base_ring(), self.characteristic(), _vars, order)

        return self.__singular


def can_convert_to_singular(R):
    """
    Return ``True`` if this ring's base field or ring can be
    represented in Singular, and the polynomial ring has at
    least one generator.

    The following base rings are supported: finite fields,
    rationals, number fields, and real and complex fields.

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular
        sage: can_convert_to_singular(PolynomialRing(QQ, names=['x']))
        True
        sage: can_convert_to_singular(PolynomialRing(ZZ, names=['x']))
        True

        sage: can_convert_to_singular(PolynomialRing(QQ, names=[]))
        False

    TESTS:

    Avoid non absolute number fields (see :issue:`23535`)::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a,b> = NumberField([x^2 - 2, x^2 - 5])                                 # needs sage.rings.number_field
        sage: can_convert_to_singular(K['s,t'])                                         # needs sage.rings.number_field
        False

    Check for :issue:`33319`::

        sage: # needs sage.rings.finite_rings
        sage: R.<x,y> = GF((2^31-1)^3)[]
        sage: R._has_singular
        True
        sage: R.<x,y> = GF((2^31+11)^2)[]
        sage: R._has_singular
        False
        sage: R.<x,y> = GF(10^20 - 11)[]
        sage: R._has_singular
        True

        sage: R.<x,y> = Zmod(10^20 + 1)[]
        sage: R._has_singular
        True
    """
    if R.ngens() == 0:
        return False

    base_ring = R.base_ring()
    if (base_ring is ZZ
        or isinstance(base_ring, (RationalField,
                                  sage.rings.abc.IntegerModRing,
                                  sage.rings.abc.RealField, sage.rings.abc.ComplexField,
                                  sage.rings.abc.RealDoubleField, sage.rings.abc.ComplexDoubleField))):
        return True
    elif isinstance(base_ring, FiniteField):
        return base_ring.characteristic() <= 2147483647
    elif isinstance(base_ring, NumberField):
        return base_ring.is_absolute()
    elif isinstance(base_ring, sage.rings.fraction_field.FractionField_generic):
        B = base_ring.base_ring()
        return (B.is_prime_field() or B is ZZ
                or (isinstance(B, FiniteField) and B.characteristic() <= 2147483647))
    elif isinstance(base_ring, RationalFunctionField):
        return base_ring.constant_field().is_prime_field()
    else:
        return False


class Polynomial_singular_repr:
    """
    Implement coercion of polynomials to Singular polynomials.

    This class is a base class for all (univariate and multivariate)
    polynomial classes which support conversion from and to
    Singular polynomials.

    Due to the incompatibility of Python extension classes and multiple inheritance,
    this just defers to module-level functions.
    """
    def _singular_(self, singular=None):
        if singular is None:
            from sage.interfaces.singular import singular
        return _singular_func(self, singular)

    def _singular_init_func(self, singular=None):
        if singular is None:
            from sage.interfaces.singular import singular
        return _singular_init_func(self, singular)


def _singular_func(self, singular=None):
    """
    Return Singular polynomial matching this polynomial.

    INPUT:

    - ``singular`` -- Singular instance to use

    EXAMPLES::

        sage: # needs sage.libs.singular
        sage: P.<a,b> = PolynomialRing(GF(7), 2)
        sage: f = (a^3 + 2*b^2*a)^7; f
        a^21 + 2*a^7*b^14
        sage: h = f._singular_(); h
        a^21+2*a^7*b^14
        sage: P(h)
        a^21 + 2*a^7*b^14
        sage: P(h^20) == f^20
        True

        sage: # needs sage.libs.singular
        sage: R.<x> = PolynomialRing(GF(7))
        sage: f = (x^3 + 2*x^2*x)^7
        sage: f
        3*x^21
        sage: h = f._singular_(); h
        3*x^21
        sage: R(h)
        3*x^21
        sage: R(h^20) == f^20
        True
    """
    if singular is None:
        from sage.interfaces.singular import singular
    self.parent()._singular_(singular).set_ring()  # this is expensive
    try:
        self.__singular._check_valid()
        if self.__singular.parent() is singular:
            return self.__singular
    except (AttributeError, ValueError):
        pass
    return _singular_init_func(self, singular)


def _singular_init_func(self, singular=None):
    """
    Return corresponding Singular polynomial but enforce that a new
    instance is created in the Singular interpreter.

    Use ``self._singular_()`` instead.
    """
    if singular is None:
        from sage.interfaces.singular import singular
    self.parent()._singular_(singular).set_ring()  # this is expensive
    self.__singular = singular(str(self))
    return self.__singular
