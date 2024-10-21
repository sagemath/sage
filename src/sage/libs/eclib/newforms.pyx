# sage_setup: distribution = sagemath-eclib
"""
Modular symbols using eclib newforms
"""

# ****************************************************************************
#       Copyright (C) 2008 Tom Boothby <boothby@u.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.eclib cimport *
from sage.libs.gmp.mpq cimport mpq_numref
from sage.libs.ntl.convert cimport mpz_to_ZZ
from sage.rings.rational_field import QQ
from sage.rings.rational cimport Rational
from sage.modular.all import Cusp


cdef class ECModularSymbol:
    r"""
    Modular symbol associated with an elliptic curve,  using John Cremona's newforms class.

    EXAMPLES::

        sage: from sage.libs.eclib.newforms import ECModularSymbol
        sage: E = EllipticCurve('11a')
        sage: M = ECModularSymbol(E,1); M
        Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

    By default, symbols are based at the cusp `\infty`, i.e. we evaluate `\{\infty,r\}`::

        sage: [M(1/i) for i in range(1,11)]
        [2/5, -8/5, -3/5, 7/5, 12/5, 12/5, 7/5, -3/5, -8/5, 2/5]

    We can also switch the base point to the cusp `0`::

        sage: [M(1/i, base_at_infinity=False) for i in range(1,11)]
        [0, -2, -1, 1, 2, 2, 1, -1, -2, 0]

    For the minus symbols this makes no difference since
    `\{0,\infty\}` is in the plus space.  Note that to evaluate minus
    symbols the space must be defined with sign 0, which makes both
    signs available::

        sage: M = ECModularSymbol(E,0); M
        Modular symbol with sign 0 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: [M(1/i, -1) for i in range(1,11)]
        [0, 0, 1, 1, 0, 0, -1, -1, 0, 0]
        sage: [M(1/i, -1, base_at_infinity=False) for i in range(1,11)]
        [0, 0, 1, 1, 0, 0, -1, -1, 0, 0]

    If the ECModularSymbol is created with sign 0 then as well as
    asking for both + and - symbols, we can also obtain both (as a
    tuple).  However it is more work to create the full modular
    symbol space::

        sage: E = EllipticCurve('11a1')
        sage: M = ECModularSymbol(E,0); M
        Modular symbol with sign 0 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: [M(1/i) for i in range(2,11)]
        [[-8/5, 0],
         [-3/5, 1],
         [7/5, 1],
         [12/5, 0],
         [12/5, 0],
         [7/5, -1],
         [-3/5, -1],
         [-8/5, 0],
         [2/5, 0]]

    The curve is automatically converted to its minimal model::

        sage: E = EllipticCurve([0,0,0,0,1/4])
        sage: ECModularSymbol(E)
        Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 over Rational Field

    Non-optimal curves are handled correctly in eclib, by comparing the ratios of real and/or imaginary periods::

        sage: from sage.libs.eclib.newforms import ECModularSymbol
        sage: E1 = EllipticCurve('11a1') # optimal
        sage: E1.period_lattice().basis()
        (1.26920930427955, 0.634604652139777 + 1.45881661693850*I)
        sage: M1 = ECModularSymbol(E1,0)
        sage: M1(0)
        [2/5, 0]
        sage: M1(1/3)
        [-3/5, 1]

    One non-optimal curve has real period 1/5 that of the optimal one, so plus symbols scale up by a factor of 5 while minus symbols are unchanged::

        sage: E2 = EllipticCurve('11a2') # not optimal
        sage: E2.period_lattice().basis()
        (0.253841860855911, 0.126920930427955 + 1.45881661693850*I)
        sage: M2 = ECModularSymbol(E2,0)
        sage: M2(0)
        [2, 0]
        sage: M2(1/3)
        [-3, 1]
        sage: all((M2(r,1)==5*M1(r,1)) for r in QQ.range_by_height(10))
        True
        sage: all((M2(r,-1)==M1(r,-1)) for r in QQ.range_by_height(10))
        True

    The other non-optimal curve has real period 5 times that of the optimal one, so plus symbols scale down by a factor of 5; again, minus symbols are unchanged::

        sage: E3 = EllipticCurve('11a3') # not optimal
        sage: E3.period_lattice().basis()
        (6.34604652139777, 3.17302326069888 + 1.45881661693850*I)
        sage: M3 = ECModularSymbol(E3,0)
        sage: M3(0)
        [2/25, 0]
        sage: M3(1/3)
        [-3/25, 1]
        sage: all((5*M3(r,1)==M1(r,1)) for r in QQ.range_by_height(10))
        True
        sage: all((M3(r,-1)==M1(r,-1)) for r in QQ.range_by_height(10))
        True

    TESTS::

        sage: ECModularSymbol.__new__(ECModularSymbol)
        Modular symbol with sign 0 over Rational Field attached to None
    """
    def __init__(self, E, sign=1, nap=1000):
        """
        Construct the modular symbol object from an elliptic curve.

        INPUT:

        - ``E``- an elliptic curve defined over Q.

        - ``sign`` -- integer; 0 or +1.  If +1, only plus modular symbols
         of this sign are available.  If 0, modular symbols of both
         signs are available but the construction is more expensive.

        - ``nap`` -- integer (default: 1000); the number of ap of E to use
         in determining the normalisation of the modular symbols.
         Note that eclib will increase this to 100*sqrt(N) if necessary.

        EXAMPLES::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E, +1)
            sage: M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E, -1)
            Traceback (most recent call last):
            ...
            ValueError: ECModularSymbol can only be created with signs +1 or 0, not -1
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E, 0)
            sage: M
            Modular symbol with sign 0 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

        TESTS:

        This one is from :issue:`8042`::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('858k2')
            sage: ECModularSymbol(E)
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + x*y = x^3 + 16353089*x - 335543012233 over Rational Field

        We allow a-invariants which are larger than 64 bits
        (:issue:`16977`)::

            sage: E = EllipticCurve([-25194941007454971, -1539281792450963687794218])  # non-minimal model of 21758k3
            sage: ECModularSymbol(E)  # long time
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + x*y = x^3 - 19440540900814*x - 32992152521343165084 over Rational Field

        The curve needs to be defined over the rationals::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: ECModularSymbol(E)
            Traceback (most recent call last):
            ...
            TypeError: the elliptic curve must have coefficients in the rational numbers
        """
        cdef ZZ_c a1, a2, a3, a4, a6, N
        cdef Curve C
        cdef Curvedata CD
        cdef CurveRed CR
        cdef int n

        if E.base_field() is not QQ:
            raise TypeError("the elliptic curve must have coefficients in the rational numbers")

        E = E.minimal_model()
        self._E = E

        # The a invariants are rational numbers with denominator 1
        mpz_to_ZZ(&a1, mpq_numref((<Rational>(E.a1())).value))
        mpz_to_ZZ(&a2, mpq_numref((<Rational>(E.a2())).value))
        mpz_to_ZZ(&a3, mpq_numref((<Rational>(E.a3())).value))
        mpz_to_ZZ(&a4, mpq_numref((<Rational>(E.a4())).value))
        mpz_to_ZZ(&a6, mpq_numref((<Rational>(E.a6())).value))

        sig_on()
        C = Curve(a1,a2,a3,a4,a6)
        CD = Curvedata(C,0)
        CR = CurveRed(CD)
        N = getconductor(CR)
        n = I2int(N)
        self.n = n
        if not (sign == 0 or sign == 1):
            sig_off()
            raise ValueError("ECModularSymbol can only be created with signs +1 or 0, not {}".format(sign))
        self.sign = sign

        self.nfs = new newforms(n, 0)
        self.nfs.createfromcurve(sign, CR, nap)
        sig_off()

    def __dealloc__(self):
        del self.nfs

    def __repr__(self):
        """
        TESTS::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E, 1); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E, 0); M
            Modular symbol with sign 0 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return "Modular symbol with sign %s over Rational Field attached to %s"%(self.sign, self._E)

    def __call__(self, r, sign=None, base_at_infinity=True):
        r"""
        Compute the value of ``self`` on `\{0,r\}` or `\{\infty, r\}` for
        rational `r`.

        INPUT:

        - ``r`` -- rational; a rational number

        - ``sign`` -- integer; either +1, -1 or 0.  If the sign of the
          space is +1, only sign +1 is allowed.  Default: ``self.sign``, or +1
          when ``self.sign==0``.

        - ``base_at_infinity`` -- boolean; if ``True``, evaluates
          `\{\infty, r\}`. Otherwise (default) evaluates `\{0,r\}`.

        OUTPUT:

        If sign=+1, the rational value of the plus modular symbol is
        returned.  If sign=-1, the rational value of the minus modular
        symbol is returned.  If sign=0, a tuple of both the plus and
        minus modular symbols is returned.  In the last two cases, the
        space must have sign 0 or an error will be raised.

        EXAMPLES::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,11)]
            [2/5, -8/5, -3/5, 7/5, 12/5, 12/5, 7/5, -3/5, -8/5, 2/5]
            sage: [M(1/i, base_at_infinity=False) for i in range(1,11)]
            [0, -2, -1, 1, 2, 2, 1, -1, -2, 0]
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,11)]
            [0, 0, 0, 0, 1, 0, 1, 1, 0, 0]
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,11)]
            [0, 0, 0, 0, 2, 0, 1, 0, -1, 0]

        When the class is created with sign 0 we can ask for +1 or -1 symbols or (by default) both::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a1')
            sage: M = ECModularSymbol(E,0)
            sage: [M(1/i, 1) for i in range(2,11)]
            [-8/5, -3/5, 7/5, 12/5, 12/5, 7/5, -3/5, -8/5, 2/5]
            sage: [M(1/i, -1) for i in range(2,11)]
            [0, 1, 1, 0, 0, -1, -1, 0, 0]
            sage: [M(1/i) for i in range(2,11)]
            [[-8/5, 0],
             [-3/5, 1],
             [7/5, 1],
             [12/5, 0],
             [12/5, 0],
             [7/5, -1],
             [-3/5, -1],
             [-8/5, 0],
             [2/5, 0]]

        When the class is created with sign +1 we can only ask for +1 symbols::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a1')
            sage: M = ECModularSymbol(E)
            sage: M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: M(0)
            2/5
            sage: M(1/3, sign=-1)
            Traceback (most recent call last):
            ...
            ValueError: impossible to evaluate a minus symbol on a plus space
            sage: M(2/3, sign=0)
            Traceback (most recent call last):
            ...
            ValueError: impossible to evaluate both symbols on a plus space

        TESTS (see :issue:`11211`)::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: M(oo)
            0
            sage: M(oo, base_at_infinity=False)
            -2/5
            sage: M(7/5)
            -13/5
            sage: M("garbage")
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'garbage' to a cusp
            sage: M(7/5)
            -13/5
        """
        cdef rational _r
        cdef rational _sp, _sm
        cdef pair[rational,rational] _spm
        cdef long n, d

        r = Cusp(r)
        d = r.denominator()
        n = r.numerator()
        if d != 0:
            n = n % d
        sig_on()
        _r = rational(n, d)
        if sign is None or sign not in [-1, 0, 1]:
            sign = self.sign
        if sign == +1:
            _sp = self.nfs.plus_modular_symbol(_r, 0, int(base_at_infinity))
        elif sign == -1:
            if self.sign != 0:
                sig_off()
                raise ValueError("impossible to evaluate a minus symbol on a plus space")
            else:
                _sm = self.nfs.minus_modular_symbol(_r, 0, int(base_at_infinity))
        else:  # sign == 0
            if self.sign != 0:
                sig_off()
                raise ValueError("impossible to evaluate both symbols on a plus space")
            else:
                _spm = self.nfs.full_modular_symbol(_r, 0, int(base_at_infinity))
                _sp = _spm.first
                _sm = _spm.second

        sig_off()

        if sign == +1:
            return Rational((rational_num(_sp), rational_den(_sp)))
        elif sign == -1:
            return Rational((rational_num(_sm), rational_den(_sm)))
        else:
            return [Rational((rational_num(_sp), rational_den(_sp))),
                    Rational((rational_num(_sm), rational_den(_sm)))]

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: M.__reduce__()
            (<class 'sage.libs.eclib.newforms.ECModularSymbol'>,
             (Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field,
              1))
        """
        return (ECModularSymbol, (self._E, self.sign))
