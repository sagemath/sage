r"""
Double Precision Real Numbers, implementation using GSL
"""

cimport libc.math

from cysignals.signals cimport sig_on, sig_off

from sage.arith.constants cimport *

from sage.libs.gsl.errno cimport gsl_set_error_handler_off
from sage.libs.gsl.math cimport *
from sage.libs.gsl.exp cimport *
from sage.libs.gsl.log cimport gsl_sf_log
from sage.libs.gsl.trig cimport *
from sage.libs.gsl.dilog cimport gsl_sf_dilog
from sage.libs.gsl.gamma cimport gsl_sf_fact, gsl_sf_gamma
from sage.libs.gsl.zeta cimport gsl_sf_zeta
from sage.libs.gsl.erf cimport gsl_sf_erf


gsl_set_error_handler_off()


cdef class RealDoubleElement_gsl(RealDoubleElement):

    def nth_root(self, int n):
        """
        Return the `n`-th root of ``self``.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        The output is a complex double if ``self`` is negative and `n` is even,
        otherwise it is a real double.

        EXAMPLES::

            sage: r = RDF(-125.0); r.nth_root(3)
            -5.000000000000001
            sage: r.nth_root(5)
            -2.6265278044037674
            sage: RDF(-2).nth_root(5)^5  # rel tol 1e-15
            -2.000000000000001
            sage: RDF(-1).nth_root(5)^5
            -1.0
            sage: RDF(3).nth_root(10)^10
            2.9999999999999982
            sage: RDF(-1).nth_root(2)
            6.123233995736757e-17 + 1.0*I
            sage: RDF(-1).nth_root(4)
            0.7071067811865476 + 0.7071067811865475*I
        """
        if n == 0:
            return RealDoubleElement(float('nan'))
        if self._value < 0:
            if GSL_IS_EVEN(n):
                from sage.rings.complex_double import CDF
                return self._complex_double_(CDF).nth_root(n)
            else:
                return - ((-self) ** (float(1)/n))

        return self ** (float(1)/n)

    cdef __pow_double(self, double exponent, double sign):
        """
        If ``sign == 1`` or ``self >= 0``, return ``self ^ exponent``.
        If ``sign == -1`` and ``self < 0``, return ``- abs(self) ^ exponent``.
        """
        cdef double v = self._value
        if v >= 0:
            if v == 1:
                return self
            elif exponent == 0:
                return self._new_c(1.0)
            elif v == 0:
                if exponent < 0:
                    raise ZeroDivisionError("0.0 cannot be raised to a negative power")
                return self
            sign = 1.0
        else:  # v < 0
            expmod2 = libc.math.fmod(exponent, 2.0)
            if expmod2 == 0.0:
                pass
            elif expmod2 == 1.0:
                sign = -1.0
            else:
                raise ValueError("negative number cannot be raised to a fractional power")
            v = -v
        return self._new_c(sign * gsl_sf_exp(gsl_sf_log(v) * exponent))

    cpdef _pow_(self, other):
        """
        Return ``self`` raised to the real double power ``other``.

        EXAMPLES::

            sage: a = RDF('1.23456')
            sage: a^a
            1.2971114817819216

        TESTS::

            sage: RDF(0) ^ RDF(0.5)
            0.0
            sage: RDF(0) ^ (1/2)
            0.0
            sage: RDF(0) ^ RDF(0)
            1.0
            sage: RDF(0) ^ RDF(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: 0.0 cannot be raised to a negative power
            sage: RDF(-1) ^ RDF(0)
            1.0
            sage: RDF(-1) ^ RDF(1)
            -1.0
            sage: RDF(-1) ^ RDF(0.5)
            Traceback (most recent call last):
            ...
            ValueError: negative number cannot be raised to a fractional power
        """
        return self.__pow_double((<RealDoubleElement>other)._value, 1)

    cpdef _pow_int(self, n):
        """
        Return ``self`` raised to the integer power ``n``.

        TESTS::

            sage: RDF(1) ^ (2^1000)
            1.0
            sage: RDF(1) ^ (2^1000 + 1)
            1.0
            sage: RDF(1) ^ (-2^1000)
            1.0
            sage: RDF(1) ^ (-2^1000 + 1)
            1.0
            sage: RDF(-1) ^ (2^1000)
            1.0
            sage: RDF(-1) ^ (2^1000 + 1)
            -1.0
            sage: RDF(-1) ^ (-2^1000)
            1.0
            sage: RDF(-1) ^ (-2^1000 + 1)
            -1.0

        ::

            sage: base = RDF(1.0000000000000002)
            sage: base._pow_int(0)
            1.0
            sage: base._pow_int(1)
            1.0000000000000002
            sage: base._pow_int(2)
            1.0000000000000004
            sage: base._pow_int(3)
            1.0000000000000007
            sage: base._pow_int(2^57)
            78962960182680.42
            sage: base._pow_int(2^57 + 1)
            78962960182680.42

        ::

            sage: base = RDF(-1.0000000000000002)
            sage: base._pow_int(0)
            1.0
            sage: base._pow_int(1)
            -1.0000000000000002
            sage: base._pow_int(2)
            1.0000000000000004
            sage: base._pow_int(3)
            -1.0000000000000007
            sage: base._pow_int(2^57)
            78962960182680.42
            sage: base._pow_int(2^57 + 1)
            -78962960182680.42
        """
        return self.__pow_double(n, -1.0 if (n & 1) else 1.0)

    cdef _pow_long(self, long n):
        """
        Compute ``self`` raised to the power ``n``.

        EXAMPLES::

            sage: RDF('1.23456') ^ 20
            67.64629770385...
            sage: RDF(3) ^ 32
            1853020188851841.0
            sage: RDF(2)^(-1024)
            5.562684646268003e-309

        TESTS::

            sage: base = RDF(1.0000000000000002)
            sage: base ^ RDF(2^31)
            1.000000476837272
            sage: base ^ (2^57)
            78962960182680.42
            sage: base ^ RDF(2^57)
            78962960182680.42
        """
        if -2048 <= n <= 2048:
            # For small exponents, it is possible that the powering
            # is exact either because the base is a power of 2
            # (e.g. 2.0^1000) or because the exact result has few
            # significant digits (e.g. 3.0^10). Here, we use the
            # square-and-multiply algorithm by GSL.
            return self._new_c(gsl_pow_int(self._value, <int>n))
        # If the exponent is sufficiently large in absolute value, the
        # result cannot be exact (except if the base is -1.0, 0.0 or
        # 1.0 but those cases are handled by __pow_double too). The
        # log-and-exp algorithm from __pow_double will be more precise
        # than square-and-multiply.

        # We do need to take care of the sign since the conversion
        # of n to double might change an odd number to an even number.
        return self.__pow_double(<double>n, -1.0 if (n & 1) else 1.0)

    cdef _log_base(self, double log_of_base):
        if self._value == 0:
            from sage.rings.real_double import RDF
            return RDF(-1)/RDF(0)
        elif self._value < 0:
            from sage.rings.real_double import RDF
            return RDF.NaN()
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) / log_of_base)
        sig_off()
        return a

    def log(self, base=None):
        """
        Return the logarithm.

        INPUT:

        - ``base`` -- integer or ``None`` (default). The base of the
          logarithm. If ``None`` is specified, the base is `e` (the so-called
          natural logarithm).

        OUTPUT:

        The logarithm of ``self``.  If ``self`` is positive, a double
        floating point number. Infinity if ``self`` is zero. A
        imaginary complex floating point number if ``self`` is
        negative.

        EXAMPLES::

            sage: RDF(2).log()
            0.6931471805599453
            sage: RDF(2).log(2)
            1.0
            sage: RDF(2).log(pi)                                                        # needs sage.symbolic
            0.6055115613982801
            sage: RDF(2).log(10)
            0.30102999566398114
            sage: RDF(2).log(1.5)
            1.7095112913514547
            sage: RDF(0).log()
            -infinity
            sage: RDF(-1).log()
            3.141592653589793*I
            sage: RDF(-1).log(2)  # rel tol 1e-15
            4.532360141827194*I

        TESTS:

        Make sure that we can take the log of small numbers accurately
        and the fix doesn't break preexisting values (:issue:`12557`)::

            sage: R = RealField(128)
            sage: def check_error(x):
            ....:   x = RDF(x)
            ....:   log_RDF = x.log()
            ....:   log_RR = R(x).log()
            ....:   diff = R(log_RDF) - log_RR
            ....:   if abs(diff) < log_RDF.ulp():
            ....:       return True
            ....:   print("logarithm check failed for %s (diff = %s ulp)"% \
            ....:       (x, diff/log_RDF.ulp()))
            ....:   return False
            sage: all( check_error(2^x) for x in range(-100,100) )
            True
            sage: all( check_error(x) for x in sxrange(0.01, 2.00, 0.01) )
            True
            sage: all( check_error(x) for x in sxrange(0.99, 1.01, 0.001) )
            True
            sage: RDF(1.000000001).log()
            1.000000082240371e-09
            sage: RDF(1e-17).log()
            -39.14394658089878
            sage: RDF(1e-50).log()
            -115.12925464970229
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(base)
        if base is None:
            return self._log_base(1)
        else:
            if isinstance(base, RealDoubleElement):
                return self._log_base(base._log_base(1))
            else:
                return self._log_base(gsl_sf_log(float(base)))

    def log2(self):
        """
        Return log to the base 2 of ``self``.

        EXAMPLES::

            sage: r = RDF(16.0)
            sage: r.log2()
            4.0

        ::

            sage: r = RDF(31.9); r.log2()
            4.9954845188775066
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(2)
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) * M_1_LN2)
        sig_off()
        return a

    def log10(self):
        """
        Return log to the base 10 of ``self``.

        EXAMPLES::

            sage: r = RDF('16.0'); r.log10()
            1.2041199826559248
            sage: r.log() / RDF(log(10))                                                # needs sage.symbolic
            1.2041199826559246
            sage: r = RDF('39.9'); r.log10()
            1.6009728956867482
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(10)
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) * M_1_LN10)
        sig_off()
        return a

    def logpi(self):
        r"""
        Return log to the base `\pi` of ``self``.

        EXAMPLES::

            sage: r = RDF(16); r.logpi()
            2.4220462455931204
            sage: r.log() / RDF(log(pi))                                                # needs sage.symbolic
            2.4220462455931204
            sage: r = RDF('39.9'); r.logpi()
            3.2203023346075152
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(M_PI)
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) * M_1_LNPI)
        sig_off()
        return a

    def exp(self):
        r"""
        Return `e^\mathtt{self}`.

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp()
            1.0

        ::

            sage: r = RDF('32.3')
            sage: a = r.exp(); a
            106588847274864.47
            sage: a.log()
            32.3

        ::

            sage: r = RDF('-32.3')
            sage: r.exp()
            9.381844588498685e-15

        ::

            sage: RDF(1000).exp()
            +infinity
        """
        sig_on()
        a = self._new_c(gsl_sf_exp(self._value))
        sig_off()
        return a

    def exp2(self):
        r"""
        Return `2^\mathtt{self}`.

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp2()
            1.0

        ::

            sage: r = RDF(32.0)
            sage: r.exp2()
            4294967295.9999967

        ::

            sage: r = RDF(-32.3)
            sage: r.exp2()
            1.8911724825302065e-10
        """
        sig_on()
        a = self._new_c(gsl_sf_exp(self._value * M_LN2))
        sig_off()
        return a

    def exp10(self):
        r"""
        Return `10^\mathtt{self}`.

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp10()
            1.0

        ::

            sage: r = RDF(32.0)
            sage: r.exp10()
            1.0000000000000069e+32

        ::

            sage: r = RDF(-32.3)
            sage: r.exp10()
            5.011872336272702e-33
        """
        sig_on()
        a = self._new_c(gsl_sf_exp(self._value * M_LN10))
        sig_off()
        return a

    def cos(self):
        """
        Return the cosine of ``self``.

        EXAMPLES::

            sage: t=RDF.pi()/2
            sage: t.cos()
            6.123233995736757e-17
        """
        return self._new_c(gsl_sf_cos(self._value))

    def sin(self):
        """
        Return the sine of ``self``.

        EXAMPLES::

            sage: RDF(2).sin()
            0.9092974268256817
        """
        return self._new_c(gsl_sf_sin(self._value))

    def dilog(self):
        r"""
        Return the dilogarithm of ``self``.

        This is defined by the
        series `\sum_n x^n/n^2` for `|x| \le 1`. When the absolute
        value of ``self`` is greater than 1, the returned value is the
        real part of (the analytic continuation to `\CC` of) the
        dilogarithm of ``self``.

        EXAMPLES::

            sage: RDF(1).dilog()  # rel tol 1.0e-13
            1.6449340668482264
            sage: RDF(2).dilog()  # rel tol 1.0e-13
            2.46740110027234
        """
        return self._new_c(gsl_sf_dilog(self._value))

    def restrict_angle(self):
        r"""
        Return a number congruent to ``self`` mod `2\pi` that lies in
        the interval `(-\pi, \pi]`.

        Specifically, it is the unique `x \in (-\pi, \pi]` such
        that ``self`` `= x + 2\pi n` for some `n \in \ZZ`.

        EXAMPLES::

            sage: RDF(pi).restrict_angle()                                              # needs sage.symbolic
            3.141592653589793
            sage: RDF(pi + 1e-10).restrict_angle()                                      # needs sage.symbolic
            -3.1415926534897936
            sage: RDF(1+10^10*pi).restrict_angle()                                      # needs sage.symbolic
            0.9999977606...
        """
        return self._new_c(gsl_sf_angle_restrict_symm(self._value))

    def tan(self):
        """
        Return the tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/3
            sage: q.tan()
            1.7320508075688767
            sage: q = RDF.pi()/6
            sage: q.tan()
            0.5773502691896256
        """
        cos = gsl_sf_cos(self._value)
        a = self._new_c(gsl_sf_sin(self._value) / cos)
        return a

    def sincos(self):
        """
        Return a pair consisting of the sine and cosine of ``self``.

        EXAMPLES::

            sage: t = RDF.pi()/6
            sage: t.sincos()
            (0.49999999999999994, 0.8660254037844387)
        """
        return self.sin(), self.cos()

    def hypot(self, other):
        r"""
        Compute the value `\sqrt{s^2 + o^2}` where `s` is ``self`` and `o`
        is ``other`` in such a way as to avoid overflow.

        EXAMPLES::

            sage: x = RDF(4e300); y = RDF(3e300)
            sage: x.hypot(y)
            5e+300
            sage: sqrt(x^2+y^2) # overflow
            +infinity
        """
        sig_on()
        a = self._new_c(gsl_sf_hypot(self._value, float(other)))
        sig_off()
        return a

    def arccos(self):
        """
        Return the inverse cosine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/3
            sage: i = q.cos()
            sage: i.arccos() == q
            True
        """
        return self._new_c(libc.math.acos(self._value))

    def arcsin(self):
        """
        Return the inverse sine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/5
            sage: i = q.sin()
            sage: i.arcsin() == q
            True
        """
        return self._new_c(libc.math.asin(self._value))

    def arctan(self):
        """
        Return the inverse tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/5
            sage: i = q.tan()
            sage: i.arctan() == q
            True
        """
        return self._new_c(libc.math.atan(self._value))

    def cosh(self):
        """
        Return the hyperbolic cosine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.cosh()
            1.0344656400955106
        """
        return self._new_c(gsl_ldexp(gsl_sf_exp(self._value) + gsl_sf_exp(-self._value), -1))  # (e^x + e^-x)/2

    def sinh(self):
        """
        Return the hyperbolic sine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.sinh()
            0.26480022760227073
        """
        return self._new_c(gsl_ldexp(gsl_sf_expm1(self._value) - gsl_sf_expm1(-self._value), -1))  # (e^x - e^-x)/2

    def tanh(self):
        """
        Return the hyperbolic tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.tanh()
            0.25597778924568454
        """
        return self.sinh() / self.cosh()

    def acosh(self):
        """
        Return the hyperbolic inverse cosine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.cosh(); i
            2.5091784786580567
            sage: abs(i.acosh()-q) < 1e-15
            True
        """
        return self._new_c(gsl_acosh(self._value))

    def arcsinh(self):
        """
        Return the hyperbolic inverse sine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.sinh(); i
            2.3012989023072947
            sage: abs(i.arcsinh()-q) < 1e-15
            True
        """
        return self._new_c(gsl_asinh(self._value))

    def arctanh(self):
        """
        Return the hyperbolic inverse tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.tanh(); i
            0.9171523356672744
            sage: i.arctanh() - q  # rel tol 1
            4.440892098500626e-16
        """
        return self._new_c(gsl_atanh(self._value))

    def sech(self):
        r"""
        Return the hyperbolic secant of ``self``.

        EXAMPLES::

            sage: RDF(pi).sech()
            0.08626673833405443
            sage: CDF(pi).sech()
            0.08626673833405443
        """
        return 1/self.cosh()

    def csch(self):
        r"""
        Return the hyperbolic cosecant of ``self``.

        EXAMPLES::

            sage: RDF(pi).csch()
            0.08658953753004694
            sage: CDF(pi).csch()  # rel tol 1e-15
            0.08658953753004696
        """
        return 1/self.sinh()

    def coth(self):
        r"""
        Return the hyperbolic cotangent of ``self``.

        EXAMPLES::

            sage: RDF(pi).coth()
            1.003741873197321
            sage: CDF(pi).coth()
            1.0037418731973213
        """
        return self.cosh() / self.sinh()

    def erf(self):
        """
        Return the value of the error function on ``self``.

        EXAMPLES::

            sage: RDF(6).erf()
            1.0
        """
        return self._new_c(gsl_sf_erf(self._value))

    @classmethod
    def _factorial(cls, int n):
        """
        Return the factorial of the integer `n` as a real number.

        EXAMPLES::

            sage: RDF.factorial(100)
            9.332621544394415e+157
        """
        if n < 0:
            raise ArithmeticError("n must be nonnegative")
        return cls(gsl_sf_fact(n))

    def gamma(self):
        """
        Return the value of the Euler gamma function on ``self``.

        EXAMPLES::

            sage: RDF(6).gamma()
            120.0
            sage: RDF(1.5).gamma()  # rel tol 1e-15
            0.8862269254527584
        """
        sig_on()
        a = self._new_c(gsl_sf_gamma(self._value))
        sig_off()
        return a

    def zeta(self):
        r"""
        Return the Riemann zeta function evaluated at this real number.

        .. NOTE::

           PARI is vastly more efficient at computing the Riemann zeta
           function. See the example below for how to use it.

        EXAMPLES::

            sage: RDF(2).zeta()  # rel tol 1e-15
            1.6449340668482269
            sage: RDF.pi()^2/6
            1.6449340668482264
            sage: RDF(-2).zeta()
            0.0
            sage: RDF(1).zeta()
            +infinity
        """
        if self._value == 1:
            return self._new_c(1)/self._new_c(0)
        return self._new_c(gsl_sf_zeta(self._value))
