"""
`L`-functions from PARI

This is a wrapper around the general PARI `L`-functions functionality.

AUTHORS:

- Frédéric Chapoton (2018) interface
"""
# ****************************************************************************
#       Copyright (C) 2018 Frédéric Chapoton <chapoton@unistra.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from operator import index as PyNumber_Index
from cypari2.gen import Gen
from sage.libs.pari import pari
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RealField
from sage.rings.complex_mpfr import ComplexField
from sage.rings.power_series_ring import PowerSeriesRing


class lfun_generic():
    r"""
    Create a PARI `L`-function (:pari:`lfun` instance).

    The arguments are::

        lfun_generic(conductor, gammaV, weight, eps, poles, residues, init)

    where

    - ``conductor`` -- integer; the conductor

    - ``gammaV`` -- list of Gamma-factor parameters, e.g. [0] for
      Riemann zeta, [0,1] for ell.curves, (see examples)

    - ``weight`` -- positive real number, usually an integer e.g. 1 for
      Riemann zeta, 2 for `H^1` of curves/`\QQ`

    - ``eps`` -- complex number; sign in functional equation

    - ``poles`` -- (default: ``[]``) list of points where `L^*(s)` has
      (simple) poles; only poles with `Re(s)>weight/2` should be
      included

    - ``residues`` -- vector of residues of `L^*(s)` in those poles or
      set ``residues='automatic'`` (default)

    - ``init`` -- list of coefficients

    RIEMANN ZETA FUNCTION:

    We compute with the Riemann Zeta function::

        sage: from sage.lfunctions.pari import lfun_generic, LFunction
        sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[1])
        sage: lf.init_coeffs([1]*2000)

    Now we can wrap this PARI `L`-function into one Sage `L`-function::

        sage: L = LFunction(lf); L
        L-series of conductor 1 and weight 1
        sage: L(1)
        Traceback (most recent call last):
        ...
        ArithmeticError: pole here
        sage: L(2)
        1.64493406684823
        sage: L.derivative(2)
        -0.937548254315844
        sage: h = RR('0.0000000000001')
        sage: (zeta(2+h) - zeta(2.))/h
        -0.937028232783632
        sage: L.taylor_series(2, k=5)
        1.64493406684823 - 0.937548254315844*z + 0.994640117149451*z^2 - 1.00002430047384*z^3 + 1.00006193307...*z^4 + O(z^5)
    """
    def __init__(self, conductor, gammaV, weight, eps, poles=[],
                 residues='automatic', prec=None, *args, **kwds):
        """
        Initialisation of a :pari:`lfun` from motivic data.

        This can happen either in one stage or in two stages: the coefficients
        can be given using the ``init`` keyword, or entered using
        the :meth:`init_coeffs`.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_generic, LFunction
            sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[1])
        """
        # before entering the coefficients, this attribute is None
        self._L = None

        self.conductor = conductor
        self.gammaV = gammaV
        self.weight = weight
        self.eps = eps
        self.poles = poles
        self.residues = residues
        self.prec = prec

        if args or kwds:
            self.init_coeffs(*args, **kwds)

    def init_coeffs(self, v, cutoff=None, w=1):
        """
        Set the coefficients `a_n` of the `L`-series.

        If `L(s)` is not equal to its dual, pass the coefficients of
        the dual as the second optional argument.

        INPUT:

        - ``v`` -- list of complex numbers or unary function

        - ``cutoff`` -- unused

        - ``w`` -- list of complex numbers or unary function

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_generic, LFunction
            sage: lf = lfun_generic(conductor=1, gammaV=[0, 1], weight=12, eps=1)
            sage: pari_coeffs = pari('k->vector(k,n,(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5)))')
            sage: lf.init_coeffs(pari_coeffs)

        Evaluate the resulting `L`-function at a point, and compare with
        the answer that one gets "by definition" (of `L`-function
        attached to a modular form)::

            sage: L = LFunction(lf)
            sage: L(14)
            0.998583063162746
            sage: a = delta_qexp(1000)
            sage: sum(a[n]/float(n)^14 for n in reversed(range(1,1000)))
            0.9985830631627461

        Illustrate that one can give a list of complex numbers for v
        (see :issue:`10937`)::

            sage: l2 = lfun_generic(conductor=1, gammaV=[0, 1], weight=12, eps=1)
            sage: l2.init_coeffs(list(delta_qexp(1000))[1:])
            sage: L2 = LFunction(l2)
            sage: L2(14)
            0.998583063162746

        TESTS:

        Verify that setting the `w` parameter does not raise an error
        (see :issue:`10937`)::

            sage: L2 = lfun_generic(conductor=1, gammaV=[0, 1], weight=12, eps=1)
            sage: L2.init_coeffs(list(delta_qexp(1000))[1:], w=[1..1000])
        """
        v = pari(v)
        if v.type() not in ('t_CLOSURE', 't_VEC'):
            raise TypeError("v (coefficients) must be a list or a function")

        # w = 0 means a*_n = a_n
        # w = 1 means a*_n = complex conjugate of a_n
        # otherwise w must be a list of coefficients
        w = pari(w)
        if w.type() not in ('t_INT', 't_CLOSURE', 't_VEC'):
            raise TypeError("w (dual coefficients) must be a list or a function or the special value 0 or 1")

        if not self.poles:
            self._L = pari.lfuncreate([v, w, self.gammaV, self.weight,
                                       self.conductor, self.eps])
        else:
            # TODO
            # poles = list(zip(poles, residues))
            self._L = pari.lfuncreate([v, w, self.gammaV, self.weight,
                                       self.conductor, self.eps,
                                       self.poles[0]])

    def __pari__(self):
        """
        Return the PARI `L`-function object.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_generic, LFunction
            sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[1])
            sage: lf.__pari__()
            Traceback (most recent call last):
            ...
            ValueError: call init_coeffs on the L-function first

            sage: lf.init_coeffs([1]*2000)
            sage: X = lf.__pari__()
            sage: X.type()
            't_VEC'
        """
        if self._L is None:
            raise ValueError("call init_coeffs on the L-function first")
        return self._L


def lfun_character(chi):
    """
    Create the `L`-function of a primitive Dirichlet character.

    If the given character is not primitive, it is replaced by its
    associated primitive character.

    OUTPUT: one :pari:`lfun` object

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_character, LFunction
        sage: chi = DirichletGroup(6).gen().primitive_character()
        sage: L = LFunction(lfun_character(chi))
        sage: L(3)
        0.884023811750080

    TESTS:

    A non-primitive one::

        sage: L = LFunction(lfun_character(DirichletGroup(6).gen()))
        sage: L(4)
        0.940025680877124

    With complex arguments::

        sage: from sage.lfunctions.pari import lfun_character, LFunction
        sage: chi = DirichletGroup(6, CC).gen().primitive_character()
        sage: L = LFunction(lfun_character(chi))
        sage: L(3)
        0.884023811750080

    Check the values::

        sage: chi = DirichletGroup(24)([1, -1, -1]); chi
        Dirichlet character modulo 24 of conductor 24
        mapping 7 |--> 1, 13 |--> -1, 17 |--> -1
        sage: Lchi = lfun_character(chi)
        sage: v = [0] + Lchi.lfunan(30).sage()
        sage: all(v[i] == chi(i) for i in (7, 13, 17))
        True
    """
    if not chi.is_primitive():
        chi = chi.primitive_character()
    G, v = chi._pari_init_()
    return pari.lfuncreate([G, v])


def lfun_hgm(motif, t):
    """
    Create the `L`-function of an hypergeometric motive.

    OUTPUT: one :pari:`lfun` object

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_hgm, LFunction
        sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
        sage: H = Hyp(gamma_list=([3,-1,-1,-1]))
        sage: L = LFunction(lfun_hgm(H, 1/5))
        sage: L(3)
        0.901925346034773
    """
    H = pari.hgminit(*motif.alpha_beta())
    lf = pari.lfunhgm(H, t)
    return pari.lfuncreate(lf)


def lfun_elliptic_curve(E):
    """
    Create the `L`-function of an elliptic curve.

    OUTPUT: one :pari:`lfun` object

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_elliptic_curve, LFunction
        sage: E = EllipticCurve('11a1')
        sage: L = LFunction(lfun_elliptic_curve(E))
        sage: L(3)
        0.752723147389051
        sage: L(1)
        0.253841860855911

    Over number fields::

        sage: K.<a> = QuadraticField(2)
        sage: E = EllipticCurve([1, a])
        sage: L = LFunction(lfun_elliptic_curve(E))
        sage: L(3)
        1.00412346717019
    """
    return pari.lfuncreate(E)


def lfun_number_field(K):
    """
    Create the Dedekind zeta function of a number field.

    OUTPUT: one :pari:`lfun` object

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_number_field, LFunction

        sage: L = LFunction(lfun_number_field(QQ))
        sage: L(3)
        1.20205690315959

        sage: K = NumberField(x**2 - 2, 'a')
        sage: L = LFunction(lfun_number_field(K))
        sage: L(3)
        1.15202784126080
        sage: L(0)
        ...0.000000000000000
    """
    return pari.lfuncreate(K)


def lfun_eta_quotient(scalings, exponents):
    """
    Return the `L`-function of an eta-quotient.

    This uses :pari:`lfunetaquo`.

    INPUT:

    - ``scalings`` -- list of integers; the scaling factors

    - ``exponents`` -- list of integers; the exponents

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_eta_quotient, LFunction
        sage: L = LFunction(lfun_eta_quotient([1], [24]))
        sage: L(1)
        0.0374412812685155

        sage: lfun_eta_quotient([6], [4])
        [[Vecsmall([7]), [Vecsmall([6]), Vecsmall([4]), 0]], 0, [0, 1], 2, 36, 1]

        sage: lfun_eta_quotient([2, 1, 4], [5, -2, -2])
        Traceback (most recent call last):
        ...
        PariError: sorry, noncuspidal eta quotient is not yet implemented
    """
    from sage.matrix.constructor import matrix
    N = len(scalings)
    if N != len(exponents):
        raise ValueError('arguments should have the same length')
    m = matrix(ZZ, N, 2, list(zip(scalings, exponents)))
    return pari.lfunetaquo(m)


def lfun_delta():
    """
    Return the `L`-function of Ramanujan's Delta modular form.

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_delta, LFunction
        sage: L = LFunction(lfun_delta())
        sage: L(1)
        0.0374412812685155
    """
    return lfun_eta_quotient([1], [24])


def lfun_quadratic_form(qf):
    """
    Return the `L`-function of a positive definite quadratic form.

    This uses :pari:`lfunqf`.

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_quadratic_form, LFunction
        sage: Q = QuadraticForm(ZZ, 2, [2, 3, 4])
        sage: L = LFunction(lfun_quadratic_form(Q))
        sage: L(3)
        0.377597233183583
    """
    if not qf.is_positive_definite():
        raise ValueError('quadratic form must be positive definite')
    return pari.lfunqf(qf.matrix())


def lfun_genus2(C):
    """
    Return the `L`-function of a curve of genus 2.

    INPUT:

    - ``C`` -- hyperelliptic curve of genus 2

    Currently, the model needs to be minimal at 2.

    This uses :pari:`lfungenus2`.

    EXAMPLES::

        sage: from sage.lfunctions.pari import lfun_genus2, LFunction
        sage: x = polygen(QQ, 'x')
        sage: C = HyperellipticCurve(x^5 + x + 1)
        sage: L = LFunction(lfun_genus2(C))
        ...
        sage: L(3)
        0.965946926261520

        sage: C = HyperellipticCurve(x^2 + x, x^3 + x^2 + 1)
        sage: L = LFunction(lfun_genus2(C))
        sage: L(2)
        0.364286342944359

    TESTS::

        sage: x = polygen(QQ, 'x')
        sage: H = HyperellipticCurve(x^15 + x + 1)
        sage: L = LFunction(lfun_genus2(H))
        Traceback (most recent call last):
        ...
        ValueError: curve must be hyperelliptic of genus 2
    """
    from sage.schemes.hyperelliptic_curves.hyperelliptic_g2 import HyperellipticCurve_g2 as hyp_g2
    if not isinstance(C, hyp_g2):
        raise ValueError('curve must be hyperelliptic of genus 2')
    P, Q = C.hyperelliptic_polynomials()
    return pari.lfungenus2(P) if not Q else pari.lfungenus2([P, Q])


class LFunction(SageObject):
    r"""
    Build the `L`-function from a PARI `L`-function.

    .. RUBRIC:: Rank 1 elliptic curve

    We compute with the `L`-series of a rank `1` curve. ::

        sage: E = EllipticCurve('37a')
        sage: L = E.lseries().dokchitser(algorithm='pari'); L
        PARI L-function associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: L(1)
        0.000000000000000
        sage: L.derivative(1)
        0.305999773834052
        sage: L.derivative(1, 2)
        0.373095594536324
        sage: L.num_coeffs()
        50
        sage: L.taylor_series(1, 4)
        0.000000000000000 + 0.305999773834052*z + 0.186547797268162*z^2 - 0.136791463097188*z^3 + O(z^4)
        sage: L.check_functional_equation()  # abs tol 4e-19
        1.08420217248550e-19

    .. RUBRIC:: Rank 2 elliptic curve

    We compute the leading coefficient and Taylor expansion of the
    `L`-series of a rank `2` elliptic curve::

        sage: E = EllipticCurve('389a')
        sage: L = E.lseries().dokchitser(algorithm='pari')
        sage: L.num_coeffs()
        163
        sage: L.derivative(1, E.rank())
        1.51863300057685
        sage: L.taylor_series(1, 4)
        ...e-19 + (...e-19)*z + 0.759316500288427*z^2 - 0.430302337583362*z^3 + O(z^4)

    .. RUBRIC:: Number field

    We compute with the Dedekind zeta function of a number field::

        sage: x = var('x')
        sage: K = NumberField(x**4 - x**2 - 1,'a')
        sage: L = K.zeta_function(algorithm='pari')
        sage: L.conductor
        400
        sage: L.num_coeffs()
        348
        sage: L(2)
        1.10398438736918
        sage: L.taylor_series(2, 3)
        1.10398438736918 - 0.215822638498759*z + 0.279836437522536*z^2 + O(z^3)

    .. RUBRIC:: Ramanujan `\Delta` L-function

    The coefficients are given by Ramanujan's tau function::

        sage: from sage.lfunctions.pari import lfun_generic, LFunction
        sage: lf = lfun_generic(conductor=1, gammaV=[0, 1], weight=12, eps=1)
        sage: tau = pari('k->vector(k,n,(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5)))')
        sage: lf.init_coeffs(tau)
        sage: L = LFunction(lf)

    Now we are ready to evaluate, etc. ::

        sage: L(1)
        0.0374412812685155
        sage: L.taylor_series(1, 3)
        0.0374412812685155 + 0.0709221123619322*z + 0.0380744761270520*z^2 + O(z^3)
    """
    def __init__(self, lfun, prec=None):
        """
        Initialization of the `L`-function from a PARI `L`-function.

        INPUT:

        - ``lfun`` -- a PARI :pari:`lfun` object or an instance of :class:`lfun_generic`
        - ``prec`` -- integer (default: 53); number of *bits* of precision

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_generic, LFunction
            sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], v=pari('k->vector(k,n,1)'))
            sage: L = LFunction(lf)
            sage: L.num_coeffs()
            4
        """
        if isinstance(lfun, lfun_generic):
            # preparation using motivic data
            self._L = lfun.__pari__()
            if prec is None:
                prec = lfun.prec
        elif isinstance(lfun, Gen):
            # already some PARI lfun
            self._L = lfun
        else:
            # create a PARI lfunction from other input data
            self._L = pari.lfuncreate(lfun)

        self._conductor = ZZ(self._L[4])  # needs check
        self._weight = ZZ(self._L[3])  # needs check
        if prec is None:
            self.prec = 53
        else:
            self.prec = PyNumber_Index(prec)
        self._RR = RealField(self.prec)
        self._CC = ComplexField(self.prec)
        # Complex field used for inputs, which ensures exact 1-to-1
        # conversion to/from PARI. Since PARI objects have a precision
        # in machine words (not bits), this is typically higher. For
        # example, the default of 53 bits of precision would become 64.
        self._CCin = ComplexField(pari.bitprecision(self._RR(1)))

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: L = LFunction(lfun_number_field(QQ)); L
            L-series of conductor 1 and weight 1
        """
        return "L-series of conductor %s and weight %s" % (self._conductor,
                                                           self._weight)

    @property
    def conductor(self):
        """
        Return the conductor.

        EXAMPLES::

            sage: from sage.lfunctions.pari import *
            sage: L = LFunction(lfun_number_field(QQ)); L.conductor
            1
            sage: E = EllipticCurve('11a')
            sage: L = LFunction(lfun_number_field(E)); L.conductor
            11
        """
        return self._conductor

    def num_coeffs(self, T=1):
        """
        Return number of coefficients `a_n` that are needed in
        order to perform most relevant `L`-function computations to
        the desired precision.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: L = E.lseries().dokchitser(algorithm='pari')
            sage: L.num_coeffs()
            27
            sage: E = EllipticCurve('5077a')
            sage: L = E.lseries().dokchitser(algorithm='pari')
            sage: L.num_coeffs()
            591

            sage: from sage.lfunctions.pari import lfun_generic, LFunction
            sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], v=pari('k->vector(k,n,1)'))
            sage: L = LFunction(lf)
            sage: L.num_coeffs()
            4
        """
        # lfuncost takes a domain
        domain = [pari(1) / 2, 1, 1]
        return ZZ(pari.lfuncost(self._L, domain)[0])

    def Lambda(self, s):
        """
        Evaluate the completed `L`-function at s.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: L = LFunction(lfun_number_field(QQ))
            sage: L.Lambda(2)
            0.523598775598299
            sage: L.Lambda(1 - 2)
            0.523598775598299
        """
        s = self._CCin(s)
        R = self._CC
        return R(pari.lfunlambda(self._L, s, precision=self.prec))

    def hardy(self, t):
        """
        Evaluate the associated Hardy function at t.

        This only works for real t.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: L = LFunction(lfun_number_field(QQ))
            sage: L.hardy(2)
            -0.962008487244041

        TESTS::

            sage: L.hardy(.4 + .3*I)
            Traceback (most recent call last):
            ...
            PariError: incorrect type in lfunhardy (t_COMPLEX)
        """
        t = self._CCin(t)
        R = self._RR
        return R(pari.lfunhardy(self._L, t, precision=self.prec))

    def derivative(self, s, D=1):
        """
        Return the derivative of the `L`-function at point s and order D.

        INPUT:

        - ``s`` -- complex number

        - ``D`` -- integer (default: 1)

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: L = LFunction(lfun_number_field(QQ))
            sage: L.derivative(2)
            -0.937548254315844
        """
        s = self._CCin(s)
        R = self._CC
        return R(pari.lfun(self._L, s, D, precision=self.prec))

    def taylor_series(self, s, k=6, var='z'):
        r"""
        Return the first `k` terms of the Taylor series expansion
        of the `L`-series about `s`.

        This is returned as a formal power series in ``var``.

        INPUT:

        - ``s`` -- complex number; point about which to expand

        - ``k`` -- integer (default: 6); series is `O(``var``^k)`

        - ``var`` -- string (default: ``'z'``); variable of power series

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: lf = lfun_number_field(QQ)
            sage: L = LFunction(lf)
            sage: L.taylor_series(2, 3)
            1.64493406684823 - 0.937548254315844*z + 0.994640117149451*z^2 + O(z^3)
            sage: E = EllipticCurve('37a')
            sage: L = E.lseries().dokchitser(algorithm='pari')
            sage: L.taylor_series(1)
            0.000000000000000 + 0.305999773834052*z + 0.186547797268162*z^2 - 0.136791463097188*z^3 + 0.0161066468496401*z^4 + 0.0185955175398802*z^5 + O(z^6)

        We compute a Taylor series where each coefficient is to high
        precision::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries().dokchitser(200,algorithm='pari')
            sage: L.taylor_series(1, 3)
            2...e-63 + (...e-63)*z + 0.75931650028842677023019260789472201907809751649492435158581*z^2 + O(z^3)

        Check that :issue:`25402` is fixed::

            sage: L = EllipticCurve("24a1").modular_form().lseries()
            sage: L.taylor_series(-1, 3)
            0.000000000000000 - 0.702565506265199*z + 0.638929001045535*z^2 + O(z^3)
        """
        pt = pari.Ser([s, 1], d=k)  # s + x + O(x^k)
        B = PowerSeriesRing(self._CC, var)
        return B(pari.lfun(self._L, pt, precision=self.prec))

    def zeros(self, maxi):
        """
        Return the zeros with imaginary part bounded by maxi.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: lf = lfun_number_field(QQ)
            sage: L = LFunction(lf)
            sage: L.zeros(20)
            [14.1347251417347]
        """
        R = self._CC
        return [R(z) for z in pari.lfunzeros(self._L, maxi)]

    def _clear_value_cache(self):
        """
        Clear the cache where values of the function are stored.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_number_field, LFunction
            sage: lf = lfun_number_field(QQ)
            sage: L = LFunction(lf)
            sage: L(4)
            1.08232323371114
            sage: L._clear_value_cache()
        """
        del self.__values

    def __call__(self, s):
        r"""
        Return the value of the `L`-function at point ``s``.

        INPUT:

        - ``s`` -- complex number

        .. NOTE::

            Evaluation of the function takes a long time, so each
            evaluation is cached. Call :meth:`_clear_value_cache` to
            clear the evaluation cache.

        EXAMPLES::

            sage: E = EllipticCurve('5077a')
            sage: L = E.lseries().dokchitser(100, algorithm='pari')
            sage: L(1)
            0.00000000000000000000000000000
            sage: L(1 + I)
            -1.3085436607849493358323930438 + 0.81298000036784359634835412129*I
        """
        s = self._CC(s)
        try:
            return self.__values[s]
        except AttributeError:
            self.__values = {}
        except KeyError:
            pass
        try:
            value = self._CC(pari.lfun(self._L, s, precision=self.prec))
        except NameError:
            raise ArithmeticError('pole here')
        else:
            return value

    def check_functional_equation(self):
        r"""
        Verify how well numerically the functional equation is satisfied.

        If what this function returns does not look like 0 at all,
        probably the functional equation is wrong, i.e. some of the
        parameters gammaV, conductor, etc., or the coefficients are wrong.

        EXAMPLES::

            sage: from sage.lfunctions.pari import lfun_generic, LFunction
            sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[1], v=pari('k->vector(k,n,1)'))
            sage: L = LFunction(lf)
            sage: L.check_functional_equation()
            4.33680868994202e-19

        If we choose the sign in functional equation for the
        `\zeta` function incorrectly, the functional equation
        does not check out::

            sage: lf = lfun_generic(conductor=1, gammaV=[0], weight=1, eps=-1, poles=[1], residues=[1])
            sage: lf.init_coeffs([1]*2000)
            sage: L = LFunction(lf)
            sage: L.check_functional_equation()
            16.0000000000000
        """
        return self._RR(2)**pari.lfuncheckfeq(self._L)
