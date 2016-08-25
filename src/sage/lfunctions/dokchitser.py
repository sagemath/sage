"""
Dokchitser's L-functions Calculator

AUTHORS:

- Tim Dokchitser (2002): original PARI code and algorithm (and the
  documentation below is based on Dokchitser's docs).

- William Stein (2006-03-08): Sage interface

TODO:

- add more examples from SAGE_EXTCODE/pari/dokchitser that illustrate
  use with Eisenstein series, number fields, etc.

- plug this code into number fields and modular forms code (elliptic
  curves are done).
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function

import copy

from sage.structure.sage_object import SageObject
from sage.rings.all import ComplexField, Integer
from sage.misc.all import verbose, sage_eval
import sage.interfaces.gp

class Dokchitser(SageObject):
    r"""
    Dokchitser's `L`-functions Calculator

    Create a Dokchitser `L`-series with

    Dokchitser(conductor, gammaV, weight, eps, poles, residues, init,
    prec)

    where

    - ``conductor`` - integer, the conductor

    - ``gammaV`` - list of Gamma-factor parameters, e.g. [0] for
      Riemann zeta, [0,1] for ell.curves, (see examples).

    - ``weight`` - positive real number, usually an integer e.g. 1 for
      Riemann zeta, 2 for `H^1` of curves/`\QQ`

    - ``eps`` - complex number; sign in functional equation

    - ``poles`` - (default: []) list of points where `L^*(s)` has
      (simple) poles; only poles with `Re(s)>weight/2` should be
      included

    - ``residues`` - vector of residues of `L^*(s)` in those poles or
      set residues='automatic' (default value)

    - ``prec`` - integer (default: 53) number of *bits* of precision

    RIEMANN ZETA FUNCTION:

    We compute with the Riemann Zeta function.

    ::

        sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
        sage: L
        Dokchitser L-series of conductor 1 and weight 1
        sage: L(1)
        Traceback (most recent call last):
        ...
        ArithmeticError
        sage: L(2)
        1.64493406684823
        sage: L(2, 1.1)
        1.64493406684823
        sage: L.derivative(2)
        -0.937548254315844
        sage: h = RR('0.0000000000001')
        sage: (zeta(2+h) - zeta(2.))/h
        -0.937028232783632
        sage: L.taylor_series(2, k=5)
        1.64493406684823 - 0.937548254315844*z + 0.994640117149451*z^2 - 1.00002430047384*z^3 + 1.00006193307...*z^4 + O(z^5)

    RANK 1 ELLIPTIC CURVE:

    We compute with the `L`-series of a rank `1`
    curve.

    ::

        sage: E = EllipticCurve('37a')
        sage: L = E.lseries().dokchitser(); L
        Dokchitser L-function associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: L(1)
        0.000000000000000
        sage: L.derivative(1)
        0.305999773834052
        sage: L.derivative(1,2)
        0.373095594536324
        sage: L.num_coeffs()
        48
        sage: L.taylor_series(1,4)
        0.000000000000000 + 0.305999773834052*z + 0.186547797268162*z^2 - 0.136791463097188*z^3 + O(z^4)
        sage: L.check_functional_equation()
        6.11218974700000e-18                            # 32-bit
        6.04442711160669e-18                            # 64-bit

    RANK 2 ELLIPTIC CURVE:

    We compute the leading coefficient and Taylor expansion of the
    `L`-series of a rank `2` curve.

    ::

        sage: E = EllipticCurve('389a')
        sage: L = E.lseries().dokchitser()
        sage: L.num_coeffs ()
        156
        sage: L.derivative(1,E.rank())
        1.51863300057685
        sage: L.taylor_series(1,4)
        -1.27685190980159e-23 + (7.23588070754027e-24)*z + 0.759316500288427*z^2 - 0.430302337583362*z^3 + O(z^4)  # 32-bit
        -2.72911738151096e-23 + (1.54658247036311e-23)*z + 0.759316500288427*z^2 - 0.430302337583362*z^3 + O(z^4)  # 64-bit

    RAMANUJAN DELTA L-FUNCTION:

    The coefficients are given by Ramanujan's tau function::

        sage: L = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
        sage: pari_precode = 'tau(n)=(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5))'
        sage: L.init_coeffs('tau(k)', pari_precode=pari_precode)

    We redefine the default bound on the coefficients: Deligne's
    estimate on tau(n) is better than the default
    coefgrow(n)=`(4n)^{11/2}` (by a factor 1024), so
    re-defining coefgrow() improves efficiency (slightly faster).

    ::

        sage: L.num_coeffs()
        12
        sage: L.set_coeff_growth('2*n^(11/2)')
        sage: L.num_coeffs()
        11

    Now we're ready to evaluate, etc.

    ::

        sage: L(1)
        0.0374412812685155
        sage: L(1, 1.1)
        0.0374412812685155
        sage: L.taylor_series(1,3)
        0.0374412812685155 + 0.0709221123619322*z + 0.0380744761270520*z^2 + O(z^3)
    """
    def __init__(self, conductor, gammaV, weight, eps, \
                       poles=[], residues='automatic', prec=53,
                       init=None):
        """
        Initialization of Dokchitser calculator EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.num_coeffs()
            4
        """
        self.conductor = conductor
        self.gammaV = gammaV
        self.weight = weight
        self.eps    = eps
        self.poles  = poles
        self.residues = residues
        self.prec   = prec
        self.__CC   = ComplexField(self.prec)
        self.__RR   = self.__CC._real_field()
        if not init is None:
            self.init_coeffs(init)
            self.__init = init
        else:
            self.__init = False

    def __reduce__(self):
        D = copy.copy(self.__dict__)
        if '_Dokchitser__gp' in D:
            del D['_Dokchitser__gp']
        return reduce_load_dokchitser, (D, )

    def _repr_(self):
        z = "Dokchitser L-series of conductor %s and weight %s"%(
                   self.conductor, self.weight)
        return z

    def __del__(self):
        self.gp().quit()

    def gp(self):
        """
        Return the gp interpreter that is used to implement this Dokchitser
        L-function.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: L = E.lseries().dokchitser()
            sage: L(2)
            0.546048036215014
            sage: L.gp()
            PARI/GP interpreter
        """
        try:
            return self.__gp
        except AttributeError:
            logfile = None
            # For debugging
            import os
            from sage.env import DOT_SAGE
            logfile = os.path.join(DOT_SAGE, 'dokchitser.log')
            g = sage.interfaces.gp.Gp(script_subdirectory='dokchitser', logfile=logfile)
            g.read('computel.gp')
            self.__gp = g
            self._gp_eval('default(realprecision, %s)'%(self.prec//3 + 2))
            self._gp_eval('conductor = %s'%self.conductor)
            self._gp_eval('gammaV = %s'%self.gammaV)
            self._gp_eval('weight = %s'%self.weight)
            self._gp_eval('sgn = %s'%self.eps)
            self._gp_eval('Lpoles = %s'%self.poles)
            self._gp_eval('Lresidues = %s'%self.residues)
            g._dokchitser = True
            return g

    def _gp_eval(self, s):
        try:
            t = self.gp().eval(s)
        except (RuntimeError, TypeError):
            raise RuntimeError("Unable to create L-series, due to precision or other limits in PARI.")
        if '***' in t:
            raise RuntimeError("Unable to create L-series, due to precision or other limits in PARI.")
        return t

    def __check_init(self):
        if not self.__init:
            raise ValueError("you must call init_coeffs on the L-function first")

    def num_coeffs(self, T=1):
        """
        Return number of coefficients `a_n` that are needed in
        order to perform most relevant `L`-function computations to
        the desired precision.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: L = E.lseries().dokchitser()
            sage: L.num_coeffs()
            26
            sage: E = EllipticCurve('5077a')
            sage: L = E.lseries().dokchitser()
            sage: L.num_coeffs()
            568
            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.num_coeffs()
            4
        """
        return Integer(self.gp().eval('cflength(%s)'%T))

    def init_coeffs(self, v, cutoff=1,
                             w=None,
                             pari_precode='',
                             max_imaginary_part=0,
                             max_asymp_coeffs=40):
        """
        Set the coefficients `a_n` of the `L`-series. If
        `L(s)` is not equal to its dual, pass the coefficients of
        the dual as the second optional argument.

        INPUT:


        -  ``v`` - list of complex numbers or string (pari
           function of k)

        -  ``cutoff`` - real number = 1 (default: 1)

        -  ``w`` - list of complex numbers or string (pari
           function of k)

        -  ``pari_precode`` - some code to execute in pari
           before calling initLdata

        -  ``max_imaginary_part`` - (default: 0): redefine if
           you want to compute L(s) for s having large imaginary part,

        -  ``max_asymp_coeffs`` - (default: 40): at most this
           many terms are generated in asymptotic series for phi(t) and
           G(s,t).


        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: pari_precode = 'tau(n)=(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5))'
            sage: L.init_coeffs('tau(k)', pari_precode=pari_precode)

        Evaluate the resulting L-function at a point, and compare with the answer that
        one gets "by definition" (of L-function attached to a modular form)::

            sage: L(14)
            0.998583063162746
            sage: a = delta_qexp(1000)
            sage: sum(a[n]/float(n)^14 for n in range(1,1000))
            0.9985830631627459

        Illustrate that one can give a list of complex numbers for v
        (see :trac:`10937`)::

            sage: L2 = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: L2.init_coeffs(list(delta_qexp(1000))[1:])
            sage: L2(14)
            0.998583063162746

        TESTS:

        Verify that setting the `w` parameter does not raise an error
        (see :trac:`10937`).  Note that the meaning of `w` does not seem to
        be documented anywhere in Dokchitser's package yet, so there is
        no claim that the example below is meaningful! ::

            sage: L2 = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: L2.init_coeffs(list(delta_qexp(1000))[1:], w=[1..1000])
        """
        if isinstance(v, tuple) and w is None:
            v, cutoff, w, pari_precode, max_imaginary_part, max_asymp_coeffs = v

        self.__init = (v, cutoff, w, pari_precode, max_imaginary_part, max_asymp_coeffs)
        gp = self.gp()
        if pari_precode != '':
            self._gp_eval(pari_precode)
        RR = self.__CC._real_field()
        cutoff = RR(cutoff)
        if isinstance(v, str):
            if w is None:
                self._gp_eval('initLdata("%s", %s)'%(v, cutoff))
                return
            self._gp_eval('initLdata("%s",%s,"%s")'%(v,cutoff,w))
            return
        if not isinstance(v, (list, tuple)):
            raise TypeError("v (=%s) must be a list, tuple, or string"%v)
        CC = self.__CC
        v = ','.join([CC(a)._pari_init_() for a in v])
        self._gp_eval('Avec = [%s]'%v)
        if w is None:
            self._gp_eval('initLdata("Avec[k]", %s)'%cutoff)
            return
        w = ','.join([CC(a)._pari_init_() for a in w])
        self._gp_eval('Bvec = [%s]'%w)
        self._gp_eval('initLdata("Avec[k]",%s,"Bvec[k]")'%cutoff)

    def __to_CC(self, s):
        s = s.replace('.E','.0E').replace(' ','')
        return self.__CC(sage_eval(s, locals={'I':self.__CC.gen(0)}))

    def _clear_value_cache(self):
        del self.__values

    def __call__(self, s, c=None):
        r"""
        INPUT:

        -  ``s`` - complex number

        .. note::

           Evaluation of the function takes a long time, so each
           evaluation is cached. Call ``self._clear_value_cache()`` to
           clear the evaluation cache.

        EXAMPLES::

            sage: E = EllipticCurve('5077a')
            sage: L = E.lseries().dokchitser(100)
            sage: L(1)
            0.00000000000000000000000000000
            sage: L(1+I)
            -1.3085436607849493358323930438 + 0.81298000036784359634835412129*I
        """
        self.__check_init()
        s = self.__CC(s)
        try:
            return self.__values[s]
        except AttributeError:
            self.__values = {}
        except KeyError:
            pass
        z = self.gp().eval('L(%s)'%s)
        if 'pole' in z:
            print(z)
            raise ArithmeticError
        elif '***' in z:
            print(z)
            raise RuntimeError
        elif 'Warning' in z:
            i = z.rfind('\n')
            msg = z[:i].replace('digits','decimal digits')
            verbose(msg, level=-1)
            ans = self.__to_CC(z[i+1:])
            self.__values[s] = ans
            return ans
        ans = self.__to_CC(z)
        self.__values[s] = ans
        return ans

    def derivative(self, s, k=1):
        r"""
        Return the `k`-th derivative of the `L`-series at
        `s`.

        .. warning::

           If `k` is greater than the order of vanishing of
           `L` at `s` you may get nonsense.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries().dokchitser()
            sage: L.derivative(1,E.rank())
            1.51863300057685
        """
        self.__check_init()
        s = self.__CC(s)
        k = Integer(k)
        z = self.gp().eval('L(%s,,%s)'%(s,k))
        if 'pole' in z:
            raise ArithmeticError(z)
        elif 'Warning' in z:
            i = z.rfind('\n')
            msg = z[:i].replace('digits','decimal digits')
            verbose(msg, level=-1)
            return self.__CC(z[i:])
        return self.__CC(z)


    def taylor_series(self, a=0, k=6, var='z'):
        r"""
        Return the first `k` terms of the Taylor series expansion
        of the `L`-series about `a`.

        This is returned as a series in ``var``, where you
        should view ``var`` as equal to `s-a`. Thus
        this function returns the formal power series whose coefficients
        are `L^{(n)}(a)/n!`.

        INPUT:


        -  ``a`` - complex number (default: 0); point about
           which to expand

        -  ``k`` - integer (default: 6), series is
           `O(``var``^k)`

        -  ``var`` - string (default: 'z'), variable of power
           series


        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.taylor_series(2, 3)
            1.64493406684823 - 0.937548254315844*z + 0.994640117149451*z^2 + O(z^3)
            sage: E = EllipticCurve('37a')
            sage: L = E.lseries().dokchitser()
            sage: L.taylor_series(1)
            0.000000000000000 + 0.305999773834052*z + 0.186547797268162*z^2 - 0.136791463097188*z^3 + 0.0161066468496401*z^4 + 0.0185955175398802*z^5 + O(z^6)

        We compute a Taylor series where each coefficient is to high
        precision.

        ::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries().dokchitser(200)
            sage: L.taylor_series(1,3)
            -9.094...e-82 + (5.1538...e-82)*z + 0.75931650028842677023019260789472201907809751649492435158581*z^2 + O(z^3)
        """
        self.__check_init()
        a = self.__CC(a)
        k = Integer(k)
        try:
            z = self.gp()('Vec(Lseries(%s,,%s))'%(a,k-1))
        except TypeError as msg:
            raise RuntimeError("%s\nUnable to compute Taylor expansion (try lowering the number of terms)"%msg)
        r = repr(z)
        if 'pole' in r:
            raise ArithmeticError(r)
        elif 'Warning' in r:
            i = r.rfind('\n')
            msg = r[:i].replace('digits','decimal digits')
            verbose(msg, level=-1)
        v = list(z)
        K = self.__CC
        v = [K(repr(x)) for x in v]
        R = self.__CC[[var]]
        return R(v,len(v))

    def check_functional_equation(self, T=1.2):
        r"""
        Verifies how well numerically the functional equation is satisfied,
        and also determines the residues if ``self.poles !=
        []`` and residues='automatic'.

        More specifically: for `T>1` (default 1.2),
        ``self.check_functional_equation(T)`` should ideally
        return 0 (to the current precision).


        -  if what this function returns does not look like 0 at all,
           probably the functional equation is wrong (i.e. some of the
           parameters gammaV, conductor etc., or the coefficients are wrong),

        -  if checkfeq(T) is to be used, more coefficients have to be
           generated (approximately T times more), e.g. call cflength(1.3),
           initLdata("a(k)",1.3), checkfeq(1.3)

        -  T=1 always (!) returns 0, so T has to be away from 1

        -  default value `T=1.2` seems to give a reasonable
           balance

        -  if you don't have to verify the functional equation or the
           L-values, call num_coeffs(1) and initLdata("a(k)",1), you need
           slightly less coefficients.


        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.check_functional_equation()
            -1.35525271600000e-20                        # 32-bit
            -2.71050543121376e-20                        # 64-bit

        If we choose the sign in functional equation for the
        `\zeta` function incorrectly, the functional equation
        doesn't check out.

        ::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=-11, poles=[1], residues=[-1], init='1')
            sage: L.check_functional_equation()
            -9.73967861488124
        """
        self.__check_init()
        z = self.gp().eval('checkfeq(%s)'%T).replace(' ','')
        return self.__CC(z)

    def set_coeff_growth(self, coefgrow):
        r"""
        You might have to redefine the coefficient growth function if the
        `a_n` of the `L`-series are not given by the
        following PARI function::

                        coefgrow(n) = if(length(Lpoles),
                                          1.5*n^(vecmax(real(Lpoles))-1),
                                          sqrt(4*n)^(weight-1));


        INPUT:


        -  ``coefgrow`` - string that evaluates to a PARI
           function of n that defines a coefgrow function.


        EXAMPLE::

            sage: L = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: pari_precode = 'tau(n)=(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5))'
            sage: L.init_coeffs('tau(k)', pari_precode=pari_precode)
            sage: L.set_coeff_growth('2*n^(11/2)')
            sage: L(1)
            0.0374412812685155
        """
        if not isinstance(coefgrow, str):
            raise TypeError("coefgrow must be a string")
        g = self.gp()
        g.eval('coefgrow(n) = %s'%(coefgrow.replace('\n',' ')))


def reduce_load_dokchitser(D):
    X = Dokchitser(1,1,1,1)
    X.__dict__ = D
    X.init_coeffs(X._Dokchitser__init)
    return X
