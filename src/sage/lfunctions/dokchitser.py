"""
Dokchitser's `L`-functions calculator

.. TODO::

    - add more examples from SAGE_EXTCODE/pari/dokchitser that illustrate
      use with Eisenstein series, number fields, etc.

    - plug this code into number fields and modular forms code (elliptic
      curves are done).

AUTHORS:

- Tim Dokchitser (2002): original PARI code and algorithm (and the
  documentation below is based on Dokchitser's docs).

- William Stein (2006-03-08): Sage interface
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import copy
import os
import re
import string

from sage.structure.sage_object import SageObject
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer import Integer
from sage.misc.verbose import verbose
import sage.interfaces.gp
from sage.env import SAGE_EXTCODE


class Dokchitser(SageObject):
    r"""
    Dokchitser's `L`-functions Calculator.

    PARI code can be found on
    `Dokchitser's homepage <https://people.maths.bris.ac.uk/~matyd/computel>`_.

    Create a Dokchitser `L`-series with

    Dokchitser(conductor, gammaV, weight, eps, poles, residues, init,
    prec)

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

    - ``prec`` -- integer (default: 53); number of *bits* of precision

    RIEMANN ZETA FUNCTION:

    We compute with the Riemann Zeta function. ::

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

    We compute with the `L`-series of a rank `1` curve. ::

        sage: E = EllipticCurve('37a')
        sage: L = E.lseries().dokchitser(algorithm='gp'); L
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
        sage: L.check_functional_equation()  # abs tol 1e-19
        6.04442711160669e-18

    RANK 2 ELLIPTIC CURVE:

    We compute the leading coefficient and Taylor expansion of the
    `L`-series of a rank `2` elliptic curve. ::

        sage: E = EllipticCurve('389a')
        sage: L = E.lseries().dokchitser(algorithm='gp')
        sage: L.num_coeffs()
        156
        sage: L.derivative(1,E.rank())
        1.51863300057685
        sage: L.taylor_series(1,4)
        -1.27685190980159e-23 + (7.23588070754027e-24)*z + 0.759316500288427*z^2 - 0.430302337583362*z^3 + O(z^4)  # 32-bit
        -2.72911738151096e-23 + (1.54658247036311e-23)*z + 0.759316500288427*z^2 - 0.430302337583362*z^3 + O(z^4)  # 64-bit

    NUMBER FIELD:

    We compute with the Dedekind zeta function of a number field. ::

        sage: x = var('x')
        sage: K = NumberField(x**4 - x**2 - 1,'a')
        sage: L = K.zeta_function(algorithm='gp')
        sage: L.conductor
        400
        sage: L.num_coeffs()
        264
        sage: L(2)
        1.10398438736918
        sage: L.taylor_series(2,3)
        1.10398438736918 - 0.215822638498759*z + 0.279836437522536*z^2 + O(z^3)

    RAMANUJAN DELTA L-FUNCTION:

    The coefficients are given by Ramanujan's tau function::

        sage: L = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
        sage: pari_precode = 'tau(n)=(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5))'
        sage: L.init_coeffs('tau(k)', pari_precode=pari_precode)

    We redefine the default bound on the coefficients: Deligne's
    estimate on tau(n) is better than the default
    coefgrow(n)= `(4n)^{11/2}` (by a factor 1024), so
    re-defining coefgrow() improves efficiency (slightly faster). ::

        sage: L.num_coeffs()
        12
        sage: L.set_coeff_growth('2*n^(11/2)')
        sage: L.num_coeffs()
        11

    Now we're ready to evaluate, etc. ::

        sage: L(1)
        0.0374412812685155
        sage: L(1, 1.1)
        0.0374412812685155
        sage: L.taylor_series(1,3)
        0.0374412812685155 + 0.0709221123619322*z + 0.0380744761270520*z^2 + O(z^3)
    """

    __gp = None
    __globals = set()
    # set of global variables defined in a run of the
    # computel.gp script that are replaced by indexed copies
    # in the computel.gp.template

    __globals_re = None
    __instance = 0  # Monotonically increasing unique instance ID
    __n_instances = 0  # Number of currently allocated instances
    __template_filename = os.path.join(SAGE_EXTCODE, 'pari', 'dokchitser',
                                       'computel.gp.template')
    __init = False

    def __new__(cls, *args, **kwargs):
        inst = super().__new__(cls, *args, **kwargs)
        inst.__instance = cls.__instance
        cls.__n_instances += 1
        cls.__instance += 1
        return inst

    def __init__(self, conductor, gammaV, weight, eps,
                 poles=None, residues='automatic', prec=53,
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
        self.eps = eps
        self.poles = poles if poles is not None else []
        self.residues = residues
        self.prec = prec
        self.__CC = ComplexField(self.prec)
        self.__RR = self.__CC._real_field()
        self.__initialized = False
        if init is not None:
            self.init_coeffs(init)

    def __reduce__(self):
        D = copy.copy(self.__dict__)
        if '_Dokchitser__gp' in D:
            del D['_Dokchitser__gp']
        return reduce_load_dokchitser, (D, )

    def _repr_(self):
        return "Dokchitser L-series of conductor %s and weight %s" % (
            self.conductor, self.weight)

    def __del__(self):
        self._teardown_gp(self.__instance)

    def gp(self):
        """
        Return the gp interpreter that is used to implement this Dokchitser
        `L`-function.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: L = E.lseries().dokchitser(algorithm='gp')
            sage: L(2)
            0.546048036215014
            sage: L.gp()
            PARI/GP interpreter
        """
        if self.__gp is None:
            self._instantiate_gp()
        elif self.__initialized:
            return self.__gp

        self.__initialized = True

        with open(self.__template_filename) as tf:
            template = string.Template(tf.read())

        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile(suffix='.gp', mode='w+t') as f:
            f.write(template.substitute(i=str(self.__instance)))
            f.flush()
            self.__gp.read(f.name)

        self._gp_eval('default(realprecision, %s)' % (self.prec // 3 + 2))
        self._gp_set_inst('conductor', self.conductor)
        self._gp_set_inst('gammaV', self.gammaV)
        self._gp_set_inst('weight', self.weight)
        self._gp_set_inst('sgn', self.eps)
        self._gp_set_inst('Lpoles', self.poles)
        self._gp_set_inst('Lresidues', self.residues)
        return self.__gp

    @classmethod
    def _instantiate_gp(cls):
        from sage.env import DOT_SAGE
        logfile = os.path.join(DOT_SAGE, 'dokchitser.log')
        cls.__gp = sage.interfaces.gp.Gp(script_subdirectory='dokchitser',
                                         logfile=logfile)
        # Read the script template and parse out all indexed global variables
        # (easy because they all end in "_$i" and there's nothing else in the
        # script that uses $)
        global_re = re.compile(r'([a-zA-Z0-9]+)_\$i')
        with open(cls.__template_filename) as f:
            for line in f:
                for m in global_re.finditer(line):
                    cls.__globals.add(m.group(1))

        cls.__globals_re = re.compile(
            '([^a-zA-Z0-9_]|^)(%s)([^a-zA-Z0-9_]|$)' % '|'.join(cls.__globals))
        return

    @classmethod
    def _teardown_gp(cls, instance=None):
        cls.__n_instances -= 1
        if cls.__n_instances == 0:
            cls.__gp.quit()
        elif instance is not None:
            # Clean up all global variables created by this instance
            for varname in cls.__globals:
                cls.__gp.eval('kill(%s_%s)' % (varname, instance))

    def _gp_call_inst(self, func, *args):
        """
        Call the specified PARI function in the GP interpreter, with the
        instance number of this `Dokchitser` instance automatically appended.

        For example, ``self._gp_call_inst('L', 1)`` is equivalent to
        ``self.gp().eval('L_N(1)')`` where ``N`` is ``self.__instance``.
        """

        cmd = '%s_%d(%s)' % (func, self.__instance,
                             ','.join(str(a) for a in args))
        return self._gp_eval(cmd)

    def _gp_set_inst(self, varname, value):
        """
        Like ``_gp_call_inst`` but sets the variable given by ``varname`` to
        the given value, appending ``_N`` to the variable name.

        If ``varname`` is a function (e.g. ``'func(n)'``) then this sets
        ``func_N(n)``).
        """

        if '(' in varname:
            funcname, args = varname.split('(', 1)
            varname = '%s_%s(%s' % (funcname, self.__instance, args)
        else:
            varname = '%s_%s' % (varname, self.__instance)

        cmd = '%s = %s' % (varname, value)
        return self._gp_eval(cmd)

    def _gp_eval(self, s):
        try:
            t = self.gp().eval(s)
        except (RuntimeError, TypeError):
            raise RuntimeError("unable to create L-series, due to precision or other limits in PARI")
        if not self.__init and '***' in t:
            # After init_coeffs is called, future calls to this method should
            # return the full output for further parsing
            raise RuntimeError("unable to create L-series, due to precision or other limits in PARI")
        t = t.replace("  *** _^_: Warning: normalizing a series with 0 leading term.\n", "")
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
            sage: L = E.lseries().dokchitser(algorithm='gp')
            sage: L.num_coeffs()
            26
            sage: E = EllipticCurve('5077a')
            sage: L = E.lseries().dokchitser(algorithm='gp')
            sage: L.num_coeffs()
            568
            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.num_coeffs()
            4

        Verify that ``num_coeffs`` works with non-real spectral
        parameters, e.g. for the `L`-function of the level 10 Maass form
        with eigenvalue 2.7341055592527126::

            sage: ev = 2.7341055592527126
            sage: L = Dokchitser(conductor=10, gammaV=[ev*i, -ev*i],weight=2,eps=1)
            sage: L.num_coeffs()
            26
        """
        return Integer(self._gp_call_inst('cflength', T))

    def init_coeffs(self, v, cutoff=1,
                    w=None,
                    pari_precode='',
                    max_imaginary_part=0,
                    max_asymp_coeffs=40):
        """
        Set the coefficients `a_n` of the `L`-series.

        If `L(s)` is not equal to its dual, pass the coefficients of
        the dual as the second optional argument.

        INPUT:

        - ``v`` -- list of complex numbers or string (pari function of k)

        - ``cutoff`` -- real number (default: 1)

        - ``w`` -- list of complex numbers or string (pari function of k)

        - ``pari_precode`` -- some code to execute in pari
          before calling initLdata

        - ``max_imaginary_part`` -- (default: 0) redefine if
          you want to compute L(s) for s having large imaginary part

        - ``max_asymp_coeffs`` -- (default: 40) at most this
          many terms are generated in asymptotic series for phi(t) and G(s,t)

        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: pari_precode = 'tau(n)=(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5))'
            sage: L.init_coeffs('tau(k)', pari_precode=pari_precode)

        Evaluate the resulting `L`-function at a point, and compare with
        the answer that one gets "by definition" (of `L`-function
        attached to a modular form)::

            sage: L(14)
            0.998583063162746
            sage: a = delta_qexp(1000)
            sage: sum(a[n]/float(n)^14 for n in reversed(range(1,1000)))
            0.9985830631627461

        Illustrate that one can give a list of complex numbers for v
        (see :issue:`10937`)::

            sage: L2 = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: L2.init_coeffs(list(delta_qexp(1000))[1:])
            sage: L2(14)
            0.998583063162746

        TESTS:

        Verify that setting the `w` parameter does not raise an error
        (see :issue:`10937`).  Note that the meaning of `w` does not seem to
        be documented anywhere in Dokchitser's package yet, so there is
        no claim that the example below is meaningful! ::

            sage: L2 = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: L2.init_coeffs(list(delta_qexp(1000))[1:], w=[1..1000])
        """
        if isinstance(v, tuple) and w is None:
            v, cutoff, w, pari_precode, max_imaginary_part, max_asymp_coeffs = v

        if pari_precode:
            # Must have __globals_re instantiated
            if self.__gp is None:
                self._instantiate_gp()

            def repl(m):
                return '%s%s_%d%s' % (m.group(1), m.group(2), self.__instance,
                                      m.group(3))

            # If any of the pre-code contains references to some of the
            # templated global variables we must replace those as well
            pari_precode = self.__globals_re.sub(repl, pari_precode)

        if pari_precode != '':
            self._gp_eval(pari_precode)
        RR = self.__CC._real_field()
        cutoff = RR(cutoff)
        if isinstance(v, str):
            if w is None:
                self._gp_call_inst('initLdata', '"%s"' % v, cutoff)
            else:
                self._gp_call_inst('initLdata', '"%s"' % v, cutoff, '"%s"' % w)
        elif not isinstance(v, (list, tuple)):
            raise TypeError("v (=%s) must be a list, tuple, or string" % v)
        else:
            CC = self.__CC
            v = ','.join(CC(a)._pari_init_() for a in v)
            self._gp_eval('Avec = [%s]' % v)
            if w is None:
                self._gp_call_inst('initLdata', '"Avec[k]"', cutoff)
            else:
                w = ','.join(CC(a)._pari_init_() for a in w)
                self._gp_eval('Bvec = [%s]' % w)
                self._gp_call_inst('initLdata', '"Avec[k]"', cutoff,
                                   '"Bvec[k]"')
        self.__init = (v, cutoff, w, pari_precode, max_imaginary_part,
                       max_asymp_coeffs)

    def _clear_value_cache(self):
        del self.__values

    def __call__(self, s, c=None):
        r"""
        INPUT:

        - ``s`` -- complex number

        - ``c`` -- internal parameter, call with `c>1` to get the same value
          with a different cutoff point (`c` close to `1`); should return the
          same answer, good to check if everything works with right precision

        .. NOTE::

           Evaluation of the function takes a long time, so each
           evaluation is cached. Call ``self._clear_value_cache()`` to
           clear the evaluation cache.

        EXAMPLES::

            sage: E = EllipticCurve('5077a')
            sage: L = E.lseries().dokchitser(100, algorithm='gp')
            sage: L(1)
            0.00000000000000000000000000000
            sage: L(1+I)
            -1.3085436607849493358323930438 + 0.81298000036784359634835412129*I
            sage: L(1+I, 1.2)
            -1.3085436607849493358323930438 + 0.81298000036784359634835412129*I

        TESTS::

            sage: L(1+I, 0)
            Traceback (most recent call last):
            ...
            RuntimeError
        """
        self.__check_init()
        s = self.__CC(s)
        try:
            return self.__values[s, c]
        except AttributeError:
            self.__values = {}
        except KeyError:
            pass
        if c is None:
            z = self._gp_call_inst('L', s)
        else:
            z = self._gp_call_inst('L', s, c)
        CC = self.__CC
        if 'pole' in z:
            print(z)
            raise ArithmeticError
        elif '***' in z:
            print(z)
            raise RuntimeError
        elif 'Warning' in z:
            i = z.rfind('\n')
            msg = z[:i].replace('digits', 'decimal digits')
            verbose(msg, level=-1)
            ans = CC(z[i + 1:])
            self.__values[s, c] = ans
            return ans
        ans = CC(z)
        self.__values[s, c] = ans
        return ans

    def derivative(self, s, k=1):
        r"""
        Return the `k`-th derivative of the `L`-series at `s`.

        .. WARNING::

           If `k` is greater than the order of vanishing of
           `L` at `s` you may get nonsense.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries().dokchitser(algorithm='gp')
            sage: L.derivative(1,E.rank())
            1.51863300057685
        """
        self.__check_init()
        s = self.__CC(s)
        k = Integer(k)
        z = self._gp_call_inst('L', s, '', k)
        if 'pole' in z:
            raise ArithmeticError(z)
        elif 'Warning' in z:
            i = z.rfind('\n')
            msg = z[:i].replace('digits', 'decimal digits')
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

        - ``a`` -- complex number (default: 0); point about which to expand

        - ``k`` -- integer (default: 6); series is `O(\texttt{var}^k)`

        - ``var`` -- string (default: ``'z'``); variable of power series

        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.taylor_series(2, 3)
            1.64493406684823 - 0.937548254315844*z + 0.994640117149451*z^2 + O(z^3)
            sage: E = EllipticCurve('37a')
            sage: L = E.lseries().dokchitser(algorithm='gp')
            sage: L.taylor_series(1)
            0.000000000000000 + 0.305999773834052*z + 0.186547797268162*z^2 - 0.136791463097188*z^3 + 0.0161066468496401*z^4 + 0.0185955175398802*z^5 + O(z^6)

        We compute a Taylor series where each coefficient is to high
        precision. ::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries().dokchitser(200, algorithm='gp')
            sage: L.taylor_series(1,3)
            ...e-82 + (...e-82)*z + 0.75931650028842677023019260789472201907809751649492435158581*z^2 + O(z^3)

        Check that :issue:`25402` is fixed::

            sage: L = EllipticCurve("24a1").modular_form().lseries()
            sage: L.taylor_series(-1, 3)
            0.000000000000000 - 0.702565506265199*z + 0.638929001045535*z^2 + O(z^3)

        Check that :issue:`25965` is fixed::

            sage: L2 = EllipticCurve("37a1").modular_form().lseries(); L2
            L-series associated to the cusp form q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: L2.taylor_series(0,4)
            0.000000000000000 - 0.357620466127498*z + 0.273373112603865*z^2 + 0.303362857047671*z^3 + O(z^4)
            sage: L2.taylor_series(0,1)
            O(z^1)
            sage: L2(0)
            0.000000000000000
        """
        self.__check_init()
        a = self.__CC(a)
        k = Integer(k)
        try:
            z = self._gp_call_inst('Lseries', a, '', k - 1)
            z = self.gp()('Vecrev(Pol(%s))' % z)
        except TypeError as msg:
            raise RuntimeError("%s\nUnable to compute Taylor expansion (try lowering the number of terms)" % msg)
        r = repr(z)
        if 'pole' in r:
            raise ArithmeticError(r)
        elif 'Warning' in r:
            i = r.rfind('\n')
            msg = r[:i].replace('digits', 'decimal digits')
            verbose(msg, level=-1)
        v = list(z)
        K = self.__CC
        v = [K(repr(x)) for x in v]
        R = self.__CC[[var]]
        return R(v, k)

    def check_functional_equation(self, T=1.2):
        r"""
        Verifies how well numerically the functional equation is satisfied,
        and also determines the residues if ``self.poles !=
        []`` and residues='automatic'.

        More specifically: for `T>1` (default: 1.2),
        ``self.check_functional_equation(T)`` should ideally
        return 0 (to the current precision).

        - if what this function returns does not look like 0 at all,
          probably the functional equation is wrong (i.e. some of the
          parameters gammaV, conductor etc., or the coefficients are wrong),

        - if checkfeq(T) is to be used, more coefficients have to be
          generated (approximately T times more), e.g. call cflength(1.3),
          initLdata("a(k)",1.3), checkfeq(1.3)

        - T=1 always (!) returns 0, so T has to be away from 1

        - default value `T=1.2` seems to give a reasonable
          balance

        - if you don't have to verify the functional equation or the
          L-values, call num_coeffs(1) and initLdata("a(k)",1), you need
          slightly less coefficients.

        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=1, poles=[1], residues=[-1], init='1')
            sage: L.check_functional_equation()  # abs tol 1e-19
            -2.71050543121376e-20

        If we choose the sign in functional equation for the
        `\zeta` function incorrectly, the functional equation
        doesn't check out. ::

            sage: L = Dokchitser(conductor=1, gammaV=[0], weight=1, eps=-11, poles=[1], residues=[-1], init='1')
            sage: L.check_functional_equation()
            -9.73967861488124
        """
        self.__check_init()
        z = self._gp_call_inst('checkfeq', T)
        z = z.replace(' ', '')
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

        - ``coefgrow`` -- string that evaluates to a PARI function of n that
          defines a coefgrow function

        EXAMPLES::

            sage: L = Dokchitser(conductor=1, gammaV=[0,1], weight=12, eps=1)
            sage: pari_precode = 'tau(n)=(5*sigma(n,3)+7*sigma(n,5))*n/12 - 35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5))'
            sage: L.init_coeffs('tau(k)', pari_precode=pari_precode)
            sage: L.set_coeff_growth('2*n^(11/2)')
            sage: L(1)
            0.0374412812685155
        """
        if not isinstance(coefgrow, str):
            raise TypeError("coefgrow must be a string")

        self._gp_set_inst('coefgrow(n)', coefgrow.replace('\n', ' '))


def reduce_load_dokchitser(D):
    X = Dokchitser(1, 1, 1, 1)
    X.__dict__ = D
    X.init_coeffs(X._Dokchitser__init)
    return X
