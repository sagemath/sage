r"""
Mathematical constants

The following standard mathematical constants are defined in Sage,
along with support for coercing them into GAP, PARI/GP, KASH,
Maxima, Mathematica, Maple, Octave, and Singular::

    sage: pi
    pi
    sage: e             # base of the natural logarithm
    e
    sage: NaN           # Not a number
    NaN
    sage: golden_ratio
    golden_ratio
    sage: log2          # natural logarithm of the real number 2
    log2
    sage: euler_gamma   # Euler's gamma constant
    euler_gamma
    sage: catalan       # the Catalan constant
    catalan
    sage: khinchin      # Khinchin's constant
    khinchin
    sage: twinprime
    twinprime
    sage: mertens
    mertens

Support for coercion into the various systems means that if, e.g.,
you want to create `\pi` in Maxima and Singular, you don't
have to figure out the special notation for each system. You just
type the following::

    sage: maxima(pi)
    %pi
    sage: singular(pi)
    pi
    sage: gap(pi)
    pi
    sage: gp(pi)
    3.141592653589793238462643383     # 32-bit
    3.1415926535897932384626433832795028842   # 64-bit
    sage: pari(pi)
    3.14159265358979
    sage: kash(pi)                    # optional - kash
    3.14159265358979323846264338328
    sage: mathematica(pi)             # optional - mathematica
    Pi
    sage: pi._maple_init_()
    'Pi'
    sage: octave(pi)                  # optional - octave
    3.14159

Arithmetic operations with constants also yield constants, which
can be coerced into other systems or evaluated.

::

    sage: a = pi + e*4/5; a
    pi + 4/5*e
    sage: maxima(a)
    %pi+(4*%e)/5
    sage: RealField(15)(a)           # 15 *bits* of precision
    5.316
    sage: gp(a)
    5.316218116357029426750873360              # 32-bit
    5.3162181163570294267508733603616328824    # 64-bit
    sage: print(mathematica(a))                  # optional - mathematica
     4 E
     --- + Pi
      5

EXAMPLES: Decimal expansions of constants

We can obtain floating point approximations to each of these
constants by coercing into the real field with given precision. For
example, to 200 binary places we have the following::

    sage: R = RealField(200); R
    Real Field with 200 bits of precision

::

    sage: R(pi)
    3.1415926535897932384626433832795028841971693993751058209749

::

    sage: R(e)
    2.7182818284590452353602874713526624977572470936999595749670

::

    sage: R(NaN)
    NaN

::

    sage: R(golden_ratio)
    1.6180339887498948482045868343656381177203091798057628621354

::

    sage: R(log2)
    0.69314718055994530941723212145817656807550013436025525412068

::

    sage: R(euler_gamma)
    0.57721566490153286060651209008240243104215933593992359880577

::

    sage: R(catalan)
    0.91596559417721901505460351493238411077414937428167213426650

::

    sage: R(khinchin)
    2.6854520010653064453097148354817956938203822939944629530512

EXAMPLES: Arithmetic with constants

::

    sage: f = I*(e+1); f
    I*e + I
    sage: f^2
    (I*e + I)^2
    sage: _.expand()
    -e^2 - 2*e - 1

::

    sage: pp = pi+pi; pp
    2*pi
    sage: R(pp)
    6.2831853071795864769252867665590057683943387987502116419499

::

    sage: s = (1 + e^pi); s
    e^pi + 1
    sage: R(s)
    24.140692632779269005729086367948547380266106242600211993445
    sage: R(s-1)
    23.140692632779269005729086367948547380266106242600211993445

::

    sage: l = (1-log2)/(1+log2); l
    -(log2 - 1)/(log2 + 1)
    sage: R(l)
    0.18123221829928249948761381864650311423330609774776013488056

::

    sage: pim = maxima(pi)
    sage: maxima.eval('fpprec : 100')
    '100'
    sage: pim.bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068b0

AUTHORS:

- Alex Clemesha (2006-01-15)

- William Stein

- Alex Clemesha, William Stein (2006-02-20): added new constants;
  removed todos

- Didier Deshommes (2007-03-27): added constants from RQDF (deprecated)

TESTS:

Coercing the sum of a bunch of the constants to many different
floating point rings::

    sage: a = pi + e + golden_ratio + log2 + euler_gamma + catalan + khinchin + twinprime + mertens; a
    mertens + twinprime + khinchin + log2 + golden_ratio + catalan + euler_gamma + pi + e
    sage: parent(a)
    Symbolic Ring
    sage: RR(a)  # abs tol 1e-13
    13.2713479401972
    sage: RealField(212)(a)
    13.2713479401972493100988191995758139408711068200030748178329712
    sage: RealField(230)(a)
    13.271347940197249310098819199575813940871106820003074817832971189555
    sage: RDF(a)  # abs tol 1e-13
    13.271347940197249
    sage: CC(a)  # abs tol 1e-13
    13.2713479401972
    sage: CDF(a)  # abs tol 1e-13
    13.271347940197249
    sage: ComplexField(230)(a)
    13.271347940197249310098819199575813940871106820003074817832971189555

Check that :issue:`8237` is fixed::

    sage: maxima('infinity').sage()
    Infinity
    sage: maxima('inf').sage()
    +Infinity
    sage: maxima('minf').sage()
    -Infinity
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008-2010 Burcin Erocal <burcin@erocal.org>
#                     2009 Mike Hansen <mhansen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################

import math
from functools import partial
from sage.rings.infinity import (infinity, minus_infinity,
                                 unsigned_infinity)
from sage.structure.richcmp import richcmp_method, op_EQ, op_GE, op_LE
from sage.symbolic.expression import register_symbol, init_pynac_I
from sage.symbolic.expression import E

constants_table = {}
constants_name_table = {}
constants_name_table[repr(infinity)] = infinity
constants_name_table[repr(unsigned_infinity)] = unsigned_infinity
constants_name_table[repr(minus_infinity)] = minus_infinity

I = init_pynac_I()

register_symbol(infinity, {'maxima': 'inf'}, 0)
register_symbol(minus_infinity, {'maxima': 'minf'}, 0)
register_symbol(unsigned_infinity, {'maxima': 'infinity'}, 0)
register_symbol(I, {'mathematica': 'I'}, 0)
register_symbol(True, {'giac': 'true',
                       'mathematica': 'True',
                       'maxima': 'true'}, 0)
register_symbol(False, {'giac': 'false',
                        'mathematica': 'False',
                        'maxima': 'false'}, 0)


def unpickle_Constant(class_name, name, conversions, latex, mathml, domain):
    """
    EXAMPLES::

        sage: from sage.symbolic.constants import unpickle_Constant
        sage: a = unpickle_Constant('Constant', 'a', {}, 'aa', '', 'positive')
        sage: a.domain()
        'positive'
        sage: latex(a)
        aa

    Note that if the name already appears in the
    ``constants_name_table``, then that will be returned instead of
    constructing a new object::

        sage: pi = unpickle_Constant('Pi', 'pi', None, None, None, None)
        sage: pi._maxima_init_()
        '%pi'
    """
    if name in constants_name_table:
        return constants_name_table[name]
    if class_name == "Constant":
        return Constant(name, conversions=conversions, latex=latex,
                        mathml=mathml, domain=domain)
    else:
        cls = globals()[class_name]
        return cls(name=name)


@richcmp_method
class Constant():
    def __init__(self, name, conversions=None, latex=None, mathml='',
                 domain='complex'):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: loads(dumps(p))
            p
        """
        self._conversions = conversions if conversions is not None else {}
        self._latex = latex if latex is not None else name
        self._mathml = mathml
        self._name = name
        self._domain = domain

        for system, value in self._conversions.items():
            setattr(self, "_%s_" % system, partial(self._generic_interface, value))
            setattr(self, "_%s_init_" % system, partial(self._generic_interface_init, value))

        from .expression import PynacConstant
        self._pynac = PynacConstant(self._name, self._latex, self._domain)
        self._serial = self._pynac.serial()
        constants_table[self._serial] = self
        constants_name_table[self._name] = self

        register_symbol(self.expression(), self._conversions)

    def __richcmp__(self, other, op):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: s = Constant('s')
            sage: p == p
            True
            sage: p == s
            False
            sage: p != s
            True
        """
        if self.__class__ == other.__class__ and self._name == other._name:
            return op in [op_EQ, op_GE, op_LE]
        else:
            return NotImplemented

    def __reduce__(self):
        """
        Add support for pickling constants.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: p.__reduce__()
            (<function unpickle_Constant at 0x...>,
             ('Constant', 'p', {}, 'p', '', 'complex'))
            sage: loads(dumps(p))
            p

            sage: pi.pyobject().__reduce__()
            (<function unpickle_Constant at 0x...>,
             ('Pi',
              'pi',
              ...,
              '\\pi',
              '<mi>&pi;</mi>',
              'positive'))
            sage: loads(dumps(pi.pyobject()))
            pi
        """
        return (unpickle_Constant, (self.__class__.__name__, self._name,
                                    self._conversions, self._latex,
                                    self._mathml, self._domain))

    def domain(self):
        """
        Return the domain of this constant.  This is either positive,
        real, or complex, and is used by Pynac to make inferences
        about expressions containing this constant.

        EXAMPLES::

            sage: p = pi.pyobject(); p
            pi
            sage: type(_)
            <class 'sage.symbolic.constants.Pi'>
            sage: p.domain()
            'positive'
        """
        return self._domain

    def expression(self):
        """
        Return an expression for this constant.

        EXAMPLES::

            sage: a = pi.pyobject()
            sage: pi2 = a.expression()
            sage: pi2
            pi
            sage: pi2 + 2
            pi + 2
            sage: pi - pi2
            0
        """
        return self._pynac.expression()

    def _symbolic_(self, SR):
        """
        Return an expression for this constant.

        INPUT:

        - ``SR`` -- a symbolic ring parent

        EXAMPLES::

            sage: SR(pi.pyobject())
            pi
            sage: pi.pyobject()._symbolic_(SR)
            pi
            sage: f(x,y) = 2
            sage: f.parent()(pi.pyobject())
            (x, y) |--> pi
        """
        return SR(self.expression())

    def name(self):
        """
        Return the name of this constant.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c')
            sage: c.name()
            'c'
        """
        return self._name

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c')
            sage: c
            c
        """
        return self._name

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c', latex=r'\xi')
            sage: latex(c)
            \xi
        """
        return self._latex

    def _mathml_(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c', mathml=r'<mi>c</mi>')
            sage: mathml(c)
            <mi>c</mi>
        """
        return self._mathml

    def _generic_interface(self, value, I):
        """
        This is a helper method used in defining the ``_X_`` methods
        where ``X`` is the name of some interface.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p', conversions=dict(maxima='%pi'))
            sage: p._maxima_(maxima)
            %pi

        The above ``_maxima_`` is constructed like ``m`` below::

            sage: from functools import partial
            sage: m = partial(p._generic_interface, '%pi')
            sage: m(maxima)
            %pi
        """
        return I(value)

    def _generic_interface_init(self, value):
        """
        This is a helper method used in defining the ``_X_init_`` methods
        where ``X`` is the name of some interface.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p', conversions=dict(maxima='%pi'))
            sage: p._maxima_init_()
            '%pi'

        The above ``_maxima_init_`` is constructed like ``mi`` below::

            sage: from functools import partial
            sage: mi = partial(p._generic_interface_init, '%pi')
            sage: mi()
            '%pi'
        """
        return value

    def _interface_(self, I):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p', conversions=dict(maxima='%pi'))
            sage: p._interface_(maxima)
            %pi
        """
        try:
            s = self._conversions[I.name()]
            return I(s)
        except KeyError:
            pass

        try:
            return getattr(self, "_%s_" % (I.name()))(I)
        except AttributeError:
            pass

        raise NotImplementedError

    def _gap_(self, gap):
        """
        Return the constant as a string in GAP. Since GAP does not
        have floating point numbers, we simply return the constant as
        a string.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: gap(p)
            p
        """
        return gap('"%s"' % self)

    def _singular_(self, singular):
        """
        Return the constant as a string in Singular. Since Singular
        does not always support floating point numbers, we simply
        return the constant as a string.  (Singular allows floating point
        numbers if the current ring has floating point coefficients,
        but not otherwise.)

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: singular(p)
            p
        """
        return singular('"%s"' % self)


class Pi(Constant):
    def __init__(self, name='pi'):
        r"""
        TESTS::

            sage: pi._latex_()
            '\\pi'
            sage: latex(pi)
            \pi
            sage: mathml(pi)
            <mi>&pi;</mi>
        """
        conversions = dict(axiom='%pi', fricas='%pi', maxima='%pi', giac='pi',
                           gp='Pi', kash='PI',
                           mathematica='Pi', matlab='pi', maple='Pi',
                           octave='pi', pari='Pi', pynac='Pi')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r"\pi", mathml="<mi>&pi;</mi>",
                          domain='positive')

    def __float__(self):
        """
        EXAMPLES::

            sage: float(pi)
            3.141592653589793
        """
        return math.pi

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: pi._mpfr_(RealField(100))
            3.1415926535897932384626433833
        """
        return R.pi()

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: pi._real_double_(RDF)
            3.141592653589793
        """
        return R.pi()

    def _sympy_(self):
        """
        Convert pi to sympy pi.

        EXAMPLES::

            sage: import sympy                                                          # needs sympy
            sage: sympy.pi == pi  # indirect doctest                                    # needs sympy
            True
        """
        import sympy
        return sympy.pi


pi = Pi().expression()
"""
The formal square root of -1.

EXAMPLES::

    sage: SR.I()
    I
    sage: SR.I()^2
    -1

Note that conversions to real fields will give TypeErrors::

    sage: float(SR.I())
    Traceback (most recent call last):
    ...
    TypeError: unable to simplify to float approximation
    sage: gp(SR.I())
    I
    sage: RR(SR.I())
    Traceback (most recent call last):
    ...
    TypeError: unable to convert '1.00000000000000*I' to a real number

Expressions involving I that are real-valued can be converted to real fields::

    sage: float(I*I)
    -1.0
    sage: RR(I*I)
    -1.00000000000000

We can convert to complex fields::

    sage: C = ComplexField(200); C
    Complex Field with 200 bits of precision
    sage: C(SR.I())
    1.0000000000000000000000000000000000000000000000000000000000*I
    sage: SR.I()._complex_mpfr_field_(ComplexField(53))
    1.00000000000000*I

    sage: SR.I()._complex_double_(CDF)
    1.0*I
    sage: CDF(SR.I())
    1.0*I

    sage: z = SR.I() + I; z
    2*I
    sage: C(z)
    2.0000000000000000000000000000000000000000000000000000000000*I
    sage: 1e8*SR.I()
    1.00000000000000e8*I

    sage: complex(SR.I())
    1j

    sage: QQbar(SR.I())
    I

    sage: abs(SR.I())
    1

    sage: SR.I().minpoly()
    x^2 + 1
    sage: maxima(2*SR.I())
    2*%i

TESTS::

    sage: repr(SR.I())
    'I'
    sage: latex(SR.I())
    i
"""

# The base of the natural logarithm, e, is not a constant in GiNaC/Sage. It is
# represented by exp(1). A dummy class to make this work with arithmetic and
# coercion is implemented in the module sage.symbolic.expression for speed.
e = E()

# Allow for backtranslation to this symbol from Mathematica (#29833).
register_symbol(e, {'mathematica': 'E'})


class NotANumber(Constant):
    """
    Not a Number
    """
    def __init__(self, name="NaN"):
        """
        EXAMPLES::

            sage: loads(dumps(NaN))
            NaN
        """
        conversions = dict(matlab='NaN')
        Constant.__init__(self, name, conversions=conversions)

    def __float__(self):
        """
        EXAMPLES::

            sage: float(NaN)
            nan
        """
        return float('nan')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: NaN._mpfr_(RealField(53))
            NaN
            sage: type(_)
            <class 'sage.rings.real_mpfr.RealNumber'>
        """
        return R('NaN')  # ??? nan in mpfr: void mpfr_set_nan (mpfr_t x)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(NaN)
            NaN
        """
        return R.NaN()

    def _sympy_(self):
        """
        Convert ``NaN`` to SymPy NaN.

        EXAMPLES::

            sage: bool(NaN._sympy_()._sage_() == NaN)                                   # needs sympy
            True
            sage: import sympy                                                          # needs sympy
            sage: sympy.nan == NaN  # this should be fixed                              # needs sympy
            False
        """
        import sympy
        return sympy.nan


NaN = NotANumber().expression()


class GoldenRatio(Constant):
    """
    The number (1+sqrt(5))/2.

    EXAMPLES::

        sage: gr = golden_ratio
        sage: RR(gr)
        1.61803398874989
        sage: R = RealField(200)
        sage: R(gr)
        1.6180339887498948482045868343656381177203091798057628621354
        sage: grm = maxima(golden_ratio);grm
        (sqrt(5)+1)/2
        sage: grm + grm
        sqrt(5)+1
        sage: float(grm + grm)
        3.23606797749979
    """
    def __init__(self, name='golden_ratio'):
        """
        EXAMPLES::

            sage: loads(dumps(golden_ratio))
            golden_ratio
        """
        conversions = dict(mathematica='(1+Sqrt[5])/2', gp='(1+sqrt(5))/2',
                           maple='(1+sqrt(5))/2', maxima='(1+sqrt(5))/2',
                           pari='(1+sqrt(5))/2', octave='(1+sqrt(5))/2',
                           kash='(1+Sqrt(5))/2', giac='(1+sqrt(5))/2')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r'\phi', domain='positive')

    def minpoly(self, bits=None, degree=None, epsilon=0):
        """
        EXAMPLES::

            sage: golden_ratio.minpoly()
            x^2 - x - 1
        """
        from sage.rings.rational_field import QQ
        x = QQ['x'].gen(0)
        return x**2 - x - 1

    def __float__(self):
        """
        EXAMPLES::

            sage: float(golden_ratio)
            1.618033988749895
            sage: golden_ratio.__float__()
            1.618033988749895
        """
        return 0.5 + math.sqrt(1.25)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(golden_ratio)
            1.618033988749895
        """
        return R('1.61803398874989484820458')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: golden_ratio._mpfr_(RealField(100))
            1.6180339887498948482045868344
            sage: RealField(100)(golden_ratio)
            1.6180339887498948482045868344
        """
        return (R(1) + R(5).sqrt()) / R(2)

    def _algebraic_(self, field):
        """
        EXAMPLES::

            sage: golden_ratio._algebraic_(QQbar)
            1.618033988749895?
            sage: QQbar(golden_ratio)
            1.618033988749895?
        """
        import sage.rings.qqbar
        return field(sage.rings.qqbar.get_AA_golden_ratio())

    def _sympy_(self):
        """
        Convert ``golden_ratio`` to SymPy GoldenRatio.

        EXAMPLES::

            sage: import sympy                                                          # needs sympy
            sage: sympy.GoldenRatio == golden_ratio  # indirect doctest                 # needs sympy
            True
        """
        import sympy
        return sympy.GoldenRatio


golden_ratio = GoldenRatio().expression()


class Log2(Constant):
    """
    The natural logarithm of the real number 2.

    EXAMPLES::

        sage: log2
        log2
        sage: float(log2)
        0.6931471805599453
        sage: RR(log2)
        0.693147180559945
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(log2)
        0.69314718055994530941723212145817656807550013436025525412068
        sage: l = (1-log2)/(1+log2); l
        -(log2 - 1)/(log2 + 1)
        sage: R(l)
        0.18123221829928249948761381864650311423330609774776013488056
        sage: maxima(log2)
        log(2)
        sage: maxima(log2).float()
        0.6931471805599453
        sage: gp(log2)
        0.6931471805599453094172321215             # 32-bit
        0.69314718055994530941723212145817656807   # 64-bit
        sage: RealField(150)(2).log()
        0.69314718055994530941723212145817656807550013
    """
    def __init__(self, name='log2'):
        """
        EXAMPLES::

            sage: loads(dumps(log2))
            log2
        """
        conversions = dict(mathematica='Log[2]', kash='Log(2)',
                           maple='log(2)', maxima='log(2)', gp='log(2)',
                           pari='log(2)', octave='log(2)')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r'\log(2)', domain='positive')

    def __float__(self):
        """
        EXAMPLES::

            sage: float(log2)
            0.6931471805599453
            sage: log2.__float__()
            0.6931471805599453
        """
        return math.log(2)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(log2)
            0.6931471805599453
        """
        return R.log2()

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: RealField(100)(log2)
            0.69314718055994530941723212146
            sage: log2._mpfr_(RealField(100))
            0.69314718055994530941723212146
        """
        return R.log2()


log2 = Log2().expression()


class EulerGamma(Constant):
    """
    The limiting difference between the harmonic series and the natural
    logarithm.

    EXAMPLES::

        sage: R = RealField()
        sage: R(euler_gamma)
        0.577215664901533
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(euler_gamma)
        0.57721566490153286060651209008240243104215933593992359880577
        sage: eg = euler_gamma + euler_gamma; eg
        2*euler_gamma
        sage: R(eg)
        1.1544313298030657212130241801648048620843186718798471976115
    """
    def __init__(self, name='euler_gamma'):
        """
        EXAMPLES::

            sage: loads(dumps(euler_gamma))
            euler_gamma
        """
        conversions = dict(kash='EulerGamma(R)', maple='gamma',
                           mathematica='EulerGamma', pari='Euler',
                           maxima='%gamma', pynac='Euler', giac='euler_gamma',
                           fricas='-digamma(1)')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r'\gamma', domain='positive')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: RealField(100)(euler_gamma)
            0.57721566490153286060651209008
            sage: euler_gamma._mpfr_(RealField(100))
            0.57721566490153286060651209008
        """
        return R.euler_constant()

    def __float__(self):
        """
        EXAMPLES::

            sage: float(euler_gamma)
            0.5772156649015329
        """
        return 0.57721566490153286060651209008

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(euler_gamma)
            0.5772156649015329
        """
        return R.euler_constant()

    def _sympy_(self):
        """
        Convert ``euler_gamma`` to SymPy EulerGamma.

        EXAMPLES::

            sage: import sympy                                                          # needs sympy
            sage: sympy.EulerGamma == euler_gamma  # indirect doctest                   # needs sympy
            True
        """
        import sympy
        return sympy.EulerGamma


euler_gamma = EulerGamma().expression()


class Catalan(Constant):
    """
    A number appearing in combinatorics defined as the Dirichlet beta
    function evaluated at the number 2.

    EXAMPLES::

        sage: catalan^2 + mertens
        mertens + catalan^2
    """
    def __init__(self, name='catalan'):
        """
        EXAMPLES::

            sage: loads(dumps(catalan))
            catalan
        """
        # kash: R is default prec
        conversions = dict(mathematica='Catalan', kash='Catalan(R)',
                           maple='Catalan', maxima='catalan',
                           pynac='Catalan')
        Constant.__init__(self, name, conversions=conversions,
                          domain='positive')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: RealField(100)(catalan)
            0.91596559417721901505460351493
            sage: catalan._mpfr_(RealField(100))
            0.91596559417721901505460351493
        """
        return R.catalan_constant()

    def _real_double_(self, R):
        """
        EXAMPLES: We coerce to the real double field::

            sage: RDF(catalan)
            0.915965594177219
        """
        return R('0.91596559417721901505460351493252')

    def __float__(self):
        """
        EXAMPLES::

            sage: float(catalan)
            0.915965594177219
        """
        return 0.91596559417721901505460351493252

    def _sympy_(self):
        """
        Convert ``catalan`` to SymPy Catalan.

        EXAMPLES::

            sage: import sympy                                                          # needs sympy
            sage: sympy.Catalan == catalan  # indirect doctest                          # needs sympy
            True
        """
        import sympy
        return sympy.Catalan


catalan = Catalan().expression()


class Khinchin(Constant):
    """
    The geometric mean of the continued fraction expansion of any
    (almost any) real number.

    EXAMPLES::

        sage: float(khinchin)
        2.6854520010653062
        sage: khinchin.n(digits=60)
        2.68545200106530644530971483548179569382038229399446295305115
        sage: m = mathematica(khinchin); m             # optional - mathematica
        Khinchin
        sage: m.N(200)                                 # optional - mathematica
        2.685452001065306445309714835481795693820382293...32852204481940961807
    """
    def __init__(self, name='khinchin'):
        """
        EXAMPLES::

            sage: loads(dumps(khinchin))
            khinchin
        """
        conversions = dict(maxima='khinchin', mathematica='Khinchin',
            pynac='Khinchin')
        Constant.__init__(self, name, conversions=conversions,
                          domain='positive')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: khinchin._mpfr_(RealField(100))
            2.6854520010653064453097148355
            sage: RealField(100)(khinchin)
            2.6854520010653064453097148355
        """
        import sage.libs.mpmath.all as a
        return a.eval_constant('khinchin', R)

    def __float__(self):
        """
        EXAMPLES::

            sage: float(khinchin)
            2.6854520010653062
        """
        return 2.6854520010653064453097148355


khinchin = Khinchin().expression()


class TwinPrime(Constant):
    r"""
    The Twin Primes constant is defined as
    `\prod 1 - 1/(p-1)^2` for primes `p > 2`.

    EXAMPLES::

        sage: float(twinprime)
        0.6601618158468696
        sage: twinprime.n(digits=60)
        0.660161815846869573927812110014555778432623360284733413319448
    """
    def __init__(self, name='twinprime'):
        """
        EXAMPLES::

            sage: loads(dumps(twinprime))
            twinprime
        """
        conversions = dict(maxima='twinprime', pynac='TwinPrime')
        Constant.__init__(self, name, conversions=conversions,
                          domain='positive')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: twinprime._mpfr_(RealField(100))
            0.66016181584686957392781211001
            sage: RealField(100)(twinprime)
            0.66016181584686957392781211001
        """
        import sage.libs.mpmath.all as a
        return a.eval_constant('twinprime', R)

    def __float__(self):
        """
        EXAMPLES::

            sage: float(twinprime)
            0.6601618158468696
        """
        return 0.66016181584686957392781211001


twinprime = TwinPrime().expression()


class Mertens(Constant):
    """
    The Mertens constant is related to the Twin Primes constant and
    appears in Mertens' second theorem.

    EXAMPLES::

        sage: float(mertens)
        0.26149721284764277
        sage: mertens.n(digits=60)
        0.261497212847642783755426838608695859051566648261199206192064
    """
    def __init__(self, name='mertens'):
        """
        EXAMPLES::

            sage: loads(dumps(mertens))
            mertens
        """
        conversions = dict(maxima='mertens', pynac='Mertens')
        Constant.__init__(self, name, conversions=conversions,
                          domain='positive')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: mertens._mpfr_(RealField(100))
            0.26149721284764278375542683861
            sage: RealField(100)(mertens)
            0.26149721284764278375542683861
        """
        import sage.libs.mpmath.all as a
        return a.eval_constant('mertens', R)

    def __float__(self):
        """
        EXAMPLES::

            sage: float(mertens)
            0.26149721284764277
        """
        return 0.26149721284764278375542683861


mertens = Mertens().expression()


class Glaisher(Constant):
    r"""
    The Glaisher-Kinkelin constant `A = \exp(\frac{1}{12}-\zeta'(-1))`.

    EXAMPLES::

        sage: float(glaisher)
        1.2824271291006226
        sage: glaisher.n(digits=60)
        1.28242712910062263687534256886979172776768892732500119206374
        sage: a = glaisher + 2
        sage: a
        glaisher + 2
        sage: parent(a)
        Symbolic Ring
    """
    def __init__(self, name='glaisher'):
        """
        EXAMPLES::

            sage: loads(dumps(glaisher))
            glaisher
        """
        conversions = dict(maxima='glaisher', pynac='Glaisher',
            mathematica='Glaisher')
        Constant.__init__(self, name, conversions=conversions,
                          domain='positive')

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: glaisher._mpfr_(RealField(100))
            1.2824271291006226368753425689
            sage: RealField(100)(glaisher)
            1.2824271291006226368753425689
        """
        import sage.libs.mpmath.all as a
        return a.eval_constant('glaisher', R)

    def __float__(self):
        """
        EXAMPLES::

            sage: float(glaisher)
            1.2824271291006226
        """
        return 1.2824271291006226368753425689


glaisher = Glaisher().expression()
