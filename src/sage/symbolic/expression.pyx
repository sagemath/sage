"""
Symbolic Expressions

RELATIONAL EXPRESSIONS:

We create a relational expression::

    sage: x = var('x')
    sage: eqn = (x-1)^2 <= x^2 - 2*x + 3
    sage: eqn.subs(x == 5)
    16 <= 18

Notice that squaring the relation squares both sides.

::

    sage: eqn^2
    (x - 1)^4 <= (x^2 - 2*x + 3)^2
    sage: eqn.expand()
    x^2 - 2*x + 1 <= x^2 - 2*x + 3

The can transform a true relational into a false one::

    sage: eqn = SR(-5) < SR(-3); eqn
    -5 < -3
    sage: bool(eqn)
    True
    sage: eqn^2
    25 < 9
    sage: bool(eqn^2)
    False

We can do arithmetic with relationals::

    sage: e = x+1 <= x-2
    sage: e + 2
    x + 3 <= x
    sage: e - 1
    x <= x - 3
    sage: e*(-1)
    -x - 1 <= -x + 2
    sage: (-2)*e
    -2*x - 2 <= -2*x + 4
    sage: e*5
    5*x + 5 <= 5*x - 10
    sage: e/5
    1/5*x + 1/5 <= 1/5*x - 2/5
    sage: 5/e
    5/(x + 1) <= 5/(x - 2)
    sage: e/(-2)
    -1/2*x - 1/2 <= -1/2*x + 1
    sage: -2/e
    -2/(x + 1) <= -2/(x - 2)

We can even add together two relations, so long as the operators are
the same::

    sage: (x^3 + x <= x - 17)  + (-x <= x - 10)
    x^3 <= 2*x - 27

Here they aren't::

    sage: (x^3 + x <= x - 17)  + (-x >= x - 10)
    Traceback (most recent call last):
    ...
    TypeError: incompatible relations


ARBITRARY SAGE ELEMENTS:

You can work symbolically with any Sage data type.  This can lead to
nonsense if the data type is strange, e.g., an element of a finite
field (at present).

We mix Singular variables with symbolic variables::

    sage: R.<u,v> = QQ[]
    sage: var('a,b,c')
    (a, b, c)
    sage: expand((u + v + a + b + c)^2)
    a^2 + 2*a*b + b^2 + 2*a*c + 2*b*c + c^2 + 2*a*u + 2*b*u + 2*c*u + u^2 + 2*a*v + 2*b*v + 2*c*v + 2*u*v + v^2

TESTS:

Test Jacobian on Pynac expressions. #5546 ::

    sage: var('x,y')
    (x, y)
    sage: f = x + y
    sage: jacobian(f, [x,y])
    [1 1]


Test if matrices work #5546 ::

    sage: var('x,y,z')
    (x, y, z)
    sage: M = matrix(2,2,[x,y,z,x])
    sage: v = vector([x,y])
    sage: M * v
    (x^2 + y^2, x*y + x*z)
    sage: v*M
    (x^2 + y*z, 2*x*y)

Test if comparison bugs from #6256 are fixed::

    sage: t = exp(sqrt(x)); u = 1/t
    sage: t*u
    1
    sage: t + u
    e^(-sqrt(x)) + e^sqrt(x)
    sage: t
    e^sqrt(x)

Test if #9947 is fixed::

    sage: real_part(1+2*(sqrt(2)+1)*(sqrt(2)-1))
    3
    sage: a=(sqrt(4*(sqrt(3) - 5)*(sqrt(3) + 5) + 48) + 4*sqrt(3))/ (sqrt(3) + 5)
    sage: a.real_part()
    4*sqrt(3)/(sqrt(3) + 5)
    sage: a.imag_part()
    sqrt(abs(4*(sqrt(3) + 5)*(sqrt(3) - 5) + 48))/(sqrt(3) + 5)
"""

###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"
include "sage/ext/python.pxi"

import operator
import ring
import sage.rings.integer
import sage.rings.rational
from sage.structure.element cimport ModuleElement, RingElement, Element
from sage.symbolic.getitem cimport OperandsWrapper
from sage.symbolic.function import get_sfunction_from_serial, SymbolicFunction
from sage.rings.rational import Rational  # Used for sqrt.
from sage.misc.derivative import multi_derivative
from sage.rings.infinity import AnInfinity, infinity, minus_infinity, unsigned_infinity
from sage.misc.decorators import rename_keyword
from sage.misc.superseded import deprecated_function_alias
from sage.structure.dynamic_class import dynamic_class

LOG_TEN_TWO_PLUS_EPSILON = 3.321928094887363 # a small overestimate of log(10,2)

cpdef bint is_Expression(x):
    """
    Returns True if *x* is a symbolic Expression.

    EXAMPLES::

        sage: from sage.symbolic.expression import is_Expression
        sage: is_Expression(x)
        True
        sage: is_Expression(2)
        False
        sage: is_Expression(SR(2))
        True
    """
    return isinstance(x, Expression)

cpdef bint is_SymbolicEquation(x):
    """
    Returns True if *x* is a symbolic equation.

    EXAMPLES:

    The following two examples are symbolic equations::

        sage: from sage.symbolic.expression import is_SymbolicEquation
        sage: is_SymbolicEquation(sin(x) == x)
        True
        sage: is_SymbolicEquation(sin(x) < x)
        True
        sage: is_SymbolicEquation(x)
        False

    This is not, since ``2==3`` evaluates to the boolean
    ``False``::

        sage: is_SymbolicEquation(2 == 3)
        False

    However here since both 2 and 3 are coerced to be symbolic, we
    obtain a symbolic equation::

        sage: is_SymbolicEquation(SR(2) == SR(3))
        True

    """
    return isinstance(x, Expression) and is_a_relational((<Expression>x)._gobj)

cdef class Expression(CommutativeRingElement):
    cpdef object pyobject(self):
        """
        Get the underlying Python object.

        OUTPUT:

        The Python object corresponding to this expression, assuming
        this expression is a single numerical value or an infinity
        representable in Python. Otherwise, a ``TypeError`` is raised.

        EXAMPLES::

            sage: var('x')
            x
            sage: b = -17/3
            sage: a = SR(b)
            sage: a.pyobject()
            -17/3
            sage: a.pyobject() is b
            True

        TESTS::

            sage: SR(oo).pyobject()
            +Infinity
            sage: SR(-oo).pyobject()
            -Infinity
            sage: SR(unsigned_infinity).pyobject()
            Infinity
            sage: SR(I*oo).pyobject()
            Traceback (most recent call last):
            ...
            TypeError: Python infinity cannot have complex phase.
        """
        cdef GConstant* c
        if is_a_constant(self._gobj):
            from sage.symbolic.constants import constants_name_table
            return constants_name_table[GEx_to_str(&self._gobj)]

        if is_a_infinity(self._gobj):
            if (ex_to_infinity(self._gobj).is_unsigned_infinity()): return unsigned_infinity
            if (ex_to_infinity(self._gobj).is_plus_infinity()):     return infinity
            if (ex_to_infinity(self._gobj).is_minus_infinity()):    return minus_infinity
            raise TypeError('Python infinity cannot have complex phase.')

        if not is_a_numeric(self._gobj):
            raise TypeError, "self must be a numeric expression"
        return py_object_from_numeric(self._gobj)

    def __init__(self, SR, x=0):
        """
        Nearly all expressions are created by calling new_Expression_from_*,
        but we need to make sure this at least doesn't leave self._gobj
        uninitialized and segfault.

        TESTS::

            sage: sage.symbolic.expression.Expression(SR)
            0
            sage: sage.symbolic.expression.Expression(SR, 5)
            5

        We test subclassing ``Expression``::

            sage: from sage.symbolic.expression import Expression
            sage: class exp_sub(Expression): pass
            sage: f = function('f')
            sage: t = f(x)
            sage: u = exp_sub(SR, t)
            sage: u.operator()
            f
        """
        self._parent = SR
        cdef Expression exp = self.coerce_in(x)
        GEx_construct_ex(&self._gobj, exp._gobj)

    def __dealloc__(self):
        """
        Delete memory occupied by this expression.
        """
        GEx_destruct(&self._gobj)

    def __getstate__(self):
        """
        Returns a tuple describing the state of this expression for pickling.

        This should return all information that will be required to unpickle
        the object. The functionality for unpickling is implemented in
        __setstate__().

        In order to pickle Expression objects, we return a tuple containing

         * 0  - as pickle version number
                in case we decide to change the pickle format in the feature
         * names of symbols of this expression
         * a string representation of self stored in a Pynac archive.

        TESTS::
            sage: var('x,y,z')
            (x, y, z)
            sage: t = 2*x*y^z+3
            sage: s = dumps(t)

            sage: t.__getstate__()
            (0,
             ['x', 'y', 'z'],
             ...)

        """
        cdef GArchive ar
        ar.archive_ex(self._gobj, "sage_ex")
        ar_str = GArchive_to_str(&ar)
        return (0, map(repr, self.variables()), ar_str)

    def _dbgprint(self):
        r"""
        Print pynac debug output to ``stderr``.

        EXAMPLES::

            sage: (1+x)._dbgprint()
            x + 1
        """
        self._gobj.dbgprint()

    def _dbgprinttree(self):
        r"""
        Print pynac expression tree debug output to ``stderr``.

        EXAMPLES:

        The expression tree is composed of Ginac primitives
        and functions, organized by the tree, with various
        other memory and hash information which will vary::

            sage: (1+x+exp(x+1))._dbgprinttree()    # not tested
            add @0x65e5960, hash=0x4727e01a, flags=0x3, nops=3
                x (symbol) @0x6209150, serial=6, hash=0x2057b15e, flags=0xf, domain=0
                1 (numeric) @0x3474cf0, hash=0x0, flags=0x7
                -----
                function exp @0x24835d0, hash=0x765c2165, flags=0xb, nops=1
                    add @0x65df570, hash=0x420740d2, flags=0xb, nops=2
                        x (symbol) @0x6209150, serial=6, hash=0x2057b15e, flags=0xf, domain=0
                        1 (numeric) @0x3474cf0, hash=0x0, flags=0x7
                        -----
                        overall_coeff
                        1 (numeric) @0x65e4df0, hash=0x7fd3, flags=0x7
                        =====
                    =====
                1 (numeric) @0x3474cf0, hash=0x0, flags=0x7
                -----
                overall_coeff
                1 (numeric) @0x663cc40, hash=0x7fd3, flags=0x7
                =====

        TESTS:

        This test is just to make sure the function is working::

            sage: (1+x+exp(x+1))._dbgprinttree()
            add @...
                x (symbol) ...
                1 (numeric) ...
                ...
                overall_coeff
                1 (numeric) ...
                =====
        """
        self._gobj.dbgprinttree();

    def __setstate__(self, state):
        """
        Initializes the state of the object from data saved in a pickle.

        During unpickling __init__ methods of classes are not called, the saved
        data is passed to the class via this function instead.

        TESTS::
            sage: var('x,y,z')
            (x, y, z)
            sage: t = 2*x*y^z+3
            sage: u = loads(dumps(t)) # indirect doctest
            sage: u
            2*x*y^z + 3
            sage: bool(t == u)
            True
            sage: u.subs(x=z)
            2*y^z*z + 3

            sage: loads(dumps(x.parent()(2)))
            2
        """
        # check input
        if state[0] != 0 or len(state) != 3:
            raise ValueError, "unknown state information"
        # set parent
        self._set_parent(ring.SR)
        # get variables
        cdef GExList sym_lst
        for name in state[1]:
            sym_lst.append_sym(\
                    ex_to_symbol((<Expression>ring.SR.symbol(name))._gobj))

        # initialize archive
        cdef GArchive ar
        GArchive_from_str(&ar, state[2], len(state[2]))

        # extract the expression from the archive
        GEx_construct_ex(&self._gobj, ar.unarchive_ex(sym_lst, <unsigned>0))

    def __copy__(self):
        """
        TESTS::

            sage: copy(x)
            x
        """
        return new_Expression_from_GEx(self._parent, self._gobj)

    def _repr_(self):
        """
        Return string representation of this symbolic expression.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: repr(x+y)
            'x + y'

        TESTS::

            # printing of modular number equal to -1 as coefficient
            sage: k.<a> = GF(9); k(2)*x
            2*x

            sage: (x+1)*Mod(6,7)
            6*x + 6

            #printing rational functions
            sage: x/y
            x/y
            sage: x/2/y
            1/2*x/y
            sage: .5*x/y
            0.500000000000000*x/y
            sage: x^(-1)
            1/x
            sage: x^(-5)
            x^(-5)
            sage: x^(-y)
            x^(-y)
            sage: 2*x^(-1)
            2/x
            sage: i*x
            I*x
            sage: -x.parent(i)
            -I
            sage: y + 3*(x^(-1))
            y + 3/x

        Printing the exp function::

            sage: x.parent(1).exp()
            e
            sage: x.exp()
            e^x

        Powers::

            sage: _ = var('A,B,n'); (A*B)^n
            (A*B)^n
            sage: (A/B)^n
            (A/B)^n
            sage: n*x^(n-1)
            n*x^(n - 1)
            sage: (A*B)^(n+1)
            (A*B)^(n + 1)
            sage: (A/B)^(n-1)
            (A/B)^(n - 1)
            sage: n*x^(n+1)
            n*x^(n + 1)
            sage: n*x^(n-1)
            n*x^(n - 1)
            sage: n*(A/B)^(n+1)
            n*(A/B)^(n + 1)
            sage: (n+A/B)^(n+1)
            (n + A/B)^(n + 1)

        Powers where the base or exponent is a Python object::

            sage: (2/3)^x
            (2/3)^x
            sage: x^CDF(1,2)
            x^(1.0 + 2.0*I)
            sage: (2/3)^(2/3)
            (2/3)^(2/3)
            sage: (-x)^(1/4)
            (-x)^(1/4)
            sage: k.<a> = GF(9)
            sage: SR(a+1)^x
            (a + 1)^x

        Check if #7876 is fixed::

            sage: (1/2-1/2*I )*sqrt(2)
            -(1/2*I - 1/2)*sqrt(2)
            sage: latex((1/2-1/2*I )*sqrt(2))
            -\left(\frac{1}{2} i - \frac{1}{2}\right) \, \sqrt{2}

        Check if :trac:`9632` is fixed::

            sage: zeta(x) + cos(x)
            cos(x) + zeta(x)
            sage: psi(1,1/3)*log(3)
            log(3)*psi(1, 1/3)
        """
        return self._parent._repr_element_(self)

    def _ascii_art_(self):
        """
        TESTS::

            sage: i = var('i')
            sage: ascii_art(sum(i^2/pi*x^i, i, 0, oo))
                          2
                         x  + x
            -------------------------------
                  3         2
            - pi*x  + 3*pi*x  - 3*pi*x + pi
            sage: ascii_art(integral(exp(x + x^2)/(x+1), x))
              /
             |
             |   2
             |  x  + x
             | e
             | ------- dx
             |  x + 1
             |
            /
        """
        from sympy import pretty, sympify
        from sage.misc.ascii_art import AsciiArt
        # FIXME:: when *sage* will use at least sympy >= 0.7.2
        # we could use a nice splitting with respect of the AsciiArt module.
        # from sage.misc.ascii_art import AsciiArt, MAX_LENGTH ## for import
        #            num_columns = MAX_LENGTH  ## option of pretty
        try:
            s = pretty(sympify(self), use_unicode=False)
        except StandardError:
            s = self
        return AsciiArt(str(s).splitlines())

    def _interface_(self, I):
        """
        EXAMPLES::

            sage: f = sin(e + 2)
            sage: f._interface_(sage.calculus.calculus.maxima)
            sin(%e+2)
        """
        if is_a_constant(self._gobj):
            return self.pyobject()._interface_(I)
        return super(Expression, self)._interface_(I)

    def _maxima_(self, session=None):
        """
        EXAMPLES::

            sage: f = sin(e + 2)
            sage: f._maxima_()
            sin(%e+2)
            sage: _.parent() is sage.calculus.calculus.maxima
            True
        """
        if session is None:
            # This chooses the Maxima interface used by calculus
            # Maybe not such a great idea because the "default" interface is another one
            from sage.calculus.calculus import maxima
            return super(Expression, self)._interface_(maxima)
        else:
            return super(Expression, self)._interface_(session)

    def _interface_init_(self, I):
        """
        EXAMPLES::

            sage: a = (pi + 2).sin()
            sage: a._maxima_init_()
            'sin((%pi)+(2))'

            sage: a = (pi + 2).sin()
            sage: a._maple_init_()
            'sin((pi)+(2))'

            sage: a = (pi + 2).sin()
            sage: a._mathematica_init_()
            'Sin[(Pi)+(2)]'

            sage: f = pi + I*e
            sage: f._pari_init_()
            '(Pi)+((exp(1))*(I))'

        TESTS:

        Check if complex numbers are converted to Maxima correctly
        :trac:`7557`::

            sage: SR(1.5*I)._maxima_init_()
            '1.5000000000000000*%i'
            sage: SR(CC.0)._maxima_init_()
            '1.0000000000000000*%i'
            sage: SR(CDF.0)._maxima_init_()
            '1.0000000000000000*%i'
        """
        from sage.symbolic.expression_conversions import InterfaceInit
        return InterfaceInit(I)(self)

    def _gap_init_(self):
        """
        Conversion of symbolic object to GAP always results in a GAP
        string.

        EXAMPLES::

            sage: gap(e + pi^2 + x^3)
            x^3 + pi^2 + e
        """
        return '"%s"'%repr(self)

    def _singular_init_(self):
        """
        Conversion of a symbolic object to Singular always results in a
        Singular string.

        EXAMPLES::

            sage: singular(e + pi^2 + x^3)
            x^3 + pi^2 + e
        """
        return '"%s"'%repr(self)

    def _magma_init_(self, magma):
        """
        Return string representation in Magma of this symbolic expression.

        Since Magma has no notation of symbolic calculus, this simply
        returns something that evaluates in Magma to a a Magma string.

        EXAMPLES::

            sage: x = var('x')
            sage: f = sin(cos(x^2) + log(x))
            sage: f._magma_init_(magma)
            '"sin(cos(x^2) + log(x))"'
            sage: magma(f)                         # optional - magma
            sin(log(x) + cos(x^2))
            sage: magma(f).Type()                  # optional - magma
            MonStgElt
        """
        return '"%s"'%repr(self)

    def _latex_(self):
        r"""
        Return string representation of this symbolic expression.

        TESTS::

            sage: var('x,y,z')
            (x, y, z)
            sage: latex(y + 3*(x^(-1)))
            y + \frac{3}{x}
            sage: latex(x^(y+z^(1/y)))
            x^{y + z^{\left(\frac{1}{y}\right)}}
            sage: latex(1/sqrt(x+y))
            \frac{1}{\sqrt{x + y}}
            sage: latex(sin(x*(z+y)^x))
            \sin\left(x {\left(y + z\right)}^{x}\right)
            sage: latex(3/2*(x+y)/z/y)
            \frac{3 \, {\left(x + y\right)}}{2 \, y z}
            sage: latex((2^(x^y)))
            2^{\left(x^{y}\right)}
            sage: latex(abs(x))
            {\left| x \right|}
            sage: latex((x*y).conjugate())
            \overline{x} \overline{y}
            sage: latex(x*(1/(x^2)+sqrt(x^7)))
            x {\left(\sqrt{x^{7}} + \frac{1}{x^{2}}\right)}

        Check spacing of coefficients of mul expressions (#3202)::

            sage: latex(2*3^x)
            2 \, 3^{x}

        Powers::

            sage: _ = var('A,B,n')
            sage: latex((n+A/B)^(n+1))
            {\left(n + \frac{A}{B}\right)}^{n + 1}
            sage: latex((A*B)^n)
            \left(A B\right)^{n}
            sage: latex((A*B)^(n-1))
            \left(A B\right)^{n - 1}

        Powers where the base or exponent is a Python object::

            sage: latex((2/3)^x)
            \left(\frac{2}{3}\right)^{x}
            sage: latex(x^CDF(1,2))
            x^{1.0 + 2.0i}
            sage: latex((2/3)^(2/3))
            \left(\frac{2}{3}\right)^{\frac{2}{3}}
            sage: latex((-x)^(1/4))
            \left(-x\right)^{\frac{1}{4}}
            sage: k.<a> = GF(9)
            sage: latex(SR(a+1)^x)
            \left(a + 1\right)^{x}

        More powers, #7406::

            sage: latex((x^pi)^e)
            {\left(x^{\pi}\right)}^{e}
            sage: latex((x^(pi+1))^e)
            {\left(x^{\pi + 1}\right)}^{e}
            sage: a,b,c = var('a b c')
            sage: latex(a^(b^c))
            a^{\left(b^{c}\right)}
            sage: latex((a^b)^c)
            {\left(a^{b}\right)}^{c}

        Separate coefficients to numerator and denominator, #7363::

            sage: latex(2/(x+1))
            \frac{2}{x + 1}
            sage: latex(1/2/(x+1))
            \frac{1}{2 \, {\left(x + 1\right)}}

        Check if rational function coefficients without a ``numerator()`` method
        are printed correctly. #8491::

            sage: latex(6.5/x)
            \frac{6.50000000000000}{x}
            sage: latex(Mod(2,7)/x)
            \frac{2}{x}

        Check if we avoid extra parenthesis in rational functions (#8688)::

            sage: latex((x+2)/(x^3+1))
            \frac{x + 2}{x^{3} + 1}
            sage: latex((x+2)*(x+1)/(x^3+1))
            \frac{{\left(x + 2\right)} {\left(x + 1\right)}}{x^{3} + 1}
            sage: latex((x+2)/(x^3+1)/(x+1))
            \frac{x + 2}{{\left(x^{3} + 1\right)} {\left(x + 1\right)}}

        Check that the sign is correct (#9086)::

            sage: latex(-1/x)
            -\frac{1}{x}
            sage: latex(1/-x)
            -\frac{1}{x}

        More tests for the sign (#9314)::

            sage: latex(-2/x)
            -\frac{2}{x}
            sage: latex(-x/y)
            -\frac{x}{y}
            sage: latex(-x*z/y)
            -\frac{x z}{y}
            sage: latex(-x/z/y)
            -\frac{x}{y z}

        Check if #9394 is fixed::

            sage: var('n')
            n
            sage: latex( e^(2*I*pi*n*x - 2*I*pi*n) )
            e^{\left(2 i \, \pi n x - 2 i \, \pi n\right)}
            sage: latex( e^(2*I*pi*n*x - (2*I+1)*pi*n) )
            e^{\left(2 i \, \pi n x - \left(2 i + 1\right) \, \pi n\right)}
            sage: x+(1-2*I)*y
            x - (2*I - 1)*y
            sage: latex(x+(1-2*I)*y)
            x - \left(2 i - 1\right) \, y

        Check if complex coefficients with denominators are displayed
        correctly #10769::

            sage: var('a x')
            (a, x)
            sage: latex(1/2*I/x)
            \frac{i}{2 \, x}
            sage: ratio = i/2* x^2/a
            sage: latex(ratio)
            \frac{i \, x^{2}}{2 \, a}

        Parenthesis in powers, #13262::

            sage: latex(1+x^(2/3)+x^(-2/3))
            x^{\frac{2}{3}} + \frac{1}{x^{\frac{2}{3}}} + 1
        """
        return self._parent._latex_element_(self)

    def _mathml_(self):
        """
        Returns a MathML representation of this object.

        EXAMPLES::

            sage: mathml(pi)
            <mi>&pi;</mi>
            sage: mathml(pi+2)
            MATHML version of the string pi + 2

        """
        from sage.misc.all import mathml
        try:
            obj = self.pyobject()
        except TypeError:
            return mathml(repr(self))
        return mathml(obj)

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: f = x^3 + 17*x -3
            sage: ZZ(f.coeff(x^3))
            1
            sage: ZZ(f.coeff(x))
            17
            sage: ZZ(f.coeff(x,0))
            -3
            sage: type(ZZ(f.coeff(x,0)))
            <type 'sage.rings.integer.Integer'>

        Coercion is done if necessary::

            sage: f = x^3 + 17/1*x
            sage: ZZ(f.coeff(x))
            17
            sage: type(ZZ(f.coeff(x)))
            <type 'sage.rings.integer.Integer'>

        If the symbolic expression is just a wrapper around an integer,
        that very same integer is returned::

            sage: n = 17; SR(n)._integer_() is n
            True
        """
        try:
            n = self.pyobject()
        except TypeError:
            raise TypeError, "unable to convert x (=%s) to an integer"%(self)
        if isinstance(n, sage.rings.integer.Integer):
            return n
        return sage.rings.integer.Integer(n)

    def __int__(self):
        """
        EXAMPLES::

            sage: int(log(8)/log(2))
            3
            sage: int(-log(8)/log(2))
            -3
            sage: int(sin(2)*100)
            90
            sage: int(-sin(2)*100)
            -90
            sage: int(SR(3^64)) == 3^64
            True
            sage: int(SR(10^100)) == 10^100
            True
            sage: int(SR(10^100-10^-100)) == 10^100 - 1
            True
            sage: int(sqrt(-3))
            Traceback (most recent call last):
            ...
            ValueError: cannot convert sqrt(-3) to int
        """
        from sage.functions.all import floor, ceil
        try:
            rif_self = sage.rings.all.RIF(self)
        except TypeError:
            raise ValueError, "cannot convert %s to int"%(self)
        if rif_self > 0 or (rif_self.contains_zero() and self > 0):
            result = floor(self)
        else:
            result = ceil(self)
        if not isinstance(result, sage.rings.integer.Integer):
            raise ValueError, "cannot convert %s to int"%(self)
        else:
            return int(result)

    def __long__(self):
        """
        EXAMPLES::

            sage: long(sin(2)*100)
            90L
        """
        return long(int(self))

    def _rational_(self):
        """
        EXAMPLES::

            sage: f = x^3 + 17/1*x - 3/8
            sage: QQ(f.coeff(x^2))
            0
            sage: QQ(f.coeff(x^3))
            1
            sage: a = QQ(f.coeff(x)); a
            17
            sage: type(a)
            <type 'sage.rings.rational.Rational'>
            sage: QQ(f.coeff(x,0))
            -3/8

        If the symbolic expression is just a wrapper around a rational,
        that very same rational is returned::

            sage: n = 17/1; SR(n)._rational_() is n
            True
        """
        try:
            n = self.pyobject()
        except TypeError:
            raise TypeError, "unable to convert %s to a rational"%self
        if isinstance(n, sage.rings.rational.Rational):
            return n
        return sage.rings.rational.Rational(n)

    cpdef _eval_self(self, R):
        """
        Evaluate this expression numerically.

        This function is used to convert symbolic expressions to ``RR``,
        ``CC``, ``float``, ``complex``, ``CIF`` and ``RIF``.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: sin(x).subs(x=5)._eval_self(RR)
            -0.958924274663138
            sage: gamma(x).subs(x=I)._eval_self(CC)
            -0.154949828301811 - 0.498015668118356*I
            sage: x._eval_self(CC)
            Traceback (most recent call last):
            ...
            TypeError: Cannot evaluate symbolic expression to a numeric value.

        Check if we can compute a real evaluation even if the expression
        contains complex coefficients::

            sage: RR((I - sqrt(2))*(I+sqrt(2)))
            -3.00000000000000
            sage: cos(I)._eval_self(RR)
            1.54308063481524
            sage: float(cos(I))
            1.5430806348152437
        """
        cdef GEx res
        try:
            res = self._gobj.evalf(0, R)
        except TypeError as err:
            # try the evaluation again with the complex field
            # corresponding to the parent R
            if R is float:
                R_complex = complex
            else:
                try:
                    R_complex = R.complex_field()
                except (TypeError, AttributeError):
                    raise err
            res = self._gobj.evalf(0, R_complex)
        if is_a_numeric(res):
            return R(py_object_from_numeric(res))
        else:
            raise TypeError, "Cannot evaluate symbolic expression to a numeric value."

    cpdef _convert(self, R):
        """
        Convert all the numeric coefficients and constants in this expression
        to the given ring `R`. This results in an expression which contains
        only variables, and functions whose arguments contain a variable.

        EXAMPLES::

            sage: f = sqrt(2) * cos(3); f
            sqrt(2)*cos(3)
            sage: f._convert(RDF)
            -1.40006081534
            sage: f._convert(float)
            -1.40006081533995

        There is nothing to convert for variables::

            sage: x._convert(CC)
            x

        Note that the output is not meant to be in the in the given ring `R`.
        Since the results of some functions will still be  floating point
        approximations::

            sage: t = log(10); t
            log(10)
            sage: t._convert(QQ)
            2.30258509299405

        ::

            sage: (0.25 / (log(5.74 /x^0.9, 10))^2 / 4)._convert(QQ)
            0.331368631904900/log(287/50/x^0.900000000000000)^2
            sage: (0.25 / (log(5.74 /x^0.9, 10))^2 / 4)._convert(CC)
            0.331368631904900/log(5.74000000000000/x^0.900000000000000)^2

        When converting to an exact domain, powers remain unevaluated::

            sage: f = sqrt(2) * cos(3); f
            sqrt(2)*cos(3)
            sage: f._convert(int)
            -0.989992496600445*sqrt(2)
        """
        cdef GEx res = self._gobj.evalf(0, R)
        return new_Expression_from_GEx(self._parent, res)

    def _mpfr_(self, R):
        """
        Return a numerical approximation of this symbolic expression in the RealField R.

        The precision of the approximation is determined by the precision of
        the input R.

        EXAMPLES::

            0.090909090909090909090909090909090909090909090909090909090909

            sage: a = sin(3); a
            sin(3)
            sage: RealField(200)(a)
            0.14112000805986722210074480280811027984693326425226558415188
            sage: a._mpfr_(RealField(100))
            0.14112000805986722210074480281
        """
        return self._eval_self(R)

    def _real_mpfi_(self, R):
        """
        Returns this expression as a real interval.

        EXAMPLES::

            sage: RIF(sqrt(2))
            1.414213562373095?
        """
        try:
            return self._eval_self(R)
        except TypeError:
            raise TypeError, "unable to simplify to a real interval approximation"

    def _complex_mpfi_(self, R):
        """
        Returns this expression as a complex interval.

        EXAMPLES::

            sage: CIF(pi)
            3.141592653589794?
        """
        try:
            return self._eval_self(R)
        except TypeError:
            raise TypeError, "unable to simplify to a complex interval approximation"

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(sin(3))
            0.14112000806
        """
        return self._eval_self(R)

    def _complex_mpfr_field_(self, R):
        """
        Return a numerical approximation to this expression in the given
        ComplexField R.

        The precision of the approximation is determined by the precision of
        the input R.

        EXAMPLES::

            sage: ComplexField(200)(SR(1/11))
            0.090909090909090909090909090909090909090909090909090909090909
            sage: zeta(x).subs(x=I)._complex_mpfr_field_(ComplexField(70))
            0.0033002236853241028742 - 0.41815544914132167669*I
            sage: gamma(x).subs(x=I)._complex_mpfr_field_(ComplexField(60))
            -0.1549498283018106... - 0.49801566811835604*I
            sage: log(x).subs(x=I)._complex_mpfr_field_(ComplexField(50))
            1.5707963267949*I

            sage: CC(sqrt(2))
            1.41421356237309
            sage: a = sqrt(-2); a
            sqrt(-2)
            sage: CC(a).imag()
            1.41421356237309
            sage: ComplexField(200)(a).imag()
            1.4142135623730950488016887242096980785696718753769480731767
            sage: ComplexField(100)((-1)^(1/10))
            0.95105651629515357211643933338 + 0.30901699437494742410229341718*I
            sage: CC(x*sin(0))
            0.000000000000000
        """
        return self._eval_self(R)

    def _complex_double_(self, R):
        """
        Return a numerical approximation to this expression in the given
        Complex Double Field R.

        EXAMPLES::

            sage: CDF(SR(1/11))
            0.0909090909091
            sage: zeta(x).subs(x=I)._complex_double_(CDF)
            0.00330022368532 - 0.418155449141*I
            sage: gamma(x).subs(x=I)._complex_double_(CDF)
            -0.154949828302 - 0.498015668118*I
            sage: log(x).subs(x=I)._complex_double_(CDF)
            1.57079632679*I
            sage: CDF((-1)^(1/3))
            0.5 + 0.866025403784*I
        """
        return self._eval_self(R)

    def __float__(self):
        """
        Return float conversion of self, assuming self is constant.
        Otherwise, raise a TypeError.

        OUTPUT:

        A ``float``. Double precision evaluation of self.

        EXAMPLES::

            sage: float(SR(12))
            12.0
            sage: float(SR(2/3))
            0.6666666666666666
            sage: float(sqrt(SR(2)))
            1.4142135623730951
            sage: float(x^2 + 1)
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to float approximation
            sage: float(SR(RIF(2)))
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to float approximation
        """
        try:
            return float(self._eval_self(float))
        except TypeError:
            raise TypeError, "unable to simplify to float approximation"

    def __complex__(self):
        """
        EXAMPLES::

            sage: complex(I)
            1j
            sage: complex(erf(3*I))
            1629.9946226015657j
        """
        try:
            return self._eval_self(complex)
        except TypeError:
            raise TypeError, "unable to simplify to complex approximation"

    def _sympy_(self):
        """
        Returns a Sympy version of this object.

        EXAMPLES::

            sage: pi._sympy_()
            pi
            sage: type(_)
            <class 'sympy.core.numbers.Pi'>

        """
        from sage.symbolic.expression_conversions import sympy
        return sympy(self)

    def _algebraic_(self, field):
        """
        Convert a symbolic expression to an algebraic number.

        EXAMPLES::

            sage: QQbar(sqrt(2) + sqrt(8))
            4.242640687119285?
            sage: AA(sqrt(2) ^ 4) == 4
            True
            sage: AA(-golden_ratio)
            -1.618033988749895?
            sage: QQbar((2*I)^(1/2))
            1 + 1*I
            sage: QQbar(e^(pi*I/3))
            0.500000000000000? + 0.866025403784439?*I

            sage: QQbar(sqrt(2))
            1.414213562373095?
            sage: AA(abs(1+I))
            1.414213562373095?
            sage: golden_ratio._algebraic_(QQbar)
            1.618033988749895?
            sage: QQbar(golden_ratio)
            1.618033988749895?

            sage: AA(x*sin(0))
            0
            sage: QQbar(x*sin(0))
            0
        """
        from sage.symbolic.expression_conversions import algebraic
        return algebraic(self, field)

    def __hash__(self):
        """
        Return hash of this expression.

        EXAMPLES::

        The hash of an object in Python or its coerced version into
        the symbolic ring is the same::

            sage: hash(SR(3/1))
            3
            sage: hash(SR(19/23))
            4
            sage: hash(19/23)
            4

        The hash for symbolic expressions are unfortunately random. Here we
        only test that the hash() function returns without error, and that
        the return type is correct::

            sage: x, y = var("x y")
            sage: t = hash(x); type(t)
            <type 'int'>
            sage: t = hash(x^y); type(t)
            <type 'int'>
            sage: type(hash(x+y))
            <type 'int'>
            sage: d = {x+y: 5}
            sage: d
            {x + y: 5}

        In this example hashing is important otherwise the answer is
        wrong::

            sage: uniq([x-x, -x+x])
            [0]

        Test if exceptions during hashing are handled properly::

            sage: t = SR(matrix(2,2,range(4)))
            sage: hash(t)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        TESTS:

        Test if hashes for fderivatives with different parameters collide.
        #6243::

            sage: f = function('f'); t = f(x,y)
            sage: u = t.derivative(x); v = t.derivative(y)
            sage: hash(u) == hash(v)
            False
            sage: d = {u: 3, v: 5}; sorted(d.values())
            [3, 5]

        More checks for fderivative hashes #6851 ::

            sage: hash(f(x).derivative(x)) == hash(f(x).derivative(x,2))
            False
            sage: d = dict( (f(x).derivative(x, i), i) for i in range(1,6) )
            sage: len(d.keys())
            5

        We create a function with 10 arguments and test if there are hash
        collisions between any of its derivatives of order at most 7. #7508::

            sage: num_vars = 10; max_order=7
            sage: X = var(' '.join(['x'+str(i) for i in range(num_vars)]))
            sage: f = function('f',*X)
            sage: hashes=set()
            sage: for length in range(1,max_order+1):  # long time (4s on sage.math, 2012)
            ...       for s in UnorderedTuples(X, length):
            ...           deriv = f.diff(*s)
            ...           h = hash(deriv)
            ...           if h in hashes:
            ...               print "deriv: %s, hash:%s"%(deriv,h)
            ...           else:
            ...               hashes.add(n)
        """
        return self._gobj.gethash()

    # Boilerplate code from sage/structure/element.pyx
    def __richcmp__(left, right, int op):
        """
        Create a formal symbolic inequality or equality.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: x + 2/3 < y^2
            x + 2/3 < y^2
            sage: x^3 -y <= y + x
            x^3 - y <= x + y
            sage: x^3 -y == y + x
            x^3 - y == x + y
            sage: x^3 - y^10 >= y + x^10
            -y^10 + x^3 >= x^10 + y
            sage: x^2 > x
            x^2 > x

        Testing trac #11309 which changes the behavior of comparison of
        comparisons::

            sage: (-x + y < 0) in [x - y < 0]
            False
            sage: (x - 1 < 0) in [x - 2 < 0]
            False
            sage: Set([-x + y < 0, x - y < 0])
            {-x + y < 0, x - y < 0}
            sage: (x < y) == (x > y)
            False
            sage: (x < 0) < (x < 1)
            False
            sage: (x < y) != (y > x)
            False
            sage: (x >= y) == (y <= x)
            True
            sage: (x > y) == (y <= x)
            False
            sage: (x < x) == (x < x)
            True
            sage: (y > y) != (y > y)
            False
            sage: (x < y) != x
            True
            sage: (x == y) == (y == x)
            True
            sage: (x != y) != (y != x)
            False
            sage: (x == y) != (x != y)
            True
            sage: (x == y) == (y != x)
            False
            sage: x == (x == x)
            False
        """
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        cdef Expression l, r

        l = left
        r = right

        # If lhs or rhs is a relation, resolve the big relation
        # immediately UNLESS the lhs and rhs are flipped versions of
        # the same relation.
        if is_a_relational(l._gobj):
            if (op != Py_EQ and op != Py_NE):
                # relations aren't <, >, <=, or >= to other things
                return False
            if is_a_relational(r._gobj):
                # both lhs and rhs are relations, so we can get to work
                if l.operator() == r.operator():
                    e2 = ( # case: (x _ y) ?= (x _ y)
                           ( l._gobj.lhs().is_equal(r._gobj.lhs()) and
                             l._gobj.rhs().is_equal(r._gobj.rhs()) ) or

                           # case: (x == y) ?= (y == x)
                           #       (x != y) ?= (y != x)
                           ( ( l.operator() == operator.eq or
                               l.operator() == operator.ne ) and
                             l._gobj.lhs().is_equal(r._gobj.rhs()) and
                             l._gobj.rhs().is_equal(r._gobj.lhs()) ))
                else:
                    e2 = ( # case: (x < y)  ?= (y > x)  (or vice versa)
                           #       (x <= y) ?= (y >= x) (or vice versa)
                           ( ( l.operator() == operator.lt and
                               r.operator() == operator.gt ) or
                             ( l.operator() == operator.gt and
                               r.operator() == operator.lt ) or
                             ( l.operator() == operator.le and
                               r.operator() == operator.ge ) or
                             ( l.operator() == operator.ge and
                               r.operator() == operator.le ) ) and
                           l._gobj.lhs().is_equal(r._gobj.rhs()) and
                           l._gobj.rhs().is_equal(r._gobj.lhs()) )

            else:
                e2 = False              # l is relational but r isn't.

            if op == Py_EQ:
                return e2
            else:                       # op == Py_NE, checked earlier.
                return not e2

        elif is_a_relational(r._gobj):  # l isn't relational but r is.
            # things aren't <, >, <=, >=, or == to relations; they
            # are, however, != to relations
            return op == Py_NE

        # neither was relational, so we can create a symbolic relation
        cdef GEx e
        if op == Py_LT:
            e = g_lt(l._gobj, r._gobj)
        elif op == Py_EQ:
            e = g_eq(l._gobj, r._gobj)
        elif op == Py_GT:
            e = g_gt(l._gobj, r._gobj)
        elif op == Py_LE:
            e = g_le(l._gobj, r._gobj)
        elif op == Py_NE:
            e = g_ne(l._gobj, r._gobj)
        elif op == Py_GE:
            e = g_ge(l._gobj, r._gobj)
        else:
            raise TypeError
        return new_Expression_from_GEx(l._parent, e)

    def assume(self):
        r"""
        Assume that this equation holds. This is relevant for symbolic
        integration, among other things.

        EXAMPLES: We call the assume method to assume that `x>2`::

            sage: (x > 2).assume()

        Bool returns True below if the inequality is *definitely* known to
        be True.

        ::

            sage: bool(x > 0)
            True
            sage: bool(x < 0)
            False

        This may or may not be True, so bool returns False::

            sage: bool(x > 3)
            False

        If you make inconsistent or meaningless assumptions,
        Sage will let you know::

            sage: forget()
            sage: assume(x<0)
            sage: assume(x>0)
            Traceback (most recent call last):
            ...
            ValueError: Assumption is inconsistent
            sage: assumptions()
            [x < 0]
            sage: forget()

        TESTS::

            sage: v,c = var('v,c')
            sage: assume(c != 0)
            sage: integral((1+v^2/c^2)^3/(1-v^2/c^2)^(3/2),v)
            83/8*v/sqrt(-v^2/c^2 + 1) - 17/8*v^3/(c^2*sqrt(-v^2/c^2 + 1)) - 1/4*v^5/(c^4*sqrt(-v^2/c^2 + 1)) - 75/8*arcsin(v/(c^2*sqrt(c^(-2))))/sqrt(c^(-2))
            sage: forget()
        """
        from sage.symbolic.assumptions import _assumptions
        from sage.calculus.calculus import maxima
        if not self.is_relational():
            raise TypeError, "self (=%s) must be a relational expression"%self
        if not self in _assumptions:
            m = self._maxima_init_assume_()
            s = maxima.assume(m)
            if str(s._sage_()[0]) in ['meaningless','inconsistent','redundant']:
                raise ValueError, "Assumption is %s"%str(s._sage_()[0])
            _assumptions.append(self)

    def forget(self):
        """
        Forget the given constraint.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: forget()
            sage: assume(x>0, y < 2)
            sage: assumptions()
            [x > 0, y < 2]
            sage: forget(y < 2)
            sage: assumptions()
            [x > 0]

        TESTS:

        Check if #7507 is fixed::

            sage: forget()
            sage: n = var('n')
            sage: foo=sin((-1)*n*pi)
            sage: foo.simplify()
            -sin(pi*n)
            sage: assume(n, 'odd')
            sage: assumptions()
            [n is odd]
            sage: foo.simplify()
            0
            sage: forget(n, 'odd')
            sage: assumptions()
            []
            sage: foo.simplify()
            -sin(pi*n)
        """
        from sage.symbolic.assumptions import _assumptions
        from sage.calculus.calculus import maxima
        if not self.is_relational():
            raise TypeError, "self (=%s) must be a relational expression"%self
        m = self._maxima_init_assume_()
        maxima.forget(m)
        try:
            _assumptions.remove(self)
        except ValueError:
            pass

    def _maxima_init_assume_(self):
        """
        Return string that when evaluated in Maxima defines the assumption
        determined by this expression.

        EXAMPLES::

            sage: f = x+2 > sqrt(3)
            sage: f._maxima_init_assume_()
            '((x)+(2))>((3/1)^(1/2))'
        """
        from sage.calculus.calculus import maxima

        l = self.lhs()._assume_str()
        r = self.rhs()._assume_str()
        op = self.operator()
        if  op is operator.eq:
            m = 'equal(%s, %s)'%(l, r)
        elif op is operator.ne:
            m = 'notequal(%s, %s)'%(l, r)
        else:
            m = '(%s)%s(%s)' % (l, maxima._relation_symbols()[op], r)
        return m

    def _assume_str(self):
        """
        TESTS::

            sage: x = var('x')
            sage: x._assume_str()
            'x'
            sage: y = function('y', x)
            sage: y._assume_str()
            'y'
            sage: abs(x)._assume_str()
            'abs(x)'
        """
        # if this is a function with a single argument which is a symbol, i.e.
        # this is of the form f(x), we pass the string 'f > 0'
        if is_a_function(self._gobj) and self.nops() == 1 and \
                is_a_symbol(self._gobj.op(0)):
                    op = self.operator()
                    # check if op is a user defined function, for builtin
                    # functions like abs() we still need to pass 'abs(x) > 0'
                    if isinstance(op, SymbolicFunction):
                        return self.operator().name()
        return self._maxima_init_()


    _is_real     = deprecated_function_alias(10859, is_real)
    _is_positive = deprecated_function_alias(10859, is_positive)
    _is_negative = deprecated_function_alias(10859, is_negative)
    _is_integer  = deprecated_function_alias(10859, is_integer)
    _is_symbol   = deprecated_function_alias(10859, is_symbol)
    _is_constant = deprecated_function_alias(10859, is_constant)
    _is_numeric  = deprecated_function_alias(10859, is_numeric)

    def is_real(self):
        """
        Returns True if this expression is known to be a real number.

        EXAMPLES::

            sage: t0 = SR.symbol("t0", domain='real')
            sage: t0.is_real()
            True
            sage: t0.is_positive()
            False
            sage: t1 = SR.symbol("t1", domain='positive')
            sage: (t0+t1).is_real()
            True
            sage: (t0+x).is_real()
            False
            sage: (t0*t1).is_real()
            True
            sage: (t0*x).is_real()
            False

        The following is real, but we can't deduce that.::

            sage: (x*x.conjugate()).is_real()
            False
        """
        return self._gobj.info(info_real)

    def is_positive(self):
        """
        Returns True if this expression is known to be positive.

        EXAMPLES::

            sage: t0 = SR.symbol("t0", domain='positive')
            sage: t0.is_positive()
            True
            sage: t0.is_negative()
            False
            sage: t0.is_real()
            True
            sage: t1 = SR.symbol("t1", domain='positive')
            sage: (t0*t1).is_positive()
            True
            sage: (t0 + t1).is_positive()
            True
            sage: (t0*x).is_positive()
            False
        """
        return self._gobj.info(info_positive)

    def is_negative(self):
        """
        Return True if this expression is known to be negative.

        EXAMPLES::

            sage: SR(-5).is_negative()
            True

        Check if we can correctly deduce negativity of mul objects::

            sage: t0 = SR.symbol("t0", domain='positive')
            sage: t0.is_negative()
            False
            sage: (-t0).is_negative()
            True
            sage: (-pi).is_negative()
            True
        """
        return self._gobj.info(info_negative)

    def is_integer(self):
        """
        Return True if this expression is known to be an integer.

        EXAMPLES::

            sage: SR(5).is_integer()
            True
        """
        return self._gobj.info(info_integer)

    def is_symbol(self):
        """
        Return True if this symbolic expression consists of only a symbol, i.e.,
        a symbolic variable.

        EXAMPLES::

            sage: x.is_symbol()
            True
            sage: var('y')
            y
            sage: y.is_symbol()
            True
            sage: (x*y).is_symbol()
            False
            sage: pi.is_symbol()
            False

        ::

            sage: ((x*y)/y).is_symbol()
            True
            sage: (x^y).is_symbol()
            False
        """
        return is_a_symbol(self._gobj)

    def is_constant(self):
        """
        Return True if this symbolic expression is a constant.

        This function is intended to provide an interface to query the internal
        representation of the expression. In this sense, the word ``constant``
        doesn't reflect the mathematical properties of the expression.
        Expressions which have no variables may return ``False``.

        EXAMPLES::

            sage: pi.is_constant()
            True
            sage: x.is_constant()
            False
            sage: SR(1).is_constant()
            False

        Note that the complex I is not a constant::

            sage: I.is_constant()
            False
            sage: I.is_numeric()
            True
        """
        return is_a_constant(self._gobj)

    def is_numeric(self):
        """
        Return True if this expression only consists of a numeric object.

        EXAMPLES::

            sage: SR(1).is_numeric()
            True
            sage: x.is_numeric()
            False
            sage: pi.is_numeric()
            False
        """
        return is_a_numeric(self._gobj)

    def is_series(self):
        """
        Return True if ``self`` is a series.

        Series are special kinds of symbolic expressions that are
        constructed via the :meth:`series` method. They usually have
        an ``Order()`` term unless the series representation is exact,
        see :meth:`is_terminating_series`.

        OUTPUT:

        Boolean. Whether ``self`` is a series symbolic
        expression. Usually, this means that it was constructed by the
        :meth:`series` method.

        Returns ``False`` if only a subexpression of the symbolic
        expression is a series.

        EXAMPLES::

            sage: SR(5).is_series()
            False
            sage: var('x')
            x
            sage: x.is_series()
            False
            sage: exp(x).is_series()
            False
            sage: exp(x).series(x,10).is_series()
            True

        Laurent series are series, too::

            sage: laurent_series = (cos(x)/x).series(x, 5)
            sage: laurent_series
            1*x^(-1) + (-1/2)*x + 1/24*x^3 + Order(x^5)
            sage: laurent_series.is_series()
            True

        Something only containing a series as a subexpression is not a
        series::

            sage: sum_expr = 1 + exp(x).series(x,5); sum_expr
            (1 + 1*x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + Order(x^5)) + 1
            sage: sum_expr.is_series()
            False
        """
        return is_a_series(self._gobj)

    def is_terminating_series(self):
        """
        Return True if ``self`` is a series without order term.

        A series is terminating if it can be represented exactly,
        without requiring an order term. See also :meth:`is_series`
        for general series.

        OUTPUT:

        Boolean. Whether ``self`` was constructed by :meth:`series`
        and has no order term.

        EXAMPLES::

            sage: (x^5+x^2+1).series(x,10)
            1 + 1*x^2 + 1*x^5
            sage: (x^5+x^2+1).series(x,10).is_series()
            True
            sage: (x^5+x^2+1).series(x,10).is_terminating_series()
            True
            sage: SR(5).is_terminating_series()
            False
            sage: var('x')
            x
            sage: x.is_terminating_series()
            False
            sage: exp(x).series(x,10).is_terminating_series()
            False
        """
        return g_is_a_terminating_series(self._gobj)

    cpdef bint is_polynomial(self, var):
        """
        Return True if self is a polynomial in the given variable.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: t = x^2 + y; t
            x^2 + y
            sage: t.is_polynomial(x)
            True
            sage: t.is_polynomial(y)
            True
            sage: t.is_polynomial(z)
            True

            sage: t = sin(x) + y; t
            y + sin(x)
            sage: t.is_polynomial(x)
            False
            sage: t.is_polynomial(y)
            True
            sage: t.is_polynomial(sin(x))
            True

        TESTS:

        Check if we can handle derivatives. #6523::

            sage: f(x) = function('f',x)
            sage: f(x).diff(x).is_zero()
            False

        Check if #11352 is fixed::

            sage: el = -1/2*(2*x^2 - sqrt(2*x - 1)*sqrt(2*x + 1) - 1)
            sage: el.is_polynomial(x)
            False
        """
        cdef Expression symbol0 = self.coerce_in(var)
        return self._gobj.is_polynomial(symbol0._gobj)

    cpdef bint is_relational(self):
        """
        Return True if self is a relational expression.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.is_relational()
            True
            sage: sin(x).is_relational()
            False
        """
        return is_a_relational(self._gobj)

    cpdef bint is_infinity(self):
        """
        Return True if self is an infinite expression.

        EXAMPLES::

            sage: SR(oo).is_infinity()
            True
            sage: x.is_infinity()
            False
        """
        return is_a_infinity(self._gobj)

    def left_hand_side(self):
        """
        If self is a relational expression, return the left hand side
        of the relation.  Otherwise, raise a ValueError.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.left_hand_side()
            (x - 1)^2
            sage: eqn.lhs()
            (x - 1)^2
            sage: eqn.left()
            (x - 1)^2
        """
        if not self.is_relational():
            raise ValueError, "self must be a relational expression"
        return new_Expression_from_GEx(self._parent, self._gobj.lhs())

    lhs = left = left_hand_side

    def right_hand_side(self):
        """
        If self is a relational expression, return the right hand side
        of the relation.  Otherwise, raise a ValueError.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 <= x^2 - 2*x + 3
            sage: eqn.right_hand_side()
            x^2 - 2*x + 3
            sage: eqn.rhs()
            x^2 - 2*x + 3
            sage: eqn.right()
            x^2 - 2*x + 3
        """
        if not self.is_relational():
            raise ValueError, "self must be a relation"
        return new_Expression_from_GEx(self._parent, self._gobj.rhs())

    rhs = right = right_hand_side

    def is_trivial_zero(self):
        """
        Check if this expression is trivially equal to zero without any
        simplification.

        This method is intended to be used in library code where trying to
        obtain a mathematically correct result by applying potentially
        expensive rewrite rules is not desirable.

        EXAMPLES::

            sage: SR(0).is_trivial_zero()
            True
            sage: SR(0.0).is_trivial_zero()
            True
            sage: SR(float(0.0)).is_trivial_zero()
            True

            sage: (SR(1)/2^1000).is_trivial_zero()
            False
            sage: SR(1./2^10000).is_trivial_zero()
            False

        The :meth:`~sage.structure.element.Element.is_zero` method
        is more capable::

            sage: t = pi + (pi - 1)*pi - pi^2
            sage: t.is_trivial_zero()
            False
            sage: t.is_zero()
            True
            sage: u = sin(x)^2 + cos(x)^2 - 1
            sage: u.is_trivial_zero()
            False
            sage: u.is_zero()
            True
        """
        return self._gobj.is_zero()

    def __nonzero__(self):
        """
        Return True unless this symbolic expression can be shown by Sage
        to be zero.  Note that deciding if an expression is zero is
        undecidable in general.

        EXAMPLES::

            sage: x = var('x')
            sage: forget()
            sage: SR(0).__nonzero__()
            False
            sage: SR(1).__nonzero__()
            True
            sage: bool(abs(x))
            True
            sage: bool(x/x - 1)
            False

        This is called by :meth:`is_zero`::

            sage: k = var('k')
            sage: pol = 1/(k-1) - 1/k - 1/k/(k-1)
            sage: pol.is_zero()
            True

            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: f.is_zero()
            True

        TESTS:

        First, a bunch of tests of nonzero (which is called by bool)
        for symbolic relations::

            sage: x = var('x')
            sage: bool((x-1)^2 == x^2 - 2*x + 1)
            True
            sage: bool(((x-1)^2 == x^2 - 2*x + 1).expand())
            True
            sage: bool(((x-1)^2 == x^2 - 2*x + 3).expand())
            False
            sage: bool(2 + x < 3 + x)
            True
            sage: bool(2 + x < 1 + x)
            False
            sage: bool(2 + x > 1 + x)
            True
            sage: bool(1 + x > 1 + x)
            False
            sage: bool(1 + x >= 1 + x)
            True
            sage: bool(1 + x < 1 + x)
            False
            sage: bool(1 + x <= 1 + x)
            True
            sage: bool(1 + x^2 != 1 + x*x)
            False
            sage: bool(1 + x^2 != 2 + x*x)
            True
            sage: bool(SR(oo) == SR(oo))
            True
            sage: bool(-SR(oo) == SR(oo))
            False
            sage: bool(-SR(oo) != SR(oo))
            True

        Next, tests to ensure assumptions are correctly used::

            sage: x, y, z = var('x, y, z')
            sage: assume(x>=y,y>=z,z>=x)
            sage: bool(x==z)
            True
            sage: bool(z<x)
            False
            sage: bool(z>y)
            False
            sage: bool(y==z)
            True
            sage: bool(y<=z)
            True
            sage: forget()
            sage: assume(x>=1,x<=1)
            sage: bool(x==1)
            True
            sage: bool(x != 1)
            False
            sage: bool(x>1)
            False
            sage: forget()
            sage: assume(x>0)
            sage: bool(x==0)
            False
            sage: bool(x != 0)
            True
            sage: bool(x == 1)
            False

        The following must be true, even though we don't
        know for sure that x isn't 1, as symbolic comparisons
        elsewhere rely on x!=y unless we are sure it is not
        true; there is no equivalent of Maxima's ``unknown``.
        Since it is False that x==1, it is True that x != 1.

        ::

            sage: bool(x != 1)
            True
            sage: forget()
            sage: assume(x>y)
            sage: bool(x==y)
            False
            sage: bool(x != y)
            True
            sage: bool(x != y) # The same comment as above applies here as well
            True
            sage: forget()

        Comparisons of infinities::

            sage: assert( (1+I)*oo == (2+2*I)*oo )
            sage: assert( SR(unsigned_infinity) == SR(unsigned_infinity) )
            sage: assert( SR(I*oo) == I*oo )
            sage: assert( SR(-oo) <= SR(oo) )
            sage: assert( SR(oo) >= SR(-oo) )
            sage: assert( SR(oo) != SR(-oo) )
            sage: assert( sqrt(2)*oo != I*oo )
        """
        if self.is_relational():
            # constants are wrappers around Sage objects, compare directly
            if is_a_constant(self._gobj.lhs()) and is_a_constant(self._gobj.rhs()):
                return self.operator()(self.lhs().pyobject(), self.rhs().pyobject())

            pynac_result = relational_to_bool(self._gobj)

            # pynac is guaranteed to give the correct answer for comparing infinities
            if is_a_infinity(self._gobj.lhs()) or is_a_infinity(self._gobj.rhs()):
                return pynac_result

            if pynac_result:
                if self.operator() == operator.ne: # this hack is necessary to catch the case where the operator is != but is False because of assumptions made
                    m = self._maxima_()
                    s = m.parent()._eval_line('is (notequal(%s,%s))'%(repr(m.lhs()),repr(m.rhs())))
                    if s == 'false':
                        return False
                    else:
                        return True
                else:
                    return True

            # If assumptions are involved, falsification is more complicated...
            need_assumptions = False
            from sage.symbolic.assumptions import assumptions
            assumption_list = assumptions()
            if assumption_list:
                vars = self.variables()
                if vars:
                    assumption_var_list = []
                    for eqn in assumption_list:
                        try:
                            assumption_var_list.append(eqn.variables())
                        except AttributeError: # if we have a GenericDeclaration
                            assumption_var_list.append((eqn._var,))
                    assumption_vars = set(sum(assumption_var_list, ()))
                    if set(vars).intersection(assumption_vars):
                        need_assumptions = True

            # Use interval fields to try and falsify the relation
            if not need_assumptions:
                res = self.test_relation()
                if res is True:
                    return True
                elif res is False:
                    return False

            # we really have to do some work here...
            # I really don't like calling Maxima to test equality.  It
            # is SUPER SUPER SLOW, and it has all the problem
            # associated with different semantics, different
            # precision, etc., that can lead to subtle bugs.  Also, a
            # lot of basic Sage objects can't be put into maxima.
            from sage.symbolic.relation import test_relation_maxima
            return test_relation_maxima(self)

        self_is_zero = self._gobj.is_zero()
        if self_is_zero:
            return False
        else:
            return not bool(self == self._parent.zero_element())

    def test_relation(self, int ntests=20, domain=None, proof=True):
        """
        Test this relation at several random values, attempting to find
        a contradiction. If this relation has no variables, it will also
        test this relation after casting into the domain.

        Because the interval fields never return false positives, we can be
        assured that if True or False is returned (and proof is False) then
        the answer is correct.

        INPUT:

        - ``ntests`` -- (default ``20``) the number of iterations to run
        - ``domain`` -- (optional) the domain from which to draw the random
          values defaults to ``CIF`` for equality testing and ``RIF`` for
          order testing
        - ``proof`` -- (default ``True``) if ``False`` and the domain is an
          interval field, regard overlapping (potentially equal) intervals as
          equal, and return ``True`` if all tests succeeded.

        OUTPUT:

        Boolean or ``NotImplemented``, meaning

        - ``True`` -- this relation holds in the domain and has no variables.

        - ``False`` -- a contradiction was found.

        - ``NotImplemented`` -- no contradiction found.

        EXAMPLES::

            sage: (3 < pi).test_relation()
            True
            sage: (0 >= pi).test_relation()
            False
            sage: (exp(pi) - pi).n()
            19.9990999791895
            sage: (exp(pi) - pi == 20).test_relation()
            False
            sage: (sin(x)^2 + cos(x)^2 == 1).test_relation()
            NotImplemented
            sage: (sin(x)^2 + cos(x)^2 == 1).test_relation(proof=False)
            True
            sage: (x == 1).test_relation()
            False
            sage: var('x,y')
            (x, y)
            sage: (x < y).test_relation()
            False

        TESTS::

            sage: all_relations = [op for name, op in sorted(operator.__dict__.items()) if len(name) == 2]
            sage: all_relations
            [<built-in function eq>, <built-in function ge>, <built-in function gt>, <built-in function le>, <built-in function lt>, <built-in function ne>]
            sage: [op(3, pi).test_relation() for op in all_relations]
            [False, False, False, True, True, True]
            sage: [op(pi, pi).test_relation() for op in all_relations]
            [True, True, False, True, False, False]

            sage: s = 'some_very_long_variable_name_which_will_definitely_collide_if_we_use_a_reasonable_length_bound_for_a_hash_that_respects_lexicographic_order'
            sage: t1, t2 = var(','.join([s+'1',s+'2']))
            sage: (t1 == t2).test_relation()
            False
            sage: (cot(pi + x) == 0).test_relation()
            NotImplemented
        """
        cdef int k, eq_count = 0
        cdef bint is_interval
        if not self.is_relational():
            raise ValueError, "self must be a relation"
        cdef operators op = relational_operator(self._gobj)
        from sage.rings.real_mpfi import is_RealIntervalField
        from sage.rings.complex_interval_field import is_ComplexIntervalField
        from sage.rings.all import RIF, CIF
        if domain is None:
            is_interval = True
            if op == equal or op == not_equal:
                domain = CIF
            else:
                domain = RIF
        else:
            is_interval = is_RealIntervalField(domain) or is_ComplexIntervalField(domain)
        zero = domain(0)
        diff = self.lhs() - self.rhs()
        vars = diff.variables()
        if op == equal:
            falsify = operator.ne
        elif op == not_equal:
            falsify = operator.eq
        elif op == less:
            falsify = operator.ge
        elif op == less_or_equal:
            falsify = operator.gt
        elif op == greater:
            falsify = operator.le
        elif op == greater_or_equal:
            falsify = operator.lt
        cdef bint equality_ok = op in [equal, less_or_equal, greater_or_equal]
        cdef int errors = 0
        val = None
        if len(vars) == 0:
            try:
                val = domain(diff)
            except (TypeError, ValueError, ArithmeticError), ex:
                pass
            else:
                if self.operator()(val, zero):
                    return True
                elif falsify(val, zero):
                    return False
                if is_interval and not proof:
                    if val.contains_zero():
                        return equality_ok
                    else:
                        return not equality_ok
        else:
            for k in range(ntests):
                try:
                    if is_interval:
                        # Let's up the prec
                        if val and k > 4 and val.contains_zero() and domain.prec() < 1000:
                            domain = domain.to_prec(int(domain.prec() * 1.5))
                        # Uniform [-1,1] isn't the best distribution to use...
                        var_dict = dict([(v, domain.random_element() * domain.random_element(-2,6).exp()) for v in vars])
                    else:
                        var_dict = dict([(v, domain.random_element()) for v in vars])
                    val = domain(diff.subs(var_dict))
                    if falsify(val, zero):
                        return False
                    if is_interval:
                        eq_count += <bint>val.contains_zero()
                except (TypeError, ValueError, ArithmeticError, AttributeError), ex:
                    errors += 1
                    if k == errors > 3 and is_ComplexIntervalField(domain):
                        domain = RIF.to_prec(domain.prec())
                    # we are plugging in random values above, don't be surprised
                    # if something goes wrong...
                    eq_count += equality_ok

        if not proof:
            if not equality_ok:
                return eq_count == 0
            elif op == equal and is_interval:
                return eq_count == ntests
            else:
                return True
        # Nothing failed, so it *may* be True, but this method doesn't wasn't
        # able to find anything.
        return NotImplemented

    def negation(self):
        """
        Returns the negated version of self, that is the relation that is
        False iff self is True.

        EXAMPLES::

            sage: (x < 5).negation()
            x >= 5
            sage: (x == sin(3)).negation()
            x != sin(3)
            sage: (2*x >= sqrt(2)).negation()
            2*x < sqrt(2)
        """
        if not self.is_relational():
            raise ValueError, "self must be a relation"
        cdef operators op = relational_operator(self._gobj)
        if op == equal:
            falsify = operator.ne
        elif op == not_equal:
            falsify = operator.eq
        elif op == less:
            falsify = operator.ge
        elif op == less_or_equal:
            falsify = operator.gt
        elif op == greater:
            falsify = operator.le
        elif op == greater_or_equal:
            falsify = operator.lt
        return falsify(self.lhs(), self.rhs())

    def contradicts(self, soln):
        """
        Returns ``True`` if this relation is violated by the given variable assignment(s).

        EXAMPLES::

            sage: (x<3).contradicts(x==0)
            False
            sage: (x<3).contradicts(x==3)
            True
            sage: (x<=3).contradicts(x==3)
            False
            sage: y = var('y')
            sage: (x<y).contradicts(x==30)
            False
            sage: (x<y).contradicts({x: 30, y: 20})
            True
        """
        return bool(self.negation().subs(soln))

    def is_unit(self):
        """
        Return True if this expression is a unit of the symbolic ring.

        EXAMPLES::

            sage: SR(1).is_unit()
            True
            sage: SR(-1).is_unit()
            True
            sage: SR(0).is_unit()
            False
        """
        if not not self:
            return True
        if self == 0:
            return False
        raise NotImplementedError

    cdef Expression coerce_in(self, z):
        """
        Quickly coerce z to be an Expression.
        """
        cdef Expression w
        try:
            w = z
            return w
        except TypeError:
            return self._parent._coerce_(z)

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add left and right.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x + y + y + x
            2*x + 2*y

            # adding relational expressions
            sage: ( (x+y) > x ) + ( x > y )
            2*x + y > x + y

            sage: ( (x+y) > x ) + x
            2*x + y > 2*x

        TESTS::

            sage: x + ( (x+y) > x )
            2*x + y > 2*x

            sage: ( x > y) + (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: (x < 1) + (y <= 2)
            x + y < 3

            sage: x + oo
            +Infinity
            sage: x - oo
            -Infinity
            sage: x + unsigned_infinity
            Infinity
            sage: x - unsigned_infinity
            Infinity

            sage: nsr = x.parent()
            sage: nsr(oo) + nsr(oo)
            +Infinity
            sage: nsr(-oo) + nsr(-oo)
            -Infinity
            sage: nsr(oo) - nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity - infinity encountered.
            sage: nsr(-oo) - nsr(-oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity - infinity encountered.

            sage: nsr(unsigned_infinity) + nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity +- infinity encountered.
            sage: nsr(unsigned_infinity) - nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity +- infinity encountered.
            sage: nsr(oo) + nsr(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity +- infinity encountered.
            sage: nsr(oo) - nsr(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity +- infinity encountered.
            sage: nsr(unsigned_infinity) + nsr(unsigned_infinity)
            Infinity
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators op
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(gadd(left._gobj.lhs(), _right._gobj.lhs()),
                               gadd(left._gobj.rhs(), _right._gobj.rhs()),
                               op)
            else:
                x = relational(gadd(left._gobj.lhs(), _right._gobj),
                               gadd(left._gobj.rhs(), _right._gobj),
                               relational_operator(left._gobj))
        elif is_a_relational(_right._gobj):
            x = relational(gadd(left._gobj, _right._gobj.lhs()),
                           gadd(left._gobj, _right._gobj.rhs()),
                           relational_operator(_right._gobj))
        else:
            x = gadd(left._gobj, _right._gobj)
        return new_Expression_from_GEx(left._parent, x)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        """
        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x - y
            x - y

            # subtracting relational expressions
            sage: ( (x+y) > x ) - ( x > y )
            y > x - y

            sage: ( (x+y) > x ) - x
            y > 0

        TESTS::

            sage: x - ( (x+y) > x )
            -y > 0

            sage: ( x > y) - (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: x - oo
            -Infinity
            sage: oo - x
            +Infinity
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(gsub(left._gobj.lhs(), _right._gobj.lhs()),
                               gsub(left._gobj.rhs(), _right._gobj.rhs()),
                               op)
            else:
                x = relational(gsub(left._gobj.lhs(), _right._gobj),
                               gsub(left._gobj.rhs(), _right._gobj),
                               relational_operator(left._gobj))
        elif is_a_relational(_right._gobj):
            x = relational(gsub(left._gobj, _right._gobj.lhs()),
                           gsub(left._gobj, _right._gobj.rhs()),
                           relational_operator(_right._gobj))
        else:
            x = gsub(left._gobj, _right._gobj)
        return new_Expression_from_GEx(left._parent, x)

    cpdef RingElement _mul_(left, RingElement right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x*y*y
            x*y^2

            # multiplying relational expressions
            sage: ( (x+y) > x ) * ( x > y )
            (x + y)*x > x*y

            sage: ( (x+y) > x ) * x
            (x + y)*x > x^2

            sage: ( (x+y) > x ) * -1
            -x - y > -x

        TESTS::

            sage: x * ( (x+y) > x )
            (x + y)*x > x^2

            sage: ( x > y) * (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: a = 1000 + 300*x + x^3 + 30*x^2
            sage: a*Mod(1,7)
            x^3 + 2*x^2 + 6*x + 6

            sage: var('z')
            z
            sage: 3*(x+y)/z
            3*(x + y)/z
            sage: (-x+z)*(3*x-3*z)
            -3*(x - z)^2

            # check if comparison of constant terms in Pynac add objects work
            sage: (y-1)*(y-2)
            (y - 1)*(y - 2)

        Check if Pynac can compute inverses of Python longs (:trac:`13107`)::

            sage: SR(4L)*SR(2L)^(-1)
            2.0

        Check for simplifications when multiplying instances of exp::

            sage: exp(x)*exp(y)
            e^(x + y)
            sage: exp(x)^2*exp(y)
            e^(2*x + y)
            sage: x^y*exp(x+y)*exp(-y)
            x^y*e^x
            sage: x^y*exp(x+y)*(x+y)*(2*x+2*y)*exp(-y)
            2*(x + y)^2*x^y*e^x
            sage: x^y*exp(x+y)*(x+y)*(2*x+2*y)*exp(-y)*exp(z)^2
            2*(x + y)^2*x^y*e^(x + 2*z)
            sage: 1/exp(x)
            e^(-x)
            sage: exp(x)/exp(y)
            e^(x - y)
            sage: A = exp(I*pi/5)
            sage: t = A*A*A*A; t
            e^(4/5*I*pi)
            sage: t*A
            -1
            sage: b = -x*A; c = b*b; c
            x^2*e^(2/5*I*pi)
            sage: u = -t*A; u
            1

        Products of non integer powers of exp are not simplified::

            sage: exp(x)^I*exp(z)^(2.5)
            (e^x)^I*(e^z)^2.50000000000000

        ::

            sage: x*oo
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity * f(x) encountered.
            sage: x*unsigned_infinity
            Traceback (most recent call last):
            ...
            ValueError: oo times number < oo not defined

            sage: SR(oo)*SR(oo)
            +Infinity
            sage: SR(-oo)*SR(oo)
            -Infinity
            sage: SR(oo)*SR(-oo)
            -Infinity
            sage: SR(unsigned_infinity)*SR(oo)
            Infinity

        Check if we are returning informative error messages in case of
        nonsensical arithmetic :trac:`13739`::

            sage: t = GF(5)(3)
            sage: u = GF(7)(4)
            sage: var('y')
            y
            sage: e = t*x + u*y
            sage: t*e
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Finite Field
            of size 7' and 'Finite Field of size 5'

        The same issue (with a different test case) was reported in
        :trac:`10960`::

            sage: K.<b> = FiniteField(9)
            sage: i*b
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Number Field
            in I with defining polynomial x^2 + 1' and 'Finite Field in b of
            size 3^2'

        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(gmul(left._gobj.lhs(), _right._gobj.lhs()),
                               gmul(left._gobj.rhs(), _right._gobj.rhs()),
                               op)
            else:
                o = relational_operator(left._gobj)
                x = relational(gmul(left._gobj.lhs(), _right._gobj),
                               gmul(left._gobj.rhs(), _right._gobj),
                               o)
        elif is_a_relational(_right._gobj):
            o = relational_operator(_right._gobj)
            x = relational(gmul(left._gobj, _right._gobj.lhs()),
                           gmul(left._gobj, _right._gobj.rhs()),
                           o)
        else:
            x = gmul(left._gobj, _right._gobj)
        return new_Expression_from_GEx(left._parent, x)

    cpdef RingElement _div_(left, RingElement right):
        """
        Divide left and right.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x/y/y
            x/y^2

            # dividing relational expressions
            sage: ( (x+y) > x ) / ( x > y )
            (x + y)/x > x/y

            sage: ( (x+y) > x ) / x
            (x + y)/x > 1

            sage: ( (x+y) > x ) / -1
            -x - y > -x

        TESTS::

            sage: x / ( (x+y) > x )
            x/(x + y) > 1

            sage: ( x > y) / (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations
            sage: x/oo
            0
            sage: oo/x
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity * f(x) encountered.

            sage: SR(oo)/SR(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0 * infinity encountered.

            sage: SR(-oo)/SR(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0 * infinity encountered.

            sage: SR(oo)/SR(-oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0 * infinity encountered.

            sage: SR(oo)/SR(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0 * infinity encountered.

            sage: SR(unsigned_infinity)/SR(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0 * infinity encountered.

            sage: SR(0)/SR(oo)
            0

            sage: SR(0)/SR(unsigned_infinity)
            0

            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Symbolic division by zero

        Check if Pynac can compute divisions of Python longs (:trac:`13107`)::

            sage: SR(1L)/SR(2L)
            0.5
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        try:
            if is_a_relational(left._gobj):
                if is_a_relational(_right._gobj):
                    op = compatible_relation(relational_operator(left._gobj),
                                             relational_operator(_right._gobj))
                    x = relational(gdiv(left._gobj.lhs(), _right._gobj.lhs()),
                                   gdiv(left._gobj.rhs(), _right._gobj.rhs()),
                                   op)
                else:
                    o = relational_operator(left._gobj)
                    x = relational(gdiv(left._gobj.lhs(), _right._gobj),
                                   gdiv(left._gobj.rhs(), _right._gobj),
                                   o)
            elif is_a_relational(_right._gobj):
                o = relational_operator(_right._gobj)
                x = relational(gdiv(left._gobj, _right._gobj.lhs()),
                               gdiv(left._gobj, _right._gobj.rhs()),
                               o)
            else:
                x = gdiv(left._gobj, _right._gobj)
            return new_Expression_from_GEx(left._parent, x)
        except Exception, msg:
            # TODO: change this to maybe cleverly do something involving Cython C++ exception handling.
            # See http://docs.cython.org/docs/wrapping_CPlusPlus.html
            if 'division by zero' in str(msg):
                raise ZeroDivisionError, "Symbolic division by zero"
            else:
                raise

    def __invert__(self):
        """
        Return the inverse of this symbolic expression.

        EXAMPLES::

            sage: ~x
            1/x
            sage: ~SR(3)
            1/3
            sage: v1=var('v1'); a = (2*erf(2*v1*arcsech(0))/v1); ~a
            1/2*v1/erf(2*v1*arcsech(0))
        """
        return 1/self

    # Boilerplate code from sage/structure/element.pyx
    def __cmp__(left, right):
        """
        Compare self and right, returning -1, 0, or 1, depending on if
        self < right, self == right, or self > right, respectively.

        Use this instead of the operators <=, <, etc. to compare symbolic
        expressions when you do not want to get a formal inequality back.

        IMPORTANT: Both self and right *must* have the same type, or
        this function won't be called.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: x.__cmp__(y)
            1
            sage: x < y
            x < y
            sage: cmp(x,y)
            1
            sage: cmp(SR(0.5), SR(0.7))
            -1
            sage: SR(0.5) < SR(0.7)
            0.500000000000000 < 0.700000000000000
            sage: cmp(SR(0.5), 0.7)
            -1
            sage: cmp(sin(SR(2)), sin(SR(1)))
            1
            sage: float(sin(SR(2)))
            0.9092974268256817
            sage: float(sin(SR(1)))
            0.8414709848078965

        TESTS:

        Check that :trac:`9880` is fixed::

            sage: b = [var('b_%s'%i) for i in range(4)]
            sage: precomp = (2^b_2 + 2)*(2^b_1 + 2^(-b_1) + 2^b_1*2^b_0 - \
                        2^b_1*2^(-b_0) - 2^(-b_1)*2^b_0 - 2^(-b_1)*2^(-b_0) + \
                        2^b_0 + 2^(-b_0) - 9) + (2^b_1 + 2^(-b_1) + \
                        2^b_1*2^b_0 - 2^b_1*2^(-b_0) - 2^(-b_1)*2^b_0 - \
                         2^(-b_1)*2^(-b_0) + 2^b_0 + 2^(-b_0) - 9)/2^b_2
            sage: repl_dict = {b_0: b_0, b_3: b_1, b_2: b_3, b_1: b_2}
            sage: P = precomp.substitute(repl_dict)
            sage: P.expand()
            -2^(-b_0)*2^(-b_2)*2^b_3 - 2^b_0*2^(-b_2)*2^b_3 -
            2^(-b_0)*2^b_2*2^b_3 + 2^b_0*2^b_2*2^b_3 - 2*2^(-b_0)*2^(-b_2)
            - 2*2^b_0*2^(-b_2) - 2*2^(-b_0)*2^b_2 + 2*2^b_0*2^b_2 +
            2^(-b_0)*2^b_3 + 2^b_0*2^b_3 + 2^(-b_2)*2^b_3 + 2^b_2*2^b_3 +
            2*2^(-b_0) + 2*2^b_0 + 2*2^(-b_2) + 2*2^b_2 - 9*2^b_3 -
            2^(-b_0)*2^(-b_2)/2^b_3 - 2^b_0*2^(-b_2)/2^b_3 -
            2^(-b_0)*2^b_2/2^b_3 + 2^b_0*2^b_2/2^b_3 + 2^(-b_0)/2^b_3 +
            2^b_0/2^b_3 + 2^(-b_2)/2^b_3 + 2^b_2/2^b_3 - 9/2^b_3 - 18

            sage: _0,b_1,b_2=var('b_0,b_1,b_2')
            sage: f = 1/27*b_2^2/(2^b_2)^2 + 1/27*b_1^2/(2^b_1)^2 + \
            1/27*b_0^2/(2^b_0)^2 + 1/27*b_2/(2^b_2)^2 - 2/81/(2^b_2)^2 + \
            1/27*b_1/(2^b_1)^2 + 8/243/(2^b_2)^2 - 1/81*b_0/(2^b_0)^2 - \
            1/27*b_1^2/((2^b_2)^2*(2^b_1)^2) - \
            1/27*b_0^2/((2^b_2)^2*(2^b_0)^2) - 20/243/(2^b_1)^2 + 1/9/2^b_0 \
            + 4/81*b_0/(2^b_0)^2 - 8/243/(2^b_2)^2 - 2/9/(2^b_2*2^b_1) - \
            2/9/(2^b_2*2^b_0) + 8/243/(2^b_1)^2 - 1/9/2^b_0 + \
            2/9/(2^b_2*2^b_1) + 2/9/(2^b_2*2^b_0) - \
            2/27*b_1*b_2/((2^b_2)^2*(2^b_1)^2) - \
            1/27*b_2^2/((2^b_2)^2*(2^b_1)^2) - \
            2/27*b_0*b_2/((2^b_2)^2*(2^b_0)^2) - \
            1/27*b_2^2/((2^b_2)^2*(2^b_0)^2) + 2/81/(2^b_1)^2 - \
            1/27*b_0^2/((2^b_1)^2*(2^b_0)^2) - \
            2/27*b_0*b_1/((2^b_1)^2*(2^b_0)^2) - \
            1/27*b_1^2/((2^b_1)^2*(2^b_0)^2) - 2/81/(2^b_0)^2 + \
            5/27*b_1/((2^b_2)^2*(2^b_1)^2) + 5/27*b_2/((2^b_2)^2*(2^b_1)^2) \
            + 5/27*b_0/((2^b_2)^2*(2^b_0)^2) + \
            5/27*b_2/((2^b_2)^2*(2^b_0)^2) + 5/27*b_0/((2^b_1)^2*(2^b_0)^2) \
            + 5/27*b_1/((2^b_1)^2*(2^b_0)^2) - 4/81/((2^b_2)^2*(2^b_1)^2) + \
            1/27*b_0^2/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + \
            2/27*b_0*b_1/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + \
            2/27*b_0*b_2/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + \
            1/27*b_1^2/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + \
            2/27*b_1*b_2/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + \
            1/27*b_2^2/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) - \
            4/81/((2^b_2)^2*(2^b_0)^2) - 4/81/((2^b_1)^2*(2^b_0)^2) - \
            11/27*b_0/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) - \
            11/27*b_1/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) - \
            11/27*b_2/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + \
            64/81/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + 35/81 \
            sage: f.nops()
            38

            sage: x,y,z = var('x y z');
            sage: print (-x+z)*(3*x-3*z)
            -3*(x - z)^2

            sage: t = var('t')
            sage: (x-t)^3
            -(t - x)^3
            sage: (-t+x)^3
            -(t - x)^3
            sage: (-x+t)^3
            (t - x)^3

        This example is from :trac:`10833`::

            sage: R.<x,c> = PolynomialRing(QQ,2)
            sage: phi(x) = x^2 + c
            sage: def iterkate(n):
            ....:     pol = x
            ....:     for i in range(1,n):
            ....:         pol = phi(pol)
            ....:     return pol
            ....:
            sage: g = expand(iterkate(7))
            sage: g.nops()
            480
        """
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compare ``left`` and ``right``.

        INPUT:

        - ``right`` -- A :class:`Expression` instance.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: a = sqrt(3)
            sage: b = x^2+1
            sage: a.__cmp__(b)   # indirect doctest
            -1
        """
        return print_order_compare(left._gobj, (<Expression>right)._gobj)

    cpdef int _cmp_add(Expression left, Expression right) except -2:
        """
        Compare ``left`` and ``right`` in the print order.

        INPUT:

        - ``right`` -- A :class:`Expression` instance.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: a = sqrt(3)
            sage: b = x^2+1
            sage: a._cmp_add(b)
            -1
            sage: b._cmp_add(a)
            1
            sage: b._cmp_add(1)
            Traceback (most recent call last):
            ...
            TypeError: Argument 'right' has incorrect type (expected
            sage.symbolic.expression.Expression, got sage.rings.integer.Integer)
        """
        return print_order_compare(left._gobj, right._gobj)

    cpdef int _cmp_mul(Expression left, Expression right) except -2:
        """
        Compare ``left`` and ``right`` in the print order for products.

        INPUT:

        - ``right`` -- A :class:`Expression` instance.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: a = sqrt(3)
            sage: b = x^2+1
            sage: a._cmp_mul(b)
            -1
            sage: b._cmp_mul(a)
            1
            sage: b._cmp_mul(1)
            Traceback (most recent call last):
            ...
            TypeError: Argument 'right' has incorrect type (expected
            sage.symbolic.expression.Expression, got sage.rings.integer.Integer)
        """
        return print_order_compare_mul(left._gobj, right._gobj)

    def __pow__(self, exp, ignored):
        """
        Return self raised to the power of exp.

        INPUT:

        - ``exp`` -- something that coerces to a symbolic expressions.
        - ``ignored`` -- the second argument that should accept a modulus
          is actually ignored.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: x.__pow__(y)
            x^y
            sage: x^(3/5)
            x^(3/5)
            sage: x^sin(x)^cos(y)
            x^(sin(x)^cos(y))

        TESTS::

            sage: (Mod(2,7)*x^2 + Mod(2,7))^7
            (2*x^2 + 2)^7

        The leading coefficient in the result above is 1 since::

            sage: t = Mod(2,7); gcd(t, t)^7
            1
            sage: gcd(t,t).parent()
            Ring of integers modulo 7

        ::

            sage: k = GF(7)
            sage: f = expand((k(1)*x^5 + k(1)*x^2 + k(2))^7); f
            x^35 + x^14 + 2

            sage: x^oo
            Traceback (most recent call last):
            ...
            ValueError: power::eval(): pow(f(x), infinity) is not defined.
            sage: SR(oo)^2
            +Infinity
            sage: SR(-oo)^2
            +Infinity
            sage: SR(-oo)^3
            -Infinity
            sage: SR(unsigned_infinity)^2
            Infinity

        Test powers of exp::

            sage: exp(2)^5
            e^10
            sage: exp(x)^5
            e^(5*x)

        Test base a Python numeric type::

            sage: int(2)^x
            2^x
            sage: float(2.3)^(x^3 - x^2 + 1/3)
            2.3^(x^3 - x^2 + 1/3)
            sage: complex(1,3)^(sqrt(2))
            (1+3j)^sqrt(2)

        Test complex numeric powers::

            sage: I^0.5
            0.707106781186548 + 0.707106781186547*I
            sage: (I + 1) ^ (0.5 + I)
            0.400667052375828 + 0.365310866736929*I
            sage: I^I
            I^I
            sage: I^x
            I^x
            sage: I^(1/2)
            sqrt(I)
            sage: I^(2/3)
            I^(2/3)
            sage: 2^(1/2)
            sqrt(2)
            sage: (2*I)^(1/2)
            sqrt(2*I)

        Test if we can take powers of elements of Q(i) #8659::

            sage: t = I.pyobject().parent()(8)
            sage: t^(1/2)
            2*sqrt(2)
            sage: (t^2)^(1/4)
            2*4^(1/4)

        Test if we can compute inverses of Python longs (:trac:`13107`)::

            sage: SR(2L)^(-1)
            0.5
        """
        cdef Expression base, nexp

        try:
            # self is an Expression and exp might not be
            base = self
            nexp = base.coerce_in(exp)
        except TypeError:
            # exp is an Expression and self might not be
            nexp = exp
            base = nexp.coerce_in(self)
        cdef GEx x
        if is_a_relational(base._gobj):
            x = relational(g_pow(base._gobj.lhs(), nexp._gobj),
                           g_pow(base._gobj.rhs(), nexp._gobj),
                           relational_operator(base._gobj))
        else:
            x = g_pow(base._gobj, nexp._gobj)
        return new_Expression_from_GEx(base._parent, x)

    def derivative(self, *args):
        """
        Returns the derivative of this expressions with respect to the
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global
        :meth:`~sage.calculus.functional.derivative` function for more
        details.

        .. seealso::

            This is implemented in the `_derivative` method (see the
            source code).

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: t = (x^2+y)^2
            sage: t.derivative(x)
            4*(x^2 + y)*x
            sage: t.derivative(x, 2)
            12*x^2 + 4*y
            sage: t.derivative(x, 2, y)
            4
            sage: t.derivative(y)
            2*x^2 + 2*y

        If the function depends on only one variable, you may omit the
        variable. Giving just a number (for the order of the derivative)
        also works::

            sage: f(x) = x^3 + sin(x)
            sage: f.derivative()
            x |--> 3*x^2 + cos(x)
            sage: f.derivative(2)
            x |--> 6*x - sin(x)

        ::

            sage: t = sin(x+y^2)*tan(x*y)
            sage: t.derivative(x)
            (tan(x*y)^2 + 1)*y*sin(y^2 + x) + cos(y^2 + x)*tan(x*y)
            sage: t.derivative(y)
            (tan(x*y)^2 + 1)*x*sin(y^2 + x) + 2*y*cos(y^2 + x)*tan(x*y)

        ::

            sage: h = sin(x)/cos(x)
            sage: derivative(h,x,x,x)
            8*sin(x)^2/cos(x)^2 + 6*sin(x)^4/cos(x)^4 + 2
            sage: derivative(h,x,3)
            8*sin(x)^2/cos(x)^2 + 6*sin(x)^4/cos(x)^4 + 2

        ::

            sage: var('x, y')
            (x, y)
            sage: u = (sin(x) + cos(y))*(cos(x) - sin(y))
            sage: derivative(u,x,y)
            -cos(x)*cos(y) + sin(x)*sin(y)
            sage: f = ((x^2+1)/(x^2-1))^(1/4)
            sage: g = derivative(f, x); g # this is a complex expression
            -1/2*((x^2 + 1)*x/(x^2 - 1)^2 - x/(x^2 - 1))/((x^2 + 1)/(x^2 - 1))^(3/4)
            sage: g.factor()
            -x/((x + 1)^2*(x - 1)^2*((x^2 + 1)/(x^2 - 1))^(3/4))

        ::

            sage: y = var('y')
            sage: f = y^(sin(x))
            sage: derivative(f, x)
            y^sin(x)*cos(x)*log(y)

        ::

            sage: g(x) = sqrt(5-2*x)
            sage: g_3 = derivative(g, x, 3); g_3(2)
            -3

        ::

            sage: f = x*e^(-x)
            sage: derivative(f, 100)
            x*e^(-x) - 100*e^(-x)

        ::

            sage: g = 1/(sqrt((x^2-1)*(x+5)^6))
            sage: derivative(g, x)
            -((x + 5)^6*x + 3*(x^2 - 1)*(x + 5)^5)/((x^2 - 1)*(x + 5)^6)^(3/2)

        TESTS::

            sage: t.derivative()
            Traceback (most recent call last):
            ...
            ValueError: No differentiation variable specified.
        """
        return multi_derivative(self, args)

    diff = differentiate = derivative

    def _derivative(self, symb=None, deg=1):
        """
        Return the deg-th (partial) derivative of self with respect to symb.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: b = (x+y)^5
            sage: b._derivative(x, 2)
            20*(x + y)^3

            sage: foo = function('foo',nargs=2)
            sage: foo(x^2,x^2)._derivative(x)
            2*x*D[0](foo)(x^2, x^2) + 2*x*D[1](foo)(x^2, x^2)

            sage: SR(1)._derivative()
            0

        If the expression is a callable symbolic expression, and no
        variables are specified, then calculate the gradient::

            sage: f(x,y)=x^2+y
            sage: f.diff() # gradient
            (x, y) |--> (2*x, 1)

        TESTS:

        Raise error if no variable is specified and there are multiple
        variables::

            sage: b._derivative()
            Traceback (most recent call last):
            ...
            ValueError: No differentiation variable specified.

        Check if #6524 is fixed::

            sage: f = function('f')
            sage: f(x)*f(x).derivative(x)*f(x).derivative(x,2)
            f(x)*D[0](f)(x)*D[0, 0](f)(x)
            sage: g = f(x).diff(x)
            sage: h = f(x).diff(x)*sin(x)
            sage: h/g
            sin(x)
        """
        if symb is None:
            # we specify a default value of None for symb and check for it here
            # to return more helpful error messages when no variable is
            # given by the multi_derivative framework
            vars = self.variables()
            if len(vars) == 1:
                symb = vars[0]
            elif len(vars) == 0:
                return self._parent(0)
            elif sage.symbolic.callable.is_CallableSymbolicExpression(self):
                return self.gradient()
            else:
                raise ValueError, "No differentiation variable specified."
        if not isinstance(deg, (int, long, sage.rings.integer.Integer)) \
                or deg < 1:
            raise TypeError, "argument deg should be an integer >= 1."
        cdef Expression symbol = self.coerce_in(symb)
        if not is_a_symbol(symbol._gobj):
            raise TypeError, "argument symb must be a symbol"
        cdef GEx x
        sig_on()
        try:
            x = self._gobj.diff(ex_to_symbol(symbol._gobj), deg)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def gradient(self, variables=None):
        r"""
        Compute the gradient of a symbolic function.

        This function returns a vector whose components are the derivatives
        of the original function with respect to the arguments of the
        original function. Alternatively, you can specify the variables as
        a list.

        EXAMPLES::

            sage: x,y = var('x y')
            sage: f = x^2+y^2
            sage: f.gradient()
            (2*x, 2*y)
            sage: g(x,y) = x^2+y^2
            sage: g.gradient()
            (x, y) |--> (2*x, 2*y)
            sage: n = var('n')
            sage: f(x,y) = x^n+y^n
            sage: f.gradient()
            (x, y) |--> (n*x^(n - 1), n*y^(n - 1))
            sage: f.gradient([y,x])
            (x, y) |--> (n*y^(n - 1), n*x^(n - 1))
        """
        from sage.modules.free_module_element import vector
        if variables is None:
            variables = self.arguments()
        return vector([self.derivative(x) for x in variables])

    def hessian(self):
        r"""
        Compute the hessian of a function. This returns a matrix components
        are the 2nd partial derivatives of the original function.

        EXAMPLES::

            sage: x,y = var('x y')
            sage: f = x^2+y^2
            sage: f.hessian()
            [2 0]
            [0 2]
            sage: g(x,y) = x^2+y^2
            sage: g.hessian()
            [(x, y) |--> 2 (x, y) |--> 0]
            [(x, y) |--> 0 (x, y) |--> 2]
        """
        from sage.matrix.constructor import matrix
        return matrix([[g.derivative(x) for x in self.arguments()]
                       for g in self.gradient()])


    def series(self, symbol, int order):
        r"""
        Return the power series expansion of self in terms of the
        given variable to the given order.

        INPUT:

        - ``symbol`` - a symbolic variable or symbolic equality
          such as ``x == 5``; if an equality is given, the
          expansion is around the value on the right hand side
          of the equality
        - ``order`` - an integer

        OUTPUT:

        A power series.

        To truncate the power series and obtain a normal expression, use the
        :meth:`truncate` command.

        EXAMPLES:

        We expand a polynomial in `x` about 0, about `1`, and also truncate
        it back to a polynomial::

            sage: var('x,y')
            (x, y)
            sage: f = (x^3 - sin(y)*x^2 - 5*x + 3); f
            x^3 - x^2*sin(y) - 5*x + 3
            sage: g = f.series(x, 4); g
            3 + (-5)*x + (-sin(y))*x^2 + 1*x^3
            sage: g.truncate()
            x^3 - x^2*sin(y) - 5*x + 3
            sage: g = f.series(x==1, 4); g
            (-sin(y) - 1) + (-2*sin(y) - 2)*(x - 1) + (-sin(y) + 3)*(x - 1)^2 + 1*(x - 1)^3
            sage: h = g.truncate(); h
            (x - 1)^3 - (x - 1)^2*(sin(y) - 3) - 2*(x - 1)*(sin(y) + 1) - sin(y) - 1
            sage: h.expand()
            x^3 - x^2*sin(y) - 5*x + 3

        We computer another series expansion of an analytic function::

            sage: f = sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x==1,3)
            (sin(1)) + (cos(1) - 2*sin(1))*(x - 1) + (-2*cos(1) + 5/2*sin(1))*(x - 1)^2 + Order((x - 1)^3)
            sage: f.series(x==1,3).truncate().expand()
            -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)

        Following the GiNaC tutorial, we use John Machin's amazing
        formula `\pi = 16 \tan^{-1}(1/5) - 4 \tan^{-1}(1/239)` to compute
        digits of `\pi`. We expand the arc tangent around 0 and insert
        the fractions 1/5 and 1/239.

        ::

            sage: x = var('x')
            sage: f = atan(x).series(x, 10); f
            1*x + (-1/3)*x^3 + 1/5*x^5 + (-1/7)*x^7 + 1/9*x^9 + Order(x^10)
            sage: float(16*f.subs(x==1/5) - 4*f.subs(x==1/239))
            3.1415926824043994

        TESTS:

        Check if #8943 is fixed::

            sage: ((1+arctan(x))**(1/x)).series(x==0, 3)
            (e) + (-1/2*e)*x + (1/8*e)*x^2 + Order(x^3)
        """
        cdef Expression symbol0 = self.coerce_in(symbol)
        cdef GEx x
        sig_on()
        try:
            x = self._gobj.series(symbol0._gobj, order, 0)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def taylor(self, *args):
        r"""
        Expands this symbolic expression in a truncated Taylor or
        Laurent series in the variable `v` around the point `a`,
        containing terms through `(x - a)^n`. Functions in more
        variables is also supported.

        INPUT:

        -  ``*args`` - the following notation is supported

           - ``x, a, n`` - variable, point, degree

           - ``(x, a), (y, b), n`` - variables with points, degree of polynomial

        EXAMPLES::

            sage: var('a, x, z')
            (a, x, z)
            sage: taylor(a*log(z), z, 2, 3)
            1/24*a*(z - 2)^3 - 1/8*a*(z - 2)^2 + 1/2*a*(z - 2) + a*log(2)

        ::

            sage: taylor(sqrt (sin(x) + a*x + 1), x, 0, 3)
            1/48*(3*a^3 + 9*a^2 + 9*a - 1)*x^3 - 1/8*(a^2 + 2*a + 1)*x^2 + 1/2*(a + 1)*x + 1

        ::

            sage: taylor (sqrt (x + 1), x, 0, 5)
            7/256*x^5 - 5/128*x^4 + 1/16*x^3 - 1/8*x^2 + 1/2*x + 1

        ::

            sage: taylor (1/log (x + 1), x, 0, 3)
            -19/720*x^3 + 1/24*x^2 - 1/12*x + 1/x + 1/2

        ::

            sage: taylor (cos(x) - sec(x), x, 0, 5)
            -1/6*x^4 - x^2

        ::

            sage: taylor ((cos(x) - sec(x))^3, x, 0, 9)
            -1/2*x^8 - x^6

        ::

            sage: taylor (1/(cos(x) - sec(x))^3, x, 0, 5)
            -15377/7983360*x^4 - 6767/604800*x^2 + 11/120/x^2 + 1/2/x^4 - 1/x^6 - 347/15120

        TESTS:

        Check that ticket #7472 is fixed (Taylor polynomial in more variables)::

            sage: x,y=var('x y'); taylor(x*y^3,(x,1),(y,1),4)
            (x - 1)*(y - 1)^3 + 3*(x - 1)*(y - 1)^2 + (y - 1)^3 + 3*(x - 1)*(y - 1) + 3*(y - 1)^2 + x + 3*y - 3
            sage: expand(_)
            x*y^3

        """
        from sage.all import SR, Integer
        A=args
        try:
            if isinstance(A[0],tuple):
                B=[]
                B.append([SR(A[i][0]) for i in range(len(A)-1)])
                B.append([A[i][1] for i in range(len(A)-1)])
            else:
                B=[A[0],SR(A[1])]
            B.append(Integer(A[len(A)-1]))
        except StandardError:
            raise NotImplementedError, "Wrong arguments passed to taylor. See taylor? for more details."
        l = self._maxima_().taylor(B)
        return self.parent()(l)



    def truncate(self):
        """
        Given a power series or expression, return the corresponding
        expression without the big oh.

        INPUT:

        - ``self`` -- a series as output by the :meth:`series` command.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: f = sin(x)/x^2
            sage: f.truncate()
            sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x,7).truncate()
            -1/5040*x^5 + 1/120*x^3 - 1/6*x + 1/x
            sage: f.series(x==1,3).truncate().expand()
            -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)
        """
        if not is_a_series(self._gobj):
            return self
        return new_Expression_from_GEx(self._parent, series_to_poly(self._gobj))

    def expand(Expression self, side=None):
        """
        Expand this symbolic expression. Products of sums and exponentiated
        sums are multiplied out, numerators of rational expressions which
        are sums are split into their respective terms, and multiplications
        are distributed over addition at all levels.

        EXAMPLES:

        We expand the expression `(x-y)^5` using both
        method and functional notation.

        ::

            sage: x,y = var('x,y')
            sage: a = (x-y)^5
            sage: a.expand()
            x^5 - 5*x^4*y + 10*x^3*y^2 - 10*x^2*y^3 + 5*x*y^4 - y^5
            sage: expand(a)
            x^5 - 5*x^4*y + 10*x^3*y^2 - 10*x^2*y^3 + 5*x*y^4 - y^5

        We expand some other expressions::

            sage: expand((x-1)^3/(y-1))
            x^3/(y - 1) - 3*x^2/(y - 1) + 3*x/(y - 1) - 1/(y - 1)
            sage: expand((x+sin((x+y)^2))^2)
            x^2 + 2*x*sin((x + y)^2) + sin((x + y)^2)^2

        We can expand individual sides of a relation::

            sage: a = (16*x-13)^2 == (3*x+5)^2/2
            sage: a.expand()
            256*x^2 - 416*x + 169 == 9/2*x^2 + 15*x + 25/2
            sage: a.expand('left')
            256*x^2 - 416*x + 169 == 1/2*(3*x + 5)^2
            sage: a.expand('right')
            (16*x - 13)^2 == 9/2*x^2 + 15*x + 25/2

        TESTS::

            sage: var('x,y')
            (x, y)
            sage: ((x + (2/3)*y)^3).expand()
            x^3 + 2*x^2*y + 4/3*x*y^2 + 8/27*y^3
            sage: expand( (x*sin(x) - cos(y)/x)^2 )
            x^2*sin(x)^2 - 2*cos(y)*sin(x) + cos(y)^2/x^2
            sage: f = (x-y)*(x+y); f
            (x + y)*(x - y)
            sage: f.expand()
            x^2 - y^2
        """
        if side is not None:
            if not is_a_relational(self._gobj):
                raise ValueError, "expansion on sides only makes sense for relations"
            if side == 'left':
                return self.operator()(self.lhs().expand(), self.rhs())
            elif side == 'right':
                return self.operator()(self.lhs(), self.rhs().expand())
            else:
                raise ValueError, "side must be 'left', 'right', or None"

        cdef GEx x
        sig_on()
        try:
            x = self._gobj.expand(0)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    expand_rational = rational_expand = expand

    def expand_trig(self, full=False, half_angles=False, plus=True, times=True):
        """
        Expands trigonometric and hyperbolic functions of sums of angles
        and of multiple angles occurring in self. For best results, self
        should already be expanded.

        INPUT:

        -  ``full`` - (default: False) To enhance user control
           of simplification, this function expands only one level at a time
           by default, expanding sums of angles or multiple angles. To obtain
           full expansion into sines and cosines immediately, set the optional
           parameter full to True.

        -  ``half_angles`` - (default: False) If True, causes
           half-angles to be simplified away.

        -  ``plus`` - (default: True) Controls the sum rule;
           expansion of sums (e.g. 'sin(x + y)') will take place only if plus
           is True.

        -  ``times`` - (default: True) Controls the product
           rule, expansion of products (e.g. sin(2\*x)) will take place only
           if times is True.


        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: sin(5*x).expand_trig()
            5*cos(x)^4*sin(x) - 10*cos(x)^2*sin(x)^3 + sin(x)^5
            sage: cos(2*x + var('y')).expand_trig()
            cos(2*x)*cos(y) - sin(2*x)*sin(y)

        We illustrate various options to this function::

            sage: f = sin(sin(3*cos(2*x))*x)
            sage: f.expand_trig()
            sin((3*cos(cos(2*x))^2*sin(cos(2*x)) - sin(cos(2*x))^3)*x)
            sage: f.expand_trig(full=True)
            sin((3*(cos(cos(x)^2)*cos(sin(x)^2) + sin(cos(x)^2)*sin(sin(x)^2))^2*(cos(sin(x)^2)*sin(cos(x)^2) - cos(cos(x)^2)*sin(sin(x)^2)) - (cos(sin(x)^2)*sin(cos(x)^2) - cos(cos(x)^2)*sin(sin(x)^2))^3)*x)
            sage: sin(2*x).expand_trig(times=False)
            sin(2*x)
            sage: sin(2*x).expand_trig(times=True)
            2*cos(x)*sin(x)
            sage: sin(2 + x).expand_trig(plus=False)
            sin(x + 2)
            sage: sin(2 + x).expand_trig(plus=True)
            cos(x)*sin(2) + cos(2)*sin(x)
            sage: sin(x/2).expand_trig(half_angles=False)
            sin(1/2*x)
            sage: sin(x/2).expand_trig(half_angles=True)
            (-1)^floor(1/2*x/pi)*sqrt(-1/2*cos(x) + 1/2)

        ALIASES:

        :meth:`trig_expand` and :meth:`expand_trig` are the same
        """
        from sage.calculus.calculus import maxima_options
        M = self._maxima_()
        P = M.parent()
        opt = maxima_options(trigexpand=full, halfangles=half_angles,
                             trigexpandplus=plus, trigexpandtimes=times)
        cmd = 'trigexpand(%s), %s'%(M.name(), opt)
        ans = P(cmd)
        return self.parent()(ans)

    trig_expand = expand_trig

    def reduce_trig(self, var=None):
        r"""
        Combines products and powers of trigonometric and hyperbolic
        sin's and cos's of x into those of multiples of x. It also
        tries to eliminate these functions when they occur in
        denominators.

        INPUT:

        - ``self`` - a symbolic expression

        - ``var`` - (default: None) the variable which is used for
          these transformations. If not specified, all variables are
          used.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: y=var('y')
            sage: f=sin(x)*cos(x)^3+sin(y)^2
            sage: f.reduce_trig()
            -1/2*cos(2*y) + 1/8*sin(4*x) + 1/4*sin(2*x) + 1/2

        To reduce only the expressions involving x we use optional parameter::

            sage: f.reduce_trig(x)
            sin(y)^2 + 1/8*sin(4*x) + 1/4*sin(2*x)

        ALIASES: :meth:`trig_reduce` and :meth:`reduce_trig` are the same
        """
        M = self._maxima_()
        P = M.parent()
        if var is None:
            cmd = 'trigreduce(%s)'%(M.name())
        else:
            cmd = 'trigreduce(%s,%s)'%(M.name(),str(var))
        ans = P(cmd)
        return self.parent()(ans)

    trig_reduce = reduce_trig

    ############################################################################
    # Pattern Matching
    ############################################################################
    def match(self, pattern):
        """
        Check if self matches the given pattern.

        INPUT:

        -  ``pattern`` -- a symbolic expression, possibly containing wildcards
           to match for

        OUTPUT:

        One of

        ``None`` if there is no match, or a dictionary mapping the
        wildcards to the matching values if a match was found. Note
        that the dictionary is empty if there were no wildcards in the
        given pattern.

        See also http://www.ginac.de/tutorial/Pattern-matching-and-advanced-substitutions.html

        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,f,g')
            (x, y, z, a, b, c, d, f, g)
            sage: w0 = SR.wild(0); w1 = SR.wild(1); w2 = SR.wild(2)
            sage: ((x+y)^a).match((x+y)^a)  # no wildcards, so empty dict
            {}
            sage: print ((x+y)^a).match((x+y)^b)
            None
            sage: t = ((x+y)^a).match(w0^w1)
            sage: t[w0], t[w1]
            (x + y, a)
            sage: print ((x+y)^a).match(w0^w0)
            None
            sage: ((x+y)^(x+y)).match(w0^w0)
            {$0: x + y}
            sage: t = ((a+b)*(a+c)).match((a+w0)*(a+w1))
            sage: t[w0], t[w1]
            (c, b)
            sage: ((a+b)*(a+c)).match((w0+b)*(w0+c))
            {$0: a}
            sage: t = ((a+b)*(a+c)).match((w0+w1)*(w0+w2))
            sage: t[w0], t[w1], t[w2]
            (a, c, b)
            sage: print ((a+b)*(a+c)).match((w0+w1)*(w1+w2))
            None
            sage: t = (a*(x+y)+a*z+b).match(a*w0+w1)
            sage: t[w0], t[w1]
            (x + y, a*z + b)
            sage: print (a+b+c+d+f+g).match(c)
            None
            sage: (a+b+c+d+f+g).has(c)
            True
            sage: (a+b+c+d+f+g).match(c+w0)
            {$0: a + b + d + f + g}
            sage: (a+b+c+d+f+g).match(c+g+w0)
            {$0: a + b + d + f}
            sage: (a+b).match(a+b+w0)
            {$0: 0}
            sage: print (a*b^2).match(a^w0*b^w1)
            None
            sage: (a*b^2).match(a*b^w1)
            {$1: 2}
            sage: (x*x.arctan2(x^2)).match(w0*w0.arctan2(w0^2))
            {$0: x}

        Beware that behind-the-scenes simplification can lead to
        surprising results in matching::

            sage: print (x+x).match(w0+w1)
            None
            sage: t = x+x; t
            2*x
            sage: t.operator()
            <built-in function mul>

        Since asking to match w0+w1 looks for an addition operator,
        there is no match.
        """
        cdef Expression p = self.coerce_in(pattern)
        cdef GExList mlst
        cdef bint res = self._gobj.match(p._gobj, mlst)
        if not res:
            return None

        cdef dict rdict = {}
        cdef GExListIter itr = mlst.begin()
        cdef GExListIter lstend = mlst.end()
        while itr.is_not_equal(lstend):
            key = new_Expression_from_GEx(self._parent, itr.obj().lhs())
            val = new_Expression_from_GEx(self._parent, itr.obj().rhs())
            rdict[key] = val
            itr.inc()
        return rdict


    def find(self, pattern):
        """
        Find all occurrences of the given pattern in this expression.

        Note that once a subexpression matches the pattern, the search doesn't
        extend to subexpressions of it.

        EXAMPLES::

            sage: var('x,y,z,a,b')
            (x, y, z, a, b)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)

            sage: (sin(x)*sin(y)).find(sin(w0))
            [sin(y), sin(x)]

            sage: ((sin(x)+sin(y))*(a+b)).expand().find(sin(w0))
            [sin(y), sin(x)]

            sage: (1+x+x^2+x^3).find(x)
            [x]
            sage: (1+x+x^2+x^3).find(x^w0)
            [x^2, x^3]

            sage: (1+x+x^2+x^3).find(y)
            []

            # subexpressions of a match are not listed
            sage: ((x^y)^z).find(w0^w1)
            [(x^y)^z]
        """
        cdef Expression p = self.coerce_in(pattern)
        cdef GExList found
        self._gobj.find(p._gobj, found)
        res = []
        cdef GExListIter itr = found.begin()
        while itr.is_not_equal(found.end()):
            res.append(new_Expression_from_GEx(self._parent, itr.obj()))
            itr.inc()
        res.sort(cmp)
        return res

    def has(self, pattern):
        """
        EXAMPLES::

            sage: var('x,y,a'); w0 = SR.wild(); w1 = SR.wild()
            (x, y, a)
            sage: (x*sin(x + y + 2*a)).has(y)
            True

        Here "x+y" is not a subexpression of "x+y+2*a" (which has the
        subexpressions "x", "y" and "2*a")::

            sage: (x*sin(x + y + 2*a)).has(x+y)
            False
            sage: (x*sin(x + y + 2*a)).has(x + y + w0)
            True

        The following fails because "2*(x+y)" automatically gets converted to
        "2*x+2*y" of which "x+y" is not a subexpression::

            sage: (x*sin(2*(x+y) + 2*a)).has(x+y)
            False

        Although x^1==x and x^0==1, neither "x" nor "1" are actually of the
        form "x^something"::

            sage: (x+1).has(x^w0)
            False

        Here is another possible pitfall, where the first expression
        matches because the term "-x" has the form "(-1)*x" in GiNaC. To check
        whether a polynomial contains a linear term you should use the
        coeff() function instead.

        ::

            sage: (4*x^2 - x + 3).has(w0*x)
            True
            sage: (4*x^2 + x + 3).has(w0*x)
            False
            sage: (4*x^2 + x + 3).has(x)
            True
            sage: (4*x^2 - x + 3).coeff(x,1)
            -1
            sage: (4*x^2 + x + 3).coeff(x,1)
            1
        """
        cdef Expression p = self.coerce_in(pattern)
        return self._gobj.has(p._gobj)

    def substitute(self, in_dict=None, **kwds):
        """
        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,f,g')
            (x, y, z, a, b, c, d, f, g)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: t = a^2 + b^2 + (x+y)^3

            # substitute with keyword arguments (works only with symbols)
            sage: t.subs(a=c)
            (x + y)^3 + b^2 + c^2

            # substitute with a dictionary argument
            sage: t.subs({a^2: c})
            (x + y)^3 + b^2 + c

            sage: t.subs({w0^2: w0^3})
            a^3 + b^3 + (x + y)^3

            # substitute with a relational expression
            sage: t.subs(w0^2 == w0^3)
            a^3 + b^3 + (x + y)^3

            sage: t.subs(w0==w0^2)
            (x^2 + y^2)^18 + a^16 + b^16

            # more than one keyword argument is accepted
            sage: t.subs(a=b, b=c)
            (x + y)^3 + b^2 + c^2

            # using keyword arguments with a dictionary is allowed
            sage: t.subs({a:b}, b=c)
            (x + y)^3 + b^2 + c^2

            # in this case keyword arguments override the dictionary
            sage: t.subs({a:b}, a=c)
            (x + y)^3 + b^2 + c^2

            sage: t.subs({a:b, b:c})
            (x + y)^3 + b^2 + c^2

        TESTS::

            sage: # no arguments return the same expression
            sage: t.subs()
            (x + y)^3 + a^2 + b^2

            # similarly for an empty dictionary argument
            sage: t.subs({})
            (x + y)^3 + a^2 + b^2

            # non keyword or dictionary argument returns error
            sage: t.subs(5)
            Traceback (most recent call last):
            ...
            TypeError: subs takes either a set of keyword arguments, a dictionary, or a symbolic relational expression

            # substitutions with infinity
            sage: (x/y).subs(y=oo)
            0
            sage: (x/y).subs(x=oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity * f(x) encountered.
            sage: (x*y).subs(x=oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity * f(x) encountered.
            sage: (x^y).subs(x=oo)
            Traceback (most recent call last):
            ...
            ValueError: power::eval(): pow(Infinity, f(x)) is not defined.
            sage: (x^y).subs(y=oo)
            Traceback (most recent call last):
            ...
            ValueError: power::eval(): pow(f(x), infinity) is not defined.
            sage: (x+y).subs(x=oo)
            +Infinity
            sage: (x-y).subs(y=oo)
            -Infinity
            sage: gamma(x).subs(x=-1)
            Infinity
            sage: 1/gamma(x).subs(x=-1)
            0

            # verify that this operation does not modify the passed dictionary (#6622)
            sage: var('v t')
            (v, t)
            sage: f = v*t
            sage: D = {v: 2}
            sage: f(D, t=3)
            6
            sage: D
            {v: 2}

        Check if #9891 is fixed::

            sage: exp(x).subs(x=log(x))
            x

        Check if :trac:`13587` is fixed::

            sage: t = tan(x)^2 - tan(x)
            sage: t.subs(x=pi/2)
            Infinity
            sage: u = gamma(x) - gamma(x-1)
            sage: u.subs(x=-1)
            Infinity
        """
        cdef dict sdict = {}
        if in_dict is not None:
            if isinstance(in_dict, Expression):
                return self._subs_expr(in_dict)
            if not isinstance(in_dict, dict):
                raise TypeError, "subs takes either a set of keyword arguments, a dictionary, or a symbolic relational expression"
            sdict.update(in_dict)

        if kwds:
            for k, v in kwds.iteritems():
                k = self._parent.var(k)
                sdict[k] = v

        cdef GExMap smap
        for k, v in sdict.iteritems():
            smap.insert(make_pair((<Expression>self.coerce_in(k))._gobj,
                                  (<Expression>self.coerce_in(v))._gobj))

        return new_Expression_from_GEx(self._parent, self._gobj.subs_map(smap))

    subs = substitute

    cpdef Expression _subs_expr(self, expr):
        """
        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,f')
            (x, y, z, a, b, c, d, f)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: (a^2 + b^2 + (x+y)^2)._subs_expr(w0^2 == w0^3)
            a^3 + b^3 + (x + y)^3
            sage: (a^4 + b^4 + (x+y)^4)._subs_expr(w0^2 == w0^3)
            a^4 + b^4 + (x + y)^4
            sage: (a^2 + b^4 + (x+y)^4)._subs_expr(w0^2 == w0^3)
            b^4 + (x + y)^4 + a^3
            sage: ((a+b+c)^2)._subs_expr(a+b == x)
            (a + b + c)^2
            sage: ((a+b+c)^2)._subs_expr(a+b+w0 == x+w0)
            (c + x)^2
            sage: (a+2*b)._subs_expr(a+b == x)
            a + 2*b
            sage: (a+2*b)._subs_expr(a+b+w0 == x+w0)
            a + 2*b
            sage: (a+2*b)._subs_expr(a+w0*b == x)
            x
            sage: (a+2*b)._subs_expr(a+b+w0*b == x+w0*b)
            a + 2*b
            sage: (4*x^3-2*x^2+5*x-1)._subs_expr(x==a)
            4*a^3 - 2*a^2 + 5*a - 1
            sage: (4*x^3-2*x^2+5*x-1)._subs_expr(x^w0==a^w0)
            4*a^3 - 2*a^2 + 5*x - 1
            sage: (4*x^3-2*x^2+5*x-1)._subs_expr(x^w0==a^(2*w0))._subs_expr(x==a)
            4*a^6 - 2*a^4 + 5*a - 1
            sage: sin(1+sin(x))._subs_expr(sin(w0)==cos(w0))
            cos(cos(x) + 1)
            sage: (sin(x)^2 + cos(x)^2)._subs_expr(sin(w0)^2+cos(w0)^2==1)
            1
            sage: (1 + sin(x)^2 + cos(x)^2)._subs_expr(sin(w0)^2+cos(w0)^2==1)
            cos(x)^2 + sin(x)^2 + 1
            sage: (17*x + sin(x)^2 + cos(x)^2)._subs_expr(w1 + sin(w0)^2+cos(w0)^2 == w1 + 1)
            17*x + 1
            sage: ((x-1)*(sin(x)^2 + cos(x)^2)^2)._subs_expr(sin(w0)^2+cos(w0)^2 == 1)
            x - 1
            """
        cdef Expression p = self.coerce_in(expr)
        return new_Expression_from_GEx(self._parent, self._gobj.subs(p._gobj))

    def substitute_expression(self, *equations):
        """
        Given a dictionary of key:value pairs, substitute all occurrences
        of key for value in self.  The substitutions can also be given
        as a number of symbolic equalities key == value; see the
        examples.

        .. warning::

           This is a formal pattern substitution, which may or may not
           have any mathematical meaning. The exact rules used at
           present in Sage are determined by Maxima's subst
           command. Sometimes patterns are not replaced even though
           one would think they should be - see examples below.

        EXAMPLES::

            sage: f = x^2 + 1
            sage: f.subs_expr(x^2 == x)
            x + 1

        ::

            sage: var('x,y,z'); f = x^3 + y^2 + z
            (x, y, z)
            sage: f.subs_expr(x^3 == y^2, z == 1)
            2*y^2 + 1

        Or the same thing giving the substitutions as a dictionary::

            sage: f.subs_expr({x^3:y^2, z:1})
            2*y^2 + 1

            sage: f = x^2 + x^4
            sage: f.subs_expr(x^2 == x)
            x^4 + x
            sage: f = cos(x^2) + sin(x^2)
            sage: f.subs_expr(x^2 == x)
            cos(x) + sin(x)

        ::

            sage: f(x,y,t) = cos(x) + sin(y) + x^2 + y^2 + t
            sage: f.subs_expr(y^2 == t)
            (x, y, t) |--> x^2 + 2*t + cos(x) + sin(y)

        The following seems really weird, but it *is* what Maple does::

            sage: f.subs_expr(x^2 + y^2 == t)
            (x, y, t) |--> x^2 + y^2 + t + cos(x) + sin(y)
            sage: maple.eval('subs(x^2 + y^2 = t, cos(x) + sin(y) + x^2 + y^2 + t)')          # optional - maple
            'cos(x)+sin(y)+x^2+y^2+t'
            sage: maxima.quit()
            sage: maxima.eval('cos(x) + sin(y) + x^2 + y^2 + t, x^2 + y^2 = t')
            'sin(y)+y^2+cos(x)+x^2+t'

        Actually Mathematica does something that makes more sense::

            sage: mathematica.eval('Cos[x] + Sin[y] + x^2 + y^2 + t /. x^2 + y^2 -> t')       # optional - mathematica
            2 t + Cos[x] + Sin[y]
        """
        if isinstance(equations[0], dict):
            eq_dict = equations[0]
            equations = [ x == eq_dict[x] for x in eq_dict.keys() ]

        if not all([is_SymbolicEquation(eq) for eq in equations]):
            raise TypeError, "each expression must be an equation"

        d = dict([(eq.lhs(), eq.rhs()) for eq in equations])
        return self.subs(d)

    subs_expr = substitute_expression

    def substitute_function(self, original, new):
        """
        Returns this symbolic expressions all occurrences of the
        function *original* replaced with the function *new*.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: foo = function('foo'); bar = function('bar')
            sage: f = foo(x) + 1/foo(pi*y)
            sage: f.substitute_function(foo, bar)
            1/bar(pi*y) + bar(x)
        """
        from sage.symbolic.expression_conversions import SubstituteFunction
        return SubstituteFunction(self, original, new)()

    def __call__(self, *args, **kwds):
        """
        Calls the :meth:`subs` on this expression.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: (x+y)(x=z^2, y=x^y)
            z^2 + x^y
        """
        return self._parent._call_element_(self, *args, **kwds)

    def variables(self):
        """
        Return sorted tuple of variables that occur in this expression.

        EXAMPLES::

            sage: (x,y,z) = var('x,y,z')
            sage: (x+y).variables()
            (x, y)
            sage: (2*x).variables()
            (x,)
            sage: (x^y).variables()
            (x, y)
            sage: sin(x+y^z).variables()
            (x, y, z)

        """
        from sage.symbolic.ring import SR
        cdef GExSet sym_set
        g_list_symbols(self._gobj, sym_set)
        res = []
        cdef GExSetIter itr = sym_set.begin()
        while itr.is_not_equal(sym_set.end()):
            res.append(new_Expression_from_GEx(SR, itr.obj()))
            itr.inc()
        res.sort(cmp=lambda x,y: -cmp(x,y))
        return tuple(res)

    def arguments(self):
        """
        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f = x + y
            sage: f.arguments()
            (x, y)

            sage: g = f.function(x)
            sage: g.arguments()
            (x,)

        """
        try:
            return self._parent.arguments()
        except AttributeError:
            return self.variables()

    args = arguments

    def number_of_arguments(self):
        """
        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f = x + y
            sage: f.number_of_arguments()
            2

            sage: g = f.function(x)
            sage: g.number_of_arguments()
            1

        ::

            sage: x,y,z = var('x,y,z')
            sage: (x+y).number_of_arguments()
            2
            sage: (x+1).number_of_arguments()
            1
            sage: (sin(x)+1).number_of_arguments()
            1
            sage: (sin(z)+x+y).number_of_arguments()
            3
            sage: (sin(x+y)).number_of_arguments()
            2

        ::

            sage: ( 2^(8/9) - 2^(1/9) )(x-1)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0
        """
        return len(self.arguments())

    def number_of_operands(self):
        """
        Returns the number of arguments of this expression.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: a.number_of_operands()
            0
            sage: (a^2 + b^2 + (x+y)^2).number_of_operands()
            3
            sage: (a^2).number_of_operands()
            2
            sage: (a*b^2*c).number_of_operands()
            3
        """
        return self._gobj.nops()

    nops = number_of_operands

    def __len__(self):
        """
        Returns the number of arguments of this expression.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: len(a)
            0
            sage: len((a^2 + b^2 + (x+y)^2))
            3
            sage: len((a^2))
            2
            sage: len(a*b^2*c)
            3
        """
        return self.number_of_operands()

    def operands(self):
        """
        Returns a list containing the operands of this expression.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: (a^2 + b^2 + (x+y)^2).operands()
            [a^2, b^2, (x + y)^2]
            sage: (a^2).operands()
            [a, 2]
            sage: (a*b^2*c).operands()
            [a, b^2, c]
        """
        from sage.symbolic.ring import SR
        return [new_Expression_from_GEx(SR, self._gobj.op(i)) \
                            for i from 0 <= i < self._gobj.nops()]

    def operator(self):
        """
        Returns the topmost operator in this expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: (x+y).operator()
            <built-in function add>
            sage: (x^y).operator()
            <built-in function pow>
            sage: (x^y * z).operator()
            <built-in function mul>
            sage: (x < y).operator()
            <built-in function lt>

            sage: abs(x).operator()
            abs
            sage: r = gamma(x).operator(); type(r)
            <class 'sage.functions.other.Function_gamma'>

            sage: psi = function('psi', nargs=1)
            sage: psi(x).operator()
            psi

            sage: r = psi(x).operator()
            sage: r == psi
            True

            sage: f = function('f', nargs=1, conjugate_func=lambda self, x: 2*x)
            sage: nf = f(x).operator()
            sage: nf(x).conjugate()
            2*x

            sage: f = function('f')
            sage: a = f(x).diff(x); a
            D[0](f)(x)
            sage: a.operator()
            D[0](f)

        TESTS::

            sage: (x <= y).operator()
            <built-in function le>
            sage: (x == y).operator()
            <built-in function eq>
            sage: (x != y).operator()
            <built-in function ne>
            sage: (x > y).operator()
            <built-in function gt>
            sage: (x >= y).operator()
            <built-in function ge>
        """
        cdef operators o
        cdef unsigned serial
        import operator
        if is_a_add(self._gobj):
            return operator.add
        elif is_a_mul(self._gobj) or is_a_ncmul(self._gobj):
            return operator.mul
        elif is_a_power(self._gobj):
            return operator.pow
        elif is_a_relational(self._gobj):
            # find the operator and return it
            o = relational_operator(self._gobj)
            if o == equal:
                return operator.eq
            elif o == not_equal:
                return operator.ne
            elif o == less:
                return operator.lt
            elif o == less_or_equal:
                return operator.le
            elif o == greater:
                return operator.gt
            elif o == greater_or_equal:
                return operator.ge
            else:
                raise RuntimeError, "operator type not known, please report this as a bug"
        elif is_a_function(self._gobj):
            # get function id
            serial = ex_to_function(self._gobj).get_serial()

            # if operator is a special function defined by us
            # find the python equivalent and return it
            res = get_sfunction_from_serial(serial)
            if res is None:
                raise RuntimeError, "cannot find SFunction in table"

            if is_a_fderivative(self._gobj):
                from sage.symbolic.pynac import paramset_from_Expression
                from sage.symbolic.operators import FDerivativeOperator
                parameter_set = paramset_from_Expression(self)
                res = FDerivativeOperator(res, parameter_set)

            return res

        # self._gobj is either a symbol, constant or numeric
        return None

    def __index__(self):
        """
        EXAMPLES::

            sage: a = range(10)
            sage: a[:SR(5)]
            [0, 1, 2, 3, 4]
        """
        return int(self._integer_())

    def iterator(self):
        """
        Return an iterator over the operands of this expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: list((x+y+z).iterator())
            [x, y, z]
            sage: list((x*y*z).iterator())
            [x, y, z]
            sage: list((x^y*z*(x+y)).iterator())
            [x + y, x^y, z]

        Note that symbols, constants and numeric objects don't have operands,
        so the iterator function raises an error in these cases::

            sage: x.iterator()
            Traceback (most recent call last):
            ...
            ValueError: expressions containing only a numeric coefficient, constant or symbol have no operands
            sage: pi.iterator()
            Traceback (most recent call last):
            ...
            ValueError: expressions containing only a numeric coefficient, constant or symbol have no operands
            sage: SR(5).iterator()
            Traceback (most recent call last):
            ...
            ValueError: expressions containing only a numeric coefficient, constant or symbol have no operands
        """
        if is_a_symbol(self._gobj) or is_a_constant(self._gobj) or \
                is_a_numeric(self._gobj):
                    raise ValueError, "expressions containing only a numeric coefficient, constant or symbol have no operands"
        return new_ExpIter_from_Expression(self)

    property op:
        def __get__(self):
            """
            Provide access to the operands of an expression through a property.


            EXAMPLES::

                sage: t = 1+x+x^2
                sage: t.op
                Operands of x^2 + x + 1
                sage: x.op
                Traceback (most recent call last):
                ...
                TypeError: expressions containing only a numeric coefficient, constant or symbol have no operands
                sage: t.op[0]
                x^2

            Indexing directly with ``t[1]`` causes problems with numpy types.

                sage: t[1]
                Traceback (most recent call last):
                ...
                TypeError: 'sage.symbolic.expression.Expression' object does not support indexing
            """
            if is_a_symbol(self._gobj) or is_a_constant(self._gobj) or \
                is_a_numeric(self._gobj):
                    raise TypeError, "expressions containing only a numeric coefficient, constant or symbol have no operands"
            cdef OperandsWrapper res = OperandsWrapper.__new__(OperandsWrapper)
            res._expr = self
            return res

    def _numerical_approx(self, prec=None, digits=None):
        """
        Return a numerical approximation this symbolic expression as
        either a real or complex number with at least the requested
        number of bits or digits of precision.

        EXAMPLES::

            sage: sin(x).subs(x=5).n()
            -0.958924274663138
            sage: sin(x).subs(x=5).n(100)
            -0.95892427466313846889315440616
            sage: sin(x).subs(x=5).n(digits=50)
            -0.95892427466313846889315440615599397335246154396460
            sage: zeta(x).subs(x=2).numerical_approx(digits=50)
            1.6449340668482264364724151666460251892189499012068

            sage: cos(3).numerical_approx(200)
            -0.98999249660044545727157279473126130239367909661558832881409
            sage: numerical_approx(cos(3),200)
            -0.98999249660044545727157279473126130239367909661558832881409
            sage: numerical_approx(cos(3), digits=10)
            -0.9899924966
            sage: (i + 1).numerical_approx(32)
            1.00000000 + 1.00000000*I
            sage: (pi + e + sqrt(2)).numerical_approx(100)
            7.2740880444219335226246195788

        TESTS:

        We test the evaluation of different infinities available in Pynac::

            sage: t = x - oo; t
            -Infinity
            sage: t.n()
            -infinity
            sage: t = x + oo; t
            +Infinity
            sage: t.n()
            +infinity
            sage: t = x - unsigned_infinity; t
            Infinity
            sage: t.n()
            +infinity

        Some expressions can't be evaluated numerically::

            sage: n(sin(x))
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression numerically
            sage: a = var('a')
            sage: (x^2 + 2*x + 2).subs(x=a).n()
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression numerically

        Make sure we've rounded up log(10,2) enough to guarantee
        sufficient precision (trac #10164)::

            sage: ks = 4*10**5, 10**6
            sage: all(len(str(e.n(digits=k)))-1 >= k for k in ks)
            True

        """
        if prec is None:
            if digits is None:
                prec = 53
            else:
                prec = int((digits+1) * LOG_TEN_TWO_PLUS_EPSILON) + 1
        from sage.rings.real_mpfr import RealField
        R = RealField(prec)
        cdef Expression x
        try:
            x = self._convert(R)
        except TypeError: # numerical approximation for real number failed
            pass          # try again with complex
            R = R.complex_field()
            x = self._convert(R)

        # we have to consider constants as well, since infinity is a constant
        # in pynac
        if is_a_numeric(x._gobj):
            res = py_object_from_numeric(x._gobj)
        elif  is_a_constant(x._gobj):
            res = x.pyobject()
        else:
            raise TypeError("cannot evaluate symbolic expression numerically")

        # Important -- the  we get might not be a valid output for numerical_approx in
        # the case when one gets infinity.
        if isinstance(res, AnInfinity):
            return res.n(prec=prec,digits=digits)
        return res

    #added this line to make doctests visible to users
    numerical_approx =_numerical_approx
    n=_numerical_approx
    N=_numerical_approx

    def round(self):
        """
        Round this expression to the nearest integer.

        EXAMPLES::

            sage: u = sqrt(43203735824841025516773866131535024)
            sage: u.round()
            207855083711803945
            sage: t = sqrt(Integer('1'*1000)).round(); print str(t)[-10:]
            3333333333
            sage: (-sqrt(110)).round()
            -10
            sage: (-sqrt(115)).round()
            -11
            sage: (sqrt(-3)).round()
            Traceback (most recent call last):
            ...
            ValueError: could not convert sqrt(-3) to a real number
        """
        try:
            return self.pyobject().round()
        except (TypeError, AttributeError):
            pass
        from sage.functions.all import floor, ceil
        try:
            rif_self = sage.rings.all.RIF(self)
        except TypeError:
            raise ValueError, "could not convert %s to a real number"%(self)
        half = 1 / sage.rings.integer.Integer(2)
        if rif_self < 0 or (rif_self.contains_zero() and self < 0):
            result = ceil(self - half)
        else:
            result = floor(self + half)
        if not isinstance(result, sage.rings.integer.Integer):
            raise ValueError, "could not convert %s to a real number"%(self)
        else:
            return result

    def function(self, *args):
        """
        Return a callable symbolic expression with the given variables.

        EXAMPLES:

        We will use several symbolic variables in the examples below::

            sage: var('x, y, z, t, a, w, n')
            (x, y, z, t, a, w, n)

        ::

            sage: u = sin(x) + x*cos(y)
            sage: g = u.function(x,y)
            sage: g(x,y)
            x*cos(y) + sin(x)
            sage: g(t,z)
            t*cos(z) + sin(t)
            sage: g(x^2, x^y)
            x^2*cos(x^y) + sin(x^2)

        ::

            sage: f = (x^2 + sin(a*w)).function(a,x,w); f
            (a, x, w) |--> x^2 + sin(a*w)
            sage: f(1,2,3)
            sin(3) + 4

        Using the :meth:`function` method we can obtain the above function
        `f`, but viewed as a function of different variables::

            sage: h = f.function(w,a); h
            (w, a) |--> x^2 + sin(a*w)

        This notation also works::

            sage: h(w,a) = f
            sage: h
            (w, a) |--> x^2 + sin(a*w)

        You can even make a symbolic expression `f` into a function
        by writing ``f(x,y) = f``::

            sage: f = x^n + y^n; f
            x^n + y^n
            sage: f(x,y) = f
            sage: f
            (x, y) |--> x^n + y^n
            sage: f(2,3)
            3^n + 2^n
        """
        # we override type checking in CallableSymbolicExpressionRing,
        # since it checks for old SymbolicVariable's
        # and do the check here instead
        from sage.symbolic.callable import CallableSymbolicExpressionRing
        from sage.symbolic.ring import is_SymbolicVariable
        for i in args:
            if not is_SymbolicVariable(i):
                break
        else:
            R = CallableSymbolicExpressionRing(args, check=False)
            return R(self)
        raise TypeError, "Must construct a function with a tuple (or list) of symbolic variables."

    ############################################################################
    # Basic arithmetic wrappers
    # which allow disabling automatic evaluation with the hold parameter
    ############################################################################
    def power(self, exp, hold=False):
        """
        Returns the current expression to the power ``exp``.

        To prevent automatic evaluation use the ``hold`` argument.

        EXAMPLES::

            sage: (x^2).power(2)
            x^4
            sage: (x^2).power(2, hold=True)
            (x^2)^2

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = (x^2).power(2, hold=True); a.simplify()
            x^4

        """
        cdef Expression nexp = self.coerce_in(exp)
        return new_Expression_from_GEx(self._parent,
                g_hold2_wrapper(g_power_construct, self._gobj, nexp._gobj,
                    hold))

    def add(self, *args, hold=False):
        """
        Return the sum of the current expression and the given arguments.

        To prevent automatic evaluation use the ``hold`` argument.

        EXAMPLES::

            sage: x.add(x)
            2*x
            sage: x.add(x, hold=True)
            x + x
            sage: x.add(x, (2+x), hold=True)
            (x + 2) + x + x
            sage: x.add(x, (2+x), x, hold=True)
            (x + 2) + x + x + x
            sage: x.add(x, (2+x), x, 2*x, hold=True)
            (x + 2) + 2*x + x + x + x

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = x.add(x, hold=True); a.simplify()
            2*x
        """
        nargs = [self.coerce_in(x) for x in args]
        cdef GExVector vec
        cdef Py_ssize_t i
        vec.push_back(self._gobj)
        for i in range(len(args)):
            vec.push_back((<Expression>nargs[i])._gobj)
        return new_Expression_from_GEx(self._parent, g_add_construct(vec, hold))

    def mul(self, *args, hold=False):
        """
        Return the product of the current expression and the given arguments.

        To prevent automatic evaluation use the ``hold`` argument.

        EXAMPLES::

            sage: x.mul(x)
            x^2
            sage: x.mul(x, hold=True)
            x*x
            sage: x.mul(x, (2+x), hold=True)
            (x + 2)*x*x
            sage: x.mul(x, (2+x), x, hold=True)
            (x + 2)*x*x*x
            sage: x.mul(x, (2+x), x, 2*x, hold=True)
            (2*x)*(x + 2)*x*x*x

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = x.mul(x, hold=True); a.simplify()
            x^2

        """
        nargs = [self.coerce_in(x) for x in args]
        cdef GExVector vec
        cdef Py_ssize_t i
        vec.push_back(self._gobj)
        for i in range(len(args)):
            vec.push_back((<Expression>nargs[i])._gobj)
        return new_Expression_from_GEx(self._parent, g_mul_construct(vec, hold))

    ############################################################################
    # Polynomial functions
    ############################################################################
    def coefficient(self, s, int n=1):
        """
        Returns the coefficient of `s^n` in this symbolic expression.

        INPUT:

        - ``s`` - expression

        - ``n`` - integer, default 1

        OUTPUT:

        A symbolic expression. The coefficient of `s^n`.

        Sometimes it may be necessary to expand or factor first, since this
        is not done automatically.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.collect(x)
            x^3*sin(x*y) + (a + y + 1/y)*x + 2*sin(x*y)/x + 100
            sage: f.coefficient(x,0)
            100
            sage: f.coefficient(x,-1)
            2*sin(x*y)
            sage: f.coefficient(x,1)
            a + y + 1/y
            sage: f.coefficient(x,2)
            0
            sage: f.coefficient(x,3)
            sin(x*y)
            sage: f.coefficient(x^3)
            sin(x*y)
            sage: f.coefficient(sin(x*y))
            x^3 + 2/x
            sage: f.collect(sin(x*y))
            a*x + x*y + (x^3 + 2/x)*sin(x*y) + x/y + 100

            sage: var('a, x, y, z')
            (a, x, y, z)
            sage: f = (a*sqrt(2))*x^2 + sin(y)*x^(1/2) + z^z
            sage: f.coefficient(sin(y))
            sqrt(x)
            sage: f.coefficient(x^2)
            sqrt(2)*a
            sage: f.coefficient(x^(1/2))
            sin(y)
            sage: f.coefficient(1)
            0
            sage: f.coefficient(x, 0)
            sqrt(x)*sin(y) + z^z

        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.coeff(ss._gobj, n))

    coeff = coefficient

    def coefficients(self, x=None):
        r"""
        Coefficients of this symbolic expression as a polynomial in x.

        INPUT:

        -  ``x`` -- optional variable.

        OUTPUT:

        A list of pairs ``(expr, n)``, where ``expr`` is a symbolic
        expression and ``n`` is a power.

        EXAMPLES::

            sage: var('x, y, a')
            (x, y, a)
            sage: p = x^3 - (x-3)*(x^2+x) + 1
            sage: p.coefficients()
            [[1, 0], [3, 1], [2, 2]]
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: p.coefficients(a)
            [[x^2 + x + 1, 0], [-2*sqrt(2)*x, 1], [2, 2]]
            sage: p.coefficients(x)
            [[2*a^2 + 1, 0], [-2*sqrt(2)*a + 1, 1], [1, 2]]

        A polynomial with wacky exponents::

            sage: p = (17/3*a)*x^(3/2) + x*y + 1/x + x^x
            sage: p.coefficients(x)
            [[1, -1], [x^x, 0], [y, 1], [17/3*a, 3/2]]
        """
        f = self._maxima_()
        maxima = f.parent()
        maxima._eval_line('load(coeflist)')
        if x is None:
            x = self.default_variable()
        x = self.parent().var(repr(x))
        G = f.coeffs(x)
        from sage.calculus.calculus import symbolic_expression_from_maxima_string
        S = symbolic_expression_from_maxima_string(repr(G))
        return S[1:]

    coeffs = coefficients

    def leading_coefficient(self, s):
        """
        Return the leading coefficient of s in self.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.leading_coefficient(x)
            sin(x*y)
            sage: f.leading_coefficient(y)
            x
            sage: f.leading_coefficient(sin(x*y))
            x^3 + 2/x
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.lcoeff(ss._gobj))

    leading_coeff = leading_coefficient

    def trailing_coefficient(self, s):
        """
        Return the trailing coefficient of s in self, i.e., the coefficient
        of the smallest power of s in self.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.trailing_coefficient(x)
            2*sin(x*y)
            sage: f.trailing_coefficient(y)
            x
            sage: f.trailing_coefficient(sin(x*y))
            a*x + x*y + x/y + 100
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.tcoeff(ss._gobj))

    trailing_coeff = trailing_coefficient

    def low_degree(self, s):
        """
        Return the exponent of the lowest nonpositive power of s in self.

        OUTPUT:

        An integer ``<= 0``.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y^10 + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + 2*sin(x*y)/x + x/y^10 + 100
            sage: f.low_degree(x)
            -1
            sage: f.low_degree(y)
            -10
            sage: f.low_degree(sin(x*y))
            0
            sage: (x^3+y).low_degree(x)
            0
        """
        cdef Expression ss = self.coerce_in(s)
        return self._gobj.ldegree(ss._gobj)

    def degree(self, s):
        """
        Return the exponent of the highest nonnegative power of s in self.

        OUTPUT:

        An integer ``>= 0``.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y^10 + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + 2*sin(x*y)/x + x/y^10 + 100
            sage: f.degree(x)
            3
            sage: f.degree(y)
            1
            sage: f.degree(sin(x*y))
            1
            sage: (x^-3+y).degree(x)
            0
        """
        cdef Expression ss = self.coerce_in(s)
        return self._gobj.degree(ss._gobj)

    def unit(self, s):
        """
        Return the unit of this expression when considered as a
        polynomial in ``s``.

        See also :meth:`content`, :meth:`primitive_part`, and
        :meth:`unit_content_primitive`.

        INPUT:

        - ``s`` -- a symbolic expression.

        OUTPUT:

        The unit part of a polynomial as a symbolic expression. It is
        defined as the sign of the leading coefficient.

        EXAMPLES::

            sage: (2*x+4).unit(x)
            1
            sage: (-2*x+1).unit(x)
            -1
            sage: (2*x+1/2).unit(x)
            1
            sage: var('y')
            y
            sage: (2*x - 4*sin(y)).unit(sin(y))
            -1
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.unit(ss._gobj))

    def content(self, s):
        """
        Return the content of this expression when considered as a
        polynomial in ``s``.

        See also :meth:`unit`, :meth:`primitive_part`, and
        :meth:`unit_content_primitive`.

        INPUT:

        - ``s`` -- a symbolic expression.

        OUTPUT:

        The content part of a polynomial as a symbolic expression. It
        is defined as the gcd of the coefficients.

        .. warning::

            The expression is considered to be a univariate polynomial
            in ``s``. The output is different from the ``content()``
            method provided by multivariate polynomial rings in Sage.

        EXAMPLES::

            sage: (2*x+4).content(x)
            2
            sage: (2*x+1).content(x)
            1
            sage: (2*x+1/2).content(x)
            1/2
            sage: var('y')
            y
            sage: (2*x + 4*sin(y)).content(sin(y))
            2
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.content(ss._gobj))

    def primitive_part(self, s):
        """
        Return the primitive polynomial of this expression when
        considered as a polynomial in ``s``.

        See also :meth:`unit`, :meth:`content`, and
        :meth:`unit_content_primitive`.

        INPUT:

        - ``s`` -- a symbolic expression.

        OUTPUT:

        The primitive polynomial as a symbolic expression. It is
        defined as the quotient by the :meth:`unit` and
        :meth:`content` parts (with respect to the variable ``s``).

        EXAMPLES::

            sage: (2*x+4).primitive_part(x)
            x + 2
            sage: (2*x+1).primitive_part(x)
            2*x + 1
            sage: (2*x+1/2).primitive_part(x)
            4*x + 1
            sage: var('y')
            y
            sage: (2*x + 4*sin(y)).primitive_part(sin(y))
            x + 2*sin(y)
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.primpart(ss._gobj))

    def unit_content_primitive(self, s):
        """
        Return the factorization into unit, content, and primitive part.

        INPUT:

        - ``s`` -- a symbolic expression, usually a symbolic
          variable. The whole symbolic expression ``self`` will be
          considered as a univariate polynomial in ``s``.

        OUTPUT:

        A triple (unit, content, primitive polynomial)` containing the
        :meth:`unit <unit>`, :meth:`content <content>`, and
        :meth:`primitive polynomial <primitive_part>`. Their product equals
        ``self``.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: ex = 9*x^3*y+3*y
            sage: ex.unit_content_primitive(x)
            (1, 3*y, 3*x^3 + 1)
            sage: ex.unit_content_primitive(y)
            (1, 9*x^3 + 3, y)
        """
        cdef Expression ss = self.coerce_in(s)
        cdef GEx unit, cont, prim
        self._gobj.unitcontprim(ss._gobj, unit, cont, prim)
        return (new_Expression_from_GEx(self._parent, unit),
                new_Expression_from_GEx(self._parent, cont),
                new_Expression_from_GEx(self._parent, prim))

    def poly(self, x=None):
        r"""
        Express this symbolic expression as a polynomial in *x*. If
        this is not a polynomial in *x*, then some coefficients may be
        functions of *x*.

        .. warning::

           This is different from :meth:`polynomial` which returns
           a Sage polynomial over a given base ring.

        EXAMPLES::

            sage: var('a, x')
            (a, x)
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: p.poly(a)
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: bool(p.poly(a) == (x-a*sqrt(2))^2 + x + 1)
            True
            sage: p.poly(x)
            2*a^2 - (2*sqrt(2)*a - 1)*x + x^2 + 1
        """
        from sage.symbolic.ring import SR
        f = self._maxima_()
        P = f.parent()
        P._eval_line('load(coeflist)')
        if x is None:
            x = self.default_variable()
        x = self._parent.var(repr(x))
        G = f.coeffs(x)
        ans = None
        for i in range(1, len(G)):
            Z = G[i]
            coeff = SR(Z[0])
            n = SR(Z[1])
            if repr(coeff) != '0':
                if repr(n) == '0':
                    xpow = SR(1)
                elif repr(n) == '1':
                    xpow = x
                else:
                    xpow = x**n
                if ans is None:
                    ans = coeff*xpow
                else:
                    ans += coeff*xpow
        return ans

    def polynomial(self, base_ring, ring=None):
        r"""
        Return this symbolic expression as an algebraic polynomial
        over the given base ring, if possible.

        The point of this function is that it converts purely symbolic
        polynomials into optimised algebraic polynomials over a given
        base ring.

        .. warning::

           This is different from :meth:`poly` which is used to rewrite
           self as a polynomial in terms of one of the variables.

        INPUT:

        -  ``base_ring`` - a ring

        EXAMPLES::

            sage: f = x^2 -2/3*x + 1
            sage: f.polynomial(QQ)
            x^2 - 2/3*x + 1
            sage: f.polynomial(GF(19))
            x^2 + 12*x + 1

        Polynomials can be useful for getting the coefficients of an
        expression::

            sage: g = 6*x^2 - 5
            sage: g.coefficients()
            [[-5, 0], [6, 2]]
            sage: g.polynomial(QQ).list()
            [-5, 0, 6]
            sage: g.polynomial(QQ).dict()
            {0: -5, 2: 6}

        ::

            sage: f = x^2*e + x + pi/e
            sage: f.polynomial(RDF)
            2.71828182846*x^2 + x + 1.15572734979
            sage: g = f.polynomial(RR); g
            2.71828182845905*x^2 + x + 1.15572734979092
            sage: g.parent()
            Univariate Polynomial Ring in x over Real Field with 53 bits of precision
            sage: f.polynomial(RealField(100))
            2.7182818284590452353602874714*x^2 + x + 1.1557273497909217179100931833
            sage: f.polynomial(CDF)
            2.71828182846*x^2 + x + 1.15572734979
            sage: f.polynomial(CC)
            2.71828182845905*x^2 + x + 1.15572734979092

        We coerce a multivariate polynomial with complex symbolic
        coefficients::

            sage: x, y, n = var('x, y, n')
            sage: f = pi^3*x - y^2*e - I; f
            pi^3*x - y^2*e - I
            sage: f.polynomial(CDF)
            (-2.71828182846)*y^2 + 31.0062766803*x - 1.0*I
            sage: f.polynomial(CC)
            (-2.71828182845905)*y^2 + 31.0062766802998*x - 1.00000000000000*I
            sage: f.polynomial(ComplexField(70))
            (-2.7182818284590452354)*y^2 + 31.006276680299820175*x - 1.0000000000000000000*I

        Another polynomial::

            sage: f = sum((e*I)^n*x^n for n in range(5)); f
            x^4*e^4 - I*x^3*e^3 - x^2*e^2 + I*x*e + 1
            sage: f.polynomial(CDF)
            54.5981500331*x^4 - 20.0855369232*I*x^3 - 7.38905609893*x^2 + 2.71828182846*I*x + 1.0
            sage: f.polynomial(CC)
            54.5981500331442*x^4 - 20.0855369231877*I*x^3 - 7.38905609893065*x^2 + 2.71828182845905*I*x + 1.00000000000000

        A multivariate polynomial over a finite field::

            sage: f = (3*x^5 - 5*y^5)^7; f
            (3*x^5 - 5*y^5)^7
            sage: g = f.polynomial(GF(7)); g
            3*x^35 + 2*y^35
            sage: parent(g)
            Multivariate Polynomial Ring in x, y over Finite Field of size 7
        """
        from sage.symbolic.expression_conversions import polynomial
        return polynomial(self, base_ring=base_ring, ring=ring)

    def _polynomial_(self, R):
        """
        Coerce this symbolic expression to a polynomial in `R`.

        EXAMPLES::

            sage: var('x,y,z,w')
            (x, y, z, w)

        ::

            sage: R = QQ[x,y,z]
            sage: R(x^2 + y)
            x^2 + y
            sage: R = QQ[w]
            sage: R(w^3 + w + 1)
            w^3 + w + 1
            sage: R = GF(7)[z]
            sage: R(z^3 + 10*z)
            z^3 + 3*z

        .. note::

           If the base ring of the polynomial ring is the symbolic ring,
           then a constant polynomial is always returned.

        ::

            sage: R = SR[x]
            sage: a = R(sqrt(2) + x^3 + y)
            sage: a
            x^3 + y + sqrt(2)
            sage: type(a)
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field'>
            sage: a.degree()
            0

        We coerce to a double precision complex polynomial ring::

            sage: f = e*x^3 + pi*y^3 + sqrt(2) + I; f
            pi*y^3 + x^3*e + sqrt(2) + I
            sage: R = CDF[x,y]
            sage: R(f)
            2.71828182846*x^3 + 3.14159265359*y^3 + 1.41421356237 + 1.0*I

        We coerce to a higher-precision polynomial ring

        ::

            sage: R = ComplexField(100)[x,y]
            sage: R(f)
            2.7182818284590452353602874714*x^3 + 3.1415926535897932384626433833*y^3 + 1.4142135623730950488016887242 + 1.0000000000000000000000000000*I

        TESTS:

        This shows that the issue at trac #5755 is fixed (attempting to
        coerce a symbolic expression to a non-symbolic polynomial ring
        caused an error::

            sage: xx = var('xx')
            sage: RDF['xx'](1.0*xx)
            xx
            sage: RDF['xx'](2.0*xx)
            2.0*xx
            sage: RR['xx'](1.0*xx)
            xx
            sage: RR['xx'](2.0*xx)
            2.00000000000000*xx

        This shows that the issue at trac #4246 is fixed (attempting to
        coerce an expression containing at least one variable that's not in
        `R` raises an error)::

            sage: x, y = var('x y')
            sage: S = PolynomialRing(Integers(4), 1, 'x')
            sage: S(x)
            x
            sage: S(y)
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Multivariate Polynomial Ring in x over Ring of integers modulo 4
            sage: S(x+y)
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Multivariate Polynomial Ring in x over Ring of integers modulo 4
            sage: (x+y)._polynomial_(S)
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Multivariate Polynomial Ring in x over Ring of integers modulo 4
        """
        from sage.symbolic.all import SR
        from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
        base_ring = R.base_ring()
        if base_ring == SR:
            if is_MPolynomialRing(R):
                return R({tuple([0]*R.ngens()):self})
            else:
                return R([self])
        return self.polynomial(None, ring=R)

    def power_series(self, base_ring):
        """
        Return algebraic power series associated to this symbolic
        expression, which must be a polynomial in one variable, with
        coefficients coercible to the base ring.

        The power series is truncated one more than the degree.

        EXAMPLES::

            sage: theta = var('theta')
            sage: f = theta^3 + (1/3)*theta - 17/3
            sage: g = f.power_series(QQ); g
            -17/3 + 1/3*theta + theta^3 + O(theta^4)
            sage: g^3
            -4913/27 + 289/9*theta - 17/9*theta^2 + 2602/27*theta^3 + O(theta^4)
            sage: g.parent()
            Power Series Ring in theta over Rational Field
        """
        v = self.variables()
        if len(v) != 1:
            raise ValueError, "self must be a polynomial in one variable but it is in the variables %s"%tuple([v])
        f = self.polynomial(base_ring)
        from sage.rings.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, names=f.parent().variable_names())
        return R(f, f.degree()+1)

    def gcd(self, b):
        """
        Return the gcd of self and b, which must be integers or polynomials over
        the rational numbers.

        TODO: I tried the massive gcd from
        http://trac.sagemath.org/sage_trac/ticket/694 on Ginac dies
        after about 10 seconds.  Singular easily does that GCD now.
        Since Ginac only handles poly gcd over QQ, we should change
        ginac itself to use Singular.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: SR(10).gcd(SR(15))
            5
            sage: (x^3 - 1).gcd(x-1)
            x - 1
            sage: (x^3 - 1).gcd(x^2+x+1)
            x^2 + x + 1
            sage: (x^3 - sage.symbolic.constants.pi).gcd(x-sage.symbolic.constants.pi)
            Traceback (most recent call last):
            ...
            ValueError: gcd: arguments must be polynomials over the rationals
            sage: gcd(x^3 - y^3, x-y)
            -x + y
            sage: gcd(x^100-y^100, x^10-y^10)
            -x^10 + y^10
            sage: gcd(expand( (x^2+17*x+3/7*y)*(x^5 - 17*y + 2/3) ), expand((x^13+17*x+3/7*y)*(x^5 - 17*y + 2/3)) )
            1/7*x^5 - 17/7*y + 2/21
        """
        cdef Expression r = self.coerce_in(b)
        cdef GEx x
        sig_on()
        try:
            x = g_gcd(self._gobj, r._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def lcm(self, b):
        """
        Return the lcm of self and b, which must be integers or
        polynomials over the rational numbers.  This is computed from
        the gcd of self and b implicitly from the relation
        self * b = gcd(self, b) * lcm(self, b).

        .. NOTE::

            In agreement with the convention in use for integers, if
            self * b == 0, then gcd(self, b) == max(self, b) and
            lcm(self, b) == 0.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: SR(10).lcm(SR(15))
            30
            sage: (x^3 - 1).lcm(x-1)
            x^3 - 1
            sage: (x^3 - 1).lcm(x^2+x+1)
            x^3 - 1
            sage: (x^3 - sage.symbolic.constants.pi).lcm(x-sage.symbolic.constants.pi)
            Traceback (most recent call last):
            ...
            ValueError: lcm: arguments must be polynomials over the rationals
            sage: lcm(x^3 - y^3, x-y)
            -x^3 + y^3
            sage: lcm(x^100-y^100, x^10-y^10)
            -x^100 + y^100
            sage: lcm(expand( (x^2+17*x+3/7*y)*(x^5 - 17*y + 2/3) ), expand((x^13+17*x+3/7*y)*(x^5 - 17*y + 2/3)) )
             1/21*(21*x^18 - 357*x^13*y + 14*x^13 + 357*x^6 + 9*x^5*y -
                     6069*x*y - 153*y^2 + 238*x + 6*y)*(21*x^7 + 357*x^6 +
                             9*x^5*y - 357*x^2*y + 14*x^2 - 6069*x*y -
                             153*y^2 + 238*x + 6*y)/(3*x^5 - 51*y + 2)

        TESTS:

        Verify that x * y = gcd(x,y) * lcm(x,y)::

            sage: x, y = var('x,y')
            sage: LRs = [(SR(10), SR(15)), (x^3-1, x-1), (x^3-y^3, x-y), (x^3-1, x^2+x+1), (SR(0), x-y)]
            sage: all((L.gcd(R) * L.lcm(R)) == L*R for L, R in LRs)
            True

        Make sure that the convention for what to do with the 0 is being respected::

            sage: gcd(x, SR(0)), lcm(x, SR(0))
            (x, 0)
            sage: gcd(SR(0), SR(0)), lcm(SR(0), SR(0))
            (0, 0)

        """
        sb = self * b
        try:
            return 0 if sb == 0 else sb / self.gcd(b)
        except ValueError:
            # make the error message refer to lcm, not gcd
            raise ValueError("lcm: arguments must be polynomials over the rationals")

    def collect(Expression self, s):
        """
        Collect the coefficients of ``s`` into a group.

        INPUT:

        - ``s`` -- the symbol whose coefficients will be collected.

        OUTPUT:

        A new expression, equivalent to the original one, with the
        coefficients of ``s`` grouped.

        .. note::

            The expression is not expanded or factored before the
            grouping takes place. For best results, call :meth:`expand`
            on the expression before :meth:`collect`.

        EXAMPLES:

        In the first term of `f`, `x` has a coefficient of `4y`. In
        the second term, `x` has a coefficient of `z`. Therefore, if
        we collect those coefficients, `x` will have a coefficient of
        `4y+z`::

            sage: x,y,z = var('x,y,z')
            sage: f = 4*x*y + x*z + 20*y^2 + 21*y*z + 4*z^2 + x^2*y^2*z^2
            sage: f.collect(x)
            x^2*y^2*z^2 + x*(4*y + z) + 20*y^2 + 21*y*z + 4*z^2

        Here we do the same thing for `y` and `z`; however, note that
        we don't factor the `y^{2}` and `z^{2}` terms before
        collecting coefficients::

            sage: f.collect(y)
            (x^2*z^2 + 20)*y^2 + (4*x + 21*z)*y + x*z + 4*z^2
            sage: f.collect(z)
            (x^2*y^2 + 4)*z^2 + 4*x*y + 20*y^2 + (x + 21*y)*z

        Sometimes, we do have to call :meth:`expand()` on the
        expression first to achieve the desired result::

            sage: f = (x + y)*(x - z)
            sage: f.collect(x)
            x^2 + x*y - x*z - y*z
            sage: f.expand().collect(x)
            x^2 + x*(y - z) - y*z

        TESTS:

        The output should be equivalent to the input::

            sage: polynomials = QQ['x']
            sage: f = SR(polynomials.random_element())
            sage: g = f.collect(x)
            sage: bool(f == g)
            True

        If ``s`` is not present in the given expression, the
        expression should not be modified. The variable `z` will not
        be present in `f` below since `f` is a random polynomial of
        maximum degree 10 in `x` and `y`::

            sage: z = var('z')
            sage: polynomials = QQ['x,y']
            sage: f = SR(polynomials.random_element(10))
            sage: g = f.collect(z)
            sage: bool(str(f) == str(g))
            True

        Check if :trac:`9046` is fixed::

            sage: var('a b x y z')
            (a, b, x, y, z)
            sage: p = -a*x^3 - a*x*y^2 + 2*b*x^2*y + 2*y^3 + x^2*z + y^2*z + x^2 + y^2 + a*x
            sage: p.collect(x)
            -a*x^3 + (2*b*y + z + 1)*x^2 + 2*y^3 + y^2*z - (a*y^2 - a)*x + y^2
        """
        cdef Expression s0 = self.coerce_in(s)
        cdef GEx x
        sig_on()
        try:
            x = self._gobj.collect(s0._gobj, False)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def collect_common_factors(self):
        """
        EXAMPLES::

            sage: var('x')
            x
            sage: (x/(x^2 + x)).collect_common_factors()
            1/(x + 1)
        """
        cdef GEx x
        sig_on()
        try:
            x = g_collect_common_factors(self._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def __abs__(self):
        """
        Return the absolute value of this expression.

        EXAMPLES::

            sage: var('x, y')
            (x, y)

        The absolute value of a symbolic expression::

            sage: abs(x^2+y^2)
            abs(x^2 + y^2)

        The absolute value of a number in the symbolic ring::

            sage: abs(SR(-5))
            5
            sage: type(abs(SR(-5)))
            <type 'sage.symbolic.expression.Expression'>

        Because this overrides a Python builtin function, we do not
        currently support a ``hold`` parameter to prevent automatic
        evaluation::

            sage: abs(SR(-5),hold=True)
            Traceback (most recent call last):
            ...
            TypeError: abs() takes no keyword arguments

        But this is possible using the method :meth:`abs`::

            sage: SR(-5).abs(hold=True)
            abs(-5)

        TESTS:

        Check if :trac:`11155` is fixed::

            sage: abs(pi+i)
            abs(pi + I)
        """
        return new_Expression_from_GEx(self._parent, g_abs(self._gobj))

    def abs(self, hold=False):
        """
        Return the absolute value of this expression.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: (x+y).abs()
            abs(x + y)

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(-5).abs(hold=True)
            abs(-5)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(-5).abs(hold=True); a.simplify()
            5

        TESTS:

        From :trac:`7557`::

            sage: var('y', domain='real')
            y
            sage: abs(exp(1.1*y*I)).simplify()
            1
            sage: var('y', domain='complex') # reset the domain for other tests
            y
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_abs, self._gobj, hold))

    def step(self, hold=False):
        """
        Return the value of the Heaviside step function, which is 0 for
        negative x, 1/2 for 0, and 1 for positive x.

        EXAMPLES::

            sage: x = var('x')
            sage: SR(1.5).step()
            1
            sage: SR(0).step()
            1/2
            sage: SR(-1/2).step()
            0
            sage: SR(float(-1)).step()
            0

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(2).step()
            1
            sage: SR(2).step(hold=True)
            step(2)

        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_step, self._gobj, hold))

    def csgn(self, hold=False):
        """
        Return the sign of self, which is -1 if self < 0, 0 if self ==
        0, and 1 if self > 0, or unevaluated when self is a nonconstant
        symbolic expression.

        If self is not real, return the complex half-plane (left or right)
        in which the number lies.  If self is pure imaginary, return the sign
        of the imaginary part of self.

        EXAMPLES::

            sage: x = var('x')
            sage: SR(-2).csgn()
            -1
            sage: SR(0.0).csgn()
            0
            sage: SR(10).csgn()
            1
            sage: x.csgn()
            csgn(x)
            sage: SR(CDF.0).csgn()
            1
            sage: SR(I).csgn()
            1
            sage: SR(-I).csgn()
            -1
            sage: SR(1+I).csgn()
            1
            sage: SR(1-I).csgn()
            1
            sage: SR(-1+I).csgn()
            -1
            sage: SR(-1-I).csgn()
            -1

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(I).csgn(hold=True)
            csgn(I)

        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_csgn, self._gobj, hold))

    def conjugate(self, hold=False):
        """
        Return the complex conjugate of this symbolic expression.

        EXAMPLES::

            sage: a = 1 + 2*I
            sage: a.conjugate()
            -2*I + 1
            sage: a = sqrt(2) + 3^(1/3)*I; a
            sqrt(2) + I*3^(1/3)
            sage: a.conjugate()
            sqrt(2) - I*3^(1/3)

            sage: SR(CDF.0).conjugate()
            -1.0*I
            sage: x.conjugate()
            conjugate(x)
            sage: SR(RDF(1.5)).conjugate()
            1.5
            sage: SR(float(1.5)).conjugate()
            1.5
            sage: SR(I).conjugate()
            -I
            sage: ( 1+I  + (2-3*I)*x).conjugate()
            (3*I + 2)*conjugate(x) - I + 1

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(I).conjugate(hold=True)
            conjugate(I)

        This also works in functional notation::

            sage: conjugate(I)
            -I
            sage: conjugate(I,hold=True)
            conjugate(I)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(I).conjugate(hold=True); a.simplify()
            -I

        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_conjugate, self._gobj, hold))

    def norm(self):
        r"""
        The complex norm of this symbolic expression, i.e.,
        the expression times its complex conjugate. If `c = a + bi` is a
        complex number, then the norm of `c` is defined as the product of
        `c` and its complex conjugate

        .. MATH::

            \text{norm}(c)
            =
            \text{norm}(a + bi)
            =
            c \cdot \overline{c}
            =
            a^2 + b^2.

        The norm of a complex number is different from its absolute value.
        The absolute value of a complex number is defined to be the square
        root of its norm. A typical use of the complex norm is in the
        integral domain `\ZZ[i]` of Gaussian integers, where the norm of
        each Gaussian integer `c = a + bi` is defined as its complex norm.

        .. SEEALSO::

            - :func:`sage.misc.functional.norm`

        EXAMPLES::

            sage: a = 1 + 2*I
            sage: a.norm()
            5
            sage: a = sqrt(2) + 3^(1/3)*I; a
            sqrt(2) + I*3^(1/3)
            sage: a.norm()
            3^(2/3) + 2
            sage: CDF(a).norm()
            4.08008382305
            sage: CDF(a.norm())
            4.08008382305
        """
        return (self*self.conjugate()).expand()

    def real_part(self, hold=False):
        """
        Return the real part of this symbolic expression.

        EXAMPLES::

            sage: x = var('x')
            sage: x.real_part()
            real_part(x)
            sage: SR(2+3*I).real_part()
            2
            sage: SR(CDF(2,3)).real_part()
            2.0
            sage: SR(CC(2,3)).real_part()
            2.00000000000000

            sage: f = log(x)
            sage: f.real_part()
            log(abs(x))

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(2).real_part()
            2
            sage: SR(2).real_part(hold=True)
            real_part(2)

        This also works using functional notation::

            sage: real_part(I,hold=True)
            real_part(I)
            sage: real_part(I)
            0

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(2).real_part(hold=True); a.simplify()
            2
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_real_part, self._gobj, hold))

    real = real_part

    def imag_part(self, hold=False):
        r"""
        Return the imaginary part of this symbolic expression.

        EXAMPLES::

            sage: sqrt(-2).imag_part()
            sqrt(2)

        We simplify `\ln(\exp(z))` to `z`.  This should only
        be for `-\pi<{\rm Im}(z)<=\pi`, but Maxima does not
        have a symbolic imaginary part function, so we cannot
        use ``assume`` to assume that first::

            sage: z = var('z')
            sage: f = log(exp(z))
            sage: f
            log(e^z)
            sage: f.simplify()
            z
            sage: forget()

        A more symbolic example::

            sage: var('a, b')
            (a, b)
            sage: f = log(a + b*I)
            sage: f.imag_part()
            arctan2(imag_part(a) + real_part(b), -imag_part(b) + real_part(a))

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: I.imag_part()
            1
            sage: I.imag_part(hold=True)
            imag_part(I)

        This also works using functional notation::

            sage: imag_part(I,hold=True)
            imag_part(I)
            sage: imag_part(I)
            1

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = I.imag_part(hold=True); a.simplify()
            1

        TESTS::

            sage: x = var('x')
            sage: x.imag_part()
            imag_part(x)
            sage: SR(2+3*I).imag_part()
            3
            sage: SR(CC(2,3)).imag_part()
            3.00000000000000
            sage: SR(CDF(2,3)).imag_part()
            3.0
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_imag_part, self._gobj, hold))

    imag = imag_part

    def sqrt(self, hold=False):
        """
        Return the square root of this expression

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: SR(2).sqrt()
            sqrt(2)
            sage: (x^2+y^2).sqrt()
            sqrt(x^2 + y^2)
            sage: (x^2).sqrt()
            sqrt(x^2)

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(4).sqrt()
            2
            sage: SR(4).sqrt(hold=True)
            sqrt(4)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(4).sqrt(hold=True); a.simplify()
            2

        To use this parameter in functional notation, you must coerce to
        the symbolic ring::

            sage: sqrt(SR(4),hold=True)
            sqrt(4)
            sage: sqrt(4,hold=True)
            Traceback (most recent call last):
            ...
            TypeError: _do_sqrt() got an unexpected keyword argument 'hold'
        """
        return new_Expression_from_GEx(self._parent,
                g_hold2_wrapper(g_power_construct, self._gobj, g_ex1_2, hold))

    def sin(self, hold=False):
        """
        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: sin(x^2 + y^2)
            sin(x^2 + y^2)
            sage: sin(sage.symbolic.constants.pi)
            0
            sage: sin(SR(1))
            sin(1)
            sage: sin(SR(RealField(150)(1)))
            0.84147098480789650665250232163029899962256306

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(0).sin()
            0
            sage: SR(0).sin(hold=True)
            sin(0)

        This also works using functional notation::

            sage: sin(0,hold=True)
            sin(0)
            sage: sin(0)
            0

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(0).sin(hold=True); a.simplify()
            0

        TESTS::

            sage: SR(oo).sin()
            Traceback (most recent call last):
            ...
            RuntimeError: sin_eval(): sin(infinity) encountered
            sage: SR(-oo).sin()
            Traceback (most recent call last):
            ...
            RuntimeError: sin_eval(): sin(infinity) encountered
            sage: SR(unsigned_infinity).sin()
            Traceback (most recent call last):
            ...
            RuntimeError: sin_eval(): sin(infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_sin, self._gobj, hold))

    def cos(self, hold=False):
        """
        Return the cosine of self.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: cos(x^2 + y^2)
            cos(x^2 + y^2)
            sage: cos(sage.symbolic.constants.pi)
            -1
            sage: cos(SR(1))
            cos(1)
            sage: cos(SR(RealField(150)(1)))
            0.54030230586813971740093660744297660373231042


        In order to get a numeric approximation use .n()::

            sage: SR(RR(1)).cos().n()
            0.540302305868140
            sage: SR(float(1)).cos().n()
            0.540302305868140

        To prevent automatic evaluation use the ``hold`` argument::

            sage: pi.cos()
            -1
            sage: pi.cos(hold=True)
            cos(pi)

        This also works using functional notation::

            sage: cos(pi,hold=True)
            cos(pi)
            sage: cos(pi)
            -1

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = pi.cos(hold=True); a.simplify()
            -1

        TESTS::

            sage: SR(oo).cos()
            Traceback (most recent call last):
            ...
            RuntimeError: cos_eval(): cos(infinity) encountered
            sage: SR(-oo).cos()
            Traceback (most recent call last):
            ...
            RuntimeError: cos_eval(): cos(infinity) encountered
            sage: SR(unsigned_infinity).cos()
            Traceback (most recent call last):
            ...
            RuntimeError: cos_eval(): cos(infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_cos, self._gobj, hold))

    def tan(self, hold=False):
        """
        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: tan(x^2 + y^2)
            tan(x^2 + y^2)
            sage: tan(sage.symbolic.constants.pi/2)
            Infinity
            sage: tan(SR(1))
            tan(1)
            sage: tan(SR(RealField(150)(1)))
            1.5574077246549022305069748074583601730872508

        To prevent automatic evaluation use the ``hold`` argument::

            sage: (pi/12).tan()
            -sqrt(3) + 2
            sage: (pi/12).tan(hold=True)
            tan(1/12*pi)

        This also works using functional notation::

            sage: tan(pi/12,hold=True)
            tan(1/12*pi)
            sage: tan(pi/12)
            -sqrt(3) + 2

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = (pi/12).tan(hold=True); a.simplify()
            -sqrt(3) + 2

        TESTS::

            sage: SR(oo).tan()
            Traceback (most recent call last):
            ...
            RuntimeError: tan_eval(): tan(infinity) encountered
            sage: SR(-oo).tan()
            Traceback (most recent call last):
            ...
            RuntimeError: tan_eval(): tan(infinity) encountered
            sage: SR(unsigned_infinity).tan()
            Traceback (most recent call last):
            ...
            RuntimeError: tan_eval(): tan(infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_tan, self._gobj, hold))

    def arcsin(self, hold=False):
        """
        Return the arcsin of x, i.e., the number y between -pi and pi
        such that sin(y) == x.

        EXAMPLES::

            sage: x.arcsin()
            arcsin(x)
            sage: SR(0.5).arcsin()
            0.523598775598299
            sage: SR(0.999).arcsin()
            1.52607123962616
            sage: SR(1/3).arcsin()
            arcsin(1/3)
            sage: SR(-1/3).arcsin()
            -arcsin(1/3)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(0).arcsin()
            0
            sage: SR(0).arcsin(hold=True)
            arcsin(0)

        This also works using functional notation::

            sage: arcsin(0,hold=True)
            arcsin(0)
            sage: arcsin(0)
            0

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(0).arcsin(hold=True); a.simplify()
            0

        TESTS::

            sage: SR(oo).arcsin()
            Traceback (most recent call last):
            ...
            RuntimeError: arcsin_eval(): arcsin(infinity) encountered
            sage: SR(-oo).arcsin()
            Traceback (most recent call last):
            ...
            RuntimeError: arcsin_eval(): arcsin(infinity) encountered
            sage: SR(unsigned_infinity).arcsin()
            Infinity
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_asin, self._gobj, hold))

    def arccos(self, hold=False):
        """
        Return the arc cosine of self.

        EXAMPLES::

            sage: x.arccos()
            arccos(x)
            sage: SR(1).arccos()
            0
            sage: SR(1/2).arccos()
            1/3*pi
            sage: SR(0.4).arccos()
            1.15927948072741
            sage: plot(lambda x: SR(x).arccos(), -1,1)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1).arccos(hold=True)
            arccos(1)

        This also works using functional notation::

            sage: arccos(1,hold=True)
            arccos(1)
            sage: arccos(1)
            0

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(1).arccos(hold=True); a.simplify()
            0

        TESTS::

            sage: SR(oo).arccos()
            Traceback (most recent call last):
            ...
            RuntimeError: arccos_eval(): arccos(infinity) encountered
            sage: SR(-oo).arccos()
            Traceback (most recent call last):
            ...
            RuntimeError: arccos_eval(): arccos(infinity) encountered
            sage: SR(unsigned_infinity).arccos()
            Infinity
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_acos, self._gobj, hold))

    def arctan(self, hold=False):
        """
        Return the arc tangent of self.

        EXAMPLES::

            sage: x = var('x')
            sage: x.arctan()
            arctan(x)
            sage: SR(1).arctan()
            1/4*pi
            sage: SR(1/2).arctan()
            arctan(1/2)
            sage: SR(0.5).arctan()
            0.463647609000806
            sage: plot(lambda x: SR(x).arctan(), -20,20)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1).arctan(hold=True)
            arctan(1)

        This also works using functional notation::

            sage: arctan(1,hold=True)
            arctan(1)
            sage: arctan(1)
            1/4*pi

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(1).arctan(hold=True); a.simplify()
            1/4*pi

        TESTS::

            sage: SR(oo).arctan()
            1/2*pi
            sage: SR(-oo).arctan()
            -1/2*pi
            sage: SR(unsigned_infinity).arctan()
            Traceback (most recent call last):
            ...
            RuntimeError: arctan_eval(): arctan(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_atan, self._gobj, hold))

    def arctan2(self, x, hold=False):
        """
        Return the inverse of the 2-variable tan function on self and x.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: x.arctan2(y)
            arctan2(x, y)
            sage: SR(1/2).arctan2(1/2)
            1/4*pi
            sage: maxima.eval('atan2(1/2,1/2)')
            '%pi/4'

            sage: SR(-0.7).arctan2(SR(-0.6))
            -2.27942259892257

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1/2).arctan2(1/2, hold=True)
            arctan2(1/2, 1/2)

        This also works using functional notation::

            sage: arctan2(1,2,hold=True)
            arctan2(1, 2)
            sage: arctan2(1,2)
            arctan(1/2)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(1/2).arctan2(1/2, hold=True); a.simplify()
            1/4*pi

        TESTS:

        We compare a bunch of different evaluation points between
        Sage and Maxima::

            sage: float(SR(0.7).arctan2(0.6))
            0.8621700546672264
            sage: maxima('atan2(0.7,0.6)')
            .862170054667226...
            sage: float(SR(0.7).arctan2(-0.6))
            2.279422598922567
            sage: maxima('atan2(0.7,-0.6)')
            2.279422598922567
            sage: float(SR(-0.7).arctan2(0.6))
            -0.8621700546672264
            sage: maxima('atan2(-0.7,0.6)')
            -.862170054667226...
            sage: float(SR(-0.7).arctan2(-0.6))
            -2.279422598922567
            sage: maxima('atan2(-0.7,-0.6)')
            -2.279422598922567
            sage: float(SR(0).arctan2(-0.6))
            3.141592653589793
            sage: maxima('atan2(0,-0.6)')
            3.141592653589793
            sage: float(SR(0).arctan2(0.6))
            0.0
            sage: maxima('atan2(0,0.6)')
            0.0
            sage: SR(0).arctan2(0) # see trac ticket #11423
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(0,0) encountered
            sage: SR(I).arctan2(1)
            arctan2(I, 1)
            sage: SR(CDF(0,1)).arctan2(1)
            arctan2(1.0*I, 1)
            sage: SR(1).arctan2(CDF(0,1))
            arctan2(1, 1.0*I)

            sage: arctan2(0,oo)
            0
            sage: SR(oo).arctan2(oo)
            1/4*pi
            sage: SR(oo).arctan2(0)
            1/2*pi
            sage: SR(-oo).arctan2(0)
            -1/2*pi
            sage: SR(-oo).arctan2(-2)
            pi
            sage: SR(unsigned_infinity).arctan2(2)
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(x, unsigned_infinity) encountered
            sage: SR(2).arctan2(oo)
            1/2*pi
            sage: SR(2).arctan2(-oo)
            -1/2*pi
            sage: SR(2).arctan2(SR(unsigned_infinity))
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(unsigned_infinity, x) encountered
        """
        cdef Expression nexp = self.coerce_in(x)
        return new_Expression_from_GEx(self._parent,
                g_hold2_wrapper(g_atan2, self._gobj, nexp._gobj, hold))

    def sinh(self, hold=False):
        r"""
        Return sinh of self.

        We have $\sinh(x) = (e^{x} - e^{-x})/2$.

        EXAMPLES::

            sage: x.sinh()
            sinh(x)
            sage: SR(1).sinh()
            sinh(1)
            sage: SR(0).sinh()
            0
            sage: SR(1.0).sinh()
            1.17520119364380
            sage: maxima('sinh(1.0)')
            1.17520119364380...

            sinh(1.0000000000000000000000000)
            sage: SR(1).sinh().n(90)
            1.1752011936438014568823819
            sage: SR(RIF(1)).sinh()
            1.175201193643802?

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arccosh(x).sinh()
            sqrt(x + 1)*sqrt(x - 1)
            sage: arccosh(x).sinh(hold=True)
            sinh(arccosh(x))

        This also works using functional notation::

            sage: sinh(arccosh(x),hold=True)
            sinh(arccosh(x))
            sage: sinh(arccosh(x))
            sqrt(x + 1)*sqrt(x - 1)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = arccosh(x).sinh(hold=True); a.simplify()
            sqrt(x + 1)*sqrt(x - 1)

        TESTS::

            sage: SR(oo).sinh()
            +Infinity
            sage: SR(-oo).sinh()
            -Infinity
            sage: SR(unsigned_infinity).sinh()
            Traceback (most recent call last):
            ...
            RuntimeError: sinh_eval(): sinh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_sinh, self._gobj, hold))

    def cosh(self, hold=False):
        r"""
        Return cosh of self.

        We have $\cosh(x) = (e^{x} + e^{-x})/2$.

        EXAMPLES::

            sage: x.cosh()
            cosh(x)
            sage: SR(1).cosh()
            cosh(1)
            sage: SR(0).cosh()
            1
            sage: SR(1.0).cosh()
            1.54308063481524
            sage: maxima('cosh(1.0)')
            1.54308063481524...
            sage: SR(1.00000000000000000000000000).cosh()
            1.5430806348152437784779056
            sage: SR(RIF(1)).cosh()
            1.543080634815244?

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arcsinh(x).cosh()
            sqrt(x^2 + 1)
            sage: arcsinh(x).cosh(hold=True)
            cosh(arcsinh(x))

        This also works using functional notation::

            sage: cosh(arcsinh(x),hold=True)
            cosh(arcsinh(x))
            sage: cosh(arcsinh(x))
            sqrt(x^2 + 1)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = arcsinh(x).cosh(hold=True); a.simplify()
            sqrt(x^2 + 1)

        TESTS::

            sage: SR(oo).cosh()
            +Infinity
            sage: SR(-oo).cosh()
            +Infinity
            sage: SR(unsigned_infinity).cosh()
            Traceback (most recent call last):
            ...
            RuntimeError: cosh_eval(): cosh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_cosh, self._gobj, hold))

    def tanh(self, hold=False):
        r"""
        Return tanh of self.

        We have $\tanh(x) = \sinh(x) / \cosh(x)$.

        EXAMPLES::

            sage: x.tanh()
            tanh(x)
            sage: SR(1).tanh()
            tanh(1)
            sage: SR(0).tanh()
            0
            sage: SR(1.0).tanh()
            0.761594155955765
            sage: maxima('tanh(1.0)')
            .7615941559557649
            sage: plot(lambda x: SR(x).tanh(), -1, 1)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arcsinh(x).tanh()
            x/sqrt(x^2 + 1)
            sage: arcsinh(x).tanh(hold=True)
            tanh(arcsinh(x))

        This also works using functional notation::

            sage: tanh(arcsinh(x),hold=True)
            tanh(arcsinh(x))
            sage: tanh(arcsinh(x))
            x/sqrt(x^2 + 1)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = arcsinh(x).tanh(hold=True); a.simplify()
            x/sqrt(x^2 + 1)

        TESTS::

            sage: SR(oo).tanh()
            1
            sage: SR(-oo).tanh()
            -1
            sage: SR(unsigned_infinity).tanh()
            Traceback (most recent call last):
            ...
            RuntimeError: tanh_eval(): tanh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_tanh, self._gobj, hold))

    def arcsinh(self, hold=False):
        """
        Return the inverse hyperbolic sine of self.

        EXAMPLES::

            sage: x.arcsinh()
            arcsinh(x)
            sage: SR(0).arcsinh()
            0
            sage: SR(1).arcsinh()
            arcsinh(1)
            sage: SR(1.0).arcsinh()
            0.881373587019543
            sage: maxima('asinh(2.0)')
            1.4436354751788...

        Sage automatically applies certain identities::

            sage: SR(3/2).arcsinh().cosh()
            1/2*sqrt(13)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(-2).arcsinh()
            -arcsinh(2)
            sage: SR(-2).arcsinh(hold=True)
            arcsinh(-2)

        This also works using functional notation::

            sage: arcsinh(-2,hold=True)
            arcsinh(-2)
            sage: arcsinh(-2)
            -arcsinh(2)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(-2).arcsinh(hold=True); a.simplify()
            -arcsinh(2)

        TESTS::

            sage: SR(oo).arcsinh()
            +Infinity
            sage: SR(-oo).arcsinh()
            -Infinity
            sage: SR(unsigned_infinity).arcsinh()
            Infinity
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_asinh, self._gobj, hold))

    def arccosh(self, hold=False):
        """
        Return the inverse hyperbolic cosine of self.

        EXAMPLES::

            sage: x.arccosh()
            arccosh(x)
            sage: SR(0).arccosh()
            1/2*I*pi
            sage: SR(1/2).arccosh()
            arccosh(1/2)
            sage: SR(CDF(1/2)).arccosh()
            1.0471975512*I
            sage: maxima('acosh(0.5)')
            1.04719755119659...*%i

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(-1).arccosh()
            I*pi
            sage: SR(-1).arccosh(hold=True)
            arccosh(-1)

        This also works using functional notation::

            sage: arccosh(-1,hold=True)
            arccosh(-1)
            sage: arccosh(-1)
            I*pi

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(-1).arccosh(hold=True); a.simplify()
            I*pi

        TESTS::

            sage: SR(oo).arccosh()
            +Infinity
            sage: SR(-oo).arccosh()
            +Infinity
            sage: SR(unsigned_infinity).arccosh()
            +Infinity
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_acosh, self._gobj, hold))

    def arctanh(self, hold=False):
        """
        Return the inverse hyperbolic tangent of self.

        EXAMPLES::

            sage: x.arctanh()
            arctanh(x)
            sage: SR(0).arctanh()
            0
            sage: SR(1/2).arctanh()
            arctanh(1/2)
            sage: SR(0.5).arctanh()
            0.549306144334055
            sage: SR(0.5).arctanh().tanh()
            0.500000000000000
            sage: maxima('atanh(0.5)')
            .5493061443340...

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(-1/2).arctanh()
            -arctanh(1/2)
            sage: SR(-1/2).arctanh(hold=True)
            arctanh(-1/2)

        This also works using functional notation::

            sage: arctanh(-1/2,hold=True)
            arctanh(-1/2)
            sage: arctanh(-1/2)
            -arctanh(1/2)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(-1/2).arctanh(hold=True); a.simplify()
            -arctanh(1/2)

        TESTS::

            sage: SR(1).arctanh()
            +Infinity
            sage: SR(-1).arctanh()
            -Infinity

            sage: SR(oo).arctanh()
            -1/2*I*pi
            sage: SR(-oo).arctanh()
            1/2*I*pi
            sage: SR(unsigned_infinity).arctanh()
            Traceback (most recent call last):
            ...
            RuntimeError: arctanh_eval(): arctanh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_atanh, self._gobj, hold))

    def exp(self, hold=False):
        """
        Return exponential function of self, i.e., e to the
        power of self.

        EXAMPLES::

            sage: x.exp()
            e^x
            sage: SR(0).exp()
            1
            sage: SR(1/2).exp()
            e^(1/2)
            sage: SR(0.5).exp()
            1.64872127070013
            sage: math.exp(0.5)
            1.6487212707001282

            sage: SR(0.5).exp().log()
            0.500000000000000
            sage: (pi*I).exp()
            -1

        To prevent automatic evaluation use the ``hold`` argument::

            sage: (pi*I).exp(hold=True)
            e^(I*pi)

        This also works using functional notation::

            sage: exp(I*pi,hold=True)
            e^(I*pi)
            sage: exp(I*pi)
            -1

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = (pi*I).exp(hold=True); a.simplify()
            -1

        TESTS:

        Test if #6377 is fixed::

            sage: SR(oo).exp()
            +Infinity
            sage: SR(-oo).exp()
            0
            sage: SR(unsigned_infinity).exp()
            Traceback (most recent call last):
            ...
            RuntimeError: exp_eval(): exp^(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_exp, self._gobj, hold))

    def log(self, b=None, hold=False):
        """
        Return the logarithm of self.

        EXAMPLES::

            sage: x, y = var('x, y')
            sage: x.log()
            log(x)
            sage: (x^y + y^x).log()
            log(x^y + y^x)
            sage: SR(0).log()
            -Infinity
            sage: SR(-1).log()
            I*pi
            sage: SR(1).log()
            0
            sage: SR(1/2).log()
            log(1/2)
            sage: SR(0.5).log()
            -0.693147180559945
            sage: SR(0.5).log().exp()
            0.500000000000000
            sage: math.log(0.5)
            -0.6931471805599453
            sage: plot(lambda x: SR(x).log(), 0.1,10)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: I.log()
            1/2*I*pi
            sage: I.log(hold=True)
            log(I)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = I.log(hold=True); a.simplify()
            1/2*I*pi

        We do not currently support a ``hold`` parameter in functional
        notation::

            sage: log(SR(-1),hold=True)
            Traceback (most recent call last):
            ...
            TypeError: log() got an unexpected keyword argument 'hold'

        TESTS::

            sage: SR(oo).log()
            +Infinity
            sage: SR(-oo).log()
            +Infinity
            sage: SR(unsigned_infinity).log()
            +Infinity
        """
        res = new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_log, self._gobj, hold))
        if b is None:
            return res
        else:
            return res/self.coerce_in(b).log(hold=hold)

    def zeta(self, hold=False):
        """
        EXAMPLES::

            sage: x, y = var('x, y')
            sage: (x/y).zeta()
            zeta(x/y)
            sage: SR(2).zeta()
            1/6*pi^2
            sage: SR(3).zeta()
            zeta(3)
            sage: SR(CDF(0,1)).zeta()
            0.00330022368532 - 0.418155449141*I
            sage: CDF(0,1).zeta()
            0.00330022368532 - 0.418155449141*I
            sage: plot(lambda x: SR(x).zeta(), -10,10).show(ymin=-3,ymax=3)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(2).zeta(hold=True)
            zeta(2)

        This also works using functional notation::

            sage: zeta(2,hold=True)
            zeta(2)
            sage: zeta(2)
            1/6*pi^2

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(2).zeta(hold=True); a.simplify()
            1/6*pi^2

        TESTS::

            sage: t = SR(1).zeta(); t
            Infinity
        """
        cdef GEx x
        sig_on()
        try:
            x = g_hold_wrapper(g_zeta, self._gobj, hold)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def factorial(self, hold=False):
        """
        Return the factorial of self.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: SR(5).factorial()
            120
            sage: x.factorial()
            factorial(x)
            sage: (x^2+y^3).factorial()
            factorial(y^3 + x^2)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(5).factorial(hold=True)
            factorial(5)

        This also works using functional notation::

            sage: factorial(5,hold=True)
            factorial(5)
            sage: factorial(5)
            120

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(5).factorial(hold=True); a.simplify()
            120
        """
        cdef GEx x
        sig_on()
        try:
            x = g_hold_wrapper(g_factorial, self._gobj, hold)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def binomial(self, k, hold=False):
        """
        Return binomial coefficient "self choose k".

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: SR(5).binomial(SR(3))
            10
            sage: x.binomial(SR(3))
            1/6*x^3 - 1/2*x^2 + 1/3*x
            sage: x.binomial(y)
            binomial(x, y)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: x.binomial(3, hold=True)
            binomial(x, 3)
            sage: SR(5).binomial(3, hold=True)
            binomial(5, 3)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(5).binomial(3, hold=True); a.simplify()
            10

        We do not currently support a ``hold`` parameter in functional
        notation::

            sage: binomial(5,3, hold=True)
            Traceback (most recent call last):
            ...
            TypeError: binomial() got an unexpected keyword argument 'hold'

        TESTS:

        Check if we handle zero correctly (#8561)::

            sage: x.binomial(0)
            1
            sage: SR(0).binomial(0)
            1
        """
        cdef Expression nexp = self.coerce_in(k)
        cdef GEx x
        sig_on()
        try:
            x = g_hold2_wrapper(g_binomial, self._gobj, nexp._gobj, hold)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def Order(self, hold=False):
        """
        Order, as in big oh notation.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: n = var('n')
            sage: t = (17*n^3).Order(); t
            Order(n^3)
            sage: t.derivative(n)
            Order(n^2)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: (17*n^3).Order(hold=True)
            Order(17*n^3)
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_Order, self._gobj, hold))

    def gamma(self, hold=False):
        """
        Return the Gamma function evaluated at self.

        EXAMPLES::

            sage: x = var('x')
            sage: x.gamma()
            gamma(x)
            sage: SR(2).gamma()
            1
            sage: SR(10).gamma()
            362880
            sage: SR(10.0r).gamma()
            362880.0
            sage: SR(CDF(1,1)).gamma()
            0.498015668118 - 0.154949828302*I

        ::

            sage: gp('gamma(1+I)')
            0.4980156681183560427136911175 - 0.1549498283018106851249551305*I # 32-bit
            0.49801566811835604271369111746219809195 - 0.15494982830181068512495513048388660520*I # 64-bit

        We plot the familiar plot of this log-convex function::

            sage: plot(gamma(x), -6,4).show(ymin=-3,ymax=3)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1/2).gamma()
            sqrt(pi)
            sage: SR(1/2).gamma(hold=True)
            gamma(1/2)

        This also works using functional notation::

            sage: gamma(1/2,hold=True)
            gamma(1/2)
            sage: gamma(1/2)
            sqrt(pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`simplify`::

            sage: a = SR(1/2).gamma(hold=True); a.simplify()
            sqrt(pi)
        """
        cdef GEx x
        sig_on()
        try:
            x = g_hold_wrapper(g_tgamma, self._gobj, hold)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def lgamma(self, hold=False):
        """
        This method is deprecated, please use the ``.log_gamma()`` function instead.

        Log gamma function evaluated at self.

        EXAMPLES::

            sage: x.lgamma()
            doctest:...: DeprecationWarning: The lgamma() function is deprecated. Use log_gamma() instead.
            See http://trac.sagemath.org/6992 for details.
            log_gamma(x)
        """
        from sage.misc.superseded import deprecation
        deprecation(6992, "The lgamma() function is deprecated. Use log_gamma() instead.")
        return self.log_gamma(hold=hold)

    def log_gamma(self, hold=False):
        """
        Return the log gamma function evaluated at self.
        This is the logarithm of gamma of self, where
        gamma is a complex function such that `gamma(n)`
        equals `factorial(n-1)`.

        EXAMPLES::

            sage: x = var('x')
            sage: x.log_gamma()
            log_gamma(x)
            sage: SR(2).log_gamma()
            0
            sage: SR(5).log_gamma()
            log(24)
            sage: a = SR(5).log_gamma(); a.n()
            3.17805383034795
            sage: SR(5-1).factorial().log()
            log(24)
            sage: set_verbose(-1); plot(lambda x: SR(x).log_gamma(), -7,8, plot_points=1000).show()
            sage: math.exp(0.5)
            1.6487212707001282
            sage: plot(lambda x: (SR(x).exp() - SR(-x).exp())/2 - SR(x).sinh(), -1, 1)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(5).log_gamma(hold=True)
            log_gamma(5)

        To evaluate again, currently we must use numerical evaluation
        via :meth:`n`::

            sage: a = SR(5).log_gamma(hold=True); a.n()
            3.17805383034795
        """
        cdef GEx x
        sig_on()
        try:
            x = g_hold_wrapper(g_lgamma, self._gobj, hold)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def default_variable(self):
        """
        Return the default variable, which is by definition the first
        variable in self, or `x` is there are no variables in self.
        The result is cached.

        EXAMPLES::

            sage: sqrt(2).default_variable()
            x
            sage: x, theta, a = var('x, theta, a')
            sage: f = x^2 + theta^3 - a^x
            sage: f.default_variable()
            a

        Note that this is the first *variable*, not the first *argument*::

            sage: f(theta, a, x) = a + theta^3
            sage: f.default_variable()
            a
            sage: f.variables()
            (a, theta)
            sage: f.arguments()
            (theta, a, x)
        """
        v = self.variables()
        if len(v) == 0:
            return self.parent().var('x')
        else:
            return v[0]

    def combine(self):
        r"""
        Returns a simplified version of this symbolic expression
        by combining all terms with the same denominator into a single
        term.

        EXAMPLES::

            sage: var('x, y, a, b, c')
            (x, y, a, b, c)
            sage: f = x*(x-1)/(x^2 - 7) + y^2/(x^2-7) + 1/(x+1) + b/a + c/a; f
            (x - 1)*x/(x^2 - 7) + y^2/(x^2 - 7) + b/a + c/a + 1/(x + 1)
            sage: f.combine()
            ((x - 1)*x + y^2)/(x^2 - 7) + (b + c)/a + 1/(x + 1)
        """
        return self.parent()(self._maxima_().combine())

    def normalize(self):
        """
        Return this expression normalized as a fraction

        .. SEEALSO:

            :meth:`numerator`, :meth:`denominator`,
            :meth:`numerator_denominator`, :meth:`combine`

        EXAMPLES::

            sage: var('x, y, a, b, c')
            (x, y, a, b, c)
            sage: g = x + y/(x + 2)
            sage: g.normalize()
            (x^2 + 2*x + y)/(x + 2)

            sage: f = x*(x-1)/(x^2 - 7) + y^2/(x^2-7) + 1/(x+1) + b/a + c/a
            sage: f.normalize()
            (a*x^3 + b*x^3 + c*x^3 + a*x*y^2 + a*x^2 + b*x^2 + c*x^2 +
                    a*y^2 - a*x - 7*b*x - 7*c*x - 7*a - 7*b - 7*c)/((x^2 -
                        7)*a*(x + 1))

        ALGORITHM: Uses GiNaC.

        """
        return new_Expression_from_GEx(self._parent, self._gobj.normal())

    def numerator(self, bint normalize = True):
        """
        Returns the numerator of this symbolic expression

        INPUT:

        - ``normalize`` -- (default: ``True``) a boolean.

        If ``normalize`` is ``True``, the expression is first normalized to
        have it as a fraction before getting the numerator.

        If ``normalize`` is ``False``, the expression is kept and if it is not
        a quotient, then this will return the expression itself.

        .. SEEALSO::

            :meth:`normalize`, :meth:`denominator`,
            :meth:`numerator_denominator`, :meth:`combine`

        EXAMPLES::

            sage: a, x, y = var('a,x,y')
            sage: f = x*(x-a)/((x^2 - y)*(x-a)); f
            x/(x^2 - y)
            sage: f.numerator()
            x
            sage: f.denominator()
            x^2 - y
            sage: f.numerator(normalize=False)
            x
            sage: f.denominator(normalize=False)
            x^2 - y

            sage: y = var('y')
            sage: g = x + y/(x + 2); g
            x + y/(x + 2)
            sage: g.numerator()
            x^2 + 2*x + y
            sage: g.denominator()
            x + 2
            sage: g.numerator(normalize=False)
            x + y/(x + 2)
            sage: g.denominator(normalize=False)
            1

        TESTS::

            sage: ((x+y)^2/(x-y)^3*x^3).numerator(normalize=False)
            (x + y)^2*x^3
            sage: ((x+y)^2*x^3).numerator(normalize=False)
            (x + y)^2*x^3
            sage: (y/x^3).numerator(normalize=False)
            y
            sage: t = y/x^3/(x+y)^(1/2); t
            y/(sqrt(x + y)*x^3)
            sage: t.numerator(normalize=False)
            y
            sage: (1/x^3).numerator(normalize=False)
            1
            sage: (x^3).numerator(normalize=False)
            x^3
            sage: (y*x^sin(x)).numerator(normalize=False)
            Traceback (most recent call last):
            ...
            TypeError: self is not a rational expression
        """
        cdef GExVector vec
        cdef GEx oper, power
        if normalize:
            return new_Expression_from_GEx(self._parent, self._gobj.numer())
        elif is_a_mul(self._gobj):
            for i from 0 <= i < self._gobj.nops():
                oper = self._gobj.op(i)
                if not is_a_power(oper):
                    vec.push_back(oper)
                else:
                    power = oper.op(1)
                    if not is_a_numeric(power):
                        raise TypeError, "self is not a rational expression"
                    elif ex_to_numeric(power).is_positive():
                        vec.push_back(oper)
            return new_Expression_from_GEx(self._parent,
                                           g_mul_construct(vec, True))
        elif is_a_power(self._gobj):
            power = self._gobj.op(1)
            if is_a_numeric(power) and ex_to_numeric(power).is_negative():
                return self._parent.one()
        return self

    def denominator(self, bint normalize=True):
        """
        Returns the denominator of this symbolic expression

        INPUT:

        - ``normalize`` -- (default: ``True``) a boolean.

        If ``normalize`` is ``True``, the expression is first normalized to
        have it as a fraction before getting the denominator.

        If ``normalize`` is ``False``, the expression is kept and if it is not
        a quotient, then this will just return 1.

        .. SEEALSO::

            :meth:`normalize`, :meth:`numerator`,
            :meth:`numerator_denominator`, :meth:`combine`

        EXAMPLES::

            sage: x, y, z, theta = var('x, y, z, theta')
            sage: f = (sqrt(x) + sqrt(y) + sqrt(z))/(x^10 - y^10 - sqrt(theta))
            sage: f.numerator()
            sqrt(x) + sqrt(y) + sqrt(z)
            sage: f.denominator()
            x^10 - y^10 - sqrt(theta)

            sage: f.numerator(normalize=False)
            (sqrt(x) + sqrt(y) + sqrt(z))
            sage: f.denominator(normalize=False)
            x^10 - y^10 - sqrt(theta)

            sage: y = var('y')
            sage: g = x + y/(x + 2); g
            x + y/(x + 2)
            sage: g.numerator(normalize=False)
            x + y/(x + 2)
            sage: g.denominator(normalize=False)
            1

        TESTS::

            sage: ((x+y)^2/(x-y)^3*x^3).denominator(normalize=False)
            (x - y)^3
            sage: ((x+y)^2*x^3).denominator(normalize=False)
            1
            sage: (y/x^3).denominator(normalize=False)
            x^3
            sage: t = y/x^3/(x+y)^(1/2); t
            y/(sqrt(x + y)*x^3)
            sage: t.denominator(normalize=False)
            sqrt(x + y)*x^3
            sage: (1/x^3).denominator(normalize=False)
            x^3
            sage: (x^3).denominator(normalize=False)
            1
            sage: (y*x^sin(x)).denominator(normalize=False)
            Traceback (most recent call last):
            ...
            TypeError: self is not a rational expression
        """
        cdef GExVector vec
        cdef GEx oper, ex, power
        if normalize:
            return new_Expression_from_GEx(self._parent, self._gobj.denom())
        elif is_a_mul(self._gobj):
            for i from 0 <= i < self._gobj.nops():
                oper = self._gobj.op(i)
                if is_a_power(oper):
                    ex = oper.op(0)
                    power = oper.op(1)
                    if not is_a_numeric(power):
                        raise TypeError, "self is not a rational expression"
                    elif ex_to_numeric(power).is_negative():
                        vec.push_back(g_pow(ex, g_abs(power)))
            return new_Expression_from_GEx(self._parent,
                                           g_mul_construct(vec, False))
        elif is_a_power(self._gobj):
            power = self._gobj.op(1)
            if is_a_numeric(power) and ex_to_numeric(power).is_negative():
                return new_Expression_from_GEx(self._parent,
                        g_pow(self._gobj.op(0), g_abs(power)))

        return self._parent.one()

    def numerator_denominator(self, bint normalize=True):
        """
        Returns the numerator and the denominator of this symbolic expression

        INPUT:

        - ``normalize`` -- (default: ``True``) a boolean.

        If ``normalize`` is ``True``, the expression is first normalized to
        have it as a fraction before getting the numerator and denominator.

        If ``normalize`` is ``False``, the expression is kept and if it is not
        a quotient, then this will return the expression itself together with
        1.

        .. SEEALSO::

            :meth:`normalize`, :meth:`numerator`, :meth:`denominator`,
            :meth:`combine`

        EXAMPLE::

            sage: x, y, a = var("x y a")
            sage: ((x+y)^2/(x-y)^3*x^3).numerator_denominator()
            ((x + y)^2*x^3, (x - y)^3)

            sage: ((x+y)^2/(x-y)^3*x^3).numerator_denominator(False)
            ((x + y)^2*x^3, (x - y)^3)

            sage: g = x + y/(x + 2)
            sage: g.numerator_denominator()
            (x^2 + 2*x + y, x + 2)
            sage: g.numerator_denominator(normalize=False)
            (x + y/(x + 2), 1)

            sage: g = x^2*(x + 2)
            sage: g.numerator_denominator()
            ((x + 2)*x^2, 1)
            sage: g.numerator_denominator(normalize=False)
            ((x + 2)*x^2, 1)

        TESTS::

            sage: ((x+y)^2/(x-y)^3*x^3).numerator_denominator(normalize=False)
            ((x + y)^2*x^3, (x - y)^3)
            sage: ((x+y)^2*x^3).numerator_denominator(normalize=False)
            ((x + y)^2*x^3, 1)
            sage: (y/x^3).numerator_denominator(normalize=False)
            (y, x^3)
            sage: t = y/x^3/(x+y)^(1/2); t
            y/(sqrt(x + y)*x^3)
            sage: t.numerator_denominator(normalize=False)
            (y, sqrt(x + y)*x^3)
            sage: (1/x^3).numerator_denominator(normalize=False)
            (1, x^3)
            sage: (x^3).numerator_denominator(normalize=False)
            (x^3, 1)
            sage: (y*x^sin(x)).numerator_denominator(normalize=False)
            Traceback (most recent call last):
            ...
            TypeError: self is not a rational expression
        """
        cdef GExVector vecnumer, vecdenom
        cdef GEx oper, ex, power
        cdef GNumeric power_num
        if normalize:
            ex = self._gobj.numer_denom()
            return (new_Expression_from_GEx(self._parent, ex.op(0)),
                    new_Expression_from_GEx(self._parent, ex.op(1)))
        elif is_a_mul(self._gobj):
            for i from 0 <= i < self._gobj.nops():
                oper = self._gobj.op(i)
                if is_a_power(oper):   # oper = ex^power
                    ex = oper.op(0)
                    power = oper.op(1)
                    if not is_a_numeric(power):
                        raise TypeError, "self is not a rational expression"
                    elif is_a_numeric(power):
                        power_num = ex_to_numeric(power)
                        if power_num.is_positive():
                            vecnumer.push_back(oper)
                        else:
                            vecdenom.push_back(g_pow(ex, g_abs(power)))
                else:
                    vecnumer.push_back(oper)
            return (new_Expression_from_GEx(self._parent,
                                            g_mul_construct(vecnumer, False)),
                    new_Expression_from_GEx(self._parent,
                                            g_mul_construct(vecdenom, False)))
        elif is_a_power(self._gobj):
            power = self._gobj.op(1)
            if is_a_numeric(power) and ex_to_numeric(power).is_positive():
                return (self, self._parent.one())
            else:
                return (self._parent.one(),
                        new_Expression_from_GEx(self._parent,
                               g_pow(self._gobj.op(0), g_abs(power))))
        else:
            return (self, self._parent.one())

    def partial_fraction(self, var=None):
        r"""
        Return the partial fraction expansion of ``self`` with
        respect to the given variable.

        INPUT:


        -  ``var`` - variable name or string (default: first
           variable)


        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: f = x^2/(x+1)^3
            sage: f.partial_fraction()
            1/(x + 1) - 2/(x + 1)^2 + 1/(x + 1)^3
            sage: f.partial_fraction()
            1/(x + 1) - 2/(x + 1)^2 + 1/(x + 1)^3

        Notice that the first variable in the expression is used by
        default::

            sage: y = var('y')
            sage: f = y^2/(y+1)^3
            sage: f.partial_fraction()
            1/(y + 1) - 2/(y + 1)^2 + 1/(y + 1)^3

            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f.partial_fraction()
            y^2/(y^3 + 3*y^2 + 3*y + 1) + 1/(x - 1)^2 + 1/(x - 1)^3

        You can explicitly specify which variable is used::

            sage: f.partial_fraction(y)
            x/(x^3 - 3*x^2 + 3*x - 1) + 1/(y + 1) - 2/(y + 1)^2 + 1/(y + 1)^3
        """
        if var is None:
            var = self.default_variable()
        return self.parent()(self._maxima_().partfrac(var))

    def maxima_methods(self):
        """
        Provides easy access to maxima methods, converting the result to a
        Sage expression automatically.

        EXAMPLES::

            sage: t = log(sqrt(2) - 1) + log(sqrt(2) + 1); t
            log(sqrt(2) + 1) + log(sqrt(2) - 1)
            sage: res = t.maxima_methods().logcontract(); res
            log((sqrt(2) + 1)*(sqrt(2) - 1))
            sage: type(res)
            <type 'sage.symbolic.expression.Expression'>
        """
        from sage.symbolic.maxima_wrapper import MaximaWrapper
        return MaximaWrapper(self)

    def rectform(self):
        r"""
        Convert this symbolic expression to rectangular form; that
        is, the form `a + bi` where `a` and `b` are real numbers and
        `i` is the imaginary unit.

        .. note::

           The name \"rectangular\" comes from the fact that, in the
           complex plane, `a` and `bi` are perpendicular.

        INPUT:

        - ``self`` -- the expression to convert.

        OUTPUT:

        A new expression, equivalent to the original, but expressed in
        the form `a + bi`.

        ALGORITHM:

        We call Maxima's ``rectform()`` and return the result unmodified.

        EXAMPLES:

        The exponential form of `\sin(x)`::

            sage: f = (e^(I*x) - e^(-I*x)) / (2*I)
            sage: f.rectform()
            sin(x)

        And `\cos(x)`::

            sage: f = (e^(I*x) + e^(-I*x)) / 2
            sage: f.rectform()
            cos(x)

        In some cases, this will simplify the given expression. For
        example, here, `e^{ik\pi}`, `\sin(k\pi)=0` should cancel
        leaving only `\cos(k\pi)` which can then be simplified::

            sage: k = var('k')
            sage: assume(k, 'integer')
            sage: f = e^(I*pi*k)
            sage: f.rectform()
            (-1)^k

        However, in general, the resulting expression may be more
        complicated than the original::

            sage: f = e^(I*x)
            sage: f.rectform()
            cos(x) + I*sin(x)

        TESTS:

        If the expression is already in rectangular form, it should be
        left alone::

            sage: a,b = var('a,b')
            sage: assume((a, 'real'), (b, 'real'))
            sage: f = a + b*I
            sage: f.rectform()
            a + I*b
            sage: forget()

        We can check with specific real numbers::

            sage: a = RR.random_element()
            sage: b = RR.random_element()
            sage: f = a + b*I
            sage: bool(f.rectform() == a + b*I)
            True

        If we decompose a complex number into its real and imaginary
        parts, they should correspond to the real and imaginary terms
        of the rectangular form::

            sage: z = CC.random_element()
            sage: a = z.real_part()
            sage: b = z.imag_part()
            sage: bool(SR(z).rectform() == a + b*I)
            True

        """
        return self.maxima_methods().rectform()

    def simplify(self):
        """
        Returns a simplified version of this symbolic expression.

        .. note::

           Currently, this just sends the expression to Maxima
           and converts it back to Sage.

        .. seealso::

           :meth:`simplify_full`, :meth:`simplify_trig`,
           :meth:`simplify_rational`, :meth:`simplify_radical`,
           :meth:`simplify_factorial`, :meth:`simplify_log`

        EXAMPLES::

            sage: a = var('a'); f = x*sin(2)/(x^a); f
            x*sin(2)/x^a
            sage: f.simplify()
            x^(-a + 1)*sin(2)
        """
        return self._parent(self._maxima_())

    def simplify_full(self):
        """
        Applies simplify_factorial, simplify_trig, simplify_rational,
        simplify_radical, simplify_log, and again simplify_rational to
        self (in that order).

        ALIAS: simplify_full and full_simplify are the same.

        EXAMPLES::

            sage: f = sin(x)^2 + cos(x)^2
            sage: f.simplify_full()
            1

        ::

            sage: f = sin(x/(x^2 + x))
            sage: f.simplify_full()
            sin(1/(x + 1))

        ::

            sage: var('n,k')
            (n, k)
            sage: f = binomial(n,k)*factorial(k)*factorial(n-k)
            sage: f.simplify_full()
            factorial(n)

        TESTS:

        There are two square roots of `$(x + 1)^2$`, so this should
        not be simplified to `$x + 1$`, :trac:`12737`::

            sage: f = sqrt((x + 1)^2)
            sage: f.simplify_full()
            sqrt(x^2 + 2*x + 1)

        The imaginary part of an expression should not change under
        simplification; :trac:`11934`::

            sage: f = sqrt(-8*(4*sqrt(2) - 7)*x^4 + 16*(3*sqrt(2) - 5)*x^3)
            sage: original = f.imag_part()
            sage: simplified = f.full_simplify().imag_part()
            sage: original - simplified
            0

        """
        x = self
        x = x.simplify_factorial()
        x = x.simplify_trig()
        x = x.simplify_rational()
        x = x.simplify_log('one')
        x = x.simplify_rational()
        return x

    full_simplify = simplify_full

    def simplify_trig(self,expand=True):
        r"""
        Optionally expands and then employs identities such as
        `\sin(x)^2 + \cos(x)^2 = 1`, `\cosh(x)^2 - \sinh(x)^2 = 1`,
        `\sin(x)\csc(x) = 1`, or `\tanh(x)=\sinh(x)/\cosh(x)`
        to simplify expressions containing tan, sec, etc., to sin,
        cos, sinh, cosh.

        INPUT:

        - ``self`` - symbolic expression

        - ``expand`` - (default:True) if True, expands trigonometric
          and hyperbolic functions of sums of angles and of multiple
          angles occurring in ``self`` first. For best results,
          ``self`` should be expanded. See also :meth:`expand_trig` to
          get more controls on this expansion.

        ALIAS: :meth:`trig_simplify` and :meth:`simplify_trig` are the same

        EXAMPLES::

            sage: f = sin(x)^2 + cos(x)^2; f
            cos(x)^2 + sin(x)^2
            sage: f.simplify()
            cos(x)^2 + sin(x)^2
            sage: f.simplify_trig()
            1
            sage: h = sin(x)*csc(x)
            sage: h.simplify_trig()
            1
            sage: k = tanh(x)*cosh(2*x)
            sage: k.simplify_trig()
            (2*sinh(x)^3 + sinh(x))/cosh(x)

        In some cases we do not want to expand::

            sage: f=tan(3*x)
            sage: f.simplify_trig()
            (4*cos(x)^2 - 1)*sin(x)/(4*cos(x)^3 - 3*cos(x))
            sage: f.simplify_trig(False)
            sin(3*x)/cos(3*x)

        """
        # much better to expand first, since it often doesn't work
        # right otherwise!
        if expand:
            return self.parent()(self._maxima_().trigexpand().trigsimp())
        else:
            return self.parent()(self._maxima_().trigsimp())

    trig_simplify = simplify_trig

    @rename_keyword(deprecation=6094, method="algorithm")
    def simplify_rational(self,algorithm='full', map=False):
        r"""
        Simplify rational expressions.

        INPUT:

        - ``self`` - symbolic expression

        - ``algorithm`` - (default: 'full') string which switches the
          algorithm for simplifications. Possible values are

          - 'simple' (simplify rational functions into quotient of two
            polynomials),

          - 'full' (apply repeatedly, if necessary)

          - 'noexpand' (convert to commmon denominator and add)

        - ``map`` - (default: False) if True, the result is an
          expression whose leading operator is the same as that of the
          expression ``self`` but whose subparts are the results of
          applying simplification rules to the corresponding subparts
          of the expressions.

        ALIAS: :meth:`rational_simplify` and :meth:`simplify_rational`
        are the same

        DETAILS: We call Maxima functions ratsimp, fullratsimp and
        xthru. If each part of the expression has to be simplified
        separately, we use Maxima function map.

        EXAMPLES::

            sage: f = sin(x/(x^2 + x))
            sage: f
            sin(x/(x^2 + x))
            sage: f.simplify_rational()
            sin(1/(x + 1))

        ::

            sage: f = ((x - 1)^(3/2) - (x + 1)*sqrt(x - 1))/sqrt((x - 1)*(x + 1)); f
            -((x + 1)*sqrt(x - 1) - (x - 1)^(3/2))/sqrt((x + 1)*(x - 1))
            sage: f.simplify_rational()
            -2*sqrt(x - 1)/sqrt(x^2 - 1)

        With ``map=True`` each term in a sum is simplified separately
        and thus the resuls are shorter for functions which are
        combination of rational and nonrational funtions. In the
        following example, we use this option if we want not to
        combine logarithm and the rational function into one
        fraction::

            sage: f=(x^2-1)/(x+1)-ln(x)/(x+2)
            sage: f.simplify_rational()
            (x^2 + x - log(x) - 2)/(x + 2)
            sage: f.simplify_rational(map=True)
            x - log(x)/(x + 2) - 1

        Here is an example from the Maxima documentation of where
        ``algorithm='simple'`` produces an (possibly useful) intermediate
        step::

            sage: y = var('y')
            sage: g = (x^(y/2) + 1)^2*(x^(y/2) - 1)^2/(x^y - 1)
            sage: g.simplify_rational(algorithm='simple')
            (x^(2*y) - 2*x^y + 1)/(x^y - 1)
            sage: g.simplify_rational()
            x^y - 1

        With option ``algorithm='noexpand'`` we only convert to common
        denominators and add. No expansion of products is performed::

            sage: f=1/(x+1)+x/(x+2)^2
            sage: f.simplify_rational()
            (2*x^2 + 5*x + 4)/(x^3 + 5*x^2 + 8*x + 4)
            sage: f.simplify_rational(algorithm='noexpand')
            ((x + 2)^2 + (x + 1)*x)/((x + 2)^2*(x + 1))
        """
        self_m = self._maxima_()
        if algorithm == 'full':
            maxima_method = 'fullratsimp'
        elif algorithm == 'simple':
            maxima_method = 'ratsimp'
        elif algorithm == 'noexpand':
            maxima_method = 'xthru'
        else:
            raise NotImplementedError, "unknown algorithm, see the help for available algorithms"
        P = self_m.parent()
        self_str=self_m.str()
        if map:
            cmd = "if atom(%s) then %s(%s) else map(%s,%s)"%(self_str,maxima_method,self_str,maxima_method,self_str)
        else:
            cmd = "%s(%s)"%(maxima_method,self_m.str())
        res = P(cmd)
        return self.parent()(res)

    rational_simplify = simplify_rational

    def simplify_factorial(self):
        """
        Simplify by combining expressions with factorials, and by
        expanding binomials into factorials.

        ALIAS: factorial_simplify and simplify_factorial are the same

        EXAMPLES:

        Some examples are relatively clear::

            sage: var('n,k')
            (n, k)
            sage: f = factorial(n+1)/factorial(n); f
            factorial(n + 1)/factorial(n)
            sage: f.simplify_factorial()
            n + 1

        ::

            sage: f = factorial(n)*(n+1); f
            (n + 1)*factorial(n)
            sage: simplify(f)
            (n + 1)*factorial(n)
            sage: f.simplify_factorial()
            factorial(n + 1)

        ::

            sage: f = binomial(n, k)*factorial(k)*factorial(n-k); f
            binomial(n, k)*factorial(k)*factorial(-k + n)
            sage: f.simplify_factorial()
            factorial(n)

        A more complicated example, which needs further processing::

            sage: f = factorial(x)/factorial(x-2)/2 + factorial(x+1)/factorial(x)/2; f
            1/2*factorial(x + 1)/factorial(x) + 1/2*factorial(x)/factorial(x - 2)
            sage: g = f.simplify_factorial(); g
            1/2*(x - 1)*x + 1/2*x + 1/2
            sage: g.simplify_rational()
            1/2*x^2 + 1/2


        TESTS:

        Check that the problem with applying full_simplify() to gamma functions (Trac 9240)
        has been fixed::

            sage: gamma(1/3)
            gamma(1/3)
            sage: gamma(1/3).full_simplify()
            gamma(1/3)
            sage: gamma(4/3)
            gamma(4/3)
            sage: gamma(4/3).full_simplify()
            1/3*gamma(1/3)

        """
        return self.parent()(self._maxima_().makefact().factcomb().minfactorial())

    factorial_simplify = simplify_factorial

    def simplify_radical(self):
        r"""
        Simplifies this symbolic expression, which can contain logs,
        exponentials, and radicals, by trying to convert it into a canonical
        form over a large class of expressions and a given ordering of
        variables.

        .. WARNING::

            As shown in the examples below, a canonical form is not always
            returned, i.e., two mathematically identical expressions might
            be simplified to different expressions.

        ALGORITHM:

        This uses the Maxima ``radcan()`` command. From the Maxima
        documentation: "All functionally equivalent forms are mapped into a
        unique form. For a somewhat larger class of expressions, produces a
        regular form. Two equivalent expressions in this class do not
        necessarily have the same appearance, but their difference can be
        simplified by radcan to zero. For some expressions radcan is quite
        time consuming. This is the cost of exploring certain relationships
        among the components of the expression for simplifications based on
        factoring and partial fraction expansions of exponents."

        .. NOTE::

            :meth:`radical_simplify`, :meth:`simplify_radical`,
            :meth:`exp_simplify`, :meth:`simplify_exp` are all the same.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)

        ::

            sage: f = log(x*y)
            sage: f.simplify_radical()
            log(x) + log(y)

        ::

            sage: f = log(8)/log(2)
            sage: f.simplify_radical()
            3

        ::

            sage: f = (log(x+x^2)-log(x))^a/log(1+x)^(a/2)
            sage: f.simplify_radical()
            log(x + 1)^(1/2*a)

        ::

            sage: f = (e^x-1)/(1+e^(x/2))
            sage: f.simplify_exp()
            e^(1/2*x) - 1

        The example below shows two expressions e1 and e2 which are
        "simplified" to different expressions, while their difference is
        "simplified" to zero, thus ``simplify_radical`` does not return a
        canonical form, except maybe for 0. ::

            sage: e1 = 1/(sqrt(5)+sqrt(2))
            sage: e2 = (sqrt(5)-sqrt(2))/3
            sage: e1.simplify_radical()
            1/(sqrt(5) + sqrt(2))
            sage: e2.simplify_radical()
            1/3*sqrt(5) - 1/3*sqrt(2)
            sage: (e1-e2).simplify_radical()
            0

        TESTS:

        This tests that :trac:`11668` has been fixed (by :trac:`12780`)::

            sage: a,b = var('a b')
            sage: A = abs((a+I*b))^2
            sage: A.simplify_radical()
            abs(a + I*b)^2
            sage: imag(A)
            0
            sage: imag(A.simplify_radical())
            0
        """
        from sage.calculus.calculus import maxima
        res = self.parent()(self._maxima_().radcan())
        return res

    radical_simplify = simplify_radical
    simplify_exp = exp_simplify = simplify_radical

    @rename_keyword(deprecation=6094, method="algorithm")
    def simplify_log(self,algorithm=None):
        r"""
        Simplifies symbolic expression, which can contain logs.

        Recursively scans the expression self, transforming
        subexpressions of the form `a1 \log(b1) + a2 \log(b2) + c` into
        `\log( b1^{a1} b2^{a2} ) + c` and simplifies inside logarithm. User
        can specify, which conditions must satisfy `a1` and `a2` to use
        this transformation in optional parameter ``algorithm``.

        INPUT:

        - ``self`` - expression to be simplified

        - ``algorithm`` - (default: None) optional, governs the condition
          on `a1` and `a2` which must be satisfied to contract expression
          `a1 \log(b1) + a2 \log(b2)`. Values are

          - None (use Maxima default, integers),

          - 'one' (1 and -1),

          - 'ratios' (integers and fractions of integers),

          - 'constants' (constants),

          - 'all' (all expressions).

          See also examples below.

        DETAILS: This uses the Maxima logcontract() command. From the
        Maxima documentation: "Recursively scans the expression expr,
        transforming subexpressions of the form a1*log(b1) +
        a2*log(b2) + c into log(ratsimp(b1^a1 * b2^a2)) + c. The user
        can control which coefficients are contracted by setting the
        option logconcoeffp to the name of a predicate function of one
        argument. E.g. if you like to generate SQRTs, you can do
        logconcoeffp:'logconfun$ logconfun(m):=featurep(m,integer) or
        ratnump(m)$ . Then logcontract(1/2*log(x)); will give
        log(sqrt(x))."

        ALIAS: :meth:`log_simplify` and :meth:`simplify_log` are the
        same

        EXAMPLES::

            sage: x,y,t=var('x y t')

        Only two first terms are contracted in the following example ,
        the logarithm with coefficient 1/2 is not contracted::

            sage: f = log(x)+2*log(y)+1/2*log(t)
            sage: f.simplify_log()
            log(x*y^2) + 1/2*log(t)

        To contract all terms in previous example use option ``algorithm``::

            sage: f.simplify_log(algorithm='ratios')
            log(sqrt(t)*x*y^2)

        This shows that the option ``algorithm`` from the previous call
        has no influence to future calls (we changed some default
        Maxima flag, and have to ensure that this flag has been
        restored)::

            sage: f.simplify_log('one')
            1/2*log(t) + log(x) + 2*log(y)

            sage: f.simplify_log('ratios')
            log(sqrt(t)*x*y^2)

            sage: f.simplify_log()
            log(x*y^2) + 1/2*log(t)

        To contract terms with no coefficient (more precisely, with
        coefficients 1 and -1) use option ``algorithm``::

            sage: f = log(x)+2*log(y)-log(t)
            sage: f.simplify_log('one')
            2*log(y) + log(x/t)

        ::

            sage: f = log(x)+log(y)-1/3*log((x+1))
            sage: f.simplify_log()
            log(x*y) - 1/3*log(x + 1)

            sage: f.simplify_log('ratios')
            log(x*y/(x + 1)^(1/3))

        `\pi` is irrational number, to contract logarithms in the following example
        we have to put ``algorithm`` to ``constants`` or ``all``::

            sage: f = log(x)+log(y)-pi*log((x+1))
            sage: f.simplify_log('constants')
            log(x*y/(x + 1)^pi)

        x*log(9) is contracted only if ``algorithm`` is ``all``::

            sage: (x*log(9)).simplify_log()
            x*log(9)
            sage: (x*log(9)).simplify_log('all')
            log(9^x)

        TESTS:

        This shows that the issue at trac #7334 is fixed. Maxima intentionally
        keeps the expression inside the log factored::

            sage: log_expr = (log(sqrt(2)-1)+log(sqrt(2)+1))
            sage: log_expr.simplify_log('all')
            log((sqrt(2) + 1)*(sqrt(2) - 1))
            sage: _.simplify_rational()
            0
            sage: log_expr.simplify_full()   # applies both simplify_log and simplify_rational
            0

        We should use the current simplification domain rather than
        set it to 'real' explicitly (:trac:`12780`)::

            sage: f = sqrt(x^2)
            sage: f.simplify_log()
            sqrt(x^2)
            sage: from sage.calculus.calculus import maxima
            sage: maxima('domain: real;')
            real
            sage: f.simplify_log()
            abs(x)
            sage: maxima('domain: complex;')
            complex

        AUTHORS:

        - Robert Marik (11-2009)
        """
        from sage.calculus.calculus import maxima
        maxima.eval('savelogexpand:logexpand$ logexpand:false$')
        if algorithm is not None:
            maxima.eval('logconcoeffp:\'logconfun$')
        if algorithm == 'ratios':
            maxima.eval('logconfun(m):= featurep(m,integer) or ratnump(m)$')
        elif algorithm == 'one':
            maxima.eval('logconfun(m):= is(m=1) or is(m=-1)$')
        elif algorithm == 'constants':
            maxima.eval('logconfun(m):= constantp(m)$')
        elif algorithm == 'all':
            maxima.eval('logconfun(m):= true$')
        elif algorithm is not None:
            raise NotImplementedError, "unknown algorithm, see the help for available algorithms"
        res = self.parent()(self._maxima_().logcontract())
        if algorithm is not None:
            maxima.eval('logconcoeffp:false$')
        maxima.eval('logexpand:savelogexpand$')
        return res

    log_simplify = simplify_log

    @rename_keyword(deprecation=6094, method="algorithm")
    def expand_log(self,algorithm='products'):
        r"""
        Simplifies symbolic expression, which can contain logs.

        Expands logarithms of powers, logarithms of products and
        logarithms of quotients.  The option ``algorithm`` specifies
        which expression types should be expanded.

        INPUT:

        - ``self`` - expression to be simplified

        - ``algorithm`` - (default: 'products') optional, governs which
          expression is expanded. Possible values are

          - 'nothing' (no expansion),

          - 'powers' (log(a^r) is expanded),

          - 'products' (like 'powers' and also log(a*b) are expanded),

          - 'all' (all possible expansion).

          See also examples below.

        DETAILS: This uses the Maxima simplifier and sets
        ``logexpand`` option for this simplifier. From the Maxima
        documentation: "Logexpand:true causes log(a^b) to become
        b*log(a). If it is set to all, log(a*b) will also simplify to
        log(a)+log(b). If it is set to super, then log(a/b) will also
        simplify to log(a)-log(b) for rational numbers a/b,
        a#1. (log(1/b), for integer b, always simplifies.) If it is
        set to false, all of these simplifications will be turned
        off. "

        ALIAS: :meth:`log_expand` and :meth:`expand_log` are the same

        EXAMPLES:

        By default powers and products (and quotients) are expanded,
        but not quotients of integers::

            sage: (log(3/4*x^pi)).log_expand()
            pi*log(x) + log(3/4)

        To expand also log(3/4) use ``algorithm='all'``::

            sage: (log(3/4*x^pi)).log_expand('all')
            pi*log(x) - log(4) + log(3)

        To expand only the power use ``algorithm='powers'``.::

            sage: (log(x^6)).log_expand('powers')
            6*log(x)

        The expression ``log((3*x)^6)`` is not expanded with
        ``algorithm='powers'``, since it is converted into product
        first::

            sage: (log((3*x)^6)).log_expand('powers')
            log(729*x^6)

        This shows that the option ``algorithm`` from the previous call
        has no influence to future calls (we changed some default
        Maxima flag, and have to ensure that this flag has been
        restored)::

            sage: (log(3/4*x^pi)).log_expand()
            pi*log(x) + log(3/4)

            sage: (log(3/4*x^pi)).log_expand('all')
            pi*log(x) - log(4) + log(3)

            sage: (log(3/4*x^pi)).log_expand()
            pi*log(x) + log(3/4)

        TESTS:

        Most of these log expansions only make sense over the
        reals. So, we should set the Maxima ``domain`` variable to
        'real' before we call out to Maxima. When we return, however, we
        should set the ``domain`` back to what it was, rather than
        assuming that it was 'complex'. See :trac:`12780`::

            sage: from sage.calculus.calculus import maxima
            sage: maxima('domain: real;')
            real
            sage: x.expand_log()
            x
            sage: maxima('domain;')
            real
            sage: maxima('domain: complex;')
            complex

        AUTHORS:

        - Robert Marik (11-2009)
        """
        from sage.calculus.calculus import maxima
        original_domain = maxima.eval('domain')
        maxima.eval('domain: real$ savelogexpand:logexpand$')
        if algorithm == 'nothing':
            maxima_method='false'
        elif algorithm == 'powers':
            maxima_method='true'
        elif algorithm == 'products':
            maxima_method='all'
        elif algorithm == 'all':
            maxima_method='super'
        else:
            raise NotImplementedError, "unknown algorithm, see the help for available algorithms"
        maxima.eval('logexpand:%s'%maxima_method)
        res = self._maxima_()
        res = res.sage()
        # Set the domain back to what it was before expand_log() was called.
        maxima.eval('domain: %s$ logexpand:savelogexpand$' % original_domain)
        return res

    log_expand = expand_log


    def factor(self, dontfactor=[]):
        """
        Factors self, containing any number of variables or functions, into
        factors irreducible over the integers.

        INPUT:


        -  ``self`` - a symbolic expression

        -  ``dontfactor`` - list (default: []), a list of
           variables with respect to which factoring is not to occur.
           Factoring also will not take place with respect to any variables
           which are less important (using the variable ordering assumed for
           CRE form) than those on the 'dontfactor' list.


        EXAMPLES::

            sage: x,y,z = var('x, y, z')
            sage: (x^3-y^3).factor()
            (x^2 + x*y + y^2)*(x - y)
            sage: factor(-8*y - 4*x + z^2*(2*y + x))
            (x + 2*y)*(z + 2)*(z - 2)
            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: F = factor(f/(36*(1 + 2*y + y^2)), dontfactor=[x]); F
            1/36*(x^2 + 2*x + 1)*(y - 1)/(y + 1)

        If you are factoring a polynomial with rational coefficients (and
        dontfactor is empty) the factorization is done using Singular
        instead of Maxima, so the following is very fast instead of
        dreadfully slow::

            sage: var('x,y')
            (x, y)
            sage: (x^99 + y^99).factor()
            (x^60 + x^57*y^3 - x^51*y^9 - x^48*y^12 + x^42*y^18 + x^39*y^21 -
            x^33*y^27 - x^30*y^30 - x^27*y^33 + x^21*y^39 + x^18*y^42 -
            x^12*y^48 - x^9*y^51 + x^3*y^57 + y^60)*(x^20 + x^19*y -
            x^17*y^3 - x^16*y^4 + x^14*y^6 + x^13*y^7 - x^11*y^9 -
            x^10*y^10 - x^9*y^11 + x^7*y^13 + x^6*y^14 - x^4*y^16 -
            x^3*y^17 + x*y^19 + y^20)*(x^10 - x^9*y + x^8*y^2 - x^7*y^3 +
            x^6*y^4 - x^5*y^5 + x^4*y^6 - x^3*y^7 + x^2*y^8 - x*y^9 +
            y^10)*(x^6 - x^3*y^3 + y^6)*(x^2 - x*y + y^2)*(x + y)
        """
        from sage.calculus.calculus import symbolic_expression_from_maxima_string, symbolic_expression_from_string
        if len(dontfactor) > 0:
            m = self._maxima_()
            name = m.name()
            cmd = 'block([dontfactor:%s],factor(%s))'%(dontfactor, name)
            return symbolic_expression_from_maxima_string(cmd)
        else:
            try:
                from sage.rings.all import QQ
                f = self.polynomial(QQ)
                w = repr(f.factor())
                return symbolic_expression_from_string(w)
            except TypeError:
                pass
            return self.parent()(self._maxima_().factor())

    def factor_list(self, dontfactor=[]):
        """
        Returns a list of the factors of self, as computed by the
        factor command.

        INPUT:

        -  ``self`` - a symbolic expression

        -  ``dontfactor`` - see docs for :meth:`factor`

        .. note::

           If you already have a factored expression and just want to
           get at the individual factors, use the `_factor_list` method
           instead.

        EXAMPLES::

            sage: var('x, y, z')
            (x, y, z)
            sage: f = x^3-y^3
            sage: f.factor()
            (x^2 + x*y + y^2)*(x - y)

        Notice that the -1 factor is separated out::

            sage: f.factor_list()
            [(x^2 + x*y + y^2, 1), (x - y, 1)]

        We factor a fairly straightforward expression::

            sage: factor(-8*y - 4*x + z^2*(2*y + x)).factor_list()
            [(x + 2*y, 1), (z + 2, 1), (z - 2, 1)]

        A more complicated example::

            sage: var('x, u, v')
            (x, u, v)
            sage: f = expand((2*u*v^2-v^2-4*u^3)^2 * (-u)^3 * (x-sin(x))^3)
            sage: f.factor()
            -(4*u^3 - 2*u*v^2 + v^2)^2*u^3*(x - sin(x))^3
            sage: g = f.factor_list(); g
            [(4*u^3 - 2*u*v^2 + v^2, 2), (u, 3), (x - sin(x), 3), (-1, 1)]

        This function also works for quotients::

            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: g = f/(36*(1 + 2*y + y^2)); g
            1/36*(x^2*y^2 + 2*x*y^2 - x^2 + y^2 - 2*x - 1)/(y^2 + 2*y + 1)
            sage: g.factor(dontfactor=[x])
            1/36*(x^2 + 2*x + 1)*(y - 1)/(y + 1)
            sage: g.factor_list(dontfactor=[x])
            [(x^2 + 2*x + 1, 1), (y + 1, -1), (y - 1, 1), (1/36, 1)]

        This example also illustrates that the exponents do not have to be
        integers::

            sage: f = x^(2*sin(x)) * (x-1)^(sqrt(2)*x); f
            (x - 1)^(sqrt(2)*x)*x^(2*sin(x))
            sage: f.factor_list()
            [(x - 1, sqrt(2)*x), (x, 2*sin(x))]
        """
        return self.factor(dontfactor=dontfactor)._factor_list()

    def _factor_list(self):
        r"""
        Turn an expression already in factored form into a list of (prime,
        power) pairs.

        This is used, e.g., internally by the :meth:`factor_list`
        command.

        EXAMPLES::

            sage: g = factor(x^3 - 1); g
            (x^2 + x + 1)*(x - 1)
            sage: v = g._factor_list(); v
            [(x^2 + x + 1, 1), (x - 1, 1)]
            sage: type(v)
            <type 'list'>
        """
        op = self.operator()
        if op is operator.mul:
            return sum([f._factor_list() for f in self.operands()], [])
        elif op is operator.pow:
            return [tuple(self.operands())]
        else:
            return [(self, 1)]

    ###################################################################
    # Units
    ###################################################################
    def convert(self, target=None):
        """
        Calls the convert function in the units package. For symbolic
        variables that are not units, this function just returns the
        variable.

        INPUT:

        - ``self`` -- the symbolic expression converting from
        - ``target`` -- (default None) the symbolic expression
          converting to

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: units.length.foot.convert()
            381/1250*meter
            sage: units.mass.kilogram.convert(units.mass.pound)
            100000000/45359237*pound

        We don't get anything new by converting an ordinary symbolic variable::

            sage: a = var('a')
            sage: a - a.convert()
            0

        Raises ValueError if self and target are not convertible::

            sage: units.mass.kilogram.convert(units.length.foot)
            Traceback (most recent call last):
            ...
            ValueError: Incompatible units
            sage: (units.length.meter^2).convert(units.length.foot)
            Traceback (most recent call last):
            ...
            ValueError: Incompatible units

        Recognizes derived unit relationships to base units and other
        derived units::

            sage: (units.length.foot/units.time.second^2).convert(units.acceleration.galileo)
            762/25*galileo
            sage: (units.mass.kilogram*units.length.meter/units.time.second^2).convert(units.force.newton)
            newton
            sage: (units.length.foot^3).convert(units.area.acre*units.length.inch)
            1/3630*(acre*inch)
            sage: (units.charge.coulomb).convert(units.current.ampere*units.time.second)
            (ampere*second)
            sage: (units.pressure.pascal*units.si_prefixes.kilo).convert(units.pressure.pounds_per_square_inch)
            1290320000000/8896443230521*pounds_per_square_inch

        For decimal answers multiply by 1.0::

            sage: (units.pressure.pascal*units.si_prefixes.kilo).convert(units.pressure.pounds_per_square_inch)*1.0
            0.145037737730209*pounds_per_square_inch

        Converting temperatures works as well::

            sage: s = 68*units.temperature.fahrenheit
            sage: s.convert(units.temperature.celsius)
            20*celsius
            sage: s.convert()
            293.150000000000*kelvin

        Trying to multiply temperatures by another unit then converting
        raises a ValueError::

            sage: wrong = 50*units.temperature.celsius*units.length.foot
            sage: wrong.convert()
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert
        """
        import units
        return units.convert(self, target)

    ###################################################################
    # solve
    ###################################################################
    def roots(self, x=None, explicit_solutions=True, multiplicities=True, ring=None):
        r"""
        Returns roots of ``self`` that can be found exactly,
        possibly with multiplicities.  Not all roots are guaranteed to
        be found.

        .. warning::

           This is *not* a numerical solver - use ``find_root`` to
           solve for self == 0 numerically on an interval.

        INPUT:

        - ``x`` - variable to view the function in terms of
          (use default variable if not given)

        - ``explicit_solutions`` - bool (default True); require that
          roots be explicit rather than implicit

        - ``multiplicities`` - bool (default True); when True, return
          multiplicities

        - ``ring`` - a ring (default None): if not None, convert
          self to a polynomial over ring and find roots over ring

        OUTPUT:

        A list of pairs ``(root, multiplicity)`` or list of roots.

        If there are infinitely many roots, e.g., a function like
        `\sin(x)`, only one is returned.

        EXAMPLES::

            sage: var('x, a')
            (x, a)

        A simple example::

            sage: ((x^2-1)^2).roots()
            [(-1, 2), (1, 2)]
            sage: ((x^2-1)^2).roots(multiplicities=False)
            [-1, 1]

        A complicated example::

            sage: f = expand((x^2 - 1)^3*(x^2 + 1)*(x-a)); f
            -a*x^8 + x^9 + 2*a*x^6 - 2*x^7 - 2*a*x^2 + 2*x^3 + a - x

        The default variable is `a`, since it is the first in
        alphabetical order::

            sage: f.roots()
            [(x, 1)]

        As a polynomial in `a`, `x` is indeed a root::

            sage: f.poly(a)
            x^9 - 2*x^7 + 2*x^3 - (x^8 - 2*x^6 + 2*x^2 - 1)*a - x
            sage: f(a=x)
            0

        The roots in terms of `x` are what we expect::

            sage: f.roots(x)
            [(a, 1), (-I, 1), (I, 1), (1, 3), (-1, 3)]

        Only one root of `\sin(x) = 0` is given::

            sage: f = sin(x)
            sage: f.roots(x)
            [(0, 1)]

        .. note::

            It is possible to solve a greater variety of equations
            using ``solve()`` and the keyword ``to_poly_solve``,
            but only at the price of possibly encountering
            approximate solutions.  See documentation for f.solve
            for more details.

        We derive the roots of a general quadratic polynomial::

            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: (a*x^2 + b*x + c).roots(x)
            [(-1/2*(b + sqrt(b^2 - 4*a*c))/a, 1), (-1/2*(b - sqrt(b^2 - 4*a*c))/a, 1)]

        By default, all the roots are required to be explicit rather than
        implicit. To get implicit roots, pass ``explicit_solutions=False``
        to ``.roots()`` ::

            sage: var('x')
            x
            sage: f = x^(1/9) + (2^(8/9) - 2^(1/9))*(x - 1) - x^(8/9)
            sage: f.roots()
            Traceback (most recent call last):
            ...
            RuntimeError: no explicit roots found
            sage: f.roots(explicit_solutions=False)
            [((2^(8/9) + x^(8/9) - 2^(1/9) - x^(1/9))/(2^(8/9) - 2^(1/9)), 1)]

        Another example, but involving a degree 5 poly whose roots don't
        get computed explicitly::

            sage: f = x^5 + x^3 + 17*x + 1
            sage: f.roots()
            Traceback (most recent call last):
            ...
            RuntimeError: no explicit roots found
            sage: f.roots(explicit_solutions=False)
            [(x^5 + x^3 + 17*x + 1, 1)]
            sage: f.roots(explicit_solutions=False, multiplicities=False)
            [x^5 + x^3 + 17*x + 1]

        Now let's find some roots over different rings::

            sage: f.roots(ring=CC)
            [(-0.0588115223184..., 1), (-1.331099917875... - 1.52241655183732*I, 1), (-1.331099917875... + 1.52241655183732*I, 1), (1.36050567903502 - 1.51880872209965*I, 1), (1.36050567903502 + 1.51880872209965*I, 1)]
            sage: (2.5*f).roots(ring=RR)
            [(-0.058811522318449..., 1)]
            sage: f.roots(ring=CC, multiplicities=False)
            [-0.05881152231844..., -1.331099917875... - 1.52241655183732*I, -1.331099917875... + 1.52241655183732*I, 1.36050567903502 - 1.51880872209965*I, 1.36050567903502 + 1.51880872209965*I]
            sage: f.roots(ring=QQ)
            []
            sage: f.roots(ring=QQbar, multiplicities=False)
            [-0.05881152231844944?, -1.331099917875796? - 1.522416551837318?*I, -1.331099917875796? + 1.522416551837318?*I, 1.360505679035020? - 1.518808722099650?*I, 1.360505679035020? + 1.518808722099650?*I]

        Root finding over finite fields::

            sage: f.roots(ring=GF(7^2, 'a'))
            [(3, 1), (4*a + 6, 2), (3*a + 3, 2)]

        TESTS::

            sage: (sqrt(3) * f).roots(ring=QQ)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(3) to a rational

        Check if #9538 is fixed::

            sage: var('f6,f5,f4,x')
            (f6, f5, f4, x)
            sage: e=15*f6*x^2 + 5*f5*x + f4
            sage: res = e.roots(x); res
            [(-1/30*(5*f5 + sqrt(25*f5^2 - 60*f4*f6))/f6, 1), (-1/30*(5*f5 - sqrt(25*f5^2 - 60*f4*f6))/f6, 1)]
            sage: e.subs(x=res[0][0]).is_zero()
            True
        """
        if x is None:
            x = self.default_variable()
        if ring is not None:
            p = self.polynomial(ring)
            return p.roots(ring=ring, multiplicities=multiplicities)

        S, mul = self.solve(x, multiplicities=True, explicit_solutions=explicit_solutions)
        if len(mul) == 0 and explicit_solutions:
            raise RuntimeError, "no explicit roots found"
        else:
            rt_muls = [(S[i].rhs(), mul[i]) for i in range(len(mul))]
        if multiplicities:
            return rt_muls
        else:
            return [ rt for rt, mul in rt_muls ]

    def solve(self, x, multiplicities=False, solution_dict=False, explicit_solutions=False, to_poly_solve=False):
        r"""
        Analytically solve the equation ``self == 0`` or an univarite
        inequality for the variable `x`.

        .. warning::

           This is not a numerical solver - use ``find_root`` to solve
           for self == 0 numerically on an interval.

        INPUT:


        -  ``x`` - variable to solve for

        -  ``multiplicities`` - bool (default: False); if True,
           return corresponding multiplicities.  This keyword is
           incompatible with ``to_poly_solve=True`` and does not make
           any sense when solving an inequality.

        -  ``solution_dict`` - bool (default: False); if True or non-zero,
           return a list of dictionaries containing solutions. Not used
           when solving an inequality.

        -  ``explicit_solutions`` - bool (default: False); require that
           all roots be explicit rather than implicit. Not used
           when solving an inequality.

        -  ``to_poly_solve`` - bool (default: False) or string; use
           Maxima's ``to_poly_solver`` package to search for more possible
           solutions, but possibly encounter approximate solutions.
           This keyword is incompatible with ``multiplicities=True``
           and is not used when solving an inequality. Setting ``to_poly_solve``
           to 'force' (string) omits Maxima's solve command (usefull when
           some solution of trigonometric equations are lost).

        EXAMPLES::

            sage: z = var('z')
            sage: (z^5 - 1).solve(z)
            [z == e^(2/5*I*pi), z == e^(4/5*I*pi), z == e^(-4/5*I*pi), z == e^(-2/5*I*pi), z == 1]

            sage: solve((z^3-1)^3, z, multiplicities=True)
            ([z == 1/2*I*sqrt(3) - 1/2, z == -1/2*I*sqrt(3) - 1/2, z == 1], [3, 3, 3])

        A simple example to show use of the keyword
        ``multiplicities``::

            sage: ((x^2-1)^2).solve(x)
            [x == -1, x == 1]
            sage: ((x^2-1)^2).solve(x,multiplicities=True)
            ([x == -1, x == 1], [2, 2])
            sage: ((x^2-1)^2).solve(x,multiplicities=True,to_poly_solve=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: to_poly_solve does not return multiplicities

        Here is how the ``explicit_solutions`` keyword functions::

            sage: solve(sin(x)==x,x)
            [x == sin(x)]
            sage: solve(sin(x)==x,x,explicit_solutions=True)
            []
            sage: solve(x*sin(x)==x^2,x)
            [x == 0, x == sin(x)]
            sage: solve(x*sin(x)==x^2,x,explicit_solutions=True)
            [x == 0]

        The following examples show use of the keyword ``to_poly_solve``::

            sage: solve(abs(1-abs(1-x)) == 10, x)
            [abs(abs(x - 1) - 1) == 10]
            sage: solve(abs(1-abs(1-x)) == 10, x, to_poly_solve=True)
            [x == -10, x == 12]

            sage: var('Q')
            Q
            sage: solve(Q*sqrt(Q^2 + 2) - 1, Q)
            [Q == 1/sqrt(Q^2 + 2)]
            sage: solve(Q*sqrt(Q^2 + 2) - 1, Q, to_poly_solve=True)
            [Q == 1/sqrt(-sqrt(2) + 1), Q == 1/sqrt(sqrt(2) + 1)]

        In some cases there may be infinitely many solutions indexed
        by a dummy variable.  If it begins with ``z``, it is implicitly
        assumed to be an integer, a real if with ``r``, and so on::

            sage: solve( sin(x)==cos(x), x, to_poly_solve=True)
            [x == 1/4*pi + pi*z...]

        An effort is made to only return solutions that satisfy the current assumptions::

            sage: solve(x^2==4, x)
            [x == -2, x == 2]
            sage: assume(x<0)
            sage: solve(x^2==4, x)
            [x == -2]
            sage: solve((x^2-4)^2 == 0, x, multiplicities=True)
            ([x == -2], [2])
            sage: solve(x^2==2, x)
            [x == -sqrt(2)]
            sage: assume(x, 'rational')
            sage: solve(x^2 == 2, x)
            []
            sage: solve(x^2==2-z, x)
            [x == -sqrt(-z + 2)]
            sage: solve((x-z)^2==2, x)
            [x == z - sqrt(2), x == z + sqrt(2)]

        There is still room for improvement::

            sage: assume(x, 'integer')
            sage: assume(z, 'integer')
            sage: solve((x-z)^2==2, x)
            [x == z - sqrt(2), x == z + sqrt(2)]

            sage: forget()

        In some cases it may be worthwhile to directly use to_poly_solve,
        if one suspects some answers are being missed::

            sage: solve(cos(x)==0,x)
            [x == 1/2*pi]
            sage: solve(cos(x)==0,x,to_poly_solve=True)
            [x == 1/2*pi]
            sage: from sage.calculus.calculus import maxima
            sage: sol = maxima(cos(x)==0).to_poly_solve(x)
            sage: sol.sage()
            [[x == -1/2*pi + 2*pi*z...], [x == 1/2*pi + 2*pi*z...]]

        If a returned unsolved expression has a denominator, but the
        original one did not, this may also be true::

            sage: solve(cos(x) * sin(x) == 1/2, x, to_poly_solve=True)
            [sin(x) == 1/2/cos(x)]
            sage: from sage.calculus.calculus import maxima
            sage: sol = maxima(cos(x) * sin(x) == 1/2).to_poly_solve(x)
            sage: sol.sage()
            [[x == 1/4*pi + pi*z...]]

        Some basic inequalities can be also solved::

            sage: x,y=var('x,y'); (ln(x)-ln(y)>0).solve(x)
            [[log(x) - log(y) > 0]]

        ::

            sage: x,y=var('x,y'); (ln(x)>ln(y)).solve(x) # not tested - output depends on system
            [[0 < y, y < x, 0 < x]]
            [[y < x, 0 < y]]

        TESTS:

        :trac:`7325` (solving inequalities)::

            sage: (x^2>1).solve(x)
            [[x < -1], [x > 1]]

        Catch error message from Maxima::

            sage: solve(acot(x),x)
            []

        ::

            sage: solve(acot(x),x,to_poly_solve=True)
            []

        :trac:`7491` fixed::

            sage: y=var('y')
            sage: solve(y==y,y)
            [y == r1]
            sage: solve(y==y,y,multiplicities=True)
            ([y == r1], [])

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: GenericDeclaration(x, 'rational').assume()
            sage: solve(x^2 == 2, x)
            []
            sage: forget()

        :trac:`8390` fixed::

            sage: solve(sin(x)==1/2,x)
            [x == 1/6*pi]

        ::

            sage: solve(sin(x)==1/2,x,to_poly_solve=True)
            [x == 1/6*pi]

        ::

            sage: solve(sin(x)==1/2,x,to_poly_solve='force')
            [x == 5/6*pi + 2*pi*z..., x == 1/6*pi + 2*pi*z...]

        :trac:`11618` fixed::

            sage: g(x)=0
            sage: solve(g(x)==0,x,solution_dict=True)
            [{x: r1}]

        :trac:`13286` fixed::

            sage: solve([x-4], [x])
            [x == 4]

        :trac:`13645`: fixed::

            sage: x.solve((1,2))
            Traceback (most recent call last):
            ...
            TypeError: 1 is not a valid variable.
        """
        import operator
        cdef Expression ex
        if is_a_relational(self._gobj):
            if self.operator() is not operator.eq:
                from sage.symbolic.relation import solve_ineq
                try:
                    return(solve_ineq(self)) # trying solve_ineq_univar
                except StandardError:
                    pass
                try:
                    return(solve_ineq([self])) # trying solve_ineq_fourier
                except StandardError:
                    raise NotImplementedError, "solving only implemented for equalities and few special inequalities, see solve_ineq"
            ex = self
        else:
            ex = (self == 0)

        if multiplicities and to_poly_solve:
            raise NotImplementedError, "to_poly_solve does not return multiplicities"

        # Take care of cases like solve([x^2-1], [x]) for consistency with
        # multiple variable input in sage.symbolic.relation.solve().
        # There *should* be only one variable in the list, since it is
        # passed from sage.symbolic.relation.solve() and multiple variables
        # there don't call this function.
        if isinstance(x, (list, tuple)):
            x = x[0]

        if x is None:
            v = ex.variables()
            if len(v) == 0:
                if multiplicities:
                    return [], []
                else:
                    return []
            x = v[0]

        if not isinstance(x, Expression):
            raise TypeError("%s is not a valid variable."%repr(x))

        m = ex._maxima_()
        P = m.parent()
        if explicit_solutions:
            P.eval('solveexplicit: true') # switches Maxima to looking for only explicit solutions
        try:
            if to_poly_solve != 'force':
                s = m.solve(x).str()
            else: # omit Maxima's solve command
                s = str([])
        except TypeError, mess: # if Maxima's solve has an error, we catch it
            if "Error executing code in Maxima" in str(mess):
                s = str([])
            else:
                raise
        if explicit_solutions:
            P.eval('solveexplicit: false') # switches Maxima back to default

        if s == 'all':
            if solution_dict:
                ans = [ {x: self.parent().var('r1')} ]
            else:
                ans = [x == self.parent().var('r1')]
            if multiplicities:
                return ans,[]
            else:
                return ans

        from sage.symbolic.relation import string_to_list_of_solutions

        X = string_to_list_of_solutions(s) # our initial list of solutions

        if multiplicities: # to_poly_solve does not return multiplicities, so in this case we end here
            if len(X) == 0:
                return X, []
            else:
                ret_multiplicities = [int(e) for e in str(P.get('multiplicities'))[1:-1].split(',')]

        ########################################################
        # Maxima's to_poly_solver package converts difficult   #
        # equations to (quasi)-polynomial systems and uses     #
        # Maxima's algsys function to try to solve them.       #
        # This allows a much larger range of solved equations, #
        # but also allows for the possibility of approximate   #
        # solutions being returned.                            #
        ########################################################
        if to_poly_solve and not multiplicities:
            if len(X)==0: # if Maxima's solve gave no solutions, only try it
                try:
                    s = m.to_poly_solve(x)
                    T = string_to_list_of_solutions(repr(s))
                    X = [t[0] for t in T]
                except StandardError: # if that gives an error, stick with no solutions
                    X = []

            for eq in X:
                # If the RHS of one solution also has the variable, or if
                # the LHS is not the variable, try another way to get solutions
                if repr(x) in map(repr, eq.rhs().variables()) or \
                        repr(x) in repr(eq.lhs()):
                    try:
                        # try to solve it using to_poly_solve
                        Y = eq._maxima_().to_poly_solve(x).sage()
                        X.remove(eq)
                        # replace with the new solutions
                        X.extend([y[0] for y in Y])
                    except TypeError, mess:
                        if "Error executing code in Maxima" in str(mess) or \
                                "unable to make sense of Maxima expression" in\
                                str(mess):
                            if explicit_solutions:
                                X.remove(eq) # this removes an implicit solution
                            else:
                                pass # we keep this implicit solution
                        else:
                            raise

        # make sure all the assumptions are satisfied
        from sage.symbolic.assumptions import assumptions
        to_check = assumptions()
        if to_check:
            for ix, soln in reversed(list(enumerate(X))):
                if soln.lhs().is_symbol():
                    if any([a.contradicts(soln) for a in to_check]):
                        del X[ix]
                        if multiplicities:
                            del ret_multiplicities[ix]
                        continue

        # if solution_dict is True:
        # Relaxed form suggested by Mike Hansen (#8553):
        if solution_dict:
            X=[dict([[sol.left(),sol.right()]]) for sol in X]

        if multiplicities:
            return X, ret_multiplicities
        else:
            return X

    def find_root(self, a, b, var=None, xtol=10e-13, rtol=4.5e-16, maxiter=100, full_output=False):
        """
        Numerically find a root of self on the closed interval [a,b] (or
        [b,a]) if possible, where self is a function in the one variable.

        INPUT:

        -  ``a, b`` - endpoints of the interval

        -  ``var`` - optional variable

        -  ``xtol, rtol`` - the routine converges when a root
           is known to lie within xtol of the value return. Should be = 0. The
           routine modifies this to take into account the relative precision
           of doubles.

        -  ``maxiter`` - integer; if convergence is not
           achieved in maxiter iterations, an error is raised. Must be = 0.

        -  ``full_output`` - bool (default: False), if True,
           also return object that contains information about convergence.


        EXAMPLES:

        Note that in this example both f(-2) and f(3) are positive,
        yet we still find a root in that interval::

            sage: f = x^2 - 1
            sage: f.find_root(-2, 3)
            1.0
            sage: f.find_root(-2, 3, x)
            1.0
            sage: z, result = f.find_root(-2, 3, full_output=True)
            sage: result.converged
            True
            sage: result.flag
            'converged'
            sage: result.function_calls
            11
            sage: result.iterations
            10
            sage: result.root
            1.0

        More examples::

            sage: (sin(x) + exp(x)).find_root(-10, 10)
            -0.588532743981862...
            sage: sin(x).find_root(-1,1)
            0.0
            sage: (1/tan(x)).find_root(3,3.5)
            3.1415926535...

        An example with a square root::

            sage: f = 1 + x + sqrt(x+2); f.find_root(-2,10)
            -1.618033988749895

        Some examples that Ted Kosan came up with::

            sage: t = var('t')
            sage: v = 0.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))
            sage: v.find_root(0, 0.002)
            0.001540327067911417...

        With this expression, we can see there is a
        zero very close to the origin::

            sage: a = .004*(8*e^(-(300*t)) - 8*e^(-(1200*t)))*(720000*e^(-(300*t)) - 11520000*e^(-(1200*t))) +.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))^2
            sage: show(plot(a, 0, .002), xmin=0, xmax=.002)

        It is easy to approximate with ``find_root``::

            sage: a.find_root(0,0.002)
            0.0004110514049349...

        Using solve takes more effort, and even then gives
        only a solution with free (integer) variables::

            sage: a.solve(t)
            []
            sage: b = a.simplify_radical(); b
            -23040*(-2.0*e^(1800*t) + 25.0*e^(900*t) - 32.0)*e^(-2400*t)
            sage: b.solve(t)
            []
            sage: b.solve(t, to_poly_solve=True)
            [t == 1/450*I*pi*z... + 1/900*log(3/4*sqrt(41) + 25/4), t == 1/450*I*pi*z... + 1/900*log(-3/4*sqrt(41) + 25/4)]
            sage: n(1/900*log(-3/4*sqrt(41) + 25/4))
            0.000411051404934985

        We illustrate that root finding is only implemented in one
        dimension::

            sage: x, y = var('x,y')
            sage: (x-y).find_root(-2,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented in 1 dimension.

        TESTS:

        Test the special case that failed for the first attempt to fix
        #3980::

            sage: t = var('t')
            sage: find_root(1/t - x,0,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented in 1 dimension.
        """
        if is_a_relational(self._gobj) and self.operator() is not operator.eq:
            raise ValueError, "Symbolic equation must be an equality."
        from sage.numerical.optimize import find_root
        if self.number_of_arguments() == 0:
            if bool(self == 0):
                return a
            else:
                raise RuntimeError, "no zero in the interval, since constant expression is not 0."
        elif self.number_of_arguments() == 1:
            f = self._fast_float_(self.default_variable())
            return find_root(f, a=a, b=b, xtol=xtol,
                             rtol=rtol,maxiter=maxiter,
                             full_output=full_output)
        else:
            raise NotImplementedError, "root finding currently only implemented in 1 dimension."

    def find_local_maximum(self, a, b, var=None, tol=1.48e-08, maxfun=500):
        r"""
        Numerically find a local maximum of the expression ``self``
        on the interval [a,b] (or [b,a]) along with the point at which the
        maximum is attained.

        See the documentation for
        :func:`find_local_minimum` for more details.

        EXAMPLES::

            sage: f = x*cos(x)
            sage: f.find_local_maximum(0,5)
            (0.5610963381910451, 0.8603335890...)
            sage: f.find_local_maximum(0,5, tol=0.1, maxfun=10)
            (0.561090323458081..., 0.857926501456...)
        """
        minval, x = (-self).find_local_minimum(a, b, var=var, tol=tol,
                                                     maxfun=maxfun)
        return -minval, x

    def find_local_minimum(self, a, b, var=None, tol=1.48e-08, maxfun=500):
        r"""
        Numerically find a local minimum of the expression ``self``
        on the interval [a,b] (or [b,a]) and the point at which it attains
        that minimum. Note that ``self`` must be a function of
        (at most) one variable.

        INPUT:

        -  ``var`` - variable (default: first variable in
           self)

        -  ``a,b`` - endpoints of interval on which to minimize
           self.

        -  ``tol`` - the convergence tolerance

        -  ``maxfun`` - maximum function evaluations


        OUTPUT:

        A tuple ``(minval, x)``, where

        - ``minval`` -- float. The minimum value that self takes on in
          the interval ``[a,b]``.

        - ``x`` -- float. The point at which self takes on the minimum
          value.

        EXAMPLES::

            sage: f = x*cos(x)
            sage: f.find_local_minimum(1, 5)
            (-3.288371395590..., 3.4256184695...)
            sage: f.find_local_minimum(1, 5, tol=1e-3)
            (-3.288371361890..., 3.4257507903...)
            sage: f.find_local_minimum(1, 5, tol=1e-2, maxfun=10)
            (-3.288370845983..., 3.4250840220...)
            sage: show(f.plot(0, 20))
            sage: f.find_local_minimum(1, 15)
            (-9.477294259479..., 9.5293344109...)

        ALGORITHM:

        Uses :func:`sage.numerical.optimize.find_local_minimum`.

        AUTHORS:

        - William Stein (2007-12-07)
        """
        from sage.numerical.optimize import find_local_minimum

        if var is None:
            var = self.default_variable()
        return find_local_minimum(self._fast_float_(var),
                                        a=a, b=b, tol=tol, maxfun=maxfun )

    find_maximum_on_interval = deprecated_function_alias(2607, find_local_maximum)
    find_minimum_on_interval = deprecated_function_alias(2607, find_local_minimum)

    ###################
    # Fast Evaluation #
    ###################
    def _fast_float_(self, *vars):
        """
        Returns an object which provides fast floating point
        evaluation of this symbolic expression.

        See :mod:`sage.ext.fast_eval` for more information.

        EXAMPLES::

            sage: f = sqrt(x+1)
            sage: ff = f._fast_float_('x')
            sage: ff(1.0)
            1.4142135623730951
            sage: type(_)
            <type 'float'>
        """
        from sage.symbolic.expression_conversions import fast_float
        return fast_float(self, *vars)

    def _fast_callable_(self, etb):
        """
        Given an ExpressionTreeBuilder *etb*, return an Expression representing
        this symbolic expression.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: x,y = var('x,y')
            sage: f = y+2*x^2
            sage: f._fast_callable_(etb)
            add(mul(ipow(v_0, 2), 2), v_1)
        """
        from sage.symbolic.expression_conversions import fast_callable
        return fast_callable(self, etb)

    def show(self):
        """
        Show this symbolic expression, i.e., typeset it nicely.

        EXAMPLES::

            sage: (x^2 + 1).show()
            x^{2}  + 1
        """
        from sage.misc.functional import _do_show
        return _do_show(self)

    def plot(self, *args, **kwds):
        """
        Plot a symbolic expression. All arguments are passed onto the standard plot command.

        EXAMPLES:

        This displays a straight line::

            sage: sin(2).plot((x,0,3))

        This draws a red oscillatory curve::

            sage: sin(x^2).plot((x,0,2*pi), rgbcolor=(1,0,0))

        Another plot using the variable theta::

            sage: var('theta')
            theta
            sage: (cos(theta) - erf(theta)).plot((theta,-2*pi,2*pi))

        A very thick green plot with a frame::

            sage: sin(x).plot((x,-4*pi, 4*pi), thickness=20, rgbcolor=(0,0.7,0)).show(frame=True)

        You can embed 2d plots in 3d space as follows::

            sage: plot(sin(x^2), (x,-pi, pi), thickness=2).plot3d(z = 1)

        A more complicated family::

            sage: G = sum([plot(sin(n*x), (x,-2*pi, 2*pi)).plot3d(z=n) for n in [0,0.1,..1]])
            sage: G.show(frame_aspect_ratio=[1,1,1/2])  # long time (5s on sage.math, 2012)

        A plot involving the floor function::

            sage: plot(1.0 - x * floor(1/x), (x,0.00001,1.0))

        Sage used to allow symbolic functions with "no arguments";
        this no longer works::

            sage: plot(2*sin, -4, 4)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and '<class 'sage.functions.trig.Function_sin'>'

        You should evaluate the function first::

            sage: plot(2*sin(x), -4, 4)

        TESTS::

            sage: f(x) = x*(1 - x)
            sage: plot(f,0,1)
        """
        from sage.symbolic.callable import is_CallableSymbolicExpression
        from sage.symbolic.ring import is_SymbolicVariable
        from sage.plot.plot import plot

        # see if the user passed a variable in.
        if kwds.has_key('param'):
            param = kwds['param']
        else:
            param = None
            for i, arg in enumerate(args):
                if is_SymbolicVariable(arg):
                    param = arg
                    args = args[:i] + args[i+1:]
                    break

        if param is None:
            if is_CallableSymbolicExpression(self):
                A = self.arguments()
                if len(A) == 0:
                    raise ValueError, "function has no input arguments"
                else:
                    param = A[0]

                f = self._plot_fast_callable(param)
            else:
                A = self.variables()
                if len(A) == 0:
                    #Here we handle the case where f is something
                    #like 2*sin, which has takes arguments which
                    #aren't explicitly given
                    n = self.number_of_arguments()
                    f = self._plot_fast_callable()
                else:
                    param = A[0]
                    try:
                        f = self._plot_fast_callable(param)
                    except NotImplementedError:
                        return self.function(param)
        else:
            try:
                f = self._plot_fast_callable(param)
            except NotImplementedError:
                return self.function(param)
        return plot(f, *args, **kwds)

    def _plot_fast_callable(self, *vars):
        """
        Internal function used for creating a fast callable version of this
        symbolic expression for plotting.

        EXAMPLES::

            sage: s = abs((1+I*x)^4); s
            abs((I*x + 1)^4)
            sage: s._plot_fast_callable(x)
            <sage.ext.interpreters.wrapper_py.Wrapper_py object at ...>
            sage: s._plot_fast_callable(x)(10)
            10201
            sage: abs((I*10+1)^4)
            10201
            sage: plot(s)
        """
        try:
            # First we try fast float.  However, this doesn't work on some
            # input where fast_callable works fine.
            return self._fast_float_(*vars)
        except (TypeError, NotImplementedError):
            # Now we try fast_callable as a fallback, since it works in some
            # cases when fast_float doesn't, e.g., when I is anywhere in the
            # expression fast_float doesn't work but fast_callable does in some
            # cases when the resulting expression is real.
            from sage.ext.fast_callable import fast_callable
            # I tried calling self._fast_callable_ but that's too complicated
            # of an interface so we just use the fast_callable function.
            return fast_callable(self, vars=vars)

    ############
    # Calculus #
    ############
    def sum(self, *args, **kwds):
        r"""
        Returns the symbolic sum
        `\sum_{v = a}^b self`

        with respect to the variable `v` with endpoints
        `a` and `b`.


        INPUT:

        -  ``v`` - a variable or variable name

        -  ``a`` - lower endpoint of the sum

        -  ``b`` - upper endpoint of the sum

        - ``algorithm`` - (default: 'maxima')  one of

                - 'maxima' - use Maxima (the default)

                - 'maple' - (optional) use Maple

                - 'mathematica' - (optional) use Mathematica

                - 'giac' - (optional) use Giac


        EXAMPLES::

            sage: k, n = var('k,n')
            sage: k.sum(k, 1, n).factor()
            1/2*(n + 1)*n

        ::

            sage: (1/k^4).sum(k, 1, oo)
            1/90*pi^4

        ::

            sage: (1/k^5).sum(k, 1, oo)
            zeta(5)

        A well known binomial identity::

            sage: assume(n>=0)
            sage: binomial(n,k).sum(k, 0, n)
            2^n

        And some truncations thereof::

            sage: binomial(n,k).sum(k,1,n)
            2^n - 1
            sage: binomial(n,k).sum(k,2,n)
            2^n - n - 1
            sage: binomial(n,k).sum(k,0,n-1)
            2^n - 1
            sage: binomial(n,k).sum(k,1,n-1)
            2^n - 2

        The binomial theorem::

            sage: x, y = var('x, y')
            sage: (binomial(n,k) * x^k * y^(n-k)).sum(k, 0, n)
            (x + y)^n

        ::

            sage: (k * binomial(n, k)).sum(k, 1, n)
            2^(n - 1)*n

        ::

            sage: ((-1)^k*binomial(n,k)).sum(k, 0, n)
            0

        ::

            sage: (2^(-k)/(k*(k+1))).sum(k, 1, oo)
            -log(2) + 1

        Summing a hypergeometric term::

            sage: (binomial(n, k) * factorial(k) / factorial(n+1+k)).sum(k, 0, n)
            1/2*sqrt(pi)/factorial(n + 1/2)

        We check a well known identity::

            sage: bool((k^3).sum(k, 1, n) == k.sum(k, 1, n)^2)
            True

        A geometric sum::

            sage: a, q = var('a, q')
            sage: (a*q^k).sum(k, 0, n)
            (a*q^(n + 1) - a)/(q - 1)

        The geometric series::

            sage: assume(abs(q) < 1)
            sage: (a*q^k).sum(k, 0, oo)
            -a/(q - 1)

        A divergent geometric series.  Don't forget
        to forget your assumptions::

            sage: forget()
            sage: assume(q > 1)
            sage: (a*q^k).sum(k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: Sum is divergent.

        This summation only Mathematica can perform::

            sage: (1/(1+k^2)).sum(k, -oo, oo, algorithm = 'mathematica')     # optional - mathematica
            pi*coth(pi)

        Use Giac to perform this summation::

            sage: (sum(1/(1+k^2), k, -oo, oo, algorithm = 'giac')).factor()       # optional - giac
            (e^(2*pi) + 1)*pi/((e^pi - 1)*(e^pi + 1))

        Use Maple as a backend for summation::

            sage: (binomial(n,k)*x^k).sum(k, 0, n, algorithm = 'maple')      # optional - maple
            (x + 1)^n

        Check that the sum in #10682 is done right::

            sage: sum(binomial(n,k)*k^2, k, 2, n)
            1/4*(n^2 + n)*2^n - n

        .. note::

           #. Sage can currently only understand a subset of the output of Maxima, Maple and
              Mathematica, so even if the chosen backend can perform the summation the
              result might not be convertable into a Sage expression.

        """
        from sage.calculus.calculus import symbolic_sum
        return symbolic_sum(self, *args, **kwds)

    def integral(self, *args, **kwds):
        """
        Compute the integral of self.  Please see
        :func:`sage.symbolic.integration.integral.integrate` for more details.

        EXAMPLES::

            sage: sin(x).integral(x,0,3)
            -cos(3) + 1
            sage: sin(x).integral(x)
            -cos(x)

        TESTS:

        We check that :trac:`12438` is resolved::

            sage: f(x) = x; f
            x |--> x
            sage: integral(f, x)
            x |--> 1/2*x^2
            sage: integral(f, x, 0, 1)
            1/2

            sage: f(x, y) = x + y
            sage: f
            (x, y) |--> x + y
            sage: integral(f, y, 0, 1)
            x |--> x + 1/2
            sage: integral(f, x, 0, 1)
            y |--> y + 1/2
            sage: _(3)
            7/2
            sage: var("z")
            z
            sage: integral(f, z, 0, 2)
            (x, y) |--> 2*x + 2*y
            sage: integral(f, z)
            (x, y) |--> (x + y)*z
        """
        from sage.symbolic.integration.integral import \
            integral, _normalize_integral_input
        from sage.symbolic.callable import \
            CallableSymbolicExpressionRing, is_CallableSymbolicExpressionRing
        R = self._parent
        if is_CallableSymbolicExpressionRing(R):
            f = ring.SR(self)
            f, v, a, b = _normalize_integral_input(f, *args)
            # Definite integral with respect to a positional variable.
            if a is not None and v in R.arguments():
                arguments = list(R.arguments())
                arguments.remove(v)
                if arguments:
                    arguments = tuple(arguments)
                    R = CallableSymbolicExpressionRing(arguments, check=False)
                else:   # all arguments are gone
                    R = ring.SR
            return R(integral(f, v, a, b, **kwds))
        return integral(self, *args, **kwds)

    integrate = integral

    def nintegral(self, *args, **kwds):
        """
        Compute the numerical integral of self.  Please see
        :obj:`sage.calculus.calculus.nintegral` for more details.

        EXAMPLES::

            sage: sin(x).nintegral(x,0,3)
            (1.989992496600..., 2.209335488557...e-14, 21, 0)
        """
        from sage.calculus.calculus import nintegral
        return nintegral(self, *args, **kwds)

    nintegrate = nintegral

    def minpoly(self, *args, **kwds):
        """
        Return the minimal polynomial of this symbolic expression.

        EXAMPLES::

            sage: golden_ratio.minpoly()
            x^2 - x - 1
        """
        try:
            obj = self.pyobject()
            return obj.minpoly()
        except AttributeError:
            pass
        except TypeError:
            pass
        from sage.calculus.calculus import minpoly
        return minpoly(self, *args, **kwds)

    def limit(self, *args, **kwds):
        """
        Return a symbolic limit.  See
        :obj:`sage.calculus.calculus.limit`

        EXAMPLES::

            sage: (sin(x)/x).limit(x=0)
            1
        """
        from sage.calculus.calculus import limit
        return limit(self, *args, **kwds)

    def laplace(self, t, s):
        """
        Return Laplace transform of self.  See
        :obj:`sage.calculus.calculus.laplace`

        EXAMPLES::

            sage: var('x,s,z')
            (x, s, z)
            sage: (z + exp(x)).laplace(x, s)
            z/s + 1/(s - 1)
        """
        from sage.calculus.calculus import laplace
        return laplace(self, t, s)

    def inverse_laplace(self, t, s):
        """
        Return inverse Laplace transform of self.  See
        :obj:`sage.calculus.calculus.inverse_laplace`

        EXAMPLES::

            sage: var('w, m')
            (w, m)
            sage: f = (1/(w^2+10)).inverse_laplace(w, m); f
            1/10*sqrt(10)*sin(sqrt(10)*m)
        """
        from sage.calculus.calculus import inverse_laplace
        return inverse_laplace(self, t, s)

    def add_to_both_sides(self, x):
        """
        Returns a relation obtained by adding *x* to both sides of
        this relation.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: eqn = x^2 + y^2 + z^2 <= 1
            sage: eqn.add_to_both_sides(-z^2)
            x^2 + y^2 <= -z^2 + 1
            sage: eqn.add_to_both_sides(I)
            x^2 + y^2 + z^2 + I <= (I + 1)
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self + x

    def subtract_from_both_sides(self, x):
        """
        Returns a relation obtained by subtracting *x* from both sides
        of this relation.

        EXAMPLES::

            sage: eqn = x*sin(x)*sqrt(3) + sqrt(2) > cos(sin(x))
            sage: eqn.subtract_from_both_sides(sqrt(2))
            sqrt(3)*x*sin(x) > -sqrt(2) + cos(sin(x))
            sage: eqn.subtract_from_both_sides(cos(sin(x)))
            sqrt(3)*x*sin(x) + sqrt(2) - cos(sin(x)) > 0
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self - x

    def multiply_both_sides(self, x, checksign=None):
        """
        Returns a relation obtained by multiplying both sides of this
        relation by *x*.

        .. note::

           The *checksign* keyword argument is currently ignored and
           is included for backward compatibility reasons only.

        EXAMPLES::

            sage: var('x,y'); f = x + 3 < y - 2
            (x, y)
            sage: f.multiply_both_sides(7)
            7*x + 21 < 7*y - 14
            sage: f.multiply_both_sides(-1/2)
            -1/2*x - 3/2 < -1/2*y + 1
            sage: f*(-2/3)
            -2/3*x - 2 < -2/3*y + 4/3
            sage: f*(-pi)
            -pi*(x + 3) < -pi*(y - 2)

        Since the direction of the inequality never changes when doing
        arithmetic with equations, you can multiply or divide the
        equation by a quantity with unknown sign::

            sage: f*(1+I)
            (I + 1)*x + 3*I + 3 < (I + 1)*y - 2*I - 2
            sage: f = sqrt(2) + x == y^3
            sage: f.multiply_both_sides(I)
            I*x + I*sqrt(2) == I*y^3
            sage: f.multiply_both_sides(-1)
            -x - sqrt(2) == -y^3

        Note that the direction of the following inequalities is
        not reversed::

            sage: (x^3 + 1 > 2*sqrt(3)) * (-1)
            -x^3 - 1 > -2*sqrt(3)
            sage: (x^3 + 1 >= 2*sqrt(3)) * (-1)
            -x^3 - 1 >= -2*sqrt(3)
            sage: (x^3 + 1 <= 2*sqrt(3)) * (-1)
            -x^3 - 1 <= -2*sqrt(3)
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self * x

    def divide_both_sides(self, x, checksign=None):
        """
        Returns a relation obtained by dividing both sides of this
        relation by *x*.

        .. note::

           The *checksign* keyword argument is currently ignored and
           is included for backward compatibility reasons only.

        EXAMPLES::

            sage: theta = var('theta')
            sage: eqn =   (x^3 + theta < sin(x*theta))
            sage: eqn.divide_both_sides(theta, checksign=False)
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: eqn.divide_both_sides(theta)
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: eqn/theta
            (x^3 + theta)/theta < sin(theta*x)/theta
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self / x


    # Functions to add later, maybe.  These were in Ginac mainly
    # implemented using a lot from cln, and I had to mostly delete
    # their implementations.   They are pretty specialized for
    # physics apps, maybe.
    # This doesn't work / isn't implemented yet / just segfaults.
    #def Li(self, x):
    #    """
    #    """
    #    cdef Expression nexp = self.coerce_in(x)
    #    return new_Expression_from_GEx(self._parent, g_Li(self._gobj, nexp._gobj))
    #def Li2(self):
    #    return new_Expression_from_GEx(self._parent, g_Li2(self._gobj))
    #def G(self, Expression y):
    #    return new_Expression_from_GEx(self._parent, g_G(self._gobj, y._gobj))
    #def G2(self, Expression s, Expression y):
    #    return new_Expression_from_GEx(self._parent, g_G2(self._gobj, s._gobj, y._gobj))
    #def SR(self, Expression p, Expression x):
    #return new_Expression_from_GEx(self._parent, g_SR(self._gobj, p._gobj, x._gobj))
    #def H(self, Expression x):
    #return new_Expression_from_GEx(self._parent, g_H(self._gobj, x._gobj))
    #def zeta2(self, Expression s):
    #    return new_Expression_from_GEx(self._parent, g_zeta2(self._gobj, s._gobj))
    #def zetaderiv(self, Expression x):
    #    return new_Expression_from_GEx(self._parent, g_zetaderiv(self._gobj, x._gobj))
    #def beta(self, Expression y):
    #    return new_Expression_from_GEx(self._parent, g_beta(self._gobj, y._gobj))
    #def psi(self):
    #    return new_Expression_from_GEx(self._parent, g_psi(self._gobj))
    #def psi2(self, Expression x):
    #    return new_Expression_from_GEx(self._parent, g_psi2(self._gobj, x._gobj))


cdef dict dynamic_class_cache = {}
cdef get_dynamic_class_for_function(unsigned serial):
    r"""
    Create a dynamic class corresponding to the function with given
    ``serial`` that includes dynamic methods defined by the function.

    Dynamic methods can be defined in a subclass ``EvaluationMethods`` in
    the function body. These will be available in symbolic expressions
    representing evaluations of the said function on some arguments.

    EXAMPLES::

        sage: from sage.symbolic.function import BuiltinFunction
        sage: class TFunc(BuiltinFunction):
        ....:     def __init__(self):
        ....:         BuiltinFunction.__init__(self, 'tfunc', nargs=1)
        ....:
        ....:     class EvaluationMethods:
        ....:         def argp1(self, x):
        ....:             '''
        ....:             Some documentation about a bogus function.
        ....:             '''
        ....:             return x+1
        ....:
        ....:         @property
        ....:         def foo(self):
        ....:             return 5
        ....:
        sage: tfunc = TFunc()
        sage: e = tfunc(x); e
        tfunc(x)
        sage: type(e)
        <class '__main__.Expression_with_dynamic_methods'>
        sage: e.argp1()
        x + 1
        sage: e.foo
        5
        sage: x.argp1()
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.symbolic.expression.Expression' object has no
        attribute 'argp1'
        sage: t = (e + 1).op[0]; t
        tfunc(x)
        sage: t
        tfunc(x)
        sage: type(t)
        <class '__main__.Expression_with_dynamic_methods'>
        sage: t.argp1()
        x + 1
        sage: import sagenb.misc.support as s
        sage: s.completions('t.argp', globals(), system='python')
        ['t.argp1']
        sage: t.argp1.__doc__.strip()
        'Some documentation about a bogus function.'

    Now with two arguments::

        sage: class TFunc2(BuiltinFunction):
        ....:     def __init__(self):
        ....:         BuiltinFunction.__init__(self, 'tfunc', nargs=2)
        ....:
        ....:     class EvaluationMethods:
        ....:         def argsum(self, x, y):
        ....:             return x + y
        ....:
        sage: tfunc2 = TFunc2()
        sage: e = tfunc2(x, 1)
        sage: e.argsum()
        x + 1
    """
    cls = dynamic_class_cache.get(serial)
    if cls is None:
        # if operator is a special function defined by us
        # find the python equivalent and return it
        func_class = get_sfunction_from_serial(serial)
        eval_methods = getattr(func_class, 'EvaluationMethods', None)
        if eval_methods is not None:
            # callable methods need to be wrapped to extract the operands
            # and pass them as arguments
            from sage.symbolic.function_factory import eval_on_operands
            from sage.structure.parent import getattr_from_other_class
            for name in dir(eval_methods):
                m = getattr_from_other_class(func_class, eval_methods, name)
                if callable(m):
                    setattr(eval_methods, name, eval_on_operands(m))
            cls = dynamic_class('Expression_with_dynamic_methods',
                    (eval_methods, Expression))
        else:
            cls = Expression

        dynamic_class_cache[serial] = cls

    return cls

cdef Expression new_Expression_from_GEx(parent, GEx juice):
    cdef Expression nex
    cdef unsigned serial
    if is_exactly_a_function(juice):
        # if the function defines any dynamic methods these are made
        # available through a dynamic class
        cls = get_dynamic_class_for_function(ex_to_function(juice).get_serial())
    else:
        cls = Expression

    nex = <Expression>PY_NEW(cls)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = parent
    return nex

cdef Expression new_Expression_from_pyobject(parent, x):
    cdef GEx exp
    GEx_construct_pyobject(exp, x)
    return new_Expression_from_GEx(parent, exp)

cdef class ExpressionIterator:
    cdef Expression _ex
    cdef int _ind
    cdef int _len
    def __iter__(self):
        """
        Return this iterator object itself.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: i = (x+y).iterator()
            sage: iter(i) is i
            True
        """
        return self

    def __next__(self):
        """
        Return the next component of the expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: i = (x+y).iterator()
            sage: i.next()
            x
        """
        cdef GEx ex
        if self._ind == self._len:
            raise StopIteration
        ex = self._ex._gobj.op(self._ind)
        self._ind+=1
        return new_Expression_from_GEx(self._ex._parent, ex)

cdef inline ExpressionIterator new_ExpIter_from_Expression(Expression ex):
    """
    Construct a new iterator over a symbolic expression.

    EXAMPLES::

        sage: x,y,z = var('x,y,z')
        sage: i = (x+y).iterator() #indirect doctest
    """
    # The const_iterator in GiNaC just keeps an integer index to the current
    # subexpression. We do the same here, to avoid the trouble of having to
    # mess with C++ class constructors/destructors.
    cdef ExpressionIterator m = <ExpressionIterator>PY_NEW(ExpressionIterator)
    m._ex = ex
    m._ind = 0
    m._len  = ex._gobj.nops()
    return m


cdef operators compatible_relation(operators lop, operators rop) except <operators>-1:
    """
    TESTS::

        sage: var('a,b,x,y')
        (a, b, x, y)
        sage: (x < a) + (y <= b)     # indirect doctest
        x + y < a + b
        sage: (x >= 4) * (y > 7)
        x*y > 28
    """
    if lop == rop:
        return lop
    elif lop == not_equal or rop == not_equal:
        raise TypeError, "incompatible relations"
    elif lop == equal:
       return rop
    elif rop == equal:
       return lop
    elif lop in [less, less_or_equal] and rop in [less, less_or_equal]:
       return less
    elif lop in [greater, greater_or_equal] and rop in [greater, greater_or_equal]:
       return greater
    else:
        raise TypeError, "incompatible relations"
