"""
Wrapper around Pynac's constants
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#                     2009 Mike Hansen <mhansen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
from __future__ import absolute_import

from sage.symbolic.expression cimport Expression, new_Expression_from_GEx
from sage.symbolic.ring import SR

from .ginac cimport *

cdef extern from "pynac/constant.h":
    pass

cdef class PynacConstant:
    def __cinit__(self, name, texname, domain_string):
        """
        Creates a constant in Pynac.

        EXAMPLES:

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f
            foo

        Note that this just creates a 'constant' object and not an
        expression.  If you want to work with this constant, you'll
        have to use the result of :meth:`expression`::

            sage: foo = f.expression()
            sage: foo + 2
            foo + 2
        """
        cdef unsigned domain
        if domain_string == 'complex':
            domain = domain_complex
        elif domain_string == 'real':
            domain = domain_real
        elif domain_string == 'positive':
            domain = domain_positive
        else:
            raise ValueError

        self._name = name

        #For the constants explicitly defined in constant.cpp in the Pynac
        #library, we use those symbols.
        if self._name == "pi":
            self.pointer = <GConstant *>&g_Pi
        elif self._name == "catalan":
            self.pointer = <GConstant *>&g_Catalan
        elif self._name == "euler_gamma":
            self.pointer = <GConstant *>&g_Euler
        else:
            GConstant_construct(&self.object, name, texname, domain)
            self.pointer = &self.object

    def __dealloc__(self):
        if self.pointer == & self.object:
            GConstant_destruct(&self.object)

    def serial(self):
        """
        Returns the underlying Pynac serial for this constant.

        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f.serial()  #random
            15
        """
        return int(self.pointer.get_serial())

    def name(self):
        """
        Returns the name of this constant.

        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f.name()
            'foo'
        """
        return self._name

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real'); f
            foo
        """
        return self.name()

    def expression(self):
        """
        Returns this constant as an Expression.

        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f + 2
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': '<type 'sage.symbolic.constants_c.PynacConstant'>' and 'Integer Ring'

            sage: foo = f.expression(); foo
            foo
            sage: foo + 2
            foo + 2
        """
        return new_Expression_from_GEx(SR, <GEx>(self.pointer[0]))

# keep exp(1) for fast access
# this is initialized in the constructor of the class E below to prevent
# circular imports while loading the library

# Note, the nearest IEEE 754 value of e is NOT the same as the correctly
# rounded decimal value of e.
# The numerical value of e                      = 2.71828182845904523536...
# Correctly rounded decimal number              = 2.7182818284590452
# Nearest IEEE 754 format number                = 2.7182818284590451
# Value returned on some or all AMD/Intel CPUs  = 2.7182818284590451
# On Sun Blade 1000 with SPARC processors       = 2.7182818284590455

cdef object exp_one

cdef class E(Expression):
    def __init__(self):
        r"""
        Dummy class to represent base of the natural logarithm.

        The base of the natural logarithm ``e`` is not a constant in GiNaC/Sage.
        It is represented by ``exp(1)``.

        This class provides a dummy object that behaves well under addition,
        multiplication, etc. and on exponentiation calls the function ``exp``.

        EXAMPLES:

        The constant defined at the top level is just ``exp(1)``::

            sage: e.operator()
            exp
            sage: e.operands()
            [1]

        Arithmetic works::

            sage: e + 2
            e + 2
            sage: 2 + e
            e + 2
            sage: 2*e
            2*e
            sage: e*2
            2*e
            sage: x*e
            x*e
            sage: var('a,b')
            (a, b)
            sage: t = e^(a+b); t
            e^(a + b)
            sage: t.operands()
            [a + b]

        Numeric evaluation, conversion to other systems, and pickling works
        as expected. Note that these are properties of the :func:`exp` function,
        not this class::

            sage: RR(e)
            2.71828182845905
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(e)
            2.7182818284590452353602874713526624977572470936999595749670
            sage: em = 1 + e^(1-e); em
            e^(-e + 1) + 1
            sage: R(em)
            1.1793740787340171819619895873183164984596816017589156131574
            sage: maxima(e).float()
            2.718281828459045
            sage: t = mathematica(e)               # optional - mathematica
            sage: t                                # optional - mathematica
            E
            sage: float(t)                         # optional - mathematica
            2.718281828459045...

            sage: loads(dumps(e))
            e

            sage: float(e)
            2.718281828459045...
            sage: e.__float__()
            2.718281828459045...
            sage: e._mpfr_(RealField(100))
            2.7182818284590452353602874714
            sage: e._real_double_(RDF)   # abs tol 5e-16
            2.718281828459045
            sage: import sympy
            sage: sympy.E == e # indirect doctest
            True

        TESTS::

            sage: t = e^a; t
            e^a
            sage: t^b
            (e^a)^b
            sage: SR(1).exp()
            e

        Testing that it works with matrices (see :trac:`4735`)::

            sage: m = matrix(QQ, 2, 2, [1,0,0,1])
            sage: e^m
            [e 0]
            [0 e]
        """
        global exp_one
        exp_one = SR.one().exp()
        Expression.__init__(self, SR, exp_one)

    def __pow__(left, right, dummy):
        """
        Call the `exp` function when taking powers of `e`.

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: t = e^a; t
            e^a
            sage: t.operator()
            exp
            sage: t.operands()
            [a]

        As opposed to::

            sage: u = SR(1).exp()^a; u
            e^a
            sage: u.operator()
            <built-in function pow>
            sage: u.operands()
            [e, a]

        It also works with matrices (see :trac:`4735`)::

            sage: m = matrix(QQ, 2, 2, [1,0,0,1])
            sage: e^m
            [e 0]
            [0 e]
            sage: A = matrix(RDF, [[1,2],[3,4]])
            sage: e^A  # rel tol 1e-15
            [51.968956198705044  74.73656456700327]
            [112.10484685050491 164.07380304920997]
        """
        if isinstance(left, E):
            if isinstance(right, E):
                return exp_one.exp()
            try:
                return right.exp()
            except AttributeError:
                return SR(right).exp()
        else:
            return SR(left)**exp_one

