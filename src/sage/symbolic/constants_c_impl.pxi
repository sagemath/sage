"""
The constant `e`
"""

# ****************************************************************************
#       Copyright (C) 2008      William Stein <wstein@gmail.com>
#       Copyright (C) 2008-2009 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009-2012 Mike Hansen <mhansen@gmail.com>
#       Copyright (C) 2010      David Kirkby
#       Copyright (C) 2012      Karl-Dieter Crisman
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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
            sage: import sympy                                                          # needs sympy
            sage: sympy.E == e  # indirect doctest                                      # needs sympy
            True

        TESTS::

            sage: t = e^a; t
            e^a
            sage: t^b
            (e^a)^b
            sage: SR(1).exp()
            e

        Testing that it works with matrices (see :issue:`4735`)::

            sage: m = matrix(QQ, 2, 2, [1,0,0,1])
            sage: e^m
            [e 0]
            [0 e]
        """
        from sage.symbolic.ring import SR

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

        This applies to the unit argument as well::

            sage: u = SR(1).exp()^a; u
            e^a
            sage: u.operator()
            exp
            sage: u.operands()
            [a]

        It also works with matrices (see :issue:`4735`)::

            sage: m = matrix(QQ, 2, 2, [1,0,0,1])
            sage: e^m
            [e 0]
            [0 e]
            sage: A = matrix(RDF, [[1,2],[3,4]])
            sage: e^A  # rel tol 5e-14
            [51.968956198705044  74.73656456700327]
            [112.10484685050491 164.07380304920997]
        """
        if isinstance(left, E):
            if isinstance(right, E):
                return exp_one.exp()
            try:
                return right.exp()
            except AttributeError:
                from sage.symbolic.ring import SR
                return SR(right).exp()
        else:
            from sage.symbolic.ring import SR
            return SR(left)**exp_one
