# distutils: sources = sage/symbolic/ginac/add.cpp sage/symbolic/ginac/archive.cpp sage/symbolic/ginac/assume.cpp sage/symbolic/ginac/basic.cpp sage/symbolic/ginac/cmatcher.cpp sage/symbolic/ginac/constant.cpp sage/symbolic/ginac/context.cpp sage/symbolic/ginac/ex.cpp sage/symbolic/ginac/expair.cpp sage/symbolic/ginac/expairseq.cpp sage/symbolic/ginac/exprseq.cpp sage/symbolic/ginac/fderivative.cpp sage/symbolic/ginac/function.cpp sage/symbolic/ginac/function_info.cpp sage/symbolic/ginac/infinity.cpp sage/symbolic/ginac/infoflagbase.cpp sage/symbolic/ginac/inifcns.cpp sage/symbolic/ginac/inifcns_comb.cpp sage/symbolic/ginac/inifcns_gamma.cpp sage/symbolic/ginac/inifcns_hyperb.cpp sage/symbolic/ginac/inifcns_hyperg.cpp sage/symbolic/ginac/inifcns_nstdsums.cpp sage/symbolic/ginac/inifcns_orthopoly.cpp sage/symbolic/ginac/inifcns_trans.cpp sage/symbolic/ginac/inifcns_trig.cpp sage/symbolic/ginac/inifcns_zeta.cpp sage/symbolic/ginac/lst.cpp sage/symbolic/ginac/matrix.cpp sage/symbolic/ginac/mpoly-ginac.cpp sage/symbolic/ginac/mpoly-singular.cpp sage/symbolic/ginac/mpoly.cpp sage/symbolic/ginac/mul.cpp sage/symbolic/ginac/normal.cpp sage/symbolic/ginac/numeric.cpp sage/symbolic/ginac/operators.cpp sage/symbolic/ginac/order.cpp sage/symbolic/ginac/power.cpp sage/symbolic/ginac/print.cpp sage/symbolic/ginac/pseries.cpp sage/symbolic/ginac/py_funcs.cpp sage/symbolic/ginac/registrar.cpp sage/symbolic/ginac/relational.cpp sage/symbolic/ginac/remember.cpp sage/symbolic/ginac/sum.cpp sage/symbolic/ginac/symbol.cpp sage/symbolic/ginac/templates.cpp sage/symbolic/ginac/upoly-ginac.cpp sage/symbolic/ginac/useries.cpp sage/symbolic/ginac/utils.cpp sage/symbolic/ginac/wildcard.cpp
# distutils: language = c++
# distutils: libraries = flint gmp SINGULAR_LIBRARIES
# distutils: extra_compile_args = -std=c++11 SINGULAR_CFLAGS
# distutils: depends = ginac/add.h ginac/archive.h ginac/assertion.h ginac/assume.h ginac/basic.h ginac/class_info.h ginac/cmatcher.h ginac/compiler.h ginac/constant.h ginac/container.h ginac/context.h ginac/ex.h ginac/ex_utils.h ginac/expair.h ginac/expairseq.h ginac/exprseq.h ginac/extern_templates.h ginac/fderivative.h ginac/flags.h ginac/function.h ginac/ginac.h ginac/infinity.h ginac/infoflagbase.h ginac/inifcns.h ginac/lst.h ginac/matrix.h ginac/mpoly.h ginac/mul.h ginac/normal.h ginac/numeric.h ginac/operators.h ginac/order.h ginac/optional.hpp ginac/power.h ginac/print.h ginac/pseries.h ginac/ptr.h ginac/py_funcs.h ginac/pynac-config.h ginac/registrar.h ginac/relational.h ginac/remember.h ginac/sum.h ginac/symbol.h ginac/templates.h ginac/tostring.h ginac/upoly.h ginac/useries-flint.h ginac/useries.h ginac/utils.h ginac/wildcard.h
# distutils: include_dirs = SINGULAR_INCDIR
# pynac/basic.h includes
#   factory/factory.h    so this ^ is needed to find it
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

This can transform a true relation into a false one::

    sage: eqn = SR(-5) < SR(-3); eqn
    -5 < -3
    sage: bool(eqn)
    True
    sage: eqn^2
    25 < 9
    sage: bool(eqn^2)
    False

We can do arithmetic with relations::

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

We can even add together two relations, as long as the operators are
the same::

    sage: (x^3 + x <= x - 17)  + (-x <= x - 10)
    x^3 <= 2*x - 27

Here they are not::

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

Test Jacobian on Pynac expressions. (:issue:`5546`) ::

    sage: var('x,y')
    (x, y)
    sage: f = x + y
    sage: jacobian(f, [x,y])
    [1 1]

Test if matrices work (:issue:`5546`) ::

    sage: var('x,y,z')
    (x, y, z)
    sage: M = matrix(2,2,[x,y,z,x])
    sage: v = vector([x,y])
    sage: M * v
    (x^2 + y^2, x*y + x*z)
    sage: v*M
    (x^2 + y*z, 2*x*y)

Test if comparison bugs from :issue:`6256` are fixed::

    sage: t = exp(sqrt(x)); u = 1/t
    sage: t*u
    1
    sage: t + u
    e^(-sqrt(x)) + e^sqrt(x)
    sage: t
    e^sqrt(x)

Test if :issue:`9947` is fixed::

    sage: r = real_part(1+2*(sqrt(2)+1)*(sqrt(2)-1)); r
    2*(sqrt(2) + 1)*(sqrt(2) - 1) + 1
    sage: r.expand()
    3
    sage: a = (sqrt(4*(sqrt(3) - 5)*(sqrt(3) + 5) + 48) + 4*sqrt(3))/ (sqrt(3) + 5)
    sage: a.real_part()
    4*sqrt(3)/(sqrt(3) + 5)
    sage: a.imag_part()
    2*sqrt(10)/(sqrt(3) + 5)

Check the fix for :issue:`25251` and :issue:`25252`::

    sage: e1 = sqrt(2)*I - sqrt(2) - 2
    sage: e2 = sqrt(2)
    sage: e1 * e2
    sqrt(2)*((I - 1)*sqrt(2) - 2)
    sage: (1 + exp(I*pi/4)) * exp(I*pi/4)
    -(1/4*I + 1/4)*sqrt(2)*(-(I + 1)*sqrt(2) - 2)

Test if :issue:`24883` is fixed::

    sage: a = exp(I*pi/4) + 1
    sage: b = 1 - exp(I*pi/4)
    sage: a*b
    1/4*((I + 1)*sqrt(2) - 2)*(-(I + 1)*sqrt(2) - 2)

Test that :issue:`20784` is fixed (equations should stay unevaluated)::

    sage: limit(1/x, x=0) == unsigned_infinity
    Infinity == Infinity
    sage: SR(unsigned_infinity) == unsigned_infinity
    Infinity == Infinity

Many tests about comparison.

Use :func:`sage.symbolic.expression.mixed_order` instead of
the operators ``<=``, ``<``, etc. to compare symbolic expressions when
you do not want to get a formal inequality::

    sage: from sage.symbolic.expression import mixed_order

    sage: a = sqrt(3)
    sage: b = x^2+1
    sage: mixed_order(a, b)   # indirect doctest
    -1

    sage: x,y = var('x,y')
    sage: x < y
    x < y
    sage: mixed_order(x, y)
    1

    sage: mixed_order(SR(0.5), SR(0.7))
    -1
    sage: SR(0.5) < SR(0.7)
    0.500000000000000 < 0.700000000000000
    sage: mixed_order(SR(0.5), 0.7)
    -1

    sage: mixed_order(sin(SR(2)), sin(SR(1)))
    1
    sage: float(sin(SR(2)))
    0.9092974268256817
    sage: float(sin(SR(1)))
    0.8414709848078965

Check that :issue:`9880` is fixed::

    sage: b = [var('b_%s'%i) for i in range(4)]
    sage: precomp = (2^b_2 + 2)*(2^b_1 + 2^(-b_1) + 2^b_1*2^b_0 - \
                2^b_1*2^(-b_0) - 2^(-b_1)*2^b_0 - 2^(-b_1)*2^(-b_0) + \
                2^b_0 + 2^(-b_0) - 9) + (2^b_1 + 2^(-b_1) + \
                2^b_1*2^b_0 - 2^b_1*2^(-b_0) - 2^(-b_1)*2^b_0 - \
                 2^(-b_1)*2^(-b_0) + 2^b_0 + 2^(-b_0) - 9)/2^b_2
    sage: repl_dict = {b_0: b_0, b_3: b_1, b_2: b_3, b_1: b_2}
    sage: P = precomp.substitute(repl_dict)
    sage: P.expand()
    2^b_0*2^b_2*2^b_3 + 2*2^b_0*2^b_2 + 2^b_0*2^b_3 + 2^b_2*2^b_3 +
    2*2^b_0 + 2*2^b_2 - 9*2^b_3 + 2^b_0*2^b_2/2^b_3 -
    2^b_0*2^b_3/2^b_2 - 2^b_2*2^b_3/2^b_0 - 2*2^b_0/2^b_2 -
    2*2^b_2/2^b_0 + 2^b_0/2^b_3 + 2^b_2/2^b_3 + 2^b_3/2^b_0 +
    2^b_3/2^b_2 + 2/2^b_0 + 2/2^b_2 - 2^b_0/(2^b_2*2^b_3) -
    2^b_2/(2^b_0*2^b_3) - 9/2^b_3 - 2^b_3/(2^b_0*2^b_2) -
    2/(2^b_0*2^b_2) + 1/(2^b_0*2^b_3) + 1/(2^b_2*2^b_3) -
    1/(2^b_0*2^b_2*2^b_3) - 18

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
    64/81/((2^b_2)^2*(2^b_1)^2*(2^b_0)^2) + 35/81
    sage: f.nops()
    38

    sage: x,y,z = var('x y z')
    sage: print((-x+z)*(3*x-3*z))
    -3*(x - z)^2

    sage: t = var('t')
    sage: (x-t)^3
    -(t - x)^3
    sage: (-t+x)^3
    -(t - x)^3
    sage: (-x+t)^3
    (t - x)^3

This example is from :issue:`10833`::

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

Check if :issue:`10849` is fixed::

    sage: t = I.parent()(-1/2)
    sage: t > 0
    False
    sage: t = I*x-1/2; t
    I*x - 1/2
    sage: t.subs(x=I*x).subs(x=0).is_positive()
    False

Check if :issue:`16397` is fixed:

    sage: mixed_order(1, sqrt(2))
    -1
    sage: mixed_order(SR(1), sqrt(2))
    -1
    sage: mixed_order(log(8), 3*log(2))
    0
    sage: bool(RLF(1) < RLF(sqrt(2)))
    True
    sage: RealSet((0, pi),[pi, pi],(pi,4))
    (0, 4)
    sage: RealSet((0, pi),[0, pi],(pi,4))
    [0, 4)
    sage: RealSet((0, pi),[0, 3.5],(pi,4))
    [0, 4)

More sanity tests::

    sage: bool(pi < pi)
    False
    sage: bool(e < e)
    False
    sage: bool(sqrt(2) < sqrt(2))
    False
    sage: bool(pi < SR.zero())
    False
"""
# ****************************************************************************
#       Copyright (C) 2008-2011 William Stein <wstein@gmail.com>
#       Copyright (C) 2008-2014 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009      Alexandru Ghitza
#       Copyright (C) 2009      Bill Cauchois
#       Copyright (C) 2009-2010 Minh Van Nguyen
#       Copyright (C) 2009-2012 Robert Bradshaw
#       Copyright (C) 2009-2013 Mike Hansen
#       Copyright (C) 2009-2015 Karl-Dieter Crisman
#       Copyright (C) 2009-2019 John H. Palmieri
#       Copyright (C) 2010      Jason Grout
#       Copyright (C) 2010      Johan Sebastian Rosenkilde Nielsen
#       Copyright (C) 2010      Mitesh Patel
#       Copyright (C) 2010      Robert Marik
#       Copyright (C) 2010-2011 Dan Drake
#       Copyright (C) 2010-2019 Jeroen Demeyer
#       Copyright (C) 2011      Florent Hivert
#       Copyright (C) 2011      Jean-Pierre Flori
#       Copyright (C) 2011      Martin Albrecht
#       Copyright (C) 2011-2012 Keshav Kini
#       Copyright (C) 2011-2014 Paul Zimmermann
#       Copyright (C) 2011-2015 Volker Braun
#       Copyright (C) 2012      Andrey Novoseltsev
#       Copyright (C) 2012      Benjamin Jones
#       Copyright (C) 2012      D. S. McNeil
#       Copyright (C) 2012      David Roe
#       Copyright (C) 2012      Titus Nicolae
#       Copyright (C) 2012-2015 Karen T. Kohl
#       Copyright (C) 2012-2021 Dima Pasechnik
#       Copyright (C) 2012-2021 Michael Orlitzky
#       Copyright (C) 2013      Jean-Baptiste Priez
#       Copyright (C) 2013-2014 Eviatar Bach
#       Copyright (C) 2014      Fredrik Johansson
#       Copyright (C) 2014      Miguel Marco
#       Copyright (C) 2014      Peter Bruin
#       Copyright (C) 2014-2016 Daniel Krenn
#       Copyright (C) 2014-2018 Clemens Heuberger
#       Copyright (C) 2014-2018 Ralf Stephan
#       Copyright (C) 2014-2020 Frédéric Chapoton
#       Copyright (C) 2014-2021 Marc Mezzarobba
#       Copyright (C) 2015      Nils Bruin
#       Copyright (C) 2015-2017 Vincent Delecroix
#       Copyright (C) 2015-2019 Benjamin Hackl
#       Copyright (C) 2016      Jori Mäntysalo
#       Copyright (C) 2016      Julian Rüth
#       Copyright (C) 2017      Emmanuel Charpentier
#       Copyright (C) 2017      Marcelo Forets
#       Copyright (C) 2017      Thierry Monteil
#       Copyright (C) 2018      Eran Assaf
#       Copyright (C) 2018      Eric Gourgoulhon
#       Copyright (C) 2018      Erik M. Bray
#       Copyright (C) 2018      Martin Rubey
#       Copyright (C) 2018      Vincent Klein
#       Copyright (C) 2019      Bruno Grenet
#       Copyright (C) 2019      Markus Wageringel
#       Copyright (C) 2020-2021 Dave Witte Morris
#       Copyright (C) 2020-2021 Matthias Koeppe
#       Copyright (C) 2021      Jonathan Kliem
#       Copyright (C) 2021      Marius Gerbershagen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from sage.ext.cplusplus cimport ccrepr, ccreadstr

from copy import copy

import operator
import sage.rings.integer
import sage.rings.rational

from cpython.object cimport Py_EQ, Py_NE, Py_LE, Py_GE, Py_LT, Py_GT

from sage.cpython.string cimport str_to_bytes, char_to_str

from sage.structure.element cimport parent
from sage.structure.element cimport RingElement, Element, Matrix
from sage.structure.element cimport Expression as Expression_abc
from sage.symbolic.complexity_measures import string_length
from sage.symbolic.function cimport SymbolicFunction
from sage.rings.rational import Rational
from sage.rings.real_mpfr cimport RealNumber
from sage.misc.derivative import multi_derivative
from sage.misc.decorators import sage_wraps
from sage.misc.latex import latex_variable_name
from sage.rings.infinity import AnInfinity, infinity, minus_infinity, unsigned_infinity
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
import sage.rings.abc
from sage.structure.dynamic_class import dynamic_class
from sage.arith.numerical_approx cimport digits_to_bits

include "pynac.pxi"
include "pynac_impl.pxi"

from sage.symbolic.symbols import symbol_table, register_symbol  # used to be defined in pynac_impl


cpdef bint is_SymbolicEquation(x) noexcept:
    """
    Return ``True`` if ``x`` is a symbolic equation.

    This function is deprecated.

    EXAMPLES:

    The following two examples are symbolic equations::

        sage: from sage.symbolic.expression import is_SymbolicEquation
        sage: is_SymbolicEquation(sin(x) == x)
        doctest:warning...
        DeprecationWarning: is_SymbolicEquation is deprecated; use
        'isinstance(x, sage.structure.element.Expression) and x.is_relational()' instead
        See https://github.com/sagemath/sage/issues/35505 for details.
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
    from sage.misc.superseded import deprecation
    deprecation(35505,
                "is_SymbolicEquation is deprecated; use "
                "'isinstance(x, sage.structure.element.Expression) and x.is_relational()' instead")
    return isinstance(x, Expression) and is_a_relational((<Expression>x)._gobj)


# Defined here but exported by sage.symbolic.ring
cpdef bint _is_SymbolicVariable(x) noexcept:
    """
    Return ``True`` if ``x`` is a variable.

    EXAMPLES::

        sage: from sage.symbolic.ring import is_SymbolicVariable
        sage: is_SymbolicVariable(x)
        doctest:warning...
        DeprecationWarning: is_SymbolicVariable is deprecated; use
        'isinstance(x, sage.structure.element.Expression) and x.is_symbol()' instead
        See https://github.com/sagemath/sage/issues/35505 for details.
        True
        sage: is_SymbolicVariable(x + 2)
        False

    TESTS::

        sage: ZZ['x']
        Univariate Polynomial Ring in x over Integer Ring
    """
    from sage.misc.superseded import deprecation
    deprecation(35505,
                "is_SymbolicVariable is deprecated; use "
                "'isinstance(x, sage.structure.element.Expression) and x.is_symbol()' instead")
    return isinstance(x, Expression) and is_a_symbol((<Expression>x)._gobj)


def _dict_update_check_duplicate(dict d1, dict d2):
    r"""
    Merge the dictionary ``d2`` into ``d1`` and check for duplicates.

    The two dictionaries must be of the form ``{expr: replacement}``. This
    function throws a :exc:`ValueError` if any expressions are substituted
    for twice.

    EXAMPLES:

    A normal merge with no conflicts::

        sage: from sage.symbolic.expression import _dict_update_check_duplicate
        sage: d1 = {'a': 1}
        sage: d2 = {'b': 2}
        sage: _dict_update_check_duplicate(d1, d2)
        sage: d1 == {'a': 1, 'b': 2}
        True

    In this case, the variable ``a`` is substituted twice resulting in
    an error::

        sage: from sage.symbolic.expression import _dict_update_check_duplicate
        sage: d1 = {'a': 1}
        sage: d2 = {'a': 2}
        sage: _dict_update_check_duplicate(d1, d2)
        Traceback (most recent call last):
        ...
        ValueError: duplicate substitution for a, got values 1 and 2

    We report only the first conflict (according to the Python sort
    order)::

        sage: from sage.symbolic.expression import _dict_update_check_duplicate
        sage: d1 = {'b': 1, 'a': 1}
        sage: d2 = {'b': 2, 'a': 2}
        sage: _dict_update_check_duplicate(d1, d2)
        Traceback (most recent call last):
        ...
        ValueError: duplicate substitution for a, got values 1 and 2
    """
    # We need to check for duplicates in a predictable order so that
    # errors are reported reliably. We only need to sort one of the
    # dictionaries to achieve that, and we suspect that d2 will
    # generally be smaller than d1, so we sort d2. This gives us a
    # list of d2's keys.
    #
    # When sorting d2, we compare the string representations of its
    # keys rather than the keys themselves. This is because comparison
    # of symbolic expressions doesn't do what the sorted() function
    # needs: `x <= y` is a symbolic inequality, and we need a
    # True/False answer. The expression 'x' <= 'y' on the other hand
    # is unambiguous.
    #
    for k in sorted(d2, key=str):
        if k in d1:
            msg = "duplicate substitution for {}, got values {} and {}"
            raise ValueError(msg.format(k, d1[k], d2[k]))

    d1.update(d2)


def _subs_make_dict(s):
    r"""
    There are a few ways we can represent a substitution. The first is
    a symbolic equation. The second is a dictionary. The third would
    be a list/tuple whose entries are expressions, dictionaries, or
    lists/tuples themselves. This function converts all such
    representations to dictionaries.

    INPUT:

    - ``s`` -- a representation of a substitution

    OUTPUT: a dictionary of substitutions

    EXAMPLES:

    An expression::

        sage: from sage.symbolic.expression import _subs_make_dict
        sage: _subs_make_dict(x == 1)
        {x: 1}

    And a dictionary (we just return it as-is)::

        sage: _subs_make_dict({x: 1})
        {x: 1}

    And finally, a tuple or a list containing one of everything::

        sage: w, x, y, z = SR.var('w, x, y, z')
        sage: actual = _subs_make_dict([w == 1, {x: 1}, [y == 1], (z == 1,)])
        sage: expected = {w: 1, y: 1, x: 1, z: 1}
        sage: actual == expected
        True

    Note that it recursively calls itself so that the following does work::

        sage: x, y, z = SR.var('x, y, z')
        sage: actual = _subs_make_dict([[x == 1], [[y == 2], [z == 3]]])
        sage: expected = {z: 3, y: 2, x: 1}
        sage: actual == expected
        True

    Check that a :exc:`TypeError` is raised if the input is not valid::

        sage: _subs_make_dict(1)
        Traceback (most recent call last):
        ...
        TypeError: not able to determine a substitution from 1
        sage: _subs_make_dict(x)
        Traceback (most recent call last):
        ...
        TypeError: not able to determine a substitution from x
        sage: _subs_make_dict(x <= 1)
        Traceback (most recent call last):
        ...
        TypeError: can only substitute equality, not inequalities; got x <= 1
    """
    if isinstance(s, dict):
        return s
    elif isinstance(s, Expression) and s.is_relational():
        if s.operator() is not operator.eq:
            msg = "can only substitute equality, not inequalities; got {}"
            raise TypeError(msg.format(s))
        return {s.lhs(): s.rhs()}
    elif isinstance(s, (tuple, list)):
        result = {}
        for d in s:
            _dict_update_check_duplicate(result, _subs_make_dict(d))
        return result
    else:
        msg = "not able to determine a substitution from {}"
        raise TypeError(msg.format(s))


def _subs_fun_make_dict(s):
    r"""
    There are a few ways we can represent a substitution. The first is
    a symbolic equation. The second is a dictionary. The third would
    be a list/tuple whose entries are expressions, dictionaries, or
    lists/tuples themselves. This function converts all such
    representations to dictionaries.

    INPUT:

    - ``s`` -- a representation of a substitution

    OUTPUT: a dictionary of substitutions

    EXAMPLES:

    An expression::

        sage: from sage.symbolic.expression import _subs_fun_make_dict
        sage: f, g = function('f g')
        sage: _subs_fun_make_dict(f(x) == g(x))
        {f: x |--> g(x)}

    And a dictionary::

        sage: _subs_fun_make_dict({f: g})
        {f: g}
        sage: _subs_fun_make_dict({f(x): g(x+1)})
        {f: x |--> g(x + 1)}

    And finally, a tuple or a list containing one of everything::

        sage: f1, g1, f2, g2, f3, g3, f4, g4 = function('f1 g1 f2 g2 f3 g3 f4 g4')
        sage: actual = _subs_fun_make_dict([f1(x) == g1(x), {f2: g2}, [f3(x) == g3(x+1)], (f4(x) == g4(x),)])
        sage: expected = {f1: (g1(x)).function(x), f2: g2, f3: (g3(x + 1)).function(x), f4: (g4(x)).function(x)}
        sage: actual == expected
        True

    Note that it recursively calls itself so that the following does work::

        sage: actual = _subs_fun_make_dict([[f1(x) == g1(x)], [[f2(x) == g2(x)], [f3(x) == g3(x)]]])
        sage: expected = {f1: (g1(x)).function(x), f2: (g2(x)).function(x), f3: (g3(x)).function(x)}
        sage: actual == expected
        True

    Check that a :exc:`TypeError` is raised if the input is not valid::

        sage: _subs_fun_make_dict(1)
        Traceback (most recent call last):
        ...
        TypeError: not able to determine a substitution from 1
        sage: _subs_fun_make_dict(f)
        Traceback (most recent call last):
        ...
        TypeError: not able to determine a substitution from f
        sage: _subs_fun_make_dict(f(x) <= g(x))
        Traceback (most recent call last):
        ...
        TypeError: can only substitute equality, not inequalities; got f(x) <= g(x)
    """
    if isinstance(s, dict):
        return dict((k, v) if not isinstance(k, Expression) else (k.operator(), v.function(*k.operands())) for k, v in s.items())
    elif isinstance(s, Expression) and s.is_relational():
        if s.operator() is not operator.eq:
            msg = "can only substitute equality, not inequalities; got {}"
            raise TypeError(msg.format(s))
        temp_fun = s.rhs().function(*s.lhs().operands())
        return {s.lhs().operator(): temp_fun}
    elif isinstance(s, (tuple, list)):
        result = {}
        for d in s:
            _dict_update_check_duplicate(result, _subs_fun_make_dict(d))
        return result
    else:
        msg = "not able to determine a substitution from {}"
        raise TypeError(msg.format(s))


cdef class Expression(Expression_abc):

    cdef GEx _gobj

    cpdef object pyobject(self):
        """
        Get the underlying Python object.

        OUTPUT:

        The Python object corresponding to this expression, assuming
        this expression is a single numerical value or an infinity
        representable in Python. Otherwise, a :exc:`TypeError` is raised.

        EXAMPLES::

            sage: var('x')
            x
            sage: b = -17.3
            sage: a = SR(b)
            sage: a.pyobject()
            -17.3000000000000
            sage: a.pyobject() is b
            True

        Integers and Rationals are converted internally though, so you
        won't get back the same object::

            sage: b = -17/3
            sage: a = SR(b)
            sage: a.pyobject()
            -17/3
            sage: a.pyobject() is b
            False

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
        if is_a_constant(self._gobj):
            from sage.symbolic.constants import constants_name_table
            return constants_name_table[ccrepr(self._gobj)]

        if is_a_infinity(self._gobj):
            if (ex_to_infinity(self._gobj).is_unsigned_infinity()):
                return unsigned_infinity
            if (ex_to_infinity(self._gobj).is_plus_infinity()):
                return infinity
            if (ex_to_infinity(self._gobj).is_minus_infinity()):
                return minus_infinity
            raise TypeError('Python infinity cannot have complex phase.')

        if not is_a_numeric(self._gobj):
            raise TypeError("self must be a numeric expression")
        return py_object_from_numeric(self._gobj)

    def __init__(self, SR, x=0):
        """
        Nearly all expressions are created by calling new_Expression_from_*,
        but we need to make sure this at least does not leave self._gobj
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
        self._gobj = GEx(exp._gobj)

    def __getstate__(self):
        """
        Return a tuple describing the state of this expression for pickling.

        This should return all information that will be required to unpickle
        the object. The functionality for unpickling is implemented in
        __setstate__().

        In order to pickle Expression objects, we return a tuple containing

         * 0  - as pickle version number
                in case we decide to change the pickle format in the feature
         * names of symbols of this expression
         * a string representation of ``self`` stored in a Pynac archive.

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
        return (0, [repr(x) for x in self.variables()], ccrepr(ar))

    def _dbgprint(self):
        r"""
        Print pynac debug output to ``stderr``.

        EXAMPLES::

            sage: (1+x)._dbgprint()
            x + 1
        """
        sig_on()
        try:
            self._gobj.dbgprint()
        finally:
            sig_off()

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

        Check that user-defined functions get the same treatment (:issue:`19194`)::

            sage: f = function('f')(x)
            sage: f._dbgprinttree()
            function f ...
                x (symbol) ...
                =====
        """
        sig_on()
        try:
            self._gobj.dbgprinttree()
        finally:
            sig_off()

    def __setstate__(self, state):
        """
        Initialize the state of the object from data saved in a pickle.

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
            raise ValueError("unknown state information")
        # set parent
        from sage.symbolic.ring import SR
        self._parent = SR
        # get variables
        cdef GExList sym_lst
        for name in state[1]:
            sym_lst.append_sym(
                ex_to_symbol((<Expression>SR.symbol(name))._gobj))

        # initialize archive
        cdef GArchive ar
        ccreadstr(ar, state[2])

        # extract the expression from the archive
        self._gobj = GEx(ar.unarchive_ex(sym_lst, <unsigned>0))

    def __copy__(self):
        """
        TESTS::

            sage: copy(x)
            x
        """
        return new_Expression_from_GEx(self._parent, self._gobj)

    def __enter__(self):
        """
        Method used by temporary variables with Python `with` to
        automatically clean up after themselves.
        """
        return self

    def __exit__(self, *args):
        """
        Method used by temporary variables with Python `with` to
        automatically clean up after themselves.

        TESTS::

            sage: symbols_copy = SR.symbols.copy()
            sage: with SR.temp_var() as t: pass
            sage: symbols_copy == SR.symbols
            True
        """
        from sage.symbolic.ring import SR
        SR.cleanup_var(self)
        return False

    def _repr_(self):
        r"""
        Return string representation of this symbolic expression.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: repr(x+y)
            'x + y'

        TESTS:

        Rational functions::

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
            1/(x^y)
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

        Check if :issue:`7876` is fixed::

            sage: (1/2-1/2*I)*sqrt(2)
            -(1/2*I - 1/2)*sqrt(2)
            sage: latex((1/2-1/2*I)*sqrt(2))
            -\left(\frac{1}{2} i - \frac{1}{2}\right) \, \sqrt{2}

        Check if :issue:`9632` is fixed::

            sage: zeta(x) + cos(x)
            cos(x) + zeta(x)
            sage: psi(1,1/3)*log(3)
            log(3)*psi(1, 1/3)
        """
        return self._parent._repr_element_(self)

    def _sympy_character_art(self, use_unicode):
        r"""
        Create character art using Sympy.

        INPUT:

        - ``use_unicode`` -- boolean; whether to allow unicode instead
          of 7-bit clean output

        OUTPUT: string

        EXAMPLES::

            sage: i = var('i')
            sage: f = integral(exp(x + x^2)/(x+1), x)
            ...
            sage: f._sympy_character_art(False)                                         # needs sympy
            '  /          \n |           \n |   2       \n |  x  + x   \n | e...'
        """
        from sympy import pretty, sympify
        # FIXME:: when *sage* will use at least sympy >= 0.7.2
        # we could use a nice splitting with respect of the AsciiArt module.
        # from sage.typeset.ascii_art import AsciiArt, MAX_LENGTH ## for import
        #            num_columns = MAX_LENGTH  ## option of pretty
        try:
            return pretty(sympify(self, evaluate=False), use_unicode=use_unicode)
        except Exception:
            return str(self)

    def _ascii_art_(self):
        """
        Ascii art magic method.

        See :mod:`sage.typeset.ascii_art` for details.

        EXAMPLES::

            sage: i = var('i')
            sage: ascii_art(sum(i^2/pi*x^i, i, 0, oo))                                  # needs sympy
                          2
                         x  + x
            -------------------------------
                  3         2
            - pi*x  + 3*pi*x  - 3*pi*x + pi
            sage: ascii_art(integral(exp(x + x^2)/(x+1), x))                            # needs sympy
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
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._sympy_character_art(False).splitlines())

    def _unicode_art_(self):
        """
        Unicode art magic method.

        See :mod:`sage.typeset.unicode_art` for details.

        EXAMPLES::

            sage: i = var('i')
            sage: unicode_art(sum(i^2/pi*x^i, i, 0, oo))                                # needs sympy
                        2
                       x  + x
            ───────────────────────────
                 3        2
            - π⋅x  + 3⋅π⋅x  - 3⋅π⋅x + π
            sage: unicode_art(integral(exp(x + x^2)/(x+1), x))                          # needs sympy
            ⌠
            ⎮   2
            ⎮  x  + x
            ⎮ ℯ
            ⎮ ─────── dx
            ⎮  x + 1
            ⌡

        TESTS:

        Check that :issue:`28891` is fixed::

            sage: unicode_art(exp(x).series(x, 4))                                      # needs sympy
                     2    3
                    x    x     ⎛ 4⎞
            1 + x + ── + ── + O⎝x ⎠
                    2    6
            sage: unicode_art(exp(x).series(x==1, 3))                                   # needs sympy
                                     2
                            ℯ⋅(x - 1)     ⎛       3       ⎞
            ℯ + ℯ⋅(x - 1) + ────────── + O⎝(x - 1) ; x → 1⎠
                                2

        Check that complex numbers are handled correctly (:issue:`28903`)::

            sage: # needs sympy
            sage: unicode_art(SR(I))
            ⅈ
            sage: unicode_art(SR(13 - I))
            13 - ⅈ
            sage: unicode_art(SR(1.3 - I))
            1.3 - ...ⅈ
            sage: unicode_art(cos(I))
            cosh(1)

            sage: unicode_art(SR(CC(1/3, 1)))                                           # needs sympy
            0.333333333333333 + 1.0⋅ⅈ
            sage: unicode_art(SR(CDF(1/3, 1)))                                          # needs sympy
            0.333333333333333 + 1.0⋅ⅈ
            sage: unicode_art(SR(RealField(100)(1/7)))                                  # needs sympy
            0.14285714285714285714285714286

            sage: K.<a> = QuadraticField(-1)                                            # needs sympy
            sage: unicode_art(SR(2 + a))                                                # needs sympy
            2 + ⅈ
            sage: unicode_art(SR(1/3 + a/2))                                            # needs sympy
            1   ⅈ
            ─ + ─
            3   2
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self._sympy_character_art(True).splitlines())

    def _interface_(self, I):
        """
        EXAMPLES::

            sage: f = sin(e + 2)
            sage: f._interface_(sage.calculus.calculus.maxima)
            sin(%e+2)
        """
        if is_a_constant(self._gobj):
            return self.pyobject()._interface_(I)
        return super()._interface_(I)

    def _maxima_(self, session=None):
        """
        EXAMPLES::

            sage: f = sin(e + 2)
            sage: f._maxima_()
            sin(%e+2)
            sage: _.parent() is sage.calculus.calculus.maxima
            True

        TESTS:

        We test that unicode characters are handled correctly
        :issue:`30122`::

            sage: var('ξ')._maxima_()
            _SAGE_VAR_ξ
        """
        if session is None:
            # This chooses the Maxima interface used by calculus
            # Maybe not such a great idea because the "default" interface is another one
            from sage.calculus.calculus import maxima
            return super()._interface_(maxima)
        else:
            return super()._interface_(session)

    def _interface_init_(self, I):
        """
        EXAMPLES::

            sage: a = (pi + 2).sin()
            sage: a._maxima_init_()
            'sin((%pi)+(2))'

            sage: a = (pi + 2).sin()
            sage: a._maple_init_()
            'sin((Pi)+(2))'

            sage: a = (pi + 2).sin()
            sage: a._mathematica_init_()
            'Sin[(Pi)+(2)]'

            sage: f = pi + I*e
            sage: f._pari_init_()
            '(Pi)+((exp(1))*(I))'

        TESTS:

        Check if complex numbers are converted to Maxima correctly
        :issue:`7557`::

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
        Convert symbolic object to GAP string.

        EXAMPLES::

            sage: gap(e + pi^2 + x^3)                                                   # needs sage.libs.gap
            x^3 + pi^2 + e
        """
        return '"%s"' % repr(self)

    def _singular_init_(self):
        """
        Conversion of a symbolic object to Singular string.

        EXAMPLES::

            sage: singular(e + pi^2 + x^3)                                              # needs sage.libs.singular
            x^3 + pi^2 + e
        """
        return '"%s"' % repr(self)

    def _magma_init_(self, magma):
        """
        Return string representation in Magma of this symbolic expression.

        Since Magma has no notation of symbolic calculus, this simply
        returns something that evaluates in Magma to a Magma string.

        EXAMPLES::

            sage: x = var('x')
            sage: f = sin(cos(x^2) + log(x))
            sage: f._magma_init_(magma)
            '"sin(cos(x^2) + log(x))"'
            sage: magma(f)                         # optional - magma
            sin(cos(x^2) + log(x))
            sage: magma(f).Type()                  # optional - magma
            MonStgElt
        """
        return '"%s"' % repr(self)

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

        Check spacing of coefficients of mul expressions (:issue:`3202` and
        :issue:`13356`)::

            sage: latex(2*3^x)
            2 \cdot 3^{x}
            sage: latex(1/2/3^x)
            \frac{1}{2 \cdot 3^{x}}
            sage: latex(1/2*3^x)
            \frac{1}{2} \cdot 3^{x}

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

        More powers (:issue:`7406`)::

            sage: latex((x^pi)^e)
            {\left(x^{\pi}\right)}^{e}
            sage: latex((x^(pi+1))^e)
            {\left(x^{\pi + 1}\right)}^{e}
            sage: a,b,c = var('a b c')
            sage: latex(a^(b^c))
            a^{\left(b^{c}\right)}
            sage: latex((a^b)^c)
            {\left(a^{b}\right)}^{c}

        Separate coefficients to numerator and denominator (:issue:`7363`)::

            sage: latex(2/(x+1))
            \frac{2}{x + 1}
            sage: latex(1/2/(x+1))
            \frac{1}{2 \, {\left(x + 1\right)}}

        Check if rational function coefficients without a ``numerator()`` method
        are printed correctly. :issue:`8491`::

            sage: latex(6.5/x)
            \frac{6.50000000000000}{x}

        Check if we avoid extra parenthesis in rational functions (:issue:`8688`)::

            sage: latex((x+2)/(x^3+1))
            \frac{x + 2}{x^{3} + 1}
            sage: latex((x+2)*(x+1)/(x^3+1))
            \frac{{\left(x + 2\right)} {\left(x + 1\right)}}{x^{3} + 1}
            sage: latex((x+2)/(x^3+1)/(x+1))
            \frac{x + 2}{{\left(x^{3} + 1\right)} {\left(x + 1\right)}}

        Check that the sign is correct (:issue:`9086`)::

            sage: latex(-1/x)
            -\frac{1}{x}
            sage: latex(1/-x)
            -\frac{1}{x}

        More tests for the sign (:issue:`9314`)::

            sage: latex(-2/x)
            -\frac{2}{x}
            sage: latex(-x/y)
            -\frac{x}{y}
            sage: latex(-x*z/y)
            -\frac{x z}{y}
            sage: latex(-x/z/y)
            -\frac{x}{y z}

        Check if :issue:`9394` is fixed::

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
        correctly (:issue:`10769`)::

            sage: var('a x')
            (a, x)
            sage: latex(1/2*I/x)
            \frac{i}{2 \, x}
            sage: ratio = i/2* x^2/a
            sage: latex(ratio)
            \frac{i \, x^{2}}{2 \, a}

        Parenthesis in powers (:issue:`13262`)::

            sage: latex(1+x^(2/3)+x^(-2/3))
            x^{\frac{2}{3}} + \frac{1}{x^{\frac{2}{3}}} + 1

        Check that pynac understands rational powers (:issue:`30446`)::

            sage: QQ((24*sqrt(3))^(100/50))==1728
            True
            sage: float((24*sqrt(3))^(100/51))
            1493.0092154...
        """
        return self._parent._latex_element_(self)

    def _mathml_(self):
        """
        Return a MathML representation of this object.

        EXAMPLES::

            sage: mathml(pi)
            <mi>&pi;</mi>
            sage: mathml(pi+2)
            MATHML version of the string pi + 2
        """
        from sage.misc.mathml import mathml
        try:
            obj = self.pyobject()
        except TypeError:
            return mathml(repr(self))
        return mathml(obj)

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: f = x^3 + 17*x -3
            sage: ZZ(f.coefficient(x^3))
            1
            sage: ZZ(f.coefficient(x))
            17
            sage: ZZ(f.coefficient(x,0))
            -3
            sage: type(ZZ(f.coefficient(x,0)))
            <class 'sage.rings.integer.Integer'>

        Coercion is done if necessary::

            sage: f = x^3 + 17/1*x
            sage: ZZ(f.coefficient(x))
            17
            sage: type(ZZ(f.coefficient(x)))
            <class 'sage.rings.integer.Integer'>

        If the symbolic expression is just a wrapper around an integer,
        that very same integer is not preserved, but a new one returned::

            sage: n = 17; SR(n)._integer_() is n
            False
        """
        try:
            n = self.pyobject()
        except TypeError:
            raise TypeError("unable to convert %r to an integer" % self)
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
        from sage.rings.real_mpfi import RIF
        try:
            rif_self = RIF(self)
        except TypeError:
            raise ValueError("cannot convert %s to int" % self)
        if rif_self > 0 or (rif_self.contains_zero() and self > 0):
            result = floor(self)
        else:
            result = ceil(self)
        if not isinstance(result, sage.rings.integer.Integer):
            raise ValueError("cannot convert %s to int" % self)
        else:
            return int(result)

    def _rational_(self):
        """
        EXAMPLES::

            sage: f = x^3 + 17/1*x - 3/8
            sage: QQ(f.coefficient(x^2))
            0
            sage: QQ(f.coefficient(x^3))
            1
            sage: a = QQ(f.coefficient(x)); a
            17
            sage: type(a)
            <class 'sage.rings.rational.Rational'>
            sage: QQ(f.coefficient(x,0))
            -3/8

        If the symbolic expression is just a wrapper around a rational,
        that very same rational is not preserved, but a new one returned::

            sage: n = 17/1; SR(n)._rational_() is n
            False
        """
        try:
            n = self.pyobject()
        except TypeError:
            raise TypeError("unable to convert %s to a rational" % self)
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
            TypeError: cannot evaluate symbolic expression to a numeric value

        Check if we can compute a real evaluation even if the expression
        contains complex coefficients::

            sage: RR((I - sqrt(2))*(I+sqrt(2)))
            -3.00000000000000
            sage: cos(I)._eval_self(RR)
            1.54308063481524
            sage: float(cos(I))  # abs tol 1e-15
            1.5430806348152437

        TESTS::

            sage: e = sqrt(2)/sqrt(abs(-(I - 1)*sqrt(2) - I - 1))
            sage: e._eval_self(float)
            0.9036020036...
        """
        try:
            res = self._convert({'parent':R})
        except TypeError as err:
            # try the evaluation again with the complex field
            # corresponding to the parent R
            if R in (float, complex):
                R_complex = complex
            else:
                try:
                    R_complex = R.complex_field()
                except (TypeError, AttributeError):
                    raise err
            res = self._convert({'parent':R_complex})

        if res.is_numeric():
            ans = res.pyobject()
            # Convert ans to R.
            if R is float and isinstance(ans, complex) and not ans.imag:
                # Python does not automatically convert "real" complex
                # numbers to floats, so we do this manually:
                ans = ans.real
            return R(ans)
        else:
            raise TypeError("cannot evaluate symbolic expression to a numeric value")

    cpdef _convert(self, kwds):
        """
        Convert all the numeric coefficients and constants in this expression
        to the given ring ``R``. This results in an expression which contains
        only variables, and functions whose arguments contain a variable.

        EXAMPLES::

            sage: f = sqrt(2) * cos(3); f
            sqrt(2)*cos(3)
            sage: f._convert({'parent':RDF})
            -1.40006081533995
            sage: f._convert({'parent':float})
            -1.4000608153399503

        There is nothing to convert for variables::

            sage: x._convert({'parent':CC})
            x

        Note that the output is not meant to be in the in the given ring ``R``.
        Since the results of some functions will still be  floating point
        approximations::

            sage: t = log(10); t
            log(10)
            sage: t._convert({'parent':ZZ})
            log(10)

        ::

            sage: (0.25 / (log(5.74 /x^0.9, 10))^2 / 4)._convert({'parent':QQ})
            1/16*log(10)^2/log(287/50/x^0.900000000000000)^2
            sage: (0.25 / (log(5.74 /x^0.9, 10))^2 / 4)._convert({'parent':CC})
            0.331368631904900/log(5.74000000000000/x^0.900000000000000)^2

        When converting to an exact domain, powers remain unevaluated::

            sage: f = sqrt(2) * cos(3); f
            sqrt(2)*cos(3)
            sage: (sqrt(2))._convert({'parent':int})
            sqrt(2)
            sage: f._convert({'parent':int})
            0

        If ``R`` has an associated complex field it is used with complex
        input::

            sage: SR(CBF(1+I))._convert({'parent':RDF})
            1.0 + 1.0*I
            sage: type(_.pyobject())
            <class 'sage.rings.complex_double.ComplexDoubleElement'>
            sage: SR(CBF(1+I))._convert({'parent':CDF})
            1.0 + 1.0*I
            sage: SR(RBF(1))._convert({'parent':RDF})
            1.0
            sage: SR(CBF(1))._convert({'parent':RDF})
            1.0
            sage: type(_.pyobject())
            <class 'sage.rings.real_double...RealDoubleElement...'>
        """
        cdef GEx res = self._gobj.evalf(0, kwds)
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
        Return this expression as a real interval.

        EXAMPLES::

            sage: RIF(sqrt(2))
            1.414213562373095?
        """
        try:
            return self._eval_self(R)
        except TypeError:
            raise TypeError("unable to simplify to a real interval approximation")

    def _complex_mpfi_(self, R):
        """
        Return this expression as a complex interval.

        EXAMPLES::

            sage: CIF(pi)
            3.141592653589794?
        """
        try:
            return self._eval_self(R)
        except TypeError:
            raise TypeError("unable to simplify to a complex interval approximation")

    def _arb_(self, R):
        r"""
        Convert this expression to a real or complex ball.

        (In spite of its name, this method also works in the complex case.)

        .. WARNING::

            The generic conversion mechanism is fragile. When rigorous results
            are essential, it is recommended to call suitable methods of real
            or complex balls instead. For example, `RBF(real(i + 1))` is better
            expressed as `(CBF(i) + 1).real()`.

        EXAMPLES::

            sage: RBF(pi, 1/1000)
            [3.14 +/- 2.60e-3]
            sage: RBF(pi/2 + 2*arctan(1))
            [3.14159265358979...]
            sage: (pi + I)._arb_(CBF)
            [3.14159265358979...] + 1.000000000000000*I
            sage: RBF(x)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a RealBall

        TESTS::

            sage: CBF(gamma(15/2, 1)).identical(CBF(15/2).gamma(1))
            True
            sage: a = RBF(abs(e+i)); (a, a.parent())
            ([2.89638673159001 +/- 3.07e-15], Real ball field with 53 bits of precision)
            sage: a = CBF(abs(e+i)); (a, a.parent())
            ([2.89638673159001 +/- 3.07e-15], Complex ball field with 53 bits of precision)
            sage: RBF(sin(7/12)^2 + real(exp(i*7/12))^2)
            [1.0000000000000 +/- 1.12e-15]
            sage: RBF(arg(sin(i+1)))
            [0.454820233309950 +/- 7.08e-16]
            sage: RBF(abs(i) + i)
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part

            sage: from sage.symbolic.function import BuiltinFunction
            sage: f = BuiltinFunction("bogus_builtin_function")
            sage: RBF(f(1))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert bogus_builtin_function(1) to a RealBall
            sage: CBF(f(1))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert bogus_builtin_function(1) to a ComplexBall

            sage: CBF(acos(float(1/4) * sqrt(int(5))))
            [0.97759655064526...]
        """
        cdef bint progress = False
        cdef int i
        # Note that we deliberately don't use _eval_self and don't try going
        # through RIF/CIF in order to avoid unsafe conversions.
        operator = self.operator()
        # Constants
        if operator is None:
            try:
                return R(self.pyobject())
            except (TypeError, ValueError):
                pass
        else:
            C = R.complex_field()
            # Intended for BuiltinFunctions with a well-defined main argument
            args = [a.pyobject() if a.is_numeric() else a
                    for a in self.operands()]
            try:
                args = operator._method_arguments(*args)
            except AttributeError:
                pass
            else:
                # If R is a real field, prefer the method of real balls if it
                # exists (sometimes leads to tighter bounds, and avoids
                # confusing inconsistencies).
                for T in ([R] if C is R else [R, C]):
                    try:
                        method = getattr(T(args[0]), operator.name())
                    except (AttributeError, TypeError, ValueError):
                        pass
                    else:
                        if callable(method):
                            return method(*args[1:])
            # Generic case: walk through the expression. In this case, we do
            # not bother trying to stay in the real field.
            try:
                for i in range(len(args)):
                    if parent(args[i]) is not C:
                        progress = True
                        args[i] = C(args[i])
            except (TypeError, ValueError):
                pass
            if progress:
                try:
                    res = self.operator()(*[C(a) for a in args])
                    return R(res)
                except (TypeError, ValueError):
                    pass
        # Typically more informative and consistent than the exceptions that
        # would propagate
        raise TypeError("unable to convert {!r} to a {!s}".format(
                self, R.element_class.__name__))

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(sin(3))
            0.1411200080598672
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
            0.09090909090909091
            sage: zeta(x).subs(x=I)._complex_double_(CDF)  # rel tol 1e-16
            0.003300223685324103 - 0.4181554491413217*I
            sage: gamma(x).subs(x=I)._complex_double_(CDF)
            -0.15494982830181067 - 0.49801566811835607*I
            sage: log(x).subs(x=I)._complex_double_(CDF)
            1.5707963267948966*I
            sage: CDF((-1)^(1/3))
            0.5000000000000001 + 0.8660254037844386*I
        """
        return self._eval_self(R)

    def __float__(self):
        """
        Return float conversion of ``self``, assuming ``self`` is constant.

        Otherwise, raise a :exc:`TypeError`.

        OUTPUT: a ``float``. Double precision evaluation of ``self``

        EXAMPLES::

            sage: float(SR(12))
            12.0
            sage: float(SR(2/3))
            0.6666666666666666
            sage: float(sqrt(SR(2)))
            1.4142135623730951
            sage: float(SR(RIF(2)))
            2.0
            sage: float(x^2 + 1)
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to float approximation

        TESTS::

            sage: float(sqrt(2)/sqrt(abs(-(I - 1)*sqrt(2) - I - 1)))
            0.9036020036...
        """
        from sage.functions.other import real, imag
        try:
            ret = float(self._eval_self(float))
        except TypeError:
            try:
                c = (self._eval_self(complex))
                if imag(c) == 0:
                    ret = real(c)
                else:
                    raise
            except TypeError:
                raise TypeError("unable to simplify to float approximation")
        return ret

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
            raise TypeError("unable to simplify to complex approximation")

    def _sympy_(self):
        """
        Return a Sympy version of this object.

        EXAMPLES::

            sage: pi._sympy_()
            pi
            sage: type(_)
            <class 'sympy.core.numbers.Pi'>
        """
        from sage.symbolic.expression_conversions import sympy_converter
        return sympy_converter(self)

    def _fricas_init_(self):
        """
        Return a FriCAS version of this object.

        EXAMPLES::

            sage: pi._fricas_()                                                 # optional - fricas
            %pi
        """
        from sage.symbolic.expression_conversions import fricas_converter
        return fricas_converter(self)

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
            sage: QQbar(SR(2*I)^(1/2))
            1 + 1*I
            sage: QQbar(e^(pi*I/3))
            0.50000000000000000? + 0.866025403784439?*I

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
        r"""
        Return hash of this expression.

        EXAMPLES:

        The hash of an object in Python or its coerced version into
        the symbolic ring is usually the same::

            sage: hash(SR(3.1)) == hash(3.1)
            True
            sage: hash(SR(19.23)) == hash(19.23)
            True
            sage: hash(SR(3/1))
            3
            sage: hash(SR(19/23)) == hash(19/23)
            True
            sage: hash(SR(2^32)) == hash(2^32)
            True
            sage: hash(SR(2^64-1)) == hash(2^64-1)
            True
            sage: hash(SR(1e100)) == hash(1e100)
            True

        The hash for symbolic expressions are unfortunately random. Here we
        only test that the hash() function returns without error, and that
        the return type is correct::

            sage: x, y = var("x y")
            sage: t = hash(x); type(t)
            <... 'int'>
            sage: t = hash(x^y); type(t)
            <... 'int'>
            sage: type(hash(x+y))
            <... 'int'>
            sage: d = {x+y: 5}
            sage: d
            {x + y: 5}

        In this example hashing is important otherwise the answer is
        wrong::

            sage: set([x-x, -x+x])
            {0}

        Test if exceptions during hashing are handled properly::

            sage: t = SR(matrix(2,2,range(4)))
            sage: hash(t)
            Traceback (most recent call last):
            ...
            RuntimeError: Python object not hashable

        TESTS:

        Test if hashes for fderivatives with different parameters collide.
        :issue:`6243`::

            sage: f = function('f'); t = f(x,y)
            sage: u = t.derivative(x); v = t.derivative(y)
            sage: hash(u) == hash(v)
            False
            sage: d = {u: 3, v: 5}; sorted(d.values())
            [3, 5]

        More checks for fderivative hashes :issue:`6851` ::

            sage: hash(f(x).derivative(x)) == hash(f(x).derivative(x,2))
            False
            sage: d = dict( (f(x).derivative(x, i), i) for i in range(1,6) )
            sage: len(d.keys())
            5

        We create a function with 10 arguments and test if there are
        hash collisions between any of its derivatives of order at
        most 7. :issue:`7508` ::

            sage: num_vars = 10; max_order=7
            sage: X = var(' '.join('x' + str(i) for i in range(num_vars)))
            sage: f = function('f')(*X)
            sage: hashes = set()
            sage: for length in range(1,max_order+1):  # long time (4s on sage.math, 2012)
            ....:     for s in UnorderedTuples(X, length):
            ....:         deriv = f.diff(*s)
            ....:         h = hash(deriv)
            ....:         if h in hashes:
            ....:             print("deriv: %s, hash:%s" % (deriv, h))
            ....:         else:
            ....:             hashes.add(n)

        Check whether `oo` keeps its hash in `SR` (:issue:`19928`)::

            sage: hash(oo) == hash(SR(oo))
            True
            sage: hash(oo) == hash((-x).subs(x=-oo))
            True
            sage: hash(-oo) == hash(SR(-oo))
            True
            sage: hash(-oo) == hash((-x).subs(x=oo))
            True
            sage: hash(unsigned_infinity) == hash(SR(unsigned_infinity))
            True

        Check a corner case for rational numbers (:issue:`28219`)::

            sage: hash(-1/3) == hash(SR(-1/3))
            True
        """
        sig_on()
        try:
            return self._gobj.gethash()
        finally:
            sig_off()

    cdef bint _rel_equal1(Expression self, Expression other) except -1:
        """
        Internal helper function.
        """
        sig_on()
        try:
            return (self._gobj.lhs().is_equal(other._gobj.lhs())
                    and self._gobj.rhs().is_equal(other._gobj.rhs()))
        finally:
            sig_off()

    cdef bint _rel_equal2(Expression self, Expression other) except -1:
        """
        Internal helper function.
        """
        sig_on()
        try:
            return (self._gobj.lhs().is_equal(other._gobj.rhs())
                    and self._gobj.rhs().is_equal(other._gobj.lhs()))
        finally:
            sig_off()

    cpdef _richcmp_(left, right, int op):
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

        Testing :issue:`11309` which changes the behavior of comparison of
        comparisons::

            sage: (-x + y < 0) in [x - y < 0]
            False
            sage: (x - 1 < 0) in [x - 2 < 0]
            False
            sage: len(Set([-x + y < 0, x - y < 0]))
            2
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
                    e2 = (  # case: (x _ y) ?= (x _ y)
                        left._rel_equal1(right) or
                        # case: (x == y) ?= (y == x)
                        #       (x != y) ?= (y != x)
                        ((l.operator() == operator.eq or
                          l.operator() == operator.ne) and
                         left._rel_equal2(right)))
                else:
                    e2 = (  # case: (x < y)  ?= (y > x)  (or vice versa)
                            #       (x <= y) ?= (y >= x) (or vice versa)
                        ((l.operator() == operator.lt and
                          r.operator() == operator.gt) or
                         (l.operator() == operator.gt and
                          r.operator() == operator.lt) or
                         (l.operator() == operator.le and
                          r.operator() == operator.ge) or
                         (l.operator() == operator.ge and
                          r.operator() == operator.le)) and
                        left._rel_equal2(right))
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
            e = (l._gobj < r._gobj)
        elif op == Py_EQ:
            e = (l._gobj == r._gobj)
        elif op == Py_GT:
            e = (l._gobj > r._gobj)
        elif op == Py_LE:
            e = (l._gobj <= r._gobj)
        elif op == Py_NE:
            e = (l._gobj != r._gobj)
        elif op == Py_GE:
            e = (l._gobj >= r._gobj)
        else:
            raise TypeError
        return new_Expression_from_GEx(l._parent, e)

    def _test_nonzero_equal(self, **options):
        r"""
        Do not perform tests for the operators ``==`` and ``!=``.

        EXAMPLES:

        The operator ``==`` has a special meaning for a symbolic expression::

            sage: x == 0
            x == 0

        In particular, it does not return a boolean, so the following check does
        not hold anymore::

            sage: (not x) == (x != 0)
            False

        TESTS::

            sage: x._test_nonzero_equal()
        """
        pass

    def assume(self):
        r"""
        Assume that this equation holds. This is relevant for symbolic
        integration, among other things.

        EXAMPLES: We call the assume method to assume that `x>2`::

            sage: (x > 2).assume()

        ``bool`` returns ``True`` below if the inequality is *definitely* known
        to be true.

        ::

            sage: bool(x > 0)
            True
            sage: bool(x < 0)
            False

        This may or may not be True, so ``bool`` returns False::

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
            -75/8*sqrt(c^2)*arcsin(sqrt(c^2)*v/c^2) + 83/8*v/sqrt(-v^2/c^2 + 1) - 17/8*v^3/(c^2*sqrt(-v^2/c^2 + 1)) - 1/4*v^5/(c^4*sqrt(-v^2/c^2 + 1))
            sage: forget()
        """
        from sage.symbolic.assumptions import _assumptions
        from sage.calculus.calculus import maxima
        if not self.is_relational():
            raise TypeError("self (=%s) must be a relational expression" % self)
        if self not in _assumptions:
            m = self._maxima_init_assume_()
            s = maxima.assume(m)
            pynac_assume_rel(self._gobj)
            if str(s._sage_()[0]) in ['meaningless','inconsistent','redundant']:
                raise ValueError("Assumption is %s" % str(s._sage_()[0]))
            _assumptions[self] = True

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

        Check if :issue:`7507` is fixed::

            sage: forget()
            sage: n = var('n')
            sage: foo = sin((-1)*n*pi)
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
            raise TypeError("self (=%s) must be a relational expression" % self)
        pynac_forget_rel(self._gobj)
        m = self._maxima_init_assume_()
        maxima.forget(m)
        try:
            del _assumptions[self]
        except KeyError:
            pass

    def _maxima_init_assume_(self):
        """
        Return string that when evaluated in Maxima defines the assumption
        determined by this expression.

        EXAMPLES::

            sage: f = x+2 > sqrt(3)
            sage: f._maxima_init_assume_()
            '((_SAGE_VAR_x)+(2))>((3)^(1/2))'
        """
        from sage.calculus.calculus import maxima

        l = self.lhs()._assume_str()
        r = self.rhs()._assume_str()
        op = self.operator()
        if op is operator.eq:
            m = 'equal(%s, %s)' % (l, r)
        elif op is operator.ne:
            m = 'notequal(%s, %s)' % (l, r)
        else:
            m = '(%s)%s(%s)' % (l, maxima._relation_symbols()[op], r)
        return m

    def _assume_str(self):
        """
        TESTS::

            sage: x = var('x')
            sage: x._assume_str()
            '_SAGE_VAR_x'
            sage: y = function('y')(x)
            sage: y._assume_str()
            'y'
            sage: abs(x)._assume_str()
            'abs(_SAGE_VAR_x)'
        """
        # if this is a function with a single argument which is a symbol, i.e.
        # this is of the form f(x), we pass the string 'f > 0'
        if is_a_function(self._gobj) and self.nops() == 1 and \
           is_a_symbol(self._gobj.op(0)):
            op = self.operator()
            # check if op is a user defined function, for builtin
            # functions like abs() we still need to pass 'abs(x) > 0'
            from sage.symbolic.function import SymbolicFunction
            if isinstance(op, SymbolicFunction):
                return self.operator().name()
        return self._maxima_init_()

    def decl_assume(self, decl):
        """
        TESTS::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'real')
            sage: x.is_real()
            False
            sage: x.decl_assume(decl._assumption)
            sage: x.is_real()
            True
        """
        pynac_assume_gdecl(self._gobj, str_to_bytes(decl))

    def decl_forget(self, decl):
        """
        TESTS::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'integer')
            sage: x.is_integer()
            False
            sage: x.decl_assume(decl._assumption)
            sage: x.is_integer()
            True
            sage: x.decl_forget(decl._assumption)
            sage: x.is_integer()
            False
        """
        pynac_forget_gdecl(self._gobj, str_to_bytes(decl))

    def has_wild(self):
        """
        Return ``True`` if this expression contains a wildcard.

        EXAMPLES::

            sage: (1 + x^2).has_wild()
            False
            sage: (SR.wild(0) + x^2).has_wild()
            True
            sage: SR.wild(0).has_wild()
            True
        """
        return haswild(self._gobj)

    def is_algebraic(self):
        """
        Return ``True`` if this expression is known to be algebraic.

        EXAMPLES::

            sage: sqrt(2).is_algebraic()
            True
            sage: (5*sqrt(2)).is_algebraic()
            True
            sage: (sqrt(2) + 2^(1/3) - 1).is_algebraic()
            True
            sage: (I*golden_ratio + sqrt(2)).is_algebraic()
            True
            sage: (sqrt(2) + pi).is_algebraic()
            False
            sage: SR(QQ(2/3)).is_algebraic()
            True
            sage: SR(1.2).is_algebraic()
            False

            sage: complex_root_of(x^3 - x^2 - x - 1, 0).is_algebraic()
            True
        """
        from sage.rings.qqbar import QQbar
        try:
            _ = QQbar(self)
        except (TypeError, ValueError, NotImplementedError):
            return False
        return True

    def is_rational_expression(self):
        """
        Return ``True`` if this expression if a rational expression, i.e.,
        a quotient of polynomials.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: ((x + y + z)/(1 + x^2)).is_rational_expression()
            True
            sage: ((1 + x + y)^10).is_rational_expression()
            True
            sage: ((1/x + z)^5 - 1).is_rational_expression()
            True
            sage: (1/(x + y)).is_rational_expression()
            True
            sage: (exp(x) + 1).is_rational_expression()
            False
            sage: (sin(x*y) + z^3).is_rational_expression()
            False
            sage: (exp(x) + exp(-x)).is_rational_expression()
            False
        """
        return all(part.is_polynomial(v)
                   for part in (self.numerator(), self.denominator())
                   for v in part.variables())

    def is_real(self):
        """
        Return ``True`` if this expression is known to be a real number.

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
            sage: t2 = SR.symbol("t2", domain='positive')
            sage: (t1**t2).is_real()
            True
            sage: (t0*x).is_real()
            False
            sage: (t0^t1).is_real()
            False
            sage: (t1^t2).is_real()
            True
            sage: gamma(pi).is_real()
            True
            sage: cosh(-3).is_real()
            True
            sage: cos(exp(-3) + log(2)).is_real()
            True
            sage: gamma(t1).is_real()
            True
            sage: (x^pi).is_real()
            False
            sage: (cos(exp(t0) + log(t1))^8).is_real()
            True
            sage: cos(I + 1).is_real()
            False
            sage: sin(2 - I).is_real()
            False
            sage: (2^t0).is_real()
            True

        The following is real, but we cannot deduce that.::

            sage: (x*x.conjugate()).is_real()
            False

        Assumption of real has the same effect as setting the domain::

            sage: forget()
            sage: assume(x, 'real')
            sage: x.is_real()
            True
            sage: cosh(x).is_real()
            True
            sage: forget()

        The real domain is also set with the integer domain::

            sage: SR.var('x', domain='integer').is_real()
            True

        TESTS:

        Check that :issue:`23093` is fixed::

            sage: sqrt(-2).is_real()
            False
        """
        return self._gobj.info(info_real)

    def is_positive(self):
        """
        Return ``True`` if this expression is known to be positive.

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

        ::

            sage: forget()
            sage: assume(x>0)
            sage: x.is_positive()
            True
            sage: cosh(x).is_positive()
            True
            sage: f = function('f')(x)
            sage: assume(f>0)
            sage: f.is_positive()
            True
            sage: forget()

        TESTS:

        Check if :issue:`18630` is fixed::

            sage: (log(1/2)).is_negative()
            True
            sage: e.is_positive()
            True
            sage: (e+1).is_positive()
            True
            sage: (2*e).is_positive()
            True
            sage: (e^3).is_positive()
            True

        ::

            sage: cosh(x).is_positive()
            False
            sage: cosh(real(x)).is_positive()
            True
            sage: (cosh(real(x))^2).is_positive()
            True
            sage: ((real(x))^2).is_positive()
            False
            sage: gamma(x^2).is_positive()
            False
            sage: gamma(x^2+1).is_positive()
            False
            sage: gamma(cosh(real(x))).is_positive()
            True
            sage: (real(x)^2).is_positive()
            False
            sage: (real(x)^2+1).is_positive()
            True
            sage: (abs(x)^2+1).is_positive()
            True
            sage: gamma(real(x)^2+1).is_positive()
            True
            sage: cos(I + 1).is_positive()
            False
            sage: sin(2 - I).is_positive()
            False

        ::

            sage: (log(1/3) * log(1/2)).is_positive()
            True
            sage: log((2**500+1)/2**500).is_positive()
            True
            sage: log(2*500/(2**500-1)).is_negative()
            True
            sage: ((-pi^(1/5))^2).is_positive()
            True
            sage: (pi^2).is_positive()
            True
            sage: ((-pi)^2).is_positive()
            True
        """
        return self._gobj.info(info_positive)

    def is_negative(self):
        """
        Return ``True`` if this expression is known to be negative.

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

        Assumptions on symbols are handled correctly::

            sage: y = var('y')
            sage: assume(y < 0)
            sage: y.is_positive()
            False
            sage: y.is_negative()
            True
            sage: forget()
        """
        return self._gobj.info(info_negative)

    def is_integer(self):
        """
        Return ``True`` if this expression is known to be an integer.

        EXAMPLES::

            sage: SR(5).is_integer()
            True

        TESTS:

        Check that integer variables are recognized (:issue:`18921`)::

            sage: _ = var('n', domain='integer')
            sage: n.is_integer()
            True

        Assumption of integer has the same effect as setting the domain::

            sage: forget()
            sage: assume(x, 'integer')
            sage: x.is_integer()
            True
            sage: forget()
        """
        return self._gobj.info(info_integer)

    def is_symbol(self):
        """
        Return ``True`` if this symbolic expression consists of only a symbol, i.e.,
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

    def _is_registered_constant_(self):
        """
        Return ``True`` if this symbolic expression is internally represented as
        a constant.

        This function is intended to provide an interface to query the internal
        representation of the expression. In this sense, the word ``constant``
        does not reflect the mathematical properties of the expression.
        Expressions which have no variables may return ``False``.

        EXAMPLES::

            sage: pi._is_registered_constant_()
            True
            sage: x._is_registered_constant_()
            False
            sage: SR(1)._is_registered_constant_()
            False

        Note that the complex I is not a constant::

            sage: SR(I)._is_registered_constant_()
            False
            sage: SR(I).is_numeric()
            True
        """
        return is_a_constant(self._gobj)

    def is_constant(self):
        """
        Return whether this symbolic expression is a constant.

        A symbolic expression is constant if it does not contain
        any variables.

        EXAMPLES::

            sage: pi.is_constant()
            True
            sage: SR(1).is_constant()
            True
            sage: SR(2).is_constant()
            True
            sage: log(2).is_constant()
            True
            sage: SR(I).is_constant()
            True
            sage: x.is_constant()
            False

        TESTS::

            sage: P.<p> = ZZ[]
            sage: SR(42).is_constant() == P(2).is_constant()
            True
        """
        return not self.variables()

    def is_numeric(self):
        """
        A Pynac numeric is an object you can do arithmetic with
        that is not a symbolic variable, function, or constant.
        Return ``True`` if this expression only consists of a numeric object.

        EXAMPLES::

            sage: SR(1).is_numeric()
            True
            sage: x.is_numeric()
            False
            sage: pi.is_numeric()
            False
            sage: sin(x).is_numeric()
            False
        """
        return is_a_numeric(self._gobj)

    def is_terminating_series(self):
        """
        Return ``True`` if ``self`` is a series without order term.

        A series is terminating if it can be represented exactly,
        without requiring an order term. You can explicitly
        request terminating series by setting the order to
        positive infinity.

        OUTPUT: boolean; whether ``self`` was constructed by :meth:`series`
        and has no order term

        EXAMPLES::

            sage: (x^5+x^2+1).series(x, +oo)
            1 + 1*x^2 + 1*x^5
            sage: (x^5+x^2+1).series(x,+oo).is_terminating_series()
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
        return False

    cpdef bint is_polynomial(self, var) noexcept:
        """
        Return ``True`` if ``self`` is a polynomial in the given variable.

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

        Check if we can handle derivatives. :issue:`6523`::

            sage: f(x) = function('f')(x)
            sage: f(x).diff(x).is_zero()
            False

        Check if :issue:`11352` is fixed::

            sage: el = -1/2*(2*x^2 - sqrt(2*x - 1)*sqrt(2*x + 1) - 1)
            sage: el.is_polynomial(x)
            False

        Check that negative exponents are handled (:issue:`15304`)::

            sage: y = var('y')
            sage: (y/x).is_polynomial(x)
            False
        """
        cdef Expression symbol0 = self.coerce_in(var)
        sig_on()
        try:
            return self._gobj.is_polynomial(symbol0._gobj)
        finally:
            sig_off()

    cpdef bint is_relational(self) noexcept:
        r"""
        Return ``True`` if ``self`` is a relational expression.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.is_relational()
            True
            sage: sin(x).is_relational()
            False
        """
        return is_a_relational(self._gobj)

    def is_exact(self):
        r"""
        Return ``True`` if this expression only contains exact numerical coefficients.

        EXAMPLES::

            sage: x, y = var('x, y')
            sage: (x+y-1).is_exact()
            True
            sage: (x+y-1.9).is_exact()
            False
            sage: x.is_exact()
            True
            sage: pi.is_exact()
            True
            sage: (sqrt(x-y) - 2*x + 1).is_exact()
            True
            sage: ((x-y)^0.5 - 2*x + 1).is_exact()
            False

        TESTS::

            sage: (sin(x*cos(2*x*pi)) - 10*y^3 - 1/(x+4)).is_exact()
            True
            sage: (sin(x*cos(2.0*x*pi)) - 10*y^3 - 1/(x+4)).is_exact()
            False
            sage: SR(42).is_exact()
            True
            sage: SR(42.01).is_exact()
            False
            sage: SR(I).is_exact()
            True
            sage: (x-I).is_exact()
            True
            sage: (x-CC(0,1)).is_exact()
            False
        """
        # generator over all numerical elements in the subexpression tree of expr
        def numelems_gen(expr):
            if expr.is_numeric():
                yield expr
            elif expr.operator() is not None:
                for op in expr.operands():
                    if op.is_numeric():
                        yield op
                    else:
                        yield from numelems_gen(op)
        # stop at the first inexact number in the subexpression tree of self,
        # and if there is no such element, then self is exact
        return all(nelem.pyobject().base_ring().is_exact()
                   for nelem in numelems_gen(self))

    cpdef bint is_infinity(self) noexcept:
        """
        Return ``True`` if ``self`` is an infinite expression.

        EXAMPLES::

            sage: SR(oo).is_infinity()
            True
            sage: x.is_infinity()
            False
        """
        return is_a_infinity(self._gobj)

    cpdef bint is_positive_infinity(self) noexcept:
        """
        Return ``True`` if ``self`` is a positive infinite expression.

        EXAMPLES::

            sage: SR(oo).is_positive_infinity()
            True
            sage: SR(-oo).is_positive_infinity()
            False
            sage: x.is_infinity()
            False
        """
        return is_a_infinity(self._gobj) and self._gobj.info(info_positive)

    cpdef bint is_negative_infinity(self) noexcept:
        """
        Return ``True`` if ``self`` is a negative infinite expression.

        EXAMPLES::

            sage: SR(oo).is_negative_infinity()
            False
            sage: SR(-oo).is_negative_infinity()
            True
            sage: x.is_negative_infinity()
            False
        """
        return is_a_infinity(self._gobj) and self._gobj.info(info_negative)

    def is_square(self):
        """
        Return ``True`` if ``self`` is the square of another symbolic expression.

        This is ``True`` for all constant, non-relational expressions
        (containing no variables or comparison), and not implemented
        otherwise.

        EXAMPLES::

            sage: SR(4).is_square()
            True
            sage: SR(5).is_square()
            True
            sage: pi.is_square()
            True
            sage: x.is_square()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square() not implemented for non-constant
            or relational elements of Symbolic Ring
            sage: r = SR(4) == SR(5)
            sage: r.is_square()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square() not implemented for non-constant
            or relational elements of Symbolic Ring
        """
        if self.is_constant() and not self.is_relational():
            # The square root of any "number" is in SR... just call
            # sqrt() on it.
            return True

        try:
            obj = self.pyobject()
        except TypeError as e:
            raise NotImplementedError("is_square() not implemented for non-constant or relational elements of Symbolic Ring")

        return obj.is_square()

    def is_callable(self):
        r"""
        Return ``True`` if ``self`` is a callable symbolic expression.

        EXAMPLES::

            sage: var('a x y z')
            (a, x, y, z)
            sage: f(x, y) = a + 2*x + 3*y + z
            sage: f.is_callable()
            True
            sage: (a+2*x).is_callable()
            False
        """
        return isinstance(self.parent(), sage.rings.abc.CallableSymbolicExpressionRing)

    def left_hand_side(self):
        """
        If ``self`` is a relational expression, return the left hand side
        of the relation.  Otherwise, raise a :exc:`ValueError`.

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
            raise ValueError("self must be a relational expression")
        return new_Expression_from_GEx(self._parent, self._gobj.lhs())

    lhs = left = left_hand_side

    def right_hand_side(self):
        """
        If ``self`` is a relational expression, return the right hand side
        of the relation.  Otherwise, raise a :exc:`ValueError`.

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
            raise ValueError("self must be a relation")
        return new_Expression_from_GEx(self._parent, self._gobj.rhs())

    rhs = right = right_hand_side

    def is_trivially_equal(self, other):
        """
        Check if this expression is trivially equal to the argument
        expression, without any simplification.

        Note that the expressions may still be subject to immediate
        evaluation.

        This method is intended to be used in library code where trying to
        obtain a mathematically correct result by applying potentially
        expensive rewrite rules is not desirable.

        EXAMPLES::

            sage: (x^2).is_trivially_equal(x^2)
            True
            sage: ((x+1)^2 - 2*x - 1).is_trivially_equal(x^2)
            False
            sage: (x*(x+1)).is_trivially_equal((x+1)*x)
            True
            sage: (x^2 + x).is_trivially_equal((x+1)*x)
            False
            sage: ((x+1)*(x+1)).is_trivially_equal((x+1)^2)
            True
            sage: (x^2 + 2*x + 1).is_trivially_equal((x+1)^2)
            False
            sage: (x^-1).is_trivially_equal(1/x)
            True
            sage: (x/x^2).is_trivially_equal(1/x)
            True
            sage: ((x^2+x) / (x+1)).is_trivially_equal(1/x)
            False

        TESTS:

        Make sure Python objects work as argument too::

            sage: x = SR(1/2)
            sage: x.is_trivially_equal(QQbar(1/2))
            True
        """
        from sage.symbolic.ring import SR
        cdef Expression _other = <Expression>(SR(other))
        sig_on()
        try:
            return self._gobj.is_equal(_other._gobj)
        finally:
            sig_off()

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
            sage: t = pi + x*pi + (pi - 1 - x)*pi - pi^2
            sage: t.is_zero()
            True
            sage: u = sin(x)^2 + cos(x)^2 - 1
            sage: u.is_trivial_zero()
            False
            sage: u.is_zero()
            True
        """
        return self._gobj.is_zero()

    def __bool__(self):
        """
        Return ``True`` unless this symbolic expression can be shown by Sage
        to be zero.  Note that deciding if an expression is zero is
        undecidable in general.

        EXAMPLES::

            sage: x = var('x')
            sage: forget()
            sage: bool(SR(0))
            False
            sage: bool(SR(1))
            True
            sage: assert(abs(x))
            sage: assert(not x/x - 1)

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
            sage: assert((x-1)^2 == x^2 - 2*x + 1)
            sage: assert(((x-1)^2 == x^2 - 2*x + 1).expand())
            sage: assert(not ((x-1)^2 == x^2 - 2*x + 3).expand())
            sage: assert(2 + x < 3 + x)
            sage: assert(not 2 + x < 1 + x)
            sage: assert(2 + x > 1 + x)
            sage: assert(not 1 + x > 1 + x)
            sage: assert(1 + x >= 1 + x)
            sage: assert(not 1 + x < 1 + x)
            sage: assert(1 + x <= 1 + x)
            sage: assert(not 1 + x^2 != 1 + x*x)
            sage: assert(1 + x^2 != 2 + x*x)
            sage: assert(SR(oo) == SR(oo))
            sage: assert(not -SR(oo) == SR(oo))
            sage: assert(-SR(oo) != SR(oo))

        Next, tests to ensure assumptions are correctly used::

            sage: x, y, z = var('x, y, z')
            sage: assume(x >= y, y >= z, z >= x)
            sage: assert(x == z)
            sage: assert(not z < x)
            sage: assert(not z > y)
            sage: assert(y == z)
            sage: assert(y <= z)
            sage: forget()
            sage: assume(x >= 1, x <= 1)
            sage: assert(x == 1)
            sage: assert(not x != 1)
            sage: assert(not x > 1)
            sage: forget()
            sage: assume(x > 0)
            sage: assert(not x == 0)
            sage: assert(x != 0)

        We cannot return undecidable or throw an exception
        at the moment so ``False`` is returned for unknown
        outcomes.

        ::

            sage: assert(not x == 1)
            sage: assert(not x != 1)
            sage: forget()
            sage: assume(x>y)
            sage: assert(not x==y)
            sage: assert(x != y) # The same comment as above applies here as well
            sage: forget()

        Comparisons of infinities::

            sage: assert( (1+I)*oo == (2+2*I)*oo )
            sage: assert( SR(unsigned_infinity) == SR(unsigned_infinity) )
            sage: assert( SR(I*oo) == I*oo )
            sage: assert( SR(-oo) <= SR(oo) )
            sage: assert( SR(oo) >= SR(-oo) )
            sage: assert( SR(oo) != SR(-oo) )
            sage: assert( sqrt(2)*oo != I*oo )

        The expression may be zero with integers but is not
        when in the complex domain (:issue:`15571`)::

            sage: a,x = var('a,x')
            sage: assume(a, 'integer')
            sage: assume(x, 'integer')
            sage: expr = a^(4*x) - (a^4)^x
            sage: expr.is_zero()
            True
            sage: forget()
            sage: assume(a, 'complex')
            sage: assume(x, 'complex')
            sage: expr.is_zero()
            False
            sage: forget()

        Check that :issue:`13326` is fixed::

            sage: assert(log(2)*Infinity == Infinity)

        More checks for comparisons with infinity (see :issue:`12967`)::

            sage: assert(SR(oo) > 5)
            sage: assert(5 < SR(oo))
            sage: assert(SR(2) < Infinity)
            sage: assert(pi < Infinity)
            sage: assert(not pi>Infinity)
            sage: assert(2*pi < Infinity)
            sage: assert(SR(pi) < SR(Infinity))
            sage: assert(sqrt(2) < oo)
            sage: assert(log(2) < oo)
            sage: assert(e < oo)
            sage: assert(e+pi < oo)
            sage: assert(e^pi < oo)
            sage: assert(not SR(2) < -oo)
            sage: assert(SR(2) > -oo)
            sage: assert(exp(2) > -oo)
            sage: assert(SR(oo) > sqrt(2))
            sage: assert(sqrt(2) < SR(oo))
            sage: assert(SR(-oo) < sqrt(2))
            sage: assert(sqrt(2) > SR(-oo))

        Check that :issue:`18360` is fixed::

            sage: f(x) = matrix()
            sage: bool(f(x) - f(x) == 0)
            True

        Check that :issue:`22857` is fixed::

            sage: a, b = var('a b', domain='positive')
            sage: bool((a-b)*b == 0)
            False

        Check that :issue:`24658` is fixed::

            sage: val = pi - 2286635172367940241408/1029347477390786609545*sqrt(2)
            sage: bool(val>0)
            False

        Check that :issue:`34341` is fixed::

            sage: bool(x^2 + 2*x + 1 != (x + 1)^2)
            False

        Check that :issue:`16031` is fixed::

            sage: expr = reduce(lambda u, v: 1/u -v, [1/pi] + list(continued_fraction(pi)[:20]))
            sage: expr.is_zero()
            False
        """
        if self.is_relational():
            # constants are wrappers around Sage objects, compare directly
            if is_a_constant(self._gobj.lhs()) and is_a_constant(self._gobj.rhs()):
                return self.operator()(self.lhs().pyobject(), self.rhs().pyobject())
            pynac_result = decide_relational(self._gobj)
            if pynac_result == relational_undecidable:
                raise ValueError('undecidable relation: ' + repr(self))

            # pynac is guaranteed to give the correct answer for comparing infinities
            if is_a_infinity(self._gobj.lhs()) or is_a_infinity(self._gobj.rhs()):
                return pynac_result == relational_true

            if pynac_result == relational_true:
                if self.operator() == operator.ne:
                    # this hack is necessary to catch the case where the
                    # operator is != but is False because of assumptions made
                    m = self._maxima_()
                    s = m.parent()._eval_line('is (notequal(%s,%s))' % (repr(m.lhs()),repr(m.rhs())))
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
                        except AttributeError:  # if we have a GenericDeclaration
                            assumption_var_list.append((eqn._var,))
                    assumption_vars = set(sum(assumption_var_list, ()))
                    if set(vars).intersection(assumption_vars):
                        need_assumptions = True

            # Use interval fields to try and falsify the relation
            if not need_assumptions:
                if pynac_result == relational_notimplemented and self.operator() == operator.ne:
                    return not (self.lhs()-self.rhs()).is_zero()
                res = self.test_relation()
                if res in (True, False):
                    return res
                res = self.operator()((self.lhs()-self.rhs()).simplify_full(), 0).test_relation()
                if res in (True, False):
                    return res

            # we really have to do some work here...
            # I really don't like calling Maxima to test equality.  It
            # is SUPER SUPER SLOW, and it has all the problem
            # associated with different semantics, different
            # precision, etc., that can lead to subtle bugs.  Also, a
            # lot of basic Sage objects can't be put into maxima.
            from sage.symbolic.relation import test_relation_maxima
            if self.variables():
                return test_relation_maxima(self)
            else:
                return False

        self_is_zero = self._gobj.is_zero()
        if self_is_zero:
            return False
        else:
            return not bool(self == self._parent.zero())

    def test_relation(self, int ntests=20, domain=None, proof=True):
        """
        Test this relation at several random values, attempting to find
        a contradiction. If this relation has no variables, it will also
        test this relation after casting into the domain.

        Because the interval fields never return false positives, we can be
        assured that if ``True`` or ``False`` is returned (and proof is
        ``False``) then the answer is correct.

        INPUT:

        - ``ntests`` -- (default: 20) the number of iterations to run
        - ``domain`` -- (optional) the domain from which to draw the random
          values defaults to ``CIF`` for equality testing and ``RIF`` for
          order testing
        - ``proof`` -- (default: ``True``) if ``False`` and the domain is an
          interval field, regard overlapping (potentially equal) intervals as
          equal, and return ``True`` if all tests succeeded.

        OUTPUT: boolean or ``NotImplemented``, meaning

        - ``True`` -- this relation holds in the domain and has no variables

        - ``False`` -- a contradiction was found

        - ``NotImplemented`` -- no contradiction found

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
            sage: (cot(-x) == -cot(x)).test_relation()
            NotImplemented

        Check that :issue:`18896` is fixed::

            sage: m = 540579833922455191419978421211010409605356811833049025*sqrt(1/2)
            sage: m1 = 382247666339265723780973363167714496025733124557617743
            sage: (m == m1).test_relation(domain=QQbar)
            False
            sage: (m == m1).test_relation()
            False

        Try the examples from :issue:`31424` and :issue:`31665`::

            sage: k = 26
            sage: bool(2/(2*pi)^(2*k) <= abs(bernoulli(2*k)/factorial(2*k)))
            True
            sage: t = log(17179815199/17179869184) + 727717503781781876485802\
            ....: 752874818120860129694543334299450155913077668355/2315841784\
            ....: 74632390847141970017375815706539969331281128078915168015826\
            ....: 259279872
            sage: v = -53985/17179869184
            sage: bool(abs(t) < 1.213*2^-56*v^4)
            True
        """
        cdef int k, eq_count = 0
        cdef bint is_interval
        if not self.is_relational():
            raise ValueError("self must be a relation")
        cdef operators op = relational_operator(self._gobj)
        if domain is None:
            is_interval = True
            if self.lhs().is_algebraic() and self.rhs().is_algebraic():
                if op == equal or op == not_equal:
                    from sage.rings.qqbar import QQbar as domain
                else:
                    from sage.rings.qqbar import AA as domain
            else:
                if op == equal or op == not_equal:
                    from sage.rings.complex_interval_field import \
                        ComplexIntervalField
                    # We don't want to be in the business of trying to
                    # ensure enough precision to solve EVERY problem,
                    # but since there are two real-life examples in
                    # Github issues 31424 and 31665 that are aided by
                    # a bump, we reluctantly enter that game.
                    domain = ComplexIntervalField(128)
                else:
                    from sage.rings.real_mpfi import RealIntervalField
                    domain = RealIntervalField(128)
        else:
            is_interval = isinstance(domain, (sage.rings.abc.RealIntervalField,
                                              sage.rings.abc.ComplexIntervalField,
                                              sage.rings.abc.AlgebraicField_common))
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
            except (TypeError, ValueError, ArithmeticError) as ex:
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
                except (TypeError, ValueError, ArithmeticError, AttributeError) as ex:
                    errors += 1
                    if k == errors > 3 and isinstance(domain, sage.rings.abc.ComplexIntervalField):
                        from sage.rings.real_mpfi import RIF
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
        Return the negated version of ``self``.

        This is the relation that is ``False`` iff ``self`` is ``True``.

        EXAMPLES::

            sage: (x < 5).negation()
            x >= 5
            sage: (x == sin(3)).negation()
            x != sin(3)
            sage: (2*x >= sqrt(2)).negation()
            2*x < sqrt(2)
        """
        if not self.is_relational():
            raise ValueError("self must be a relation")
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
        Return ``True`` if this relation is violated by the given variable assignment(s).

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
        Return ``True`` if this expression is a unit of the symbolic ring.

        Note that a proof may be attempted to get the result. To avoid
        this use ``(ex-1).is_trivial_zero()``.

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
        try:
            return <Expression?>z
        except TypeError:
            return self._parent.coerce(z)

    cpdef _add_(left, right):
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
                x = relational(left._gobj.lhs() + _right._gobj.lhs(),
                               left._gobj.rhs() + _right._gobj.rhs(),
                               op)
            else:
                x = relational(left._gobj.lhs() + _right._gobj,
                               left._gobj.rhs() + _right._gobj,
                               relational_operator(left._gobj))
        elif is_a_relational(_right._gobj):
            x = relational(left._gobj + _right._gobj.lhs(),
                           left._gobj + _right._gobj.rhs(),
                           relational_operator(_right._gobj))
        else:
            x = left._gobj + _right._gobj
        return new_Expression_from_GEx(left._parent, x)

    cpdef _sub_(left, right):
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
                x = relational(left._gobj.lhs() - _right._gobj.lhs(),
                               left._gobj.rhs() - _right._gobj.rhs(),
                               op)
            else:
                x = relational(left._gobj.lhs() - _right._gobj,
                               left._gobj.rhs() - _right._gobj,
                               relational_operator(left._gobj))
        elif is_a_relational(_right._gobj):
            x = relational(left._gobj - _right._gobj.lhs(),
                           left._gobj - _right._gobj.rhs(),
                           relational_operator(_right._gobj))
        else:
            x = left._gobj - _right._gobj
        return new_Expression_from_GEx(left._parent, x)

    cpdef _mul_(left, right):
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

            sage: var('z')
            z
            sage: 3*(x+y)/z
            3*(x + y)/z
            sage: (-x+z)*(3*x-3*z)
            -3*(x - z)^2

            # check if comparison of constant terms in Pynac add objects work
            sage: (y-1)*(y-2)
            (y - 1)*(y - 2)

        Check if Pynac can compute inverses of Python longs (:issue:`13107`)::

            sage: SR(4L)*SR(2L)^(-1)
            2

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
            sage: A = exp(I*pi/5, hold=True)
            sage: t = A*A; t
            1/4*sqrt(5) + 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4
            sage: A^5
            -1
            sage: b = -x*A; c = b*b; c
            1/4*x^2*(sqrt(5) + I*sqrt(2*sqrt(5) + 10) - 1)
            sage: exp(x)^I*exp(z)^(2.5)
            e^(I*x + 2.50000000000000*z)

        ::

            sage: x*oo
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity * f(x) encountered.
            sage: x*unsigned_infinity
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: infinity * f(x) encountered.

            sage: SR(oo)*SR(oo)
            +Infinity
            sage: SR(-oo)*SR(oo)
            -Infinity
            sage: SR(oo)*SR(-oo)
            -Infinity
            sage: SR(unsigned_infinity)*SR(oo)
            Infinity

        Check if we are returning informative error messages in case of
        nonsensical arithmetic (:issue:`10960`:, :issue:`13739` and
        :issue:`24072`)::

            sage: GF(5)(3) * SR.var('x')
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Finite Field of size 5' and 'Symbolic Ring'

            sage: b = polygen(FiniteField(9), 'b')
            sage: SR('I') * b
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations

        Check that :issue:`18360` is fixed::

            sage: f(x) = matrix()
            sage: f(x)*1
            []

        Check that floating point numbers +/- 1.0 are treated
        differently from integers +/- 1 (:issue:`12257`)::

            sage: (1*x).operator()
            sage: (1.0*x).operator()
            <function mul_vararg...
            sage: 1.0 * pi
            1.00000000000000*pi
            sage: 1.000000*(x+2)
            1.00000000000000*x + 2.00000000000000
            sage: -1.0*x
            -1.00000000000000*x
            sage: -1.0/x
            -1.00000000000000/x
            sage: (-1.0*x)*(1.0/x)
            -1.00000000000000
            sage: sin(1.0*pi)
            sin(1.00000000000000*pi)

        Check that infinities multiply correctly (:issue:`23427`)::

            sage: SR(-oo) * SR(-oo)
            +Infinity
            sage: SR(-oo) * SR(oo)
            -Infinity
            sage: SR(-oo) * SR(unsigned_infinity)
            Infinity
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(left._gobj.lhs() * _right._gobj.lhs(),
                               left._gobj.rhs() * _right._gobj.rhs(),
                               op)
            else:
                o = relational_operator(left._gobj)
                x = relational(left._gobj.lhs() * _right._gobj,
                               left._gobj.rhs() * _right._gobj,
                               o)
        elif is_a_relational(_right._gobj):
            o = relational_operator(_right._gobj)
            x = relational(left._gobj * _right._gobj.lhs(),
                           left._gobj * _right._gobj.rhs(),
                           o)
        else:
            x = left._gobj * _right._gobj
        return new_Expression_from_GEx(left._parent, x)

    cpdef _div_(left, right):
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
            ZeroDivisionError: symbolic division by zero

        Check if Pynac can compute divisions of Python longs (:issue:`13107`)::

            sage: SR(1L)/SR(2L)
            1/2
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        try:
            if is_a_relational(left._gobj):
                if is_a_relational(_right._gobj):
                    op = compatible_relation(relational_operator(left._gobj),
                                             relational_operator(_right._gobj))
                    x = relational(left._gobj.lhs() / _right._gobj.lhs(),
                                   left._gobj.rhs() / _right._gobj.rhs(),
                                   op)
                else:
                    o = relational_operator(left._gobj)
                    x = relational(left._gobj.lhs() / _right._gobj,
                                   left._gobj.rhs() / _right._gobj,
                                   o)
            elif is_a_relational(_right._gobj):
                o = relational_operator(_right._gobj)
                x = relational(left._gobj / _right._gobj.lhs(),
                               left._gobj / _right._gobj.rhs(),
                               o)
            else:
                x = left._gobj / _right._gobj
            return new_Expression_from_GEx(left._parent, x)
        except Exception as msg:
            # TODO: change this to maybe cleverly do something involving Cython C++ exception handling.
            # See http://docs.cython.org/docs/wrapping_CPlusPlus.html
            if 'division by zero' in str(msg):
                raise ZeroDivisionError("symbolic division by zero")
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
            sage: v1=var('v1'); a = (2*erf(2*v1*arcsech(1/2))/v1); ~a
            1/2*v1/erf(2*v1*arcsech(1/2))
        """
        return 1/self

    cpdef int _cmp_add(Expression left, Expression right) except -2:
        """
        Compare ``left`` and ``right`` in the print order.

        INPUT:

        - ``right`` -- a :class:`Expression` instance

        OUTPUT: boolean

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

        - ``right`` -- a :class:`Expression` instance

        OUTPUT: boolean

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

    cpdef _pow_(self, other):
        r"""
        Return ``self`` raised to the power ``other``.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: x._pow_(y)
            x^y
            sage: x^(3/5)
            x^(3/5)
            sage: x^sin(x)^cos(y)
            x^(sin(x)^cos(y))

        Immediate simplifications are applied::

            sage: x = SR.symbol('x', domain='real')
            sage: (x^3)^(1/3)
            x
            sage: (x^4)^(1/4)
            abs(x)
            sage: (x^8)^(1/4)
            x^2
            sage: (x^-4)^(1/4)
            1/abs(x)
            sage: (x^-8)^(1/4)
            x^(-2)
            sage: forget()

        TESTS::

            sage: (Mod(2,7)*x^2 + Mod(2,7))^7
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Ring of integers modulo 7' and 'Symbolic Ring'

        The leading coefficient in the result above is 1 since::

            sage: t = Mod(2,7); gcd(t, t)^7
            1
            sage: gcd(t,t).parent()
            Ring of integers modulo 7

        ::

            sage: k = GF(7)
            sage: f = expand((k(1)*x^5 + k(1)*x^2 + k(2))^7); f  # known bug
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
            sage: ex = exp(sqrt(x))^x; ex
            (e^sqrt(x))^x
            sage: latex(ex)
            \left(e^{\sqrt{x}}\right)^{x}

        Test simplification of powers involving the reciprocal
        logarithm of the (positive) base::

            sage: 2^(1/log(2))
            e
            sage: 2^(x/log(2))
            e^x
            sage: 2^(-x^2/2/log(2))
            1/e^(1/2*x^2)
            sage: x^(x/log(x))
            x^(x/log(x))
            sage: assume(x > 0)
            sage: x^(x/log(x))
            e^x
            sage: forget()

        Test base a Python numeric type::

            sage: int(2)^x
            2^x
            sage: float(2.3)^(x^3 - x^2 + 1/3)
            2.3^(x^3 - x^2 + 1/3)
            sage: complex(1,3)^(sqrt(2))
            (1+3j)^sqrt(2)

        Test complex numeric powers::

            sage: symI = SR(I)
            sage: symI^0.5
            0.707106781186548 + 0.707106781186547*I
            sage: (symI + 1) ^ (0.5 + symI)
            0.400667052375828 + 0.365310866736929*I
            sage: symI^symI
            I^I
            sage: symI^x
            I^x
            sage: symI^(1/2)
            sqrt(I)
            sage: symI^(2/3)
            I^(2/3)
            sage: 2^(1/2)
            sqrt(2)
            sage: (2*symI)^(1/2)
            sqrt(2*I)

        Test if we can take powers of elements of `\QQ(i)` (:issue:`8659`)::

            sage: t = QuadraticField(-1, 'I')(8)
            sage: t^(1/2)
            2*sqrt(2)
            sage: (t^2)^(1/4)
            2*4^(1/4)

        Test if we can compute inverses of Python longs (:issue:`13107`)::

            sage: SR(2L)^(-1)
            1/2

        Symbolic powers with ``None`` shouldn't crash (:issue:`17523`)::

            sage: None^pi
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for ** or pow(): 'NoneType' and 'sage.symbolic.expression.Expression'
            sage: sin(x)^None
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for ** or pow(): 'sage.symbolic.expression.Expression' and 'NoneType'

        Check that :issue:`18088` is fixed::

            sage: SR(0)^SR(0)
            1

        Check that floating point numbers +/- 1.0 are treated
        differently from integers +/- 1 (:issue:`12257`)::

            sage: (x^1).operator()
            sage: (x^1.0).operator()
            <built-in function pow>
            sage: x^1.0
            x^1.00000000000000
            sage: x^-1.0
            x^(-1.00000000000000)
            sage: 0^1.0
            0.000000000000000
            sage: exp(x)^1.0
            e^(1.00000000000000*x)

        Check that :issue:`23921` is resolved::

            sage: assume(SR.an_element() > 0)
            sage: A.<n> = AsymptoticRing('(SR_+)^n * n^SR', SR)
            sage: elem = SR(2)^n
            sage: (elem, elem.parent())
            (2^n, Asymptotic Ring <SR^n * n^SR> over Symbolic Ring)

        Check that pynac understands rational powers (:issue:`30446`,
        :issue:`28620`, :issue:`30304`, and :issue:`30786`)::

            sage: QQ((24*sqrt(3))^(100/50))==1728
            True
            sage: float((24*sqrt(3))^(100/51))
            1493.0092154...
            sage: t = ((1/10)*I/pi)^(3/2)
            sage: t^2
            -1/1000*I/pi^3
            sage: (2*pi)^QQ(2)
            4*pi^2
            sage: exp(-3*ln(-9*x)/3)
            -1/9/x

        Check that :issue:`31137` is also fixed::

            sage: _ = var('A, L, G, R, f, k, n, q, u, beta, gamma', domain='positive')
            sage: a = I*R^2*f^3*k*q*A*u
            sage: b = 2*pi*L*R^2*G*f^4*k^2*q - 2*pi*L*R^2*G*f^4*q - 2*pi*L*R^2*beta^2*G*q
            sage: c = (2*I*pi*L*R^2*beta*gamma*q + 2*I*pi*L*R*(beta + q))*G*f^3
            sage: d = 2*(pi*(beta^2 + 1)*L*R^2*q + pi*L*R*beta*gamma*q + pi*L*beta)*G*f^2
            sage: e = (-2*I*pi*L*R^2*beta*gamma*q - 2*I*pi*(beta^2*q + beta)*L*R)*G*f
            sage: expr = a / ((b + c + d + e)*n)
            sage: R = ((sqrt(expr.real()^2 + expr.imag()^2).factor())^2).factor()
            sage: Rs = R.subs(f = 2*beta)
            sage: len(str(Rs))
            520
        """
        cdef Expression nexp = <Expression>other
        cdef GEx x
        if is_a_relational(self._gobj):
            x = relational(g_pow(self._gobj.lhs(), nexp._gobj),
                           g_pow(self._gobj.rhs(), nexp._gobj),
                           relational_operator(self._gobj))
        else:
            x = g_pow(self._gobj, nexp._gobj)
        return new_Expression_from_GEx(self._parent, x)

    cpdef _pow_int(self, other):
        """
        TESTS::

            sage: var('x')
            x
            sage: cos(x)._pow_int(-3)
            cos(x)^(-3)
            sage: SR(-oo)._pow_int(2)
            +Infinity
            sage: pi._pow_int(4)
            pi^4
        """
        return self._pow_(self._parent(other))

    def derivative(self, *args):
        """
        Return the derivative of this expressions with respect to the
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global
        :meth:`~sage.calculus.functional.derivative` function for more
        details.

        .. SEEALSO::

            This is implemented in the ``_derivative`` method (see the
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

        Some expressions can't be cleanly differentiated by the
        chain rule::

            sage: _ = var('x', domain='real')
            sage: _ = var('w z')
            sage: (x^z).conjugate().diff(x)
            conjugate(x^(z - 1))*conjugate(z)
            sage: (w^z).conjugate().diff(w)
            w^(z - 1)*z*D[0](conjugate)(w^z)
            sage: atanh(x).real_part().diff(x)
            -1/(x^2 - 1)
            sage: atanh(x).imag_part().diff(x)
            0
            sage: atanh(w).real_part().diff(w)
            -D[0](real_part)(arctanh(w))/(w^2 - 1)
            sage: atanh(w).imag_part().diff(w)
            -D[0](imag_part)(arctanh(w))/(w^2 - 1)
            sage: abs(log(x)).diff(x)
            1/2*(conjugate(log(x))/x + log(x)/x)/abs(log(x))
            sage: abs(log(z)).diff(z)
            1/2*(conjugate(log(z))/z + log(z)/conjugate(z))/abs(log(z))
            sage: forget()

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
        Return the ``deg``-th (partial) derivative of ``self`` with respect to ``symb``.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: b = (x+y)^5
            sage: b._derivative(x, 2)
            20*(x + y)^3

            sage: foo = function('foo', nargs=2)
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

        Check if :issue:`6524` is fixed::

            sage: f = function('f')
            sage: f(x)*f(x).derivative(x)*f(x).derivative(x,2)
            f(x)*diff(f(x), x)*diff(f(x), x, x)
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
            elif self.is_callable():
                return self.gradient()
            else:
                raise ValueError("No differentiation variable specified.")
        if not isinstance(deg, (int, sage.rings.integer.Integer)) \
                or deg < 1:
            raise TypeError("argument deg should be an integer >= 1.")
        cdef Expression symbol = self.coerce_in(symb)
        if not is_a_symbol(symbol._gobj):
            raise TypeError("argument symb must be a symbol")
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

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.scalarfield.DiffScalarField.gradient`
            of scalar fields on Euclidean spaces (and more generally
            pseudo-Riemannian manifolds), in particular for computing the
            gradient in curvilinear coordinates.
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

    def series(self, symbol, order=None):
        r"""
        Return the power series expansion of ``self`` in terms of the
        given variable to the given order.

        INPUT:

        - ``symbol`` -- a symbolic variable or symbolic equality
          such as ``x == 5``; if an equality is given, the
          expansion is around the value on the right hand side
          of the equality
        - ``order`` -- integer; if nothing given, it is set
          to the global default (``20``), which can be changed
          using :func:`set_series_precision`

        OUTPUT: a power series

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
            3 + (-5)*x + (-sin(y))*x^2 + 1*x^3 + Order(x^4)
            sage: g.truncate()
            x^3 - x^2*sin(y) - 5*x + 3
            sage: g = f.series(x==1, 4); g
            (-sin(y) - 1) + (-2*sin(y) - 2)*(x - 1) + (-sin(y) + 3)*(x - 1)^2
            + 1*(x - 1)^3 + Order((x - 1)^4)
            sage: h = g.truncate(); h
            (x - 1)^3 - (x - 1)^2*(sin(y) - 3) - 2*(x - 1)*(sin(y) + 1) - sin(y) - 1
            sage: h.expand()
            x^3 - x^2*sin(y) - 5*x + 3

        We computer another series expansion of an analytic function::

            sage: f = sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x)
            1*x^(-1) + (-1/6)*x + ... + Order(x^20)
            sage: f.series(x==1,3)
            (sin(1)) + (cos(1) - 2*sin(1))*(x - 1) + (-2*cos(1) + 5/2*sin(1))*(x - 1)^2
            + Order((x - 1)^3)
            sage: f.series(x==1,3).truncate().expand()
            -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)

        Expressions formed by combining series can be expanded
        by applying series again::

            sage: (1/(1-x)).series(x, 3)+(1/(1+x)).series(x,3)
            (1 + 1*x + 1*x^2 + Order(x^3)) + (1 + (-1)*x + 1*x^2 + Order(x^3))
            sage: _.series(x,3)
            2 + 2*x^2 + Order(x^3)
            sage: (1/(1-x)).series(x, 3)*(1/(1+x)).series(x,3)
            (1 + 1*x + 1*x^2 + Order(x^3))*(1 + (-1)*x + 1*x^2 + Order(x^3))
            sage: _.series(x,3)
            1 + 1*x^2 + Order(x^3)

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

        Check if :issue:`8943` is fixed::

            sage: ((1+arctan(x))**(1/x)).series(x==0, 3)
            (e) + (-1/2*e)*x + (1/8*e)*x^2 + Order(x^3)

        Order may be negative::

            sage: f = sin(x)^(-2); f.series(x, -1)
            1*x^(-2) + Order(1/x)

        Check if changing global series precision does it right::

            sage: set_series_precision(3)
            sage: (1/(1-2*x)).series(x)
            1 + 2*x + 4*x^2 + Order(x^3)
            sage: set_series_precision(20)

        Check that :issue:`31645` is fixed::

            sage: (x^(-1) + 1).series(x,1)
            1*x^(-1) + 1 + Order(x)

        Check that :issue:`32115` is fixed::

            sage: exp(log(1+x)*(1/x)).series(x)
            (e) + (-1/2*e)*x + (11/24*e)*x^2 + (-7/16*e)*x^3 + (2447/5760*e)*x^4 + ...

        Check that :issue:`32640` is fixed::

            sage: ((1 - x)^-x).series(x, 8)
            1 + 1*x^2 + 1/2*x^3 + 5/6*x^4 + 3/4*x^5 + 33/40*x^6 + 5/6*x^7 + Order(x^8)
        """
        cdef Expression symbol0 = self.coerce_in(symbol)
        cdef GEx x
        cdef SymbolicSeries nex
        cdef int cprec
        cdef int options
        if order is infinity:
            options = 4
            order = None
        else:
            options = 2

        if order is None:
            from sage.misc.defaults import series_precision
            cprec = series_precision()
        else:
            cprec = order
        sig_on()
        try:
            x = self._gobj.expand(0).series(symbol0._gobj, cprec, options)
            nex = SymbolicSeries.__new__(SymbolicSeries)
            nex._parent = self._parent
            nex._gobj = GEx(x)
        finally:
            sig_off()
        return nex

    def residue(self, symbol):
        """
        Calculate the residue of ``self`` with respect to ``symbol``.

        INPUT:

        - ``symbol`` -- a symbolic variable or symbolic equality such
          as ``x == 5``. If an equality is given, the expansion is
          around the value on the right hand side of the equality,
          otherwise at ``0``.

        OUTPUT: the residue of ``self``

        Say, ``symbol`` is ``x == a``, then this function calculates
        the residue of ``self`` at `x=a`, i.e., the coefficient of
        `1/(x-a)` of the series expansion of ``self`` around `a`.

        EXAMPLES::

            sage: (1/x).residue(x == 0)
            1
            sage: (1/x).residue(x == oo)
            -1
            sage: (1/x^2).residue(x == 0)
            0
            sage: (1/sin(x)).residue(x == 0)
            1
            sage: var('q, n, z')
            (q, n, z)
            sage: (-z^(-n-1)/(1-z/q)^2).residue(z == q).simplify_full()
            (n + 1)/q^n
            sage: var('s')
            s
            sage: zeta(s).residue(s == 1)
            1

        We can also compute the residue at more general places,
        given that the pole is recognized::

            sage: k = var('k', domain='integer')
            sage: (gamma(1+x)/(1 - exp(-x))).residue(x==2*I*pi*k)
            gamma(2*I*pi*k + 1)
            sage: csc(x).residue(x==2*pi*k)
            1

        TESTS::

            sage: (exp(x)/sin(x)^4).residue(x == 0)
            5/6

        Check that :issue:`18372` is resolved::

            sage: (1/(x^2 - x - 1)).residue(x == 1/2*sqrt(5) + 1/2)
            1/5*sqrt(5)

        Check that :issue:`20084` is fixed::

            sage: (1/(1 - 2^-x)).residue(x == 2*pi*I/log(2))
            1/log(2)
        """
        if symbol.is_relational():
            x = symbol.lhs()
            a = symbol.rhs()
        else:
            x = symbol
            a = 0
        if a == infinity:
            return (-self.subs({x: 1/x}) / x**2).residue(x == 0)
        return self.subs({x: x + a}).series(x == 0, 0).coefficient(x, -1)

    def taylor(self, *args):
        r"""
        Expand this symbolic expression in a truncated Taylor or
        Laurent series in the variable `v` around the point `a`,
        containing terms through `(x - a)^n`. Functions in more
        variables is also supported.

        INPUT:

        - ``*args`` -- the following notation is supported

           - ``x``, ``a``, ``n`` -- variable, point, degree

           - ``(x, a)``, ``(y, b)``, ``n`` -- variables with points, degree of polynomial

        EXAMPLES::

            sage: var('a, x, z')
            (a, x, z)
            sage: taylor(a*log(z), z, 2, 3)
            1/24*a*(z - 2)^3 - 1/8*a*(z - 2)^2 + 1/2*a*(z - 2) + a*log(2)

        ::

            sage: taylor(sqrt(sin(x) + a*x + 1), x, 0, 3)
            1/48*(3*a^3 + 9*a^2 + 9*a - 1)*x^3 - 1/8*(a^2 + 2*a + 1)*x^2 + 1/2*(a + 1)*x + 1

        ::

            sage: taylor(sqrt(x + 1), x, 0, 5)
            7/256*x^5 - 5/128*x^4 + 1/16*x^3 - 1/8*x^2 + 1/2*x + 1

        ::

            sage: taylor(1/log(x + 1), x, 0, 3)
            -19/720*x^3 + 1/24*x^2 - 1/12*x + 1/x + 1/2

        ::

            sage: taylor(cos(x) - sec(x), x, 0, 5)
            -1/6*x^4 - x^2

        ::

            sage: taylor((cos(x) - sec(x))^3, x, 0, 9)
            -1/2*x^8 - x^6

        ::

            sage: taylor(1/(cos(x) - sec(x))^3, x, 0, 5)
            -15377/7983360*x^4 - 6767/604800*x^2 + 11/120/x^2 + 1/2/x^4 - 1/x^6 - 347/15120

        TESTS:

        Check that issue :issue:`7472` is fixed (Taylor polynomial in
        more variables)::

            sage: x,y = var('x y'); taylor(x*y^3,(x,1),(y,1),4)
            (x - 1)*(y - 1)^3 + 3*(x - 1)*(y - 1)^2 + (y - 1)^3 + 3*(x - 1)*(y - 1) + 3*(y - 1)^2 + x + 3*y - 3
            sage: expand(_)
            x*y^3
        """
        from sage.rings.integer import Integer
        from sage.symbolic.ring import SR
        A = args
        try:
            if isinstance(A[0],tuple):
                B=[]
                B.append([SR(A[i][0]) for i in range(len(A)-1)])
                B.append([A[i][1] for i in range(len(A)-1)])
            else:
                B=[A[0],SR(A[1])]
            B.append(Integer(A[len(A)-1]))
        except Exception:
            raise NotImplementedError("Wrong arguments passed to taylor. See taylor? for more details.")
        l = self._maxima_().taylor(B)
        return self.parent()(l)

    def truncate(self):
        """
        Given a power series or expression, return the corresponding
        expression without the big oh.

        INPUT:

        - ``self`` -- a series as output by the :meth:`series` command

        OUTPUT: a symbolic expression

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
        return self

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
            x^2 + 2*x*sin(x^2 + 2*x*y + y^2) + sin(x^2 + 2*x*y + y^2)^2

        Observe that :meth:`expand` also expands function arguments::

            sage: f(x) = function('f')(x)
            sage: fx = f(x*(x+1)); fx
            f((x + 1)*x)
            sage: fx.expand()
            f(x^2 + x)

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

            sage: a,b,c = var('a,b,c')
            sage: x,y = var('x,y', domain='real')
            sage: p,q = var('p,q', domain='positive')
            sage: (c/2*(5*(3*a*b*x*y*p*q)^2)^(7/2*c)).expand()
            1/2*45^(7/2*c)*(a^2*b^2*x^2*y^2)^(7/2*c)*c*p^(7*c)*q^(7*c)
            sage: ((-(-a*x*p)^3*(b*y*p)^3)^(c/2)).expand()
            (a^3*b^3*x^3*y^3)^(1/2*c)*p^(3*c)
            sage: x,y,p,q = var('x,y,p,q', domain='complex')

        Check that :issue:`18568` is fixed::

            sage: ((x+sqrt(2)*x)^2).expand()
            2*sqrt(2)*x^2 + 3*x^2

        Check that :issue:`21360` is fixed::

            sage: ((x^(x/2) + 1)^2).expand()
            2*x^(1/2*x) + x^x + 1
            sage: ((x^(1/2*x))^2).expand()
            x^x
            sage: ((x^(2*x))^2).expand()
            x^(4*x)

        Check that exactness is preserved::

            sage: ((x+1.001)^2).expand()
            x^2 + 2.00200000000000*x + 1.00200100000000
            sage: ((x+1.001)^3).expand()
            x^3 + 3.00300000000000*x^2 + 3.00600300000000*x + 1.00300300100000

        Check that :issue:`21302` is fixed::

            sage: ((x+1)^-2).expand()
            1/(x^2 + 2*x + 1)
            sage: (((x-1)/(x+1))^2).expand()
            x^2/(x^2 + 2*x + 1) - 2*x/(x^2 + 2*x + 1) + 1/(x^2 + 2*x + 1)

        Check that :issue:`30688` is fixed::

            sage: assume(x < 0)
            sage: sqrt(-x).expand()
            sqrt(-x)
            sage: ((-x)^(3/4)).expand()
            (-x)^(3/4)
            sage: forget()

        Check that :issue:`31077` and :issue:`31585` are fixed (also see :issue:`31679`)::

            sage: a,b,c,d = var("a b c d")
            sage: ((a + b + c)^30 * (3*b + d - 5/d)^3).expand().subs(a=0,b=2,c=-1)
            d^3 + 18*d^2 + 93*d - 465/d + 450/d^2 - 125/d^3 + 36

        Check that :issue:`31411` is fixed::

            sage: q, j = var("q, j")
            sage: A = q^(2/3) + q^(2/5)
            sage: B = product(1 - q^j, j, 1, 31) * q^(1/24)
            sage: bool((A * B).expand() == (A * B.expand()).expand())
            True
        """
        if side is not None:
            if not is_a_relational(self._gobj):
                raise ValueError("expansion on sides only makes sense for relations")
            if side == 'left':
                return self.operator()(self.lhs().expand(), self.rhs())
            elif side == 'right':
                return self.operator()(self.lhs(), self.rhs().expand())
            else:
                raise ValueError("side must be 'left', 'right', or None")

        cdef GEx x
        sig_on()
        try:
            x = self._gobj.expand(0)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    expand_rational = rational_expand = expand

    def expand_trig(self, full=False, half_angles=False, plus=True, times=True):
        r"""
        Expand trigonometric and hyperbolic functions of sums of angles
        and of multiple angles occurring in ``self``.

        For best results, ``self`` should already be expanded.

        INPUT:

        - ``full`` -- boolean (default: ``False``); to enhance user control
          of simplification, this function expands only one level at a time
          by default, expanding sums of angles or multiple angles. To obtain
          full expansion into sines and cosines immediately, set the optional
          parameter full to ``True``.

        - ``half_angles`` -- boolean (default: ``False``); if ``True``, causes
          half-angles to be simplified away

        - ``plus`` -- boolean (default: ``True``); controls the sum rule.
          Expansion of sums (e.g. `\sin(x + y)`) will take place only if
          ``plus`` is ``True``.

        - ``times`` -- boolean (default: ``True``); controls the product
          rule, expansion of products (e.g. `\sin(2 x)`) will take place only
          if ``times`` is ``True``.

        OUTPUT: a symbolic expression

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
            sin((3*(cos(cos(x)^2)*cos(sin(x)^2)
                  + sin(cos(x)^2)*sin(sin(x)^2))^2*(cos(sin(x)^2)*sin(cos(x)^2)
                                                     - cos(cos(x)^2)*sin(sin(x)^2))
                  - (cos(sin(x)^2)*sin(cos(x)^2) - cos(cos(x)^2)*sin(sin(x)^2))^3)*x)
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

        If the expression contains terms which are factored, we expand first::

            sage: (x, k1, k2) = var('x, k1, k2')
            sage: cos((k1-k2)*x).expand().expand_trig()
            cos(k1*x)*cos(k2*x) + sin(k1*x)*sin(k2*x)

        ALIAS:

        :meth:`trig_expand` and :meth:`expand_trig` are the same
        """
        from sage.calculus.calculus import maxima_options
        M = self._maxima_()
        P = M.parent()
        opt = maxima_options(trigexpand=full, halfangles=half_angles,
                             trigexpandplus=plus, trigexpandtimes=times)
        cmd = 'trigexpand(%s), %s' % (M.name(), opt)
        ans = P(cmd)
        return self.parent()(ans)

    trig_expand = expand_trig

    def reduce_trig(self, var=None):
        r"""
        Combine products and powers of trigonometric and hyperbolic
        sin's and cos's of x into those of multiples of x. It also
        tries to eliminate these functions when they occur in
        denominators.

        INPUT:

        - ``self`` -- a symbolic expression

        - ``var`` -- (default: ``None``) the variable which is used for
          these transformations. If not specified, all variables are
          used.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: y = var('y')
            sage: f = sin(x)*cos(x)^3+sin(y)^2
            sage: f.reduce_trig()
            -1/2*cos(2*y) + 1/8*sin(4*x) + 1/4*sin(2*x) + 1/2

        To reduce only the expressions involving x we use optional parameter::

            sage: f.reduce_trig(x)
            sin(y)^2 + 1/8*sin(4*x) + 1/4*sin(2*x)

        ALIAS: :meth:`trig_reduce` and :meth:`reduce_trig` are the same
        """
        M = self._maxima_()
        P = M.parent()
        if var is None:
            cmd = 'trigreduce(%s)' % (M.name())
        else:
            cmd = 'trigreduce(%s,%s)' % (M.name(),'_SAGE_VAR_' + str(var))
        ans = P(cmd)
        return self.parent()(ans)

    trig_reduce = reduce_trig

    ############################################################################
    # Pattern Matching
    ############################################################################
    def match(self, pattern):
        """
        Check if ``self`` matches the given pattern.

        INPUT:

        - ``pattern`` -- a symbolic expression, possibly containing wildcards
          to match for

        OUTPUT:

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
            sage: print(((x+y)^a).match((x+y)^b))
            None
            sage: t = ((x+y)^a).match(w0^w1)
            sage: t[w0], t[w1]
            (x + y, a)
            sage: print(((x+y)^a).match(w0^w0))
            None
            sage: ((x+y)^(x+y)).match(w0^w0)
            {$0: x + y}
            sage: t = ((a+b)*(a+c)).match((a+w0)*(a+w1))
            sage: set([t[w0], t[w1]]) == set([b, c])
            True
            sage: ((a+b)*(a+c)).match((w0+b)*(w0+c))
            {$0: a}
            sage: t = ((a+b)*(a+c)).match((w0+w1)*(w0+w2))
            sage: t[w0]
            a
            sage: set([t[w1], t[w2]]) == set([b, c])
            True
            sage: t = ((a+b)*(a+c)).match((w0+w1)*(w1+w2))
            sage: t[w1]
            a
            sage: set([t[w0], t[w2]]) == set([b, c])
            True
            sage: t = (a*(x+y)+a*z+b).match(a*w0+w1)
            sage: s = set([t[w0], t[w1]])
            sage: s == set([x+y, a*z+b]) or s == set([z, a*(x+y)+b])
            True
            sage: print((a+b+c+d+f+g).match(c))
            None
            sage: (a+b+c+d+f+g).has(c)
            True
            sage: (a+b+c+d+f+g).match(c+w0)
            {$0: a + b + d + f + g}
            sage: (a+b+c+d+f+g).match(c+g+w0)
            {$0: a + b + d + f}
            sage: (a+b).match(a+b+w0) # known bug
            {$0: 0}
            sage: print((a*b^2).match(a^w0*b^w1))
            None
            sage: (a*b^2).match(a*b^w1)
            {$1: 2}
            sage: (x*x.arctan2(x^2)).match(w0*w0.arctan2(w0^2))
            {$0: x}

        Beware that behind-the-scenes simplification can lead to
        surprising results in matching::

            sage: print((x+x).match(w0+w1))
            None
            sage: t = x+x; t
            2*x
            sage: t.operator()
            <function mul_vararg ...>

        Since asking to match w0+w1 looks for an addition operator,
        there is no match.
        """
        cdef Expression p = self.coerce_in(pattern)
        cdef GExList mlst
        cdef bint res
        sig_on()
        try:
            res = self._gobj.match(p._gobj, mlst)
        finally:
            sig_off()
        if not res:
            return None

        cdef dict rdict = {}
        cdef GExListIter itr = mlst.begin()
        cdef GExListIter lstend = mlst.end()
        while itr != lstend:
            key = new_Expression_from_GEx(self._parent, itr.obj().lhs())
            val = new_Expression_from_GEx(self._parent, itr.obj().rhs())
            rdict[key] = val
            itr.inc()
        return rdict

    def find(self, pattern):
        """
        Find all occurrences of the given pattern in this expression.

        Note that once a subexpression matches the pattern, the search does
        not extend to subexpressions of it.

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
        sig_on()
        try:
            self._gobj.find(p._gobj, found)
        finally:
            sig_off()
        res = []
        cdef GExListIter itr = found.begin()
        while itr != found.end():
            res.append(new_Expression_from_GEx(self._parent, itr.obj()))
            itr.inc()
        res = print_sorted(res)
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
            sage: (4*x^2 - x + 3).coefficient(x,1)
            -1
            sage: (4*x^2 + x + 3).coefficient(x,1)
            1
        """
        cdef Expression p = self.coerce_in(pattern)
        return self._gobj.has(p._gobj)

    def subs(self, *args, **kwds):
        """
        Substitute the given subexpressions in this expression.

        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,f,g')
            (x, y, z, a, b, c, d, f, g)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: t = a^2 + b^2 + (x+y)^3

        Substitute with keyword arguments (works only with symbols)::

            sage: t.subs(a=c)
            (x + y)^3 + b^2 + c^2
            sage: t.subs(b=19, x=z)
            (y + z)^3 + a^2 + 361

        Substitute with a dictionary argument::

            sage: t.subs({a^2: c})
            (x + y)^3 + b^2 + c

            sage: t.subs({w0^2: w0^3})
            a^3 + b^3 + (x + y)^3

        Substitute with one or more relational expressions::

            sage: t.subs(w0^2 == w0^3)
            a^3 + b^3 + (x + y)^3

            sage: t.subs(w0 == w0^2)
            a^8 + b^8 + (x^2 + y^2)^6

            sage: t.subs(a == b, b == c)
            (x + y)^3 + b^2 + c^2

        Any number of arguments is accepted::

            sage: t.subs(a=b, b=c)
            (x + y)^3 + b^2 + c^2

            sage: t.subs({a:b}, b=c)
            (x + y)^3 + b^2 + c^2

            sage: t.subs([x == 3, y == 2], a == 2, {b:3})
            138

        It can even accept lists of lists::

            sage: eqn1 = (a*x + b*y == 0)
            sage: eqn2 = (1 + y == 0)
            sage: soln = solve([eqn1, eqn2], [x, y])
            sage: soln
            [[x == b/a, y == -1]]
            sage: f = x + y
            sage: f.subs(soln)
            b/a - 1

        Duplicate assignments will throw an error::

            sage: t.subs({a:b}, a=c)
            Traceback (most recent call last):
            ...
            ValueError: duplicate substitution for a, got values b and c

            sage: t.subs([x == 1], a = 1, b = 2, x = 2)
            Traceback (most recent call last):
            ...
            ValueError: duplicate substitution for x, got values 1 and 2

        All substitutions are performed at the same time::

             sage: t.subs({a:b, b:c})
             (x + y)^3 + b^2 + c^2

        Substitutions are done term by term, in other words Sage is not
        able to identify partial sums in a substitution (see :issue:`18396`)::

            sage: f = x + x^2 + x^4
            sage: f.subs(x = y)
            y^4 + y^2 + y
            sage: f.subs(x^2 == y)             # one term is fine
            x^4 + x + y
            sage: f.subs(x + x^2 == y)         # partial sum does not work
            x^4 + x^2 + x
            sage: f.subs(x + x^2 + x^4 == y)   # whole sum is fine
            y

        Note that it is the very same behavior as in Maxima::

            sage: E = 'x^4 + x^2 + x'
            sage: subs = [('x','y'), ('x^2','y'), ('x^2+x','y'), ('x^4+x^2+x','y')]

            sage: cmd = '{}, {}={}'
            sage: for s1,s2 in subs:
            ....:     maxima.eval(cmd.format(E, s1, s2))
            'y^4+y^2+y'
            'y+x^4+x'
            'x^4+x^2+x'
            'y'

        Or as in Maple::

            sage: cmd = 'subs({}={}, {})'              # optional - maple
            sage: for s1,s2 in subs:                   # optional - maple
            ....:     maple.eval(cmd.format(s1,s2, E))
            'y^4+y^2+y'
            'x^4+x+y'
            'x^4+x^2+x'
            'y'

        But Mathematica does something different on the third example::

            sage: cmd = '{} /. {} -> {}'                    # optional - mathematica
            sage: for s1,s2 in subs:                        # optional - mathematica
            ....:     mathematica.eval(cmd.format(E,s1,s2))
                 2    4
            y + y  + y
                 4
            x + x  + y
             4
            x  + y
            y

        The same, with formatting more suitable for cut and paste::

            sage: for s1,s2 in subs:                        # optional - mathematica
            ....:     mathematica(cmd.format(E,s1,s2))
            y + y^2 + y^4
            x + x^4 + y
            x^4 + y
            y

        .. WARNING::

            Unexpected results may occur if the left-hand side of some substitution
            is not just a single variable (or is a "wildcard" variable). For example,
            the result of ``cos(cos(cos(x))).subs({cos(x) : x})`` is ``x``, because
            the substitution is applied repeatedly. Such repeated substitutions (and
            pattern-matching code that may be somewhat unpredictable) are disabled
            only in the basic case where the left-hand side of every substitution is
            a variable. In particular, although the result of
            ``(x^2).subs({x : sqrt(x)})`` is ``x``, the result of
            ``(x^2).subs({x : sqrt(x), y^2 : y})`` is ``sqrt(x)``, because repeated
            substitution is enabled by the presence of the expression ``y^2`` in the
            left-hand side of one of the substitutions, even though that particular
            substitution does not get applied.

        TESTS:

        No arguments return the same expression::

            sage: t = a^2 + b^2 + (x+y)^3
            sage: t.subs()
            (x + y)^3 + a^2 + b^2

        Similarly for a empty dictionary, empty tuples and empty lists::

            sage: t.subs({}, (), [], ())
            (x + y)^3 + a^2 + b^2

        Invalid argument returns error::

            sage: t.subs(5)
            Traceback (most recent call last):
            ...
            TypeError: not able to determine a substitution from 5

        Substitutions with infinity::

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

        Verify that this operation does not modify the passed
        dictionary (:issue:`6622`)::

            sage: var('v t')
            (v, t)
            sage: f = v*t
            sage: D = {v: 2}
            sage: f(D, t=3)
            6
            sage: D
            {v: 2}

        Check if :issue:`9891` is fixed::

            sage: exp(x).subs(x=log(x))
            x

        Check if :issue:`13587` is fixed::

            sage: t = tan(x)^2 - tan(x)
            sage: t.subs(x=pi/2)
            Infinity
            sage: u = gamma(x) - gamma(x-1)
            sage: u.subs(x=-1)
            Infinity

        More checks for ``subs``::

            sage: var('x,y,z'); f = x^3 + y^2 + z
            (x, y, z)
            sage: f.subs(x^3 == y^2, z == 1)
            2*y^2 + 1
            sage: f.subs({x^3:y^2, z:1})
            2*y^2 + 1
            sage: f = x^2 + x^4
            sage: f.subs(x^2 == x)
            x^4 + x
            sage: f = cos(x^2) + sin(x^2)
            sage: f.subs(x^2 == x)
            cos(x) + sin(x)
            sage: f(x,y,t) = cos(x) + sin(y) + x^2 + y^2 + t
            sage: f.subs(y^2 == t)
            (x, y, t) |--> x^2 + 2*t + cos(x) + sin(y)
            sage: f.subs(x^2 + y^2 == t)
            (x, y, t) |--> x^2 + y^2 + t + cos(x) + sin(y)

        Check that inverses in sums are recognized::

            sage: (1 + 1/x).subs({x: 1/x})
            x + 1
            sage: (x + 1/x^2).subs({x: 1/x})
            x^2 + 1/x
            sage: (sqrt(x) + 1/sqrt(x)).subs({x: 1/x})
            sqrt(x) + 1/sqrt(x)

        Check that :issue:`30378` is fixed::

            sage: (x^2).subs({x: sqrt(x)})
            x
            sage: f(x) = x^2
            sage: f(sqrt(x))
            x
            sage: a = var("a")
            sage: f = function("f")
            sage: integrate(f(x), x, 0, a).subs(a=cos(a))
            integrate(f(x), x, 0, cos(a))

        Check that :issue:`31554` is fixed::

            sage: a,b,c,d,x,y = var("a b c d x y")
            sage: with hold:
            ....:     print((2*x^0*a + b*y^1).subs({x:c, y:c*d}))
            b*c*d + 2*a

        Check that :issue:`31530` is fixed::

            sage: a, b = var("a b")
            sage: (a + b*x).series(x, 2).subs(a=a, b=b)
            (a) + (b)*x + Order(x^2)

        Check that :issue:`31585` is fixed::

            sage: m = -2^31
            sage: (-x).subs(x=m)
            2147483648
            sage: abs(x).subs(x=m)
            2147483648
            sage: (2*x).subs(x=m)
            -4294967296
            sage: (m*x + 1)*x
            -(2147483648*x - 1)*x
            sage: m = -2^63
            sage: (-x).subs(x=m)
            9223372036854775808
            sage: abs(x).subs(x=m)
            9223372036854775808
            sage: (2*x).subs(x=m)
            -18446744073709551616
            sage: (m*x + 1)*x
            -(9223372036854775808*x - 1)*x
        """
        cdef dict sdict = {}
        cdef GEx res

        if args and args[0] is None:
            # this is needed because sometimes this function get called as
            # expr.substitute(None, **kwds). This is because its signature used
            # to be (in_dict=None, **kwds) instead of (*args, **kwds)
            # (see issue #12834)
            args = args[1:]

        for a in args:
            _dict_update_check_duplicate(sdict, _subs_make_dict(a))

        if kwds:
            # Ensure that the keys are symbolic variables.
            varkwds = {self._parent.var(k): v for k,v in kwds.items()}
            # Check for duplicate
            _dict_update_check_duplicate(sdict, varkwds)

        cdef GExMap smap
        for k, v in sdict.items():
            smap.insert(make_pair((<Expression>self.coerce_in(k))._gobj,
                                  (<Expression>self.coerce_in(v))._gobj))
        res = self._gobj.subs_map(smap, 0)
        return new_Expression_from_GEx(self._parent, res)

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
        cdef GEx res
        sig_on()
        try:
            res = self._gobj.subs(p._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, res)

    def substitute_function(self, *args, **kwds):
        """
        Substitute the given functions by their replacements in this expression.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: foo = function('foo'); bar = function('bar')
            sage: f = foo(x) + 1/foo(pi*y)

        Substitute with a dictionary::

            sage: f.substitute_function({foo: bar})
            1/bar(pi*y) + bar(x)
            sage: f.substitute_function({foo(x): bar(x)})
            1/bar(pi*y) + bar(x)

        If the function expression to be substituted includes its arguments, the right hand side can be an arbitrary symbolic expression::

            sage: f.substitute_function({foo(x): x^2})
            x^2 + 1/(pi^2*y^2)

        Substitute with keyword arguments (works only if no function arguments are given)::

            sage: f.substitute_function(foo=bar)
            1/bar(pi*y) + bar(x)

        Substitute with a relational expression::

            sage: f.substitute_function(foo(x)==bar(x))
            1/bar(pi*y) + bar(x)
            sage: f.substitute_function(foo(x)==bar(x+1))
            1/bar(pi*y + 1) + bar(x + 1)

        All substitutions are performed at the same time::

            sage: g = foo(x) + 1/bar(pi*y)
            sage: g.substitute_function({foo: bar, bar: foo})
            1/foo(pi*y) + bar(x)

        Any number of arguments is accepted::

            sage: g.substitute_function({foo: bar}, bar(x) == x^2)
            1/(pi^2*y^2) + bar(x)

        As well as lists of substitutions::

            sage: g.substitute_function([foo(x) == 1, bar(x) == x])
            1/(pi*y) + 1

        Alternative syntax::

            sage: g.substitute_function(foo, bar)
            1/bar(pi*y) + bar(x)

        Duplicate assignments will throw an error::

            sage: g.substitute_function({foo:bar}, foo(x) == x^2)
            Traceback (most recent call last):
            ...
            ValueError: duplicate substitution for foo, got values bar and x |--> x^2

            sage: g.substitute_function([foo(x) == x^2], foo = bar)
            Traceback (most recent call last):
            ...
            ValueError: duplicate substitution for foo, got values x |--> x^2 and bar

        TESTS:

        Make sure :issue:`17849` is fixed::

            sage: ex = sin(x) + atan2(0,0,hold=True)
            sage: ex.substitute_function(sin,cos)
            arctan2(0, 0) + cos(x)
            sage: ex = sin(x) + hypergeometric([1, 1], [2], -1)
            sage: ex.substitute_function(sin,cos)
            cos(x) + hypergeometric((1, 1), (2,), -1)
        """
        from sage.symbolic.expression_conversions import SubstituteFunction
        from sage.calculus.calculus import _find_func, _is_function

        sdict = {}

        if len(args) == 2 and _is_function(args[0]):
            sdict = {args[0]: args[1]}
        else:
            for a in args:
                _dict_update_check_duplicate(sdict, _subs_fun_make_dict(a))

        if kwds:
            # Ensure that the keys are functions.
            funkwds = {_find_func(k): v for k,v in kwds.items()}
            # Check for duplicate
            _dict_update_check_duplicate(sdict, funkwds)

        return SubstituteFunction(self, sdict)()

    def exponentialize(self):
        r"""
        Return this symbolic expression with all circular and hyperbolic
        functions replaced by their respective exponential
        expressions.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: sin(x).exponentialize()
            -1/2*I*e^(I*x) + 1/2*I*e^(-I*x)
            sage: sec(x).exponentialize()
            2/(e^(I*x) + e^(-I*x))
            sage: tan(x).exponentialize()
            (-I*e^(I*x) + I*e^(-I*x))/(e^(I*x) + e^(-I*x))
            sage: sinh(x).exponentialize()
            -1/2*e^(-x) + 1/2*e^x
            sage: sech(x).exponentialize()
            2/(e^(-x) + e^x)
            sage: tanh(x).exponentialize()
            -(e^(-x) - e^x)/(e^(-x) + e^x)

        TESTS:

        Check that ``u(x).exponentialize().demoivre(force=True)``
        is identity::

            sage: x = SR.var("x")
            sage: all([bool(u(x).exponentialize().demoivre(force=True) == u(x))
            ....:      for u in (sin, cos, tan, csc, sec, cot,
            ....:                sinh, cosh, tanh, csch, sech, coth)])
            True

        Check that differentiation and exponentialization commute::

            sage: x = SR.var("x")
            sage: all([bool(u(x).diff(x).exponentialize() ==
            ....:           u(x).exponentialize().diff(x))
            ....:      for u in (sin, cos, tan, csc, sec, cot,
            ....:                sinh, cosh, tanh, csch, sech, coth)])
            True
        """
        from sage.symbolic.expression_conversions import Exponentialize
        return Exponentialize(self)()

    def demoivre(self, force=False):
        r"""
        Return this symbolic expression with complex exponentials
        (optionally all exponentials) replaced by (at least partially)
        trigonometric/hyperbolic expressions.

        EXAMPLES::

            sage: x, a, b = SR.var("x, a, b")
            sage: exp(a + I*b).demoivre()
            (cos(b) + I*sin(b))*e^a
            sage: exp(I*x).demoivre()
            cos(x) + I*sin(x)
            sage: exp(x).demoivre()
            e^x
            sage: exp(x).demoivre(force=True)
            cosh(x) + sinh(x)

        TESTS:

        Check that de Moivre transformation correctly commutes
        with differentiation::

            sage: x = SR.var("x")
            sage: f = function("f")
            sage: bool(f(exp(I*x)).diff(x).demoivre() ==
            ....:      f(exp(I*x)).demoivre().diff(x))
            True
        """
        from sage.symbolic.expression_conversions import DeMoivre
        return DeMoivre(self, force)()

    def half_angle(self):
        """
        Replace all occurrences of trigonometric (or hyperbolic)
        functions by rational fractions of the (hyperbolic) tangent
        of half the original argument.

        This can help highlight the algebraic structure of an expression,
        which can be useful e.g. for integration.

        This method has no direct relation with the ``half_angles``
        argument of the :meth:`trig_expand` method.

        EXAMPLES::

            sage: x, t = var("x, t")
            sage: cos(x).half_angle().subs(tan(x/2) == t)
            -(t^2 - 1)/(t^2 + 1)

        Note that this structure highlighting works better after expansion::

            sage: x, t = var("x, t")
            sage: a = (cos(3*x)/(4-cos(x)))
            sage: b = a.trig_expand()
            sage: a.half_angle().subs(tan(x/2) == t).simplify_full()
            (2*(t^2 + 1)*cos(3/2*x)^2 - t^2 - 1)/(5*t^2 + 3)
            sage: b.half_angle().subs(tan(x/2) == t).simplify_full()
            -(t^6 - 15*t^4 + 15*t^2 - 1)/(5*t^6 + 13*t^4 + 11*t^2 + 3)

        TESTS::

            sage: all((u(x) == u(x).half_angle()).subs(x == 2*x).trig_simplify()
            ....:     for u in (sin, cos, tan, csc, sec, cot,
            ....:               sinh, cosh, tanh, csch, sech, coth))
            True
        """
        from sage.symbolic.expression_conversions import HalfAngle
        return HalfAngle(self)()

    def substitution_delayed(self, pattern, replacement):
        """
        Replace all occurrences of pattern by the result of replacement.

        In contrast to :meth:`subs`, the pattern may contains wildcards
        and the replacement can depend on the particular term matched by the
        pattern.

        INPUT:

        - ``pattern`` -- an :class:`Expression`, usually containing wildcards

        - ``replacement`` -- a function; its argument is a dictionary
          mapping the wildcard occurring in ``pattern`` to the actual
          values.  If it returns ``None``, this occurrence of ``pattern`` is
          not replaced. Otherwise, it is replaced by the output of
          ``replacement``.

        OUTPUT: an :class:`Expression`

        EXAMPLES::

            sage: var('x y')
            (x, y)
            sage: w0 = SR.wild(0)
            sage: sqrt(1 + 2*x + x^2).substitution_delayed(
            ....:     sqrt(w0), lambda d: sqrt(factor(d[w0]))
            ....: )
            sqrt((x + 1)^2)
            sage: def r(d):
            ....:    if x not in d[w0].variables():
            ....:        return cos(d[w0])
            sage: (sin(x^2 + x) + sin(y^2 + y)).substitution_delayed(sin(w0), r)
            cos(y^2 + y) + sin(x^2 + x)

        .. SEEALSO::

            :meth:`match`
        """
        result = self
        for matched in self.find(pattern):
            r = replacement(matched.match(pattern))
            if r is not None:
                result = result.subs({matched: r})
        return result

    def __call__(self, *args, **kwds):
        """
        Call the :meth:`subs` on this expression.

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
        while itr != sym_set.end():
            res.append(new_Expression_from_GEx(SR, itr.obj()))
            itr.inc()
        res = print_sorted(res)[::-1]
        return tuple(res)

    def free_variables(self):
        """
        Return sorted tuple of unbound variables that occur in this
        expression.

        EXAMPLES::

            sage: (x,y,z) = var('x,y,z')
            sage: (x+y).free_variables()
            (x, y)
            sage: (2*x).free_variables()
            (x,)
            sage: (x^y).free_variables()
            (x, y)
            sage: sin(x+y^z).free_variables()
            (x, y, z)
            sage: _ = function('f')
            sage: e = limit( f(x,y), x=0 ); e
            limit(f(x, y), x, 0)
            sage: e.free_variables()
            (y,)
        """
        from sage.symbolic.ring import SR
        cdef GSymbolSet sym_set
        sym_set = self._gobj.free_symbols()
        res = []
        cdef GSymbolSetIter itr = sym_set.begin()
        while itr != sym_set.end():
            res.append(new_Expression_from_GEx(SR, GEx(itr.obj())))
            itr.inc()
        res = print_sorted(res)[::-1]
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
        """
        return len(self.arguments())

    def number_of_operands(self):
        """
        Return the number of operands of this expression.

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

    def _unpack_operands(self):
        """
        Unpack the operands of this expression converting each to a Python
        object if possible.

        This corresponds to the conversion performed when arguments of a
        function are unpacked as they are being passed to custom methods of
        a symbolic function.

        EXAMPLES::

            sage: t = SR._force_pyobject((1, 2, x, x+1, x+2))
            sage: t._unpack_operands()
            (1, 2, x, x + 1, x + 2)
            sage: type(t._unpack_operands())
            <... 'tuple'>
            sage: list(map(type, t._unpack_operands()))
            [<class 'sage.rings.integer.Integer'>, <class 'sage.rings.integer.Integer'>, <class 'sage.symbolic.expression.Expression'>, <class 'sage.symbolic.expression.Expression'>, <class 'sage.symbolic.expression.Expression'>]
            sage: u = SR._force_pyobject((t, x^2))
            sage: u._unpack_operands()
            ((1, 2, x, x + 1, x + 2), x^2)
            sage: type(u._unpack_operands()[0])
            <... 'tuple'>
        """
        return unpack_operands(self)

    def operands(self):
        """
        Return a list containing the operands of this expression.

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
        return [new_Expression_from_GEx(SR, self._gobj.op(i))
                for i in range(self._gobj.nops())]

    def operator(self):
        """
        Return the topmost operator in this expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: (x+y).operator()
            <function add_vararg ...>
            sage: (x^y).operator()
            <built-in function pow>
            sage: (x^y * z).operator()
            <function mul_vararg ...>
            sage: (x < y).operator()
            <built-in function lt>

            sage: abs(x).operator()
            abs
            sage: r = gamma(x).operator(); type(r)
            <class 'sage.functions.gamma.Function_gamma'>

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
            diff(f(x), x)
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
            sage: SR._force_pyobject( (x, x + 1, x + 2) ).operator()
            <... 'tuple'>
            sage: exp(x).series(x,3).operator()
            <function add_vararg ...>
        """
        cdef operators o
        cdef unsigned serial
        from sage.symbolic.operators import add_vararg, mul_vararg
        if is_a_add(self._gobj):
            return add_vararg
        elif is_a_mul(self._gobj):
            return mul_vararg
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
                raise RuntimeError("operator type not known, please report this as a bug")
        elif is_a_function(self._gobj):
            # get function id
            serial = ex_to_function(self._gobj).get_serial()

            # if operator is a special function defined by us
            # find the python equivalent and return it
            res = get_sfunction_from_serial(serial)
            if res is None:
                raise RuntimeError("cannot find SymbolicFunction in table")

            if is_a_fderivative(self._gobj):
                parameter_set = paramset_from_Expression(self)
                from sage.symbolic.operators import FDerivativeOperator
                res = FDerivativeOperator(res, parameter_set)

            return res
        elif is_exactly_a_exprseq(self._gobj):
            return tuple
        elif is_a_series(self._gobj):
            from sage.symbolic.operators import add_vararg
            return add_vararg

        # self._gobj is either a symbol, constant or numeric
        return None

    def __index__(self):
        """
        EXAMPLES::

            sage: a = list(range(10))
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

        Note that symbols, constants and numeric objects do not have operands,
        so the iterator function raises an error in these cases::

            sage: x.iterator()
            Traceback (most recent call last):
            ...
            ValueError: expressions containing only a numeric coefficient,
            constant or symbol have no operands
            sage: pi.iterator()
            Traceback (most recent call last):
            ...
            ValueError: expressions containing only a numeric coefficient,
            constant or symbol have no operands
            sage: SR(5).iterator()
            Traceback (most recent call last):
            ...
            ValueError: expressions containing only a numeric coefficient,
            constant or symbol have no operands
        """
        if (is_a_symbol(self._gobj) or is_a_constant(self._gobj) or
                is_a_numeric(self._gobj)):
            raise ValueError("expressions containing only a numeric coefficient, constant or symbol have no operands")
        return new_ExpIter_from_Expression(self)

    @property
    def op(self):
        """
        Provide access to the operands of an expression through a property.

        EXAMPLES::

            sage: t = 1+x+x^2
            sage: t.op
            Operands of x^2 + x + 1
            sage: x.op
            Traceback (most recent call last):
            ...
            TypeError: expressions containing only a numeric coefficient,
            constant or symbol have no operands
            sage: t.op[0]
            x^2

        Indexing directly with ``t[1]`` causes problems with numpy types.

            sage: t[1]
            Traceback (most recent call last):
            ...
            TypeError: 'sage.symbolic.expression.Expression' object ...
        """
        if (is_a_symbol(self._gobj) or is_a_constant(self._gobj) or
                is_a_numeric(self._gobj)):
            raise TypeError("expressions containing only a numeric coefficient, constant or symbol have no operands")
        cdef OperandsWrapper res = OperandsWrapper.__new__(OperandsWrapper)
        res._expr = self
        return res

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        """
        Return a numerical approximation of ``self`` with ``prec`` bits
        (or decimal ``digits``) of precision.

        No guarantee is made about the accuracy of the result.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- which algorithm to use to compute this
          approximation

        If neither ``prec`` nor ``digits`` is given, the default
        precision is 53 bits (roughly 16 digits).

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
            Traceback (most recent call last):
            ...
            ValueError: can only convert signed infinity to RR

        Some expressions cannot be evaluated numerically::

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
        sufficient precision (:issue:`10164`)::

            sage: ks = 4*10**5, 10**6
            sage: all(len(str(e.n(digits=k)))-1 >= k for k in ks)
            True

        Symbolic sums with definite endpoints are expanded (:issue:`9424`)::

            sage: (k,n) = var('k,n')
            sage: f(n) = sum(abs(-k*k+n),k,1,n)
            sage: ex = f(n=8); ex
            sum(abs(-k^2 + 8), k, 1, 8)
            sage: ex.n()
            162.000000000000
            sage: (ex+1).n()
            163.000000000000

        Check if :issue:`24418` is fixed::

            sage: numerical_approx(2^(450232897/4888643760))
            1.06591892580915
        """
        if prec is None:
            prec = digits_to_bits(digits)

        from sage.symbolic.expression_conversions import ExpressionTreeWalker

        class DefiniteSumExpander(ExpressionTreeWalker):
            def composition(self, ex, operator):
                if hasattr(operator, 'name') and operator.name() == 'sum' and (
                        is_a_numeric((<Expression>ex.operands()[2])._gobj)
                        and is_a_numeric((<Expression>ex.operands()[3])._gobj)):
                    from sage.calculus.calculus import symbolic_sum
                    return symbolic_sum(*(ex.operands()))
                return super().composition(ex, operator)

        s = DefiniteSumExpander(self)
        cdef Expression x = self._parent(s())
        from sage.rings.real_mpfr import RealField
        R = RealField(prec)
        kwds = {'parent': R, 'algorithm': algorithm}
        try:
            x = x._convert(kwds)
        except TypeError:  # numerical approximation for real number failed
            pass           # try again with complex
            kwds['parent'] = R.complex_field()
            x = x._convert(kwds)

        # we have to consider constants as well, since infinity is a constant
        # in pynac
        if is_a_numeric(x._gobj):
            res = py_object_from_numeric(x._gobj)
        elif is_a_constant(x._gobj):
            res = x.pyobject()
        else:
            raise TypeError("cannot evaluate symbolic expression numerically")

        # Important -- the  we get might not be a valid output for numerical_approx in
        # the case when one gets infinity.
        if isinstance(res, AnInfinity):
            return res.n(prec=prec, digits=digits)
        return res

    def round(self):
        """
        Round this expression to the nearest integer.

        EXAMPLES::

            sage: u = sqrt(43203735824841025516773866131535024)
            sage: u.round()
            207855083711803945
            sage: t = sqrt(Integer('1'*1000)).round(); print(str(t)[-10:])
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
        from sage.rings.real_mpfi import RIF
        try:
            rif_self = RIF(self)
        except TypeError:
            raise ValueError("could not convert %s to a real number" % self)
        half = 1 / sage.rings.integer.Integer(2)
        if rif_self < 0 or (rif_self.contains_zero() and self < 0):
            result = ceil(self - half)
        else:
            result = floor(self + half)
        if not isinstance(result, sage.rings.integer.Integer):
            raise ValueError("could not convert %s to a real number" % self)
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
        for i in args:
            if not (isinstance(i, Expression) and i.is_symbol()):
                break
        else:
            from sage.symbolic.callable import CallableSymbolicExpressionRing
            R = CallableSymbolicExpressionRing(args, check=False)
            return R(self)
        raise TypeError(f"must construct a function with symbolic variables as arguments, got {args}.")

    ############################################################################
    # Basic arithmetic wrappers
    # which allow disabling automatic evaluation with the hold parameter
    ############################################################################
    def power(self, exp, hold=False):
        """
        Return the current expression to the power ``exp``.

        To prevent automatic evaluation use the ``hold`` argument.

        EXAMPLES::

            sage: (x^2).power(2)
            x^4
            sage: (x^2).power(2, hold=True)
            (x^2)^2

        To then evaluate again, we use :meth:`unhold`::

            sage: a = (x^2).power(2, hold=True); a.unhold()
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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = x.add(x, hold=True); a.unhold()
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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = x.mul(x, hold=True); a.unhold()
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
    def coefficient(self, s, n=1):
        """
        Return the coefficient of `s^n` in this symbolic expression.

        INPUT:

        - ``s`` -- expression

        - ``n`` -- expression, default 1

        OUTPUT: a symbolic expression. The coefficient of `s^n`

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
            z^z

        Any coefficient can be queried::

            sage: (x^2 + 3*x^pi).coefficient(x, pi)
            3
            sage: (2^x + 5*x^x).coefficient(x, x)
            5

        TESTS:

        Check if :issue:`9505` is fixed::

            sage: var('x,y,z')
            (x, y, z)
            sage: f = x*y*z^2
            sage: f.coefficient(x*y)
            z^2
            sage: f.coefficient(x*y, 2)
            Traceback (most recent call last):
            ...
            TypeError: n != 1 only allowed for s being a variable

        Check that :issue:`19996` is fixed::

            sage: (x^(1/2)).coefficient(x, QQ(1)/3)
            0
            sage: (x^(1/2)).coefficient(x, 1/3)
            0
        """
        cdef Expression ss = self.coerce_in(s)
        cdef Expression nn = self.coerce_in(n)
        cdef GEx r
        if n != 1 and not is_a_symbol(ss._gobj):
            raise TypeError("n != 1 only allowed for s being a variable")

        # the following is a temporary fix for GiNaC bug #9505
        if is_a_mul(ss._gobj):  # necessarily n=1 here
            res = self
            for i in range(ss._gobj.nops()):
                res = res.coefficient(new_Expression_from_GEx(self._parent, ss._gobj.op(i)))
            return res
        sig_on()
        try:
            r = self._gobj.coeff(ss._gobj, nn._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    def coefficients(self, x=None, sparse=True):
        r"""
        Return the coefficients of this symbolic expression as a polynomial in x.

        INPUT:

        - ``x`` -- (optional) variable

        OUTPUT: depending on the value of ``sparse``,

        - A list of pairs ``(expr, n)``, where ``expr`` is a symbolic
          expression and ``n`` is a power (``sparse=True``, default)

        - A list of expressions where the ``n``-th element is the coefficient of
          ``x^n`` when ``self`` is seen as polynomial in ``x`` (``sparse=False``).

        EXAMPLES::

            sage: var('x, y, a')
            (x, y, a)
            sage: p = x^3 - (x-3)*(x^2+x) + 1
            sage: p.coefficients()
            [[1, 0], [3, 1], [2, 2]]
            sage: p.coefficients(sparse=False)
            [1, 3, 2]
            sage: p = x - x^3 + 5/7*x^5
            sage: p.coefficients()
            [[1, 1], [-1, 3], [5/7, 5]]
            sage: p.coefficients(sparse=False)
            [0, 1, 0, -1, 0, 5/7]
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: p.coefficients(a)
            [[x^2 + x + 1, 0], [-2*sqrt(2)*x, 1], [2, 2]]
            sage: p.coefficients(a, sparse=False)
            [x^2 + x + 1, -2*sqrt(2)*x, 2]
            sage: p.coefficients(x)
            [[2*a^2 + 1, 0], [-2*sqrt(2)*a + 1, 1], [1, 2]]
            sage: p.coefficients(x, sparse=False)
            [2*a^2 + 1, -2*sqrt(2)*a + 1, 1]

        TESTS:

        The behaviour is undefined with noninteger or negative exponents::

            sage: p = (17/3*a)*x^(3/2) + x*y + 1/x + 2*x^x + 5*x^y
            sage: rset = set([(1, -1), (y, 1), (17/3*a, 3/2), (2, x), (5, y)])
            sage: all((pair[0],pair[1]) in rset for pair in p.coefficients(x))
            True
            sage: p.coefficients(x, sparse=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot return dense coefficient list with noninteger exponents

        Series coefficients are now handled correctly (:issue:`17399`)::


            sage: s = (1/(1-x)).series(x,6); s
            1 + 1*x + 1*x^2 + 1*x^3 + 1*x^4 + 1*x^5 + Order(x^6)
            sage: s.coefficients()
            [[1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5]]
            sage: s.coefficients(x, sparse=False)
            [1, 1, 1, 1, 1, 1]
            sage: x,y = var("x,y")
            sage: s = (1/(1-y*x-x)).series(x,3); s
            1 + (y + 1)*x + ((y + 1)^2)*x^2 + Order(x^3)
            sage: s.coefficients(x, sparse=False)
            [1, y + 1, (y + 1)^2]

        We can find coefficients of symbolic functions, :issue:`12255`::

            sage: g = function('g')(var('t'))
            sage: f = 3*g + g**2 + t
            sage: f.coefficients(g)
            [[t, 0], [3, 1], [1, 2]]

        Handle bound variable strictly as part of a constant::

            sage: (sin(1+x)*sin(1+x^2)).coefficients(x)
            [[sin(x^2 + 1)*sin(x + 1), 0]]
            sage: (sin(1+x)*sin(1+x^2)*x).coefficients(x)
            [[sin(x^2 + 1)*sin(x + 1), 1]]

        Check that :issue:`23545` is fixed::

            sage: (x^2/(1+x)).coefficients()
            [[x^2/(x + 1), 0]]
            sage: (1+x+exp(x^2/(1+x))).coefficients()
            [[e^(x^2/(x + 1)) + 1, 0], [1, 1]]
            sage: (1/x).coefficients()
            [[1, -1]]
            sage: ((1+x)^pi).coefficients()
            [[(x + 1)^pi, 0]]
        """
        cdef vector[pair[GEx,GEx]] vec
        cdef Expression xx
        if x is None:
            x = self.default_variable()
        xx = self.coerce_in(x)
        sig_on()
        try:
            self._gobj.coefficients(xx._gobj, vec)
        finally:
            sig_off()
        l = []
        for p in vec:
            l.append([new_Expression_from_GEx(self._parent, p.first),
                new_Expression_from_GEx(self._parent, p.second)])
        if sparse is True:
            return l

        from sage.rings.integer_ring import ZZ
        if any(not c[1] in ZZ for c in l):
            raise ValueError("cannot return dense coefficient list with noninteger exponents")
        if not l:
            l = [[0, 0]]
        val = l[0][1]
        if val < 0:
            raise ValueError("cannot return dense coefficient list with negative valuation")
        deg = l[-1][1]
        ret = [ZZ(0)] * int(deg+1)
        for c in l:
            ret[c[1]] = c[0]
        return ret

    def list(self, x=None):
        r"""
        Return the coefficients of this symbolic expression as a polynomial in x.

        INPUT:

        - ``x`` -- (optional) variable

        OUTPUT:

        A list of expressions where the ``n``-th element is the coefficient of
        ``x^n`` when ``self`` is seen as polynomial in ``x``.

        EXAMPLES::

            sage: var('x, y, a')
            (x, y, a)
            sage: (x^5).list()
            [0, 0, 0, 0, 0, 1]
            sage: p = x - x^3 + 5/7*x^5
            sage: p.list()
            [0, 1, 0, -1, 0, 5/7]
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: p.list(a)
            [x^2 + x + 1, -2*sqrt(2)*x, 2]
            sage: s = (1/(1-x)).series(x,6); s
            1 + 1*x + 1*x^2 + 1*x^3 + 1*x^4 + 1*x^5 + Order(x^6)
            sage: s.list()
            [1, 1, 1, 1, 1, 1]
        """
        return self.coefficients(x=x, sparse=False)

    def leading_coefficient(self, s):
        """
        Return the leading coefficient of ``s`` in ``self``.

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
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.lcoeff(ss._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    leading_coeff = leading_coefficient

    def trailing_coefficient(self, s):
        """
        Return the trailing coefficient of s in ``self``, i.e., the coefficient
        of the smallest power of s in ``self``.

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
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.tcoeff(ss._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    trailing_coeff = trailing_coefficient

    def low_degree(self, s):
        """
        Return the exponent of the lowest power of ``s`` in ``self``.

        OUTPUT: integer

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
            sage: (x+x**2).low_degree(x)
            1
        """
        cdef Expression ss = self.coerce_in(s)
        sig_on()
        try:
            return new_Expression_from_GEx(self._parent,
                                     GEx(self._gobj.ldegree(ss._gobj)))
        finally:
            sig_off()

    def degree(self, s):
        """
        Return the exponent of the highest power of ``s`` in ``self``.

        OUTPUT: integer

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
            sage: (1/x+1/x**2).degree(x)
            -1
        """
        cdef Expression ss = self.coerce_in(s)
        sig_on()
        try:
            return new_Expression_from_GEx(self._parent,
                                     GEx(self._gobj.degree(ss._gobj)))
        finally:
            sig_off()

    def unit(self, s):
        """
        Return the unit of this expression when considered as a
        polynomial in ``s``.

        See also :meth:`content`, :meth:`primitive_part`, and
        :meth:`unit_content_primitive`.

        INPUT:

        - ``s`` -- a symbolic expression

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
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.unit(ss._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    def content(self, s):
        """
        Return the content of this expression when considered as a
        polynomial in ``s``.

        See also :meth:`unit`, :meth:`primitive_part`, and
        :meth:`unit_content_primitive`.

        INPUT:

        - ``s`` -- a symbolic expression

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
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.content(ss._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    def primitive_part(self, s):
        """
        Return the primitive polynomial of this expression when
        considered as a polynomial in ``s``.

        See also :meth:`unit`, :meth:`content`, and
        :meth:`unit_content_primitive`.

        INPUT:

        - ``s`` -- a symbolic expression

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
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.primpart(ss._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

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
        sig_on()
        try:
            self._gobj.unitcontprim(ss._gobj, unit, cont, prim)
        finally:
            sig_off()
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
        if x is None:
            x = self.default_variable()
        G = self.coefficients(x)
        ans = None
        for Z in G:
            coeff = SR(Z[0])
            n = SR(Z[1])
            if not coeff.is_trivial_zero():
                if n.is_trivial_zero():
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

    def polynomial(self, base_ring=None, ring=None):
        r"""
        Return this symbolic expression as an algebraic polynomial
        over the given base ring, if possible.

        The point of this function is that it converts purely symbolic
        polynomials into optimised algebraic polynomials over a given
        base ring.

        You can specify either the base ring (``base_ring``) you want
        the output polynomial to be over, or you can specify the full
        polynomial ring (``ring``) you want the output polynomial to
        be an element of.

        INPUT:

        - ``base_ring`` -- (optional) the base ring for the polynomial

        - ``ring`` -- (optional) the parent for the polynomial

        .. warning::

           This is different from :meth:`poly` which is used to rewrite
           ``self`` as a polynomial in terms of one of the variables.

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
            sage: g.polynomial(QQ).monomial_coefficients()
            {0: -5, 2: 6}

        ::

            sage: f = x^2*e + x + pi/e
            sage: f.polynomial(RDF)  # abs tol 5e-16
            2.718281828459045*x^2 + x + 1.1557273497909217
            sage: g = f.polynomial(RR); g
            2.71828182845905*x^2 + x + 1.15572734979092
            sage: g.parent()
            Univariate Polynomial Ring in x over Real Field with 53 bits of precision
            sage: f.polynomial(RealField(100))
            2.7182818284590452353602874714*x^2 + x + 1.1557273497909217179100931833
            sage: f.polynomial(CDF)  # abs tol 5e-16
            2.718281828459045*x^2 + x + 1.1557273497909217
            sage: f.polynomial(CC)
            2.71828182845905*x^2 + x + 1.15572734979092

        We coerce a multivariate polynomial with complex symbolic
        coefficients::

            sage: x, y, n = var('x, y, n')
            sage: f = pi^3*x - y^2*e - I; f
            pi^3*x - y^2*e - I
            sage: f.polynomial(CDF)  # abs tol 1e-15
            (-2.718281828459045)*y^2 + 31.006276680299816*x - 1.0*I
            sage: f.polynomial(CC)
            (-2.71828182845905)*y^2 + 31.0062766802998*x - 1.00000000000000*I
            sage: f.polynomial(ComplexField(70))
            (-2.7182818284590452354)*y^2 + 31.006276680299820175*x - 1.0000000000000000000*I

        Another polynomial::

            sage: f = sum((e*I)^n*x^n for n in range(5)); f
            x^4*e^4 - I*x^3*e^3 - x^2*e^2 + I*x*e + 1
            sage: f.polynomial(CDF)   # abs tol 5e-16
            54.598150033144236*x^4 - 20.085536923187668*I*x^3 - 7.38905609893065*x^2
             + 2.718281828459045*I*x + 1.0
            sage: f.polynomial(CC)
            54.5981500331442*x^4 - 20.0855369231877*I*x^3 - 7.38905609893065*x^2
             + 2.71828182845905*I*x + 1.00000000000000

        A multivariate polynomial over a finite field::

            sage: f = (3*x^5 - 5*y^5)^7; f
            (3*x^5 - 5*y^5)^7
            sage: g = f.polynomial(GF(7)); g
            3*x^35 + 2*y^35
            sage: parent(g)
            Multivariate Polynomial Ring in x, y over Finite Field of size 7

        We check to make sure constants are converted appropriately::

            sage: (pi*x).polynomial(SR)
            pi*x

        Using the ``ring`` parameter, you can also create polynomials
        rings over the symbolic ring where only certain variables are
        considered generators of the polynomial ring and the others
        are considered "constants"::

            sage: a, x, y = var('a,x,y')
            sage: f = a*x^10*y+3*x
            sage: B = f.polynomial(ring=SR['x,y'])
            sage: B.coefficients()
            [a, 3]
        """
        from sage.symbolic.expression_conversions import polynomial
        return polynomial(self, base_ring=base_ring, ring=ring)

    def laurent_polynomial(self, base_ring=None, ring=None):
        r"""
        Return this symbolic expression as a Laurent polynomial
        over the given base ring, if possible.

        INPUT:

        - ``base_ring`` -- (optional) the base ring for the polynomial

        - ``ring`` -- (optional) the parent for the polynomial

        You can specify either the base ring (``base_ring``) you want
        the output Laurent polynomial to be over, or you can specify the full
        laurent polynomial ring (``ring``) you want the output laurent
        polynomial to be an element of.

        EXAMPLES::

            sage: f = x^2 -2/3/x + 1
            sage: f.laurent_polynomial(QQ)
            -2/3*x^-1 + 1 + x^2
            sage: f.laurent_polynomial(GF(19))
            12*x^-1 + 1 + x^2
        """
        from sage.symbolic.expression_conversions import laurent_polynomial
        return laurent_polynomial(self, base_ring=base_ring, ring=ring)

    def _polynomial_(self, R):
        """
        Coerce this symbolic expression to a polynomial in `R`.

        EXAMPLES::

            sage: var('x,y,z,w')
            (x, y, z, w)

        ::

            sage: R = QQ['x,y,z']
            sage: R(x^2 + y)
            x^2 + y
            sage: R = QQ['w']
            sage: R(w^3 + w + 1)
            w^3 + w + 1
            sage: R = GF(7)['z']
            sage: R(z^3 + 10*z)
            z^3 + 3*z

        .. NOTE::

           If the base ring of the polynomial ring is the symbolic ring,
           then a constant polynomial is always returned.

        ::

            sage: R = SR['x']
            sage: a = R(sqrt(2) + x^3 + y)
            sage: a
            x^3 + y + sqrt(2)
            sage: type(a)
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_field_with_category.element_class'>
            sage: a.degree()
            0

        We coerce to a double precision complex polynomial ring::

            sage: f = e*x^3 + pi*y^3 + sqrt(2) + I; f
            pi*y^3 + x^3*e + sqrt(2) + I
            sage: R = CDF['x,y']
            sage: R(f)  # abs tol 1e-15
            2.718281828459045*x^3 + 3.141592653589793*y^3 + 1.414213562373095 + 1.0*I

        We coerce to a higher-precision polynomial ring::

            sage: R = ComplexField(100)['x,y']
            sage: R(f)
            2.7182818284590452353602874714*x^3 + 3.1415926535897932384626433833*y^3
             + 1.4142135623730950488016887242 + 1.0000000000000000000000000000*I

        TESTS:

        This shows that the issue at :issue:`5755` is fixed (attempting to
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

        This shows that the issue at :issue:`4246` is fixed (attempting to
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
        from sage.symbolic.ring import SR
        from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
        base_ring = R.base_ring()
        if base_ring == SR:
            if isinstance(R, MPolynomialRing_base):
                return R({tuple([0]*R.ngens()):self})
            else:
                return R([self])
        return self.polynomial(None, ring=R)

    def fraction(self, base_ring):
        """
        Return this expression as element of the algebraic fraction
        field over the base ring given.

        EXAMPLES::

            sage: fr = (1/x).fraction(ZZ); fr
            1/x
            sage: parent(fr)
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: parent(((pi+sqrt(2)/x).fraction(SR)))
            Fraction Field of Univariate Polynomial Ring in x over Symbolic Ring
            sage: parent(((pi+sqrt(2))/x).fraction(SR))
            Fraction Field of Univariate Polynomial Ring in x over Symbolic Ring
            sage: y = var('y')
            sage: fr = ((3*x^5 - 5*y^5)^7/(x*y)).fraction(GF(7)); fr
            (3*x^35 + 2*y^35)/(x*y)
            sage: parent(fr)
            Fraction Field of Multivariate Polynomial Ring in x, y over Finite Field of size 7

        TESTS:

        Check that :issue:`17736` is fixed::

            sage: a,b,c = var('a,b,c')
            sage: fr = (1/a).fraction(QQ); fr
            1/a
            sage: parent(fr)
            Fraction Field of Univariate Polynomial Ring in a over Rational Field
            sage: parent((b/(a+sin(c))).fraction(SR))
            Fraction Field of Multivariate Polynomial Ring in a, b over Symbolic Ring
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.fraction_field import FractionField
        from sage.symbolic.ring import SR

        nu = SR(self.numerator()).polynomial(base_ring)
        de = SR(self.denominator()).polynomial(base_ring)
        vars = sorted(set(nu.variables() + de.variables()), key=repr)
        R = FractionField(PolynomialRing(base_ring, vars))
        return R(self.numerator())/R(self.denominator())

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
            raise ValueError("self must be a polynomial in one variable but it is in the variables %s" % tuple([v]))
        f = self.polynomial(base_ring)
        from sage.rings.power_series_ring import PowerSeriesRing
        R = PowerSeriesRing(base_ring, names=f.parent().variable_names())
        return R(f, f.degree()+1)

    def gcd(self, b):
        r"""
        Return the symbolic gcd of ``self`` and ``b``.

        Note that the polynomial GCD is unique up to the multiplication
        by an invertible constant. The following examples make sure all
        results are caught.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: SR(10).gcd(SR(15))
            5
            sage: (x^3 - 1).gcd(x-1) / (x-1) in QQ
            True
            sage: (x^3 - 1).gcd(x^2+x+1) / (x^2+x+1) in QQ
            True
            sage: (x^3 - x^2*pi + x^2 - pi^2).gcd(x-pi) / (x-pi) in QQ
            True
            sage: gcd(sin(x)^2 + sin(x), sin(x)^2 - 1) / (sin(x) + 1) in QQ
            True
            sage: gcd(x^3 - y^3, x-y) / (x-y) in QQ
            True
            sage: gcd(x^100-y^100, x^10-y^10) / (x^10-y^10) in QQ
            True
            sage: r = gcd(expand( (x^2+17*x+3/7*y)*(x^5 - 17*y + 2/3) ), expand((x^13+17*x+3/7*y)*(x^5 - 17*y + 2/3)) )
            sage: r / (x^5 - 17*y + 2/3) in QQ
            True

        Embedded Sage objects of all kinds get basic support. Note that
        full algebraic GCD is not implemented yet::

            sage: gcd(I - I*x, x^2 - 1)
            x - 1
            sage: gcd(I + I*x, x^2 - 1)
            x + 1
            sage: alg = SR(QQbar(sqrt(2) + I*sqrt(3)))
            sage: gcd(alg + alg*x, x^2 - 1)  # known bug (Issue #28489)
            x + 1
            sage: gcd(alg - alg*x, x^2 - 1)  # known bug (Issue #28489)
            x - 1
            sage: sqrt2 = SR(QQbar(sqrt(2)))
            sage: gcd(sqrt2 + x, x^2 - 2)    # known bug
            1

        TESTS:

        Check if :issue:`10284` is fixed::

            sage: u = var('u')
            sage: v = var('v')
            sage: w = var('w')
            sage: x = var('x')
            sage: y = var('y')
            sage: z = var('z')
            sage: e = 792*z^8*w^4*x^3*y^4*u^7 + 24*z^4*w^4*x^2*y^3*u^4 + \
                    264*z^8*w^3*x^2*y^7*u^5 + 198*z^4*w^5*x^5*y*u^6  + 110*z^2*w^3*x^5*y^4*u^6 \
                    - 120*z^8*w*x^4*u^6 - 480*z^5*w*x^4*y^6*u^8 - 720*z^7*x^3*y^3*u^7 + \
                    165*z^4*w^2*x^4*y*u^5 + 450*z^8*w^6*x^2*y*u^8 + 40*z^2*w^3*x^3*y^3*u^6 - \
                    288*z^7*w^2*x^3*y^6*u^6  + 250*z^6*w^4*x^2*y^4*u^8 + \
                    576*z^7*w^7*x^2*y^4*u^8  - 80*z^6*w^2*x^5*y^3*u^7 - 144*z^8*w^4*x^5*u^7 + \
                    120*z^4*w*x^2*y^6*u^6 + 320*z^5*w^5*x^2*y^7*u^8 + 192*z^7*w^6*x*y^7*u^6 - \
                    12*z^4*w^3*x^3*y^5*u^6  - 36*z^4*w^4*x^4*y^2*u^8 + 72*z^4*w^5*x^3*u^6  - \
                    20*z^2*w^2*x^4*y^5*u^8 + 660*z^8*w*x^2*y^4*u^6 + 66*z^4*w^4*x^4*y^4*u^4 + \
                    440*z^6*w^2*x^3*y^7*u^7  - 30*z^4*w*x^3*y^2*u^7 - 48*z^8*w^3*x^4*y^3*u^5 + \
                    72*z^6*w^2*x*y^6*u^4 - 864*z^7*w^3*x^4*y^3*u^8 + 480*z^7*w^4*x*y^4*u^7 + \
                    60*z^4*w^2*x^2*u^5 + 375*z^8*w^3*x*y*u^7 + 150*z^8*w^5*x*y^4*u^6 + \
                    180*z^6*x*y^3*u^5 + 216*z^6*w^3*x^2*y^3*u^6;
            sage: d = e.diff(x)
            sage: gcd(d,e) / (u^4*z^2) in QQ
            True

        Check that :issue:`23793` is fixed::

            sage: gcd(I + I*x, x^2 - 1)
            x + 1

        Check that arguments are expanded before GCD (:issue:`23845`)::

            sage: P = (x+1)^2 + 1
            sage: gcd(P, P.expand())
            x^2 + 2*x + 2
        """
        cdef Expression r = self.coerce_in(b)
        cdef GEx x = g_gcd(self._gobj, r._gobj)
        return new_Expression_from_GEx(self._parent, x)

    def gosper_sum(self, *args):
        """
        Return the summation of this hypergeometric expression using
        Gosper's algorithm.

        INPUT:

        - a symbolic expression that may contain rational functions,
          powers, factorials, gamma function terms, binomial
          coefficients, and Pochhammer symbols that are rational-linear
          in their arguments

        - the main variable and, optionally, summation limits

        EXAMPLES::

            sage: a,b,k,m,n = var('a b k m n')
            sage: SR(1).gosper_sum(n)
            n
            sage: SR(1).gosper_sum(n,5,8)
            4
            sage: n.gosper_sum(n)
            1/2*(n - 1)*n
            sage: n.gosper_sum(n,0,5)
            15
            sage: n.gosper_sum(n,0,m)
            1/2*(m + 1)*m
            sage: n.gosper_sum(n,a,b)
            -1/2*(a + b)*(a - b - 1)

        ::

            sage: (factorial(m + n)/factorial(n)).gosper_sum(n)
            n*factorial(m + n)/((m + 1)*factorial(n))
            sage: (binomial(m + n, n)).gosper_sum(n)
            n*binomial(m + n, n)/(m + 1)
            sage: (binomial(m + n, n)).gosper_sum(n, 0, a)
            (a + m + 1)*binomial(a + m, a)/(m + 1)
            sage: (binomial(m + n, n)).gosper_sum(n, 0, 5)
            1/120*(m + 6)*(m + 5)*(m + 4)*(m + 3)*(m + 2)
            sage: (rising_factorial(a,n)/rising_factorial(b,n)).gosper_sum(n)
            (b + n - 1)*gamma(a + n)*gamma(b)/((a - b + 1)*gamma(a)*gamma(b + n))
            sage: factorial(n).gosper_term(n)
            Traceback (most recent call last):
            ...
            ValueError: expression not Gosper-summable
        """
        cdef Expression s, a, b
        cdef GEx x
        cdef int i = 1
        if len(args) > 1:
            n, l1, l2 = args
            s = self.coerce_in(n)
            a = self.coerce_in(l1)
            b = self.coerce_in(l2)
            sig_on()
            try:
                x = g_gosper_sum_definite(self._gobj, s._gobj,
                        a._gobj, b._gobj, &i)
            finally:
                sig_off()
            if i == 0:
                raise ValueError("expression not Gosper-summable")
            return new_Expression_from_GEx(self._parent, x)
        else:
            s = self.coerce_in(args[0])
            sig_on()
            try:
                x = g_gosper_sum_indefinite(self._gobj, s._gobj, &i)
            finally:
                sig_off()
            if i == 0:
                raise ValueError("expression not Gosper-summable")
            return new_Expression_from_GEx(self._parent, x)

    def gosper_term(self, n):
        r"""
        Return Gosper's hypergeometric term for ``self``.

        Suppose ``f``=``self`` is a hypergeometric term such that:

        .. math::

            s_n = \sum_{k=0}^{n-1} f_k

        and `f_k` doesn't depend on `n`. Return a hypergeometric
        term `g_n` such that `g_{n+1} - g_n = f_n`.

        EXAMPLES::

            sage: _ = var('n')
            sage: SR(1).gosper_term(n)
            n
            sage: n.gosper_term(n)
            1/2*(n^2 - n)/n
            sage: (n*factorial(n)).gosper_term(n)
            1/n
            sage: factorial(n).gosper_term(n)
            Traceback (most recent call last):
            ...
            ValueError: expression not Gosper-summable
        """
        cdef Expression s
        s = self.coerce_in(n)
        sig_on()
        try:
            x = g_gosper_term(self._gobj, s._gobj)
        except ValueError:
            raise ValueError("expression not Gosper-summable")
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def WZ_certificate(self, n, k):
        r"""
        Return the Wilf-Zeilberger certificate for this hypergeometric
        summand in ``n``, ``k``.

        To prove the identity `\sum_k F(n,k)=\textrm{const}` it suffices
        to show that `F(n+1,k)-F(n,k)=G(n,k+1)-G(n,k),` with `G=RF` and
        `R` the WZ certificate.

        EXAMPLES:

        To show that `\sum_k \binom{n}{k} = 2^n` do::

            sage: _ = var('k n')
            sage: F(n,k) = binomial(n,k) / 2^n
            sage: c = F(n,k).WZ_certificate(n,k); c
            1/2*k/(k - n - 1)
            sage: G(n,k) = c * F(n,k); G
            (n, k) |--> 1/2*k*binomial(n, k)/(2^n*(k - n - 1))
            sage: (F(n+1,k) - F(n,k) - G(n,k+1) + G(n,k)).simplify_full()
            0
        """
        cdef Expression s, f
        s = self.coerce_in(k)
        f = self.subs(n == n+1) - self
        sig_on()
        try:
            x = g_gosper_term(f._gobj, s._gobj)
        finally:
            sig_off()
        g = new_Expression_from_GEx(self._parent, x)
        return (f*g / self).to_gamma().gamma_normalize().simplify_full().factor()

    def lcm(self, b):
        """
        Return the lcm of ``self`` and ``b``.

        The lcm is computed from the gcd of ``self`` and ``b`` implicitly from
        the relation ``self * b = gcd(self, b) * lcm(self, b)``.

        .. NOTE::

            In agreement with the convention in use for integers, if
            ``self * b == 0``, then ``gcd(self, b) == max(self, b)`` and
            ``lcm(self, b) == 0``.

        .. NOTE::

            Since the polynomial lcm is computed from the gcd, and the
            polynomial gcd is unique up to a constant factor (which can
            be negative), the polynomial lcm is unique up to a factor of -1.

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
            (pi - x^3)*(pi - x)
            sage: lcm(x^3 - y^3, x-y) / (x^3 - y^3) in [1,-1]
            True
            sage: lcm(x^100-y^100, x^10-y^10) / (x^100 - y^100) in [1,-1]
            True
            sage: a = expand( (x^2+17*x+3/7*y)*(x^5 - 17*y + 2/3) )
            sage: b = expand((x^13+17*x+3/7*y)*(x^5 - 17*y + 2/3) )
            sage: gcd(a,b) * lcm(a,b) / (a * b) in [1,-1]
            True

        The result is not automatically simplified::

            sage: ex = lcm(sin(x)^2 - 1, sin(x)^2 + sin(x)); ex
            (sin(x)^2 + sin(x))*(sin(x)^2 - 1)/(sin(x) + 1)
            sage: ex.simplify_full()
            sin(x)^3 - sin(x)

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
            return 0 if sb.is_trivial_zero() else sb / self.gcd(b)
        except ValueError:
            # make the error message refer to lcm, not gcd
            raise ValueError("lcm: arguments must be polynomials over the rationals")

    def resultant(self, other, var):
        """
        Compute the resultant of this polynomial expression and the first
        argument with respect to the variable given as the second
        argument.

        EXAMPLES::

            sage: _ = var('a b n k u x y')
            sage: x.resultant(y, x)
            y
            sage: (x+y).resultant(x-y, x)
            -2*y
            sage: r = (x^4*y^2+x^2*y-y).resultant(x*y-y*a-x*b+a*b+u,x)
            sage: r.coefficient(a^4)
            b^4*y^2 - 4*b^3*y^3 + 6*b^2*y^4 - 4*b*y^5 + y^6
            sage: x.resultant(sin(x), x)
            Traceback (most recent call last):
            ...
            RuntimeError: resultant(): arguments must be polynomials
        """
        cdef Expression o = self.coerce_in(other)
        cdef Expression v = self.coerce_in(var)
        cdef GEx x
        sig_on()
        try:
            x = g_resultant(self._gobj, o._gobj, v._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def collect(Expression self, s):
        """
        Collect the coefficients of ``s`` into a group.

        INPUT:

        - ``s`` -- the symbol whose coefficients will be collected

        OUTPUT:

        A new expression, equivalent to the original one, with the
        coefficients of ``s`` grouped.

        .. NOTE::

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
        we do not factor the `y^{2}` and `z^{2}` terms before
        collecting coefficients::

            sage: f.collect(y)
            (x^2*z^2 + 20)*y^2 + (4*x + 21*z)*y + x*z + 4*z^2
            sage: f.collect(z)
            (x^2*y^2 + 4)*z^2 + 4*x*y + 20*y^2 + (x + 21*y)*z

        The terms are collected, whether the expression
        is expanded or not::

            sage: f = (x + y)*(x - z)
            sage: f.collect(x)
            x^2 + x*(y - z) - y*z
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

        Check if :issue:`9046` is fixed::

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

    def horner(self, x):
        """
        Rewrite this expression as a polynomial in Horner form in ``x``.

        EXAMPLES::

            sage: add((i+1)*x^i for i in range(5)).horner(x)
            (((5*x + 4)*x + 3)*x + 2)*x + 1

            sage: x, y, z = SR.var('x,y,z')
            sage: (x^5 + y*cos(x) + z^3 + (x + y)^2 + y^x).horner(x)
            z^3 + ((x^3 + 1)*x + 2*y)*x + y^2 + y*cos(x) + y^x

            sage: expr = sin(5*x).expand_trig(); expr
            5*cos(x)^4*sin(x) - 10*cos(x)^2*sin(x)^3 + sin(x)^5
            sage: expr.horner(sin(x))
            (5*cos(x)^4 - (10*cos(x)^2 - sin(x)^2)*sin(x)^2)*sin(x)
            sage: expr.horner(cos(x))
            sin(x)^5 + 5*(cos(x)^2*sin(x) - 2*sin(x)^3)*cos(x)^2

        TESTS::

            sage: SR(0).horner(x), SR(1).horner(x), x.horner(x)
            (0, 1, x)
            sage: (x^(1/3)).horner(x)
            Traceback (most recent call last):
            ...
            ValueError: cannot return dense coefficient list with noninteger exponents
        """
        coef = self.coefficients(x, sparse=False)
        res = coef[-1]
        for c in reversed(coef[:-1]):
            res = res*x + c
        return res

    def _evaluate_polynomial(self, pol):
        """
        Evaluate a univariate polynomial on this expression.

        EXAMPLES::

            sage: pol = QQ['s'](range(5))
            sage: pol(sin(x))
            4*sin(x)^4 + 3*sin(x)^3 + 2*sin(x)^2 + sin(x)

        TESTS::

            sage: SR(0.1)._evaluate_polynomial(pol)
            0.123400000000000

            sage: Pol.<x> = SR[]
            sage: pol = x^2 - 1
            sage: pol(1)
            0
            sage: pol(1).parent()
            Symbolic Ring
        """
        cdef Expression zero
        if not isinstance(pol[0], Expression):  # avoid infinite recursion
            try:
                x = self.pyobject()
                y = pol(x)  # may fail if x is a symbolic constant
            except TypeError:
                pass
            else:
                return new_Expression_from_pyobject(self._parent, y)
        zero = self._parent.zero()
        return zero.add(*(pol[i]*self**i for i in range(pol.degree() + 1)))

    def collect_common_factors(self):
        """
        This function does not perform a full factorization but only
        looks for factors which are already explicitly present.

        Polynomials can often be brought into a more compact form by
        collecting common factors from the terms of sums. This is
        accomplished by this function.

        EXAMPLES::

            sage: var('x')
            x
            sage: (x/(x^2 + x)).collect_common_factors()
            1/(x + 1)

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: (a*x+a*y).collect_common_factors()
            a*(x + y)
            sage: (a*x^2+2*a*x*y+a*y^2).collect_common_factors()
            (x^2 + 2*x*y + y^2)*a
            sage: (a*(b*(a+c)*x+b*((a+c)*x+(a+c)*y)*y)).collect_common_factors()
            ((x + y)*y + x)*(a + c)*a*b
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
            <class 'sage.symbolic.expression.Expression'>

        Because this overrides a Python builtin function, we do not
        currently support a ``hold`` parameter to prevent automatic
        evaluation::

            sage: abs(SR(-5),hold=True)
            Traceback (most recent call last):
            ...
            TypeError: ...abs() takes no keyword arguments

        But this is possible using the method :meth:`abs`::

            sage: SR(-5).abs(hold=True)
            abs(-5)

        TESTS:

        Check if :issue:`11155` is fixed::

            sage: abs(pi+i)
            abs(pi + I)
        """
        cdef GEx r
        sig_on()
        try:
            r = g_abs(self._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(-5).abs(hold=True); a.unhold()
            5

        TESTS:

        From :issue:`7557`::

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
        Return the value of the unit step function, which is 0 for
        negative x, 1 for 0, and 1 for positive x.

        .. SEEALSO::

            :class:`sage.functions.generalized.FunctionUnitStep`

        EXAMPLES::

            sage: x = var('x')
            sage: SR(1.5).step()
            1
            sage: SR(0).step()
            1
            sage: SR(-1/2).step()
            0
            sage: SR(float(-1)).step()
            0

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(2).step()
            1
            sage: SR(2).step(hold=True)
            unit_step(2)
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_step, self._gobj, hold))

    def csgn(self, hold=False):
        """
        Return the sign of ``self``, which is -1 if ``self < 0``, 0 if
        ``self == 0``, and 1 if ``self > 0``, or unevaluated when ``self`` is a
        nonconstant symbolic expression.

        If ``self`` is not real, return the complex half-plane (left or right)
        in which the number lies.  If ``self`` is pure imaginary, return the sign
        of the imaginary part of ``self``.

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(I).conjugate(hold=True); a.unhold()
            -I
        """
        return new_Expression_from_GEx(self._parent,
                g_hold_wrapper(g_conjugate, self._gobj, hold))

    def norm(self):
        r"""
        Return the complex norm of this symbolic expression, i.e.,
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

            :func:`sage.misc.functional.norm`

        EXAMPLES::

            sage: a = 1 + 2*I
            sage: a.norm()
            5
            sage: a = sqrt(2) + 3^(1/3)*I; a
            sqrt(2) + I*3^(1/3)
            sage: a.norm()
            3^(2/3) + 2
            sage: CDF(a).norm()
            4.080083823051...
            sage: CDF(a.norm())
            4.080083823051904
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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(2).real_part(hold=True); a.unhold()
            2

        TESTS:

        Check that :issue:`12807` is fixed::

            sage: (6*exp(i*pi/3)-6*exp(i*2*pi/3)).real_part()
            6

        Check that :issue:`28357` is fixed::

            sage: m = var('m')
            sage: assume(m, 'integer')
            sage: (I^m).real_part()
            cos(1/2*pi*m)
            sage: (I^m).imag_part()
            sin(1/2*pi*m)
            sage: forget()

        Check that :issue:`29400` is fixed::

            sage: cot(1 + i).imag().n() - (1/tan(1 + i)).imag().n()  # abs tol 10^-12
            0.00000000000000
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

            sage: SR(I).imag_part()
            1
            sage: SR(I).imag_part(hold=True)
            imag_part(I)

        This also works using functional notation::

            sage: imag_part(I, hold=True)
            imag_part(I)
            sage: imag_part(SR(I))
            1

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(I).imag_part(hold=True); a.unhold()
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
        Return the square root of this expression.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: SR(2).sqrt()
            sqrt(2)
            sage: (x^2+y^2).sqrt()
            sqrt(x^2 + y^2)
            sage: (x^2).sqrt()
            sqrt(x^2)

        Immediate simplifications are applied::

            sage: sqrt(x^2)
            sqrt(x^2)
            sage: x = SR.symbol('x', domain='real')
            sage: sqrt(x^2)
            abs(x)
            sage: forget()
            sage: assume(x<0)
            sage: sqrt(x^2)
            -x
            sage: sqrt(x^4)
            x^2
            sage: forget()
            sage: x = SR.symbol('x', domain='real')
            sage: sqrt(x^4)
            x^2
            sage: sqrt(sin(x)^2)
            abs(sin(x))
            sage: sqrt((x+1)^2)
            abs(x + 1)
            sage: forget()
            sage: assume(x<0)
            sage: sqrt((x-1)^2)
            -x + 1
            sage: forget()

        Using the ``hold`` parameter it is possible to prevent automatic
        evaluation::

            sage: SR(4).sqrt()
            2
            sage: SR(4).sqrt(hold=True)
            sqrt(4)

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(4).sqrt(hold=True); a.unhold()
            2

        To use this parameter in functional notation, you must coerce to
        the symbolic ring::

            sage: sqrt(SR(4),hold=True)
            sqrt(4)
            sage: sqrt(4,hold=True)
            Traceback (most recent call last):
            ...
            TypeError: ..._do_sqrt() got an unexpected keyword argument 'hold'
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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(0).sin(hold=True); a.unhold()
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
        Return the cosine of ``self``.

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = pi.cos(hold=True); a.unhold()
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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = (pi/12).tan(hold=True); a.unhold()
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
            1/6*pi
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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(0).arcsin(hold=True); a.unhold()
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
        Return the arc cosine of ``self``.

        EXAMPLES::

            sage: x.arccos()
            arccos(x)
            sage: SR(1).arccos()
            0
            sage: SR(1/2).arccos()
            1/3*pi
            sage: SR(0.4).arccos()
            1.15927948072741
            sage: plot(lambda x: SR(x).arccos(), -1,1)                                  # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1).arccos(hold=True)
            arccos(1)

        This also works using functional notation::

            sage: arccos(1, hold=True)
            arccos(1)
            sage: arccos(1)
            0

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(1).arccos(hold=True); a.unhold()
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
        Return the arc tangent of ``self``.

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
            sage: plot(lambda x: SR(x).arctan(), -20,20)                                # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1).arctan(hold=True)
            arctan(1)

        This also works using functional notation::

            sage: arctan(1, hold=True)
            arctan(1)
            sage: arctan(1)
            1/4*pi

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(1).arctan(hold=True); a.unhold()
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
        Return the inverse of the 2-variable tan function on ``self`` and ``x``.

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(1/2).arctan2(1/2, hold=True); a.unhold()
            1/4*pi

        TESTS:

        We compare a bunch of different evaluation points between
        Sage and Maxima::

            sage: float(SR(0.7).arctan2(0.6))
            0.8621700546672264
            sage: maxima('atan2(0.7,0.6)')
            0.862170054667226...
            sage: float(SR(0.7).arctan2(-0.6))
            2.279422598922567
            sage: maxima('atan2(0.7,-0.6)')
            2.279422598922567
            sage: float(SR(-0.7).arctan2(0.6))
            -0.8621700546672264
            sage: maxima('atan2(-0.7,0.6)')
            -0.862170054667226...
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
            sage: SR(0).arctan2(0) # see github issue #21614
            NaN
            sage: SR(I).arctan2(1)
            arctan2(I, 1)
            sage: SR(CDF(0,1)).arctan2(1)
            Traceback (most recent call last):
            ...
            ValueError: power::eval(): division by zero
            sage: SR(1).arctan2(CDF(0,1))
            Traceback (most recent call last):
            ...
            ValueError: power::eval(): division by zero

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
        Return sinh of ``self``.

        We have `\sinh(x) = (e^{x} - e^{-x})/2`.

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

        To then evaluate again, we use :meth:`unhold`::

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
        Return cosh of ``self``.

        We have `\cosh(x) = (e^{x} + e^{-x})/2`.

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = arcsinh(x).cosh(hold=True); a.unhold()
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
        Return tanh of ``self``.

        We have `\tanh(x) = \sinh(x) / \cosh(x)`.

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
            0.7615941559557649
            sage: plot(lambda x: SR(x).tanh(), -1, 1)                                   # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arcsinh(x).tanh()
            x/sqrt(x^2 + 1)
            sage: arcsinh(x).tanh(hold=True)
            tanh(arcsinh(x))

        This also works using functional notation::

            sage: tanh(arcsinh(x), hold=True)
            tanh(arcsinh(x))
            sage: tanh(arcsinh(x))
            x/sqrt(x^2 + 1)

        To then evaluate again, we use :meth:`unhold`::

            sage: a = arcsinh(x).tanh(hold=True); a.unhold()
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
        Return the inverse hyperbolic sine of ``self``.

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(-2).arcsinh(hold=True); a.unhold()
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
        Return the inverse hyperbolic cosine of ``self``.

        EXAMPLES::

            sage: x.arccosh()
            arccosh(x)
            sage: SR(0).arccosh()
            1/2*I*pi
            sage: SR(1/2).arccosh()
            arccosh(1/2)
            sage: SR(CDF(1/2)).arccosh() #  rel tol 1e-15
            1.0471975511965976*I
            sage: z = maxima('acosh(0.5)')
            sage: z.real(), z.imag()  # abs tol 1e-15
            (0.0, 1.047197551196598)

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(-1).arccosh(hold=True); a.unhold()
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
        Return the inverse hyperbolic tangent of ``self``.

        EXAMPLES::

            sage: x.arctanh()
            arctanh(x)
            sage: SR(0).arctanh()
            0
            sage: SR(1/2).arctanh()
            1/2*log(3)
            sage: SR(0.5).arctanh()
            0.549306144334055
            sage: SR(0.5).arctanh().tanh()
            0.500000000000000
            sage: maxima('atanh(0.5)')  # abs tol 2e-16
            0.5493061443340548

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(-1/2).arctanh()
            -1/2*log(3)
            sage: SR(-1/2).arctanh(hold=True)
            arctanh(-1/2)

        This also works using functional notation::

            sage: arctanh(-1/2,hold=True)
            arctanh(-1/2)
            sage: arctanh(-1/2)
            -1/2*log(3)

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(-1/2).arctanh(hold=True); a.unhold()
            -1/2*log(3)

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
        r"""
        Return exponential function of ``self``, i.e., `e` to the
        power of ``self``.

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = (pi*I).exp(hold=True); a.unhold()
            -1

        TESTS:

        Test if :issue:`6377` is fixed::

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
        Return the logarithm of ``self``.

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
            sage: plot(lambda x: SR(x).log(), 0.1,10)                                   # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        To prevent automatic evaluation use the ``hold`` argument::

            sage: I.log()
            1/2*I*pi
            sage: I.log(hold=True)
            log(I)

        To then evaluate again, we use :meth:`unhold`::

            sage: a = I.log(hold=True); a.unhold()
            1/2*I*pi

        The ``hold`` parameter also works in functional notation::

            sage: log(-1, hold=True)
            log(-1)
            sage: log(-1)
            I*pi

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
            sage: SR(CDF(0,1)).zeta()  # abs tol 1e-16                                  # needs sage.libs.pari
            0.003300223685324103 - 0.4181554491413217*I
            sage: CDF(0,1).zeta()  # abs tol 1e-16                                      # needs sage.libs.pari
            0.003300223685324103 - 0.4181554491413217*I
            sage: plot(lambda x: SR(x).zeta(), -10,10).show(ymin=-3, ymax=3)            # needs sage.plot

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(2).zeta(hold=True)
            zeta(2)

        This also works using functional notation::

            sage: zeta(2, hold=True)
            zeta(2)
            sage: zeta(2)
            1/6*pi^2

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(2).zeta(hold=True); a.unhold()
            1/6*pi^2

        TESTS::

            sage: t = SR(1).zeta(); t
            Infinity
        """
        cdef GEx x = g_hold_wrapper(g_zeta, self._gobj, hold)
        return new_Expression_from_GEx(self._parent, x)

    def factorial(self, hold=False):
        """
        Return the factorial of ``self``.

        OUTPUT: a symbolic expression

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

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(5).factorial(hold=True); a.unhold()
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

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: SR(5).binomial(SR(3))
            10
            sage: x.binomial(SR(3))
            1/6*(x - 1)*(x - 2)*x
            sage: x.binomial(y)
            binomial(x, y)

        To prevent automatic evaluation use the ``hold`` argument::

            sage: x.binomial(3, hold=True)
            binomial(x, 3)
            sage: SR(5).binomial(3, hold=True)
            binomial(5, 3)

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(5).binomial(3, hold=True); a.unhold()
            10

        The ``hold`` parameter is also supported in functional notation::

            sage: binomial(5,3, hold=True)
            binomial(5, 3)

        TESTS:

        Check if we handle zero correctly (:issue:`8561`)::

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
        Return the order of the expression, as in big oh notation.

        OUTPUT: a symbolic expression

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

    def gamma(self, *, hold=False):
        """
        Return the Gamma function evaluated at ``self``.

        EXAMPLES::

            sage: x = var('x')
            sage: x.gamma()
            gamma(x)
            sage: SR(2).gamma()
            1
            sage: SR(10).gamma()
            362880
            sage: SR(10.0r).gamma()  # For ARM: rel tol 2e-15
            362880.0
            sage: SR(CDF(1,1)).gamma()
            0.49801566811835607 - 0.15494982830181067*I

        ::

            sage: gp('gamma(1+I)')
            0.49801566811835604271369111746219809195 - 0.15494982830181068512495513048388660520*I

        We plot the familiar plot of this log-convex function::

            sage: plot(gamma(x), -6, 4).show(ymin=-3, ymax=3)                           # needs sage.plot

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(1/2).gamma()
            sqrt(pi)
            sage: SR(1/2).gamma(hold=True)
            gamma(1/2)

        This also works using functional notation::

            sage: gamma(1/2, hold=True)
            gamma(1/2)
            sage: gamma(1/2)
            sqrt(pi)

        To then evaluate again, we use :meth:`unhold`::

            sage: a = SR(1/2).gamma(hold=True); a.unhold()
            sqrt(pi)

        TESTS:

        Check that no confusion with the incomplete gamma function is
        possible::

            sage: x, y = SR.var('x,y')
            sage: x.gamma(y)
            Traceback (most recent call last):
            ...
            TypeError: ...gamma() takes exactly 0 positional arguments (1 given)
        """
        cdef GEx x
        sig_on()
        try:
            x = g_hold_wrapper(g_gamma, self._gobj, hold)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def log_gamma(self, hold=False):
        """
        Return the log gamma function evaluated at ``self``.
        This is the logarithm of gamma of ``self``, where
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
            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(-1)
            sage: plot(lambda x: SR(x).log_gamma(), -7,8, plot_points=1000).show()      # needs sage.plot
            sage: math.exp(0.5)
            1.6487212707001282
            sage: plot(lambda x: (SR(x).exp() - SR(-x).exp())/2 - SR(x).sinh(), -1, 1)  # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        To prevent automatic evaluation use the ``hold`` argument::

            sage: SR(5).log_gamma(hold=True)
            log_gamma(5)

        To evaluate again, currently we must use numerical evaluation
        via :meth:`n`::

            sage: a = SR(5).log_gamma(hold=True); a.n()
            3.17805383034795
        """
        # Note: g_lgamma calls back into Python, must not wrap in sig_on/sig_off
        cdef GEx x = g_hold_wrapper(g_lgamma, self._gobj, hold)
        return new_Expression_from_GEx(self._parent, x)

    def default_variable(self):
        """
        Return the default variable, which is by definition the first
        variable in ``self``, or `x` if there are no variables in ``self``.
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

    def combine(self, bint deep=False):
        r"""
        Return a simplified version of this symbolic expression
        by combining all toplevel terms with the same denominator into
        a single term.

        Please use the keyword ``deep=True`` to apply the process
        recursively.

        EXAMPLES::

            sage: var('x, y, a, b, c')
            (x, y, a, b, c)
            sage: f = x*(x-1)/(x^2 - 7) + y^2/(x^2-7) + 1/(x+1) + b/a + c/a; f
            (x - 1)*x/(x^2 - 7) + y^2/(x^2 - 7) + b/a + c/a + 1/(x + 1)
            sage: f.combine()
            ((x - 1)*x + y^2)/(x^2 - 7) + (b + c)/a + 1/(x + 1)
            sage: (1/x + 1/x^2 + (x+1)/x).combine()
            (x + 2)/x + 1/x^2
            sage: ex = 1/x + ((x + 1)/x - 1/x)/x^2 + (x+1)/x; ex
            (x + 1)/x + 1/x + ((x + 1)/x - 1/x)/x^2
            sage: ex.combine()
            (x + 2)/x + ((x + 1)/x - 1/x)/x^2
            sage: ex.combine(deep=True)
            (x + 2)/x + 1/x^2
            sage: (1+sin((x + 1)/x - 1/x)).combine(deep=True)
            sin(1) + 1
        """
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.combine_fractions(deep)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    def normalize(self):
        """
        Return this expression normalized as a fraction.

        .. SEEALSO::

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

        TESTS:

        Check that :issue:`19775` is fixed::

            sage: a,b,c,d,e,y = var('a,b,c,d,e,y')
            sage: ((x - 2*y)^4/(x^2 - 4*y^2)^2).normalize()
            (x - 2*y)^2/(x + 2*y)^2
            sage: f = ((x - 2*y)^4/(x^2 - 4*y^2)^2 + 1)*(y + a)*(2*y + x) / (4*y^2 + x^2)
            sage: f.normalize()
            2*(a + y)/(x + 2*y)
            sage: (c/a - b*c^2/(a^2*(b*c/a-d)) + c*d/(a*(b*c/a-d))).normalize()
            0
            sage: (e + c/a - b*c^2/(a^2*(b*c/a-d)) + c*d/(a*(b*c/a-d))).normalize()
            e

        Check that :issue:`23861` is fixed::

            sage: (x^(2*pi) + x^(-2*pi) - 2).normalize()
            (x^(4*pi) - 2*x^(2*pi) + 1)/x^(2*pi)
            sage: (e^2 + e^(-2) - 2).normalize()
            (e^4 - 2*e^2 + 1)/e^2
            sage: (e^(2*pi) - e^(-2*pi)).normalize()
            (e^(4*pi) - 1)/e^(2*pi)

        ALGORITHM: Uses GiNaC.
        """
        cdef GEx r
        sig_on()
        try:
            r = self._gobj.normal(0, False, True)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, r)

    def numerator(self, bint normalize=True):
        """
        Return the numerator of this symbolic expression.

        INPUT:

        - ``normalize`` -- boolean (default: ``True``)

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
            sage: n = var('n'); assume(n,'integer'); assume(n>0); (e^(2*n)/(e^(2*n) - 1)).numerator()
            e^(2*n)
        """
        cdef GExVector vec
        cdef GEx oper, power, ex
        if normalize:
            sig_on()
            try:
                ex = self._gobj.numer()
            finally:
                sig_off()
            return new_Expression_from_GEx(self._parent, ex)
        elif is_a_mul(self._gobj):
            for i in range(self._gobj.nops()):
                oper = self._gobj.op(i)
                if not is_a_power(oper):
                    vec.push_back(oper)
                else:
                    power = oper.op(1)
                    if not is_a_numeric(power):
                        raise TypeError("self is not a rational expression")
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
        Return the denominator of this symbolic expression.

        INPUT:

        - ``normalize`` -- boolean (default: ``True``)

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
            sig_on()
            try:
                ex = self._gobj.denom()
            finally:
                sig_off()
            return new_Expression_from_GEx(self._parent, ex)
        elif is_a_mul(self._gobj):
            for i in range(self._gobj.nops()):
                oper = self._gobj.op(i)
                if is_a_power(oper):
                    ex = oper.op(0)
                    power = oper.op(1)
                    if not is_a_numeric(power):
                        raise TypeError("self is not a rational expression")
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
        Return the numerator and the denominator of this symbolic expression.

        INPUT:

        - ``normalize`` -- boolean (default: ``True``)

        If ``normalize`` is ``True``, the expression is first normalized to
        have it as a fraction before getting the numerator and denominator.

        If ``normalize`` is ``False``, the expression is kept and if it is not
        a quotient, then this will return the expression itself together with
        1.

        .. SEEALSO::

            :meth:`normalize`, :meth:`numerator`, :meth:`denominator`,
            :meth:`combine`

        EXAMPLES::

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
            sig_on()
            try:
                ex = self._gobj.numer_denom()
            finally:
                sig_off()
            return (new_Expression_from_GEx(self._parent, ex.op(0)),
                    new_Expression_from_GEx(self._parent, ex.op(1)))
        elif is_a_mul(self._gobj):
            for i in range(self._gobj.nops()):
                oper = self._gobj.op(i)
                if is_a_power(oper):   # oper = ex^power
                    ex = oper.op(0)
                    power = oper.op(1)
                    if not is_a_numeric(power):
                        raise TypeError("self is not a rational expression")
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

        - ``var`` -- variable name or string (default: first variable)

        OUTPUT: a symbolic expression

        .. SEEALSO:: :meth:`partial_fraction_decomposition`

        EXAMPLES::

            sage: f = x^2/(x+1)^3
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

    def partial_fraction_decomposition(self, var=None):
        r"""
        Return the partial fraction decomposition of ``self`` with
        respect to the given variable.

        INPUT:

        - ``var`` -- variable name or string (default: first variable)

        OUTPUT: list of symbolic expressions

        .. SEEALSO:: :meth:`partial_fraction`

        EXAMPLES::

            sage: f = x^2/(x+1)^3
            sage: f.partial_fraction_decomposition()
            [1/(x + 1), -2/(x + 1)^2, (x + 1)^(-3)]
            sage: (4+f).partial_fraction_decomposition()
            [1/(x + 1), -2/(x + 1)^2, (x + 1)^(-3), 4]

        Notice that the first variable in the expression is used by
        default::

            sage: y = var('y')
            sage: f = y^2/(y+1)^3
            sage: f.partial_fraction_decomposition()
            [1/(y + 1), -2/(y + 1)^2, (y + 1)^(-3)]

            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f.partial_fraction_decomposition()
            [y^2/(y^3 + 3*y^2 + 3*y + 1), (x - 1)^(-2), (x - 1)^(-3)]

        You can explicitly specify which variable is used::

            sage: f.partial_fraction_decomposition(y)
            [1/(y + 1), -2/(y + 1)^2, (y + 1)^(-3), x/(x^3 - 3*x^2 + 3*x - 1)]
        """
        if var is None:
            var = self.default_variable()
        return [self.parent()(ex)
                for ex in self._maxima_().partfrac(var).args()]

    def maxima_methods(self):
        """
        Provide easy access to maxima methods, converting the result to a
        Sage expression automatically.

        EXAMPLES::

            sage: t = log(sqrt(2) - 1) + log(sqrt(2) + 1); t
            log(sqrt(2) + 1) + log(sqrt(2) - 1)
            sage: res = t.maxima_methods().logcontract(); res
            log((sqrt(2) + 1)*(sqrt(2) - 1))
            sage: type(res)
            <class 'sage.symbolic.expression.Expression'>
        """
        from sage.symbolic.maxima_wrapper import MaximaWrapper
        return MaximaWrapper(self)

    def rectform(self):
        r"""
        Convert this symbolic expression to rectangular form; that
        is, the form `a + bi` where `a` and `b` are real numbers and
        `i` is the imaginary unit.

        .. NOTE::

           The name \"rectangular\" comes from the fact that, in the
           complex plane, `a` and `bi` are perpendicular.

        INPUT:

        - ``self`` -- the expression to convert

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
            sage: f = SR(a + b*I)
            sage: abs(f.rectform() - (a + b*I))  # abs tol 1e-16
            0.0

        If we decompose a complex number into its real and imaginary
        parts, they should correspond to the real and imaginary terms
        of the rectangular form::

            sage: z = CC.random_element()
            sage: a = z.real_part()
            sage: b = z.imag_part()
            sage: abs(SR(z).rectform() - (a + b*I))  # abs tol 1e-16
            0.0
        """
        return self.maxima_methods().rectform()

    def unhold(self, exclude=None):
        """
        Evaluates any held operations (with the ``hold`` keyword) in the
        expression

        INPUT:

        - ``self`` -- an expression with held operations
        - ``exclude`` -- (default: ``None``) a list of operators to exclude from
          evaluation. Excluding arithmetic operators does not yet work (see
          :issue:`10169`).

        OUTPUT:

        A new expression with held operations, except those in ``exclude``,
        evaluated

        EXAMPLES::

            sage: a = exp(I * pi, hold=True)
            sage: a
            e^(I*pi)
            sage: a.unhold()
            -1
            sage: b = x.add(x, hold=True)
            sage: b
            x + x
            sage: b.unhold()
            2*x
            sage: (a + b).unhold()
            2*x - 1
            sage: c = (x.mul(x, hold=True)).add(x.mul(x, hold=True), hold=True)
            sage: c
            x*x + x*x
            sage: c.unhold()
            2*x^2
            sage: sin(tan(0, hold=True), hold=True).unhold()
            0
            sage: sin(tan(0, hold=True), hold=True).unhold(exclude=[sin])
            sin(0)
            sage: (e^sgn(0, hold=True)).unhold()
            1
            sage: (e^sgn(0, hold=True)).unhold(exclude=[exp])
            e^0
            sage: log(3).unhold()
            log(3)
        """
        if self.operator():
            from sage.symbolic.expression_conversions import HoldRemover
            h = HoldRemover(self, exclude)
            return h()
        else:
            return self

    def simplify(self, algorithm='maxima', **kwds):
        """
        Return a simplified version of this symbolic expression.

        INPUT:

        - ``algorithm`` -- one of :

            - ``maxima`` : (default) sends the expression to
              ``maxima`` and converts it back to Sage

            - ``sympy`` : converts the expression to ``sympy``,
              simplifies it (passing any optional keyword(s)), and
              converts the result to Sage

            - ``giac`` : converts the expression to ``giac``,
              simplifies it, and converts the result to Sage

            - ``fricas`` : converts the expression to ``fricas``,
              simplifies it, and converts the result to Sage

        .. SEEALSO::

           :meth:`simplify_full`, :meth:`simplify_trig`,
           :meth:`simplify_rational`, :meth:`simplify_rectform`
           :meth:`simplify_factorial`, :meth:`simplify_log`,
           :meth:`simplify_real`, :meth:`simplify_hypergeometric`,
           :meth:`canonicalize_radical`

        EXAMPLES::

            sage: a = var('a'); f = x*sin(2)/(x^a); f
            x*sin(2)/x^a
            sage: f.simplify()
            x^(-a + 1)*sin(2)

        Some simplifications are quite algorithm-specific::

            sage: x, t = var("x, t")
            sage: ex = cos(t).exponentialize()
            sage: ex = ex.subs((sin(t).exponentialize()==x).solve(t)[0])
            sage: ex
            1/2*I*x + 1/2*I*sqrt(x^2 - 1) + 1/2/(I*x + I*sqrt(x^2 - 1))
            sage: ex.simplify()
            1/2*I*x + 1/2*I*sqrt(x^2 - 1) + 1/(2*I*x + 2*I*sqrt(x^2 - 1))
            sage: ex.simplify(algorithm='sympy')
            I*(x^2 + sqrt(x^2 - 1)*x - 1)/(x + sqrt(x^2 - 1))
            sage: ex.simplify(algorithm='giac')  # needs giac
            I*sqrt(x^2 - 1)
            sage: ex.simplify(algorithm='fricas')  # optional - fricas
            (I*x^2 + I*sqrt(x^2 - 1)*x - I)/(x + sqrt(x^2 - 1))

        TESTS:

        Check that :issue:`14637` is fixed::

            sage: assume(x > 0, x < pi/2)
            sage: acos(cos(x)).simplify()
            x
            sage: forget()

        Check that simplifying with sympy works correctly::

            sage: expr = (-1/5*(2*sqrt(6)*(sqrt(5) - 5) + 11*sqrt(5) - 11)/(2*sqrt(6)*sqrt(5) - 11))
            sage: expr.simplify(algorithm='sympy')
            1/5*sqrt(5) - 1/5
        """
        if algorithm == "maxima":
            return self._parent(self._maxima_())
        if algorithm == "sympy":
            return self._sympy_().simplify(**kwds)._sage_()
        if algorithm == "giac":
            return self._giac_().simplify()._sage_()
        if algorithm == "fricas":
            return self._fricas_().simplify()._sage_()
        raise ValueError(f"algorithm {algorithm} unknown to simplify")

    def simplify_full(self):
        """
        Apply :meth:`simplify_factorial`, :meth:`simplify_rectform`,
        :meth:`simplify_trig`, :meth:`simplify_rational`, and
        then :meth:`expand_sum` to ``self`` (in that order).

        ALIAS: ``simplify_full`` and ``full_simplify`` are the same.

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

        There are two square roots of `(x + 1)^2`, so this should
        not be simplified to `x + 1`, see :issue:`12737`::

            sage: f = sqrt((x + 1)^2)
            sage: f.simplify_full()
            sqrt(x^2 + 2*x + 1)

        The imaginary part of an expression should not change under
        simplification; :issue:`11934`::

            sage: f = sqrt(-8*(4*sqrt(2) - 7)*x^4 + 16*(3*sqrt(2) - 5)*x^3)
            sage: original = f.imag_part()
            sage: simplified = f.full_simplify().imag_part()
            sage: original - simplified
            0

        The invalid simplification from :issue:`12322` should not occur
        after :issue:`12737`::

            sage: t = var('t')
            sage: assume(t, 'complex')
            sage: assumptions()
            [t is complex]
            sage: f = (1/2)*log(2*t) + (1/2)*log(1/t)
            sage: f.simplify_full()
            1/2*log(2*t) - 1/2*log(t)
            sage: forget()

        Complex logs are not contracted, :issue:`17556`::

            sage: x,y = SR.var('x,y')
            sage: assume(y, 'complex')
            sage: f = log(x*y) - (log(x) + log(y))
            sage: f.simplify_full()
            log(x*y) - log(x) - log(y)
            sage: forget()

        The simplifications from :meth:`simplify_rectform` are
        performed, :issue:`17556`::

            sage: f = ( e^(I*x) - e^(-I*x) ) / ( I*e^(I*x) + I*e^(-I*x) )
            sage: f.simplify_full()
            sin(x)/cos(x)

        Check that :issue:`20846` is fixed::

            sage: ((1/6*pi^2).series(x)).simplify_full()
            1/6*pi^2
        """
        x = self
        x = x.simplify_factorial()
        x = x.simplify_rectform()
        x = x.simplify_trig()
        x = x.simplify_rational()
        x = x.expand_sum()
        return x

    full_simplify = simplify_full

    def simplify_hypergeometric(self, algorithm='maxima'):
        """
        Simplify an expression containing hypergeometric or confluent
        hypergeometric functions.

        INPUT:

        - ``algorithm`` -- (default: ``'maxima'``) the algorithm to use for
          for simplification. Implemented are ``'maxima'``, which uses Maxima's
          ``hgfred`` function, and ``'sage'``, which uses an algorithm
          implemented in the hypergeometric module

        ALIAS: :meth:`hypergeometric_simplify` and
        :meth:`simplify_hypergeometric` are the same

        EXAMPLES::

            sage: hypergeometric((5, 4), (4, 1, 2, 3),
            ....:                x).simplify_hypergeometric()
            1/144*x^2*hypergeometric((), (3, 4), x) +...
            1/3*x*hypergeometric((), (2, 3), x) + hypergeometric((), (1, 2), x)
            sage: (2*hypergeometric((), (), x)).simplify_hypergeometric()
            2*e^x
            sage: (nest(lambda y: hypergeometric([y], [1], x), 3, 1)  # not tested, unstable
            ....:  .simplify_hypergeometric())
            laguerre(-laguerre(-e^x, x), x)
            sage: (nest(lambda y: hypergeometric([y], [1], x), 3, 1)  # not tested, unstable
            ....:  .simplify_hypergeometric(algorithm='sage'))
            hypergeometric((hypergeometric((e^x,), (1,), x),), (1,), x)
            sage: hypergeometric_M(1, 3, x).simplify_hypergeometric()
            -2*((x + 1)*e^(-x) - 1)*e^x/x^2
            sage: (2 * hypergeometric_U(1, 3, x)).simplify_hypergeometric()
            2*(x + 1)/x^2
        """
        from sage.functions.hypergeometric import (hypergeometric,
                                                   hypergeometric_M,
                                                   hypergeometric_U,
                                                   closed_form)
        from sage.calculus.calculus import maxima
        try:
            op = self.operator()
        except RuntimeError:
            return self
        ops = self.operands()

        if op == hypergeometric_M or op == hypergeometric_U:
            return self.generalized().simplify_hypergeometric(algorithm)

        if algorithm not in ('maxima', 'sage'):
            raise NotImplementedError(
                    "unknown algorithm: '{}'".format(algorithm))

        def simplify(o):
            return o.simplify_hypergeometric(algorithm)

        if op == hypergeometric:
            a = [simplify(o) for o in ops[0].operands()]
            b = [simplify(o) for o in ops[1].operands()]
            t = simplify(ops[2])

            if algorithm == 'maxima':
                R = self.parent()
                return R(maxima.hgfred(a, b, t))
            elif algorithm == 'sage':
                return closed_form(hypergeometric(a, b, t))

        if not op:
            return self

        return op(*(simplify(o) for o in ops))

    hypergeometric_simplify = simplify_hypergeometric

    def simplify_rectform(self, complexity_measure=string_length):
        r"""
        Attempt to simplify this expression by expressing it in the
        form `a + bi` where both `a` and `b` are real. This
        transformation is generally not a simplification, so we use
        the given ``complexity_measure`` to discard
        non-simplifications.

        INPUT:

        - ``self`` -- the expression to simplify

        - ``complexity_measure`` -- (default:
          ``sage.symbolic.complexity_measures.string_length``) a
          function taking a symbolic expression as an argument and
          returning a measure of that expressions complexity. If
          ``None`` is supplied, the simplification will be performed
          regardless of the result.

        OUTPUT:

        If the transformation produces a simpler expression (according
        to ``complexity_measure``) then that simpler expression is
        returned. Otherwise, the original expression is returned.

        ALGORITHM:

        We first call :meth:`rectform()` on the given
        expression. Then, the supplied complexity measure is used to
        determine whether or not the result is simpler than the
        original expression.

        EXAMPLES:

        The exponential form of `\tan(x)`::

            sage: f = ( e^(I*x) - e^(-I*x) ) / ( I*e^(I*x) + I*e^(-I*x) )
            sage: f.simplify_rectform()
            sin(x)/cos(x)

        This should not be expanded with Euler's formula since the
        resulting expression is longer when considered as a string,
        and the default ``complexity_measure`` uses string length to
        determine which expression is simpler::

            sage: f = e^(I*x)
            sage: f.simplify_rectform()
            e^(I*x)

        However, if we pass ``None`` as our complexity measure, it
        is::

            sage: f = e^(I*x)
            sage: f.simplify_rectform(complexity_measure = None)
            cos(x) + I*sin(x)

        TESTS:

        When given ``None``, we should always call :meth:`rectform()`
        and return the result::

            sage: polynomials = QQ['x']
            sage: f = SR(polynomials.random_element())
            sage: g = f.simplify_rectform(complexity_measure = None)
            sage: bool(g == f.rectform())
            True
        """
        simplified_expr = self.rectform()

        if complexity_measure is None:
            return simplified_expr

        if complexity_measure(simplified_expr) < complexity_measure(self):
            return simplified_expr
        else:
            return self

    def simplify_real(self):
        r"""
        Simplify the given expression over the real numbers. This allows
        the simplification of `\sqrt{x^{2}}` into `\left|x\right|` and
        the contraction of `\log(x) + \log(y)` into `\log(xy)`.

        INPUT:

        - ``self`` -- the expression to convert

        OUTPUT:

        A new expression, equivalent to the original one under the
        assumption that the variables involved are real.

        EXAMPLES::

            sage: f = sqrt(x^2)
            sage: f.simplify_real()
            abs(x)

        ::

            sage: y = SR.var('y')
            sage: f = log(x) + 2*log(y)
            sage: f.simplify_real()
            log(x*y^2)

        TESTS:

        We set the Maxima ``domain`` variable to 'real' before we call
        out to Maxima. When we return, however, we should set the
        ``domain`` back to what it was, rather than assuming that it
        was 'complex'::

            sage: from sage.calculus.calculus import maxima
            sage: maxima('domain: real;')
            real
            sage: x.simplify_real()
            x
            sage: maxima('domain;')
            real
            sage: maxima('domain: complex;')
            complex

        We forget the assumptions that our variables are real after
        simplification; make sure we don't forget an assumption that
        existed before we were called::

            sage: assume(x, 'real')
            sage: x.simplify_real()
            x
            sage: assumptions()
            [x is real]
            sage: forget()

        We also want to be sure that we don't forget assumptions on
        other variables::

            sage: x,y,z = SR.var('x,y,z')
            sage: assume(y, 'integer')
            sage: assume(z, 'antisymmetric')
            sage: x.simplify_real()
            x
            sage: assumptions()
            [y is integer, z is antisymmetric]
            sage: forget()

        No new assumptions should exist after the call::

            sage: assumptions()
            []
            sage: x.simplify_real()
            x
            sage: assumptions()
            []
        """
        from sage.symbolic.assumptions import assume, assumptions, forget
        from sage.calculus.calculus import maxima
        original_domain = maxima.eval('domain')
        original_assumptions = assumptions()

        maxima.eval('domain: real$')

        # We might as well go all the way and tell Maxima to assume
        # that all variables are real. Since we're setting the
        # simplification domain (and it's indiscriminate), you'd
        # better not call this unless your variables really are real
        # anyway.
        for v in self.variables():
            assume(v, 'real')

        # This will round trip through Maxima, essentially performing
        # self.simplify() in the process.
        result = self.simplify_log()

        # Set the domain back to what it was before we were called.
        maxima.eval('domain: %s$' % original_domain)

        # Forget all assumptions, and restore the ones that existed
        # when we were called. This is much simpler than the bookkeeping
        # necessary otherwise.
        forget()
        for assumption in original_assumptions:
            assume(assumption)

        return result

    def simplify_trig(self, expand=True):
        r"""
        Optionally expand and then employ identities such as
        `\sin(x)^2 + \cos(x)^2 = 1`, `\cosh(x)^2 - \sinh(x)^2 = 1`,
        `\sin(x)\csc(x) = 1`, or `\tanh(x)=\sinh(x)/\cosh(x)`
        to simplify expressions containing tan, sec, etc., to sin,
        cos, sinh, cosh.

        INPUT:

        - ``self`` -- symbolic expression

        - ``expand`` -- boolean (default: ``True``); if ``True``, expands
          trigonometric and hyperbolic functions of sums of angles and of
          multiple angles occurring in ``self`` first. For best results,
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

            sage: f = tan(3*x)
            sage: f.simplify_trig()
            -(4*cos(x)^2 - 1)*sin(x)/(4*cos(x)*sin(x)^2 - cos(x))
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

    def simplify_rational(self, algorithm='full', map=False):
        r"""
        Simplify rational expressions.

        INPUT:

        - ``self`` -- symbolic expression

        - ``algorithm`` -- (default: ``'full'``) string which switches the
          algorithm for simplifications. Possible values are

          - ``'simple'`` (simplify rational functions into quotient of two
            polynomials)

          - ``'full'`` (apply repeatedly, if necessary)

          - ``'noexpand'`` (convert to common denominator and add)

        - ``map`` -- boolean (default: ``False``); if ``True``, the result is
          an expression whose leading operator is the same as that of the
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
        and thus the results are shorter for functions which are
        combination of rational and nonrational functions. In the
        following example, we use this option if we want not to
        combine logarithm and the rational function into one
        fraction::

            sage: f = (x^2-1)/(x+1)-ln(x)/(x+2)
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

            sage: f = 1/(x+1)+x/(x+2)^2
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
            raise NotImplementedError("unknown algorithm, see the help for available algorithms")
        P = self_m.parent()
        self_str = self_m.str()
        if map:
            cmd = "if atom(%s) then %s(%s) else map(%s,%s)" % (self_str, maxima_method, self_str, maxima_method, self_str)
        else:
            cmd = "%s(%s)" % (maxima_method, self_m.str())
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

        Check that the problem with applying ``full_simplify()`` to gamma
        functions (:issue:`9240`) has been fixed::

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

    def to_gamma(self):
        """
        Convert factorial, binomial, and Pochhammer symbol
        expressions to their gamma function equivalents.

        EXAMPLES::

            sage: m,n = var('m n', domain='integer')
            sage: factorial(n).to_gamma()
            gamma(n + 1)
            sage: binomial(m,n).to_gamma()
            gamma(m + 1)/(gamma(m - n + 1)*gamma(n + 1))
        """
        cdef GEx x
        sig_on()
        try:
            x = to_gamma(self._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def gamma_normalize(self):
        """
        Return the expression with any gamma functions that have
        a common base converted to that base.

        Additionally the expression is normalized so any fractions
        can be simplified through cancellation.

        EXAMPLES::

            sage: m,n = var('m n', domain='integer')
            sage: (gamma(n+2)/gamma(n)).gamma_normalize()
            (n + 1)*n
            sage: (gamma(n+2)*gamma(n)).gamma_normalize()
            (n + 1)*n*gamma(n)^2
            sage: (gamma(n+2)*gamma(m-1)/gamma(n)/gamma(m+1)).gamma_normalize()
            (n + 1)*n/((m - 1)*m)

        Check that :issue:`22826` is fixed::

            sage: _ = var('n')
            sage: (n-1).gcd(n+1)
            1
            sage: ex = (n-1)^2*gamma(2*n+5)/gamma(n+3) + gamma(2*n+3)/gamma(n+1)
            sage: ex.gamma_normalize()
            (4*n^3 - 2*n^2 - 7*n + 7)*gamma(2*n + 3)/((n + 1)*gamma(n + 1))
        """
        cdef GEx x
        sig_on()
        try:
            x = gamma_normalize(self._gobj)
        finally:
            sig_off()
        return new_Expression_from_GEx(self._parent, x)

    def expand_sum(self):
        r"""
        For every symbolic sum in the given expression, try to expand it,
        symbolically or numerically.

        While symbolic sum expressions with constant limits are evaluated
        immediately on the command line, unevaluated sums of this kind can
        result from, e.g., substitution of limit variables.

        INPUT:

        - ``self`` -- symbolic expression

        EXAMPLES::

            sage: (k,n) = var('k,n')
            sage: ex = sum(abs(-k*k+n),k,1,n)(n=8); ex
            sum(abs(-k^2 + 8), k, 1, 8)
            sage: ex.expand_sum()
            162
            sage: f(x,k) = sum((2/n)*(sin(n*x)*(-1)^(n+1)), n, 1, k)
            sage: f(x,2)
            -2*sum((-1)^n*sin(n*x)/n, n, 1, 2)
            sage: f(x,2).expand_sum()
            -sin(2*x) + 2*sin(x)

        We can use this to do floating-point approximation as well::

            sage: (k,n) = var('k,n')
            sage: f(n)=sum(sqrt(abs(-k*k+n)),k,1,n)
            sage: f(n=8)
            sum(sqrt(abs(-k^2 + 8)), k, 1, 8)
            sage: f(8).expand_sum()
            sqrt(41) + sqrt(17) + 2*sqrt(14) + 3*sqrt(7) + 2*sqrt(2) + 3
            sage: f(8).expand_sum().n()
            31.7752256945384

        See :issue:`9424` for making the following no longer raise
        an error::

            sage: f(8).n()
            31.7752256945384
        """
        return self.parent()(self._maxima_().simplify_sum())

    def canonicalize_radical(self):
        r"""
        Choose a canonical branch of the given expression.

        The square root, cube root, natural log, etc. functions are
        multi-valued. The ``canonicalize_radical()`` method will
        choose *one* of these values based on a heuristic.

        For example, ``sqrt(x^2)`` has two values: ``x``, and
        ``-x``. The ``canonicalize_radical()`` function will choose
        *one* of them, consistently, based on the behavior of the
        expression as ``x`` tends to positive infinity. The solution
        chosen is the one which exhibits this same behavior. Since
        ``sqrt(x^2)`` approaches positive infinity as ``x`` does, the
        solution chosen is ``x`` (which also tends to positive
        infinity).

        .. WARNING::

            As shown in the examples below, a canonical form is not always
            returned, i.e., two mathematically identical expressions might
            be converted to different expressions.

            Assumptions are not taken into account during the
            transformation. This may result in a branch choice
            inconsistent with your assumptions.

        ALGORITHM:

        This uses the Maxima ``radcan()`` command. From the Maxima
        documentation:

        .. pull-quote::

            Simplifies an expression, which can contain logs,
            exponentials, and radicals, by converting it into a form
            which is canonical over a large class of expressions and a
            given ordering of variables; that is, all functionally
            equivalent forms are mapped into a unique form. For a
            somewhat larger class of expressions, radcan produces a
            regular form. Two equivalent expressions in this class do
            not necessarily have the same appearance, but their
            difference can be simplified by radcan to zero.

            For some expressions radcan is quite time consuming. This
            is the cost of exploring certain relationships among the
            components of the expression for simplifications based on
            factoring and partial fraction expansions of exponents.

        EXAMPLES:

        ``canonicalize_radical()`` can perform some of the same
        manipulations as :meth:`log_expand`::

            sage: y = SR.symbol('y')
            sage: f = log(x*y)
            sage: f.log_expand()
            log(x) + log(y)
            sage: f.canonicalize_radical()
            log(x) + log(y)

        And also handles some exponential functions::

            sage: f = (e^x-1)/(1+e^(x/2))
            sage: f.canonicalize_radical()
            e^(1/2*x) - 1

        It can also be used to change the base of a logarithm when the
        arguments to ``log()`` are positive real numbers::

            sage: f = log(8)/log(2)
            sage: f.canonicalize_radical()
            3

        ::

            sage: a = SR.symbol('a')
            sage: f = (log(x+x^2)-log(x))^a/log(1+x)^(a/2)
            sage: f.canonicalize_radical()
            log(x + 1)^(1/2*a)

        The simplest example of counter-intuitive behavior is what
        happens when we take the square root of a square::

            sage: sqrt(x^2).canonicalize_radical()
            x

        If you don't want this kind of "simplification," don't use
        ``canonicalize_radical()``.

        This behavior can also be triggered when the expression under
        the radical is not given explicitly as a square::

            sage: sqrt(x^2 - 2*x + 1).canonicalize_radical()
            x - 1

        Another place where this can become confusing is with
        logarithms of complex numbers. Suppose ``x`` is complex with
        ``x == r*e^(I*t)`` (``r`` real). Then ``log(x)`` is
        ``log(r) + I*(t + 2*k*pi)`` for some integer ``k``.

        Calling ``canonicalize_radical()`` will choose a branch,
        eliminating the solutions for all choices of ``k`` but
        one. Simplified by hand, the expression below is
        ``(1/2)*log(2) + I*pi*k`` for integer ``k``. However,
        ``canonicalize_radical()`` will take each log expression, and
        choose one particular solution, dropping the other. When the
        results are subtracted, we're left with no imaginary part::

            sage: f = (1/2)*log(2*x) + (1/2)*log(1/x)
            sage: f.canonicalize_radical()
            1/2*log(2)

        Naturally the result is wrong for some choices of ``x``::

            sage: f(x = -1)
            I*pi + 1/2*log(2)

        The example below shows two expressions e1 and e2 which are
        "simplified" to different expressions, while their difference
        is "simplified" to zero; thus ``canonicalize_radical()`` does
        not return a canonical form::

            sage: e1 = 1/(sqrt(5)+sqrt(2))
            sage: e2 = (sqrt(5)-sqrt(2))/3
            sage: e1.canonicalize_radical()
            1/(sqrt(5) + sqrt(2))
            sage: e2.canonicalize_radical()
            1/3*sqrt(5) - 1/3*sqrt(2)
            sage: (e1-e2).canonicalize_radical()
            0

        The issue reported in :issue:`3520` is a case where
        ``canonicalize_radical()`` causes a numerical integral to be
        calculated incorrectly::

            sage: f1 = sqrt(25 - x) * sqrt( 1 + 1/(4*(25-x)) )
            sage: f2 = f1.canonicalize_radical()
            sage: numerical_integral(f1.real(), 0, 1)[0] # abs tol 1e-10
            4.974852579915647
            sage: numerical_integral(f2.real(), 0, 1)[0] # abs tol 1e-10
            -4.974852579915647

        TESTS:

        This tests that :issue:`11668` has been fixed (by :issue:`12780`)::

            sage: a,b = var('a b', domain='real')
            sage: A = abs((a+I*b))^2
            sage: imag(A)
            0
            sage: A.canonicalize_radical() # not implemented
            a^2 + b^2
            sage: imag(A.canonicalize_radical())
            0
        """
        from sage.calculus.calculus import maxima
        return self.parent()(self._maxima_().radcan())

    def simplify_log(self, algorithm=None):
        r"""
        Simplify a (real) symbolic expression that contains logarithms.

        The given expression is scanned recursively, transforming
        subexpressions of the form `a \log(b) + c \log(d)` into
        `\log(b^{a} d^{c})` before simplifying within the ``log()``.

        The user can specify conditions that `a` and `c` must satisfy
        before this transformation will be performed using the optional
        parameter ``algorithm``.

        .. WARNING::

            This is only safe to call if every variable in the given
            expression is assumed to be real. The simplification it performs
            is in general not valid over the complex numbers. For example::

                sage: x,y = SR.var('x,y')
                sage: f = log(x*y) - (log(x) + log(y))
                sage: f(x=-1, y=i)
                -2*I*pi
                sage: f.simplify_log()
                0

        INPUT:

        - ``self`` -- expression to be simplified

        - ``algorithm`` -- (default: ``None``) governs the condition
          on `a` and `c` which must be satisfied to contract expression
          `a \log(b) + c \log(d)`. Values are

          - ``None`` (use Maxima default, integers),

          - ``'one'`` (1 and -1),

          - ``'ratios'`` (rational numbers),

          - ``'constants'`` (constants),

          - ``'all'`` (all expressions).

        ALGORITHM:

        This uses the Maxima ``logcontract()`` command.

        ALIAS:

        :meth:`log_simplify` and :meth:`simplify_log` are the same.

        EXAMPLES::

            sage: x,y,t = var('x y t')

        Only two first terms are contracted in the following example;
        the logarithm with coefficient `\frac{1}{2}` is not contracted::

            sage: f = log(x)+2*log(y)+1/2*log(t)
            sage: f.simplify_log()
            log(x*y^2) + 1/2*log(t)

        To contract all terms in the previous example, we use the
        ``'ratios'`` ``algorithm``::

            sage: f.simplify_log(algorithm='ratios')
            log(sqrt(t)*x*y^2)

        To contract terms with no coefficient (more precisely, with
        coefficients `1` and `-1`), we use the ``'one'``
        ``algorithm``::

            sage: f = log(x)+2*log(y)-log(t)
            sage: f.simplify_log('one')
            2*log(y) + log(x/t)

        ::

            sage: f = log(x)+log(y)-1/3*log((x+1))
            sage: f.simplify_log()
            log(x*y) - 1/3*log(x + 1)

            sage: f.simplify_log('ratios')
            log(x*y/(x + 1)^(1/3))

        `\pi` is an irrational number; to contract logarithms in the
        following example we have to set ``algorithm`` to ``'constants'``
        or ``'all'``::

            sage: f = log(x)+log(y)-pi*log((x+1))
            sage: f.simplify_log('constants')
            log(x*y/(x + 1)^pi)

        ``x*log(9)`` is contracted only if ``algorithm`` is ``'all'``::

            sage: (x*log(9)).simplify_log()
            2*x*log(3)
            sage: (x*log(9)).simplify_log('all')
            log(3^(2*x))

        TESTS:

        Ensure that the option ``algorithm`` from one call has no
        influence upon future calls (a Maxima flag was set, and we have
        to ensure that its value has been restored)::

            sage: f = log(x)+2*log(y)+1/2*log(t)
            sage: f.simplify_log('one')
            1/2*log(t) + log(x) + 2*log(y)

            sage: f.simplify_log('ratios')
            log(sqrt(t)*x*y^2)

            sage: f.simplify_log()
            log(x*y^2) + 1/2*log(t)

        This shows that the issue at :issue:`7334` is fixed. Maxima
        intentionally keeps the expression inside the log factored::

            sage: log_expr = (log(sqrt(2)-1)+log(sqrt(2)+1))
            sage: log_expr.simplify_log('all')
            log((sqrt(2) + 1)*(sqrt(2) - 1))
            sage: _.simplify_rational()
            0

        We should use the current simplification domain rather than
        set it to 'real' explicitly (:issue:`12780`)::

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
            raise NotImplementedError("unknown algorithm, see the help for available algorithms")
        res = self.parent()(self._maxima_().logcontract())
        if algorithm is not None:
            maxima.eval('logconcoeffp:false$')
        maxima.eval('logexpand:savelogexpand$')
        return res

    log_simplify = simplify_log

    def expand_log(self, algorithm='products'):
        r"""
        Simplify symbolic expression, which can contain logs.

        Expands logarithms of powers, logarithms of products and
        logarithms of quotients.  The option ``algorithm`` specifies
        which expression types should be expanded.

        INPUT:

        - ``self`` -- expression to be simplified

        - ``algorithm`` -- (default: ``'products'``) governs which
          expression is expanded. Possible values are

          - ``'nothing'`` (no expansion),

          - ``'powers'`` (log(a^r) is expanded),

          - ``'products'`` (like 'powers' and also log(a*b) are expanded),

          - ``'all'`` (all possible expansion).

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
            pi*log(x) + log(3) - 2*log(2)

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
            pi*log(x) + log(3) - 2*log(2)

            sage: (log(3/4*x^pi)).log_expand()
            pi*log(x) + log(3/4)

        TESTS:

        Most of these log expansions only make sense over the
        reals. So, we should set the Maxima ``domain`` variable to
        'real' before we call out to Maxima. When we return, however, we
        should set the ``domain`` back to what it was, rather than
        assuming that it was 'complex'. See :issue:`12780`::

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
            raise NotImplementedError("unknown algorithm, see the help for available algorithms")
        maxima.eval('logexpand:%s' % maxima_method)
        res = self._maxima_()
        res = res.sage()
        # Set the domain back to what it was before expand_log() was called.
        maxima.eval('domain: %s$ logexpand:savelogexpand$' % original_domain)
        return res

    log_expand = expand_log

    def distribute(self, recursive=True):
        """
        Distribute some indexed operators over similar operators in
        order to allow further groupings or simplifications.

        Implemented cases (so far):

        - Symbolic sum of a sum ==> sum of symbolic sums

        - Integral (definite or not) of a sum ==> sum of integrals.

        - Symbolic product of a product ==> product of symbolic products.

        INPUT:

        - ``recursive`` -- (default: ``True``) the distribution proceeds
          along the subtrees of the expression

        TESTS:

            sage: var("j,k,p,q", domain='integer')
            (j, k, p, q)
            sage: X,Y,Z,f,g = function("X,Y,Z,f,g")
            sage: var("x,a,b")
            (x, a, b)
            sage: sum(X(j)+Y(j),j,1,p)
            sum(X(j) + Y(j), j, 1, p)
            sage: sum(X(j)+Y(j),j,1,p).distribute()
            sum(X(j), j, 1, p) + sum(Y(j), j, 1, p)
            sage: integrate(f(x)+g(x),x)
            integrate(f(x) + g(x), x)
            sage: integrate(f(x)+g(x),x).distribute()
            integrate(f(x), x) + integrate(g(x), x)
            sage: result = integrate(f(x)+g(x),x,a,b)
            ...
            sage: result
            integrate(f(x) + g(x), x, a, b)
            sage: result = integrate(f(x)+g(x),x,a,b).distribute()
            ...
            sage: result
            integrate(f(x), x, a, b) + integrate(g(x), x, a, b)
            sage: sum(X(j)+sum(Y(k)+Z(k),k,1,q),j,1,p)
            sum(X(j) + sum(Y(k) + Z(k), k, 1, q), j, 1, p)
            sage: sum(X(j)+sum(Y(k)+Z(k),k,1,q),j,1,p).distribute()
            sum(sum(Y(k), k, 1, q) + sum(Z(k), k, 1, q), j, 1, p) + sum(X(j), j, 1, p)
            sage: sum(X(j)+sum(Y(k)+Z(k),k,1,q),j,1,p).distribute(recursive=False)
            sum(X(j), j, 1, p) + sum(sum(Y(k) + Z(k), k, 1, q), j, 1, p)
            sage: maxima("product(X(j)*Y(j),j,1,p)").sage()
            product(X(j)*Y(j), j, 1, p)
            sage: maxima("product(X(j)*Y(j),j,1,p)").sage().distribute()
            product(X(j), j, 1, p)*product(Y(j), j, 1, p)


        AUTHORS:

        - Emmanuel Charpentier, Ralf Stephan (05-2017)
        """
        from sage.functions.other import symbolic_sum as opsum, \
            symbolic_product as opprod
        from sage.symbolic.integration.integral \
            import indefinite_integral as opii, definite_integral as opdi
        from sage.symbolic.operators import add_vararg as opadd, \
            mul_vararg as opmul
        from sage.misc.misc_c import prod
        from sage.symbolic.ring import SR

        def treat_term(op, term, args):
            l = copy(args)
            l.insert(0, term)
            return op(*l)

        if self.parent() is not SR:
            return self

        op = self.operator()
        if op is None:
            return self

        if op in {opsum, opdi, opii}:
            sa = self.operands()[0].expand()
            op1 = sa.operator()
            if op1 is opadd:
                la = self.operands()[1:]
                aa = sa.operands()
                if recursive:
                    return sum(treat_term(op, t.distribute(), la) for t in aa)
                return sum(treat_term(op, t, la) for t in aa)
            return self
        if op is opprod:
            sa = self.operands()[0].expand()
            op1 = sa.operator()
            if op1 is opmul:
                la = self.operands()[1:]
                aa = sa.operands()
                if recursive:
                    return prod(treat_term(op, t.distribute(), la) for t in aa)
                return prod(treat_term(op, t, la) for t in aa)
            return self

        if recursive:
            done = [t.distribute() for t in self.operands()]
            return op(*done)
        return self

    def factor(self, dontfactor=None):
        """
        Factor the expression, containing any number of variables or functions, into
        factors irreducible over the integers.

        INPUT:

        - ``self`` -- a symbolic expression

        - ``dontfactor`` -- list (default: ``[]``); a list of
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

        TESTS:

        Check that :issue:`21529` is fixed::

            sage: f(x) = function('f')(x)
            sage: (f(x).diff(x)^2-1).factor()
            (diff(f(x), x) + 1)*(diff(f(x), x) - 1)

        Check that :issue:`27304` is fixed::

            sage: factor(2*exp(x) + exp(-x))
            (2*e^(2*x) + 1)*e^(-x)
            sage: factor(x*exp(-x) + exp(-x))
            (x + 1)*e^(-x)
            sage: factor(x + sqrt(x))
            x + sqrt(x)
            sage: factor((x + sqrt(x))/(x - sqrt(x)))
            (x + sqrt(x))/(x - sqrt(x))

        Check that :issue:`33640` is fixed::

            sage: ((x + 1)^2 - 2*x - 1).factor()
            x^2
        """
        from sage.calculus.calculus import symbolic_expression_from_maxima_string
        cdef GEx x
        cdef bint b
        if dontfactor or not self.is_rational_expression():
            if dontfactor is None:
                dontfactor = []
            m = self._maxima_()
            name = m.name()
            varstr = ','.join('_SAGE_VAR_' + str(v) for v in dontfactor)
            cmd = 'block([dontfactor:[%s]],factor(%s))' % (varstr, name)
            return symbolic_expression_from_maxima_string(cmd)
        sig_on()
        try:
            b = g_factor(self._gobj, x)
        finally:
            sig_off()
        if b:
            return new_Expression_from_GEx(self._parent, x)
        else:
            return self

    def factor_list(self, dontfactor=None):
        """
        Return a list of the factors of self, as computed by the
        factor command.

        INPUT:

        - ``self`` -- a symbolic expression

        - ``dontfactor`` -- see docs for :meth:`factor`

        .. NOTE::

           If you already have a factored expression and just want to
           get at the individual factors, use the ``_factor_list`` method
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
            <... 'list'>
        """
        op = self.operator()
        from sage.symbolic.operators import mul_vararg
        if op is mul_vararg:
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
        Call the convert function in the units package. For symbolic
        variables that are not units, this function just returns the
        variable.

        INPUT:

        - ``self`` -- the symbolic expression converting from
        - ``target`` -- (default: ``None``) the symbolic expression
          converting to

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: units.length.foot.convert()
            381/1250*meter
            sage: units.mass.kilogram.convert(units.mass.pound)
            100000000/45359237*pound

        We do not get anything new by converting an ordinary symbolic variable::

            sage: a = var('a')
            sage: a - a.convert()
            0

        Raises :exc:`ValueError` if ``self`` and ``target`` are not convertible::

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
            ValueError: cannot convert
        """
        from sage.symbolic import units
        return units.convert(self, target)

    ###################################################################
    # solve
    ###################################################################
    def roots(self, x=None, explicit_solutions=True, multiplicities=True, ring=None):
        r"""
        Return roots of ``self`` that can be found exactly,
        possibly with multiplicities.  Not all roots are guaranteed to
        be found.

        .. warning::

           This is *not* a numerical solver - use :meth:`find_root` to
           solve for ``self == 0`` numerically on an interval.

        INPUT:

        - ``x`` -- variable to view the function in terms of
          (use default variable if not given)

        - ``explicit_solutions`` -- boolean (default: ``True``); require that
          roots be explicit rather than implicit

        - ``multiplicities`` -- boolean (default: ``True``); when ``True``, return
          multiplicities

        - ``ring`` -- a ring (default: ``None``); if not ``None``, convert
          ``self`` to a polynomial over ``ring`` and find roots over ``ring``

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

        .. NOTE::

            It is possible to solve a greater variety of equations
            using :func:`solve` and the keyword ``to_poly_solve``,
            but only at the price of possibly encountering
            approximate solutions.  See documentation for :meth:`solve`
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

        Another example, but involving a degree 5 poly whose roots do not
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

        Now let us find some roots over different rings::

            sage: f.roots(ring=CC)
            [(-0.0588115223184..., 1),
             (-1.331099917875... - 1.52241655183732*I, 1),
             (-1.331099917875... + 1.52241655183732*I, 1),
             (1.36050567903502 - 1.51880872209965*I, 1),
             (1.36050567903502 + 1.51880872209965*I, 1)]
            sage: (2.5*f).roots(ring=RR)
            [(-0.058811522318449..., 1)]
            sage: f.roots(ring=CC, multiplicities=False)
            [-0.05881152231844...,
             -1.331099917875... - 1.52241655183732*I,
             -1.331099917875... + 1.52241655183732*I,
             1.36050567903502 - 1.51880872209965*I,
             1.36050567903502 + 1.51880872209965*I]
            sage: f.roots(ring=QQ)
            []
            sage: f.roots(ring=QQbar, multiplicities=False)
            [-0.05881152231844944?,
             -1.331099917875796? - 1.522416551837318?*I,
             -1.331099917875796? + 1.522416551837318?*I,
             1.360505679035020? - 1.518808722099650?*I,
             1.360505679035020? + 1.518808722099650?*I]

        Root finding over finite fields::

            sage: f.roots(ring=GF(7^2, 'a'))                                            # needs sage.rings.finite_rings
            [(3, 1), (4*a + 6, 2), (3*a + 3, 2)]

        TESTS::

            sage: (sqrt(3) * f).roots(ring=QQ)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(3) to a rational

        Check if :issue:`9538` is fixed::

            sage: var('f6,f5,f4,x')
            (f6, f5, f4, x)
            sage: e = 15*f6*x^2 + 5*f5*x + f4
            sage: res = e.roots(x); res
            [(-1/30*(5*f5 + sqrt(25*f5^2 - 60*f4*f6))/f6, 1),
             (-1/30*(5*f5 - sqrt(25*f5^2 - 60*f4*f6))/f6, 1)]
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
            raise RuntimeError("no explicit roots found")
        else:
            rt_muls = [(S[i].rhs(), mul[i]) for i in range(len(mul))]
        if multiplicities:
            return rt_muls
        else:
            return [rt for rt, mul in rt_muls]

    def solve(self, x, multiplicities=False, solution_dict=False, explicit_solutions=False, to_poly_solve=False, algorithm=None, domain=None):
        r"""
        Analytically solve the equation ``self == 0`` or a univariate
        inequality for the variable `x`.

        .. warning::

           This is not a numerical solver -- use :meth:`find_root` to solve
           for ``self == 0`` numerically on an interval.

        INPUT:

        - ``x`` -- variable(s) to solve for

        - ``multiplicities`` -- boolean (default: ``False``); if ``True``,
          return corresponding multiplicities.  This keyword is
          incompatible with ``to_poly_solve=True`` and does not make
          any sense when solving an inequality.

        - ``solution_dict`` -- boolean (default: ``False``); if ``True`` or nonzero,
          return a list of dictionaries containing solutions. Not used
          when solving an inequality.

        - ``explicit_solutions`` -- boolean (default: ``False``); require that
          all roots be explicit rather than implicit. Not used
          when solving an inequality.

        - ``to_poly_solve`` -- boolean (default: ``False``) or string; use
          Maxima's ``to_poly_solver`` package to search for more possible
          solutions, but possibly encounter approximate solutions.
          This keyword is incompatible with ``multiplicities=True``
          and is not used when solving an inequality. Setting ``to_poly_solve``
          to ``'force'`` omits Maxima's solve command (useful when
          some solutions of trigonometric equations are lost).

        EXAMPLES::

            sage: z = var('z')
            sage: (z^5 - 1).solve(z)
            [z == 1/4*sqrt(5) + 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4,
             z == -1/4*sqrt(5) + 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4,
             z == -1/4*sqrt(5) - 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4,
             z == 1/4*sqrt(5) - 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4,
             z == 1]

            sage: solve((z^3-1)^3, z, multiplicities=True)
            ([z == 1/2*I*sqrt(3) - 1/2, z == -1/2*I*sqrt(3) - 1/2, z == 1], [3, 3, 3])

        TESTS:

        Check that :issue:`20755` is indeed fixed::

            sage: w = x^4 - (1+3*i)*x^3 - (2-4*i)*x^2 + (6-2*i)*x - 4 - 4*i
            sage: w.solve(x,multiplicities=True)
            ([x == -1/2*sqrt(2*I) + 3/2*I - 1/2,
              x == 1/2*sqrt(2*I) + 3/2*I - 1/2,
              x == (-I + 1), x == (I + 1)],
             [1, 1, 1, 1])

        See :func:`sage.symbolic.relation.solve` or the output of ``solve?``
        for extensive documentation.
        """
        from sage.symbolic.relation import solve
        return solve(self, x, multiplicities=multiplicities,
                     solution_dict=solution_dict,
                     explicit_solutions=explicit_solutions,
                     to_poly_solve=to_poly_solve,
                     algorithm=algorithm,
                     domain=domain)

    def solve_diophantine(self, x=None, solution_dict=False):
        """
        Solve a polynomial equation in the integers (a so called Diophantine).

        If the argument is just a polynomial expression, equate to zero.
        If ``solution_dict=True``, return a list of dictionaries instead of
        a list of tuples.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: solve_diophantine(3*x == 4)                                           # needs sympy
            []
            sage: solve_diophantine(x^2 - 9)                                            # needs sympy
            [-3, 3]
            sage: sorted(solve_diophantine(x^2 + y^2 == 25))                            # needs sympy
            [(-5, 0), (-4, -3), (-4, 3), (-3, -4), (-3, 4), (0, -5)...

        The function is used when ``solve()`` is called with all variables
        assumed integer::

            sage: assume(x, 'integer')
            sage: assume(y, 'integer')
            sage: sorted(solve(x*y == 1, (x,y)))                                        # needs sympy
            [(-1, -1), (1, 1)]

        You can also pick specific variables, and get the solution as
        a dictionary::

            sage: # needs sympy
            sage: solve_diophantine(x*y == 10, x)
            [-10, -5, -2, -1, 1, 2, 5, 10]
            sage: sorted(solve_diophantine(x*y - y == 10, (x,y)))
            [(-9, -1), (-4, -2), (-1, -5), (0, -10), (2, 10), (3, 5), (6, 2), (11, 1)]
            sage: res = solve_diophantine(x*y - y == 10, solution_dict=True)
            sage: sol = [{y: -5, x: -1}, {y: -10, x: 0}, {y: -1, x: -9}, {y: -2, x: -4},
            ....:        {y: 10, x: 2}, {y: 1, x: 11}, {y: 2, x: 6}, {y: 5, x: 3}]
            sage: all(solution in res
            ....:     for solution in sol) and bool(len(res) == len(sol))
            True

        If the solution is parametrized the parameter(s) are not defined,
        but you can substitute them with specific integer values::

            sage: # needs sympy
            sage: x,y,z = var('x,y,z')
            sage: sol = solve_diophantine(x^2-y == 0); sol
            (t, t^2)
            sage: [(sol[0].subs(t=t),sol[1].subs(t=t)) for t in range(-3,4)]
            [(-3, 9), (-2, 4), (-1, 1), (0, 0), (1, 1), (2, 4), (3, 9)]
            sage: sol = solve_diophantine(x^2 + y^2 == z^2); sol
            (2*p*q, p^2 - q^2, p^2 + q^2)
            sage: [(sol[0].subs(p=p,q=q), sol[1].subs(p=p,q=q), sol[2].subs(p=p,q=q))
            ....:  for p in range(1,4) for q in range(1,4)]
            [(2, 0, 2), (4, -3, 5), (6, -8, 10), (4, 3, 5), (8, 0, 8),
             (12, -5, 13), (6, 8, 10), (12, 5, 13), (18, 0, 18)]

        Solve Brahmagupta-Pell equations::

            sage: sol = sorted(solve_diophantine(x^2 - 2*y^2 == 1), key=str); sol       # needs sympy
            [(-sqrt(2)*(2*sqrt(2) + 3)^t + sqrt(2)*(-2*sqrt(2) + 3)^t
               - 3/2*(2*sqrt(2) + 3)^t - 3/2*(-2*sqrt(2) + 3)^t,...
            sage: [(sol[1][0].subs(t=t).simplify_full(),                                # needs sympy
            ....:   sol[1][1].subs(t=t).simplify_full()) for t in range(-1,5)]
            [(1, 0), (3, -2), (17, -12), (99, -70), (577, -408), (3363, -2378)]

        TESTS::

            sage: solve_diophantine(x^2 - y, x, y)                                      # needs sympy
            Traceback (most recent call last):
            ...
            AttributeError: please use a tuple or list for several variables.

        .. SEEALSO::

            http://docs.sympy.org/latest/modules/solvers/diophantine.html
        """
        from sage.symbolic.ring import SR
        from sympy.solvers.diophantine import diophantine

        if not isinstance(solution_dict, bool):
            raise AttributeError("please use a tuple or list for several variables.")
        if is_a_relational(self._gobj) and self.operator() is operator.eq:
            ex = self.lhs() - self.rhs()
        else:
            ex = self
        sympy_ex = ex._sympy_()
        solutions = diophantine(sympy_ex)
        if isinstance(solutions, set):
            solutions = list(solutions)

        if len(solutions) == 0:
            return []
        if not isinstance(solutions[0], tuple):
            solutions = [SR(sol) for sol in solutions]
        else:
            solutions = [tuple(SR(s) for s in sol) for sol in solutions]
        if x is None:
            wanted_vars = ex.variables()
            var_idx = list(range(len(ex.variables())))
        else:
            if isinstance(x, (list, tuple)):
                wanted_vars = x
            else:
                wanted_vars = [x]
            var_idx = [ex.variables().index(v) for v in wanted_vars]

        if solution_dict is False:
            if len(wanted_vars) == 1:
                ret = sorted([sol[var_idx[0]] for sol in solutions])
            else:
                ret = [tuple([sol[i] for i in var_idx]) for sol in solutions]
        else:
            ret = [dict([[ex.variables()[i],sol[i]] for i in var_idx]) for sol in solutions]

        if len(ret) == 1:
            ret = ret[0]
        return ret

    def find_root(self, a, b, var=None, xtol=10e-13, rtol=2.0**-50,
                  maxiter=100, full_output=False, imaginary_tolerance=1e-8):
        """
        Numerically find a root of ``self`` on the closed interval [a,b] (or
        [b,a]) if possible, where ``self`` is a function in one variable.
        Note: this function only works in fixed (machine) precision, it is not
        possible to get arbitrary precision approximations with it.

        INPUT:

        - ``a``, ``b`` -- endpoints of the interval

        - ``var`` -- (optional) variable

        - ``xtol, rtol`` -- the routine converges when a root is known to lie
          within ``xtol`` of the value return. Should be nonnegative. The
          routine modifies this to take into account the relative precision of
          doubles.

        - ``maxiter`` -- integer; if convergence is not achieved in maxiter
          iterations, an error is raised. Must be nonnegative.

        - ``full_output`` -- boolean (default: ``False``); if ``True``,
          also return object that contains information about convergence

        - ``imaginary_tolerance`` -- (default: ``1e-8``) if an imaginary
          number arises (due, for example, to numerical issues), this
          tolerance specifies how large it has to be in magnitude before
          we raise an error. In other words, imaginary parts smaller than
          this are ignored when we are expecting a real answer.

        EXAMPLES:

        Note that in this example both f(-2) and f(3) are positive,
        yet we still find a root in that interval::

            sage: # needs scipy
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

            sage: (sin(x) + exp(x)).find_root(-10, 10)                                  # needs scipy
            -0.588532743981862...
            sage: sin(x).find_root(-1,1)                                                # needs scipy
            0.0

        This example was fixed along with :issue:`4942` -
        there was an error in the example
        pi is a root for tan(x), but an asymptote to 1/tan(x)
        added an example to show handling of both cases::

            sage: (tan(x)).find_root(3,3.5)                                             # needs scipy
            3.1415926535...
            sage: (1/tan(x)).find_root(3, 3.5)                                          # needs scipy
            Traceback (most recent call last):
            ...
            NotImplementedError: Brent's method failed to find a zero for f on the interval

        An example with a square root::

            sage: f = 1 + x + sqrt(x+2); f.find_root(-2,10)                             # needs scipy
            -1.618033988749895

        Some examples that Ted Kosan came up with::

            sage: t = var('t')
            sage: v = 0.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))
            sage: v.find_root(0, 0.002)                                                 # needs scipy
            0.001540327067911417...

        With this expression, we can see there is a
        zero very close to the origin::

            sage: a = .004*(8*e^(-(300*t)) - 8*e^(-(1200*t)))*(720000*e^(-(300*t)) - 11520000*e^(-(1200*t))) +.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))^2
            sage: show(plot(a, 0, .002), xmin=0, xmax=.002)                             # needs sage.plot

        It is easy to approximate with ``find_root``::

            sage: a.find_root(0,0.002)                                                  # needs scipy
            0.0004110514049349...

        Using solve takes more effort, and even then gives
        only a solution with free (integer) variables::

            sage: a.solve(t)
            []
            sage: b = a.canonicalize_radical(); b
            (46080.0*e^(1800*t) - 576000.0*e^(900*t) + 737280.0)*e^(-2400*t)
            sage: b.solve(t)
            []
            sage: b.solve(t, to_poly_solve=True)
            [t == 1/450*I*pi*z... + 1/900*log(-3/4*sqrt(41) + 25/4),
             t == 1/450*I*pi*z... + 1/900*log(3/4*sqrt(41) + 25/4)]
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
        :issue:`3980`::

            sage: t = var('t')
            sage: find_root(1/t - x,0,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented
            in 1 dimension.

        Ensure that complex expressions do not cause a problem if they
        appear only as intermediate results as in :issue:`24536`::

            sage: x = SR.symbol('x', domain='real')
            sage: f = (sqrt(x) - I).abs()
            sage: f.find_root(-2, 2, rtol=1e-6)  # abs tol 1e-6                         # needs scipy
            -1.0000000049668551
        """
        if is_a_relational(self._gobj) and self.operator() is not operator.eq:
            raise ValueError("Symbolic equation must be an equality.")
        from sage.numerical.optimize import find_root
        if self.number_of_arguments() == 0:
            if self.is_trivial_zero():
                return a
            else:
                raise RuntimeError("no zero in the interval, since constant expression is not 0.")
        elif self.number_of_arguments() == 1:
            # We use the plotting fast_callable() and its wrapper
            # because numerical root finding is a lot like plotting in
            # the sense that we usually want to tolerate intermediate
            # complex expressions like abs(x + I), so long as the end
            # result has a "negligible" complex part.
            from sage.ext.fast_callable import FastCallableFloatWrapper
            f = self._plot_fast_callable(self.default_variable())
            ff = FastCallableFloatWrapper(f, imag_tol=imaginary_tolerance)
            return find_root(ff, a=a, b=b, xtol=xtol,
                             rtol=rtol,maxiter=maxiter,
                             full_output=full_output)
        else:
            raise NotImplementedError("root finding currently only implemented in 1 dimension.")

    def find_local_maximum(self, a, b, var=None, tol=1.48e-08,
                           maxfun=500, imaginary_tolerance=1e-8):
        r"""
        Numerically find a local maximum of the expression ``self``
        on the interval [a,b] (or [b,a]) along with the point at which the
        maximum is attained.

        See the documentation for
        :func:`find_local_minimum` for more details.

        EXAMPLES::

            sage: f = x*cos(x)
            sage: f.find_local_maximum(0,5)                                             # needs scipy
            (0.5610963381910451, 0.8603335890...)
            sage: f.find_local_maximum(0,5, tol=0.1, maxfun=10)                         # needs scipy
            (0.561090323458081..., 0.857926501456...)
        """
        minval, x = (-self).find_local_minimum(a, b, var, tol,
                                               maxfun, imaginary_tolerance)
        return -minval, x

    def find_local_minimum(self, a, b, var=None, tol=1.48e-08,
                           maxfun=500, imaginary_tolerance=1e-8):
        r"""
        Numerically find a local minimum of the expression ``self``
        on the interval [a,b] (or [b,a]) and the point at which it attains
        that minimum. Note that ``self`` must be a function of
        (at most) one variable.

        INPUT:

        - ``a`` -- real number; left endpoint of interval on which to minimize

        - ``b`` -- real number; right endpoint of interval on which to minimize

        - ``var`` -- variable (default: first variable in ``self``); the
          variable in ``self`` to maximize over

        - ``tol`` -- positive real (default: 1.48e-08); the convergence
          tolerance

        - ``maxfun`` -- natural number (default: 500); maximum function
          evaluations

        - ``imaginary_tolerance`` -- (default: ``1e-8``) if an imaginary
          number arises (due, for example, to numerical issues), this
          tolerance specifies how large it has to be in magnitude before
          we raise an error. In other words, imaginary parts smaller than
          this are ignored when we are expecting a real answer.

        OUTPUT:

        A tuple ``(minval, x)``, where

        - ``minval`` -- float; the minimum value that ``self`` takes on in
          the interval ``[a,b]``

        - ``x`` -- float; the point at which ``self`` takes on the minimum
          value

        EXAMPLES::

            sage: # needs scipy
            sage: f = x*cos(x)
            sage: f.find_local_minimum(1, 5)
            (-3.288371395590..., 3.4256184695...)
            sage: f.find_local_minimum(1, 5, tol=1e-3)
            (-3.288371361890..., 3.4257507903...)
            sage: f.find_local_minimum(1, 5, tol=1e-2, maxfun=10)
            (-3.288370845983..., 3.4250840220...)
            sage: show(f.plot(0, 20))                                                   # needs sage.plot
            sage: f.find_local_minimum(1, 15)
            (-9.477294259479..., 9.5293344109...)

        TESTS:

        Ensure that complex expressions do not cause a problem if they
        appear only as intermediate results as in :issue:`24536`::

            sage: x = SR.symbol('x', domain='real')
            sage: f = (x + I).abs()
            sage: f.find_local_minimum(-1,1)  # abs tol 1e-7                            # needs scipy
            (1.0, 1.6937685757340167e-08)

        ALGORITHM:

        Uses :func:`sage.numerical.optimize.find_local_minimum`.

        AUTHORS:

        - William Stein (2007-12-07)
        """
        from sage.numerical.optimize import find_local_minimum
        from sage.ext.fast_callable import fast_callable

        if var is None:
            var = self.default_variable()

        # We use the plotting fast_callable() and its wrapper
        # because numerical minimization is a lot like plotting in
        # the sense that we usually want to tolerate intermediate
        # complex expressions like abs(x + I), so long as the end
        # result has a "negligible" complex part.
        from sage.ext.fast_callable import FastCallableFloatWrapper
        f = self._plot_fast_callable(var)
        ff = FastCallableFloatWrapper(f, imag_tol=imaginary_tolerance)

        return find_local_minimum(ff, a=a, b=b, tol=tol, maxfun=maxfun)

    ###################
    # Fast Evaluation #
    ###################

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
        r"""
        Pretty-print this symbolic expression.

        This typesets it nicely and prints it immediately.

        OUTPUT:

        This method does not return anything. Like ``print``, output
        is sent directly to the screen.

        Note that the output depends on the display preferences. For details,
        see :func:`~sage.repl.rich_output.pretty_print.pretty_print`.

        EXAMPLES::

            sage: (x^2 + 1).show()
            x^2 + 1

        TESTS::

            sage: dm = get_display_manager()
            sage: dm.preferences.text = 'ascii_art'

        EXAMPLES::

            sage: %display ascii_art  # not tested
            sage: (x^2 + 1).show()
             2
            x  + 1

        TESTS:

        After the previous example, we need to reset the text display
        preferences::

            sage: dm.preferences.text = None
        """
        from sage.repl.rich_output.pretty_print import pretty_print
        pretty_print(self)

    def plot(self, *args, **kwds):
        """
        Plot a symbolic expression. All arguments are passed onto the standard :func:`plot` command.

        EXAMPLES:

        This displays a straight line::

            sage: sin(2).plot((x,0,3))                                                  # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        This draws a red oscillatory curve::

            sage: sin(x^2).plot((x,0,2*pi), rgbcolor=(1,0,0))                           # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Another plot using the variable theta::

            sage: var('theta')
            theta
            sage: (cos(theta) - erf(theta)).plot((theta,-2*pi,2*pi))                    # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        A very thick green plot with a frame::

            sage: sin(x).plot((x, -4*pi, 4*pi),                                         # needs sage.plot
            ....:             thickness=20, rgbcolor=(0,0.7,0)).show(frame=True)

        You can embed 2d plots in 3d space as follows::

            sage: plot(sin(x^2), (x, -pi, pi), thickness=2).plot3d(z=1)         # long time, needs sage.plot
            Graphics3d Object

        A more complicated family::

            sage: G = sum(plot(sin(n*x), (x, -2*pi, 2*pi)).plot3d(z=n)                  # needs sage.plot
            ....:         for n in [0,0.1,..1])
            sage: G.show(frame_aspect_ratio=[1,1,1/2])  # long time (5s on sage.math, 2012), needs sage.plot

        A plot involving the floor function::

            sage: plot(1.0 - x * floor(1/x), (x,0.00001,1.0))                           # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Sage used to allow symbolic functions with "no arguments";
        this no longer works::

            sage: plot(2*sin, -4, 4)                                                    # needs sage.plot
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
            'Integer Ring' and '<class 'sage.functions.trig.Function_sin'>'

        You should evaluate the function first::

            sage: plot(2*sin(x), -4, 4)                                                 # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        TESTS::

            sage: f(x) = x*(1 - x)
            sage: plot(f, 0, 1)                                                         # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.plot import plot

        # see if the user passed a variable in.
        if 'param' in kwds:
            param = kwds['param']
        else:
            param = None
            for i, arg in enumerate(args):
                if isinstance(arg, Expression) and arg.is_symbol():
                    param = arg
                    args = args[:i] + args[i+1:]
                    break

        if param is None:
            if self.is_callable():
                A = self.arguments()
                if len(A) == 0:
                    raise ValueError("function has no input arguments")
                else:
                    param = A[0]

                f = self._plot_fast_callable(param)
            else:
                A = self.variables()
                if len(A) == 0:
                    # Here we handle the case where f is something
                    # like ``sin``, which has takes arguments which
                    # aren't explicitly given
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

            sage: x = var('x', domain='real')
            sage: s = abs((1+I*x)^4)
            sage: f = s._plot_fast_callable(x)
            sage: abs(f(10) - abs((I*10+1)^4)) < 1e-11
            True
            sage: plot(s)                                                               # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Check that :issue:`15030` is fixed::

            sage: abs(log(x))._plot_fast_callable(x)(-0.2) # abs tol 1e-10
            3.52985761682672
            sage: f = function('f', evalf_func=lambda self,x,parent: I*x)
            sage: plot(abs(f(x)), 0,5)                                                  # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.ext.fast_callable import fast_callable
        from sage.rings.complex_double import CDF
        return fast_callable(self, vars=vars, expect_one_var=True, domain=CDF)

    ############
    # Calculus #
    ############
    def sum(self, *args, **kwds):
        r"""
        Return the symbolic sum `\sum_{v = a}^b` ``self``.

        with respect to the variable `v` with endpoints
        `a` and `b`.

        INPUT:

        - ``v`` -- a variable or variable name

        - ``a`` -- lower endpoint of the sum

        - ``b`` -- upper endpoint of the sum

        - ``algorithm`` -- (default: ``'maxima'``)  one of

          - ``'maxima'`` -- use Maxima (the default)
          - ``'maple'`` -- (optional) use Maple
          - ``'mathematica'`` -- (optional) use Mathematica
          - ``'giac'`` -- (optional) use Giac
          - ``'sympy'`` -- use SymPy

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

        A divergent geometric series.  Do not forget
        to `forget` your assumptions::

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

            sage: # needs giac
            sage: (sum(1/(1+k^2), k, -oo, oo, algorithm = 'giac')).factor()
            pi*(e^(2*pi) + 1)/((e^pi + 1)*(e^pi - 1))

        Use Maple as a backend for summation::

            sage: (binomial(n,k)*x^k).sum(k, 0, n, algorithm = 'maple')      # optional - maple
            (x + 1)^n

        .. NOTE::

           #. Sage can currently only understand a subset of the output of Maxima, Maple and
              Mathematica, so even if the chosen backend can perform the summation the
              result might not be convertible into a usable Sage expression.

        TESTS:

        Check that the sum in :issue:`10682` is done right::

            sage: sum(binomial(n,k)*k^2, k, 2, n)
            1/4*(n^2 + n)*2^n - n

        This sum used to give a wrong result (:issue:`9635`) but
        now gives correct results with all relevant assumptions::

            sage: (n,k,j)=var('n,k,j')
            sage: sum(binomial(n,k)*binomial(k-1,j)*(-1)**(k-1-j),k,j+1,n)
            -(-1)^j*sum((-1)^k*binomial(k - 1, j)*binomial(n, k), k, j + 1, n)
            sage: assume(j>-1)
            sage: sum(binomial(n,k)*binomial(k-1,j)*(-1)**(k-1-j),k,j+1,n)
            1
            sage: forget()
            sage: assume(n>=j)
            sage: sum(binomial(n,k)*binomial(k-1,j)*(-1)**(k-1-j),k,j+1,n)
            -(-1)^j*sum((-1)^k*binomial(k - 1, j)*binomial(n, k), k, j + 1, n)
            sage: forget()
            sage: assume(j==-1)
            sage: sum(binomial(n,k)*binomial(k-1,j)*(-1)**(k-1-j),k,j+1,n)
            1
            sage: forget()
            sage: assume(j<-1)
            sage: sum(binomial(n,k)*binomial(k-1,j)*(-1)**(k-1-j),k,j+1,n)
            -(-1)^j*sum((-1)^k*binomial(k - 1, j)*binomial(n, k), k, j + 1, n)
            sage: forget()

        Check that :issue:`16176` is fixed::

            sage: n = var('n')
            sage: sum(log(1-1/n^2),n,2,oo)
            -log(2)

        Check that :issue:`21801` is fixed::

            sage: n = SR.var('n')
            sage: sum(1/((n+1)*(2*n-1)), n, 0, oo)
            2/3*log(2) - 2/3
            sage: _.n()
            -0.204568546293370
            sage: f(n) = (-1)^(n+1)/(3*n+6*(-1)^n)
            sage: sum(f(2*n)+f(2*n+1), n, 0, oo)
            1/3*log(2) - 1/3
        """
        from sage.calculus.calculus import symbolic_sum
        return symbolic_sum(self, *args, **kwds)

    def prod(self, *args, **kwds):
        r"""
        Return the symbolic product `\prod_{v = a}^b` ``self``.

        This is the product respect to the variable `v` with endpoints `a` and `b`.

        INPUT:

        - ``expression`` -- a symbolic expression

        - ``v`` -- a variable or variable name

        - ``a`` -- lower endpoint of the product

        - ``b`` -- upper endpoint of the product

        - ``algorithm`` -- (default: ``'maxima'``)  one of

          - ``'maxima'`` -- use Maxima (the default)
          - ``'giac'`` -- (optional) use Giac
          - ``'sympy'`` -- use SymPy

        - ``hold`` -- boolean (default: ``False``); if ``True``, don't evaluate

        TESTS:

            sage: i, k, n = var('i,k,n')
            sage: k.prod(k, 1, n)
            factorial(n)
            sage: (x + i*(i+1)/2).prod(i,1,4)
            x^4 + 20*x^3 + 127*x^2 + 288*x + 180
            sage: (i^2).prod(i,1,7)
            25401600
            sage: f = function('f')
            sage: f(i).prod(i,1,7)
            f(7)*f(6)*f(5)*f(4)*f(3)*f(2)*f(1)
            sage: f(i).prod(i,1,n)
            product(f(i), i, 1, n)
            sage: assume(k>0)
            sage: (x^k).integrate(x,0,1).prod(k,1,n)
            1/factorial(n + 1)
            sage: f(i).prod(i,1,n).log().log_expand()
            sum(log(f(i)), i, 1, n)
        """
        from sage.calculus.calculus import symbolic_product
        return symbolic_product(self, *args, **kwds)

    def integral(self, *args, **kwds):
        """
        Compute the integral of ``self``.

        Please see :func:`sage.symbolic.integration.integral.integrate` for more details.

        EXAMPLES::

            sage: sin(x).integral(x,0,3)
            -cos(3) + 1
            sage: sin(x).integral(x)
            -cos(x)

        TESTS:

        We check that :issue:`12438` is resolved::

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

        We check that :issue:`13097` is resolved (sage doesn't
        crash). If giac is available, you may even get a usable
        answer::

            sage: f = ln(1+4/5*sin(x))
            sage: integrate(f, x, -3.1415, 3.1415)  # random
            integrate(log(4/5*sin(x) + 1), x, -3.14150000000000,
            3.14150000000000)
            sage: # needs sage.libs.giac
            sage: integrate(f, x, -3.1415, 3.1415)  # tol 10e-6
            -1.40205228301000
        """
        from sage.symbolic.integration.integral import \
            integral, _normalize_integral_input
        R = self._parent
        if isinstance(R, sage.rings.abc.CallableSymbolicExpressionRing):
            from sage.symbolic.ring import SR
            f = SR(self)
            f, v, a, b = _normalize_integral_input(f, *args)
            # Definite integral with respect to a positional variable.
            if a is not None and v in R.arguments():
                arguments = list(R.arguments())
                arguments.remove(v)
                if arguments:
                    arguments = tuple(arguments)
                    from sage.symbolic.callable import CallableSymbolicExpressionRing
                    R = CallableSymbolicExpressionRing(arguments, check=False)
                else:   # all arguments are gone
                    R = SR
            return R(integral(f, v, a, b, **kwds))
        return integral(self, *args, **kwds)

    integrate = integral

    def nintegral(self, *args, **kwds):
        """
        Compute the numerical integral of ``self``.

        Please see :obj:`sage.calculus.calculus.nintegral` for more details.

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
        Return a symbolic limit.

        See :obj:`sage.calculus.calculus.limit`

        EXAMPLES::

            sage: (sin(x)/x).limit(x=0)
            1
        """
        from sage.calculus.calculus import limit
        return limit(self, *args, **kwds)

    def laplace(self, t, s):
        """
        Return Laplace transform of ``self``.

        See :obj:`sage.calculus.calculus.laplace`

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
        Return inverse Laplace transform of ``self``.

        See :obj:`sage.calculus.calculus.inverse_laplace`

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
        Return a relation obtained by adding ``x`` to both sides of
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
            raise TypeError("this expression must be a relation")
        return self + x

    def subtract_from_both_sides(self, x):
        """
        Return a relation obtained by subtracting ``x`` from both sides
        of this relation.

        EXAMPLES::

            sage: eqn = x*sin(x)*sqrt(3) + sqrt(2) > cos(sin(x))
            sage: eqn.subtract_from_both_sides(sqrt(2))
            sqrt(3)*x*sin(x) > -sqrt(2) + cos(sin(x))
            sage: eqn.subtract_from_both_sides(cos(sin(x)))
            sqrt(3)*x*sin(x) + sqrt(2) - cos(sin(x)) > 0
        """
        if not is_a_relational(self._gobj):
            raise TypeError("this expression must be a relation")
        return self - x

    def multiply_both_sides(self, x, checksign=None):
        """
        Return a relation obtained by multiplying both sides of this
        relation by ``x``.

        .. NOTE::

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
            raise TypeError("this expression must be a relation")
        return self * x

    def divide_both_sides(self, x, checksign=None):
        """
        Return a relation obtained by dividing both sides of this
        relation by ``x``.

        .. NOTE::

           The ``checksign`` keyword argument is currently ignored and
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
            raise TypeError("this expression must be a relation")
        return self / x

    def implicit_derivative(self, Y, X, n=1):
        """
        Return the `n`-th derivative of `Y` with respect to `X` given
        implicitly by this expression.

        INPUT:

        - ``Y`` -- the dependent variable of the implicit expression

        - ``X`` -- the independent variable with respect to which the
          derivative is taken

        - ``n`` -- (default: 1) the order of the derivative

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: f = cos(x)*sin(y)
            sage: f.implicit_derivative(y, x)
            sin(x)*sin(y)/(cos(x)*cos(y))
            sage: g = x*y^2
            sage: g.implicit_derivative(y, x, 3)
            -1/4*(y + 2*y/x)/x^2 + 1/4*(2*y^2/x - y^2/x^2)/(x*y) - 3/4*y/x^3

        It is an error to not include an independent variable term
        in the expression::

            sage: (cos(x)*sin(x)).implicit_derivative(y, x)
            Traceback (most recent call last):
            ...
            ValueError: Expression cos(x)*sin(x) contains no y terms


        TESTS:

        Check that the symbols registry is not polluted::

            sage: var('x,y')
            (x, y)
            sage: psr = copy(SR.symbols)
            sage: (x^6*y^5).implicit_derivative(y, x, 3)
            -792/125*y/x^3 + 12/25*(15*x^4*y^5 + 28*x^3*y^5)/(x^6*y^4) - 36/125*(20*x^5*y^4 + 43*x^4*y^4)/(x^7*y^3)
            sage: psr == SR.symbols
            True
        """
        from sage.symbolic.ring import SR
        from sage.symbolic.function import SymbolicFunction

        if not self.has(Y):
            raise ValueError("Expression {} contains no {} terms".format(self, Y))
        with SR.temp_var() as x:
            with SR.temp_var() as yy:
                y = SymbolicFunction('y', 1)(x)
                f = SymbolicFunction('f', 2)(x, yy)
                Fx = f.diff(x)
                Fy = f.diff(yy)
                G = -(Fx/Fy)
                G = G.subs({yy: y})
                di = {y.diff(x): -self.diff(X)/self.diff(Y)}
                R = G
                S = G.diff(x, n - 1)
                for i in range(n + 1):
                    di[y.diff(x, i + 1).subs({x: x})] = R
                    S = S.subs(di)
                    R = G.diff(x, i)
                    for j in range(n + 1 - i):
                        di[f.diff(x, i, yy, j).subs({x: x, yy: y})] = self.diff(X, i, Y, j)
                        S = S.subs(di)
                return S


cpdef _repr_Expression(x):
    r"""
    Return the string representation of the expression ``x``.

    EXAMPLES::

        sage: SR._repr_element_(x+2)
        'x + 2'
    """
    return ccrepr((<Expression>x)._gobj)


cpdef _latex_Expression(x):
    r"""
    Return the standard LaTeX version of the expression `x`.

    EXAMPLES::

        sage: latex(sin(x+2))
        \sin\left(x + 2\right)
        sage: latex(var('theta') + 2)
        \theta + 2
    """
    return char_to_str(GEx_to_str_latex(&(<Expression>x)._gobj))


def solve_diophantine(f, *args, **kwds):
    """
    Solve a Diophantine equation.

    The argument, if not given as symbolic equation, is set equal to zero.
    It can be given in any form that can be converted to symbolic. Please
    see :meth:`Expression.solve_diophantine` for a detailed
    synopsis.

    EXAMPLES::

        sage: R.<a,b> = PolynomialRing(ZZ); R
        Multivariate Polynomial Ring in a, b over Integer Ring
        sage: solve_diophantine(a^2 - 3*b^2 + 1)
        []
        sage: sorted(solve_diophantine(a^2 - 3*b^2 + 2), key=str)
        [(-1/2*sqrt(3)*(sqrt(3) + 2)^t + 1/2*sqrt(3)*(-sqrt(3) + 2)^t - 1/2*(sqrt(3) + 2)^t - 1/2*(-sqrt(3) + 2)^t,
          -1/6*sqrt(3)*(sqrt(3) + 2)^t + 1/6*sqrt(3)*(-sqrt(3) + 2)^t - 1/2*(sqrt(3) + 2)^t - 1/2*(-sqrt(3) + 2)^t),
          (1/2*sqrt(3)*(sqrt(3) + 2)^t - 1/2*sqrt(3)*(-sqrt(3) + 2)^t + 1/2*(sqrt(3) + 2)^t + 1/2*(-sqrt(3) + 2)^t,
           1/6*sqrt(3)*(sqrt(3) + 2)^t - 1/6*sqrt(3)*(-sqrt(3) + 2)^t + 1/2*(sqrt(3) + 2)^t + 1/2*(-sqrt(3) + 2)^t)]
    """
    from sage.symbolic.ring import SR

    if not isinstance(f, Expression):
        f = SR(f)
    return f.solve_diophantine(*args, **kwds)


def _eval_on_operands(f):
    """
    Given a function ``f``, return a new function which takes a symbolic
    expression as first argument and prepends the operands of that
    expression to the arguments of ``f``.

    EXAMPLES::

        sage: def f(ex, x, y):
        ....:     '''
        ....:     Some documentation.
        ....:     '''
        ....:     return x + 2*y
        ....:
        sage: f(None, x, 1)
        x + 2
        sage: from sage.symbolic.expression import _eval_on_operands
        sage: g = _eval_on_operands(f)
        sage: var('a,b')
        (a, b)
        sage: g(a + b)
        a + 2*b
        sage: print(g.__doc__.strip())
        Some documentation.
    """
    @sage_wraps(f)
    def new_f(ex, *args, **kwds):
        new_args = list(ex._unpack_operands())
        new_args.extend(args)
        return f(ex, *new_args, **kwds)
    return new_f


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
        ....:     class EvaluationMethods():
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
        attribute 'argp1'...
        sage: t = (e + 1).op[0]; t
        tfunc(x)
        sage: t
        tfunc(x)
        sage: type(t)
        <class '__main__.Expression_with_dynamic_methods'>
        sage: t.argp1()
        x + 1
        sage: import sage.interfaces.tab_completion as s
        sage: s.completions('t.argp', globals())
        ['t.argp1']
        sage: t.argp1.__doc__.strip()
        'Some documentation about a bogus function.'

    Now with two arguments::

        sage: class TFunc2(BuiltinFunction):
        ....:     def __init__(self):
        ....:         BuiltinFunction.__init__(self, 'tfunc', nargs=2)
        ....:
        ....:     class EvaluationMethods():
        ....:         def argsum(self, x, y):
        ....:             return x + y
        ....:
        sage: tfunc2 = TFunc2()
        sage: e = tfunc2(x, 1)
        sage: e.argsum()
        x + 1
    """
    cls = dynamic_class_cache.get(serial)
    if cls is not None:
        return cls

    func_class = get_sfunction_from_serial(serial)
    try:
        eval_methods = func_class.EvaluationMethods
    except AttributeError:
        cls = Expression
    else:
        cls = dynamic_class('Expression_with_dynamic_methods',
                            (Expression,),
                            eval_methods, prepend_cls_bases=False)
        # Fix methods from eval_methods, wrapping them to extract
        # the operands and pass them as arguments
        for name, meth in eval_methods.__dict__.items():
            if not callable(meth):
                continue
            meth = _eval_on_operands(meth)
            setattr(cls, name, meth)

    dynamic_class_cache[serial] = cls
    return cls


cdef Expression new_Expression_from_GEx(parent, GEx juice):
    cdef type cls
    cdef Expression nex
    if is_exactly_a_function(juice):
        # if the function defines any dynamic methods these are made
        # available through a dynamic class
        cls = <type>get_dynamic_class_for_function(ex_to_function(juice).get_serial())
    else:
        cls = Expression

    nex = <Expression>cls.__new__(cls)
    nex._gobj = GEx(juice)
    nex._parent = parent
    return nex


cpdef new_Expression(parent, x):
    r"""
    Convert ``x`` into the symbolic expression ring ``parent``.

    This is the element constructor.

    EXAMPLES::

        sage: a = SR(-3/4); a
        -3/4
        sage: type(a)
        <class 'sage.symbolic.expression.Expression'>
        sage: a.parent()
        Symbolic Ring
        sage: K.<a> = QuadraticField(-3)                                                # needs sage.rings.number_field
        sage: a + sin(x)                                                                # needs sage.rings.number_field
        I*sqrt(3) + sin(x)
        sage: x = var('x'); y0,y1 = PolynomialRing(ZZ,2,'y').gens()
        sage: x+y0/y1
        x + y0/y1
        sage: x.subs(x=y0/y1)
        y0/y1
        sage: x + int(1)
        x + 1
    """
    cdef GEx exp
    if isinstance(x, Expression):
        return new_Expression_from_GEx(parent, (<Expression>x)._gobj)
    if hasattr(x, '_symbolic_'):
        return x._symbolic_(parent)
    elif isinstance(x, str):
        try:
            from sage.calculus.calculus import symbolic_expression_from_string
            return parent(symbolic_expression_from_string(x))
        except SyntaxError as err:
            msg, s, pos = err.args
            raise TypeError("%s: %s !!! %s" % (msg, s[:pos], s[pos:]))

    from sage.rings.infinity import (infinity, minus_infinity,
                                     unsigned_infinity)
    from sage.structure.factorization import Factorization
    from sage.categories.sets_cat import Sets

    if isinstance(x, RealNumber):
        if x.is_NaN():
            from sage.symbolic.constants import NaN
            return NaN
        exp = x
    elif isinstance(x, (float, complex)):
        if not (x == x):
            from sage.symbolic.constants import NaN
            return NaN
        exp = x
    elif isinstance(x, int):
        try:
            exp = GEx(<long>x)
        except OverflowError:
            exp = x
    elif x is infinity:
        return new_Expression_from_GEx(parent, g_Infinity)
    elif x is minus_infinity:
        return new_Expression_from_GEx(parent, g_mInfinity)
    elif x is unsigned_infinity:
        return new_Expression_from_GEx(parent, g_UnsignedInfinity)
    elif isinstance(x, (RingElement, Matrix)):
        if x.parent().characteristic():
            raise TypeError('positive characteristic not allowed in symbolic computations')
        exp = x
    elif isinstance(x, Factorization):
        from sage.misc.misc_c import prod
        from sage.symbolic.ring import SR

        return prod([SR(p)**e for p,e in x], SR(x.unit()))
    elif x in Sets():
        from sage.rings.integer_ring import ZZ
        from sage.rings.qqbar import AA
        from sage.rings.rational_field import QQ
        from sage.rings.semirings.non_negative_integer_semiring import NN
        from sage.sets.real_set import RealSet
        if (x.is_finite() or x in (NN, ZZ, QQ, AA)
                or isinstance(x, RealSet)):
            exp = x
        else:
            raise TypeError(f"unable to convert {x!r} to a symbolic expression")
    else:
        raise TypeError(f"unable to convert {x!r} to a symbolic expression")

    return new_Expression_from_GEx(parent, exp)


cpdef new_Expression_from_pyobject(parent, x, bint force=True, bint recursive=True):
    r"""
    Wrap the given Python object in a symbolic expression even if it
    cannot be coerced to the Symbolic Ring.

    INPUT:

    - ``parent`` -- a symbolic ring

    - ``x`` -- a Python object

    - ``force`` -- boolean (default: ``True``); if ``True``, the Python object
      is taken as is without attempting coercion or list traversal

    - ``recursive`` -- boolean (default: ``True``); disables recursive
      traversal of lists

    EXAMPLES::

        sage: t = SR._force_pyobject(QQ); t   # indirect doctest
        Rational Field
        sage: type(t)
        <class 'sage.symbolic.expression.Expression'>

        sage: from sage.symbolic.expression import new_Expression_from_pyobject
        sage: t = new_Expression_from_pyobject(SR, 17); t
        17
        sage: type(t)
        <class 'sage.symbolic.expression.Expression'>

        sage: t2 = new_Expression_from_pyobject(SR, t, False); t2
        17
        sage: t2 is t
        True

        sage: tt = new_Expression_from_pyobject(SR, t, True); tt
        17
        sage: tt is t
        False
    """
    cdef GEx exp
    cdef GExprSeq ex_seq
    cdef GExVector ex_v
    if force:
        exp = x

    else:
        # first check if we can do it the nice way
        if isinstance(x, Expression):
            return x
        try:
            return parent.coerce(x)
        except TypeError:
            pass

        # tuples can be packed into exprseq
        if isinstance(x, (tuple, list)):
            from sage.symbolic.ring import SR
            for e in x:
                obj = SR._force_pyobject(e, force=(not recursive))
                ex_v.push_back((<Expression>obj)._gobj)

            ex_seq = GExprSeq(ex_v)

            exp = GEx(ex_seq)
        else:
            exp = x

    return new_Expression_from_GEx(parent, exp)


cpdef new_Expression_wild(parent, unsigned int n=0):
    r"""
    Return the n-th wild-card for pattern matching and substitution.

    INPUT:

    - ``parent`` -- a symbolic ring

    - ``n`` -- nonnegative integer

    OUTPUT: n-th wildcard expression

    EXAMPLES::

        sage: x,y = var('x,y')
        sage: w0 = SR.wild(0); w1 = SR.wild(1)
        sage: pattern = sin(x)*w0*w1^2; pattern
        $1^2*$0*sin(x)
        sage: f = atan(sin(x)*3*x^2); f
        arctan(3*x^2*sin(x))
        sage: f.has(pattern)
        True
        sage: f.subs(pattern == x^2)
        arctan(x^2)
    """
    return new_Expression_from_GEx(parent, g_wild(n))


cpdef new_Expression_symbol(parent, name=None, latex_name=None, domain=None):
    r"""
    Look up or create a symbol.

    EXAMPLES::

        sage: t0 = SR.symbol("t0")
        sage: t0.conjugate()
        conjugate(t0)

        sage: t1 = SR.symbol("t1", domain='real')
        sage: t1.conjugate()
        t1

        sage: t0.abs()
        abs(t0)

        sage: t0_2 = SR.symbol("t0", domain='positive')
        sage: t0_2.abs()
        t0
        sage: bool(t0_2 == t0)
        True
        sage: t0.conjugate()
        t0

        sage: SR.symbol() # temporary variable
        symbol...
    """
    cdef GSymbol symb
    cdef Expression e

    # check if there is already a symbol with same name
    e = parent.symbols.get(name)

    # fast path to get an already existing variable
    if e is not None:
        if domain is None:
            if latex_name is None:
                return e

        # get symbol
        symb = ex_to_symbol(e._gobj)
        if latex_name is not None:
            symb.set_texname(str_to_bytes(latex_name))
        if domain is not None:
            symb.set_domain(sage_domain_to_ginac_domain(domain))
        e._gobj = GEx(symb)
        if domain is not None:
            send_sage_domain_to_maxima(e, domain)

        return e

    else:  # initialize a new symbol
        # Construct expression
        e = <Expression>Expression.__new__(Expression)
        e._parent = parent

        if name is None:  # Check if we need a temporary anonymous new symbol
            symb = ginac_new_symbol()
            name = symb.get_name().decode('ascii')
            if domain is not None:
                symb.set_domain(sage_domain_to_ginac_domain(domain))
        else:
            if latex_name is None:
                latex_name = latex_variable_name(name)
            if domain is not None:
                ginac_domain = sage_domain_to_ginac_domain(domain)
            else:
                ginac_domain = domain_complex
            symb = ginac_symbol(str_to_bytes(name),
                                str_to_bytes(latex_name), ginac_domain)

        e._gobj = GEx(symb)
        parent.symbols[name] = e
        if domain is not None:
            send_sage_domain_to_maxima(e, domain)

    return e


cdef unsigned sage_domain_to_ginac_domain(object domain) except? 3474701533:
    """
    TESTS::

        sage: var('x', domain='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo': domain must be one of 'complex', 'real', 'positive' or 'integer'
    """
    # convert the domain argument to something easy to parse
    if domain is RR or domain == 'real':
        return domain_real
    elif domain == 'positive':
        return domain_positive
    elif isinstance(domain, sage.rings.abc.ComplexField) or domain == 'complex':
        return domain_complex
    elif domain is ZZ or domain == 'integer':
        return domain_integer
    else:
        raise ValueError(repr(domain)+": domain must be one of 'complex', 'real', 'positive' or 'integer'")

cdef void send_sage_domain_to_maxima(Expression v, object domain) except +:
    from sage.symbolic.assumptions import assume
    # convert the domain argument to something easy to parse
    if domain is RR or domain == 'real':
        assume(v, 'real')
    elif domain == 'positive':
        assume(v>0)
    elif isinstance(domain, sage.rings.abc.ComplexField) or domain == 'complex':
        assume(v, 'complex')
    elif domain is ZZ or domain == 'integer':
        assume(v, 'integer')
    else:
        raise ValueError(repr(domain)+": domain must be one of 'complex', 'real', 'positive' or 'integer'")


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
            sage: next(i)
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
    cdef ExpressionIterator m = <ExpressionIterator>ExpressionIterator.__new__(ExpressionIterator)
    m._ex = ex
    m._ind = 0
    m._len = ex._gobj.nops()
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
        raise TypeError("incompatible relations")
    elif lop == equal:
        return rop
    elif rop == equal:
        return lop
    elif lop in [less, less_or_equal] and rop in [less, less_or_equal]:
        return less
    elif lop in [greater, greater_or_equal] and rop in [greater, greater_or_equal]:
        return greater
    else:
        raise TypeError("incompatible relations")

cdef class hold_class:
    """
    Instances of this class can be used with Python `with`.

    EXAMPLES::

        sage: with hold:
        ....:     tan(1/12*pi)
        ....:
        tan(1/12*pi)
        sage: tan(1/12*pi)
        -sqrt(3) + 2
        sage: with hold:
        ....:     2^5
        ....:
        32
        sage: with hold:
        ....:     SR(2)^5
        ....:
        2^5
        sage: with hold:
        ....:     t=tan(1/12*pi)
        ....:
        sage: t
        tan(1/12*pi)
        sage: t.unhold()
        -sqrt(3) + 2
    """
    def __enter__(self):
        """
        EXAMPLES::

            sage: hold.__enter__()
            sage: SR(2)^5
            2^5
            sage: hold.__exit__()
            sage: SR(2)^5
            32
        """
        g_set_state(b'hold', True)

    def __exit__(self, *args):
        """
        EXAMPLES::

            sage: hold.__enter__()
            sage: SR(2)^5
            2^5
            sage: hold.__exit__()
            sage: SR(2)^5
            32
        """
        g_set_state(b'hold', False)

    def start(self):
        """
        Start a hold context.

        EXAMPLES::

            sage: hold.start()
            sage: SR(2)^5
            2^5
            sage: hold.stop()
            sage: SR(2)^5
            32
        """
        self.__enter__()

    def stop(self):
        """
        Stop any hold context.

        EXAMPLES::

            sage: hold.start()
            sage: SR(2)^5
            2^5
            sage: hold.stop()
            sage: SR(2)^5
            32
        """
        self.__exit__()

hold = hold_class()


include "comparison_impl.pxi"
include "constants_c_impl.pxi"
include "getitem_impl.pxi"
include "pynac_constant_impl.pxi"
include "pynac_function_impl.pxi"
include "series_impl.pxi"
include "substitution_map_impl.pxi"
