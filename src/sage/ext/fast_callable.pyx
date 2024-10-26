r"""
Fast Expression Evaluation

For many applications such as numerical integration, differential
equation approximation, plotting a 3d surface, optimization problems,
Monte-Carlo simulations, etc., one wishes to pass around and evaluate
a single algebraic expression many, many times at various floating
point values.  Other applications may need to evaluate an expression
many times in interval arithmetic, or in a finite field.  Doing this
via recursive calls over a python representation of the object (even
if Maxima or other outside packages are not involved) is extremely
inefficient.

This module provides a function, :func:`fast_callable`, to
transform such expressions into a form where they can be evaluated
quickly::

    sage: # needs sage.symbolic
    sage: f = sin(x) + 3*x^2
    sage: ff = fast_callable(f, vars=[x])
    sage: ff(3.5)
    36.3992167723104
    sage: ff(RIF(3.5))
    36.39921677231038?

By default, :func:`fast_callable` only removes some interpretive
overhead from the evaluation, but all of the individual arithmetic
operations are done using standard Sage arithmetic.  This is still a
huge win over :mod:`sage.calculus`, which evidently has a lot of overhead.
Compare the cost of evaluating Wilkinson's polynomial (in unexpanded
form) at `x = 30`::

    sage: # needs sage.symbolic
    sage: wilk = prod((x-i) for i in [1 .. 20]); wilk
    (x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5)*(x - 6)*(x - 7)*(x - 8)*(x - 9)*(x - 10)*(x - 11)*(x - 12)*(x - 13)*(x - 14)*(x - 15)*(x - 16)*(x - 17)*(x - 18)*(x - 19)*(x - 20)
    sage: timeit('wilk.subs(x=30)')         # random    # long time
    625 loops, best of 3: 1.43 ms per loop
    sage: fc_wilk = fast_callable(wilk, vars=[x])
    sage: timeit('fc_wilk(30)')             # random    # long time
    625 loops, best of 3: 9.72 us per loop

You can specify a particular domain for the evaluation using
``domain=``::

    sage: fc_wilk_zz = fast_callable(wilk, vars=[x], domain=ZZ)                         # needs sage.symbolic

The meaning of ``domain=D`` is that each intermediate and final result
is converted to type ``D``.  For instance, the previous example of
``sin(x) + 3*x^2`` with ``domain=D`` would be equivalent to
``D(D(sin(D(x))) + D(D(3)*D(D(x)^2)))``.  (This example also
demonstrates the one exception to the general rule: if an exponent is an
integral constant, then it is not wrapped with ``D()``.)

At first glance, this seems like a very bad idea if you want to
compute quickly.  And it is a bad idea, for types where we don't
have a special interpreter.  It's not too bad of a slowdown, though.
To mitigate the costs, we check whether the value already has
the correct parent before we call ``D``.

We don't yet have a special interpreter with domain ``ZZ``, so we can see
how that compares to the generic ``fc_wilk`` example above::

    sage: timeit('fc_wilk_zz(30)')          # random    # long time                     # needs sage.symbolic
    625 loops, best of 3: 15.4 us per loop

However, for other types, using ``domain=D`` will get a large speedup,
because we have special-purpose interpreters for those types.  One
example is ``RDF``.  Since with ``domain=RDF`` we know that every single
operation will be floating-point, we can just execute the
floating-point operations directly and skip all the Python object
creations that you would get from actually using ``RDF`` objects::

    sage: fc_wilk_rdf = fast_callable(wilk, vars=[x], domain=RDF)                       # needs sage.symbolic
    sage: timeit('fc_wilk_rdf(30.0)')       # random    # long time                     # needs sage.symbolic
    625 loops, best of 3: 7 us per loop

The domain does not need to be a Sage type; for instance, ``domain=float``
also works.  (We actually use the same fast interpreter for ``domain=float``
and ``domain=RDF``; the only difference is that when ``domain=RDF`` is used,
the return value is an ``RDF`` element, and when ``domain=float`` is used,
the return value is a Python :class:`float`.) ::

    sage: fc_wilk_float = fast_callable(wilk, vars=[x], domain=float)                   # needs sage.symbolic
    sage: timeit('fc_wilk_float(30.0)')     # random    # long time                     # needs sage.symbolic
    625 loops, best of 3: 5.04 us per loop

We also have support for ``RR``::

    sage: fc_wilk_rr = fast_callable(wilk, vars=[x], domain=RR)                         # needs sage.symbolic
    sage: timeit('fc_wilk_rr(30.0)')        # random    # long time                     # needs sage.symbolic
    625 loops, best of 3: 13 us per loop

For ``CC``::

    sage: fc_wilk_cc = fast_callable(wilk, vars=[x], domain=CC)                         # needs sage.symbolic
    sage: timeit('fc_wilk_cc(30.0)')        # random    # long time                     # needs sage.symbolic
    625 loops, best of 3: 23 us per loop

And support for ``CDF``::

    sage: fc_wilk_cdf = fast_callable(wilk, vars=[x], domain=CDF)                       # needs sage.symbolic
    sage: timeit('fc_wilk_cdf(30.0)')       # random    # long time                     # needs sage.symbolic
    625 loops, best of 3: 10.2 us per loop

Currently, :func:`fast_callable` can accept two kinds of objects:
polynomials (univariate and multivariate) and symbolic expressions
(elements of the Symbolic Ring).  For polynomials, you can omit the
``vars`` argument; the variables will default to the ring generators (in
the order used when creating the ring). ::

    sage: K.<x,y,z> = QQ[]
    sage: p = 10*y + 100*z + x
    sage: fp = fast_callable(p)
    sage: fp(1,2,3)
    321

But you can also specify the variable names to override the default
ordering (you can include extra variable names here, too). ::

    sage: fp = fast_callable(p, vars=('x','w','z','y'))

For symbolic expressions, you need to specify the variable names, so
that :func:`fast_callable` knows what order to use. ::

    sage: # needs sage.symbolic
    sage: var('y,z,x')
    (y, z, x)
    sage: f = 10*y + 100*z + x
    sage: ff = fast_callable(f, vars=(x,y,z))
    sage: ff(1,2,3)
    321

You can also specify extra variable names::

    sage: ff = fast_callable(f, vars=('x','w','z','y'))                                 # needs sage.symbolic
    sage: ff(1,2,3,4)                                                                   # needs sage.symbolic
    341

This should be enough for normal use of :func:`fast_callable`; let's
discuss some more advanced topics.

Sometimes it may be useful to create a fast version of an expression
without going through symbolic expressions or polynomials; perhaps
because you want to describe to :func:`fast_callable` an expression
with common subexpressions.

Internally, :func:`fast_callable` works in two stages: it constructs
an expression tree from its argument, and then it builds a
fast evaluator from that expression tree.  You can bypass the first phase
by building your own expression tree and passing that directly to
:func:`fast_callable`, using an :class:`ExpressionTreeBuilder`. ::

    sage: from sage.ext.fast_callable import ExpressionTreeBuilder
    sage: etb = ExpressionTreeBuilder(vars=('x','y','z'))

An :class:`ExpressionTreeBuilder` has three interesting methods:
:meth:`constant`, :meth:`var`, and :meth:`call`.
All of these methods return :class:`ExpressionTree` objects.

The :meth:`var` method takes a string, and returns an expression tree
for the corresponding variable. ::

    sage: x = etb.var('x')
    sage: y = etb.var('y')
    sage: z = etb.var('y')

Expression trees support Python's numeric operators, so you can easily
build expression trees representing arithmetic expressions. ::

    sage: v1 = (x+y)*(y+z) + (y//z)

The :meth:`constant` method takes a Sage value, and returns an
expression tree representing that value. ::

    sage: v2 = etb.constant(3.14159) * x + etb.constant(1729) * y

The :meth:`call` method takes a sage/Python function and zero or more
expression trees, and returns an expression tree representing
the function call. ::

    sage: v3 = etb.call(sin, v1+v2)
    sage: v3                                                                            # needs sage.symbolic
    sin(add(add(mul(add(v_0, v_1), add(v_1, v_1)), floordiv(v_1, v_1)),
            add(mul(3.14159000000000, v_0), mul(1729, v_1))))

Many sage/Python built-in functions are specially handled; for instance,
when evaluating an expression involving :func:`sin` over ``RDF``,
the C math library function :func:`sin` is called.  Arbitrary functions
are allowed, but will be much slower since they will call back to
Python code on every call; for example, the following will work. ::

    sage: def my_sqrt(x): return pow(x, 0.5)
    sage: e = etb.call(my_sqrt, v1); e
    {my_sqrt}(add(mul(add(v_0, v_1), add(v_1, v_1)), floordiv(v_1, v_1)))
    sage: fast_callable(e)(1, 2, 3)
    3.60555127546399

To provide :func:`fast_callable` for your own class (so that
``fast_callable(x)`` works when ``x`` is an instance of your
class), implement a method ``_fast_callable_(self, etb)`` for your class.
This method takes an :class:`ExpressionTreeBuilder`, and returns an
expression tree built up using the methods described above.

EXAMPLES::

    sage: var('x')                                                                      # needs sage.symbolic
    x
    sage: f = fast_callable(sqrt(x^7+1), vars=[x], domain=float)                        # needs sage.symbolic

::

    sage: f(1)                                                                          # needs sage.symbolic
    1.4142135623730951
    sage: f.op_list()                                                                   # needs sage.symbolic
    [('load_arg', 0), ('ipow', 7), ('load_const', 1.0), 'add', 'sqrt', 'return']

To interpret that last line, we load argument 0 (``x`` in this case) onto
the stack, push the constant `7.0` onto the stack, call the :func:`pow`
functionn(which takes 2 arguments from the stack), push the constant `1.0`,
add the top two arguments of the stack, and then call :func:`sqrt`.

Here we take :func:`sin` of the first argument and add it to ``f``::

    sage: from sage.ext.fast_callable import ExpressionTreeBuilder
    sage: etb = ExpressionTreeBuilder('x')
    sage: x = etb.var('x')
    sage: f = etb.call(sqrt, x^7 + 1)
    sage: g = etb.call(sin, x)
    sage: fast_callable(f+g).op_list()
    [('load_arg', 0), ('ipow', 7), ('load_const', 1), 'add',
     ('py_call', <function sqrt at ...>, 1), ('load_arg', 0), ('py_call', sin, 1),
     'add', 'return']


AUTHOR:

- Carl Witty (2009-02): initial version (heavily inspired by
  Robert Bradshaw's fast_eval.pyx)

.. TODO::

    The following bits of text were written for the module docstring.
    They are not true yet, but I hope they will be true someday, at
    which point I will move them into the main text.

    The final interesting method of :class:`ExpressionTreeBuilder` is
    :meth:`choice`.  This produces conditional expressions, like the C
    ``COND ? T : F`` expression or the Python ``T if COND else F``.
    This lets you define piecewise functions using :func:`fast_callable`. ::

        sage: v4 = etb.choice(v3 >= etb.constant(0), v1, v2)  # not tested

    The arguments are ``(COND, T, F)`` (the same order as in C), so the
    above means that if ``v3`` evaluates to a nonnegative number,
    then ``v4`` will evaluate to the result of ``v1``;
    otherwise, ``v4`` will evaluate to the result of ``v2``.

    Let's see an example where we see that :func:`fast_callable` does not
    evaluate common subexpressions more than once.  We'll make a
    :func:`fast_callable` expression that gives the result
    of 16 iterations of the Mandelbrot function. ::

        sage: etb = ExpressionTreeBuilder('c')
        sage: z = etb.constant(0)
        sage: c = etb.var('c')
        sage: for i in range(16):
        ....:     z = z*z + c
        sage: mand = fast_callable(z, domain=CDF)                                       # needs sage.rings.complex_double

    Now ``ff`` does 32 complex arithmetic operations on each call
    (16 additions and 16 multiplications).  However, if ``z*z`` produced
    code that evaluated ``z`` twice, then this would do many
    thousands of arithmetic operations instead.

    Note that the handling for common subexpressions only checks whether
    expression trees are the same Python object; for instance, the following
    code will evaluate ``x+1`` twice::

        sage: etb = ExpressionTreeBuilder('x')
        sage: x = etb.var('x')
        sage: (x+1)*(x+1)
        mul(add(v_0, 1), add(v_0, 1))

    but this code will only evaluate ``x+1`` once::

        sage: v = x+1; v*v
        mul(add(v_0, 1), add(v_0, 1))
"""


#*****************************************************************************
#       Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2009 Carl Witty <Carl.Witty@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator
from copy import copy
from math import isnan, nan

import sage.rings.abc
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.element cimport parent
from sage.structure.element cimport Expression as Expression_abc


def fast_callable(x, domain=None, vars=None,
                  expect_one_var=False):
    r"""
    Given an expression ``x``, compile it into a form that can be quickly
    evaluated, given values for the variables in ``x``.

    Currently, ``x`` can be an expression object, an element of ``SR``, or a
    (univariate or multivariate) polynomial; this list will probably
    be extended soon.

    By default, ``x`` is evaluated the same way that a Python function
    would evaluate it -- addition maps to ``PyNumber_Add``, etc.  However,
    you can specify ``domain=D`` where ``D`` is some Sage parent or Python
    type; in this case, all arithmetic is done in that domain.  If we
    have a special-purpose interpreter for that parent (like ``RDF`` or
    :class:`float`), ``domain=...`` will trigger the use of that interpreter.

    If ``vars`` is ``None`` and ``x`` is a polynomial, then we will use the
    generators of ``parent(x)`` as the variables; otherwise, ``vars`` must be
    specified (unless ``x`` is a symbolic expression with only one variable,
    and ``expect_one_var`` is ``True``, in which case we will use that variable).

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: var('x')
        x
        sage: expr = sin(x) + 3*x^2
        sage: f = fast_callable(expr, vars=[x])
        sage: f(2)
        sin(2) + 12
        sage: f(2.0)
        12.9092974268257

    We have special fast interpreters for ``domain=float`` and ``domain=RDF``.
    (Actually it's the same interpreter; only the return type varies.)
    Note that the float interpreter is not actually more accurate than
    the ``RDF`` interpreter; elements of ``RDF`` just don't display all
    their digits. We have special fast interpreter for ``domain=CDF``::

        sage: # needs sage.symbolic
        sage: f_float = fast_callable(expr, vars=[x], domain=float)
        sage: f_float(2)
        12.909297426825681
        sage: f_rdf = fast_callable(expr, vars=[x], domain=RDF)
        sage: f_rdf(2)
        12.909297426825681
        sage: f_cdf = fast_callable(expr, vars=[x], domain=CDF)
        sage: f_cdf(2)
        12.909297426825681
        sage: f_cdf(2+I)
        10.40311925062204 + 11.510943740958707*I
        sage: f = fast_callable(expr, vars=('z','x','y'))
        sage: f(1, 2, 3)
        sin(2) + 12

        sage: K.<x> = QQ[]
        sage: p = -1/4*x^6 + 1/2*x^5 - x^4 - 12*x^3 + 1/2*x^2 - 1/95*x - 1/2
        sage: fp = fast_callable(p, domain=RDF)
        sage: fp.op_list()
        [('load_arg', 0), ('load_const', -0.25), 'mul', ('load_const', 0.5), 'add',
         ('load_arg', 0), 'mul', ('load_const', -1.0), 'add', ('load_arg', 0), 'mul',
         ('load_const', -12.0), 'add', ('load_arg', 0), 'mul', ('load_const', 0.5),
         'add', ('load_arg', 0), 'mul', ('load_const', -0.010526315789473684),
         'add', ('load_arg', 0), 'mul', ('load_const', -0.5), 'add', 'return']
        sage: fp(3.14159)
        -552.4182988917153

        sage: K.<x,y,z> = QQ[]
        sage: p = x*y^2 + 1/3*y^2 - x*z - y*z
        sage: fp = fast_callable(p, domain=RDF)
        sage: fp.op_list()
        [('load_const', 0.0), ('load_const', 1.0), ('load_arg', 0), ('ipow', 1),
         ('load_arg', 1), ('ipow', 2), 'mul', 'mul', 'add',
         ('load_const', 0.3333333333333333), ('load_arg', 1), ('ipow', 2), 'mul', 'add',
         ('load_const', -1.0), ('load_arg', 0), ('ipow', 1), ('load_arg', 2),
         ('ipow', 1), 'mul', 'mul', 'add', ('load_const', -1.0), ('load_arg', 1),
         ('ipow', 1), ('load_arg', 2), ('ipow', 1), 'mul', 'mul', 'add', 'return']
        sage: fp(e, pi, sqrt(2))   # abs tol 3e-14                                      # needs sage.symbolic
        21.831120464939584
        sage: symbolic_result = p(e, pi, sqrt(2)); symbolic_result                      # needs sage.symbolic
        pi^2*e + 1/3*pi^2 - sqrt(2)*pi - sqrt(2)*e
        sage: n(symbolic_result)                                                        # needs sage.symbolic
        21.8311204649396

    ::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=('x','y'), domain=float)
        sage: x = etb.var('x')
        sage: y = etb.var('y')
        sage: expr = etb.call(sin, x^2 + y); expr
        sin(add(ipow(v_0, 2), v_1))
        sage: fc = fast_callable(expr, domain=float)
        sage: fc(5, 7)
        0.5514266812416906

    Check that :func:`fast_callable` also works for symbolic functions with
    evaluation functions::

        sage: # needs sage.symbolic
        sage: def evalf_func(self, x, y, parent):
        ....:     return parent(x*y) if parent is not None else x*y
        sage: x,y = var('x,y')
        sage: f = function('f', evalf_func=evalf_func)
        sage: fc = fast_callable(f(x, y), vars=[x, y])
        sage: fc(3, 4)
        f(3, 4)

    And also when there are complex values involved::

        sage: # needs sage.symbolic
        sage: def evalf_func(self, x, y, parent):
        ....:     return parent(I*x*y) if parent is not None else I*x*y
        sage: g = function('g', evalf_func=evalf_func)
        sage: fc = fast_callable(g(x, y), vars=[x, y])
        sage: fc(3, 4)
        g(3, 4)
        sage: fc2 = fast_callable(g(x, y), domain=complex, vars=[x, y])
        sage: fc2(3, 4)
        12j
        sage: fc3 = fast_callable(g(x, y), domain=float, vars=[x, y])
        sage: fc3(3, 4)
        Traceback (most recent call last):
            ...
        TypeError: unable to simplify to float approximation
    """
    cdef Expression et
    if isinstance(x, Expression):
        et = x
        vars = et._etb._vars
    else:
        if not vars:
            # fast_float passes empty list/tuple
            vars = None

        if isinstance(x, Expression_abc) and x.is_callable():
            if vars is None:
                vars = x.arguments()
            if expect_one_var and len(vars) != 1:
                raise ValueError(f"passed expect_one_var=True, but the callable expression takes {len(vars)} arguments")
        elif isinstance(x, Expression_abc):
            if vars is None:
                vars = x.variables()
                if expect_one_var and len(vars) <= 1:
                    if not vars:
                        vars = ['EXTRA_VAR0']
                else:
                    raise ValueError("list of variables must be specified for symbolic expressions")

            def to_var(var):
                if isinstance(var, Expression_abc) and var.is_symbol():
                    return var
                from sage.symbolic.ring import SR
                return SR.var(var)
            vars = [to_var(var) for var in vars]
            # Convert to a callable symbolic expression
            x = x.function(*vars)

        if vars is None:
            from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
            from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
            if isinstance(x.parent(), PolynomialRing_general) or isinstance(x.parent(), MPolynomialRing_base):
                vars = x.parent().variable_names()
            else:
                # constant
                vars = ()

        etb = ExpressionTreeBuilder(vars=vars, domain=domain)
        et = x._fast_callable_(etb)

    builder, str = _builder_and_stream(vars=vars, domain=domain)

    generate_code(et, str)
    str.instr('return')
    return builder(str.get_current())


def _builder_and_stream(vars, domain):
    r"""
    Return a builder and a stream.

    This is an internal function used only once, by :func:`fast_callable`.

    INPUT:

    - ``vars`` -- a sequence of variable names

    - ``domain`` -- a Sage parent or Python type or ``None``; if non-``None``,
      all arithmetic is done in that domain

    OUTPUT: a :class:`Wrapper`, an class:`InstructionStream`

    EXAMPLES::

        sage: from sage.ext.fast_callable import _builder_and_stream
        sage: _builder_and_stream(["x", "y"], ZZ)
        (<class '...interpreters.wrapper_el.Wrapper_el'>,
         <sage.ext.fast_callable.InstructionStream object at 0x...>)
        sage: _builder_and_stream(["x", "y"], RR)                                       # needs sage.rings.real_mpfr
        (<class '...interpreters.wrapper_rr.Wrapper_rr'>,
         <sage.ext.fast_callable.InstructionStream object at 0x...>)

    Modularized test with sagemath-categories after :issue:`35095`, which has
    (a basic version of) ``RDF`` but not the specialized interpreter for it.
    In this case, the function falls back to using the :class:`Element`
    interpreter::

        sage: domain = RDF
        sage: from sage.structure.element import Element as domain
        sage: _builder_and_stream(["x", "y"], domain)
        (<class '...interpreters.wrapper_el.Wrapper_el'>,
         <sage.ext.fast_callable.InstructionStream object at 0x...>)
    """
    if isinstance(domain, sage.rings.abc.RealField):
        try:
            from sage.ext.interpreters.wrapper_rr import metadata, Wrapper_rr as builder
        except ImportError:
            pass
        else:
            return builder, InstructionStream(metadata, len(vars), domain)

    if isinstance(domain, sage.rings.abc.ComplexField):
        try:
            from sage.ext.interpreters.wrapper_cc import metadata, Wrapper_cc as builder
        except ImportError:
            pass
        else:
            return builder, InstructionStream(metadata, len(vars), domain)

    if isinstance(domain, sage.rings.abc.RealDoubleField) or domain is float:
        try:
            from sage.ext.interpreters.wrapper_rdf import metadata, Wrapper_rdf as builder
        except ImportError:
            pass
        else:
            return builder, InstructionStream(metadata, len(vars), domain)

    if isinstance(domain, sage.rings.abc.ComplexDoubleField):
        try:
            from sage.ext.interpreters.wrapper_cdf import metadata, Wrapper_cdf as builder
        except ImportError:
            pass
        else:
            return builder, InstructionStream(metadata, len(vars), domain)

    if domain is None:
        from sage.ext.interpreters.wrapper_py import metadata, Wrapper_py as builder
        return builder, InstructionStream(metadata, len(vars))

    from sage.ext.interpreters.wrapper_el import metadata, Wrapper_el as builder
    return builder, InstructionStream(metadata, len(vars), domain)


def function_name(fn):
    r"""
    Given a function, return a string giving a name for the function.

    For functions we recognize, we use our standard opcode name for the
    function (so :func:`operator.add` becomes ``'add'``, and
    :func:`sage.functions.trig.sin` becomes ``'sin'``).

    For functions we don't recognize, we try to come up with a name,
    but the name will be wrapped in braces; this is a signal that
    we'll definitely use a slow Python call to call this function.
    (We may use a slow Python call even for functions we do recognize,
    if we're targeting an interpreter without an opcode for the function.)

    Only used when printing Expressions.

    EXAMPLES::

        sage: from sage.ext.fast_callable import function_name
        sage: function_name(operator.pow)
        'pow'
        sage: function_name(cos)                                                        # needs sage.symbolic
        'cos'
        sage: function_name(factorial)
        '{factorial}'
    """
    from sage.structure.dynamic_class import DynamicMetaclass
    if isinstance(type(fn), DynamicMetaclass):
        return "{%r}" % fn
    builtins = get_builtin_functions()
    if fn in builtins:
        return builtins[fn]
    try:
        return "{%s}" % fn.__name__
    except AttributeError:
        return "{%r}" % fn


cdef class ExpressionTreeBuilder:
    r"""
    A class with helper methods for building Expressions.

    An instance of this class is passed to :meth:`_fast_callable_` methods;
    you can also instantiate it yourself to create your own expressions
    for :func:`fast_callable`, bypassing :meth:`_fast_callable_`.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder('x')
        sage: x = etb.var('x')
        sage: (x+3)*5
        mul(add(v_0, 3), 5)
    """

    cdef readonly object _domain
    cdef readonly object _vars

    def __init__(self, vars, domain=None):
        r"""
        Initialize an instance of :class:`ExpressionTreeBuilder`.

        Takes a list or tuple of variable names to use, and also an optional
        ``domain``.  If a ``domain`` is given, then creating an
        :class:`ExpressionConstant` node with the :meth:`__call__`, make, or
        constant methods will convert the value into the given domain.

        Note that this is the only effect of the domain parameter.  It
        is quite possible to use different domains for
        :class:`ExpressionTreeBuilder` and for :func:`fast_callable`; in that case,
        constants will be converted twice (once when building the
        :class:`Expression`, and once when generating code).

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder('x')
            sage: etb(3^50)
            717897987691852588770249
            sage: etb = ExpressionTreeBuilder('x', domain=RR)
            sage: etb(3^50)
            7.17897987691853e23
        """

        if isinstance(vars, tuple):
            vars = list(vars)
        elif not isinstance(vars, list):
            vars = [vars]

        vars = [self._clean_var(v) for v in vars]

        self._domain = domain
        self._vars = vars

    def __call__(self, x):
        r"""
        Try to convert the given value to an :class:`Expression`.

        If it is already an Expression, just return it.  If it has a
        :meth:`_fast_callable_` method, then call the method with ``self`` as
        an argument.  Otherwise, use ``self.constant()`` to turn it into a constant.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder('x')
            sage: v = etb(3); v, type(v)
            (3, <class 'sage.ext.fast_callable.ExpressionConstant'>)
            sage: v = etb(polygen(QQ)); v, type(v)
            (v_0, <class 'sage.ext.fast_callable.ExpressionVariable'>)
            sage: v is etb(v)
            True
        """
        if isinstance(x, Expression):
            return x

        try:
            fc = x._fast_callable_
        except AttributeError:
            return self.constant(x)

        return fc(self)

    def _clean_var(self, v):
        r"""
        Give a variable name, given a variable.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder('x')
            sage: var('x')                                                              # needs sage.symbolic
            x
            sage: etb._clean_var(x)                                                     # needs sage.symbolic
            'x'
            sage: x = polygen(RR); x
            x
            sage: etb._clean_var(x)
            'x'
        """
        # There should be a better way to do this.  (Maybe there is.)
        if not isinstance(v, str):
            v = str(v)
            if '*' in v:
                v = v[v.index('*')+1:]
        return v

    def constant(self, c):
        r"""
        Turn the argument into an :class:`ExpressionConstant`, converting it to
        our domain if we have one.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder('x')
            sage: etb.constant(pi)                                                      # needs sage.symbolic
            pi
            sage: etb = ExpressionTreeBuilder('x', domain=RealField(200))               # needs sage.rings.real_mpfr
            sage: etb.constant(pi)                                                      # needs sage.rings.real_mpfr sage.symbolic
            3.1415926535897932384626433832795028841971693993751058209749
        """
        if self._domain is not None:
            c = self._domain(c)
        return ExpressionConstant(self, c)

    def var(self, v):
        r"""
        Turn the argument into an :class:`ExpressionVariable`.  Look it up in
        the list of variables.  (Variables are matched by name.)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: var('a,b,some_really_long_name')
            (a, b, some_really_long_name)
            sage: x = polygen(QQ)
            sage: etb = ExpressionTreeBuilder(vars=('a','b',some_really_long_name, x))
            sage: etb.var(some_really_long_name)
            v_2
            sage: etb.var('some_really_long_name')
            v_2
            sage: etb.var(x)
            v_3
            sage: etb.var('y')
            Traceback (most recent call last):
            ...
            ValueError: Variable 'y' not found...
        """
        var_name = self._clean_var(v)
        try:
            ind = self._vars.index(var_name)
        except ValueError:
            raise ValueError(f"Variable '{var_name}' not found in {self._vars}")
        return ExpressionVariable(self, ind)

    def _var_number(self, n):
        r"""
        Given an integer, return the variable with that index.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=('a','b','c','d'))
            sage: etb._var_number(0)
            v_0
            sage: etb._var_number(3)
            v_3
            sage: etb._var_number(4)
            Traceback (most recent call last):
            ...
            ValueError: Variable number 4 out of range
        """
        if 0 <= n < len(self._vars):
            return ExpressionVariable(self, n)
        raise ValueError("Variable number %d out of range" % n)

    def call(self, fn, *args):
        r"""
        Construct a call node, given a function and a list of arguments.

        The arguments will be converted to Expressions using
        :meth:`ExpressionTreeBuilder.__call__`.

        As a special case, notice if the function is :func:`operator.pow` and
        the second argument is integral, and construct an :class:`ExpressionIPow`
        instead.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: etb.call(cos, x)
            cos(v_0)
            sage: etb.call(sin, 1)
            sin(1)
            sage: etb.call(sin, etb(1))
            sin(1)
            sage: etb.call(factorial, x+57)
            {factorial}(add(v_0, 57))
            sage: etb.call(operator.pow, x, 543)
            ipow(v_0, 543)
        """
        if fn is operator.pow:
            base, exponent = args
            return self(base)**exponent
        else:
            return ExpressionCall(self, fn, [self(a) for a in args])

    def choice(self, cond, iftrue, iffalse):
        r"""
        Construct a choice node (a conditional expression), given the
        condition, and the values for the true and false cases.

        (It's possible to create choice nodes, but they don't work yet.)

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: etb.choice(etb.call(operator.eq, x, 0), 0, 1/x)                       # needs sage.symbolic
            (0 if {eq}(v_0, 0) else div(1, v_0))
        """
        return ExpressionChoice(self,
                                cond,
                                self(iftrue),
                                self(iffalse))


# Cache these values, to make expression building a tiny bit faster
# (by skipping the hash-table lookup in the operator module).
cdef op_add = operator.add
cdef op_sub = operator.sub
cdef op_mul = operator.mul
cdef op_div
try:
    op_div = operator.div
except AttributeError:
    op_div = object()  # Unique object not equal to anything else
cdef op_truediv = operator.truediv
cdef op_floordiv = operator.floordiv
cdef op_pow = operator.pow
cdef op_neg = operator.neg
cdef op_abs = operator.abs
cdef op_inv = operator.inv


cdef class Expression:
    r"""
    Represent an expression for :func:`fast_callable`.

    Supports the standard Python arithmetic operators; if arithmetic
    is attempted between an Expression and a non-Expression, the
    non-Expression is converted to an expression (using the
    :meth:`__call__` method of the Expression's :class:`ExpressionTreeBuilder`).

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=(x,))
        sage: x = etb.var(x)
        sage: etb(x)
        v_0
        sage: etb(3)
        3
        sage: etb.call(sin, x)
        sin(v_0)
        sage: (x+1)/(x-1)
        div(add(v_0, 1), sub(v_0, 1))
        sage: x//5
        floordiv(v_0, 5)
        sage: -abs(~x)
        neg(abs(inv(v_0)))
    """

    cdef ExpressionTreeBuilder _etb

    def __init__(self, etb):
        r"""
        Initialize an Expression.  Sets the _etb member.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: v = etb(3); v  # indirect doctest                                     # needs sage.symbolic
            3
            sage: v._get_etb() is etb                                                   # needs sage.symbolic
            True
        """
        self._etb = etb

    def _get_etb(self):
        r"""
        Return the ExpressionTreeBuilder used to build a given expression.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: v = etb(3); v                                                         # needs sage.symbolic
            3
            sage: v._get_etb() is etb                                                   # needs sage.symbolic
            True
        """
        return self._etb

    def __add__(s, o):
        r"""
        Compute a sum of two Expressions.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x+x
            add(v_0, v_0)
            sage: x+1
            add(v_0, 1)
            sage: 1+x
            add(1, v_0)
            sage: x.__add__(1)
            add(v_0, 1)
            sage: x.__radd__(1)
            add(1, v_0)
        """
        return _expression_binop_helper(s, o, op_add)

    def __sub__(s, o):
        r"""
        Compute a difference of two Expressions.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x-x
            sub(v_0, v_0)
            sage: x-1
            sub(v_0, 1)
            sage: 1-x
            sub(1, v_0)
            sage: x.__sub__(1)
            sub(v_0, 1)
            sage: x.__rsub__(1)
            sub(1, v_0)
        """
        return _expression_binop_helper(s, o, op_sub)

    def __mul__(s, o):
        r"""
        Compute a product of two Expressions.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x*x
            mul(v_0, v_0)
            sage: x*1
            mul(v_0, 1)
            sage: 1*x
            mul(1, v_0)
            sage: x.__mul__(1)
            mul(v_0, 1)
            sage: x.__rmul__(1)
            mul(1, v_0)
        """
        return _expression_binop_helper(s, o, op_mul)

    def __truediv__(s, o):
        r"""
        Compute a quotient of two Expressions.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x/x
            div(v_0, v_0)
            sage: x/1
            div(v_0, 1)
            sage: 1/x
            div(1, v_0)
            sage: x.__truediv__(1)
            div(v_0, 1)
            sage: x.__rtruediv__(1)
            div(1, v_0)
        """
        return _expression_binop_helper(s, o, op_truediv)

    def __floordiv__(s, o):
        r"""
        Compute the floordiv (the floor of the quotient) of two Expressions.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x//x
            floordiv(v_0, v_0)
            sage: x//1
            floordiv(v_0, 1)
            sage: 1//x
            floordiv(1, v_0)
            sage: x.__floordiv__(1)
            floordiv(v_0, 1)
            sage: x.__rfloordiv__(1)
            floordiv(1, v_0)
        """
        return _expression_binop_helper(s, o, op_floordiv)

    def __pow__(s, o, dummy):
        r"""
        Compute a power expression from two Expressions.

        If the second Expression is a constant integer, then return
        an :class:`ExpressionIPow` instead of an :class:`ExpressionCall`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x^x
            pow(v_0, v_0)
            sage: x^1
            ipow(v_0, 1)
            sage: x.__pow__(1)
            ipow(v_0, 1)
            sage: x.__pow__(1.0)
            pow(v_0, 1.00000000000000)
            sage: x.__rpow__(1)
            pow(1, v_0)
        """
        # XXX There is a performance regression from the original
        # fast_float here; it would replace small integer powers with
        # multiplication.  We can't do this safely until we support
        # common subexpression elimination (or at least the dup instruction).
        # (Plus, we should consider how strict a semantics we want;
        # probably this sort of optimization should be controlled by a
        # flag.)

        cdef Expression es
        if isinstance(o, (int, Integer)):
            es = s
            return ExpressionIPow(es._etb, s, o)
        else:
            # I really don't like this, but I can't think of a better way
            if isinstance(o, Expression_abc) and o in ZZ:
                es = s
                return ExpressionIPow(es._etb, s, ZZ(o))
            else:
                return _expression_binop_helper(s, o, op_pow)

    def __neg__(self):
        r"""
        Compute the negation of an Expression.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: -x
            neg(v_0)
            sage: x.__neg__()
            neg(v_0)
        """
        return ExpressionCall(self._etb, op_neg, [self])

    def __abs__(self):
        r"""
        Compute the absolute value of an Expression.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: abs(x)
            abs(v_0)
            sage: x.abs()
            abs(v_0)
            sage: x.__abs__()
            abs(v_0)
        """
        return ExpressionCall(self._etb, op_abs, [self])

    def abs(self):
        r"""
        Compute the absolute value of an Expression.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: abs(x)
            abs(v_0)
            sage: x.abs()
            abs(v_0)
            sage: x.__abs__()
            abs(v_0)
        """
        return ExpressionCall(self._etb, op_abs, [self])

    def __invert__(self):
        r"""
        Compute the inverse of an Expression.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: ~x
            inv(v_0)
            sage: x.__invert__()
            inv(v_0)
        """
        return ExpressionCall(self._etb, op_inv, [self])


cdef class ExpressionConstant(Expression):
    r"""
    An Expression that represents an arbitrary constant.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=(x,))                                    # needs sage.symbolic
        sage: type(etb(3))                                                              # needs sage.symbolic
        <class 'sage.ext.fast_callable.ExpressionConstant'>
    """

    cdef object _value

    def __init__(self, etb, c):
        r"""
        Initialize an :class:`ExpressionConstant`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder, ExpressionConstant
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: etb(3)
            3
            sage: v = ExpressionConstant(etb, 3); v
            3
            sage: v._get_etb() is etb
            True
            sage: v.value()
            3
            sage: v.value() == 3
            True
        """
        Expression.__init__(self, etb)
        self._value = c

    def value(self):
        r"""
        Return the constant value of an :class:`ExpressionConstant`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: etb(3).value()                                                        # needs sage.symbolic
            3
        """
        return self._value

    def __repr__(self):
        r"""
        Give a string representing this :class:`ExpressionConstant`.
        (We use the repr of its value.)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: v = etb.constant(pi)
            sage: v
            pi
            sage: repr(v)
            'pi'
            sage: v.__repr__()
            'pi'
        """
        return repr(self._value)

cdef class ExpressionVariable(Expression):
    r"""
    An Expression that represents a variable.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=(x,))                                    # needs sage.symbolic
        sage: type(etb.var(x))                                                          # needs sage.symbolic
        <class 'sage.ext.fast_callable.ExpressionVariable'>
    """
    cdef int _variable_index

    def __init__(self, etb, int n):
        r"""
        Initialize an ExpressionVariable.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder, ExpressionVariable
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: etb(x)
            v_0
            sage: v = ExpressionVariable(etb, 0); v
            v_0
            sage: v._get_etb() is etb
            True
            sage: v.variable_index()
            0
        """
        Expression.__init__(self, etb)
        self._variable_index = n

    def variable_index(self):
        r"""
        Return the variable index of an :class:`ExpressionVariable`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: etb(x).variable_index()                                               # needs sage.symbolic
            0
        """
        return self._variable_index

    def __repr__(self):
        r"""
        Give a string representing this :class:`ExpressionVariable`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: v = etb._var_number(0)
            sage: v
            v_0
            sage: repr(v)
            'v_0'
            sage: v.__repr__()
            'v_0'
        """
        # Should we look up the variable name in self._etb, instead?
        # I think not.. I like the emphasis that we're totally removed
        # from the original expression when we have an Expression.
        return "v_%d" % self._variable_index


cdef class ExpressionCall(Expression):
    r"""
    An Expression that represents a function call.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=(x,))                                    # needs sage.symbolic
        sage: type(etb.call(sin, x))                                                    # needs sage.symbolic
        <class 'sage.ext.fast_callable.ExpressionCall'>
    """
    cdef object _function
    cdef object _arguments

    def __init__(self, etb, fn, args):
        r"""
        Initialize an :class:`ExpressionCall`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder, ExpressionCall
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: etb.call(factorial, x)
            {factorial}(v_0)
            sage: v = ExpressionCall(etb, factorial, [x]); v
            {factorial}(v_0)
            sage: v._get_etb() is etb
            True
            sage: v.function()
            factorial
            sage: v.arguments()
            [v_0]
        """
        Expression.__init__(self, etb)
        self._function = fn
        self._arguments = args

    def function(self):
        r"""
        Return the function from this :class:`ExpressionCall`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: etb.call(sin, x).function()                                           # needs sage.symbolic
            sin
        """
        return self._function

    def arguments(self):
        r"""
        Return the arguments from this :class:`ExpressionCall`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: etb.call(sin, x).arguments()                                          # needs sage.symbolic
            [v_0]
        """
        return copy(self._arguments)

    def __repr__(self):
        r"""
        Give a string representing this :class:`ExpressionCall`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb.var(x)
            sage: etb.call(operator.add, x, 1)
            add(v_0, 1)
            sage: etb.call(factorial, x)
            {factorial}(v_0)
            sage: v = etb.call(sin, x)
            sage: v
            sin(v_0)
            sage: repr(v)
            'sin(v_0)'
            sage: v.__repr__()
            'sin(v_0)'
        """
        fn = function_name(self._function)
        return '%s(%s)' % (fn, ', '.join(repr(a) for a in self._arguments))


cdef class ExpressionIPow(Expression):
    r"""
    A power Expression with an integer exponent.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=(x,))                                    # needs sage.symbolic
        sage: type(etb.var('x')^17)                                                     # needs sage.symbolic
        <class 'sage.ext.fast_callable.ExpressionIPow'>
    """
    cdef object _base
    cdef object _exponent

    def __init__(self, etb, base, exponent):
        r"""
        Initialize an ExpressionIPow.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder, ExpressionIPow
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: x^(-12)
            ipow(v_0, -12)
            sage: v = ExpressionIPow(etb, x, 55); v
            ipow(v_0, 55)
            sage: v._get_etb() is etb
            True
            sage: v.base()
            v_0
            sage: v.exponent()
            55
        """
        Expression.__init__(self, etb)
        self._base = base
        self._exponent = exponent

    def base(self):
        r"""
        Return the base from this :class:`ExpressionIPow`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: (etb(33)^42).base()                                                   # needs sage.symbolic
            33
        """
        return self._base

    def exponent(self):
        r"""
        Return the exponent from this :class:`ExpressionIPow`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: (etb(x)^(-1)).exponent()                                              # needs sage.symbolic
            -1
        """
        return self._exponent

    def __repr__(self):
        r"""
        Give a string representing this :class:`ExpressionIPow`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb.var(x)
            sage: x^3
            ipow(v_0, 3)
            sage: x^(-2)
            ipow(v_0, -2)
            sage: v = (x+1)^3
            sage: v
            ipow(add(v_0, 1), 3)
            sage: repr(v)
            'ipow(add(v_0, 1), 3)'
            sage: v.__repr__()
            'ipow(add(v_0, 1), 3)'
        """
        return 'ipow(%s, %d)' % (repr(self._base), self._exponent)


cdef class ExpressionChoice(Expression):
    r"""
    A conditional expression.

    (It's possible to create choice nodes, but they don't work yet.)

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=(x,))                                    # needs sage.symbolic
        sage: etb.choice(etb.call(operator.eq, x, 0), 0, 1/x)                           # needs sage.symbolic
        (0 if {eq}(v_0, 0) else div(1, v_0))
    """

    cdef object _cond
    cdef object _iftrue
    cdef object _iffalse

    def __init__(self, etb, cond, iftrue, iffalse):
        r"""
        Initialize an :class:`ExpressionChoice`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder, ExpressionChoice
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: etb.choice(x, ~x, 0)
            (inv(v_0) if v_0 else 0)
            sage: v = ExpressionChoice(etb, x, ~x, etb(0)); v
            (inv(v_0) if v_0 else 0)
            sage: v._get_etb() is etb
            True
            sage: v.condition()
            v_0
            sage: v.if_true()
            inv(v_0)
            sage: v.if_false()
            0
        """
        Expression.__init__(self, etb)
        self._cond = cond
        self._iftrue = iftrue
        self._iffalse = iffalse

    def condition(self):
        r"""
        Return the condition of an :class:`ExpressionChoice`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: x = etb(x)                                                            # needs sage.symbolic
            sage: etb.choice(x, ~x, 0).condition()                                      # needs sage.symbolic
            v_0
        """
        return self._cond

    def if_true(self):
        r"""
        Return the true branch of an :class:`ExpressionChoice`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: x = etb(x)                                                            # needs sage.symbolic
            sage: etb.choice(x, ~x, 0).if_true()                                        # needs sage.symbolic
            inv(v_0)
        """
        return self._iftrue

    def if_false(self):
        r"""
        Return the false branch of an :class:`ExpressionChoice`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))                                # needs sage.symbolic
            sage: x = etb(x)                                                            # needs sage.symbolic
            sage: etb.choice(x, ~x, 0).if_false()                                       # needs sage.symbolic
            0
        """
        return self._iffalse

    def __repr__(self):
        r"""
        Give a string representation for this :class:`ExpressionChoice`.
        (Based on the syntax for Python conditional expressions.)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=(x,))
            sage: x = etb(x)
            sage: v = etb.choice(x, ~x, 0)
            sage: v
            (inv(v_0) if v_0 else 0)
            sage: repr(v)
            '(inv(v_0) if v_0 else 0)'
            sage: v.__repr__()
            '(inv(v_0) if v_0 else 0)'
        """
        return '(%s if %s else %s)' % (repr(self._iftrue),
                                       repr(self._cond),
                                       repr(self._iffalse))


cpdef _expression_binop_helper(s, o, op):
    r"""
    Make an :class:`Expression` for (``s`` ``op`` ``o``).

    Either ``s`` or ``o`` (or both) must already be an :class:`Expression`.

    EXAMPLES::

        sage: from sage.ext.fast_callable import _expression_binop_helper, ExpressionTreeBuilder
        sage: var('x,y')                                                                # needs sage.symbolic
        (x, y)
        sage: etb = ExpressionTreeBuilder(vars=(x,y))                                   # needs sage.symbolic
        sage: x = etb(x)                                                                # needs sage.symbolic

    Now ``x`` is an :class:`Expression`, but ``y`` is not.  Still, all the following
    cases work::

        sage: _expression_binop_helper(x, x, operator.add)                              # needs sage.symbolic
        add(v_0, v_0)
        sage: _expression_binop_helper(x, y, operator.add)                              # needs sage.symbolic
        add(v_0, v_1)
        sage: _expression_binop_helper(y, x, operator.add)                              # needs sage.symbolic
        add(v_1, v_0)
    """
    # The Cython way of handling operator overloading on cdef classes
    # (which is inherited from Python) is quite annoying.  Inside the
    # code for a binary operator, you know that either the first or
    # second argument (or both) is a member of your class, but you
    # don't know which.

    # If there is an arithmetic operator between an Expression and
    # a non-Expression, I want to convert the non-Expression into
    # an Expression.  But to do that, I need the ExpressionTreeBuilder
    # from the Expression.

    cdef Expression self
    cdef Expression other

    if not isinstance(o, Expression):
        self = s
        other = self._etb(o)
    elif not isinstance(s, Expression):
        other = o
        self = other._etb(s)
    else:
        self = s
        other = o
        assert self._etb is other._etb

    return ExpressionCall(self._etb, op, [self, other])


class IntegerPowerFunction():
    r"""
    This class represents the function `x^n` for an arbitrary integral power `n`.

    That is, ``IntegerPowerFunction(2)`` is the squaring function;
    ``IntegerPowerFunction(-1)`` is the reciprocal function.

    EXAMPLES::

        sage: from sage.ext.fast_callable import IntegerPowerFunction
        sage: square = IntegerPowerFunction(2)
        sage: square
        (^2)
        sage: square(pi)                                                                # needs sage.symbolic
        pi^2
        sage: square(I)                                                                 # needs sage.symbolic
        -1
        sage: square(RIF(-1, 1)).str(style='brackets')                                  # needs sage.rings.real_interval_field
        '[0.0000000000000000 .. 1.0000000000000000]'
        sage: IntegerPowerFunction(-1)
        (^(-1))
        sage: IntegerPowerFunction(-1)(22/7)
        7/22
        sage: v = Integers(123456789)(54321)
        sage: v^9876543210
        79745229
        sage: IntegerPowerFunction(9876543210)(v)
        79745229
    """

    def __init__(self, n):
        r"""
        Initialize an :class:`IntegerPowerFunction`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import IntegerPowerFunction
            sage: cube = IntegerPowerFunction(3)
            sage: cube
            (^3)
            sage: cube(AA(7)^(1/3))                                                     # needs sage.rings.number_field
            7.000000000000000?
            sage: cube.exponent
            3
        """
        self.exponent = n

    def __repr__(self):
        r"""
        Return a string representing this :class:`IntegerPowerFunction`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import IntegerPowerFunction
            sage: square = IntegerPowerFunction(2)
            sage: square
            (^2)
            sage: repr(square)
            '(^2)'
            sage: square.__repr__()
            '(^2)'
            sage: repr(IntegerPowerFunction(-57))
            '(^(-57))'
        """
        if self.exponent >= 0:
            return "(^%s)" % self.exponent
        else:
            return "(^(%s))" % self.exponent

    def __call__(self, x):
        r"""
        Call this :class:`IntegerPowerFunction`, to compute a power of its argument.

        EXAMPLES::

            sage: from sage.ext.fast_callable import IntegerPowerFunction
            sage: square = IntegerPowerFunction(2)
            sage: square.__call__(5)
            25
            sage: square(5)
            25
        """
        return x**self.exponent


cdef dict builtin_functions = None
cpdef dict get_builtin_functions():
    r"""
    Return a dictionary from Sage and Python functions to opcode names.

    The result is cached.

    The dictionary is used in :class:`ExpressionCall`.

    EXAMPLES::

        sage: from sage.ext.fast_callable import get_builtin_functions
        sage: builtins = get_builtin_functions()
        sage: sorted(set(builtins.values()))                                            # needs sage.symbolic
        ['abs', 'acos', 'acosh', 'add', 'asin', 'asinh', 'atan', 'atanh', 'ceil',
         'cos', 'cosh', 'cot', 'csc', 'div', 'exp', 'floor', 'floordiv', 'inv', 'log',
         'mul', 'neg', 'pow', 'sec', 'sin', 'sinh', 'sqrt', 'sub', 'tan', 'tanh']
        sage: builtins[sin]                                                             # needs sage.symbolic
        'sin'
        sage: builtins[ln]
        'log'
    """
    # We delay building builtin_functions to break a circular import
    # between sage.functions and this file.
    global builtin_functions
    if builtin_functions is not None:
        return builtin_functions
    builtin_functions = {
        op_add: 'add',
        op_sub: 'sub',
        op_mul: 'mul',
        op_div: 'div',
        op_truediv: 'div',
        op_floordiv: 'floordiv',
        op_abs: 'abs',
        op_neg: 'neg',
        op_inv: 'inv',
        op_pow: 'pow',
    }

    # not handled: atan2, log2, log10
    import sage.functions.all as func_all
    for fn in ('sqrt', 'ceil', 'floor',
               'sin', 'cos', 'tan', 'sec', 'csc', 'cot',
               'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh',
               'asinh', 'acosh', 'atanh', 'exp', 'log'):
        builtin_functions[getattr(func_all, fn)] = fn
    builtin_functions[func_all.abs_symbolic] = 'abs'
    builtin_functions[func_all.ln] = 'log'
    return builtin_functions


cdef class InstructionStream  # forward declaration


cpdef generate_code(Expression expr, InstructionStream stream):
    r"""
    Generate code from an :class:`Expression` tree; write the result into an
    :class:`InstructionStream`.

    In :func:`fast_callable`, first we create an :class:`Expression`, either directly
    with an :class:`ExpressionTreeBuilder` or with :meth:`_fast_callable_` methods.
    Then we optimize the :class:`Expression` in tree form.  (Unfortunately,
    this step is currently missing -- we do no optimizations.)

    Then we linearize the :class:`Expression` into a sequence of instructions,
    by walking the :class:`Expression` and sending the corresponding stack
    instructions to an :class:`InstructionStream`.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder, generate_code, InstructionStream

        sage: # needs sage.symbolic
        sage: etb = ExpressionTreeBuilder('x')
        sage: x = etb.var('x')
        sage: expr = (x+pi) * (x+1)
        sage: from sage.ext.interpreters.wrapper_py import metadata, Wrapper_py
        sage: instr_stream = InstructionStream(metadata, 1)
        sage: generate_code(expr, instr_stream)
        sage: instr_stream.instr('return')
        sage: v = Wrapper_py(instr_stream.get_current())
        sage: type(v)
        <class '...interpreters.wrapper_py.Wrapper_py'>
        sage: v(7)
        8*pi + 56

    TESTS::

        sage: def my_sin(x): return sin(x)
        sage: def my_norm(x, y): return x*x + y*y
        sage: def my_sqrt(x):
        ....:     if x < 0: raise ValueError("sqrt of negative number")
        ....:     return sqrt(x, extend=False)

        sage: # needs sage.symbolic
        sage: fc = fast_callable(expr, domain=RealField(130))
        sage: fc(0)
        3.1415926535897932384626433832795028842
        sage: fc(1)
        8.2831853071795864769252867665590057684
        sage: fc = fast_callable(expr, domain=RDF)
        sage: fc(0)
        3.141592653589793
        sage: fc(1)
        8.283185307179586
        sage: fc.op_list()
        [('load_arg', 0), ('load_const', pi), 'add', ('load_arg', 0), ('load_const', 1), 'add', 'mul', 'return']
        sage: fc = fast_callable(etb.call(sin, x) + etb.call(sqrt, x), domain=RDF)
        sage: fc(1)
        1.8414709848078965
        sage: fc.op_list()
        [('load_arg', 0), 'sin', ('load_arg', 0), 'sqrt', 'add', 'return']
        sage: fc = fast_callable(etb.call(sin, x) + etb.call(sqrt, x))
        sage: fc(1)
        sin(1) + 1
        sage: fc.op_list()
        [('load_arg', 0), ('py_call', sin, 1), ('load_arg', 0), ('py_call', <function sqrt at ...>, 1), 'add', 'return']
        sage: fc = fast_callable(etb.call(my_sin, x), domain=RDF)
        sage: fc(3)
        0.1411200080598672

        sage: # needs sage.rings.real_mpfr sage.symbolic
        sage: fc = fast_callable(etb.call(my_sin, x), domain=RealField(100))
        sage: fc(3)
        0.14112000805986722210074480281
        sage: fc.op_list()
        [('load_arg', 0), ('py_call', <function my_sin at 0x...>, 1), 'return']

        sage: # needs sage.symbolic
        sage: fc = fast_callable(etb.call(my_sqrt, x), domain=RDF)
        sage: fc(3)
        1.7320508075688772
        sage: parent(fc(3))
        Real Double Field
        sage: fc(-3)
        Traceback (most recent call last):
        ...
        ValueError: sqrt of negative number
        sage: fc = fast_callable(etb.call(my_sqrt, x), domain=RR)
        sage: fc(3)
        1.73205080756888
        sage: fc(-3)
        Traceback (most recent call last):
        ...
        ValueError: sqrt of negative number

        sage: etb2 = ExpressionTreeBuilder(('y','z'))
        sage: y = etb2.var('y')
        sage: z = etb2.var('z')
        sage: fc = fast_callable(etb2.call(sqrt, etb2.call(my_norm, y, z)), domain=RDF)
        sage: fc(3, 4)
        5.0
        sage: fc.op_list()
        [('load_arg', 0), ('load_arg', 1), ('py_call', <function my_norm at 0x...>, 2), 'sqrt', 'return']
        sage: fc.python_calls()
        [<function my_norm at 0x...>]
        sage: fc = fast_callable(etb2.call(sqrt, etb2.call(my_norm, y, z)), domain=RR)  # needs sage.rings.real_mpfr
        sage: fc(3, 4)                                                                  # needs sage.rings.real_mpfr
        5.00000000000000
        sage: fc = fast_callable(etb2.call(my_norm, y, z), domain=ZZ)
        sage: fc(3, 4)
        25
        sage: fc.op_list()
        [('load_arg', 0), ('load_arg', 1), ('py_call', <function my_norm at 0x...>, 2), 'return']

        sage: # needs sage.symbolic
        sage: fc = fast_callable(expr)
        sage: fc(3.0r)
        4.0*pi + 12.0
        sage: fc = fast_callable(x+3, domain=ZZ)
        sage: fc(4)
        7
        sage: fc = fast_callable(x/3, domain=ZZ)
        sage: fc(4)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
        sage: fc(6)
        2
        sage: fc = fast_callable(etb.call(sin, x), domain=ZZ)
        sage: fc(0)
        0
        sage: fc(3)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert sin(3) to an integer

    ::

        sage: # needs sage.symbolic
        sage: fc = fast_callable(etb(x)^100)
        sage: fc(pi)
        pi^100
        sage: fc = fast_callable(etb(x)^100, domain=ZZ)
        sage: fc(2)
        1267650600228229401496703205376
        sage: fc = fast_callable(etb(x)^100, domain=RIF)
        sage: fc(RIF(-2))
        1.2676506002282295?e30
        sage: fc = fast_callable(etb(x)^100, domain=RDF)
        sage: fc.op_list()
        [('load_arg', 0), ('ipow', 100), 'return']
        sage: fc(1.1)
        13780.61233982...
        sage: fc = fast_callable(etb(x)^100, domain=RR)                                 # needs sage.rings.real_mpfr
        sage: fc.op_list()                                                              # needs sage.rings.real_mpfr
        [('load_arg', 0), ('ipow', 100), 'return']
        sage: fc(1.1)                                                                   # needs sage.rings.real_mpfr
        13780.6123398224
        sage: fc = fast_callable(etb(x)^(-100), domain=RDF)
        sage: fc.op_list()
        [('load_arg', 0), ('ipow', -100), 'return']
        sage: fc(1.1)
        7.25657159014...e-05
        sage: fc = fast_callable(etb(x)^(-100), domain=RR)                              # needs sage.rings.real_mpfr
        sage: fc(1.1)                                                                   # needs sage.rings.real_mpfr
        0.0000725657159014814
        sage: expo = 2^32
        sage: base = (1.0).nextabove()
        sage: fc = fast_callable(etb(x)^expo, domain=RDF)
        sage: fc.op_list()
        [('load_arg', 0), ('py_call', (^4294967296), 1), 'return']
        sage: fc(base)        # rel tol 1e-15
        1.0000009536747712
        sage: RDF(base)^expo
        1.0000009536747712
        sage: fc = fast_callable(etb(x)^expo, domain=RR)                                # needs sage.rings.real_mpfr
        sage: fc.op_list()                                                              # needs sage.rings.real_mpfr
        [('load_arg', 0), ('py_call', (^4294967296), 1), 'return']
        sage: fc(base)                                                                  # needs sage.rings.real_mpfr
        1.00000095367477
        sage: base^expo
        1.00000095367477

    Make sure we do not overflow the stack with highly nested expressions
    (:issue:`11766`)::

        sage: # needs sage.rings.real_mpfr
        sage: R.<x> = CC[]
        sage: f = R(list(range(100000)))
        sage: ff = fast_callable(f)
        sage: f(0.5)
        2.00000000000000
        sage: ff(0.5)
        2.00000000000000
        sage: f(0.9), ff(0.9)
        (90.0000000000000, 90.0000000000000)
    """
    cdef ExpressionConstant econst
    cdef ExpressionVariable evar
    cdef ExpressionCall ecall
    cdef ExpressionChoice echoice

    # Maintain our own stack to avoid crashing on deeply-nested expressions.
    cdef list todo = [expr]
    do_call = Expression(None)
    while len(todo):
        expr = todo.pop()
        if isinstance(expr, ExpressionConstant):
            econst = expr
            stream.load_const(econst._value)
        elif isinstance(expr, ExpressionVariable):
            evar = expr
            stream.load_arg(evar._variable_index)
        elif isinstance(expr, ExpressionCall):
            ecall = expr
            todo.append(expr)
            todo.append(do_call)
            for arg in reversed(ecall._arguments):
                todo.append(arg)
            continue
        elif expr is do_call:
            # arguments already evaluated, make the call
            ecall = todo.pop()
            fn = ecall._function
            opname = get_builtin_functions().get(fn)
            if opname is not None:
                if stream.has_instr(opname):
                    stream.instr0(opname, ())
                    continue
            if stream.has_instr('py_call'):
                stream.instr('py_call', fn, len(ecall._arguments))
            else:
                raise ValueError("Unhandled function %s in generate_code" % fn)
        elif isinstance(expr, ExpressionIPow):
            base = expr.base()
            exponent = expr.exponent()
            metadata = stream.get_metadata()
            ipow_range = metadata.ipow_range
            if ipow_range is True:
                use_ipow = True
            elif isinstance(ipow_range, tuple):
                a,b = ipow_range
                use_ipow = (a <= exponent <= b)
            else:
                use_ipow = False
            generate_code(base, stream)
            if use_ipow:
                stream.instr('ipow', exponent)
            else:
                stream.instr('py_call', IntegerPowerFunction(exponent), 1)
        else:
            raise ValueError("Unhandled expression kind %s in generate_code" % type(expr))


cdef class InterpreterMetadata  # forward declaration


cdef class InstructionStream:
    r"""
    An :class:`InstructionStream` takes a sequence of instructions (passed in by
    a series of method calls) and computes the data structures needed
    by the interpreter.  This is the stage where we switch from operating
    on :class:`Expression` trees to a linear representation.  If we had a peephole
    optimizer (we don't) it would go here.

    Currently, this class is not very general; it only works for
    interpreters with a fixed set of memory chunks (with fixed names).
    Basically, it only works for stack-based expression interpreters.
    It should be generalized, so that the interpreter metadata includes
    a description of the memory chunks involved and the instruction stream
    can handle any interpreter.

    Once you're done adding instructions, you call :meth:`get_current` to retrieve
    the information needed by the interpreter (as a Python dictionary).
    """

    cdef InterpreterMetadata _metadata
    cdef list _instrs
    cdef list _bytecode
    cdef list _constants
    cdef object _constant_locs
    cdef object _py_constants
    cdef object _py_constant_locs
    cdef int _stack_cur_size
    cdef int _stack_max_size
    cdef int _n_args
    cdef object _domain

    def __init__(self, metadata, n_args, domain=None):
        r"""
        Initialize an :class:`InstructionStream`.

        INPUT:

        - ``metadata`` -- the ``metadata_by_opname`` from a wrapper module

        - ``n_args`` -- the number of arguments accessible by the generated code
          (this is just passed to the wrapper class)

        - ``domain`` -- the domain of interpretation (this is just passed to the
          wrapper class)

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: instr_stream.get_current()
            {'args': 1,
             'code': [],
             'constants': [],
             'domain': None,
             'py_constants': [],
             'stack': 0}
            sage: md = instr_stream.get_metadata()
            sage: type(md)
            <class 'sage.ext.fast_callable.InterpreterMetadata'>
            sage: md.by_opname['py_call']
            (CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs']), 3)
            sage: md.by_opcode[3]
            ('py_call', CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs']))
        """
        self._metadata = metadata
        self._instrs = []
        self._bytecode = []
        self._constants = []
        self._constant_locs = {}
        self._py_constants = []
        self._py_constant_locs = {}
        self._stack_cur_size = 0
        self._stack_max_size = 0
        self._domain = domain
        self._n_args = n_args

    def load_const(self, c):
        r"""
        Add a ``'load_const'`` instruction to this :class:`InstructionStream`.

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream, op_list
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: instr_stream.load_const(5)
            sage: instr_stream.current_op_list()
            [('load_const', 5)]
            sage: instr_stream.load_const(7)
            sage: instr_stream.load_const(5)
            sage: instr_stream.current_op_list()
            [('load_const', 5), ('load_const', 7), ('load_const', 5)]

        Note that constants are shared: even though we load 5 twice, it
        only appears once in the constant table. ::

            sage: instr_stream.get_current()['constants']
            [5, 7]
        """
        self.instr('load_const', c)

    def load_arg(self, n):
        r"""
        Add a ``'load_arg'`` instruction to this :class:`InstructionStream`.

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 12)
            sage: instr_stream.load_arg(5)
            sage: instr_stream.current_op_list()
            [('load_arg', 5)]
            sage: instr_stream.load_arg(3)
            sage: instr_stream.current_op_list()
            [('load_arg', 5), ('load_arg', 3)]
        """
        self.instr('load_arg', n)

    cpdef bint has_instr(self, opname) noexcept:
        r"""
        Check whether this :class:`InstructionStream` knows how to generate code
        for a given instruction.

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: instr_stream.has_instr('return')
            True
            sage: instr_stream.has_instr('factorial')
            False
            sage: instr_stream.has_instr('abs')
            True
        """
        return (opname in self._metadata.by_opname)

    def instr(self, opname, *args):
        r"""
        Generate code in this :class:`InstructionStream` for the given instruction
        and arguments.

        The opname is used to look up a :class:`CompilerInstrSpec`; the
        :class:`CompilerInstrSpec` describes how to interpret the arguments.
        (This is documented in the class docstring for :class:`CompilerInstrSpec`.)

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: instr_stream.instr('load_arg', 0)
            sage: instr_stream.instr('sin')                                             # needs sage.symbolic
            sage: instr_stream.instr('py_call', math.sin, 1)
            sage: instr_stream.instr('abs')
            sage: instr_stream.instr('factorial')
            Traceback (most recent call last):
            ...
            KeyError: 'factorial'
            sage: instr_stream.instr('return')
            sage: instr_stream.current_op_list()                                        # needs sage.symbolic
            [('load_arg', 0), 'sin', ('py_call', <built-in function sin>, 1), 'abs', 'return']
        """
        self.instr0(opname, args)

    cdef instr0(self, opname, tuple args):
        """
        Cdef version of instr. (Can't cpdef because of star args.)
        """
        cdef int i

        spec, opcode = self._metadata.by_opname[opname]
        assert len(spec.parameters) == len(args)

        cdef int n_inputs = spec.n_inputs
        cdef int n_outputs = spec.n_outputs

        self._bytecode.append(opcode)
        for i in range(len(args)):
            if spec.parameters[i] == 'constants':
                # XXX bad for strict-mode floating-point constants
                # (doesn't handle signed 0, NaN)
                arg = args[i]
                if (arg,parent(arg)) in self._constant_locs:
                    self._bytecode.append(self._constant_locs[(arg,parent(arg))])
                else:
                    loc = len(self._constants)
                    self._constants.append(arg)
                    self._constant_locs[(arg,parent(arg))] = loc
                    self._bytecode.append(loc)
            elif spec.parameters[i] == 'args':
                self._bytecode.append(args[i])
            elif spec.parameters[i] == 'code':
                self._bytecode.append(args[i])
            elif spec.parameters[i] == 'n_inputs':
                self._bytecode.append(args[i])
                n_inputs = args[i]
            elif spec.parameters[i] == 'n_outputs':
                self._bytecode.append(args[i])
                n_outputs = args[i]
            elif spec.parameters[i] == 'py_constants':
                arg = args[i]
                if arg in self._py_constant_locs:
                    self._bytecode.append(self._py_constant_locs[arg])
                else:
                    loc = len(self._py_constants)
                    self._py_constants.append(arg)
                    self._py_constant_locs[arg] = loc
                    self._bytecode.append(loc)
            else:
                raise ValueError

        self._stack_cur_size -= n_inputs
        self._stack_cur_size += n_outputs
        self._stack_max_size = max(self._stack_max_size, self._stack_cur_size)

    def get_metadata(self):
        r"""
        Return the interpreter metadata being used by the current :class:`InstructionStream`.

        The code generator sometimes uses this to decide which code
        to generate.

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: md = instr_stream.get_metadata()
            sage: type(md)
            <class 'sage.ext.fast_callable.InterpreterMetadata'>
        """
        return self._metadata

    def current_op_list(self):
        r"""
        Return the list of instructions that have been added to this
        :class:`InstructionStream` so far.

        It's OK to call this, then add more instructions.

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: instr_stream.instr('load_arg', 0)
            sage: instr_stream.instr('py_call', math.sin, 1)
            sage: instr_stream.instr('abs')
            sage: instr_stream.instr('return')
            sage: instr_stream.current_op_list()
            [('load_arg', 0), ('py_call', <built-in function sin>, 1), 'abs', 'return']
        """
        return op_list(self.get_current(), self._metadata)

    def get_current(self):
        r"""
        Return the current state of the :class:`InstructionStream`, as a dictionary
        suitable for passing to a wrapper class.

        NOTE: The dictionary includes internal data structures of the
        :class:`InstructionStream`; you must not modify it.

        EXAMPLES::

            sage: from sage.ext.interpreters.wrapper_el import metadata
            sage: from sage.ext.interpreters.wrapper_rdf import metadata                # needs sage.modules
            sage: from sage.ext.fast_callable import InstructionStream
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: instr_stream.get_current()
            {'args': 1,
             'code': [],
             'constants': [],
             'domain': None,
             'py_constants': [],
             'stack': 0}
            sage: instr_stream.instr('load_arg', 0)
            sage: instr_stream.instr('py_call', math.sin, 1)
            sage: instr_stream.instr('abs')
            sage: instr_stream.instr('return')
            sage: instr_stream.current_op_list()
            [('load_arg', 0), ('py_call', <built-in function sin>, 1), 'abs', 'return']
            sage: instr_stream.get_current()
            {'args': 1,
             'code': [0, 0, 3, 0, 1, 12, 2],
             'constants': [],
             'domain': None,
             'py_constants': [<built-in function sin>],
             'stack': 1}
        """
        d = {'args': self._n_args,
             'constants': self._constants,
             'py_constants': self._py_constants,
             'stack': self._stack_max_size,
             'code': self._bytecode,
             'domain': self._domain}
        return d


cdef class InterpreterMetadata():
    r"""
    The interpreter metadata for a :func:`fast_callable` interpreter.

    Currently consists of a dictionary mapping instruction names to
    (:class:`CompilerInstrSpec`, opcode) pairs, a list mapping opcodes to
    (instruction name, :class:`CompilerInstrSpec`) pairs, and a range of exponents
    for which the ``'ipow'`` instruction can be used.  This range can be
    ``False`` (if the ``'ipow'`` instruction should never be used), a pair of
    two integers `(a, b)`, if ``'ipow'`` should be used for `a \le n \le b`, or
    ``True``, if ``'ipow'`` should always be used.  When ``'ipow'`` cannot be
    used, then we fall back on calling :class:`IntegerPowerFunction`.

    See the class docstring for :class:`CompilerInstrSpec` for more information.

    NOTE: You must not modify the metadata.
    """
    cdef public dict by_opname
    cdef public list by_opcode
    cdef public ipow_range

    def __init__(self, by_opname, by_opcode, ipow_range):
        r"""
        Initialize an InterpreterMetadata object.

        EXAMPLES::

            sage: from sage.ext.fast_callable import InterpreterMetadata
            sage: metadata = InterpreterMetadata(by_opname={'opname dict goes here': True},
            ....:                                by_opcode=['opcode list goes here'],
            ....:                                ipow_range=(2, 57))
            sage: metadata.by_opname
            {'opname dict goes here': True}
            sage: metadata.by_opcode
            ['opcode list goes here']
            sage: metadata.ipow_range
            (2, 57)
        """
        self.by_opname = by_opname
        self.by_opcode = by_opcode
        self.ipow_range = ipow_range


class CompilerInstrSpec():
    r"""
    Describe a single instruction to the :func:`fast_callable` code generator.

    An instruction has a number of stack inputs, a number of stack
    outputs, and a parameter list describing extra arguments that
    must be passed to the :meth:`InstructionStream.instr` method (that end up
    as extra words in the code).

    The parameter list is a list of strings.  Each string is one of
    the following:

    - ``'args'`` -- the instruction argument refers to an input argument of the
      wrapper class; it is just appended to the code

    - ``'constants'``, ``'py_constants'`` -- the instruction argument is a value; the
      value is added to the corresponding list (if it's not already there) and
      the index is appended to the code.

    - ``'n_inputs'``, ``'n_outputs'`` -- the instruction actually takes a variable
      number of inputs or outputs (the ``n_inputs`` and ``n_outputs`` attributes of
      this instruction are ignored). The instruction argument specifies the
      number of inputs or outputs (respectively); it is just appended to the
      code.
    """

    def __init__(self, n_inputs, n_outputs, parameters):
        r"""
        Initialize a :class:`CompilerInstrSpec`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import CompilerInstrSpec
            sage: CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])
            CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])
        """
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
        self.parameters = parameters

    def __repr__(self):
        r"""
        Give a string representation for this :class:`CompilerInstrSpec`.

        EXAMPLES::

            sage: from sage.ext.fast_callable import CompilerInstrSpec
            sage: v = CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])
            sage: v
            CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])
            sage: repr(v)
            "CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])"
            sage: v.__repr__()
            "CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])"
        """
        return "CompilerInstrSpec(%d, %d, %s)" % (self.n_inputs, self.n_outputs, self.parameters)


def op_list(args, metadata):
    r"""
    Given a dictionary with the result of calling :meth:`get_current` on an
    :class:`InstructionStream`, and the corresponding interpreter metadata,
    return a list of the instructions, in a simple somewhat
    human-readable format.

    For debugging only.  (That is, it's probably not a good idea to
    try to programmatically manipulate the result of this function;
    the expected use is just to print the returned list to the
    screen.)

    There's probably no reason to call this directly; if you
    have a wrapper object, call :func:`op_list` on it; if you have an
    :class:`InstructionStream` object, call :meth:`current_op_list` on it.

    EXAMPLES::

        sage: from sage.ext.interpreters.wrapper_el import metadata
        sage: from sage.ext.interpreters.wrapper_rdf import metadata                    # needs sage.modules
        sage: from sage.ext.fast_callable import InstructionStream, op_list
        sage: instr_stream = InstructionStream(metadata, 1)
        sage: instr_stream.instr('load_arg', 0)
        sage: instr_stream.instr('abs')
        sage: instr_stream.instr('return')
        sage: instr_stream.current_op_list()
        [('load_arg', 0), 'abs', 'return']
        sage: op_list(instr_stream.get_current(), metadata)
        [('load_arg', 0), 'abs', 'return']
    """
    ops = []
    code = args['code']
    while len(code):
        opcode = code[0]
        code = code[1:]
        (opname, instr) = metadata.by_opcode[opcode]
        if len(instr.parameters):
            op = [opname]
            for p in instr.parameters:
                p_loc = code[0]
                code = code[1:]
                if p in ('args', 'code', 'n_inputs', 'n_outputs'):
                    op.append(p_loc)
                else:
                    op.append(args[p][p_loc])
            ops.append(tuple(op))
        else:
            ops.append(opname)
    return ops


cdef class Wrapper:
    r"""
    The parent class for all :func:`fast_callable` wrappers.

    Implements shared behavior (currently only debugging).
    """

    def __init__(self, args, metadata):
        r"""
        Initialize a Wrapper object.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder, generate_code, InstructionStream
            sage: etb = ExpressionTreeBuilder('x')
            sage: x = etb.var('x')
            sage: expr = (x+pi) * (x+1)
            sage: from sage.ext.interpreters.wrapper_py import metadata, Wrapper_py
            sage: instr_stream = InstructionStream(metadata, 1)
            sage: generate_code(expr, instr_stream)
            sage: instr_stream.instr('return')
            sage: v = Wrapper_py(instr_stream.get_current())
            sage: v.get_orig_args()
            {'args': 1,
             'code': [0, 0, 1, 0, 4, 0, 0, 1, 1, 4, 6, 2],
             'constants': [pi, 1],
             'domain': None,
             'py_constants': [],
             'stack': 3}
            sage: v.op_list()
            [('load_arg', 0), ('load_const', pi), 'add', ('load_arg', 0), ('load_const', 1), 'add', 'mul', 'return']
        """

        # We only keep the original arguments for debugging (op_list(), etc.);
        # is it worth the memory cost?  (Note that we may be holding on to
        # large objects that could otherwise be garbage collected, for
        # instance.)
        self._orig_args = args
        self._metadata = metadata

    def get_orig_args(self):
        r"""
        Get the original arguments used when initializing this wrapper.

        (Probably only useful when writing doctests.)

        EXAMPLES::

            sage: fast_callable(sin(x)/x, vars=[x], domain=RDF).get_orig_args()         # needs sage.symbolic
            {'args': 1,
             'code': [0, 0, 16, 0, 0, 8, 2],
             'constants': [],
             'domain': Real Double Field,
             'py_constants': [],
             'stack': 2}
        """
        return self._orig_args

    def op_list(self):
        r"""
        Return the list of instructions in this wrapper.

        EXAMPLES::

            sage: fast_callable(cos(x) * x, vars=[x], domain=RDF).op_list()             # needs sage.symbolic
            [('load_arg', 0), ('load_arg', 0), 'cos', 'mul', 'return']
        """
        return op_list(self._orig_args, self._metadata)

    def python_calls(self):
        r"""
        List the Python functions that are called in this wrapper.

        (Python function calls are slow, so ideally this list would
        be empty.  If it is not empty, then perhaps there is an
        optimization opportunity where a Sage developer could speed
        this up by adding a new instruction to the interpreter.)

        EXAMPLES::

            sage: fast_callable(abs(sin(x)), vars=[x], domain=RDF).python_calls()       # needs sage.symbolic
            []
            sage: fast_callable(abs(sin(factorial(x))), vars=[x]).python_calls()        # needs sage.symbolic
            [factorial, sin]
        """
        ops = self.op_list()
        py_calls = []
        for op in ops:
            if isinstance(op, tuple) and op[0] == 'py_call':
                py_calls.append(op[1])
        return py_calls


class FastCallableFloatWrapper:
    r"""
    A class to alter the return types of the fast-callable functions.

    When applying numerical routines (including plotting) to symbolic
    expressions and functions, we generally first convert them to a
    faster form with :func:`fast_callable`.  That function takes a
    ``domain`` parameter that forces the end (and all intermediate)
    results of evaluation to a specific type.  Though usually always
    want the end result to be of type :class:`float`, correctly choosing
    the ``domain`` presents some problems:

      * :class:`float` is a bad choice because it's common for real
        functions to have complex terms in them. Moreover precision
        issues can produce terms like ``1.0 + 1e-12*I`` that are hard
        to avoid if calling ``real()`` on everything is infeasible.

      * :class:`complex` has essentially the same problem as :class:`float`.
        There are several symbolic functions like :func:`min_symbolic`,
        :func:`max_symbolic`, and :func:`floor` that are unable to
        operate on complex numbers.

      * ``None`` leaves the types of the inputs/outputs alone, but due
        to the lack of a specialized interpreter, slows down evaluation
        by an unacceptable amount.

      * ``CDF`` has none of the other issues, because ``CDF`` has its
        own specialized interpreter, a lexicographic ordering (for
        min/max), and supports :func:`floor`. However, most numerical
        functions cannot handle complex numbers, so using ``CDF``
        would require us to wrap every evaluation in a
        ``CDF``-to-:class:`float` conversion routine. That would slow
        things down less than a domain of ``None`` would, but is
        unattractive mainly because of how invasive it would be to
        "fix" the output everywhere.

    Creating a new fast-callable interpreter that has different input
    and output types solves most of the problems with a ``CDF``
    domain, but :func:`fast_callable` and the interpreter classes in
    :mod:`sage.ext.interpreters` are not really written with that in
    mind. The ``domain`` parameter to :func:`fast_callable`, for
    example, is expecting a single Sage ring that corresponds to one
    interpreter. You can make it accept, for example, a string like
    "CDF-to-float", but the hacks required to make that work feel
    wrong.

    Thus we arrive at this solution: a class to wrap the result of
    :func:`fast_callable`. Whenever we need to support intermediate
    complex terms in a numerical routine, we can set ``domain=CDF``
    while creating its fast-callable incarnation, and then wrap the
    result in this class. The :meth:`__call__` method of this class then
    ensures that the ``CDF`` output is converted to a :class:`float` if
    its imaginary part is within an acceptable tolerance.

    EXAMPLES:

    An error is thrown if the answer is complex::

        sage: # needs sage.symbolic
        sage: from sage.ext.fast_callable import FastCallableFloatWrapper
        sage: f = sqrt(x)
        sage: ff = fast_callable(f, vars=[x], domain=CDF)
        sage: fff = FastCallableFloatWrapper(ff, imag_tol=1e-8)
        sage: fff(1)
        1.0
        sage: fff(-1)
        Traceback (most recent call last):
        ...
        ValueError: complex fast-callable function result
        1.0*I for arguments (-1,)
    """
    def __init__(self, ff, imag_tol):
        r"""
        Construct a :class:`FastCallableFloatWrapper`.

        INPUT:

          - ``ff`` -- a fast-callable wrapper over ``CDF``; an instance of
            :class:`sage.ext.interpreters.Wrapper_cdf`, usually constructed
            with :func:`fast_callable`.

          - ``imag_tol`` -- float; how big of an imaginary part we're willing
            to ignore before raising an error

        OUTPUT:

        An instance of :class:`FastCallableFloatWrapper` that can be
        called just like ``ff``, but that always returns a :class:`float`
        if no error is raised. A :exc:`ValueError` is raised if the
        imaginary part of the result exceeds ``imag_tol``.

        EXAMPLES:

        The wrapper will ignore an imaginary part smaller in magnitude
        than ``imag_tol``, but not one larger::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import FastCallableFloatWrapper
            sage: f = x
            sage: ff = fast_callable(f, vars=[x], domain=CDF)
            sage: fff = FastCallableFloatWrapper(ff, imag_tol=1e-8)
            sage: fff(I*1e-9)
            0.0
            sage: fff = FastCallableFloatWrapper(ff, imag_tol=1e-12)
            sage: fff(I*1e-9)
            Traceback (most recent call last):
            ...
            ValueError: complex fast-callable function result 1e-09*I for
            arguments (1.00000000000000e-9*I,)
        """
        self._ff = ff
        self._imag_tol = imag_tol

    def __call__(self, *args):
        r"""
        Evaluate the underlying fast-callable and convert the result to :class:`float`.

        TESTS:

        Evaluation either returns a :class:`float`, or raises a :exc:`ValueError`::

            sage: # needs sage.symbolic
            sage: from sage.ext.fast_callable import FastCallableFloatWrapper
            sage: f = x
            sage: ff = fast_callable(f, vars=[x], domain=CDF)
            sage: fff = FastCallableFloatWrapper(ff, imag_tol=0.1)
            sage: try:
            ....:     result = fff(CDF.random_element())
            ....: except ValueError:
            ....:     result = float(0)
            sage: type(result) is float
            True
        """
        z = self._ff(*args)

        if abs(z.imag()) < self._imag_tol:
            return float(z.real())
        elif isnan(z.real()) and isnan(z.imag()):
            # A fast-callable with domain=float or domain=RDF will
            # return NaN if that's what the underlying function itself
            # returns. When performing the same computation over CDF,
            # however, we get something special; witness:
            #
            #     sage: RDF(0)*math.inf
            #     NaN
            #     sage: CDF(0)*math.inf
            #     NaN + NaN*I
            #
            # Thus for backwards-compatibility, we handle this special
            # case to return NaN to anyone expecting it.
            return nan
        else:
            raise ValueError(f"complex fast-callable function result {z} " +
                             f"for arguments {args}")
