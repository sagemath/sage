r"""
Utilities for Calculus

This module defines helper functions which are used for simplifications
and display of symbolic expressions.

AUTHORS:

- Michal Bejger (2015) : class :class:`ExpressionNice`
- Eric Gourgoulhon (2015) : simplification functions
- Travis Scrimshaw (2016): review tweaks

"""

#******************************************************************************
#
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import division

from sage.symbolic.expression import Expression

def simplify_sqrt_real(expr):
    r"""
    Simplify ``sqrt`` in symbolic expressions in the real domain.

    EXAMPLES:

    Simplifications of basic expressions::

        sage: from sage.manifolds.utilities import simplify_sqrt_real
        sage: simplify_sqrt_real( sqrt(x^2) )
        abs(x)
        sage: assume(x<0)
        sage: simplify_sqrt_real( sqrt(x^2) )
        -x
        sage: simplify_sqrt_real( sqrt(x^2-2*x+1) )
        -x + 1
        sage: simplify_sqrt_real( sqrt(x^2) + sqrt(x^2-2*x+1) )
        -2*x + 1

    This improves over Sage's
    :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`,
    which yields incorrect results when ``x < 0``::

        sage: forget()  # removes the assumption x<0
        sage: sqrt(x^2).canonicalize_radical()
        x
        sage: assume(x<0)
        sage: sqrt(x^2).canonicalize_radical() # wrong output
        x
        sage: sqrt(x^2-2*x+1).canonicalize_radical() # wrong output
        x - 1
        sage: ( sqrt(x^2) + sqrt(x^2-2*x+1) ).canonicalize_radical() # wrong output
        2*x - 1

    Simplification of nested ``sqrt``'s::

        sage: forget()  # removes the assumption x<0
        sage: simplify_sqrt_real( sqrt(1 + sqrt(x^2)) )
        sqrt(abs(x) + 1)
        sage: assume(x<0)
        sage: simplify_sqrt_real( sqrt(1 + sqrt(x^2)) )
        sqrt(-x + 1)
        sage: simplify_sqrt_real( sqrt(x^2 + sqrt(4*x^2) + 1) )
        -x + 1

    Again, :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
    fails on the last one::

        sage: (sqrt(x^2 + sqrt(4*x^2) + 1)).canonicalize_radical()  # wrong output
        x + 1

    """
    from sage.symbolic.ring import SR
    from sage.functions.other import sqrt
    # 1/ Search for the sqrt's in expr
    sexpr = str(expr)
    if 'sqrt(' not in sexpr:  # no sqrt to simplify
        return expr
    if 'D[' in sexpr:
        return expr    #!# the code below is not capable of simplifying
                       # expressions with symbolic derivatives denoted by Pynac
                       # symbols of the type D[0]
    # Lists to store the positions of all the top-level sqrt's in sexpr:
    pos_sqrts = []  # position of first character, i.e. 's' of 'sqrt(...)'
    pos_after = []  # position of character immediatelty after 'sqrt(...)'
    the_sqrts = []  # the sqrt sub-expressions in sexpr, i.e. 'sqrt(...)'
    pos_max = len(sexpr) - 6
    pos = 0
    while pos < pos_max:
        if sexpr[pos:pos+5] == 'sqrt(':
            pos_sqrts.append(pos)
            parenth = 1
            scan = pos+5
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1
            the_sqrts.append( sexpr[pos:scan] )
            pos_after.append(scan)
            pos = scan
        else:
            pos += 1
    # 2/ Search for sub-sqrt's:
    for i in range(len(the_sqrts)):
        argum = the_sqrts[i][5:-1]  # the sqrt argument
        if 'sqrt(' in argum:
            simpl = simplify_sqrt_real(SR(argum))
            the_sqrts[i] = 'sqrt(' + str(simpl) + ')'
    # 3/ Simplifications of the sqrt's
    new_expr = ""    # will contain the result
    pos0 = 0
    for i, pos in enumerate(pos_sqrts):
        # radcan is called on each sqrt:
        x = SR(the_sqrts[i])
        argum = x.operands()[0] # the argument of sqrt
        den = argum.denominator()
        if not (den == 1):  # the argument of sqrt is a fraction
            # NB: after #19312 (integrated in Sage 6.10.beta7), the above
            # cannot be written as
            #    if den != 1!:
            num = argum.numerator()
            if num < 0 or den < 0:
                x = sqrt(-num) / sqrt(-den)  # new equivalent expression for x
        simpl = SR(x._maxima_().radcan())
        if str(simpl)[:5] == 'sqrt(' or str(simpl)[:7] == '1/sqrt(':
            # no further simplification seems possible:
            ssimpl = str(simpl)
        else:
            # the absolute value of radcan's output is taken, the call to
            # simplify() taking into account possible assumptions regarding the
            # sign of simpl:
            ssimpl = str(abs(simpl).simplify())
        # search for abs(1/sqrt(...)) term to simplify it into 1/sqrt(...):
        pstart = ssimpl.find('abs(1/sqrt(')
        if pstart != -1:
            ssimpl = ssimpl[:pstart] + ssimpl[pstart+3:] # getting rid of 'abs'
        new_expr += sexpr[pos0:pos] + '(' + ssimpl + ')'
        pos0 = pos_after[i]
    new_expr += sexpr[pos0:]
    return SR(new_expr)

def simplify_abs_trig(expr):
    r"""
    Simplify ``abs(sin(...))`` in symbolic expressions.

    EXAMPLES::

        sage: forget()  # for doctests only
        sage: M = Manifold(3, 'M', structure='topological')
        sage: X.<x,y,z> = M.chart(r'x y:(0,pi) z:(-pi/3,0)')
        sage: X.coord_range()
        x: (-oo, +oo); y: (0, pi); z: (-1/3*pi, 0)

    Since ``x`` spans all `\RR`, no simplification of ``abs(sin(x))``
    occurs, while ``abs(sin(y))`` and ``abs(sin(3*z))`` are correctly
    simplified, given that `y \in (0,\pi)` and `z \in (-\pi/3,0)`::

        sage: from sage.manifolds.utilities import simplify_abs_trig
        sage: simplify_abs_trig( abs(sin(x)) + abs(sin(y)) + abs(sin(3*z)) )
        abs(sin(x)) + sin(y) - sin(3*z)

    Note that neither Sage's function
    :meth:`~sage.symbolic.expression.Expression.simplify_trig` nor
    :meth:`~sage.symbolic.expression.Expression.simplify_full`
    works in this case::

        sage: s = abs(sin(x)) + abs(sin(y)) + abs(sin(3*z))
        sage: s.simplify_trig()
        abs(4*cos(z)^2 - 1)*abs(sin(z)) + abs(sin(x)) + abs(sin(y))
        sage: s.simplify_full()
        abs(4*cos(z)^2 - 1)*abs(sin(z)) + abs(sin(x)) + abs(sin(y))

    despite the following assumptions hold::

        sage: assumptions()
        [x is real, y is real, y > 0, y < pi, z is real, z > -1/3*pi, z < 0]

    Additional checks are::

        sage: simplify_abs_trig( abs(sin(y/2)) )  # shall simplify
        sin(1/2*y)
        sage: simplify_abs_trig( abs(sin(2*y)) )  # must not simplify
        abs(sin(2*y))
        sage: simplify_abs_trig( abs(sin(z/2)) )  # shall simplify
        -sin(1/2*z)
        sage: simplify_abs_trig( abs(sin(4*z)) )  # must not simplify
        abs(sin(4*z))

    """
    from sage.symbolic.ring import SR
    from sage.symbolic.constants import pi
    sexpr = str(expr)
    if 'abs(sin(' not in sexpr:  # nothing to simplify
        return expr
    tp = []
    val = []
    for pos in range(len(sexpr)):
        if sexpr[pos:pos+8] == 'abs(sin(':
            # finding the end of abs argument:
            scan = pos+4 # start of abs
            parenth = 1
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1
            pos_abs_end = scan
            # finding the end of sin argument:
            scan = pos+8 # start of sin
            parenth = 1
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1
            pos_sin_end = scan
            # if the abs contains only the sinus, the simplification can be tried:
            if pos_sin_end == pos_abs_end-1:
                tp.append(pos)
                val.append( sexpr[pos:pos_abs_end] )
    simp = []
    for v in val:
        # argument of the sinus:
        sx = v[8:-2]
        x = SR(sx)
        if x>=0 and x<=pi:
            simp.append('sin(' + sx + ')')
        elif x>=-pi and x<=0:
            simp.append('(-sin(' + sx + '))')
        else:
            simp.append(v)  # no simplification is applicable
    nexpr = ""
    pos0 = 0
    for i, pos in enumerate(tp):
        nexpr += sexpr[pos0:pos] + simp[i]
        pos0 = pos + len(val[i])
    nexpr += sexpr[pos0:]
    return SR(nexpr)


def simplify_chain_real(expr):
    r"""
    Apply a chain of simplifications to a symbolic expression, assuming the
    real domain.

    This is the simplification chain used in calculus involving coordinate
    functions on real manifolds, as implemented in
    :class:`~sage.manifolds.coord_func_symb.CoordFunctionSymb`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sage.symbolic.expression.Expression.simplify_factorial`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :func:`simplify_sqrt_real`
    #. :func:`simplify_abs_trig`
    #. :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_log`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`

    EXAMPLES:

    We consider variables that are coordinates of a chart on a real manifold::

        sage: forget()  # for doctest only
        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart('x:(0,1) y')

    The following assumptions then hold::

        sage: assumptions()
        [x is real, x > 0, x < 1, y is real]

    and we have::

        sage: from sage.manifolds.utilities import simplify_chain_real
        sage: s = sqrt(y^2)
        sage: simplify_chain_real(s)
        abs(y)

    The above result is correct since ``y`` is real. It is obtained by
    :meth:`~sage.symbolic.expression.Expression.simplify_real` as well,
    but not by :meth:`~sage.symbolic.expression.Expression.simplify_full`::

        sage: s.simplify_real()
        abs(y)
        sage: s.simplify_full()
        sqrt(y^2)

    Furthermore, we have::

        sage: s = sqrt(x^2-2*x+1)
        sage: simplify_chain_real(s)
        -x + 1

    which is correct since `x \in (0,1)`. On this example, neither
    :meth:`~sage.symbolic.expression.Expression.simplify_real`
    nor :meth:`~sage.symbolic.expression.Expression.simplify_full`,
    nor :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
    give satisfactory results::

        sage: s.simplify_real()  # unsimplified output
        sqrt(x^2 - 2*x + 1)
        sage: s.simplify_full()  # unsimplified output
        sqrt(x^2 - 2*x + 1)
        sage: s.canonicalize_radical()  # wrong output since x in (0,1)
        x - 1

    Other simplifications::

        sage: s = abs(sin(pi*x))
        sage: simplify_chain_real(s)  # correct output since x in (0,1)
        sin(pi*x)
        sage: s.simplify_real()  # unsimplified output
        abs(sin(pi*x))
        sage: s.simplify_full()  # unsimplified output
        abs(sin(pi*x))

    ::

        sage: s = cos(y)^2 + sin(y)^2
        sage: simplify_chain_real(s)
        1
        sage: s.simplify_real()  # unsimplified output
        cos(y)^2 + sin(y)^2
        sage: s.simplify_full()  # OK
        1

    """
    expr = expr.simplify_factorial()
    expr = expr.simplify_trig()
    expr = expr.simplify_rational()
    expr = simplify_sqrt_real(expr)
    expr = simplify_abs_trig(expr)
    expr = expr.canonicalize_radical()
    expr = expr.simplify_log('one')
    expr = expr.simplify_rational()
    expr = expr.simplify_trig()
    return expr

def simplify_chain_generic(expr):
    r"""
    Apply a chain of simplifications to a symbolic expression.

    This is the simplification chain used in calculus involving coordinate
    functions on manifolds over fields different from `\RR`, as implemented in
    :class:`~sage.manifolds.coord_func_symb.CoordFunctionSymb`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sage.symbolic.expression.Expression.simplify_factorial`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rectform`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :meth:`~sage.symbolic.expression.Expression.expand_sum`

    NB: for the time being, this is identical to
    :meth:`~sage.symbolic.expression.Expression.simplify_full`.

    EXAMPLES:

    We consider variables that are coordinates of a chart on a complex
    manifold::

        sage: forget()  # for doctest only
        sage: M = Manifold(2, 'M', structure='topological', field='complex')
        sage: X.<x,y> = M.chart()

    Then neither ``x`` nor ``y`` is assumed to be real::

        sage: assumptions()
        []

    Accordingly, ``simplify_chain_generic`` does not simplify
    ``sqrt(x^2)`` to ``abs(x)``::

        sage: from sage.manifolds.utilities import simplify_chain_generic
        sage: s = sqrt(x^2)
        sage: simplify_chain_generic(s)
        sqrt(x^2)

    This contrasts with the behavior of
    :func:`~sage.manifolds.utilities.simplify_chain_real`.

    Other simplifications::

        sage: s = (x+y)^2 - x^2 -2*x*y - y^2
        sage: simplify_chain_generic(s)
        0
        sage: s = (x^2 - 2*x + 1) / (x^2 -1)
        sage: simplify_chain_generic(s)
        (x - 1)/(x + 1)
        sage: s = cos(2*x) - 2*cos(x)^2 + 1
        sage: simplify_chain_generic(s)
        0

    """
    expr = expr.simplify_factorial()
    expr = expr.simplify_rectform()
    expr = expr.simplify_trig()
    expr = expr.simplify_rational()
    expr = expr.expand_sum()
    return expr



#******************************************************************************

class ExpressionNice(Expression):
    r"""
    Subclass of :class:`~sage.symbolic.expression.Expression` for a
    "human-friendly" display of partial derivatives and the possibility to
    shorten the display by skipping the arguments of symbolic functions.

    INPUT:

    - ``ex`` -- symbolic expression

    EXAMPLES:

    An expression formed with callable symbolic expressions::

        sage: var('x y z')
        (x, y, z)
        sage: f = function('f')(x, y)
        sage: g = f.diff(y).diff(x)
        sage: h = function('h')(y, z)
        sage: k = h.diff(z)
        sage: fun = x*g + y*(k-z)^2

    The standard Pynac display of partial derivatives::

        sage: fun
        y*(z - D[1](h)(y, z))^2 + x*D[0, 1](f)(x, y)
        sage: latex(fun)
        y {\left(z - D[1]\left(h\right)\left(y, z\right)\right)}^{2}
         + x D[0, 1]\left(f\right)\left(x, y\right)

    With :class:`ExpressionNice`, the Pynac notation ``D[...]`` is replaced
    by textbook-like notation::

        sage: from sage.manifolds.utilities import ExpressionNice
        sage: ExpressionNice(fun)
        y*(z - d(h)/dz)^2 + x*d^2(f)/dxdy
        sage: latex(ExpressionNice(fun))
        y {\left(z - \frac{\partial\,h}{\partial z}\right)}^{2}
         + x \frac{\partial^2\,f}{\partial x\partial y}

    An example when function variables are themselves functions::

        sage: f = function('f')(x, y)
        sage: g = function('g')(x, f)  # the second variable is the function f
        sage: fun = (g.diff(x))*x - x^2*f.diff(x,y)
        sage: fun
        -x^2*D[0, 1](f)(x, y) + (D[0](f)(x, y)*D[1](g)(x, f(x, y)) + D[0](g)(x, f(x, y)))*x
        sage: ExpressionNice(fun)
        -x^2*d^2(f)/dxdy + (d(f)/dx*d(g)/d(f(x, y)) + d(g)/dx)*x
        sage: latex(ExpressionNice(fun))
        -x^{2} \frac{\partial^2\,f}{\partial x\partial y}
         + {\left(\frac{\partial\,f}{\partial x}
           \frac{\partial\,g}{\partial \left( f\left(x, y\right) \right)}
         + \frac{\partial\,g}{\partial x}\right)} x

    Note that ``D[1](g)(x, f(x,y))`` is rendered as ``d(g)/d(f(x, y))``.

    An example with multiple differentiations::

        sage: fun = f.diff(x,x,y,y,x)*x
        sage: fun
        x*D[0, 0, 0, 1, 1](f)(x, y)
        sage: ExpressionNice(fun)
        x*d^5(f)/dx^3dy^2
        sage: latex(ExpressionNice(fun))
        x \frac{\partial^5\,f}{\partial x ^ 3\partial y ^ 2}

    Parentheses are added around powers of partial derivatives to avoid any
    confusion::

        sage: fun = f.diff(y)^2
        sage: fun
        D[1](f)(x, y)^2
        sage: ExpressionNice(fun)
        (d(f)/dy)^2
        sage: latex(ExpressionNice(fun))
        \left(\frac{\partial\,f}{\partial y}\right)^{2}

    The explicit mention of function arguments can be omitted for the sake of
    brevity::

        sage: fun = fun*f
        sage: ExpressionNice(fun)
        f(x, y)*(d(f)/dy)^2
        sage: Manifold.options.omit_function_arguments=True
        sage: ExpressionNice(fun)
        f*(d(f)/dy)^2
        sage: latex(ExpressionNice(fun))
        f \left(\frac{\partial\,f}{\partial y}\right)^{2}
        sage: Manifold.options._reset()
        sage: ExpressionNice(fun)
        f(x, y)*(d(f)/dy)^2
        sage: latex(ExpressionNice(fun))
        f\left(x, y\right) \left(\frac{\partial\,f}{\partial y}\right)^{2}

    """
    def __init__(self, ex):
        r"""
        Initialize ``self``.

        TESTS::

            sage: f = function('f')(x)
            sage: df = f.diff(x)
            sage: df
            D[0](f)(x)
            sage: from sage.manifolds.utilities import ExpressionNice
            sage: df_nice = ExpressionNice(df)
            sage: df_nice
            d(f)/dx

        """
        from sage.symbolic.ring import SR
        self._parent = SR
        Expression.__init__(self, SR, x=ex)

    def _repr_(self):
        r"""
        String representation of the object.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: f = function('f')(x, y)
            sage: g = f.diff(y).diff(x)
            sage: h = function('h')(y, z)
            sage: k = h.diff(z)
            sage: fun = x*g + y*(k-z)^2
            sage: fun
            y*(z - D[1](h)(y, z))^2 + x*D[0, 1](f)(x, y)
            sage: from sage.manifolds.utilities import ExpressionNice
            sage: ExpressionNice(fun)
            y*(z - d(h)/dz)^2 + x*d^2(f)/dxdy

        """
        d = self._parent._repr_element_(self)

        import re

        # find all occurences of diff
        list_d = []
        _list_derivatives(self, list_d)

        # process the list
        for m in list_d:
            funcname = m[1]
            diffargs = m[3]
            numargs = len(diffargs)

            if numargs > 1:
                numargs = "^" + str(numargs)
            else:
                numargs = ""

            variables = m[4]
            strv = list(str(v) for v in variables)

            # checking if the variable is composite
            for i in range(len(strv)):
                if bool(re.search(r'[+|-|/|*|^|(|)]', strv[i])):
                    strv[i] = "(" + strv[i] + ")"

            # dictionary to group multiple occurences of differentiation: d/dxdx -> d/dx^2 etc.
            occ = dict((i, strv[i] + "^" + str(diffargs.count(i))
                       if (diffargs.count(i)>1) else strv[i])
                       for i in diffargs)

            res = "d" + str(numargs) + "(" + str(funcname) + ")/d" + "d".join(
                               [i for i in occ.values()])

            # str representation of the operator
            s = self._parent._repr_element_(m[0])

            # if diff operator is raised to some power (m[5]), put brackets around
            if m[5]:
                res = "(" + res + ")^" + str(m[5])
                o = s + "^" + str(m[5])
            else:
                o = s

            d = d.replace(o, res)

        from sage.manifolds.manifold import TopologicalManifold
        if TopologicalManifold.options.omit_function_arguments:
            list_f = []
            _list_functions(self, list_f)

            for m in list_f:
                d = d.replace(m[1] + m[2], m[1])

        return d

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        EXAMPLE::

            sage: var('x y z')
            (x, y, z)
            sage: f = function('f')(x, y)
            sage: g = f.diff(y).diff(x)
            sage: h = function('h')(y, z)
            sage: k = h.diff(z)
            sage: fun = x*g + y*(k-z)^2
            sage: fun
            y*(z - D[1](h)(y, z))^2 + x*D[0, 1](f)(x, y)
            sage: from sage.manifolds.utilities import ExpressionNice
            sage: ExpressionNice(fun)
            y*(z - d(h)/dz)^2 + x*d^2(f)/dxdy
            sage: latex(ExpressionNice(fun))
            y {\left(z - \frac{\partial\,h}{\partial z}\right)}^{2} + x \frac{\partial^2\,f}{\partial x\partial y}

        Testing the behavior if no latex_name of the function is given::

            sage: f = function('f_x')(x, y)
            sage: fun = f.diff(y)
            sage: latex(ExpressionNice(fun))
            \frac{\partial\,f_{x}}{\partial y}

        If latex_name, it should be used in LaTeX output:

            sage: f = function('f_x', latex_name=r"{\cal F}")(x,y)
            sage: fun = f.diff(y)
            sage: latex(ExpressionNice(fun))
            \frac{\partial\,{\cal F}}{\partial y}

        """
        d = self._parent._latex_element_(self)

        import re

        # find all occurences of diff
        list_d = []
        _list_derivatives(self, list_d)

        for m in list_d:
            if str(m[1]) == str(m[2]):
                funcname = str(m[1])
            else:
                funcname = str(m[2])

            diffargs = m[3]
            numargs = len(diffargs)

            if numargs > 1:
                numargs = "^" + str(numargs)
            else:
                numargs = ""

            variables = m[4]

            from sage.misc.latex import latex
            strv = [str(v) for v in variables]
            latv = [latex(v) for v in variables]

            # checking if the variable is composite
            for i, val in enumerate(strv):
                if bool(re.search(r'[+|-|/|*|^|(|)]', val)):
                    latv[i] = "\left(" + latv[i] + "\\right)"

            # dictionary to group multiple occurences of differentiation: d/dxdx -> d/dx^2 etc.
            occ = {i: (latv[i] + "^" + latex(diffargs.count(i))
                       if diffargs.count(i) > 1 else latv[i])
                   for i in diffargs}

            res = "\\frac{\partial" + numargs + "\," + funcname + \
                  "}{\partial " + "\partial ".join(i for i in occ.values()) + "}"

            # representation of the operator
            s = self._parent._latex_element_(m[0])

            # if diff operator is raised to some power (m[5]), put brackets around
            if m[5]:
                res = "\left(" + res + "\\right)^{" + str(m[5]) + "}"
                o = s + "^{" + str(m[5]) + "}"
            else:
                o = s

            d = d.replace(o, res)

        from sage.manifolds.manifold import TopologicalManifold
        if TopologicalManifold.options.omit_function_arguments:
            list_f = []
            _list_functions(self, list_f)

            for m in list_f:
                d = d.replace(str(m[3]) + str(m[4]), str(m[3]))

        return d


def _list_derivatives(ex, list_d, exponent=0):
    r"""
    Function to find the occurrences of ``FDerivativeOperator`` in a symbolic
    expression; inspired by
    http://ask.sagemath.org/question/10256/how-can-extract-different-terms-from-a-symbolic-expression/?answer=26136#post-id-26136

    INPUT:

    - ``ex`` -- symbolic expression to be analyzed
    - ``exponent`` -- (optional) exponent of ``FDerivativeOperator``,
      passed to a next level in the expression tree

    OUTPUT:

    - ``list_d`` -- tuple containing the details of ``FDerivativeOperator``
      found, in the following order:

      1. operator
      2. function name
      3. LaTeX function name
      4. parameter set
      5. operands
      6. exponent (if found, else 0)

    TESTS::

        sage: f = function('f_x', latex_name=r"{\cal F}")(x)
        sage: df = f.diff(x)^2
        sage: from sage.manifolds.utilities import _list_derivatives
        sage: list_d = []
        sage: _list_derivatives(df, list_d)
        sage: list_d
        [(D[0](f_x)(x), 'f_x', {\cal F}, [0], [x], 2)]

    """
    op = ex.operator()
    operands = ex.operands()

    import operator
    from sage.misc.latex import latex, latex_variable_name
    from sage.symbolic.operators import FDerivativeOperator

    if op:
        if op is operator.pow:
            if isinstance(operands[0].operator(), FDerivativeOperator):
                exponent = operands[1]

        if isinstance(op, FDerivativeOperator):
            parameter_set = op.parameter_set()
            function = repr(op.function())
            latex_function = latex(op.function())

            # case when no latex_name given
            if function == latex_function:
                latex_function = latex_variable_name(str(op.function()))

            list_d.append((ex, function, latex_function, parameter_set,
                           operands, exponent))

        for operand in operands:
            _list_derivatives(operand, list_d, exponent)


def _list_functions(ex, list_f):
    r"""
    Function to find the occurrences of symbolic functions in a symbolic
    expression.

    INPUT:

    - ``ex`` -- symbolic expression to be analyzed

    OUTPUT:

    - ``list_f`` -- tuple containing the details of a symbolic function found,
      in the following order:

      1. operator
      2. function name
      3. arguments
      4. LaTeX version of function name
      5. LaTeX version of arguments

    TESTS::

        sage: var('x y z')
        (x, y, z)
        sage: f = function('f', latex_name=r"{\cal F}")(x, y)
        sage: g = function('g_x')(x, y)
        sage: d = sin(x)*g.diff(x)*x*f - x^2*f.diff(x,y)/g
        sage: from sage.manifolds.utilities import _list_functions
        sage: list_f = []
        sage: _list_functions(d, list_f)
        sage: list_f
        [(f, 'f', '(x, y)', {\cal F}, \left(x, y\right)),
         (g_x, 'g_x', '(x, y)', 'g_{x}', \left(x, y\right))]

    """
    op = ex.operator()
    operands = ex.operands()

    from sage.misc.latex import latex, latex_variable_name

    if op:
        # FIXME: This hack is needed because the NewSymbolicFunction is
        #   a class defined inside of the *function* function_factory().
        if str(type(op)) == "<class 'sage.symbolic.function_factory.NewSymbolicFunction'>":
            repr_function = repr(op)
            latex_function = latex(op)

            # case when no latex_name given
            if repr_function == latex_function:
                latex_function = latex_variable_name(str(op))

            repr_args = repr(ex.arguments())
            # remove comma in case of singleton
            if len(ex.arguments()) == 1:
                repr_args = repr_args.replace(",","")

            latex_args = latex(ex.arguments())

            list_f.append((op, repr_function, repr_args, latex_function, latex_args))

        for operand in operands:
            _list_functions(operand, list_f)

#******************************************************************************

def set_axes_labels(graph, xlabel, ylabel, zlabel, **kwds):
    r"""
    Set axes labels for a 3D graphics object ``graph``.

    This is a workaround for the lack of axes labels in 3D plots.
    This sets the labels as :func:`~sage.plot.plot3d.shapes2.text3d`
    objects at locations determined from the bounding box of the
    graphic object ``graph``.

    INPUT:

    - ``graph`` -- :class:`~sage.plot.plot3d.base.Graphics3d`;
      a 3D graphic object
    - ``xlabel`` -- string for the x-axis label
    - ``ylabel`` -- string for the y-axis label
    - ``zlabel`` -- string for the z-axis label
    - ``**kwds`` -- options (e.g. color) for text3d

    OUTPUT:

    - the 3D graphic object with text3d labels added

    EXAMPLES::

        sage: g = sphere()
        sage: g.all
        [Graphics3d Object]
        sage: from sage.manifolds.utilities import set_axes_labels
        sage: ga = set_axes_labels(g, 'X', 'Y', 'Z', color='red')
        sage: ga.all  # the 3D frame has now axes labels
        [Graphics3d Object, Graphics3d Object,
         Graphics3d Object, Graphics3d Object]

    """
    from sage.plot.plot3d.shapes2 import text3d
    xmin, ymin, zmin = graph.bounding_box()[0]
    xmax, ymax, zmax = graph.bounding_box()[1]
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    x1 = xmin + dx / 2
    y1 = ymin + dy / 2
    z1 = zmin + dz / 2
    xmin1 = xmin - dx / 20
    xmax1 = xmax + dx / 20
    ymin1 = ymin - dy / 20
    zmin1 = zmin - dz / 20
    graph += text3d('  ' + xlabel, (x1, ymin1, zmin1), **kwds)
    graph += text3d('  ' + ylabel, (xmax1, y1, zmin1), **kwds)
    graph += text3d('  ' + zlabel, (xmin1, ymin1, z1), **kwds)
    return graph

