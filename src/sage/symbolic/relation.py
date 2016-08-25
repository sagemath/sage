r"""
Symbolic Equations and Inequalities

Sage can solve symbolic equations and inequalities. For
example, we derive the quadratic formula as follows::

    sage: a,b,c = var('a,b,c')
    sage: qe = (a*x^2 + b*x + c == 0)
    sage: qe
    a*x^2 + b*x + c == 0
    sage: print(solve(qe, x))
    [
    x == -1/2*(b + sqrt(b^2 - 4*a*c))/a,
    x == -1/2*(b - sqrt(b^2 - 4*a*c))/a
    ]


The operator, left hand side, and right hand side
--------------------------------------------------

Operators::

    sage: eqn = x^3 + 2/3 >= x - pi
    sage: eqn.operator()
    <built-in function ge>
    sage: (x^3 + 2/3 < x - pi).operator()
    <built-in function lt>
    sage: (x^3 + 2/3 == x - pi).operator()
    <built-in function eq>

Left hand side::

    sage: eqn = x^3 + 2/3 >= x - pi
    sage: eqn.lhs()
    x^3 + 2/3
    sage: eqn.left()
    x^3 + 2/3
    sage: eqn.left_hand_side()
    x^3 + 2/3

Right hand side::

    sage: (x + sqrt(2) >= sqrt(3) + 5/2).right()
    sqrt(3) + 5/2
    sage: (x + sqrt(2) >= sqrt(3) + 5/2).rhs()
    sqrt(3) + 5/2
    sage: (x + sqrt(2) >= sqrt(3) + 5/2).right_hand_side()
    sqrt(3) + 5/2


Arithmetic
----------
Add two symbolic equations::

    sage: var('a,b')
    (a, b)
    sage: m = 144 == -10 * a + b
    sage: n = 136 == 10 * a + b
    sage: m + n
    280 == 2*b
    sage: int(-144) + m
    0 == -10*a + b - 144

Subtract two symbolic equations::

    sage: var('a,b')
    (a, b)
    sage: m = 144 == 20 * a + b
    sage: n = 136 == 10 * a + b
    sage: m - n
    8 == 10*a
    sage: int(144) - m
    0 == -20*a - b + 144

Multiply two symbolic equations::

    sage: x = var('x')
    sage: m = x == 5*x + 1
    sage: n = sin(x) == sin(x+2*pi, hold=True)
    sage: m * n
    x*sin(x) == (5*x + 1)*sin(2*pi + x)
    sage: m = 2*x == 3*x^2 - 5
    sage: int(-1) * m
    -2*x == -3*x^2 + 5

Divide two symbolic equations::

    sage: x = var('x')
    sage: m = x == 5*x + 1
    sage: n = sin(x) == sin(x+2*pi, hold=True)
    sage: m/n
    x/sin(x) == (5*x + 1)/sin(2*pi + x)
    sage: m = x != 5*x + 1
    sage: n = sin(x) != sin(x+2*pi, hold=True)
    sage: m/n
    x/sin(x) != (5*x + 1)/sin(2*pi + x)

Substitution
------------

Substitution into relations::

    sage: x, a = var('x, a')
    sage: eq = (x^3 + a == sin(x/a)); eq
    x^3 + a == sin(x/a)
    sage: eq.substitute(x=5*x)
    125*x^3 + a == sin(5*x/a)
    sage: eq.substitute(a=1)
    x^3 + 1 == sin(x)
    sage: eq.substitute(a=x)
    x^3 + x == sin(1)
    sage: eq.substitute(a=x, x=1)
    x + 1 == sin(1/x)
    sage: eq.substitute({a:x, x:1})
    x + 1 == sin(1/x)

Solving
-------

We can solve equations::

    sage: x = var('x')
    sage: S = solve(x^3 - 1 == 0, x)
    sage: S
    [x == 1/2*I*sqrt(3) - 1/2, x == -1/2*I*sqrt(3) - 1/2, x == 1]
    sage: S[0]
    x == 1/2*I*sqrt(3) - 1/2
    sage: S[0].right()
    1/2*I*sqrt(3) - 1/2
    sage: S = solve(x^3 - 1 == 0, x, solution_dict=True)
    sage: S
    [{x: 1/2*I*sqrt(3) - 1/2}, {x: -1/2*I*sqrt(3) - 1/2}, {x: 1}]
    sage: z = 5
    sage: solve(z^2 == sqrt(3),z)
    Traceback (most recent call last):
    ...
    TypeError: 5 is not a valid variable.

We illustrate finding multiplicities of solutions::

    sage: f = (x-1)^5*(x^2+1)
    sage: solve(f == 0, x)
    [x == -I, x == I, x == 1]
    sage: solve(f == 0, x, multiplicities=True)
    ([x == -I, x == I, x == 1], [1, 1, 5])

We can also solve many inequalities::

    sage: solve(1/(x-1)<=8,x)
    [[x < 1], [x >= (9/8)]]

We can numerically find roots of equations::

    sage: (x == sin(x)).find_root(-2,2)
    0.0
    sage: (x^5 + 3*x + 2 == 0).find_root(-2,2,x)
    -0.6328345202421523
    sage: (cos(x) == sin(x)).find_root(10,20)
    19.634954084936208

We illustrate some valid error conditions::

    sage: (cos(x) != sin(x)).find_root(10,20)
    Traceback (most recent call last):
    ...
    ValueError: Symbolic equation must be an equality.
    sage: (SR(3)==SR(2)).find_root(-1,1)
    Traceback (most recent call last):
    ...
    RuntimeError: no zero in the interval, since constant expression is not 0.

There must be at most one variable::

    sage: x, y = var('x,y')
    sage: (x == y).find_root(-2,2)
    Traceback (most recent call last):
    ...
    NotImplementedError: root finding currently only implemented in 1 dimension.

Assumptions
-----------

Forgetting assumptions::

    sage: var('x,y')
    (x, y)
    sage: forget() #Clear assumptions
    sage: assume(x>0, y < 2)
    sage: assumptions()
    [x > 0, y < 2]
    sage: (y < 2).forget()
    sage: assumptions()
    [x > 0]
    sage: forget()
    sage: assumptions()
    []


Miscellaneous
-------------

Conversion to Maxima::

    sage: x = var('x')
    sage: eq = (x^(3/5) >= pi^2 + e^i)
    sage: eq._maxima_init_()
    '(_SAGE_VAR_x)^(3/5) >= ((%pi)^(2))+(exp(0+%i*1))'
    sage: e1 = x^3 + x == sin(2*x)
    sage: z = e1._maxima_()
    sage: z.parent() is sage.calculus.calculus.maxima
    True
    sage: z = e1._maxima_(maxima)
    sage: z.parent() is maxima
    True
    sage: z = maxima(e1)
    sage: z.parent() is maxima
    True

Conversion to Maple::

    sage: x = var('x')
    sage: eq = (x == 2)
    sage: eq._maple_init_()
    'x = 2'

Comparison::

    sage: x = var('x')
    sage: (x>0) == (x>0)
    True
    sage: (x>0) == (x>1)
    False
    sage: (x>0) != (x>1)
    True

Variables appearing in the relation::

    sage: var('x,y,z,w')
    (x, y, z, w)
    sage: f =  (x+y+w) == (x^2 - y^2 - z^3);   f
    w + x + y == -z^3 + x^2 - y^2
    sage: f.variables()
    (w, x, y, z)

LaTeX output::

    sage: latex(x^(3/5) >= pi)
    x^{\frac{3}{5}} \geq \pi

When working with the symbolic complex number `I`, notice that comparisons do not
automatically simplify even in trivial situations::

    sage: I^2 == -1
    -1 == -1
    sage: I^2 < 0
    -1 < 0
    sage: (I+1)^4 > 0
    -4 > 0

Nevertheless, if you force the comparison, you get the right answer (:trac:`7160`)::

    sage: bool(I^2 == -1)
    True
    sage: bool(I^2 < 0)
    True
    sage: bool((I+1)^4 > 0)
    False

More Examples
-------------

::

    sage: x,y,a = var('x,y,a')
    sage: f = x^2 + y^2 == 1
    sage: f.solve(x)
    [x == -sqrt(-y^2 + 1), x == sqrt(-y^2 + 1)]

::

    sage: f = x^5 + a
    sage: solve(f==0,x)
    [x == 1/4*(-a)^(1/5)*(sqrt(5) + I*sqrt(2*sqrt(5) + 10) - 1), x == -1/4*(-a)^(1/5)*(sqrt(5) - I*sqrt(-2*sqrt(5) + 10) + 1), x == -1/4*(-a)^(1/5)*(sqrt(5) + I*sqrt(-2*sqrt(5) + 10) + 1), x == 1/4*(-a)^(1/5)*(sqrt(5) - I*sqrt(2*sqrt(5) + 10) - 1), x == (-a)^(1/5)]

You can also do arithmetic with inequalities, as illustrated
below::

    sage: var('x y')
    (x, y)
    sage: f = x + 3 == y - 2
    sage: f
    x + 3 == y - 2
    sage: g = f - 3; g
    x == y - 5
    sage: h =  x^3 + sqrt(2) == x*y*sin(x)
    sage: h
    x^3 + sqrt(2) == x*y*sin(x)
    sage: h - sqrt(2)
    x^3 == x*y*sin(x) - sqrt(2)
    sage: h + f
    x^3 + x + sqrt(2) + 3 == x*y*sin(x) + y - 2
    sage: f = x + 3 < y - 2
    sage: g = 2 < x+10
    sage: f - g
    x + 1 < -x + y - 12
    sage: f + g
    x + 5 < x + y + 8
    sage: f*(-1)
    -x - 3 < -y + 2

TESTS:

We test serializing symbolic equations::

    sage: eqn = x^3 + 2/3 >= x
    sage: loads(dumps(eqn))
    x^3 + 2/3 >= x
    sage: loads(dumps(eqn)) == eqn
    True

AUTHORS:

- Bobby Moretti: initial version (based on a trick that Robert
  Bradshaw suggested).

- William Stein: second version

- William Stein (2007-07-16): added arithmetic with symbolic equations

"""
from __future__ import print_function

import operator


def test_relation_maxima(relation):
    """
    Return True if this (in)equality is definitely true. Return False
    if it is false or the algorithm for testing (in)equality is
    inconclusive.

    EXAMPLES::

        sage: from sage.symbolic.relation import test_relation_maxima
        sage: k = var('k')
        sage: pol = 1/(k-1) - 1/k -1/k/(k-1);
        sage: test_relation_maxima(pol == 0)
        True
        sage: f = sin(x)^2 + cos(x)^2 - 1
        sage: test_relation_maxima(f == 0)
        True
        sage: test_relation_maxima( x == x )
        True
        sage: test_relation_maxima( x != x )
        False
        sage: test_relation_maxima( x > x )
        False
        sage: test_relation_maxima( x^2 > x )
        False
        sage: test_relation_maxima( x + 2 > x )
        True
        sage: test_relation_maxima( x - 2 > x )
        False

    Here are some examples involving assumptions::

        sage: x, y, z = var('x, y, z')
        sage: assume(x>=y,y>=z,z>=x)
        sage: test_relation_maxima(x==z)
        True
        sage: test_relation_maxima(z<x)
        False
        sage: test_relation_maxima(z>y)
        False
        sage: test_relation_maxima(y==z)
        True
        sage: forget()
        sage: assume(x>=1,x<=1)
        sage: test_relation_maxima(x==1)
        True
        sage: test_relation_maxima(x>1)
        False
        sage: test_relation_maxima(x>=1)
        True
        sage: test_relation_maxima(x!=1)
        False
        sage: forget()
        sage: assume(x>0)
        sage: test_relation_maxima(x==0)
        False
        sage: test_relation_maxima(x>-1)
        True
        sage: test_relation_maxima(x!=0)
        True
        sage: test_relation_maxima(x!=1)
        False
        sage: forget()

    TESTS:

    Ensure that ``canonicalize_radical()`` and ``simplify_log`` are not
    used inappropriately, :trac:`17389`. Either one would simplify ``f``
    to zero below::

        sage: x,y = SR.var('x,y')
        sage: assume(y, 'complex')
        sage: f = log(x*y) - (log(x) + log(y))
        sage: f(x=-1, y=i)
        -2*I*pi
        sage: test_relation_maxima(f == 0)
        False
        sage: forget()

    Ensure that the ``sqrt(x^2)`` -> ``abs(x)`` simplification is not
    performed when testing equality::

        sage: assume(x, 'complex')
        sage: f = sqrt(x^2) - abs(x)
        sage: test_relation_maxima(f == 0)
        False
        sage: forget()

    If assumptions are made, ``simplify_rectform()`` is used::

        sage: assume(x, 'real')
        sage: f1 = ( e^(I*x) - e^(-I*x) ) / ( I*e^(I*x) + I*e^(-I*x) )
        sage: f2 = sin(x)/cos(x)
        sage: test_relation_maxima(f1 - f2 == 0)
        True
        sage: forget()

    But not if ``x`` itself is complex::

        sage: assume(x, 'complex')
        sage: f1 = ( e^(I*x) - e^(-I*x) ) / ( I*e^(I*x) + I*e^(-I*x) )
        sage: f2 = sin(x)/cos(x)
        sage: test_relation_maxima(f1 - f2 == 0)
        False
        sage: forget()

    If assumptions are made, then ``simplify_factorial()`` is used::

        sage: n,k = SR.var('n,k')
        sage: assume(n, 'integer')
        sage: assume(k, 'integer')
        sage: f1 = factorial(n+1)/factorial(n)
        sage: f2 = n + 1
        sage: test_relation_maxima(f1 - f2 == 0)
        True
        sage: forget()

    """
    m = relation._maxima_()

    #Handle some basic cases first
    if repr(m) in ['0=0']:
        return True
    elif repr(m) in ['0#0', '1#1']:
        return False

    if relation.operator() == operator.eq: # operator is equality
        try:
            s = m.parent()._eval_line('is (equal(%s,%s))'%(repr(m.lhs()),repr(m.rhs())))
        except TypeError:
            raise ValueError("unable to evaluate the predicate '%s'" % repr(relation))

    elif relation.operator() == operator.ne: # operator is not equal
        try:
            s = m.parent()._eval_line('is (notequal(%s,%s))'%(repr(m.lhs()),repr(m.rhs())))
        except TypeError:
            raise ValueError("unable to evaluate the predicate '%s'" % repr(relation))

    else: # operator is < or > or <= or >=, which Maxima handles fine
        try:
            s = m.parent()._eval_line('is (%s)'%repr(m))
        except TypeError:
            raise ValueError("unable to evaluate the predicate '%s'" % repr(relation))

    if s == 'true':
        return True
    elif s == 'false':
        return False # if neither of these, s=='unknown' and we try a few other tricks

    if relation.operator() != operator.eq:
        return False

    difference = relation.lhs() - relation.rhs()
    if repr(difference) == '0':
        return True

    # Try to apply some simplifications to see if left - right == 0.
    #
    # TODO: If simplify_log() is ever removed from simplify_full(), we
    # can replace all of these individual simplifications with a
    # single call to simplify_full(). That would work in cases where
    # two simplifications are needed consecutively; the current
    # approach does not.
    #
    simp_list = [difference.simplify_factorial(),
                 difference.simplify_rational(),
                 difference.simplify_rectform(),
                 difference.simplify_trig()]
    for f in simp_list:
        try:
            if repr( f() ).strip() == "0":
                return True
                break
        except Exception:
            pass
    return False


def string_to_list_of_solutions(s):
    r"""
    Used internally by the symbolic solve command to convert the output
    of Maxima's solve command to a list of solutions in Sage's symbolic
    package.

    EXAMPLES:

    We derive the (monic) quadratic formula::

        sage: var('x,a,b')
        (x, a, b)
        sage: solve(x^2 + a*x + b == 0, x)
        [x == -1/2*a - 1/2*sqrt(a^2 - 4*b), x == -1/2*a + 1/2*sqrt(a^2 - 4*b)]

    Behind the scenes when the above is evaluated the function
    :func:`string_to_list_of_solutions` is called with input the
    string `s` below::

        sage: s = '[x=-(sqrt(a^2-4*b)+a)/2,x=(sqrt(a^2-4*b)-a)/2]'
        sage: sage.symbolic.relation.string_to_list_of_solutions(s)
         [x == -1/2*a - 1/2*sqrt(a^2 - 4*b), x == -1/2*a + 1/2*sqrt(a^2 - 4*b)]
    """
    from sage.categories.all import Objects
    from sage.structure.sequence import Sequence
    from sage.calculus.calculus import symbolic_expression_from_maxima_string
    v = symbolic_expression_from_maxima_string(s, equals_sub=True)
    return Sequence(v, universe=Objects(), cr_str=True)

###########
# Solving #
###########

def solve(f, *args, **kwds):
    r"""
    Algebraically solve an equation or system of equations (over the
    complex numbers) for given variables. Inequalities and systems
    of inequalities are also supported.

    INPUT:

    -  ``f`` - equation or system of equations (given by a
       list or tuple)

    -  ``*args`` - variables to solve for.

    -  ``solution_dict`` - bool (default: False); if True or non-zero,
       return a list of dictionaries containing the solutions. If there
       are no solutions, return an empty list (rather than a list containing
       an empty dictionary). Likewise, if there's only a single solution,
       return a list containing one dictionary with that solution.

    There are a few optional keywords if you are trying to solve a single
    equation.  They may only be used in that context.

    -  ``multiplicities`` - bool (default: False); if True,
       return corresponding multiplicities.  This keyword is
       incompatible with ``to_poly_solve=True`` and does not make
       any sense when solving inequalities.

    -  ``explicit_solutions`` - bool (default: False); require that
       all roots be explicit rather than implicit. Not used
       when solving inequalities.

    -  ``to_poly_solve`` - bool (default: False) or string; use
       Maxima's ``to_poly_solver`` package to search for more possible
       solutions, but possibly encounter approximate solutions.
       This keyword is incompatible with ``multiplicities=True``
       and is not used when solving inequalities. Setting ``to_poly_solve``
       to 'force' (string) omits Maxima's solve command (useful when
       some solutions of trigonometric equations are lost).


    EXAMPLES::

        sage: x, y = var('x, y')
        sage: solve([x+y==6, x-y==4], x, y)
        [[x == 5, y == 1]]
        sage: solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y)
        [[x == -1/2*I*sqrt(3) - 1/2, y == -sqrt(-1/2*I*sqrt(3) + 3/2)],
         [x == -1/2*I*sqrt(3) - 1/2, y == sqrt(-1/2*I*sqrt(3) + 3/2)],
         [x == 1/2*I*sqrt(3) - 1/2, y == -sqrt(1/2*I*sqrt(3) + 3/2)],
         [x == 1/2*I*sqrt(3) - 1/2, y == sqrt(1/2*I*sqrt(3) + 3/2)],
         [x == 0, y == -1],
         [x == 0, y == 1]]
        sage: solve([sqrt(x) + sqrt(y) == 5, x + y == 10], x, y)
        [[x == -5/2*I*sqrt(5) + 5, y == 5/2*I*sqrt(5) + 5], [x == 5/2*I*sqrt(5) + 5, y == -5/2*I*sqrt(5) + 5]]
        sage: solutions=solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y, solution_dict=True)
        sage: for solution in solutions: print("{} , {}".format(solution[x].n(digits=3), solution[y].n(digits=3)))
        -0.500 - 0.866*I , -1.27 + 0.341*I
        -0.500 - 0.866*I , 1.27 - 0.341*I
        -0.500 + 0.866*I , -1.27 - 0.341*I
        -0.500 + 0.866*I , 1.27 + 0.341*I
        0.000 , -1.00
        0.000 , 1.00

    Whenever possible, answers will be symbolic, but with systems of
    equations, at times approximations will be given, due to the
    underlying algorithm in Maxima::

        sage: sols = solve([x^3==y,y^2==x],[x,y]); sols[-1], sols[0]
        ([x == 0, y == 0], [x == (0.3090169943749475 + 0.9510565162951535*I), y == (-0.8090169943749475 - 0.5877852522924731*I)])
        sage: sols[0][0].rhs().pyobject().parent()
        Complex Double Field

    If ``f`` is only one equation or expression, we use the solve method
    for symbolic expressions, which defaults to exact answers only::

        sage: solve([y^6==y],y)
        [y == 1/4*sqrt(5) + 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4, y == -1/4*sqrt(5) + 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4, y == -1/4*sqrt(5) - 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4, y == 1/4*sqrt(5) - 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4, y == 1, y == 0]
        sage: solve( [y^6 == y], y)==solve( y^6 == y, y)
        True

    Here we demonstrate very basic use of the optional keywords for
    a single expression to be solved::

        sage: ((x^2-1)^2).solve(x)
        [x == -1, x == 1]
        sage: ((x^2-1)^2).solve(x,multiplicities=True)
        ([x == -1, x == 1], [2, 2])
        sage: solve(sin(x)==x,x)
        [x == sin(x)]
        sage: solve(sin(x)==x,x,explicit_solutions=True)
        []
        sage: solve(abs(1-abs(1-x)) == 10, x)
        [abs(abs(x - 1) - 1) == 10]
        sage: solve(abs(1-abs(1-x)) == 10, x, to_poly_solve=True)
        [x == -10, x == 12]

    .. note::

        For more details about solving a single equation, see
        the documentation for the single-expression
        :meth:`~sage.symbolic.expression.Expression.solve`.

    ::

        sage: from sage.symbolic.expression import Expression
        sage: Expression.solve(x^2==1,x)
        [x == -1, x == 1]

    We must solve with respect to actual variables::

        sage: z = 5
        sage: solve([8*z + y == 3, -z +7*y == 0],y,z)
        Traceback (most recent call last):
        ...
        TypeError: 5 is not a valid variable.

    If we ask for dictionaries containing the solutions, we get them::

        sage: solve([x^2-1],x,solution_dict=True)
        [{x: -1}, {x: 1}]
        sage: solve([x^2-4*x+4],x,solution_dict=True)
        [{x: 2}]
        sage: res = solve([x^2 == y, y == 4],x,y,solution_dict=True)
        sage: for soln in res: print("x: %s, y: %s" % (soln[x], soln[y]))
        x: 2, y: 4
        x: -2, y: 4

    If there is a parameter in the answer, that will show up as
    a new variable.  In the following example, ``r1`` is a real free
    variable (because of the ``r``)::

        sage: solve([x+y == 3, 2*x+2*y == 6],x,y)
        [[x == -r1 + 3, y == r1]]

    Especially with trigonometric functions, the dummy variable may
    be implicitly an integer (hence the ``z``)::

        sage: solve([cos(x)*sin(x) == 1/2, x+y == 0],x,y)
        [[x == 1/4*pi + pi*z79, y == -1/4*pi - pi*z79]]

    Expressions which are not equations are assumed to be set equal
    to zero, as with `x` in the following example::

        sage: solve([x, y == 2],x,y)
        [[x == 0, y == 2]]

    If ``True`` appears in the list of equations it is
    ignored, and if ``False`` appears in the list then no
    solutions are returned. E.g., note that the first
    ``3==3`` evaluates to ``True``, not to a
    symbolic equation.

    ::

        sage: solve([3==3, 1.00000000000000*x^3 == 0], x)
        [x == 0]
        sage: solve([1.00000000000000*x^3 == 0], x)
        [x == 0]

    Here, the first equation evaluates to ``False``, so
    there are no solutions::

        sage: solve([1==3, 1.00000000000000*x^3 == 0], x)
        []

    Completely symbolic solutions are supported::

        sage: var('s,j,b,m,g')
        (s, j, b, m, g)
        sage: sys = [ m*(1-s) - b*s*j, b*s*j-g*j ];
        sage: solve(sys,s,j)
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
        sage: solve(sys,(s,j))
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
        sage: solve(sys,[s,j])
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]

    Inequalities can be also solved::

        sage: solve(x^2>8,x)
        [[x < -2*sqrt(2)], [x > 2*sqrt(2)]]

    We use ``use_grobner`` in Maxima if no solution is obtained from
    Maxima's ``to_poly_solve``::

       sage: x,y=var('x y'); c1(x,y)=(x-5)^2+y^2-16; c2(x,y)=(y-3)^2+x^2-9
       sage: solve([c1(x,y),c2(x,y)],[x,y])
       [[x == -9/68*sqrt(55) + 135/68, y == -15/68*sqrt(11)*sqrt(5) + 123/68], [x == 9/68*sqrt(55) + 135/68, y == 15/68*sqrt(11)*sqrt(5) + 123/68]]

    TESTS::

        sage: solve([sin(x)==x,y^2==x],x,y)
        [sin(x) == x, y^2 == x]
        sage: solve(0==1,x)
        Traceback (most recent call last):
        ...
        TypeError:  The first argument must be a symbolic expression or a list of symbolic expressions.

    Test if the empty list is returned, too, when (a list of)
    dictionaries (is) are requested (#8553)::

        sage: solve([SR(0)==1],x)
        []
        sage: solve([SR(0)==1],x,solution_dict=True)
        []
        sage: solve([x==1,x==-1],x)
        []
        sage: solve([x==1,x==-1],x,solution_dict=True)
        []
        sage: solve((x==1,x==-1),x,solution_dict=0)
        []

    Relaxed form, suggested by Mike Hansen (#8553)::

        sage: solve([x^2-1],x,solution_dict=-1)
        [{x: -1}, {x: 1}]
        sage: solve([x^2-1],x,solution_dict=1)
        [{x: -1}, {x: 1}]
        sage: solve((x==1,x==-1),x,solution_dict=-1)
        []
        sage: solve((x==1,x==-1),x,solution_dict=1)
        []

    This inequality holds for any real ``x`` (:trac:`8078`)::

        sage: solve(x^4+2>0,x)
        [x < +Infinity]

    Test for user friendly input handling :trac:`13645`::

        sage: poly.<a,b> = PolynomialRing(RR)
        sage: solve([a+b+a*b == 1], a)
        Traceback (most recent call last):
        ...
        TypeError: The first argument to solve() should be a symbolic expression or a list of symbolic expressions, cannot handle <type 'bool'>
        sage: solve([a, b], (1, a))
        Traceback (most recent call last):
        ...
        TypeError: 1 is not a valid variable.
        sage: solve([x == 1], (1, a))
        Traceback (most recent call last):
        ...
        TypeError: (1, a) are not valid variables.

    Test that the original version of a system in the French Sage book
    now works (:trac:`14306`)::

        sage: var('y,z')
        (y, z)
        sage: solve([x^2 * y * z == 18, x * y^3 * z == 24, x * y * z^4 == 6], x, y, z)
        [[x == 3, y == 2, z == 1], [x == (1.337215067... - 2.685489874...*I), y == (-1.700434271... + 1.052864325...*I), z == (0.9324722294... - 0.3612416661...*I)], ...]
    """
    from sage.symbolic.expression import is_Expression
    if is_Expression(f): # f is a single expression
        ans = f.solve(*args,**kwds)
        return ans

    if not isinstance(f, (list, tuple)):
        raise TypeError("The first argument must be a symbolic expression or a list of symbolic expressions.")

    if len(f)==1:
        # f is a list with a single element
        if is_Expression(f[0]):
            # if its a symbolic expression call solve method of this expression
            return f[0].solve(*args,**kwds)
        # otherwise complain
        raise TypeError("The first argument to solve() should be a symbolic "
                        "expression or a list of symbolic expressions, "
                        "cannot handle %s"%repr(type(f[0])))

    # f is a list of such expressions or equations
    from sage.symbolic.ring import is_SymbolicVariable

    if len(args)==0:
        raise TypeError("Please input variables to solve for.")
    if is_SymbolicVariable(args[0]):
        variables = args
    else:
        variables = tuple(args[0])

    for v in variables:
        if not is_SymbolicVariable(v):
            raise TypeError("%s is not a valid variable."%repr(v))

    try:
        f = [s for s in f if s is not True]
    except TypeError:
        raise ValueError("Unable to solve %s for %s"%(f, args))

    if any(s is False for s in f):
        return []

    from sage.calculus.calculus import maxima
    m = maxima(f)

    try:
        s = m.solve(variables)
    except Exception: # if Maxima gave an error, try its to_poly_solve
        try:
            s = m.to_poly_solve(variables)
        except TypeError as mess: # if that gives an error, raise an error.
            if "Error executing code in Maxima" in str(mess):
                raise ValueError("Sage is unable to determine whether the system %s can be solved for %s"%(f,args))
            else:
                raise

    if len(s)==0: # if Maxima's solve gave no solutions, try its to_poly_solve
        try:
            s = m.to_poly_solve(variables)
        except Exception: # if that gives an error, stick with no solutions
            s = []

    if len(s)==0: # if to_poly_solve gave no solutions, try use_grobner
        try:
            s = m.to_poly_solve(variables,'use_grobner=true')
        except Exception: # if that gives an error, stick with no solutions
            s = []

    sol_list = string_to_list_of_solutions(repr(s))

    # Relaxed form suggested by Mike Hansen (#8553):
    if kwds.get('solution_dict', False):
        if len(sol_list)==0: # fixes IndexError on empty solution list (#8553)
            return []
        if isinstance(sol_list[0], list):
            sol_dict=[dict([[eq.left(),eq.right()] for eq in solution])
                    for solution in sol_list]
        else:
            sol_dict=[{eq.left():eq.right()} for eq in sol_list]

        return sol_dict
    else:
        return sol_list

def solve_mod(eqns, modulus, solution_dict = False):
    r"""
    Return all solutions to an equation or list of equations modulo the
    given integer modulus. Each equation must involve only polynomials
    in 1 or many variables.

    By default the solutions are returned as `n`-tuples, where `n`
    is the number of variables appearing anywhere in the given
    equations. The variables are in alphabetical order.

    INPUT:


    -  ``eqns`` - equation or list of equations

    -  ``modulus`` - an integer

    -  ``solution_dict`` - bool (default: False); if True or non-zero,
       return a list of dictionaries containing the solutions. If there
       are no solutions, return an empty list (rather than a list containing
       an empty dictionary). Likewise, if there's only a single solution,
       return a list containing one dictionary with that solution.


    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(4, 2), (4, 6), (4, 9), (4, 13)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5::

        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    We can solve with respect to a bigger modulus if it consists only of small prime factors::

        sage: [d] = solve_mod([5*x + y == 3, 2*x - 3*y == 9], 3*5*7*11*19*23*29, solution_dict = True)
        sage: d[x]
        12915279
        sage: d[y]
        8610183

    For cases where there are relatively few solutions and the prime
    factors are small, this can be efficient even if the modulus itself
    is large::

        sage: sorted(solve_mod([x^2 == 41], 10^20))
        [(4538602480526452429,), (11445932736758703821,), (38554067263241296179,),
        (45461397519473547571,), (54538602480526452429,), (61445932736758703821,),
        (88554067263241296179,), (95461397519473547571,)]

    We solve a simple equation modulo 2::

        sage: x,y = var('x,y')
        sage: solve_mod([x == y], 2)
        [(0, 0), (1, 1)]

    .. warning::

       The current implementation splits the modulus into prime
       powers, then naively enumerates all possible solutions
       (starting modulo primes and then working up through prime
       powers), and finally combines the solution using the Chinese
       Remainder Theorem.  The interface is good, but the algorithm is
       very inefficient if the modulus has some larger prime factors! Sage
       *does* have the ability to do something much faster in certain
       cases at least by using Groebner basis, linear algebra
       techniques, etc. But for a lot of toy problems this function as
       is might be useful. At least it establishes an interface.


    TESTS:

    Make sure that we short-circuit in at least some cases::

        sage: solve_mod([2*x==1], 2*next_prime(10^50))
        []

    Try multi-equation cases::

        sage: x, y, z = var("x y z")
        sage: solve_mod([2*x^2 + x*y, -x*y+2*y^2+x-2*y, -2*x^2+2*x*y-y^2-x-y], 12)
        [(0, 0), (4, 4), (0, 3), (4, 7)]
        sage: eqs = [-y^2+z^2, -x^2+y^2-3*z^2-z-1, -y*z-z^2-x-y+2, -x^2-12*z^2-y+z]
        sage: solve_mod(eqs, 11)
        [(8, 5, 6)]

    Confirm that modulus 1 now behaves as it should::

        sage: x, y = var("x y")
        sage: solve_mod([x==1], 1)
        [(0,)]
        sage: solve_mod([2*x^2+x*y, -x*y+2*y^2+x-2*y, -2*x^2+2*x*y-y^2-x-y], 1)
        [(0, 0)]


    """
    from sage.rings.all import Integer, Integers, crt_basis
    from sage.symbolic.expression import is_Expression
    from sage.misc.all import cartesian_product_iterator
    from sage.modules.all import vector
    from sage.matrix.all import matrix

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    eqns = [eq if is_Expression(eq) else (eq.lhs()-eq.rhs()) for eq in eqns]
    modulus = Integer(modulus)
    if modulus < 1:
        raise ValueError("the modulus must be a positive integer")
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort(key=repr)

    if modulus == 1: # degenerate case
        ans = [tuple(Integers(1)(0) for v in vars)]
        return ans

    factors = modulus.factor()
    crt_basis = vector(Integers(modulus), crt_basis([p**i for p,i in factors]))
    solutions = []

    has_solution = True
    for p,i in factors:
        solution =_solve_mod_prime_power(eqns, p, i, vars)
        if len(solution) > 0:
            solutions.append(solution)
        else:
            has_solution = False
            break


    ans = []
    if has_solution:
        for solution in cartesian_product_iterator(solutions):
            solution_mat = matrix(Integers(modulus), solution)
            ans.append(tuple(c.dot_product(crt_basis) for c in solution_mat.columns()))

    # if solution_dict == True:
    # Relaxed form suggested by Mike Hansen (#8553):
    if solution_dict:
        sol_dict = [dict(zip(vars, solution)) for solution in ans]
        return sol_dict
    else:
        return ans

def _solve_mod_prime_power(eqns, p, m, vars):
    r"""
    Internal help function for solve_mod, does little checking since it expects
    solve_mod to do that

    Return all solutions to an equation or list of equations modulo p^m.
    Each equation must involve only polynomials
    in 1 or many variables.

    The solutions are returned as `n`-tuples, where `n`
    is the number of variables in vars.

    INPUT:


    -  ``eqns`` - equation or list of equations

    -  ``p`` - a prime

    -  ``i`` - an integer > 0

    -  ``vars`` - a list of variables to solve for


    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(4, 2), (4, 6), (4, 9), (4, 13)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5::

        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    We solve a simple equation modulo 2::

        sage: x,y = var('x,y')
        sage: solve_mod([x == y], 2)
        [(0, 0), (1, 1)]


    .. warning::

       Currently this constructs possible solutions by building up
       from the smallest prime factor of the modulus.  The interface
       is good, but the algorithm is horrible if the modulus isn't the
       product of many small primes! Sage *does* have the ability to
       do something much faster in certain cases at least by using the
       Chinese Remainder Theorem, Groebner basis, linear algebra
       techniques, etc. But for a lot of toy problems this function as
       is might be useful. At the very least, it establishes an
       interface.

    TESTS:

    Confirm we can reproduce the first few terms of OEIS A187719::

        sage: from sage.symbolic.relation import _solve_mod_prime_power
        sage: [sorted(_solve_mod_prime_power([x^2==41], 10, i, [x]))[0][0] for i in [1..13]]
        [1, 21, 71, 1179, 2429, 47571, 1296179, 8703821, 26452429, 526452429,
        13241296179, 19473547571, 2263241296179]

    """
    from sage.rings.all import Integers, PolynomialRing
    from sage.modules.all import vector
    from sage.misc.all import cartesian_product_iterator

    mrunning = 1
    ans = []
    for mi in xrange(m):
        mrunning *= p
        R = Integers(mrunning)
        S = PolynomialRing(R, len(vars), vars)
        eqns_mod = [S(eq) for eq in eqns]
        if mi == 0:
            possibles = cartesian_product_iterator([xrange(len(R)) for _ in xrange(len(vars))])
        else:
            shifts = cartesian_product_iterator([xrange(p) for _ in xrange(len(vars))])
            pairs = cartesian_product_iterator([shifts, ans])
            possibles = (tuple(vector(t)+vector(shift)*(mrunning//p)) for shift, t in pairs)
        ans = list(t for t in possibles if all(e(*t) == 0 for e in eqns_mod))
        if not ans: return ans

    return ans

def solve_ineq_univar(ineq):
    """
    Function solves rational inequality in one variable.

    INPUT:

    - ``ineq`` - inequality in one variable

    OUTPUT:

    - ``list`` -- output is list of solutions as a list of simple inequalities
      output [A,B,C] means (A or B or C) each A, B, C is again a list and
      if A=[a,b], then A means (a and b). The list is empty if there is no
      solution.

    EXAMPLES::

        sage: from sage.symbolic.relation import solve_ineq_univar
        sage: solve_ineq_univar(x-1/x>0)
        [[x > -1, x < 0], [x > 1]]

        sage: solve_ineq_univar(x^2-1/x>0)
        [[x < 0], [x > 1]]

        sage: solve_ineq_univar((x^3-1)*x<=0)
        [[x >= 0, x <= 1]]

    ALGORITHM:

    Calls Maxima command solve_rat_ineq

    AUTHORS:

    - Robert Marik (01-2010)
    """
    ineqvar = ineq.variables()
    if len(ineqvar) != 1:
        raise NotImplementedError("The command solve_ineq_univar accepts univariate inequalities only. Your variables are " + ineqvar)
    ineq0 = ineq._maxima_()
    ineq0.parent().eval("if solve_rat_ineq_loaded#true then (solve_rat_ineq_loaded:true,load(\"solve_rat_ineq.mac\")) ")
    sol = ineq0.solve_rat_ineq().sage()
    if repr(sol)=="all":
        from sage.rings.infinity import Infinity
        sol = [ineqvar[0]<Infinity]
    return sol

def solve_ineq_fourier(ineq,vars=None):
    """
    Solves sytem of inequalities using Maxima and fourier elimination

    Can be used for system of linear inequalities and for some types
    of nonlinear inequalities. For examples see the section EXAMPLES
    below and http://maxima.cvs.sourceforge.net/viewvc/maxima/maxima/share/contrib/fourier_elim/rtest_fourier_elim.mac


    INPUT:

    - ``ineq`` - list with system of inequalities

    - ``vars`` - optionally list with variables for fourier elimination.

    OUTPUT:

    - ``list`` - output is list of solutions as a list of simple inequalities
      output [A,B,C] means (A or B or C) each A, B, C is again a list and
      if A=[a,b], then A means (a and b). The list is empty if there is no
      solution.

    EXAMPLES::

        sage: from sage.symbolic.relation import solve_ineq_fourier
        sage: y=var('y')
        sage: solve_ineq_fourier([x+y<9,x-y>4],[x,y])
        [[y + 4 < x, x < -y + 9, y < (5/2)]]
        sage: solve_ineq_fourier([x+y<9,x-y>4],[y,x])
        [[y < min(x - 4, -x + 9)]]

        sage: solve_ineq_fourier([x^2>=0])
        [[x < +Infinity]]

        sage: solve_ineq_fourier([log(x)>log(y)],[x,y])
        [[y < x, 0 < y]]
        sage: solve_ineq_fourier([log(x)>log(y)],[y,x])
        [[0 < y, y < x, 0 < x]]

    Note that different systems will find default variables in different
    orders, so the following is not tested::

        sage: solve_ineq_fourier([log(x)>log(y)])  # random (one of the following appears)
        [[0 < y, y < x, 0 < x]]
        [[y < x, 0 < y]]

    ALGORITHM:

    Calls Maxima command fourier_elim

    AUTHORS:

    - Robert Marik (01-2010)
    """
    if vars is None:
        setvars = set([])
        for i in (ineq):
            setvars = setvars.union(set(i.variables()))
            vars =[i for i in setvars]
    ineq0 = [i._maxima_() for i in ineq]
    ineq0[0].parent().eval("if fourier_elim_loaded#true then (fourier_elim_loaded:true,load(\"fourier_elim\"))")
    sol = ineq0[0].parent().fourier_elim(ineq0,vars)
    ineq0[0].parent().eval("or_to_list(x):=\
        if not atom(x) and op(x)=\"or\" then args(x) \
        else [x]")
    sol = sol.or_to_list().sage()
    if repr(sol) == "[emptyset]":
        sol = []
    if repr(sol) == "[universalset]":
        from sage.rings.infinity import Infinity
        sol = [[i<Infinity for i in vars]]
    return sol

def solve_ineq(ineq, vars=None):
    """
    Solves inequalities and systems of inequalities using Maxima.
    Switches between rational inequalities
    (sage.symbolic.relation.solve_ineq_rational)
    and fourier elimination (sage.symbolic.relation.solve_ineq_fouried).
    See the documentation of these functions for more details.

    INPUT:

    - ``ineq`` - one inequality or a list of inequalities

      Case1: If ``ineq`` is one equality, then it should be rational
      expression in one varible. This input is passed to
      sage.symbolic.relation.solve_ineq_univar function.

      Case2: If ``ineq`` is a list involving one or more
      inequalities, than the input is passed to
      sage.symbolic.relation.solve_ineq_fourier function. This
      function can be used for system of linear inequalities and
      for some types of nonlinear inequalities. See
      http://maxima.cvs.sourceforge.net/viewvc/maxima/maxima/share/contrib/fourier_elim/rtest_fourier_elim.mac
      for a big gallery of problems covered by this algorithm.

    - ``vars`` - optional parameter with list of variables. This list
      is used only if fourier elimination is used. If omitted or if
      rational inequality is solved, then variables are determined
      automatically.

    OUTPUT:

    - ``list`` -- output is list of solutions as a list of simple inequalities
      output [A,B,C] means (A or B or C) each A, B, C is again a list and
      if A=[a,b], then A means (a and b).

    EXAMPLES::

        sage: from sage.symbolic.relation import solve_ineq

    Inequalities in one variable. The variable is detected automatically::

        sage: solve_ineq(x^2-1>3)
        [[x < -2], [x > 2]]

        sage: solve_ineq(1/(x-1)<=8)
        [[x < 1], [x >= (9/8)]]

    System of inequalities with automatically detected inequalities::

        sage: y=var('y')
        sage: solve_ineq([x-y<0,x+y-3<0],[y,x])
        [[x < y, y < -x + 3, x < (3/2)]]
        sage: solve_ineq([x-y<0,x+y-3<0],[x,y])
        [[x < min(-y + 3, y)]]

    Note that although Sage will detect the variables automatically,
    the order it puts them in may depend on the system, so the following
    command is only guaranteed to give you one of the above answers::

        sage: solve_ineq([x-y<0,x+y-3<0])  # random
        [[x < y, y < -x + 3, x < (3/2)]]

    ALGORITHM:

    Calls solve_ineq_fourier if inequalities are list and
    solve_ineq_univar of the inequality is symbolic expression. See
    the description of these commands for more details related to the
    set of inequalities which can be solved. The list is empty if
    there is no solution.

    AUTHORS:

    - Robert Marik (01-2010)
    """
    if isinstance(ineq,list):
        return(solve_ineq_fourier(ineq, vars))
    else:
        return(solve_ineq_univar(ineq))
