# sage_setup: distribution = sagemath-symbolics
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

You can even substitute multivariable and matrix
expressions::

    sage: x,y = var('x, y')
    sage: M = Matrix([[x+1,y],[x^2,y^3]]); M
    [x + 1     y]
    [  x^2   y^3]
    sage: M.substitute({x:0,y:1})
    [1 1]
    [0 1]

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

We can also solve equations involving matrices. The following
example defines a multivariable function ``f(x,y)``, then solves
for where the partial derivatives with respect to ``x``
and ``y`` are zero. Then it substitutes one of the solutions
into the Hessian matrix ``H`` for ``f``::

    sage: f(x,y) = x^2*y+y^2+y
    sage: solutions = solve(list(f.diff()),[x,y],solution_dict=True)
    sage: solutions == [{x: -I, y: 0}, {x: I, y: 0}, {x: 0, y: -1/2}]
    True
    sage: H = f.diff(2) # Hessian matrix
    sage: H.subs(solutions[2])
    [(x, y) |--> -1  (x, y) |--> 0]
    [ (x, y) |--> 0  (x, y) |--> 2]
    sage: H(x,y).subs(solutions[2])
    [-1  0]
    [ 0  2]

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

    sage: SR(I)^2 == -1
    -1 == -1
    sage: SR(I)^2 < 0
    -1 < 0
    sage: (SR(I)+1)^4 > 0
    -4 > 0

Nevertheless, if you force the comparison, you get the right answer (:issue:`7160`)::

    sage: bool(SR(I)^2 == -1)
    True
    sage: bool(SR(I)^2 < 0)
    True
    sage: bool((SR(I)+1)^4 > 0)
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
from itertools import product
import operator


def test_relation_maxima(relation):
    """
    Return ``True`` if this (in)equality is definitely true. Return ``False``
    if it is false or the algorithm for testing (in)equality is inconclusive.

    EXAMPLES::

        sage: from sage.symbolic.relation import test_relation_maxima
        sage: k = var('k')
        sage: pol = 1/(k-1) - 1/k -1/k/(k-1)
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
    used inappropriately, :issue:`17389`. Either one would simplify ``f``
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

    In case an equation is to be solved for non-integers, ''assume()''
    is used::

        sage: k = var('k')
        sage: assume(k,'noninteger')
        sage: solve([k^3==1],k)
        [k == 1/2*I*sqrt(3) - 1/2, k == -1/2*I*sqrt(3) - 1/2]
        sage: assumptions()
        [k is noninteger]
    """
    m = relation._maxima_()

    # Handle some basic cases first
    if repr(m) in ['0=0']:
        return True
    elif repr(m) in ['0#0', '1#1']:
        return False

    if relation.operator() == operator.eq:  # operator is equality
        try:
            s = m.parent()._eval_line('is (equal(%s,%s))' % (repr(m.lhs()),
                                                             repr(m.rhs())))
        except TypeError:
            raise ValueError("unable to evaluate the predicate '%s'" % repr(relation))

    elif relation.operator() == operator.ne:  # operator is not equal
        try:
            s = m.parent()._eval_line('is (notequal(%s,%s))' % (repr(m.lhs()),
                                                                repr(m.rhs())))
        except TypeError:
            raise ValueError("unable to evaluate the predicate '%s'" % repr(relation))

    else:  # operator is < or > or <= or >=, which Maxima handles fine
        try:
            s = m.parent()._eval_line('is (%s)' % repr(m))
        except TypeError:
            raise ValueError("unable to evaluate the predicate '%s'" % repr(relation))

    if s == 'true':
        return True
    elif s == 'false':
        return False  # if neither of these, s=='unknown' and we try a few other tricks

    if relation.operator() != operator.eq:
        return False

    difference = relation.lhs() - relation.rhs()
    if difference.is_trivial_zero():
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
            if f().is_trivial_zero():
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
    from sage.categories.objects import Objects
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

    - ``f`` -- equation or system of equations (given by a list or tuple)

    - ``*args`` -- variables to solve for

    - ``solution_dict`` -- boolean (default: ``False``); if ``True`` or nonzero,
      return a list of dictionaries containing the solutions. If there
      are no solutions, return an empty list (rather than a list containing
      an empty dictionary). Likewise, if there's only a single solution,
      return a list containing one dictionary with that solution.

    There are a few optional keywords if you are trying to solve a single
    equation.  They may only be used in that context.

    - ``multiplicities`` -- boolean (default: ``False``); if True,
      return corresponding multiplicities.  This keyword is
      incompatible with ``to_poly_solve=True`` and does not make
      any sense when solving inequalities.

    - ``explicit_solutions`` -- boolean (default: ``False``); require that
      all roots be explicit rather than implicit. Not used
      when solving inequalities.

    - ``to_poly_solve`` -- boolean (default: ``False``) or string; use
      Maxima's ``to_poly_solver`` package to search for more possible
      solutions, but possibly encounter approximate solutions.
      This keyword is incompatible with ``multiplicities=True``
      and is not used when solving inequalities. Setting ``to_poly_solve``
      to 'force' (string) omits Maxima's solve command (useful when
      some solutions of trigonometric equations are lost).

    - ``algorithm`` -- string (default: ``'maxima'``); to use SymPy's
      solvers set this to 'sympy'. Note that SymPy is always used
      for diophantine equations. Another choice is 'giac'.

    - ``domain`` -- string (default: ``'complex'``); setting this to 'real'
      changes the way SymPy solves single equations; inequalities
      are always solved in the real domain.

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
        sage: solutions = solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y, solution_dict=True)
        sage: for solution in solutions: print("{} , {}".format(solution[x].n(digits=3), solution[y].n(digits=3)))
        -0.500 - 0.866*I , -1.27 + 0.341*I
        -0.500 - 0.866*I , 1.27 - 0.341*I
        -0.500 + 0.866*I , -1.27 - 0.341*I
        -0.500 + 0.866*I , 1.27 + 0.341*I
        0.000 , -1.00
        0.000 , 1.00

    Whenever possible, answers will be symbolic, but with systems of
    equations, at times approximations will be given by Maxima, due to the
    underlying algorithm::

        sage: sols = solve([x^3==y,y^2==x], [x,y]); sols[-1], sols[0] # abs tol 1e-15
        ([x == 0, y == 0],
         [x == (0.3090169943749475 + 0.9510565162951535*I),
          y == (-0.8090169943749475 - 0.5877852522924731*I)])
        sage: sols[0][0].rhs().pyobject().parent()
        Complex Double Field

        sage: solve([y^6==y],y)
        [y == 1/4*sqrt(5) + 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4,
         y == -1/4*sqrt(5) + 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4,
         y == -1/4*sqrt(5) - 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4,
         y == 1/4*sqrt(5) - 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4,
         y == 1,
         y == 0]
        sage: solve( [y^6 == y], y)==solve( y^6 == y, y)
        True

    Here we demonstrate very basic use of the optional keywords::

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
    a new variable.  In the following example, ``r1`` is an arbitrary
    constant (because of the ``r``)::

        sage: forget()
        sage: x, y = var('x,y')
        sage: solve([x+y == 3, 2*x+2*y == 6],x,y)
        [[x == -r1 + 3, y == r1]]

        sage: var('b, c')
        (b, c)
        sage: solve((b-1)*(c-1), [b,c])
        [[b == 1, c == r...], [b == r..., c == 1]]

    Especially with trigonometric functions, the dummy variable may
    be implicitly an integer (hence the ``z``)::

        sage: solve( sin(x)==cos(x), x, to_poly_solve=True)
        [x == 1/4*pi + pi*z...]
        sage: solve([cos(x)*sin(x) == 1/2, x+y == 0],x,y)
        [[x == 1/4*pi + pi*z..., y == -1/4*pi - pi*z...]]

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
        sage: sys = [ m*(1-s) - b*s*j, b*s*j-g*j ]
        sage: solve(sys,s,j)
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
        sage: solve(sys,(s,j))
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]
        sage: solve(sys,[s,j])
        [[s == 1, j == 0], [s == g/b, j == (b - g)*m/(b*g)]]

        sage: z = var('z')
        sage: solve((x-z)^2==2, x)
        [x == z - sqrt(2), x == z + sqrt(2)]

    Inequalities can be also solved::

        sage: solve(x^2>8,x)
        [[x < -2*sqrt(2)], [x > 2*sqrt(2)]]
        sage: x,y = var('x,y'); (ln(x)-ln(y)>0).solve(x)
        [[log(x) - log(y) > 0]]
        sage: x,y = var('x,y'); (ln(x)>ln(y)).solve(x)  # random
        [[0 < y, y < x, 0 < x]]
        [[y < x, 0 < y]]

    A simple example to show the use of the keyword
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

    The following examples show the use of the keyword ``to_poly_solve``::

        sage: solve(abs(1-abs(1-x)) == 10, x)
        [abs(abs(x - 1) - 1) == 10]
        sage: solve(abs(1-abs(1-x)) == 10, x, to_poly_solve=True)
        [x == -10, x == 12]

        sage: var('Q')
        Q
        sage: solve(Q*sqrt(Q^2 + 2) - 1, Q)
        [Q == 1/sqrt(Q^2 + 2)]

    The following example is a regression in Maxima 5.39.0.
    It used to be possible to get one more solution here,
    namely ``1/sqrt(sqrt(2) + 1)``, see
    https://sourceforge.net/p/maxima/bugs/3276/::

        sage: solve(Q*sqrt(Q^2 + 2) - 1, Q, to_poly_solve=True)
        [Q == -sqrt(-sqrt(2) - 1), Q == sqrt(sqrt(2) + 1)*(sqrt(2) - 1)]

    An effort is made to only return solutions that satisfy
    the current assumptions::

        sage: solve(x^2==4, x)
        [x == -2, x == 2]
        sage: assume(x<0)
        sage: solve(x^2==4, x)
        [x == -2]
        sage: solve((x^2-4)^2 == 0, x, multiplicities=True)
        ([x == -2], [2])
        sage: solve(x^2==2, x)
        [x == -sqrt(2)]
        sage: z = var('z')
        sage: solve(x^2==2-z, x)
        [x == -sqrt(-z + 2)]
        sage: assume(x, 'rational')
        sage: solve(x^2 == 2, x)
        []

    In some cases it may be worthwhile to directly use ``to_poly_solve``
    if one suspects some answers are being missed::

        sage: forget()
        sage: solve(cos(x)==0, x)
        [x == 1/2*pi]
        sage: solve(cos(x)==0, x, to_poly_solve=True)
        [x == 1/2*pi]
        sage: solve(cos(x)==0, x, to_poly_solve='force')
        [x == 1/2*pi + pi*z...]

    The same may also apply if a returned unsolved expression has a
    denominator, but the original one did not::

        sage: solve(cos(x) * sin(x) == 1/2, x, to_poly_solve=True)
        [sin(x) == 1/2/cos(x)]
        sage: solve(cos(x) * sin(x) == 1/2, x, to_poly_solve=True, explicit_solutions=True)
        [x == 1/4*pi + pi*z...]
        sage: solve(cos(x) * sin(x) == 1/2, x, to_poly_solve='force')
        [x == 1/4*pi + pi*z...]

    We use ``use_grobner`` in Maxima if no solution is obtained from
    Maxima's ``to_poly_solve``::

        sage: x,y = var('x y')
        sage: c1(x,y) = (x-5)^2+y^2-16
        sage: c2(x,y) = (y-3)^2+x^2-9
        sage: solve([c1(x,y),c2(x,y)],[x,y])
        [[x == -9/68*sqrt(55) + 135/68, y == -15/68*sqrt(55) + 123/68],
         [x == 9/68*sqrt(55) + 135/68, y == 15/68*sqrt(55) + 123/68]]

    We use SymPy for Diophantine equations, see
    ``Expression.solve_diophantine``::

        sage: assume(x, 'integer')
        sage: assume(z, 'integer')
        sage: solve((x-z)^2==2, x)
        []

        sage: forget()

    The following shows some more of SymPy's capabilities that cannot be
    handled by Maxima::

        sage: _ = var('t')
        sage: r = solve([x^2 - y^2/exp(x), y-1], x, y, algorithm='sympy')
        sage: (r[0][x], r[0][y])
        (2*lambert_w(-1/2), 1)
        sage: solve(-2*x**3 + 4*x**2 - 2*x + 6 > 0, x, algorithm='sympy')
        [x < 1/3*(1/2)^(1/3)*(9*sqrt(77) + 79)^(1/3) + 2/3*(1/2)^(2/3)/(9*sqrt(77) + 79)^(1/3) + 2/3]
        sage: solve(sqrt(2*x^2 - 7) - (3 - x),x,algorithm='sympy')
        [x == -8, x == 2]
        sage: solve(sqrt(2*x + 9) - sqrt(x + 1) - sqrt(x + 4),x,algorithm='sympy')
        [x == 0]
        sage: r = solve([x + y + z + t, -z - t], x, y, z, t, algorithm='sympy')
        sage: (r[0][x], r[0][z])
        (-y, -t)
        sage: r = solve([x^2+y+z, y+x^2+z, x+y+z^2], x, y,z, algorithm='sympy')
        sage: (r[0][x], r[0][y])
        (z, -(z + 1)*z)
        sage: (r[1][x], r[1][y])
        (-z + 1, -z^2 + z - 1)
        sage: solve(abs(x + 3) - 2*abs(x - 3),x,algorithm='sympy',domain='real')
        [x == 1, x == 9]

    We cannot translate all results from SymPy but we can at least
    print them::

        sage: solve(sinh(x) - 2*cosh(x),x,algorithm='sympy')
        [ImageSet(Lambda(_n, I*(2*_n*pi + pi/2) + log(sqrt(3))), Integers),
         ImageSet(Lambda(_n, I*(2*_n*pi - pi/2) + log(sqrt(3))), Integers)]
        sage: solve(2*sin(x) - 2*sin(2*x), x,algorithm='sympy')
        [ImageSet(Lambda(_n, 2*_n*pi), Integers),
         ImageSet(Lambda(_n, 2*_n*pi + pi), Integers),
         ImageSet(Lambda(_n, 2*_n*pi + 5*pi/3), Integers),
         ImageSet(Lambda(_n, 2*_n*pi + pi/3), Integers)]

        sage: solve(x^5 + 3*x^3 + 7, x, algorithm='sympy')[0] # known bug
        complex_root_of(x^5 + 3*x^3 + 7, 0)

    A basic interface to Giac is provided::

        sage: solve([(2/3)^x-2], [x], algorithm='giac')
        ...[[-log(2)/(log(3) - log(2))]]

        sage: f = (sin(x) - 8*cos(x)*sin(x))*(sin(x)^2 + cos(x)) - (2*cos(x)*sin(x) - sin(x))*(-2*sin(x)^2 + 2*cos(x)^2 - cos(x))
        sage: solve(f, x, algorithm='giac')
        ...[-2*arctan(sqrt(2)), 0, 2*arctan(sqrt(2)), pi]

        sage: x, y = SR.var('x,y')
        sage: solve([x+y-4,x*y-3],[x,y],algorithm='giac')
        [[1, 3], [3, 1]]

    TESTS::

        sage: solve([sin(x)==x,y^2==x],x,y)
        [sin(x) == x, y^2 == x]
        sage: solve(0==1,x)
        Traceback (most recent call last):
        ...
        TypeError:  The first argument must be a symbolic expression or a list of symbolic expressions.

    Test if the empty list is returned, too, when (a list of)
    dictionaries (is) are requested (:issue:`8553`)::

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

    Relaxed form, suggested by Mike Hansen (:issue:`8553`)::

        sage: solve([x^2-1],x,solution_dict=-1)
        [{x: -1}, {x: 1}]
        sage: solve([x^2-1],x,solution_dict=1)
        [{x: -1}, {x: 1}]
        sage: solve((x==1,x==-1),x,solution_dict=-1)
        []
        sage: solve((x==1,x==-1),x,solution_dict=1)
        []

    This inequality holds for any real ``x`` (:issue:`8078`)::

        sage: solve(x^4+2>0,x)
        [x < +Infinity]

    Test for user friendly input handling :issue:`13645`::

        sage: poly.<a,b> = PolynomialRing(RR)
        sage: solve([a+b+a*b == 1], a)
        Traceback (most recent call last):
        ...
        TypeError: a is not a valid variable.
        sage: a,b = var('a,b')
        sage: solve([a+b+a*b == 1], a)
        [a == -(b - 1)/(b + 1)]
        sage: solve([a, b], (1, a))
        Traceback (most recent call last):
        ...
        TypeError: 1 is not a valid variable.
        sage: solve([x == 1], (1, a))
        Traceback (most recent call last):
        ...
        TypeError: 1 is not a valid variable.
        sage: x.solve((1,2))
        Traceback (most recent call last):
        ...
        TypeError: 1 is not a valid variable.

    Test that the original version of a system in the French Sage book
    now works (:issue:`14306`)::

        sage: var('y,z')
        (y, z)
        sage: solve([x^2 * y * z == 18, x * y^3 * z == 24, x * y * z^4 == 6], x, y, z)
        [[x == 3, y == 2, z == 1],
         [x == (1.337215067... - 2.685489874...*I),
          y == (-1.700434271... + 1.052864325...*I),
          z == (0.9324722294... - 0.3612416661...*I)],
         ...]

    :issue:`13286` fixed::

        sage: solve([x-4], [x])
        [x == 4]

    Test for a list of non-symbolic expressions as first argument
    (:issue:`31714`)::

        sage: solve([1], x)
        Traceback (most recent call last):
        ...
        TypeError: The first argument to solve() should be a symbolic expression
        or a list of symbolic expressions.
    """
    from sage.structure.element import Expression
    explicit_solutions = kwds.get('explicit_solutions', None)
    multiplicities = kwds.get('multiplicities', None)
    to_poly_solve = kwds.get('to_poly_solve', None)
    solution_dict = kwds.get('solution_dict', False)
    algorithm = kwds.get('algorithm', None)
    domain = kwds.get('domain', None)

    if len(args) > 1:
        x = args
    else:
        x = args[0]
    if isinstance(x, (list, tuple)):
        for i in x:
            if not isinstance(i, Expression):
                raise TypeError("%s is not a valid variable." % repr(i))
    elif x is None:
        vars = f.variables()
        if len(vars) == 0:
            if multiplicities:
                return [], []
            else:
                return []
        x = vars[0]
    elif not isinstance(x, Expression):
        raise TypeError("%s is not a valid variable." % repr(x))

    if isinstance(f, (list, tuple)) and len(f) == 1:
        # f is a list with a single element
        if isinstance(f[0], Expression):
            f = f[0]
        else:
            raise TypeError("The first argument to solve() should be a "
                            "symbolic expression or a list of symbolic "
                            "expressions.")

    if isinstance(f, Expression):  # f is a single expression
        return _solve_expression(f, x, explicit_solutions, multiplicities, to_poly_solve, solution_dict, algorithm, domain)

    if not isinstance(f, (list, tuple)):
        raise TypeError("The first argument must be a symbolic expression or a list of symbolic expressions.")

    # f is a list of such expressions or equations

    if not args:
        raise TypeError("Please input variables to solve for.")
    if isinstance(x, Expression) and x.is_symbol():
        variables = args
    else:
        variables = tuple(x)

    for v in variables:
        if not (isinstance(v, Expression) and v.is_symbol()):
            raise TypeError("%s is not a valid variable." % repr(v))

    try:
        f = [s for s in f if s is not True]
    except TypeError:
        raise ValueError("Unable to solve %s for %s" % (f, args))

    if any(s is False for s in f):
        return []

    if algorithm == 'sympy':
        from sympy import solve as ssolve
        from sage.interfaces.sympy import sympy_set_to_list
        if isinstance(f, Expression):  # f is a single expression
            sympy_f = f._sympy_()
        else:
            sympy_f = [s._sympy_() for s in f]
        if isinstance(f, Expression) and f.is_symbol():
            sympy_vars = (x._sympy_(),)
        else:
            sympy_vars = tuple([v._sympy_() for v in x])
        if len(sympy_vars) > 1 or not isinstance(f, Expression):
            ret = ssolve(sympy_f, sympy_vars, dict=True)
            if isinstance(ret, dict):
                if solution_dict:
                    l = []
                    for d in ret:
                        r = {}
                        for (v, ex) in d.items():
                            r[v._sage_()] = ex._sage_()
                        l.append(r)
                    return l
                else:
                    return [[v._sage_() == ex._sage_()
                             for v, ex in d.items()]
                            for d in ret]
            elif isinstance(ret, list):
                l = []
                for sol in ret:
                    r = {}
                    for (v, ex) in sol.items():
                        r[v._sage_()] = ex._sage_()
                    l.append(r)
                return l
            else:
                return sympy_set_to_list(ret, sympy_vars)

    if algorithm == 'giac':
        return _giac_solver(f, x, solution_dict)

    from sage.calculus.calculus import maxima
    m = maxima(f)

    try:
        s = m.solve(variables)
    except Exception:  # if Maxima gave an error, try its to_poly_solve
        try:
            s = m.to_poly_solve(variables)
        except TypeError as mess:  # if that gives an error, raise an error.
            if "Error executing code in Maxima" in str(mess):
                raise ValueError("Sage is unable to determine whether the system %s can be solved for %s" % (f, args))
            else:
                raise

    if len(s) == 0:  # if Maxima's solve gave no solutions, try its to_poly_solve
        try:
            s = m.to_poly_solve(variables)
        except Exception:  # if that gives an error, stick with no solutions
            s = []

    if len(s) == 0:  # if to_poly_solve gave no solutions, try use_grobner
        try:
            s = m.to_poly_solve(variables, 'use_grobner=true')
        except Exception:  # if that gives an error, stick with no solutions
            s = []

    sol_list = string_to_list_of_solutions(repr(s))

    # Relaxed form suggested by Mike Hansen (#8553):
    if kwds.get('solution_dict', None):
        if not sol_list:  # fixes IndexError on empty solution list (#8553)
            return []
        if isinstance(sol_list[0], list):
            sol_dict = [{eq.left(): eq.right() for eq in solution}
                        for solution in sol_list]
        else:
            sol_dict = [{eq.left(): eq.right()} for eq in sol_list]

        return sol_dict
    else:
        return sol_list


def _solve_expression(f, x, explicit_solutions, multiplicities,
                      to_poly_solve, solution_dict, algorithm, domain):
    """
    Solve an expression ``f``. For more information, see :func:`solve`.

    .. NOTE::

        This is an auxiliary function only meant to be called
        from :func:`solve`.

    TESTS:

    :issue:`7325` (solving inequalities)::

        sage: (x^2>1).solve(x)
        [[x < -1], [x > 1]]

    Catch error message from Maxima::

        sage: solve(acot(x),x)
        Traceback (most recent call last):
        ...
        TypeError: ECL says: cot: argument 0 isn't in the domain of cot.

    ::

        sage: solve(acot(x),x,to_poly_solve=True)
        Traceback (most recent call last):
        ...
        TypeError: ECL says: cot: argument 0 isn't in the domain of cot.

    :issue:`7491` fixed::

        sage: y = var('y')
        sage: solve(y==y,y)
        [y == r1]
        sage: solve(y==y,y,multiplicities=True)
        ([y == r1], [])

        sage: from sage.symbolic.assumptions import GenericDeclaration
        sage: GenericDeclaration(x, 'rational').assume()
        sage: solve(x^2 == 2, x)
        []
        sage: forget()

    :issue:`8390` fixed::

        sage: solve(sin(x)==1/2,x)
        [x == 1/6*pi]

    ::

        sage: solve(sin(x)==1/2,x,to_poly_solve=True)
        [x == 1/6*pi]

    ::

        sage: solve(sin(x)==1/2, x, to_poly_solve='force')
        [x == 5/6*pi + 2*pi*z..., x == 1/6*pi + 2*pi*z...]

    :issue:`11618` fixed::

        sage: g(x)=0
        sage: solve(g(x)==0,x,solution_dict=True)
        [{x: r1}]

    :issue:`17128`: fixed::

        sage: var('x,y')
        (x, y)
        sage: f = x+y
        sage: sol = f.solve([x, y], solution_dict=True)
        sage: sol[0].get(x) + sol[0].get(y)
        0

    :issue:`16651` fixed::

        sage: (x^7-x-1).solve(x, to_poly_solve=True)     # abs tol 1e-6
        [x == 1.11277569705,
         x == (-0.363623519329 - 0.952561195261*I),
         x == (0.617093477784 - 0.900864951949*I),
         x == (-0.809857800594 - 0.262869645851*I),
         x == (-0.809857800594 + 0.262869645851*I),
         x == (0.617093477784 + 0.900864951949*I),
         x == (-0.363623519329 + 0.952561195261*I)]

    :issue:`31452` fixed::

        sage: solve([x==3], [x], solution_dict=True)
        [{x: 3}]
        sage: solve([x==3], [x], solution_dict=True, algorithm='sympy')
        [{x: 3}]
    """
    from sage.structure.element import Expression

    if f.is_relational():
        if f.operator() is not operator.eq:
            if algorithm == 'sympy':
                from sympy import S, solveset
                from sage.interfaces.sympy import sympy_set_to_list
                if isinstance(x, Expression) and x.is_symbol():
                    sympy_vars = (x._sympy_(),)
                else:
                    sympy_vars = tuple([v._sympy_() for v in x])
                ret = solveset(f._sympy_(), sympy_vars[0], S.Reals)
                return sympy_set_to_list(ret, sympy_vars)
            elif algorithm == 'giac':
                return _giac_solver(f, x, solution_dict)
            else:
                try:
                    return solve_ineq(f)  # trying solve_ineq_univar
                except Exception:
                    pass
                try:
                    return solve_ineq([f])  # trying solve_ineq_fourier
                except Exception:
                    raise NotImplementedError("solving only implemented for equalities and few special inequalities, see solve_ineq")
        ex = f
    else:
        ex = (f == 0)

    if multiplicities and to_poly_solve:
        raise NotImplementedError("to_poly_solve does not return multiplicities")
    # check if all variables are assumed integer;
    # if so, we have a Diophantine

    def has_integer_assumption(v):
        from sage.symbolic.assumptions import assumptions, GenericDeclaration
        alist = assumptions()
        return any(isinstance(a, GenericDeclaration) and a.has(v) and
                   a._assumption in ['even', 'odd', 'integer', 'integervalued']
                   for a in alist)
    if len(ex.variables()) and all(has_integer_assumption(var) for var in ex.variables()):
        return f.solve_diophantine(x, solution_dict=solution_dict)

    if algorithm == 'sympy':
        from sympy import S, solveset
        from sage.interfaces.sympy import sympy_set_to_list
        if isinstance(x, Expression) and x.is_symbol():
            sympy_vars = (x._sympy_(),)
        else:
            sympy_vars = tuple([v._sympy_() for v in x])
        if domain == 'real':
            ret = solveset(ex._sympy_(), sympy_vars[0], S.Reals)
        else:
            ret = solveset(ex._sympy_(), sympy_vars[0])
        ret = sympy_set_to_list(ret, sympy_vars)
        if solution_dict:
            ret = [{sol.left(): sol.right()} for sol in ret]
        return ret

    if algorithm == 'giac':
        return _giac_solver(f, x, solution_dict)

    # from here on, maxima is used for solution
    m = ex._maxima_()
    P = m.parent()
    if explicit_solutions:
        P.eval('solveexplicit: true')  # switches Maxima to looking for only explicit solutions
    try:
        if to_poly_solve != 'force':
            s = m.solve(x).str()
        else:  # omit Maxima's solve command
            s = str([])
    except TypeError as mess:  # if Maxima's solve has an error, we catch it
        if "Error executing code in Maxima" in str(mess):
            s = str([])
        else:
            raise
    if explicit_solutions:
        P.eval('solveexplicit: false')  # switches Maxima back to default

    if s == 'all':
        if solution_dict:
            ans = [{x: f.parent().var('r1')}]
        else:
            ans = [x == f.parent().var('r1')]
        if multiplicities:
            return ans, []
        else:
            return ans

    X = string_to_list_of_solutions(s)  # our initial list of solutions

    if multiplicities:  # to_poly_solve does not return multiplicities, so in this case we end here
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
    if to_poly_solve:
        if len(X) == 0:
            # Maxima's solve gave no solutions
            solutions_so_far = [ex]
            ignore_exceptions = True
        else:
            solutions_so_far = X
            ignore_exceptions = False
        X = []
        for eq in solutions_so_far:
            if eq.lhs().is_symbol() and (eq.lhs() == x) and (x not in eq.rhs().variables()):
                X.append(eq)
                continue
            try:
                m = eq._maxima_()
                s = m.to_poly_solve(x, options='algexact:true')
                T = string_to_list_of_solutions(repr(s))
                X.extend([t[0] for t in T])
            except TypeError as mess:
                if ignore_exceptions:
                    continue
                elif "Error executing code in Maxima" in str(mess) or \
                     "unable to make sense of Maxima expression" in \
                     str(mess):
                    if not explicit_solutions:
                        X.append(eq)  # we keep this implicit solution
                else:
                    raise

    # make sure all the assumptions are satisfied
    from sage.symbolic.assumptions import assumptions
    to_check = assumptions()
    if to_check:
        for ix, soln in reversed(list(enumerate(X))):
            if soln.lhs().is_symbol():
                if any(a.contradicts(soln) for a in to_check):
                    del X[ix]
                    if multiplicities:
                        del ret_multiplicities[ix]
                    continue

    if solution_dict:
        if isinstance(x, (list, tuple)) and len(x) > 1:
            X = [{sol.left(): sol.right() for sol in b} for b in X]
        else:
            X = [{sol.left(): sol.right()} for sol in X]

    if multiplicities:
        return X, ret_multiplicities
    else:
        return X


def _giac_solver(f, x, solution_dict=False):
    """
    Solve a system of equations using libgiac.

    INPUT:

    - ``f`` -- equation or list of equations
    - ``x`` -- variable or list of variables
    - ``solution_dict`` -- boolean (default: ``False``)

    EXAMPLES::

        sage: solve([(2/3)^x-2], [x], algorithm='giac')
        ...[[-log(2)/(log(3) - log(2))]]
        sage: solve([(2/3)^x-2], [x], algorithm='giac', solution_dict=True)
        ...[{x: -log(2)/(log(3) - log(2))}]

        sage: f = (sin(x) - 8*cos(x)*sin(x))*(sin(x)^2 + cos(x)) - (2*cos(x)*sin(x) - sin(x))*(-2*sin(x)^2 + 2*cos(x)^2 - cos(x))
        sage: solve(f, x, algorithm='giac')
        ...[-2*arctan(sqrt(2)), 0, 2*arctan(sqrt(2)), pi]
        sage: solve(f, x, algorithm='giac', solution_dict=True)
        ...[{x: -2*arctan(sqrt(2))}, {x: 0}, {x: 2*arctan(sqrt(2))}, {x: pi}]

        sage: x, y = SR.var('x,y')
        sage: solve([x+y-7,x*y-10],[x,y],algorithm='giac')
        [[2, 5], [5, 2]]
    """
    from sage.libs.giac.giac import libgiac
    giac_f = libgiac(f)
    giac_vars = libgiac(x)
    ret = giac_f.solve(giac_vars)
    sols = ret.sage()
    if solution_dict:
        if not sols:
            return []
        if isinstance(sols[0], list):
            return [dict(zip(x, solution)) for solution in sols]
        return [{x: sx} for sx in sols]
    return sols


def solve_mod(eqns, modulus, solution_dict=False):
    r"""
    Return all solutions to an equation or list of equations modulo the
    given integer modulus. Each equation must involve only polynomials
    in 1 or many variables.

    By default the solutions are returned as `n`-tuples, where `n`
    is the number of variables appearing anywhere in the given
    equations. The variables are in alphabetical order.

    INPUT:

    - ``eqns`` -- equation or list of equations

    - ``modulus`` -- integer

    - ``solution_dict`` -- boolean (default: ``False``); if ``True`` or nonzero,
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
    from sage.rings.finite_rings.integer_mod_ring import Integers
    from sage.rings.integer import Integer
    from sage.rings.integer_ring import crt_basis
    from sage.structure.element import Expression
    from sage.modules.free_module_element import vector
    from sage.matrix.constructor import matrix

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    eqns = [eq if isinstance(eq, Expression) else (eq.lhs() - eq.rhs()) for eq in eqns]
    modulus = Integer(modulus)
    if modulus < 1:
        raise ValueError("the modulus must be a positive integer")
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort(key=repr)

    if modulus == 1:  # degenerate case
        ans = [tuple(Integers(1)(0) for v in vars)]
        return ans

    factors = modulus.factor()
    crt_basis = vector(Integers(modulus), crt_basis([p**i for p, i in factors]))
    solutions = []

    has_solution = True
    for p, i in factors:
        solution = _solve_mod_prime_power(eqns, p, i, vars)
        if len(solution) > 0:
            solutions.append(solution)
        else:
            has_solution = False
            break

    ans = []
    if has_solution:
        for solution in product(*solutions):
            solution_mat = matrix(Integers(modulus), solution)
            ans.append(tuple(c.dot_product(crt_basis) for c in solution_mat.columns()))

    # if solution_dict == True:
    # Relaxed form suggested by Mike Hansen (#8553):
    if solution_dict:
        return [dict(zip(vars, solution)) for solution in ans]
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

    - ``eqns`` -- equation or list of equations

    - ``p`` -- a prime

    - ``i`` -- integer > 0

    - ``vars`` -- list of variables to solve for

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
       is good, but the algorithm is horrible if the modulus is not the
       product of many small primes! Sage *does* have the ability to
       do something much faster in certain cases at least by using the
       Chinese Remainder Theorem, Groebner basis, linear algebra
       techniques, etc. But for a lot of toy problems this function as
       is might be useful. At the very least, it establishes an
       interface.

    TESTS:

    Confirm we can reproduce the first few terms of :oeis:`A187719`::

        sage: from sage.symbolic.relation import _solve_mod_prime_power
        sage: [sorted(_solve_mod_prime_power([x^2==41], 10, i, [x]))[0][0] for i in [1..13]]
        [1, 21, 71, 1179, 2429, 47571, 1296179, 8703821, 26452429, 526452429,
        13241296179, 19473547571, 2263241296179]
    """
    from sage.rings.finite_rings.integer_mod_ring import Integers
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.modules.free_module_element import vector

    mrunning = 1
    ans = []
    for mi in range(m):
        mrunning *= p
        R = Integers(mrunning)
        S = PolynomialRing(R, len(vars), vars)
        eqns_mod = [S(eq) for eq in eqns]
        if mi == 0:
            possibles = product(*[range(len(R)) for _ in range(len(vars))])
        else:
            shifts = product(*[range(p) for _ in range(len(vars))])
            pairs = product(shifts, ans)
            possibles = (tuple(vector(t) + vector(shift) * (mrunning // p))
                         for shift, t in pairs)
        ans = list(t for t in possibles if all(e(*t) == 0 for e in eqns_mod))
        if not ans:
            return ans

    return ans


def solve_ineq_univar(ineq):
    """
    Function solves rational inequality in one variable.

    INPUT:

    - ``ineq`` -- inequality in one variable

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

    Calls Maxima command ``solve_rat_ineq``

    AUTHORS:

    - Robert Marik (01-2010)
    """
    ineqvar = ineq.variables()
    if len(ineqvar) != 1:
        raise NotImplementedError("The command solve_ineq_univar accepts univariate inequalities only. Your variables are " + ineqvar)
    ineq0 = ineq._maxima_()
    ineq0.parent().eval("if solve_rat_ineq_loaded#true then (solve_rat_ineq_loaded:true,load(\"solve_rat_ineq.mac\")) ")
    sol = ineq0.solve_rat_ineq().sage()
    if repr(sol) == "all":
        from sage.rings.infinity import Infinity
        sol = [ineqvar[0] < Infinity]
    return sol


def solve_ineq_fourier(ineq, vars=None):
    """
    Solve system of inequalities using Maxima and Fourier elimination.

    Can be used for system of linear inequalities and for some types
    of nonlinear inequalities. For examples, see the example section
    below and http://maxima.cvs.sourceforge.net/viewvc/maxima/maxima/share/contrib/fourier_elim/rtest_fourier_elim.mac


    INPUT:

    - ``ineq`` -- list with system of inequalities

    - ``vars`` -- optionally list with variables for Fourier elimination

    OUTPUT:

    - ``list`` -- output is list of solutions as a list of simple inequalities
      output [A,B,C] means (A or B or C) each A, B, C is again a list and
      if A=[a,b], then A means (a and b). The list is empty if there is no
      solution.

    EXAMPLES::

        sage: from sage.symbolic.relation import solve_ineq_fourier
        sage: y = var('y')
        sage: solve_ineq_fourier([x+y<9,x-y>4],[x,y])
        [[y + 4 < x, x < -y + 9, y < (5/2)]]
        sage: solve_ineq_fourier([x+y<9,x-y>4],[y,x])[0][0](x=42)
        y < -33

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

    Calls Maxima command ``fourier_elim``

    AUTHORS:

    - Robert Marik (01-2010)
    """
    if vars is None:
        setvars = set([])
        for i in (ineq):
            setvars = setvars.union(set(i.variables()))
            vars = list(setvars)
    ineq0 = [i._maxima_() for i in ineq]
    ineq0[0].parent().eval("if fourier_elim_loaded#true then (fourier_elim_loaded:true,load(\"fourier_elim\"))")
    sol = ineq0[0].parent().fourier_elim(ineq0, vars)
    ineq0[0].parent().eval("or_to_list(x):=\
        if not atom(x) and op(x)=\"or\" then args(x) \
        else [x]")
    sol = sol.or_to_list().sage()
    if repr(sol) == "[emptyset]":
        sol = []
    if repr(sol) == "[universalset]":
        from sage.rings.infinity import Infinity
        sol = [[i < Infinity for i in vars]]
    return sol


def solve_ineq(ineq, vars=None):
    """
    Solve inequalities and systems of inequalities using Maxima.
    Switches between rational inequalities
    (sage.symbolic.relation.solve_ineq_rational)
    and Fourier elimination (sage.symbolic.relation.solve_ineq_fouried).
    See the documentation of these functions for more details.

    INPUT:

    - ``ineq`` -- one inequality or a list of inequalities

      Case1: If ``ineq`` is one equality, then it should be rational
      expression in one variable. This input is passed to
      sage.symbolic.relation.solve_ineq_univar function.

      Case2: If ``ineq`` is a list involving one or more
      inequalities, than the input is passed to
      sage.symbolic.relation.solve_ineq_fourier function. This
      function can be used for system of linear inequalities and
      for some types of nonlinear inequalities. See
      http://maxima.cvs.sourceforge.net/viewvc/maxima/maxima/share/contrib/fourier_elim/rtest_fourier_elim.mac
      for a big gallery of problems covered by this algorithm.

    - ``vars`` -- (optional) parameter with list of variables. This list
      is used only if Fourier elimination is used. If omitted or if
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

        sage: y = var('y')
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

    Calls ``solve_ineq_fourier`` if inequalities are list and
    ``solve_ineq_univar`` of the inequality is symbolic expression. See
    the description of these commands for more details related to the
    set of inequalities which can be solved. The list is empty if
    there is no solution.

    AUTHORS:

    - Robert Marik (01-2010)
    """
    if isinstance(ineq, list):
        return solve_ineq_fourier(ineq, vars)
    return solve_ineq_univar(ineq)
