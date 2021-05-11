r"""
Propositional Calculus

Formulas consist of the following operators:

* ``&`` -- and
* ``|`` -- or
* ``~`` -- not
* ``^`` -- xor
* ``->`` -- if-then
* ``<->`` -- if and only if

Operators can be applied to variables that consist of a leading letter and
trailing underscores and alphanumerics.  Parentheses may be used to explicitly
show order of operation.

AUTHORS:

- Chris Gorecki (2006): initial version, propcalc, boolformula,
  logictable, logicparser, booleval

- Michael Greenberg -- boolopt

- Paul Scurek (2013-08-05): updated docstring formatting

- Paul Scurek (2013-08-12): added :func:`~sage.logic.propcalc.get_formulas()`,
  :func:`~sage.logic.propcalc.consistent()`,
  :func:`~sage.logic.propcalc.valid_consequence()`

EXAMPLES:

We can create boolean formulas in different ways::

    sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
    sage: g = propcalc.formula("boolean<->algebra")
    sage: (f&~g).ifthen(f)
    ((a&((b|c)^a->c)<->b)&(~(boolean<->algebra)))->(a&((b|c)^a->c)<->b)

We can create a truth table from a formula::

    sage: f.truthtable()
    a      b      c      value
    False  False  False  True
    False  False  True   True
    False  True   False  False
    False  True   True   False
    True   False  False  True
    True   False  True   False
    True   True   False  True
    True   True   True   True
    sage: f.truthtable(end=3)
    a      b      c      value
    False  False  False  True
    False  False  True   True
    False  True   False  False
    sage: f.truthtable(start=4)
    a      b      c      value
    True   False  False  True
    True   False  True   False
    True   True   False  True
    True   True   True   True
    sage: propcalc.formula("a").truthtable()
    a      value
    False  False
    True   True

Now we can evaluate the formula for a given set of input::

    sage: f.evaluate({'a':True, 'b':False, 'c':True})
    False
    sage: f.evaluate({'a':False, 'b':False, 'c':True})
    True

And we can convert a boolean formula to conjunctive normal form::

    sage: f.convert_cnf_table()
    sage: f
    (a|~b|c)&(a|~b|~c)&(~a|b|~c)
    sage: f.convert_cnf_recur()
    sage: f
    (a|~b|c)&(a|~b|~c)&(~a|b|~c)

Or determine if an expression is satisfiable, a contradiction, or a tautology::

    sage: f = propcalc.formula("a|b")
    sage: f.is_satisfiable()
    True
    sage: f = f & ~f
    sage: f.is_satisfiable()
    False
    sage: f.is_contradiction()
    True
    sage: f = f | ~f
    sage: f.is_tautology()
    True

The equality operator compares semantic equivalence::

    sage: f = propcalc.formula("(a|b)&c")
    sage: g = propcalc.formula("c&(b|a)")
    sage: f == g
    True
    sage: g = propcalc.formula("a|b&c")
    sage: f == g
    False

It is an error to create a formula with bad syntax::

    sage: propcalc.formula("")
    Traceback (most recent call last):
    ...
    SyntaxError: malformed statement
    sage: propcalc.formula("a&b~(c|(d)")
    Traceback (most recent call last):
    ...
    SyntaxError: malformed statement
    sage: propcalc.formula("a&&b")
    Traceback (most recent call last):
    ...
    SyntaxError: malformed statement
    sage: propcalc.formula("a&b a")
    Traceback (most recent call last):
    ...
    SyntaxError: malformed statement

    It is also an error to not abide by the naming conventions.
    sage: propcalc.formula("~a&9b")
    Traceback (most recent call last):
    ...
    NameError: invalid variable name 9b: identifiers must begin with a letter and contain only alphanumerics and underscores
"""
# *****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Chris Gorecki <chris.k.gorecki@gmail.com>
#       Copyright (C) 2013 Paul Scurek <scurek86@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

### TODO:
### converts (cnf) returns w/o change

from . import boolformula
from . import logicparser
from sage.misc.superseded import deprecated_function_alias


def formula(s):
    r"""
    Return an instance of :class:`BooleanFormula`.

    INPUT:

    - ``s`` -- a string that contains a logical expression

    OUTPUT:

    An instance of :class:`BooleanFormula`.

    EXAMPLES:

    This example illustrates ways to create a boolean formula::

        sage: f = propcalc.formula("a&~b|c")
        sage: g = propcalc.formula("a^c<->b")
        sage: f&g|f
        ((a&~b|c)&(a^c<->b))|(a&~b|c)

    We now demonstrate some possible errors::

        sage: propcalc.formula("((a&b)")
        Traceback (most recent call last):
        ...
        SyntaxError: malformed statement
        sage: propcalc.formula("_a&b")
        Traceback (most recent call last):
        ...
        NameError: invalid variable name _a: identifiers must begin with a letter and contain only alphanumerics and underscores
    """
    try:
        parse_tree, vars_order = logicparser.parse(s)
        f = boolformula.BooleanFormula(s, parse_tree, vars_order)
        f.truthtable(0, 1)
    except (KeyError, RuntimeError, IndexError, SyntaxError):
        msg = "malformed statement"
        raise SyntaxError(msg)
    return f

def get_formulas(*statements):
    r"""
    Convert statements and parse trees into instances of
    :class:`BooleanFormula`.

    INPUT:

    - ``*statements`` -- strings or lists; a list must be a
      full syntax parse tree of a formula, and a string must
      be a string representation of a formula

    OUTPUT:

    The converted formulas in a list.

    EXAMPLES:

    This example illustrates converting strings into boolean formulas.

    ::

        sage: f = "a&(~c<->d)"
        sage: g = "d|~~b"
        sage: h = "~(a->c)<->(d|~c)"
        sage: propcalc.get_formulas(f, g, h)
        [a&(~c<->d), d|~~b, ~(a->c)<->(d|~c)]

    ::

        sage: A, B, C = propcalc.get_formulas("(a&b)->~c", "c", "~(a&b)")
        sage: A
        (a&b)->~c
        sage: B
        c
        sage: C
        ~(a&b)

    We can also convert parse trees into formulas.

    ::

        sage: t = ['a']
        sage: u = ['~', ['|', ['&', 'a', 'b'], ['~', 'c']]]
        sage: v = "b->(~c<->d)"
        sage: formulas= propcalc.get_formulas(t, u, v)
        sage: formulas[0]
        a
        sage: formulas[1]
        ~((a&b)|~c)
        sage: formulas[2]
        b->(~c<->d)

    AUTHORS:

    - Paul Scurek (2013-08-12)
    """
    formulas = []

    for statement in statements:
        try:
            if isinstance(statement, str):
                formulas.append(formula(statement))
            elif isinstance(statement, list):
                formulas.append(formula(logicparser.recover_formula(statement)))
            else:
                raise TypeError
        except (SyntaxError, NameError):
            raise SyntaxError("malformed statement")
        except TypeError:
            raise TypeError
    return formulas

def consistent(*formulas):
    r"""
    Determine if the formulas are logically consistent.

    INPUT:

    - ``*formulas`` -- instances of :class:`BooleanFormula`

    OUTPUT:

    A boolean value to be determined as follows:

    - ``True`` - if the formulas are logically consistent

    - ``False`` - if the formulas are not logically consistent

    EXAMPLES:

    This example illustrates determining if formulas are logically consistent.

    ::

        sage: f, g, h, i = propcalc.get_formulas("a<->b", "~b->~c", "d|g", "c&a")
        sage: propcalc.consistent(f, g, h, i)
        True

    ::

        sage: j, k, l, m = propcalc.get_formulas("a<->b", "~b->~c", "d|g", "c&~a")
        sage: propcalc.consistent(j, k ,l, m)
        False

    AUTHORS:

    - Paul Scurek (2013-08-12)
    """
    # make sure only instances of :class:`BooleanFormula` were passed as arguments
    for formula in formulas[1:]:
        if not isinstance(formula, boolformula.BooleanFormula):
            raise TypeError("consistent() takes BooleanFormula() class instances as arguments")

    # conjoin all of the formulas with &
    conjunction = formulas[0]
    for formula in formulas[1:]:
        conjunction = conjunction & formula

    # if conjunction is a contradiction, the formulas are inconsistent
    return not conjunction.is_contradiction()

# define function ``valid_consequence`` for backward compatibility
valid_consequence = deprecated_function_alias(28052, boolformula.is_consequence)
