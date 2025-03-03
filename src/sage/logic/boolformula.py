r"""
Boolean Formulas

Formulas consist of the operators ``&``, ``|``, ``~``, ``^``, ``->``, ``<->``,
corresponding to ``and``, ``or``, ``not``, ``xor``, ``if...then``, ``if and
only if``.  Operators can be applied to variables that consist of a leading
letter and trailing underscores and alphanumerics.  Parentheses may be used
to explicitly show order of operation.

EXAMPLES:

Create boolean formulas and combine them with
:meth:`~sage.logic.boolformula.BooleanFormula.ifthen()` method::

    sage: import sage.logic.propcalc as propcalc
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

Now we can evaluate the formula for a given set of inputs::

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

It is also an error to not abide by the naming conventions::

    sage: propcalc.formula("~a&9b")
    Traceback (most recent call last):
    ...
    NameError: invalid variable name 9b: identifiers must begin with a letter and contain only alphanumerics and underscores

AUTHORS:

- Chris Gorecki (2006): initial version

- Paul Scurek (2013-08-03): added polish_notation, full_tree,
  updated docstring formatting

- Paul Scurek (2013-08-08): added
  :meth:`~sage.logic.boolformula.BooleanFormula.implies()`
"""
# *****************************************************************************
#       Copyright (C) 2006 William Stein <wstein.gmail.com>
#       Copyright (C) 2006 Chris Gorecki <chris.k.gorecki@gmail.com>
#       Copyright (C) 2013 Paul Scurek <scurek86@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import booleval
from . import logictable
from . import logicparser
# import boolopt
from sage.misc.flatten import flatten

latex_operators = [('&', '\\wedge '),
                   ('|', '\\vee '),
                   ('~', '\\neg '),
                   ('^', '\\oplus '),
                   ('<->', '\\leftrightarrow '),
                   ('->', '\\rightarrow ')]


class BooleanFormula:
    """
    Boolean formulas.

    INPUT:

    - ``self`` -- calling object

    - ``exp`` -- string; this contains the boolean expression
      to be manipulated

    - ``tree`` -- list; this contains the parse tree of the expression

    - ``vo`` -- list; this contains the variables in the expression, in the
      order that they appear; each variable only occurs once in the list
    """
    __expression = ""
    __tree = []
    __vars_order = []

    def __init__(self, exp, tree, vo):
        r"""
        Initialize the data fields.

        EXAMPLES:

        This example illustrates the creation of a statement::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b|~(c|a)")
            sage: s
            a&b|~(c|a)
        """
        self.__expression = exp.replace(' ', '')
        self.__tree = tree
        self.__vars_order = vo

    def __repr__(self):
        r"""
        Return a string representation of this statement.

        EXAMPLES::

            sage: import sage.logic.propcalc as propcalc
            sage: propcalc.formula("man->monkey&human")
            man->monkey&human
        """
        return self.__expression

    def _latex_(self):
        r"""
        Return a LaTeX representation of this statement.

        OUTPUT: string containing the latex code for the statement

        EXAMPLES::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("man->monkey&human")
            sage: latex(s)
            man\rightarrow monkey\wedge human

            sage: f = propcalc.formula("a & ((~b | c) ^ a -> c) <-> ~b")
            sage: latex(f)
            a\wedge ((\neg b\vee c)\oplus a\rightarrow c)\leftrightarrow \neg b
        """
        latex_expression = self.__expression
        for old, new in latex_operators:
            latex_expression = latex_expression.replace(old, new)
        return latex_expression

    def polish_notation(self):
        r"""
        Convert the calling boolean formula into polish notation.

        OUTPUT: string representation of the formula in polish notation

        EXAMPLES:

        This example illustrates converting a formula to polish notation::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("~~a|(c->b)")
            sage: f.polish_notation()
            '|~~a->cb'

            sage: g = propcalc.formula("(a|~b)->c")
            sage: g.polish_notation()
            '->|a~bc'

        AUTHORS:

        - Paul Scurek (2013-08-03)
        """
        return ''.join(flatten(logicparser.polish_parse(repr(self))))

    def tree(self):
        r"""
        Return the parse tree of this boolean expression.

        OUTPUT: the parse tree as a nested list

        EXAMPLES:

        This example illustrates how to find the parse tree of a boolean
        formula::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("man -> monkey & human")
            sage: s.tree()
            ['->', 'man', ['&', 'monkey', 'human']]

        ::

            sage: f = propcalc.formula("a & ((~b | c) ^ a -> c) <-> ~b")
            sage: f.tree()
            ['<->',
             ['&', 'a', ['->', ['^', ['|', ['~', 'b', None], 'c'], 'a'], 'c']],
             ['~', 'b', None]]

        .. NOTE::

            This function is used by other functions in the logic module
            that perform semantic operations on a boolean formula.
        """
        return self.__tree

    def full_tree(self):
        r"""
        Return a full syntax parse tree of the calling formula.

        OUTPUT: the full syntax parse tree as a nested list

        EXAMPLES:

        This example shows how to find the full syntax parse tree
        of a formula::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a->(b&c)")
            sage: s.full_tree()
            ['->', 'a', ['&', 'b', 'c']]

            sage: t = propcalc.formula("a & ((~b | c) ^ a -> c) <-> ~b")
            sage: t.full_tree()
            ['<->', ['&', 'a', ['->', ['^', ['|', ['~', 'b'], 'c'], 'a'], 'c']], ['~', 'b']]

            sage: f = propcalc.formula("~~(a&~b)")
            sage: f.full_tree()
            ['~', ['~', ['&', 'a', ['~', 'b']]]]

        .. NOTE::

            This function is used by other functions in the logic module
            that perform syntactic operations on a boolean formula.

        AUTHORS:

        - Paul Scurek (2013-08-03)
        """
        return logicparser.polish_parse(repr(self))

    def __or__(self, other):
        r"""
        Overload the ``|`` operator to 'or' two statements together.

        INPUT:

        - ``other`` -- boolean formula; this is the statement
          on the right side of the operator

        OUTPUT:

        A boolean formula of the form ``self | other``.

        EXAMPLES:

        This example illustrates combining two formulas with ``|``::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s | f
            (a&b)|(c^d)
        """
        return self.add_statement(other, '|')

    def __and__(self, other):
        r"""
        Overload the ``&`` operator to 'and' two statements together.

        INPUT:

        - ``other`` -- boolean formula; this is the formula on
          the right side of the operator

        OUTPUT: a boolean formula of the form ``self & other``

        EXAMPLES:

        This example shows how to combine two formulas with ``&``::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s & f
            (a&b)&(c^d)
        """
        return self.add_statement(other, '&')

    def __xor__(self, other):
        r"""
        Overload the ``^`` operator to 'xor' two statements together.

        INPUT:

        - ``other`` -- boolean formula; this is the formula on
          the right side of the operator

        OUTPUT: a boolean formula of the form ``self ^ other``

        EXAMPLES:

        This example illustrates how to combine two formulas with ``^``::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s ^ f
            (a&b)^(c^d)
        """
        return self.add_statement(other, '^')

    def __pow__(self, other):
        r"""
        Overload the ``^`` operator to 'xor' two statements together.

        INPUT:

        - ``other`` -- boolean formula; this is the formula on
          the right side of the operator

        OUTPUT: a boolean formula of the form ``self ^ other``

        EXAMPLES:

        This example shows how to combine two formulas with ``^``::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s ^ f
            (a&b)^(c^d)

        .. TODO::

            This function seems to be identical to ``__xor__``.
            Thus, this function should be replaced with ``__xor__`` everywhere
            that it appears in the logic module. Then it can be deleted
            altogether.
        """
        return self.add_statement(other, '^')

    def __invert__(self):
        r"""
        Overload the ``~`` operator to 'not' a statement.

        OUTPUT: a boolean formula of the form ``~self``

        EXAMPLES:

        This example shows how to negate a boolean formula::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: ~s
            ~(a&b)
        """
        exp = '~(' + self.__expression + ')'
        parse_tree, vars_order = logicparser.parse(exp)
        return BooleanFormula(exp, parse_tree, vars_order)

    def ifthen(self, other):
        r"""
        Combine two formulas with the ``->`` operator.

        INPUT:

        - ``other`` -- boolean formula; this is the formula
          on the right side of the operator

        OUTPUT:

        A boolean formula of the form ``self -> other``.

        EXAMPLES:

        This example illustrates how to combine two formulas with '->'::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s.ifthen(f)
            (a&b)->(c^d)
        """
        return self.add_statement(other, '->')

    def iff(self, other):
        r"""
        Combine two formulas with the ``<->`` operator.

        INPUT:

        - ``other`` -- boolean formula; this is the formula
          on the right side of the operator

        OUTPUT:

        A boolean formula of the form ``self <-> other``.

        EXAMPLES:

        This example illustrates how to combine two formulas with '<->'::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s.iff(f)
            (a&b)<->(c^d)
        """
        return self.add_statement(other, '<->')

    def __eq__(self, other):
        r"""
        Overload the ``==`` operator to determine logical equivalence.

        INPUT:

        - ``other`` -- boolean formula; this is the formula
          on the right side of the comparator

        OUTPUT: a boolean value to be determined as follows:

        - ``True`` if ``self`` and ``other`` are logically equivalent

        - ``False`` if ``self`` and ``other`` are not logically equivalent

        EXAMPLES:

        This example shows how to determine logical equivalence::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("(a|b)&c")
            sage: g = propcalc.formula("c&(b|a)")
            sage: f == g
            True

        ::

            sage: g = propcalc.formula("a|b&c")
            sage: f == g
            False
        """
        return self.equivalent(other)

    def truthtable(self, start=0, end=-1):
        r"""
        Return a truth table for the calling formula.

        INPUT:

        - ``start`` -- (default: 0) an integer; this is the first
          row of the truth table to be created

        - ``end`` -- (default: -1) an integer; this is the last
          row of the truth table to be created

        OUTPUT: the truth table as a 2-D array

        EXAMPLES:

        This example illustrates the creation of a truth table::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b|~(c|a)")
            sage: s.truthtable()
            a      b      c      value
            False  False  False  True
            False  False  True   False
            False  True   False  True
            False  True   True   False
            True   False  False  False
            True   False  True   False
            True   True   False  True
            True   True   True   True

        We can now create a truthtable of rows 1 to 4, inclusive::

            sage: s.truthtable(1, 5)
            a      b      c      value
            False  False  True   False
            False  True   False  True
            False  True   True   False
            True   False  False  False

        .. NOTE::

            Each row of the table corresponds to a binary number, with
            each variable associated to a column of the number, and taking on
            a true value if that column has a value of 1.  Please see the
            logictable module for details.  The function returns a table that
            start inclusive and end exclusive so ``truthtable(0, 2)`` will
            include row 0, but not row 2.

            When sent with no start or end parameters, this is an
            exponential time function requiring `O(2^n)` time, where
            `n` is the number of variables in the expression.
        """
        maximum = 2 ** len(self.__vars_order)
        if end < 0:
            end = maximum
        end = min(end, maximum)
        start = max(start, 0)
        start = min(start, maximum)
        keys, table = [], []
        vars = {}
        for var in self.__vars_order:
            vars[var] = False
            keys.insert(0, var)
        keys = list(keys)
        for i in range(start, end):
            j = 0
            row = []
            for key in keys:
                bit = self.get_bit(i, j)
                vars[key] = bit
                j += 1
                row.insert(0, bit)
            row.append(booleval.eval_formula(self.__tree, vars))
            table.append(row)
        keys.reverse()
        table = logictable.Truthtable(table, keys)
        return table

    def evaluate(self, var_values):
        r"""
        Evaluate a formula for the given input values.

        INPUT:

        - ``var_values`` -- dictionary; this contains the
          pairs of variables and their boolean values

        OUTPUT: the result of the evaluation as a boolean

        EXAMPLES:

        This example illustrates the evaluation of a boolean formula::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&b|c")
            sage: f.evaluate({'a':False, 'b':False, 'c':True})
            True
            sage: f.evaluate({'a':True, 'b':False, 'c':False})
            False
        """
        return booleval.eval_formula(self.__tree, var_values)

    def is_satisfiable(self):
        r"""
        Determine if the formula is ``True`` for some assignment of values.

        OUTPUT: a boolean value to be determined as follows:

        - ``True`` if there is an assignment of values that makes the
          formula ``True``.

        - ``False`` if the formula cannot be made ``True`` by any assignment
          of values.

        EXAMPLES:

        This example illustrates how to check a formula for satisfiability::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a|b")
            sage: f.is_satisfiable()
            True

            sage: g = f & (~f)
            sage: g.is_satisfiable()
            False
        """
        table = self.truthtable().get_table_list()
        return any(row[-1] is True for row in table[1:])

    def is_tautology(self):
        r"""
        Determine if the formula is always ``True``.

        OUTPUT: a boolean value to be determined as follows:

        - ``True`` if the formula is a tautology.

        - ``False`` if the formula is not a tautology.

        EXAMPLES:

        This example illustrates how to check if a formula is a tautology::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a|~a")
            sage: f.is_tautology()
            True

            sage: f = propcalc.formula("a&~a")
            sage: f.is_tautology()
            False

            sage: f = propcalc.formula("a&b")
            sage: f.is_tautology()
            False
        """
        return not (~self).is_satisfiable()

    def is_contradiction(self):
        r"""
        Determine if the formula is always ``False``.

        OUTPUT: a boolean value to be determined as follows:

        - ``True`` if the formula is a contradiction.

        - ``False`` if the formula is not a contradiction.

        EXAMPLES:

        This example illustrates how to check if a formula is a contradiction.

        ::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&~a")
            sage: f.is_contradiction()
            True

            sage: f = propcalc.formula("a|~a")
            sage: f.is_contradiction()
            False

            sage: f = propcalc.formula("a|b")
            sage: f.is_contradiction()
            False
        """
        return not self.is_satisfiable()

    def is_consequence(self, *hypotheses):
        r"""
        Determine if ``self`` (the desired conclusion) is a logical consequence of the
        hypotheses. The function call ``is_consequence(conclusion, *hypotheses)`` is a
        synonym for ``conclusion.is_consequence(*hypotheses)``.

        INPUT:

        - ``*hypotheses`` -- instances of :class:`BooleanFormula`

        OUTPUT: a boolean value to be determined as follows:

        - ``True`` -- if ``self`` (the desired conclusion) is a logical consequence
          of the set of hypotheses

        - ``False`` -- if ``self`` (the desired conclusion) is not a logical consequence
          of the set of hypotheses

        EXAMPLES::

            sage: from sage.logic.propcalc import formula
            sage: formula("a | b").is_consequence(formula("b"))
            True
            sage: formula("a & b").is_consequence(formula("b"))
            False
            sage: formula("b").is_consequence(formula("a"), formula("a -> b"))
            True
            sage: formula("b -> a").is_consequence(formula("a -> b"))
            False
            sage: formula("~b -> ~a").is_consequence(formula("a -> b"))
            True

        ::

            sage: f, g, h = propcalc.get_formulas("a & ~b", "c -> b", "c | e")
            sage: propcalc.formula("a & e").is_consequence(f, g, h)
            True
            sage: i = propcalc.formula("a & ~e")
            sage: i.is_consequence(f, g, h)
            False
            sage: from sage.logic.boolformula import is_consequence
            sage: is_consequence(i, f, g, h)
            False
            sage: is_consequence(propcalc.formula("((p <-> q) & r) -> ~c"), f, g, h)
            True

        Only a tautology is a logical consequence of an empty set of formulas::

            sage: propcalc.formula("a | ~a").is_consequence()
            True
            sage: propcalc.formula("a | b").is_consequence()
            False

        TESTS:

        Arguments must be instances of :class:`BooleanFormula` (not strings, for example)::

            sage: propcalc.formula("a | b").is_consequence("a | b")
            Traceback (most recent call last):
            ...
            TypeError: is_consequence only takes instances of BooleanFormula() class as input

        AUTHORS:

        - Paul Scurek (2013-08-12)
        """
        # make sure every argument is an instance of :class:`BooleanFormula`
        for formula in (self,) + hypotheses:
            if not isinstance(formula, BooleanFormula):
                raise TypeError("is_consequence only takes instances of BooleanFormula() class as input")

        if not hypotheses:
            # if there are no hypotheses, then we just want to know whether self is a tautology
            return self.is_tautology()
        else:
            # conjoin all of the hypotheses into a single Boolean formula
            conjunction = hypotheses[0]
            for hypothesis in hypotheses[1:]:
                conjunction = conjunction & hypothesis

            return conjunction.implies(self)

    def implies(self, other):
        r"""
        Determine if calling formula implies other formula.

        INPUT:

        - ``self`` -- calling object

        - ``other`` -- instance of :class:`BooleanFormula`

        OUTPUT: a boolean value to be determined as follows:

        - ``True`` -- if ``self`` implies ``other``

        - ``False`` -- if ``self does not imply ``other``

        EXAMPLES:

        This example illustrates determining if one formula implies another::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a<->b")
            sage: g = propcalc.formula("b->a")
            sage: f.implies(g)
            True

        ::

            sage: h = propcalc.formula("a->(a|~b)")
            sage: i = propcalc.formula("a")
            sage: h.implies(i)
            False

        AUTHORS:

        - Paul Scurek (2013-08-08)
        """
        # input validation
        if not isinstance(other, BooleanFormula):
            raise TypeError("implies() takes an instance of the BooleanFormula() class as input")

        conditional = self.ifthen(other)
        return (conditional).is_tautology()

    def equivalent(self, other):
        r"""
        Determine if two formulas are semantically equivalent.

        INPUT:

        - ``self`` -- calling object

        - ``other`` -- instance of BooleanFormula class

        OUTPUT: a boolean value to be determined as follows:

        ``True`` -- if the two formulas are logically equivalent

        ``False`` -- if the two formulas are not logically equivalent

        EXAMPLES:

        This example shows how to check for logical equivalence::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("(a|b)&c")
            sage: g = propcalc.formula("c&(a|b)")
            sage: f.equivalent(g)
            True

            sage: g = propcalc.formula("a|b&c")
            sage: f.equivalent(g)
            False
        """
        return self.iff(other).is_tautology()

    def convert_cnf_table(self):
        r"""
        Convert boolean formula to conjunctive normal form.

        OUTPUT: an instance of :class:`BooleanFormula` in conjunctive normal form

        EXAMPLES:

        This example illustrates how to convert a formula to cnf::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a ^ b <-> c")
            sage: s.convert_cnf()
            sage: s
            (a|b|~c)&(a|~b|c)&(~a|b|c)&(~a|~b|~c)

        We now show that :meth:`convert_cnf` and :meth:`convert_cnf_table`
        are aliases::

            sage: t = propcalc.formula("a ^ b <-> c")
            sage: t.convert_cnf_table(); t
            (a|b|~c)&(a|~b|c)&(~a|b|c)&(~a|~b|~c)
            sage: t == s
            True

        .. NOTE::

            This method creates the cnf parse tree by examining the logic
            table of the formula.  Creating the table requires `O(2^n)` time
            where `n` is the number of variables in the formula.
        """
        str = ''
        t = self.truthtable()
        table = t.get_table_list()
        vars = table[0]
        for row in table[1:]:
            if row[-1] is False:
                str += '('
                for i in range(len(row) - 1):
                    if row[i] is True:
                        str += '~'
                    str += vars[i]
                    str += '|'
                str = str[:-1] + ')&'
        self.__expression = str[:-1]
        # in case of tautology
        if len(self.__expression) == 0:
            self.__expression = '(' + self.__vars_order[0] + '|~' + self.__vars_order[0] + ')'
        self.__tree, self.__vars_order = logicparser.parse(self.__expression)

    convert_cnf = convert_cnf_table

    def convert_cnf_recur(self):
        r"""
        Convert boolean formula to conjunctive normal form.

        OUTPUT: an instance of :class:`BooleanFormula` in conjunctive normal form

        EXAMPLES:

        This example hows how to convert a formula to conjunctive normal form::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a^b<->c")
            sage: s.convert_cnf_recur()
            sage: s
            (~a|a|c)&(~b|a|c)&(~a|b|c)&(~b|b|c)&(~c|a|b)&(~c|~a|~b)

        .. NOTE::

            This function works by applying a set of rules that are
            guaranteed to convert the formula.  Worst case the converted
            expression has an `O(2^n)` increase in size (and time as well), but
            if the formula is already in CNF (or close to) it is only `O(n)`.

            This function can require an exponential blow up in space from the
            original expression.  This in turn can require large amounts of
            time. Unless a formula is already in (or close to) being in cnf
            :meth:`convert_cnf()` is typically preferred, but results can vary.
        """
        self.__tree = logicparser.apply_func(self.__tree, self.reduce_op)
        self.__tree = logicparser.apply_func(self.__tree, self.dist_not)
        self.__tree = logicparser.apply_func(self.__tree, self.dist_ors)
        self.convert_expression()

    def satformat(self):
        r"""
        Return the satformat representation of a boolean formula.

        OUTPUT: the satformat of the formula as a string

        EXAMPLES:

        This example illustrates how to find the satformat of a formula::

            sage: import sage.logic.propcalc as propcalc
            sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
            sage: f.convert_cnf()
            sage: f
            (a|~b|c)&(a|~b|~c)&(~a|b|~c)
            sage: f.satformat()
            'p cnf 3 0\n1 -2 3 0 1 -2 -3 \n0 -1 2 -3'

        .. NOTE::

            See www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps for a
            description of satformat.

            If the instance of boolean formula has not been converted to
            CNF form by a call to :meth:`convert_cnf()` or
            :meth:`convert_cnf_recur()`, then :meth:`satformat()` will call
            :meth:`convert_cnf()`. Please see the notes for
            :meth:`convert_cnf()` and :meth:`convert_cnf_recur()` for
            performance issues.
        """
        self.convert_cnf_table()
        s = ''
        vars_num = {}
        i = 0
        clauses = 0
        for e in self.__vars_order:
            vars_num[e] = str(i + 1)
            i += 1
        i = 0
        w = 1
        while i < len(self.__expression):
            c = self.__expression[i]
            if c == ')':
                clauses += 1
            if c in '()|':
                i += 1
                continue
            if c == '~':
                s += '-'
            elif c == '&':
                s += '0 '
            else:
                varname = ''
                if self.__expression[i] not in '|) ':
                    varname += self.__expression[i]
                    i += 1
                s += vars_num[varname] + ' '
            if len(s) >= (w * 15) and s[-1] != '-':
                s += '\n'
                w += 1
            i += 1
        s = 'p cnf ' + str(len(self.__vars_order)) + ' ' + str(clauses) + '\n' + s
        return s[:-1]

#    def simplify(self):
#        r"""
#        This function uses the propcalc package to simplify an expression to
#        its minimal form.
#
#        OUTPUT: a simplified expression
#
#        EXAMPLES::

#            sage: import sage.logic.propcalc as propcalc
#            sage: f = propcalc.formula("a&((b|c)^a->c)<->b")
#            sage: f.truthtable()
#            a      b      c      value
#            False  False  False  True
#            False  False  True   True
#            False  True   False  False
#            False  True   True   False
#            True   False  False  True
#            True   False  True   False
#            True   True   False  True
#            True   True   True   True
#            sage: f.simplify()
#            (~a&~b)|(a&~b&~c)|(a&b)
#            sage: f.truthtable()
#            a      b      c      value
#            False  False  False  True
#            False  False  True   True
#            False  True   False  False
#            False  True   True   False
#            True   False  False  True
#            True   False  True   False
#            True   True   False  True
#            True   True   True   True
#
#        .. NOTE::
#
#            If the instance of boolean formula has not been converted to
#            cnf form by a call to convert_cnf() or convert_cnf_recur()
#            satformat() will call convert_cnf().  Please see the notes for
#            convert_cnf() and convert_cnf_recur() for performance issues.
#        """
#        exp = ''
#        self.__tree = logicparser.apply_func(self.__tree, self.reduce_op)
#        plf = logicparser.apply_func(self.__tree, self.convert_opt)
#        wff = boolopt.PLFtoWFF()(plf) # convert to positive-normal form
#        wtd = boolopt.WFFtoDNF()
#        dnf = wtd(wff)
#        dnf = wtd.clean(dnf)
#        if dnf == [] or dnf == [[]]:
#            exp = self.__vars_order[0] + '&~' + self.__vars_order[0] + ' '
#        opt = boolopt.optimize(dnf)
#        if exp == '' and (opt == [] or opt == [[]]):
#            exp = self.__vars_order[0] + '|~' + self.__vars_order[0] + ' '
#        if exp == '':
#            for con in opt:
#                s = '('
#                for prop in con:
#                    if prop[0] == 'notprop':
#                       s += '~'
#                    s += prop[1] + '&'
#                exp += s[:-1] + ')|'
#        self.__expression = exp[:-1]
#        self.__tree, self.__vars_order = logicparser.parse(self.__expression)
#        return BooleanFormula(self.__expression, self.__tree, self.__vars_order)

    def convert_opt(self, tree):
        r"""
        Convert a parse tree to the tuple form used by :meth:`bool_opt()`.

        INPUT:

        - ``tree`` -- list; this is a branch of a
          parse tree and can only contain the '&', '|'
          and '~' operators along with variables

        OUTPUT: a 3-tuple

        EXAMPLES:

        This example illustrates the conversion of a formula into its
        corresponding tuple::

            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("a&(b|~c)")
            sage: tree = ['&', 'a', ['|', 'b', ['~', 'c', None]]]
            sage: logicparser.apply_func(tree, s.convert_opt)
            ('and', ('prop', 'a'), ('or', ('prop', 'b'), ('not', ('prop', 'c'))))

        .. NOTE::

            This function only works on one branch of the parse tree. To
            apply the function to every branch of a parse tree, pass the
            function as an argument in
            :func:`~sage.logic.logicparser.apply_func()` in
            :mod:`~sage.logic.logicparser`.
        """
        if not isinstance(tree[1], tuple) and tree[1] is not None:
            lval = ('prop', tree[1])
        else:
            lval = tree[1]
        if not isinstance(tree[2], tuple) and tree[2] is not None:
            rval = ('prop', tree[2])
        else:
            rval = tree[2]
        if tree[0] == '~':
            return ('not', lval)
        if tree[0] == '&':
            op = 'and'
        if tree[0] == '|':
            op = 'or'
        return (op, lval, rval)

    def add_statement(self, other, op):
        r"""
        Combine two formulas with the given operator.

        INPUT:

        - ``other`` -- instance of :class:`BooleanFormula`; this
          is the formula on the right of the operator

        - ``op`` -- string; this is the operator used to
          combine the two formulas

        OUTPUT: the result as an instance of :class:`BooleanFormula`

        EXAMPLES:

        This example shows how to create a new formula from two others::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: f = propcalc.formula("c^d")
            sage: s.add_statement(f, '|')
            (a&b)|(c^d)

            sage: s.add_statement(f, '->')
            (a&b)->(c^d)
        """
        exp = '(' + self.__expression + ')' + op + '(' + other.__expression + ')'
        parse_tree, vars_order = logicparser.parse(exp)
        return BooleanFormula(exp, parse_tree, vars_order)

    def get_bit(self, x, c):
        r"""
        Determine if bit ``c`` of the number ``x`` is 1.

        INPUT:

        - ``x`` -- integer; this is the number from
          which to take the bit

        - ``c`` -- integer; this is the but number to
          be taken, where 0 is the low order bit

        OUTPUT: a boolean to be determined as follows:

        - ``True`` if bit ``c`` of ``x`` is 1.

        - ``False`` if bit c of x is not 1.

        EXAMPLES:

        This example illustrates the use of :meth:`get_bit`::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a&b")
            sage: s.get_bit(2, 1)
            True
            sage: s.get_bit(8, 0)
            False

        It is not an error to have a bit out of range::

            sage: s.get_bit(64, 7)
            False

        Nor is it an error to use a negative number::

            sage: s.get_bit(-1, 3)
            False
            sage: s.get_bit(64, -1)
            True
            sage: s.get_bit(64, -2)
            False

        .. NOTE::

            The 0 bit is the low order bit.  Errors should be handled
            gracefully by a return of ``False``, and negative numbers ``x``
            always return ``False`` while a negative ``c`` will index from the
            high order bit.
        """
        bits = []
        while x > 0:
            b = bool(x % 2)
            x = x // 2
            bits.append(b)
        if c > len(bits) - 1:
            return False
        else:
            return bits[c]

    def reduce_op(self, tree):
        r"""
        Convert if-and-only-if, if-then, and xor operations to operations
        only involving and/or operations.

        INPUT:

        - ``tree`` -- list; this represents a branch
          of a parse tree

        OUTPUT:

        A new list with no '^', '->', or '<->' as first element of list.

        EXAMPLES:

        This example illustrates the use of :meth:`reduce_op` with
        :func:`apply_func`::

            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("a->b^c")
            sage: tree = ['->', 'a', ['^', 'b', 'c']]
            sage: logicparser.apply_func(tree, s.reduce_op)
            ['|', ['~', 'a', None], ['&', ['|', 'b', 'c'], ['~', ['&', 'b', 'c'], None]]]

        .. NOTE::

            This function only operates on a single branch of a parse tree.
            To apply the function to an entire parse tree, pass the function
            as an argument to :func:`~sage.logic.logicparser.apply_func()`
            in :mod:`~sage.logic.logicparser`.
        """
        if tree[0] == '<->':
            # parse tree for (~tree[1]|tree[2])&(~tree[2]|tree[1])
            new_tree = ['&', ['|', ['~', tree[1], None], tree[2]],
                        ['|', ['~', tree[2], None], tree[1]]]
        elif tree[0] == '^':
            # parse tree for (tree[1]|tree[2])&~(tree[1]&tree[2])
            new_tree = ['&', ['|', tree[1], tree[2]],
                        ['~', ['&', tree[1], tree[2]], None]]
        elif tree[0] == '->':
            # parse tree for ~tree[1]|tree[2]
            new_tree = ['|', ['~', tree[1], None], tree[2]]
        else:
            new_tree = tree
        return new_tree

    def dist_not(self, tree):
        r"""
        Distribute '~' operators over '&' and '|' operators.

        INPUT:

        - ``tree`` a list; this represents a branch
          of a parse tree

        OUTPUT: a new list

        EXAMPLES:

        This example illustrates the distribution of '~' over '&'::

            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("~(a&b)")
            sage: tree = ['~', ['&', 'a', 'b'], None]
            sage: logicparser.apply_func(tree, s.dist_not) #long time
            ['|', ['~', 'a', None], ['~', 'b', None]]

        .. NOTE::

            This function only operates on a single branch of a parse tree.
            To apply the function to an entire parse tree, pass the function
            as an argument to :func:`~sage.logic.logicparser.apply_func()`
            in :mod:`~sage.logic.logicparser`.
        """
        if tree[0] == '~' and isinstance(tree[1], list):
            op = tree[1][0]
            if op != '~':
                op = '|' if op == '&' else '&'
                new_tree = [op, ['~', tree[1][1], None], ['~', tree[1][2], None]]
                return logicparser.apply_func(new_tree, self.dist_not)
            else:
                # cancel double negative
                return tree[1][1]
        else:
            return tree

    def dist_ors(self, tree):
        r"""
        Distribute '|' over '&'.

        INPUT:

        - ``tree`` -- list; this represents a branch of
          a parse tree

        OUTPUT: a new list

        EXAMPLES:

        This example illustrates the distribution of '|' over '&'::

            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("(a&b)|(a&c)")
            sage: tree = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
            sage: logicparser.apply_func(tree, s.dist_ors) #long time
            ['&', ['&', ['|', 'a', 'a'], ['|', 'b', 'a']], ['&', ['|', 'a', 'c'], ['|', 'b', 'c']]]

        .. NOTE::

            This function only operates on a single branch of a parse tree.
            To apply the function to an entire parse tree, pass the function
            as an argument to :func:`~sage.logic.logicparser.apply_func()`
            in :mod:`~sage.logic.logicparser`.
        """
        if tree[0] == '|' and isinstance(tree[2], list) and tree[2][0] == '&':
            new_tree = ['&', ['|', tree[1], tree[2][1]],
                        ['|', tree[1], tree[2][2]]]
            return logicparser.apply_func(new_tree, self.dist_ors)
        if tree[0] == '|' and isinstance(tree[1], list) and tree[1][0] == '&':
            new_tree = ['&', ['|', tree[1][1], tree[2]],
                        ['|', tree[1][2], tree[2]]]
            return logicparser.apply_func(new_tree, self.dist_ors)
        return tree

    def to_infix(self, tree):
        r"""
        Convert a parse tree from prefix to infix form.

        INPUT:

        - ``tree`` -- list; this represents a branch
          of a parse tree

        OUTPUT: a new list

        EXAMPLES:

        This example shows how to convert a parse tree from prefix to
        infix form::

            sage: import sage.logic.propcalc as propcalc, sage.logic.logicparser as logicparser
            sage: s = propcalc.formula("(a&b)|(a&c)")
            sage: tree = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
            sage: logicparser.apply_func(tree, s.to_infix)
            [['a', '&', 'b'], '|', ['a', '&', 'c']]

        .. NOTE::

            This function only operates on a single branch of a parse tree.
            To apply the function to an entire parse tree, pass the function
            as an argument to :func:`~sage.logic.logicparser.apply_func()`
            in :mod:`~sage.logic.logicparser`.
        """
        if tree[0] != '~':
            return [tree[1], tree[0], tree[2]]
        return tree

    def convert_expression(self):
        r"""
        Convert the string representation of a formula to conjunctive
        normal form.

        EXAMPLES::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a^b<->c")
            sage: s.convert_expression(); s
            a^b<->c
        """
        ttree = self.__tree[:]
        ttree = logicparser.apply_func(ttree, self.to_infix)
        self.__expression = ''
        str_tree = str(ttree)
        open_flag = False
        i = 0
        for c in str_tree:
            if i < len(str_tree) - 1:
                op = self.get_next_op(str_tree[i:])
            if op == '|' and not open_flag:
                self.__expression += '('
                open_flag = True
            if i < len(str_tree) - 2 and str_tree[i + 1] == '&' and open_flag:
                open_flag = False
                self.__expression += ')'
            if str_tree[i:i + 4] == 'None':
                i += 4
            if i < len(str_tree) and str_tree[i] not in ' \',[]':
                self.__expression += str_tree[i]
            i += 1
        if open_flag is True:
            self.__expression += ')'

    def get_next_op(self, str):
        r"""
        Return the next operator in a string.

        INPUT:

        - ``str`` -- string; this contains a logical
          expression

        OUTPUT: the next operator as a string

        EXAMPLES:

        This example illustrates how to find the next operator in a formula::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("f&p")
            sage: s.get_next_op("abra|cadabra")
            '|'

        .. NOTE::

            The parameter ``str`` is not necessarily the string
            representation of the calling object.
        """
        i = 0
        while i < len(str) - 1 and str[i] != '&' and str[i] != '|':
            i += 1
        return str[i]

    def length(self):
        r"""
        Return the length of ``self``.

        OUTPUT:

        The length of the Boolean formula. This is the number of operators plus
        the number of variables (counting multiplicity). Parentheses are ignored.

        EXAMPLES::

            sage: import sage.logic.propcalc as propcalc
            sage: s = propcalc.formula("a")
            sage: s.length()
            1
            sage: s = propcalc.formula("(a)")
            sage: s.length()
            1
            sage: s = propcalc.formula("~a")
            sage: s.length()
            2
            sage: s = propcalc.formula("a -> b")
            sage: s.length()
            3
            sage: s = propcalc.formula("alpha -> beta")
            sage: s.length()
            3
            sage: s = propcalc.formula("a -> a")
            sage: s.length()
            3
            sage: s = propcalc.formula("~(a -> b)")
            sage: s.length()
            4
            sage: s = propcalc.formula("((a&b)|(a&c))->~d")
            sage: s.length()
            10

        TESTS::

            sage: s = propcalc.formula("(((alpha) -> ((beta))))")
            sage: s.length()
            3
        """
        return len(flatten(self.full_tree()))

    # For backward compatibility, we allow `self.length()` to be called as
    # `len(self)`, but this may be deprecated in the future (see :issue:`32148`):
    __len__ = length


# allow is_consequence to be called as a function (not only as a method of BooleanFormula)
is_consequence = BooleanFormula.is_consequence
