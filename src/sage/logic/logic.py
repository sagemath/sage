r"""
Symbolic Logic Expressions

An expression is created from a string that consists of the
operators ``!``, ``&``, ``|``, ``->``, ``<->``, which correspond to the
logical functions not, and, or, if then, if and only if, respectively.
Variable names must start with a letter and contain only
alpha-numerics and the underscore character.

AUTHORS:

- Chris Gorecki (2007): initial version

- William Stein (2007-08-31): integration into Sage 2.8.4

- Paul Scurek (2013-08-03): updated docstring formatting
"""
# ****************************************************************************
#       Copyright (C) 2007 Chris Gorecki <chris.k.gorecki@gmail.com>
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Paul Scurek <scurek86@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import string

# constants
tok_list = ['OPAREN', 'CPAREN', 'AND', 'OR', 'NOT', 'IFTHEN', 'IFF']
bin_list = ['AND', 'OR', 'IFTHEN', 'IFF']
operators = '()&|!<->'
# variables
vars = {}
vars_order = []


class SymbolicLogic:
    """
    EXAMPLES:

    This example illustrates how to create a boolean formula and print
    its table::

        sage: log = SymbolicLogic()
        sage: s = log.statement("a&b|!(c|a)")
        sage: t = log.truthtable(s)
        sage: log.print_table(t)
        a     | b     | c     | value |
        --------------------------------
        False | False | False | True  |
        False | False | True  | False |
        False | True  | False | True  |
        False | True  | True  | False |
        True  | False | False | False |
        True  | False | True  | False |
        True  | True  | False | True  |
        True  | True  | True  | True  |
    """
    def statement(self, s):
        r"""
        Return a token list to be used by other functions in the class.

        INPUT:

        - ``s`` -- string containing the logic expression to be manipulated

        - ``global vars`` -- dictionary with variable names as keys and the
          variables' current boolean values as dictionary values

        - ``global vars_order`` -- list of the variables in the order
          that they are found

        OUTPUT: list of length three containing the following in this order:

        1. a list of tokens
        2. a dictionary of variable/value pairs
        3. a list of the variables in the order they were found

        EXAMPLES:

        This example illustrates the creation of a statement::

            sage: log = SymbolicLogic()
            sage: s = log.statement("a&b|!(c|a)")
            sage: s2 = log.statement("!((!(a&b)))")

        It is an error to use invalid variable names::

            sage: s = log.statement("3fe & @q")
            Invalid variable name:  3fe
            Invalid variable name:  @q

        It is also an error to use invalid syntax::

            sage: s = log.statement("a&&b")
            Malformed Statement
            sage: s = log.statement("a&((b)")
            Malformed Statement
        """
        global vars, vars_order
        toks, vars, vars_order = ['OPAREN'], {}, []
        tokenize(s, toks)
        statement = [toks, vars, vars_order]
        try:                # verify the syntax
            eval(toks)
        except (KeyError, RuntimeError):
            print('Malformed Statement')
            return []
        return statement

    def truthtable(self, statement, start=0, end=-1):
        r"""
        Return a truth table.

        INPUT:

        - ``statement`` -- list; it contains the tokens and the two global
          variables vars and vars_order

        - ``start`` -- integer (default: 0); this represents the row of
          the truth table from which to start

        - ``end`` -- integer (default: -1); this represents the last row
          of the truth table to be created

        OUTPUT:

        The truth table as a 2d array with the creating formula tacked
        to the front.

        EXAMPLES:

        This example illustrates the creation of a statement::

            sage: log = SymbolicLogic()
            sage: s = log.statement("a&b|!(c|a)")
            sage: t = log.truthtable(s) #creates the whole truth table

        We can now create truthtable of rows 1 to 5::

            sage: s2 = log.truthtable(s, 1, 5); s2
            [[['OPAREN',
               'a',
               'AND',
               'b',
               'OR',
               'NOT',
               'OPAREN',
               'c',
               'OR',
               'a',
               'CPAREN',
               'CPAREN'],
              {'a': 'True', 'b': 'False', 'c': 'False'},
              ['a', 'b', 'c']],
             ['False', 'False', 'True', 'False'],
             ['False', 'True', 'False', 'True'],
             ['False', 'True', 'True', 'False'],
             ['True', 'False', 'False', 'False']]

        .. NOTE::

            When sent with no start or end parameters this is an
            exponential time function requiring `O(2^n)` time, where
            `n` is the number of variables in the logic expression

        TESTS:

        Verify that :issue:`32676` is fixed::

            sage: s = log.statement("a&b|!(c|a)")
            sage: copy_s2 = copy(s[2])
            sage: t = log.truthtable(s)
            sage: s[2] == copy_s2
            True
        """
        global vars, vars_order
        toks, vars, vars_order = statement
        if end == -1:
            end = 2 ** len(vars)
        table = [statement]
        keys = vars_order
        for i in range(start, end):
            j = 0
            row = []
            for key in reversed(keys):
                bit = get_bit(i, j)
                vars[key] = bit
                j += 1
                row.insert(0, bit)
            row.append(eval(toks))
            table.append(row)
        return table

    def print_table(self, table):
        r"""
        Return a truthtable corresponding to the given statement.

        INPUT:

        - ``table`` -- object created by :meth:`truthtable()` method; it
          contains the variable values and the evaluation of the statement

        OUTPUT: a formatted version of the truth table

        EXAMPLES:

        This example illustrates the creation of a statement and
        its truth table::

            sage: log = SymbolicLogic()
            sage: s = log.statement("a&b|!(c|a)")
            sage: t = log.truthtable(s) #creates the whole truth table
            sage: log.print_table(t)
            a     | b     | c     | value |
            --------------------------------
            False | False | False | True  |
            False | False | True  | False |
            False | True  | False | True  |
            False | True  | True  | False |
            True  | False | False | False |
            True  | False | True  | False |
            True  | True  | False | True  |
            True  | True  | True  | True  |

        We can also print a shortened table::

            sage: t = log.truthtable(s, 1, 5)
            sage: log.print_table(t)
            a     | b     | c     | value |
            --------------------------------
            False | False | True  | False |
            False | True  | False | True  |
            False | True  | True  | False |
            True  | False | False | False |

        TESTS:

        Verify that :issue:`32676` is fixed::

            sage: table = log.truthtable(log.statement("A->B"))
            sage: table_copy = table.copy()
            sage: log.print_table(table)
            ...
            sage: table_copy == table
            True
        """
        statement = table[0]
        table = table[1:]
        vars_order = statement[2].copy()
        vars_len = []
        line = s = ""
        vars_order.append('value')
        for var in vars_order:
            vars_len.append(len(var))
            s = var + ' '
            while len(s) < len('False '):
                s += ' '
            s += '| '
            line += s
        print(line)
        print(len(line) * '-')
        for row in table:
            line = s = ""
            i = 0
            for e in row:
                j = 2 if e == 'True' else 1
                s = e + ' ' * j
                if i < len(vars_len):
                    while len(s) <= vars_len[i]:
                        s += ' '
                s += '| '
                line += s
                i += 1
            print(line)
        print("")

    def combine(self, statement1, statement2):
        r"""
        Return a new statement which contains the
        two statements or'd together.

        INPUT:

        - ``statement1`` -- the first statement
        - ``statement2`` -- the second statement

        OUTPUT: a new statement which or'd the given statements together

        EXAMPLES::

            sage: log = SymbolicLogic()
            sage: s1 = log.statement("(a&b)")
            sage: s2 = log.statement("b")
            sage: log.combine(s1,s2)
            [['OPAREN',
              'OPAREN',
              'OPAREN',
              'a',
              'AND',
              'b',
              'CPAREN',
              'CPAREN',
              'OR',
              'OPAREN',
              'b',
              'CPAREN',
              'CPAREN'],
             {'a': 'False', 'b': 'False'},
             ['a', 'b', 'b']]
        """
        toks = ['OPAREN'] + statement1[0] + ['OR'] + statement2[0] + ['CPAREN']
        variables = dict(statement1[1])
        variables.update(statement2[1])
        var_order = statement1[2] + statement2[2]
        return [toks, variables, var_order]

    # TODO: implement the simplify function which calls
    # a c++ implementation of the ESPRESSO algorithm
    # to simplify the truthtable: probably Minilog
    def simplify(self, table):
        """
        Call a C++ implementation of the ESPRESSO algorithm to simplify the
        given truth table.

        .. TODO::

            Implement this method.

        EXAMPLES::

            sage: log = SymbolicLogic()
            sage: s = log.statement("a&b|!(c|a)")
            sage: t = log.truthtable(s)
            sage: log.simplify(t)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def prove(self, statement):
        """
        A function to test to see if the statement is a tautology or
        contradiction by calling a C++ library.

        .. TODO::

            Implement this method.

        EXAMPLES::

            sage: log = SymbolicLogic()
            sage: s = log.statement("a&b|!(c|a)")
            sage: log.prove(s)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


def get_bit(x, c):
    r"""
    Determine if bit ``c`` of the number ``x`` is 1.

    INPUT:

    - ``x`` -- integer; this is the number from which to take the bit

    - ``c`` -- integer; this is the bit number to be taken

    OUTPUT: a boolean value to be determined as follows:

    - ``True`` if bit ``c`` of ``x`` is 1.

    - ``False`` if bit ``c`` of ``x`` is not 1.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    EXAMPLES::

        sage: from sage.logic.logic import get_bit
        sage: get_bit(int(2), int(1))
        'True'
        sage: get_bit(int(8), int(0))
        'False'
    """
    bits = []
    while x > 0:
        b = 'False' if x % 2 == 0 else 'True'
        x = x // 2
        bits.append(b)
    if c > len(bits) - 1:
        return 'False'
    else:
        return bits[c]


def eval(toks):
    r"""
    Evaluate the expression contained in ``toks``.

    INPUT:

    - ``toks`` -- list of tokens; this represents a boolean expression

    OUTPUT: a boolean value to be determined as follows:

    - ``True`` if expression evaluates to ``True``.

    - ``False`` if expression evaluates to ``False``.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.
        The evaluations rely on setting the values of the variables in the
        global dictionary vars.

    TESTS::

        sage: log = SymbolicLogic()
        sage: s = log.statement("a&b|!(c|a)")
        sage: sage.logic.logic.eval(s[0])
        'True'
    """
    stack = []
    for tok in toks:
        stack.append(tok)
        if tok == 'CPAREN':
            lrtoks = []
            while tok != 'OPAREN':
                tok = stack.pop()
                lrtoks.insert(0, tok)
            stack.append(eval_ltor_toks(lrtoks[1:-1]))
    if len(stack) > 1:
        raise RuntimeError
    return stack[0]


def eval_ltor_toks(lrtoks):
    r"""
    Evaluates the expression contained in ``lrtoks``.

    INPUT:

    - ``lrtoks`` -- list of tokens; this represents a part of a boolean
      formula that contains no inner parentheses

    OUTPUT: a boolean value to be determined as follows:

    - ``True`` if expression evaluates to ``True``.

    - ``False`` if expression evaluates to ``False``.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.
        The evaluations rely on setting the values of the variables in the
        global dictionary vars.

    TESTS::

        sage: log = SymbolicLogic()
        sage: s = log.statement("a&b|!c")
        sage: ltor = s[0][1:-1]; ltor
        ['a', 'AND', 'b', 'OR', 'NOT', 'c']
        sage: sage.logic.logic.eval_ltor_toks(ltor)
        'True'
    """
    reduce_monos(lrtoks)        # monotonic ! operators go first
    reduce_bins(lrtoks)         # then the binary operators
    if len(lrtoks) > 1:
        raise RuntimeError
    return lrtoks[0]


def reduce_bins(lrtoks):
    r"""
    Evaluate ``lrtoks`` to a single boolean value.

    INPUT:

    - ``lrtoks`` -- list of tokens; this represents a part of a boolean
      formula that contains no inner parentheses or monotonic operators

    OUTPUT:

    ``None``; the pointer to lrtoks is now a list containing
    ``True`` or ``False``.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: log = SymbolicLogic()
        sage: s = log.statement("a&b|c")
        sage: lrtoks = s[0][1:-1]; lrtoks
        ['a', 'AND', 'b', 'OR', 'c']
        sage: sage.logic.logic.reduce_bins(lrtoks); lrtoks
        ['False']
    """
    i = 0
    while i < len(lrtoks):
        if lrtoks[i] in bin_list:
            args = [lrtoks[i - 1], lrtoks[i], lrtoks[i + 1]]
            lrtoks[i - 1] = eval_bin_op(args)
            del lrtoks[i]
            del lrtoks[i]
            reduce_bins(lrtoks)
        i += 1


def reduce_monos(lrtoks):
    r"""
    Replace monotonic operator/variable pairs with a boolean value.

    INPUT:

    - ``lrtoks`` -- list of tokens; this represents a part of a boolean
      expression that contains now inner parentheses

    OUTPUT:

    ``None``; the pointer to ``lrtoks`` is now a list containing
    monotonic operators.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: log = SymbolicLogic()
        sage: s = log.statement("!a&!b")
        sage: lrtoks = s[0][1:-1]; lrtoks
        ['NOT', 'a', 'AND', 'NOT', 'b']
        sage: sage.logic.logic.reduce_monos(lrtoks); lrtoks
        ['True', 'AND', 'True']
    """
    i = 0
    while i < len(lrtoks):
        if lrtoks[i] == 'NOT':
            args = [lrtoks[i], lrtoks[i + 1]]
            lrtoks[i] = eval_mon_op(args)
            del lrtoks[i + 1]
        i += 1


def eval_mon_op(args):
    r"""
    Return a boolean value based on the truth table of the operator
    in ``args``.

    INPUT:

    - ``args`` -- list of length 2; this contains the token 'NOT' and
      then a variable name

    OUTPUT: a boolean value to be determined as follows:

    - ``True`` if the variable in ``args`` is ``False``.

    - ``False`` if the variable in ``args`` is ``True``.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: log = SymbolicLogic()
        sage: s = log.statement("!(a&b)|!a"); s
        [['OPAREN', 'NOT', 'OPAREN', 'a', 'AND', 'b', 'CPAREN', 'OR',
         'NOT', 'a', 'CPAREN'],
         {'a': 'False', 'b': 'False'},
         ['a', 'b']]
        sage: sage.logic.logic.eval_mon_op(['NOT', 'a'])
        'True'
    """
    val = vars[args[1]] if args[1] != 'True' and args[1] != 'False' else args[1]

    return 'False' if val == 'True' else 'True'


def eval_bin_op(args):
    r"""
    Return a boolean value based on the truth table of the operator
    in ``args``.

    INPUT:

    - ``args`` -- list of length 3; this contains a variable name,
      then a binary operator, and then a variable name, in that order

    OUTPUT:

    A boolean value; this is the evaluation of the operator based on the
    truth values of the variables.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: log = SymbolicLogic()
        sage: s = log.statement("!(a&b)"); s
        [['OPAREN', 'NOT', 'OPAREN', 'a', 'AND', 'b', 'CPAREN', 'CPAREN'],
         {'a': 'False', 'b': 'False'},
         ['a', 'b']]
        sage: sage.logic.logic.eval_bin_op(['a', 'AND', 'b'])
        'False'
    """
    if args[0] == 'False':
        lval = 'False'
    elif args[0] == 'True':
        lval = 'True'
    else:
        lval = vars[args[0]]

    if args[2] == 'False':
        rval = 'False'
    elif args[2] == 'True':
        rval = 'True'
    else:
        rval = vars[args[2]]

    if args[1] == 'AND':
        return eval_and_op(lval, rval)
    elif args[1] == 'OR':
        return eval_or_op(lval, rval)
    elif args[1] == 'IFTHEN':
        return eval_ifthen_op(lval, rval)
    elif args[1] == 'IFF':
        return eval_iff_op(lval, rval)


def eval_and_op(lval, rval):
    r"""
    Apply the 'and' operator to ``lval`` and ``rval``.

    INPUT:

    - ``lval`` -- string; this represents the value of the variable
      appearing to the left of the 'and' operator

    - ``rval`` -- string; this represents the value of the variable
      appearing to the right of the 'and' operator

    OUTPUT: the result of applying 'and' to ``lval`` and ``rval`` as a string

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: sage.logic.logic.eval_and_op('False', 'False')
        'False'
        sage: sage.logic.logic.eval_and_op('False', 'True')
        'False'
        sage: sage.logic.logic.eval_and_op('True', 'False')
        'False'
        sage: sage.logic.logic.eval_and_op('True', 'True')
        'True'
    """
    return 'True' if (lval == 'True' == rval) else 'False'


def eval_or_op(lval, rval):
    r"""
    Apply the 'or' operator to ``lval`` and ``rval``.

    INPUT:

    - ``lval`` -- string; this represents the value of the variable
      appearing to the left of the 'or' operator

    - ``rval`` -- string; this represents the value of the variable
      appearing to the right of the 'or' operator

    OUTPUT: string representing the result of applying 'or' to ``lval`` and ``rval``

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: sage.logic.logic.eval_or_op('False', 'False')
        'False'
        sage: sage.logic.logic.eval_or_op('False', 'True')
        'True'
        sage: sage.logic.logic.eval_or_op('True', 'False')
        'True'
        sage: sage.logic.logic.eval_or_op('True', 'True')
        'True'
    """
    return 'True' if (lval == 'True' or rval == 'True') else 'False'


def eval_ifthen_op(lval, rval):
    r"""
    Apply the 'if then' operator to ``lval`` and ``rval``.

    INPUT:

    - ``lval`` -- string; this represents the value of the variable
      appearing to the left of the 'if then' operator

    - ``rval`` -- string; this represents the value of the variable
      appearing to the right of the 'if then' operator

    OUTPUT:

    A string representing the result of applying 'if then' to
    ``lval`` and ``rval``.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: sage.logic.logic.eval_ifthen_op('False', 'False')
        'True'
        sage: sage.logic.logic.eval_ifthen_op('False', 'True')
        'True'
        sage: sage.logic.logic.eval_ifthen_op('True', 'False')
        'False'
        sage: sage.logic.logic.eval_ifthen_op('True', 'True')
        'True'
    """
    return 'False' if (lval == 'True' and rval == 'False') else 'True'


def eval_iff_op(lval, rval):
    r"""
    Apply the 'if and only if' operator to ``lval`` and ``rval``.

    INPUT:

    - ``lval`` -- string; this represents the value of the variable
      appearing to the left of the 'if and only if' operator

    - ``rval`` -- string; this represents the value of the variable
      appearing to the right of the 'if and only if' operator

    OUTPUT:

    A string representing the result of applying 'if and only if'
    to ``lval`` and ``rval``.

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    TESTS::

        sage: sage.logic.logic.eval_iff_op('False', 'False')
        'True'
        sage: sage.logic.logic.eval_iff_op('False', 'True')
        'False'
        sage: sage.logic.logic.eval_iff_op('True', 'False')
        'False'
        sage: sage.logic.logic.eval_iff_op('True', 'True')
        'True'
    """
    return 'True' if (lval == rval) else 'False'


def tokenize(s, toks):
    r"""
    Tokenize ``s`` and place the tokens of ``s`` in ``toks``.

    INPUT:

    - ``s`` -- string; this contains a boolean expression

    - ``toks`` -- list; this will be populated with the tokens of ``s``

    OUTPUT: none; the tokens of ``s`` are placed in ``toks``

    .. NOTE::

        This function is for internal use by the :class:`SymbolicLogic` class.

    EXAMPLES::

        sage: from sage.logic.logic import tokenize
        sage: toks = []
        sage: tokenize("(a&b)|c", toks)
        sage: toks
        ['OPAREN', 'a', 'AND', 'b', 'CPAREN', 'OR', 'c', 'CPAREN']
    """
    i = 0
    while i < len(s):
        tok = ""
        skip = valid = 1
        if s[i] == '(':
            tok = tok_list[0]
        elif s[i] == ')':
            tok = tok_list[1]
        elif s[i] == '&':
            tok = tok_list[2]
        elif s[i] == '|':
            tok = tok_list[3]
        elif s[i] == '!':
            tok = tok_list[4]
        elif s[i:i + 2] == '->':
            tok = tok_list[5]
            skip = 2
        elif s[i:i + 3] == '<->':
            tok = tok_list[6]
            skip = 3

        if tok:
            toks.append(tok)
            i += skip
            continue
        else:
            # token is a variable name
            if s[i] == ' ':
                i += 1
                continue

            while i < len(s) and s[i] not in operators and s[i] != ' ':
                tok += s[i]
                i += 1

            if tok:
                if tok[0] not in string.ascii_letters:
                    valid = 0
                for c in tok:
                    if not (c in string.ascii_letters
                            or c in string.digits or c == '_'):
                        valid = 0

            if valid == 1:
                toks.append(tok)
                vars[tok] = 'False'
                if tok not in vars_order:
                    vars_order.append(tok)
            else:
                print('Invalid variable name: ', tok)
                toks = []

    toks.append('CPAREN')
