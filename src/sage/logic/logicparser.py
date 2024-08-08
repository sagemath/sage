r"""
Module that creates and modifies parse trees of well formed boolean formulas.

A parse tree of a boolean formula is a nested list, where each branch is either
a single variable, or a formula composed of either two variables and a binary
operator or one variable and a unary operator. The function parse produces
a parse tree that is simplified for the purposes of more efficient truth value
evaluation. The function :func:`~sage.logic.logicparser.polish_parse()`
produces the full parse tree of a boolean formula which is used in functions
related to proof and inference.  That is,
:func:`~sage.logic.logicparser.parse()` is meant to be used with functions
in the logic module that perform semantic operations on a boolean formula,
and :func:`~sage.logic.logicparser.polish_parse()` is to be used with
functions that perform syntactic operations on a boolean formula.

AUTHORS:

- Chris Gorecki (2007): initial version

- Paul Scurek (2013-08-01): added polish_parse, cleaned up python code,
  updated docstring formatting

- Paul Scurek (2013-08-06): added
  :func:`~sage.logic.logicparser.recover_formula()`,
  :func:`~sage.logic.logicparser.recover_formula_internal()`,
  :func:`~sage.logic.logicparser.prefix_to_infix()`,
  :func:`~sage.logic.logicparser.to_infix_internal()`

- Paul Scurek (2013-08-08): added get_trees, error handling in polish_parse,
  recover_formula_internal, and tree_parse

EXAMPLES:

Find the parse tree and variables of a string representation of a boolean
formula::

    sage: import sage.logic.logicparser as logicparser
    sage: s = 'a|b&c'
    sage: t = logicparser.parse(s)
    sage: t
    (['|', 'a', ['&', 'b', 'c']], ['a', 'b', 'c'])

Find the full syntax parse tree of a string representation of a boolean
formula::

    sage: import sage.logic.logicparser as logicparser
    sage: s = '(a&b)->~~c'
    sage: logicparser.polish_parse(s)
    ['->', ['&', 'a', 'b'], ['~', ['~', 'c']]]

Find the tokens and distinct variables of a boolean formula::

    sage: import sage.logic.logicparser as logicparser
    sage: s = '~(a|~b)<->(c->c)'
    sage: logicparser.tokenize(s)
    (['(', '~', '(', 'a', '|', '~', 'b', ')', '<->', '(', 'c', '->', 'c', ')', ')'], ['a', 'b', 'c'])

Find the parse tree of a boolean formula from a list of the formula's tokens::

    sage: import sage.logic.logicparser as logicparser
    sage: t = ['(', 'a', '->', '~', 'c', ')']
    sage: logicparser.tree_parse(t)
    ['->', 'a', ['~', 'c', None]]
    sage: r = ['(', '~', '~', 'a', '|', 'b', ')']
    sage: logicparser.tree_parse(r)
    ['|', 'a', 'b']

Find the full syntax parse tree of a boolean formula from a list of tokens::

    sage: import sage.logic.logicparser as logicparser
    sage: t = ['(', 'a', '->', '~', 'c', ')']
    sage: logicparser.tree_parse(t, polish = True)
    ['->', 'a', ['~', 'c']]
    sage: r = ['(', '~', '~', 'a', '|', 'b', ')']
    sage: logicparser.tree_parse(r, polish = True)
    ['|', ['~', ['~', 'a']], 'b']
"""
# ****************************************************************************
#       Copyright (C) 2007 Chris Gorecki <chris.k.gorecki@gmail.com>
#       Copyright (C) 2013 Paul Scurek <scurek86@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import string


__symbols = '()&|~<->^'
__op_list = ['~', '&', '|', '^', '->', '<->']


def parse(s):
    r"""
    Return a parse tree from a boolean formula ``s``.

    INPUT:

    - ``s`` -- string containing a boolean formula

    OUTPUT:

    A list containing the parse tree and a list containing the
    variables in a boolean formula in this order:

    1. the list containing the parse tree
    2. the list containing the variables

    EXAMPLES:

    This example illustrates how to produce the parse tree of a boolean
    formula ``s``::

        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|b&c'
        sage: t = logicparser.parse(s)
        sage: t
        (['|', 'a', ['&', 'b', 'c']], ['a', 'b', 'c'])
    """
    toks, vars_order = tokenize(s)
    tree = tree_parse(toks)
    # special case of tree == single variable
    if isinstance(tree, str):
        return ['&', tree, tree], vars_order
    return tree, vars_order


def polish_parse(s):
    r"""
    Return the full syntax parse tree from a boolean formula ``s``.

    INPUT:

    - ``s`` -- string containing a boolean expression

    OUTPUT: the full syntax parse tree as a nested list

    EXAMPLES:

    This example illustrates how to find the full syntax parse tree
    of a boolean formula::

        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|~~b'
        sage: t = logicparser.polish_parse(s)
        sage: t
        ['|', 'a', ['~', ['~', 'b']]]

    AUTHORS:

    - Paul Scurek (2013-08-03)
    """
    if (s.count('(') != s.count(')')) or not s:
        raise SyntaxError("malformed statement")

    toks, vars_order = tokenize(s)
    tree = tree_parse(toks, polish=True)
    # special case where the formula s is a single variable
    if isinstance(tree, str):
        return vars_order
    return tree


def get_trees(*statements):
    r"""
    Return the full syntax parse trees of the statements.

    INPUT:

    - ``*statements`` -- strings or :class:`BooleanFormula` instances

    OUTPUT: the parse trees in a list

    EXAMPLES:

    This example illustrates finding the parse trees of multiple formulas.

    ::

        sage: import sage.logic.propcalc as propcalc
        sage: import sage.logic.logicparser as logicparser
        sage: f = propcalc.formula("((a|b)&~~c)")
        sage: g = "a<->(~(c))"
        sage: h = "~b"
        sage: logicparser.get_trees(f, g, h)
        [['&', ['|', 'a', 'b'], ['~', ['~', 'c']]],
        ['<->', 'a', ['~', 'c']],
        ['~', 'b']]

    ::

        sage: i = "(~q->p)"
        sage: j = propcalc.formula("a")
        sage: logicparser.get_trees(i, j)
        [['->', ['~', 'q'], 'p'], ['a']]

    ::

        sage: k = "p"
        sage: logicparser.get_trees(k)
        [['p']]

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    trees = []
    from . import boolformula
    for statement in statements:
        if not isinstance(statement, boolformula.BooleanFormula):
            try:
                trees.append(polish_parse(statement))
            except (NameError, SyntaxError):
                raise SyntaxError("malformed statement")
        else:
            trees.append(statement.full_tree())
    return trees


def recover_formula(prefix_tree):
    r"""
    Recover the formula from a parse tree in prefix form.

    INPUT:

    - ``prefix_tree`` -- list; this is a full syntax parse
      tree in prefix form

    OUTPUT: the formula as a string

    EXAMPLES:

    This example illustrates the recovery of a formula from a parse tree::

        sage: import sage.logic.propcalc as propcalc
        sage: import sage.logic.logicparser as logicparser
        sage: t = ['->', ['&', 'a', ['~', ['~', 'c']]], ['~', ['|', ['~', 'c'], 'd']]]
        sage: logicparser.recover_formula(t)
        '(a&~~c)->~(~c|d)'

        sage: f = propcalc.formula("a&(~~c|d)")
        sage: logicparser.recover_formula(f.full_tree())
        'a&(~~c|d)'

        sage: r = ['~', 'a']
        sage: logicparser.recover_formula(r)
        '~a'

        sage: s = ['d']
        sage: logicparser.recover_formula(s)
        'd'

    .. NOTE::

        The function :func:`~sage.logic.logicparser.polish_parse()` may be
        passed as an argument, but :func:`~sage.logic.logicparser.tree_parse()`
        may not unless the parameter ``polish`` is set to ``True``.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    formula = ''
    if not isinstance(prefix_tree, list):
        raise TypeError("the input must be a parse tree as a list")

    formula = apply_func(prefix_tree, recover_formula_internal)
    if prefix_tree[0] == '~' or len(prefix_tree) == 1:
        return formula
    return formula[1:-1]


def recover_formula_internal(prefix_tree):
    r"""
    Recover the formula from a parse tree in prefix form.

    INPUT:

    - ``prefix_tree`` -- list; this is a simple tree
      with at most one operator in prefix form

    OUTPUT: the formula as a string

    EXAMPLES:

    This example illustrates recovering the formula from a parse tree::

        sage: import sage.logic.logicparser as logicparser
        sage: import sage.logic.propcalc as propcalc
        sage: t = ['->', 'a', 'b']
        sage: logicparser.recover_formula_internal(t)
        '(a->b)'

        sage: r = ['~', 'c']
        sage: logicparser.recover_formula_internal(r)
        '~c'

        sage: s = ['d']
        sage: logicparser.recover_formula_internal(s)
        'd'

    We can pass :func:`~sage.logic.logicparser.recover_formula_internal()`
    as an argument in :func:`~sage.logic.logicparser.apply_func()`::

        sage: f = propcalc.formula("~(d|c)<->(a&~~~c)")
        sage: logicparser.apply_func(f.full_tree(), logicparser.recover_formula_internal)
        '(~(d|c)<->(a&~~~c))'

    .. NOTE::

        This function is for internal use by :mod:`~sage.logic.logicparser`.
        The function recovers the formula of a simple parse tree in prefix
        form. A simple parse tree contains at most one operator.

        The function :func:`~sage.logic.logicparser.polish_parse()` may be
        passed as an argument, but :func:`~sage.logic.logicparser.tree_parse()`
        may not unless the parameter ``polish`` is set to ``True``.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    from .propcalc import formula as propcalc_formula
    if len(prefix_tree) == 3:
        bool_formula = '(' + prefix_tree[1] + prefix_tree[0] + prefix_tree[2] + ')'
    else:
        bool_formula = ''.join(prefix_tree)

    try:
        bool_formula = propcalc_formula(bool_formula)
    except (SyntaxError, NameError):
        raise SyntaxError

    return repr(bool_formula)


def prefix_to_infix(prefix_tree):
    r"""
    Convert a parse tree from prefix form to infix form.

    INPUT:

    - ``prefix_tree`` -- list; this is a full syntax parse
      tree in prefix form

    OUTPUT: list containing the tree in infix form

    EXAMPLES:

    This example illustrates converting a prefix tree to an infix tree::

        sage: import sage.logic.logicparser as logicparser
        sage: import sage.logic.propcalc as propcalc
        sage: t = ['|', ['~', 'a'], ['&', 'b', 'c']]
        sage: logicparser.prefix_to_infix(t)
        [['~', 'a'], '|', ['b', '&', 'c']]

    ::

        sage: f = propcalc.formula("(a&~b)<->~~~(c|d)")
        sage: logicparser.prefix_to_infix(f.full_tree())
        [['a', '&', ['~', 'b']], '<->', ['~', ['~', ['~', ['c', '|', 'd']]]]]

    .. NOTE::

        The function :func:`~sage.logic.logicparser.polish_parse()` may be
        passed as an argument, but :func:`~sage.logic.logicparser.tree_parse()`
        may not unless the parameter ``polish`` is set to ``True``.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    if not isinstance(prefix_tree, list):
        raise TypeError("the input must be a parse tree as a list")
    return apply_func(prefix_tree, to_infix_internal)


def to_infix_internal(prefix_tree):
    r"""
    Convert a simple parse tree from prefix form to infix form.

    INPUT:

    - ``prefix_tree`` -- list; this is a simple parse tree
      in prefix form with at most one operator

    OUTPUT: the tree in infix form as a list

    EXAMPLES:

    This example illustrates converting a simple tree from prefix
    to infix form::

        sage: import sage.logic.logicparser as logicparser
        sage: import sage.logic.propcalc as propcalc
        sage: t = ['|', 'a', 'b']
        sage: logicparser.to_infix_internal(t)
        ['a', '|', 'b']

    We can pass :func:`~sage.logic.logicparser.to_infix_internal()` as an
    argument in :func:`~sage.logic.logicparser.apply_func()`::

        sage: f = propcalc.formula("(a&~b)<->~~~(c|d)")
        sage: logicparser.apply_func(f.full_tree(), logicparser.to_infix_internal)
        [['a', '&', ['~', 'b']], '<->', ['~', ['~', ['~', ['c', '|', 'd']]]]]

    .. NOTE::

        This function is for internal use by :mod:`~sage.logic.logicparser`.
        It converts a simple parse tree from prefix form to infix form. A
        simple parse tree contains at most one operator.

        The function :func:`polish_parse` may be passed as an argument,
        but :func:`~sage.logic.logicparser.tree_parse()` may not unless the
        parameter ``polish`` is set to ``True``.

    AUTHORS:

    - Paul Scurek (2013-08-06)
    """
    if prefix_tree[0] != '~' and len(prefix_tree) == 3:
        return [prefix_tree[1], prefix_tree[0], prefix_tree[2]]
    return prefix_tree


def tokenize(s):
    r"""
    Return the tokens and the distinct variables appearing in a boolean
    formula ``s``.

    INPUT:

    - ``s`` -- string representation of a boolean formula

    OUTPUT:

    The tokens and variables as an ordered pair of lists in the following
    order:

    1. A list containing the tokens of ``s``, in the order they appear in ``s``
    2. A list containing the distinct variables in ``s``, in the order
       they appear in ``s``

    EXAMPLES:

    This example illustrates how to tokenize a string representation of a
    boolean formula::

        sage: import sage.logic.logicparser as logicparser
        sage: s = 'a|b&c'
        sage: t = logicparser.tokenize(s)
        sage: t
        (['(', 'a', '|', 'b', '&', 'c', ')'], ['a', 'b', 'c'])
    """
    i = 0
    toks = ['(']
    vars_order = []

    while i < len(s):
        tok = ""
        skip = valid = 1
        if s[i] in '()~&|^':
            tok = s[i]
        elif s[i:i + 2] == '->':
            tok = '->'
            skip = 2
        elif s[i:i + 3] == '<->':
            tok = '<->'
            skip = 3
        # check to see if '-', '<' or '>' are used incorrectly
        elif s[i] in '<->':
            raise SyntaxError("'{}' can only be used as part of the operators '<->' or '->'.".format(s[i]))

        if len(tok) > 0:
            toks.append(tok)
            i += skip
            continue

        # token is a variable name
        if s[i] == ' ':
            i += 1
            continue

        while i < len(s) and s[i] not in __symbols and s[i] != ' ':
            tok += s[i]
            i += 1

        if len(tok) > 0:
            if tok[0] not in string.ascii_letters:
                valid = 0
            for c in tok:
                if c not in string.ascii_letters and c not in string.digits and c != '_':
                    valid = 0

        if valid == 1:
            toks.append(tok)
            if tok not in vars_order:
                vars_order.append(tok)
        else:
            msg = "invalid variable name " + tok
            msg += ": identifiers must begin with a letter and contain only "
            msg += "alphanumerics and underscores"
            raise NameError(msg)

    toks.append(')')
    return toks, vars_order


def tree_parse(toks, polish=False):
    r"""
    Return a parse tree from the tokens in ``toks``.

    INPUT:

    - ``toks`` -- list of tokens from a boolean formula

    - ``polish`` -- boolean (default: ``False``); when ``True``,
      :func:`~sage.logic.logicparser.tree_parse()` will return
      the full syntax parse tree

    OUTPUT:

    A parse tree in the form of a nested list that depends on ``polish``
    as follows:

    - If ``False``, then return a simplified parse tree.

    - If ``True``, then return the full syntax parse tree.

    EXAMPLES:

    This example illustrates the use of
    :func:`~sage.logic.logicparser.tree_parse()` when ``polish`` is ``False``::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['(', 'a', '|', 'b', '&', 'c', ')']
        sage: logicparser.tree_parse(t)
        ['|', 'a', ['&', 'b', 'c']]

    We now demonstrate the use of :func:`~sage.logic.logicparser.tree_parse()`
    when ``polish`` is ``True``::

        sage: t = ['(', 'a', '->', '~', '~', 'b', ')']
        sage: logicparser.tree_parse(t)
        ['->', 'a', 'b']
        sage: t = ['(', 'a', '->', '~', '~', 'b', ')']
        sage: logicparser.tree_parse(t, polish = True)
        ['->', 'a', ['~', ['~', 'b']]]
    """
    if toks[1] in ['|', '&', '->', '<->', '^']:
        raise SyntaxError

    stack = []
    for tok in toks:
        stack.append(tok)
        if tok == ')':
            lrtoks = []
            while tok != '(':
                tok = stack.pop()
                lrtoks.insert(0, tok)
            branch = parse_ltor(lrtoks[1:-1], polish=polish)
            stack.append(branch)
    return stack[0]


def parse_ltor(toks, n=0, polish=False):
    r"""
    Return a parse tree from ``toks``, where each token in ``toks`` is atomic.

    INPUT:

    - ``toks`` -- list of tokens; each token is atomic

    - ``n`` -- integer (default: 0) representing which order of
      operations are occurring

    - ``polish`` -- boolean (default: ``False``); when ``True``, double
      negations are not cancelled and negated statements are turned into
      list of length two

    OUTPUT: the parse tree as a nested list that depends on ``polish`` as follows:

    - If ``False``, then return a simplified parse tree.

    - If ``True``, then return the full syntax parse tree.

    EXAMPLES:

    This example illustrates the use of
    :func:`~sage.logic.logicparser.parse_ltor()` when ``polish`` is ``False``::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['a', '|', 'b', '&', 'c']
        sage: logicparser.parse_ltor(t)
        ['|', 'a', ['&', 'b', 'c']]

    ::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['a', '->', '~', '~', 'b']
        sage: logicparser.parse_ltor(t)
        ['->', 'a', 'b']

    We now repeat the previous example, but with ``polish`` being ``True``::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['a', '->', '~', '~', 'b']
        sage: logicparser.parse_ltor(t, polish = True)
        ['->', 'a', ['~', ['~', 'b']]]
    """
    i = 0
    for tok in toks:
        if tok == __op_list[n]:
            if tok == '~':
                if not polish:
                    # cancel double negations
                    if toks[i] == '~' and toks[i + 1] == '~':
                        del toks[i]
                        del toks[i]
                        return parse_ltor(toks, n)
                    args = [toks[i], toks[i + 1], None]
                    toks[i] = args
                    del toks[i + 1]
                    return parse_ltor(toks, n)
                # This executes when creating the full syntax parse tree
                else:
                    j = i
                    while toks[j] == '~':
                        j += 1
                    while j > i:
                        args = [toks[j - 1], toks[j]]
                        toks[j - 1] = args
                        del toks[j]
                        j -= 1
                    return parse_ltor(toks, n=n, polish=polish)
            else:
                args = [toks[i - 1], toks[i], toks[i + 1]]
                toks[i - 1] = [args[1], args[0], args[2]]
                del toks[i]
                del toks[i]
                return parse_ltor(toks, n)
        i += 1
    if n + 1 < len(__op_list):
        return parse_ltor(toks, n + 1)
    if len(toks) > 1:
        raise SyntaxError
    return toks[0]


def apply_func(tree, func):
    r"""
    Apply ``func`` to each node of ``tree``, and return a new parse tree.

    INPUT:

    - ``tree`` -- a parse tree of a boolean formula

    - ``func`` -- a function to be applied to each node of tree; this may
      be a function that comes from elsewhere in the logic module

    OUTPUT: the new parse tree in the form of a nested list

    EXAMPLES:

    This example uses :func:`~sage.logic.logicparser.apply_func()` where
    ``func`` switches two entries of tree::

        sage: import sage.logic.logicparser as logicparser
        sage: t = ['|', ['&', 'a', 'b'], ['&', 'a', 'c']]
        sage: f = lambda t: [t[0], t[2], t[1]]
        sage: logicparser.apply_func(t, f)
        ['|', ['&', 'c', 'a'], ['&', 'b', 'a']]
    """
    # used when full syntax parse tree is passed as argument
    if len(tree) == 1:
        return func(tree)
    # used when full syntax parse tree is passed as argument
    elif len(tree) == 2:
        rval = apply_func(tree[1], func)
        return func([tree[0], rval])
    elif isinstance(tree[1], list) and isinstance(tree[2], list):
        lval = apply_func(tree[1], func)
        rval = apply_func(tree[2], func)
    elif isinstance(tree[1], list):
        lval = apply_func(tree[1], func)
        rval = tree[2]
    elif isinstance(tree[2], list):
        lval = tree[1]
        rval = apply_func(tree[2], func)
    else:
        lval = tree[1]
        rval = tree[2]
    return func([tree[0], lval, rval])
