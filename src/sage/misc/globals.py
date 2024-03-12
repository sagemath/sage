# sage_setup: distribution = sagemath-objects
r"""
Managing the main global namespace

TESTS:

The following test, verifying that :issue:`16181` has been resolved, needs
to stay at the beginning of this file so that its context is not
poisoned by other tests::

    sage: from sage.misc.globals import inject_variable
    sage: inject_variable('a', 0)
    sage: a
    0

Check the fix from :issue:`8323`::

    sage: 'name' in globals()
    False
    sage: 'func' in globals()
    False
"""

# ****************************************************************************
#       Copyright (C) 2010-2017 Nicolas M. Thiery
#                     2014      Darij Grinberg
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import warnings


def get_main_globals():
    """
    Return the main global namespace.

    EXAMPLES::

        sage: from sage.misc.globals import get_main_globals
        sage: G = get_main_globals()
        sage: bla = 1
        sage: G['bla']
        1
        sage: bla = 2
        sage: G['bla']
        2
        sage: G['ble'] = 5
        sage: ble
        5

    This is analogous to :func:`globals`, except that it can be called
    from any function, even if it is in a Python module::

        sage: def f():
        ....:     G = get_main_globals()
        ....:     assert G['bli'] == 14
        ....:     G['blo'] = 42
        sage: bli = 14
        sage: f()
        sage: blo
        42

    ALGORITHM:

    The main global namespace is discovered by going up the frame
    stack until the frame for the :mod:`__main__` module is found.
    Should this frame not be found (this should not occur in normal
    operation), an exception "ValueError: call stack is not deep
    enough" will be raised by ``_getframe``.

    See :meth:`inject_variable_test` for a real test that this works
    within deeply nested calls in a function defined in a Python
    module.
    """
    import sys
    depth = 0
    while True:
        G = sys._getframe(depth).f_globals
        if G.get("__name__", None) == "__main__":
            break
        depth += 1
    return G


def inject_variable(name, value, warn=True):
    """
    Inject a variable into the main global namespace.

    INPUT:

    - ``name``  -- a string
    - ``value`` -- anything
    - ``warn`` -- a boolean (default: :obj:`False`)

    EXAMPLES::

        sage: from sage.misc.globals import inject_variable
        sage: inject_variable("a", 314)
        sage: a
        314

    A warning is issued the first time an existing value is overwritten::

        sage: inject_variable("a", 271)
        doctest:...: RuntimeWarning: redefining global value `a`
        sage: a
        271
        sage: inject_variable("a", 272)
        sage: a
        272

    That's because warn seem to not reissue twice the same warning:

        sage: from warnings import warn
        sage: warn("blah")
        doctest:...: UserWarning: blah
        sage: warn("blah")

    Warnings can be disabled::

        sage: b = 3
        sage: inject_variable("b", 42, warn=False)
        sage: b
        42

    Use with care!
    """
    assert isinstance(name, str)
    # Using globals() does not work, even in Cython, because
    # inject_variable is called not only from the interpreter, but
    # also from functions in various modules.
    G = get_main_globals()
    if name in G and warn:
        warnings.warn("redefining global value `%s`" % name,
                      RuntimeWarning, stacklevel=2)
    G[name] = value


def inject_variable_test(name, value, depth):
    """
    A function for testing deep calls to inject_variable

    EXAMPLES::

        sage: from sage.misc.globals import inject_variable_test
        sage: inject_variable_test("a0", 314, 0)
        sage: a0
        314
        sage: inject_variable_test("a1", 314, 1)
        sage: a1
        314
        sage: inject_variable_test("a2", 314, 2)
        sage: a2
        314
        sage: inject_variable_test("a2", 271, 2)
        doctest:...: RuntimeWarning: redefining global value `a2`
        sage: a2
        271
    """
    if depth == 0:
        inject_variable(name, value)
    else:
        inject_variable_test(name, value, depth - 1)
