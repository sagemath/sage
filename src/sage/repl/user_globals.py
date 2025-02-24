# sage_setup: distribution = sagemath-repl
# sage.doctest: needs sage.modules
r"""
User-interface globals

These represent "globals" in the top-level user interface (command
line, notebook). When starting up Sage, ``user_globals`` becomes
identical to the ``globals()`` of the top-level frame. Later on,
functions can inject globals.

.. NOTE::

    Injecting globals is a feature which has some legitimate uses, but
    beware of over-using it. A natural use case for injecting globals is
    for "compile-and-import" functions, such as :func:`cython`.

EXAMPLES:

This is how a typical user interface initializes the globals::

    sage: ui_globals = globals()  # or wherever the user interface stores its globals
    sage: from sage import all_cmdline
    sage: from sage.repl.user_globals import initialize_globals
    sage: _ = initialize_globals(all_cmdline, ui_globals)

Now everything which was imported in ``all_cmdline`` is available as a
global::

    sage: from sage.repl.user_globals import get_global, set_global
    sage: get_global("Matrix")
    <cyfunction matrix at ...>

This is exactly the same::

    sage: ui_globals["Matrix"]
    <cyfunction matrix at ...>

We inject a global::

    sage: set_global("myvar", "Hello World!")
    sage: get_global("myvar")
    'Hello World!'

If we set up things correctly, this new variable is now actually
available as global::

    sage: myvar
    'Hello World!'

AUTHORS:

- Jeroen Demeyer (2015-03-30): initial version (:issue:`12446`)
"""

# ****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


user_globals = None


def _check():
    """
    Raise :exc:`RuntimeError` if ``user_globals`` has not been initialized.

    EXAMPLES::

        sage: import sage.repl.user_globals
        sage: sage.repl.user_globals._check()
        sage: sage.repl.user_globals.user_globals = None
        sage: sage.repl.user_globals._check()
        Traceback (most recent call last):
        ...
        RuntimeError: the user-space globals dictionary has not been initialized...
    """
    if user_globals is None:
        raise RuntimeError(
            "the user-space globals dictionary has not been initialized. "
            "Use initialize_globals() or set_globals() or use a different "
            "function which doesn't need these globals")


def get_globals():
    """
    Return the dictionary of all user globals.

    EXAMPLES::

        sage: from sage.repl.user_globals import get_globals, initialize_globals
        sage: initialize_globals(sage.all)
        sage: get_globals()["Matrix"]
        <cyfunction matrix at ...>
    """
    _check()
    return user_globals


def set_globals(g):
    """
    Set the dictionary of all user globals to ``g``.

    INPUT:

    - ``g`` -- dictionary; typically, this will be some dictionary
      given by the user interface or just ``globals()``

    EXAMPLES::

        sage: from sage.repl.user_globals import get_globals, set_globals
        sage: my_dict = dict()
        sage: set_globals(my_dict)
        sage: my_dict is get_globals()
        True
    """
    global user_globals
    user_globals = g


def initialize_globals(all, g=None):
    """
    Set the user globals dictionary to ``g`` and assign everything
    which was imported in module ``all`` as global.

    INPUT:

    - ``all`` -- a module whose globals will be injected

    - ``g`` -- dictionary; see :func:`set_globals`. If this is
      ``None``, keep the current globals dictionary.

    EXAMPLES::

        sage: my_globs = {"foo": "bar"}
        sage: from sage.repl.user_globals import initialize_globals
        sage: initialize_globals(sage.all, my_globs)
        sage: my_globs["foo"]
        'bar'
        sage: my_globs["Matrix"]
        <cyfunction matrix at ...>

    Remove ``Matrix`` from the globals and initialize again without
    changing the dictionary::

        sage: del my_globs["Matrix"]
        sage: initialize_globals(sage.all)
        sage: my_globs["Matrix"]
        <cyfunction matrix at ...>
    """
    if g is not None:
        set_globals(g)
    for key in dir(all):
        if key[0] != '_':
            user_globals[key] = getattr(all, key)
    from sage.misc.lazy_import import clean_namespace
    clean_namespace(user_globals)


def get_global(name):
    """
    Return the value of global variable ``name``.

    Raise :exc:`NameError` if there is no such global variable.

    INPUT:

    - ``name`` -- string representing a variable name

    OUTPUT: the value of variable ``name``

    EXAMPLES::

        sage: from sage.repl.user_globals import get_global
        sage: the_answer = 42
        sage: get_global("the_answer")
        42
        sage: get_global("the_question")
        Traceback (most recent call last):
        ...
        NameError: name 'the_question' is not defined
    """
    _check()
    try:
        return user_globals[name]
    except KeyError:
        raise NameError("name {!r} is not defined".format(name))


def set_global(name, value):
    """
    Assign ``value`` to global variable ``name``. This is equivalent
    to executing ``name = value`` in the global namespace.

    INPUT:

    - ``name`` -- string representing a variable name

    - ``value`` -- a value to assign to the variable

    EXAMPLES::

        sage: from sage.repl.user_globals import set_global
        sage: set_global("the_answer", 42)
        sage: the_answer
        42
    """
    _check()
    user_globals[name] = value
