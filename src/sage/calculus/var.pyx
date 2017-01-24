# cython: old_style_globals=True
"""
Symbolic variables
"""

from sage.symbolic.function_factory import function as new_function
from sage.symbolic.ring import SR


def var(*args, **kwds):
    r"""
    Create a symbolic variable with the name *s*.

    INPUT:

    - ``args`` -- A single string ``var('x y')``, a list of strings
      ``var(['x','y'])``, or multiple strings ``var('x', 'y')``. A
      single string can be either a single variable name, or a space
      or comma separated list of variable names. In a list or tuple of
      strings, each entry is one variable. If multiple arguments are
      specified, each argument is taken to be one variable. Spaces
      before or after variable names are ignored.

    - ``kwds`` -- keyword arguments can be given to specify domain and
      custom latex_name for variables. See EXAMPLES for usage.

    .. NOTE::

       The new variable is both returned and automatically injected
       into the global namespace. If you need a symbolic variable in
       library code, you must use either ``SR.var()``
       or ``SR.symbol()``.

    OUTPUT:

    If a single symbolic variable was created, the variable
    itself. Otherwise, a tuple of symbolic variables. The variable
    names are checked to be valid Python identifiers and a
    ``ValueError`` is raised otherwise.

    EXAMPLES:

    Here are the different ways to define three variables ``x``, ``y``,
    and ``z`` in a single line::

        sage: var('x y z')
        (x, y, z)
        sage: var('x, y, z')
        (x, y, z)
        sage: var(['x', 'y', 'z'])
        (x, y, z)
        sage: var('x', 'y', 'z')
        (x, y, z)
        sage: var('x'), var('y'), var(z)
        (x, y, z)

    We define some symbolic variables::

        sage: var('n xx yy zz')
        (n, xx, yy, zz)

    Then we make an algebraic expression out of them::

        sage: f = xx^n + yy^n + zz^n; f
        xx^n + yy^n + zz^n

    By default, var returns a complex variable. To define real or positive
    variables we can specify the domain as::

        sage: x = var('x', domain=RR); x; x.conjugate()
        x
        x
        sage: y = var('y', domain='real'); y.conjugate()
        y
        sage: y = var('y', domain='positive'); y.abs()
        y

    Custom latex expression can be assigned to variable::

        sage: x = var('sui', latex_name="s_{u,i}"); x._latex_()
        '{s_{u,i}}'

    In notebook, we can also colorize latex expression::

        sage: x = var('sui', latex_name="\\color{red}{s_{u,i}}"); x._latex_()
        '{\\color{red}{s_{u,i}}}'

    We can substitute a new variable name for n::

        sage: f(n = var('sigma'))
        xx^sigma + yy^sigma + zz^sigma

    If you make an important built-in variable into a symbolic variable,
    you can get back the original value using restore::

        sage: var('QQ RR')
        (QQ, RR)
        sage: QQ
        QQ
        sage: restore('QQ')
        sage: QQ
        Rational Field

    We make two new variables separated by commas::

        sage: var('theta, gamma')
        (theta, gamma)
        sage: theta^2 + gamma^3
        gamma^3 + theta^2

    The new variables are of type Expression, and belong
    to the symbolic expression ring::

        sage: type(theta)
        <type 'sage.symbolic.expression.Expression'>
        sage: parent(theta)
        Symbolic Ring

    TESTS::

        sage: var('q',ns=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: The new (Pynac) symbolics are now the only symbolics; please do not use keyword `ns` any longer.
        sage: q
        Traceback (most recent call last):
        ...
        NameError: name 'q' is not defined
        sage: var('q',ns=1)
        doctest:...: DeprecationWarning: The new (Pynac) symbolics are now the only symbolics; please do not use keyword 'ns' any longer.
        See http://trac.sagemath.org/6559 for details.
        q
    """
    if len(args)==1:
        name = args[0]
    else:
        name = args
    G = globals()  # this is the reason the code must be in Cython.
    if 'ns' in kwds:
        if kwds['ns']:
            from sage.misc.superseded import deprecation
            deprecation(6559, "The new (Pynac) symbolics are now the only symbolics; please do not use keyword 'ns' any longer.")
        else:
            raise NotImplementedError("The new (Pynac) symbolics are now the only symbolics; please do not use keyword `ns` any longer.")
        kwds.pop('ns')
    v = SR.var(name, **kwds)
    if isinstance(v, tuple):
        for x in v:
            G[repr(x)] = x
    else:
        G[repr(v)] = v
    return v

def function(s, *args, **kwds):
    r"""
    Create a formal symbolic function with the name *s*.

    INPUT:

    - ``s`` - a string, either a single variable name, or a space or
      comma separated list of variable names.

    - ``**kwds`` - keyword arguments. Either one of the following two
        keywords can be used to customize latex representation of
        symbolic functions:

            (1) latex_name=LaTeX
                where ``LaTeX`` is any valid latex expression.
                Ex: f = function('f', latex_name="\\mathcal{F}")
                See EXAMPLES for more.

            (2) print_latex_func=my_latex_print
                where ``my_latex_print`` is any callable function
                that returns a valid latex expression.
                Ex: f = function('f', print_latex_func=my_latex_print)
                See EXAMPLES for an explicit usage.

    .. NOTE::

       The new function is both returned and automatically injected
       into the global namespace.  If you use this function in library
       code, it is better to use sage.symbolic.function_factory.function,
       since it won't touch the global namespace.

    EXAMPLES:

    We create a formal function called supersin ::

        sage: function('supersin')
        supersin

    We can immediately use supersin in symbolic expressions::

        sage: y, z, A = var('y z A')
        sage: supersin(y+z) + A^3
        A^3 + supersin(y + z)

    We can define other functions in terms of supersin::

        sage: g(x,y) = supersin(x)^2 + sin(y/2)
        sage: g
        (x, y) |--> supersin(x)^2 + sin(1/2*y)
        sage: g.diff(y)
        (x, y) |--> 1/2*cos(1/2*y)
        sage: k = g.diff(x); k
        (x, y) |--> 2*supersin(x)*diff(supersin(x), x)

    Custom typesetting of symbolic functions in LaTeX, either using latex_name
    keyword::

        sage: function('riemann', latex_name="\\mathcal{R}")
        riemann
        sage: latex(riemann(x))
        \mathcal{R}\left(x\right)

    or passing a custom callable function that returns a latex expression::

        sage: mu,nu = var('mu,nu')
        sage: def my_latex_print(self, *args): return "\\psi_{%s}"%(', '.join(map(latex, args)))
        sage: function('psi', print_latex_func=my_latex_print)
        psi
        sage: latex(psi(mu,nu))
        \psi_{\mu, \nu}

    In Sage 4.0, you must now use the :meth:`substitute_function`
    method to replace functions::

        sage: k.substitute_function(supersin, sin)
        2*cos(x)*sin(x)
        
    TESTS:

    Make sure that :trac:`15860` is fixed and whitespaces are removed::
    
        sage: function('A, B')
        (A, B)
        sage: B
        B   
    """
    if len(args) > 0:
        from sage.misc.superseded import deprecation
        deprecation(17447, "Calling function('f',x) is deprecated. Use function('f')(x) instead.")
        return function(s, **kwds)(*args)

    G = globals()  # this is the reason the code must be in Cython.
    v = new_function(s, **kwds)
    if isinstance(v, tuple):
        for x in v:
            G[repr(x)] = x
    else:
        G[repr(v)] = v
    return v


def clear_vars():
    """
    Delete all 1-letter symbolic variables that are predefined at
    startup of Sage.  Any one-letter global variables that are not
    symbolic variables are not cleared.

    EXAMPLES::

        sage: var('x y z')
        (x, y, z)
        sage: (x+y)^z
        (x + y)^z
        sage: k = 15
        sage: clear_vars()
        sage: (x+y)^z
        Traceback (most recent call last):
        ...
        NameError: name 'x' is not defined
        sage: expand((e + i)^2)
        e^2 + 2*I*e - 1
        sage: k
        15
    """
    G = globals()
    from sage.symbolic.ring import is_SymbolicVariable
    for i in list(range(65, 65 + 26)) + list(range(97, 97 + 26)):
        if chr(i) in G and is_SymbolicVariable(G[chr(i)]):
            # We check to see if there is a corresponding pyobject
            # associated with the expression.  This will work for
            # constants which we want to keep, but will fail for
            # variables that we want to delete.
            try:
                G[chr(i)].pyobject()
            except TypeError:
                del G[chr(i)]
