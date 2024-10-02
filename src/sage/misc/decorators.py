# sage_setup: distribution = sagemath-objects
"""
Decorators

Python decorators for use in Sage.

AUTHORS:

- Tim Dumol (5 Dec 2009) -- initial version.
- Johan S. R. Nielsen (2010) -- collect decorators from various modules.
- Johan S. R. Nielsen (8 apr 2011) -- improve introspection on decorators.
- Simon King (2011-05-26) -- improve introspection of sage_wraps. Put this
  file into the reference manual.
- Julian Rueth (2014-03-19): added ``decorator_keywords`` decorator
"""
# *****************************************************************************
#       Copyright (C) 2009 Tim Dumol
#                     2010,2011 Johan S. R. Nielsen
#                     2011 Simon King <simon.king@uni-jena.de>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from functools import (partial, update_wrapper, WRAPPER_ASSIGNMENTS,
                       WRAPPER_UPDATES)
from copy import copy

from sage.misc.sageinspect import (sage_getsource, sage_getsourcelines,
                                   sage_getargspec)

from inspect import FullArgSpec


def sage_wraps(wrapped, assigned=WRAPPER_ASSIGNMENTS, updated=WRAPPER_UPDATES):
    r"""
    Decorator factory which should be used in decorators for making sure that
    meta-information on the decorated callables are retained through the
    decorator, such that the introspection functions of
    ``sage.misc.sageinspect`` retrieves them correctly. This includes
    documentation string, source, and argument specification. This is an
    extension of the Python standard library decorator functools.wraps.

    That the argument specification is retained from the decorated functions
    implies, that if one uses ``sage_wraps`` in a decorator which intentionally
    changes the argument specification, one should add this information to
    the special attribute ``_sage_argspec_`` of the wrapping function (for an
    example, see e.g. ``@options`` decorator in this module).

    EXAMPLES:

    Demonstrate that documentation string and source are retained from the
    decorated function::

        sage: def square(f):
        ....:     @sage_wraps(f)
        ....:     def new_f(x):
        ....:         return f(x)*f(x)
        ....:     return new_f
        sage: @square
        ....: def g(x):
        ....:     "My little function"
        ....:     return x
        sage: g(2)
        4
        sage: g(x)                                                                      # needs sage.symbolic
        x^2
        sage: g.__doc__
        'My little function'
        sage: from sage.misc.sageinspect import sage_getsource, sage_getsourcelines, sage_getfile
        sage: sage_getsource(g)
        '@square...def g(x)...'

    Demonstrate that the argument description are retained from the
    decorated function through the special method (when left
    unchanged) (see :issue:`9976`)::

        sage: def diff_arg_dec(f):
        ....:     @sage_wraps(f)
        ....:     def new_f(y, some_def_arg=2):
        ....:         return f(y+some_def_arg)
        ....:     return new_f
        sage: @diff_arg_dec
        ....: def g(x):
        ....:     return x
        sage: g(1)
        3
        sage: g(1, some_def_arg=4)
        5
        sage: from sage.misc.sageinspect import sage_getargspec
        sage: sage_getargspec(g)
        FullArgSpec(args=['x'], varargs=None, varkw=None, defaults=None,
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    Demonstrate that it correctly gets the source lines and the source
    file, which is essential for interactive code edition; note that we
    do not test the line numbers, as they may easily change::

        sage: P.<x,y> = QQ[]
        sage: I = P*[x,y]
        sage: sage_getfile(I.interreduced_basis)       # known bug
        '.../sage/rings/polynomial/multi_polynomial_ideal.py'
        sage: sage_getsourcelines(I.interreduced_basis)                                 # needs sage.libs.singular
        (['    @handle_AA_and_QQbar\n',
          '    @singular_gb_standard_options\n',
          '    @libsingular_gb_standard_options\n',
          '    def interreduced_basis(self):\n',
          ...
          '        return self.basis.reduced()\n'], ...)

    The ``f`` attribute of the decorated function refers to the
    original function::

        sage: foo = object()
        sage: @sage_wraps(foo)
        ....: def func():
        ....:     pass
        sage: wrapped = sage_wraps(foo)(func)
        sage: wrapped.f is foo
        True

    Demonstrate that sage_wraps works for non-function callables
    (:issue:`9919`)::

        sage: def square_for_met(f):
        ....:   @sage_wraps(f)
        ....:   def new_f(self, x):
        ....:       return f(self,x)*f(self,x)
        ....:   return new_f
        sage: class T:
        ....:   @square_for_met
        ....:   def g(self, x):
        ....:       "My little method"
        ....:       return x
        sage: t = T()
        sage: t.g(2)
        4
        sage: t.g.__doc__
        'My little method'

    The bug described in :issue:`11734` is fixed::

        sage: def square(f):
        ....:     @sage_wraps(f)
        ....:     def new_f(x):
        ....:         return f(x)*f(x)
        ....:     return new_f
        sage: f = lambda x:x^2
        sage: g = square(f)
        sage: g(3)  # this line used to fail for some people if these command were manually entered on the sage prompt
        81
    """
    # TRAC 9919: Workaround for bug in @update_wrapper when used with
    # non-function callables.
    assigned = set(assigned).intersection(set(dir(wrapped)))
    # end workaround

    def f(wrapper, assigned=assigned, updated=updated):
        update_wrapper(wrapper, wrapped, assigned=assigned, updated=updated)
        # For backwards-compatibility with old versions of sage_wraps
        wrapper.f = wrapped
        # For forwards-compatibility with functools.wraps on Python 3
        wrapper.__wrapped__ = wrapped
        wrapper._sage_src_ = lambda: sage_getsource(wrapped)
        wrapper._sage_src_lines_ = lambda: sage_getsourcelines(wrapped)
        # Getting the signature right in documentation by Sphinx (Issue 9976)
        # The attribute _sage_argspec_() is read by Sphinx if present and used
        # as the argspec of the function instead of using reflection.
        wrapper._sage_argspec_ = lambda: sage_getargspec(wrapped)
        return wrapper
    return f


# Infix operator decorator
class infix_operator:
    """
    A decorator for functions which allows for a hack that makes
    the function behave like an infix operator.

    This decorator exists as a convenience for interactive use.

    EXAMPLES:

    An infix dot product operator::

        sage: @infix_operator('multiply')
        ....: def dot(a, b):
        ....:     '''Dot product.'''
        ....:     return a.dot_product(b)
        sage: u = vector([1, 2, 3])                                                     # needs sage.modules
        sage: v = vector([5, 4, 3])                                                     # needs sage.modules
        sage: u *dot* v                                                                 # needs sage.modules
        22

    An infix element-wise addition operator::

        sage: # needs sage.modules
        sage: @infix_operator('add')
        ....: def eadd(a, b):
        ....:   return a.parent([i + j for i, j in zip(a, b)])
        sage: u = vector([1, 2, 3])
        sage: v = vector([5, 4, 3])
        sage: u +eadd+ v
        (6, 6, 6)
        sage: 2*u +eadd+ v
        (7, 8, 9)

    A hack to simulate a postfix operator::

        sage: @infix_operator('or')
        ....: def thendo(a, b):
        ....:     return b(a)
        sage: x |thendo| cos |thendo| (lambda x: x^2)                                   # needs sage.symbolic
        cos(x)^2
    """

    operators = {
        'add': {'left': '__add__', 'right': '__radd__'},
        'multiply': {'left': '__mul__', 'right': '__rmul__'},
        'or': {'left': '__or__', 'right': '__ror__'},
    }

    def __init__(self, precedence):
        """
        INPUT:

        - ``precedence`` -- one of ``'add'``, ``'multiply'``, or ``'or'``
          indicating the new operator's precedence in the order of operations
        """
        self.precedence = precedence

    def __call__(self, func):
        """Returns a function which acts as an inline operator."""

        left_meth = self.operators[self.precedence]['left']
        right_meth = self.operators[self.precedence]['right']
        wrapper_name = func.__name__
        wrapper_members = {
            'function': staticmethod(func),
            left_meth: _infix_wrapper._left,
            right_meth: _infix_wrapper._right,
            '_sage_src_': lambda: sage_getsource(func)
        }
        for attr in WRAPPER_ASSIGNMENTS:
            try:
                wrapper_members[attr] = getattr(func, attr)
            except AttributeError:
                pass

        wrapper = type(wrapper_name, (_infix_wrapper,), wrapper_members)

        wrapper_inst = wrapper()
        wrapper_inst.__dict__.update(getattr(func, '__dict__', {}))
        return wrapper_inst


class _infix_wrapper:
    function = None

    def __init__(self, left=None, right=None):
        """
        Initialize the actual infix object, with possibly a specified left
        and/or right operand.
        """
        self.left = left
        self.right = right

    def __call__(self, *args, **kwds):
        """Call the passed function."""
        return self.function(*args, **kwds)

    def _left(self, right):
        """The function for the operation on the left (e.g., __add__)."""
        if self.left is None:
            if self.right is None:
                new = copy(self)
                new.right = right
                return new
            else:
                raise SyntaxError("Infix operator already has its "
                                  "right argument")
        else:
            return self.function(self.left, right)

    def _right(self, left):
        """The function for the operation on the right (e.g., __radd__)."""
        if self.right is None:
            if self.left is None:
                new = copy(self)
                new.left = left
                return new
            else:
                raise SyntaxError("Infix operator already has its "
                                  "left argument")
        else:
            return self.function(left, self.right)


def decorator_defaults(func):
    """
    This function allows a decorator to have default arguments.

    Normally, a decorator can be called with or without arguments.
    However, the two cases call for different types of return values.
    If a decorator is called with no parentheses, it should be run
    directly on the function.  However, if a decorator is called with
    parentheses (i.e., arguments), then it should return a function
    that is then in turn called with the defined function as an
    argument.

    This decorator allows us to have these default arguments without
    worrying about the return type.

    EXAMPLES::

        sage: from sage.misc.decorators import decorator_defaults
        sage: @decorator_defaults
        ....: def my_decorator(f,*args,**kwds):
        ....:   print(kwds)
        ....:   print(args)
        ....:   print(f.__name__)

        sage: @my_decorator
        ....: def my_fun(a,b):
        ....:   return a,b
        {}
        ()
        my_fun
        sage: @my_decorator(3,4,c=1,d=2)
        ....: def my_fun(a,b):
        ....:   return a,b
        {'c': 1, 'd': 2}
        (3, 4)
        my_fun
    """
    @sage_wraps(func)
    def my_wrap(*args, **kwds):
        if len(kwds) == 0 and len(args) == 1:
            # call without parentheses
            return func(*args)
        else:
            return lambda f: func(f, *args, **kwds)
    return my_wrap


class suboptions:
    def __init__(self, name, **options):
        """
        A decorator for functions which collects all keywords
        starting with ``name+'_'`` and collects them into a dictionary
        which will be passed on to the wrapped function as a
        dictionary called ``name_options``.

        The keyword arguments passed into the constructor are taken
        to be default for the ``name_options`` dictionary.

        EXAMPLES::

            sage: from sage.misc.decorators import suboptions
            sage: s = suboptions('arrow', size=2)
            sage: s.name
            'arrow_'
            sage: s.options
            {'size': 2}
        """
        self.name = name + "_"
        self.options = options

    def __call__(self, func):
        """
        Return a wrapper around ``func``.

        EXAMPLES::

            sage: from sage.misc.decorators import suboptions
            sage: def f(*args, **kwds): print(sorted(kwds.items()))
            sage: f = suboptions('arrow', size=2)(f)
            sage: f(size=2)
            [('arrow_options', {'size': 2}), ('size', 2)]
            sage: f(arrow_size=3)
            [('arrow_options', {'size': 3})]
            sage: f(arrow_options={'size':4})
            [('arrow_options', {'size': 4})]
            sage: f(arrow_options={'size':4}, arrow_size=5)
            [('arrow_options', {'size': 5})]

         Demonstrate that the introspected argument specification of the
         wrapped function is updated (see :issue:`9976`)::

            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(f)
            FullArgSpec(args=['arrow_size'], varargs='args', varkw='kwds', defaults=(2,),
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
        """
        @sage_wraps(func)
        def wrapper(*args, **kwds):
            suboptions = copy(self.options)
            suboptions.update(kwds.pop(self.name+"options", {}))

            # Collect all the relevant keywords in kwds
            # and put them in suboptions
            for key, value in list(kwds.items()):
                if key.startswith(self.name):
                    suboptions[key[len(self.name):]] = value
                    del kwds[key]

            kwds[self.name + "options"] = suboptions

            return func(*args, **kwds)

        # Add the options specified by @options to the signature of the wrapped
        # function in the Sphinx-generated documentation (Issue 9976), using the
        # special attribute _sage_argspec_ (see e.g. sage.misc.sageinspect)
        def argspec():
            argspec = sage_getargspec(func)

            def listForNone(l):
                return l if l is not None else []
            newArgs = [self.name + opt for opt in self.options.keys()]
            args = (argspec.args if argspec.args is not None else []) + newArgs
            defaults = (argspec.defaults if argspec.defaults is not None else ()) \
                + tuple(self.options.values())
            # Note: argspec.defaults is not always a tuple for some reason
            return FullArgSpec(args, argspec.varargs, argspec.varkw, defaults,
                               kwonlyargs=[], kwonlydefaults=None, annotations={})
        wrapper._sage_argspec_ = argspec

        return wrapper


class options:
    def __init__(self, **options):
        """
        A decorator for functions which allows for default options to be
        set and reset by the end user.  Additionally, if one needs to, one
        can get at the original keyword arguments passed into the
        decorator.

        TESTS::

            sage: from sage.misc.decorators import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: o.options
            {'rgbcolor': (0, 0, 1)}
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: o.original_opts
            True
            sage: loads(dumps(o)).options
            {'rgbcolor': (0, 0, 1)}

        Demonstrate that the introspected argument specification of the wrapped
        function is updated (see :issue:`9976`)::

            sage: from sage.misc.decorators import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, sorted(kwds.items())))
            sage: f1 = o(f)
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(f1)
            FullArgSpec(args=['rgbcolor'], varargs='args', varkw='kwds', defaults=((0, 0, 1),),
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
        """
        self.options = options
        self.original_opts = options.pop('__original_opts', False)

    def __call__(self, func):
        """
        EXAMPLES::

            sage: from sage.misc.decorators import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, sorted(kwds.items())))
            sage: f1 = o(f)
            sage: f1()
            () [('rgbcolor', (0, 0, 1))]
            sage: f1(rgbcolor=1)
            () [('rgbcolor', 1)]
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: f2 = o(f)
            sage: f2(alpha=1)
            () [('__original_opts', {'alpha': 1}), ('alpha', 1), ('rgbcolor', (0, 0, 1))]
        """
        @sage_wraps(func)
        def wrapper(*args, **kwds):
            options = copy(wrapper.options)
            if self.original_opts:
                options['__original_opts'] = kwds
            options.update(kwds)
            return func(*args, **options)

        # Add the options specified by @options to the signature of the wrapped
        # function in the Sphinx-generated documentation (Issue 9976), using the
        # special attribute _sage_argspec_ (see e.g. sage.misc.sageinspect)
        def argspec():
            argspec = sage_getargspec(func)
            args = ((argspec.args if argspec.args is not None else []) +
                    list(self.options))
            defaults = (argspec.defaults or ()) + tuple(self.options.values())
            # Note: argspec.defaults is not always a tuple for some reason
            return FullArgSpec(args, argspec.varargs, argspec.varkw, defaults,
                               kwonlyargs=[], kwonlydefaults=None, annotations={})

        wrapper._sage_argspec_ = argspec

        def defaults():
            """
            Return the default options.

            EXAMPLES::

                sage: from sage.misc.decorators import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds):
                ....:     print("{} {}".format(args, sorted(kwds.items())))
                sage: f = o(f)
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.defaults()
                {'rgbcolor': (0, 0, 1)}
            """
            return copy(self.options)

        def reset():
            """
            Reset the options to the defaults.

            EXAMPLES::

                sage: from sage.misc.decorators import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds):
                ....:     print("{} {}".format(args, sorted(kwds.items())))
                sage: f = o(f)
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.options
                {'rgbcolor': (1, 1, 1)}
                sage: f()
                () [('rgbcolor', (1, 1, 1))]
                sage: f.reset()
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f()
                () [('rgbcolor', (0, 0, 1))]
            """
            wrapper.options = copy(self.options)

        wrapper.options = copy(self.options)
        wrapper.reset = reset
        wrapper.reset.__doc__ = """
        Reset the options to the defaults.

        Defaults:
        %s
        """ % self.options

        wrapper.defaults = defaults
        wrapper.defaults.__doc__ = """
        Return the default options.

        Defaults:
        %s
        """ % self.options

        return wrapper


class rename_keyword:
    def __init__(self, deprecated=None, deprecation=None, **renames):
        """
        A decorator which renames keyword arguments and optionally
        deprecates the new keyword.

        INPUT:

        - ``deprecation`` -- integer; the github issue number where the
          deprecation was introduced

        - the rest of the arguments is a list of keyword arguments in the
          form ``renamed_option='existing_option'``.  This will have the
          effect of renaming ``renamed_option`` so that the function only
          sees ``existing_option``.  If both ``renamed_option`` and
          ``existing_option`` are passed to the function, ``existing_option``
          will override the ``renamed_option`` value.

        EXAMPLES::

            sage: from sage.misc.decorators import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: r.renames
            {'color': 'rgbcolor'}
            sage: loads(dumps(r)).renames
            {'color': 'rgbcolor'}

        To deprecate an old keyword::

            sage: r = rename_keyword(deprecation=13109, color='rgbcolor')
        """
        assert deprecated is None, 'Use @rename_keyword(deprecation=<issue_number>, ...)'
        self.renames = renames
        self.deprecation = deprecation

    def __call__(self, func):
        """
        Rename keywords.

        EXAMPLES::

            sage: from sage.misc.decorators import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, kwds))
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(rgbcolor=1)
            () {'rgbcolor': 1}
            sage: f(color=1)
            () {'rgbcolor': 1}

        We can also deprecate the renamed keyword::

            sage: r = rename_keyword(deprecation=13109, deprecated_option='new_option')
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, kwds))
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(new_option=1)
            () {'new_option': 1}
            sage: f(deprecated_option=1)
            doctest:...: DeprecationWarning: use the option 'new_option' instead of 'deprecated_option'
            See https://github.com/sagemath/sage/issues/13109 for details.
            () {'new_option': 1}
        """
        @sage_wraps(func)
        def wrapper(*args, **kwds):
            for old_name, new_name in self.renames.items():
                if old_name in kwds and new_name not in kwds:
                    if self.deprecation is not None:
                        from sage.misc.superseded import deprecation
                        deprecation(self.deprecation, "use the option "
                                    "%r instead of %r" % (new_name, old_name))
                    kwds[new_name] = kwds[old_name]
                    del kwds[old_name]
            return func(*args, **kwds)

        return wrapper


class specialize:
    r"""
    A decorator generator that returns a decorator that in turn
    returns a specialized function for function ``f``. In other words,
    it returns a function that acts like ``f`` with arguments
    ``*args`` and ``**kwargs`` supplied.

    INPUT:

    - ``*args``, ``**kwargs`` -- arguments to specialize the function for

    OUTPUT: a decorator that accepts a function ``f`` and specializes it
    with ``*args`` and ``**kwargs``

    EXAMPLES::

        sage: f = specialize(5)(lambda x, y: x+y)
        sage: f(10)
        15
        sage: f(5)
        10
        sage: @specialize("Bon Voyage")
        ....: def greet(greeting, name):
        ....:     print("{0}, {1}!".format(greeting, name))
        sage: greet("Monsieur Jean Valjean")
        Bon Voyage, Monsieur Jean Valjean!
        sage: greet(name = 'Javert')
        Bon Voyage, Javert!
    """
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, f):
        return sage_wraps(f)(partial(f, *self.args, **self.kwargs))


def decorator_keywords(func):
    r"""
    A decorator for decorators with optional keyword arguments.

    EXAMPLES::

        sage: from sage.misc.decorators import decorator_keywords
        sage: @decorator_keywords
        ....: def preprocess(f=None, processor=None):
        ....:     def wrapper(*args, **kwargs):
        ....:         if processor is not None:
        ....:             args, kwargs = processor(*args, **kwargs)
        ....:         return f(*args, **kwargs)
        ....:     return wrapper

    This decorator can be called with and without arguments::

        sage: @preprocess
        ....: def foo(x): return x
        sage: foo(None)
        sage: foo(1)
        1

        sage: def normalize(x): return ((0,),{}) if x is None else ((x,),{})
        sage: @preprocess(processor=normalize)
        ....: def foo(x): return x
        sage: foo(None)
        0
        sage: foo(1)
        1
    """
    @sage_wraps(func)
    def wrapped(f=None, **kwargs):
        if f is None:
            return sage_wraps(func)(lambda f: func(f, **kwargs))
        else:
            return func(f, **kwargs)
    return wrapped
