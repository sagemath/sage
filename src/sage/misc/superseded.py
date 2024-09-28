# sage_setup: distribution = sagemath-objects
"""
Handling Superseded Functionality

The main mechanism in Sage to deal with superseded functionality is to
add a deprecation warning. This will be shown once, the first time
that the deprecated function is called.

Note that all doctests in the following use the github issue number
:issue:`13109`, which is where this mandatory argument to
:func:`deprecation` was introduced.

Functions and classes
---------------------
"""


########################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

from warnings import warn
import inspect

from sage.misc.lazy_attribute import lazy_attribute


def _check_issue_number(issue_number):
    """
    Check that the argument is likely to be a valid github issue number.

    INPUT:

    - ``issue_number`` -- anything

    OUTPUT:

    This function returns nothing. A :exc:`ValueError` or :exc:`TypeError`
    is raised if the argument cannot be a valid issue number.

    EXAMPLES::

        sage: from sage.misc.superseded import _check_issue_number
        sage: _check_issue_number(1)
        sage: _check_issue_number(0)
        Traceback (most recent call last):
        ...
        ValueError: 0 is not a valid github issue number
        sage: _check_issue_number(int(10))
        sage: _check_issue_number(10.0)
        Traceback (most recent call last):
        ...
        TypeError: 10.0000000000000 is not a valid github issue number
        sage: _check_issue_number('10')
        Traceback (most recent call last):
        ...
        TypeError: '10' is not a valid github issue number
    """
    try:
        issue_number = issue_number.__index__()
    except Exception:
        raise TypeError('%r is not a valid github issue number' % issue_number)
    if issue_number <= 0:
        raise ValueError('%r is not a valid github issue number' % issue_number)


def deprecation(issue_number, message, stacklevel=4):
    r"""
    Issue a deprecation warning.

    INPUT:

    - ``issue_number`` -- integer; the github issue number where the
      deprecation is introduced

    - ``message`` -- string; an explanation why things are deprecated
      and by what it should be replaced

    - ``stack_level`` -- integer (default: `4`); this is passed on to
      :func:`warnings.warn`

    EXAMPLES::

        sage: def foo():
        ....:  sage.misc.superseded.deprecation(13109, 'the function foo is replaced by bar')
        sage: foo()
        doctest:...: DeprecationWarning: the function foo is replaced by bar
        See https://github.com/sagemath/sage/issues/13109 for details.

    .. SEEALSO::

        :func:`experimental`,
        :func:`warning`.
    """
    warning(issue_number, message, DeprecationWarning, stacklevel)


def deprecation_cython(issue_number, message, stacklevel=3):
    r"""
    Issue a deprecation warning -- for use in cython functions.

    TESTS:

    We check that `deprecation_cython` in a cython function generates a warning
    with the same callsite reference as `deprecation` in a python function, whereas
    `deprecation` in a cython function does not::

        sage: # needs sage.misc.cython
        sage: cython(
        ....: '''
        ....: from sage.misc.superseded import deprecation_cython, deprecation
        ....: def foo1():
        ....:     deprecation_cython(100, "boo")
        ....: def foo2():
        ....:     deprecation(100, "boo")
        ....: ''')
        sage: def foo3():
        ....:     deprecation(100, "boo")
        sage: if True:  # Execute the three "with" blocks as one doctest
        ....:     with warnings.catch_warnings(record=True) as w1:
        ....:        warnings.simplefilter("always")
        ....:        foo1()
        ....:     with warnings.catch_warnings(record=True) as w2:
        ....:        warnings.simplefilter("always")
        ....:        foo2()
        ....:     with warnings.catch_warnings(record=True) as w3:
        ....:        warnings.simplefilter("always")
        ....:        foo3()
        sage: w1[0].filename == w3[0].filename
        True
        sage: w2[0].filename == w3[0].filename
        False
    """
    warning(issue_number, message, DeprecationWarning, stacklevel)


def warning(issue_number, message, warning_class=Warning, stacklevel=3):
    r"""
    Issue a warning.

    INPUT:

    - ``issue_number`` -- integer; the github issue number where the
      deprecation is introduced

    - ``message`` -- string; an explanation what is going on

    - ``warning_class`` -- (default: ``Warning``) a class inherited
      from a Python :class:`~exceptions.Warning`

    - ``stack_level`` -- integer (default: `3`); this is passed on to
      :func:`warnings.warn`

    EXAMPLES::

        sage: def foo():
        ....:     sage.misc.superseded.warning(
        ....:         99999,
        ....:         'The syntax will change in future.',
        ....:         FutureWarning)
        sage: foo()
        doctest:...: FutureWarning: The syntax will change in future.
        See https://github.com/sagemath/sage/issues/99999 for details.

    .. SEEALSO::

        :func:`deprecation`,
        :func:`experimental`,
        :class:`exceptions.Warning`.
    """
    _check_issue_number(issue_number)
    message += '\nSee https://github.com/sagemath/sage/issues/' + str(issue_number) + ' for details.'

    # Stack level 3 to get the line number of the code which called
    # the deprecated function which called this function.
    warn(message, warning_class, stacklevel)


def experimental_warning(issue_number, message, stacklevel=4):
    r"""
    Issue a warning that the functionality or class is experimental
    and might change in future.

    INPUT:

    - ``issue_number`` -- integer; the github issue number where the
      experimental functionality was introduced

    - ``message`` -- string; an explanation what is going on

    - ``stack_level`` -- integer (default: `4`); this is passed on to
      :func:`warnings.warn`

    EXAMPLES::

        sage: def foo():
        ....:    sage.misc.superseded.experimental_warning(
        ....:        66666, 'This function is experimental and '
        ....:               'might change in future.')
        sage: foo()
        doctest:...: FutureWarning: This function is experimental and
        might change in future.
        See https://github.com/sagemath/sage/issues/66666 for details.

    .. SEEALSO::

        :class:`mark_as_experimental`,
        :func:`warning`,
        :func:`deprecation`.
    """
    warning(issue_number, message, FutureWarning, stacklevel)


class experimental:
    def __init__(self, issue_number, stacklevel=4):
        """
        A decorator which warns about the experimental/unstable status of
        the decorated class/method/function.

        INPUT:

        - ``issue_number`` -- integer; the github issue number where this
          code was introduced

        - ``stack_level`` -- integer (default: `4`); this is passed on to
          :func:`warnings.warn`

        EXAMPLES::

            sage: @sage.misc.superseded.experimental(issue_number=79997)
            ....: def foo(*args, **kwargs):
            ....:     print("{} {}".format(args, kwargs))
            sage: foo(7, what='Hello')
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/79997 for details.
            (7,) {'what': 'Hello'}

        ::

            sage: class bird(SageObject):
            ....:     @sage.misc.superseded.experimental(issue_number=99999)
            ....:     def __init__(self, *args, **kwargs):
            ....:         print("piep {} {}".format(args, kwargs))
            sage: _ = bird(99)
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/99999 for details.
            piep (99,) {}

        TESTS:

        The following test works together with the doc-test for
        :meth:`__experimental_self_test` to demonstrate that warnings are issued only
        once, even in doc-tests (see :issue:`20601`).
        ::

            sage: from sage.misc.superseded import __experimental_self_test
            sage: _ = __experimental_self_test("A")
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/88888 for details.
            I'm A

        .. SEEALSO::

            :func:`experimental`,
            :func:`warning`,
            :func:`deprecation`.
        """
        self.issue_number = issue_number
        self.stacklevel = stacklevel

    def __call__(self, func):
        """
        Print experimental warning.

        INPUT:

        - ``func`` -- the function to decorate

        OUTPUT: the wrapper to this function

        TESTS::

            sage: def foo(*args, **kwargs):
            ....:     print("{} {}".format(args, kwargs))
            sage: from sage.misc.superseded import experimental
            sage: ex_foo = experimental(issue_number=99399)(foo)
            sage: ex_foo(3, what='Hello')
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/99399 for details.
            (3,) {'what': 'Hello'}
        """
        from sage.misc.decorators import sage_wraps

        @sage_wraps(func)
        def wrapper(*args, **kwds):
            if not wrapper._already_issued:
                experimental_warning(self.issue_number,
                            'This class/method/function is marked as '
                            'experimental. It, its functionality or its '
                            'interface might change without a '
                            'formal deprecation.',
                            self.stacklevel)
                wrapper._already_issued = True
            return func(*args, **kwds)
        wrapper._already_issued = False

        return wrapper


class __experimental_self_test:
    r"""
    This is a class only to demonstrate with a doc-test that the @experimental
    decorator only issues a warning message once (see :issue:`20601`).

    The test below does not issue a warning message because that warning has
    already been issued by a previous doc-test in the @experimental code. Note
    that this behaviour cannot be demonstrated within a single documentation
    string: Sphinx will itself supress multiple issued warnings.

    TESTS::

        sage: from sage.misc.superseded import __experimental_self_test
        sage: _ = __experimental_self_test("B")
        I'm B
    """
    @experimental(issue_number=88888)
    def __init__(self, x):
        print("I'm " + x)


class DeprecatedFunctionAlias:
    """
    A wrapper around methods or functions which automatically prints a
    deprecation message. See :func:`deprecated_function_alias`.

    AUTHORS:

    - Florent Hivert (2009-11-23), with the help of Mike Hansen.
    - Luca De Feo (2011-07-11), printing the full module path when different from old path
    """
    def __init__(self, issue_number, func, module, instance=None, unbound=None):
        r"""
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g = deprecated_function_alias(13109, number_of_partitions)            # needs sage.combinat
            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g.__doc__                                                             # needs sage.combinat
            'Deprecated: Use :func:`number_of_partitions` instead.\nSee :issue:`13109` for details.\n\n'
        """
        _check_issue_number(issue_number)
        try:
            self.__dict__.update(func.__dict__)
        except AttributeError:
            pass  # Cython classes don't have __dict__
        self.func = func
        self.issue_number = issue_number
        self.instance = instance  # for use with methods
        self.unbound = unbound
        self.__module__ = module
        if isinstance(func, type(deprecation)):
            sphinxrole = "func"
        else:
            sphinxrole = "meth"
        doc = 'Deprecated: '
        doc += 'Use :' + sphinxrole + ':`' + self.func.__name__ + '` instead.\n'
        doc += 'See :issue:`' + str(self.issue_number) + '` for details.\n\n'
        self.__doc__ = doc

    @lazy_attribute
    def __name__(self):
        r"""
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g = deprecated_function_alias(13109, number_of_partitions)            # needs sage.combinat
            sage: g.__name__                                                            # needs sage.combinat
            'g'

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class cls():
            ....:    def new_meth(self): return 42
            ....:    old_meth = deprecated_function_alias(13109, new_meth)
            sage: cls.old_meth.__name__
            'old_meth'
            sage: cls().old_meth.__name__
            'old_meth'

            sage: cython('\n'.join([                                                    # needs sage.misc.cython
            ....:     r"from sage.misc.superseded import deprecated_function_alias",
            ....:     r"cdef class cython_cls():",
            ....:     r"    def new_cython_meth(self):",
            ....:     r"        return 1",
            ....:     r"    old_cython_meth = deprecated_function_alias(13109, new_cython_meth)"
            ....: ]))
            sage: cython_cls().old_cython_meth.__name__                                 # needs sage.misc.cython
            'old_cython_meth'
        """
        # first look through variables in stack frames
        for frame in inspect.stack():
            for name, obj in frame[0].f_globals.items():
                if obj is self:
                    return name
        # then search object that contains self as method
        import gc
        import copy
        gc.collect()

        def is_class(gc_ref):
            if not isinstance(gc_ref, dict):
                return False
            is_python_class = '__module__' in gc_ref or '__package__' in gc_ref
            is_cython_class = '__new__' in gc_ref
            return is_python_class or is_cython_class
        search_for = self if (self.unbound is None) else self.unbound
        for ref in gc.get_referrers(search_for):
            if is_class(ref) and ref is not self.__dict__:
                ref_copy = copy.copy(ref)
                for key, val in ref_copy.items():
                    if val is search_for:
                        return key
        raise AttributeError("The name of this deprecated function cannot be determined")

    def __call__(self, *args, **kwds):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: def bla(): return 42
            sage: blo = deprecated_function_alias(13109, bla)
            sage: blo()
            doctest:...: DeprecationWarning: blo is deprecated. Please use bla instead.
            See https://github.com/sagemath/sage/issues/13109 for details.
            42
        """
        if self.instance is None and self.__module__ != self.func.__module__:
            other = self.func.__module__ + "." + self.func.__name__
        else:
            other = self.func.__name__

        deprecation(self.issue_number,
                    "{} is deprecated. Please use {} instead.".format(
                        self.__name__, other))
        if self.instance is None:
            return self.func(*args, **kwds)
        else:
            return self.func(self.instance, *args, **kwds)

    def __get__(self, inst, cls=None):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class cls():
            ....:    def new_meth(self): return 42
            ....:    old_meth = deprecated_function_alias(13109, new_meth)
            sage: obj = cls()
            sage: obj.old_meth.instance is obj
            True

        :issue:`19125`::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class A:
            ....:    def __init__(self, x):
            ....:        self.x = x
            ....:    def f(self, y):
            ....:        return self.x+y
            ....:    g = deprecated_function_alias(42, f)
            sage: a1 = A(1)
            sage: a2 = A(2)
            sage: a1.g(a2.g(0))
            doctest:...: DeprecationWarning: g is deprecated. Please use f instead.
            See https://github.com/sagemath/sage/issues/42 for details.
            3
            sage: a1.f(a2.f(0))
            3
        """
        if inst is None:
            return self  # Unbound method lookup on class
        else:
            # Return a bound method wrapper
            return DeprecatedFunctionAlias(self.issue_number, self.func,
                                           self.__module__, instance=inst,
                                           unbound=self)


def deprecated_function_alias(issue_number, func):
    """
    Create an aliased version of a function or a method which raises a
    deprecation warning message.

    If f is a function or a method, write
    ``g = deprecated_function_alias(issue_number, f)``
    to make a deprecated aliased version of f.

    INPUT:

    - ``issue_number`` -- integer; the github issue number where the
      deprecation is introduced

    - ``func`` -- the function or method to be aliased

    EXAMPLES::

        sage: from sage.misc.superseded import deprecated_function_alias
        sage: g = deprecated_function_alias(13109, number_of_partitions)                # needs sage.combinat sage.libs.flint
        sage: g(5)                                                                      # needs sage.combinat sage.libs.flint
        doctest:...: DeprecationWarning: g is deprecated.
        Please use sage.combinat.partition.number_of_partitions instead.
        See https://github.com/sagemath/sage/issues/13109 for details.
        7

    This also works for methods::

        sage: class cls():
        ....:    def new_meth(self): return 42
        ....:    old_meth = deprecated_function_alias(13109, new_meth)
        sage: cls().old_meth()
        doctest:...: DeprecationWarning: old_meth is deprecated. Please use new_meth instead.
        See https://github.com/sagemath/sage/issues/13109 for details.
        42

    :issue:`11585`::

        sage: def a(): pass
        sage: b = deprecated_function_alias(13109, a)
        sage: b()
        doctest:...: DeprecationWarning: b is deprecated. Please use a instead.
        See https://github.com/sagemath/sage/issues/13109 for details.

    AUTHORS:

     - Florent Hivert (2009-11-23), with the help of Mike Hansen.
     - Luca De Feo (2011-07-11), printing the full module path when different from old path
    """
    module_name = None
    frame0 = inspect.currentframe()
    if frame0:
        frame1 = frame0.f_back
        if frame1:
            module_name = inspect.getmodulename(frame1.f_code.co_filename)
    if module_name is None:
        module_name = '__main__'
    return DeprecatedFunctionAlias(issue_number, func, module_name)
