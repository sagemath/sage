# sage_setup: distribution = sagemath-repl
r"""
Fixtures to help testing functionality

Utilities which modify or replace code to help with doctesting functionality.
Wrappers, proxies and mockups are typical examples of fixtures.

AUTHORS:

- Martin von Gagern (2014-12-15): AttributeAccessTracerProxy and trace_method
- Martin von Gagern (2015-01-02): Factor out TracerHelper and reproducible_repr

EXAMPLES:

You can use :func:`trace_method` to see how a method
communicates with its surroundings::

    sage: class Foo():
    ....:     def f(self):
    ....:         self.y = self.g(self.x)
    ....:     def g(self, arg):
    ....:         return arg + 1
    ....:
    sage: foo = Foo()
    sage: foo.x = 3
    sage: from sage.doctest.fixtures import trace_method
    sage: trace_method(foo, "f")
    sage: foo.f()
    enter f()
      read x = 3
      call g(3) -> 4
      write y = 4
    exit f -> None
"""

#*****************************************************************************
#       Copyright (C) 2014-2015 Martin von Gagern <Martin.vGagern@gmx.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import wraps


def reproducible_repr(val):
    r"""
    String representation of an object in a reproducible way.

    .. NOTE::

        This function is deprecated, in most cases it suffices to use
        the automatic sorting of dictionary keys and set items by a displayhook.
        See :func:`sage.doctest.forker.init_sage`.
        If used in a format string, use :func:`IPython.lib.pretty.pretty`.
        In the rare cases where the ordering of the elements is not reliable
        or transitive, ``sorted`` with a sane key can be used instead.

    This tries to ensure that the returned string does not depend on
    factors outside the control of the doctest.
    One example is the order of elements in a hash-based structure.
    For most objects, this is simply the ``repr`` of the object.

    All types for which special handling had been implemented are
    covered by the examples below. If a doctest requires special
    handling for additional types, this function may be extended
    appropriately. It is an error if an argument to this function has
    a non-reproducible ``repr`` implementation and is not explicitly
    mentioned in an example case below.

    INPUT:

    - ``val`` -- an object to be represented

    OUTPUT:

    A string representation of that object, similar to what ``repr``
    returns but for certain cases with more guarantees to ensure
    exactly the same result for semantically equivalent objects.

    EXAMPLES::

        sage: # not tested (test fails because of deprecation warning)
        sage: from sage.doctest.fixtures import reproducible_repr
        sage: print(reproducible_repr(set(["a", "c", "b", "d"])))
        set(['a', 'b', 'c', 'd'])
        sage: print(reproducible_repr(frozenset(["a", "c", "b", "d"])))
        frozenset(['a', 'b', 'c', 'd'])
        sage: print(reproducible_repr([1, frozenset("cab"), set("bar"), 0]))
        [1, frozenset(['a', 'b', 'c']), set(['a', 'b', 'r']), 0]
        sage: print(reproducible_repr({3.0: "three", "2": "two", 1: "one"}))            # needs sage.rings.real_mpfr
        {'2': 'two', 1: 'one', 3.00000000000000: 'three'}
        sage: print(reproducible_repr("foo\nbar"))  # demonstrate default case
        'foo\nbar'

    TESTS:

    Ensures deprecation warning is printed out::

        sage: from sage.doctest.fixtures import reproducible_repr
        sage: print(reproducible_repr(set(["a", "c", "b", "d"])))
        doctest:warning...
        DeprecationWarning: reproducible_repr is deprecated, see its documentation for details
        See https://github.com/sagemath/sage/issues/39420 for details.
        set(['a', 'b', 'c', 'd'])
    """
    from sage.misc.superseded import deprecation
    deprecation(39420, 'reproducible_repr is deprecated, see its documentation for details')

    def sorted_pairs(iterable, pairs=False):
        # We don't know whether container data structures will have
        # homogeneous types, and if not, whether comparisons will work
        # in a sane way. For this reason, we sort by representation first.
        res = sorted((reproducible_repr(item), item) for item in iterable)
        if not pairs:
            res = [r for r, i in res]
        return res

    if isinstance(val, frozenset):
        itms = sorted_pairs(val)
        return "frozenset([{}])".format(", ".join(itms))
    if isinstance(val, set):
        itms = sorted_pairs(val)
        return "set([{}])".format(", ".join(itms))
    if isinstance(val, dict):
        keys = sorted_pairs(val.keys(), True)
        itms = ["{}: {}".format(r, reproducible_repr(val[k])) for r, k in keys]
        return ("{{{}}}".format(", ".join(itms)))
    if isinstance(val, list):
        itms = map(reproducible_repr, val)
        return ("[{}]".format(", ".join(itms)))
    return repr(val)


class AttributeAccessTracerHelper:

    def __init__(self, delegate, prefix="  ", reads=True):
        r"""
        Helper to print proxied access to attributes.

        This class does the actual printing of access traces
        for objects proxied by :class:`AttributeAccessTracerProxy`.
        The fact that it's not a proxy at the same time
        helps avoiding complicated attribute access syntax.

        INPUT:

        - ``delegate`` -- the actual object to be proxied

        - ``prefix`` -- (default: ``"  "``) string to prepend to each printed output

        - ``reads`` -- (default: ``True``) whether to trace read access as well

        EXAMPLES::

            sage: class Foo():
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import AttributeAccessTracerHelper
            sage: pat = AttributeAccessTracerHelper(foo)
            sage: pat.set("x", 2)
              write x = 2
            sage: pat.get("x")
              read x = 2
            2
            sage: pat.get("f")(3)
              call f(3) -> 4
            4
        """
        self.delegate = delegate
        self.prefix = prefix
        self.reads = reads

    def get(self, name):
        r"""
        Read an attribute from the wrapped delegate object.

        If that value is a method (i.e. a callable object which is not
        contained in the dictionary of the object itself but instead
        inherited from some class) then it is replaced by a wrapper
        function to report arguments and return value.
        Otherwise an attribute read access is reported.

        EXAMPLES::

            sage: class Foo():
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: foo.x = 2
            sage: from sage.doctest.fixtures import AttributeAccessTracerHelper
            sage: pat = AttributeAccessTracerHelper(foo)
            sage: pat.get("x")
              read x = 2
            2
            sage: pat.get("f")(3)
              call f(3) -> 4
            4
        """
        val = getattr(self.delegate, name)
        from IPython.lib.pretty import pretty
        if callable(val) and name not in self.delegate.__dict__:
            @wraps(val)
            def wrapper(*args, **kwds):
                arglst = [pretty(arg) for arg in args]
                arglst.extend("{}={}".format(k, pretty(v))
                              for k, v in sorted(kwds.items()))
                res = val(*args, **kwds)
                print("{}call {}({}) -> {}"
                      .format(self.prefix, name, ", ".join(arglst),
                              pretty(res)))
                return res
            return wrapper
        else:
            if self.reads:
                print("{}read {} = {}".format(self.prefix, name,
                                              pretty(val)))
            return val

    def set(self, name, val):
        r"""
        Write an attribute to the wrapped delegate object.

        The name and new value are also reported in the output.

        EXAMPLES::

            sage: class Foo():
            ....:     pass
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import AttributeAccessTracerHelper
            sage: pat = AttributeAccessTracerHelper(foo)
            sage: pat.set("x", 2)
              write x = 2
            sage: foo.x
            2
        """
        from IPython.lib.pretty import pretty
        print("{}write {} = {}".format(self.prefix, name,
                                       pretty(val)))
        setattr(self.delegate, name, val)


class AttributeAccessTracerProxy:

    def __init__(self, delegate, **kwds):
        r"""
        Proxy object which prints all attribute and method access to an object.

        The implementation is kept lean since all access to attributes of
        the proxy itself requires complicated syntax.
        For this reason, the actual handling of attribute access
        is delegated to a :class:`AttributeAccessTracerHelper`.

        INPUT:

        - ``delegate`` -- the actual object to be proxied

        - ``prefix`` -- (default: ``"  "``)
          string to prepend to each printed output

        - ``reads`` -- (default: ``True``)
          whether to trace read access as well

        EXAMPLES::

            sage: class Foo():
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import AttributeAccessTracerProxy
            sage: pat = AttributeAccessTracerProxy(foo)
            sage: pat.x = 2
              write x = 2
            sage: pat.x
              read x = 2
            2
            sage: pat.f(3)
              call f(3) -> 4
            4

        .. automethod:: __getattribute__
        .. automethod:: __setattr__
        """
        helper = AttributeAccessTracerHelper(delegate, **kwds)
        object.__setattr__(self, "helper", helper)

    def __getattribute__(self, name):
        r"""
        Read an attribute from the wrapped delegate object.

        If that value is a method (i.e. a callable object which is not
        contained in the dictionary of the object itself but instead
        inherited from some class) then it is replaced by a wrapper
        function to report arguments and return value.
        Otherwise an attribute read access is reported.

        EXAMPLES::

            sage: class Foo():
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: foo.x = 2
            sage: from sage.doctest.fixtures import AttributeAccessTracerProxy
            sage: pat = AttributeAccessTracerProxy(foo)
            sage: pat.x
              read x = 2
            2
            sage: pat.f(3)
              call f(3) -> 4
            4
        """
        helper = object.__getattribute__(self, "helper")
        return helper.get(name)

    def __setattr__(self, name, val):
        r"""
        Write an attribute to the wrapped delegate object.

        The name and new value are also reported in the output.

        EXAMPLES::

            sage: class Foo():
            ....:     pass
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import AttributeAccessTracerProxy
            sage: pat = AttributeAccessTracerProxy(foo)
            sage: pat.x = 2
              write x = 2
            sage: foo.x
            2
        """
        helper = object.__getattribute__(self, "helper")
        return helper.set(name, val)


def trace_method(obj, meth, **kwds):
    r"""
    Trace the doings of a given method.
    It prints method entry with arguments,
    access to members and other methods during method execution
    as well as method exit with return value.

    INPUT:

    - ``obj`` -- the object containing the method

    - ``meth`` -- the name of the method to be traced

    - ``prefix`` -- (default: ``"  "``)
      string to prepend to each printed output

    - ``reads`` -- (default: ``True``)
      whether to trace read access as well

    EXAMPLES::

        sage: class Foo():
        ....:     def f(self, arg=None):
        ....:         self.y = self.g(self.x)
        ....:         if arg: return arg*arg
        ....:     def g(self, arg):
        ....:         return arg + 1
        ....:
        sage: foo = Foo()
        sage: foo.x = 3
        sage: from sage.doctest.fixtures import trace_method
        sage: trace_method(foo, "f")
        sage: foo.f()
        enter f()
          read x = 3
          call g(3) -> 4
          write y = 4
        exit f -> None
        sage: foo.f(3)
        enter f(3)
          read x = 3
          call g(3) -> 4
          write y = 4
        exit f -> 9
        9
    """
    from sage.cpython.getattr import raw_getattr
    f = raw_getattr(obj, meth)
    t = AttributeAccessTracerProxy(obj, **kwds)

    @wraps(f)
    def g(*args, **kwds):
        from IPython.lib.pretty import pretty
        arglst = [pretty(arg) for arg in args]
        arglst.extend("{}={}".format(k, pretty(v))
                      for k, v in sorted(kwds.items()))
        print("enter {}({})".format(meth, ", ".join(arglst)))
        res = f(t, *args, **kwds)
        print("exit {} -> {}".format(meth, pretty(res)))
        return res
    setattr(obj, meth, g)
