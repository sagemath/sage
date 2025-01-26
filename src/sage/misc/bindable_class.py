# sage_setup: distribution = sagemath-objects
"""
Bindable classes
"""
# ****************************************************************************
#       Copyright (C) 2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import functools
from sage.misc.nested_class import NestedClassMetaclass
from sage.misc.classcall_metaclass import ClasscallMetaclass


class BindableClass(metaclass=ClasscallMetaclass):
    """
    Bindable classes.

    This class implements a binding behavior for nested classes that
    derive from it. Namely, if a nested class ``Outer.Inner`` derives
    from ``BindableClass``, and if ``outer`` is an instance of
    ``Outer``, then ``outer.Inner(...)`` is equivalent to
    ``Outer.Inner(outer, ...)``.

    EXAMPLES:

    Let us consider the following class ``Outer`` with a nested class ``Inner``::

        sage: from sage.misc.nested_class import NestedClassMetaclass
        sage: class Outer(metaclass=NestedClassMetaclass):
        ....:     class Inner:
        ....:         def __init__(self, *args):
        ....:             print(args)
        ....:     def f(self, *args):
        ....:         print("{} {}".format(self, args))
        ....:     @staticmethod
        ....:     def f_static(*args):
        ....:         print(args)

        sage: outer = Outer()

    By default, when ``Inner`` is a class nested in ``Outer``,
    accessing ``outer.Inner`` returns the ``Inner`` class as is::

        sage: outer.Inner is Outer.Inner
        True

    In particular, ``outer`` is completely ignored in the following call::

        sage: x = outer.Inner(1,2,3)
        (1, 2, 3)

    This is similar to what happens with a static method::

        sage: outer.f_static(1,2,3)
        (1, 2, 3)

    In some cases, we would want instead ``Inner`` to receive ``outer``
    as parameter, like in a usual method call::

        sage: outer.f(1,2,3)
        <__main__.Outer object at ...> (1, 2, 3)

    To this end, ``outer.f`` returns a *bound method*::

        sage: outer.f
        <bound method Outer.f of <__main__.Outer object at ...>>

    so that ``outer.f(1,2,3)`` is equivalent to::

        sage: Outer.f(outer, 1,2,3)
        <__main__.Outer object at ...> (1, 2, 3)

    :class:`BindableClass` gives this binding behavior to all its subclasses::

        sage: from sage.misc.bindable_class import BindableClass
        sage: class Outer(metaclass=NestedClassMetaclass):
        ....:     class Inner(BindableClass):
        ....:         " some documentation "
        ....:         def __init__(self, outer, *args):
        ....:             print("{} {}".format(outer, args))

    Calling ``Outer.Inner`` returns the (unbound) class as usual::

        sage: Outer.Inner
        <class '__main__.Outer.Inner'>

    However, ``outer.Inner(1,2,3)`` is equivalent to ``Outer.Inner(outer, 1,2,3)``::

        sage: outer = Outer()
        sage: x = outer.Inner(1,2,3)
        <__main__.Outer object at ...> (1, 2, 3)

    To achieve this, ``outer.Inner`` returns (some sort of) bound class::

        sage: outer.Inner
        <bound class '__main__.Outer.Inner' of <__main__.Outer object at ...>>

    .. NOTE::

        This is not actually a class, but an instance of
        :class:`functools.partial`::

            sage: type(outer.Inner).mro()
            [<class 'sage.misc.bindable_class.BoundClass'>,
             <... 'functools.partial'>,
             <... 'object'>]

        Still, documentation works as usual::

            sage: outer.Inner.__doc__
            '...some documentation '

    TESTS::

        sage: from sage.misc.bindable_class import Outer
        sage: TestSuite(Outer.Inner).run()
        sage: outer = Outer()
        sage: TestSuite(outer.Inner).run(skip=["_test_pickling"])
    """
    @staticmethod
    def __classget__(cls, instance, owner):
        """
        Bind ``cls`` to ``instance``, returning a ``BoundClass``.

        INPUT:

        - ``instance`` -- an object of the outer class or ``None``

        For technical details, see the Section :python:`Implementing Descriptor
        <reference/datamodel.html#implementing-descriptors>` in the Python
        reference manual.

        EXAMPLES::

            sage: from sage.misc.bindable_class import Outer
            sage: Outer.Inner
            <class 'sage.misc.bindable_class.Outer.Inner'>
            sage: Outer().Inner
            <bound class 'sage.misc.bindable_class.Outer.Inner' of <sage.misc.bindable_class.Outer object at ...>>
        """
        if instance is None:
            return cls
        return BoundClass(cls, instance)
        # We probably do not need to use sage_wraps, since
        # sageinspect already supports partial functions
        # return sage_wraps(cls)(BoundClass(cls, instance))


class BoundClass(functools.partial):
    """
    TESTS::

        sage: from sage.misc.sageinspect import *
        sage: from sage.misc.bindable_class import Outer
        sage: x = Outer()
        sage: c = x.Inner; c
        <bound class 'sage.misc.bindable_class.Outer.Inner' of <sage.misc.bindable_class.Outer object at ...>>

    Introspection works, at least partially::

        sage: sage_getdoc(c).strip()
        'Some documentation for Outer.Inner'
        sage: sage_getfile(c)
        '.../sage/misc/bindable_class.py'

        sage: c = x.Inner2
        sage: sage_getdoc(c).strip()
        'Some documentation for Inner2'
        sage: sage_getsourcelines(c)
        (['class Inner2(BindableClass):...], ...)

    .. warning::

        Since ``c`` is not a class (as tested by :func:`inspect.isclass`),
        and has a ``__call__`` method, IPython's introspection
        (with ``c?``) insists on showing not only its
        documentation but also its class documentation and call
        documentation (see :meth:`IPython.OInspect.Inspector.pdoc`)
        if available.

        Until a better approach is found, we reset the documentation
        of :class:`BoundClass` below, and make an exception for
        :meth:`__init__` to the strict rule that every method should
        be doctested::

            sage: c.__class__.__doc__
            sage: c.__class__.__init__.__doc__

    Make sure classes which inherit from :class:`functools.partial` have the correct
    syntax, see :issue:`14748`::

        sage: import warnings
        sage: warnings.simplefilter('error', DeprecationWarning)
        sage: import functools
        sage: def f(x, y): return x^y
        sage: g = functools.partial(f, 2, 3)
        sage: g()
        8

    The following has correct syntax and no :class:`DeprecationWarning`::

        sage: class mynewpartial(functools.partial):
        ....:     def __init__(self, f, i, j):
        ....:         functools.partial.__init__(self)
        sage: g = mynewpartial(f, 2, 3)
        sage: g()
        8
    """
    __doc__ = None  # See warning above

    def __init__(self, *args):
        super().__init__()
        self.__doc__ = self.func.__doc__

    def __repr__(self):
        """
        TESTS:

            sage: from sage.misc.bindable_class import Outer
            sage: x = Outer(); x
            <sage.misc.bindable_class.Outer object at ...>
            sage: x.Inner
            <bound class 'sage.misc.bindable_class.Outer.Inner' of <sage.misc.bindable_class.Outer object at ...>>
        """
        return "<bound %s of %s>" % (repr(self.func)[1:-1], self.args[0])


##############################################################################
# Test classes
##############################################################################

class Inner2(BindableClass):
    """
    Some documentation for Inner2
    """


# We need NestedClassMetaclass to work around a Python pickling bug
class Outer(metaclass=NestedClassMetaclass):
    """
    A class with a bindable nested class, for testing purposes
    """
    class Inner(BindableClass):
        """
        Some documentation for Outer.Inner
        """

    Inner2 = Inner2
