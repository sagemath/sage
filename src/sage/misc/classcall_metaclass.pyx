# sage_setup: distribution = sagemath-objects
r"""
Special Methods for Classes

AUTHORS:

- Nicolas M. Thiery (2009-2011) implementation of
  ``__classcall__``, ``__classget__``, ``__classcontains__``;
- Florent Hivert (2010-2012): implementation of ``__classcall_private__``,
  documentation, Cythonization and optimization.
"""

# ****************************************************************************
#       Copyright (C) 2009      Nicolas M. Thiery <nthiery at users.sf.net>
#       Copyright (C) 2010-2012 Florent Hivert <Florent.Hivert at lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport *
from cpython.type cimport type as pytype

__all__ = ['ClasscallMetaclass', 'typecall', 'timeCall']


cdef class ClasscallMetaclass(NestedClassMetaclass):
    r"""
    A metaclass providing support for special methods for classes.

    From the Section :python:`Special method names
    <reference/datamodel.html#special-method-names>` of the Python Reference
    Manual:

        \`a class ``cls`` can implement certain operations on its instances
        that are invoked by special syntax (such as arithmetic operations or
        subscripting and slicing) by defining methods with special
        names\'.

    The purpose of this metaclass is to allow for the class ``cls`` to
    implement analogues of those special methods for the operations on the
    class itself.

    Currently, the following special methods are supported:

     - ``.__classcall__`` (and ``.__classcall_private__``) for
       customizing ``cls(...)`` (analogue of ``.__call__``).

     - ``.__classcontains__`` for customizing membership testing
       ``x in cls`` (analogue of ``.__contains__``).

     - ``.__classget__`` for customizing the binding behavior in
       ``foo.cls`` (analogue of ``.__get__``).

    See the documentation of :meth:`__call__` and of :meth:`__get__`
    and :meth:`__contains__` for the description of the respective
    protocols.

    .. WARNING::

        For technical reasons, ``__classcall__``,
        ``__classcall_private__``, ``__classcontains__``, and
        ``__classget__`` must be defined as :func:`staticmethod`'s,
        even though they receive the class itself as their first
        argument.

    .. WARNING::

        For efficiency reasons, the resolution for the special methods
        is done once for all, upon creation of the class. Thus, later
        dynamic changes to those methods are ignored. But see also
        :meth:`_set_classcall`.

    ``ClasscallMetaclass`` is an extension of the base :class:`type`.

    .. TODO:: find a good name for this metaclass.

    TESTS::

        sage: PerfectMatchings(2).list()                                                # needs sage.combinat
        [[(1, 2)]]

    .. NOTE::

        If a class is put in this metaclass it automatically becomes a
        new-style class::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Foo(metaclass=ClasscallMetaclass): pass
            sage: x = Foo(); x
            <__main__.Foo object at 0x...>
            sage: issubclass(Foo, object)
            True
            sage: isinstance(Foo, type)
            True
    """

    def __cinit__(self, *args, **opts):
        r"""
        TESTS::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class FOO(metaclass=ClasscallMetaclass): pass
            sage: isinstance(FOO, ClasscallMetaclass)  # indirect doctest
            True
        """
        if '__classcall_private__' in self.__dict__:
            self.classcall = self.__classcall_private__
        elif hasattr(self, "__classcall__"):
            self.classcall = self.__classcall__
        else:
            self.classcall = None

        self.classcontains = getattr(self, "__classcontains__", None)
        self.classget = getattr(self, "__classget__", None)

    def _set_classcall(cls, function):
        r"""
        Change dynamically the classcall function for this class.

        EXAMPLES::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class FOO(metaclass=ClasscallMetaclass): pass
            sage: FOO()
            <__main__.FOO object at ...>

        For efficiency reason, the resolution of the ``__classcall__``
        method is done once for all, upon creation of the class. Thus,
        later dynamic changes to this method are ignored by FOO::

            sage: FOO.__classcall__ = ConstantFunction(1)
            sage: FOO()
            <__main__.FOO object at ...>

        but not by subclasses created later on::

            sage: class BAR(FOO): pass
            sage: BAR()
            1

        To update the ``classcall`` special function for FOO, one
        should use this setter::

            sage: FOO._set_classcall(ConstantFunction(2))
            sage: FOO()
            2

        Note that it has no influence on subclasses::

            sage: class BAR(FOO): pass
            sage: BAR()
            1
        """
        cls.classcall = function

    def __call__(cls, *args, **kwds):
        r"""
        This method implements ``cls(<some arguments>)``.

        Let ``cls`` be a class in :class:`ClasscallMetaclass`, and
        consider a call of the form::

            cls(<some arguments>)

        - If ``cls`` defines a method ``__classcall_private__``, then
          this results in a call to::

            cls.__classcall_private__(cls, <some arguments>)

        - Otherwise, if ``cls`` has a method ``__classcall__``, then instead
          the following is called::

            cls.__classcall__(cls, <some arguments>)

        - If neither of these two methods are implemented, then the standard
          ``type.__call__(cls, <some arguments>)`` is called, which in turn
          uses :meth:`~object.__new__` and :meth:`~object.__init__` as usual
          (see Section :python:`Basic Customization
          <reference/datamodel.html#basic-customization>` in the Python
          Reference Manual).

        .. warning:: for technical reasons, ``__classcall__`` must be
            defined as a :func:`staticmethod`, even though it receives
            the class itself as its first argument.

        EXAMPLES::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Foo(metaclass=ClasscallMetaclass):
            ....:     @staticmethod
            ....:     def __classcall__(cls):
            ....:         print("calling classcall")
            ....:         return type.__call__(cls)
            ....:     def __new__(cls):
            ....:         print("calling new")
            ....:         return super(Foo, cls).__new__(cls)
            ....:     def __init__(self):
            ....:         print("calling init")
            sage: Foo()
            calling classcall
            calling new
            calling init
            <__main__.Foo object at ...>

        This behavior is inherited::

            sage: class Bar(Foo): pass
            sage: Bar()
            calling classcall
            calling new
            calling init
            <__main__.Bar object at ...>

        We now show the usage of ``__classcall_private__``::

            sage: class FooNoInherits(object, metaclass=ClasscallMetaclass):
            ....:     @staticmethod
            ....:     def __classcall_private__(cls):
            ....:         print("calling private classcall")
            ....:         return type.__call__(cls)
            sage: FooNoInherits()
            calling private classcall
            <__main__.FooNoInherits object at ...>

        Here the behavior is not inherited::

            sage: class BarNoInherits(FooNoInherits): pass
            sage: BarNoInherits()
            <__main__.BarNoInherits object at ...>

        We now show the usage of both::

            sage: class Foo2(object, metaclass=ClasscallMetaclass):
            ....:     @staticmethod
            ....:     def __classcall_private__(cls):
            ....:         print("calling private classcall")
            ....:         return type.__call__(cls)
            ....:     @staticmethod
            ....:     def __classcall__(cls):
            ....:         print("calling classcall with %s" % cls)
            ....:         return type.__call__(cls)
            ...
            sage: Foo2()
            calling private classcall
            <__main__.Foo2 object at ...>

            sage: class Bar2(Foo2): pass
            sage: Bar2()
            calling classcall with <class '__main__.Bar2'>
            <__main__.Bar2 object at ...>


        .. rubric:: Discussion

        Typical applications include the implementation of factories or of
        unique representation (see :class:`UniqueRepresentation`). Such
        features are traditionally implemented by either using a wrapper
        function, or fiddling with :meth:`~object.__new__`.

        The benefit, compared with fiddling directly with
        :meth:`~object.__new__` is a clear separation of the three distinct
        roles:

        - ``cls.__classcall__``: what ``cls(<...>)`` does
        - ``cls.__new__``: memory allocation for a *new* instance
        - ``cls.__init__``: initialization of a newly created instance

        The benefit, compared with using a wrapper function, is that the
        user interface has a single handle for the class::

            sage: x = Partition([3,2,2])                                                # needs sage.combinat
            sage: isinstance(x, Partition)      # not implemented                       # needs sage.combinat

        instead of::

            sage: isinstance(x, sage.combinat.partition.Partition)                      # needs sage.combinat
            True

        Another difference is that ``__classcall__`` is inherited by
        subclasses, which may be desirable, or not. If not, one should
        instead define the method ``__classcall_private__`` which will
        not be called for subclasses. Specifically, if a class ``cls``
        defines both methods ``__classcall__`` and
        ``__classcall_private__`` then, for any subclass ``sub`` of ``cls``:

        - ``cls(<args>)`` will call ``cls.__classcall_private__(cls, <args>)``
        - ``sub(<args>)`` will call ``cls.__classcall__(sub, <args>)``


        TESTS:

        We check that the idiom ``method_name in cls.__dict__`` works
        for extension types::

           sage: "_sage_" in SageObject.__dict__, "_sage_" in Parent.__dict__
           (True, False)

        We check for memory leaks::

            sage: class NOCALL(object, metaclass=ClasscallMetaclass):
            ....:    pass
            sage: sys.getrefcount(NOCALL())
            1

        We check that exceptions are correctly handled::

            sage: class Exc(object, metaclass=ClasscallMetaclass):
            ....:     @staticmethod
            ....:     def __classcall__(cls):
            ....:         raise ValueError("Calling classcall")
            sage: Exc()
            Traceback (most recent call last):
            ...
            ValueError: Calling classcall
        """
        if cls.classcall is not None:
            return cls.classcall(cls, *args, **kwds)
        else:
            # Fast version of type.__call__(cls, *args, **kwds)
            return (<PyTypeObject*>type).tp_call(cls, args, kwds)

    def __get__(cls, instance, owner):
        r"""
        This method implements instance binding behavior for nested classes.

        Suppose that a class ``Outer`` contains a nested class ``cls`` which
        is an instance of this metaclass. For any object ``obj`` of ``cls``,
        this method implements a instance binding behavior for ``obj.cls`` by
        delegating it to ``cls.__classget__(Outer, obj, owner)`` if available.
        Otherwise, ``obj.cls`` results in ``cls``, as usual.

        Similarly, a class binding as in ``Outer.cls`` is delegated
        to ``cls.__classget__(Outer, None, owner)`` if available and
        to ``cls`` if not.

        .. warning:: for technical reasons, ``__classget__`` must be
            defined as a :func:`staticmethod`, even though it receives
            the class itself as its first argument.

        For technical details, and in particular the description of the
        ``owner`` argument, see the Section :python:`Implementing Descriptor
        <reference/datamodel.html#implementing-descriptors>` in the Python
        reference manual.

        EXAMPLES:

        We show how to implement a nested class ``Outer.Inner`` with a
        binding behavior, as if it was a method of ``Outer``: namely,
        for ``obj`` an instance of ``Outer``, calling
        ``obj.Inner(...)`` is equivalent to ``Outer.Inner(obj, ...)``::

            sage: import functools
            sage: from sage.misc.nested_class import NestedClassMetaclass
            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Outer(metaclass=NestedClassMetaclass):
            ....:     class Inner(metaclass=ClasscallMetaclass):
            ....:         @staticmethod
            ....:         def __classget__(cls, instance, owner):
            ....:             print("calling __classget__(%s, %s, %s)" % (
            ....:                        cls, instance, owner))
            ....:             if instance is None:
            ....:                 return cls
            ....:             return functools.partial(cls, instance)
            ....:         def __init__(self, instance):
            ....:             self.instance = instance
            sage: obj = Outer()
            sage: bar = obj.Inner()
            calling __classget__(<class '__main__.Outer.Inner'>, <__main__.Outer object at 0x...>, <class '__main__.Outer'>)
            sage: bar.instance == obj
            True

        Calling ``Outer.Inner`` returns the (unbinded) class as usual::

            sage: Inner = Outer.Inner
            calling __classget__(<class '__main__.Outer.Inner'>, None, <class '__main__.Outer'>)
            sage: Inner
            <class '__main__.Outer.Inner'>
            sage: type(bar) is Inner
            True

        .. warning:: Inner has to be a new style class (i.e. a subclass of object).

        .. warning::

            Calling ``obj.Inner`` no longer returns a class::

                sage: bind = obj.Inner
                calling __classget__(<class '__main__.Outer.Inner'>, <__main__.Outer object at 0x...>, <class '__main__.Outer'>)
                sage: bind
                functools.partial(<class '__main__.Outer.Inner'>, <__main__.Outer object at 0x...>)
        """
        if cls.classget:
            return cls.classget(cls, instance, owner)
        else:
            return cls

    def __contains__(cls, x):
        r"""
        This method implements membership testing for a class.

        Let ``cls`` be a class in :class:`ClasscallMetaclass`, and consider
        a call of the form::

            x in cls

        If ``cls`` defines a method ``__classcontains__``, then this
        results in a call to::

           cls.__classcontains__(cls, x)

        .. warning:: for technical reasons, ``__classcontains__`` must
            be defined as a :func:`staticmethod`, even though it
            receives the class itself as its first argument.

        EXAMPLES:

        We construct a class which implements membership testing, and
        which contains ``1`` and no other x::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Foo(metaclass=ClasscallMetaclass):
            ....:     @staticmethod
            ....:     def __classcontains__(cls, x):
            ....:         return x == 1
            sage: 1 in Foo
            True
            sage: 2 in Foo
            False

        We now check that for a class without ``__classcontains__``
        method, we emulate the usual error message::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Bar(metaclass=ClasscallMetaclass): pass
            sage: 1 in Bar
            Traceback (most recent call last):
            ...
            TypeError: argument of type 'type' is not iterable
        """
        if cls.classcontains:
            return cls.classcontains(cls, x)
        else:
            return x in object


def typecall(pytype cls, *args, **kwds):
    r"""
    Object construction.

    This is a faster equivalent to ``type.__call__(cls, <some arguments>)``.

    INPUT:

    - ``cls`` -- the class used for constructing the instance; it must be
      a builtin type or a new style class (inheriting from :class:`object`)

    EXAMPLES::

        sage: from sage.misc.classcall_metaclass import typecall
        sage: class Foo(): pass
        sage: typecall(Foo)
        <__main__.Foo object at 0x...>
        sage: typecall(list)
        []
        sage: typecall(Integer, 2)
        2
    """
    return (<PyTypeObject*>type).tp_call(cls, args, kwds)


# Class for timing::

class CRef():
    def __init__(self, i):
        """
        TESTS::

            sage: from sage.misc.classcall_metaclass import CRef
            sage: P = CRef(2); P.i
            3
        """
        self.i = i+1


class C2(object, metaclass=ClasscallMetaclass):
    def __init__(self, i):
        """
        TESTS::

            sage: from sage.misc.classcall_metaclass import C2
            sage: P = C2(2); P.i
            3
        """
        self.i = i+1


class C3(object, metaclass = ClasscallMetaclass):
    def __init__(self, i):
        """
        TESTS::

            sage: from sage.misc.classcall_metaclass import C3
            sage: P = C3(2); P.i
            3
        """
        self.i = i+1


class C2C(object, metaclass=ClasscallMetaclass):
    @staticmethod
    def __classcall__(cls, i):
        """
        TESTS::

            sage: from sage.misc.classcall_metaclass import C2C
            sage: C2C(2)
            3
        """
        return i+1


def timeCall(T, int n, *args):
    r"""
    We illustrate some timing when using the classcall mechanism.

    EXAMPLES::

        sage: from sage.misc.classcall_metaclass import (
        ....:     ClasscallMetaclass, CRef, C2, C3, C2C, timeCall)
        sage: timeCall(object, 1000)

    For reference let construct basic objects and a basic Python class::

        sage: %timeit timeCall(object, 1000)   # not tested
        625 loops, best of 3: 41.4 µs per loop

        sage: i1 = int(1); i3 = int(3) # don't use Sage's Integer
        sage: class PRef():
        ....:     def __init__(self, i):
        ....:         self.i = i+i1

    For a Python class, compared to the reference class there is a 10%
    overhead in using :class:`ClasscallMetaclass` if there is no classcall
    defined::

        sage: class P(metaclass=ClasscallMetaclass):
        ....:     def __init__(self, i):
        ....:         self.i = i+i1

        sage: %timeit timeCall(PRef, 1000, i3)   # not tested
        625 loops, best of 3: 420 µs per loop
        sage: %timeit timeCall(P, 1000, i3)      # not tested
        625 loops, best of 3: 458 µs per loop

    For a Cython class (not cdef since they doesn't allows metaclasses), the
    overhead is a little larger::

        sage: %timeit timeCall(CRef, 1000, i3)   # not tested
        625 loops, best of 3: 266 µs per loop
        sage: %timeit timeCall(C2, 1000, i3)     # not tested
        625 loops, best of 3: 298 µs per loop

    Let's now compare when there is a classcall defined::

        sage: class PC(object, metaclass=ClasscallMetaclass):
        ....:     @staticmethod
        ....:     def __classcall__(cls, i):
        ....:         return i+i1
        sage: %timeit timeCall(C2C, 1000, i3)   # not tested
        625 loops, best of 3: 148 µs per loop
        sage: %timeit timeCall(PC, 1000, i3)    # not tested
        625 loops, best of 3: 289 µs per loop

    The overhead of the indirection ( ``C(...) ->
    ClasscallMetaclass.__call__(...) -> C.__classcall__(...)``) is
    unfortunately quite large in this case (two method calls instead of
    one). In reasonable usecases, the overhead should be mostly hidden by the
    computations inside the classcall::

        sage: %timeit timeCall(C2C.__classcall__, 1000, C2C, i3)  # not tested
        625 loops, best of 3: 33 µs per loop
        sage: %timeit timeCall(PC.__classcall__, 1000, PC, i3)    # not tested
        625 loops, best of 3: 131 µs per loop

    Finally, there is no significant difference between Cython's V2 and V3
    syntax for metaclass::

        sage: %timeit timeCall(C2, 1000, i3)   # not tested
        625 loops, best of 3: 330 µs per loop
        sage: %timeit timeCall(C3, 1000, i3)   # not tested
        625 loops, best of 3: 328 µs per loop
    """
    cdef int i
    for 0<=i<n:
        T(*args)
