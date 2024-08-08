# sage_setup: distribution = sagemath-objects
"""
Lazy attributes

AUTHORS:

- Nicolas Thiery (2008): Initial version
- Nils Bruin (2013-05): Cython version
"""

# ****************************************************************************
#       Copyright (C) 2008 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef class _lazy_attribute():
    """
    Cython base class for lazy attributes.

    EXAMPLES:

    Only Python subclasses of this class are supposed to be instantiated::

        sage: from sage.misc.lazy_attribute import _lazy_attribute
        sage: _lazy_attribute(lambda x:1)
        Traceback (most recent call last):
        ...
        NotImplementedError: Only instantiate wrapper python class
    """

    cdef public f
    cdef public str __name__

    def __init__(self, f):
        r"""
        Constructor for lazy attributes.

        EXAMPLES::

            sage: def f(x):
            ....:     "doc of f"
            ....:     return 1
            ....:
            sage: x = lazy_attribute(f); x
            <sage.misc.lazy_attribute.lazy_attribute object at ...>
            sage: x.__doc__
            'doc of f'
            sage: x.__name__
            'f'
            sage: x.__module__
            '__main__'

        TESTS:

        We check that :issue:`9251` is solved::

            sage: Parent.element_class
            <sage.misc.lazy_attribute.lazy_attribute object at 0x...>
            sage: Parent.element_class.__doc__[91:147]
            'The (default) class for the elements of this parent.\n\n   '
            sage: Parent.element_class.__name__
            'element_class'
            sage: Parent.element_class.__module__
            'sage.misc.lazy_attribute'
        """
        raise NotImplementedError("Only instantiate wrapper python class")

    def _sage_src_lines_(self):
        r"""
        Return the source code location for the wrapped function.

        EXAMPLES::

            sage: from sage.misc.sageinspect import sage_getsourcelines
            sage: g = lazy_attribute(sage.misc.banner.banner)
            sage: (src, lines) = sage_getsourcelines(g)
            sage: src[0]
            'def banner():\n'
            sage: lines
            87
        """
        from sage.misc.sageinspect import sage_getsourcelines
        return sage_getsourcelines(self.f)

    def __get__(self, a, cls):
        """
        Implement the attribute access protocol.

        EXAMPLES::

            sage: class A: pass
            sage: def f(x): return 1
            ...
            sage: f = lazy_attribute(f)
            sage: f.__get__(A(), A)
            1
        """
        cdef CM
        cdef result
        if a is None: # when doing cls.x for cls a class and x a lazy attribute
            return self
        try:
            # _cached_methods is supposed to be a public Cython attribute.
            # Apparently, these are *not* subject to name mangling.
            CM = getattr(a, '_cached_methods')
            if CM is None:
                CM = {}
                setattr(a, '_cached_methods', CM)
        except AttributeError as msg:
            CM = None
        if CM is not None:
            try:
                return CM[self.__name__]
            except KeyError:
                pass
        result = self.f(a)
        if result is NotImplemented:
            # Workaround: we make sure that cls is the class
            # where the lazy attribute self is actually defined.
            # This avoids running into an infinite loop
            # See About descriptor specifications
            for supercls in cls.__mro__:
                if self.__name__ in supercls.__dict__ and self is supercls.__dict__[self.__name__]:
                    cls = supercls
            return getattr(super(cls, a),self.__name__)
        try:
            setattr(a, self.__name__, result)
        except AttributeError:
            if CM is not None:
                CM[self.__name__] = result
                return result
            raise
        return result


class lazy_attribute(_lazy_attribute):
    r"""
    A lazy attribute for an object is like a usual attribute, except
    that, instead of being computed when the object is constructed
    (i.e. in ``__init__``), it is computed on the fly the first time it
    is accessed.

    For constant values attached to an object, lazy attributes provide
    a shorter syntax and automatic caching (unlike methods), while
    playing well with inheritance (like methods): a subclass can
    easily override a given attribute; you don't need to call the
    super class constructor, etc.

    Technically, a :class:`lazy_attribute` is a non-data descriptor (see
    Invoking Descriptors in the Python reference manual).

    EXAMPLES:

    We create a class whose instances have a lazy attribute ``x``::

        sage: class A():
        ....:     def __init__(self):
        ....:         self.a=2 # just to have some data to calculate from
        ....:
        ....:     @lazy_attribute
        ....:     def x(self):
        ....:         print("calculating x in A")
        ....:         return self.a + 1
        ....:

    For an instance ``a`` of ``A``, ``a.x`` is calculated the first time it
    is accessed, and then stored as a usual attribute::

        sage: a = A()
        sage: a.x
        calculating x in A
        3
        sage: a.x
        3

    .. rubric:: Implementation details

    We redo the same example, but opening the hood to see what happens to
    the internal dictionary of the object::

        sage: a = A()
        sage: a.__dict__
        {'a': 2}
        sage: a.x
        calculating x in A
        3
        sage: a.__dict__
        {'a': 2, 'x': 3}
        sage: a.x
        3
        sage: timeit('a.x') # random
        625 loops, best of 3: 89.6 ns per loop

    This shows that, after the first calculation, the attribute ``x``
    becomes a usual attribute; in particular, there is no time penalty
    to access it.

    A lazy attribute may be set as usual, even before its first access,
    in which case the lazy calculation is completely ignored::

        sage: a = A()
        sage: a.x = 4
        sage: a.x
        4
        sage: a.__dict__
        {'a': 2, 'x': 4}

    Class binding results in the lazy attribute itself::

        sage: A.x
        <sage.misc.lazy_attribute.lazy_attribute object at ...>

    .. rubric:: Conditional definitions

    The function calculating the attribute may return NotImplemented
    to declare that, after all, it is not able to do it. In that case,
    the attribute lookup proceeds in the super class hierarchy::

        sage: class B(A):
        ....:     @lazy_attribute
        ....:     def x(self):
        ....:         if hasattr(self, "y"):
        ....:             print("calculating x from y in B")
        ....:             return self.y
        ....:         else:
        ....:             print("y not there; B does not define x")
        ....:             return NotImplemented
        ....:
        sage: b = B()
        sage: b.x
        y not there; B does not define x
        calculating x in A
        3
        sage: b = B()
        sage: b.y = 1
        sage: b.x
        calculating x from y in B
        1

    .. rubric:: Attribute existence testing

    Testing for the existence of an attribute with hasattr currently
    always triggers its full calculation, which may not be desirable
    when the calculation is expensive::

        sage: a = A()
        sage: hasattr(a, "x")
        calculating x in A
        True

    It would be great if we could take over the control somehow, if at
    all possible without a special implementation of hasattr, so as to
    allow for something like::

        sage: class A ():
        ....:     @lazy_attribute
        ....:     def x(self, existence_only=False):
        ....:         if existence_only:
        ....:             print("testing for x existence")
        ....:             return True
        ....:         else:
        ....:             print("calculating x in A")
        ....:             return 3
        ....:
        sage: a = A()
        sage: hasattr(a, "x") # todo: not implemented
        testing for x existence
        sage: a.x
        calculating x in A
        3
        sage: a.x
        3

    Here is a full featured example, with both conditional definition
    and existence testing::

        sage: class B(A):
        ....:     @lazy_attribute
        ....:     def x(self, existence_only=False):
        ....:         if hasattr(self, "y"):
        ....:             if existence_only:
        ....:                 print("testing for x existence in B")
        ....:                 return True
        ....:             else:
        ....:                 print("calculating x from y in B")
        ....:                 return self.y
        ....:         else:
        ....:             print("y not there; B does not define x")
        ....:             return NotImplemented
        ....:
        sage: b = B()
        sage: hasattr(b, "x") # todo: not implemented
        y not there; B does not define x
        testing for x existence
        True
        sage: b.x
        y not there; B does not define x
        calculating x in A
        3
        sage: b = B()
        sage: b.y = 1
        sage: hasattr(b, "x") # todo: not implemented
        testing for x existence in B
        True
        sage: b.x
        calculating x from y in B
        1


    .. rubric:: lazy attributes and introspection

    .. TODO::

        Make the following work nicely::

            sage: b.x?                # todo: not implemented
            sage: b.x??               # todo: not implemented

    Right now, the first one includes the doc of this class, and the
    second one brings up the code of this class, both being not very
    useful.

    TESTS:

    .. rubric:: Lazy attributes and Cython

    This attempts to check that lazy attributes work with built-in
    functions like cpdef methods::

        sage: class A:
        ....:     def __len__(x):
        ....:         return int(5)
        ....:     len = lazy_attribute(len)
        ....:
        sage: A().len
        5

    Since :issue:`11115`, extension classes derived from
    :class:`~sage.structure.parent.Parent` can inherit a lazy attribute,
    such as ``element_class``::

        sage: cython_code = ["from sage.structure.parent cimport Parent",
        ....: "from sage.structure.element cimport Element",
        ....: "cdef class MyElement(Element): pass",
        ....: "cdef class MyParent(Parent):",
        ....: "    Element = MyElement"]
        sage: cython('\n'.join(cython_code))                                            # needs sage.misc.cython
        sage: P = MyParent(category=Rings())                                            # needs sage.misc.cython
        sage: P.element_class    # indirect doctest                                     # needs sage.misc.cython
        <class '...MyElement'>

    .. rubric:: About descriptor specifications

    The specifications of descriptors (see 3.4.2.3 Invoking
    Descriptors in the Python reference manual) are incomplete
    w.r.t. inheritance, and maybe even ill-implemented. We illustrate
    this on a simple class hierarchy, with an instrumented descriptor::

        sage: class descriptor():
        ....:     def __get__(self, obj, cls):
        ....:         print(cls)
        ....:         return 1
        sage: class A():
        ....:     x = descriptor()
        sage: class B(A):
        ....:     pass
        ....:

    This is fine::

        sage: A.x
        <class '__main__.A'>
        1

    The behaviour for the following case is not specified (see Instance Binding)
    when ``x`` is not in the dictionary of ``B`` but in that of some super
    category::

        sage: B().x
        <class '__main__.B'>
        1

    It would seem more natural (and practical!) to get ``A`` rather than ``B``.

    From the specifications for Super Binding, it would be expected to
    get ``A`` and not ``B`` as cls parameter::

        sage: super(B, B()).x
        <class '__main__.B'>
        1

    Due to this, the natural implementation runs into an infinite loop
    in the following example::

        sage: class A():
        ....:     @lazy_attribute
        ....:     def unimplemented_A(self):
        ....:         return NotImplemented
        ....:     @lazy_attribute
        ....:     def unimplemented_AB(self):
        ....:         return NotImplemented
        ....:     @lazy_attribute
        ....:     def unimplemented_B_implemented_A(self):
        ....:         return 1
        ....:
        sage: class B(A):
        ....:     @lazy_attribute
        ....:     def unimplemented_B(self):
        ....:         return NotImplemented
        ....:     @lazy_attribute
        ....:     def unimplemented_AB(self):
        ....:         return NotImplemented
        ....:     @lazy_attribute
        ....:     def unimplemented_B_implemented_A(self):
        ....:         return NotImplemented
        ....:
        sage: class C(B):
        ....:     pass
        ....:

    This is the simplest case where, without workaround, we get an
    infinite loop::

        sage: hasattr(B(), "unimplemented_A") # todo: not implemented
        False

    .. TODO::

        Improve the error message::

            sage: B().unimplemented_A # todo: not implemented
            Traceback (most recent call last):
            ...
            AttributeError: 'super' object has no attribute 'unimplemented_A'...

    We now make some systematic checks::

        sage: B().unimplemented_A
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'unimplemented_A'...
        sage: B().unimplemented_B
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'unimplemented_B'...
        sage: B().unimplemented_AB
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'unimplemented_AB'...
        sage: B().unimplemented_B_implemented_A
        1

        sage: C().unimplemented_A()
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'unimplemented_A'...
        sage: C().unimplemented_B()
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'unimplemented_B'...
        sage: C().unimplemented_AB()
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'unimplemented_AB'...
        sage: C().unimplemented_B_implemented_A # todo: not implemented
        1
    """
    def __init__(self, f):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: def f(x):
            ....:     "doc of f"
            ....:     return 1
            ....:
            sage: x = lazy_attribute(f); x
            <sage.misc.lazy_attribute.lazy_attribute object at ...>
            sage: x.__doc__
            'doc of f'
            sage: x.__name__
            'f'
            sage: x.__module__
            '__main__'
        """
        self.f = f
        if hasattr(f, "__doc__"):
            self.__doc__ = f.__doc__
        elif hasattr(f, "__doc__"): # Needed to handle Cython methods
            self.__doc__ = f.__doc__
        if hasattr(f, "__name__"):
            self.__name__ = f.__name__
        elif hasattr(f, "__name__"): # Needed to handle Cython methods
            self.__name__ = f.__name__
        if hasattr(f, "__module__"):
            self.__module__ = f.__module__


class lazy_class_attribute(lazy_attribute):
    """
    A lazy class attribute for a class is like a usual class attribute,
    except that, instead of being computed when the class is constructed, it
    is computed on the fly the first time it is accessed, either through the
    class itself or trough on of its objects.

    This is very similar to :class:`lazy_attribute` except that the attribute
    is a class attribute. More precisely, once computed, the lazy class
    attribute is stored in the class rather than in the object. The lazy class
    attribute is only computed once for all the objects::

        sage: class Cl():
        ....:     @lazy_class_attribute
        ....:     def x(cls):
        ....:          print("computing x")
        ....:          return 1
        sage: Cl.x
        computing x
        1
        sage: Cl.x
        1

    As for a any usual class attribute it is also possible to access it from
    an object::

        sage: b = Cl()
        sage: b.x
        1

    First access from an object also properly triggers the computation::

        sage: class Cl1():
        ....:     @lazy_class_attribute
        ....:     def x(cls):
        ....:          print("computing x")
        ....:          return 1
        sage: Cl1().x
        computing x
        1
        sage: Cl1().x
        1

    .. WARNING::

        The behavior of lazy class attributes with respect to inheritance is
        not specified. It currently depends on the evaluation order::

            sage: class A():
            ....:     @lazy_class_attribute
            ....:     def x(cls):
            ....:          print("computing x")
            ....:          return str(cls)
            ....:     @lazy_class_attribute
            ....:     def y(cls):
            ....:          print("computing y")
            ....:          return str(cls)
            sage: class B(A):
            ....:     pass

            sage: A.x
            computing x
            "<class '__main__.A'>"
            sage: B.x
            "<class '__main__.A'>"

            sage: B.y
            computing y
            "<class '__main__.B'>"
            sage: A.y
            computing y
            "<class '__main__.A'>"
            sage: B.y
            "<class '__main__.B'>"

    TESTS::

        sage: "x" in b.__dict__
        False
    """
    def __get__(self, _, cls):
        """
        Implement the attribute access protocol.

        EXAMPLES::

            sage: class A: pass
            sage: def f(x): return 1
            ....:
            sage: f = lazy_class_attribute(f)
            sage: f.__get__(A(), A)
            1
        """
        result = self.f(cls)
        if result is NotImplemented:
            return getattr(super(cls, cls), self.__name__)
        setattr(cls, self.__name__, result)
        return result
