# sage_setup: distribution = sagemath-objects
"""
Fast methods via Cython

This module provides extension classes with useful methods of cython speed,
that python classes can inherit.

.. NOTE::

    This module provides a cython base class :class:`WithEqualityById`
    implementing unique instance behaviour, and a cython base class
    :class:`FastHashable_class`, which has a quite fast hash
    whose value can be freely chosen at initialisation time.

AUTHOR:

- Simon King (2013-02): Original version
- Simon King (2013-10): Add :class:`Singleton`
"""

# ****************************************************************************
#       Copyright (C) 2013 Simon A. King <simon.king at uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.constant_function import ConstantFunction

from cpython.object cimport Py_EQ, Py_NE


cdef class WithEqualityById:
    """
    Provide hash and equality test based on identity.

    .. NOTE::

        This class provides the unique representation behaviour of
        :class:`~sage.structure.unique_representation.UniqueRepresentation`,
        together with :class:`~sage.structure.unique_representation.CachedRepresentation`.

    EXAMPLES:

    Any instance of :class:`~sage.structure.unique_representation.UniqueRepresentation`
    inherits from :class:`WithEqualityById`.
    ::

        sage: class MyParent(Parent):
        ....:   def __init__(self, x):
        ....:       self.x = x
        ....:   def __hash__(self):
        ....:       return hash(self.x)
        sage: class MyUniqueParent(UniqueRepresentation, MyParent): pass
        sage: issubclass(MyUniqueParent, sage.misc.fast_methods.WithEqualityById)
        True

    Inheriting from :class:`WithEqualityById` provides unique representation
    behaviour::

        sage: a = MyUniqueParent(1)
        sage: b = MyUniqueParent(2)
        sage: c = MyUniqueParent(1)
        sage: a is c
        True
        sage: d = MyUniqueParent(-1)
        sage: a == d
        False

    The hash inherited from ``MyParent`` is replaced by a hash that coincides
    with :class:`object`'s hash::

        sage: hash(a) == hash(a.x)
        False
        sage: hash(a) == object.__hash__(a)
        True

    .. WARNING::

        It is possible to inherit from
        :class:`~sage.structure.unique_representation.UniqueRepresentation`
        and then overload equality test in a way that destroys the unique
        representation property. We strongly recommend against it!  You should
        use :class:`~sage.structure.unique_representation.CachedRepresentation`
        instead.

    ::

        sage: class MyNonUniqueParent(MyUniqueParent):
        ....:   def __eq__(self, other):
        ....:       return self.x^2 == other.x^2
        sage: a = MyNonUniqueParent(1)
        sage: d = MyNonUniqueParent(-1)
        sage: a is MyNonUniqueParent(1)
        True
        sage: a == d
        True
        sage: a is d
        False
    """
    def __hash__(self):
        """
        The hash provided by this class coincides with that of ``<class 'object'>``.

        TESTS::

            sage: class MyParent(Parent):
            ....:   def __init__(self, x):
            ....:       self.x = x
            ....:   def __hash__(self):
            ....:       return hash(self.x)
            sage: class MyUniqueParent(UniqueRepresentation, MyParent): pass
            sage: issubclass(MyUniqueParent, sage.misc.fast_methods.WithEqualityById)
            True
            sage: a = MyUniqueParent(1)
            sage: hash(a) == hash(a.x)
            False
            sage: hash(a) == object.__hash__(a)
            True

            sage: from sage.misc.fast_methods import WithEqualityById
            sage: o1 = WithEqualityById()
            sage: o2 = WithEqualityById()
            sage: hash(o1) == hash(o2)
            False
        """
        # This is the default hash function in Python's object.c:
        return hash_by_id(<void *>self)

    def __richcmp__(self, other, int op):
        """
        Equality test provided by this class is by identity.

        TESTS::

            sage: class MyParent(Parent):
            ....:   def __init__(self, x):
            ....:       self.x = x
            ....:   def __hash__(self):
            ....:       return hash(self.x)
            sage: class MyUniqueParent(UniqueRepresentation, MyParent): pass
            sage: issubclass(MyUniqueParent, sage.misc.fast_methods.WithEqualityById)
            True
            sage: a = MyUniqueParent(1)
            sage: b = MyUniqueParent(-1)

        Equality test takes into account identity::

            sage: a == b
            False

        When comparing with an object which is not an instance of
        ``WithEqualityById``, the other object determines the
        comparison::

            sage: class AlwaysEqual:
            ....:     def __eq__(self, other):
            ....:         return True
            sage: AlwaysEqual() == a
            True
            sage: a == AlwaysEqual()
            True

        Check that :issue:`19628` is fixed::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazyQQ = LazyImport('sage.rings.rational_field', 'QQ')
            sage: PolynomialRing(lazyQQ, 'ijk') is PolynomialRing(QQ, 'ijk')
            True
            sage: PolynomialRing(QQ, 'ijkl') is PolynomialRing(lazyQQ, 'ijkl')
            True
        """
        # This only makes sense if "other" is also of type WithEqualityById
        if type(self) is not type(other):
            if not isinstance(other, WithEqualityById):
                return NotImplemented

        if op == Py_EQ:
            return self is other
        elif op == Py_NE:
            return self is not other
        return NotImplemented


cdef class FastHashable_class:
    """
    A class that has a fast hash method, returning a pre-assigned value.

    .. NOTE::

        This is for internal use only. The class has a cdef attribute
        ``_hash``, that needs to be assigned (for example, by calling
        the init method, or by a direct assignment using
        cython). This is slower than using :func:`provide_hash_by_id`,
        but has the advantage that the hash can be prescribed, by
        assigning a cdef attribute ``_hash``.

    TESTS::

        sage: from sage.misc.fast_methods import FastHashable_class
        sage: H = FastHashable_class(123)
        sage: hash(H)
        123
    """
    def __init__(self, h):
        """
        TESTS::

            sage: from sage.misc.fast_methods import FastHashable_class
            sage: H = FastHashable_class(123)
            sage: hash(H)   # indirect doctest
            123
        """
        self._hash = h

    def __hash__(self):
        """
        TESTS::

            sage: from sage.misc.fast_methods import FastHashable_class
            sage: H = FastHashable_class(123)
            sage: hash(H)   # indirect doctest
            123
        """
        return self._hash


class Singleton(WithEqualityById, metaclass=ClasscallMetaclass):
    """
    A base class for singletons.

    A singleton is a class that allows to create not more than a
    single instance. This instance can also belong to a subclass, but
    it is not possible to have several subclasses of a singleton all
    having distinct unique instances.

    In order to create a singleton, just add :class:`Singleton`
    to the list of base classes::

        sage: from sage.misc.fast_methods import Singleton
        sage: class C(Singleton, SageObject):
        ....:     def __init__(self):
        ....:         print("creating singleton")
        sage: c = C()
        creating singleton
        sage: c2 = C()
        sage: c is c2
        True

    The unique instance of a singleton stays in memory as long as the
    singleton itself does.

    Pickling, copying, hashing, and comparison are provided for by
    :class:`Singleton` according to the singleton paradigm. Note
    that pickling fails if the class is replaced by a sub-sub-class
    after creation of the instance::

        sage: class D(C):
        ....:     pass
        sage: import __main__      # This is only needed ...
        sage: __main__.C = C       # ... in doctests
        sage: __main__.D = D       # same here, only in doctests
        sage: orig = type(c)
        sage: c.__class__ = D
        sage: orig == type(c)
        False
        sage: loads(dumps(c))
        Traceback (most recent call last):
        ...
        AssertionError: <class '__main__.D'> is not a direct subclass of <class 'sage.misc.fast_methods.Singleton'>
    """
    @staticmethod
    def __classcall__(cls):
        """
        Create an instance ``O`` of the given class ``cls``, and make it
        so that in future both ``cls.__call__`` and ``O.__class__.__call__``
        are constant functions returning ``O``.

        EXAMPLES::

            sage: from sage.misc.fast_methods import Singleton
            sage: class C(Singleton, Parent):
            ....:     def __init__(self):
            ....:         print("creating singleton")
            ....:         Parent.__init__(self, base=ZZ, category=Rings())
            sage: c = C()
            creating singleton
            sage: import __main__      # This is only needed ...
            sage: __main__.C = C       # ... in doctests
            sage: loads(dumps(c)) is copy(c) is C()  # indirect doctest
            True
        """
        assert cls.mro()[1] == Singleton, "{} is not a direct subclass of {}".format(cls, Singleton)
        res = typecall(cls)
        cf = ConstantFunction(res)
        cls._set_classcall(cf)
        res.__class__._set_classcall(cf)
        return res

    def __copy__(self):
        """
        There is a unique instance of a singleton, hence, copying
        returns ``self``.

        EXAMPLES::

            sage: from sage.misc.fast_methods import Singleton
            sage: class C(Singleton, Parent):
            ....:     def __init__(self):
            ....:         print("creating singleton")
            ....:         Parent.__init__(self, base=ZZ, category=Rings())
            sage: c = C()
            creating singleton
            sage: import __main__      # This is only needed ...
            sage: __main__.C = C       # ... in doctests
            sage: loads(dumps(c)) is copy(c) is C()  # indirect doctest
            True
        """
        return self

    def __reduce__(self):
        """
        There is a unique instance of a singleton, hence, pickling
        returns ``self``.

        EXAMPLES::

            sage: from sage.misc.fast_methods import Singleton
            sage: class C(Singleton, Parent):
            ....:     def __init__(self):
            ....:         print("creating singleton")
            ....:         Parent.__init__(self, base=ZZ, category=Rings())
            ....:
            sage: c = C()
            creating singleton
            sage: import __main__      # This is only needed ...
            sage: __main__.C = C       # ... in doctests
            sage: loads(dumps(c)) is copy(c) is C()  # indirect doctest
            True

        The pickle data mainly consist of the class of the unique instance,
        which may be a subclass of the original class used to create the
        instance.If the class is replaced by a sub-sub-class after creation
        of the instance, pickling fails. See the doctest
        in :class:`Singleton`.
        """
        return self.__class__, ()
