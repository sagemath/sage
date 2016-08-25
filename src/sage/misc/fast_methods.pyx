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

#******************************************************************************
#  Copyright (C) 2013 Simon A. King <simon.king at uni-jena.de>
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
#                  http://www.gnu.org/licenses/
#******************************************************************************
from __future__ import print_function

from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.constant_function import ConstantFunction
from sage.misc.lazy_attribute import lazy_class_attribute

from cpython.bool cimport *
from cpython.ref cimport *

cdef extern from "Python.h":
    cdef size_t SIZEOF_VOID_P

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
        ...     def __init__(self, x):
        ...         self.x = x
        ...     def __cmp__(self,other):
        ...         return cmp(self.x^2,other.x^2)
        ...     def __hash__(self):
        ...         return hash(self.x)
        sage: class MyUniqueParent(UniqueRepresentation, MyParent): pass
        sage: issubclass(MyUniqueParent, sage.misc.fast_methods.WithEqualityById)
        True

    Inheriting from :class:`WithEqualityById` provides unique representation
    behaviour. In particular, the comparison inherited from ``MyParent``
    is overloaded::

        sage: a = MyUniqueParent(1)
        sage: b = MyUniqueParent(2)
        sage: c = MyUniqueParent(1)
        sage: a is c
        True
        sage: d = MyUniqueParent(-1)
        sage: a == d
        False

    Note, however, that Python distinguishes between "comparison by cmp"
    and "comparison by binary relations"::

        sage: cmp(a,d)
        0

    The comparison inherited from ``MyParent`` will be used in those cases
    in which identity does not give sufficient information to find the relation::

        sage: a < b
        True
        sage: b > d
        True

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
        ...     def __eq__(self, other):
        ...         return self.x^2 == other.x^2
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
        The hash provided by this class coincides with that of ``<type 'object'>``.

        TESTS::

            sage: class MyParent(Parent):
            ...     def __init__(self, x):
            ...         self.x = x
            ...     def __cmp__(self,other):
            ...         return cmp(self.x^2,other.x^2)
            ...     def __hash__(self):
            ...         return hash(self.x)
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

    def __richcmp__(self, other, int m):
        """
        Equality test provided by this class is by identity.

        TESTS::

            sage: class MyParent(Parent):
            ...     def __init__(self, x):
            ...         self.x = x
            ...     def __cmp__(self,other):
            ...         return cmp(self.x^2,other.x^2)
            ...     def __hash__(self):
            ...         return hash(self.x)
            sage: class MyUniqueParent(UniqueRepresentation, MyParent): pass
            sage: issubclass(MyUniqueParent, sage.misc.fast_methods.WithEqualityById)
            True
            sage: a = MyUniqueParent(1)
            sage: b = MyUniqueParent(-1)

        Comparison with ``cmp`` is still using what is inherited
        from ``MyParent``::

            sage: cmp(a,b)
            0

        However, equality test takes into account identity::

            sage: a == b
            False

        In cases in which rich comparison by identity gives no final answer,
        the comparison inherited from ``MyParent`` is consulted again::

            sage: a <= b and b >= a
            True
            sage: a < b
            False

        """
        cdef object out
        if self is other:
            if m == 2: # ==
                return True
            elif m == 3: # !=
                return False
            else:
                # <= or >= or NotImplemented
                return m==1 or m==5 or NotImplemented
        else:
            if m == 2:
                return False
            elif m == 3:
                return True
            else:
                return NotImplemented



cdef class FastHashable_class:
    """
    A class that has a fast hash method, returning a pre-assigned value.

    NOTE:

    This is for internal use only. The class has a cdef attribute ``_hash``,
    that needs to be assigned (for example, by calling the init method, or by
    a direct assignement using cython). This is slower than using
    :func:`provide_hash_by_id`, but has the advantage that the hash can be
    prescribed, by assigning a cdef attribute ``_hash``.

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

class Singleton(WithEqualityById):
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
        AssertionError: (("<class '__main__.D'> is not a direct
        subclass of <class 'sage.misc.fast_methods.Singleton'>",),
        <class '__main__.D'>, ())
    """
    __metaclass__ = ClasscallMetaclass

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

