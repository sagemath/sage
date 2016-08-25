"""
Python 2 and 3 Compatibility
"""
from __future__ import absolute_import

from .six import *

def with_metaclass(meta, *bases):
    """
    Create a base class with a metaclass.

    TESTS::

        sage: from sage.misc.six import with_metaclass
        sage: class Meta(type): pass
        sage: class X(with_metaclass(Meta)): pass
        sage: type(X) is Meta
        True
        sage: issubclass(X, object)
        True
        sage: class Base(object): pass
        sage: class X(with_metaclass(Meta, Base)): pass
        sage: type(X) is Meta
        True
        sage: issubclass(X, Base)
        True
        sage: class Base2(object): pass
        sage: class X(with_metaclass(Meta, Base, Base2)): pass
        sage: type(X) is Meta
        True
        sage: issubclass(X, Base)
        True
        sage: issubclass(X, Base2)
        True
        sage: X.__mro__ == (X, Base, Base2, object)
        True

    Check that :trac:`18503` is fixed, i.e. that with_metaclass
    works with cdef'ed metaclasses::

        sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
        sage: class X(with_metaclass(ClasscallMetaclass)): pass
        sage: type(X) is ClasscallMetaclass
        True
        sage: X.__mro__ == (X, object)
        True
    """
    # This requires a bit of explanation: the basic idea is to make a dummy
    # metaclass for one level of class instantiation that replaces itself with
    # the actual metaclass.
    class metaclass(type):
        def __new__(cls, name, _, __dict__):
            return meta.__new__(meta, name, bases, __dict__)
    return type.__new__(metaclass, 'temporary_class', (), {})
