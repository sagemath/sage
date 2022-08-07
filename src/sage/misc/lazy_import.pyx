r"""
Lazy imports

This module allows one to lazily import objects into a namespace,
where the actual import is delayed until the object is actually called
or inspected. This is useful for modules that are expensive to import
or may cause circular references, though there is some overhead in its
use.

EXAMPLES::

    sage: lazy_import('sage.rings.all', 'ZZ')
    sage: type(ZZ)
    <class 'sage.misc.lazy_import.LazyImport'>
    sage: ZZ(4.0)
    4

By default, a warning is issued if a lazy import module is resolved
during Sage's startup. In case a lazy import's sole purpose is to
break a circular reference and it is known to be resolved at startup
time, one can use the ``at_startup`` option::

    sage: lazy_import('sage.rings.all', 'ZZ', at_startup=True)

This option can also be used as an intermediate step toward not
importing by default a module that is used in several places, some of
which can already afford to lazy import the module but not all.

A lazy import that is marked as "at_startup" will print a message if
it is actually resolved after the startup, so that the developer knows
that (s)he can remove the flag::

    sage: ZZ
    doctest:warning...
    UserWarning: Option ``at_startup=True`` for lazy import ZZ not needed anymore
    Integer Ring

.. SEEALSO:: :func:`lazy_import`, :class:`LazyImport`

AUTHOR:

 - Robert Bradshaw
"""

# ****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Keep OLD division semantics for Python 2 compatibility, such that
# lazy imports support old and true division.

cimport cython
from cpython.object cimport PyObject_RichCompare
from cpython.number cimport PyNumber_TrueDivide, PyNumber_Power, PyNumber_Index

cdef extern from *:
    int likely(int) nogil  # Defined by Cython

import os
import pickle
from warnings import warn
import inspect
from . import sageinspect

from sage.features import FeatureNotPresentError

cdef inline obj(x):
    if type(x) is LazyImport:
        return (<LazyImport>x).get_object()
    else:
        return x


# boolean to determine whether Sage is still starting up
cdef bint startup_guard = True

cdef bint finish_startup_called = False


cpdef finish_startup():
    """
    Finish the startup phase.

    This function must be called exactly once at the end of the Sage
    import process (:mod:`~sage.all`).

    TESTS::

        sage: from sage.misc.lazy_import import finish_startup
        sage: finish_startup()
        Traceback (most recent call last):
        ...
        AssertionError: finish_startup() must be called exactly once
    """
    global startup_guard, finish_startup_called
    assert startup_guard, 'finish_startup() must be called exactly once'
    startup_guard = False
    finish_startup_called = True


cpdef ensure_startup_finished():
    """
    Make sure that the startup phase is finished.

    In contrast to :func:`finish_startup`, this function can
    be called repeatedly.

    TESTS::

        sage: from sage.misc.lazy_import import ensure_startup_finished
        sage: ensure_startup_finished()
    """
    global startup_guard
    startup_guard = False


cpdef bint is_during_startup():
    """
    Return whether Sage is currently starting up.

    OUTPUT:

    Boolean

    TESTS::

        sage: from sage.misc.lazy_import import is_during_startup
        sage: is_during_startup()
        False
    """
    global startup_guard
    return startup_guard


cpdef test_fake_startup():
    """
    For testing purposes only.

    Switch the startup lazy import guard back on.

    EXAMPLES::

        sage: sage.misc.lazy_import.test_fake_startup()
        sage: lazy_import('sage.rings.all', 'ZZ', 'my_ZZ')
        sage: my_ZZ(123)
        doctest:warning...
        UserWarning: Resolving lazy import ZZ during startup
        123
        sage: sage.misc.lazy_import.finish_startup()
    """
    global startup_guard, finish_startup_called
    startup_guard = True
    finish_startup_called = False


@cython.final
cdef class LazyImport():
    """
    EXAMPLES::

        sage: from sage.misc.lazy_import import LazyImport
        sage: my_integer = LazyImport('sage.rings.all', 'Integer')
        sage: my_integer(4)
        4
        sage: my_integer('101', base=2)
        5
        sage: my_integer(3/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    cdef readonly _object  # The actual object if imported, None otherwise
    cdef _module
    cdef _name
    cdef _as_name
    cdef _namespace
    cdef bint _at_startup
    cdef _deprecation
    cdef _feature

    def __init__(self, module, name, as_name=None, at_startup=False, namespace=None,
                 deprecation=None, feature=None):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: type(lazy_ZZ)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: lazy_ZZ._get_object() is ZZ
            True
            sage: type(lazy_ZZ)
            <class 'sage.misc.lazy_import.LazyImport'>
        """
        self._object = None
        self._module = module
        self._name = name
        self._as_name = name if as_name is None else as_name
        self._namespace = namespace
        self._at_startup = at_startup
        self._deprecation = deprecation
        self._feature = feature

    cdef inline get_object(self):
        """
        Faster, Cython-only partially-inlined version of ``_get_object``.
        """
        if likely(self._object is not None):
            return self._object
        return self._get_object()

    cpdef _get_object(self):
        """
        Return the wrapped object, importing it if necessary.

        OUTPUT:

        - the wrapped object

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer_ring = LazyImport('sage.rings.all', 'ZZ')
            sage: my_integer_ring._object is None
            True
            sage: my_integer_ring._get_object()
            Integer Ring
            sage: my_integer_ring._object is None
            False
            sage: my_rats = LazyImport('sage.rings.rational_field', 'QQ', at_startup=True)
            sage: my_rats
            doctest:warning...
            UserWarning: Option ``at_startup=True`` for lazy import QQ not needed anymore
            Rational Field
        """
        if self._object is not None:
            return self._object

        if startup_guard and not self._at_startup:
            warn(f"Resolving lazy import {self._name} during startup")
        elif self._at_startup and not startup_guard:
            if finish_startup_called:
                warn(f"Option ``at_startup=True`` for lazy import {self._name} not needed anymore")

        try:
            self._object = getattr(__import__(self._module, {}, {}, [self._name]), self._name)
        except ImportError as e:
            if self._feature:
                raise FeatureNotPresentError(self._feature, reason=f'Importing {self._name} failed: {e}')
            raise

        if self._deprecation is not None:
            from sage.misc.superseded import deprecation_cython as deprecation
            try:
                trac_number, message = self._deprecation
            except TypeError:
                trac_number = self._deprecation
                import_command = f'from {self._module} import {self._name}'
                if self._as_name != self._name:
                    import_command += f' as {self._as_name}'
                message = f'\nImporting {self._as_name} from here is deprecated; please use "{import_command}" instead.'
            deprecation(trac_number, message)
        # Replace the lazy import in the namespace by the actual object
        name = self._as_name
        if self._namespace is not None:
            if self._namespace.get(name) is self:
                self._namespace[name] = self._object
        return self._object

    def _get_deprecation_ticket(self):
        """
        Return the ticket number of the deprecation, or 0 if this lazy
        import is not deprecated.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: H = LazyImport('sage.categories.homsets', 'Homsets')
            sage: H._get_deprecation_ticket()
            0
            sage: H = LazyImport('sage.categories.homsets', 'Homsets', deprecation=10668)
            sage: H._get_deprecation_ticket()
            10668
            sage: H = LazyImport('sage.categories.homsets', 'Homsets', deprecation=(10668, "this is deprecated"))
            sage: H._get_deprecation_ticket()
            10668
        """
        if self._deprecation is None:
            return 0
        try:
            return self._deprecation[0]
        except TypeError:
            return self._deprecation

    def _instancedoc_(self):
        """
        Return the docstring of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime.__doc__ is is_prime.__doc__
            True

        TESTS:

        Check that :trac:`19475` is fixed::

            sage: 'A subset of the real line' in RealSet.__doc__
            True
        """
        return sageinspect.sage_getdoc_original(self.get_object())

    def _sage_src_(self):
        """
        Return the source of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: 'def is_prime(' in my_isprime._sage_src_()
            True
        """
        return sageinspect.sage_getsource(self.get_object())

    def _sage_argspec_(self):
        """
        Return the argspec of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: rm = LazyImport('sage.all', 'random_matrix')
            sage: rm._sage_argspec_()
            ArgSpec(args=['ring', 'nrows', 'ncols', 'algorithm', 'implementation'], varargs='args', keywords='kwds', defaults=(None, 'randomize', None))
        """
        return sageinspect.sage_getargspec(self.get_object())

    def __getattr__(self, attr):
        """
        Attribute lookup on self defers to attribute lookup on the
        wrapped object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer = LazyImport('sage.rings.all', 'Integer')
            sage: my_integer.sqrt is Integer.sqrt
            True
        """
        return getattr(self.get_object(), attr)

    # We need to wrap all the slot methods, as they are not forwarded
    # via getattr.

    def __dir__(self):
        """
        Tab completion on self defers to completion on the wrapped
        object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: dir(lazy_ZZ) == dir(ZZ)
            True
        """
        return dir(self.get_object())

    def __call__(self, *args, **kwds):
        """
        Calling self calls the wrapped object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: is_prime(12) == my_isprime(12)
            True
            sage: is_prime(13) == my_isprime(13)
            True
        """
        return self.get_object()(*args, **kwds)

    def __repr__(self):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: repr(lazy_ZZ) == repr(ZZ)
            True
        """
        try:
            obj = self.get_object()
            return repr(obj)
        except FeatureNotPresentError as e:
            return "Failed lazy import:\n" + str(e)

    def __str__(self):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: str(lazy_ZZ) == str(ZZ)
            True
        """
        return str(self.get_object())

    def __bool__(self):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: bool(lazy_ZZ) == bool(ZZ)
            True
        """
        return bool(self.get_object())

    def __hash__(self):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: hash(lazy_ZZ) == hash(ZZ)
            True
        """
        return hash(self.get_object())

    def __richcmp__(left, right, int op):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: lazy_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: lazy_ZZ == ZZ
            True
            sage: lazy_ZZ == RR
            False
        """
        return PyObject_RichCompare(obj(left), obj(right), op)

    def __len__(self):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: len(version_info)
            5
        """
        return len(self.get_object())

    def __get__(self, instance, owner):
        """
        EXAMPLES:

        Here we show how to take a function in a module, and lazy
        import it as a method of a class. For the sake of this
        example, we add manually a function in sage.all::

            sage: def my_method(self): return self
            sage: sage.all.my_method = my_method

        Now we lazy import it as a method of a new class ``Foo``::

            sage: from sage.misc.lazy_import import LazyImport
            sage: class Foo():
            ....:     my_method = LazyImport('sage.all', 'my_method')

        Now we can use it as a usual method::

            sage: Foo().my_method()
            <__main__.Foo object at ...>
            sage: Foo.my_method
            <function my_method at 0x...>
            sage: Foo().my_method
            <bound method my_method of <__main__.Foo object at ...>>

        When a :class:`LazyImport` method is a method (or attribute)
        of a class, then extra work must be done to replace this
        :class:`LazyImport` object with the actual object. See the
        documentation of :meth:`_get_object` for an explanation of
        this.

        .. NOTE::

           For a :class:`LazyImport` object that appears in a class
           namespace, we need to do something special. Indeed, the
           class namespace dictionary at the time of the class
           definition is not the one that actually gets used. Thus,
           ``__get__`` needs to manually modify the class dict::

               sage: class Foo():
               ....:     lazy_import('sage.all', 'plot')
               sage: class Bar(Foo):
               ....:     pass
               sage: type(Foo.__dict__['plot'])
               <class 'sage.misc.lazy_import.LazyImport'>

           We access the ``plot`` method::

               sage: Bar.plot
               <function plot at 0x...>

           Now ``plot`` has been replaced in the dictionary of ``Foo``::

               sage: type(Foo.__dict__['plot'])
               <... 'function'>
        """
        # Don't use the namespace of the class definition
        self._namespace = None
        obj = self.get_object()

        name = self._as_name
        for cls in inspect.getmro(owner):
            if cls.__dict__.get(name) is self:
                setattr(cls, name, obj)
                break

        # Check whether the imported object is itself a descriptor
        try:
            get = obj.__get__
        except AttributeError:
            return obj
        else:
            return get(instance, owner)

    def __getitem__(self, key):
        """
        TESTS::

            sage: import sys
            sage: py_version = sys.version_info[0]
            sage: lazy_import('sys', 'version_info')
            sage: version_info[0] == py_version
            True
        """
        return self.get_object()[key]

    def __setitem__(self, key, value):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: sage.all.foo = list(range(10))
            sage: lazy_foo = LazyImport('sage.all', 'foo')
            sage: lazy_foo[1] = 100
            sage: print(lazy_foo)
            [0, 100, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: sage.all.foo
            [0, 100, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        self.get_object()[key] = value

    def __delitem__(self, key):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: sage.all.foo = list(range(10))
            sage: lazy_foo = LazyImport('sage.all', 'foo')
            sage: del lazy_foo[1]
            sage: print(lazy_foo)
            [0, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: print(sage.all.foo)
            [0, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        del self.get_object()[key]

    def __iter__(self):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: iter(version_info)
            <...iterator object at ...>
        """
        return iter(self.get_object())

    def __contains__(self, item):
        """
        TESTS::

            sage: import sys
            sage: py_version = sys.version_info[0]
            sage: lazy_import('sys', 'version_info')
            sage: py_version in version_info
            True

            sage: lazy_import('sys', 'version_info')
            sage: 2000 not in version_info
            True
        """
        return item in self.get_object()

    def __add__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo + 1
            11
        """
        return obj(left) + obj(right)

    def __sub__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo - 1
            9
        """
        return obj(left) - obj(right)

    def __mul__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo * 2
            20
        """
        return obj(left) * obj(right)

    def __matmul__(left, right):
        """
        TESTS::

            sage: from sympy import Matrix
            sage: sage.all.foo = Matrix([[1,1],[0,1]])
            sage: lazy_import('sage.all', 'foo')
            sage: foo.__matmul__(foo)
            Matrix([
            [1, 2],
            [0, 1]])
        """
        return obj(left) @ obj(right)

    def __floordiv__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo  // 3
            3
        """
        return obj(left) // obj(right)

    def __truediv__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: operator.truediv(foo, 3)
            10/3
        """
        return PyNumber_TrueDivide(obj(left), obj(right))

    def __pow__(left, right, mod):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo ** 2
            100
        """
        return PyNumber_Power(obj(left), obj(right), obj(mod))

    def __mod__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo % 7
            3
        """
        return obj(left) % obj(right)

    def __lshift__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo << 3
            80
        """
        return obj(left) << obj(right)

    def __rshift__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo >> 2
            2
        """
        return obj(left) >> obj(right)

    def __and__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo & 7
            2
        """
        return obj(left) & obj(right)

    def __or__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo | 7
            15
        """
        return obj(left) | obj(right)

    def __xor__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: foo ^^ 7
            13
        """
        return obj(left) ^ obj(right)

    def __neg__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: -foo
            -10
        """
        return -self.get_object()

    def __pos__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: +foo
            10
        """
        return +self.get_object()

    def __abs__(self):
        """
        TESTS::

            sage: sage.all.foo = -1000
            sage: lazy_import('sage.all', 'foo')
            sage: abs(foo)
            1000
        """
        return abs(self.get_object())

    def __invert__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: ~foo
            1/10
        """
        return ~self.get_object()

    def __complex__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: complex(foo)
            (10+0j)
        """
        return complex(self.get_object())

    def __int__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: int(foo)
            10
        """
        return int(self.get_object())

    def __float__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: float(foo)
            10.0
        """
        return float(self.get_object())

    def __oct__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: oct(foo)
            '0o12'
        """
        return self.get_object().__oct__()

    def __hex__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: hex(foo)
            '0xa'
        """
        return self.get_object().__hex__()

    def __index__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: list(range(100))[foo]
            10
        """
        return PyNumber_Index(self.get_object())

    def __copy__(self):
        """
        Support copy()

        TESTS::


            sage: from sage.misc.lazy_import import LazyImport
            sage: sage.all.foo = [[1,2], 3]
            sage: lazy_foo = LazyImport('sage.all', 'foo')
            sage: a = copy(lazy_foo)
            sage: a is sage.all.foo        # copy
            False
            sage: a[0] is sage.all.foo[0]  # copy but not deep
            True
            sage: type(lazy_foo) is LazyImport
            True
        """
        import copy
        return copy.copy(self.get_object())

    def __deepcopy__(self, memo=None):
        """
        Support copy()

        TESTS::

            sage: from sage.misc.lazy_import import LazyImport
            sage: sage.all.foo = [[1,2], 3]
            sage: lazy_foo = LazyImport('sage.all', 'foo')
            sage: a = deepcopy(lazy_foo)
            sage: a is sage.all.foo        # copy
            False
            sage: a[0] is sage.all.foo[0]  # deep copy
            False
            sage: type(lazy_foo) is LazyImport
            True
        """
        import copy
        return copy.deepcopy(self.get_object())

    def __instancecheck__(self, x):
        """
        Support ``isinstance()``.

        EXAMPLES::

            sage: lazy_import('sage.rings.rational_field', 'RationalField')
            sage: isinstance(QQ, RationalField)
            True

        No object is an instance of a class that cannot be imported::

            sage: lazy_import('sage.xxxxx_does_not_exist', 'DoesNotExist')
            sage: isinstance(QQ, DoesNotExist)
            False
        """
        try:
            return isinstance(x, self.get_object())
        except ImportError:
            return False

    def __subclasscheck__(self, x):
        """
        Support ``issubclass()``.

        EXAMPLES::

            sage: lazy_import('sage.structure.parent', 'Parent')
            sage: issubclass(RationalField, Parent)
            True

        No class is a subclass of a class that cannot be imported::

            sage: lazy_import('sage.xxxxx_does_not_exist', 'DoesNotExist')
            sage: issubclass(RationalField, DoesNotExist)
            False
        """
        try:
            return issubclass(x, self.get_object())
        except ImportError:
            return False


def lazy_import(module, names, as_=None, *,
    at_startup=False, namespace=None, deprecation=None, feature=None):
    """
    Create a lazy import object and inject it into the caller's global
    namespace. For the purposes of introspection and calling, this is
    like performing a lazy "from module import name" where the import
    is delayed until the object actually is used or inspected.

    INPUT:

    - ``module`` -- a string representing the module to import

    - ``names`` -- a string or list of strings representing the names to
      import from module

    - ``as_`` -- (optional) a string or list of strings representing the
      names of the objects in the importing module. This is analogous to
      ``from ... import ... as ...``.

    - ``at_startup`` -- a boolean (default: ``False``);
      whether the lazy import is supposed to be resolved at startup time

    - ``namespace`` -- the namespace where importing the names; by default,
      import the names to current namespace

    - ``deprecation`` -- (optional) if not ``None``, a deprecation warning
      will be issued when the object is actually imported;
      ``deprecation`` should be either a trac number (integer) or a
      pair ``(trac_number, message)``

    - ``feature`` -- a python module (optional), if it cannot be imported
      an appropriate error is raised

    .. SEEALSO:: :mod:`sage.misc.lazy_import`, :class:`LazyImport`

    EXAMPLES::

        sage: lazy_import('sage.rings.all', 'ZZ')
        sage: type(ZZ)
        <class 'sage.misc.lazy_import.LazyImport'>
        sage: ZZ(4.0)
        4
        sage: lazy_import('sage.rings.all', 'RDF', 'my_RDF')
        sage: my_RDF._get_object() is RDF
        True
        sage: my_RDF(1/2)
        0.5

        sage: lazy_import('sage.all', ['QQ', 'RR'], ['my_QQ', 'my_RR'])
        sage: my_QQ._get_object() is QQ
        True
        sage: my_RR._get_object() is RR
        True

    Upon the first use, the object is injected directly into
    the calling namespace::

        sage: lazy_import('sage.all', 'ZZ', 'my_ZZ')
        sage: my_ZZ is ZZ
        False
        sage: my_ZZ(37)
        37
        sage: my_ZZ is ZZ
        True

    We check that :func:`lazy_import` also works for methods::

        sage: class Foo():
        ....:     lazy_import('sage.all', 'plot')
        sage: class Bar(Foo):
        ....:     pass
        sage: type(Foo.__dict__['plot'])
        <class 'sage.misc.lazy_import.LazyImport'>
        sage: 'EXAMPLES' in Bar.plot.__doc__
        True
        sage: type(Foo.__dict__['plot'])
        <... 'function'>

    If deprecated then a deprecation warning is issued::

        sage: lazy_import('sage.all', 'Qp', 'my_Qp', deprecation=14275)
        sage: my_Qp(5)
        doctest:...: DeprecationWarning:
        Importing my_Qp from here is deprecated; please use "from sage.all import Qp as my_Qp" instead.
        See http://trac.sagemath.org/14275 for details.
        5-adic Field with capped relative precision 20

    An example of deprecation with a message::

        sage: lazy_import('sage.all', 'Qp', 'my_Qp_msg', deprecation=(14275, "This is an example."))
        sage: my_Qp_msg(5)
        doctest:...: DeprecationWarning: This is an example.
        See http://trac.sagemath.org/14275 for details.
        5-adic Field with capped relative precision 20

    An example of an import relying on a feature::

        sage: from sage.features import PythonModule
        sage: lazy_import('ppl', 'equation', feature=PythonModule('ppl', spkg='pplpy'))
        sage: equation
        <built-in function equation>
        sage: lazy_import('PyNormaliz', 'NmzListConeProperties', feature=PythonModule('PyNormaliz', spkg='pynormaliz'))  # optional - pynormaliz
        sage: NmzListConeProperties  # optional - pynormaliz
        <built-in function NmzListConeProperties>
        sage: lazy_import('foo', 'not_there', feature=PythonModule('foo', spkg='non-existing-package'))
        sage: not_there
        Failed lazy import:
        foo is not available.
        Importing not_there failed: No module named 'foo'...
        No equivalent system packages for ... are known to Sage...
    """
    if as_ is None:
        as_ = names
    if isinstance(names, basestring):
        names = [names]
        as_ = [as_]
    else:
        names = list(names)
        as_ = list(as_)
    if namespace is None:
        namespace = inspect.currentframe().f_locals
    if "*" in names:
        ix = names.index("*")
        all = get_star_imports(module)
        names[ix:ix+1] = all
        as_[ix:ix+1] = all
    for name, alias in zip(names, as_):
        namespace[alias] = LazyImport(module, name, alias, at_startup, namespace, deprecation, feature)


star_imports = None

def save_cache_file():
    """
    Used to save the cached * import names.

    TESTS::

        sage: import sage.misc.lazy_import
        sage: sage.misc.lazy_import.save_cache_file()
    """
    from sage.misc.temporary_file import atomic_write
    from .lazy_import_cache import get_cache_file

    global star_imports
    if star_imports is None:
        star_imports = {}
    cache_file = get_cache_file()
    cache_dir = os.path.dirname(cache_file)

    os.makedirs(cache_dir, exist_ok=True)
    with atomic_write(cache_file, binary=True) as f:
        pickle.dump(star_imports, f)

def get_star_imports(module_name):
    """
    Lookup the list of names in a module that would be imported with "import \\*"
    either via a cache or actually importing.

    EXAMPLES::

        sage: from sage.misc.lazy_import import get_star_imports
        sage: 'get_star_imports' in get_star_imports('sage.misc.lazy_import')
        True
        sage: 'EllipticCurve' in get_star_imports('sage.schemes.all')
        True

    TESTS::

        sage: import os, tempfile
        sage: fd, cache_file = tempfile.mkstemp()
        sage: os.write(fd, b'invalid')
        7
        sage: os.close(fd)
        sage: import sage.misc.lazy_import as lazy
        sage: import sage.misc.lazy_import_cache as cache
        sage: cache.get_cache_file = (lambda: cache_file)
        sage: lazy.star_imports = None
        sage: lazy.get_star_imports('sage.schemes.all')
        doctest:...: UserWarning: star_imports cache is corrupted
        [...]
        sage: os.remove(cache_file)
    """
    global star_imports
    if star_imports is None:
        from .lazy_import_cache import get_cache_file
        star_imports = {}
        try:
            with open(get_cache_file(), "rb") as cache_file:
                star_imports = pickle.load(cache_file)
        except IOError:        # file does not exist
            pass
        except Exception:  # unpickling failed
            import warnings
            warnings.warn('star_imports cache is corrupted')
    try:
        return star_imports[module_name]
    except KeyError:
        module = __import__(module_name, {}, {}, ["*"])
        if hasattr(module, "__all__"):
            all = module.__all__
        else:
            all = [key for key in dir(module) if key[0] != "_"]
        star_imports[module_name] = all
        return all


# Add support for _instancedoc_
from sage.misc.instancedoc import instancedoc
instancedoc(LazyImport)
