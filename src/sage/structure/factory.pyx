# sage_setup: distribution = sagemath-objects
r"""
Factory for cached representations

.. SEEALSO::

    :mod:`sage.structure.unique_representation`

Using a :class:`UniqueFactory` is one way of implementing a *cached
representation behaviour*. In spite of its name, using a
:class:`UniqueFactory` is not enough to ensure the *unique representation
behaviour*. See :mod:`~sage.structure.unique_representation` for a
detailed explanation.

With a :class:`UniqueFactory`, one can preprocess the given arguments. There
is special support for specifying a subset of the arguments that serve as the
unique key, so that still *all* given arguments are used to create a new
instance, but only the specified subset is used to look up in the
cache. Typically, this is used to construct objects that accept an optional
``check=[True|False]`` argument, but whose result should be unique
regardless of said optional argument. (This use case should be handled with
care, though: Any checking which isn't done in the ``create_key`` or
``create_key_and_extra_args`` method will be done only when a new object is
generated, but not when a cached object is retrieved from cache.
Consequently, if the factory is once called with ``check=False``, a
subsequent call with ``check=True`` cannot be expected to perform all checks
unless these checks are all in the ``create_key`` or
``create_key_and_extra_args`` method.)

For a class derived from
:class:`~sage.structure.unique_representation.CachedRepresentation`, argument
preprocessing can be obtained by providing a custom static ``__classcall__``
or ``__classcall_private__`` method, but this seems less transparent. When
argument preprocessing is not needed or the preprocess is not very
sophisticated, then generally
:class:`~sage.structure.unique_representation.CachedRepresentation` is much
easier to use than a factory.

AUTHORS:

- Robert Bradshaw (2008): initial version.
- Simon King (2013): extended documentation.
- Julian Rueth (2014-05-09): use ``_cache_key`` if parameters are unhashable
"""

#*****************************************************************************
#  Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#                2014 Julian Rueth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import types

from sage.structure.sage_object cimport SageObject

cdef sage_version
from sage.version import version as sage_version
from sage.cpython.string cimport bytes_to_str

sage_version = sage_version.split('.')
for i in range(len(sage_version)):
    try:
        sage_version[i] = int(sage_version[i])
    except ValueError:
        pass
sage_version = tuple(sage_version)

cimport sage.misc.weak_dict
from sage.misc.cachefunc cimport cache_key as _cache_key


cdef class UniqueFactory(SageObject):
    """
    This class is intended to make it easy to cache objects.

    It is based on the idea that the object is uniquely defined by a set of
    defining data (the key). There is also the possibility of some
    non-defining data (extra args) which will be used in initial creation,
    but not affect the caching.

    .. WARNING::

        This class only provides *cached representation behaviour*. Hence,
        using :class:`UniqueFactory`, it is still possible to create distinct
        objects that evaluate equal. Unique representation behaviour can be
        added, for example, by additionally inheriting from
        :class:`sage.misc.fast_methods.WithEqualityById`.

    The objects created are cached (using weakrefs) based on their key and
    returned directly rather than re-created if requested again. Pickling is
    taken care of by the factory, and will return the same object for the same
    version of Sage, and distinct (but hopefully equal) objects for different
    versions of Sage.

    .. WARNING::

        The objects returned by a :class:`UniqueFactory` must be instances of
        new style classes (hence, they must be instances of :class:`object`)
        that must not only allow a weak reference, but must accept general
        attribute assignment. Otherwise, pickling won't work.

    USAGE:

    A *unique factory* provides a way to create objects from parameters
    (the type of these objects can depend on the parameters, and is often
    determined only at runtime) and to cache them by a certain key
    derived from these parameters, so that when the factory is being
    called again with the same parameters (or just with parameters which
    yield the same key), the object is being returned from cache rather
    than constructed anew.

    An implementation of a unique factory consists of a factory class and
    an instance of this factory class.

    The factory class has to be a class inheriting from ``UniqueFactory``.
    Typically it only needs to implement :meth:`create_key` (a method that
    creates a key from the given parameters, under which key the object
    will be stored in the cache) and :meth:`create_object` (a method that
    returns the actual object from the key). Sometimes, one would also
    implement :meth:`create_key_and_extra_args` (this differs from
    :meth:`create_key` in allowing to also create some additional
    arguments from the given parameters, which arguments then get passed
    to :meth:`create_object` and thus can have an effect on the initial
    creation of the object, but do *not* affect the key) or
    :meth:`other_keys`. Other methods are not supposed to be overloaded.

    The factory class itself cannot be called to create objects. Instead,
    an instance of the factory class has to be created first. For
    technical reasons, this instance has to be provided with a name that
    allows Sage to find its definition. Specifically, the name of the
    factory instance (or the full path to it, if it is not in the global
    namespace) has to be passed to the factory class as a string variable.
    So, if our factory class has been called ``A`` and is located in
    ``sage/spam/battletoads.py``, then we need to define an instance (say,
    ``B``) of ``A`` by writing ``B = A("sage.spam.battletoads.B")``
    (or ``B = A("B")`` if this ``B`` will be imported into global
    namespace). This instance can then be used to create objects (by
    calling ``B(*parameters)``).

    Notice that the objects created by the factory don't inherit from the
    factory class. They *do* know about the factory that created them (this
    information, along with the keys under which this factory caches them,
    is stored in the ``_factory_data`` attributes of the objects), but not
    via inheritance.

    EXAMPLES:

    The below examples are rather artificial and illustrate particular
    aspects. For a "real-life" usage case of ``UniqueFactory``, see
    the finite field factory in :mod:`sage.rings.finite_rings.finite_field_constructor`.

    In many cases, a factory class is implemented by providing the two
    methods :meth:`create_key` and :meth:`create_object`. In our example,
    we want to demonstrate how to use "extra arguments" to choose a specific
    implementation, with preference given to an instance found in the cache,
    even if its implementation is different. Hence, we implement
    :meth:`create_key_and_extra_args` rather than :meth:`create_key`, putting
    the chosen implementation into the extra arguments. Then, in the
    :meth:`create_object` method, we create and return instances of the
    specified implementation.
    ::

        sage: from sage.structure.factory import UniqueFactory
        sage: class MyFactory(UniqueFactory):
        ....:     def create_key_and_extra_args(self, *args, **kwds):
        ....:         return args, {'impl':kwds.get('impl', None)}
        ....:     def create_object(self, version, key, **extra_args):
        ....:         impl = extra_args['impl']
        ....:         if impl == 'C':
        ....:             return C(*key)
        ....:         if impl == 'D':
        ....:             return D(*key)
        ....:         return E(*key)
        ....:

    Now we can create a factory instance. It is supposed to be found under the
    name ``'F'`` in the ``__main__`` module. Note that in an interactive
    session, ``F`` would automatically be in the ``__main__`` module. Hence,
    the second and third of the following four lines are only needed in
    doctests.  ::

        sage: F = MyFactory("__main__.F")
        sage: import __main__
        sage: __main__.F = F
        sage: loads(dumps(F)) is F
        True

    Now we create three classes ``C``, ``D`` and ``E``. The first is a Cython
    extension-type class that does not allow weak references nor attribute
    assignment. The second is a Python class that is not derived from
    :class:`object`. The third allows attribute assignment and is derived
    from :class:`object`.  ::

        sage: cython("cdef class C: pass")                                              # needs sage.misc.cython
        sage: class D:
        ....:     def __init__(self, *args):
        ....:         self.t = args
        ....:     def __repr__(self):
        ....:         return "D%s"%repr(self.t)
        ....:
        sage: class E(D, object): pass

    Again, being in a doctest, we need to put the class ``D`` into the
    ``__main__`` module, so that Python can find it::

        sage: import __main__
        sage: __main__.D = D

    It is impossible to create an instance of ``C`` with our factory, since it
    does not allow weak references::

        sage: F(1, impl='C')                                                            # needs sage.misc.cython
        Traceback (most recent call last):
        ...
        TypeError: cannot create weak reference to '....C' object

    Let us try again, with a Cython class that does allow weak
    references. Now, creation of an instance using the factory works::

        sage: cython(                                                                   # needs sage.misc.cython
        ....: '''
        ....: cdef class C:
        ....:     cdef __weakref__
        ....: ''')
        ....:
        sage: c = F(1, impl='C')                                                        # needs sage.misc.cython
        sage: isinstance(c, C)                                                          # needs sage.misc.cython
        True

    The cache is used when calling the factory again---even if it is suggested
    to use a different implementation. This is because the implementation is
    only considered an "extra argument" that does not count for the key.
    ::

        sage: c is F(1, impl='C') is F(1, impl='D') is F(1)                             # needs sage.misc.cython
        True

    However, pickling and unpickling does not use the cache. This is because
    the factory has tried to assign an attribute to the instance that provides
    information on the key used to create the instance, but failed::

        sage: loads(dumps(c)) is c                                                      # needs sage.misc.cython
        False
        sage: hasattr(c, '_factory_data')                                               # needs sage.misc.cython
        False

    We have already seen that our factory will only take the requested
    implementation into account if the arguments used as key have not been
    used yet. So, we use other arguments to create an instance of class
    ``D``::

        sage: d = F(2, impl='D')
        sage: isinstance(d, D)
        True

    The factory only knows about the pickling protocol used by new style
    classes. Hence, again, pickling and unpickling fails to use the cache,
    even though the "factory data" are now available (this is not the case
    on Python 3 which *only* has new style classes)::

        sage: loads(dumps(d)) is d
        True
        sage: d._factory_data
        (<__main__.MyFactory object at ...>,
         (...),
         (2,),
         {'impl': 'D'})

    Only when we have a new style class that can be weak referenced and allows
    for attribute assignment, everything works::

        sage: e = F(3)
        sage: isinstance(e, E)
        True
        sage: loads(dumps(e)) is e
        True
        sage: e._factory_data
        (<__main__.MyFactory object at ...>,
         (...),
         (3,),
         {'impl': None})
    """
    cdef readonly _name
    cdef readonly _cache

    def __init__(self, name):
        """
        INPUT:

        - ``name`` -- string; a name in the global namespace referring to
          ``self`` or a fully qualified path name to ``self``, which is used to
          locate the factory on unpickling

        EXAMPLES::

            sage: from sage.structure.factory import UniqueFactory
            sage: fake_factory = UniqueFactory('ZZ')
            sage: loads(dumps(fake_factory))
            Integer Ring
            sage: fake_factory = UniqueFactory('sage.rings.rational_field.QQ')
            sage: loads(dumps(fake_factory))
            Rational Field
        """
        self._name = name
        self._cache = sage.misc.weak_dict.WeakValueDictionary()

    def __reduce__(self):
        """
        EXAMPLES::

            sage: A = FiniteField(127)
            sage: A is loads(dumps(A)) # indirect doctest
            True

            sage: # needs sage.rings.finite_rings
            sage: B = FiniteField(3^3,'b')
            sage: B is loads(dumps(B))
            True
            sage: C = FiniteField(2^16,'c')
            sage: C is loads(dumps(C))
            True
            sage: D = FiniteField(3^20,'d')
            sage: D is loads(dumps(D))
            True

        TESTS::

            sage: loads(dumps(FiniteField)) is FiniteField
            True
            sage: from sage.structure.test_factory import test_factory
            sage: loads(dumps(test_factory)) is test_factory
            True
        """
        return lookup_global, (self._name,)

    def __call__(self, *args, **kwds):
        """
        This is the method invoked to create objects. It first creates a key
        from the given parameters, then if an object with that key already
        exists returns it, and otherwise creates one and stores a weak reference
        to it in its dictionary.

        Do not override this method. Instead, override ``create_key`` and
        ``create_object`` and put the docstring in the body of the class.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: _ = test_factory(1,2,3); _
            Making object (1, 2, 3)
            <sage.structure.test_factory.A object at ...>

        It already created one, so don't re-create::

            sage: test_factory(1,2,3)
            <sage.structure.test_factory.A object at ...>
            sage: test_factory(1,2,3) is test_factory(1,2,3)
            True

        Of course, with a different key, a new object will be created::

            sage: test_factory(1,2,3) is test_factory(1,2,4)
            Making object (1, 2, 4)
            False
        """
        key, kwds = self.create_key_and_extra_args(*args, **kwds)
        version = self.get_version(sage_version)
        return self.get_object(version, key, kwds)

    cpdef get_object(self, version, key, extra_args):
        """
        Return the object corresponding to ``key``, creating it with
        ``extra_args`` if necessary (for example, it isn't in the cache
        or it is unpickling from an older version of Sage).

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: a = test_factory.get_object(3.0, 'a', {}); a
            Making object a
            <sage.structure.test_factory.A object at ...>
            sage: test_factory.get_object(3.0, 'a', {}) is test_factory.get_object(3.0, 'a', {})
            True
            sage: test_factory.get_object(3.0, 'a', {}) is test_factory.get_object(3.1, 'a', {})
            Making object a
            False
            sage: test_factory.get_object(3.0, 'a', {}) is test_factory.get_object(3.0, 'b', {})
            Making object b
            False

        TESTS:

        Check that :issue:`16317` has been fixed, i.e., caching works for
        unhashable objects::

            sage: K.<u> = Qq(4)                                                         # needs sage.rings.padics
            sage: d = test_factory.get_object(3.0, (K(1), 'c'), {})                     # needs sage.rings.padics
            Making object (1 + O(2^20), 'c')
            sage: d is test_factory.get_object(3.0, (K(1), 'c'), {})                    # needs sage.rings.padics
            True
        """
        cache_key = key
        try:
            try:
                return self._cache[version, cache_key]
            except TypeError: # key is unhashable
                cache_key = _cache_key(cache_key)
                return self._cache[version, cache_key]
        except KeyError:
            pass
        obj = self.create_object(version, key, **extra_args)
        self._cache[version, cache_key] = obj
        try:
            for key in self.other_keys(key, obj):
                try:
                    self._cache[version, key] = obj
                except TypeError: # key is unhashable
                    self._cache[version, _cache_key(key)] = obj
            obj._factory_data = self, version, key, extra_args

            # Install a custom __reduce__ method on the instance "obj"
            # that we just created. We only do this if the class of
            # "obj" has a generic __reduce__ method, which is either
            # object.__reduce__ or __reduce_cython__, the
            # auto-generated pickling function for Cython.
            f = obj.__class__.__reduce__
            if f.__objclass__ is object or f.__name__ == "__reduce_cython__":
                obj.__reduce_ex__ = types.MethodType(generic_factory_reduce, obj)
        except AttributeError:
            pass
        return obj

    cpdef get_version(self, sage_version):
        """
        This is provided to allow more or less granular control over
        pickle versioning. Objects pickled in the same version of Sage
        will unpickle to the same rather than simply equal objects. This
        can provide significant gains as arithmetic must be performed on
        objects with identical parents. However, if there has been an
        incompatible change (e.g. in element representation) we want the
        version number to change so coercion is forced between the two
        parents.

        Defaults to the Sage version that is passed in, but coarser
        granularity can be provided.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.get_version((3,1,0))
            (3, 1, 0)
        """
        return sage_version

    def create_key_and_extra_args(self, *args, **kwds):
        r"""
        Return a tuple containing the key (uniquely defining data)
        and any extra arguments (empty by default).

        Defaults to :meth:`create_key`.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.create_key_and_extra_args(1, 2, key=5)
            ((1, 2), {})
            sage: GF.create_key_and_extra_args(3)
            ((3, ('x',), None, 'modn', 3, 1, True, None, None, None, True, False), {})
        """
        return self.create_key(*args, **kwds), {}

    def create_key(self, *args, **kwds):
        """
        Given the parameters (arguments and keywords), create a key
        that uniquely determines this object.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.create_key(1, 2, key=5)
            (1, 2)
        """
        raise NotImplementedError

    def create_object(self, version, key, **extra_args):
        """
        Create the object from the key and extra arguments. This is only
        called if the object was not found in the cache.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.create_object(0, (1,2,3))
            Making object (1, 2, 3)
            <sage.structure.test_factory.A object at ...>
            sage: test_factory('a')
            Making object ('a',)
            <sage.structure.test_factory.A object at ...>
            sage: test_factory('a') # NOT called again
            <sage.structure.test_factory.A object at ...>
        """
        raise NotImplementedError

    cpdef other_keys(self, key, obj):
        """
        Sometimes during object creation, certain defaults are chosen which
        may result in a new (more specific) key. This allows the more specific
        key to be regarded as equivalent to the original key returned by
        :meth:`create_key` for the purpose of lookup in the cache, and is used
        for pickling.

        EXAMPLES:

        The ``GF`` factory used to have a custom :meth:`other_keys`
        method, but this was removed in :issue:`16934`::

            sage: # needs sage.libs.linbox sage.rings.finite_rings
            sage: key, _ = GF.create_key_and_extra_args(27, 'k'); key
            (27, ('k',), x^3 + 2*x + 1, 'givaro', 3, 3, True, None, 'poly', True, True, True)
            sage: K = GF.create_object(0, key); K
            Finite Field in k of size 3^3
            sage: GF.other_keys(key, K)
            []

            sage: K = GF(7^40, 'a')                                                     # needs sage.rings.finite_rings
            sage: loads(dumps(K)) is K                                                  # needs sage.rings.finite_rings
            True
        """
        return []

    cpdef reduce_data(self, obj):
        """
        The results of this function can be returned from
        :meth:`__reduce__`. This is here so the factory internals can
        change without having to re-write :meth:`__reduce__` methods
        that use it.

        EXAMPLES::

            sage: # needs sage.modules
            sage: from sage.modules.free_module import FreeModuleFactory_with_standard_basis as F
            sage: V = F(ZZ, 5)
            sage: factory, data = F.reduce_data(V)
            sage: factory(*data)
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
            sage: factory(*data) is V
            True

            sage: from sage.structure.test_factory import test_factory
            sage: a = test_factory(1, 2)
            Making object (1, 2)
            sage: test_factory.reduce_data(a)
            (<built-in function generic_factory_unpickle>,
             (<sage.structure.test_factory.UniqueFactoryTester object at ...>,
              (...),
              (1, 2),
              {}))

        Note that the ellipsis ``(...)`` here stands for the Sage
        version.
        """
        return generic_factory_unpickle, obj._factory_data

# This is used to handle old UniqueFactory pickles
factory_unpickles = {}


def register_factory_unpickle(name, callable):
    """
    Register a callable to handle the unpickling from an old
    :class:`UniqueFactory` object.

    :class:`UniqueFactory` pickles use a global name through
    :func:`generic_factory_unpickle()`, so the usual
    :func:`~sage.misc.persist.register_unpickle_override()`
    cannot be used here.

    .. SEEALSO::

        :func:`generic_factory_unpickle()`

    TESTS:

    This is similar to the example given in
    :func:`generic_factory_unpickle()`, but here we will use a function to
    explicitly return a polynomial ring.

    First, we create the factory. In a doctest, it is needed to explicitly put
    it into ``__main__``, so that it can be located when pickling. Also, it is
    needed that we work with a new-style class::

        sage: from sage.structure.factory import UniqueFactory, register_factory_unpickle
        sage: import __main__
        sage: class OldStuff():
        ....:     def __init__(self, n, **extras):
        ....:         self.n = n
        ....:     def __repr__(self):
        ....:         return "Rotten old thing of level {}".format(self.n)
        sage: __main__.OldStuff = OldStuff
        sage: class MyFactory(UniqueFactory):
        ....:     def create_object(self, version, key, **extras):
        ....:         return OldStuff(key[0])
        ....:     def create_key(self, *args):
        ....:         return args
        sage: G = MyFactory('__main__.G')
        sage: __main__.G = G
        sage: a = G(3); a
        Rotten old thing of level 3
        sage: loads(dumps(a)) is a
        True

    Now, we create a pickle (the string returned by ``dumps(a)``)::

        sage: s = dumps(a)

    We create the function which will handle the unpickling::

        sage: def foo(n, **kwds):
        ....:     return PolynomialRing(QQ, n, 'x')
        sage: register_factory_unpickle('__main__.G', foo)

    The old pickle correctly unpickles as an explicit polynomial ring::

        sage: loads(s)
        Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
    """
    #global factory_unpickles
    factory_unpickles[name] = callable


def generic_factory_unpickle(factory, *args):
    """
    Method used for unpickling the object.

    The unpickling mechanism needs a plain Python function to call.
    It takes a factory as the first argument, passes the rest of the
    arguments onto the factory's :meth:`UniqueFactory.get_object` method.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.modules.free_module import FreeModuleFactory_with_standard_basis as F
        sage: V = F(ZZ, 5)
        sage: func, data = F.reduce_data(V)
        sage: func is sage.structure.factory.generic_factory_unpickle
        True
        sage: sage.structure.factory.generic_factory_unpickle(*data) is V
        True

    TESTS:

    The following was enabled in :issue:`16349`. Suppose we have defined
    (somewhere in the library of an old Sage version) a unique factory; in our
    example below, it returns polynomial rings. Now suppose that we want to
    replace the factory by something else, say, a class that provides the
    unique parent behaviour using
    :class:`~sage.structure.unique_representation.UniqueRepresentation`. We
    show here how to make it possible to unpickle a pickle created with the
    factory, automatically turning it into an instance of the new class.

    First, we create the factory. In a doctest, it is needed to explicitly put
    it into ``__main__``, so that it can be located when pickling. Also, it is
    needed that we work with a new-style class::

        sage: from sage.structure.factory import UniqueFactory
        sage: import __main__
        sage: class OldStuff():
        ....:     def __init__(self, n, **extras):
        ....:         self.n = n
        ....:     def __repr__(self):
        ....:         return "Rotten old thing of level {}".format(self.n)
        sage: __main__.OldStuff = OldStuff
        sage: class MyFactory(UniqueFactory):
        ....:     def create_object(self, version, key, **extras):
        ....:         return OldStuff(key[0])
        ....:     def create_key(self, *args):
        ....:         return args
        sage: F = MyFactory('__main__.F')
        sage: __main__.F = F
        sage: a = F(3); a
        Rotten old thing of level 3
        sage: loads(dumps(a)) is a
        True

    Now, we create a pickle (the string returned by ``dumps(a)``)::

        sage: s = dumps(a)

    We create a new class, derived from
    :class:`~sage.structure.unique_representation.UniqueRepresentation`, that
    shall replace the old factory. In particular, the class has to have the
    same name as the old factory, and has to be put into the same module
    (here: ``__main__``). We turn it into a sub-class of the old class, but
    this is just to save the effort of writing a new init method::

        sage: from sage.structure.unique_representation import UniqueRepresentation
        sage: class F(UniqueRepresentation, OldStuff):
        ....:     def __repr__(self):
        ....:         return "Shiny new thing of level {}".format(self.n)
        sage: __main__.F = F

    The old pickle correctly unpickles as an instance of the new class, which
    is of course different from the instance of the old class, but exhibits
    unique object behaviour as well::

        sage: b = loads(s); b
        Shiny new thing of level 3
        sage: a is b
        False
        sage: loads(dumps(b)) is b
        True
    """
    cdef UniqueFactory F
    if factory is not None:
        try:
            F = factory
            return F.get_object(*args)
        except TypeError:
            pass
    # See Issue #16349: When replacing a UniqueFactory by something else (e.g.,
    # a UniqueRepresentation), then we get the object by calling.
    #
    # The first argument of a UniqueFactory pickle is a version number. We
    # strip this.
    return factory(*args[1], **args[2])


def generic_factory_reduce(self, proto):
    """
    Used to provide a ``__reduce__`` method if one does not already exist.

    EXAMPLES::

        sage: V = QQ^6                                                                  # needs sage.modules
        sage: sage.structure.factory.generic_factory_reduce(V, 1) == V.__reduce_ex__(1)             # needs sage.modules
        True
    """
    if self._factory_data is None:
        raise NotImplementedError("__reduce__ not implemented for %s" % type(self))
    else:
        return self._factory_data[0].reduce_data(self)


def lookup_global(name):
    """
    Used in unpickling the factory itself.

    EXAMPLES::

        sage: from sage.structure.factory import lookup_global
        sage: lookup_global('ZZ')
        Integer Ring
        sage: lookup_global('sage.rings.integer_ring.ZZ')
        Integer Ring
    """
    name = bytes_to_str(name, encoding='ASCII')
    try:
        return factory_unpickles[name]
    except KeyError:
        pass

    if '.' in name:
        module, name = name.rsplit('.', 1)
        all = __import__(module, fromlist=[name])
    else:
        try:
            import sage.all as all
        except ImportError:
            try:
                import sage.all__sagemath_modules as all
            except ImportError:
                try:
                    import sage.all__sagemath_pari as all
                except ImportError:
                    import sage.all__sagemath_categories as all
    return getattr(all, name)


# Old imports required for unpickling old pickles
from sage.structure.test_factory import test_factory
