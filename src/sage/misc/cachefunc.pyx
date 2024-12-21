# sage_setup: distribution = sagemath-objects
r"""
Cached Functions and Methods

AUTHORS:

- William Stein: initial version, (inspired by conversation with Justin Walker)
- Mike Hansen: added doctests and made it work with class methods.
- Willem Jan Palenstijn: add CachedMethodCaller for binding cached methods to
  instances.
- Tom Boothby: added DiskCachedFunction.
- Simon King: improved performance, more doctests, cython version,
  CachedMethodCallerNoArgs, weak cached function, cached special methods.
- Julian Rueth (2014-03-19, 2014-05-09, 2014-05-12): added ``key`` parameter, allow caching
  for unhashable elements, added ``do_pickle`` parameter

EXAMPLES:

By :issue:`11115`, cached functions and methods are now also
available in Cython code. The following examples cover various ways
of usage.

Python functions::

    sage: @cached_function
    ....: def test_pfunc(x):
    ....:     '''
    ....:     Some documentation
    ....:     '''
    ....:     return -x
    sage: test_pfunc(5) is test_pfunc(5)
    True

In some cases, one would only want to keep the result in cache as long
as there is any other reference to the result. By :issue:`12215`, this is
enabled for :class:`~sage.structure.unique_representation.UniqueRepresentation`,
which is used to create unique parents: If an algebraic structure, such
as a finite field, is only temporarily used, then it will not stay in
cache forever. That behaviour is implemented using ``weak_cached_function``,
that behaves the same as ``cached_function``, except that it uses a
:class:`~sage.misc.weak_dict.CachedWeakValueDictionary` for storing the results.

Cython cdef functions do not allow arbitrary decorators.
However, one can wrap a Cython function and turn it into
a cached function, by :issue:`11115`. We need to provide
the name that the wrapped method or function should have,
since otherwise the name of the original function would
be used::

    sage: # needs sage.misc.cython
    sage: cython('''cpdef test_funct(x): return -x''')
    sage: wrapped_funct = cached_function(test_funct, name='wrapped_funct')
    sage: wrapped_funct
    Cached version of <cyfunction test_funct at ...>
    sage: wrapped_funct.__name__
    'wrapped_funct'
    sage: wrapped_funct(5)
    -5
    sage: wrapped_funct(5) is wrapped_funct(5)
    True

We can proceed similarly for cached methods of Cython classes,
provided that they allow attribute assignment or have a public
attribute ``_cached_methods`` of type ``<dict>``. Since
:issue:`11115`, this is the case for all classes inheriting from
:class:`~sage.structure.parent.Parent`. See below for a more explicit
example. By :issue:`12951`, cached methods of extension classes can
be defined by simply using the decorator. However, an indirect
approach is still needed for cpdef methods::

    sage: # needs sage.misc.cython
    sage: cython_code = ['cpdef test_meth(self, x):',
    ....: '    "some doc for a wrapped cython method"',
    ....: '    return -x',
    ....: 'from sage.misc.cachefunc import cached_method',
    ....: 'from sage.structure.parent cimport Parent',
    ....: 'cdef class MyClass(Parent):',
    ....: '    @cached_method',
    ....: '    def direct_method(self, x):',
    ....: '        "Some doc for direct method"',
    ....: '        return 2*x',
    ....: '    wrapped_method = cached_method(test_meth,name="wrapped_method")']
    sage: cython(os.linesep.join(cython_code))
    sage: O = MyClass()
    sage: O.direct_method
    Cached version of <cyfunction MyClass.direct_method at ...>
    sage: O.wrapped_method
    Cached version of <cyfunction test_meth at ...>
    sage: O.wrapped_method.__name__
    'wrapped_method'
    sage: O.wrapped_method(5)
    -5
    sage: O.wrapped_method(5) is O.wrapped_method(5)
    True
    sage: O.direct_method(5)
    10
    sage: O.direct_method(5) is O.direct_method(5)
    True

In some cases, one would only want to keep the result in cache as long
as there is any other reference to the result. By :issue:`12215`, this is
enabled for :class:`~sage.structure.unique_representation.UniqueRepresentation`,
which is used to create unique parents: If an algebraic structure, such
as a finite field, is only temporarily used, then it will not stay in
cache forever. That behaviour is implemented using ``weak_cached_function``,
that behaves the same as ``cached_function``, except that it uses a
:class:`~sage.misc.weak_dict.CachedWeakValueDictionary` for storing the results.

By :issue:`11115`, even if a parent does not allow attribute
assignment, it can inherit a cached method from the parent class of a
category (previously, the cache would have been broken)::

    sage: cython_code = ["from sage.misc.cachefunc import cached_method",
    ....: "from sage.misc.cachefunc import cached_in_parent_method",
    ....: "from sage.categories.category import Category",
    ....: "from sage.categories.objects import Objects",
    ....: "class MyCategory(Category):",
    ....: "    @cached_method",
    ....: "    def super_categories(self):",
    ....: "        return [Objects()]",
    ....: "    class ElementMethods:",
    ....: "        @cached_method",
    ....: "        def element_cache_test(self):",
    ....: "            return -self",
    ....: "        @cached_in_parent_method",
    ....: "        def element_via_parent_test(self):",
    ....: "            return -self",
    ....: "    class ParentMethods:",
    ....: "        @cached_method",
    ....: "        def one(self):",
    ....: "            return self.element_class(self,1)",
    ....: "        @cached_method",
    ....: "        def invert(self, x):",
    ....: "            return -x"]
    sage: cython('\n'.join(cython_code))                                                # needs sage.misc.cython
    sage: C = MyCategory()                                                              # needs sage.misc.cython

In order to keep the memory footprint of elements small, it was
decided to not support the same freedom of using cached methods
for elements: If an instance of a class derived from
:class:`~sage.structure.element.Element` does not allow attribute
assignment, then a cached method inherited from the category of
its parent will break, as in the class ``MyBrokenElement`` below.

However, there is a class :class:`~sage.structure.element.ElementWithCachedMethod`
that has generally a slower attribute access, but fully supports
cached methods. We remark, however, that cached methods are
*much* faster if attribute access works. So, we expect that
:class:`~sage.structure.element.ElementWithCachedMethod` will
hardly by used.
::

    sage: # needs sage.misc.cython
    sage: cython_code = ["from sage.structure.element cimport Element, ElementWithCachedMethod", "from cpython.object cimport PyObject_RichCompare",
    ....: "cdef class MyBrokenElement(Element):",
    ....: "    cdef public object x",
    ....: "    def __init__(self, P, x):",
    ....: "        self.x=x",
    ....: "        Element.__init__(self,P)",
    ....: "    def __neg__(self):",
    ....: "        return MyBrokenElement(self.parent(),-self.x)",
    ....: "    def _repr_(self):",
    ....: "        return '<%s>'%self.x",
    ....: "    def __hash__(self):",
    ....: "        return hash(self.x)",
    ....: "    cpdef _richcmp_(left, right, int op):",
    ....: "        return PyObject_RichCompare(left.x, right.x, op)",
    ....: "    def raw_test(self):",
    ....: "        return -self",
    ....: "cdef class MyElement(ElementWithCachedMethod):",
    ....: "    cdef public object x",
    ....: "    def __init__(self, P, x):",
    ....: "        self.x=x",
    ....: "        ElementWithCachedMethod.__init__(self,P)",
    ....: "    def __neg__(self):",
    ....: "        return MyElement(self.parent(),-self.x)",
    ....: "    def _repr_(self):",
    ....: "        return '<%s>'%self.x",
    ....: "    def __hash__(self):",
    ....: "        return hash(self.x)",
    ....: "    cpdef _richcmp_(left, right, int op):",
    ....: "        return PyObject_RichCompare(left.x, right.x, op)",
    ....: "    def raw_test(self):",
    ....: "        return -self",
    ....: "from sage.structure.parent cimport Parent",
    ....: "cdef class MyParent(Parent):",
    ....: "    Element = MyElement"]
    sage: cython('\n'.join(cython_code))
    sage: P = MyParent(category=C)
    sage: ebroken = MyBrokenElement(P, 5)
    sage: e = MyElement(P, 5)

The cached methods inherited by the parent works::

    sage: # needs sage.misc.cython
    sage: P.one()
    <1>
    sage: P.one() is P.one()
    True
    sage: P.invert(e)
    <-5>
    sage: P.invert(e) is P.invert(e)
    True

The cached methods inherited by ``MyElement`` works::

    sage: # needs sage.misc.cython
    sage: e.element_cache_test()
    <-5>
    sage: e.element_cache_test() is e.element_cache_test()
    True
    sage: e.element_via_parent_test()
    <-5>
    sage: e.element_via_parent_test() is e.element_via_parent_test()
    True

The other element class can only inherit a ``cached_in_parent_method``, since
the cache is stored in the parent. In fact, equal elements share the cache,
even if they are of different types::

    sage: e == ebroken                                                                  # needs sage.misc.cython
    True
    sage: type(e) == type(ebroken)                                                      # needs sage.misc.cython
    False
    sage: ebroken.element_via_parent_test() is e.element_via_parent_test()              # needs sage.misc.cython
    True

However, the cache of the other inherited method breaks, although the method
as such works::

    sage: ebroken.element_cache_test()                                                  # needs sage.misc.cython
    <-5>
    sage: ebroken.element_cache_test() is ebroken.element_cache_test()                  # needs sage.misc.cython
    False

The cache can be emptied::

    sage: # needs sage.misc.cython
    sage: a = test_pfunc(5)
    sage: test_pfunc.clear_cache()
    sage: a is test_pfunc(5)
    False
    sage: a = P.one()
    sage: P.one.clear_cache()
    sage: a is P.one()
    False

Since ``e`` and ``ebroken`` share the cache, when we empty it for one element
it is empty for the other as well::

    sage: b = ebroken.element_via_parent_test()                                         # needs sage.misc.cython
    sage: e.element_via_parent_test.clear_cache()                                       # needs sage.misc.cython
    sage: b is ebroken.element_via_parent_test()                                        # needs sage.misc.cython
    False

Introspection works::

    sage: # needs sage.misc.cython
    sage: from sage.misc.edit_module import file_and_line
    sage: from sage.misc.sageinspect import sage_getdoc, sage_getfile, sage_getsource
    sage: print(sage_getdoc(test_pfunc))
       Some documentation
    sage: print(sage_getdoc(O.wrapped_method))
    some doc for a wrapped cython method
    <BLANKLINE>
    sage: print(sage_getdoc(O.direct_method))
    Some doc for direct method
    <BLANKLINE>
    sage: print(sage_getsource(O.wrapped_method))
    cpdef test_meth(self, x):
        "some doc for a wrapped cython method"
        return -x
    sage: print(sage_getsource(O.direct_method))
    @cached_method
    def direct_method(self, x):
        "Some doc for direct method"
        return 2*x

It is a very common special case to cache a method that has no
arguments. In that special case, the time needed to access the cache
can be drastically reduced by using a special implementation. The
cached method decorator automatically determines which implementation
ought to be chosen. A typical example is
:meth:`sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal.gens`
(no arguments) versus
:meth:`sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal.groebner_basis`
(several arguments)::

    sage: P.<a,b,c,d> = QQ[]
    sage: I = P * [a, b]
    sage: I.gens()
    [a, b]
    sage: I.gens() is I.gens()
    True
    sage: I.groebner_basis()                                                            # needs sage.libs.singular
    [a, b]
    sage: I.groebner_basis() is I.groebner_basis()                                      # needs sage.libs.singular
    True
    sage: type(I.gens)
    <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
    sage: type(I.groebner_basis)
    <class 'sage.misc.cachefunc.CachedMethodCaller'>

By :issue:`12951`, the cached_method decorator is also supported on non-c(p)def
methods of extension classes, as long as they either support attribute assignment
or have a public attribute of type ``<dict>`` called ``_cached_methods``. The
latter is easy::

    sage: # needs sage.misc.cython
    sage: cython_code = [
    ....: "from sage.misc.cachefunc import cached_method",
    ....: "cdef class MyClass:",
    ....: "    cdef public dict _cached_methods",
    ....: "    @cached_method",
    ....: "    def f(self, a, b):",
    ....: "        return a*b"]
    sage: cython(os.linesep.join(cython_code))
    sage: P = MyClass()
    sage: P.f(2, 3)
    6
    sage: P.f(2, 3) is P.f(2, 3)
    True

Providing attribute access is a bit more tricky, since it is needed that
an attribute inherited by the instance from its class can be overridden
on the instance. That is why providing a ``__getattr__`` would not be
enough in the following example::

    sage: # needs sage.misc.cython
    sage: cython_code = [
    ....: "from sage.misc.cachefunc import cached_method",
    ....: "cdef class MyOtherClass:",
    ....: "    cdef dict D",
    ....: "    def __init__(self):",
    ....: "        self.D = {}",
    ....: "    def __setattr__(self, n, v):",
    ....: "        self.D[n] = v",
    ....: "    def __getattribute__(self, n):",
    ....: "        try:",
    ....: "            return self.D[n]",
    ....: "        except KeyError:",
    ....: "            pass",
    ....: "        return getattr(type(self),n).__get__(self)",
    ....: "    @cached_method",
    ....: "    def f(self, a, b):",
    ....: "        return a+b"]
    sage: cython(os.linesep.join(cython_code))
    sage: Q = MyOtherClass()
    sage: Q.f(2, 3)
    5
    sage: Q.f(2, 3) is Q.f(2, 3)
    True

Note that supporting attribute access is somehow faster than the
easier method::

    sage: timeit("a = P.f(2,3)")   # random                                             # needs sage.misc.cython
    625 loops, best of 3: 1.3 µs per loop
    sage: timeit("a = Q.f(2,3)")   # random                                             # needs sage.misc.cython
    625 loops, best of 3: 931 ns per loop

Some immutable objects (such as `p`-adic numbers) cannot implement a
reasonable hash function because their ``==`` operator has been
modified to return ``True`` for objects which might behave differently
in some computations::

    sage: # needs sage.rings.padics
    sage: K.<a> = Qq(9)
    sage: b = a.add_bigoh(1)
    sage: c = a + 3
    sage: b
    a + O(3)
    sage: c
    a + 3 + O(3^20)
    sage: b == c
    True
    sage: b == a
    True
    sage: c == a
    False

If such objects defined a non-trivial hash function, this would break
caching in many places. However, such objects should still be usable
in caches. This can be achieved by defining an appropriate method
``_cache_key``::

    sage: # needs sage.rings.padics
    sage: hash(b)
    Traceback (most recent call last):
    ...
    TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'
    sage: from sage.misc.cachefunc import cached_method
    sage: @cached_method
    ....: def f(x): return x == a
    sage: f(b)
    True
    sage: f(c)  # if b and c were hashable, this would return True
    False
    sage: b._cache_key()
    (..., ((0, 1),), 0, 1)
    sage: c._cache_key()
    (..., ((0, 1), (1,)), 0, 20)

.. NOTE::

    This attribute will only be accessed if the object itself
    is not hashable.

An implementation must make sure that for elements ``a`` and ``b``,
if ``a != b``, then also ``a._cache_key() != b._cache_key()``.
In practice this means that the ``_cache_key`` should always include
the parent as its first argument::

    sage: S.<a> = Qq(4)                                                                 # needs sage.rings.padics
    sage: d = a.add_bigoh(1)                                                            # needs sage.rings.padics
    sage: b._cache_key() == d._cache_key()  # this would be True if the parents were not included
    False
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                          Mike Hansen <mhansen@gmail.com>
#                     2011 Simon King <simon.king@uni-jena.de>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#                     2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef extern from "methodobject.h":
    cdef int METH_NOARGS, METH_O
    cdef int PyCFunction_GetFlags(object op) except -1

import os
from sage.misc.sageinspect import sage_getfile_relative, sage_getsourcelines, sage_getargspec
from inspect import isfunction

from sage.misc.weak_dict cimport CachedWeakValueDictionary
from sage.misc.decorators import decorator_keywords

cdef frozenset special_method_names = frozenset(
    ['__abs__', '__add__',
     '__and__', '__call__', '__cmp__', '__coerce__', '__complex__', '__contains__', '__del__',
     '__delattr__', '__delete__', '__delitem__', '__delslice__', '__dir__', '__div__',
     '__eq__', '__float__', '__floordiv__', '__format__', '__ge__', '__get__', '__getattr__',
     '__getattribute__', '__getitem__', '__getslice__', '__gt__', '__hash__', '__hex__',
     '__iadd__', '__iand__', '__idiv__', '__ifloordiv__', '__ilshift__', '__imod__', '__imul__',
     '__index__', '__init__', '__instancecheck__', '__int__', '__invert__', '__ior__', '__ipow__',
     '__irshift__', '__isub__', '__iter__', '__itruediv__', '__ixor__', '__le__', '__len__',
     '__length_hint__', '__long__', '__lshift__', '__lt__', '__missing__', '__mod__', '__mul__',
     '__ne__', '__neg__', '__new__', '__oct__', '__or__', '__pos__', '__pow__',
     '__radd__', '__rand__', '__rdiv__', '__repr__', '__reversed__', '__rfloordiv__', '__rlshift__',
     '__rmod__', '__rmul__', '__ror__', '__rpow__', '__rrshift__', '__rshift__', '__rsub__',
     '__rtruediv__', '__rxor__', '__set__', '__setattr__', '__setitem__', '__setslice__', '__sizeof__',
     '__str__', '__sub__', '__subclasscheck__', '__truediv__', '__unicode__', '__xor__', 'next'])


def _cached_function_unpickle(module, name, cache=None):
    """
    Unpickle the cache function ``name`` defined in ``module``.

    This function loads ``name`` from ``module`` (it does not restore the code
    of the actual function when it was pickled.) The cache is restored from
    ``cache`` if present.

    INPUT:

    - ``module`` -- the name of the module to import the function from
    - ``name`` -- the name of the cached function
    - ``cache`` -- list of cached key value pairs

    TESTS::

        sage: type(hilbert_class_polynomial)                                            # needs sage.schemes
        <class 'sage.misc.cachefunc.CachedFunction'>
        sage: loads(dumps(hilbert_class_polynomial)) is hilbert_class_polynomial  #indirect doctest                     # needs sage.schemes
        True

    Verify that the ``cache`` parameter works::

        sage: @cached_function(do_pickle=True)
        ....: def f(n): return n
        sage: import __main__
        sage: __main__.f = f
        sage: f(0)
        0
        sage: ((0,),()) in f.cache
        True

        sage: s = dumps(f)
        sage: f.clear_cache()
        sage: ((0,),()) in f.cache
        False
        sage: f = loads(s)
        sage: ((0,),()) in f.cache
        True
        sage: f(0)
        0
    """
    ret = getattr(__import__(module, fromlist=['']), name)
    if cache is not None:
        ret.cache.update(cache)
    return ret


cdef class NonpicklingDict(dict):
    r"""
    A special dict which does not pickle its contents.

    EXAMPLES::

        sage: from sage.misc.cachefunc import NonpicklingDict
        sage: d = NonpicklingDict()
        sage: d[0] = 0
        sage: loads(dumps(d))
        {}
    """
    def __reduce__(self):
        r"""
        Return data required to unpickle this dictionary.

        EXAMPLES::

            sage: from sage.misc.cachefunc import NonpicklingDict
            sage: d = NonpicklingDict()
            sage: d[0] = 0
            sage: d.__reduce__()
            (<class 'sage.misc.cachefunc.NonpicklingDict'>, ())
        """
        return NonpicklingDict, ()


cdef unhashable_key = object()

cpdef inline dict_key(o):
    """
    Return a key to cache object ``o`` in a dict.

    This is different from ``cache_key`` since the ``cache_key`` might
    get confused with the key of a hashable object. Therefore, such keys
    include ``unhashable_key`` which acts as a unique marker which is
    certainly not stored in the dictionary otherwise.

    EXAMPLES::

        sage: from sage.misc.cachefunc import dict_key
        sage: dict_key(42)
        42
        sage: K.<u> = Qq(9)                                                             # needs sage.rings.padics
        sage: dict_key(u)                                                               # needs sage.rings.padics
        (<object object at ...>, (..., 20))
    """
    try:
        hash(o)
    except TypeError:
        o = (unhashable_key, cache_key_unhashable(o))
    return o


cpdef inline cache_key(o):
    r"""
    Helper function to return a hashable key for ``o`` which can be used for
    caching.

    This function is intended for objects which are not hashable such as
    `p`-adic numbers. The difference from calling an object's ``_cache_key``
    method directly, is that it also works for tuples and unpacks them
    recursively (if necessary, i.e., if they are not hashable).

    EXAMPLES::

        sage: from sage.misc.cachefunc import cache_key
        sage: K.<u> = Qq(9)                                                             # needs sage.rings.padics
        sage: a = K(1); a                                                               # needs sage.rings.padics
        1 + O(3^20)
        sage: cache_key(a)                                                              # needs sage.rings.padics
        (..., ((1,),), 0, 20)

    This function works if ``o`` is a tuple. In this case it unpacks its
    entries recursively::

        sage: o = (1, 2, (3, a))                                                        # needs sage.rings.padics
        sage: cache_key(o)                                                              # needs sage.rings.padics
        (1, 2, (3, (..., ((1,),), 0, 20)))

    Note that tuples are only partially unpacked if some of its entries are
    hashable::

        sage: o = (1/2, a)                                                              # needs sage.rings.padics
        sage: cache_key(o)                                                              # needs sage.rings.padics
        (1/2, (..., ((1,),), 0, 20))
    """
    try:
        hash(o)
    except TypeError:
        o = cache_key_unhashable(o)
    return o


cdef cache_key_unhashable(o):
    """
    Return a key for caching an item which is unhashable.
    """
    if isinstance(o, tuple):
        return tuple(cache_key(item) for item in o)
    try:
        k = o._cache_key()
    except AttributeError:
        raise TypeError("unhashable type: {!r}".format(type(o).__name__))
    return cache_key(k)


cdef class CachedFunction():
    """
    Create a cached version of a function, which only recomputes
    values it hasn't already computed. Synonyme: ``cached_function``

    INPUT:

    - ``f`` -- a function
    - ``name`` -- (optional string) name that the cached version
      of ``f`` should be provided with
    - ``key`` -- (optional callable) takes the input and returns a
      key for the cache, typically one would use this to normalize input
    - ``do_pickle`` -- (optional boolean) whether or not the contents of the
      cache should be included when pickling this function; the default is not
      to include them.

    If ``f`` is a function, do either ``g = CachedFunction(f)``
    or ``g = cached_function(f)`` to make a cached version of ``f``,
    or put ``@cached_function`` right before the definition of ``f``
    (i.e., use Python decorators)::

        @cached_function
        def f(...):
            ....

    The inputs to the function must be hashable or they must define
    :meth:`sage.structure.sage_object.SageObject._cache_key`.

    EXAMPLES::

        sage: @cached_function
        ....: def mul(x, y=2):
        ....:     return x*y
        sage: mul(3)
        6

    We demonstrate that the result is cached, and that, moreover,
    the cache takes into account the various ways of providing
    default arguments::

        sage: mul(3) is mul(3,2)
        True
        sage: mul(3,y=2) is mul(3,2)
        True

    The user can clear the cache::

        sage: a = mul(4)
        sage: mul.clear_cache()
        sage: a is mul(4)
        False

    It is also possible to explicitly override the cache with
    a different value::

        sage: mul.set_cache('foo',5)
        sage: mul(5,2)
        'foo'

    The parameter ``key`` can be used to ignore parameters for
    caching. In this example we ignore the parameter ``algorithm``::

        sage: @cached_function(key=lambda x,y,algorithm: (x,y))
        ....: def mul(x, y, algorithm='default'):
        ....:     return x*y
        sage: mul(1,1,algorithm='default') is mul(1,1,algorithm='algorithm') is mul(1,1) is mul(1,1,'default')
        True
    """
    def __init__(self, f, *, classmethod=False, name=None, key=None, do_pickle=None):
        """
        Create a cached version of a function, which only recomputes
        values it hasn't already computed. A custom name can be
        provided by an optional argument "name".

        If ``f`` is a function, do either ``g = CachedFunction(f)``
        to make a cached version of ``f``, or put ``@CachedFunction``
        right before the definition of ``f`` (i.e., use Python decorators)::

            @CachedFunction
            def f(...):
                ....

        The inputs to the function must be hashable or they must define
        :meth:`sage.structure.sage_object.SageObject._cache_key`.

        TESTS::

            sage: # needs sage.combinat
            sage: g = CachedFunction(number_of_partitions)
            sage: g.__name__
            'number_of_partitions'
            sage: 'partitions' in sage.misc.sageinspect.sage_getdoc(g)
            True
            sage: g(5)                                                                  # needs sage.libs.flint
            7
            sage: g.cache                                                               # needs sage.libs.flint
            {((5, 'default'), ()): 7}

            sage: def f(t=1): print(t)
            sage: h = CachedFunction(f)
            sage: w = walltime()
            sage: h(); h(1); h(t=1)
            1
            sage: walltime(w) < 2
            True

        By default, the contents of the cache are not pickled::

            sage: @cached_function
            ....: def f(n): return None
            sage: import __main__
            sage: __main__.f = f
            sage: for i in range(100): f(i)
            sage: len(f.cache)
            100

            sage: s = dumps(f)
            sage: f.clear_cache()
            sage: f = loads(s)
            sage: len(f.cache)
            0

        If ``do_pickle`` is set, then the cache is pickled::

            sage: @cached_function(do_pickle=True)
            ....: def f(n): return None
            sage: __main__.f = f
            sage: for i in range(100): f(i)
            sage: len(f.cache)
            100

            sage: s = dumps(f)
            sage: f.clear_cache()
            sage: f = loads(s)
            sage: len(f.cache)
            100
        """
        self.is_classmethod = classmethod
        self._common_init(f, None, name=name, key=key, do_pickle=do_pickle)
        self.cache = {} if do_pickle else NonpicklingDict()

    def _common_init(self, f, argument_fixer, name=None, key=None, do_pickle=None):
        """
        Perform initialization common to CachedFunction and CachedMethodCaller.

        TESTS::

            sage: @cached_function
            ....: def test_cache(x):
            ....:     return -x
            sage: test_cache.__name__  # indirect doctest
            'test_cache'
        """
        self.f = f
        self.key = key
        self.do_pickle = do_pickle
        if name is not None:
            self.__name__ = name
        else:
            self.__name__ = f.__name__
        try:
            self.__cached_module__ = f.__module__
        except AttributeError:
            self.__cached_module__ = f.__objclass__.__module__
        if argument_fixer is not None:
            # it is None unless the argument fixer
            # was known previously. See #15038.
            self._argument_fixer = argument_fixer

    @property
    def __module__(self):
        return self.__cached_module__

    cdef get_key_args_kwds(self, tuple args, dict kwds):
        """
        Return the key in the cache to be used when ``args`` and
        ``kwds`` are passed in as parameters.

        See ``get_key`` for the Python interface and tests.
        """
        # The key for "no arguments" is cached in empty_key
        if not args and not kwds:
            if self.empty_key is None:
                self.empty_key = self.fix_args_kwds(args, kwds)
            return self.empty_key

        return self.fix_args_kwds(args, kwds)

    cdef int argfix_init(self) except -1:
        """
        TESTS::

            sage: @cached_function
            ....: def test_cache(x):
            ....:     return -x
            sage: test_cache(1)
            -1
        """
        self._argument_fixer = ArgumentFixer(self.f,
                                             classmethod=self.is_classmethod)

    cdef fix_args_kwds(self, tuple args, dict kwds):
        r"""
        Normalize parameters to obtain a key for the cache.

        TESTS::

            sage: @cached_function(key=lambda x,y,algorithm: (x,y))
            ....: def mul(x, y, algorithm='default'):
            ....:     return x*y
            sage: mul.get_key(1,1,"default")  # indirect doctest
            (1, 1)
        """
        if self._argument_fixer is None:
            self.argfix_init()

        k = self._argument_fixer.fix_to_pos_args_kwds(args, kwds)
        if self.key is None:
            return k
        else:
            return self.key(*k[0], **dict(k[1]))

    def __reduce__(self):
        """
        Pickling of cached functions.

        TESTS::

            sage: type(hilbert_class_polynomial)                                        # needs sage.schemes
            <class 'sage.misc.cachefunc.CachedFunction'>
            sage: loads(dumps(hilbert_class_polynomial)) is hilbert_class_polynomial  #indirect doctest                 # needs sage.schemes
            True
        """
        return _cached_function_unpickle, (self.__cached_module__, self.__name__, self.cache)

    #########
    #  Introspection
    #
    # We provide some methods explicitly, and
    # forward other questions to the cached function.

    def _instancedoc_(self):
        """
        Provide documentation for the cached function.

        A cached function shall inherit the documentation
        from the function that is wrapped, not from the
        documentation of the wrapper.

        TESTS::

            sage: P.<x,y> = QQ[]
            sage: I = P * [x,y]
            sage: from sage.misc.sageinspect import sage_getdoc
            sage: print(sage_getdoc(I.groebner_basis))  # indirect doctest
            WARNING: the enclosing module is marked...
               Return the reduced Groebner basis of this ideal.
            ...

        Test that :issue:`15184` is fixed::

            sage: from sage.misc.sageinspect import sage_getfile
            sage: type(I.groebner_basis)
            <class 'sage.misc.cachefunc.CachedMethodCaller'>
            sage: os.path.exists(sage_getfile(I.groebner_basis))
            True

        Test that :issue:`18064` is fixed::

            sage: @cached_function
            ....: def f():
            ....:     return 3
            sage: f.__doc__
            'File: ... (starting at line 1)\n'
        """
        from sage.misc.sageinspect import _extract_embedded_position
        f = self.f
        doc = f.__doc__ or ''
        if not doc or _extract_embedded_position(doc) is None:
            try:
                sourcelines = sage_getsourcelines(f)
                filename = sage_getfile_relative(f)
                # this is a rather expensive way of getting the line number, because
                # retrieving the source requires reading the source file and in many
                # cases this is not required (in cython it's embedded in the docstring,
                # on code objects you'll find it in co_filename and co_firstlineno)
                # however, this hasn't been factored out yet in sageinspect
                # and the logic in sage_getsourcelines is rather intricate.
                file_info = "File: {} (starting at line {})".format(filename, sourcelines[1]) + os.linesep

                doc = file_info + doc
            except IOError:
                pass
        return doc

    def _sage_src_(self):
        """
        Return the source code for the wrapped function.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getsource
            sage: g = CachedFunction(number_of_partitions)                              # needs sage.combinat
            sage: 'flint' in sage_getsource(g)  # indirect doctest                      # needs sage.combinat
            True
        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.f)

    def _sage_src_lines_(self):
        r"""
        Return the list of source lines and the first line number
        of the wrapped function.

        TESTS::

            sage: P.<x,y> = QQ[]
            sage: I = P*[x,y]
            sage: from sage.misc.sageinspect import sage_getsourcelines
            sage: l = '        elif algorithm.startswith("macaulay2:"):\n'
            sage: l in sage_getsourcelines(I.groebner_basis)[0]  # indirect doctest
            True
        """
        from sage.misc.sageinspect import sage_getsourcelines
        return sage_getsourcelines(self.f)

    def _sage_argspec_(self):
        """
        Return the argspec of the wrapped function or method.

        This was implemented in :issue:`11115`.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: I = P*[x,y]
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(I.groebner_basis)   # indirect doctest
            FullArgSpec(args=['self', 'algorithm', 'deg_bound', 'mult_bound', 'prot'],
                        varargs='args', varkw='kwds', defaults=('', None, None, False),
                        kwonlyargs=[], kwonlydefaults=None, annotations={})
        """
        return sage_getargspec(self.f)

    def __call__(self, *args, **kwds):
        """
        Return value from cache or call the wrapped function,
        caching the output.

        TESTS::

            sage: # needs sage.combinat sage.libs.flint
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.cache
            {((5, 'default'), ()): 7}
            sage: a = g(10^5)   # indirect doctest
            sage: a == number_of_partitions(10^5)
            True
            sage: a is g(10^5)
            True
            sage: a is number_of_partitions(10^5)
            True

        Check that :issue:`16316` has been fixed, i.e., caching works for
        immutable unhashable objects which define
        :meth:`sage.structure.sage_object.SageObject._cache_key`::

            sage: # needs sage.rings.padics
            sage: @cached_function
            ....: def f(x): return x+x
            sage: K.<u> = Qq(4)
            sage: x = K(1,1); x
            1 + O(2)
            sage: y = K(1,2); y
            1 + O(2^2)
            sage: x == y
            True
            sage: f(x) is f(x)
            True
            sage: f(y) is not f(x)
            True
        """
        k = self.get_key_args_kwds(args, kwds)

        try:
            try:
                return self.cache[k]
            except TypeError:  # k is not hashable
                k = dict_key(k)
                return self.cache[k]
        except KeyError:
            w = self.f(*args, **kwds)
            self.cache[k] = w
            return w

    def cached(self, *args, **kwds):
        """
        Return the result from the cache if available. If the value is
        not cached, raise :exc:`KeyError`.

        EXAMPLES::

            sage: @cached_function
            ....: def f(x):
            ....:     return x
            sage: f.cached(5)
            Traceback (most recent call last):
            ...
            KeyError: ((5,), ())
            sage: f(5)
            5
            sage: f.cached(5)
            5
        """
        k = self.get_key_args_kwds(args, kwds)

        try:
            return self.cache[k]
        except TypeError:  # k is not hashable
            k = dict_key(k)
            return self.cache[k]

    def is_in_cache(self, *args, **kwds):
        """
        Check if the argument list is in the cache.

        EXAMPLES::

            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, z, y=0):
            ....:         return self._x*z+y
            sage: a = Foo(2)
            sage: a.f.is_in_cache(3)
            False
            sage: a.f(3)
            6
            sage: a.f.is_in_cache(3,y=0)
            True

        TESTS:

        Check that :issue:`16316` has been fixed, i.e., caching works for
        immutable unhashable objects which define
        :meth:`sage.structure.sage_object.SageObject._cache_key`::

            sage: # needs sage.rings.padics
            sage: @cached_function
            ....: def f(x): return x
            sage: K.<u> = Qq(4)
            sage: x = K(1,1); x
            1 + O(2)
            sage: f.is_in_cache(x)
            False
            sage: f(x)
            1 + O(2)
            sage: f.is_in_cache(x)
            True
        """
        k = self.get_key_args_kwds(args, kwds)
        try:
            return k in self.cache
        except TypeError:  # k is not hashable
            k = dict_key(k)
            return k in self.cache

    def set_cache(self, value, *args, **kwds):
        """
        Set the value for those args and keyword args
        Mind the unintuitive syntax (value first).
        Any idea on how to improve that welcome!

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.flint
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.cache
            {((5, 'default'), ()): 7}
            sage: g.set_cache(17, 5)
            sage: g.cache
            {((5, 'default'), ()): 17}
            sage: g(5)
            17

        TESTS:

        Check that :issue:`16316` has been fixed, i.e., caching works for
        immutable unhashable objects which define
        :meth:`sage.structure.sage_object.SageObject._cache_key`::

            sage: # needs sage.rings.padics
            sage: @cached_function
            ....: def f(x): return x
            sage: K.<u> = Qq(4)
            sage: x = K(1,1); x
            1 + O(2)
            sage: f.set_cache(x, x)
            sage: f.is_in_cache(x)
            True

        DEVELOPER NOTE:

        Is there a way to use the following intuitive syntax?

        ::

            sage: g(5) = 19    # todo: not implemented
            sage: g(5)         # todo: not implemented
            19
        """
        k = self.get_key_args_kwds(args, kwds)
        try:
            self.cache[k] = value
        except TypeError:  # k is not hashable
            k = dict_key(k)
            self.cache[k] = value

    def get_key(self, *args, **kwds):
        """
        Return the key in the cache to be used when ``args``
        and ``kwds`` are passed in as parameters.

        EXAMPLES::

            sage: @cached_function
            ....: def foo(x):
            ....:     return x^2
            sage: foo(2)
            4
            sage: foo.get_key(2)
            ((2,), ())
            sage: foo.get_key(x=3)
            ((3,), ())

        Examples for cached methods::

            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, y, z=0):
            ....:         return self._x * y + z
            sage: a = Foo(2)
            sage: z = a.f(37)
            sage: k = a.f.get_key(37); k
            ((37, 0), ())
            sage: a.f.cache[k] is z
            True

        Note that the method does not test whether there are
        too many arguments, or wrong argument names::

            sage: a.f.get_key(1,2,3,x=4,y=5,z=6)
            ((1, 2, 3), (('x', 4), ('y', 5), ('z', 6)))

        It does, however, take into account the different
        ways of providing named arguments, possibly with a
        default value::

            sage: a.f.get_key(5)
            ((5, 0), ())
            sage: a.f.get_key(y=5)
            ((5, 0), ())
            sage: a.f.get_key(5,0)
            ((5, 0), ())
            sage: a.f.get_key(5,z=0)
            ((5, 0), ())
            sage: a.f.get_key(y=5,z=0)
            ((5, 0), ())
        """
        return self.get_key_args_kwds(args, kwds)

    def __repr__(self):
        """
        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)                              # needs sage.combinat
            sage: g     # indirect doctest                                              # needs sage.combinat
            Cached version of <function number_of_partitions at 0x...>
        """
        try:
            return "Cached version of {}".format(self.f)
        except AttributeError:
            return "Cached version of a method (pending reassignment)"

    def clear_cache(self):
        """
        Clear the cache dictionary.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)                                                              # needs sage.libs.flint
            sage: g.cache                                                               # needs sage.libs.flint
            {((5, 'default'), ()): 7}
            sage: g.clear_cache()
            sage: g.cache
            {}
        """
        self.cache.clear()

    def precompute(self, arglist, num_processes=1):
        """
        Cache values for a number of inputs.  Do the computation
        in parallel, and only bother to compute values that we
        haven't already cached.

        INPUT:

        - ``arglist`` -- list (or iterables) of arguments for which
          the method shall be precomputed

        - ``num_processes`` -- number of processes used by
          :func:`~sage.parallel.decorate.parallel`

        EXAMPLES::

            sage: @cached_function
            ....: def oddprime_factors(n):
            ....:     l = [p for p,e in factor(n) if p != 2]
            ....:     return len(l)
            sage: oddprime_factors.precompute(range(1,100), 4)
            sage: oddprime_factors.cache[(25,),()]
            1
        """
        from sage.parallel.decorate import parallel, normalize_input
        P = parallel(num_processes)(self.f)
        cdef list arglist2 = []
        for a in arglist:
            ak = normalize_input(a)
            if self.get_key_args_kwds(ak[0], ak[1]) not in self.cache:
                arglist2.append(ak)
        for ((args, kwargs), val) in P(arglist2):
            self.set_cache(val, *args, **kwargs)


cached_function = decorator_keywords(CachedFunction)


cdef class WeakCachedFunction(CachedFunction):
    """
    A version of :class:`CachedFunction` using weak references on the
    values.

    If ``f`` is a function, do either ``g = weak_cached_function(f)`` to make
    a cached version of ``f``, or put ``@weak_cached_function`` right before
    the definition of ``f`` (i.e., use Python decorators)::

        @weak_cached_function
        def f(...):
            ...

    As an exception meant to improve performance, the most recently
    computed values are strongly referenced. The number of strongly
    cached values can be controlled by the ``cache`` keyword.

    EXAMPLES::

        sage: from sage.misc.cachefunc import weak_cached_function
        sage: class A: pass
        sage: @weak_cached_function(cache=0)
        ....: def f():
        ....:     print("doing a computation")
        ....:     return A()
        sage: a = f()
        doing a computation

    The result is cached::

        sage: b = f()
        sage: a is b
        True

    However, if there are no strong references left, the result is
    deleted, and thus a new computation takes place::

        sage: del a
        sage: del b
        sage: a = f()
        doing a computation

    Above, we used the ``cache=0`` keyword. With a larger value, the
    most recently computed values are cached anyway, even if they are
    not referenced::

        sage: @weak_cached_function(cache=3)
        ....: def f(x):
        ....:     print("doing a computation for x={}".format(x))
        ....:     return A()
        sage: a = f(1); del a
        doing a computation for x=1
        sage: a = f(2), f(1); del a
        doing a computation for x=2
        sage: a = f(3), f(1); del a
        doing a computation for x=3
        sage: a = f(4), f(1); del a
        doing a computation for x=4
        doing a computation for x=1
        sage: a = f(5), f(1); del a
        doing a computation for x=5

    The parameter ``key`` can be used to ignore parameters for
    caching. In this example we ignore the parameter ``algorithm``::

        sage: @weak_cached_function(key=lambda x,algorithm: x)
        ....: def mod_ring(x, algorithm='default'):
        ....:     return IntegerModRing(x)
        sage: mod_ring(1,algorithm='default') is mod_ring(1,algorithm='algorithm') is mod_ring(1) is mod_ring(1,'default')
        True

    TESTS:

    Check that :issue:`16316` has been fixed, i.e., caching works for
    immutable unhashable objects which define
    :meth:`sage.structure.sage_object.SageObject._cache_key`::

        sage: # needs sage.rings.padics
        sage: from sage.misc.cachefunc import weak_cached_function
        sage: @weak_cached_function
        ....: def f(x): return x+x
        sage: K.<u> = Qq(4)
        sage: R.<t> = K[]
        sage: x = t + K(1,1); x
        (1 + O(2^20))*t + 1 + O(2)
        sage: y = t + K(1,2); y
        (1 + O(2^20))*t + 1 + O(2^2)
        sage: x == y
        True
        sage: f(x) is f(x)
        True
        sage: f(y) is not f(x)
        True

    Examples and tests for ``is_in_cache``::

        sage: from sage.misc.cachefunc import weak_cached_function
        sage: class A:
        ....:     def __init__(self, x):
        ....:         self.x = x
        sage: @weak_cached_function(cache=0)
        ....: def f(n):
        ....:    return A(n)
        sage: a = f(5)

    The key 5 is in the cache, as long as there is a strong
    reference to the corresponding value::

        sage: f.is_in_cache(5)
        True

    However, if there are no strong references left, the cached
    item is removed from the cache::

        sage: del a
        sage: f.is_in_cache(5)
        False

    Check that :issue:`16316` has been fixed, i.e., caching works for
    immutable unhashable objects which define
    :meth:`sage.structure.sage_object.SageObject._cache_key`::

        sage: # needs sage.rings.padics
        sage: from sage.misc.cachefunc import weak_cached_function
        sage: @weak_cached_function
        ....: def f(x): return x
        sage: K.<u> = Qq(4)
        sage: R.<t> = K[]
        sage: f.is_in_cache(t)
        False
        sage: f(t)
        (1 + O(2^20))*t
        sage: f.is_in_cache(t)
        True

    Examples and tests for ``set_cache``::

        sage: from sage.misc.cachefunc import weak_cached_function
        sage: @weak_cached_function
        ....: def f(n):
        ....:     raise RuntimeError
        sage: f.set_cache(ZZ, 5)
        sage: f(5)
        Integer Ring

    Check that :issue:`16316` has been fixed, i.e., caching works for
    immutable unhashable objects which define
    :meth:`sage.structure.sage_object.SageObject._cache_key`::

        sage: # needs sage.rings.padics
        sage: from sage.misc.cachefunc import weak_cached_function
        sage: @weak_cached_function
        ....: def f(x): return x
        sage: K.<u> = Qq(4)
        sage: R.<t> = K[]
        sage: f.set_cache(t,t)
        sage: f.is_in_cache(t)
        True
    """
    def __init__(self, f, *, classmethod=False, name=None, key=None, **kwds):
        """
        The inputs to the function must be hashable or they must define
        :meth:`sage.structure.sage_object.SageObject._cache_key`.
        The outputs to the function must be weakly referenceable.

        TESTS::

            sage: from sage.misc.cachefunc import weak_cached_function
            sage: class A: pass
            sage: @weak_cached_function
            ....: def f():
            ....:     return A()
            sage: f
            Cached version of <function f at ...>

        We demonstrate that pickling works, provided the uncached function
        is available::

            sage: import __main__
            sage: __main__.f = f
            sage: loads(dumps(f))
            Cached version of <function f at ...>
            sage: str(f.cache)
            '<CachedWeakValueDictionary at 0x...>'
        """
        self._common_init(f, None, name=name, key=key)
        self.cache = CachedWeakValueDictionary(**kwds)


weak_cached_function = decorator_keywords(WeakCachedFunction)


class CachedMethodPickle():
    """
    This class helps to unpickle cached methods.

    .. NOTE::

        Since :issue:`8611`, a cached method is an attribute
        of the instance (provided that it has a ``__dict__``).
        Hence, when pickling the instance, it would be attempted
        to pickle that attribute as well, but this is a problem,
        since functions cannot be pickled, currently. Therefore,
        we replace the actual cached method by a place holder,
        that kills itself as soon as any attribute is requested.
        Then, the original cached attribute is reinstated. But the
        cached values are in fact saved (if ``do_pickle`` is set.)

    EXAMPLES::

        sage: R.<x, y, z> = PolynomialRing(QQ, 3)
        sage: I = R * (x^3 + y^3 + z^3, x^4 - y^4)
        sage: I.groebner_basis()                                                        # needs sage.libs.singular
        [y^5*z^3 - 1/4*x^2*z^6 + 1/2*x*y*z^6 + 1/4*y^2*z^6,
         x^2*y*z^3 - x*y^2*z^3 + 2*y^3*z^3 + z^6,
         x*y^3 + y^4 + x*z^3, x^3 + y^3 + z^3]
        sage: I.groebner_basis
        Cached version of <function ...groebner_basis at 0x...>

    We now pickle and unpickle the ideal. The cached method
    ``groebner_basis`` is replaced by a placeholder::

        sage: J = loads(dumps(I))
        sage: J.groebner_basis
        Pickle of the cached method "groebner_basis"

    But as soon as any other attribute is requested from the
    placeholder, it replaces itself by the cached method, and
    the entries of the cache are actually preserved::

        sage: J.groebner_basis.is_in_cache()                                            # needs sage.libs.singular
        True
        sage: J.groebner_basis
        Cached version of <function ...groebner_basis at 0x...>
        sage: J.groebner_basis() == I.groebner_basis()                                  # needs sage.libs.singular
        True

    TESTS:

    Since :issue:`11115`, there is a special implementation for
    cached methods that don't take arguments::

        sage: class A:
        ....:     @cached_method(do_pickle=True)
        ....:     def f(self): return 1
        ....:     @cached_method(do_pickle=True)
        ....:     def g(self, x): return x

        sage: import __main__
        sage: __main__.A = A
        sage: a = A()
        sage: type(a.f)
        <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
        sage: type(a.g)
        <class 'sage.misc.cachefunc.CachedMethodCaller'>

    We demonstrate that both implementations can be pickled and
    preserve the cache. For that purpose, we assign nonsense to the
    cache. Of course, it is a very bad idea to override the cache in
    that way.  So, please don't try this at home::

        sage: a.f.set_cache(0)
        sage: a.f()
        0
        sage: a.g.set_cache(0,x=1)
        sage: a.g(1)
        0
        sage: b = loads(dumps(a))
        sage: b.f()
        0
        sage: b.g(1)
        0

    Anyway, the cache will be automatically reconstructed after
    clearing it::

        sage: a.f.clear_cache()
        sage: a.f()
        1

        sage: a.g.clear_cache()
        sage: a.g(1)
        1

    AUTHOR:

    - Simon King (2011-01)
    """
    def __init__(self, inst, name, cache=None):
        """
        INPUT:

        - ``inst`` -- some instance
        - ``name`` -- string; usually the name of an attribute
          of ``inst`` to which ``self`` is assigned

        TESTS::

            sage: from sage.misc.cachefunc import CachedMethodPickle
            sage: P = CachedMethodPickle(1, 'foo')
            sage: P
            Pickle of the cached method "foo"
        """
        self._instance = inst
        self._name = name
        self._cache = cache

    def __repr__(self):
        """
        TESTS::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R * (x^3 + y^3 + z^3, x^4 - y^4)
            sage: G = I.groebner_basis()                                                # needs sage.libs.singular
            sage: J = loads(dumps(I))
            sage: J.groebner_basis  # indirect doctest
            Pickle of the cached method "groebner_basis"
        """
        return 'Pickle of the cached method "{}"'.format(self._name)

    def __reduce__(self):
        """
        This class is a pickle. However, sometimes, pickles
        need to be pickled another time.

        TESTS::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R * (x^3 + y^3 + z^3, x^4 - y^4)
            sage: I.groebner_basis()                                                    # needs sage.libs.singular
            [y^5*z^3 - 1/4*x^2*z^6 + 1/2*x*y*z^6 + 1/4*y^2*z^6,
             x^2*y*z^3 - x*y^2*z^3 + 2*y^3*z^3 + z^6,
             x*y^3 + y^4 + x*z^3, x^3 + y^3 + z^3]
            sage: J = loads(dumps(I))
            sage: J.groebner_basis
            Pickle of the cached method "groebner_basis"

        When we now pickle ``J``, the pickle of the cached method
        needs to be taken care of::

            sage: K = loads(dumps(J))  # indirect doctest
            sage: K.groebner_basis
            Pickle of the cached method "groebner_basis"
            sage: K.groebner_basis.cache                                                # needs sage.libs.singular
            {(('', None, None, False), ()):
            [y^5*z^3 - 1/4*x^2*z^6 + 1/2*x*y*z^6 + 1/4*y^2*z^6,
             x^2*y*z^3 - x*y^2*z^3 + 2*y^3*z^3 + z^6,
             x*y^3 + y^4 + x*z^3, x^3 + y^3 + z^3]}
        """
        return CachedMethodPickle, (self._instance, self._name, self._cache)

    def __call__(self, *args, **kwds):
        """
        The purpose of this call method is to kill ``self`` and to
        replace it by an actual :class:`CachedMethodCaller`. The last
        thing that ``self`` does before disappearing is to call the
        :class:`CachedMethodCaller` and return the result.

        EXAMPLES::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens
            Cached version of <function ...gens at 0x...>
            sage: J = loads(dumps(I))
            sage: J.gens
            Pickle of the cached method "gens"
            sage: J.gens()   # indirect doctest
            [a, b]
            sage: J.gens
            Cached version of <function ...gens at 0x...>
        """
        self._instance.__dict__.__delitem__(self._name)
        CM = getattr(self._instance, self._name)
        if self._cache is not None:
            if isinstance(CM, CachedMethodCallerNoArgs):
                CM.cache = self._cache
            else:
                for k, v in self._cache:
                    CM.cache[k] = v
        return CM(*args, **kwds)

    def __getattr__(self, s):
        """
        TESTS::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R * (x^3 + y^3 + z^3, x^4 - y^4)
            sage: G = I.groebner_basis()                                                # needs sage.libs.singular
            sage: J = loads(dumps(I))
            sage: J.groebner_basis
            Pickle of the cached method "groebner_basis"

        If an attribute of name ``s`` is requested (say,
        ``is_in_cache``), the attribute ``self._name`` of
        ``self._instance`` is deleted. Then, the attribute
        of name ``s`` of the attribute ``self._name`` of
        ``self._instance`` is requested. Since ``self._name``
        is a cached method defined for the class of
        ``self._instance``, retrieving the just-deleted
        attribute ``self._name`` succeeds.

        In that way, the unpickling of the cached method is
        finally accomplished::

            sage: J.groebner_basis.is_in_cache()  # indirect doctest                    # needs sage.libs.singular
            True
            sage: J.groebner_basis
            Cached version of <function ...groebner_basis at 0x...>
        """
        self._instance.__dict__.__delitem__(self._name)
        CM = getattr(self._instance, self._name)
        if self._cache is not None:
            if isinstance(CM, CachedMethodCallerNoArgs):
                CM.cache = self._cache
            else:
                for k, v in self._cache:
                    CM.cache[k] = v
        return getattr(CM, s)


cdef object _COMPUTING = object()


cdef class CachedMethodCaller(CachedFunction):
    """
    Utility class that is used by :class:`CachedMethod` to bind a
    cached method to an instance.

    .. NOTE::

        Since :issue:`11115`, there is a special implementation
        :class:`CachedMethodCallerNoArgs` for methods that do not take
        arguments.

    EXAMPLES::

        sage: class A:
        ....:    @cached_method
        ....:    def bar(self, x):
        ....:        return x^2
        sage: a = A()
        sage: a.bar
        Cached version of <function ...bar at 0x...>
        sage: type(a.bar)
        <class 'sage.misc.cachefunc.CachedMethodCaller'>
        sage: a.bar(2) is a.bar(x=2)
        True

    TESTS:

    As of :issue:`15692` the contents of the cache are not pickled anymore::

        sage: import __main__
        sage: __main__.A = A
        sage: len(a.bar.cache)
        1
        sage: b = loads(dumps(a))
        sage: len(b.bar.cache)
        0

    The parameter ``do_pickle`` can be used to change this behaviour::

        sage: class A:
        ....:    @cached_method(do_pickle=True)
        ....:    def bar(self, x):
        ....:        return x^2

        sage: __main__.A = A
        sage: a = A()
        sage: a.bar(2)
        4
        sage: len(a.bar.cache)
        1
        sage: b = loads(dumps(a))
        sage: len(b.bar.cache)
        1
    """
    def __init__(self, CachedMethod cachedmethod, inst, *, cache=None, name=None, key=None, do_pickle=None):
        """
        EXAMPLES::

            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, *args):
            ....:         return self._x^2
            sage: a = Foo(2)
            sage: a.f.cache
            {}
            sage: a.f()
            4
            sage: a.f.cache
            {((), ()): 4}
        """
        # initialize CachedFunction. Since the cached method is actually bound
        # to an instance, it now makes sense to initialise the ArgumentFixer
        # and reuse it for all bound cached method callers of the unbound
        # cached method.
        if cachedmethod._cachedfunc._argument_fixer is None:
            cachedmethod._cachedfunc.argfix_init()
        self._common_init(cachedmethod._cachedfunc.f,
                          cachedmethod._cachedfunc._argument_fixer,
                          name=name,
                          key=key,
                          do_pickle=do_pickle)
        if cache is None:
            self.cache = NonpicklingDict() if do_pickle else {}
        else:
            self.cache = cache
        self._instance = inst
        self._cachedmethod = cachedmethod

    def __reduce__(self):
        """
        The pickle of a :class:`CachedMethodCaller` unpickles
        to a :class:`CachedMethodPickle`, that is able to replace
        itself by a copy of the original :class:`CachedMethodCaller`.

        TESTS::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R * (x^3 + y^3 + z^3, x^4 - y^4)
            sage: G = I.groebner_basis()                                                # needs sage.libs.singular
            sage: J = loads(dumps(I))  # indirect doctest
            sage: J.groebner_basis
            Pickle of the cached method "groebner_basis"
            sage: J.groebner_basis.is_in_cache()                                        # needs sage.libs.singular
            True
            sage: J.groebner_basis
            Cached version of <function ...groebner_basis at 0x...>
        """
        if isinstance(self._cachedmethod, CachedInParentMethod) or hasattr(self._instance, self._cachedmethod._cache_name):
            return CachedMethodPickle, (self._instance, self.__name__)
        else:
            return CachedMethodPickle, (self._instance, self.__name__, self.cache)

    def _instance_call(self, *args, **kwds):
        """
        Call the cached method without using the cache.

        EXAMPLES::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P * [a, b]
            sage: I.groebner_basis()                                                    # needs sage.libs.singular
            [a, b]
            sage: I.groebner_basis._instance_call() is I.groebner_basis()               # needs sage.libs.singular
            False
            sage: I.groebner_basis._instance_call() == I.groebner_basis()               # needs sage.libs.singular
            True

        ::

            sage: class Foo():
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, n=2):
            ....:         return self._x^n
            sage: a = Foo(2)
            sage: a.f()
            4

        Usually, a cached method is indeed cached::

            sage: a.f() is a.f()
            True

        However, when it becomes necessary, one can call it without
        using the cache::

            sage: a.f._instance_call() is a.f()
            False
            sage: a.f._instance_call() == a.f()
            True
        """
        return self.f(self._instance, *args, **kwds)

    cdef fix_args_kwds(self, tuple args, dict kwds):
        r"""
        Normalize parameters to obtain a key for the cache.

        TESTS::

            sage: class A():
            ....:     def _f_normalize(self, x, algorithm): return x
            ....:     @cached_method(key=_f_normalize)
            ....:     def f(self, x, algorithm='default'): return x
            sage: a = A()
            sage: a.f(1, algorithm='default') is a.f(1) is a.f(1, algorithm='algorithm')
            True
        """
        if self._argument_fixer is None:
            self.argfix_init()

        k = self._argument_fixer.fix_to_pos_args_kwds(args, kwds)
        if self.key is None:
            return k
        else:
            return self.key(self._instance, *k[0], **dict(k[1]))

    def __call__(self, *args, **kwds):
        """
        Call the cached method.

        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class Foo:
            ....:     @cached_method
            ....:     def f(self, x, y=1):
            ....:         return x+y
            ....:     g = deprecated_function_alias(57, f)
            sage: a = Foo()
            sage: a.f(1)  #indirect doctest
            2

        The result is cached, taking into account
        the three ways of providing (named) arguments::

            sage: a.f(5) is a.f(5,1)
            True
            sage: a.f(5) is a.f(5,y=1)
            True
            sage: a.f(5) is a.f(y=1,x=5)
            True

        The method can be called as a bound function using the same cache::

            sage: a.f(5) is Foo.f(a, 5)
            True
            sage: a.f(5) is Foo.f(a,5,1)
            True
            sage: a.f(5) is Foo.f(a, 5,y=1)
            True
            sage: a.f(5) is Foo.f(a, y=1,x=5)
            True

        Cached methods are compatible with
        :meth:`sage.misc.superseded.deprecated_function_alias`::

            sage: a.g(5) is a.f(5)
            doctest:...: DeprecationWarning: g is deprecated. Please use f instead.
            See https://github.com/sagemath/sage/issues/57 for details.
            True
            sage: Foo.g(a, 5) is a.f(5)
            True
            sage: Foo.g(a, y=1,x=5) is a.f(5)
            True

        We test that :issue:`5843` is fixed::

            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, y):
            ....:         return self._x
            sage: a = Foo(2)
            sage: b = Foo(3)
            sage: a.f(b.f)
            2

        Check that :issue:`16316` has been fixed, i.e., caching works for
        immutable unhashable objects which define
        :meth:`sage.structure.sage_object.SageObject._cache_key`::

            sage: # needs sage.rings.padics
            sage: K.<u> = Qq(4)
            sage: class A():
            ....:   @cached_method
            ....:   def f(self, x): return x+x
            sage: a = A()
            sage: x = K(1,1); x
            1 + O(2)
            sage: y = K(1,2); y
            1 + O(2^2)
            sage: x == y
            True
            sage: a.f(x) is a.f(x)
            True
            sage: a.f(y) is not a.f(x)
            True
        """
        if self._instance is None:
            # cached method bound to a class
            instance = args[0]
            args = args[1:]
            return self._cachedmethod.__get__(instance)(*args, **kwds)

        k = self.get_key_args_kwds(args, kwds)

        cdef dict cache = <dict>self.cache
        try:
            try:
                return cache[k]
            except TypeError:  # k is not hashable
                k = dict_key(k)
                return cache[k]
        except KeyError:
            w = self._instance_call(*args, **kwds)
            cache[k] = w
            return w

    def cached(self, *args, **kwds):
        """
        Return the result from the cache if available. If the value is
        not cached, raise :exc:`KeyError`.

        EXAMPLES::

            sage: class CachedMethodTest():
            ....:     @cached_method
            ....:     def f(self, x):
            ....:         return x
            sage: o = CachedMethodTest()
            sage: CachedMethodTest.f.cached(o, 5)
            Traceback (most recent call last):
            ...
            KeyError: ((5,), ())
            sage: o.f.cached(5)
            Traceback (most recent call last):
            ...
            KeyError: ((5,), ())
            sage: o.f(5)
            5
            sage: CachedMethodTest.f.cached(o, 5)
            5
            sage: o.f.cached(5)
            5
        """
        if self._instance is None:
            # cached method bound to a class
            instance = args[0]
            args = args[1:]
            return self._cachedmethod.__get__(instance).cached(*args, **kwds)

        k = self.get_key_args_kwds(args, kwds)

        try:
            return self.cache[k]
        except TypeError:  # k is not hashable
            k = dict_key(k)
            return self.cache[k]

    def __get__(self, inst, cls):
        r"""
        Get a :class:`CachedMethodCaller` bound to a specific
        instance of the class of the cached method.

        NOTE:

        :class:`CachedMethodCaller` has a separate ``__get__``
        since the categories framework creates and caches the
        return value of ``CachedMethod.__get__`` with
        ``inst==None``.

        This getter attempts to assign a bound method as an
        attribute to the given instance. If this is not
        possible (for example, for some extension classes),
        it is attempted to find an attribute ``_cached_methods``,
        and store/retrieve the bound method there. In that
        way, cached methods can be implemented for extension
        classes deriving from :class:`~sage.structure.parent.Parent`
        and :class:`~sage.structure.element.Element`.

        TESTS:

        Due to the separate ``__get__`` method, it is possible
        to define a cached method in one class and use it as
        an attribute of another class. ::

            sage: class Foo:
            ....:     @cached_method
            ....:     def f(self, y):
            ....:         return y - 1
            sage: class Bar:
            ....:     f = Foo.f
            sage: b1 = Bar()
            sage: b2 = Bar()

        The :class:`CachedMethod` is replaced by an instance
        of :class:`CachedMethodCaller` that (by :issue:`8611`)
        is set as an attribute. Hence, we have::

            sage: b1.f is b1.f
            True

        Any instance of ``Bar`` gets its own instance of
        :class:`CachedMethodCaller`::

            sage: b1.f is b2.f
            False

        The method caller knows the instance that it belongs
        to::

            sage: Foo.f._instance is None
            True
            sage: b1.f._instance is b1
            True
            sage: b2.f._instance is b2
            True

        An extension class can inherit a cached method from the
        parent or element class of a category (:issue:`11115`).
        See :class:`CachedMethodCaller` for examples.

        Verify that :issue:`16337` has been resolved::

            sage: class Foo:
            ....:     @cached_method(key=lambda self,y: y+1)
            ....:     def f(self, y):
            ....:         return y - 1
            sage: class Bar:
            ....:     f = Foo.f

            sage: b = Bar()
            sage: b.f(0)
            -1
            sage: b.f.cache
            {1: -1}
        """
        # This is for Parents or Elements that do not allow attribute assignment
        try:
            return (<dict>inst._cached_methods)[self._cachedmethod._cachedfunc.__name__]
        except (AttributeError, TypeError, KeyError):
            pass

        cls = type(self)
        Caller = cls(self._cachedmethod, inst,
                     cache=self._cachedmethod._get_instance_cache(inst),
                     name=self._cachedmethod._cachedfunc.__name__,
                     key=self.key, do_pickle=self.do_pickle)

        try:
            setattr(inst, self._cachedmethod._cachedfunc.__name__, Caller)
            return Caller
        except AttributeError:
            pass
        try:
            if inst._cached_methods is None:
                inst._cached_methods = {self._cachedmethod._cachedfunc.__name__: Caller}
            else:
                (<dict>inst._cached_methods)[self._cachedmethod._cachedfunc.__name__] = Caller
        except AttributeError:
            pass
        return Caller

    def precompute(self, arglist, num_processes=1):
        """
        Cache values for a number of inputs.  Do the computation
        in parallel, and only bother to compute values that we
        haven't already cached.

        INPUT:

        - ``arglist`` -- list (or iterables) of arguments for which
          the method shall be precomputed

        - ``num_processes`` -- number of processes used by
          :func:`~sage.parallel.decorate.parallel`

        EXAMPLES::

            sage: class Foo():
            ....:     @cached_method
            ....:     def f(self, i):
            ....:         return i^2
            sage: foo = Foo()
            sage: foo.f(3)
            9
            sage: foo.f(1)
            1
            sage: foo.f.precompute(range(2), 2)
            sage: foo.f.cache == {((0,), ()): 0, ((1,), ()): 1, ((3,), ()): 9}
            True
        """
        from sage.parallel.decorate import parallel, normalize_input
        P = parallel(num_processes)(self._instance_call)
        cdef list arglist2 = []
        for a in arglist:
            ak = normalize_input(a)
            if self.get_key_args_kwds(ak[0], ak[1]) not in self.cache:
                arglist2.append(ak)
        for ((args, kwargs), val) in P(arglist2):
            self.set_cache(val, *args, **kwargs)


cdef class CachedMethodCallerNoArgs(CachedFunction):
    """
    Utility class that is used by :class:`CachedMethod` to bind a
    cached method to an instance, in the case of a method that does
    not accept any arguments except ``self``.

    .. NOTE::

        The return value ``None`` would not be cached. So, if you have
        a method that does not accept arguments and may return ``None``
        after a lengthy computation, then ``@cached_method`` should not
        be used.

    EXAMPLES::

        sage: P.<a,b,c,d> = QQ[]
        sage: I = P*[a,b]
        sage: I.gens
        Cached version of <function ...gens at 0x...>
        sage: type(I.gens)
        <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
        sage: I.gens is I.gens
        True
        sage: I.gens() is I.gens()
        True

    TESTS:

    As of :issue:`15692` the contents of the cache are not pickled anymore::

        sage: class A:
        ....:    @cached_method
        ....:    def bar(self):
        ....:        return 4
        sage: import __main__
        sage: __main__.A = A
        sage: a = A()
        sage: a.bar()
        4
        sage: a.bar.cache
        4
        sage: b = loads(dumps(a))
        sage: b.bar.cache

    The parameter ``do_pickle`` can be used to change this behaviour::

        sage: class A:
        ....:    @cached_method(do_pickle=True)
        ....:    def bar(self):
        ....:        return 4

        sage: __main__.A = A
        sage: a = A()
        sage: a.bar()
        4
        sage: a.bar.cache
        4
        sage: b = loads(dumps(a))
        sage: b.bar.cache
        4

    AUTHOR:

    - Simon King (2011-04)
    """
    def __init__(self, inst, f, cache=None, name=None, do_pickle=None):
        """
        EXAMPLES::

            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self):
            ....:         return self._x^2
            sage: a = Foo(2)
            sage: print(a.f.cache)
            None
            sage: a.f()
            4
            sage: a.f.cache
            4
        """
        # initialize CachedFunction
        if isinstance(f, str):
            try:
                F = getattr(inst.__class__, f)
            except AttributeError:
                F = getattr(inst, f)
            if isinstance(F, CachedFunction):
                f = F.f
            else:
                f = F
        self._common_init(f, None, name=name, do_pickle=do_pickle)
        # This is for unpickling a CachedMethodCallerNoArgs out
        # of an old CachedMethodCaller:
        cachename = '_cache__' + self.__name__
        if hasattr(inst, cachename):
            # This is for data that are pickled in an old format
            CACHE = getattr(inst, cachename)
            if len(CACHE) > 1:
                raise TypeError("Apparently you are opening a pickle in which '{}' was a method accepting arguments".format(name))
            if len(CACHE) == 1:
                self.cache = CACHE.values()[0]
            else:
                self.cache = cache
            delattr(inst, cachename)
        else:
            self.cache = cache  # None means: the underlying method will be called
        self._instance = inst

    def __reduce__(self):
        """
        Since functions cannot be pickled, the cached method caller
        is pickled by a :class:`CachedMethodPickle`, that replaces
        itself by an actual :class:`CachedMethodCallerNoArgs` as soon
        as it is asked to do anything.

        TESTS::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens()
            [a, b]
            sage: I.gens
            Cached version of <function ...gens at 0x...>
            sage: J = loads(dumps(I))
            sage: J.gens
            Pickle of the cached method "gens"
            sage: J.gens.cache  # the cache is dropped because gens is not marked with do_pickle=True
            sage: J.gens
            Cached version of <function ...gens at 0x...>
        """
        if self.do_pickle:
            return CachedMethodPickle, (self._instance, self.__name__, self.cache)
        return CachedMethodPickle, (self._instance, self.__name__)

    def _instance_call(self):
        """
        Call the cached method without using the cache.

        EXAMPLES::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens()
            [a, b]
            sage: I.gens._instance_call() is I.gens()
            False
            sage: I.gens._instance_call() == I.gens()
            True
        """
        return self.f(self._instance)

    def __call__(self):
        """
        Call the cached method.

        EXAMPLES::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens()    # indirect doctest
            [a, b]
            sage: I.gens() is I.gens()
            True
        """
        if self.cache is None:
            f = self.f
            self.cache = f(self._instance)
        return self.cache

    def set_cache(self, value):
        """
        Override the cache with a specific value.

        .. NOTE::

            ``None`` is not suitable for a cached value. It would be
            interpreted as an empty cache, forcing a new computation.

        EXAMPLES::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens()
            [a, b]
            sage: I.gens.set_cache('bar')
            sage: I.gens()
            'bar'

        The cache can be emptied and thus the original value will
        be reconstructed::

            sage: I.gens.clear_cache()
            sage: I.gens()
            [a, b]

        The attempt to assign ``None`` to the cache fails::

            sage: I.gens.set_cache(None)
            sage: I.gens()
            [a, b]
        """
        self.cache = value

    def clear_cache(self):
        r"""
        Clear the cache dictionary.

        EXAMPLES::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens()
            [a, b]
            sage: I.gens.set_cache('bar')
            sage: I.gens()
            'bar'

        The cache can be emptied and thus the original value will
        be reconstructed::

            sage: I.gens.clear_cache()
            sage: I.gens()
            [a, b]
        """
        self.cache = None

    def is_in_cache(self):
        """
        Answers whether the return value is already in the cache.

        .. NOTE::

            Recall that a cached method without arguments cannot cache
            the return value ``None``.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: I = P*[x,y]
            sage: I.gens.is_in_cache()
            False
            sage: I.gens()
            [x, y]
            sage: I.gens.is_in_cache()
            True
        """
        return self.cache is not None

    def __get__(self, inst, cls):
        """
        Get a :class:`CachedMethodCallerNoArgs` bound to a specific
        instance of the class of the cached method.

        NOTE:

        :class:`CachedMethodCallerNoArgs` has a separate ``__get__``
        since the categories framework creates and caches the
        return value of ``CachedMethod.__get__`` with
        ``inst==None``.

        This getter attempts to assign a bound method as an
        attribute to the given instance. If this is not
        possible (for example, for some extension classes),
        it is attempted to find an attribute ``_cached_methods``,
        and store/retrieve the bound method there. In that
        way, cached methods can be implemented for extension
        classes deriving from :class:`~sage.structure.parent.Parent`
        and :class:`~sage.structure.element.Element`.

        TESTS:

        Due to the separate ``__get__`` method, it is possible
        to define a cached method in one class and use it as
        an attribute of another class. ::

            sage: class Foo:
            ....:     def __init__(self, n):
            ....:         self.__n = n
            ....:     @cached_method
            ....:     def f(self):
            ....:         return self.__n^2
            sage: class Bar:
            ....:     f = Foo.f
            sage: b1 = Bar()
            sage: b2 = Bar()

        The :class:`CachedMethod` is replaced by an instance of
        :class:`CachedMethodCallerNoArgs` that is set as an
        attribute. Hence, we have::

            sage: b1.f is b1.f
            True
            sage: type(b1.f)
            <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>

        Any instance of ``Bar`` gets its own instance of
        :class:`CachedMethodCaller`::

            sage: b1.f is b2.f
            False

        The method caller knows the instance that it belongs
        to::

            sage: Foo.f._instance is None
            True
            sage: b1.f._instance is b1
            True
            sage: b2.f._instance is b2
            True
        """
        # This is for Parents or Elements that do not allow attribute assignment
        try:
            return (<dict>inst._cached_methods)[self.__name__]
        except (AttributeError, TypeError, KeyError):
            pass
        Caller = CachedMethodCallerNoArgs(inst, self.f, name=self.__name__, do_pickle=self.do_pickle)
        try:
            setattr(inst, self.__name__, Caller)
            return Caller
        except AttributeError:
            pass
        try:
            if inst._cached_methods is None:
                inst._cached_methods = {self.__name__: Caller}
            else:
                (<dict>inst._cached_methods)[self.__name__] = Caller
        except AttributeError:
            pass
        return Caller


cdef class GloballyCachedMethodCaller(CachedMethodCaller):
    """
    Implementation of cached methods in case that the cache is not
    stored in the instance, but in some global object. In particular,
    it is used to implement :class:`CachedInParentMethod`.

    The only difference is that the instance is used as part of the
    key.
    """
    cdef get_key_args_kwds(self, tuple args, dict kwds):
        """
        Return the key in the cache to be used when ``args`` and
        ``kwds`` are passed in as parameters.

        EXAMPLES::

            sage: class MyParent(Parent):
            ....:     pass
            sage: class MyElement():           # indirect doctest
            ....:     def __init__(self, x):
            ....:         self.x = x
            ....:     def parent(self):
            ....:         return MyParent()
            ....:     @cached_in_parent_method
            ....:     def f(self):
            ....:         return self.x^2
            sage: a = MyElement(2)
            sage: a.f.get_key()
            (<__main__.MyElement object at 0x...>, ((), ()))
            sage: a.f.get_key()[0] is a
            True
        """
        k = self.fix_args_kwds(args, kwds)
        return (self._instance, k)


cdef class CachedMethod():
    """
    A decorator that creates a cached version of an instance
    method of a class.

    .. NOTE::

        For proper behavior, the method must be a pure function (no side
        effects). Arguments to the method must be hashable or transformed into
        something hashable using ``key`` or they must define
        :meth:`sage.structure.sage_object.SageObject._cache_key`.

    EXAMPLES::

        sage: class Foo():
        ....:     @cached_method
        ....:     def f(self, t, x=2):
        ....:         print('computing')
        ....:         return t**x
        sage: a = Foo()

    The example shows that the actual computation
    takes place only once, and that the result is
    identical for equivalent input::

        sage: res = a.f(3, 2); res
        computing
        9
        sage: a.f(t = 3, x = 2) is res
        True
        sage: a.f(3) is res
        True

    Note, however, that the :class:`CachedMethod` is replaced by a
    :class:`CachedMethodCaller` or :class:`CachedMethodCallerNoArgs`
    as soon as it is bound to an instance or class::

        sage: P.<a,b,c,d> = QQ[]
        sage: I = P*[a,b]
        sage: type(I.__class__.gens)
        <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>

    So, you would hardly ever see an instance of this class alive.

    The parameter ``key`` can be used to pass a function which creates a
    custom cache key for inputs. In the following example, this parameter is
    used to ignore the ``algorithm`` keyword for caching::

        sage: class A():
        ....:     def _f_normalize(self, x, algorithm): return x
        ....:     @cached_method(key=_f_normalize)
        ....:     def f(self, x, algorithm='default'): return x
        sage: a = A()
        sage: a.f(1, algorithm='default') is a.f(1) is a.f(1, algorithm='algorithm')
        True

    The parameter ``do_pickle`` can be used to enable pickling of the cache.
    Usually the cache is not stored when pickling::

        sage: class A():
        ....:     @cached_method
        ....:     def f(self, x): return None
        sage: import __main__
        sage: __main__.A = A
        sage: a = A()
        sage: a.f(1)
        sage: len(a.f.cache)
        1
        sage: b = loads(dumps(a))
        sage: len(b.f.cache)
        0

    When ``do_pickle`` is set, the pickle contains the contents of the cache::

        sage: class A():
        ....:     @cached_method(do_pickle=True)
        ....:     def f(self, x): return None
        sage: __main__.A = A
        sage: a = A()
        sage: a.f(1)
        sage: len(a.f.cache)
        1
        sage: b = loads(dumps(a))
        sage: len(b.f.cache)
        1

    Cached methods cannot be copied like usual methods, see :issue:`12603`.
    Copying them can lead to very surprising results::

        sage: class A:
        ....:     @cached_method
        ....:     def f(self):
        ....:         return 1
        sage: class B:
        ....:     g=A.f
        ....:     def f(self):
        ....:         return 2

        sage: b=B()
        sage: b.f()
        2
        sage: b.g()
        1
        sage: b.f()
        1
    """
    def __init__(self, f, name=None, key=None, do_pickle=None):
        """
        EXAMPLES::

            sage: class Foo():
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, n):
            ....:         return self._x^n
            ....:     @cached_method
            ....:     def f0(self):
            ....:         return self._x^2
            sage: a = Foo(2)
            sage: a.f(2)
            4
            sage: a.f0()
            4

        The computations in method ``f`` are tried to store in a
        dictionary assigned to the instance ``a``::

            sage: hasattr(a, '_cache__f')
            True
            sage: a._cache__f
            {((2,), ()): 4}

        As a shortcut, useful to speed up internal computations,
        the same dictionary is also available as an attribute
        of the ``CachedMethodCaller``::

            sage: type(a.f)
            <class 'sage.misc.cachefunc.CachedMethodCaller'>
            sage: a.f.cache is a._cache__f
            True

        Note that if the instance ``a`` would not accept attribute
        assignment, the computations would still be cached in
        ``a.f.cache``, and they would in fact be preserved when
        pickling.

        The cached method ``f0`` accepts no arguments, which allows
        for an improved way of caching: By an attribute of the cached
        method itself. This cache is *only* available in that way, i.e.,
        it is not additionally stored as an attribute of ``a``::

            sage: type(a.f0)
            <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
            sage: a.f0.cache
            4
            sage: sorted(n for n in dir(a) if not n.startswith('__'))
            ['_cache__f', '_x', 'f', 'f0']

        The cached method has its name and module set::

            sage: f = Foo.__dict__["f"]
            sage: f.__name__
            'f'
            sage: f.__module__
            '__main__'
        """
        self._cache_name = '_cache__' + (name or f.__name__)
        self._cachedfunc = CachedFunction(f, classmethod=True, name=name, key=key, do_pickle=do_pickle)
        self.__name__ = self._cachedfunc.__name__
        self.__cached_module__ = self._cachedfunc.__module__

    @property
    def __module__(self):
        return self.__cached_module__

    def __call__(self, inst, *args, **kwds):
        """
        Call the cached method as a function on an instance.

        INPUT:

        - ``inst`` -- an instance on which the method is to be called
        - further positional or named arguments

        EXAMPLES::


            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class Foo():
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, n=2):
            ....:         return self._x^n
            ....:     g = deprecated_function_alias(57, f)
            sage: a = Foo(2)
            sage: Foo.__dict__['f'](a)
            4

        This uses the cache as usual::

            sage: Foo.__dict__['f'](a) is a.f()
            True

        This feature makes cached methods compatible with
        :meth:`sage.misc.superseded.deprecated_function_alias`::

            sage: a.g() is a.f()
            doctest:...: DeprecationWarning: g is deprecated. Please use f instead.
            See https://github.com/sagemath/sage/issues/57 for details.
            True
            sage: Foo.g(a) is a.f()
            True
        """
        return self.__get__(inst)(*args, **kwds)

    cpdef _get_instance_cache(self, inst):
        """
        Return the cache dictionary.

        TESTS::

            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     @cached_method
            ....:     def f(self, n=2):
            ....:         return self._x^n
            sage: a = Foo(2)
            sage: a.f()
            4

        Note that we cannot provide a direct test, since ``a.f`` is
        an instance of :class:`CachedMethodCaller`.  But during its
        initialisation, this method was called in order to provide the
        cached method caller with its cache, and, if possible, assign
        it to an attribute of ``a``.  So, the following is an indirect
        doctest::

            sage: a.f.cache    # indirect doctest
            {((2,), ()): 4}
            sage: a._cache__f
            {((2,), ()): 4}
        """
        default = {} if self._cachedfunc.do_pickle else NonpicklingDict()
        try:
            return inst.__dict__.setdefault(self._cache_name, default)
        except AttributeError:
            return default

    def __get__(self, object inst, cls):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached method.

        TESTS::

            sage: class Foo:
            ....:     @cached_method
            ....:     def f(self):
            ....:         return 1
            ....:     @cached_method
            ....:     def g(self, x, n=3):
            ....:         return x^n
            sage: a = Foo()
            sage: type(a.f)
            <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
            sage: type(a.g)
            <class 'sage.misc.cachefunc.CachedMethodCaller'>

        By :issue:`8611`, it is attempted to set the
        CachedMethodCaller as an attribute of the instance ``a``,
        replacing the original cached attribute. Therefore, the
        ``__get__`` method will be used only once, which saves much
        time. Hence, we have::

            sage: a.f is a.f
            True
            sage: a.g is a.g
            True

        Verify that :issue:`16337` has been resolved::

            sage: class Foo:
            ....:     @cached_method(key=lambda self, x:x+1)
            ....:     def f(self, x=0):
            ....:         return x

            sage: a = Foo()
            sage: a.f(0)
            0
            sage: a.f.cache
            {1: 0}
        """
        # This is for Parents or Elements that do not allow attribute assignment:
        cdef str name
        try:
            name = self._cachedfunc.__name__
        except AttributeError:
            name = self.__name__
        try:
            return (<dict>inst._cached_methods)[name]
        except (AttributeError, TypeError, KeyError):
            pass
        # Apparently we need to construct the caller.
        # Since we have an optimized version for functions that do not accept arguments,
        # we need to analyse the argspec
        f = self._cachedfunc.f
        if self.nargs == 0:
            if isinstance(f, object) and not isfunction(f):
                try:
                    if METH_NOARGS&PyCFunction_GetFlags(f.__get__(inst, cls)):
                        self.nargs = 1
                except Exception:
                    pass
            if self.nargs == 0:
                args, varargs, keywords, defaults, kwonlyargs, kwonlydefaults, annotations = sage_getargspec(f)
                if varargs is None and keywords is None and len(args)<=1:
                    self.nargs = 1
                else:
                    self.nargs = 2  # don't need the exact number
        if self.nargs == 1:
            Caller = CachedMethodCallerNoArgs(inst, f, name=name, do_pickle=self._cachedfunc.do_pickle)
        else:
            Caller = CachedMethodCaller(self, inst,
                                        cache=self._get_instance_cache(inst),
                                        name=name,
                                        key=self._cachedfunc.key,
                                        do_pickle=self._cachedfunc.do_pickle)
        try:
            setattr(inst, name, Caller)
            return Caller
        except AttributeError:
            pass
        try:
            if inst._cached_methods is None:
                inst._cached_methods = {name: Caller}
            else:
                (<dict>inst._cached_methods)[name] = Caller
        except AttributeError:
            pass
        return Caller


cdef class CachedSpecialMethod(CachedMethod):
    """
    Cached version of *special* python methods.

    IMPLEMENTATION:

    For new style classes ``C``, it is not possible to override a special
    method, such as ``__hash__``, in the ``__dict__`` of an instance ``c`` of
    ``C``, because Python will for efficiency reasons always use what is
    provided by the class, not by the instance.

    By consequence, if ``__hash__`` would be wrapped by using
    :class:`CachedMethod`, then ``hash(c)`` will access ``C.__hash__`` and bind
    it to ``c``, which means that the ``__get__`` method of
    :class:`CachedMethod` will be called. But there, we assume that Python has
    already inspected ``__dict__``, and thus a :class:`CachedMethodCaller`
    will be created over and over again.

    Here, the ``__get__`` method will explicitly access the ``__dict__``, so that
    ``hash(c)`` will rely on a single :class:`CachedMethodCaller` stored in
    the ``__dict__``.

    EXAMPLES::

        sage: class C:
        ....:     @cached_method
        ....:     def __hash__(self):
        ....:         print("compute hash")
        ....:         return int(5)
        ....:
        sage: c = C()
        sage: type(C.__hash__)
        <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>

    The hash is computed only once, subsequent calls will use the value from
    the cache. This was implemented in :issue:`12601`.

    ::

        sage: hash(c)       # indirect doctest
        compute hash
        5
        sage: hash(c)
        5
    """
    def __get__(self, object inst, cls):
        """
        Bind a :class:`CachedMethodCaller` to a specific instance, using ``__dict__``.

        EXAMPLES::

            sage: class C:
            ....:     @cached_method
            ....:     def __hash__(self):
            ....:         print("compute hash")
            ....:         return int(5)
            sage: c = C()
            sage: type(C.__hash__)
            <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
            sage: hash(c)       # indirect doctest
            compute hash
            5
            sage: hash(c)
            5

        Verify that :issue:`16337` has been resolved::

            sage: class Foo:
            ....:     @cached_method(key=lambda self, x:x+1)
            ....:     def __hash__(self, x=0):
            ....:         return x

            sage: a = Foo()
            sage: a.__hash__(0)
            0
            sage: a.__hash__.cache
            {1: 0}
        """
        # This is for Parents or Elements that do not allow attribute assignment:
        cdef str name
        try:
            name = self._cachedfunc.__name__
        except AttributeError:
            name = self.__name__
        cdef dict D = None
        if inst is not None:
            try:
                D = inst.__dict__
            except (TypeError, AttributeError):
                try:
                    D = inst._cached_methods
                except (TypeError, AttributeError):
                    raise TypeError("For a cached special method, either attribute assignment or a public '_cached_methods' attribute of type <dict> is needed")
            if D is None:
                # This can only happen in the case of _cached_methods
                D = inst._cached_methods = {}
            else:
                try:
                    return D[name]
                except KeyError:
                    pass
        # Apparently we need to construct the caller.
        # Since we have an optimized version for functions that do not accept arguments,
        # we need to analyse the argspec
        f = self._cachedfunc.f
        if self.nargs == 0:
            args, varargs, keywords, defaults, kwonlyargs, kwonlydefaults, annotations = sage_getargspec(f)
            if varargs is None and keywords is None and len(args)<=1:
                self.nargs = 1
                Caller = CachedMethodCallerNoArgs(inst, f, name=name, do_pickle=self._cachedfunc.do_pickle)
            else:
                self.nargs = 2  # don't need the exact number
                Caller = CachedMethodCaller(self, inst,
                                            cache=self._get_instance_cache(inst),
                                            name=name,
                                            key=self._cachedfunc.key,
                                            do_pickle=self._cachedfunc.do_pickle)
        elif self.nargs == 1:
            Caller = CachedMethodCallerNoArgs(inst, f, name=name, do_pickle=self._cachedfunc.do_pickle)
        else:
            Caller = CachedMethodCaller(self, inst,
                                        cache=self._get_instance_cache(inst),
                                        name=name,
                                        key=self._cachedfunc.key,
                                        do_pickle=self._cachedfunc.do_pickle)
        if inst is not None:
            try:
                setattr(inst, name, Caller)
                return Caller
            except AttributeError:
                pass
            D[name] = Caller
        return Caller


@decorator_keywords
def cached_method(f, name=None, key=None, do_pickle=None):
    """
    A decorator for cached methods.

    EXAMPLES:

    In the following examples, one can see how a cached method works
    in application. Below, we demonstrate what is done behind the scenes::

        sage: class C:
        ....:     @cached_method
        ....:     def __hash__(self):
        ....:         print("compute hash")
        ....:         return int(5)
        ....:     @cached_method
        ....:     def f(self, x):
        ....:         print("computing cached method")
        ....:         return x*2
        sage: c = C()
        sage: type(C.__hash__)
        <class 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
        sage: hash(c)
        compute hash
        5

    When calling a cached method for the second time with the same arguments,
    the value is gotten from the cache, so that a new computation is not
    needed::

        sage: hash(c)
        5
        sage: c.f(4)
        computing cached method
        8
        sage: c.f(4) is c.f(4)
        True

    Different instances have distinct caches::

        sage: d = C()
        sage: d.f(4) is c.f(4)
        computing cached method
        False
        sage: d.f.clear_cache()
        sage: c.f(4)
        8
        sage: d.f(4)
        computing cached method
        8

    Using cached methods for the hash and other special methods was
    implemented in :issue:`12601`, by means of :class:`CachedSpecialMethod`. We
    show that it is used behind the scenes::

        sage: cached_method(c.__hash__)
        <sage.misc.cachefunc.CachedSpecialMethod object at ...>
        sage: cached_method(c.f)
        <sage.misc.cachefunc.CachedMethod object at ...>

    The parameter ``do_pickle`` can be used if the contents of the cache should be
    stored in a pickle of the cached method. This can be dangerous with special
    methods such as ``__hash__``::

        sage: class C:
        ....:     @cached_method(do_pickle=True)
        ....:     def __hash__(self):
        ....:         return id(self)

        sage: import __main__
        sage: __main__.C = C
        sage: c = C()
        sage: hash(c)  # random output
        sage: d = loads(dumps(c))
        sage: hash(d) == hash(c)
        True

    However, the contents of a method's cache are not pickled unless ``do_pickle``
    is set::

        sage: class C:
        ....:     @cached_method
        ....:     def __hash__(self):
        ....:         return id(self)

        sage: __main__.C = C
        sage: c = C()
        sage: hash(c)  # random output
        sage: d = loads(dumps(c))
        sage: hash(d) == hash(c)
        False
    """
    cdef str fname = name or f.__name__
    if fname in special_method_names:
        return CachedSpecialMethod(f, name, key=key, do_pickle=do_pickle)
    return CachedMethod(f, name, key=key, do_pickle=do_pickle)


cdef class CachedInParentMethod(CachedMethod):
    r"""
    A decorator that creates a cached version of an instance
    method of a class.

    In contrast to :class:`CachedMethod`,
    the cache dictionary is an attribute of the parent of
    the instance to which the method belongs.

    ASSUMPTION:

    This way of caching works only if

    - the instances *have* a parent, and
    - the instances are hashable (they are part of the cache key) or they
      define :meth:`sage.structure.sage_object.SageObject._cache_key`

    NOTE:

    For proper behavior, the method must be a pure function (no side effects).
    If this decorator is used on a method, it will have identical output on
    equal elements. This is since the element is part of the hash key.
    Arguments to the method must be hashable or define
    :meth:`sage.structure.sage_object.SageObject._cache_key`.  The instance it
    is assigned to must be hashable.

    Examples can be found at :mod:`~sage.misc.cachefunc`.
    """

    def __init__(self, f, name=None, key=None, do_pickle=None):
        """
        Construct a new method with cache stored in the parent of the instance.

        See also ``cached_method`` and ``cached_function``.

        EXAMPLES::

            sage: class MyParent(Parent):
            ....:     pass
            sage: class Foo:                         # indirect doctest
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:         self._parent = MyParent()
            ....:     def parent(self):
            ....:         return self._parent
            ....:     @cached_in_parent_method
            ....:     def f(self):
            ....:         return self._x^2
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: b = Foo(3)
            sage: b.f()
            9
            sage: hasattr(a.parent(), '_cache__element_f')
            True

        For speeding up internal computations, this dictionary
        is also accessible as an attribute of the CachedMethodCaller
        (by :issue:`8611`)::

            sage: a.parent()._cache__element_f is a.f.cache
            True

        TESTS:

        Test that ``key`` works::

            sage: class A():
            ....:     def __init__(self):
            ....:         self._parent = MyParent()
            ....:     def parent(self): return self._parent
            ....:     def _f_normalize(self, x, algorithm): return x
            ....:     @cached_in_parent_method(key=_f_normalize)
            ....:     def f(self, x, algorithm='default'): return x
            sage: a = A()
            sage: a.f(1, algorithm='default') is a.f(1) is a.f(1, algorithm='algorithm')
            True

        Test that ``do_pickle`` works. Usually the contents of the cache are not
        pickled::

            sage: class A():
            ....:     def __init__(self):
            ....:         self._parent = MyParent()
            ....:     def parent(self): return self._parent
            ....:     @cached_in_parent_method
            ....:     def f(self, x): return x
            sage: import __main__
            sage: __main__.A = A
            sage: __main__.MyParent = MyParent
            sage: a = A()
            sage: a.f(1)
            1
            sage: len(a.f.cache)
            1
            sage: b = loads(dumps(a))
            sage: len(b.f.cache)
            0

        Pickling can be enabled with ``do_pickle``::

            sage: class A():
            ....:     def __init__(self):
            ....:         self._parent = MyParent()
            ....:     def parent(self): return self._parent
            ....:     @cached_in_parent_method(do_pickle=True)
            ....:     def f(self, x): return x
            sage: __main__.A = A
            sage: a = A()
            sage: a.f(1)
            1
            sage: len(a.f.cache)
            1
            sage: b = loads(dumps(a))
            sage: len(b.f.cache)
            1
        """
        self._cache_name = '_cache__' + 'element_' + (name or f.__name__)
        self._cachedfunc = CachedFunction(f, classmethod=True, name=name, key=key, do_pickle=do_pickle)

    cpdef _get_instance_cache(self, inst):
        """
        Return the cache dictionary, which is stored in the parent.

        EXAMPLES::

            sage: class MyParent(Parent):
            ....:     pass
            sage: class Foo:
            ....:     def __init__(self, x):
            ....:         self._x = x
            ....:     _parent = MyParent()
            ....:     def parent(self):
            ....:         return self._parent
            ....:     def __eq__(self, other):
            ....:         return self._x^2 == other._x^2
            ....:     def __hash__(self):
            ....:         return hash(self._x^2)
            ....:     def __repr__(self):
            ....:         return 'My %s'%self._x
            ....:     @cached_in_parent_method
            ....:     def f(self):
            ....:         return self._x^3
            sage: a = Foo(2)
            sage: a.f()
            8
            sage: a.f.cache   #indirect doctest
            {(My 2, ((), ())): 8}

        Since the key for the cache depends on equality of
        the instances, we obtain *identical* result for
        *equal* instance - even though in this particular
        example the result is wrong::

            sage: b = Foo(-2)
            sage: a is not b
            True
            sage: a == b
            True
            sage: b.f() is a.f()
            True

        Non-equal instances do not share the result of
        the cached method, but they do share the cache::

            sage: c = Foo(3)
            sage: c.f()
            27
            sage: c.f.cache is a.f.cache  #indirect doctest
            True

        Note that the cache is also available as an
        attribute of the cached method, which speeds
        up internal computations::

            sage: a.f.cache is b.f.cache is c.f._cachedmethod._get_instance_cache(c)
            True
        """
        default = {} if self._cachedfunc.do_pickle else NonpicklingDict()
        if inst is None:
            return default
        try:
            P = inst.parent()
            return P.__dict__.setdefault(self._cache_name, default)
        except AttributeError:
            pass
        if not hasattr(P, '_cached_methods'):
            raise TypeError("The parent of this element does not allow attribute assignment\n" +
                            "    and does not descend from the Parent base class.\n" +
                            "    Cannot use CachedInParentMethod.")
        if P._cached_methods is None:
            P._cached_methods = {}
        return (<dict>P._cached_methods).setdefault(self._cache_name, default)

    def __get__(self, inst, cls):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached-in-parent method.
        """
        Caller = GloballyCachedMethodCaller(self, inst, cache=self._get_instance_cache(inst), key=self._cachedfunc.key, do_pickle=self._cachedfunc.do_pickle)
        try:
            setattr(inst, self._cachedfunc.__name__, Caller)
        except AttributeError:
            pass
        return Caller


cached_in_parent_method = decorator_keywords(CachedInParentMethod)


class FileCache():
    """
    :class:`FileCache` is a dictionary-like class which stores keys
    and values on disk.  The keys take the form of a tuple ``(A,K)``

    - ``A`` -- tuple of objects ``t`` where each ``t`` is an
      exact object which is uniquely identified by a short string

    - ``K`` -- tuple of tuples ``(s,v)`` where ``s`` is a valid
      variable name and ``v`` is an exact object which is uniquely
      identified by a short string with letters [a-zA-Z0-9-._]

    The primary use case is the :class:`DiskCachedFunction`.  If
    ``memory_cache == True``, we maintain a cache of objects seen
    during this session in memory -- but we don't load them from
    disk until necessary.  The keys and values are stored in a
    pair of files:

    - ``prefix-argstring.key.sobj`` contains the ``key`` only,
    - ``prefix-argstring.sobj`` contains the tuple ``(key,val)``

    where ``self[key] == val``.

    .. NOTE::

        We assume that each :class:`FileCache` lives in its own directory.
        Use **extreme** caution if you wish to break that assumption.
    """
    def __init__(self, dir, prefix='', memory_cache=False):
        """
        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = True)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((),())]
            1
        """
        if not dir or dir[-1] != '/':
            dir += '/'
        self._dir = dir
        os.makedirs(dir, exist_ok=True)

        self._prefix = prefix + '-'

        if memory_cache:
            self._cache = {}
        else:
            self._cache = None

    def file_list(self):
        """
        Return the list of files corresponding to ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = True, prefix='t')
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: for f in sorted(FC.file_list()): print(f[len(dir):])
            t-.key.sobj
            t-.sobj
            t-1_2.key.sobj
            t-1_2.sobj
            t-a-1.1.key.sobj
            t-a-1.1.sobj
        """
        files = []
        prefix = self._prefix
        dir = self._dir
        l = len(prefix)
        for f in os.listdir(dir):
            if f[:l] == prefix:
                files.append(dir + f)
        return files

    def items(self):
        """
        Return a list of tuples ``(k,v)`` where ``self[k] = v``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: I = FC.items()
            sage: I.sort(); I
            [(((), ()), 1), (((1,), (('a', 1),)), 3), (((1, 2), ()), 2)]
        """
        return [(k, self[k]) for k in self]

    def values(self):
        """
        Return a list of values that are stored in ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: FC[((),(('a',1),))] = 4
            sage: v = FC.values()
            sage: v.sort(); v
            [1, 2, 3, 4]
        """
        return [self[k] for k in self]

    def __iter__(self):
        """
        Return a list of keys of ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: for k in sorted(FC): print(k)
            ((), ())
            ((1,), (('a', 1),))
            ((1, 2), ())
        """
        return iter(self.keys())

    def keys(self):
        """
        Return a list of keys ``k`` where ``self[k]`` is defined.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: K = FC.keys()
            sage: K.sort(); K
            [((), ()), ((1,), (('a', 1),)), ((1, 2), ())]
        """
        cdef list K = []
        from sage.misc.persist import load
        for f in self.file_list():
            if f[-9:] == '.key.sobj':
                K.append(load(f))
        return K

    def clear(self):
        """
        Clear all key, value pairs from ``self`` and unlink the associated files
        from the file cache.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC1 = FileCache(dir, memory_cache=False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache=False, prefix='foo')
            sage: k1 = ((), (('a', 1),))
            sage: t1 = randint(0, 1000)
            sage: k2 = ((), (('b', 1),))
            sage: t2 = randint(0, 1000)
            sage: FC1[k1] = t1
            sage: FC2[k2] = t2
            sage: FC2.clear()
            sage: k1 in FC1
            False
            sage: k2 in FC1
            False
        """
        for k in self:
            del self[k]

    def _filename(self, key):
        """
        Compute the filename associated with a certain key.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False, prefix='foo')
            sage: N = FC._filename(((1,2), (('a',1),('b',2))))
            sage: print(N[len(dir):])
            foo-a-1_b-2.1_2
            sage: N = FC._filename(((), (('a',1),('b',2))))
            sage: print(N[len(dir):])
            foo-a-1_b-2
            sage: N = FC._filename(((1,2), ()))
            sage: print(N[len(dir):])
            foo-1_2
        """
        a, k = key
        kwdstr = '_'.join('%s-%s' % x for x in k)
        argstr = '_'.join('%s' % x for x in a)
        if kwdstr and argstr:
            keystr = kwdstr + '.' + argstr
        else:
            keystr = kwdstr + argstr
        return self._dir + self._prefix + keystr

    def __contains__(self, key):
        """
        Return ``True`` if ``self[key]`` is defined and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: FC[k] = True
            sage: k in FC
            True
            sage: ((),()) in FC
            False
        """
        return os.path.exists(self._filename(key) + '.key.sobj')

    def __getitem__(self, key):
        """
        Return the value set by ``self[key] = val``, in this session
        or a previous one.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC1 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: t = randint(0, 1000)
            sage: FC1[k] = t
            sage: FC2[k] == FC1[k] == t
            True
            sage: FC1[(1,2),(('a',4),('b',2))]
            Traceback (most recent call last):
            ...
            KeyError: ((1, 2), (('a', 4), ('b', 2)))
        """
        from sage.misc.persist import load

        cache = self._cache
        if cache is not None:
            if key in cache:
                return cache[key]

        f = self._filename(key) + '.sobj'
        try:
            k, v = load(f)
        except IOError:
            raise KeyError(key)
        if k != key:
            raise RuntimeError("cache corrupted")

        if cache is not None:
            cache[key] = v
        return v

    def __setitem__(self, key, value):
        """
        Set ``self[key] = value`` and stores both key and value on
        disk.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC1 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: t = randint(0, 1000)
            sage: FC1[k] = t
            sage: FC2[k] == t
            True
            sage: FC1[k] = 2000
            sage: FC2[k]!= t
            True
        """
        from sage.misc.persist import save

        f = self._filename(key)

        save(key, f+'.key.sobj')
        save((key, value), f + '.sobj')
        if self._cache is not None:
            self._cache[key] = value

    def __delitem__(self, key):
        """
        Delete the ``key, value`` pair from ``self`` and unlink the associated
        files from the file cache.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC1 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: t = randint(0, 1000)
            sage: FC1[k] = t
            sage: del FC2[k]
            sage: k in FC1
            False
       """
        f = self._filename(key)
        cache = self._cache
        if cache is not None and key in cache:
            del self._cache[key]
        if os.path.exists(f + '.sobj'):
            os.remove(f + '.sobj')
        if os.path.exists(f + '.key.sobj'):
            os.remove(f + '.key.sobj')


class DiskCachedFunction(CachedFunction):
    """
    Works similar to CachedFunction, but instead, we keep the
    cache on disk (optionally, we keep it in memory too).

    EXAMPLES::

        sage: from sage.misc.cachefunc import DiskCachedFunction
        sage: dir = tmp_dir()
        sage: factor = DiskCachedFunction(factor, dir, memory_cache=True)
        sage: f = factor(2775); f
        3 * 5^2 * 37
        sage: f is factor(2775)
        True
    """
    def __init__(self, f, dir, memory_cache=False, key=None):
        """
        EXAMPLES::

            sage: from sage.misc.cachefunc import DiskCachedFunction
            sage: def foo(x): sleep(x)
            sage: dir = tmp_dir()
            sage: bar = DiskCachedFunction(foo, dir, memory_cache = False)
            sage: w = walltime()
            sage: for i in range(10): bar(1)
            sage: walltime(w) < 2
            True
        """
        CachedFunction.__init__(self, f, key=key)
        prefix = f.__name__
        self.cache = FileCache(dir, prefix=prefix, memory_cache = memory_cache)


class disk_cached_function:
    """
    Decorator for :class:`DiskCachedFunction`.

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: @disk_cached_function(dir)
        ....: def foo(x): return next_prime(2^x)%x
        sage: x = foo(200); x                                                           # needs sage.libs.pari
        11
        sage: @disk_cached_function(dir)
        ....: def foo(x): return 1/x
        sage: foo(200)                                                                  # needs sage.libs.pari
        11
        sage: foo.clear_cache()
        sage: foo(200)
        1/200
    """
    def __init__(self, dir, memory_cache=False, key=None):
        """
        EXAMPLES::

            sage: dir = tmp_dir()
            sage: @disk_cached_function(dir, memory_cache=True)
            ....: def foo(x): return next_prime(2^x)
            sage: x = foo(200)                                                          # needs sage.libs.pari
            sage: x is foo(200)                                                         # needs sage.libs.pari
            True
            sage: @disk_cached_function(dir, memory_cache=False)
            ....: def foo(x): return next_prime(2^x)
            sage: x is foo(200)                                                         # needs sage.libs.pari
            False
        """
        self._dir = dir
        self._memory_cache = memory_cache
        self._key = key

    def __call__(self, f):
        """
        EXAMPLES::

            sage: dir = tmp_dir()
            sage: @disk_cached_function(dir)
            ....: def foo(x): return ModularSymbols(x)
            sage: foo(389)                                                              # needs sage.modular
            Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2
             with sign 0 over Rational Field
        """
        return DiskCachedFunction(f, self._dir, memory_cache=self._memory_cache, key=self._key)


# Add support for _instancedoc_
from sage.misc.instancedoc import instancedoc
instancedoc(CachedFunction)
instancedoc(WeakCachedFunction)
instancedoc(CachedMethodCaller)
instancedoc(CachedMethodCallerNoArgs)
instancedoc(GloballyCachedMethodCaller)
instancedoc(DiskCachedFunction)
