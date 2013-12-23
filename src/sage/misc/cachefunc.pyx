r"""
Cached Functions and Methods

AUTHORS:

- William Stein (inspired by conversation with Justin Walker).
- Mike Hansen (added doctests and made it work with class methods).
- Willem Jan Palenstijn (add CachedMethodCaller for binding cached
  methods to instances).
- Tom Boothby (added DiskCachedFunction).
- Simon King (improved performance, more doctests, cython version,
  CachedMethodCallerNoArgs, weak cached function, cached special methods).

EXAMPLES:

By trac ticket #11115, cached functions and methods are now also
available in Cython code. The following examples cover various ways
of usage.

Python functions::

    sage: @cached_function
    ... def test_pfunc(x):
    ...     '''
    ...     Some documentation
    ...     '''
    ...     return -x
    ...
    sage: test_pfunc(5) is test_pfunc(5)
    True

In some cases, one would only want to keep the result in cache as long
as there is any other reference to the result. By :trac:`12215`, this is
enabled for :class:`~sage.structure.unique_representation.UniqueRepresentation`,
which is used to create unique parents: If an algebraic structure, such
as a finite field, is only temporarily used, then it will not stay in
cache forever. That behaviour is implemented using ``weak_cached_function``,
that behaves the same as ``cached_function``, except that it uses a
:class:`~sage.misc.weak_dict.WeakValueDictionary` for storing the results.
::

    sage: from sage.misc.cachefunc import weak_cached_function
    sage: class A: pass
    sage: @weak_cached_function
    ... def f():
    ...    print "doing a computation"
    ...    return A()
    ...
    sage: a = f()
    doing a computation

The result is cached::

    sage: b = f()
    sage: a is b
    True

However, if there are no strong references left, the result
may be garbage collected, and thus a new computation would
take place::

    sage: del a
    sage: del b
    sage: import gc
    sage: n = gc.collect()
    sage: a = f()
    doing a computation

Cython cdef functions do not allow arbitrary decorators.
However, one can wrap a Cython function and turn it into
a cached function, by :trac:`11115`. We need to provide
the name that the wrapped method or function should have,
since otherwise the name of the original function would
be used::

    sage: cython('''cpdef test_funct(x): return -x''')
    sage: wrapped_funct = cached_function(test_funct, name='wrapped_funct')
    sage: wrapped_funct
    Cached version of <built-in function test_funct>
    sage: wrapped_funct.__name__
    'wrapped_funct'
    sage: wrapped_funct(5)
    -5
    sage: wrapped_funct(5) is wrapped_funct(5)
    True

We can proceed similarly for cached methods of Cython classes,
provided that they allow attribute assignment or have a public
attribute ``__cached_methods`` of type ``<dict>``. Since
:trac:`11115`, this is the case for all classes inheriting from
:class:`~sage.structure.parent.Parent`. See below for a more explicit
example. By :trac:`12951`, cached methods of extension classes can
be defined by simply using the decorater. However, an indirect
approach is still needed for cpdef methods::

    sage: cython_code = ['cpdef test_meth(self,x):',
    ... '    "some doc for a wrapped cython method"',
    ... '    return -x',
    ... 'from sage.all import cached_method',
    ... 'from sage.structure.parent cimport Parent',
    ... 'cdef class MyClass(Parent):',
    ... '    @cached_method',
    ... '    def direct_method(self, x):',
    ... '        "Some doc for direct method"',
    ... '        return 2*x',
    ... '    wrapped_method = cached_method(test_meth,name="wrapped_method")']
    sage: cython(os.linesep.join(cython_code))
    sage: O = MyClass()
    sage: O.direct_method
    Cached version of <method 'direct_method' of '...MyClass' objects>
    sage: O.wrapped_method
    Cached version of <built-in function test_meth>
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
as there is any other reference to the result. By :trac:`12215`, this is
enabled for :class:`~sage.structure.unique_representation.UniqueRepresentation`,
which is used to create unique parents: If an algebraic structure, such
as a finite field, is only temporarily used, then it will not stay in
cache forever. That behaviour is implemented using ``weak_cached_function``,
that behaves the same as ``cached_function``, except that it uses a
:class:`~sage.misc.weak_dict.WeakValueDictionary` for storing the results.
::

    sage: from sage.misc.cachefunc import weak_cached_function
    sage: class A: pass
    sage: @weak_cached_function
    ... def f():
    ...    print "doing a computation"
    ...    return A()
    ...
    sage: a = f()
    doing a computation

The result is cached::

    sage: b = f()
    sage: a is b
    True

However, if there are no strong references left, the result
may be garbage collected, and thus a new computation would
take place::

    sage: del a
    sage: del b
    sage: import gc
    sage: n = gc.collect()
    sage: a = f()
    doing a computation

By :trac:`11115`, even if a parent does not allow attribute
assignment, it can inherit a cached method from the parent class of a
category (previously, the cache would have been broken)::

    sage: cython_code = ["from sage.all import cached_method, cached_in_parent_method, Category, Objects",
    ... "class MyCategory(Category):",
    ... "    @cached_method",
    ... "    def super_categories(self):",
    ... "        return [Objects()]",
    ... "    class ElementMethods:",
    ... "        @cached_method",
    ... "        def element_cache_test(self):",
    ... "            return -self",
    ... "        @cached_in_parent_method",
    ... "        def element_via_parent_test(self):",
    ... "            return -self",
    ... "    class ParentMethods:",
    ... "        @cached_method",
    ... "        def one(self):",
    ... "            return self.element_class(self,1)",
    ... "        @cached_method",
    ... "        def invert(self, x):",
    ... "            return -x"]
    sage: cython('\n'.join(cython_code))
    sage: C = MyCategory()

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

    sage: cython_code = ["from sage.structure.element cimport Element, ElementWithCachedMethod",
    ... "cdef class MyBrokenElement(Element):",
    ... "    cdef public object x",
    ... "    def __init__(self,P,x):",
    ... "        self.x=x",
    ... "        Element.__init__(self,P)",
    ... "    def __neg__(self):",
    ... "        return MyBrokenElement(self.parent(),-self.x)",
    ... "    def _repr_(self):",
    ... "        return '<%s>'%self.x",
    ... "    def __hash__(self):",
    ... "        return hash(self.x)",
    ... "    def __cmp__(left, right):",
    ... "        return (<Element>left)._cmp(right)",
    ... "    def __richcmp__(left, right, op):",
    ... "        return (<Element>left)._richcmp(right,op)",
    ... "    cdef int _cmp_c_impl(left, Element right) except -2:",
    ... "        return cmp(left.x,right.x)",
    ... "    def raw_test(self):",
    ... "        return -self",
    ... "cdef class MyElement(ElementWithCachedMethod):",
    ... "    cdef public object x",
    ... "    def __init__(self,P,x):",
    ... "        self.x=x",
    ... "        ElementWithCachedMethod.__init__(self,P)",
    ... "    def __neg__(self):",
    ... "        return MyElement(self.parent(),-self.x)",
    ... "    def _repr_(self):",
    ... "        return '<%s>'%self.x",
    ... "    def __hash__(self):",
    ... "        return hash(self.x)",
    ... "    def __cmp__(left, right):",
    ... "        return (<Element>left)._cmp(right)",
    ... "    def __richcmp__(left, right, op):",
    ... "        return (<Element>left)._richcmp(right,op)",
    ... "    cdef int _cmp_c_impl(left, Element right) except -2:",
    ... "        return cmp(left.x,right.x)",
    ... "    def raw_test(self):",
    ... "        return -self",
    ... "from sage.structure.parent cimport Parent",
    ... "cdef class MyParent(Parent):",
    ... "    Element = MyElement"]
    sage: cython('\n'.join(cython_code))
    sage: P = MyParent(category=C)
    sage: ebroken = MyBrokenElement(P,5)
    sage: e = MyElement(P,5)

The cached methods inherited by the parent works::

    sage: P.one()
    <1>
    sage: P.one() is P.one()
    True
    sage: P.invert(e)
    <-5>
    sage: P.invert(e) is P.invert(e)
    True

The cached methods inherited by ``MyElement`` works::

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

    sage: e == ebroken
    True
    sage: type(e) == type(ebroken)
    False
    sage: ebroken.element_via_parent_test() is e.element_via_parent_test()
    True

However, the cache of the other inherited method breaks, although the method
as such works::

    sage: ebroken.element_cache_test()
    <-5>
    sage: ebroken.element_cache_test() is ebroken.element_cache_test()
    False

The cache can be emptied::

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

    sage: b = ebroken.element_via_parent_test()
    sage: e.element_via_parent_test.clear_cache()
    sage: b is ebroken.element_via_parent_test()
    False

Introspection works::

    sage: from sage.misc.edit_module import file_and_line
    sage: from sage.misc.sageinspect import sage_getdoc, sage_getfile, sage_getsource
    sage: print sage_getdoc(test_pfunc)
       Some documentation
    sage: print sage_getdoc(O.wrapped_method)
    some doc for a wrapped cython method
    <BLANKLINE>
    sage: print sage_getdoc(O.direct_method)
    Some doc for direct method
    <BLANKLINE>
    sage: print sage_getsource(O.wrapped_method)
    cpdef test_meth(self,x):
        "some doc for a wrapped cython method"
        return -x
    sage: print sage_getsource(O.direct_method)
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
    sage: I = P*[a,b]
    sage: I.gens()
    [a, b]
    sage: I.gens() is I.gens()
    True
    sage: I.groebner_basis()
    [a, b]
    sage: I.groebner_basis() is I.groebner_basis()
    True
    sage: type(I.gens)
    <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
    sage: type(I.groebner_basis)
    <type 'sage.misc.cachefunc.CachedMethodCaller'>

By :trac:`12951`, the cached_method decorator is also supported on non-c(p)def
methods of extension classes, as long as they either support attribute assignment
or have a public attribute of type ``<dict>`` called ``__cached_methods``. The
latter is easy::

    sage: cython_code = [
    ... "from sage.misc.cachefunc import cached_method",
    ... "cdef class MyClass:",
    ... "    cdef public dict __cached_methods",
    ... "    @cached_method",
    ... "    def f(self, a,b):",
    ... "        return a*b"]
    sage: cython(os.linesep.join(cython_code))
    sage: P = MyClass()
    sage: P.f(2,3)
    6
    sage: P.f(2,3) is P.f(2,3)
    True

Providing attribute access is a bit more tricky, since it is needed that
an attribute inherited by the instance from its class can be overridden
on the instance. That is why providing a ``__getattr__`` would not be
enough in the following example::

    sage: cython_code = [
    ... "from sage.misc.cachefunc import cached_method",
    ... "cdef class MyOtherClass:",
    ... "    cdef dict D",
    ... "    def __init__(self):",
    ... "        self.D = {}",
    ... "    def __setattr__(self, n,v):",
    ... "        self.D[n] = v",
    ... "    def __getattribute__(self, n):",
    ... "        try:",
    ... "            return self.D[n]",
    ... "        except KeyError:",
    ... "            pass",
    ... "        return getattr(type(self),n).__get__(self)",
    ... "    @cached_method",
    ... "    def f(self, a,b):",
    ... "        return a+b"]
    sage: cython(os.linesep.join(cython_code))
    sage: Q = MyOtherClass()
    sage: Q.f(2,3)
    5
    sage: Q.f(2,3) is Q.f(2,3)
    True

Note that supporting attribute access is somehow faster than the
easier method::

    sage: timeit("a = P.f(2,3)")   # random
    625 loops, best of 3: 1.3 µs per loop
    sage: timeit("a = Q.f(2,3)")   # random
    625 loops, best of 3: 931 ns per loop

"""
########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                          Mike Hansen <mhansen@gmail.com>
#                     2011 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from function_mangling import ArgumentFixer
import os
from sage.misc.sageinspect import sage_getfile, sage_getsourcelines, sage_getargspec

import sage.misc.weak_dict
from sage.misc.weak_dict import WeakValueDictionary

cdef frozenset special_method_names = frozenset(['__abs__', '__add__',
            '__and__', '__call__', '__cmp__', '__coerce__', '__complex__', '__contains__', '__del__',
            '__delattr__', '__delete__', '__delitem__', '__delslice__', '__dir__', '__div__',
            '__eq__', '__float__', '__floordiv__', '__format__', '__ge__', '__get__', '__getattr__',
            '__getattribute__', '__getitem__', '__getslice__', '__gt__', '__hash__', '__hex__',
            '__iadd__', '__iand__', '__idiv__', '__ifloordiv__', '__ilshift__', '__imod__', '__imul__',
            '__index__', '__init__', '__instancecheck__', '__int__', '__invert__', '__ior__', '__ipow__',
            '__irshift__', '__isub__', '__iter__', '__itruediv__', '__ixor__', '__le__', '__len__',
            '__length_hint__', '__long__', '__lshift__', '__lt__', '__missing__', '__mod__', '__mul__',
            '__ne__', '__neg__', '__new__', '__nonzero__', '__oct__', '__or__', '__pos__', '__pow__',
            '__radd__', '__rand__', '__rdiv__', '__repr__', '__reversed__', '__rfloordiv__', '__rlshift__',
            '__rmod__', '__rmul__', '__ror__', '__rpow__', '__rrshift__', '__rshift__', '__rsub__',
            '__rtruediv__', '__rxor__', '__set__', '__setattr__', '__setitem__', '__setslice__', '__sizeof__',
            '__str__', '__sub__', '__subclasscheck__', '__truediv__', '__unicode__', '__xor__', 'next'])

def _cached_function_unpickle(module,name):
    """
    Unpickling of cached functions.

    NOTE:

    Pickling and unpickling of cached functions is by importing
    from the module in which the function is defined.

    INPUT:

    - ``module``: A string, describing the module to import the
      function from.
    - ``name``: A string, name of the to-be-imported cached function.

    EXAMPLE::

        sage: type(cunningham_prime_factors)
        <type 'sage.misc.cachefunc.CachedFunction'>
        sage: loads(dumps(cunningham_prime_factors)) is cunningham_prime_factors #indirect doctest
        True

    """
    return getattr(__import__(module, fromlist=['']),name)

cdef class CachedFunction(object):
    """
    Create a cached version of a function, which only recomputes
    values it hasn't already computed. Synonyme: ``cached_function``

    INPUT:

    - ``f`` -- a function
    - ``name`` (optional string) -- name that the cached version
      of ``f`` should be provided with.

    If ``f`` is a function, do either ``g = CachedFunction(f)``
    or ``g = cached_function(f)`` to make a cached version of ``f``,
    or put ``@cached_function`` right before the definition of ``f``
    (i.e., use Python decorators)::

        @cached_function
        def f(...):
            ....

    The inputs to the function must be hashable.

    EXAMPLES::

        sage: @cached_function
        ... def mul(x, y=2):
        ...     return x*y
        ...
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

    """
    def __init__(self, f, classmethod=False, name=None):
        """
        Create a cached version of a function, which only recomputes
        values it hasn't already computed. A custom name can be
        provided by an optional argument "name".

        If f is a function, do either g = CachedFunction(f) to make
        a cached version of f, or put @CachedFunction right before
        the definition of f (i.e., use Python decorators, but then
        the optional argument ``name`` can not be used)::

            @CachedFunction
            def f(...):
                ....

        The inputs to the function must be hashable.

        TESTS::

            sage: g = CachedFunction(number_of_partitions)
            sage: g.__name__
            'number_of_partitions'
            sage: 'partitions' in sage.misc.sageinspect.sage_getdoc(g)
            True
            sage: g(5)
            7
            sage: g.cache
            {((5, None, 'default'), ()): 7}
            sage: def f(t=1): print(t)
            sage: h = CachedFunction(f)
            sage: w = walltime()
            sage: h(); h(1); h(t=1)
            1
            sage: walltime(w) < 2
            True

        """
        self.is_classmethod = classmethod
        self._common_init(f, None, name=name)
        self.cache = {}

    def _common_init(self, f, argument_fixer, name=None):
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
        if name is not None:
            self.__name__ = name
        elif hasattr(f, "func_name"):
            self.__name__ = f.func_name
        else:
            self.__name__ = f.__name__
        try:
            self.__module__ = f.__module__
        except AttributeError:
            self.__module__ = f.__objclass__.__module__
        if argument_fixer is not None: # it is None unless the argument fixer
                                       # was known previously. See #15038.
            self._argument_fixer = argument_fixer
            self._fix_to_pos = argument_fixer.fix_to_pos

    cdef argfix_init(self):
        """
        Perform initialization common to CachedFunction and CachedMethodCaller.

        TESTS::

            sage: @cached_function
            ....: def test_cache(x):
            ....:     return -x
            sage: test_cache(1)
            -1
            sage: test_cache._fix_to_pos is not None  # indirect doctest
            True

        """
        A = ArgumentFixer(self.f,classmethod=self.is_classmethod)
        self._argument_fixer = A
        self._fix_to_pos = A.fix_to_pos

    def __reduce__(self):
        """
        Pickling of cached functions.

        TEST::

            sage: type(cunningham_prime_factors)
            <type 'sage.misc.cachefunc.CachedFunction'>
            sage: loads(dumps(cunningham_prime_factors)) is cunningham_prime_factors #indirect doctest
            True

        """
        return _cached_function_unpickle, (self.__module__, self.__name__)

    #########
    ## Introspection
    ##
    ## We provide some methods explicitly, and
    ## forward other questions to the cached function.

    def _sage_doc_(self):
        """
        Provide documentation for the cached function.

        A cached function shall inherit the documentation
        from the function that is wrapped, not from the
        documentation of the wrapper.

        TEST::

            sage: P.<x,y> = QQ[]
            sage: I = P*[x,y]
            sage: from sage.misc.sageinspect import sage_getdoc
            sage: print sage_getdoc(I.groebner_basis) # indirect doctest
               Return the reduced Groebner basis of this ideal.
            ...
               ALGORITHM: Uses Singular, Magma (if available), Macaulay2 (if
               available), or a toy implementation.

        """
        f = self.f
        if hasattr(f, "func_doc"):
            try:
                sourcelines = sage_getsourcelines(f)
            except IOError:
                sourcelines = None
            if sourcelines is not None:
                from sage.env import SAGE_SRC, SAGE_LIB
                filename = sage_getfile(f)
                # The following is a heuristics to get
                # the file name of the cached function
                # or method
                if filename.startswith(SAGE_SRC):
                    filename = filename[len(SAGE_SRC):]
                elif filename.startswith(SAGE_LIB):
                    filename = filename[len(SAGE_LIB):]
                file_info = "File: %s (starting at line %d)\n"%(filename,sourcelines[1])
                doc = file_info+(f.func_doc or '')
            else:
                doc = f.func_doc
        else:
            doc = f.__doc__
        return doc

    def _sage_src_(self):
        """
        Returns the source code for the wrapped function.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getsource
            sage: g = CachedFunction(number_of_partitions)
            sage: 'bober' in sage_getsource(g)  # indirect doctest
            True

        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.f)

    def _sage_src_lines_(self):
        r"""
        Returns the list of source lines and the first line number
        of the wrapped function.

        TEST::

            sage: P.<x,y> = QQ[]
            sage: I = P*[x,y]
            sage: from sage.misc.sageinspect import sage_getsourcelines
            sage: l = "        elif algorithm == 'macaulay2:gb':\n"
            sage: l in sage_getsourcelines(I.groebner_basis)[0] # indirect doctest
            True

        """
        from sage.misc.sageinspect import sage_getsourcelines
        return sage_getsourcelines(self.f)

    def _sage_argspec_(self):
        """
        Return the argspec of the wrapped function or method.

        This was implemented in trac ticket #11115.

        EXAMPLE::

            sage: P.<x,y> = QQ[]
            sage: I = P*[x,y]
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(I.groebner_basis)   # indirect doctest
            ArgSpec(args=['self', 'algorithm', 'deg_bound', 'mult_bound', 'prot'],
            varargs='args', keywords='kwds', defaults=('', None, None,
            False))

        """
        return sage_getargspec(self.f)

    def __call__(self, *args, **kwds):
        """
        Return value from cache or call the wrapped function,
        caching the output.

        TESTS::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: a = g(10^5)   # indirect doctest
            sage: a == number_of_partitions(10^5)
            True
            sage: a is g(10^5)
            True
            sage: a is number_of_partitions(10^5)
            True

        """
        # We shortcut a common case of no arguments
        if args or kwds:
            if self._argument_fixer is None:
                self.argfix_init()
            k = self._fix_to_pos(*args, **kwds)
        else:
            if self._default_key is not None:
                k = self._default_key
            else:
                if self._argument_fixer is None:
                    self.argfix_init()
                k = self._default_key = self._fix_to_pos()

        try:
            return (<dict>self.cache)[k]
        except KeyError:
            w = self.f(*args, **kwds)
            self.cache[k] = w
            return w

    cpdef get_cache(self):
        """
        Returns the cache dictionary.

        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}

        """
        return self.cache

    def is_in_cache(self, *args, **kwds):
        """
        Checks if the argument list is in the cache.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, z, y=0):
            ...           return self._x*z+y
            ...
            sage: a = Foo(2)
            sage: a.f.is_in_cache(3)
            False
            sage: a.f(3)
            6
            sage: a.f.is_in_cache(3,y=0)
            True
        """
        if self._argument_fixer is None:
            self.argfix_init()
        return self._fix_to_pos(*args, **kwds) in (<dict>self.cache)

    def set_cache(self, value, *args, **kwds):
        """
        Set the value for those args and keyword args
        Mind the unintuitive syntax (value first).
        Any idea on how to improve that welcome!

        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: g.set_cache(17, 5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 17}
            sage: g(5)
            17

        DEVELOPER NOTE:

        Is there a way to use the following intuitive syntax?

        ::

            sage: g(5) = 19    # todo: not implemented
            sage: g(5)         # todo: not implemented
            19
        """
        if self._argument_fixer is None:
            self.argfix_init()
        (<dict>self.cache)[self._fix_to_pos(*args, **kwds)] = value

    def get_key(self, *args, **kwds):
        """
        Returns the key in the cache to be used when args
        and kwds are passed in as parameters.

        EXAMPLES::

            sage: @cached_function
            ... def foo(x):
            ...    return x^2
            ...
            sage: foo(2)
            4
            sage: foo.get_key(2)
            ((2,), ())
            sage: foo.get_key(x=3)
            ((3,), ())
        """
        if self._argument_fixer is None:
            self.argfix_init()
        return self._fix_to_pos(*args, **kwds)

    def __repr__(self):
        """
        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: g     # indirect doctest
            Cached version of <function number_of_partitions at 0x...>
        """
        try:
            return "Cached version of %s"%self.f
        except AttributeError:
            return "Cached version of a method (pending reassignment)"

    cpdef clear_cache(self):
        """
        Clear the cache dictionary.

        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: g.clear_cache()
            sage: g.get_cache()
            {}
        """
        cdef object cache = self.cache
        for key in cache.keys():
            del cache[key]

    def precompute(self, arglist, num_processes=1):
        """
        Cache values for a number of inputs.  Do the computation
        in parallel, and only bother to compute values that we
        haven't already cached.

        EXAMPLES::

            sage: @cached_function
            ... def oddprime_factors(n):
            ...     l = [p for p,e in factor(n) if p != 2]
            ...     return len(l)
            sage: oddprime_factors.precompute(range(1,100), 4)
            sage: oddprime_factors.cache[(25,),()]
            1
        """
        from sage.parallel.decorate import parallel, normalize_input
        P = parallel(num_processes)(self.f)
        has_key = self.cache.has_key
        if self._argument_fixer is None:
            self.argfix_init()
        get_key = self._fix_to_pos
        new = lambda x: not has_key(get_key(*x[0],**x[1]))
        arglist = filter(new, map(normalize_input, arglist))
        for ((args,kwargs), val) in P(arglist):
            self.set_cache(val, *args, **kwargs)


cached_function = CachedFunction

cdef class WeakCachedFunction(CachedFunction):
    """
    A version of :class:`CachedFunction` using weak references on the values.

    If f is a function, do either ``g = weak_cached_function(f)`` to make
    a cached version of f, or put ``@weak_cached_function`` right before
    the definition of f (i.e., use Python decorators, but then
    the optional argument ``name`` can not be used)::

        @weak_cached_function
        def f(...):
            ...

    EXAMPLES::

        sage: from sage.misc.cachefunc import weak_cached_function
        sage: class A: pass
        sage: @weak_cached_function
        ... def f():
        ...    print "doing a computation"
        ...    return A()
        ...
        sage: a = f()
        doing a computation

    The result is cached::

        sage: b = f()
        sage: a is b
        True

    However, if there are no strong references left, the result
    may be garbage collected, and thus a new computation would
    take place::

        sage: del a
        sage: del b
        sage: import gc
        sage: n = gc.collect()
        sage: a = f()
        doing a computation

    """
    def __init__(self, f, classmethod=False, name=None):
        """
        The inputs to the function must be hashable.
        The outputs to the function must be weakly referenceable.

        TESTS::

            sage: from sage.misc.cachefunc import weak_cached_function
            sage: class A: pass
            sage: @weak_cached_function
            ... def f():
            ...    return A()
            ...
            sage: f
            Cached version of <function f at ...>

        We demonstrate that pickling works, provided the uncached function
        is available::

            sage: import __main__
            sage: __main__.f = f
            sage: loads(dumps(f))
            Cached version of <function f at ...>
            sage: str(f.cache)
            '<WeakValueDictionary at 0x...>'

        """
        self._common_init(f, None, name=name)
        self.cache = WeakValueDictionary()
    def __call__(self, *args, **kwds):
        """
        Return value from cache or call the wrapped function,
        caching the output.

        TESTS::

            sage: from sage.misc.cachefunc import weak_cached_function
            sage: class A: pass
            sage: @weak_cached_function
            ... def f():
            ...    print "doing a computation"
            ...    return A()
            ...
            sage: a = f()    # indirect doctest
            doing a computation

        The result is cached::

            sage: b = f()
            sage: a is b
            True

        However, if there are no strong references left, the result
        may be garbage collected, and thus a new computation would
        take place::

            sage: del a
            sage: del b
            sage: import gc
            sage: n = gc.collect()
            sage: a = f()
            doing a computation

        """
        # We shortcut a common case of no arguments
        if args or kwds:
            if self._argument_fixer is None:
                self.argfix_init()
            k = self._fix_to_pos(*args, **kwds)
        else:
            if self._default_key is not None:
                k = self._default_key
            else:
                if self._argument_fixer is None:
                    self.argfix_init()
                k = self._default_key = self._fix_to_pos()

        try:
            return self.cache[k]
        except KeyError:
            w = self.f(*args, **kwds)
            self.cache[k] = w
            return w
    def is_in_cache(self, *args, **kwds):
        """
        Checks if the argument list is in the cache.

        EXAMPLES::

            sage: from sage.misc.cachefunc import weak_cached_function
            sage: class A:
            ...     def __init__(self, x):
            ...         self.x = x
            ...
            sage: @weak_cached_function
            ... def f(n):
            ...    return A(n)
            ...
            sage: a = f(5)

        The key 5 is in the cache, as long as there is a strong
        reference to the corresponding value::

            sage: f.is_in_cache(5)
            True

        However, if there are no strong references left, the cached
        item is removed from cache after garbage collection::

            sage: del a
            sage: import gc
            sage: n = gc.collect()
            sage: f.is_in_cache(5)
            False

        """
        if self._argument_fixer is None:
            self.argfix_init()
        return self._fix_to_pos(*args, **kwds) in self.cache

    def set_cache(self, value, *args, **kwds):
        """
        Set the value for those args and keyword args
        Mind the unintuitive syntax (value first).
        Any idea on how to improve that welcome!

        It is required that the given value is weak
        referenceable. The item will be removed from
        cache if the value is garbage collected.

        EXAMPLES::

            sage: from sage.misc.cachefunc import weak_cached_function
            sage: @weak_cached_function
            ... def f(n):
            ...     raise RuntimeError
            ...
            sage: f.set_cache(ZZ, 5)
            sage: f(5)
            Integer Ring

        """
        if self._argument_fixer is None:
            self.argfix_init()
        self.cache[self._fix_to_pos(*args, **kwds)] = value


weak_cached_function = WeakCachedFunction

class CachedMethodPickle(object):
    """
    This class helps to unpickle cached methods.

    NOTE:

    Since trac ticket #8611, a cached method is an attribute
    of the instance (provided that it has a ``__dict__``).
    Hence, when pickling the instance, it would be attempted
    to pickle that attribute as well, but this is a problem,
    since functions can not be pickled, currently. Therefore,
    we replace the actual cached method by a place holder,
    that kills itself as soon as any attribute is requested.
    Then, the original cached attribute is reinstated. But the
    cached values are in fact saved.

    EXAMPLE::

        sage: R.<x, y, z> = PolynomialRing(QQ, 3)
        sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
        sage: I.groebner_basis()
        [y^5*z^3 - 1/4*x^2*z^6 + 1/2*x*y*z^6 + 1/4*y^2*z^6,
         x^2*y*z^3 - x*y^2*z^3 + 2*y^3*z^3 + z^6,
         x*y^3 + y^4 + x*z^3, x^3 + y^3 + z^3]
        sage: I.groebner_basis
        Cached version of <function groebner_basis at 0x...>

    We now pickle and unpickle the ideal. The cached method
    ``groebner_basis`` is replaced by a placeholder::

        sage: J = loads(dumps(I))
        sage: J.groebner_basis
        Pickle of the cached method "groebner_basis"

    But as soon as any other attribute is requested from the
    placeholder, it replaces itself by the cached method, and
    the entries of the cache are actually preserved::

        sage: J.groebner_basis.is_in_cache()
        True
        sage: J.groebner_basis
        Cached version of <function groebner_basis at 0x...>
        sage: J.groebner_basis() == I.groebner_basis()
        True

    TESTS:

    Since Trac Ticket #11115, there is a special implementation for
    cached methods that don't take arguments::

        sage: P.<a,b,c,d> = QQ[]
        sage: I = P*[a,b]
        sage: type(I.gens)
        <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
        sage: type(I.groebner_basis)
        <type 'sage.misc.cachefunc.CachedMethodCaller'>

    We demonstrate that both implementations can be pickled and
    preserve the cache. For that purpose, we assign nonsense to the
    cache. Of course, it is a very bad idea to override the cache in
    that way.  So, please don't try this at home::

        sage: I.groebner_basis.set_cache('foo',algorithm='singular')
        sage: I.groebner_basis(algorithm='singular')
        'foo'
        sage: I.gens.set_cache('bar')
        sage: I.gens()
        'bar'
        sage: J = loads(dumps(I))
        sage: J.gens()
        'bar'
        sage: J.groebner_basis(algorithm='singular')
        'foo'

    Anyway, the cache will be automatically reconstructed after
    clearing it::

        sage: J.gens.clear_cache()
        sage: J.gens()
        [a, b]
        sage: J.groebner_basis.clear_cache()
        sage: J.groebner_basis(algorithm='singular')
        [a, b]

    AUTHOR:

    - Simon King (2011-01)
    """
    def __init__(self, inst, name, cache=None):
        """
        INPUT:

        - ``inst`` - some instance.
        - ``name`` (string) - usually the name of an attribute
          of ``inst`` to which ``self`` is assigned.

        TEST::

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
        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: G = I.groebner_basis()
            sage: J = loads(dumps(I))
            sage: J.groebner_basis  #indirect doctest
            Pickle of the cached method "groebner_basis"
        """
        return 'Pickle of the cached method "%s"'%self._name

    def __reduce__(self):
        """
        This class is a pickle. However, sometimes, pickles
        need to be pickled another time.

        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: I.groebner_basis()
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
            sage: K.groebner_basis.cache
            {(('', None, None, False), ()):
            [y^5*z^3 - 1/4*x^2*z^6 + 1/2*x*y*z^6 + 1/4*y^2*z^6,
             x^2*y*z^3 - x*y^2*z^3 + 2*y^3*z^3 + z^6,
             x*y^3 + y^4 + x*z^3, x^3 + y^3 + z^3]}
        """
        return CachedMethodPickle,(self._instance,self._name,self._cache)

    def __call__(self,*args,**kwds):
        """
        The purpose of this call method is to kill ``self`` and to
        replace it by an actual :class:`CachedMethodCaller`. The last
        thing that ``self`` does before disappearing is to call the
        :class:`CachedMethodCaller` and return the result.

        EXAMPLE::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens
            Cached version of <function gens at 0x...>
            sage: J = loads(dumps(I))
            sage: J.gens
            Pickle of the cached method "gens"
            sage: J.gens()   # indirect doctest
            [a, b]
            sage: J.gens
            Cached version of <function gens at 0x...>

        """
        self._instance.__dict__.__delitem__(self._name)
        CM = getattr(self._instance,self._name)
        if self._cache is not None:
            if isinstance(CM, CachedMethodCallerNoArgs):
                CM.cache = self._cache
            else:
                for k,v in self._cache:
                    CM.cache[k] = v
        return CM(*args,**kwds)

    def __getattr__(self,s):
        """
        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: G = I.groebner_basis()
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

            sage: J.groebner_basis.is_in_cache()  #indirect doctest
            True
            sage: J.groebner_basis
            Cached version of <function groebner_basis at 0x...>

        """
        self._instance.__dict__.__delitem__(self._name)
        CM = getattr(self._instance,self._name)
        if self._cache is not None:
            if isinstance(CM, CachedMethodCallerNoArgs):
                CM.cache = self._cache
            else:
                for k,v in self._cache:
                    CM.cache[k] = v
        return getattr(CM,s)

cdef class CachedMethodCaller(CachedFunction):
    """
    Utility class that is used by :class:`CachedMethod` to bind a
    cached method to an instance.

    NOTE:

    Since Trac Ticket #11115, there is a special implementation
    :class:`CachedMethodCallerNoArgs` for methods that do not take
    arguments.

    EXAMPLE::

        sage: class A:
        ...    @cached_method
        ...    def bar(self,x):
        ...        return x^2
        sage: a = A()
        sage: a.bar
        Cached version of <function bar at 0x...>
        sage: type(a.bar)
        <type 'sage.misc.cachefunc.CachedMethodCaller'>
        sage: a.bar(2) is a.bar(x=2)
        True

    """
    def __init__(self, CachedMethod cachedmethod, inst, cache=None, inst_in_key=False, name=None):
        """
        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self,*args):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f.get_cache()
            {}
            sage: a.f()
            4
            sage: a.f.get_cache()
            {((), ()): 4}
        """
        # initialize CachedFunction. Since the cached method is actually bound
        # to an instance, it now makes sense to initialise the ArgumentFixer
        # and re-use it for all bound cached method callers of the unbound
        # cached method.
        if cachedmethod._cachedfunc._argument_fixer is None:
            cachedmethod._cachedfunc.argfix_init()
        self._common_init(cachedmethod._cachedfunc.f, cachedmethod._cachedfunc._argument_fixer, name=name)
        self.cache = {} if cache is None else cache
        self._instance = inst
        self._inst_in_key = inst_in_key
        self._cachedmethod = cachedmethod

    def __reduce__(self):
        """
        The pickle of a :class:`CachedMethodCaller` unpickles
        to a :class:`CachedMethodPickle`, that is able to replace
        itself by a copy of the original :class:`CachedMethodCaller`.

        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: G = I.groebner_basis()
            sage: J = loads(dumps(I))  #indirect doctest
            sage: J.groebner_basis
            Pickle of the cached method "groebner_basis"
            sage: J.groebner_basis.is_in_cache()
            True
            sage: J.groebner_basis
            Cached version of <function groebner_basis at 0x...>

        """
#        if hasattr(self._instance,self._cachedmethod._cache_name):
#            return CachedMethodPickle,(self._instance,self.__name__)
        if isinstance(self._cachedmethod, CachedInParentMethod) or hasattr(self._instance,self._cachedmethod._cache_name):
            return CachedMethodPickle,(self._instance,self.__name__)
        return CachedMethodPickle,(self._instance,self.__name__,self.cache.items())

    def _instance_call(self, *args, **kwds):
        """
        Call the cached method without using the cache.

        EXAMPLE::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.groebner_basis()
            [a, b]
            sage: I.groebner_basis._instance_call() is I.groebner_basis()
            False
            sage: I.groebner_basis._instance_call() == I.groebner_basis()
            True

        """
        return self._cachedmethod._instance_call(self._instance, *args, **kwds)

    def __call__(self, *args, **kwds):
        """
        Call the cached method.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, x,y=1):
            ...           return x+y
            ...
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

        We test that #5843 is fixed::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, y):
            ...           return self._x
            ...
            sage: a = Foo(2)
            sage: b = Foo(3)
            sage: a.f(b.f)
            2
        """
        # We shortcut a common case of no arguments
        # and we avoid calling another python function,
        # although that means to duplicate code.
        cdef int lenargs
        cdef int nargs
        cdef tuple k
        cdef dict cache = self.cache
        if kwds:
            if self._argument_fixer is None:
                self.argfix_init()
            if self._inst_in_key:
                k = (self._instance,self._fix_to_pos(*args, **kwds))
            else:
                k = self._fix_to_pos(*args, **kwds)
        else:
            if args:
                lenargs = len(args)
                nargs = self._argument_fixer._nargs
                if self._inst_in_key:
                    if lenargs>=nargs:
                        k = (self._instance,(args,()))
                    else:
                        k = (self._instance,(<tuple>args+(<tuple>self._argument_fixer._default_tuple)[-nargs+lenargs:],()))
                else:
                    if lenargs>=nargs:
                        k = (args,())
                    else:
                        k = (<tuple>args+(<tuple>self._argument_fixer._default_tuple)[-nargs+lenargs:],())
            elif self._default_key is not None:
                k = self._default_key
            else:
                if self._argument_fixer is None:
                    self.argfix_init()
                if self._inst_in_key:
                    k = self._default_key = (self._instance,self._fix_to_pos())
                else:
                    k = self._default_key = self._fix_to_pos()
        try:
            return cache[k]
        except KeyError:
            w = self._cachedmethod._instance_call(self._instance, *args, **kwds)
            cache[k] = w
            return w

    def get_key(self, *args, **kwds):
        """
        Convert arguments to the key for this instance's cache.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, y, z=0):
            ...           return self._x * y + z
            ...
            sage: a = Foo(2)
            sage: z = a.f(37)
            sage: k = a.f.get_key(37); k
            ((37, 0), ())
            sage: a.f.get_cache()[k] is z
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
        if self._argument_fixer is None:
            self.argfix_init()
        if self._inst_in_key:
            return (self._instance,self._fix_to_pos(*args,**kwds))
        return self._fix_to_pos(*args,**kwds)

    def __get__(self, inst, cls): #cls=None):
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
        it is attempted to find an attribute ``__cached_methods``,
        and store/retrieve the bound method there. In that
        way, cached methods can be implemented for extension
        classes deriving from :class:`~sage.structure.parent.Parent`
        and :class:`~sage.structure.element.Element`.

        TESTS:

        Due to the separate ``__get__`` method, it is possible
        to define a cached method in one class and use it as
        an attribute of another class.

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, y):
            ...           return y - 1
            sage: class Bar:
            ...       f = Foo.f
            sage: b1 = Bar()
            sage: b2 = Bar()

        The :class:`CachedMethod` is replaced by an instance
        of :class:`CachedMethodCaller` that (by trac ticket
        #8611) is set as an attribute. Hence, we have::

            sage: b1.f is b1.f
            True

        Any instance of ``Bar`` gets its own instance of
        :class:`CachedMethodCaller``::

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
        parent or element class of a category (trac ticket #11115).
        See :class:`CachedMethodCaller` for examples.

        """
        # This is for Parents or Elements that do not allow attribute assignment
        try:
            return (<dict>inst.__cached_methods)[self._cachedmethod._cachedfunc.__name__]
        except (AttributeError,TypeError,KeyError):
            pass
        Caller = CachedMethodCaller(self._cachedmethod, inst, cache=self._cachedmethod._get_instance_cache(inst), inst_in_key=self._inst_in_key, name=self._cachedmethod._cachedfunc.__name__)
        try:
            setattr(inst,self._cachedmethod._cachedfunc.__name__, Caller)
            return Caller
        except AttributeError,msg:
            pass
        try:
            if inst.__cached_methods is None:
                inst.__cached_methods = {self._cachedmethod._cachedfunc.__name__ : Caller}
            else:
                (<dict>inst.__cached_methods)[self._cachedmethod._cachedfunc.__name__] = Caller
        except AttributeError,msg:
            pass
        return Caller

cdef class CachedMethodCallerNoArgs(CachedFunction):
    """
    Utility class that is used by :class:`CachedMethod` to bind a
    cached method to an instance, in the case of a method that does
    not accept any arguments except ``self``.

    NOTE:

    The return value ``None`` would not be cached. So, if you have
    a method that does not accept arguments and may return ``None``
    after a lengthy computation, then ``@cached_method`` should not
    be used.

    EXAMPLE::

        sage: P.<a,b,c,d> = QQ[]
        sage: I = P*[a,b]
        sage: I.gens
        Cached version of <function gens at 0x...>
        sage: type(I.gens)
        <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
        sage: I.gens is I.gens
        True
        sage: I.gens() is I.gens()
        True

    AUTHOR:

    - Simon King (2011-04)
    """
    def __init__(self, inst, f, cache=None, name=None):
        """
        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: print a.f.get_cache()
            None
            sage: a.f()
            4
            sage: a.f.get_cache()
            4

        """
        # initialize CachedFunction
        if isinstance(f,basestring):
            try:
                F = getattr(inst.__class__,f)
            except AttributeError:
                 F = getattr(inst,f)
            if isinstance(F,CachedFunction):
                f = F.f
            else:
                f = F
        self._common_init(f, None, name=name)
        # This is for unpickling a CachedMethodCallerNoArgs out
        # of an old CachedMethodCaller:
        cachename = '_cache__' + self.__name__
        if hasattr(inst, cachename):
            # This is for data that are pickled in an old format
            CACHE = getattr(inst, cachename)
            if len(CACHE)>1:
                raise TypeError, "Apparently you are opening a pickle in which '%s' was a method accepting arguments"%name
            if len(CACHE)==1:
                self.cache = CACHE.values()[0]
            else:
                self.cache = cache
            delattr(inst, cachename)
        else:
            self.cache = cache  # None means: the underlying method will be called
        self._instance = inst

    def __reduce__(self):
        """
        Since functions can not be pickled, the cached method caller
        is pickled by a :class:`CachedMethodPickle`, that replaces
        itself by an actual :class:`CachedMethodCallerNoArgs` as soon
        as it is asked to do anything.

        TEST::

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]
            sage: I.gens()
            [a, b]
            sage: I.gens
            Cached version of <function gens at 0x...>
            sage: J = loads(dumps(I))
            sage: J.gens
            Pickle of the cached method "gens"
            sage: J.gens.cache
            [a, b]
            sage: J.gens
            Cached version of <function gens at 0x...>

        """
        return CachedMethodPickle,(self._instance,self.__name__,self.cache)

    def _instance_call(self):
        """
        Call the cached method without using the cache.

        EXAMPLE::

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

        EXAMPLE::

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

        NOTE:

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

    cpdef clear_cache(self):
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

        NOTE:

        Recall that a cached method without arguments can not cache
        the return value ``None``.

        EXAMPLE::

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

    def __get__(self, inst, cls): #cls=None):
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
        it is attempted to find an attribute ``__cached_methods``,
        and store/retrieve the bound method there. In that
        way, cached methods can be implemented for extension
        classes deriving from :class:`~sage.structure.parent.Parent`
        and :class:`~sage.structure.element.Element`.

        TESTS:

        Due to the separate ``__get__`` method, it is possible
        to define a cached method in one class and use it as
        an attribute of another class.

            sage: class Foo:
            ...       def __init__(self, n):
            ...           self.__n = n
            ...       @cached_method
            ...       def f(self):
            ...           return self.__n^2
            ...
            sage: class Bar:
            ...       f = Foo.f
            ...
            sage: b1 = Bar()
            sage: b2 = Bar()

        The :class:`CachedMethod` is replaced by an instance of
        :class:`CachedMethodCallerNoArgs` that is set as an
        attribute. Hence, we have::

            sage: b1.f is b1.f
            True
            sage: type(b1.f)
            <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>

        Any instance of ``Bar`` gets its own instance of
        :class:`CachedMethodCaller``::

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
            return (<dict>inst.__cached_methods)[self.__name__]
        except (AttributeError,TypeError,KeyError),msg:
            pass
        Caller = CachedMethodCallerNoArgs(inst, self.f, name=self.__name__)
        try:
            setattr(inst,self.__name__, Caller)
            return Caller
        except AttributeError:
            pass
        try:
            if inst.__cached_methods is None:
                inst.__cached_methods = {self.__name__ : Caller}
            else:
                (<dict>inst.__cached_methods)[self.__name__] = Caller
        except AttributeError,msg:
            pass
        return Caller

cdef class CachedMethod(object):
    """
    A decorator that creates a cached version of an instance
    method of a class.

    NOTE:

    For proper behavior, the method must be a pure function
    (no side effects). Arguments to the method must be hashable.

    EXAMPLES::

        sage: class Foo(object):
        ...       @cached_method
        ...       def f(self, t, x=2):
        ...           print 'computing'
        ...           return t**x
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
        <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>

    So, you would hardly ever see an instance of this class alive.
    """
    def __init__(self, f, name=None):
        """
        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self,n):
            ...           return self._x^n
            ...       @cached_method
            ...       def f0(self):
            ...           return self._x^2
            ...
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
            <type 'sage.misc.cachefunc.CachedMethodCaller'>
            sage: a.f.cache is a._cache__f
            True

        Note that if the instance ``a`` would not accept attribute
        assignment, the computations would still be cached in
        ``a.f.cache``, and they would in fact be preserved when
        pickling.

        The cached method ``f0`` accepts no arguments, which allows
        for an improved way of caching: By an attribute of the cached
        method itsel. This cache is *only* available in that way, i.e.,
        it is not additionally stored as an attribute of ``a``::

            sage: type(a.f0)
            <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
            sage: a.f0.cache
            4
            sage: sorted(dir(a))
            ['__doc__', '__init__', '__module__', '_cache__f', '_x', 'f', 'f0']

        """
        self._cache_name = '_cache__' + (name or f.__name__)
        self._cachedfunc = CachedFunction(f, classmethod=True, name=name)

    def _instance_call(self, inst, *args, **kwds):
        """
        Call the cached method *without* using the cache.

        INPUT:

        - ``inst`` - an instance at which the method is to be called
        - Further positional or named arguments.

        EXAMPLES::

            sage: class Foo(object):
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self,n=2):
            ...           return self._x^n
            ...
            sage: a = Foo(2)
            sage: a.f()
            4

        Usually, a cached method is indeed cached::

            sage: a.f() is a.f()
            True

        However, when it becomes necessary, one can call it without
        using the cache. Note that ``a.f`` is an instance of
        :class:`CachedMethodCaller`.  But its
        :meth:`CachedMethodCaller._instance_call` relies on this
        method, so, we have an indirect doctest::

            sage: a.f._instance_call() is a.f() # indirect doctest
            False
            sage: a.f._instance_call() == a.f()
            True

        """
        return self._cachedfunc.f(inst, *args, **kwds)

    cpdef dict _get_instance_cache(self, inst):
        """
        Returns the cache dictionary.

        TESTS::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self,n=2):
            ...           return self._x^n
            ...
            sage: a = Foo(2)
            sage: a.f()
            4

        Note that we can not provide a direct test, since ``a.f`` is
        an instance of :class:`CachedMethodCaller`.  But during its
        initialisation, this method was called in order to provide the
        cached method caller with its cache, and, if possible, assign
        it to an attribute of ``a``.  So, the following is an indirect
        doctest::

            sage: a.f.get_cache()    # indirect doctest
            {((2,), ()): 4}
            sage: a._cache__f
            {((2,), ()): 4}

        """
        try:
            return inst.__dict__.setdefault(self._cache_name, {})
        except AttributeError:
            return {}

    def __get__(self, object inst, cls): #cls=None):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached method.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self):
            ...           return 1
            ...       @cached_method
            ...       def g(self, x,n=3):
            ...           return x^n
            ...
            sage: a = Foo()
            sage: type(a.f)
            <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
            sage: type(a.g)
            <type 'sage.misc.cachefunc.CachedMethodCaller'>

        By trac ticket #8611, it is attempted to set the
        CachedMethodCaller as an attribute of the instance ``a``,
        replacing the original cached attribute. Therefore, the
        ``__get__`` method will be used only once, which saves much
        time. Hence, we have::

            sage: a.f is a.f
            True
            sage: a.g is a.g
            True

        """
        # This is for Parents or Elements that do not allow attribute assignment:
        cdef str name
        try:
            name = self._cachedfunc.__name__
        except AttributeError:
            name = self.__name__
        try:
            return (<dict>inst.__cached_methods)[name]
        except (AttributeError,TypeError,KeyError),msg:
            pass
        # Apparently we need to construct the caller.
        # Since we have an optimized version for functions that do not accept arguments,
        # we need to analyse the argspec
        f = (<CachedFunction>self._cachedfunc).f
        if self.nargs==0:
            args, varargs, keywords, defaults = sage_getargspec(f)
            if varargs is None and keywords is None and len(args)<=1:
                self.nargs = 1
                Caller = CachedMethodCallerNoArgs(inst, f, name=name)
            else:
                self.nargs = 2 # don't need the exact number
                Caller = CachedMethodCaller(self, inst,
                                            cache=self._get_instance_cache(inst),
                                            name=name)
        elif self.nargs==1:
            Caller = CachedMethodCallerNoArgs(inst, f, name=name)
        else:
            Caller = CachedMethodCaller(self, inst,
                                        cache=self._get_instance_cache(inst),
                                        name=name)
        try:
            setattr(inst, name, Caller)
            return Caller
        except AttributeError:
            pass
        try:
            if inst.__cached_methods is None:
                inst.__cached_methods = {name : Caller}
            else:
                (<dict>inst.__cached_methods)[name] = Caller
        except AttributeError:
            pass
        return Caller

        # Note: a simpler approach to this would be
        # def caller(*args, **kwds):
        #     return self._instance_call(inst, *args, **kwds)
        # return caller
        # The disadvantage to this is that it does not provide
        # is_in_cache(), set_cache(), clear_cache(), ... methods.

cdef class CachedSpecialMethod(CachedMethod):
    """
    Cached version of *special* python methods.

    IMPLEMENTATION:

    For new style classes ``C``, it is not possible to override a special
    method, such as ``__hash__``, in the ``__dict__`` of an instance ``c`` of
    ``C``, because Python will for efficiency reasons always use what is
    provided by the class, not by the instance.

    By consequence, if ``__hash__`` would be wrapped by using
    :class:`CachedMethod`, then ``hash(c)`` will access `C.__hash__` and bind
    it to ``c``, which means that the ``__get__`` method of
    :class:`CachedMethod` will be called. But there, we assume that Python has
    already inspected ``__dict__``, and thus a :class:`CachedMethodCaller`
    will be created over and over again.

    Here, the `__get__` method will explicitly access the `__dict__`, so that
    ``hash(c)`` will rely on a single :class:`CachedMethodCaller` stored in
    the ``__dict__``.

    EXAMPLES::

        sage: class C:
        ....:     @cached_method
        ....:     def __hash__(self):
        ....:         print "compute hash"
        ....:         return int(5)
        ....:
        sage: c = C()
        sage: type(C.__hash__)
        <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>

    The hash is computed only once, subsequent calls will use the value from
    the cache. This was implemented in :trac:`12601`.

        sage: hash(c)       # indirect doctest
        compute hash
        5
        sage: hash(c)
        5

    """
    def __get__(self, object inst, cls):
        """
        Bind a :class:`CachedMethodCaller` to a specific instance, using `__dict__`.

        EXAMPLES::

            sage: class C:
            ....:     @cached_method
            ....:     def __hash__(self):
            ....:         print "compute hash"
            ....:         return int(5)
            ....:
            sage: c = C()
            sage: type(C.__hash__)
            <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
            sage: hash(c)       # indirect doctest
            compute hash
            5
            sage: hash(c)
            5
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
                    D = inst.__cached_methods
                except (TypeError, AttributeError):
                    raise TypeError("For a cached special method, either attribute assignment or a public '__cached_methods' attribute of type <dict> is needed")
            if D is None:
                # This can only happen in the case of __cached_methods
                D = inst.__cached_methods = {}
            else:
                try:
                    return D[name]
                except KeyError:
                    pass
        # Apparently we need to construct the caller.
        # Since we have an optimized version for functions that do not accept arguments,
        # we need to analyse the argspec
        f = (<CachedFunction>self._cachedfunc).f
        if self.nargs==0:
            args, varargs, keywords, defaults = sage_getargspec(f)
            if varargs is None and keywords is None and len(args)<=1:
                self.nargs = 1
                Caller = CachedMethodCallerNoArgs(inst, f, name=name)
            else:
                self.nargs = 2 # don't need the exact number
                Caller = CachedMethodCaller(self, inst,
                                            cache=self._get_instance_cache(inst),
                                            name=name)
        elif self.nargs==1:
            Caller = CachedMethodCallerNoArgs(inst, f, name=name)
        else:
            Caller = CachedMethodCaller(self, inst,
                                        cache=self._get_instance_cache(inst),
                                        name=name)
        if inst is not None:
            try:
                setattr(inst,name, Caller)
                return Caller
            except AttributeError:
                pass
            D[name] = Caller
        return Caller

def cached_method(f, name=None):
    """
    
    EXAMPLES:

    In the following examples, one can see how a cached method works in applicationy.
    Below, we demonstrate what is done behind the scenes::

        sage: class C:
        ....:     @cached_method
        ....:     def __hash__(self):
        ....:         print "compute hash"
        ....:         return int(5)
        ....:     @cached_method
        ....:     def f(self, x):
        ....:         print "computing cached method"
        ....:         return x*2
        sage: c = C()
        sage: type(C.__hash__)
        <type 'sage.misc.cachefunc.CachedMethodCallerNoArgs'>
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

    Using cached methods for the hash and other special methods was
    implemented in :trac:`12601`, by means of :class:`CachedSpecialMethod`. We
    show that it is used behind the scenes::

        sage: cached_method(c.__hash__)
        <sage.misc.cachefunc.CachedSpecialMethod object at ...>
        sage: cached_method(c.f)
        <sage.misc.cachefunc.CachedMethod object at ...>

    """
    cdef str fname = name or f.__name__
    if fname in special_method_names:
        return CachedSpecialMethod(f, name)
    return CachedMethod(f, name)

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
    - the instances are hashable (they are part of the cache key).

    NOTE:

    For proper behavior, the method must be a pure function (no side
    effects). If this decorator is used on a method, it will have
    identical output on equal elements. This is since the element is
    part of the hash key. Arguments to the method and the instance
    it is assigned to must be hashable.

    Examples can be found at :mod:`~sage.misc.cachefunc`.

    """

    def __init__(self, f, name=None):
        """
        Constructs a new method with cache stored in the parent of the instance.

        See also ``cached_method`` and ``cached_function``.

        EXAMPLES::

            sage: class MyParent(Parent):
            ...       pass
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       _parent = MyParent()
            ...       def parent(self):
            ...           return self._parent
            ...       @cached_in_parent_method  #indirect doctest
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: hasattr(a.parent(), '_cache__element_f')
            True

        For speeding up internal computations, this dictionary
        is also accessible as an attribute of the CachedMethodCaller
        (by trac ticket #8611)::

            sage: a.parent()._cache__element_f is a.f.cache
            True
        """
        self._cache_name = '_cache__' + 'element_' + (name or f.__name__)
        self._cachedfunc = CachedFunction(f, classmethod=True, name=name)

    cpdef dict _get_instance_cache(self, inst):
        """
        Returns the cache dictionary, which is stored in the parent.

        EXAMPLES::

            sage: class MyParent(Parent):
            ...       pass
            ...
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       _parent = MyParent()
            ...       def parent(self):
            ...           return self._parent
            ...       def __eq__(self, other):
            ...           return self._x^2 == other._x^2
            ...       def __hash__(self):
            ...           return hash(self._x^2)
            ...       def __repr__(self):
            ...           return 'My %s'%self._x
            ...       @cached_in_parent_method
            ...       def f(self):
            ...           return self._x^3
            ...
            sage: a = Foo(2)
            sage: a.f()
            8
            sage: a.f.get_cache()   #indirect doctest
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
            sage: c.f.get_cache() is a.f.get_cache() #indirect doctest
            True

        Note that the cache is also available as an
        attribute of the cached method, which speeds
        up internal computations::

            sage: a.f.cache is b.f.get_cache() is c.f._cachedmethod._get_instance_cache(c)
            True

        """
        if inst is None:
            return {}
        try:
            P = inst.parent()
            return P.__dict__.setdefault(self._cache_name, {})
        except AttributeError:
            pass
        if not hasattr(P,'__cached_methods'):
            raise TypeError, "The parent of this element does not allow attribute assignment\n    and does not descend from the Parent base class.\n    Can not use CachedInParentMethod."
        if P.__cached_methods is None:
            P.__cached_methods = {}
        return (<dict>P.__cached_methods).setdefault(self._cache_name, {})

    def __get__(self, inst, cls): #cls=None):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached-in-parent method.
        """
        Caller = CachedMethodCaller(self, inst, cache=self._get_instance_cache(inst), inst_in_key=True)
        try:
            setattr(inst,self._cachedfunc.__name__, Caller)
        except AttributeError:
            pass
        return Caller

cached_in_parent_method = CachedInParentMethod

class FileCache:
    """
    FileCache is a dictionary-like class which stores keys and
    values on disk.  The keys take the form of a tuple (A,K)

    - A is a tuple of objects t where each t is an exact
      object which is uniquely identified by a short string.

    - K is a tuple of tuples (s,v) where s is a valid
      variable name and v is an exact object which is uniquely
      identified by a short string with letters [a-zA-Z0-9-._]

    The primary use case is the DiskCachedFunction.  If
    ``memory_cache == True``, we maintain a cache of objects seen
    during this session in memory -- but we don't load them from
    disk until necessary.  The keys and values are stored in a
    pair of files:

    - ``prefix-argstring.key.sobj`` contains the ``key`` only,
    - ``prefix-argstring.sobj`` contains the tuple ``(key,val)``

    where ``self[key] == val``.

    NOTE:

    We assume that each FileCache lives in its own directory.
    Use **extreme** caution if you wish to break that assumption.
    """
    def __init__(self, dir, prefix = '', memory_cache = False):
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
        from sage.misc.misc import sage_makedirs
        if len(dir) == 0 or dir[-1] != '/':
            dir += '/'
        self._dir = dir
        sage_makedirs(dir)

        self._prefix = prefix + '-'

        if memory_cache:
            self._cache = {}
        else:
            self._cache = None

    def file_list(self):
        """
        Returns the list of files corresponding to self.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = True, prefix='t')
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: for f in sorted(FC.file_list()): print f[len(dir):]
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
                files.append( dir + f )
        return files

    def items(self):
        """
        Returns a list of tuples ``(k,v)`` where ``self[k] = v``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: I = FC.items()
            sage: I.sort(); print I
            [(((), ()), 1), (((1,), (('a', 1),)), 3), (((1, 2), ()), 2)]
        """
        return [(k,self[k]) for k in self]

    def values(self):
        """
        Returns a list of values that are stored in ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: FC[((),(('a',1),))] = 4
            sage: v = FC.values()
            sage: v.sort(); print v
            [1, 2, 3, 4]
        """
        return [self[k] for k in self]

    def __iter__(self):
        """
        Returns a list of keys of ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: for k in sorted(FC): print k
            ((), ())
            ((1,), (('a', 1),))
            ((1, 2), ())
        """
        return self.keys().__iter__()

    def keys(self):
        """
        Returns a list of keys ``k`` where ``self[k]`` is defined.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: K = FC.keys()
            sage: K.sort(); print K
            [((), ()), ((1,), (('a', 1),)), ((1, 2), ())]
        """
        cdef list K = []
        from sage.structure.sage_object import load
        for f in self.file_list():
            if f[-9:] == '.key.sobj':
                K.append(load(f))
        return K

    def _filename(self, key):
        """
        Computes the filename associated with a certain key.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False, prefix='foo')
            sage: N = FC._filename(((1,2), (('a',1),('b',2))))
            sage: print N[len(dir):]
            foo-a-1_b-2.1_2
            sage: N = FC._filename(((), (('a',1),('b',2))))
            sage: print N[len(dir):]
            foo-a-1_b-2
            sage: N = FC._filename(((1,2), ()))
            sage: print N[len(dir):]
            foo-1_2
        """
        a,k = key
        kwdstr = '_'.join(['%s-%s'%x for x in k])
        argstr = '_'.join(['%s'%x for x in a])
        if kwdstr and argstr:
            keystr = kwdstr + '.' + argstr
        else:
            keystr = kwdstr + argstr
        return self._dir + self._prefix + keystr

    def has_key(self, key):
        """
        Returns ``True`` if ``self[key]`` is defined and False otherwise.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: FC[k] = True
            sage: FC.has_key(k)
            True
            sage: FC.has_key(((),()))
            False
        """
        return os.path.exists(self._filename(key) + '.key.sobj')

    def __getitem__(self, key):
        """
        Returns the value set by ``self[key] = val``, in this session
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
        from sage.structure.sage_object import load

        cache = self._cache
        if cache is not None:
            if cache.has_key(key):
                return cache[key]

        f = self._filename(key) + '.sobj'
        try:
            k,v = load(f)
        except IOError:
            raise KeyError, key
        if k != key:
            raise RuntimeError, "cache corrupted"

        if cache is not None:
            cache[key] = v
        return v

    def __setitem__(self, key, value):
        """
        Sets ``self[key] = value`` and stores both key and value on
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
        from sage.structure.sage_object import save

        f = self._filename(key)

        save(key, f+'.key.sobj')
        save((key,value), f + '.sobj')
        if self._cache is not None:
            self._cache[key] = value

    def __delitem__(self, key):
        """
        Delete the key,value pair from self and unlink the associated
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
            sage: FC1.has_key(k)
            False
       """
        f = self._filename(key)
        cache = self._cache
        if cache is not None and cache.has_key(key):
            del self._cache[key]
        if os.path.exists(f + '.sobj'):
            os.remove(f + '.sobj')
        if  os.path.exists(f + '.key.sobj'):
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
    def __init__(self, f, dir, memory_cache=False):
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
        CachedFunction.__init__(self, f)
        prefix = f.__name__
        self.cache = FileCache(dir, prefix=prefix, memory_cache = memory_cache)


class disk_cached_function:
    """
    Decorator for :class:`DiskCachedFunction`.

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: @disk_cached_function(dir)
        ... def foo(x): return next_prime(2^x)%x
        sage: x = foo(200);x
        11
        sage: @disk_cached_function(dir)
        ... def foo(x): return 1/x
        sage: foo(200)
        11
        sage: foo.clear_cache()
        sage: foo(200)
        1/200
    """
    def __init__(self, dir, memory_cache = False):
        """
        EXAMPLES::

            sage: dir = tmp_dir()
            sage: @disk_cached_function(dir, memory_cache=True)
            ... def foo(x): return next_prime(2^x)
            sage: x = foo(200)
            sage: x is foo(200)
            True
            sage: @disk_cached_function(dir, memory_cache=False)
            ... def foo(x): return next_prime(2^x)
            sage: x is foo(200)
            False
        """
        self._dir = dir
        self._memory_cache = memory_cache

    def __call__(self, f):
        """
        EXAMPLES::

            sage: dir = tmp_dir()
            sage: @disk_cached_function(dir)
            ... def foo(x): return ModularSymbols(x)
            sage: foo(389)
            Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2 with sign 0 over Rational Field
        """
        return DiskCachedFunction(f, self._dir, memory_cache = self._memory_cache)

class ClearCacheOnPickle(object):
    r"""
    This class implements an appropriate ``__getstate__`` method that
    clears the cache of the methods (see @cached_method) before
    passing them on to the caller, typically the pickle and copy modules.

    The implemented ``__getstate__`` method calls the ``__getstate__``
    methods of classes later in the method resolution
    order. Therefore, classes which want this behaviour should inherit
    first from this one.

    EXAMPLE:

    In the following example, we create a Python class that inherits
    from multivariate polynomial ideals, but does not pickle cached
    values.  We provide the definition in Cython, however, since
    interactive Cython definitions provide introspection by trac
    ticket #9976, whereas Python definitions don't.
    ::

        sage: P.<a,b,c,d> = QQ[]
        sage: I = P*[a,b]
        sage: classdef = ['from sage.misc.cachefunc import ClearCacheOnPickle',
        ...    'from sage.all import QQ',
        ...    'P = QQ["a","b","c","d"]; I = P*[P.gen(0),P.gen(1)]',
        ...    'class MyClass(ClearCacheOnPickle,I.__class__):',
        ...    '    def __init__(self,ring,gens):',
        ...    '        I.__class__.__init__(self,ring,gens)',
        ...    '    def __getnewargs__(self):',
        ...    '        return (self._Ideal_generic__ring,self._Ideal_generic__gens)']
        sage: cython('\n'.join(classdef))

    We destroy the cache of two methods of ``I`` on purpose
    (demonstrating that the two different implementations of cached
    methods are correctly dealt with).  Pickling ``I`` preserves the
    cache::

        sage: I.gens.set_cache('bar')
        sage: I.groebner_basis.set_cache('foo',algorithm='singular')
        sage: J = loads(dumps(I))
        sage: J.gens()
        'bar'
        sage: J.groebner_basis(algorithm='singular')
        'foo'

    However, if we have an ideal that additionally descends from
    :class:`ClearCacheOnPickle`, the carefully corrupted cache is not
    pickled::

        sage: A = MyClass(P,[a,b])
        sage: A
        Ideal (a, b) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
        sage: A.gens.set_cache('foo')
        sage: A.groebner_basis.set_cache('bar',algorithm='singular')
        sage: A.gens()
        'foo'
        sage: A.groebner_basis(algorithm='singular')
        'bar'
        sage: B = loads(dumps(A))
        sage: B.gens()
        [a, b]
        sage: B.groebner_basis(algorithm='singular')
        [a, b]
        sage: A.gens()
        'foo'

    """
    def __getstate__(self):
        r"""
        The idea is to remove that might provide a cache to some cached method
        from the return value of the ``__getstate__`` method.

        EXAMPLE:

        In the following example, we create a Python class that
        inherits from multivariate polynomial ideals, but clears the
        cache as well.

            sage: P.<a,b,c,d> = QQ[]
            sage: I = P*[a,b]

        We destroy the cache of two methods if ``I`` on purpose
        (demonstrating that the two different implementations of cached
        methods are correctly dealt with).  Pickling ``I`` preserves the
        cache::

            sage: I.gens.set_cache('bar')
            sage: I.groebner_basis.set_cache('foo',algorithm='singular')
            sage: J = loads(dumps(I))
            sage: J.gens()
            'bar'
            sage: J.groebner_basis(algorithm='singular')
            'foo'

        However, if we do the same with a class that additionally
        descends from :class:`ClearCacheOnPickle`, the cache is not
        pickled. We provide the definition in Cython, however, since
        interactive Cython definitions provide introspection by trac
        ticket #9976, whereas Python definitions don't.
        ::

            sage: classdef = ['from sage.misc.cachefunc import ClearCacheOnPickle',
            ...    'from sage.all import QQ',
            ...    'from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal',
            ...    'class MyClass(ClearCacheOnPickle,MPolynomialIdeal):',
            ...    '    def __init__(self,ring,gens):',
            ...    '        MPolynomialIdeal.__init__(self,ring,gens)',
            ...    '    def __getnewargs__(self):',
            ...    '        return (self._Ideal_generic__ring,self._Ideal_generic__gens)']
            sage: cython('\n'.join(classdef))
            sage: A = MyClass(P,[a,b])
            sage: A
            Ideal (a, b) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
            sage: A.gens.set_cache('foo')
            sage: A.groebner_basis.set_cache('bar',algorithm='singular')
            sage: A.gens()
            'foo'
            sage: A.groebner_basis(algorithm='singular')
            'bar'
            sage: B = loads(dumps(A))
            sage: B.gens()
            [a, b]
            sage: B.groebner_basis(algorithm='singular')
            [a, b]
            sage: A.gens()
            'foo'

        And here is why the example works::

            sage: ST = I.__getstate__(); ST[0],sorted(ST[1].items())
            (Monoid of ideals of Multivariate Polynomial Ring in a, b, c, d over Rational Field, [('_Ideal_generic__gens', (a, b)), ('_Ideal_generic__ring', Multivariate Polynomial Ring in a, b, c, d over Rational Field), ('_cache__groebner_basis', {(('singular', None, None, False), ()): 'foo'}), ('_gb_by_ordering', {}), ('gens', Cached version of <function gens at 0x...>), ('groebner_basis', Cached version of <function groebner_basis at 0x...>)])
            sage: ST = A.__getstate__(); ST[0],sorted(ST[1].items())
            (Monoid of ideals of Multivariate Polynomial Ring in a, b, c, d over Rational Field, [('_Ideal_generic__gens', (a, b)), ('_Ideal_generic__ring', Multivariate Polynomial Ring in a, b, c, d over Rational Field), ('_gb_by_ordering', {})])

        """
        OrigState = super(ClearCacheOnPickle, self).__getstate__()
        def clear_list(T):
            L = []
            for x in T:
                if isinstance(x,list):
                    L.append(clear_list(x))
                elif isinstance(x,tuple):
                    L.append(clear_tuple(x))
                elif isinstance(x,dict):
                    L.append(clear_dict(x))
                elif not isinstance(x,CachedFunction):
                    L.append(x)
            return L
        def clear_tuple(T):
            return tuple(clear_list(T))
        def clear_dict(T):
            D = {}
            for key,value in T.iteritems():
                if not ((type(key) == str and key[0:8] == '_cache__') or
                            isinstance(value,CachedFunction)):
                    if isinstance(value,list):
                        D[key] = clear_list(value)
                    elif isinstance(value,tuple):
                        D[key] = clear_tuple(value)
                    elif isinstance(value,dict):
                        D[key] = clear_dict(value)
                    else:
                        D[key] = value
            return D
        if isinstance(OrigState,tuple):
            return clear_tuple(OrigState)
        if isinstance(OrigState,list):
            return clear_list(OrigState)
        if isinstance(OrigState,dict):
            return clear_dict(OrigState)
        return OrigState
