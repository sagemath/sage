# sage_setup: distribution = sagemath-objects
r"""
Unique Representation

Abstract classes for cached and unique representation behavior.

.. SEEALSO::

   :class:`sage.structure.factory.UniqueFactory`

AUTHORS:

- Nicolas M. Thiery (2008): Original version.
- Simon A. King (2013-02): Separate cached and unique representation.
- Simon A. King (2013-08): Extended documentation.


What is a cached representation?
================================

Instances of a class have a *cached representation behavior* when several
instances constructed with the same arguments share the same memory
representation. For example, calling twice::

    sage: G = SymmetricGroup(6)                                                         # needs sage.groups
    sage: H = SymmetricGroup(6)                                                         # needs sage.groups

to create the symmetric group on six elements gives back the same
object::

    sage: G is H                                                                        # needs sage.groups
    True

This is a standard design pattern. Besides saving memory, it allows for
sharing cached data (say representation theoretical information about a
group). And of course a look-up in the cache is faster than the creation of a
new object.

Implementing a cached representation
------------------------------------

Sage provides two standard ways to create a cached representation:
:class:`CachedRepresentation` and
:class:`~sage.structure.factory.UniqueFactory`. Note that, in spite of its
name, :class:`~sage.structure.factory.UniqueFactory` does not ensure *unique*
representation behaviour, which will be explained below.

Using :class:`CachedRepresentation`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is often very easy to use :class:`CachedRepresentation`: One simply writes
a Python class and adds :class:`CachedRepresentation` to the list of base
classes. If one does so, then the arguments used to create an instance of this
class will by default also be used as keys for the cache::

    sage: from sage.structure.unique_representation import CachedRepresentation
    sage: class C(CachedRepresentation):
    ....:     def __init__(self, a, b=0):
    ....:         self.a = a
    ....:         self.b = b
    ....:     def __repr__(self):
    ....:         return "C(%s, %s)"%(self.a, self.b)
    sage: a = C(1)
    sage: a is C(1)
    True

In addition, pickling just works, provided that Python is able to look up the
class. Hence, in the following two lines, we explicitly put the class into the
``__main__`` module. This is needed in doctests, but not in an interactive
session::

    sage: import __main__
    sage: __main__.C = C
    sage: loads(dumps(a)) is a
    True

Often, this very easy approach is sufficient for applications. However, there
are some pitfalls. Since the arguments are used for caching, all arguments
must be hashable, i.e., must be valid as dictionary keys::

    sage: C((1,2))
    C((1, 2), 0)
    sage: C([1,2])
    Traceback (most recent call last):
    ...
    TypeError: unhashable type: 'list'

In addition, equivalent ways of providing the arguments are *not*
automatically normalised when forming the cache key, and hence different but
equivalent arguments may yield distinct instances::

    sage: C(1) is C(1,0)
    False
    sage: C(1) is C(a=1)
    False
    sage: repr(C(1)) == repr(C(a=1))
    True

It should also be noted that the arguments are compared by equality, not by
identity. This is often desired, but can imply subtle problems. For example,
since ``C(1)`` already is in the cache, and since the unit elements in
different finite fields are all equal to the integer one, we find::

    sage: GF(5)(1) == 1 == GF(3)(1)
    True
    sage: C(1) is C(GF(3)(1)) is C(GF(5)(1))
    True

But ``C(2)`` is not in the cache, and the number two is not equal in different
finite fields (i. e., ``GF(5)(2) == GF(3)(2)`` returns as ``False``), even
though it is equal to the number two in the ring of integers (
``GF(5)(2) == 2 == GF(3)(2)`` returns as ``True``; equality is not transitive
when comparing elements of *distinct* algebraic structures!!). Hence, we
have::

    sage: GF(5)(2) == GF(3)(2)
    False
    sage: C(GF(3)(2)) is C(GF(5)(2))
    False

Normalising the arguments
.........................

:class:`CachedRepresentation` uses the metaclass
:class:`~sage.misc.classcall_metaclass.ClasscallMetaclass`. Its
``__classcall__`` method is a
:class:`~sage.misc.cachefunc.WeakCachedFunction`.  This function creates an
instance of the given class using the given arguments, unless it finds the
result in the cache. This has the following implications:

- The arguments must be valid dictionary keys (i.e., they must be hashable;
  see above).
- It is a weak cache, hence, if the user does not keep a reference to the
  resulting instance, then it may be removed from the cache during garbage
  collection.
- It is possible to preprocess the input arguments by implementing a
  ``__classcall__`` or a ``__classcall_private__`` method, but in order to
  benefit from caching, :meth:`CachedRepresentation.__classcall__` should at
  some point be called.

.. NOTE::

    For technical reasons, it is needed that ``__classcall__`` respectively
    ``__classcall_private__`` are "static methods", i.e., they are callable
    objects that do not bind to an instance or class. For example, a
    :class:`~sage.misc.cachefunc.cached_function` can be used here, because it
    is callable, but does not bind to an instance or class, because it has no
    ``__get__()`` method. A usual Python function, however, has a
    ``__get__()`` method and would thus under normal circumstances bind to an
    instance or class, and thus the instance or class would be passed to the
    function as the first argument. To prevent a callable object from being
    bound to the instance or class, one can prepend the ``@staticmethod``
    decorator to the definition; see :class:`staticmethod`.

    For more on Python's ``__get__()`` method, see:
    https://docs.python.org/2/howto/descriptor.html

.. WARNING::

    If there is preprocessing, then the preprocessed arguments
    passed to :meth:`CachedRepresentation.__classcall__` must be invariant
    under the preprocessing. That is to say, preprocessing the input
    arguments twice must have the same effect as preprocessing the input
    arguments only once. That is to say, the preprocessing must be idempotent.

The reason for this warning lies in the way pickling is implemented. If the
preprocessed arguments are passed to
:meth:`CachedRepresentation.__classcall__`, then the resulting instance will
store the *preprocessed* arguments in some attribute, and will use them for
pickling. If the pickle is unpickled, then preprocessing is applied to the
preprocessed arguments---and this second round of preprocessing must not
change the arguments further, since otherwise a different instance would be
created.

We illustrate the warning by an example. Imagine that one has instances that
are created with an integer-valued argument, but only depend on the *square*
of the argument. It would be a mistake to square the given argument during
preprocessing::

    sage: class WrongUsage(CachedRepresentation):
    ....:     @staticmethod
    ....:     def __classcall__(cls, n):
    ....:         return super().__classcall__(cls, n^2)
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:     def __repr__(self):
    ....:         return "Something(%d)"%self.n
    sage: import __main__
    sage: __main__.WrongUsage = WrongUsage # This is only needed in doctests
    sage: w = WrongUsage(3); w
    Something(9)
    sage: w._reduction
    (<class '__main__.WrongUsage'>, (9,), {})

Indeed, the reduction data are obtained from the preprocessed argument. By
consequence, if the resulting instance is pickled and unpickled, the argument
gets squared *again*::

    sage: loads(dumps(w))
    Something(81)

Instead, the preprocessing should only take the absolute value of the given
argument, while the squaring should happen inside of the ``__init__`` method,
where it won't mess with the cache::

    sage: class BetterUsage(CachedRepresentation):
    ....:     @staticmethod
    ....:     def __classcall__(cls, n):
    ....:         return super().__classcall__(cls, abs(n))
    ....:     def __init__(self, n):
    ....:         self.n = n^2
    ....:     def __repr__(self):
    ....:         return "SomethingElse(%d)"%self.n
    sage: __main__.BetterUsage = BetterUsage # This is only needed in doctests
    sage: b = BetterUsage(3); b
    SomethingElse(9)
    sage: loads(dumps(b)) is b
    True
    sage: b is BetterUsage(-3)
    True

In our next example, we create a cached representation class ``C`` that
returns an instance of a sub-class ``C1`` or ``C2`` depending on the given
arguments. This is implemented in a static ``__classcall_private__`` method of
``C``, letting it choose the sub-class according to the given arguments. Since
a ``__classcall_private__`` method will be ignored on sub-classes, the caching
of :class:`CachedRepresentation` is available to both ``C1`` and ``C2``. But
for illustration, we overload the static ``__classcall__`` method on ``C2``,
doing some argument preprocessing. We also create a sub-class ``C2b`` of
``C2``, demonstrating that the ``__classcall__`` method is used on the
sub-class (in contrast to a ``__classcall_private__`` method!).  ::

    sage: class C(CachedRepresentation):
    ....:     @staticmethod
    ....:     def __classcall_private__(cls, n, implementation=0):
    ....:         if not implementation:
    ....:             return C.__classcall__(cls, n)
    ....:         if implementation==1:
    ....:             return C1(n)
    ....:         if implementation>1:
    ....:             return C2(n,implementation)
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:     def __repr__(self):
    ....:         return "C(%d, 0)"%self.n
    sage: class C1(C):
    ....:     def __repr__(self):
    ....:         return "C1(%d)"%self.n
    sage: class C2(C):
    ....:     @staticmethod
    ....:     def __classcall__(cls, n, implementation=0):
    ....:         if implementation:
    ....:             return super().__classcall__(cls, (n,)*implementation)
    ....:         return super().__classcall__(cls, n)
    ....:     def __init__(self, t):
    ....:         self.t = t
    ....:     def __repr__(self):
    ....:         return "C2(%s)"%repr(self.t)
    sage: class C2b(C2):
    ....:     def __repr__(self):
    ....:         return "C2b(%s)"%repr(self.t)
    sage: __main__.C2 = C2      # not needed in an interactive session
    sage: __main__.C2b = C2b

In the above example, ``C`` drops the argument ``implementation`` if it
evaluates to ``False``, and since the cached ``__classcall__`` is called in
this case, we have::

    sage: C(1)
    C(1, 0)
    sage: C(1) is C(1,0)
    True
    sage: C(1) is C(1,0) is C(1,None) is C(1,[])
    True

(Note that we were able to bypass the issue of arguments having to be
hashable by catching the empty list ``[]`` during preprocessing in the
``__classcall_private__`` method. Similarly, unhashable arguments can
be made hashable -- e. g., lists normalized to tuples -- in the
``__classcall_private__`` method before they are further delegated to
``__classcall__``. See
:class:`~sage.combinat.crystals.elementary_crystals.TCrystal` for an
example.)

If we call ``C1`` directly or if we provide ``implementation=1`` to ``C``, we
obtain an instance of ``C1``. Since it uses the ``__classcall__`` method
inherited from :class:`CachedRepresentation`, the resulting instances are
cached::

    sage: C1(2)
    C1(2)
    sage: C(2, implementation=1)
    C1(2)
    sage: C(2, implementation=1) is C1(2)
    True

The class ``C2`` preprocesses the input arguments. Instances can, again, be
obtained directly or by calling ``C``::

    sage: C(1, implementation=3)
    C2((1, 1, 1))
    sage: C(1, implementation=3) is C2(1,3)
    True

The argument preprocessing of ``C2`` is inherited by ``C2b``, since
``__classcall__`` and not ``__classcall_private__`` is used. Pickling works,
since the preprocessing of arguments is idempotent::

    sage: c2b = C2b(2,3); c2b
    C2b((2, 2, 2))
    sage: loads(dumps(c2b)) is c2b
    True

Using :class:`~sage.structure.factory.UniqueFactory`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For creating a cached representation using a factory, one has to

- create a class *separately* from the factory. This class **must** inherit
  from :class:`object`. Its instances **must** allow attribute assignment.
- write a method ``create_key`` (or ``create_key_and_extra_args``) that
  creates the cache key from the given arguments.
- write a method ``create_object`` that creates an instance of the class
  from a given cache key.
- create an instance of the factory with a name that allows to conclude where
  it is defined.

An example::

    sage: class C():
    ....:     def __init__(self, t):
    ....:         self.t = t
    ....:     def __repr__(self):
    ....:         return "C%s"%repr(self.t)
    sage: from sage.structure.factory import UniqueFactory
    sage: class MyFactory(UniqueFactory):
    ....:     def create_key(self, n, m=None):
    ....:         if isinstance(n, (tuple,list)) and m is None:
    ....:             return tuple(n)
    ....:         return (n,)*m
    ....:     def create_object(self, version, key, **extra_args):
    ....:         # We ignore version and extra_args
    ....:         return C(key)

Now, we define an instance of the factory, stating that it can be found under
the name ``'F'`` in the ``__main__`` module. By consequence, pickling works::

    sage: F = MyFactory("__main__.F")
    sage: __main__.F = F                # not needed in an interactive session
    sage: loads(dumps(F)) is F
    True

We can now create *cached* instances of ``C`` by calling the factory. The
cache only takes into account the key computed with the method ``create_key``
that we provided. Hence, different given arguments may result in the same
instance. Note that, again, the cache is weak, hence, the instance might be
removed from the cache during garbage collection, unless an external reference
is preserved.
::

    sage: a = F(1, 2); a
    C(1, 1)
    sage: a is F((1,1))
    True

**If** the class of the returned instances is a sub-class of :class:`object`,
and **if** the resulting instance allows attribute assignment, then pickling
of the resulting instances is automatically provided for, and respects the
cache.  ::

    sage: loads(dumps(a)) is a
    True

This is because an attribute is stored that explains how the instance was
created::

    sage: a._factory_data
    (<__main__.MyFactory object at ...>, (...), (1, 1), {})

.. NOTE::

    If a class is used that does not inherit from :class:`object` then unique
    pickling is *not* provided.

Caching is only available if the factory is called. If an instance of the
class is directly created, then the cache is not used::

    sage: C((1,1))
    C(1, 1)
    sage: C((1,1)) is a
    False

Comparing the two ways of implementing a cached representation
--------------------------------------------------------------

In this sub-section, we discuss advantages and disadvantages of the two ways
of implementing a cached representation, depending on the type of application.

Simplicity and transparency
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many cases, turning a class into a cached representation requires nothing
more than adding :class:`CachedRepresentation` to the list of base classes of
this class. This is, of course, a very easy and convenient way. Writing a
factory would involve a lot more work.

If preprocessing of the arguments is needed, then we have seen how to do this
by a ``__classcall_private__`` or ``__classcall__`` method. But these are
double underscore methods and hence, for example, invisible in the
automatically created reference manual. Moreover, preprocessing *and* caching
are implemented in the same method, which might be confusing. In a unique
factory, these two tasks are cleanly implemented in two separate methods.
With a factory, it is possible to create the resulting instance by arguments
that are different from the key used for caching. This is significantly
restricted with CachedRepresentation due to the requirement that argument
preprocessing be idempotent.

Hence, if advanced preprocessing is needed, then
:class:`~sage.structure.factory.UniqueFactory` might be easier and more
transparent to use than :class:`CachedRepresentation`.

Class inheritance
^^^^^^^^^^^^^^^^^

Using :class:`CachedRepresentation` has the advantage that one has a class and
creates cached instances of this class by the usual Python syntax::

    sage: G = SymmetricGroup(6)                                                                     # needs sage.groups
    sage: issubclass(SymmetricGroup, sage.structure.unique_representation.CachedRepresentation)     # needs sage.groups
    True
    sage: isinstance(G, SymmetricGroup)                                                             # needs sage.groups
    True

In contrast, a factory is just a callable object that returns something that
has absolutely nothing to do with the factory, and may in fact return
instances of quite different classes::

    sage: isinstance(GF, sage.structure.factory.UniqueFactory)
    True
    sage: K5 = GF(5)
    sage: type(K5)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>

    sage: # needs sage.rings.finite_rings
    sage: K25 = GF(25, 'x')
    sage: type(K25)                                                                     # needs sage.libs.linbox
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
    sage: Kp = GF(next_prime_power(1000000)^2, 'x')
    sage: type(Kp)
    <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>

This can be confusing to the user. Namely, the user might determine the class
of an instance and try to create further instances by calling the class rather
than the factory---which is a mistake since it works around the cache (and
also since the class might be more restrictive than the factory -- i. e., the
type of ``K5`` in the above doctest cannot be called on a prime power which
is not a prime). This mistake can more easily be avoided by using
:class:`CachedRepresentation`.

We have seen above that one can easily create new cached-representation
classes by subclassing an existing cached-representation class, even making
use of an existing argument preprocess. This would be much more complicated
with a factory. Namely, one would need to rewrite old factories making them
aware of the new classes, and/or write new factories for the new classes.

Python versus extension classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`CachedRepresentation` uses a metaclass, namely
:class:`~sage.misc.classcall_metaclass.ClasscallMetaclass`. Hence, it can
currently not be a Cython extension class. Moreover, it is supposed to be used
by providing it as a base class. But in typical applications, one also has
another base class, say, :class:`~sage.structure.parent.Parent`. Hence, one
would like to create a class with at least two base classes, which is
currently impossible in Cython extension classes.

In other words, when using :class:`CachedRepresentation`, one must work with
Python classes. These can be defined in Cython code (``.pyx`` files) and can
thus benefit from Cython's speed inside of their methods, but they must not be
``cdef class`` and can thus not use ``cdef`` attributes or methods.

Such restrictions do not exist when using a factory. However, if attribute
assignment does not work, then the automatic pickling provided by
:class:`~sage.structure.factory.UniqueFactory` will not be available.

What is a unique representation?
================================

Instances of a class have a *unique instance behavior* when instances of this
class evaluate equal if and only if they are identical. Sage provides the base
class :class:`~sage.misc.fast_methods.WithEqualityById`, which provides
comparison by identity and a hash that is determined by the memory address of
the instance. Both the equality test and the hash are implemented in Cython
and are very fast, even when one has a Python class inheriting from
:class:`~sage.misc.fast_methods.WithEqualityById`.

In many applications, one wants to combine unique instance and cached
representation behaviour. This is called *unique representation* behaviour.
We have seen above that symmetric groups have a *cached* representation
behaviour. However, they do not show the *unique* representation behaviour,
since they are equal to groups created in a totally different way, namely to
subgroups::

    sage: # needs sage.groups
    sage: G = SymmetricGroup(6)
    sage: G3 = G.subgroup([G((1,2,3,4,5,6)), G((1,2))])
    sage: G is G3
    False
    sage: type(G) == type(G3)
    False
    sage: G == G3
    True

The unique representation behaviour can conveniently be implemented with a
class that inherits from :class:`UniqueRepresentation`: By adding
:class:`UniqueRepresentation` to the base classes, the class will
simultaneously inherit from :class:`CachedRepresentation` and from
:class:`~sage.misc.fast_methods.WithEqualityById`.

If the class is also a :class:`~sage.structure.parent.Parent`,
then we says it satisfies the *unique parent* condition.

For example, a symmetric function algebra is uniquely determined by the base
ring. Thus, it is reasonable to use :class:`UniqueRepresentation` in this
case::

    sage: isinstance(SymmetricFunctions(CC), SymmetricFunctions)                        # needs sage.combinat
    True
    sage: issubclass(SymmetricFunctions, UniqueRepresentation)                          # needs sage.combinat
    True

:class:`UniqueRepresentation` differs from :class:`CachedRepresentation` only
by adding :class:`~sage.misc.fast_methods.WithEqualityById` as a base
class. Hence, the above examples of argument preprocessing work for
:class:`UniqueRepresentation` as well.

Note that a cached representation created with
:class:`~sage.structure.factory.UniqueFactory` does *not* automatically
provide unique representation behaviour, in spite of its name! Hence, for
unique representation behaviour, one has to implement hash and equality test
accordingly, for example by inheriting from
:class:`~sage.misc.fast_methods.WithEqualityById`.
"""
# ****************************************************************************
#  Copyright (C) 2008 Nicolas M. Thiery <nthiery at users.sf.net>
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import weak_cached_function
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.fast_methods import WithEqualityById


class WithPicklingByInitArgs(metaclass=ClasscallMetaclass):
    r"""
    Classes derived from :class:`WithPicklingByInitArgs` store the arguments
    passed to :meth:`__init__` to implement pickling.

    This class is for objects that are semantically immutable and determined
    by the class and the arguments passed to :meth:`__init__`.
    The class also provides implementations of :meth:`__copy__` and
    :func:`__deepcopy__`, which simply return the object.
    """

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        Construct a new object of this class and store the arguments passed to ``__init__``.

        TESTS::

            sage: from sage.structure.unique_representation import WithPicklingByInitArgs
            sage: class MyClass(WithPicklingByInitArgs):
            ....:     def __init__(self, value):
            ....:         self.value = value
            ....:     def __eq__(self, other):
            ....:         if type(self) != type(other):
            ....:             return False
            ....:         return self.value == other.value
            sage: import __main__
            sage: __main__.MyClass = MyClass  # This is only needed in doctests
            sage: x = MyClass(1)
            sage: x == loads(dumps(x))
            True
            sage: y = MyClass(1)
            sage: x is y                # No Cached/UniqueRepresentation behavior
            False
        """
        instance = typecall(cls, *args, **options)
        assert isinstance(instance, cls)
        if instance.__class__.__reduce__ == WithPicklingByInitArgs.__reduce__:
            instance._reduction = (cls, args, options)
        return instance

    def __reduce__(self):
        """
        Return the arguments that have been passed to
        :meth:`__new__<object.__new__>` to construct this object,
        as per the pickle protocol.

        See also :class:`CachedRepresentation` and
        :class:`UniqueRepresentation` for a discussion.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: x.__reduce__()          # indirect doctest
            (<function unreduce at ...>, (<class 'sage.structure.unique_representation.UniqueRepresentation'>, (), {}))
        """
        return (unreduce, self._reduction)

    def __copy__(self):
        """
        Return ``self``, as a semantic copy of ``self``.

        This assumes that the object is semantically immutable.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: x is copy(x)    # indirect doctest
            True
        """
        return self

    def __deepcopy__(self, memo):
        """
        Return ``self``, as a semantic deep copy of ``self``.

        This assumes that the object is semantically immutable.

        EXAMPLES::

            sage: from copy import deepcopy
            sage: x = UniqueRepresentation()
            sage: x is deepcopy(x)      # indirect doctest
            True
        """
        return self


def unreduce(cls, args, keywords):
    """
    Calls a class on the given arguments::

        sage: sage.structure.unique_representation.unreduce(Integer, (1,), {})
        1

    .. TODO::

        should reuse something preexisting ...

    """
    return cls(*args, **keywords)


class CachedRepresentation(WithPicklingByInitArgs):
    """
    Classes derived from CachedRepresentation inherit a weak cache for their
    instances.

    .. NOTE::

        If this class is used as a base class, then instances are (weakly)
        cached, according to the arguments used to create the instance.
        Pickling is provided, of course by using the cache.

    .. NOTE::

        Using this class, one can have arbitrary hash and comparison.
        Hence, *unique* representation behaviour is *not* provided.

    .. SEEALSO::

        :class:`UniqueRepresentation`, :mod:`~sage.structure.unique_representation`

    EXAMPLES:

    Providing a class with a weak cache for the instances is easy: Just
    inherit from :class:`CachedRepresentation`::

        sage: from sage.structure.unique_representation import CachedRepresentation
        sage: class MyClass(CachedRepresentation):
        ....:     # all the rest as usual
        ....:     pass

    We start with a simple class whose constructor takes a single
    value as argument (TODO: find a more meaningful example)::

        sage: class MyClass(CachedRepresentation):
        ....:     def __init__(self, value):
        ....:         self.value = value
        ....:     def __eq__(self, other):
        ....:         if type(self) != type(other):
        ....:             return False
        ....:         return self.value == other.value

    Two coexisting instances of ``MyClass`` created with the same argument data
    are guaranteed to share the same identity. Since :issue:`12215`, this is
    only the case if there is some strong reference to the returned instance,
    since otherwise it may be garbage collected::

        sage: x = MyClass(1)
        sage: y = MyClass(1)
        sage: x is y               # There is a strong reference
        True
        sage: z = MyClass(2)
        sage: x is z
        False

    In particular, modifying any one of them modifies the other
    (reference effect)::

        sage: x.value = 3
        sage: x.value, y.value
        (3, 3)
        sage: y.value = 1
        sage: x.value, y.value
        (1, 1)

    The arguments can consist of any combination of positional or keyword
    arguments, as taken by a usual :meth:`__init__ <object.__init__>`
    function. However, all values passed in should be hashable::

        sage: MyClass(value = [1,2,3])
        Traceback (most recent call last):
        ...
        TypeError: unhashable type: 'list'

    .. rubric:: Argument preprocessing

    Sometimes, one wants to do some preprocessing on the arguments, to
    put them in some canonical form. The following example illustrates
    how to achieve this; it takes as argument any iterable, and
    canonicalizes it into a tuple (which is hashable!)::

        sage: class MyClass2(CachedRepresentation):
        ....:     @staticmethod
        ....:     def __classcall__(cls, iterable):
        ....:         t = tuple(iterable)
        ....:         return super().__classcall__(cls, t)
        ....:
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: x = MyClass2([1,2,3])
        sage: y = MyClass2(tuple([1,2,3]))
        sage: z = MyClass2(i for i in [1,2,3])
        sage: x.value
        (1, 2, 3)
        sage: x is y, y is z
        (True, True)

    A similar situation arises when the constructor accepts default
    values for some of its parameters. Alas, the obvious
    implementation does not work::

        sage: class MyClass3(CachedRepresentation):
        ....:     def __init__(self, value=3):
        ....:         self.value = value
        sage: MyClass3(3) is MyClass3()
        False

    Instead, one should do::

        sage: class MyClass3(UniqueRepresentation):
        ....:     @staticmethod
        ....:     def __classcall__(cls, value=3):
        ....:         return super().__classcall__(cls, value)
        ....:
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: MyClass3(3) is MyClass3()
        True

    A bit of explanation is in order. First, the call ``MyClass2([1,2,3])``
    triggers a call to ``MyClass2.__classcall__(MyClass2, [1,2,3])``. This is
    an extension of the standard Python behavior, needed by
    :class:`CachedRepresentation`, and implemented by the
    :class:`~sage.misc.classcall_metaclass.ClasscallMetaclass`. Then,
    ``MyClass2.__classcall__`` does the desired transformations on the
    arguments. Finally, it uses ``super`` to call the default implementation
    of ``__classcall__`` provided by :class:`CachedRepresentation`. This one
    in turn handles the caching and, if needed, constructs and initializes a
    new object in the class using :meth:`__new__<object.__new__>` and
    :meth:`__init__<object.__init__>` as usual.

    Constraints:

    - :meth:`__classcall__` is a staticmethod (like, implicitly,
      :meth:`__new__<object.__new__>`)
    - the preprocessing on the arguments should be idempotent. That is, if
      ``MyClass2.__classcall__(<arguments>)`` calls
      ``CachedRepresentation.__classcall__(<preprocessed_arguments>)``, then
      ``MyClass2.__classcall__(<preprocessed_arguments>)`` should also result
      in a call to ``CachedRepresentation.__classcall__(<preprocessed_arguments>)``.
    - ``MyClass2.__classcall__`` should return the result of
      :meth:`CachedRepresentation.__classcall__` without modifying it.

    Other than that ``MyClass2.__classcall__`` may play any tricks, like
    acting as a factory and returning objects from other classes.

    .. WARNING::

        It is possible, but strongly discouraged, to let the ``__classcall__``
        method of a class ``C`` return objects that are not instances of
        ``C``. Of course, instances of a *subclass* of ``C`` are fine. Compare
        the examples in :mod:`~sage.structure.unique_representation`.

    We illustrate what is meant by an "idempotent" preprocessing. Imagine
    that one has instances that are created with an integer-valued argument,
    but only depend on the *square* of the argument. It would be a mistake to
    square the given argument during preprocessing::

        sage: class WrongUsage(CachedRepresentation):
        ....:     @staticmethod
        ....:     def __classcall__(cls, n):
        ....:         return super().__classcall__(cls, n^2)
        ....:     def __init__(self, n):
        ....:         self.n = n
        ....:     def __repr__(self):
        ....:         return "Something(%d)"%self.n
        sage: import __main__
        sage: __main__.WrongUsage = WrongUsage # This is only needed in doctests
        sage: w = WrongUsage(3); w
        Something(9)
        sage: w._reduction
        (<class '__main__.WrongUsage'>, (9,), {})

    Indeed, the reduction data are obtained from the preprocessed
    arguments. By consequence, if the resulting instance is pickled and
    unpickled, the argument gets squared *again*::

        sage: loads(dumps(w))
        Something(81)

    Instead, the preprocessing should only take the absolute value of the
    given argument, while the squaring should happen inside of the
    ``__init__`` method, where it won't mess with the cache::

        sage: class BetterUsage(CachedRepresentation):
        ....:     @staticmethod
        ....:     def __classcall__(cls, n):
        ....:         return super().__classcall__(cls, abs(n))
        ....:     def __init__(self, n):
        ....:         self.n = n^2
        ....:     def __repr__(self):
        ....:         return "SomethingElse(%d)"%self.n
        sage: __main__.BetterUsage = BetterUsage # This is only needed in doctests
        sage: b = BetterUsage(3); b
        SomethingElse(9)
        sage: loads(dumps(b)) is b
        True
        sage: b is BetterUsage(-3)
        True

    .. rubric:: Cached representation and mutability

    :class:`CachedRepresentation` is primarily intended for implementing
    objects which are (at least semantically) immutable. This is in
    particular assumed by the default implementations of ``copy`` and
    ``deepcopy``::

        sage: copy(x) is x
        True
        sage: from copy import deepcopy
        sage: deepcopy(x) is x
        True

    However, in contrast to :class:`UniqueRepresentation`, using
    :class:`CachedRepresentation` allows for a comparison that is not by
    identity::

        sage: t = MyClass(3)
        sage: z = MyClass(2)
        sage: t.value = 2

    Now ``t`` and ``z`` are non-identical, but equal::

        sage: t.value == z.value
        True
        sage: t == z
        True
        sage: t is z
        False

    .. rubric:: More on cached representation and identity

    :class:`CachedRepresentation` is implemented by means of a cache.
    This cache uses weak references in general, but strong references to
    the most recently created objects. Hence, when all other references
    to, say, ``MyClass(1)`` have been deleted, the instance is
    eventually deleted from memory (after enough other objects have been
    created to remove the strong reference to ``MyClass(1)``). A later
    call to ``MyClass(1)`` reconstructs the instance from scratch::

        sage: class SomeClass(UniqueRepresentation):
        ....:     def __init__(self, i):
        ....:         print("creating new instance for argument %s" % i)
        ....:         self.i = i
        ....:     def __del__(self):
        ....:         print("deleting instance for argument %s" % self.i)
        sage: class OtherClass(UniqueRepresentation):
        ....:     def __init__(self, i):
        ....:         pass
        sage: O = SomeClass(1)
        creating new instance for argument 1
        sage: O is SomeClass(1)
        True
        sage: O is SomeClass(2)
        creating new instance for argument 2
        False
        sage: L = [OtherClass(i) for i in range(200)]
        deleting instance for argument 2
        sage: del O
        deleting instance for argument 1
        sage: O = SomeClass(1)
        creating new instance for argument 1
        sage: del O
        sage: del L
        sage: L = [OtherClass(i) for i in range(200)]
        deleting instance for argument 1

    .. rubric:: Cached representation and pickling

    The default Python pickling implementation (by reconstructing an object
    from its class and dictionary, see "The pickle protocol" in the Python
    Library Reference) does not preserve cached representation, as Python has
    no chance to know whether and where the same object already exists.

    :class:`CachedRepresentation` tries to ensure appropriate pickling by
    implementing a :meth:`__reduce__ <object.__reduce__>` method returning the
    arguments passed to the constructor::

        sage: import __main__             # Fake MyClass being defined in a python module
        sage: __main__.MyClass = MyClass
        sage: x = MyClass(1)
        sage: loads(dumps(x)) is x
        True

    :class:`CachedRepresentation` uses the :meth:`__reduce__
    <object.__reduce__>` pickle protocol rather than :meth:`__getnewargs__
    <object.__getnewargs__>` because the latter does not handle keyword
    arguments::

        sage: x = MyClass(value = 1)
        sage: x.__reduce__()
        (<function unreduce at ...>, (<class '__main__.MyClass'>, (), {'value': 1}))
        sage: x is loads(dumps(x))
        True

    .. NOTE::

        The default implementation of :meth:`__reduce__ <object.__reduce__>`
        in :class:`CachedRepresentation` requires to store the constructor's
        arguments in the instance dictionary upon construction::

            sage: x.__dict__
            {'_reduction': (<class '__main__.MyClass'>, (), {'value': 1}), 'value': 1}

        It is often easy in a derived subclass to reconstruct the constructor's
        arguments from the instance data structure. When this is the case,
        :meth:`__reduce__ <object.__reduce__>` should be overridden; automagically
        the arguments won't be stored anymore::

            sage: class MyClass3(UniqueRepresentation):
            ....:     def __init__(self, value):
            ....:         self.value = value
            ....:
            ....:     def __reduce__(self):
            ....:         return (MyClass3, (self.value,))
            sage: import __main__; __main__.MyClass3 = MyClass3  # Fake MyClass3 being defined in a python module
            sage: x = MyClass3(1)
            sage: loads(dumps(x)) is x
            True
            sage: x.__dict__
            {'value': 1}

    .. rubric:: Migrating classes to ``CachedRepresentation`` and unpickling

    We check that, when migrating a class to :class:`CachedRepresentation`,
    older pickles can still be reasonably unpickled. Let us create a
    (new style) class, and pickle one of its instances::

        sage: class MyClass4():
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: pickle = dumps(MyClass4(1))

    It can be unpickled::

        sage: y = loads(pickle)
        sage: y.value
        1

    Now, we upgrade the class to derive from :class:`UniqueRepresentation`,
    which inherits from :class:`CachedRepresentation`::

        sage: class MyClass4(UniqueRepresentation, object):
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: __main__.MyClass4 = MyClass4

    The pickle can still be unpickled::

        sage: y = loads(pickle)
        sage: y.value
        1

    Note however that, for the reasons explained above, unique
    representation is not guaranteed in this case::

        sage: y is MyClass4(1)
        False

    .. TODO::

        Illustrate how this can be fixed on a case by case basis.

    Now, we redo the same test for a class deriving from SageObject::

        sage: class MyClass4(SageObject):
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: pickle = dumps(MyClass4(1))

        sage: class MyClass4(UniqueRepresentation, SageObject):
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: __main__.MyClass4 = MyClass4
        sage: y = loads(pickle)
        sage: y.value
        1

    Caveat: unpickling instances of a formerly old-style class is not supported yet by default::

        sage: class MyClass4:
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: pickle = dumps(MyClass4(1))

        sage: class MyClass4(UniqueRepresentation, SageObject):
        ....:     def __init__(self, value):
        ....:         self.value = value
        sage: __main__.MyClass4 = MyClass4
        sage: y = loads(pickle)  # todo: not implemented
        sage: y.value            # todo: not implemented
        1

    .. rubric:: Rationale for the current implementation

    :class:`CachedRepresentation` and derived classes use the
    :class:`~sage.misc.classcall_metaclass.ClasscallMetaclass`
    of the standard Python type. The following example explains why.

    We define a variant of ``MyClass`` where the calls to
    :meth:`__init__<object.__init__>` are traced::

        sage: class MyClass(CachedRepresentation):
        ....:     def __init__(self, value):
        ....:         print("initializing object")
        ....:         self.value = value

    Let us create an object twice::

        sage: x = MyClass(1)
        initializing object
        sage: z = MyClass(1)

    As desired the :meth:`__init__<object.__init__>` method was only called
    the first time, which is an important feature.

    As far as we can tell, this is not achievable while just using
    :meth:`__new__<object.__new__>` and :meth:`__init__<object.__init__>` (as
    defined by type; see Section :python:`Basic Customization
    <reference/datamodel.html#basic-customization>` in the Python Reference
    Manual). Indeed, :meth:`__init__<object.__init__>` is called
    systematically on the result of :meth:`__new__<object.__new__>` whenever
    the result is an instance of the class.

    Another difficulty is that argument preprocessing (as in the example
    above) cannot be handled by :meth:`__new__<object.__new__>`, since the
    unprocessed arguments will be passed down to
    :meth:`__init__<object.__init__>`.
    """

    @weak_cached_function(cache=128)  # automatically a staticmethod
    def __classcall__(cls, *args, **options):
        """
        Construct a new object of this class or reuse an existing one.

        See also :class:`CachedRepresentation` and
        :class:`UniqueRepresentation` for a discussion.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: y = UniqueRepresentation()
            sage: x is y   # indirect doctest
            True
        """
        return super().__classcall__(cls, *args, **options)

    @classmethod
    def _clear_cache_(cls):
        """
        Remove all instances of this class from the cache.

        EXAMPLES:

        If ``cls`` overloads :meth:`~sage.structure.unique_representation.CachedRepresentation.__classcall__`,
        clearing the cache still works, because ``cls.mro()``
        is searched until a ``__classcall__`` with an attribute
        ``cache`` is found::

            sage: class A(UniqueRepresentation):
            ....:     def __init__(self, x):
            ....:         pass
            sage: class B(A):
            ....:     @staticmethod
            ....:     def __classcall__(cls, *args, **kwds):
            ....:          return super().__classcall__(cls, *args, **kwds)
            sage: class C(B): pass
            sage: a = A(1)
            sage: b = B(2)
            sage: c = C(3)
            sage: a is A(1)
            True
            sage: b is B(2)
            True
            sage: c is C(3)
            True
            sage: B._clear_cache_()

        Now, all instances of (sub-classes of) ``B`` have disappeared
        from the cache::

            sage: a is A(1)
            True
            sage: b is B(2)
            False
            sage: c is C(3)
            False

        Here is a similar example, using a private classcall in the class
        ``B``, which is not inherited by ``C``::

            sage: class A(UniqueRepresentation):
            ....:     def __init__(self, x):
            ....:         pass
            sage: class B(A):
            ....:     @staticmethod
            ....:     def __classcall_private__(cls, *args, **kwds):
            ....:         print("Private B")
            ....:         return super().__classcall__(cls, *args, **kwds)
            sage: class C(B): pass
            sage: a = A(1)
            sage: b = B(2)
            Private B
            sage: c = C(3)
            sage: a is A(1)
            True
            sage: b is B(2)
            Private B
            True
            sage: c is C(3)
            True
            sage: B._clear_cache_()

        Again, all instances of (sub-classes of) ``B`` have disappeared
        from the cache::

            sage: a is A(1)
            True
            sage: b is B(2)
            Private B
            False
            sage: c is C(3)
            False
        """
        cache = None
        # not clear why the loop below is taking the last C
        for C in cls.mro():
            try:
                cache = C.__classcall__.cache
            except AttributeError:
                pass
        del_list = [k for k in cache if issubclass(k[0][0], cls)]
        for k in del_list:
            del cache[k]


class UniqueRepresentation(WithEqualityById, CachedRepresentation):
    r"""
    Classes derived from ``UniqueRepresentation`` inherit a unique
    representation behavior for their instances.

    .. SEEALSO::

        :mod:`~sage.structure.unique_representation`

    EXAMPLES:

    The short story: to construct a class whose instances have a
    unique representation behavior one just has to do::

        sage: class MyClass(UniqueRepresentation):
        ....:     # all the rest as usual
        ....:     pass

    Everything below is for the curious or for advanced usage.

    .. rubric:: What is unique representation?

    Instances of a class have a *unique representation behavior* when
    instances evaluate equal if and only if they are identical (i.e., share
    the same memory representation), if and only if they were created using
    equal arguments. For example, calling twice::

        sage: f = SymmetricFunctions(QQ)                                                # needs sage.combinat sage.modules
        sage: g = SymmetricFunctions(QQ)                                                # needs sage.combinat sage.modules

    to create the symmetric function algebra over `\QQ` actually gives back the
    same object::

        sage: f == g                                                                    # needs sage.combinat sage.modules
        True
        sage: f is g                                                                    # needs sage.combinat sage.modules
        True

    This is a standard design pattern. It allows for sharing cached data (say
    representation theoretical information about a group) as well as for very
    fast hashing and equality testing. This behaviour is typically desirable
    for parents and categories. It can also be useful for intensive
    computations where one wants to cache all the operations on a small set of
    elements (say the multiplication table of a small group), and access this
    cache as quickly as possible.

    :class:`UniqueRepresentation` is very easy to use: a class just needs to
    derive from it, or make sure some of its super classes does. Also, it
    groups together the class and the factory in a single gadget::

        sage: isinstance(SymmetricFunctions(CC), SymmetricFunctions)                    # needs sage.combinat sage.modules
        True
        sage: issubclass(SymmetricFunctions, UniqueRepresentation)                      # needs sage.combinat sage.modules
        True

    This nice behaviour is not available when one just uses a factory::

        sage: isinstance(GF(7), GF)
        Traceback (most recent call last):
        ...
        TypeError: isinstance() arg 2 must be a type...

        sage: isinstance(GF, sage.structure.factory.UniqueFactory)
        True

    In addition, :class:`~sage.structure.factory.UniqueFactory` only provides
    the *cached* representation behaviour, but not the *unique* representation
    behaviour---the examples in :mod:`~sage.structure.unique_representation`
    explain this difference.

    On the other hand, the :class:`UniqueRepresentation` class is more
    intrusive, as it imposes a behavior (and a metaclass) on all the
    subclasses. In particular, the unique representation behaviour is imposed
    on *all* subclasses (unless the ``__classcall__`` method is overloaded and
    not called in the subclass, which is not recommended). Its implementation
    is also more technical, which leads to some subtleties.

    EXAMPLES:

    We start with a simple class whose constructor takes a single value as
    argument. This pattern is similar to what is done in
    :class:`sage.combinat.sf.sf.SymmetricFunctions`::

        sage: class MyClass(UniqueRepresentation):
        ....:     def __init__(self, value):
        ....:         self.value = value

    Two coexisting instances of ``MyClass`` created with the same argument
    data are guaranteed to share the same identity. Since :issue:`12215`, this
    is only the case if there is some strong reference to the returned
    instance, since otherwise it may be garbage collected::

        sage: x = MyClass(1)
        sage: y = MyClass(1)
        sage: x is y               # There is a strong reference
        True
        sage: z = MyClass(2)
        sage: x is z
        False

    In particular, modifying any one of them modifies the other
    (reference effect)::

        sage: x.value = 3
        sage: x.value, y.value
        (3, 3)
        sage: y.value = 1
        sage: x.value, y.value
        (1, 1)

    When comparing two instances of a unique representation with ``==``
    or ``!=`` comparison by identity is used::

        sage: x == y
        True
        sage: x is y
        True
        sage: z = MyClass(2)
        sage: x == z
        False
        sage: x is z
        False
        sage: x != y
        False
        sage: x != z
        True

    A hash function equivalent to :meth:`object.__hash__` is used, which is
    compatible with comparison by identity. However this means that the hash
    function may change in between Sage sessions, or even within the same Sage
    session.
    ::

        sage: hash(x) == object.__hash__(x)
        True

    .. WARNING::

        It is possible to inherit from
        :class:`~sage.structure.unique_representation.UniqueRepresentation`
        and then overload comparison in a way that destroys the unique
        representation property. We strongly recommend against it!  You should
        use :class:`~sage.structure.unique_representation.CachedRepresentation`
        instead.

    .. rubric:: Mixing super types and super classes

    TESTS:

    For the record, this test did fail with previous implementation
    attempts::

        sage: class bla(UniqueRepresentation, SageObject):
        ....:     pass
        sage: b = bla()
    """
