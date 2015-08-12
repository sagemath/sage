"""
Dynamic classes

.. rubric:: Why dynamic classes?

The short answer:

- Multiple inheritance is a powerful tool for constructing new classes
  by combining preexisting building blocks.
- There is a combinatorial explosion in the number of potentially
  useful classes that can be produced this way.
- The implementation of standard mathematical constructions calls for
  producing such combinations automatically.
- Dynamic classes, i.e. classes created on the fly by the Python
  interpreter, are a natural mean to achieve this.

The long answer:

Say we want to construct a new class ``MyPermutation`` for
permutations in a given set `S` (in Sage, `S` will be modelled by a
parent, but we won't discuss this point here).  First, we have to
choose a data structure for the permutations, typically among the
following:

- Stored by cycle type
- Stored by code
- Stored in list notation
  - C arrays of short ints (for small permutations)
  - python lists of ints (for huge permutations)
  - ...
- Stored by reduced word
- Stored as a function
- ...

Luckily, the Sage library provides (or will provide) classes
implementing each of those data structures. Those classes all share a
common interface (or possibly a common abstract base class). So we can
just derive our class from the chosen one::

    class MyPermutation(PermutationCycleType):
        ...

Then we may want to further choose a specific memory behavior (unique
representation, copy-on-write) which (hopefuly) can again be achieved
by inheritance::

    class MyPermutation(UniqueRepresentation, PermutationCycleType):
         ...

Finaly, we may want to endow the permutations in `S` with further
operations coming from the (algebraic) structure of `S`:

- group operations
- or just monoid operations (for a subset of permutations not stable by inverse)
- poset operations (for left/right/Bruhat order)
- word operations (searching for substrings, patterns, ...)

Or any combination thereof. Now, our class typically looks like::

    class MyPermutation(UniqueRepresentation, PermutationCycleType, PosetElement, GroupElement):
         ...

Note the combinatorial explosion in the potential number of classes
which can be created this way.


In practice, such classes will be used in mathematical constructions
like::

    SymmetricGroup(5).subset(... TODO: find a good example in the context above ...)

In such a construction, the structure of the result, and therefore the
operations on its elements can only be determined at execution
time. Let us take another standard construction::

    A = cartesian_product( B, C )

Depending on the structure of `B` and `C`, and possibly on further
options passed down by the user, `A` may be:

- an enumerated set
- a group
- an algebra
- a poset
- ...

Or any combination thereof.

Hardcoding classes for all potential combinations would be at best
tedious. Furthermore, this would require a cumbersome mechanism to
lookup the appropriate class depending on the desired combination.

Instead, one may use the ability of Python to create new classes
dynamicaly::

    type("class name", tuple of base classes, dictionary of methods)

This paradigm is powerful, but there are some technicalities to
address. The purpose of this library is to standardize its use within
Sage, and in particular to ensure that the constructed classes are
reused whenever possible (unique representation), and can be pickled.

.. rubric:: Combining dynamic classes and Cython classes

Cython classes cannot inherit from a dynamic class (there might be
some partial support for this in the future). On the other hand, such
an inheritance can be partially emulated using :meth:`__getattr__`. See
``sage.categories.examples.semigroups_cython`` for an example.

"""

#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import weak_cached_function
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.inherit_comparison import InheritComparisonMetaclass, InheritComparisonClasscallMetaclass

def dynamic_class(name, bases, cls=None, reduction=None, doccls=None,
                  prepend_cls_bases=True, cache=True):
    r"""
    INPUT:

    - ``name`` -- a string
    - ``bases`` -- a tuple of classes
    - ``cls`` -- a class or ``None``
    - ``reduction`` -- a tuple or ``None``
    - ``doccls`` -- a class or ``None``
    - ``prepend_cls_bases`` -- a boolean (default: ``True``)
    - ``cache`` -- a boolean or ``"ignore_reduction"`` (default: ``True``)

    Constructs dynamically a new class ``C`` with name ``name``, and
    bases ``bases``. If ``cls`` is provided, then its methods will be
    inserted into ``C``, and its bases will be prepended to ``bases``
    (unless ``prepend_cls_bases`` is ``False``).

    The module, documentation and source instrospection is taken from
    ``doccls``, or ``cls`` if ``doccls`` is ``None``, or ``bases[0]``
    if both are ``None`` (therefore ``bases`` should be non empty if
    ``cls` is ``None``).

    The constructed class can safely be pickled (assuming the
    arguments themselves can).

    Unless ``cache`` is ``False``, the result is cached, ensuring unique
    representation of dynamic classes.

    See :mod:`sage.structure.dynamic_class` for a discussion of the
    dynamic classes paradigm, and its relevance to Sage.

    EXAMPLES:

    To setup the stage, we create a class Foo with some methods,
    cached methods, and lazy attributes, and a class Bar::

        sage: from sage.misc.lazy_attribute import lazy_attribute
        sage: from sage.misc.cachefunc import cached_function
        sage: from sage.structure.dynamic_class import dynamic_class
        sage: class Foo(object):
        ...       "The Foo class"
        ...       def __init__(self, x):
        ...           self._x = x
        ...       @cached_method
        ...       def f(self):
        ...           return self._x^2
        ...       def g(self):
        ...           return self._x^2
        ...       @lazy_attribute
        ...       def x(self):
        ...           return self._x
        ...
        sage: class Bar:
        ...       def bar(self):
        ...           return self._x^2
        ...

    We now create a class FooBar which is a copy of Foo, except that it
    also inherits from Bar::

        sage: FooBar = dynamic_class("FooBar", (Bar,), Foo)
        sage: x = FooBar(3)
        sage: x.f()
        9
        sage: x.f() is x.f()
        True
        sage: x.x
        3
        sage: x.bar()
        9
        sage: FooBar.__name__
        'FooBar'
        sage: FooBar.__module__
        '__main__'

        sage: Foo.__bases__
        (<type 'object'>,)
        sage: FooBar.__bases__
        (<type 'object'>, <class __main__.Bar at ...>)
        sage: Foo.mro()
        [<class '__main__.Foo'>, <type 'object'>]
        sage: FooBar.mro()
        [<class '__main__.FooBar'>, <type 'object'>, <class __main__.Bar at ...>]

    .. RUBRIC:: Pickling

    Dynamic classes are pickled by construction. Namely, upon
    unpickling, the class will be reconstructed by recalling
    dynamic_class with the same arguments::

        sage: type(FooBar).__reduce__(FooBar)
        (<function dynamic_class at ...>, ('FooBar', (<class __main__.Bar at ...>,), <class '__main__.Foo'>, None, None))

    Technically, this is achieved by using a metaclass, since the
    Python pickling protocol for classes is to pickle by name::

        sage: type(FooBar)
        <class 'sage.structure.dynamic_class.DynamicMetaclass'>


    The following (meaningless) example illustrates how to customize
    the result of the reduction::

        sage: BarFoo = dynamic_class("BarFoo", (Foo,), Bar, reduction = (str, (3,)))
        sage: type(BarFoo).__reduce__(BarFoo)
        (<type 'str'>, (3,))
        sage: loads(dumps(BarFoo))
        '3'

    .. RUBRIC:: Caching

    By default, the built class is cached::

         sage: dynamic_class("FooBar", (Bar,), Foo) is FooBar
         True
         sage: dynamic_class("FooBar", (Bar,), Foo, cache=True) is FooBar
         True

    and the result depends on the reduction::

        sage: dynamic_class("BarFoo", (Foo,), Bar, reduction = (str, (3,))) is BarFoo
        True
        sage: dynamic_class("BarFoo", (Foo,), Bar, reduction = (str, (2,))) is BarFoo
        False

    With ``cache=False``, a new class is created each time::

         sage: FooBar1 = dynamic_class("FooBar", (Bar,), Foo, cache=False); FooBar1
         <class '__main__.FooBar'>
         sage: FooBar2 = dynamic_class("FooBar", (Bar,), Foo, cache=False); FooBar2
         <class '__main__.FooBar'>
         sage: FooBar1 is FooBar
         False
         sage: FooBar2 is FooBar1
         False

    With ``cache="ignore_reduction"``, the class does not depend on
    the reduction::

        sage: BarFoo = dynamic_class("BarFoo", (Foo,), Bar, reduction = (str, (3,)), cache="ignore_reduction")
        sage: dynamic_class("BarFoo", (Foo,), Bar, reduction = (str, (2,)), cache="ignore_reduction") is BarFoo
        True

    In particular, the reduction used is that provided upon creating the
    first class::

        sage: dynamic_class("BarFoo", (Foo,), Bar, reduction = (str, (2,)), cache="ignore_reduction")._reduction
        (<type 'str'>, (3,))

    .. WARNING::

        The behaviour upon creating several dynamic classes from the
        same data but with different values for ``cache`` option is
        currently left unspecified. In other words, for a given
        application, it is recommended to consistently use the same
        value for that option.

    TESTS::

        sage: import __main__
        sage: __main__.Foo = Foo
        sage: __main__.Bar = Bar
        sage: x = FooBar(3)
        sage: x.__dict__      # Breaks without the __dict__ deletion in dynamic_class_internal
        {'_x': 3}

        sage: type(FooBar).__reduce__(FooBar)
        (<function dynamic_class at ...>, ('FooBar', (<class __main__.Bar at ...>,), <class '__main__.Foo'>, None, None))
        sage: import cPickle
        sage: cPickle.loads(cPickle.dumps(FooBar)) == FooBar
        True

    We check that instrospection works reasonably::

        sage: sage.misc.sageinspect.sage_getdoc(FooBar)
        'The Foo class\n'

    Finally, we check that classes derived from UniqueRepresentation
    are handled gracefuly (despite them also using a metaclass)::

        sage: FooUnique = dynamic_class("Foo", (Bar, UniqueRepresentation))
        sage: loads(dumps(FooUnique)) is FooUnique
        True

    """
    bases = tuple(bases)
    #assert(len(bases) > 0 )
    assert(isinstance(name, str))
    #    assert(cls is None or issubtype(type(cls), type) or type(cls) is classobj)
    if cache is True:
        return dynamic_class_internal(name, bases, cls, reduction, doccls, prepend_cls_bases)
    elif cache is False:
        # bypass the cached method
        return dynamic_class_internal.f(name, bases, cls, reduction, doccls, prepend_cls_bases)
    else: # cache = "ignore_reduction"
        result = dynamic_class_internal(name, bases, cls, False, doccls, prepend_cls_bases)
        if result._reduction is False:
            result._reduction = reduction
        return result


@weak_cached_function
def dynamic_class_internal(name, bases, cls=None, reduction=None, doccls=None, prepend_cls_bases=True):
    r"""
    See sage.structure.dynamic_class.dynamic_class? for indirect doctests.

    TESTS::

        sage: Foo1 = sage.structure.dynamic_class.dynamic_class_internal("Foo", (object,))
        sage: Foo2 = sage.structure.dynamic_class.dynamic_class_internal("Foo", (object,), doccls = sage.structure.dynamic_class.TestClass)
        sage: Foo3 = sage.structure.dynamic_class.dynamic_class_internal("Foo", (object,), cls    = sage.structure.dynamic_class.TestClass)
        sage: all(Foo.__name__  == 'Foo'    for Foo in [Foo1, Foo2, Foo3])
        True
        sage: all(Foo.__bases__ == (object,) for Foo in [Foo1, Foo2, Foo3])
        True
        sage: Foo1.__module__ == object.__module__
        True
        sage: Foo2.__module__ == sage.structure.dynamic_class.TestClass.__module__
        True
        sage: Foo3.__module__ == sage.structure.dynamic_class.TestClass.__module__
        True
        sage: Foo1.__doc__ == object.__doc__
        True
        sage: Foo2.__doc__ == sage.structure.dynamic_class.TestClass.__doc__
        True
        sage: Foo3.__doc__ == sage.structure.dynamic_class.TestClass.__doc__
        True

    We check that instrospection works reasonably::

        sage: from sage.misc.sageinspect import sage_getfile, sage_getsourcelines
        sage: sage_getfile(Foo2)
        '.../sage/structure/dynamic_class.py'
        sage: sage_getfile(Foo3)
        '.../sage/structure/dynamic_class.py'
        sage: sage_getsourcelines(Foo2)
        (['class TestClass:...'], ...)
        sage: sage_getsourcelines(Foo3)
        (['class TestClass:...'], ...)
        sage: sage_getsourcelines(Foo2())
        (['class TestClass:...'], ...)
        sage: sage_getsourcelines(Foo3())
        (['class TestClass:...'], ...)
        sage: sage_getsourcelines(Foo3().bla)
        (['    def bla():...'], ...)

    """
    if reduction is None:
        reduction = (dynamic_class, (name, bases, cls, reduction, doccls))
    if cls is not None:
        methods = dict(cls.__dict__)
        # Anything else that should not be kept?
        if "__dict__" in methods:
            methods.__delitem__("__dict__")
        if prepend_cls_bases:
            bases = cls.__bases__ + bases
    else:
        methods = {}
    if doccls is None:
        if cls is not None:
            doccls = cls
        else:
            assert bases != ()
            doccls = bases[0]
    methods['_reduction'] = reduction
    if "_sage_src_lines_" not in methods:
        from sage.misc.sageinspect import sage_getsourcelines
        @staticmethod
        def _sage_src_lines():
            return sage_getsourcelines(doccls)
        methods['_sage_src_lines_'] = _sage_src_lines
    methods['__doc__'] = doccls.__doc__
    methods['__module__'] = doccls.__module__

    metaclass = DynamicMetaclass
    # The metaclass of a class must derive from the metaclasses of its
    # bases. The following handles the case where one of the base
    # classes is a known Sage metaclass.  This approach won't scale
    # well if we start using metaclasses seriously in Sage.
    for base in bases:
        if isinstance(base, InheritComparisonClasscallMetaclass):
            metaclass = DynamicInheritComparisonClasscallMetaclass
        elif isinstance(base, ClasscallMetaclass):
            metaclass = DynamicClasscallMetaclass
        elif isinstance(base, InheritComparisonMetaclass):
            metaclass = DynamicInheritComparisonMetaclass
    return metaclass(name, bases, methods)

class DynamicMetaclass(type):
    """
    A metaclass implementing an appropriate reduce-by-construction method
    """
    def __reduce__(self):
        """
        See :func:`sage.structure.dynamic_class.dynamic_class` for
        non-trivial tests.

        TESTS::

            sage: class Foo: pass
            sage: class DocClass: pass
            sage: C = sage.structure.dynamic_class.dynamic_class_internal("bla", (object,), Foo, doccls = DocClass)
            sage: type(C).__reduce__(C)
            (<function dynamic_class at ...>,
             ('bla', (<type 'object'>,), <class __main__.Foo at ...>, None, <class __main__.DocClass at ...>))
            sage: C = sage.structure.dynamic_class.dynamic_class_internal("bla", (object,), Foo, doccls = DocClass, reduction = "blah")
            sage: type(C).__reduce__(C)
            'blah'
        """
        return self._reduction

class DynamicClasscallMetaclass(DynamicMetaclass, ClasscallMetaclass):
    pass

class DynamicInheritComparisonMetaclass(DynamicMetaclass, InheritComparisonMetaclass):
    pass

class DynamicInheritComparisonClasscallMetaclass(DynamicMetaclass, InheritComparisonClasscallMetaclass):
    pass

# This registers the appropriate reduction methods (see Trac #5985)
import copy_reg
for M in [DynamicMetaclass,
          DynamicClasscallMetaclass,
          DynamicInheritComparisonMetaclass,
          DynamicInheritComparisonClasscallMetaclass]:
    copy_reg.pickle(M, M.__reduce__)


class TestClass:
    """
    A class used for checking that introspection works
    """
    def bla():
        """
        bla ...
        """
        pass
