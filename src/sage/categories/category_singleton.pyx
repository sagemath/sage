# sage_setup: distribution = sagemath-objects
r"""
Singleton categories
"""
# ****************************************************************************
#  Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.constant_function import ConstantFunction
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category import Category
from sage.structure.category_object cimport CategoryObject
from sage.structure.dynamic_class import DynamicMetaclass

from cpython.type cimport PyType_IsSubtype

# This helper class is used to implement Category_singleton.__contains__
# In particular, the docstring is what appears upon C.__contains__?
# for C a singleton category like Fields().

cdef class Category_contains_method_by_parent_class:
    """
    Return whether ``x`` is an object in this category.

    More specifically, returns ``True`` if and only if ``x`` has a
    category which is a subcategory of this one.

    EXAMPLES::

        sage: ZZ in Sets()
        True
    """

    def __init__(self, category):
        """
        TESTS::

            sage: from sage.categories.category_singleton import Category_contains_method_by_parent_class
            sage: Category_contains_method_by_parent_class(Rings())
            <sage.categories.category_singleton.Category_contains_method_by_parent_class object at ...
        """
        self._parent_class_of_category = <type> category.parent_class

    def __call__(self, x):
        """
        EXAMPLES::

            sage: from sage.categories.category_singleton import Category_contains_method_by_parent_class
            sage: in_Fields = Category_contains_method_by_parent_class(Fields())
            sage: in_Fields(QQ)
            True
            sage: in_Fields(ZZ)
            False
            sage: in_Fields(1)              # Not a CategoryObject
            False
            sage: in_Fields(int(1))         # Not a SageObject
            False

        TESTS:

            The following used to segfault in a preliminary version of the
            code::

                sage: None in Rings()
                False
        """
        if x is None:
            return False
        cdef CategoryObject y
        try:
            y = x
            return PyType_IsSubtype(<type>((y._category or y.category()).parent_class), self._parent_class_of_category)
        except AttributeError:
            return False
        except TypeError: # this is for objects that aren't CategoryObjects
            try:
                return PyType_IsSubtype(<type>(x.category().parent_class), self._parent_class_of_category)
            except AttributeError:
                return False


class Category_singleton(Category):
    """
    A base class for implementing singleton category.

    A *singleton* category is a category whose class takes no
    parameters like ``Fields()`` or ``Rings()``. See also the
    :wikipedia:`Singleton design pattern <Singleton_pattern>`.

    This is a subclass of :class:`Category`, with a couple
    optimizations for singleton categories.

    The main purpose is to make the idioms::

        sage: QQ in Fields()
        True
        sage: ZZ in Fields()
        False

    as fast as possible, and in particular competitive to calling a
    constant Python method, in order to foster its systematic use
    throughout the Sage library. Such tests are time critical, in
    particular when creating a lot of polynomial rings over small
    fields like in the elliptic curve code.

    EXAMPLES::

        sage: from sage.categories.category_singleton import Category_singleton
        sage: class MyRings(Category):
        ....:     def super_categories(self): return Rings().super_categories()
        sage: class MyRingsSingleton(Category_singleton):
        ....:     def super_categories(self): return Rings().super_categories()

    We create three rings. One of them is contained in the usual category of
    rings, one in the category of "my rings" and the third in the category of
    "my rings singleton"::

        sage: R = QQ['x,y']
        sage: R1 = Parent(category = MyRings())
        sage: R2 = Parent(category = MyRingsSingleton())
        sage: R in MyRings()
        False
        sage: R1 in MyRings()
        True
        sage: R1 in MyRingsSingleton()
        False
        sage: R2 in MyRings()
        False
        sage: R2 in MyRingsSingleton()
        True

    One sees that containment tests for the singleton class is a lot faster
    than for a usual class::

        sage: # not tested
        sage: timeit("R in MyRings()", number=10000)
        10000 loops, best of 3: 7.12 µs per loop
        sage: timeit("R1 in MyRings()", number=10000)
        10000 loops, best of 3: 6.98 µs per loop
        sage: timeit("R in MyRingsSingleton()", number=10000)
        10000 loops, best of 3: 3.08 µs per loop
        sage: timeit("R2 in MyRingsSingleton()", number=10000)
        10000 loops, best of 3: 2.99 µs per loop

    So this is an improvement, but not yet competitive with a pure
    Cython method::

        sage: timeit("R.is_ring()", number=10000)                     # not tested
        10000 loops, best of 3: 383 ns per loop

    However, it is competitive with a Python method. Actually it is faster,
    if one stores the category in a variable::

        sage: _Rings = Rings()
        sage: R3 = Parent(category = _Rings)
        sage: R3.is_ring.__module__
        'sage.categories.rings'
        sage: timeit("R3.is_ring()", number=10000)                    # not tested
        10000 loops, best of 3: 2.64 µs per loop
        sage: timeit("R3 in Rings()", number=10000)                   # not tested
        10000 loops, best of 3: 3.01 µs per loop
        sage: timeit("R3 in _Rings", number=10000)                    # not tested
        10000 loops, best of 3: 652 ns per loop

    This might not be easy to further optimize, since the time is
    consumed in many different spots::

        sage: timeit("MyRingsSingleton.__classcall__()", number=10000)# not tested
        10000 loops, best of 3: 306 ns per loop

        sage: X = MyRingsSingleton()
        sage: timeit("R in X  ", number=10000)                        # not tested
        10000 loops, best of 3: 699 ns per loop

        sage: c = MyRingsSingleton().__contains__
        sage: timeit("c(R)", number = 10000)                          # not tested
        10000 loops, best of 3: 661 ns per loop

    .. WARNING::

        A singleton concrete class `A` should not have a subclass `B`
        (necessarily concrete). Otherwise, creating an instance `a` of
        `A` and an instance `b` of `B` would break the singleton
        principle: `A` would have two instances `a` and `b`.

        With the current implementation only direct subclasses of
        :class:`Category_singleton` are supported::

            sage: class MyRingsSingleton(Category_singleton):
            ....:     def super_categories(self): return Rings().super_categories()
            sage: class Disaster(MyRingsSingleton): pass
            sage: Disaster()
            Traceback (most recent call last):
            ...
            AssertionError: <class '__main__.Disaster'> is not a direct subclass of <class 'sage.categories.category_singleton.Category_singleton'>

        However, it is acceptable for a direct subclass `R` of
        :class:`Category_singleton` to create its unique instance as
        an instance of a subclass of itself (in which case, its the
        subclass of `R` which is concrete, not `R` itself). This is
        used for example to plug in extra category code via a dynamic
        subclass::

            sage: from sage.categories.category_singleton import Category_singleton
            sage: class R(Category_singleton):
            ....:     def super_categories(self): return [Sets()]
            sage: R()
            Category of r
            sage: R().__class__
            <class '__main__.R_with_category'>
            sage: R().__class__.mro()
            [<class '__main__.R_with_category'>,
             <class '__main__.R'>,
             <class 'sage.categories.category_singleton.Category_singleton'>,
             <class 'sage.categories.category.Category'>,
             <class 'sage.structure.unique_representation.UniqueRepresentation'>,
             <class 'sage.misc.fast_methods.WithEqualityById'>,
             <class 'sage.structure.unique_representation.CachedRepresentation'>,
             <class 'sage.structure.unique_representation.WithPicklingByInitArgs'>,
             <class 'sage.structure.sage_object.SageObject'>,
             <class '__main__.R.subcategory_class'>,
             <class 'sage.categories.sets_cat.Sets.subcategory_class'>,
             <class 'sage.categories.sets_with_partial_maps.SetsWithPartialMaps.subcategory_class'>,
             <class 'sage.categories.objects.Objects.subcategory_class'>,
             <class 'object'>]
            sage: R() is R()
            True
            sage: R() is R().__class__()
            True

        In that case, ``R`` is an abstract class and has a single
        concrete subclass, so this does not break the Singleton design
        pattern.

        .. SEEALSO:: :meth:`Category.__classcall__`, :meth:`Category.__init__`

    TESTS::

        sage: import __main__
        sage: __main__.MyRings = MyRings
        sage: __main__.MyRingsSingleton = MyRingsSingleton
        sage: TestSuite(MyRingsSingleton()).run(skip=["_test_category"])

    .. NOTE::

        The ``_test_category`` test is failing because
        ``MyRingsSingleton()`` is not a subcategory of the join of its
        super categories::

            sage: C = MyRingsSingleton()
            sage: C.super_categories()
            [Category of rngs, Category of semirings]
            sage: Rngs() & Semirings()
            Category of rings
            sage: C.is_subcategory(Rings())
            False

        Oh well; it's not really relevant for those tests.
    """

    # That is just an optimized constant cached_method
    @staticmethod
    def __classcall__(object cls, *args):
        """
        Return ``cls()`` and cache the result in ``cls``.

        INPUT:

        - ``*args`` -- some constant arguments

        Most of the time, ``args`` is meant to be empty. However some
        singleton categories, in particular axiom categories of
        singleton categories, may require a constant argument.
        ``*args`` is passed down to :meth:`__init__`, and ignored upon
        later calls.

        .. SEEALSO:: :class:`sage.categories.category_with_axiomCategoryWithAxiom_singleton`

        TESTS::

            sage: from sage.categories.category_singleton import Category_singleton
            sage: class MyRingsSingleton(Category_singleton):
            ....:     def super_categories(self): return Rings().super_categories()
            sage: MyRingsSingleton()
            Category of my rings singleton

        Instantiating :class:`Category_singleton` triggers an assertion error::

            sage: Category_singleton()
            Traceback (most recent call last):
            ...
            AssertionError: <class 'sage.categories.category_singleton.Category_singleton'> is not a direct subclass of <class 'sage.categories.category_singleton.Category_singleton'>

        Instantiating a subclass of a subclass of :class:`Category_singleton`
        also triggers an assertion error::

            sage: class MyStuff(Category_singleton):
            ....:     def super_categories(self): return [Sets()]
            sage: class MySubStuff(MyStuff): pass
            sage: MySubStuff()
            Traceback (most recent call last):
            ...
            AssertionError: <class '__main__.MySubStuff'> is not a direct subclass of <class 'sage.categories.category_singleton.Category_singleton'>

        even if ``MyStuff`` has already been instantiated::

            sage: MyStuff()
            Category of my stuff
            sage: MySubStuff()
            Traceback (most recent call last):
            ...
            AssertionError: <class '__main__.MySubStuff'> is not a direct subclass of <class 'sage.categories.category_singleton.Category_singleton'>
        """
        if isinstance(cls, DynamicMetaclass):  # cls is something like Rings_with_category
            cls = cls.__base__
        # TODO: find a better way to check that cls is an abstract class
        from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
        assert (cls.__mro__[1] is Category_singleton or cls.__mro__[1] is CategoryWithAxiom_singleton), \
            "{} is not a direct subclass of {}".format(cls, Category_singleton)
        obj = super().__classcall__(cls, *args)
        cls._set_classcall(ConstantFunction(obj))
        obj.__class__._set_classcall(ConstantFunction(obj))
        return obj

    @lazy_class_attribute
    def __contains__(cls):
        """
        TESTS::

            sage: from sage.categories.category_singleton import Category_singleton
            sage: class MyRingsSingleton(Category_singleton):
            ....:     def super_categories(self): return Rings().super_categories()
            sage: ZZ in MyRingsSingleton()
            False
            sage: Parent(category=MyRingsSingleton()) in MyRingsSingleton()
            True
        """
        return Category_contains_method_by_parent_class(cls())
