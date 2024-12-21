# sage_setup: distribution = sagemath-objects
r"""
Covariant Functorial Constructions

A *functorial construction* is a collection of functors
`(F_{Cat})_{Cat}` (indexed by a collection of categories) which associate
to a sequence of parents `(A, B, ...)` in a category `Cat` a parent
`F_{Cat}(A, B, ...)`. Typical examples of functorial
constructions are :obj:`cartesian_product` and :obj:`tensor_product`.

The category of `F_{Cat}(A, B, ...)`, which only depends on `Cat`, is
called the (functorial) construction category.

A functorial construction is *(category)-covariant* if for every
categories `Cat` and `SuperCat`, the category of `F_{Cat}(A, B, ...)`
is a subcategory of the category of `F_{SuperCat}(A, B, ...)` whenever
`Cat` is a subcategory of `SuperCat`. A functorial construction is
*(category)-regressive* if the category of `F_{Cat}(A, B, ...)` is a
subcategory of `Cat`.

The goal of this module is to provide generic support for covariant
functorial constructions. In particular, given some parents `A`, `B`,
..., in respective categories `Cat_A`, `Cat_B`, ..., it provides tools
for calculating the best known category for the parent
`F(A,B,...)`. For examples, knowing that Cartesian products of
semigroups (resp. monoids, groups) have a semigroup (resp. monoid,
group) structure, and given a group `B` and two monoids `A` and `C` it
can calculate that `A \times B \times C` is naturally endowed with a
monoid structure.

See :class:`CovariantFunctorialConstruction`,
:class:`CovariantConstructionCategory` and
:class:`RegressiveCovariantConstructionCategory` for more details.

AUTHORS:

 - Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.misc.cachefunc import cached_function, cached_method
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category import Category
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.dynamic_class import DynamicMetaclass


class CovariantFunctorialConstruction(UniqueRepresentation, SageObject):
    r"""
    An abstract class for construction functors `F` (eg `F` = Cartesian
    product, tensor product, `\QQ`-algebra, ...) such that:

     - Each category `Cat` (eg `Cat=` ``Groups()``) can provide a category
       `F_{Cat}` for parents constructed via this functor (e.g.
       `F_{Cat} =` ``CartesianProductsOf(Groups())``).

     - For every category `Cat`, `F_{Cat}` is a subcategory of
       `F_{SuperCat}` for every super category `SuperCat` of
       `Cat` (the functorial construction is (category)-covariant).

     - For parents `A`, `B`, ..., respectively in the categories
       `Cat_A`, `Cat_B`, ..., the category of `F(A,B,...)` is
       `F_{Cat}` where `Cat` is the meet of the categories `Cat_A`,
       `Cat_B`, ...,.

    This covers two slightly different use cases:

     - In the first use case, one uses directly the construction
       functor to create new parents::

         sage: tensor()  # todo: not implemented (add an example)

       or even new elements, which indirectly constructs the
       corresponding parent::

         sage: tensor(...) # todo: not implemented

     - In the second use case, one implements a parent, and then put
       it in the category `F_{Cat}` to specify supplementary
       mathematical information about that parent.

    The main purpose of this class is to handle automatically the
    trivial part of the category hierarchy. For example,
    ``CartesianProductsOf(Groups())`` is set automatically as a
    subcategory of ``CartesianProductsOf(Monoids())``.

    In practice, each subclass of this class should provide the
    following attributes:

    - ``_functor_category`` -- string which should match the name of
      the nested category class to be used in each category to
      specify information and generic operations for elements of this
      category

    - ``_functor_name`` -- string which specifies the name of the
      functor, and also (when relevant) of the method on parents and
      elements used for calling the construction

    TODO: What syntax do we want for `F_{Cat}`? For example, for the
    tensor product construction, which one do we want among (see
    chat on IRC, on 07/12/2009):

     - ``tensor(Cat)``
     - ``tensor((Cat, Cat))``
     - ``tensor.of((Cat, Cat))``
     - ``tensor.category_from_categories((Cat, Cat, Cat))``
     - ``Cat.TensorProducts()``

    The syntax ``Cat.TensorProducts()`` does not supports well multivariate
    constructions like ``tensor.of([Algebras(), HopfAlgebras(), ...])``.
    Also it forces every category to be (somehow) aware of all the
    tensorial construction that could apply to it, even those which
    are only induced from super categories.

    Note: for each functorial construction, there probably is one (or several)
    largest categories on which it applies. For example, the
    :func:`~sage.categories.cartesian_product.CartesianProducts` construction makes
    only sense for concrete categories, that is subcategories of
    ``Sets()``. Maybe we want to model this one way or the other.
    """

    def category_from_parents(self, parents):
        """
        Return the category of `F(A,B,...)` for `A,B,...` parents.

        INPUT:

         - ``self`` -- a functor `F`
         - ``parents`` -- a list (or iterable) of parents

        EXAMPLES::

            sage: E = CombinatorialFreeModule(QQ, ["a", "b", "c"])                      # needs sage.modules
            sage: tensor.category_from_parents((E, E, E))                               # needs sage.modules
            Category of tensor products of
             finite dimensional vector spaces with basis over Rational Field
        """
        from sage.structure.parent import Parent
        assert all(isinstance(parent, Parent) for parent in parents)
        # Should we pass a set of categories to reduce the cache size?
        # But then this would impose that, for any constructor, the
        # category of the result does not depend on the order/repetition
        # of the categories of the parents
        return self.category_from_categories(tuple({parent.category() for parent in parents}))

    @cached_method
    def category_from_categories(self, categories):
        """
        Return the category of `F(A,B,...)` for `A,B,...` parents in
        the given categories.

        INPUT:

         - ``self`` -- a functor `F`
         - ``categories`` -- a non empty tuple of categories

        EXAMPLES::

            sage: Cat1 = Rings()
            sage: Cat2 = Groups()
            sage: cartesian_product.category_from_categories((Cat1, Cat1, Cat1))
            Join of Category of rings and ...
                and Category of Cartesian products of monoids
                and Category of Cartesian products of commutative additive groups

            sage: cartesian_product.category_from_categories((Cat1, Cat2))
            Category of Cartesian products of monoids
        """
        assert len(categories) > 0
        return self.category_from_category(Category.meet(categories))

    def category_from_category(self, category):
        """
        Return the category of `F(A,B,...)` for `A,B,...` parents in
        ``category``.

        INPUT:

         - ``self`` -- a functor `F`
         - ``category`` -- a category

        EXAMPLES::

            sage: tensor.category_from_category(ModulesWithBasis(QQ))
            Category of tensor products of vector spaces with basis over Rational Field

        # TODO: add support for parametrized functors
        """
        return getattr(category, self._functor_category)()

    def _repr_(self):
        """
        EXAMPLES::

            sage: tensor                       # indirect doctest
            The tensor functorial construction
        """
        return "The %s functorial construction" % self._functor_name

    def __call__(self, args, **kwargs):
        """
        Functorial construction application.

        INPUT:

         - ``self`` -- a covariant functorial construction `F`
         - ``args`` -- a tuple (or iterable) of parents or elements

        Returns `F(args)`

        EXAMPLES::

            sage: E = CombinatorialFreeModule(QQ, ["a", "b", "c"]); E.rename('E')       # needs sage.modules
            sage: tensor((E, E, E))                                                     # needs sage.modules
            E # E # E
        """
        args = tuple(args)  # a bit brute force; let's see if this becomes a bottleneck later
        assert all(hasattr(arg, self._functor_name) for arg in args)
        assert len(args) > 0
        return getattr(args[0], self._functor_name)(*args[1:], **kwargs)


class FunctorialConstructionCategory(Category): # Should this be CategoryWithBase?
    """
    Abstract class for categories `F_{Cat}` obtained through a
    functorial construction
    """

    @lazy_class_attribute
    def _base_category_class(cls):
        """
        Recover the class of the base category.

        OUTPUT: a *tuple* whose single entry is the base category class

        .. WARNING::

            This is only used for functorial construction categories
            that are not implemented as nested classes, and won't work
            otherwise.

        .. SEEALSO:: :meth:`__classcall__`

        EXAMPLES::

            sage: GradedModules._base_category_class
            (<class 'sage.categories.modules.Modules'>,)
            sage: GradedAlgebrasWithBasis._base_category_class
            (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>,)

        The reason for wrapping the base category class in a tuple is
        that, often, the base category class implements a
        :meth:`__classget__` method which would get in the way upon
        attribute access::

            sage: F = GradedAlgebrasWithBasis
            sage: F._foo = F._base_category_class[0]
            sage: F._foo
            Traceback (most recent call last):
            ...
            AssertionError: base category class for <...AlgebrasWithBasis'>
             mismatch; expected <...Algebras'>,
             got <...GradedAlgebrasWithBasis'>

        We note that because ``Algebras.WithBasis`` is not lazily imported
        on startup (see :issue:`22955`), the test fails at a different
        point in the code. However, if this import becomes lazy again, then
        the following error will be generated and can replace the above::

            sage: F._foo  # not tested
            Traceback (most recent call last):
            ...
            ValueError: could not infer axiom for the nested class
            <...AlgebrasWithBasis'> of <...GradedAlgebrasWithBasis'>

        .. TODO::

            The logic is very similar to that implemented in
            :class:`CategoryWithAxiom._base_category_class`. Find a
            way to refactor this to avoid the duplication.
        """
        module_name = cls.__module__.replace(cls._functor_category.lower() + "_","")
        import sys
        name = cls.__name__.replace(cls._functor_category, "")
        __import__(module_name)
        module = sys.modules[module_name]
        return (module.__dict__[name],)

    @staticmethod
    def __classcall__(cls, category=None, *args):
        """
        Make ``XXXCat(**)`` a shorthand for ``Cat(**).XXX()``.

        EXAMPLES::

            sage: GradedModules(ZZ)   # indirect doctest
            Category of graded modules over Integer Ring
            sage: Modules(ZZ).Graded()
            Category of graded modules over Integer Ring
            sage: Modules.Graded(ZZ)
            Category of graded modules over Integer Ring
            sage: GradedModules(ZZ) is Modules(ZZ).Graded()
            True

        .. SEEALSO:: :meth:`_base_category_class`

        .. TODO::

            The logic is very similar to that implemented in
            :class:`CategoryWithAxiom.__classcall__`. Find a way to
            refactor this to avoid the duplication.
        """
        base_category_class = cls._base_category_class[0]
        if isinstance(category, base_category_class):
            return super().__classcall__(cls, category, *args)
        else:
            return cls.category_of(base_category_class(category, *args))

    @staticmethod
    def __classget__(cls, base_category, base_category_class):
        r"""
        Special binding for covariant constructions.

        This implements a hack allowing e.g. ``category.Subquotients``
        to recover the default ``Subquotients`` method defined in
        ``Category``, even if it has been overridden by a
        ``Subquotients`` class.

        EXAMPLES::

            sage: Sets.Subquotients
            <class 'sage.categories.sets_cat.Sets.Subquotients'>
            sage: Sets().Subquotients
            Cached version of <function ...Subquotients at ...>

        This method also initializes the attribute
        ``_base_category_class`` if not already set::

            sage: Sets.Subquotients._base_category_class
            (<class 'sage.categories.sets_cat.Sets'>,)

        It also forces the resolution of lazy imports (see :issue:`15648`)::

            sage: type(Algebras.__dict__["Graded"])
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: Algebras.Graded
            <class 'sage.categories.graded_algebras.GradedAlgebras'>
            sage: type(Algebras.__dict__["Graded"])
            <class 'sage.misc.classcall_metaclass.ClasscallMetaclass'>

        .. TODO::

            The logic is very similar to that implemented in
            :class:`CategoryWithAxiom.__classget__`. Find a way to
            refactor this to avoid the duplication.
        """
        if base_category is not None:
            assert base_category.__class__ is base_category_class
            assert isinstance(base_category_class, DynamicMetaclass)
        if isinstance(base_category_class, DynamicMetaclass):
            base_category_class = base_category_class.__base__
        if "_base_category_class" not in cls.__dict__:
            cls._base_category_class = (base_category_class,)
        else:
            assert cls._base_category_class[0] is base_category_class, \
                "base category class for {} mismatch; expected {}, got {}".format(
                 cls, cls._base_category_class[0], base_category_class)

        # Workaround #15648: if Sets.Subquotients is a LazyImport object,
        # this forces the substitution of the object back into Sets
        # to avoid resolving the lazy import over and over
        if isinstance(base_category_class.__dict__[cls._functor_category], LazyImport):
            setattr(base_category_class, cls._functor_category, cls)
        if base_category is None:
            return cls
        return getattr(super(base_category.__class__.__base__, base_category),
                       cls._functor_category)

    @classmethod
    @cached_function
    def category_of(cls, category, *args):
        """
        Return the image category of the functor `F_{Cat}`.

        This is the main entry point for constructing the category
        `F_{Cat}` of parents `F(A,B,...)` constructed from parents
        `A,B,...` in `Cat`.

        INPUT:

        - ``cls`` -- the category class for the functorial construction `F`
        - ``category`` -- a category `Cat`
        - ``*args`` -- further arguments for the functor

        EXAMPLES::

            sage: C = sage.categories.tensor.TensorProductsCategory
            sage: C.category_of(ModulesWithBasis(QQ))
            Category of tensor products of vector spaces with basis over Rational Field

            sage: C = sage.categories.algebra_functor.AlgebrasCategory
            sage: C.category_of(FiniteMonoids(), QQ)
            Join of Category of finite dimensional algebras with basis over Rational Field
                and Category of monoid algebras over Rational Field
                and Category of finite set algebras over Rational Field
        """
        # TODO: find a better test
        # the purpose is to test whether ``category`` implements that functor
        functor_category = getattr(category.__class__, cls._functor_category)
        if isinstance(functor_category, type) and issubclass(functor_category, Category):
            return functor_category(category, *args)
        else:
            return cls.default_super_categories(category, *args)

    def __init__(self, category, *args):
        r"""
        TESTS::

            sage: from sage.categories.covariant_functorial_construction import CovariantConstructionCategory
            sage: class FooBars(CovariantConstructionCategory):
            ....:     _functor_category = "FooBars"
            ....:     _base_category_class = (Category,)
            sage: Category.FooBars = lambda self: FooBars.category_of(self)
            sage: C = FooBars(ModulesWithBasis(ZZ))
            sage: C
            Category of foo bars of modules with basis over Integer Ring
            sage: C.base_category()
            Category of modules with basis over Integer Ring
            sage: latex(C)
            \mathbf{FooBars}(\mathbf{ModulesWithBasis}_{\Bold{Z}})
            sage: import __main__; __main__.FooBars = FooBars # Fake FooBars being defined in a python module
            sage: TestSuite(C).run()
        """
        assert isinstance(category, Category)
        self._base_category = category
        self._args = args
        super().__init__(*args)

    def base_category(self):
        """
        Return the base category of the category ``self``.

        For any category ``B`` = `F_{Cat}` obtained through a functorial
        construction `F`, the call ``B.base_category()`` returns the
        category `Cat`.

        EXAMPLES::

            sage: Semigroups().Quotients().base_category()
            Category of semigroups
        """
        return self._base_category

    def _make_named_class_key(self, name):
        return getattr(self._base_category, name)

    def extra_super_categories(self):
        """
        Return the extra super categories of a construction category.

        Default implementation which returns ``[]``.

        EXAMPLES::

            sage: Sets().Subquotients().extra_super_categories()
            []
            sage: Semigroups().Quotients().extra_super_categories()
            []
        """
        return []

    def super_categories(self):
        """
        Return the super categories of a construction category.

        EXAMPLES::

            sage: Sets().Subquotients().super_categories()
            [Category of sets]
            sage: Semigroups().Quotients().super_categories()
            [Category of subquotients of semigroups, Category of quotients of sets]
        """
        return Category.join([self.__class__.default_super_categories(self.base_category(), *self._args)] +
                             self.extra_super_categories(),
                             as_list=True)

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Semigroups().Subquotients()  # indirect doctest
            Category of subquotients of semigroups
        """
        return "%s of %s" % (Category._repr_object_names(self), self.base_category()._repr_object_names())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(Semigroups().Subquotients())   # indirect doctest
            \mathbf{Subquotients}(\mathbf{Semigroups})
            sage: latex(ModulesWithBasis(QQ).TensorProducts())
            \mathbf{TensorProducts}(\mathbf{WithBasis}_{\Bold{Q}})
            sage: latex(Semigroups().Algebras(QQ))
            \mathbf{Algebras}(\mathbf{Semigroups})
        """
        from sage.misc.latex import latex
        return "\\mathbf{%s}(%s)" % (self._short_name(), latex(self.base_category()))


class CovariantConstructionCategory(FunctorialConstructionCategory):
    """
    Abstract class for categories `F_{Cat}` obtained through a
    covariant functorial construction
    """

    @classmethod
    def default_super_categories(cls, category, *args):
        r"""
        Return the default super categories of `F_{Cat}(A,B,...)` for
        `A,B,...` parents in `Cat`.

        INPUT:

        - ``cls`` -- the category class for the functor `F`
        - ``category`` -- a category `Cat`
        - ``*args`` -- further arguments for the functor

        OUTPUT: a (join) category

        The default implementation is to return the join of the
        categories of `F(A,B,...)` for `A,B,...` in turn in each of
        the super categories of ``category``.

        This is implemented as a class method, in order to be able to
        reconstruct the functorial category associated to each of the
        super categories of ``category``.

        EXAMPLES:

        Bialgebras are both algebras and coalgebras::

            sage: Bialgebras(QQ).super_categories()
            [Category of algebras over Rational Field,
             Category of coalgebras over Rational Field]

        Hence tensor products of bialgebras are tensor products of
        algebras and tensor products of coalgebras::

            sage: Bialgebras(QQ).TensorProducts().super_categories()
            [Category of tensor products of algebras over Rational Field,
             Category of tensor products of coalgebras over Rational Field]

        Here is how :meth:`default_super_categories` was called internally::

            sage: C = sage.categories.tensor.TensorProductsCategory
            sage: C.default_super_categories(Bialgebras(QQ))
            Join of Category of tensor products of algebras over Rational Field
                and Category of tensor products of coalgebras over Rational Field

        We now show a similar example, with the ``Algebra`` functor
        which takes a parameter `\QQ`::

            sage: FiniteMonoids().super_categories()
            [Category of monoids, Category of finite semigroups]
            sage: sorted(FiniteMonoids().Algebras(QQ).super_categories(), key=str)
            [Category of finite dimensional algebras with basis over Rational Field,
             Category of finite set algebras over Rational Field,
             Category of monoid algebras over Rational Field]

        Note that neither the category of *finite* semigroup algebras
        nor that of monoid algebras appear in the result; this is
        because there is currently nothing specific implemented about them.

        Here is how :meth:`default_super_categories` was called internally::

            sage: C = sage.categories.algebra_functor.AlgebrasCategory
            sage: C.default_super_categories(FiniteMonoids(), QQ)
            Join of Category of finite dimensional algebras with basis over Rational Field
                and Category of monoid algebras over Rational Field
                and Category of finite set algebras over Rational Field
        """
        return Category.join([getattr(cat, cls._functor_category)(*args)
                              for cat in category._super_categories
                              if hasattr(cat, cls._functor_category)])

    def is_construction_defined_by_base(self):
        r"""
        Return whether the construction is defined by the base of ``self``.

        EXAMPLES:

        The graded functorial construction is defined by the modules
        category. Hence this method returns ``True`` for graded
        modules and ``False`` for other graded xxx categories::

            sage: Modules(ZZ).Graded().is_construction_defined_by_base()
            True
            sage: Algebras(QQ).Graded().is_construction_defined_by_base()
            False
            sage: Modules(ZZ).WithBasis().Graded().is_construction_defined_by_base()
            False

        This is implemented as follows: given the base category `A`
        and the construction `F` of ``self``, that is ``self=A.F()``,
        check whether no super category of `A` has `F` defined.

        .. NOTE::

            Recall that, when `A` does not implement the construction
            ``F``, a join category is returned. Therefore, in such
            cases, this method is not available::

                sage: Bialgebras(QQ).Graded().is_construction_defined_by_base()
                Traceback (most recent call last):
                ...
                AttributeError: 'JoinCategory_with_category' object has
                no attribute 'is_construction_defined_by_base'
        """
        base = self.base_category()
        f = self._functor_category
        return not any(hasattr(C, f) for C in base.super_categories())

    def additional_structure(self):
        r"""
        Return the additional structure defined by ``self``.

        By default, a functorial construction category ``A.F()``
        defines additional structure if and only if `A` is the
        category defining `F`. The rationale is that, for a
        subcategory `B` of `A`, the fact that `B.F()` morphisms shall
        preserve the `F`-specific structure is already imposed by
        `A.F()`.

        .. SEEALSO::

            - :meth:`Category.additional_structure`.
            - :meth:`is_construction_defined_by_base`.

        EXAMPLES::

            sage: Modules(ZZ).Graded().additional_structure()
            Category of graded modules over Integer Ring
            sage: Algebras(ZZ).Graded().additional_structure()

        TESTS::

            sage: Modules(ZZ).Graded().additional_structure.__module__
            'sage.categories.covariant_functorial_construction'
        """
        if self.is_construction_defined_by_base():
            return self
        else:
            return None


class RegressiveCovariantConstructionCategory(CovariantConstructionCategory):
    """
    Abstract class for categories `F_{Cat}` obtained through a
    regressive covariant functorial construction
    """

    @classmethod
    def default_super_categories(cls, category, *args):
        """
        Return the default super categories of `F_{Cat}(A,B,...)` for
        `A,B,...` parents in `Cat`.

        INPUT:

        - ``cls`` -- the category class for the functor `F`
        - ``category`` -- a category `Cat`
        - ``*args`` -- further arguments for the functor

        OUTPUT: a join category

        This implements the property that an induced subcategory is a
        subcategory.

        EXAMPLES:

        A subquotient of a monoid is a monoid, and a subquotient of
        semigroup::

            sage: Monoids().Subquotients().super_categories()
            [Category of monoids, Category of subquotients of semigroups]

        TESTS::

            sage: C = Monoids().Subquotients()
            sage: C.__class__.default_super_categories(C.base_category(), *C._args)
            Category of unital subquotients of semigroups
        """
        return Category.join([category,
                              super().default_super_categories(category, *args)])
