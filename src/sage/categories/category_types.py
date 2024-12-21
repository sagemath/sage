# sage_setup: distribution = sagemath-objects
"""
Specific category classes

This is placed in a separate file from categories.py to avoid circular imports
(as morphisms must be very low in the hierarchy with the new coercion model).
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category import Category, CategoryWithParameters, JoinCategory
from sage.misc.lazy_import import lazy_import
from sage.misc.unknown import Unknown

lazy_import('sage.categories.objects', 'Objects')
lazy_import('sage.misc.latex', 'latex')

lazy_import('sage.categories.chain_complexes', 'ChainComplexes',
            deprecation=29917)

####################################################################
#   Different types of categories
####################################################################

#############################################################
# Category of elements of some object
#############################################################


class Elements(Category):
    """
    The category of all elements of a given parent.

    EXAMPLES::

        sage: a = IntegerRing()(5)
        sage: C = a.category(); C
        Category of elements of Integer Ring
        sage: a in C
        True
        sage: 2/3 in C
        False
        sage: loads(C.dumps()) == C
        True
    """
    def __init__(self, object):
        """
        EXAMPLES::

            sage: TestSuite(Elements(ZZ)).run()
        """
        Category.__init__(self)
        self.__object = object

    @classmethod
    def an_instance(cls):
        """
        Return an instance of this class.

        EXAMPLES::

            sage: Elements.an_instance()
            Category of elements of Rational Field
        """
        from sage.rings.rational_field import QQ
        return cls(QQ)

    def _call_(self, x):
        """
        EXAMPLES::

            sage: # needs sage.modules
            sage: V = VectorSpace(QQ, 3)
            sage: x = V.0
            sage: C = x.category()
            sage: C
            Category of elements of Vector space of dimension 3 over Rational Field
            sage: w = C([1, 2, 3]); w  # indirect doctest
            (1, 2, 3)
            sage: w.category()
            Category of elements of Vector space of dimension 3 over Rational Field
        """
        return self.__object(x)

    def super_categories(self):
        """
        EXAMPLES::

            sage: Elements(ZZ).super_categories()
            [Category of objects]

        .. TODO::

            Check that this is what we want.
        """
        return [Objects()]

    def object(self):
        """
        EXAMPLES::

            sage: Elements(ZZ).object()
            Integer Ring
        """
        return self.__object

    def __reduce__(self):
        """
        EXAMPLES::

            sage: C = Elements(ZZ)
            sage: loads(dumps(C)) == C
            True
        """
        return Elements, (self.__object, )

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Elements(ZZ)._repr_object_names()
            'elements of Integer Ring'
        """
        return "elements of %s" % self.object()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: V = VectorSpace(QQ, 3)                                                # needs sage.modules
            sage: x = V.0                                                               # needs sage.modules
            sage: latex(x.category())  # indirect doctest                               # needs sage.modules
            \mathbf{Elt}_{\Bold{Q}^{3}}
        """
        return "\\mathbf{Elt}_{%s}" % latex(self.__object)


#############################################################
# Category of objects over some base object
#############################################################
class Category_over_base(CategoryWithParameters):
    r"""
    A base class for categories over some base object.

    INPUT:

    - ``base`` -- a category `C` or an object of such a category

    Assumption: the classes for the parents, elements, morphisms, of
    ``self`` should only depend on `C`. See :issue:`11935` for details.

    EXAMPLES::

        sage: Algebras(GF(2)).element_class is Algebras(GF(3)).element_class
        True

        sage: C = GF(2).category()
        sage: Algebras(GF(2)).parent_class is Algebras(C).parent_class
        True

        sage: C = ZZ.category()
        sage: Algebras(ZZ).element_class is Algebras(C).element_class
        True
    """

    def __init__(self, base, name=None):
        r"""
        Initialize ``self``.

        The ``name`` parameter is ignored.

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: C = Schemes(S); C
            Category of schemes over Integer Ring
            sage: C.__class__.__init__ == sage.categories.category_types.Category_over_base.__init__
            True
            sage: C.base() is S
            True
            sage: TestSuite(C).run()
        """
        self.__base = base
        Category.__init__(self)

    def _test_category_over_bases(self, **options):
        """
        Run generic tests on this category with parameters.

        .. SEEALSO:: :class:`TestSuite`.

        EXAMPLES::

            sage: Modules(QQ)._test_category_over_bases()
        """
        tester = self._tester(**options)
        from sage.categories.category_singleton import Category_singleton

        from .bimodules import Bimodules
        from .schemes import Schemes
        for cat in self.super_categories():
            tester.assertTrue(isinstance(cat, (Category_singleton, Category_over_base,
                                            Bimodules, Schemes)),
                           "The super categories of a category over base should"
                           " be a category over base (or the related Bimodules)"
                           " or a singleton category")

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        Since :issue:`11935`, the element and parent classes of a
        category over base only depend on the category of the base (or
        the base itself if it is a category).

        .. SEEALSO::

            - :meth:`CategoryWithParameters`
            - :meth:`CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Modules(ZZ)._make_named_class_key('element_class')
            Join of Category of Dedekind domains
             and Category of euclidean domains
             and Category of noetherian rings
             and Category of infinite enumerated sets
             and Category of metric spaces
            sage: Modules(QQ)._make_named_class_key('parent_class')
            Join of Category of number fields
             and Category of quotient fields
             and Category of metric spaces
            sage: Schemes(Spec(ZZ))._make_named_class_key('parent_class')
            Category of schemes
            sage: ModularAbelianVarieties(QQ)._make_named_class_key('parent_class')
            Join of Category of number fields
             and Category of quotient fields
             and Category of metric spaces
            sage: Algebras(Fields())._make_named_class_key('morphism_class')
            Category of fields
        """
        return getattr(self.__base if isinstance(self.__base, Category) else self.__base.category(), name)

    @classmethod
    def an_instance(cls):
        """
        Return an instance of this class.

        EXAMPLES::

            sage: Algebras.an_instance()
            Category of algebras over Rational Field
        """
        from sage.rings.rational_field import QQ
        return cls(QQ)

    def base(self):
        """
        Return the base over which elements of this category are
        defined.

        EXAMPLES::

            sage: C = Algebras(QQ)
            sage: C.base()
            Rational Field
        """
        return self.__base

    def _repr_object_names(self):
        r"""
        Return the name of the objects of this category.

        .. SEEALSO:: :meth:`Category._repr_object_names`

        EXAMPLES::

            sage: Algebras(QQ)._repr_object_names()
            'algebras over Rational Field'
            sage: Algebras(Fields())._repr_object_names()
            'algebras over fields'
            sage: Algebras(GF(2).category())._repr_object_names()
            'algebras over (finite enumerated fields and subquotients of monoids and quotients of semigroups)'
        """
        base = self.__base
        if isinstance(base, Category):
            if isinstance(base, JoinCategory):
                name = '('+' and '.join(C._repr_object_names() for C in base.super_categories())+')'
            else:
                name = base._repr_object_names()
        else:
            name = base
        return Category._repr_object_names(self) + " over %s" % name

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(ModulesWithBasis(ZZ))
            \mathbf{ModulesWithBasis}_{\Bold{Z}}
        """
        return "\\mathbf{%s}_{%s}" % (self._label, latex(self.__base))

#    def construction(self):
#        return (self.__class__, self.__base)

# How to deal with HomsetWithBase
#     def _homset(self, X, Y):
#         """
#         Given two objects X and Y in this category, returns the
#         collection of the morphisms of this category between X and Y
#         """
#         assert(X in self and Y in self)
#         from sage.categories.homset import Homset, HomsetWithBase
#         if X._base is not X and X._base is not None: # does this ever fail?
#             return HomsetWithBase(X, Y, self)
#         else:
#             return Homset(X, Y, self)

#############################################################
# Category of objects over some base ring
#############################################################


class AbelianCategory(Category):
    def is_abelian(self):
        """
        Return ``True`` as ``self`` is an abelian category.

        EXAMPLES::

            sage: CommutativeAdditiveGroups().is_abelian()
            True
        """
        return True


class Category_over_base_ring(Category_over_base):
    def __init__(self, base, name=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = Algebras(GF(2)); C
            Category of algebras over Finite Field of size 2
            sage: TestSuite(C).run()
        """
        from sage.categories.rings import Rings
        if not (base in Rings() or
                isinstance(base, Category) and base.is_subcategory(Rings())):
            raise ValueError("base must be a ring or a subcategory of Rings()")
        Category_over_base.__init__(self, base, name)

    def base_ring(self):
        """
        Return the base ring over which elements of this category are
        defined.

        EXAMPLES::

            sage: C = Algebras(GF(2))
            sage: C.base_ring()
            Finite Field of size 2
        """
        return self.base()

    def _subcategory_hook_(self, C):
        """
        A quick test whether a category ``C`` may be a subcategory of
        this category.

        INPUT:

        - ``C`` -- a category (type not tested)

        OUTPUT:

        A boolean if it is certain that ``C`` is (or is not) a
        subcategory of ``self``. :obj:`~sage.misc.unknown.Unknown`
        otherwise.

        EXAMPLES:

        The answer is ``False`` if the subcategory class of ``C`` is
        not a subclass of the subcategory class of ``self``::

            sage: Algebras(QQ)._subcategory_hook_(VectorSpaces(QQ))
            False
            sage: VectorSpaces(QQ)._subcategory_hook_(Algebras(ZZ))
            False

        .. WARNING::

            This test currently includes some false negatives::

                sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(Fields().Finite()))
                False
                sage: Modules(Rings())._subcategory_hook_(Modules(GroupAlgebras(Rings())))
                False

        The answer is ``Unknown`` if ``C`` is not a category over base ring::

            sage: VectorSpaces(QQ)._subcategory_hook_(VectorSpaces(QQ) & Rings())
            Unknown

            sage: # needs sage.combinat sage.modules
            sage: Sym = SymmetricFunctions(QQ)
            sage: from sage.combinat.sf.sfa import SymmetricFunctionsBases
            sage: Modules(QQ)._subcategory_hook_(SymmetricFunctionsBases(Sym))
            Unknown
            sage: SymmetricFunctionsBases(Sym).is_subcategory(Modules(QQ))
            True

        Case 1: the two bases are categories; then the base of ``C``
        shall be a subcategory of the base of ``self``::

            sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(Fields()))
            True
            sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(Fields().Finite()))  # todo: not implemented
            True
            sage: VectorSpaces(Fields().Finite())._subcategory_hook_(Algebras(Fields()))
            False

        Case 2: the base of ``self`` is a category; then the base of
        ``C`` shall be a parent in this category::

            sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(QQ))                 # todo: not implemented
            True
            sage: VectorSpaces(Fields().Finite())._subcategory_hook_(Algebras(QQ))
            False

        Case 3: the two bases are parents; then they should coincide::

            sage: VectorSpaces(QQ)._subcategory_hook_(Algebras(QQ))
            True
            sage: VectorSpaces(CC)._subcategory_hook_(Algebras(QQ))         # base ring in different categories         # needs sage.rings.real_mpfr
            False
            sage: VectorSpaces(GF(2))._subcategory_hook_(Algebras(GF(3)))   # base ring in the same category
            False

        Note; we need both previous tests since the distinction is
        made respectively using the parent class or the base ring::

            sage: issubclass(Algebras(QQ).parent_class,                                 # needs sage.modules
            ....:            VectorSpaces(CC).parent_class)
            False
            sage: issubclass(Algebras(GF(2)).parent_class,
            ....:            VectorSpaces(GF(3)).parent_class)
            True

        Check that :issue:`16618` is fixed: this ``_subcategory_hook_``
        method is only valid for :class:`Category_over_base_ring`, not
        :class:`Category_over_base`::

            sage: # needs sage.groups
            sage: from sage.categories.category_types import Category_over_base
            sage: D = Modules(Rings())
            sage: class Cs(Category_over_base):
            ....:    def super_categories(self):
            ....:        return [D]
            sage: C = Cs(SymmetricGroup(3))
            sage: C.is_subcategory(D)
            True
            sage: D._subcategory_hook_(C)
            Unknown
            sage: import __main__
            sage: __main__.Cs = Cs # Fake Cs being defined in a python module
            sage: TestSuite(C).run()
        """
        if not issubclass(C.parent_class, self.parent_class):
            return False
        if not isinstance(C, Category_over_base_ring):
            return Unknown
        base_ring = self.base_ring()
        if C.base_ring() is base_ring:
            return True
        if isinstance(base_ring, Category):
            if isinstance(C.base(), Category):
                return C.base().is_subcategory(base_ring)
            # otherwise, C.base() is a parent
            return C.base() in base_ring
        return False

    def __contains__(self, x):
        """
        Return whether ``x`` is an object of this category.

        In most cases, ``x`` is an object in this category, if and
        only if the category of ``x`` is a subcategory of ``self``.
        Exception: ``x`` is also an object in this category if ``x``
        is in a category over a base ring category ``C``, and ``self``
        is a category over a base ring in ``C``.

        This method implements this exception.

        EXAMPLES::

            sage: QQ['x'] in Algebras(QQ)
            True
            sage: ZZ['x'] in Algebras(ZZ)
            True

        We also would want the following to hold::

            sage: QQ['x'] in Algebras(Fields()) # todo: not implemented
            True
        """
        try:
            # The issubclass test handles extension types or when the
            # category is not fully initialized
            if isinstance(x, self.parent_class) or \
               issubclass(x.category().parent_class, self.parent_class):
                if isinstance(self.base(), Category):
                    return True
                else:
                    return x.base_ring() is self.base_ring()
            else:
                return super().__contains__(x)
        except AttributeError:
            return False


#############################################################
# Category of objects in some ambient object
#############################################################
class Category_in_ambient(Category):
    def __init__(self, ambient, name=None):
        """
        Initialize ``self``.

        The parameter ``name`` is ignored.

        EXAMPLES::

            sage: C = Ideals(IntegerRing())
            sage: TestSuite(C).run()
        """
        self.__ambient = ambient
        Category.__init__(self)

    def ambient(self):
        """
        Return the ambient object in which objects of this category are
        embedded.

        EXAMPLES::

            sage: C = Ideals(IntegerRing())
            sage: C.ambient()
            Integer Ring
        """
        return self.__ambient

    def _repr_(self):
        """
        EXAMPLES::

            sage: Ideals(IntegerRing())
            Category of ring ideals in Integer Ring
        """
        return Category._repr_(self) + " in %s" % self.__ambient

#    def construction(self):
#        return (self.__class__, self.__ambient)


class Category_module(AbelianCategory, Category_over_base_ring):
    pass


class Category_ideal(Category_in_ambient):

    @classmethod
    def an_instance(cls):
        """
        Return an instance of this class.

        EXAMPLES::

            sage: AlgebraIdeals.an_instance()
            Category of algebra ideals in Univariate Polynomial Ring in x over Rational Field
        """
        from sage.rings.rational_field import QQ
        return cls(QQ['x'])

    def ring(self):
        """
        Return the ambient ring used to describe objects ``self``.

        EXAMPLES::

            sage: C = Ideals(IntegerRing())
            sage: C.ring()
            Integer Ring
        """
        return self.ambient()

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: C = Ideals(IntegerRing())
            sage: IntegerRing().zero_ideal() in C
            True
        """
        if super().__contains__(x):
            return True
        from sage.rings.ideal import Ideal_generic
        return isinstance(x, Ideal_generic) and x.ring() == self.ring()

    def __call__(self, v):
        """
        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Ig = [x, y]
            sage: I = R.ideal(Ig)
            sage: C = Ideals(R)
            sage: C(Ig)
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: I == C(I)
            True
        """
        if v in self:
            return v
        return self.ring().ideal(v)
