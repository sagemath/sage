"""
Specific category classes

This is placed in a separate file from categories.py to avoid circular imports
(as morphisms must be very low in the hierarchy with the new coercion model).
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex
from sage.misc.unknown import Unknown
from .category import JoinCategory, Category, CategoryWithParameters
from sage.misc.lazy_import import lazy_import
lazy_import('sage.categories.objects', 'Objects')

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
        Category.__init__(self)
        self.__object = object

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

        EXAMPLES::

            sage: Elements(ZZ)
            Category of elements of Integer Ring
        """
        from sage.rings.rational_field import QQ
        return cls(QQ)

    def _call_(self, x):
        """
        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: x = V.0
            sage: C = x.category()
            sage: C
            Category of elements of Vector space of dimension 3 over Rational Field
            sage: w = C([1,2,3]); w # indirect doctest
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

        TODO:

        check that this is what we want.
        """
        return [Objects()]

    def object(self):
        return self.__object

    def __reduce__(self):
        return Elements, (self.__object, )

    def _repr_object_names(self):
        return "elements of %s"%self.object()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: x = V.0
            sage: latex(x.category()) # indirect doctest
            \mathbf{Elt}_{\Bold{Q}^{3}}
        """
        return "\\mathbf{Elt}_{%s}"%latex(self.__object)


#############################################################
# Category of objects over some base object
#############################################################
class Category_over_base(CategoryWithParameters):
    r"""
    A base class for categories over some base object

    INPUT:

    - ``base`` -- a category `C` or an object of such a category

    Assumption: the classes for the parents, elements, morphisms, of
    ``self`` should only depend on `C`. See :trac:`11935` for details.

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
        Category.__init__(self, name)

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
            tester.assert_(isinstance(cat, (Category_singleton, Category_over_base,
                                            Bimodules, Schemes)),
                           "The super categories of a category over base should"
                           " be a category over base (or the related Bimodules)"
                           " or a singleton category")

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        Since :trac:`11935`, the element and parent classes of a
        category over base only depend on the category of the base (or
        the base itself if it is a category).

        .. SEEALSO::

            - :meth:`CategoryWithParameters`
            - :meth:`CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Modules(ZZ)._make_named_class_key('element_class')
            Join of Category of euclidean domains
             and Category of infinite enumerated sets
             and Category of metric spaces
            sage: Modules(QQ)._make_named_class_key('parent_class')
            Join of Category of quotient fields and Category of metric spaces
            sage: Schemes(Spec(ZZ))._make_named_class_key('parent_class')
            Category of schemes
            sage: ModularAbelianVarieties(QQ)._make_named_class_key('parent_class')
            Join of Category of quotient fields and Category of metric spaces
            sage: Algebras(Fields())._make_named_class_key('morphism_class')
            Category of fields
        """
        if isinstance(self.__base, Category):
            return self.__base
        return self.__base.category()

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

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
            'algebras over (finite fields and subquotients of monoids and quotients of semigroups)'
        """
        base = self.__base
        if isinstance(base, Category):
            if isinstance(base, JoinCategory):
                name = '('+' and '.join(C._repr_object_names() for C in base.super_categories())+')'
            else:
                name = base._repr_object_names()
        else:
            name = base
        return Category._repr_object_names(self) + " over %s"%name

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(ModulesWithBasis(ZZ))
            \mathbf{ModulesWithBasis}_{\Bold{Z}}
        """
        return "\\mathbf{%s}_{%s}"%(self._label, latex(self.__base))

    def _subcategory_hook_(self, C):
        """
        A quick test whether a category ``C`` may be a subcategory of
        this category.

        INPUT:

        - ``C`` -- a category (type not tested)

        OUTPUT:

        A boolean if it is certain that ``C`` is (or is not) a
        subcategory of self. :obj:`~sage.misc.unknown.Unknown`
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

        The answer is unknown if ``C`` is not a category over base::

            sage: VectorSpaces(QQ)._subcategory_hook_(VectorSpaces(QQ) & Rings())
            Unknown

        Otherwise, the answer is ``True`` in the three following
        cases, and ``False`` otherwise.

        Case 1: the two bases are categories, and the base of ``C`` is
        a subcategory of the base of ``self``::

            sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(Fields()))
            True
            sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(Fields().Finite())) # todo: not implemented
            True
            sage: VectorSpaces(Fields().Finite())._subcategory_hook_(Algebras(Fields()))
            False

        Case 2: the base of ``self`` is a category, and the base of
        ``C`` is a parent in this category::

            sage: VectorSpaces(Fields())._subcategory_hook_(Algebras(QQ))                # todo: not implemented
            True
            sage: VectorSpaces(Fields().Finite())._subcategory_hook_(Algebras(QQ))
            False

        Case 3: the two bases are parents and coincide::

            sage: VectorSpaces(QQ)._subcategory_hook_(Algebras(QQ))
            True
            sage: VectorSpaces(CC)._subcategory_hook_(Algebras(QQ))       # base ring in different categoriess
            False
            sage: VectorSpaces(GF(2))._subcategory_hook_(Algebras(GF(3))) # base ring in the same category
            Unknown

        In this last example, it would be better to return ``False`` ;
        however this is only guaranteed correct if e.g. ``self`` and
        ``C`` are both categories over a base ring. Here is an example
        where this is not the case::

            sage: Sym = SymmetricFunctions(QQ)
            sage: from sage.combinat.sf.sfa import SymmetricFunctionsBases
            sage: SymmetricFunctionsBases(Sym).is_subcategory(Modules(QQ))
            True
            sage: Modules(QQ)._subcategory_hook_(SymmetricFunctionsBases(Sym))
            Unknown

        Maybe such situations should be forbidden.

        """
        if not issubclass(C.parent_class, self.parent_class):
            return False
        if not isinstance(C, Category_over_base):
            return Unknown
        if C.base() is self.__base:
            return True
        if isinstance(self.__base, Category):
            if isinstance(C.base(), Category):
                return C.base().is_subcategory(self.__base)
            # else C.base() is a parent
            return C.base() in self.__base
        return Unknown

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
        if not (base in Rings or
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
                return super(Category_over_base_ring, self).__contains__(x)
        except AttributeError:
            return False


#############################################################
# Category of objects in some ambient object
#############################################################
class Category_in_ambient(Category):
    def __init__(self, ambient, name=None):
        self.__ambient = ambient
        Category.__init__(self, name)

    def ambient(self):
        """
        Return the ambient object in which objects of this category are
        embedded.
        """
        return self.__ambient

    def _repr_(self):
        return Category._repr_(self) + " in %s"%self.__ambient

#    def construction(self):
#        return (self.__class__, self.__ambient)

class Category_module(AbelianCategory, Category_over_base_ring):
    pass

class Category_ideal(Category_in_ambient):

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

        EXAMPLES::

            sage: AlgebraIdeals.an_instance()
            Category of algebra ideals in Univariate Polynomial Ring in x over Rational Field
        """
        from sage.rings.rational_field import QQ
        return cls(QQ['x'])

    def ring(self):
        return self.ambient()

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: C = Ideals(IntegerRing())
            sage: IntegerRing().zero_ideal() in C
            True
        """
        if super(Category_ideal, self).__contains__(x):
            return True
        from sage.rings.ideal import is_Ideal
        if is_Ideal(x) and x.ring() == self.ring():
            return True
        return False

    def __call__(self, v):
        if v in self:
            return v
        return self.ring().ideal(v)

# TODO: make this into a better category
#############################################################
# ChainComplex
#############################################################
class ChainComplexes(Category_module):
    """
    The category of all chain complexes over a base ring.

    EXAMPLES::

        sage: ChainComplexes(RationalField())
        Category of chain complexes over Rational Field

        sage: ChainComplexes(Integers(9))
        Category of chain complexes over Ring of integers modulo 9

     TESTS::

        sage: TestSuite(ChainComplexes(RationalField())).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: ChainComplexes(Integers(9)).super_categories()
            [Category of modules over Ring of integers modulo 9]
        """
        from sage.categories.all import Fields, Modules, VectorSpaces
        base_ring = self.base_ring()
        if base_ring in Fields():
            return [VectorSpaces(base_ring)]
        return [Modules(base_ring)]

