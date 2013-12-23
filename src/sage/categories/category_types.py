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

from sage.misc.latex import latex
from sage.misc.unknown import Unknown
from category import Category, CategoryWithParameters
from objects import Objects

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
# Category of sequences of elements of objects
#############################################################
class Sequences(Category):
    """
    The category of sequences of elements of a given object.

    This category is deprecated.

    EXAMPLES::

        sage: v = Sequence([1,2,3]); v
        [1, 2, 3]
        sage: C = v.category(); C
        Category of sequences in Integer Ring
        sage: loads(C.dumps()) == C
        True
        sage: Sequences(ZZ) is C
        True

        True
        sage: Sequences(ZZ).category()
        Category of objects
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

    def super_categories(self):
        """
        EXAMPLES::

            sage: Sequences(ZZ).super_categories()
            [Category of objects]
        """
        return [Objects()]  # Almost FiniteEnumeratedSets() except for possible repetitions

    def _call_(self, x):
        """
        EXAMPLES::

            sage: v = Sequence([1,2,3]); v
            [1, 2, 3]
            sage: C = v.category(); C
            Category of sequences in Integer Ring
            sage: w = C([2/1, 3/1, GF(2)(5)]); w # indirect doctest
            [2, 3, 1]
            sage: w.category()
            Category of sequences in Integer Ring
        """
        from sage.structure.sequence import Sequence
        return Sequence(x, self.__object)

    def object(self):
        return self.__object

    def __reduce__(self):
        return Sequences, (self.__object, )

    def _repr_object_names(self):
        return "sequences in %s"%self.object()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: v = Sequence([1,2,3])
            sage: latex(v.category()) # indirect doctest
            \mathbf{Seq}_{\Bold{Z}}
        """

        return "\\mathbf{Seq}_{%s}"%latex(self.__object)

#############################################################
# Category of objects over some base object
#############################################################
class Category_over_base(CategoryWithParameters):
    def __init__(self, base, name=None):
        self.__base = base
        Category.__init__(self, name)

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        Since :trac:`11935`, the element and parent classes of a
        category over base only depend on the categories of the base
        ring.

        .. SEEALSO::

            - :meth:`CategoryWithParameters`
            - :meth:`CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Modules(ZZ)._make_named_class_key('element_class')
            Category of euclidean domains
            sage: Modules(QQ)._make_named_class_key('parent_class')
            Category of quotient fields
            sage: Schemes(Spec(ZZ))._make_named_class_key('parent_class')
            Category of Schemes
            sage: ModularAbelianVarieties(QQ)._make_named_class_key('parent_class')
            Category of quotient fields
        """
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
        return Category._repr_object_names(self) + " over %s"%self.__base

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(ModulesWithBasis(QQ))
            \mathbf{ModulesWithBasis}_{\Bold{Q}}
        """
        return "\\mathbf{%s}_{%s}"%(self._label, latex(self.__base))

    def _subcategory_hook_(self, C):
        """
        A quick test whether a category ``C`` may be subcategory of
        this category.

        INPUT:

        - ``C`` -- a category (type not tested)

        OUTPUT:

        A boolean if it is certain that ``C`` is (or is not) a
        subcategory of self. :obj:`~sage.misc.unknown.Unknown`
        otherwise.

        EXAMPLES::

            sage: Algebras(QQ)._subcategory_hook_(HopfAlgebras(QQ))
            True
            sage: Algebras(QQ)._subcategory_hook_(HopfAlgebras(ZZ))
            False
            sage: VectorSpaces(QQ)._subcategory_hook_(VectorSpaces(QQ).hom_category())
            True
            sage: VectorSpaces(QQ)._subcategory_hook_(Category.join([VectorSpaces(QQ).hom_category(),Rings()]))
            Unknown

        """
        if not issubclass(C.parent_class, self.parent_class):
            return False
        try:
            if C.base() is self.__base:
                return True
        except AttributeError:
            pass
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
        from sage.categories.rings import Rings
        assert base in Rings(), "base must be a ring"
        Category_over_base.__init__(self, base, name)

    def base_ring(self):
        """
        Return the base ring over which elements of this category are
        defined.
        """
        return self.base()

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

#############################################################
# TODO: make those two into real categories (with super_category, ...)

# SimplicialComplex
#############################################################
class SimplicialComplexes(Category):
    """
    The category of simplicial complexes.

    EXAMPLES::

        sage: SimplicialComplexes()
        Category of simplicial complexes

    TESTS::

        sage: TestSuite(SimplicialComplexes()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: SimplicialComplexes().super_categories()
            [Category of objects]
        """
        return [Objects()] # anything better?

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
            [Category of modules with basis over Ring of integers modulo 9]
        """
        from sage.categories.all import Fields, FreeModules, VectorSpaces
        base_ring = self.base_ring()
        if base_ring in Fields():
            return [VectorSpaces(base_ring)]
        return [FreeModules(base_ring)]

