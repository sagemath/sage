# sage_setup: distribution = sagemath-objects
r"""
Homsets

The class :class:`Hom` is the base class used to represent sets of morphisms
between objects of a given category.
:class:`Hom` objects are usually "weakly" cached upon creation so that they
don't have to be generated over and over but can be garbage collected together
with the corresponding objects when these are not strongly ref'ed anymore.

EXAMPLES:

In the following, the :class:`Hom` object is indeed cached::

    sage: K = GF(17)
    sage: H = Hom(ZZ, K)
    sage: H
    Set of Homomorphisms from Integer Ring to Finite Field of size 17
    sage: H is Hom(ZZ, K)
    True

Nonetheless, garbage collection occurs when the original references are
overwritten::

    sage: # needs sage.libs.pari
    sage: for p in prime_range(200):
    ....:     K = GF(p)
    ....:     H = Hom(ZZ, K)
    sage: import gc
    sage: _ = gc.collect()
    sage: from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn as FF
    sage: L = [x for x in gc.get_objects() if isinstance(x, FF)]
    sage: len(L)
    1
    sage: L
    [Finite Field of size 199]

AUTHORS:

- David Kohel and William Stein

- David Joyner (2005-12-17): added examples

- William Stein (2006-01-14): Changed from Homspace to Homset.

- Nicolas M. Thiery (2008-12-): Updated for the new category framework

- Simon King (2011-12): Use a weak cache for homsets

- Simon King (2013-02): added examples
"""

# ****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>,
#                     William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories import morphism
from sage.categories.category import Category, JoinCategory
from sage.misc.fast_methods import WithEqualityById
from sage.misc.lazy_attribute import lazy_attribute

###################################
# Use the weak "triple" dictionary
# introduced in github issue #715
# with weak values, as introduced in
# github issue #14159
from sage.structure.coerce_dict import TripleDict
from sage.structure.dynamic_class import dynamic_class
from sage.structure.parent import Parent, Set_generic
from sage.structure.unique_representation import UniqueRepresentation

_cache = TripleDict(weak_values=True)


def Hom(X, Y, category=None, check=True):
    """
    Create the space of homomorphisms from X to Y in the category ``category``.

    INPUT:

    - ``X`` -- an object of a category

    - ``Y`` -- an object of a category

    - ``category`` -- a category in which the morphisms must be
      (default: the meet of the categories of ``X`` and ``Y``);
      both ``X`` and ``Y`` must belong to that category

    - ``check`` -- boolean (default: ``True``); whether to check the
      input, and in particular that ``X`` and ``Y`` belong to
      ``category``.

    OUTPUT: a homset in category

    EXAMPLES::

        sage: V = VectorSpace(QQ, 3)                                                    # needs sage.modules
        sage: Hom(V, V)                                                                 # needs sage.modules
        Set of Morphisms (Linear Transformations) from
        Vector space of dimension 3 over Rational Field to
        Vector space of dimension 3 over Rational Field
        sage: G = AlternatingGroup(3)                                                   # needs sage.groups
        sage: Hom(G, G)                                                                 # needs sage.groups
        Set of Morphisms
         from Alternating group of order 3!/2 as a permutation group
           to Alternating group of order 3!/2 as a permutation group
           in Category of finite enumerated permutation groups
        sage: Hom(ZZ, QQ, Sets())
        Set of Morphisms from Integer Ring to Rational Field in Category of sets

        sage: Hom(FreeModule(ZZ, 1), FreeModule(QQ, 1))                                 # needs sage.modules
        Set of Morphisms
         from Ambient free module of rank 1 over the principal ideal domain Integer Ring
           to Vector space of dimension 1 over Rational Field
           in Category of commutative additive groups
        sage: Hom(FreeModule(QQ, 1), FreeModule(ZZ, 1))                                 # needs sage.modules
        Set of Morphisms
         from Vector space of dimension 1 over Rational Field
           to Ambient free module of rank 1 over the principal ideal domain Integer Ring
           in Category of commutative additive groups

    Here, we test against a memory leak that has been fixed at :issue:`11521` by
    using a weak cache::

        sage: # needs sage.libs.pari
        sage: for p in prime_range(10^3):
        ....:  K = GF(p)
        ....:  a = K(0)
        sage: import gc
        sage: gc.collect()       # random
        624
        sage: from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn as FF
        sage: L = [x for x in gc.get_objects() if isinstance(x, FF)]
        sage: len(L), L[0]
        (1, Finite Field of size 997)

    To illustrate the choice of the category, we consider the
    following parents as running examples::

        sage: X = ZZ; X
        Integer Ring
        sage: Y = SymmetricGroup(3); Y                                                  # needs sage.groups
        Symmetric group of order 3! as a permutation group

    By default, the smallest category containing both ``X`` and ``Y``,
    is used::

        sage: Hom(X, Y)                                                                 # needs sage.groups
        Set of Morphisms from Integer Ring
         to Symmetric group of order 3! as a permutation group
         in Category of enumerated monoids

    Otherwise, if ``category`` is specified, then ``category`` is used,
    after checking that ``X`` and ``Y`` are indeed in ``category``::

        sage: Hom(X, Y, Magmas())                                                       # needs sage.groups
        Set of Morphisms
         from Integer Ring
         to Symmetric group of order 3! as a permutation group
         in Category of magmas

        sage: Hom(X, Y, Groups())                                                       # needs sage.groups
        Traceback (most recent call last):
        ...
        ValueError: Integer Ring is not in Category of groups

    A parent (or a parent class of a category) may specify how to
    construct certain homsets by implementing a method ``_Hom_(self,
    codomain, category)``. This method should either construct the
    requested homset or raise a :exc:`TypeError`. This hook is currently
    mostly used to create homsets in some specific subclass of
    :class:`Homset` (e.g. :class:`sage.rings.homset.RingHomset`)::

        sage: Hom(QQ,QQ).__class__
        <class 'sage.rings.homset.RingHomset_generic_with_category'>

    Do not call this hook directly to create homsets, as it does not
    handle unique representation::

        sage: Hom(QQ,QQ) == QQ._Hom_(QQ, category=QQ.category())
        True
        sage: Hom(QQ,QQ) is QQ._Hom_(QQ, category=QQ.category())
        False

    TESTS:

    Homset are unique parents::

        sage: k = GF(5)
        sage: H1 = Hom(k, k)
        sage: H2 = Hom(k, k)
        sage: H1 is H2
        True

    Moreover, if no category is provided, then the result is identical
    with the result for the meet of the categories of the domain and
    the codomain::

        sage: Hom(QQ, ZZ) is Hom(QQ,ZZ, Category.meet([QQ.category(), ZZ.category()]))
        True

    Some doc tests in :mod:`sage.rings` (need to) break the unique
    parent assumption. But if domain or codomain are not unique
    parents, then the homset will not fit. That is to say, the hom set
    found in the cache will have a (co)domain that is equal to, but
    not identical with, the given (co)domain.

    By :issue:`9138`, we abandon the uniqueness of homsets, if the
    domain or codomain break uniqueness::

        sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
        sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ, 3, order='degrevlex')
        sage: Q.<x,y,z>=MPolynomialRing_polydict_domain(QQ, 3, order='degrevlex')
        sage: P == Q
        True
        sage: P is Q
        False

    Hence, ``P`` and ``Q`` are not unique parents. By consequence, the
    following homsets aren't either::

        sage: H1 = Hom(QQ,P)
        sage: H2 = Hom(QQ,Q)
        sage: H1 == H2
        True
        sage: H1 is H2
        False

    It is always the most recently constructed homset that remains in
    the cache::

        sage: H2 is Hom(QQ,Q)
        True

    Variation on the theme::

        sage: # needs sage.modules
        sage: U1 = FreeModule(ZZ, 2)
        sage: U2 = FreeModule(ZZ, 2, inner_product_matrix=matrix([[1,0], [0,-1]]))
        sage: U1 == U2, U1 is U2
        (False, False)
        sage: V = ZZ^3
        sage: H1 = Hom(U1, V); H2 = Hom(U2, V)
        sage: H1 == H2, H1 is H2
        (False, False)
        sage: H1 = Hom(V, U1); H2 = Hom(V, U2)
        sage: H1 == H2, H1 is H2
        (False, False)

    Since :issue:`11900`, the meet of the categories of the given arguments is
    used to determine the default category of the homset. This can also be a
    join category, as in the following example::

        sage: PA = Parent(category=Algebras(QQ))
        sage: PJ = Parent(category=Rings() & Modules(QQ))
        sage: Hom(PA, PJ)
        Set of Homomorphisms
         from <sage.structure.parent.Parent object at ...>
           to <sage.structure.parent.Parent object at ...>
        sage: Hom(PA, PJ).category()
        Category of homsets of
         unital magmas and right modules over Rational Field
          and left modules over Rational Field
        sage: Hom(PA, PJ, Rngs())
        Set of Morphisms
         from <sage.structure.parent.Parent object at ...>
           to <sage.structure.parent.Parent object at ...> in Category of rngs

    .. TODO::

        - Design decision: how much of the homset comes from the
          category of ``X`` and ``Y``, and how much from the specific
          ``X`` and ``Y``.  In particular, do we need several parent
          classes depending on ``X`` and ``Y``, or does the difference
          only lie in the elements (i.e.  the morphism), and of course
          how the parent calls their constructors.
        - Specify the protocol for the ``_Hom_`` hook in case of ambiguity
          (e.g. if both a parent and some category thereof provide one).

    TESTS:

    Facade parents over plain Python types are supported::

        sage: from sage.sets.pythonclass import Set_PythonType
        sage: R = Set_PythonType(int)
        sage: S = Set_PythonType(float)
        sage: Hom(R, S)
        Set of Morphisms from Set of Python objects of class 'int' to Set of Python objects of class 'float' in Category of sets

    Checks that the domain and codomain are in the specified
    category. Case of a non parent::

        sage: # needs sage.graphs
        sage: S = SimplicialComplex([[1,2], [1,4]]); S.rename('S')
        sage: Hom(S, S, SimplicialComplexes())
        Set of Morphisms from S to S in Category of finite simplicial complexes
        sage: Hom(Set(), S, Sets())
        Set of Morphisms from {} to S in Category of sets
        sage: Hom(S, Set(), Sets())
        Set of Morphisms from S to {} in Category of sets
        sage: H = Hom(S, S, ChainComplexes(QQ))
        Traceback (most recent call last):
        ...
        ValueError: S is not in Category of chain complexes over Rational Field

    Those checks are done with the natural idiom ``X in category``,
    and not ``X.category().is_subcategory(category)`` as it used to be
    before :issue:`16275` (see :issue:`15801` for a real use case)::

        sage: # needs sage.graphs
        sage: class PermissiveCategory(Category):
        ....:     def super_categories(self): return [Objects()]
        ....:     def __contains__(self, X): return True
        sage: C = PermissiveCategory(); C.rename('Permissive category')
        sage: S.category().is_subcategory(C)
        False
        sage: S in C
        True
        sage: Hom(S, S, C)
        Set of Morphisms from S to S in Permissive category

    With ``check=False``, uninitialized parents, as can appear upon
    unpickling, are supported. Case of a parent::

        sage: cls = type(Set())
        sage: S = unpickle_newobj(cls, ())  # A non parent
        sage: H = Hom(S, S, SimplicialComplexes(), check=False)
        sage: H = Hom(S, S, Sets(),                check=False)
        sage: H = Hom(S, S, ChainComplexes(QQ),    check=False)

    Case of a non parent::

        sage: # needs sage.graphs
        sage: cls = type(SimplicialComplex([[1,2], [1,4]]))
        sage: S = unpickle_newobj(cls, ())
        sage: H = Hom(S, S, Sets(),                check=False)
        sage: H = Hom(S, S, Groups(),              check=False)
        sage: H = Hom(S, S, SimplicialComplexes(), check=False)

    Typical example where unpickling involves calling Hom on an
    uninitialized parent::

        sage: P.<x,y> = QQ['x,y']
        sage: Q = P.quotient([x^2 - 1, y^2 - 1])
        sage: q = Q.an_element()                                                        # needs sage.libs.singular
        sage: explain_pickle(dumps(Q))                                                  # needs sage.libs.singular
        pg_...
        ... = pg_dynamic_class('QuotientRing_generic_with_category', (pg_QuotientRing_generic, pg_getattr(..., 'parent_class')), None, None, pg_QuotientRing_generic)
        si... = unpickle_newobj(..., ())
        ...
        si... = pg_unpickle_MPolynomialRing_libsingular(..., ('x', 'y'), ...)
        si... = ... pg_Hom(si..., si..., ...) ...
        sage: Q == loads(dumps(Q))
        True

    Check that the ``_Hom_`` method of the ``category`` input is used::

        sage: from sage.categories.category_types import Category_over_base_ring
        sage: class ModulesWithHom(Category_over_base_ring):
        ....:     def super_categories(self):
        ....:         return [Modules(self.base_ring())]
        ....:     class ParentMethods:
        ....:         def _Hom_(self, Y, category=None):
        ....:             print("Modules")
        ....:             raise TypeError
        sage: class AlgebrasWithHom(Category_over_base_ring):
        ....:     def super_categories(self):
        ....:         return [Algebras(self.base_ring()), ModulesWithHom(self.base_ring())]
        ....:     class ParentMethods:
        ....:         def _Hom_(self, Y, category=None):
        ....:             R = self.base_ring()
        ....:             if category is not None and category.is_subcategory(Algebras(R)):
        ....:                 print("Algebras")
        ....:             raise TypeError
        sage: from sage.structure.element import Element
        sage: class Foo(Parent):
        ....:     def _coerce_map_from_base_ring(self):
        ....:         return self._generic_coerce_map(self.base_ring())
        ....:     class Element(Element):
        ....:         pass
        sage: X = Foo(base=QQ, category=AlgebrasWithHom(QQ))
        sage: H = Hom(X, X, ModulesWithHom(QQ))
        Modules
    """
    # This should use cache_function instead
    # However some special handling is currently needed for
    # domains/docomains that break the unique parent condition. Also,
    # at some point, it somehow broke the coercion (see e.g. sage -t
    # sage.rings.real_mpfr). To be investigated.
    global _cache
    key = (X,Y,category)
    try:
        H = _cache[key]
    except KeyError:
        H = None
    if H is not None:
        # Return H unless the domain or codomain breaks the unique parent condition
        if H.domain() is X and H.codomain() is Y:
            return H

    # Determines the category
    if category is None:
        category = X.category()._meet_(Y.category())
        # Recurse to make sure that Hom(X, Y) and Hom(X, Y, category) are identical
        # No need to check the input again
        H = Hom(X, Y, category, check=False)
    else:
        if check:
            if not isinstance(category, Category):
                raise TypeError("Argument category (= {}) must be a category.".format(category))
            for O in [X, Y]:
                try:
                    category_mismatch = O not in category
                except Exception:
                    # An error should not happen, this here is just to be on
                    # the safe side.
                    category_mismatch = True
                # A category mismatch does not necessarily mean that an error
                # should be raised. Instead, it could be the case that we are
                # unpickling an old pickle (that doesn't set the "check"
                # argument to False). In this case, it could be that the
                # (co)domain is not properly initialised, which we are
                # checking now. See trac #16275 and #14793.
                if category_mismatch and O._is_category_initialized():
                    # At this point, we can be rather sure that O is properly
                    # initialised, and thus its string representation is
                    # available for the following error message. It simply
                    # belongs to the wrong category.
                    raise ValueError("{} is not in {}".format(O, category))

        # Construct H
        try: # _Hom_ hook from the parent
            H = X._Hom_(Y, category)
        except (AttributeError, TypeError):
            # Workaround in case the above fails, but the category
            # also provides a _Hom_ hook.
            # FIXME:
            # - If X._Hom_ actually comes from category and fails, it
            #   will be called twice.
            # - This is bound to fail if X is an extension type and
            #   does not actually inherit from category.parent_class
            # For join categories, we check all of the direct super
            #   categories as the parent_class of the join category is
            #   not (necessarily) inherited and join categories do not
            #   implement a _Hom_ (see trac #23418).
            if not isinstance(category, JoinCategory):
                cats = [category]
            else:
                cats = category.super_categories()
            H = None
            for C in cats:
                try:
                    H = C.parent_class._Hom_(X, Y, category=category)
                    break
                except (AttributeError, TypeError):
                    pass
            if H is None:
                # By default, construct a plain homset.
                H = Homset(X, Y, category=category, check=check)
    _cache[key] = H
    if isinstance(X, UniqueRepresentation) and isinstance(Y, UniqueRepresentation):
        if not isinstance(H, WithEqualityById):
            try:
                H.__class__ = dynamic_class(H.__class__.__name__+"_with_equality_by_id", (WithEqualityById, H.__class__), doccls=H.__class__)
            except Exception:
                pass
    return H


def hom(X, Y, f):
    """
    Return ``Hom(X,Y)(f)``, where ``f`` is data that defines an element of
    ``Hom(X,Y)``.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: phi = hom(R, QQ, [2])
        sage: phi(x^2 + 3)
        7
    """
    return Hom(X,Y)(f)


def End(X, category=None):
    r"""
    Create the set of endomorphisms of ``X`` in the category category.

    INPUT:

    - ``X`` -- anything

    - ``category`` -- (optional) category in which to coerce ``X``

    OUTPUT: a set of endomorphisms in category

    EXAMPLES::

        sage: V = VectorSpace(QQ, 3)                                                    # needs sage.modules
        sage: End(V)                                                                    # needs sage.modules
        Set of Morphisms (Linear Transformations)
         from Vector space of dimension 3 over Rational Field
         to Vector space of dimension 3 over Rational Field

    ::

        sage: # needs sage.groups
        sage: G = AlternatingGroup(3)
        sage: S = End(G); S
        Set of Morphisms
         from Alternating group of order 3!/2 as a permutation group
         to Alternating group of order 3!/2 as a permutation group
         in Category of finite enumerated permutation groups
        sage: S.domain()
        Alternating group of order 3!/2 as a permutation group

    To avoid creating superfluous categories, a homset in a category
    ``Cs()`` is in the homset category of the lowest full super category
    ``Bs()`` of ``Cs()`` that implements ``Bs.Homsets`` (or the join
    thereof if there are several). For example, finite groups form a
    full subcategory of unital magmas: any unital magma morphism
    between two finite groups is a finite group morphism. Since finite
    groups currently implement nothing more than unital magmas about
    their homsets, we have::

        sage: # needs sage.groups
        sage: G = GL(3, 3)
        sage: G.category()
        Category of finite groups
        sage: H = Hom(G, G)
        sage: H.homset_category()
        Category of finite groups
        sage: H.category()
        Category of endsets of unital magmas

    Similarly, a ring morphism just needs to preserve addition,
    multiplication, zero, and one. Accordingly, and since the category
    of rings implements nothing specific about its homsets, a ring
    homset is currently constructed in the category of homsets of
    unital magmas and unital additive magmas::

        sage: H = Hom(ZZ,ZZ,Rings())
        sage: H.category()
        Category of endsets of unital magmas and additive unital additive magmas
    """
    return Hom(X,X, category)


def end(X, f):
    """
    Return ``End(X)(f)``, where ``f`` is data that defines an element of
    ``End(X)``.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: phi = end(R, [x + 1])
        sage: phi
        Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
          Defn: x |--> x + 1
        sage: phi(x^2 + 5)
        x^2 + 2*x + 6
    """
    return End(X)(f)


class Homset(Set_generic):
    """
    The class for collections of morphisms in a category.

    EXAMPLES::

        sage: H = Hom(QQ^2, QQ^3)                                                       # needs sage.modules
        sage: loads(H.dumps()) is H                                                     # needs sage.modules
        True

    Homsets of unique parents are unique as well::

        sage: H = End(AffineSpace(2, names='x,y'))
        sage: loads(dumps(AffineSpace(2, names='x,y'))) is AffineSpace(2, names='x,y')
        True
        sage: loads(dumps(H)) is H
        True

    Conversely, homsets of non-unique parents are non-unique::

        sage: P11 = ProductProjectiveSpaces(QQ, [1, 1])
        sage: H = End(P11)
        sage: loads(dumps(P11)) is ProductProjectiveSpaces(QQ, [1, 1])
        False
        sage: loads(dumps(P11)) == ProductProjectiveSpaces(QQ, [1, 1])
        True
        sage: loads(dumps(H)) is H
        False
        sage: loads(dumps(H)) == H
        True
    """
    def __init__(self, X, Y, category=None, base=None, check=True):
        r"""
        TESTS::

            sage: X = ZZ['x']; X.rename('X')
            sage: Y = ZZ['y']; Y.rename('Y')
            sage: f = X.hom([0], Y)
            sage: class MyHomset(Homset):
            ....:     def _an_element_(self):
            ....:         return sage.categories.morphism.SetMorphism(self, f)
            sage: import __main__; __main__.MyHomset = MyHomset  # fakes MyHomset being defined in a Python module
            sage: H = MyHomset(X, Y, category=Monoids(), base = ZZ)
            sage: H
            Set of Morphisms from X to Y in Category of monoids
            sage: TestSuite(H).run()

            sage: H = MyHomset(X, Y, category=1, base = ZZ)
            Traceback (most recent call last):
            ...
            TypeError: category (=1) must be a category

            sage: H = MyHomset(X, Y, category=1, base = ZZ, check = False)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.integer.Integer' object has no attribute 'Homsets'...
            sage: P.<t> = ZZ[]
            sage: f = P.hom([1/2*t])
            sage: f.parent().domain()
            Univariate Polynomial Ring in t over Integer Ring
            sage: f.domain() is f.parent().domain()
            True

        Test that ``base_ring`` is initialized properly::

            sage: R = QQ['x']
            sage: Hom(R, R).base_ring()
            Rational Field
            sage: Hom(R, R, category=Sets()).base_ring()
            sage: Hom(R, R, category=Modules(QQ)).base_ring()
            Rational Field
            sage: Hom(QQ^3, QQ^3, category=Modules(QQ)).base_ring()                     # needs sage.modules
            Rational Field

        For whatever it's worth, the ``base`` arguments takes precedence::

            sage: MyHomset(ZZ^3, ZZ^3, base=QQ).base_ring()                             # needs sage.modules
            Rational Field
        """
        self._domain = X
        self._codomain = Y
        if category is None:
            category = X.category()
        self.__category = category
        if check:
            if not isinstance(category, Category):
                raise TypeError("category (=%s) must be a category" % category)
            # if not X in category:
            #     raise TypeError("X (=%s) must be in category (=%s)" % (X, category))
            # if not Y in category:
            #     raise TypeError("Y (=%s) must be in category (=%s)" % (Y, category))

        if base is None and hasattr(category, "WithBasis"):
            # The above is a lame but fast check that category is a
            # subcategory of Modules(...). That will do until
            # CategoryObject.base_ring will be gone and not prevent
            # anymore from implementing base_ring in Modules.Homsets.ParentMethods.
            # See also #15801.
            base = X.base_ring()

        Parent.__init__(self, base=base,
                        category=category.Endsets() if X is Y else category.Homsets())

    def __reduce__(self):
        """
        Implement pickling by construction for Homsets.

        Homsets are unpickled using the function
        :func:`~sage.categories.homset.Hom` which is cached:
        ``Hom(domain, codomain, category, check=False)``.

        .. NOTE::

            It can happen, that ``Hom(X,X)`` is called during
            unpickling with an uninitialized instance ``X`` of a Python
            class. In some of these cases, testing that ``X in
            category`` can trigger ``X.category()``. This in turn can
            raise a error, or return a too large category (``Sets()``,
            for example) and (worse!) assign this larger category to
            the ``X._category`` cdef attribute, so that it would
            subsequently seem that ``X``'s category was initialised.

            Beside speed considerations, this is the main rationale
            for disabling checks upon unpickling.

            .. SEEALSO:: :issue:`14793`, :issue:`16275`

        EXAMPLES::

            sage: H = Hom(QQ^2, QQ^3)                                                   # needs sage.modules
            sage: H.__reduce__()                                                        # needs sage.modules
            (<function Hom at ...>,
             (Vector space of dimension 2 over Rational Field,
              Vector space of dimension 3 over Rational Field,
              Category of finite dimensional vector spaces with basis over
                 (number fields and quotient fields and metric spaces),
              False))

        TESTS::

            sage: loads(H.dumps()) is H                                                 # needs sage.modules
            True

        Homsets of non-unique parents are non-unique as well::

            sage: # needs sage.groups
            sage: G = PermutationGroup([[(1, 2, 3), (4, 5)], [(3, 4)]])
            sage: G is loads(dumps(G))
            False
            sage: H = Hom(G, G)
            sage: H is loads(dumps(H))
            False
            sage: H == loads(dumps(H))
            True
        """
        return Hom, (self._domain, self._codomain, self.__category, False)

    def _repr_(self):
        """
        TESTS::

            sage: Hom(ZZ^2, QQ, category=Sets())._repr_()                               # needs sage.modules
            'Set of Morphisms from Ambient free module of rank 2 over the principal ideal domain Integer Ring to Rational Field in Category of sets'
        """
        return "Set of Morphisms from {} to {} in {}".format(self._domain,
            self._codomain, self.__category)

    def __hash__(self):
        """
        The hash is obtained from domain, codomain and base.

        TESTS::

            sage: hash(Hom(ZZ, QQ)) == hash((ZZ, QQ, ZZ))
            True
            sage: hash(Hom(QQ, ZZ)) == hash((QQ, ZZ, QQ))
            True

            sage: E = EllipticCurve('37a')                                              # needs sage.schemes
            sage: H = E(0).parent(); H                                                  # needs sage.schemes
            Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: hash(H) == hash((H.domain(), H.codomain(), H.base()))                 # needs sage.schemes
            True
        """
        return hash((self._domain, self._codomain, self.base()))

    def __bool__(self) -> bool:
        """
        TESTS::

            sage: bool(Hom(ZZ, QQ))
            True
        """
        return True

    def homset_category(self):
        """
        Return the category that this is a Hom in, i.e., this is typically
        the category of the domain or codomain object.

        EXAMPLES::

            sage: H = Hom(AlternatingGroup(4), AlternatingGroup(7))                     # needs sage.groups
            sage: H.homset_category()                                                   # needs sage.groups
            Category of finite enumerated permutation groups
        """
        return self.__category

    def _element_constructor_(self, x, check=None, **options):
        r"""
        Construct a morphism in this homset from ``x`` if possible.

        EXAMPLES::

            sage: H = Hom(SymmetricGroup(4), SymmetricGroup(7))                         # needs sage.groups
            sage: phi = Hom(SymmetricGroup(5), SymmetricGroup(6)).natural_map()         # needs sage.groups
            sage: phi                                                                   # needs sage.groups
            Coercion morphism:
              From: Symmetric group of order 5! as a permutation group
              To:   Symmetric group of order 6! as a permutation group

        When converting `\phi` into `H`, some coerce maps are applied. Note
        that (in contrast to what is stated in the following string
        representation) it is safe to use the resulting map, since a composite
        map prevents the codomains of all constituent maps from garbage
        collection, if there is a strong reference to its domain (which is the
        case here)::

            sage: H(phi)                                                                # needs sage.groups
            Composite map:
              From: Symmetric group of order 4! as a permutation group
              To:   Symmetric group of order 7! as a permutation group
              Defn:   (map internal to coercion system -- copy before use)
                    Coercion map:
                      From: Symmetric group of order 4! as a permutation group
                      To:   Symmetric group of order 5! as a permutation group
                    then
                      Coercion morphism:
                      From: Symmetric group of order 5! as a permutation group
                      To:   Symmetric group of order 6! as a permutation group
                    then
                      (map internal to coercion system -- copy before use)
                    Coercion map:
                      From: Symmetric group of order 6! as a permutation group
                      To:   Symmetric group of order 7! as a permutation group

        Also note that making a copy of the resulting map will automatically
        make strengthened copies of the composed maps::

            sage: copy(H(phi))                                                          # needs sage.groups
            Composite map:
              From: Symmetric group of order 4! as a permutation group
              To:   Symmetric group of order 7! as a permutation group
              Defn:   Coercion map:
                      From: Symmetric group of order 4! as a permutation group
                      To:   Symmetric group of order 5! as a permutation group
                    then
                      Coercion morphism:
                      From: Symmetric group of order 5! as a permutation group
                      To:   Symmetric group of order 6! as a permutation group
                    then
                      Coercion map:
                      From: Symmetric group of order 6! as a permutation group
                      To:   Symmetric group of order 7! as a permutation group
            sage: H = Hom(ZZ, ZZ, Sets())
            sage: f = H( lambda x: x + 1 )
            sage: f.parent()
            Set of Morphisms from Integer Ring to Integer Ring in Category of sets
            sage: f.domain()
            Integer Ring
            sage: f.codomain()
            Integer Ring
            sage: f(1), f(2), f(3)
            (2, 3, 4)

            sage: H = Hom(Set([1,2,3]), Set([1,2,3]))
            sage: f = H(lambda x: 4 - x)
            sage: f.parent()
            Set of Morphisms from {1, 2, 3} to {1, 2, 3} in Category of finite enumerated sets
            sage: f(1), f(2), f(3) # todo: not implemented

            sage: H = Hom(ZZ, QQ, Sets())
            sage: f = H( ConstantFunction(2/3) )
            sage: f.parent()
            Set of Morphisms from Integer Ring to Rational Field in Category of sets
            sage: f(1), f(2), f(3)
            (2/3, 2/3, 2/3)

        By :issue:`14711`, conversion and coerce maps should be copied
        before using them outside of the coercion system::

            sage: H = Hom(ZZ,QQ['t'], CommutativeAdditiveGroups())
            sage: P.<t> = ZZ[]
            sage: f = P.hom([2*t])
            sage: phi = H._generic_convert_map(f.parent()); phi
            Conversion map:
              From: Set of Homomorphisms from Univariate Polynomial Ring in t over Integer Ring to Univariate Polynomial Ring in t over Integer Ring
              To:   Set of Morphisms from Integer Ring to Univariate Polynomial Ring in t over Rational Field in Category of commutative additive groups
            sage: H._generic_convert_map(f.parent())(f)
            Composite map:
              From: Integer Ring
              To:   Univariate Polynomial Ring in t over Rational Field
              Defn:   (map internal to coercion system -- copy before use)
                    Polynomial base injection morphism:
                      From: Integer Ring
                      To:   Univariate Polynomial Ring in t over Integer Ring
                    then
                      Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
                      Defn: t |--> 2*t
                    then
                      (map internal to coercion system -- copy before use)
                    Ring morphism:
                      From: Univariate Polynomial Ring in t over Integer Ring
                      To:   Univariate Polynomial Ring in t over Rational Field
            sage: copy(H._generic_convert_map(f.parent())(f))
            Composite map:
              From: Integer Ring
              To:   Univariate Polynomial Ring in t over Rational Field
              Defn:   Polynomial base injection morphism:
                      From: Integer Ring
                      To:   Univariate Polynomial Ring in t over Integer Ring
                    then
                      Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
                      Defn: t |--> 2*t
                    then
                      Ring morphism:
                      From: Univariate Polynomial Ring in t over Integer Ring
                      To:   Univariate Polynomial Ring in t over Rational Field
                      Defn: Induced from base ring by
                            Natural morphism:
                              From: Integer Ring
                              To:   Rational Field

        TESTS::

            sage: # needs sage.groups
            sage: G.<x,y,z> = FreeGroup()
            sage: H = Hom(G, G)
            sage: H(H.identity())
            Identity endomorphism of Free Group on generators {x, y, z}
            sage: H()
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 0 to an element of
             Set of Morphisms from Free Group on generators {x, y, z}
             to Free Group on generators {x, y, z} in Category of infinite groups
            sage: H("whatever")
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'whatever' to an element of
             Set of Morphisms from Free Group on generators {x, y, z}
             to Free Group on generators {x, y, z} in Category of infinite groups
            sage: HH = Hom(H, H)
            sage: HH(HH.identity(), foo='bar')
            Traceback (most recent call last):
            ...
            NotImplementedError: no keywords are implemented for
             constructing elements of ...

        AUTHORS:

        - Robert Bradshaw, with changes by Nicolas M. Thiery
        """
        if options:
            # TODO: this is specific for ModulesWithBasis; generalize
            # this to allow homsets and categories to provide more
            # morphism constructors (on_algebra_generators, ...)
            try:
                call_with_keywords = self.__call_on_basis__
            except AttributeError:
                if 'base_map' in options:
                    raise NotImplementedError("base_map not supported for this Homset; you may need to specify a category")
                raise NotImplementedError("no keywords are implemented for constructing elements of {}".format(self))
            options.setdefault("category", self.homset_category())
            return call_with_keywords(**options)

        if isinstance(x, morphism.Morphism):
            if x.domain() != self.domain():
                mor = x.domain()._internal_coerce_map_from(self.domain())
                if mor is None:
                    raise TypeError("Incompatible domains: x (=%s) cannot be an element of %s" % (x,self))
                x = x * mor
            if x.codomain() != self.codomain():
                mor = self.codomain()._internal_coerce_map_from(x.codomain())
                if mor is None:
                    raise TypeError("Incompatible codomains: x (=%s) cannot be an element of %s" % (x,self))
                x = mor * x
            return x

        if callable(x):
            return self.element_class_set_morphism(self, x)

        raise TypeError("unable to convert {!r} to an element of {}".format(x, self))

    @lazy_attribute
    def _abstract_element_class(self):
        """
        An abstract class for the elements of this homset.

        This class is built from the element class of the homset
        category and the morphism class of the category.  This makes
        it possible for a category to provide code for its morphisms
        and for morphisms of all its subcategories, full or not.

        .. NOTE::

            The element class of ``C.Homsets()`` will be inherited by
            morphisms in *full* subcategories of ``C``, while the morphism
            class of ``C`` will be inherited by *all* subcategories of
            ``C``. Hence, if some feature of a morphism depends on the
            algebraic properties of the homsets, it should be implemented by
            ``C.Homsets.ElementMethods``, but if it depends only on the
            algebraic properties of domain and codomain, it should be
            implemented in ``C.MorphismMethods``.

            At this point, the homset element classes take precedence over the
            morphism classes. But this may be subject to change.


        .. TODO::

            - Make sure this class is shared whenever possible.
            - Flatten join category classes

        .. SEEALSO::

            - :meth:`Parent._abstract_element_class`

        EXAMPLES:

        Let's take a homset of finite commutative groups as example; at
        this point this is the simplest one to create (gosh)::

            sage: # needs sage.groups
            sage: cat = Groups().Finite().Commutative()
            sage: C3 = PermutationGroup([(1,2,3)])
            sage: C3._refine_category_(cat)
            sage: C2 = PermutationGroup([(1,2)])
            sage: C2._refine_category_(cat)
            sage: H = Hom(C3, C2, cat)
            sage: H.homset_category()
            Category of finite commutative groups
            sage: H.category()
            Category of homsets of unital magmas
            sage: cls = H._abstract_element_class; cls
            <class 'sage.categories.homsets.GroupHomset_libgap_with_category._abstract_element_class'>
            sage: cls.__bases__ == (H.category().element_class, H.homset_category().morphism_class)
            True

        A morphism of finite commutative semigroups is also a morphism
        of semigroups, of magmas, ...; it thus inherits code from all
        those categories::

            sage: # needs sage.groups
            sage: issubclass(cls, Semigroups().Finite().morphism_class)
            True
            sage: issubclass(cls, Semigroups().morphism_class)
            True
            sage: issubclass(cls, Magmas().Commutative().morphism_class)
            True
            sage: issubclass(cls, Magmas().morphism_class)
            True
            sage: issubclass(cls, Sets().morphism_class)
            True

        Recall that FiniteMonoids() is a full subcategory of
        ``Monoids()``, but not of ``FiniteSemigroups()``. Thus::

            sage: issubclass(cls, Monoids().Finite().Homsets().element_class)           # needs sage.groups
            True
            sage: issubclass(cls, Semigroups().Finite().Homsets().element_class)        # needs sage.groups
            False
        """
        class_name = "%s._abstract_element_class" % self.__class__.__name__
        return dynamic_class(class_name, (self.category().element_class, self.homset_category().morphism_class))

    @lazy_attribute
    def element_class_set_morphism(self):
        """
        A base class for elements of this homset which are
        also ``SetMorphism``, i.e. implemented by mean of a
        Python function.

        This is currently plain ``SetMorphism``, without inheritance
        from categories.

        .. TODO::

            Refactor during the upcoming homset cleanup.

        EXAMPLES::

            sage: H = Hom(ZZ, ZZ)
            sage: H.element_class_set_morphism
            <class 'sage.categories.morphism.SetMorphism'>
        """
        return self.__make_element_class__(morphism.SetMorphism)

    def __eq__(self, other):
        """
        For two homsets, it is tested whether the domain, the codomain and
        the category coincide.

        EXAMPLES::

            sage: H1 = Hom(ZZ,QQ, CommutativeAdditiveGroups())
            sage: H2 = Hom(ZZ,QQ)
            sage: H1 == H2
            False
            sage: H1 == loads(dumps(H1))
            True
        """
        if not isinstance(other, Homset):
            return False
        return (self._domain == other._domain
                and self._codomain == other._codomain
                and self.__category == other.__category)

    def __ne__(self, other):
        """
        Check for not-equality of ``self`` and ``other``.

        EXAMPLES::

            sage: H1 = Hom(ZZ,QQ, CommutativeAdditiveGroups())
            sage: H2 = Hom(ZZ,QQ)
            sage: H3 = Hom(ZZ['t'],QQ, CommutativeAdditiveGroups())
            sage: H4 = Hom(ZZ,QQ['t'], CommutativeAdditiveGroups())
            sage: H1 != H2
            True
            sage: H1 != loads(dumps(H1))
            False
            sage: H1 != H3 != H4 != H1
            True
        """
        return not (self == other)

    def __contains__(self, x):
        """
        Test whether the parent of the argument is ``self``.

        TESTS::

            sage: P.<t> = ZZ[]
            sage: f = P.hom([1/2*t])
            sage: f in Hom(ZZ['t'],QQ['t'])  # indirect doctest
            True
            sage: f in Hom(ZZ['t'],QQ['t'], CommutativeAdditiveGroups())  # indirect doctest
            False
        """
        try:
            return x.parent() == self
        except AttributeError:
            pass
        return False

    def natural_map(self):
        """
        Return the "natural map" of this homset.

        .. NOTE::

            By default, a formal coercion morphism is returned.

        EXAMPLES::

            sage: H = Hom(ZZ['t'],QQ['t'], CommutativeAdditiveGroups())
            sage: H.natural_map()
            Coercion morphism:
              From: Univariate Polynomial Ring in t over Integer Ring
              To:   Univariate Polynomial Ring in t over Rational Field
            sage: H = Hom(QQ['t'], GF(3)['t'])
            sage: H.natural_map()
            Traceback (most recent call last):
            ...
            TypeError: natural coercion morphism
            from Univariate Polynomial Ring in t over Rational Field
            to Univariate Polynomial Ring in t over Finite Field of size 3 not defined
        """
        return morphism.FormalCoercionMorphism(self)   # good default in many cases

    def identity(self):
        """
        The identity map of this homset.

        .. NOTE::

            Of course, this only exists for sets of endomorphisms.

        EXAMPLES::

            sage: H = Hom(QQ,QQ)
            sage: H.identity()
            Identity endomorphism of Rational Field
            sage: H = Hom(ZZ,QQ)
            sage: H.identity()
            Traceback (most recent call last):
            ...
            TypeError: identity map only defined for endomorphisms; try natural_map() instead
            sage: H.natural_map()
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        if self.is_endomorphism_set():
            return morphism.IdentityMorphism(self)
        raise TypeError("identity map only defined for endomorphisms; try natural_map() instead")

    def one(self):
        """
        The identity map of this homset.

        .. NOTE::

            Of course, this only exists for sets of endomorphisms.

        EXAMPLES::

            sage: K = GaussianIntegers()                                                # needs sage.rings.number_field
            sage: End(K).one()                                                          # needs sage.rings.number_field
            Identity endomorphism of Gaussian Integers generated by I
             in Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        return self.identity()

    def domain(self):
        """
        Return the domain of this homset.

        EXAMPLES::

            sage: P.<t> = ZZ[]
            sage: f = P.hom([1/2*t])
            sage: f.parent().domain()
            Univariate Polynomial Ring in t over Integer Ring
            sage: f.domain() is f.parent().domain()
            True
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of this homset.

        EXAMPLES::

            sage: P.<t> = ZZ[]
            sage: f = P.hom([1/2*t])
            sage: f.parent().codomain()
            Univariate Polynomial Ring in t over Rational Field
            sage: f.codomain() is f.parent().codomain()
            True
        """
        return self._codomain

    def reversed(self):
        """
        Return the corresponding homset, but with the domain and codomain
        reversed.

        EXAMPLES::

            sage: # needs sage.modules
            sage: H = Hom(ZZ^2, ZZ^3); H
            Set of Morphisms from Ambient free module of rank 2 over
             the principal ideal domain Integer Ring to Ambient free module
             of rank 3 over the principal ideal domain Integer Ring in
             Category of finite dimensional modules with basis over (Dedekind
             domains and euclidean domains and noetherian rings
             and infinite enumerated sets and metric spaces)
            sage: type(H)
            <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
            sage: H.reversed()
            Set of Morphisms from Ambient free module of rank 3 over
             the principal ideal domain Integer Ring to Ambient free module
             of rank 2 over the principal ideal domain Integer Ring in
             Category of finite dimensional modules with basis over (Dedekind
             domains and euclidean domains and noetherian rings
             and infinite enumerated sets and metric spaces)
            sage: type(H.reversed())
            <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        """
        return Hom(self.codomain(), self.domain(),
                   category=self.homset_category())


# Really needed???
class HomsetWithBase(Homset):
    def __init__(self, X, Y, category=None, check=True, base=None):
        r"""
        TESTS::

            sage: X = ZZ['x']; X.rename('X')
            sage: Y = ZZ['y']; Y.rename('Y')
            sage: f = X.hom([0], Y)
            sage: class MyHomset(HomsetWithBase):
            ....:     def _an_element_(self):
            ....:         return sage.categories.morphism.SetMorphism(self, f)
            sage: import __main__; __main__.MyHomset = MyHomset # fakes MyHomset being defined in a Python module
            sage: H = MyHomset(X, Y, category=Monoids())
            sage: H
            Set of Morphisms from X to Y in Category of monoids
            sage: H.base()
            Integer Ring
            sage: TestSuite(H).run()
        """
        if base is None:
            base = X.base_ring()
        Homset.__init__(self, X, Y, check=check, category=category, base=base)


def is_Homset(x):
    """
    Return ``True`` if ``x`` is a set of homomorphisms in a category.

    EXAMPLES::

        sage: from sage.categories.homset import is_Homset
        sage: P.<t> = ZZ[]
        sage: f = P.hom([1/2*t])
        sage: is_Homset(f)
        doctest:warning...
        DeprecationWarning: the function is_Homset is deprecated;
        use 'isinstance(..., Homset)' instead
        See https://github.com/sagemath/sage/issues/37922 for details.
        False
        sage: is_Homset(f.category())
        False
        sage: is_Homset(f.parent())
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37922, "the function is_Homset is deprecated; use 'isinstance(..., Homset)' instead")
    return isinstance(x, Homset)


def is_Endset(x):
    """
    Return ``True`` if ``x`` is a set of endomorphisms in a category.

    EXAMPLES::

        sage: from sage.categories.homset import is_Endset
        sage: P.<t> = ZZ[]
        sage: f = P.hom([1/2*t])
        sage: is_Endset(f.parent())
        doctest:warning...
        DeprecationWarning: the function is_Endset is deprecated;
        use 'isinstance(..., Homset) and ....is_endomorphism_set()' instead
        See https://github.com/sagemath/sage/issues/37922 for details.
        False
        sage: g = P.hom([2*t])
        sage: is_Endset(g.parent())
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37922, "the function is_Endset is deprecated; use 'isinstance(..., Homset) and ....is_endomorphism_set()' instead")
    return isinstance(x, Homset) and x.is_endomorphism_set()
