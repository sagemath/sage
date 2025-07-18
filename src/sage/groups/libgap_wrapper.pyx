"""
LibGAP-based Groups

This module provides helper class for wrapping GAP groups via
:mod:`~sage.libs.gap.libgap`. See :mod:`~sage.groups.free_group` for an
example how they are used.

The parent class keeps track of the GAP element object, to use it
in your Python parent you have to derive both from the suitable group
parent and :class:`ParentLibGAP` ::

    sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
    sage: from sage.groups.group import Group
    sage: class FooElement(ElementLibGAP):
    ....:     pass
    sage: class FooGroup(Group, ParentLibGAP):
    ....:     Element = FooElement
    ....:     def __init__(self):
    ....:         lg = libgap(libgap.CyclicGroup(3))    # dummy
    ....:         ParentLibGAP.__init__(self, lg)
    ....:         Group.__init__(self)

Note how we call the constructor of both superclasses to initialize
``Group`` and ``ParentLibGAP`` separately. The parent class implements
its output via LibGAP::

    sage: FooGroup()
    <pc group of size 3 with 1 generator>
    sage: type(FooGroup().gap())
    <class 'sage.libs.gap.element.GapElement'>

The element class is a subclass of
:class:`~sage.structure.element.MultiplicativeGroupElement`. To use
it, you just inherit from :class:`ElementLibGAP` ::

    sage: element = FooGroup().an_element()
    sage: element
    f1

The element class implements group operations and printing via LibGAP::

    sage: element._repr_()
    'f1'
    sage: element * element
    f1^2

AUTHORS:

- Volker Braun
"""

##############################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.libs.gap.libgap import libgap
from sage.libs.gap.element cimport GapElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.structure.richcmp cimport richcmp


class ParentLibGAP(SageObject):
    """
    A class for parents to keep track of the GAP parent.

    This is not a complete group in Sage, this class is only a base
    class that you can use to implement your own groups with
    LibGAP. See :mod:`~sage.groups.libgap_group` for a minimal example
    of a group that is actually usable.

    Your implementation definitely needs to supply

    * ``__reduce__()``: serialize the LibGAP group. Since GAP does not
      support Python pickles natively, you need to figure out yourself
      how you can recreate the group from a pickle.

    INPUT:

    - ``libgap_parent`` -- the libgap element that is the parent in
      GAP

    - ``ambient`` -- a derived class of :class:`ParentLibGAP` or
      ``None`` (default); the ambient class if ``libgap_parent`` has
      been defined as a subgroup

    EXAMPLES::

        sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
        sage: from sage.groups.group import Group
        sage: class FooElement(ElementLibGAP):
        ....:     pass
        sage: class FooGroup(Group, ParentLibGAP):
        ....:     Element = FooElement
        ....:     def __init__(self):
        ....:         lg = libgap(libgap.CyclicGroup(3))    # dummy
        ....:         ParentLibGAP.__init__(self, lg)
        ....:         Group.__init__(self)
        sage: FooGroup()
        <pc group of size 3 with 1 generator>
    """

    def __init__(self, libgap_parent, ambient=None):
        """
        The Python constructor.

        TESTS::

            sage: G = FreeGroup(3)
            sage: TestSuite(G).run()

        We check that :issue:`19270` is fixed::

            sage: G = GL(2,5)
            sage: g = G( matrix([[1,0],[0,4]]))
            sage: H = G.subgroup([g])
            sage: g in H
            True
        """
        assert isinstance(libgap_parent, GapElement)
        self._libgap = libgap_parent
        self._ambient = ambient
        if ambient is not None:
            phi = self.hom(lambda x: ambient(x.gap()), codomain=ambient)  # the .gap() avoids an infinite recursion
            ambient.register_coercion(phi)

    def ambient(self):
        """
        Return the ambient group of a subgroup.

        OUTPUT:

        A group containing ``self``. If ``self`` has not been defined
        as a subgroup, we just return ``self``.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G.ambient() is G
            True
        """
        if self._ambient is None:
            return self
        else:
            return self._ambient

    def is_subgroup(self):
        """
        Return whether the group was defined as a subgroup of a bigger
        group.

        You can access the containing group with :meth:`ambient`.

        OUTPUT: boolean

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G.is_subgroup()
            False
        """
        return self._ambient is not None

    def _Hom_(self, G, category=None, check=True):
        r"""
        Return the set of group homomorphisms from ``self`` to ``G``.

        INPUT:

        - ``G`` -- group; the codomain
        - ``cat`` -- category

        OUTPUT: the set of homomorphisms from ``self`` to ``G``

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: F.Hom(F)
            Set of Morphisms from Free Group on generators {a, b}
             to Free Group on generators {a, b} in Category of infinite groups
        """
        from sage.groups.libgap_morphism import GroupHomset_libgap
        return GroupHomset_libgap(self, G, category=category, check=check)

    def _subgroup_constructor(self, libgap_subgroup):
        """
        Return the class of a subgroup.

        You should override this with a derived class. Its constructor
        must accept the same arguments as :meth:`__init__`.

        OUTPUT: a new instance of a group (derived class of :class:`ParentLibGAP`)

        TESTS::

            sage: F.<a,b> = FreeGroup()
            sage: G = F.subgroup([a^2*b]);  G
            Group([ a^2*b ])
            sage: F._subgroup_constructor(G.gap())._repr_()
            'Group([ a^2*b ])'
        """
        from sage.groups.libgap_group import GroupLibGAP
        return GroupLibGAP(libgap_subgroup, ambient=self)

    def subgroup(self, generators):
        """
        Return the subgroup generated.

        INPUT:

        - ``generators`` -- list/tuple/iterable of group elements

        OUTPUT: the subgroup generated by ``generators``

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F.subgroup([a^2*b]);  G
            Group([ a^2*b ])
            sage: G.gens()
            (a^2*b,)

        We check that coercions between the subgroup and its ambient group work::

            sage: F.0 * G.0
            a^3*b

        Checking that :issue:`19270` is fixed::

            sage: gens = [w.matrix() for w in WeylGroup(['B', 3])]
            sage: G = MatrixGroup(gens)
            sage: import itertools
            sage: diagonals = itertools.product((1,-1), repeat=3)
            sage: subgroup_gens = [diagonal_matrix(L) for L in diagonals]
            sage: G.subgroup(subgroup_gens)
            Subgroup with 8 generators of Matrix group over Rational Field with 48 generators

        TESTS:

        Check that :issue:`19010` is fixed::

            sage: G = WeylGroup(['B',3])
            sage: H = G.subgroup([G[14], G[17]])
            sage: all(g*h in G and h*g in G for g in G for h in H)
            True
        """
        generators = [ g if isinstance(g, GapElement) else self(g).gap()
                       for g in generators ]
        G = self.gap()
        H = G.Subgroup(generators)
        return self._subgroup_constructor(H)

    def gap(self):
        """
        Return the GAP representation of ``self``.

        OUTPUT: a :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES::

            sage: G = FreeGroup(3);  G
            Free Group on generators {x0, x1, x2}
            sage: G.gap()
            <free group on the generators [ x0, x1, x2 ]>
            sage: G.gap().parent()
            C library interface to GAP
            sage: type(G.gap())
            <class 'sage.libs.gap.element.GapElement'>

        This can be useful, for example, to call GAP functions that
        are not wrapped in Sage::

            sage: G = FreeGroup(3)
            sage: H = G.gap()
            sage: H.DirectProduct(H)
            <fp group on the generators [ f1, f2, f3, f4, f5, f6 ]>
            sage: H.DirectProduct(H).RelatorsOfFpGroup()
            [ f1^-1*f4^-1*f1*f4, f1^-1*f5^-1*f1*f5, f1^-1*f6^-1*f1*f6, f2^-1*f4^-1*f2*f4,
              f2^-1*f5^-1*f2*f5, f2^-1*f6^-1*f2*f6, f3^-1*f4^-1*f3*f4, f3^-1*f5^-1*f3*f5,
              f3^-1*f6^-1*f3*f6 ]

        We can also convert directly to libgap::

            sage: libgap(GL(2, ZZ))
            GL(2,Integers)
        """
        return self._libgap

    _libgap_ = _gap_ = gap

    def ngens(self):
        """
        Return the number of generators of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: G = FreeGroup(2)
            sage: G.ngens()
            2

        TESTS::

            sage: type(G.ngens())
            <class 'sage.rings.integer.Integer'>
        """
        return Integer(len(self.gens()))

    def _repr_(self):
        """
        Return a string representation.

        TESTS::

            sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
            sage: G.<a,b> =FreeGroup()
            sage: ParentLibGAP._repr_(G)
            '<free group on the generators [ a, b ]>'
        """
        return self._libgap._repr_()

    def minimal_normal_subgroups(self):
        """
        Return the nontrivial minimal normal subgroups of ``self``.

        EXAMPLES::

            sage: SL(2,GF(49)).minimal_normal_subgroups()
            [Subgroup with 1 generators (
             [6 0]
             [0 6]
             ) of Special Linear Group of degree 2 over Finite Field in z2 of size 7^2]
        """
        return [self._subgroup_constructor(gap_subgroup)
                for gap_subgroup in self._libgap.MinimalNormalSubgroups()]

    def maximal_normal_subgroups(self):
        """
        Return the maximal proper normal subgroups of ``self``.

        This raises an error if `G/[G, G]` is infinite, yielding infinitely
        many maximal normal subgroups.

        EXAMPLES::

            sage: SL(2,GF(49)).minimal_normal_subgroups()
            [Subgroup with 1 generators (
             [6 0]
             [0 6]
             ) of Special Linear Group of degree 2 over Finite Field in z2 of size 7^2]
        """
        return [self._subgroup_constructor(gap_subgroup)
                for gap_subgroup in self._libgap.MaximalNormalSubgroups()]

    @cached_method
    def gens(self) -> tuple:
        """
        Return the generators of the group.

        EXAMPLES::

            sage: G = FreeGroup(2)
            sage: G.gens()
            (x0, x1)
            sage: H = FreeGroup('a, b, c')
            sage: H.gens()
            (a, b, c)

        :meth:`generators` is an alias for :meth:`gens` ::

            sage: G = FreeGroup('a, b')
            sage: G.generators()
            (a, b)
            sage: H = FreeGroup(3, 'x')
            sage: H.generators()
            (x0, x1, x2)
        """
        return tuple(self.element_class(self, g) for g in self._libgap.GeneratorsOfGroup())

    generators = gens

    def gen(self, i):
        """
        Return the `i`-th generator of ``self``.

        .. WARNING::

            Indexing starts at `0` as usual in Sage/Python. Not as in
            GAP, where indexing starts at `1`.

        INPUT:

        - ``i`` -- integer between `0` (inclusive) and :meth:`ngens`
          (exclusive); the index of the generator

        OUTPUT: the `i`-th generator of the group

        EXAMPLES::

            sage: G = FreeGroup('a, b')
            sage: G.gen(0)
            a
            sage: G.gen(1)
            b
        """
        if not (0 <= i < self.ngens()):
            raise ValueError('i must be in range(ngens)')
        return self.gens()[i]

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G.one()
            1
            sage: G.one() == G([])
            True
            sage: G.one().Tietze()
            ()
        """
        return self.element_class(self, self.gap().Identity())

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: G._an_element_()
            a*b
        """
        from sage.misc.misc_c import prod
        gens = self.gens()
        if gens:
            return prod(gens)
        else:
            return self.one()


cdef class ElementLibGAP(MultiplicativeGroupElement):
    """
    A class for LibGAP-based Sage group elements.

    INPUT:

    - ``parent`` -- the Sage parent

    - ``libgap_element`` -- the libgap element that is being wrapped

    EXAMPLES::

        sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
        sage: from sage.groups.group import Group
        sage: class FooElement(ElementLibGAP):
        ....:     pass
        sage: class FooGroup(Group, ParentLibGAP):
        ....:     Element = FooElement
        ....:     def __init__(self):
        ....:         lg = libgap(libgap.CyclicGroup(3))    # dummy
        ....:         ParentLibGAP.__init__(self, lg)
        ....:         Group.__init__(self)
        sage: FooGroup()
        <pc group of size 3 with 1 generator>
        sage: FooGroup().gens()
        (f1,)
    """

    def __init__(self, parent, libgap_element):
        """
        The Python constructor.

        TESTS::

            sage: G = FreeGroup(2)
            sage: g = G.an_element()
            sage: TestSuite(g).run()
        """
        MultiplicativeGroupElement.__init__(self, parent)
        assert isinstance(parent, ParentLibGAP)
        if isinstance(libgap_element, GapElement):
            self._libgap = libgap_element
        else:
            if libgap_element == 1:
                self._libgap = self.parent().gap().Identity()
            else:
                raise TypeError('need a libgap group element or "1" in constructor')

    cpdef GapElement gap(self):
        """
        Return a LibGAP representation of the element.

        OUTPUT: a :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: x
            a*b*a^-1*b^-1
            sage: xg = x.gap()
            sage: xg
            a*b*a^-1*b^-1
            sage: type(xg)
            <class 'sage.libs.gap.element.GapElement'>

        TESTS::

            sage: libgap(FreeGroup('a, b').an_element())
            a*b
            sage: type(libgap(FreeGroup('a, b').an_element()))
            <class 'sage.libs.gap.element.GapElement'>
        """
        return self._libgap

    _libgap_ = _gap_ = gap

    def _test_libgap_conversion(self, **options):
        r"""
        TESTS::

            sage: FreeGroup(2).an_element()._test_libgap_conversion()
        """
        tester = self._tester(**options)
        tester.assertIs(libgap(self), self.gap())

    def _test_libgap_reconstruction(self, **options):
        r"""
        TESTS::

            sage: FreeGroup(2).an_element()._test_libgap_reconstruction()
        """
        tester = self._tester(**options)
        P = self.parent()
        tester.assertEqual(self, P.element_class(P, libgap(self)))

    def is_one(self):
        """
        Test whether the group element is the trivial element.

        OUTPUT: boolean

        EXAMPLES::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: x.is_one()
            False
            sage: (x * ~x).is_one()
            True
        """
        return self == self.parent().one()

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: a._repr_()
            'a'
            sage: type(a)
            <class 'sage.groups.free_group.FreeGroup_class_with_category.element_class'>

            sage: x = G([1, 2, -1, -2])
            sage: x._repr_()
            'a*b*a^-1*b^-1'
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: y._repr_()
            'b^3*a*b^-3'

            sage: G.one()
            1
        """
        if self.is_one():
            return '1'
        else:
            return self._libgap._repr_()

    def _latex_(self):
        r"""
        Return a LaTeX representation.

        OUTPUT: string; a valid LaTeX math command sequence

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.FreeGroup('a', 'b'))
            sage: g = G.gen(0) * G.gen(1)
            sage: latex(g)
            \text{\texttt{a*b}}
        """
        from sage.misc.latex import latex
        return latex(self._repr_())

    cpdef _mul_(left, right):
        """
        Multiplication of group elements.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x*y    # indirect doctest
            a*b*a^-1*b^2*a*b^-3
            sage: y*x    # indirect doctest
            b^3*a*b^-3*a*b*a^-1*b^-1
            sage: x*y == x._mul_(y)
            True
            sage: y*x == y._mul_(x)
            True
        """
        P = (<ElementLibGAP> left)._parent
        return P.element_class(P, (<ElementLibGAP> left)._libgap *
                                  (<ElementLibGAP> right)._libgap)

    cpdef _richcmp_(left, right, int op):
        """
        This method implements comparison.

        TESTS::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: G_gap = G.gap()
            sage: G_gap == G_gap    # indirect doctest
            True
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x == x*y*y^(-1)     # indirect doctest
            True
            sage: x < y
            True
        """
        return richcmp((<ElementLibGAP>left)._libgap,
                       (<ElementLibGAP>right)._libgap, op)

    cpdef _div_(left, right):
        """
        Division of group elements.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x/y # indirect doctest
            a*b*a^-1*b^2*a^-1*b^-3
            sage: y/x # indirect doctest
            b^3*a*b^-2*a*b^-1*a^-1
            sage: x/y == x.__truediv__(y)
            True
            sage: x/y == y.__truediv__(x)
            False
        """
        P = (<ElementLibGAP> left)._parent
        return P.element_class(P, (<ElementLibGAP> left)._libgap /
                                  (<ElementLibGAP> right)._libgap)

    def __pow__(self, n, dummy):
        """
        Implement exponentiation.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: y^(2) # indirect doctest
            b^3*a^2*b^-3
            sage: x^(-3) # indirect doctest
            (b*a*b^-1*a^-1)^3
            sage: y^3  ==  y.__pow__(3)
            True
        """
        if n not in IntegerRing():
            raise TypeError("exponent must be an integer")
        P = (<ElementLibGAP> self)._parent
        return P.element_class(P, (<ElementLibGAP> self)._libgap ** n)

    def __invert__(self):
        """
        Return the inverse of ``self``.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x.__invert__()
            b*a*b^-1*a^-1
            sage: y.__invert__()
            b^3*a^-1*b^-3
            sage: ~x
            b*a*b^-1*a^-1
            sage: x.inverse()
            b*a*b^-1*a^-1
        """
        P = (<ElementLibGAP> self)._parent
        return P.element_class(P, (<ElementLibGAP> self)._libgap.Inverse())

    inverse = __invert__

    def order(self):
        r"""
        Return the multiplicative order.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.GL(2, 3))
            sage: a,b = G.gens()
            sage: print(a.order())
            2
            sage: print(a.multiplicative_order())
            2

            sage: z = Mod(0, 3)
            sage: o = Mod(1, 3)
            sage: G(libgap([[o,o],[z,o]])).order()
            3
        """
        return self._libgap.Order().sage()

    multiplicative_order = order

    def is_conjugate(self, other):
        r"""
        Return whether the elements ``self`` and ``other`` are conjugate.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.GL(2, 3))
            sage: a,b = G.gens()
            sage: a.is_conjugate(b)
            False
            sage: a.is_conjugate((a*b^2) * a * ~(a*b^2))
            True
        """
        return libgap.IsConjugate(self.parent(), self, other).sage()

    def normalizer(self):
        r"""
        Return the normalizer of the cyclic group generated by this element.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.GL(3,3))
            sage: a,b = G.gens()
            sage: H = a.normalizer()
            sage: H
            <group of 3x3 matrices over GF(3)>
            sage: H.cardinality()
            96
            sage: all(g*a == a*g for g in H)
            True
        """
        from sage.groups.libgap_group import GroupLibGAP
        P = self.parent()
        return GroupLibGAP(libgap.Normalizer(P, self), ambient=P)

    def nth_roots(self, n):
        r"""
        Return the set of ``n``-th roots of this group element.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.GL(3,3))
            sage: a,b = G.gens()
            sage: g = a*b**2*a*~b
            sage: r = g.nth_roots(4)
            sage: r
            [[ [ Z(3), Z(3), Z(3)^0 ], [ Z(3)^0, Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3), 0*Z(3) ] ],
             [ [ Z(3)^0, Z(3)^0, Z(3) ], [ Z(3), Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0, 0*Z(3) ] ]]
            sage: r[0]**4 == r[1]**4 == g
            True
        """
        P = self.parent()
        return [P.element_class(P, g) for g in libgap.NthRootsInGroup(P.gap(), self._libgap, n)]
