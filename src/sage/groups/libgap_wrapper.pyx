"""
LibGAP-based Groups

This module provides helper class for wrapping GAP groups via
:mod:`~sage.libs.gap.libgap`. See :mod:`~sage.groups.free_group` for an
example how they are used.

The parent class keeps track of the libGAP element object, to use it
in your Python parent you have to derive both from the suitable group
parent and :class:`ParentLibGAP` ::

    sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
    sage: from sage.groups.group import Group
    sage: class FooElement(ElementLibGAP):
    ...       pass
    sage: class FooGroup(Group, ParentLibGAP):
    ...       Element = FooElement
    ...       def __init__(self):
    ...           lg = libgap(libgap.CyclicGroup(3))    # dummy
    ...           ParentLibGAP.__init__(self, lg)
    ...           Group.__init__(self)

Note how we call the constructor of both superclasses to initialize
``Group`` and ``ParentLibGAP`` separately. The parent class implements
its output via LibGAP::

    sage: FooGroup()
    <pc group of size 3 with 1 generators>
    sage: type(FooGroup().gap())
    <type 'sage.libs.gap.element.GapElement'>

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
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.libs.gap.element cimport GapElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.structure.element cimport Element


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
      GAP.

    - ``ambient`` -- A derived class of :class:`ParentLibGAP` or
      ``None`` (default). The ambient class if ``libgap_parent`` has
      been defined as a subgroup.

    EXAMPLES::

        sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
        sage: from sage.groups.group import Group
        sage: class FooElement(ElementLibGAP):
        ...       pass
        sage: class FooGroup(Group, ParentLibGAP):
        ...       Element = FooElement
        ...       def __init__(self):
        ...           lg = libgap(libgap.CyclicGroup(3))    # dummy
        ...           ParentLibGAP.__init__(self, lg)
        ...           Group.__init__(self)
        sage: FooGroup()
        <pc group of size 3 with 1 generators>
    """

    def __init__(self, libgap_parent, ambient=None):
        """
        The Python constructor.

        TESTS::

            sage: G = FreeGroup(3)
            sage: TestSuite(G).run()
        """
        assert isinstance(libgap_parent, GapElement)
        self._libgap = libgap_parent
        self._ambient = ambient

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

        You can access the contaning group with :meth:`ambient`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G.is_subgroup()
            False
        """
        return self._ambient is not None

    def _subgroup_constructor(self, libgap_subgroup):
        """
        Return the class of a subgroup.

        You should override this with a derived class. Its constructor
        must accept the same arguments as :meth:`__init__`.

        OUTPUT:

        A new instance of a group (derived class of
        :class:`ParentLibGAP`).

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

        - ``generators`` -- a list/tuple/iterable of group elements.

        OUTPUT:

        The subgroup generated by ``generators``.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F.subgroup([a^2*b]);  G
            Group([ a^2*b ])
            sage: G.gens()
            (a^2*b,)
        """
        generators = [ g if isinstance(g, GapElement) else g.gap()
                       for g in generators ]
        G = self.gap()
        H = G.Subgroup(generators)
        return self._subgroup_constructor(H)

    def gap(self):
        """
        Returns the gap representation of self

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES::

            sage: G = FreeGroup(3);  G
            Free Group on generators {x0, x1, x2}
            sage: G.gap()
            <free group on the generators [ x0, x1, x2 ]>
            sage: G.gap().parent()
            C library interface to GAP
            sage: type(G.gap())
            <type 'sage.libs.gap.element.GapElement'>

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

    @cached_method
    def _gap_gens(self):
        """
        Return the generators as a LibGAP object

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES:

            sage: G = FreeGroup(2)
            sage: G._gap_gens()
            [ x0, x1 ]
            sage: type(_)
            <type 'sage.libs.gap.element.GapElement_List'>
        """
        return self._libgap.GeneratorsOfGroup()

    @cached_method
    def ngens(self):
        """
        Return the number of generators of self.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: G = FreeGroup(2)
            sage: G.ngens()
            2

        TESTS::

            sage: type(G.ngens())
            <type 'sage.rings.integer.Integer'>
        """
        return self._gap_gens().Length().sage()

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        TESTS::

            sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
            sage: G.<a,b> =FreeGroup()
            sage: ParentLibGAP._repr_(G)
            '<free group on the generators [ a, b ]>'
        """
        return self._libgap._repr_()

    def gen(self, i):
        """
        Return the `i`-th generator of self.

        .. warning::

            Indexing starts at `0` as usual in Sage/Python. Not as in
            GAP, where indexing starts at `1`.

        INPUT:

        - ``i`` -- integer between `0` (inclusive) and :meth:`ngens`
          (exclusive). The index of the generator.

        OUTPUT:

        The `i`-th generator of the group.

        EXAMPLES::

            sage: G = FreeGroup('a, b')
            sage: G.gen(0)
            a
            sage: G.gen(1)
            b
        """
        if not (0 <= i < self.ngens()):
            raise ValueError('i must be in range(ngens)')
        gap = self._gap_gens()[i]
        return self.element_class(self, gap)

    @cached_method
    def gens(self):
        """
        Returns the generators of the group.

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
        return tuple( self.gen(i) for i in range(self.ngens()) )

    generators = gens

    @cached_method
    def one(self):
        """
        Returns the identity element of self

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
        Returns an element of self.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: G._an_element_()
            a*b
        """
        from sage.misc.all import prod
        return prod(self.gens())



cdef class ElementLibGAP(MultiplicativeGroupElement):
    """
    A class for LibGAP-based Sage group elements

    INPUT:

    - ``parent`` -- the Sage parent

    - ``libgap_element`` -- the libgap element that is being wrapped

    EXAMPLES::

        sage: from sage.groups.libgap_wrapper import ElementLibGAP, ParentLibGAP
        sage: from sage.groups.group import Group
        sage: class FooElement(ElementLibGAP):
        ...       pass
        sage: class FooGroup(Group, ParentLibGAP):
        ...       Element = FooElement
        ...       def __init__(self):
        ...           lg = libgap(libgap.CyclicGroup(3))    # dummy
        ...           ParentLibGAP.__init__(self, lg)
        ...           Group.__init__(self)
        sage: FooGroup()
        <pc group of size 3 with 1 generators>
        sage: FooGroup().gens()
        (f1,)
    """

    def __init__(self, parent, libgap_element):
        """
        The Python constructor

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
        Returns a LibGAP representation of the element

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement`

        EXAMPLES::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: x = G([1, 2, -1, -2])
            sage: x
            a*b*a^-1*b^-1
            sage: xg = x.gap()
            sage: xg
            a*b*a^-1*b^-1
            sage: type(xg)
            <type 'sage.libs.gap.element.GapElement'>
        """
        return self._libgap

    _gap_ = gap

    def is_one(self):
        """
        Test whether the group element is the trivial element.

        OUTPUT:

        Boolean.

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

        OUTPUT:

        String.

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
        Return a LaTeX representation

        OUTPUT:

        String. A valid LaTeX math command sequence.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.FreeGroup('a', 'b'))
            sage: g = G.gen(0) * G.gen(1)
            sage: g._latex_()
            "ab%\n"
        """
        try:
            return self.gap().LaTeX()
        except ValueError:
            from sage.misc.latex import latex
            return latex(self._repr_())

    cpdef _mul_(left, right):
        """
        Multiplication of group elements

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
        P = left.parent()
        return P.element_class(P, left.gap() * right.gap())

    cpdef int _cmp_(left, right) except -2:
        """
        This method implements comparison.

        TESTS::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: G_gap = G.gap()
            sage: G_gap == G_gap    # indirect doctest
            True
            sage: cmp(G.gap(), G.gap())   # indirect doctest
            0
            sage: x = G([1, 2, -1, -2])
            sage: y = G([2, 2, 2, 1, -2, -2, -2])
            sage: x == x*y*y^(-1)     # indirect doctest
            True
            sage: cmp(x,y)
            -1
            sage: x < y
            True
        """
        return cmp((<ElementLibGAP>left)._libgap,
                   (<ElementLibGAP>right)._libgap)

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
            sage: x/y == x.__div__(y)
            True
            sage: x/y == y.__div__(x)
            False
        """
        P = left.parent()
        return P.element_class(P, left.gap() / right.gap())

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
        P = self.parent()
        return P.element_class(P, self.gap() ** n)

    def __invert__(self):
        """
        Return the inverse of self.

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
        P = self.parent()
        return P.element_class(P, self.gap().Inverse())

    inverse = __invert__

