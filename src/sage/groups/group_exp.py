r"""
Functor that converts a commutative additive group into a multiplicative group.

AUTHORS:

- Mark Shimozono (2013): initial version
"""
#*****************************************************************************
#       Copyright (C) 2013 <mshimo at math.vt.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.functor import Functor
from sage.categories.groups import Groups
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class GroupExp(Functor):
    r"""
    A functor that converts a commutative additive group into an isomorphic
    multiplicative group.

    More precisely, given a commutative additive group `G`, define the exponential
    of `G` to be the isomorphic group with elements denoted
    `e^g` for every `g \in G` and but with product in multiplicative notation

    .. MATH::

        e^g e^h = e^{g+h} \qquad\text{for all $g,h \in G$.}

    The class :class:`GroupExp` implements the sage functor which sends a commutative
    additive group `G` to its exponential.

    The creation of an instance of the functor :class:`GroupExp` requires no input::

        sage: E = GroupExp(); E
        Functor from Category of commutative additive groups to Category of groups

    The :class:`GroupExp` functor (denoted `E` in the examples) can be applied to two kinds of input.
    The first is a commutative additive group. The output is its exponential.
    This is accomplished by :meth:`_apply_functor`::

        sage: EZ = E(ZZ); EZ
        Multiplicative form of Integer Ring

    Elements of the exponentiated group can be created and manipulated as follows::

        sage: x = EZ(-3); x
        -3
        sage: x.parent()
        Multiplicative form of Integer Ring
        sage: EZ(-1) * EZ(6) == EZ(5)
        True
        sage: EZ(3)^(-1)
        -3
        sage: EZ.one()
        0

    The second kind of input the :class:`GroupExp` functor accepts, is a homomorphism of commutative additive groups.
    The output is the multiplicative form of the homomorphism. This is achieved by :meth:`_apply_functor_to_morphism`::

        sage: L = RootSystem(['A',2]).ambient_space()
        sage: EL = E(L)
        sage: W = L.weyl_group(prefix='s')
        sage: s2 = W.simple_reflection(2)
        sage: def my_action(mu):
        ....:     return s2.action(mu)
        sage: from sage.categories.morphism import SetMorphism
        sage: from sage.categories.homset import Hom
        sage: f = SetMorphism(Hom(L, L, CommutativeAdditiveGroups()), my_action)
        sage: F = E(f); F
        Generic endomorphism of
         Multiplicative form of Ambient space of the Root system of type ['A', 2]
        sage: v = L.an_element(); v
        (2, 2, 3)
        sage: y = F(EL(v)); y
        (2, 3, 2)
        sage: y.parent()
        Multiplicative form of Ambient space of the Root system of type ['A', 2]
    """
    def __init__(self):
        r"""
        Initialize the :class:`GroupExp` functor.

        EXAMPLES::

            sage: F = GroupExp()
            sage: F.domain()
            Category of commutative additive groups
            sage: F.codomain()
            Category of groups
        """
        Functor.__init__(self, CommutativeAdditiveGroups(), Groups())

    def _apply_functor(self, x):
        r"""
        Given a commutative additive group, return the isomorphic
        multiplicative group.

        INPUT:

        - A commutative additive group `x`

        OUTPUT: an isomorphic group whose operation is multiplication rather
        than addition

        In the following example, ``self`` is the functor ``GroupExp()``,
        ``x`` is the additive group ``QQ^2``, and the output group is stored as ``EQ2``.

        EXAMPLES::

            sage: EQ2 = GroupExp()(QQ^2)
            sage: x = EQ2(vector(QQ,(-2,1))); x
            (-2, 1)
            sage: x^(-1)
            (2, -1)
            sage: x*x
            (-4, 2)
            sage: EQ2(vector(QQ,(-1,1)))*EQ2(vector(QQ,(3,4))) == EQ2(vector(QQ,(2,5)))
            True
            sage: EQ2.one()
            (0, 0)
        """
        return GroupExp_Class(x)

    def _apply_functor_to_morphism(self, f):
        r"""
        Given a morphism of commutative additive groups, return the corresponding morphism
        of multiplicative groups.

        INPUT:

        - A homomorphism `f` of commutative additive groups

        OUTPUT: the above homomorphism, but between the corresponding
        multiplicative groups

        In the following example, ``self`` is the functor :class:`GroupExp` and `f`
        is an endomorphism of the additive group of integers.

        EXAMPLES::

            sage: def double(x):
            ....:     return x + x
            sage: from sage.categories.morphism import SetMorphism
            sage: from sage.categories.homset import Hom
            sage: f = SetMorphism(Hom(ZZ, ZZ, CommutativeAdditiveGroups()), double)
            sage: E = GroupExp()
            sage: EZ = E._apply_functor(ZZ)
            sage: F = E._apply_functor_to_morphism(f)
            sage: F.domain() == EZ
            True
            sage: F.codomain() == EZ
            True
            sage: F(EZ(3)) == EZ(3) * EZ(3)
            True
        """
        new_domain = self._apply_functor(f.domain())
        new_codomain = self._apply_functor(f.codomain())
        new_f = lambda a: new_codomain(f(a.value))
        return SetMorphism(Hom(new_domain, new_codomain, Groups()), new_f)


class GroupExpElement(ElementWrapper, MultiplicativeGroupElement):
    r"""
    An element in the exponential of a commutative additive group.

    INPUT:

    - ``self`` -- the exponentiated group element being created
    - ``parent`` -- the exponential group (parent of ``self``)
    - ``x`` -- the commutative additive group element being wrapped to form ``self``

    EXAMPLES::

        sage: from sage.groups.group_exp import GroupExpElement
        sage: G = QQ^2
        sage: EG = GroupExp()(G)
        sage: z = GroupExpElement(EG, vector(QQ, (1, -3))); z
        (1, -3)
        sage: z.parent()
        Multiplicative form of Vector space of dimension 2 over Rational Field
        sage: EG(vector(QQ, (1, -3))) == z
        True
    """
    def __init__(self, parent, x):
        r"""
        EXAMPLES::

            sage: G = QQ^2
            sage: EG = GroupExp()(G)
            sage: x = EG.an_element(); x
            (1, 0)
            sage: TestSuite(x).run(skip="_test_category")

        See the documentation of :meth:`sage.structure.element_wrapper.ElementWrapper.__init__`
        for the reason behind skipping the category test.
        """
        if x not in parent._G:
            raise ValueError("%s is not an element of %s" % (x, parent._G))
        ElementWrapper.__init__(self, parent, x)

    def __invert__(self):
        r"""
        Invert the element ``self``.

        EXAMPLES::

            sage: EZ = GroupExp()(ZZ)
            sage: EZ(-3).inverse()  # indirect doctest
            3
        """
        return GroupExpElement(self.parent(), -self.value)

    def __mul__(self, x):
        r"""
        Multiply ``self`` by `x`.

        EXAMPLES::

            sage: G = GroupExp()(ZZ)
            sage: x = G(2)
            sage: x.__mul__(G(3))
            5
            sage: G.product(G(2), G(3))
            5
        """
        return GroupExpElement(self.parent(), self.value + x.value)


class GroupExp_Class(UniqueRepresentation, Parent):
    r"""
    The multiplicative form of a commutative additive group.

    INPUT:

    - ``G`` -- a commutative additive group

    OUTPUT: the multiplicative form of `G`

    EXAMPLES::

        sage: GroupExp()(QQ)
        Multiplicative form of Rational Field
    """
    def __init__(self, G):
        r"""

        EXAMPLES::

            sage: EG = GroupExp()(QQ^2)
            sage: TestSuite(EG).run(skip="_test_elements")
        """
        if G not in CommutativeAdditiveGroups():
            raise TypeError("%s must be a commutative additive group" % G)
        self._G = G
        Parent.__init__(self, category=Groups())

    def _repr_(self):
        r"""
        Return a string describing the multiplicative form of a commutative additive group.

        EXAMPLES::

            sage: GroupExp()(ZZ) # indirect doctest
            Multiplicative form of Integer Ring
        """
        return "Multiplicative form of %s" % self._G

    def _element_constructor_(self, x):
        r"""
        Construct the multiplicative group element, which wraps the additive
        group element `x`.

        EXAMPLES::

            sage: G = GroupExp()(ZZ)
            sage: G(4)  # indirect doctest
            4
        """
        return GroupExpElement(self, x)

    def one(self):
        r"""
        Return the identity element of the multiplicative group.

        EXAMPLES::

            sage: G = GroupExp()(ZZ^2)
            sage: G.one()
            (0, 0)
            sage: x = G.an_element(); x
            (1, 0)
            sage: x == x * G.one()
            True
        """
        return GroupExpElement(self, self._G.zero())

    def an_element(self):
        r"""
        Return an element of the multiplicative group.

        EXAMPLES::

            sage: L = RootSystem(['A',2]).weight_lattice()
            sage: EL = GroupExp()(L)
            sage: x = EL.an_element(); x
            2*Lambda[1] + 2*Lambda[2]
            sage: x.parent()
            Multiplicative form of Weight lattice of the Root system of type ['A', 2]
        """
        return GroupExpElement(self, self._G.an_element())

    def product(self, x, y):
        r"""
        Return the product of `x` and `y` in the multiplicative group.

        EXAMPLES::

            sage: G = GroupExp()(ZZ)
            sage: G.product(G(2),G(7))
            9
            sage: x = G(2)
            sage: x.__mul__(G(7))
            9
        """
        return GroupExpElement(self, x.value + y.value)

    def group_generators(self):
        r"""
        Return generators of ``self``.

        EXAMPLES::

            sage: GroupExp()(ZZ).group_generators()
            (1,)
        """
        if hasattr(self._G, 'gens'):
            additive_generators = self._G.gens()
        else:
            raise AttributeError("Additive group has no method 'gens'")
        return tuple([self(x) for x in additive_generators])


GroupExp_Class.Element = GroupExpElement
