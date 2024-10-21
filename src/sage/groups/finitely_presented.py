# sage_setup: distribution = sagemath-groups
"""
Finitely Presented Groups

Finitely presented groups are constructed as quotients of
:mod:`~sage.groups.free_group`::

    sage: F.<a,b,c> = FreeGroup()
    sage: G = F / [a^2, b^2, c^2, a*b*c*a*b*c]
    sage: G
    Finitely presented group < a, b, c | a^2, b^2, c^2, (a*b*c)^2 >

One can create their elements by multiplying the generators or by
specifying a Tietze list (see
:meth:`~sage.groups.finitely_presented.FinitelyPresentedGroupElement.Tietze`)
as in the case of free groups::

    sage: G.gen(0) * G.gen(1)
    a*b
    sage: G([1,2,-1])
    a*b*a^-1
    sage: a.parent()
    Free Group on generators {a, b, c}
    sage: G.inject_variables()
    Defining a, b, c
    sage: a.parent()
    Finitely presented group < a, b, c | a^2, b^2, c^2, (a*b*c)^2 >

Notice that, even if they are represented in the same way, the
elements of a finitely presented group and the elements of the
corresponding free group are not the same thing.  However, they can be
converted from one parent to the other::

    sage: F.<a,b,c> = FreeGroup()
    sage: G = F / [a^2,b^2,c^2,a*b*c*a*b*c]
    sage: F([1])
    a
    sage: G([1])
    a
    sage: F([1]) is G([1])
    False
    sage: F([1]) == G([1])
    False
    sage: G(a*b/c)
    a*b*c^-1
    sage: F(G(a*b/c))
    a*b*c^-1

Finitely presented groups are implemented via GAP. You can use the
:meth:`~sage.groups.libgap_wrapper.ParentLibGAP.gap` method to access
the underlying LibGAP object::

    sage: G = FreeGroup(2)
    sage: G.inject_variables()
    Defining x0, x1
    sage: H = G / (x0^2, (x0*x1)^2, x1^2)
    sage: H.gap()
    <fp group on the generators [ x0, x1 ]>

This can be useful, for example, to use GAP functions that are not yet
wrapped in Sage::

    sage: H.gap().LowerCentralSeries()
    [ Group(<fp, no generators known>), Group(<fp, no generators known>) ]

The same holds for the group elements::

    sage: G = FreeGroup(2)
    sage: H = G / (G([1, 1]), G([2, 2, 2]), G([1, 2, -1, -2]));  H
    Finitely presented group < x0, x1 | x0^2, x1^3, x0*x1*x0^-1*x1^-1 >
    sage: a = H([1])
    sage: a
    x0
    sage: a.gap()
    x0
    sage: a.gap().Order()
    2
    sage: type(_)    # note that the above output is not a Sage integer
    <class 'sage.libs.gap.element.GapElement_Integer'>

You can use call syntax to replace the generators with a set of
arbitrary ring elements. For example, take the free abelian group
obtained by modding out the commutator subgroup of the free group::

    sage: G = FreeGroup(2)
    sage: G_ab = G / [G([1, 2, -1, -2])];  G_ab
    Finitely presented group < x0, x1 | x0*x1*x0^-1*x1^-1 >
    sage: a,b = G_ab.gens()
    sage: g =  a * b
    sage: M1 = matrix([[1,0],[0,2]])
    sage: M2 = matrix([[0,1],[1,0]])
    sage: g(3, 5)
    15
    sage: g(M1, M1)
    [1 0]
    [0 4]
    sage: M1*M2 == M2*M1   # matrices do not commute
    False
    sage: g(M1, M2)
    Traceback (most recent call last):
    ...
    ValueError: the values do not satisfy all relations of the group

.. WARNING::

    Some methods are not guaranteed to finish since the word problem
    for finitely presented groups is, in general, undecidable. In
    those cases the process may run until the available memory is
    exhausted.

REFERENCES:

- :wikipedia:`Presentation_of_a_group`

- :wikipedia:`Word_problem_for_groups`

AUTHOR:

- Miguel Angel Marco Buzunariz
"""

# ****************************************************************************
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import GCD as gcd
from sage.categories.morphism import SetMorphism
from sage.groups.free_group import FreeGroup
from sage.groups.free_group import FreeGroupElement
from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.set import Set
from sage.structure.richcmp import richcmp, richcmp_method
from sage.structure.unique_representation import CachedRepresentation


class GroupMorphismWithGensImages(SetMorphism):
    r"""
    Class used for morphisms from finitely presented groups to
    other groups. It just adds the images of the generators at the
    end of the representation.

    EXAMPLES::

        sage: F = FreeGroup(3)
        sage: G = F / [F([1, 2, 3, 1, 2, 3]), F([1, 1, 1])]
        sage: H = AlternatingGroup(3)
        sage: HS = G.Hom(H)
        sage: from sage.groups.finitely_presented import GroupMorphismWithGensImages
        sage: GroupMorphismWithGensImages(HS, lambda a: H.one())
        Generic morphism:
        From: Finitely presented group < x0, x1, x2 | (x0*x1*x2)^2, x0^3 >
        To:   Alternating group of order 3!/2 as a permutation group
        Defn: x0 |--> ()
              x1 |--> ()
              x2 |--> ()
    """
    def _repr_defn(self):
        r"""
        Return the part of the representation that includes the images of the generators.

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [F([1,2,3,1,2,3]),F([1,1,1])]
            sage: H = AlternatingGroup(3)
            sage: HS = G.Hom(H)
            sage: from sage.groups.finitely_presented import GroupMorphismWithGensImages
            sage: f = GroupMorphismWithGensImages(HS, lambda a: H.one())
            sage: f._repr_defn()
            'x0 |--> ()\nx1 |--> ()\nx2 |--> ()'
        """
        return '\n'.join(f'{i} |--> {self(i)}' for i in self.domain().gens())


class FinitelyPresentedGroupElement(FreeGroupElement):
    """
    A wrapper of GAP's Finitely Presented Group elements.

    The elements are created by passing the Tietze list that determines them.

    EXAMPLES::

        sage: G = FreeGroup('a, b')
        sage: H = G / [G([1]), G([2, 2, 2])]
        sage: H([1, 2, 1, -1])
        a*b
        sage: H([1, 2, 1, -2])
        a*b*a*b^-1
        sage: x = H([1, 2, -1, -2])
        sage: x
        a*b*a^-1*b^-1
        sage: y = H([2, 2, 2, 1, -2, -2, -2])
        sage: y
        b^3*a*b^-3
        sage: x*y
        a*b*a^-1*b^2*a*b^-3
        sage: x^(-1)
        b*a*b^-1*a^-1
    """

    def __init__(self, parent, x, check=True):
        """
        The Python constructor.

        See :class:`FinitelyPresentedGroupElement` for details.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: H = G / [G([1]), G([2, 2, 2])]
            sage: H([1, 2, 1, -1])
            a*b

            sage: TestSuite(G).run()
            sage: TestSuite(H).run()
            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2, 2, 2]))
            sage: x = H([1, 2, -1, -2])
            sage: TestSuite(x).run()
            sage: TestSuite(G.one()).run()
        """
        if not isinstance(x, GapElement):
            F = parent.free_group()
            free_element = F(x)
            fp_family = parent.gap().Identity().FamilyObj()
            x = libgap.ElementOfFpGroup(fp_family, free_element.gap())
        ElementLibGAP.__init__(self, parent, x)

    def __reduce__(self):
        """
        Used in pickling.

        TESTS::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a*b, a^2]
            sage: G.inject_variables()
            Defining a, b
            sage: a.__reduce__()
            (Finitely presented group < a, b | a*b, a^2 >, ((1,),))
            sage: (a*b*a^-1).__reduce__()
            (Finitely presented group < a, b | a*b, a^2 >, ((1, 2, -1),))

            sage: F.<a,b,c> = FreeGroup('a, b, c')
            sage: G = F.quotient([a*b*c/(b*c*a), a*b*c/(c*a*b)])
            sage: G.inject_variables()
            Defining a, b, c
            sage: x = a*b*c
            sage: x.__reduce__()
            (Finitely presented group < a, b, c | a*b*c*a^-1*c^-1*b^-1, a*b*c*b^-1*a^-1*c^-1 >,
             ((1, 2, 3),))
        """
        return (self.parent(), (self.Tietze(),))

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / [a^2, b^3]
            sage: H.gen(0)
            a
            sage: H.gen(0)._repr_()
            'a'
            sage: H.one()
            1
        """
        # computing that an element is actually one can be very expensive
        if self.Tietze() == ():
            return '1'
        else:
            return self.gap()._repr_()

    @cached_method
    def Tietze(self):
        """
        Return the Tietze list of the element.

        The Tietze list of a word is a list of integers that represent
        the letters in the word.  A positive integer `i` represents
        the letter corresponding to the `i`-th generator of the group.
        Negative integers represent the inverses of generators.

        OUTPUT: tuple of integers

        EXAMPLES::

            sage: G = FreeGroup('a, b')
            sage: H = G / (G([1]), G([2, 2, 2]))
            sage: H.inject_variables()
            Defining a, b
            sage: a.Tietze()
            (1,)
            sage: x = a^2*b^(-3)*a^(-2)
            sage: x.Tietze()
            (1, 1, -2, -2, -2, -1, -1)
        """
        tl = self.gap().UnderlyingElement().TietzeWordAbstractWord()
        return tuple(tl.sage())

    def __call__(self, *values, **kwds):
        """
        Replace the generators of the free group with ``values``.

        INPUT:

        - ``*values`` -- list/tuple/iterable of the same length as
          the number of generators

        - ``check=True`` -- boolean keyword (default: ``True``); whether to
          verify that ``values`` satisfy the relations in the finitely
          presented group

        OUTPUT: the product of ``values`` in the order and with exponents
        specified by ``self``

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / [a/b];  H
            Finitely presented group < a, b | a*b^-1 >
            sage: H.simplified()
            Finitely presented group < a |  >

        The generator `b` can be eliminated using the relation `a=b`. Any
        values that you plug into a word must satisfy this relation::

            sage: A, B = H.gens()
            sage: w = A^2 * B
            sage: w(2,2)
            8
            sage: w(3,3)
            27
            sage: w(1,2)
            Traceback (most recent call last):
            ...
            ValueError: the values do not satisfy all relations of the group
            sage: w(1, 2, check=False)    # result depends on presentation of the group element
            2
        """
        values = list(values)
        if kwds.get('check', True):
            for rel in self.parent().relations():
                rel = rel(values)
                if rel != 1:
                    raise ValueError('the values do not satisfy all relations of the group')
        return super().__call__(values)


class RewritingSystem():
    """
    A class that wraps GAP's rewriting systems.

    A rewriting system is a set of rules that allow to transform
    one word in the group to an equivalent one.

    If the rewriting system is confluent, then the transformed
    word is a unique reduced form of the element of the group.

    .. WARNING::

        Note that the process of making a rewriting system confluent
        might not end.

    INPUT:

    - ``G`` -- a group

    REFERENCES:

    - :wikipedia:`Knuth-Bendix_completion_algorithm`

    EXAMPLES::

        sage: F.<a,b> = FreeGroup()
        sage: G = F / [a*b/a/b]
        sage: k = G.rewriting_system()
        sage: k
        Rewriting system of Finitely presented group < a, b | a*b*a^-1*b^-1 >
        with rules:
            a*b*a^-1*b^-1    --->    1

        sage: k.reduce(a*b*a*b)
        (a*b)^2
        sage: k.make_confluent()
        sage: k
        Rewriting system of Finitely presented group < a, b | a*b*a^-1*b^-1 >
        with rules:
            b^-1*a^-1    --->    a^-1*b^-1
            b^-1*a    --->    a*b^-1
            b*a^-1    --->    a^-1*b
            b*a    --->    a*b

        sage: k.reduce(a*b*a*b)
        a^2*b^2

    .. TODO::

        - Include support for different orderings (currently only shortlex
          is used).

        - Include the GAP package kbmag for more functionalities, including
          automatic structures and faster compiled functions.

    AUTHORS:

    - Miguel Angel Marco Buzunariz (2013-12-16)
    """
    def __init__(self, G):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F.<a,b,c> = FreeGroup()
            sage: G = F / [a^2, b^3, c^5]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b, c | a^2, b^3, c^5 >
            with rules:
                a^2    --->    1
                b^3    --->    1
                c^5    --->    1
        """
        self._free_group = G.free_group()
        self._fp_group = G
        self._fp_group_gap = G.gap()
        self._monoid_isomorphism = self._fp_group_gap.IsomorphismFpMonoid()
        self._monoid = self._monoid_isomorphism.Image()
        self._gap = self._monoid.KnuthBendixRewritingSystem()

    def __repr__(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: F.<a> = FreeGroup()
            sage: G = F / [a^2]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a | a^2 >
            with rules:
                a^2    --->    1
        """
        ret = "Rewriting system of {}\nwith rules:".format(self._fp_group)
        for i in sorted(self.rules().items()):  # Make sure they are sorted to the repr is unique
            ret += "\n    {}    --->    {}".format(i[0], i[1])
        return ret

    def free_group(self):
        """
        The free group after which the rewriting system is defined.

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [ [1,2,3], [-1,-2,-3] ]
            sage: k = G.rewriting_system()
            sage: k.free_group()
            Free Group on generators {x0, x1, x2}
        """
        return self._free_group

    def finitely_presented_group(self):
        """
        The finitely presented group where the rewriting system is defined.

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [ [1,2,3], [-1,-2,-3], [1,1], [2,2] ]
            sage: k = G.rewriting_system()
            sage: k.make_confluent()
            sage: k
            Rewriting system of Finitely presented group < x0, x1, x2 | x0*x1*x2, x0^-1*x1^-1*x2^-1, x0^2, x1^2 >
            with rules:
                x0^-1    --->    x0
                x1^-1    --->    x1
                x2^-1    --->    x2
                x0^2    --->    1
                x0*x1    --->    x2
                x0*x2    --->    x1
                x1*x0    --->    x2
                x1^2    --->    1
                x1*x2    --->    x0
                x2*x0    --->    x1
                x2*x1    --->    x0
                x2^2    --->    1
            sage: k.finitely_presented_group()
            Finitely presented group < x0, x1, x2 | x0*x1*x2, x0^-1*x1^-1*x2^-1, x0^2, x1^2 >
        """
        return self._fp_group

    def reduce(self, element):
        """
        Apply the rules in the rewriting system to the element, to obtain
        a reduced form.

        If the rewriting system is confluent, this reduced form is unique
        for all words representing the same element.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F/[a^2, b^3, (a*b/a)^3, b*a*b*a]
            sage: k = G.rewriting_system()
            sage: k.reduce(b^4)
            b
            sage: k.reduce(a*b*a)
            a*b*a
        """
        eg = self._fp_group(element).gap()
        egim = self._monoid_isomorphism.Image(eg)
        red = self.gap().ReducedForm(egim.UnderlyingElement())
        redfpmon = self._monoid.One().FamilyObj().ElementOfFpMonoid(red)
        reducfpgr = self._monoid_isomorphism.PreImagesRepresentative(redfpmon)
        tz = reducfpgr.UnderlyingElement().TietzeWordAbstractWord(self._free_group.gap().GeneratorsOfGroup())
        return self._fp_group(tz.sage())

    def gap(self):
        """
        The gap representation of the rewriting system.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F/[a*a,b*b]
            sage: k = G.rewriting_system()
            sage: k.gap()
            Knuth Bendix Rewriting System for Monoid( [ a, A, b, B ] ) with rules
            [ [ a*A, <identity ...> ], [ A*a, <identity ...> ],
              [ b*B, <identity ...> ], [ B*b, <identity ...> ],
              [ a^2, <identity ...> ], [ b^2, <identity ...> ] ]
        """
        return self._gap

    def rules(self):
        """
        Return the rules that form the rewriting system.

        OUTPUT:

        A dictionary containing the rules of the rewriting system.
        Each key is a word in the free group, and its corresponding
        value is the word to which it is reduced.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a*a*a,b*b*a*a]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^3, b^2*a^2 >
            with rules:
                a^3    --->    1
                b^2*a^2    --->    1

            sage: k.rules()
            {a^3: 1, b^2*a^2: 1}
            sage: k.make_confluent()
            sage: sorted(k.rules().items())
            [(a^-2, a), (a^-1*b^-1, a*b), (a^-1*b, b^-1), (a^2, a^-1),
             (a*b^-1, b), (b^-1*a^-1, a*b), (b^-1*a, b), (b^-2, a^-1),
             (b*a^-1, b^-1), (b*a, a*b), (b^2, a)]
        """
        dic = {}
        grules = self.gap().Rules()
        for i in grules:
            a, b = i
            afpmon = self._monoid.One().FamilyObj().ElementOfFpMonoid(a)
            afg = self._monoid_isomorphism.PreImagesRepresentative(afpmon)
            atz = afg.UnderlyingElement().TietzeWordAbstractWord(self._free_group.gap().GeneratorsOfGroup())
            af = self._free_group(atz.sage())
            if len(af.Tietze()) != 0:
                bfpmon = self._monoid.One().FamilyObj().ElementOfFpMonoid(b)
                bfg = self._monoid_isomorphism.PreImagesRepresentative(bfpmon)
                btz = bfg.UnderlyingElement().TietzeWordAbstractWord(self._free_group.gap().GeneratorsOfGroup())
                bf = self._free_group(btz.sage())
                dic[af] = bf
        return dic

    def is_confluent(self):
        """
        Return ``True`` if the system is confluent and ``False`` otherwise.

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [F([1,2,1,2,1,3,-1]),F([2,2,2,1,1,2]),F([1,2,3])]
            sage: k = G.rewriting_system()
            sage: k.is_confluent()
            False
            sage: k
            Rewriting system of Finitely presented group < x0, x1, x2 | (x0*x1)^2*x0*x2*x0^-1, x1^3*x0^2*x1, x0*x1*x2 >
            with rules:
                x0*x1*x2    --->    1
                x1^3*x0^2*x1    --->    1
                (x0*x1)^2*x0*x2*x0^-1    --->    1

            sage: k.make_confluent()
            sage: k.is_confluent()
            True
            sage: k
            Rewriting system of Finitely presented group < x0, x1, x2 | (x0*x1)^2*x0*x2*x0^-1, x1^3*x0^2*x1, x0*x1*x2 >
            with rules:
                x0^-1    --->    x0
                x1^-1    --->    x1
                x0^2    --->    1
                x0*x1    --->    x2^-1
                x0*x2^-1    --->    x1
                x1*x0    --->    x2
                x1^2    --->    1
                x1*x2^-1    --->    x0*x2
                x1*x2    --->    x0
                x2^-1*x0    --->    x0*x2
                x2^-1*x1    --->    x0
                x2^-2    --->    x2
                x2*x0    --->    x1
                x2*x1    --->    x0*x2
                x2^2    --->    x2^-1
        """
        return self._gap.IsConfluent().sage()

    def make_confluent(self):
        """
        Apply the Knuth-Bendix algorithm to try to transform the rewriting
        system into a confluent one.

        Note that this method does not return any object, just changes the
        rewriting system internally.

        .. WARNING::

            This algorithm is not granted to finish. Although it may be useful
            in some occasions to run it, interrupt it manually after some time
            and use then the transformed rewriting system. Even if it is not
            confluent, it could be used to reduce some words.

        ALGORITHM:

        Uses GAP's ``MakeConfluent``.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a^2,b^3,(a*b/a)^3,b*a*b*a]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^2, b^3, a*b^3*a^-1, (b*a)^2 >
            with rules:
                a^2    --->    1
                b^3    --->    1
                (b*a)^2    --->    1
                a*b^3*a^-1    --->    1

            sage: k.make_confluent()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^2, b^3, a*b^3*a^-1, (b*a)^2 >
            with rules:
                a^-1    --->    a
                a^2    --->    1
                b^-1*a    --->    a*b
                b^-2    --->    b
                b*a    --->    a*b^-1
                b^2    --->    b^-1
        """
        try:
            self._gap.MakeConfluent()
        except ValueError:
            raise ValueError('could not make the system confluent')


@richcmp_method
class FinitelyPresentedGroup(GroupMixinLibGAP, CachedRepresentation, Group, ParentLibGAP):
    """
    A class that wraps GAP's Finitely Presented Groups.

    .. WARNING::

        You should use
        :meth:`~sage.groups.free_group.FreeGroup_class.quotient` to
        construct finitely presented groups as quotients of free
        groups. Any class inheriting this one should define
        ``__reduce__ = CachedRepresentation.__reduce__``
        after importing ``CachedRepresentation``.

    EXAMPLES::

        sage: G.<a,b> = FreeGroup()
        sage: H = G / [a, b^3]
        sage: H
        Finitely presented group < a, b | a, b^3 >
        sage: H.gens()
        (a, b)

        sage: F.<a,b> = FreeGroup('a, b')
        sage: J = F / (F([1]), F([2, 2, 2]))
        sage: J is H
        True

        sage: G = FreeGroup(2)
        sage: H = G / (G([1, 1]), G([2, 2, 2]))
        sage: H.gens()
        (x0, x1)
        sage: H.gen(0)
        x0
        sage: H.ngens()
        2
        sage: H.gap()
        <fp group on the generators [ x0, x1 ]>
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement'>
    """
    Element = FinitelyPresentedGroupElement

    def __init__(self, free_group, relations, category=None, libgap_fpgroup=None):
        """
        The Python constructor.

        TESTS::

            sage: G = FreeGroup('a, b')
            sage: H = G / (G([1]), G([2])^3)
            sage: H
            Finitely presented group < a, b | a, b^3 >

            sage: F = FreeGroup('a, b')
            sage: J = F / (F([1]), F([2, 2, 2]))
            sage: J is H
            True

            sage: A5 = libgap(AlternatingGroup(5))
            sage: A5gapfp = A5.IsomorphismFpGroup().Range()
            sage: A5gapfp
            <fp group of size 60 on the generators [ A_5.1, A_5.2 ]>
            sage: A5sage = A5gapfp.sage(); A5sage;
            Finitely presented group < A_5.1, A_5.2 | A_5.1^5*A_5.2^-5, A_5.1^5*(A_5.2^-1*A_5.1^-1)^2, (A_5.1^-2*A_5.2^2)^2 >
            sage: A5sage.inject_variables()
            Traceback (most recent call last):
            ...
            ValueError: variable names have not yet been set using self._assign_names(...)

        Check that pickling works::

            sage: G = FreeGroup(2) / [2 * (1, 2, -1, -2)]
            sage: loads(dumps(G))
            Finitely presented group < x0, x1 | (x0*x1*x0^-1*x1^-1)^2 >
            sage: G.__reduce__()[1][1]
            (Free Group on generators {x0, x1}, ((x0*x1*x0^-1*x1^-1)^2,))

            sage: TestSuite(H).run()
            sage: TestSuite(J).run()
        """
        from sage.groups.free_group import is_FreeGroup
        assert is_FreeGroup(free_group)
        assert isinstance(relations, tuple)
        self._free_group = free_group
        self._relations = relations
        try:
            self._assign_names(free_group.variable_names())
        except ValueError:
            pass
        if libgap_fpgroup is None:
            libgap_fpgroup = free_group.gap() / libgap([rel.gap() for rel in relations])
        ParentLibGAP.__init__(self, libgap_fpgroup)
        Group.__init__(self, category=category)

    def __hash__(self):
        """
        Make hashable.

        EXAMPLES::

            sage: G = FreeGroup(2) / [(1, 2, 2, 1)]
            sage: G.__hash__() == hash((G.free_group(), G.relations()))
            True
        """
        return hash((self._free_group, self._relations))

    def __richcmp__(self, other, op):
        """
        Rich comparison of ``self`` and ``other``.

        EXAMPLES::

            sage: G1 = FreeGroup(2) / [(1, 2, 2, 1, 2, 1)]
            sage: G2 = libgap(G1).sage()
            sage: G1 == G2
            True
            sage: G1 is G2
            False
        """
        if not isinstance(other, self.__class__):
            from sage.structure.richcmp import op_NE
            return (op == op_NE)
        self_data = (self._free_group, self._relations)
        other_data = (other._free_group, other._relations)
        return richcmp(self_data, other_data, op)

    def _repr_(self) -> str:
        """
        Return a string representation.

        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2])^3)
            sage: H  # indirect doctest
            Finitely presented group < a, b | a, b^3 >
            sage: H._repr_()
            'Finitely presented group < a, b | a, b^3 >'
        """
        gens = ', '.join(self._free_group._gen_names)
        rels = ', '.join(str(r) for r in self.relations())
        return 'Finitely presented group ' + '< ' + gens + ' | ' + rels + ' >'

    def _latex_(self):
        """
        Return a LaTeX representation.

        OUTPUT: string; a valid LaTeX math command sequence

        TESTS::

            sage: F = FreeGroup(4)
            sage: F.inject_variables()
            Defining x0, x1, x2, x3
            sage: G = F.quotient([x0*x2, x3*x1*x3, x2*x1*x2])
            sage: G._latex_()
            '\\langle x_{0}, x_{1}, x_{2}, x_{3} \\mid x_{0}\\cdot x_{2} , x_{3}\\cdot x_{1}\\cdot x_{3} , x_{2}\\cdot x_{1}\\cdot x_{2}\\rangle'
        """
        r = '\\langle '
        for i in range(self.ngens()):
            r = r+self.gen(i)._latex_()
            if i < self.ngens()-1:
                r = r+', '
        r = r+' \\mid '
        for i in range(len(self._relations)):
            r = r+(self._relations)[i]._latex_()
            if i < len(self.relations())-1:
                r = r+' , '
        r = r+'\\rangle'
        return r

    def free_group(self):
        """
        Return the free group (without relations).

        OUTPUT: a :func:`~sage.groups.free_group.FreeGroup`

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.free_group()
            Free Group on generators {a, b, c}
            sage: H.free_group() is G
            True
        """
        return self._free_group

    def relations(self):
        """
        Return the relations of the group.

        OUTPUT: the relations as a tuple of elements of :meth:`free_group`

        EXAMPLES::

            sage: F = FreeGroup(5, 'x')
            sage: F.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: G = F.quotient([x0*x2, x3*x1*x3, x2*x1*x2])
            sage: G.relations()
            (x0*x2, x3*x1*x3, x2*x1*x2)
            sage: all(rel in F for rel in G.relations())
            True
        """
        return self._relations

    @cached_method
    def cardinality(self, limit=4096000):
        """
        Compute the cardinality of ``self``.

        INPUT:

        - ``limit`` -- integer (default: 4096000); the maximal number
          of cosets before the computation is aborted

        OUTPUT: integer or ``Infinity``; the number of elements in the group

        EXAMPLES::

            sage: G.<a,b> = FreeGroup('a, b')
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.cardinality()
            6

            sage: F.<a,b,c> = FreeGroup()
            sage: J = F / (F([1]), F([2, 2, 2]))
            sage: J.cardinality()
            +Infinity

        ALGORITHM:

            Uses GAP.

        .. WARNING::

            This is in general not a decidable problem, so it is not
            guaranteed to give an answer. If the group is infinite, or
            too big, you should be prepared for a long computation
            that consumes all the memory without finishing if you do
            not set a sensible ``limit``.
        """
        with libgap.global_context('CosetTableDefaultMaxLimit', limit):
            if not libgap.IsFinite(self.gap()):
                from sage.rings.infinity import Infinity
                return Infinity
            try:
                size = self.gap().Size()
            except ValueError:
                raise ValueError('Coset enumeration ran out of memory, is the group finite?')
        return size.sage()

    order = cardinality

    def as_permutation_group(self, limit=4096000):
        """
        Return an isomorphic permutation group.

        The generators of the resulting group correspond to the images
        by the isomorphism of the generators of the given group.

        INPUT:

        - ``limit`` -- integer (default: 4096000); the maximal number
          of cosets before the computation is aborted

        OUTPUT:

        A Sage
        :func:`~sage.groups.perm_gps.permgroup.PermutationGroup`. If
        the number of cosets exceeds the given ``limit``, a
        :exc:`ValueError` is returned.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.as_permutation_group()
            Permutation Group with generators [(1,2)(3,5)(4,6), (1,3,4)(2,5,6)]

            sage: G.<a,b> = FreeGroup()
            sage: H = G / [a^3*b]
            sage: H.as_permutation_group(limit=1000)
            Traceback (most recent call last):
            ...
            ValueError: Coset enumeration exceeded limit, is the group finite?

        ALGORITHM:

        Uses GAP's coset enumeration on the trivial subgroup.

        .. WARNING::

            This is in general not a decidable problem (in fact, it is
            not even possible to check if the group is finite or
            not). If the group is infinite, or too big, you should be
            prepared for a long computation that consumes all the
            memory without finishing if you do not set a sensible
            ``limit``.
        """
        with libgap.global_context('CosetTableDefaultMaxLimit', limit):
            try:
                trivial_subgroup = self.gap().TrivialSubgroup()
                coset_table = self.gap().CosetTable(trivial_subgroup).sage()
            except ValueError:
                raise ValueError('Coset enumeration exceeded limit, is the group finite?')
        from sage.combinat.permutation import Permutation
        from sage.groups.perm_gps.permgroup import PermutationGroup
        return PermutationGroup([
            Permutation(coset_table[2*i]) for i in range(len(coset_table)//2)])

    def direct_product(self, H, reduced=False, new_names=True):
        r"""
        Return the direct product of ``self`` with finitely presented
        group ``H``.

        Calls GAP function ``DirectProduct``, which returns the direct
        product of a list of groups of any representation.

        From [Joh1990]_ (p. 45, proposition 4): If `G`, `H` are groups
        presented by `\langle X \mid R \rangle` and `\langle Y \mid S \rangle`
        respectively, then their direct product has the presentation
        `\langle X, Y \mid R, S, [X, Y] \rangle` where `[X, Y]` denotes the
        set of commutators `\{ x^{-1} y^{-1} x y \mid x \in X, y \in Y \}`.

        INPUT:

        - ``H`` -- a finitely presented group

        - ``reduced`` -- boolean (default: ``False``); if ``True``, then
          attempt to reduce the presentation of the product group

        - ``new_names`` -- boolean (default: ``True``); if ``True``, then
          lexicographical variable names are assigned to the generators of
          the group to be returned. If ``False``, the group to be returned
          keeps the generator names of the two groups forming the direct
          product. Note that one cannot ask to reduce the output and ask
          to keep the old variable names, as they may change meaning
          in the output group if its presentation is reduced.

        OUTPUT: the direct product of ``self`` with ``H`` as a finitely
        presented group

        EXAMPLES::

            sage: G = FreeGroup()
            sage: C12 =  ( G / [G([1,1,1,1])] ).direct_product( G / [G([1,1,1])]); C12
            Finitely presented group < a, b | a^4, b^3, a^-1*b^-1*a*b >
            sage: C12.order(), C12.as_permutation_group().is_cyclic()
            (12, True)
            sage: klein = ( G / [G([1,1])] ).direct_product( G / [G([1,1])]); klein
            Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
            sage: klein.order(), klein.as_permutation_group().is_cyclic()
            (4, False)

        We can keep the variable names from ``self`` and ``H`` to examine how
        new relations are formed::

            sage: F = FreeGroup("a"); G = FreeGroup("g")
            sage: X = G / [G.0^12]; A = F / [F.0^6]
            sage: X.direct_product(A, new_names=False)
            Finitely presented group < g, a | g^12, a^6, g^-1*a^-1*g*a >
            sage: A.direct_product(X, new_names=False)
            Finitely presented group < a, g | a^6, g^12, a^-1*g^-1*a*g >

        Or we can attempt to reduce the output group presentation::

            sage: F = FreeGroup("a"); G = FreeGroup("g")
            sage: X = G / [G.0]; A = F / [F.0]
            sage: X.direct_product(A, new_names=True)
            Finitely presented group < a, b | a, b, a^-1*b^-1*a*b >
            sage: X.direct_product(A, reduced=True, new_names=True)
            Finitely presented group <  |  >

        But we cannot do both::

            sage: K = FreeGroup(['a','b'])
            sage: D = K / [K.0^5, K.1^8]
            sage: D.direct_product(D, reduced=True, new_names=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce output and keep old variable names

        TESTS::

            sage: G = FreeGroup()
            sage: Dp = (G / [G([1,1])]).direct_product( G / [G([1,1,1,1,1,1])] )
            sage: Dp.as_permutation_group().is_isomorphic(PermutationGroup(['(1,2)','(3,4,5,6,7,8)']))
            True
            sage: C7 = G / [G.0**7]; C6 =  G / [G.0**6]
            sage: C14 = G / [G.0**14]; C3 =  G / [G.0**3]
            sage: C7.direct_product(C6).is_isomorphic(C14.direct_product(C3))
            #I  Forcing finiteness test
            True
            sage: F = FreeGroup(2); D = F / [F([1,1,1,1,1]),F([2,2]),F([1,2])**2]
            sage: D.direct_product(D).as_permutation_group().is_isomorphic(
            ....: direct_product_permgroups([DihedralGroup(5),DihedralGroup(5)]))
            True

        AUTHORS:

        - Davis Shurbert (2013-07-20): initial version
        """
        from sage.groups.free_group import FreeGroup, _lexi_gen

        if not isinstance(H, FinitelyPresentedGroup):
            raise TypeError("input must be a finitely presented group")
        if reduced and not new_names:
            raise ValueError("cannot reduce output and keep old variable names")

        fp_product = libgap.DirectProduct([self.gap(), H.gap()])
        GAP_gens = fp_product.FreeGeneratorsOfFpGroup()
        if new_names:
            name_itr = _lexi_gen()  # Python generator for lexicographical variable names
            gen_names = [next(name_itr) for i in GAP_gens]
        else:
            gen_names = [str(g) for g in self.gens()] + [str(g) for g in H.gens()]
        # Build the direct product in Sage for better variable names
        ret_F = FreeGroup(gen_names)
        ret_rls = tuple([ret_F(rel_word.TietzeWordAbstractWord(GAP_gens).sage())
                         for rel_word in fp_product.RelatorsOfFpGroup()])
        ret_fpg = FinitelyPresentedGroup(ret_F, ret_rls)
        if reduced:
            ret_fpg = ret_fpg.simplified()
        return ret_fpg

    def semidirect_product(self, H, hom, check=True, reduced=False):
        r"""
        The semidirect product of ``self`` with ``H`` via ``hom``.

        If there exists a homomorphism `\phi` from a group `G` to the
        automorphism group of a group `H`, then we can define the semidirect
        product of `G` with `H` via `\phi` as the Cartesian product of `G`
        and `H` with the operation

        .. MATH::

                (g_1, h_1)(g_2, h_2) = (g_1 g_2, \phi(g_2)(h_1) h_2).

        INPUT:

        - ``H`` -- finitely presented group which is implicitly acted on
          by ``self`` and can be naturally embedded as a normal subgroup
          of the semidirect product

        - ``hom`` -- homomorphism from ``self`` to the automorphism group
          of ``H``. Given as a pair, with generators of ``self`` in the
          first slot and the images of the corresponding generators in the
          second. These images must be automorphisms of ``H``, given again
          as a pair of generators and images.

        - ``check`` -- boolean (default: ``True``); if ``False`` the defining
          homomorphism and automorphism images are not tested for validity.
          This test can be costly with large groups, so it can be bypassed
          if the user is confident that his morphisms are valid.

        - ``reduced`` -- boolean (default: ``False``); if ``True`` then the
          method attempts to reduce the presentation of the output group

        OUTPUT:

        The semidirect product of ``self`` with ``H`` via ``hom`` as a
        finitely presented group. See
        :meth:`PermutationGroup_generic.semidirect_product
        <sage.groups.perm_gps.permgroup.PermutationGroup_generic.semidirect_product>`
        for a more in depth explanation of a semidirect product.

        AUTHORS:

        - Davis Shurbert (8-1-2013)

        EXAMPLES:

        Group of order 12 as two isomorphic semidirect products::

            sage: D4 = groups.presentation.Dihedral(4)
            sage: C3 = groups.presentation.Cyclic(3)
            sage: alpha1 = ([C3.gen(0)],[C3.gen(0)])
            sage: alpha2 = ([C3.gen(0)],[C3([1,1])])
            sage: S1 = D4.semidirect_product(C3, ([D4.gen(1), D4.gen(0)],[alpha1,alpha2]))
            sage: C2 = groups.presentation.Cyclic(2)
            sage: Q = groups.presentation.DiCyclic(3)
            sage: a = Q([1]); b = Q([-2])
            sage: alpha = (Q.gens(), [a,b])
            sage: S2 = C2.semidirect_product(Q, ([C2.0],[alpha]))
            sage: S1.is_isomorphic(S2)
            #I  Forcing finiteness test
            True

        Dihedral groups can be constructed as semidirect products
        of cyclic groups::

            sage: C2 = groups.presentation.Cyclic(2)
            sage: C8 = groups.presentation.Cyclic(8)
            sage: hom = (C2.gens(), [ ([C8([1])], [C8([-1])]) ])
            sage: D = C2.semidirect_product(C8, hom)
            sage: D.as_permutation_group().is_isomorphic(DihedralGroup(8))
            True

        You can attempt to reduce the presentation of the output group::

            sage: D = C2.semidirect_product(C8, hom); D
            Finitely presented group < a, b | a^2, b^8, a^-1*b*a*b >
            sage: D = C2.semidirect_product(C8, hom, reduced=True); D
            Finitely presented group < a, b | a^2, a*b*a*b, b^8 >

            sage: C3 = groups.presentation.Cyclic(3)
            sage: C4 = groups.presentation.Cyclic(4)
            sage: hom = (C3.gens(), [(C4.gens(), C4.gens())])
            sage: C3.semidirect_product(C4, hom)
            Finitely presented group < a, b | a^3, b^4, a^-1*b*a*b^-1 >
            sage: D = C3.semidirect_product(C4, hom, reduced=True); D
            Finitely presented group < a, b | a^3, b^4, a^-1*b*a*b^-1 >
            sage: D.as_permutation_group().is_cyclic()
            True

        You can turn off the checks for the validity of the input morphisms.
        This check is expensive but behavior is unpredictable if inputs are
        invalid and are not caught by these tests::

            sage: C5 = groups.presentation.Cyclic(5)
            sage: C12 = groups.presentation.Cyclic(12)
            sage: hom = (C5.gens(), [(C12.gens(), C12.gens())])
            sage: sp = C5.semidirect_product(C12, hom, check=False); sp
            Finitely presented group < a, b | a^5, b^12, a^-1*b*a*b^-1 >
            sage: sp.as_permutation_group().is_cyclic(), sp.order()
            (True, 60)

        TESTS:

        The following was fixed in Gap-4.7.2::

            sage: C5.semidirect_product(C12, hom) == sp
            True

        A more complicated semidirect product::

            sage: C = groups.presentation.Cyclic(7)
            sage: D = groups.presentation.Dihedral(5)
            sage: id1 = ([C.0], [(D.gens(),D.gens())])
            sage: Se1 =  C.semidirect_product(D, id1)
            sage: id2 = (D.gens(), [(C.gens(),C.gens()),(C.gens(),C.gens())])
            sage: Se2 =  D.semidirect_product(C ,id2)
            sage: Dp1 = C.direct_product(D)
            sage: Dp1.is_isomorphic(Se1), Dp1.is_isomorphic(Se2)
            #I  Forcing finiteness test
            #I  Forcing finiteness test
            (True, True)

        Most checks for validity of input are left to GAP to handle::

            sage: bad_aut = ([C.0], [(D.gens(),[D.0, D.0])])
            sage: C.semidirect_product(D, bad_aut)
            Traceback (most recent call last):
            ...
            ValueError: images of input homomorphism must be automorphisms
            sage: bad_hom = ([D.0, D.1], [(C.gens(),C.gens())])
            sage: D.semidirect_product(C, bad_hom)
            Traceback (most recent call last):
            ...
            GAPError: Error, <gens> and <imgs> must be lists of same length
        """
        from sage.groups.free_group import FreeGroup, _lexi_gen

        if not isinstance(H, FinitelyPresentedGroup):
            raise TypeError("input must be a finitely presented group")

        GAP_self = self.gap()
        GAP_H = H.gap()
        auto_grp = libgap.AutomorphismGroup(H.gap())
        self_gens = [h.gap() for h in hom[0]]
        # construct image automorphisms in GAP
        GAP_aut_imgs = [libgap.GroupHomomorphismByImages(GAP_H, GAP_H, [g.gap() for g in gns],
                        [i.gap() for i in img]) for (gns, img) in hom[1]]

        # check for automorphism validity in images of operation defining homomorphism,
        # and construct the defining homomorphism.
        if check:
            if not all(a in libgap.List(libgap.AutomorphismGroup(GAP_H))
                       for a in GAP_aut_imgs):
                raise ValueError("images of input homomorphism must be automorphisms")
            GAP_def_hom = libgap.GroupHomomorphismByImages(GAP_self, auto_grp, self_gens, GAP_aut_imgs)
        else:
            GAP_def_hom = GAP_self.GroupHomomorphismByImagesNC(auto_grp, self_gens, GAP_aut_imgs)

        prod = libgap.SemidirectProduct(GAP_self, GAP_def_hom, GAP_H)
        # Convert pc group to fp group
        if prod.IsPcGroup():
            prod = libgap.Image(libgap.IsomorphismFpGroupByPcgs(prod.FamilyPcgs(), 'x'))
        if not prod.IsFpGroup():
            raise NotImplementedError("unable to convert GAP output to equivalent Sage fp group")

        # Convert GAP group object to Sage via Tietze
        # lists for readability of variable names
        GAP_gens = prod.FreeGeneratorsOfFpGroup()
        name_itr = _lexi_gen()  # Python generator for lexicographical variable names
        ret_F = FreeGroup([next(name_itr) for i in GAP_gens])
        ret_rls = tuple([ret_F(rel_word.TietzeWordAbstractWord(GAP_gens).sage())
                         for rel_word in prod.RelatorsOfFpGroup()])
        ret_fpg = FinitelyPresentedGroup(ret_F, ret_rls)
        if reduced:
            ret_fpg = ret_fpg.simplified()
        return ret_fpg

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2, 2, 2]))
            sage: H([1, 2, 1, -1]) # indirect doctest
            a*b
            sage: H([1, 2, 1, -2]) # indirect doctest
            a*b*a*b^-1
        """
        if len(args) != 1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        if x == 1:
            return self.one()
        try:
            P = x.parent()
        except AttributeError:
            return self.element_class(self, x, **kwds)
        if P is self._free_group:
            return self.element_class(self, x.Tietze(), **kwds)
        return self.element_class(self, x, **kwds)

    @cached_method
    def abelian_invariants(self):
        r"""
        Return the abelian invariants of ``self``.

        The abelian invariants are given by a list of integers
        `(i_1, \ldots, i_j)`, such that the abelianization of the group is
        isomorphic to `\ZZ / (i_1) \times \cdots \times \ZZ / (i_j)`.

        EXAMPLES::

            sage: G = FreeGroup(4, 'g')
            sage: G.inject_variables()
            Defining g0, g1, g2, g3
            sage: H = G.quotient([g1^2, g2*g1*g2^(-1)*g1^(-1), g1*g3^(-2), g0^4])
            sage: H.abelian_invariants()
            (0, 4, 4)

        ALGORITHM:

        Uses GAP.
        """
        invariants = self.gap().AbelianInvariants()
        return tuple(i.sage() for i in invariants)

    @cached_method
    def abelianization_map(self):
        r"""
        Return the abelianization map of ``self``.

        OUTPUT: the abelianization map of ``self`` as a homomorphism of
        finitely presented groups

        EXAMPLES::

            sage: G = FreeGroup(4, 'g')
            sage: G.inject_variables(verbose=False)
            sage: H = G.quotient([g1^2, g2*g1*g2^(-1)*g1^(-1), g1*g3^(-2), g0^4])
            sage: H.abelianization_map()
            Group morphism:
                From: Finitely presented group  < g0, g1, g2, g3 | g1^2, g2*g1*g2^-1*g1^-1, g1*g3^-2, g0^4 >
                To:   Finitely presented group  < f2, f3, f4 | f2^-1*f3^-1*f2*f3, f2^-1*f4^-1*f2*f4, f3^-1*f4^-1*f3*f4, f2^4, f3^4 >
            sage: g = FreeGroup(0) / []
            sage: g.abelianization_map()
            Group endomorphism of Finitely presented group  <  |  >
        """
        if not self.generators():
            return self.hom(codomain=self, im_gens=[])
        hom_ab_libgap = libgap(self).MaximalAbelianQuotient()
        ab_libgap = hom_ab_libgap.Range()
        hom_ab_fp = ab_libgap.IsomorphismFpGroup()
        ab_libgap_fp = hom_ab_fp.Range()
        hom_simply = ab_libgap_fp.IsomorphismSimplifiedFpGroup()
        ab = hom_simply.Range().sage()
        images = []
        for f in self.gens():
            f0 = hom_ab_libgap.Image(f)
            f1 = hom_ab_fp.Image(f0)
            f2 = hom_simply.Image(f1)
            L = f2.UnderlyingElement().LetterRepAssocWord()
            images.append(ab([int(j) for j in L]))
        return self.hom(codomain=ab, im_gens=images, check=False)

    @cached_method
    def abelianization_to_algebra(self, ring=QQ):
        r"""
        Return the group algebra of the abelianization of ``self``
        together with the monomials representing the generators of ``self``.

        INPUT:

        - ``ring`` -- (default: ``QQ``) the base ring for
          the group algebra of ``self``

        OUTPUT:

        - ``ab`` -- the abelianization  of ``self`` as a finitely presented group
          with a minimal number `n` of generators
        - ``R`` -- a Laurent polynomial ring with `n` variables with base ring ``ring``
        - ``ideal`` -- list of generators of an ideal ``I`` in ``R`` such that ``R/I``
          is the group algebra of the abelianization over ``ring``
        - ``image`` -- list  with the images of the generators of ``self`` in ``R/I``

        EXAMPLES::

            sage: G = FreeGroup(4, 'g')
            sage: G.inject_variables()
            Defining g0, g1, g2, g3
            sage: H = G.quotient([g1^2, g2*g1*g2^(-1)*g1^(-1), g1*g3^(-2), g0^4])
            sage: H.abelianization_to_algebra()
            (Finitely presented group  < f2, f3, f4 | f2^-1*f3^-1*f2*f3, f2^-1*f4^-1*f2*f4,
                                                      f3^-1*f4^-1*f3*f4, f2^4, f3^4 >,
             Multivariate Laurent Polynomial Ring in f2, f3, f4 over Rational Field,
             [f2^4 - 1, f3^4 - 1], [f2^-1*f3^-2, f3^-2, f4, f3])
            sage: g=FreeGroup(0) / []
            sage: g.abelianization_to_algebra()
            (Finitely presented group  <  |  >, Rational Field, [], [])
        """
        if not self.generators():
            return self, ring, [], []
        hom_ab = self.abelianization_map()
        ab = hom_ab.codomain()
        R = LaurentPolynomialRing(ring, ab.gens())
        ideal = []
        for a in ab.relations():
            a_T = a.Tietze()
            a_S = Set(a_T)
            if a_S.cardinality() == 1:
                j = a_T[0]
                m = len(a_T)
                ideal.append(R.gen(j - 1) ** m - 1)
        images0 = [hom_ab(g).Tietze() for g in self.gens()]
        images = []
        for L in images0:
            p = R.one()
            for a in L:
                if a > 0:
                    p *= R.gen(a - 1)
                elif a < 0:
                    p /= R.gen(-a - 1)
            images.append(p)
        return ab, R, ideal, images

    def simplification_isomorphism(self):
        """
        Return an isomorphism from ``self`` to a finitely presented group with
        a (hopefully) simpler presentation.

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G / [a*b*c, a*b^2, c*b/c^2]
            sage: I = H.simplification_isomorphism()
            sage: I
            Group morphism:
              From: Finitely presented group < a, b, c | a*b*c, a*b^2, c*b*c^-2 >
              To:   Finitely presented group < b |  >
            sage: I(a)
            b^-2
            sage: I(b)
            b
            sage: I(c)
            b

        TESTS::

            sage: F = FreeGroup(1)
            sage: G = F.quotient([F.0])
            sage: h = G.simplification_isomorphism(); h
            Group morphism:
              From: Finitely presented group < x | x >
              To:   Finitely presented group <  |  >
            sage: h(G.gen(0))
            1

        ALGORITHM:

        Uses GAP.
        """
        II = self.gap().IsomorphismSimplifiedFpGroup()
        cod = II.Range().sage()
        phi = [cod(II.ImageElm(x)) for x in self.gap().GeneratorsOfGroup()]
        return self.hom(codomain=cod, im_gens=phi, check=False)
        # II = self.gap().IsomorphismSimplifiedFpGroup()
        # codomain = II.Range().sage()
        # phi = lambda x: codomain(II.ImageElm(x.gap()))
        # HS = self.Hom(codomain)
        # return GroupMorphismWithGensImages(HS, phi)

    def simplified(self):
        """
        Return an isomorphic group with a (hopefully) simpler presentation.

        OUTPUT:

        A new finitely presented group. Use
        :meth:`simplification_isomorphism` if you want to know the
        isomorphism.

        EXAMPLES::

            sage: G.<x,y> = FreeGroup()
            sage: H = G /  [x ^5, y ^4, y*x*y^3*x ^3]
            sage: H
            Finitely presented group < x, y | x^5, y^4, y*x*y^3*x^3 >
            sage: H.simplified()
            Finitely presented group < x, y | y^4, y*x*y^-1*x^-2, x^5 >

        A more complicate example::

            sage: G.<e0, e1, e2, e3, e4, e5, e6, e7, e8, e9> = FreeGroup()
            sage: rels = [e6, e5, e3, e9, e4*e7^-1*e6, e9*e7^-1*e0,
            ....:         e0*e1^-1*e2, e5*e1^-1*e8, e4*e3^-1*e8, e2]
            sage: H = G.quotient(rels);  H
            Finitely presented group < e0, e1, e2, e3, e4, e5, e6, e7, e8, e9 |
            e6, e5, e3, e9, e4*e7^-1*e6, e9*e7^-1*e0, e0*e1^-1*e2, e5*e1^-1*e8, e4*e3^-1*e8, e2 >
            sage: H.simplified()
            Finitely presented group < e0 | e0^2 >
        """
        return self.simplification_isomorphism().codomain()

    def sorted_presentation(self):
        """
        Return the same presentation with the relations sorted to ensure
        equality.

        OUTPUT: a new finitely presented group with the relations sorted

        EXAMPLES::

            sage: G = FreeGroup(2) / [(1, 2, -1, -2), ()]; G
            Finitely presented group < x0, x1 | x0*x1*x0^-1*x1^-1, 1 >
            sage: G.sorted_presentation()
            Finitely presented group < x0, x1 | 1, x1^-1*x0^-1*x1*x0 >
        """
        F = FreeGroup(self.ngens())
        L0 = [r.Tietze() for r in self.relations()]
        L1 = []
        for rel in L0:
            C = [rel]
            C.extend(rel[j + 1:] + rel[:j + 1] for j in range(len(rel) - 1))
            C1 = [tuple(-j for j in reversed(l)) for l in C]
            C += C1
            C.sort()
            L1.append(C[0])
        L1.sort()
        return F / L1

    def epimorphisms(self, H):
        r"""
        Return the epimorphisms from ``self`` to `H`, up to automorphism of `H`.

        INPUT:

        - ``H`` -- another group

        EXAMPLES::

            sage: F = FreeGroup(3)
            sage: G = F / [F([1, 2, 3, 1, 2, 3]), F([1, 1, 1])]
            sage: H = AlternatingGroup(3)
            sage: for quo in G.epimorphisms(H):
            ....:   for a in G.gens():
            ....:       print(a, "|-->", quo(a))
            ....:   print("-----")
            x0 |--> ()
            x1 |--> (1,3,2)
            x2 |--> (1,2,3)
            -----
            x0 |--> (1,3,2)
            x1 |--> ()
            x2 |--> (1,2,3)
            -----
            x0 |--> (1,3,2)
            x1 |--> (1,2,3)
            x2 |--> ()
            -----
            x0 |--> (1,2,3)
            x1 |--> (1,2,3)
            x2 |--> (1,2,3)
            -----

        ALGORITHM:

        Uses libgap's GQuotients function.
        """
        # from sage.misc.misc_c import prod
        # HomSpace = self.Hom(H)
        Gg = libgap(self)
        Hg = libgap(H)
        gquotients = Gg.GQuotients(Hg)
        res = []
        # the following closure is needed to attach a specific value of quo to
        # each function in the different morphisms
        # fmap = lambda tup: (lambda a: H(prod(tup[abs(i)-1]**sign(i) for i in a.Tietze())))
        for quo in gquotients:
            # tup = tuple(H(quo.ImageElm(i.gap()).sage()) for i in self.gens())
            # fhom = GroupMorphismWithGensImages(HomSpace, fmap(tup))
            fhom = self.hom(codomain=H, im_gens=[H(quo.ImageElm(a.gap())) for a in self.gens()])
            res.append(fhom)
        return res

    def alexander_matrix(self, im_gens=None):
        """
        Return the Alexander matrix of the group.

        This matrix is given by the fox derivatives of the relations
        with respect to the generators.

        - ``im_gens`` -- (optional) the images of the generators

        OUTPUT:

        A matrix with coefficients in the group algebra. If ``im_gens`` is
        given, the coefficients will live in the same algebra as the given
        values. The result depends on the (fixed) choice of presentation.

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G.quotient([a*b/a/b, a*c/a/c, c*b/c/b])
            sage: H.alexander_matrix()
            [     1 - a*b*a^-1 a - a*b*a^-1*b^-1                 0]
            [     1 - a*c*a^-1                 0 a - a*c*a^-1*c^-1]
            [                0 c - c*b*c^-1*b^-1      1 - c*b*c^-1]

        If we introduce the images of the generators, we obtain the
        result in the corresponding algebra.

        ::

            sage: G.<a,b,c,d,e> = FreeGroup()
            sage: H = G.quotient([a*b/a/b, a*c/a/c, a*d/a/d, b*c*d/(c*d*b), b*c*d/(d*b*c)])
            sage: H.alexander_matrix()
            [              1 - a*b*a^-1          a - a*b*a^-1*b^-1                          0                          0                          0]
            [              1 - a*c*a^-1                          0          a - a*c*a^-1*c^-1                          0                          0]
            [              1 - a*d*a^-1                          0                          0          a - a*d*a^-1*d^-1                          0]
            [                         0             1 - b*c*d*b^-1   b - b*c*d*b^-1*d^-1*c^-1      b*c - b*c*d*b^-1*d^-1                          0]
            [                         0        1 - b*c*d*c^-1*b^-1             b - b*c*d*c^-1 b*c - b*c*d*c^-1*b^-1*d^-1                          0]
            sage: R.<t1,t2,t3,t4> = LaurentPolynomialRing(ZZ)
            sage: H.alexander_matrix([t1,t2,t3,t4])
            [    -t2 + 1      t1 - 1           0           0           0]
            [    -t3 + 1           0      t1 - 1           0           0]
            [    -t4 + 1           0           0      t1 - 1           0]
            [          0  -t3*t4 + 1      t2 - 1  t2*t3 - t3           0]
            [          0     -t4 + 1 -t2*t4 + t2   t2*t3 - 1           0]
        """
        rel = self.relations()
        gen = self._free_group.gens()
        return matrix(len(rel), len(gen),
                      lambda i, j: rel[i].fox_derivative(gen[j], im_gens))

    @cached_method
    def abelian_alexander_matrix(self, ring=QQ, simplified=True):
        """
        Return the Alexander matrix of the group with values in the group
        algebra of the abelianized.

        INPUT:

        - ``ring`` -- (default: ``QQ``) the base ring of the
          group algebra
        - ``simplified`` -- boolean (default: ``False``); if set to
          ``True`` use Gauss elimination and erase rows and columns

        OUTPUT:

        - ``A`` -- a matrix with coefficients in ``R``
        - ``ideal`` -- an list of generators of an ideal ``I`` of
          ``R = A.base_ring()`` such that ``R/I`` is the group algebra of the
          abelianization of ``self``

        EXAMPLES::

            sage: G.<a,b,c> = FreeGroup()
            sage: H = G.quotient([a*b/a/b, a*c/a/c, c*b/c/b])
            sage: A, ideal = H.abelian_alexander_matrix()
            sage: A
            [-f2 + 1  f1 - 1       0]
            [-f3 + 1       0  f1 - 1]
            [      0  f3 - 1 -f2 + 1]
            sage: A.base_ring()
            Multivariate Laurent Polynomial Ring in f1, f2, f3 over Rational Field
            sage: ideal
            []
            sage: G = FreeGroup(3)/[(2, 1, 1), (1, 2, 2, 3, 3)]
            sage: A, ideal = G.abelian_alexander_matrix(simplified=True); A
            [-f3^2 - f3^4 - f3^6         f3^3 + f3^6]
            sage: g = FreeGroup(1) / []
            sage: g.abelian_alexander_matrix()
            ([], [])
            sage: g.abelian_alexander_matrix()[0].base_ring()
            Univariate Laurent Polynomial Ring in f1 over Rational Field
            sage: g = FreeGroup(0) / []
            sage: A, ideal = g.abelian_alexander_matrix(); A
            []
            sage: A.base_ring()
            Rational Field
        """
        ab, R, ideal, images = self.abelianization_to_algebra(ring=ring)
        A = self.alexander_matrix(im_gens=images)
        if A.base_ring() != R:
            A = A.change_ring(R)
        if simplified:
            n, m = A.dimensions()
            if n == 0 or m == 0:
                return A, ideal
            simpli = True
            while simpli:
                i = 0
                j = 0
                unidad = False
                while not unidad and i < n and j < m:
                    p = A[i, j]
                    unidad = p.is_unit()
                    if unidad:
                        A.swap_rows(0, i)
                        A.swap_columns(0, j)
                        for k in range(1, n):
                            A.add_multiple_of_row(k, 0, -A[k, 0] * p ** -1)
                        A = A.delete_rows([0]).delete_columns([0])
                        n, m = A.dimensions()
                    else:
                        if j < m - 1:
                            j += 1
                        else:
                            i += 1
                            j = 0
                simpli = unidad
        return A, ideal

    def characteristic_varieties(self, ring=QQ, matrix_ideal=None, groebner=False):
        r"""
        Return the characteristic varieties of the group ``self``.

        There are several definitions of the characteristic varieties of a
        group `G`, see e.g. [CS1999a]_. Let `\Lambda` be the group algebra of
        `G/G'` and `\mathbb{T}` its associated algebraic variety (a torus).
        Each element `\xi\in\mathbb{T}` defines a local system of coefficients
        and the `k`-th characteristic variety is

        .. MATH::

            V_k(G) = \{\xi\in\mathbb{T}\mid \dim H^1(G;\xi)\geq k\}.

        These varieties are defined by ideals in `\Lambda`.

        INPUT:

        - ``ring`` -- (default: ``QQ``) the base ring of the group algebra
        - ``groebner`` -- boolean (default: ``False``); if set to
          ``True`` the minimal associated primes of the ideals and their
          groebner bases are computed; ignored if the base ring
          is not a field

        OUTPUT:

        A dictionary with keys the indices of the varieties. If ``groebner`` is ``False``
        the values are the ideals defining the characteristic varieties.
        If it is ``True``, lists for Gröbner bases for the ideal of each irreducible
        component, stopping when the first time a characteristic variety is empty.

        EXAMPLES::

            sage: L = [2*(i, j) + 2* (-i, -j) for i, j in ((1, 2), (2, 3), (3, 1))]
            sage: G = FreeGroup(3) / L
            sage: G.characteristic_varieties(groebner=True)
            {0: [(0,)],
             1: [(f1 - 1, f2 - 1, f3 - 1), (f1*f3 + 1, f2 - 1), (f1*f2 + 1, f3 - 1), (f2*f3 + 1, f1 - 1),
                 (f2*f3 + 1, f1 - f2), (f2*f3 + 1, f1 - f3), (f1*f3 + 1, f2 - f3)],
             2: [(f1 - 1, f2 - 1, f3 - 1), (f1 + 1, f2 - 1, f3 - 1), (f1 - 1, f2 - 1, f3 + 1),
                 (f3^2 + 1, f1 - f3, f2 - f3), (f1 - 1, f2 + 1, f3 - 1)],
             3: [(f1 - 1, f2 - 1, f3 - 1)],
             4: []}
            sage: G = FreeGroup(2)/[2*(1,2,-1,-2)]
            sage: G.characteristic_varieties()
            {0: Ideal (0) of Multivariate Laurent Polynomial Ring in f1, f2 over Rational Field,
             1: Ideal (f2 - 1, f1 - 1) of Multivariate Laurent Polynomial Ring in f1, f2 over Rational Field,
             2: Ideal (f2 - 1, f1 - 1) of Multivariate Laurent Polynomial Ring in f1, f2 over Rational Field,
             3: Ideal (1) of Multivariate Laurent Polynomial Ring in f1, f2 over Rational Field}
            sage: G.characteristic_varieties(ring=ZZ)
            {0: Ideal (0) of Multivariate Laurent Polynomial Ring in f1, f2 over Integer Ring,
             1: Ideal (2*f2 - 2, 2*f1 - 2) of Multivariate Laurent Polynomial Ring in f1, f2 over Integer Ring,
             2: Ideal (f2 - 1, f1 - 1) of Multivariate Laurent Polynomial Ring in f1, f2 over Integer Ring,
             3: Ideal (1) of Multivariate Laurent Polynomial Ring in f1, f2 over Integer Ring}
            sage: G = FreeGroup(2)/[(1,2,1,-2,-1,-2)]
            sage: G.characteristic_varieties()
            {0: Ideal (0) of Univariate Laurent Polynomial Ring in f2 over Rational Field,
             1: Ideal (-1 + 2*f2 - 2*f2^2 + f2^3) of Univariate Laurent Polynomial Ring in f2 over Rational Field,
             2: Ideal (1) of Univariate Laurent Polynomial Ring in f2 over Rational Field}
            sage: G.characteristic_varieties(groebner=True)
            {0: [0], 1: [-1 + f2, 1 - f2 + f2^2], 2: []}
            sage: G = FreeGroup(2)/[3 * (1, ), 2 * (2, )]
            sage: G.characteristic_varieties(groebner=True)
            {0: [-1 + F1, 1 + F1, 1 - F1 + F1^2, 1 + F1 + F1^2], 1: [1 - F1 + F1^2],  2: []}
            sage: G = FreeGroup(2)/[2 * (2, )]
            sage: G.characteristic_varieties(groebner=True)
            {0: [(f1 + 1,), (f1 - 1,)], 1: [(f1 + 1,), (f1 - 1, f2 - 1)], 2: []}
            sage: G = (FreeGroup(0) / [])
            sage: G.characteristic_varieties()
            {0: Principal ideal (0) of Rational Field,
             1: Principal ideal (1) of Rational Field}
            sage: G.characteristic_varieties(groebner=True)
            {0: [(0,)], 1: [(1,)]}
        """
        if self.ngens() == 0:
            if groebner:
                return {j: [(ring(j),)] for j in (0, 1)}
            return {j: ring.ideal(j) for j in (0, 1)}
        A, rels = self.abelian_alexander_matrix(ring=ring, simplified=True)
        R = A.base_ring()
        eval_1 = {x: ring(1) for x in R.gens()}
        A_scalar = A.apply_map(lambda p: p.subs(eval_1))
        n = A.ncols()
        n1 = n - A_scalar.rank()
        ideal_1 = R.ideal([x - 1 for x in R.gens()])
        S = R.polynomial_ring()
        K = R.base_ring()
        id_rels = R.ideal(rels)
        res = {}
        bound = n + 1
        for j in range(bound + 1):
            J = id_rels + A.fitting_ideal(j)
            # J = R.ideal(id_rels.gens() + A.fitting_ideal(j).gens())
            if j <= n1:
                J1 = K.ideal([K(p.subs(eval_1)) for p in J.gens()])
                if J1:
                    J *= ideal_1
            res[j] = R.ideal(J.gens_reduced())
            if R(1) in res[j].gens():
                bound = j
                break
        if not groebner or not ring.is_field():
            return res
        if R.ngens() == 1:
            res = {j: gcd(S(p) for p in res[j].gens()) for j in range(bound + 1)}
            char_var = {}
            strict = True
            j = 0
            while strict and j <= bound:
                if res[j] == 0:
                    char_var[j] = [R(0)]
                else:
                    fct = [q[0] for q in R(res[j]).factor()]
                    if fct:
                        char_var[j] = fct
                    else:
                        char_var[j] = []
                        strict = False
                j += 1
            return char_var
        char_var = {}
        strict = True
        j = 0
        while strict and j <= bound:
            LJ = res[j].minimal_associated_primes()
            fct = [id.groebner_basis() for id in LJ]
            char_var[j] = fct
            if not fct:
                strict = False
            j += 1
        return char_var

    def rewriting_system(self):
        """
        Return the rewriting system corresponding to the finitely presented
        group. This rewriting system can be used to reduce words with respect
        to the relations.

        If the rewriting system is transformed into a confluent one, the
        reduction process will give as a result the (unique) reduced form
        of an element.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: G = F / [a^2,b^3,(a*b/a)^3,b*a*b*a]
            sage: k = G.rewriting_system()
            sage: k
            Rewriting system of Finitely presented group < a, b | a^2, b^3, a*b^3*a^-1, b*a*b*a >
            with rules:
                a^2    --->    1
                b^3    --->    1
                b*a*b*a    --->    1
                a*b^3*a^-1    --->    1

            sage: G([1,1,2,2,2])
            a^2*b^3
            sage: k.reduce(G([1,1,2,2,2]))
            1
            sage: k.reduce(G([2,2,1]))
            b^2*a
            sage: k.make_confluent()
            sage: k.reduce(G([2,2,1]))
            a*b
        """
        return RewritingSystem(self)

    from sage.groups.generic import structure_description
