r"""
Finitely generated abelian groups with GAP.

This module provides a python wrapper for abelian groups in GAP.

EXAMPLES::

    sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
    sage: AbelianGroupGap([3,5])
    Abelian group with gap, generator orders (3, 5)

For infinite abelian groups we use the GAP package ``Polycyclic``::

    sage: AbelianGroupGap([3,0])   # optional - gap_packages
    Abelian group with gap, generator orders (3, 0)

AUTHORS:

- Simon Brandhorst (2018-01-17): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP
from sage.groups.group import AbelianGroup as AbelianGroupBase
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.groups import Groups

class AbelianGroupElement_gap(ElementLibGAP):
    r"""
    An element of an abelian group via libgap.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: G = AbelianGroupGap([3,6])
        sage: G.gens()
        (f1, f2)
    """
    def __init__(self, parent, x, check=True):
        """
        The Python constructor.

        See :class:`AbelianGroupElement_gap` for details.

        INPUT:

        - ``parent`` -- an instance of :class:`AbelianGroup_gap`
        - ``x`` -- an instance of :class:`sage.libs.gap.element.GapElement`
        - ``check`` -- boolean (default: ``True``); check
          if ``x`` is an element  of the group

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,6])
            sage: g = G.an_element()
            sage: TestSuite(g).run()
        """
        if check and x not in parent.gap():
            raise ValueError("%s is not in the group %s" % (x, parent))
        ElementLibGAP.__init__(self, parent, x)

    def __hash__(self):
        r"""
        Return the hash of this element.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,4])
            sage: g = G.an_element()
            sage: g.__hash__()    # random
            1693277541873681615
        """
        return hash(self.parent()) ^ hash(self.exponents())

    def __reduce__(self):
        r"""
        Implement pickling.

        OUTPUT:

        - a tuple ``f`` such that this element is ``f[0](*f[1])``

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,4])
            sage: g = G.an_element()
            sage: g == loads(dumps(g))
            True
            sage: g.__reduce__()
            (Abelian group with gap, generator orders (3, 2, 4), ((1, 1, 1),))
        """
        return self.parent(), (self.exponents(),)

    def exponents(self):
        r"""
        Return the tuple of exponents of this element.

        OUTPUT:

        - a tuple of integers

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4,7,9])
            sage: gens = G.gens()
            sage: g = gens[0]^2 * gens[1]^4 * gens[2]^8
            sage: g.exponents()
            (2, 4, 8)
            sage: S = G.subgroup(G.gens()[:1])
            sage: s = S.gens()[0]
            sage: s
            f1
            sage: s.exponents()
            (1,)

        It can handle quite large groups too::

            sage: G = AbelianGroupGap([2^10, 5^10])
            sage: f1, f2 = G.gens()
            sage: g = f1^123*f2^789
            sage: g.exponents()
            (123, 789)

        .. WARNING::

            Crashes for very large groups.

        .. TODO::

            Make exponents work for very large groups.
            This could be done by using Pcgs in gap.
        """

        P = self.parent()
        # better than Factorization as this does not create the
        # whole group in memory.
        f = P.gap().EpimorphismFromFreeGroup()
        x = f.PreImagesRepresentative(self.gap())
        L = x.ExtRepOfObj().sage()
        Lgens = L[::2]
        Lexpo = L[1::2]
        exp = []
        orders = P.gens_orders()
        i = 0
        for k in range(len(P.gens())):
            if k + 1 not in Lgens:
                exp.append(0)
            else:
                i = Lgens.index(k + 1)
                exp.append(Lexpo[i] % orders[k])
        return tuple(exp)

    def order(self):
        r"""
        Return the order of this element.

        OUTPUT:

        - an integer or infinity

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4])
            sage: g = G.gens()[0]
            sage: g.order()
            4
            sage: G = AbelianGroupGap([0])          # optional - gap_packages
            sage: g = G.gens()[0]                   # optional - gap_packages
            sage: g.order()                         # optional - gap_packages
            +Infinity
        """
        return self.gap().Order().sage()

class AbelianGroupElement_polycyclic(AbelianGroupElement_gap):
    r"""
    An element of an abelian group using the GAP package ``Polycyclic``.

    TESTS::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: G = AbelianGroupGap([4,7,0])          # optional - gap_packages
        sage: TestSuite(G.an_element()).run()       # optional - gap_packages
    """
    def exponents(self):
        r"""
        Return the tuple of exponents of ``self``.

        OUTPUT:

        - a tuple of integers

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4,7,0])          # optional - gap_packages
            sage: gens = G.gens()                       # optional - gap_packages
            sage: g = gens[0]^2 * gens[1]^4 * gens[2]^8 # optional - gap_packages
            sage: g.exponents()                         # optional - gap_packages
            (2, 4, 8)

        Efficiently handles very large groups::

            sage: G = AbelianGroupGap([2^30,5^30,0])    # optional - gap_packages
            sage: f1, f2, f3 = G.gens()                 # optional - gap_packages
            sage: (f1^12345*f2^123456789).exponents()   # optional - gap_packages
            (12345, 123456789, 0)
        """
        return tuple(self.gap().Exponents().sage())

class AbelianGroup_gap(UniqueRepresentation, GroupMixinLibGAP, ParentLibGAP, AbelianGroupBase):
    r"""
    Finitely generated abelian groups implemented in GAP.

    Needs the gap package ``Polycyclic`` in case the group is infinite.

    INPUT:

    - ``G`` -- a GAP group
    - ``category`` -- a category
    - ``ambient`` -- (optional) an :class:`AbelianGroupGap`

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: G = AbelianGroupGap([3, 2, 5])
        sage: G
        Abelian group with gap, generator orders (3, 2, 5)
    """
    def __init__(self, G, category, ambient=None):
        r"""
        Create an instance of this class.

        See :class:`AbelianGroup_gap` for details

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,5])
            sage: TestSuite(G).run()
        """
        AbelianGroupBase.__init__(self, category=category)
        ParentLibGAP.__init__(self, G, ambient=ambient)

    Element = AbelianGroupElement_gap

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` canonically coerces to ``self``.

        INPUT:

        - ``S`` -- anything

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: G._coerce_map_from_(S)
            True
            sage: S._coerce_map_from_(G)
            False
        """
        if isinstance(S, AbelianGroup_gap):
            return S.is_subgroup_of(self)
        return super(AbelianGroup_gap, self)._coerce_map_from_(S)

    def _element_constructor_(self, x, check=True):
        r"""
        Defines coercions and conversions.

        INPUT:

        - ``x`` -- an element of this group, a GAP element

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3])
            sage: A = AbelianGroup([2,3])
            sage: a = A.an_element()
            sage: a
            f0*f1
            sage: G(a)
            f1*f2
            sage: A = AdditiveAbelianGroup([2,3])
            sage: a = A.an_element()
            sage: a
            (1, 0)
            sage: G(a)
            f1

        For general ``fgp_modules`` conversion is implemented if our
        group is in Smith form::

            sage: G = AbelianGroupGap([6])
            sage: A = ZZ^2
            sage: e0,e1 = A.gens()
            sage: A = A / A.submodule([2*e0, 3*e1])
            sage: a = 2 * A.an_element()
            sage: a
            (2)
            sage: G(a)
            f2

        TESTS:

        Document that :trac:`31428` is fixed::

            sage: A = AbelianGroupGap([])
            sage: A([]) == A.one()
            True
        """
        if isinstance(x, AbelianGroupElement_gap):
            try:
                if x in self._cover:
                    x = self.gap().NaturalHomomorphism().Image(x.gap())
                else:
                    x = x.gap()
            except AttributeError:
                x = x.gap()
        elif x == 1 or x == ():
            x = self.gap().Identity()
        elif not isinstance(x, GapElement):
            from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
            from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroupElement
            from sage.modules.fg_pid.fgp_element import FGP_Element
            if isinstance(x, AbelianGroupElement):
                exp = x.exponents()
            elif isinstance(x, AdditiveAbelianGroupElement):
                exp = x._hermite_lift()
            elif isinstance(x, FGP_Element):
                exp = x.vector()
            else:
                from sage.modules.free_module_element import vector
                exp = vector(ZZ, x)
            # turn the exponents into a gap element
            gens_gap = self.gens()
            orders = self.gens_orders()
            if len(exp) != len(gens_gap):
                raise ValueError("input does not match the number of generators")
            x = self.one()
            for g, e, m in zip(gens_gap, exp, orders):
                if m != 0:
                    e = e % m
                x *= g**e
            x = x.gap()
        return self.element_class(self, x, check=check)

    def all_subgroups(self):
        r"""
        Return the list of all subgroups of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2, 3])
            sage: G.all_subgroups()
            [Subgroup of Abelian group with gap, generator orders (2, 3) generated by (),
             Subgroup of Abelian group with gap, generator orders (2, 3) generated by (f1,),
             Subgroup of Abelian group with gap, generator orders (2, 3) generated by (f2,),
             Subgroup of Abelian group with gap, generator orders (2, 3) generated by (f2, f1)]
        """
        subgroups_gap = self.gap().AllSubgroups()
        subgroups_sage = []
        for G in subgroups_gap:
            S = self.subgroup(G.GeneratorsOfGroup())
            subgroups_sage.append(S)
        return subgroups_sage

    def automorphism_group(self):
        r"""
        Return the group of automorphisms of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2, 3])
            sage: G.aut()
            Full group of automorphisms of Abelian group with gap, generator orders (2, 3)
        """
        from sage.groups.abelian_gps.abelian_aut import AbelianGroupAutomorphismGroup
        return AbelianGroupAutomorphismGroup(self)

    aut = automorphism_group

    def is_trivial(self):
        r"""
        Return ``True`` if this group is the trivial group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([])
            sage: G
            Abelian group with gap, generator orders ()
            sage: G.is_trivial()
            True
            sage: AbelianGroupGap([1]).is_trivial()
            True
            sage: AbelianGroupGap([1,1,1]).is_trivial()
            True
            sage: AbelianGroupGap([2]).is_trivial()
            False
            sage: AbelianGroupGap([2,1]).is_trivial()
            False
        """
        return 1 == self.order()

    def identity(self):
        r"""
        Return the identity element of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4,10])
            sage: G.identity()
            1
        """
        return self.one()

    @cached_method
    def elementary_divisors(self):
        r"""
        Return the elementary divisors of this group.

        See :meth:`sage.groups.abelian_gps.abelian_group_gap.elementary_divisors`.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: G.elementary_divisors()
            (2, 60)
        """
        ediv = self.gap().AbelianInvariants().sage()
        from sage.matrix.constructor import diagonal_matrix
        ed = diagonal_matrix(ZZ, ediv).elementary_divisors()
        return tuple(d for d in ed if d != 1)

    @cached_method
    def exponent(self):
        r"""
        Return the exponent of this abelian group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,7])
            sage: G
            Abelian group with gap, generator orders (2, 3, 7)
            sage: G = AbelianGroupGap([2,4,6])
            sage: G
            Abelian group with gap, generator orders (2, 4, 6)
            sage: G.exponent()
            12
        """
        return self.gap().Exponent().sage()

    @cached_method
    def gens_orders(self):
        r"""
        Return the orders of the generators.

        Use :meth:`elementary_divisors` if you are looking for an
        invariant of the group.

        OUTPUT:

        - a tuple of integers

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: Z2xZ3 = AbelianGroupGap([2,3])
            sage: Z2xZ3.gens_orders()
            (2, 3)
            sage: Z2xZ3.elementary_divisors()
            (6,)
            sage: Z6 = AbelianGroupGap([6])
            sage: Z6.gens_orders()
            (6,)
            sage: Z6.elementary_divisors()
            (6,)
            sage: Z2xZ3.is_isomorphic(Z6)
            True
            sage: Z2xZ3 is Z6
            False
        """
        from sage.rings.infinity import Infinity
        orders = []
        for g in self.gens():
            order = g.order()
            if order == Infinity:
                order = 0
            orders.append(order)
        return tuple(orders)

    def is_subgroup_of(self, G):
        r"""
        Return if ``self`` is a subgroup of ``G`` considered in
        the same ambient group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S1 = G.subgroup(gen)
            sage: S1.is_subgroup_of(G)
            True
            sage: S2 = G.subgroup(G.gens()[1:])
            sage: S2.is_subgroup_of(S1)
            False
        """
        if not isinstance(G, AbelianGroup_gap):
            raise ValueError("input must be an instance of AbelianGroup_gap")
        if not self.ambient() is G.ambient():
            return False
        return G.gap().IsSubsemigroup(self).sage()

    def _subgroup_constructor(self, libgap_subgroup):
        r"""
        Create a subgroup from the input.

        See :class:`~sage.groups.libgap_wrapper`. Override this in derived
        classes.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,3,4,5])
            sage: S = A.subgroup(A.gens()[:1])
            sage: T = A._subgroup_constructor(S.gap())
            sage: S is T
            True
        """
        ambient = self.ambient()
        generators = libgap_subgroup.GeneratorsOfGroup()
        generators = tuple([ambient(g) for g in generators])
        return AbelianGroupSubgroup_gap(ambient, generators)

    def quotient(self, N):
        r"""
        Return the quotient of this group by the normal subgroup `N`.

        INPUT:

        - ``N`` -- a subgroup
        - ``check`` -- bool (default: ``True``) check if `N` is normal

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,3,4,5])
            sage: S = A.subgroup(A.gens()[:1])
            sage: A.quotient(S)
            Quotient abelian group with generator orders (1, 3, 4, 5)
        """
        if not isinstance(N, AbelianGroup_gap):
            raise TypeError("not an abelian group")
        if not N.is_subgroup_of(self):
            raise ValueError("not a subgroup")
        return AbelianGroupQuotient_gap(self, N)

    def subgroup(self, gens):
        r"""
        Return the subgroup of this group generated by ``gens``.

        INPUT:

        - ``gens`` -- a list of elements coercible into this group

        OUTPUT:

        - a subgroup

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: S
            Subgroup of Abelian group with gap, generator orders (2, 3, 4, 5)
             generated by (f1, f2)
            sage: g = G.an_element()
            sage: s = S.an_element()
            sage: g * s
            f2^2*f3*f5
            sage: G = AbelianGroupGap([3,4,0,2])     # optional - gap_packages
            sage: gen = G.gens()[:2]                 # optional - gap_packages
            sage: S = G.subgroup(gen)                # optional - gap_packages
            sage: g = G.an_element()                 # optional - gap_packages
            sage: s = S.an_element()                 # optional - gap_packages
            sage: g * s                              # optional - gap_packages
            g1^2*g2^2*g3*g4

        TESTS::

            sage: h = G.gens()[3]
            sage: h in S
            False
        """
        gens = tuple([self(g) for g in gens])
        return AbelianGroupSubgroup_gap(self.ambient(), gens)

class AbelianGroupGap(AbelianGroup_gap):
    r"""
    Abelian groups implemented using GAP.

    INPUT:

    - ``generator_orders`` -- a list of nonnegative integers where `0`
      gives a factor isomorphic to `\ZZ`

    OUTPUT:

    - an abelian group

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: AbelianGroupGap([3,6])
        Abelian group with gap, generator orders (3, 6)
        sage: AbelianGroupGap([3,6,5])
        Abelian group with gap, generator orders (3, 6, 5)
        sage: AbelianGroupGap([3,6,0])      # optional - gap_packages
        Abelian group with gap, generator orders (3, 6, 0)

    .. WARNING::

        Needs the GAP package ``Polycyclic`` in case the group is infinite.
    """
    @staticmethod
    def __classcall_private__(cls, generator_orders):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A1 = AbelianGroupGap((2,3,4))
            sage: A2 = AbelianGroupGap([4/2,3,4])
            sage: A1 is A2
            True
        """
        generator_orders = tuple([ZZ(e) for e in generator_orders])
        if any(e < 0 for e in generator_orders):
            return ValueError("generator orders must be nonnegative")
        return super(AbelianGroupGap, cls).__classcall__(cls, generator_orders)

    def __init__(self, generator_orders):
        r"""
        Constructor.

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroup((2,3,4))
            sage: TestSuite(A).run()
        """
        category = Groups().Commutative()
        if 0 in generator_orders:
            if not libgap.LoadPackage("Polycyclic"):
                raise ImportError("unable to import polycyclic package")
            G = libgap.eval("AbelianPcpGroup(%s)" % list(generator_orders))
            category = category.Infinite()
            self.Element = AbelianGroupElement_polycyclic
        else:
            G = libgap.AbelianGroup(generator_orders)
            category = category.Finite().Enumerated()
        AbelianGroup_gap.__init__(self, G, category=category)

    def _latex_(self):
        r"""
        Return a latex representation of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,6])
            sage: G._latex_()
            \text{\texttt{Abelian group with gap, generator orders }} \left(2, 6\right)
        """
        from sage.misc.latex import latex
        base = r"\text{\texttt{Abelian group with gap, generator orders }}"
        return base + latex(self.gens_orders())

    def _repr_(self):
        r"""
        Return a string representation of this group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,6])
            sage: G._repr_()
            'Abelian group with gap, generator orders (2, 6)'
        """
        return "Abelian group with gap, generator orders " + str(self.gens_orders())

    def __reduce__(self):
        r"""
        Implements pickling.

        We have to work around the fact that gap does not provide pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,5])
            sage: G == loads(dumps(G))
            True
            sage: G is loads(dumps(G))
            True
        """
        return AbelianGroupGap, (self.gens_orders(),)

class AbelianGroupSubgroup_gap(AbelianGroup_gap):
    r"""
    Subgroups of abelian groups with GAP.

    INPUT:

    - ``ambient`` -- the ambient group
    - ``gens`` -- generators of the subgroup

    .. NOTE::

        Do not construct this class directly. Instead use
        :meth:`~sage.groups.abelian_groups.AbelianGroupGap.subgroup`.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: G = AbelianGroupGap([2,3,4,5])
        sage: gen = G.gens()[:2]
        sage: S = G.subgroup(gen)
    """
    def __init__(self, ambient, gens):
        r"""
        Initialize this subgroup.

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap, AbelianGroupSubgroup_gap
            sage: G = AbelianGroupGap([])
            sage: gen = G.gens()
            sage: A = AbelianGroupSubgroup_gap(G, gen)
            sage: TestSuite(A).run()

        Check that we are in the correct category::

            sage: G = AbelianGroupGap([2,3,0])      # optional - gap_packages
            sage: g = G.gens()                      # optional - gap_packages
            sage: H1 = G.subgroup([g[0],g[1]])      # optional - gap_packages
            sage: H1 in Groups().Finite()           # optional - gap_packages
            True
            sage: H2 = G.subgroup([g[0],g[2]])      # optional - gap_packages
            sage: H2 in Groups().Infinite()         # optional - gap_packages
            True
        """
        gens_gap = tuple([g.gap() for g in gens])
        G = ambient.gap().Subgroup(gens_gap)
        from sage.rings.infinity import Infinity
        category = Groups().Commutative()
        if G.Size().sage() < Infinity:
            category = category.Finite()
        else:
            category = category.Infinite()
        category = category.Subobjects()
        AbelianGroup_gap.__init__(self, G, ambient=ambient, category=category)

    def _repr_(self):
        r"""
        Return a string representation of this subgroup.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: S
            Subgroup of Abelian group with gap, generator orders (2, 3, 4, 5)
             generated by (f1, f2)
        """
        return "Subgroup of %s generated by %s"%(self.ambient(),self.gens())

    def __reduce__(self):
        r"""
        Implements pickling.

        We have to work around the fact that gap does not provide pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: S == loads(dumps(S))
            True
            sage: S is loads(dumps(S))
            True
        """
        amb = self.ambient()
        # avoid infinite loop
        gens = tuple([amb(g) for g in self.gens()])
        return amb.subgroup, (gens,)

    def lift(self, x):
        """
        Coerce to the ambient group.

        The terminology comes from the category framework and the more general notion of a subquotient.

        INPUT:

        - ``x`` -- an element of this subgroup

        OUTPUT:

        The corresponding element of the ambient group

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4])
            sage: g = G.gen(0)
            sage: H = G.subgroup([g^2])
            sage: h = H.gen(0); h
            f2
            sage: h.parent()
            Subgroup of Abelian group with gap, generator orders (4,) generated by (f2,)
            sage: H.lift(h)
            f2
            sage: H.lift(h).parent()
            Abelian group with gap, generator orders (4,)
        """
        return self.ambient()(x)

    def retract(self, x):
        """
        Convert an element of the ambient group into this subgroup.

        The terminology comes from the category framework and the more general notion of a subquotient.

        INPUT:

        - ``x`` -- an element of the ambient group that actually lies in this subgroup.

        OUTPUT:

        The corresponding element of this subgroup

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([4])
            sage: g = G.gen(0)
            sage: H = G.subgroup([g^2])
            sage: H.retract(g^2)
            f2
            sage: H.retract(g^2).parent()
            Subgroup of Abelian group with gap, generator orders (4,) generated by (f2,)
        """
        return self(x)

class AbelianGroupQuotient_gap(AbelianGroup_gap):
    r"""
    Quotients of abelian groups by a subgroup.

    .. NOTE::

        Do not call this directly. Instead use :meth:`quotient`.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([4,3])
        sage: N = A.subgroup([A.gen(0)^2])
        sage: Q1 = A.quotient(N)
        sage: Q1
        Quotient abelian group with generator orders (2, 3)
        sage: Q2 = Q1.quotient(Q1.subgroup(Q1.gens()[:1]))
        sage: Q2
        Quotient abelian group with generator orders (1, 3)
    """
    def __init__(self, G, N):
        r"""
        Constructor.

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: Q = G.quotient(S)
            sage: TestSuite(Q).run()
        """
        self._cover = G
        self._relations = N
        libgap_group = G.gap().FactorGroup(N.gap())
        category = G.category()
        AbelianGroup_gap.__init__(self, libgap_group, category=category, ambient=None)

    def _repr_(self):
        r"""
        Return a string representation of this subgroup.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: G.quotient(S)
            Quotient abelian group with generator orders (1, 1, 4, 5)
        """
        return "Quotient abelian group with generator orders " + str(
               self.gens_orders())

    def __reduce__(self):
        r"""
        Implements pickling.

        We have to work around the fact that gap does not provide pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([4])
            sage: N = A.subgroup([A.gen(0)^2])
            sage: Q = A.quotient(N)
            sage: Q is loads(dumps(Q))
            True
        """
        G = self._cover
        N = self._relations
        return G.quotient, (N, )

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` canonically coerces to ``self``.

        INPUT:

        - ``S`` -- anything

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: Q = G.quotient(S)
            sage: Q._coerce_map_from_(G)
            True
            sage: G._coerce_map_from_(Q)
            False
        """
        if isinstance(S, AbelianGroup_gap):
            return self._cover._coerce_map_from_(S)


    def cover(self):
        r"""
        Return the covering group of this quotient group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: Q = G.quotient(S)
            sage: Q.cover() is G
            True
        """
        return self._cover

    def relations(self):
        r"""
        Return the relations of this quotient group.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([2,3,4,5])
            sage: gen = G.gens()[:2]
            sage: S = G.subgroup(gen)
            sage: Q = G.quotient(S)
            sage: Q.relations() is S
            True
        """
        return self._relations

    @cached_method
    def natural_homomorphism(self):
        r"""
        Return the defining homomorphism into ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([4])
            sage: N = A.subgroup([A.gen(0)^2])
            sage: Q = A.quotient(N)
            sage: Q.natural_homomorphism()
            Group morphism:
            From: Abelian group with gap, generator orders (4,)
            To:   Quotient abelian group with generator orders (2,)
        """
        phi = self.gap().NaturalHomomorphism()
        Hom = self._cover.Hom(self)
        return Hom(phi)

    def lift(self, x):
        r"""
        Lift an element to the cover.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([4])
            sage: N = A.subgroup([A.gen(0)^2])
            sage: Q = A.quotient(N)
            sage: Q.lift(Q.0)
            f1
        """
        x = self(x)
        return self.natural_homomorphism().lift(x)
