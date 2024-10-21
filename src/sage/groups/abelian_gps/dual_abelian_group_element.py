# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.rings.number_field
"""
Elements (characters) of the dual group of a finite Abelian group

To obtain the dual group of a finite Abelian group, use the
:meth:`~sage.groups.abelian_gps.abelian_group.dual_group` method::

    sage: F = AbelianGroup([2,3,5,7,8], names='abcde')
    sage: F
    Multiplicative Abelian group isomorphic to C2 x C3 x C5 x C7 x C8

    sage: Fd = F.dual_group(names='ABCDE'); Fd
    Dual of Abelian Group isomorphic to Z/2Z x Z/3Z x Z/5Z x Z/7Z x Z/8Z
    over Cyclotomic Field of order 840 and degree 192

The elements of the dual group can be evaluated on elements of the original group::

    sage: a,b,c,d,e = F.gens()
    sage: A,B,C,D,E = Fd.gens()
    sage: A*B^2*D^7
    A*B^2
    sage: A(a)
    -1
    sage: B(b)
    zeta840^140 - 1
    sage: CC(_)     # abs tol 1e-8
    -0.499999999999995 + 0.866025403784447*I
    sage: A(a*b)
    -1
    sage: (A*B*C^2*D^20*E^65).exponents()
    (1, 1, 2, 6, 1)
    sage: B^(-1)
    B^2

AUTHORS:

- David Joyner (2006-07); based on abelian_group_element.py.

- David Joyner (2006-10); modifications suggested by William Stein.

- Volker Braun (2012-11) port to new Parent base. Use tuples for immutables.
  Default to cyclotomic base ring.
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Joyner<wdjoyner@gmail.com>
#       Copyright (C) 2012 Volker Braun<vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.arith.functions import lcm as LCM
from sage.groups.abelian_gps.element_base import AbelianGroupElementBase


def is_DualAbelianGroupElement(x) -> bool:
    """
    Test whether ``x`` is a dual Abelian group element.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.groups.abelian_gps.dual_abelian_group import is_DualAbelianGroupElement
        sage: F = AbelianGroup(5, [5,5,7,8,9], names=list("abcde")).dual_group()
        sage: is_DualAbelianGroupElement(F)
        doctest:warning...
        DeprecationWarning: The function is_DualAbelianGroupElement is deprecated;
        use 'isinstance(..., DualAbelianGroupElement)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        False
        sage: is_DualAbelianGroupElement(F.an_element())
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_DualAbelianGroupElement is deprecated; "
                "use 'isinstance(..., DualAbelianGroupElement)' instead.")
    return isinstance(x, DualAbelianGroupElement)


class DualAbelianGroupElement(AbelianGroupElementBase):
    """
    Base class for abelian group elements
    """

    def __call__(self, g):
        """
        Evaluate ``self`` on a group element ``g``.

        OUTPUT: an element in
        :meth:`~sage.groups.abelian_gps.dual_abelian_group.DualAbelianGroup_class.base_ring`

        EXAMPLES::

            sage: F = AbelianGroup(5, [2,3,5,7,8], names='abcde')
            sage: a,b,c,d,e = F.gens()
            sage: Fd = F.dual_group(names='ABCDE')
            sage: A,B,C,D,E = Fd.gens()
            sage: A*B^2*D^7
            A*B^2
            sage: A(a)
            -1
            sage: B(b)
            zeta840^140 - 1
            sage: CC(B(b))    # abs tol 1e-8
            -0.499999999999995 + 0.866025403784447*I
            sage: A(a*b)
            -1

        TESTS::

            sage: F = AbelianGroup(1, [7], names='a')
            sage: a, = F.gens()
            sage: Fd = F.dual_group(names='A', base_ring=GF(29))
            sage: A, = Fd.gens()
            sage: A(a)                                                                  # needs sage.libs.pari
            16
        """
        F = self.parent().base_ring()
        expsX = self.exponents()
        expsg = g.exponents()
        order = self.parent().gens_orders()
        N = LCM(order)
        order_not = [N / o for o in order]
        zeta = F.zeta(N)
        return F.prod(zeta**(expsX[i] * expsg[i] * order_not[i])
                      for i in range(len(expsX)))

    def word_problem(self, words):
        """
        This is a rather hackish method and is included for completeness.

        The word problem for an instance of DualAbelianGroup as it can
        for an AbelianGroup. The reason why is that word problem
        for an instance of AbelianGroup simply calls GAP (which
        has abelian groups implemented) and invokes "EpimorphismFromFreeGroup"
        and "PreImagesRepresentative". GAP does not have duals of
        abelian groups implemented. So, by using the same name
        for the generators, the method below converts the problem for
        the dual group to the corresponding problem on the group
        itself and uses GAP to solve that.

        EXAMPLES::

            sage: G = AbelianGroup(5,[3, 5, 5, 7, 8], names='abcde')
            sage: Gd = G.dual_group(names='abcde')
            sage: a,b,c,d,e = Gd.gens()
            sage: u = a^3*b*c*d^2*e^5
            sage: v = a^2*b*c^2*d^3*e^3
            sage: w = a^7*b^3*c^5*d^4*e^4
            sage: x = a^3*b^2*c^2*d^3*e^5
            sage: y = a^2*b^4*c^2*d^4*e^5
            sage: e.word_problem([u,v,w,x,y])                                           # needs sage.libs.gap
            [[b^2*c^2*d^3*e^5, 245]]
        """
        from sage.libs.gap.libgap import libgap
        A = libgap.AbelianGroup(self.parent().gens_orders())
        gens = A.GeneratorsOfGroup()
        gap_g = libgap.Product([gi**Li for gi, Li in zip(gens, self.list())])
        gensH = [libgap.Product([gi**Li for gi, Li in zip(gens, w.list())])
                 for w in words]
        H = libgap.Group(gensH)

        hom = H.EpimorphismFromFreeGroup()
        ans = hom.PreImagesRepresentative(gap_g)

        resu = ans.ExtRepOfObj().sage()  # (indice, power, indice, power, etc)
        indices = resu[0::2]
        powers = resu[1::2]
        return [[words[indi - 1], powi] for indi, powi in zip(indices, powers)]
