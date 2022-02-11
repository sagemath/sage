"""
Elements (characters) of the dual group of a finite Abelian group

To obtain the dual group of a finite Abelian group, use the
:meth:`~sage.groups.abelian_gps.abelian_group.dual_group` method::

    sage: F = AbelianGroup([2,3,5,7,8], names="abcde")
    sage: F
    Multiplicative Abelian group isomorphic to C2 x C3 x C5 x C7 x C8

    sage: Fd = F.dual_group(names="ABCDE")
    sage: Fd
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

import operator

from sage.arith.all import LCM
from sage.misc.misc_c import prod
from sage.groups.abelian_gps.element_base import AbelianGroupElementBase
from functools import reduce


def add_strings(x, z=0):
    """
    This was in sage.misc.misc but commented out. Needed to add
    lists of strings in the word_problem method below.

    Return the sum of the elements of x.  If x is empty,
    return z.

    INPUT:

    - ``x`` -- iterable

    - ``z`` -- the ``0`` that will be returned if ``x`` is empty.

    OUTPUT:

    The sum of the elements of ``x``.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.dual_abelian_group_element import add_strings
        sage: add_strings([], z='empty')
        'empty'
        sage: add_strings(['a', 'b', 'c'])
        'abc'
    """
    if len(x) == 0:
        return z
    if not isinstance(x, list):
        m = iter(x)
        y = next(m)
        return reduce(operator.add, m, y)
    else:
        return reduce(operator.add, x[1:], x[0])


def is_DualAbelianGroupElement(x):
    """
    Test whether ``x`` is a dual Abelian group element.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.dual_abelian_group import is_DualAbelianGroupElement
        sage: F = AbelianGroup(5,[5,5,7,8,9],names = list("abcde")).dual_group()
        sage: is_DualAbelianGroupElement(F)
        False
        sage: is_DualAbelianGroupElement(F.an_element())
        True
    """
    return isinstance(x, DualAbelianGroupElement)


class DualAbelianGroupElement(AbelianGroupElementBase):
    """
    Base class for abelian group elements
    """

    def __call__(self, g):
        """
        Evaluate ``self`` on a group element ``g``.

        OUTPUT:

        An element in
        :meth:`~sage.groups.abelian_gps.dual_abelian_group.DualAbelianGroup_class.base_ring`.

        EXAMPLES::

            sage: F = AbelianGroup(5, [2,3,5,7,8], names="abcde")
            sage: a,b,c,d,e = F.gens()
            sage: Fd = F.dual_group(names="ABCDE")
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

            sage: F = AbelianGroup(1, [7], names="a")
            sage: a, = F.gens()
            sage: Fd = F.dual_group(names="A", base_ring=GF(29))
            sage: A, = Fd.gens()
            sage: A(a)
            16
        """
        F = self.parent().base_ring()
        expsX = self.exponents()
        expsg = g.exponents()
        order = self.parent().gens_orders()
        N = LCM(order)
        order_not = [N / o for o in order]
        zeta = F.zeta(N)
        return F.prod(zeta ** (expsX[i] * expsg[i] * order_not[i])
                      for i in range(len(expsX)))

    def word_problem(self, words, display=True):
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

            sage: G = AbelianGroup(5,[3, 5, 5, 7, 8],names="abcde")
            sage: Gd = G.dual_group(names="abcde")
            sage: a,b,c,d,e = Gd.gens()
            sage: u = a^3*b*c*d^2*e^5
            sage: v = a^2*b*c^2*d^3*e^3
            sage: w = a^7*b^3*c^5*d^4*e^4
            sage: x = a^3*b^2*c^2*d^3*e^5
            sage: y = a^2*b^4*c^2*d^4*e^5
            sage: e.word_problem([u,v,w,x,y],display=False)
            [[b^2*c^2*d^3*e^5, 245]]

        The command e.word_problem([u,v,w,x,y],display=True) returns
        the same list but also prints $e = (b^2*c^2*d^3*e^5)^245$.
        """
        ## First convert the problem to one using AbelianGroups
        import copy
        from sage.groups.abelian_gps.abelian_group import AbelianGroup
        from sage.interfaces.all import gap
        M = self.parent()
        G = M.group()
        gens = M.variable_names()
        g = prod([G.gen(i)**(self.list()[i]) for i in range(G.ngens())])
        gap.eval("l:=One(Rationals)")            ## trick needed for LL line below to keep Sage from parsing
        s1 = "gens := GeneratorsOfGroup(%s)"%G._gap_init_()
        gap.eval(s1)
        for i in range(len(gens)):
           cmd = ("%s := gens["+str(i+1)+"]")%gens[i]
           gap.eval(cmd)
        s2 = "g0:=%s; gensH:=%s"%(str(g),words)
        gap.eval(s2)
        s3 = 'G:=Group(gens); H:=Group(gensH)'
        gap.eval(s3)
        phi = gap.eval("hom:=EpimorphismFromFreeGroup(H)")
        l1 = gap.eval("ans:=PreImagesRepresentative(hom,g0)")
        l2 = copy.copy(l1)
        l4 = []
        l3 = l1.split("*")
        for i in range(1,len(words)+1):
            l2 = l2.replace("x"+str(i),"("+str(words[i-1])+")")
        l3 = eval(gap.eval("L3:=ExtRepOfObj(ans)"))
        nn = eval(gap.eval("n:=Int(Length(L3)/2)"))
        LL1 = eval(gap.eval("L4:=List([l..n],i->L3[2*i])"))         ## note the l not 1
        LL2 = eval(gap.eval("L5:=List([l..n],i->L3[2*i-1])"))       ## note the l not 1
        if display:
            s = str(g)+" = "+add_strings(["("+str(words[LL2[i]-1])+")^"+str(LL1[i])+"*" for i in range(nn)])
            m = len(s)
            print("      ", s[:m-1], "\n")
        return [[words[LL2[i]-1],LL1[i]] for i in range(nn)]
