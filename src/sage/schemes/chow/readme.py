# -*- coding: utf-8 -*-
r"""
AUTHORS:

- Manfred Lehn

- Christoph Sorger

The following are doctests of the package not supposed to figure in the
documentation. Mostly corner cases are considered to ease refactoring...

Sheaves
-------

TESTS::

    sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')

    sage: E = Sheaf(5, 0, [])  # Incorrect input in first argument.
    Traceback (most recent call last):
    ...
    TypeError: ChowScheme expected, got 5

    sage: E = Sheaf(P2, 'x', [])  # Incorrect input in second argument.
    Traceback (most recent call last):
    ...
    TypeError: Expect an integer, got x

    sage: P3 = Proj(3)
    sage: [P3.o(n).euler_characteristic() for n in range(-4,5)]
    [-1, 0, 0, 0, 1, 4, 10, 20, 35]

Grass
-----

TESTS::

    sage: G = Grass(3, 1, chern_class='h')  # Quotients of rank 1 of C^3
    sage: G.gens()
    (h,)
    sage: G.rels()
    [h^3]
    sage: G.sheaves["universal_sub"]
    Bundle(Grass(3, 1), 2, [1, -h, h^2])
    sage: G.sheaves["universal_quotient"]
    Bundle(Grass(3, 1), 1, [1, h])

    sage: G = Grass(1, 3, chern_class='h')  # Subspaces of rank 1 of C^3
    sage: G.gens()
    (h,)
    sage: G.rels()
    [h^3]
    sage: G.sheaves["universal_sub"]
    Bundle(Grass(1, 3), 1, [1, -h])
    sage: G.sheaves["universal_quotient"]
    Bundle(Grass(1, 3), 2, [1, h, h^2])

    sage: G52 = Grass(5, 2)
    sage: G35 = Grass(3, 5)
    sage: G52 == G35
    True
    sage: G25 = Grass(2, 5)
    sage: G53 = Grass(5, 3)
    sage: G53 == G25
    True
    sage: G53.point_class() == G53.chowring()(str(G52.point_class()))
    True
    sage: G53.sheaves["universal_sub"]
    Bundle(Grass(5, 3), 2, [1, -c1, c1^2 - c2])
    sage: G53.sheaves["universal_quotient"]
    Bundle(Grass(5, 3), 3, [1, c1, c2, -c1^3 + 2*c1*c2])
    sage: G52.sheaves["universal_sub"]
    Bundle(Grass(5, 2), 3, [1, -c1, c1^2 - c2, -c1^3 + 2*c1*c2])
    sage: G52.sheaves["universal_quotient"]
    Bundle(Grass(5, 2), 2, [1, c1, c2])

    sage: P = G53.chowring().cover_ring()
    sage: I = G53.chowring().defining_ideal()
    sage: I.gens()
    [c2^4, c1*c2^3, c1^2*c2^2 - 2*c2^3, c1^3*c2 - 3/2*c1*c2^2, c1^4 - c1^2*c2 - c2^2]
    sage: J = G52.chowring().defining_ideal()
    sage: J.gens()
    [c2^4, c1*c2^3, c1^2*c2^2 - c2^3, c1^3*c2 - 2*c1*c2^2, c1^4 - 3*c1^2*c2 + c2^2]
    sage: IJ = P.ideal([P(str(x)) for x in J.gens()])
    sage: P.quotient(I) == P.quotient(IJ)
    False
    sage: G53.betti_numbers() == G52.betti_numbers()
    True
    sage: G53.chowring().basis()
    [c2^3, c1*c2^2, c2^2, c1^2*c2, c1*c2, c2, c1^3, c1^2, c1, 1]
    sage: G52.chowring().basis()
    [c2^3, c1*c2^2, c2^2, c1^2*c2, c1*c2, c2, c1^3, c1^2, c1, 1]
    sage: G53.chowring().dual_basis()
    [1, 1/2*c1, 2*c1^2 - 3*c2, -c1^2 + 2*c2, -3*c1^3 + 5*c1*c2, 2*c1^2*c2 - 3*c2^2, 2*c1^3 - 3*c1*c2, -c1^2*c2 + 2*c2^2, 1/2*c1*c2^2, c2^3]
    sage: G52.chowring().dual_basis()
    [1, c1, -c1^2 + 2*c2, c1^2 - c2, -2*c1^3 + 5*c1*c2, -c1^2*c2 + 2*c2^2, c1^3 - 2*c1*c2, c1^2*c2 - c2^2, c1*c2^2, c2^3]

    sage: G31 = Grass(3, 1)
    sage: G32 = Grass(3, 2)
    sage: G31 == G32
    True

    sage: G41 = Grass(4, 1)
    sage: G43 = Grass(4, 3)
    sage: G41 == G43
    True


Twisted cubics
--------------

    Absoluter Fall::
    sage: P = PointChowScheme
    sage: W = Bundle(P, 4, [1])
    sage: I = incidence_variety(W)
    sage: IE = I.sheaves["K1"] * I.sheaves["L2"].dual()
    sage: IF = I.sheaves["K2"] * I.sheaves["L2"].dual() * I.sheaves["K1"].determinant()
    sage: X = variety_of_nets_of_quadrics(W)
    sage: AI = I.chowring()
    sage: e1 = IE.chern_classes()[1]
    sage: e2 = IE.chern_classes()[2]
    sage: e3 = IE.chern_classes()[3]
    sage: f2 = IF.chern_classes()[2]
    sage: f = I.hom([e1, e2, e3, f2], X)
    sage: g = Blowup(f, domain_name="E", codomain_name="H")
    sage: H = g.codomain()
    sage: H.betti_numbers()
    [1, 2, 6, 10, 16, 19, 22, 19, 16, 10, 6, 2, 1]
    sage: H.euler_number()
    130
    sage: TH = H.tangent_bundle()
    sage: top = TH.chern_classes()[H.dimension()]
    sage: top.integral() == H.euler_number()
    True

Relative case over P(W^*)::

    sage: P = Grass(1, 5, 'w')
    sage: W = P.sheaves["universal_quotient"]
    sage: f = map_incidence_to_nets_of_quadrics(W, domain_name='I', codomain_name='X')
    sage: g = Blowup(f, domain_name='Exc', codomain_name='H')
    sage: Exc, H = g.domain(), g.codomain()
    sage: I, X = Exc.base_chowscheme(), H.base_chowscheme()
    sage: K1, K2 = I.sheaves["K1"], I.sheaves["K2"]
    sage: L1, L2 = I.sheaves["L1"], I.sheaves["L2"]
    sage: Q = L1 + K2.dual()
    sage: Exc_Q = Exc.base_morphism().upperstar(Q)
    sage: Exc_K1 = Exc.base_morphism().upperstar(K1)
    sage: AExc = Exc_Q.symm(2) * Exc_K1.determinant() * Exc.o(-1)  # S^2Q(beta-eta)
    sage: AExc.rank()
    6
    sage: E, F, XW = X.sheaves["E"], X.sheaves["F"], X.base_morphism().upperstar(W)
    sage: AX = (F * XW.symm(2)) - (E * XW.symm(3)) + XW.symm(5); AX.rank()
    16
    sage: AX.chern_classes()[X.dimension()].integral()  # E-S Prop. 7.9
    256676750
    sage: AH = H.base_morphism().upperstar(AX)
    sage: gAExc = g.lowerstar(AExc, normal_bundle=Exc.o(-1))
    sage: V = AH - gAExc
    sage: top = V.chern_classes()[H.dimension()]
    sage: top.integral()  # The number of twisted cubics on a general quintic threefold
    317206375

"""

# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************
