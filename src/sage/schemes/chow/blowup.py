# -*- coding: utf-8 -*-
r"""
The blowup of a ChowScheme Y along a ChowScheme X

Suppose given two projective non singular varieties `X` and `Y` and a morphism
of ChowSchemes `f:X\rightarrow Y`. Then the :class:`Blowup` computes the Blowup
of `Y` along `X`.

EXAMPLE (the Veronese embedding)::

    sage: P2.<h> = Proj(2, 'h')
    sage: P5.<k> = Proj(5, 'k')
    sage: f = P2.hom([2*h], P5)
    sage: g = Blowup(f)

Note that :class:`Blowup` returns a morphism as follows:

.. math::

    \begin{array}{ccc}
      \widetilde{X}&\xrightarrow{g}&\widetilde{Y}\\
      \downarrow&&\downarrow\scriptstyle{\sigma}{}\\
      X & \xrightarrow{f} & Y
   \end{array}

Hence in order to get `B = \widetilde{Y}` and to check for example for
generators, relations, tangent bundle or betti numbers of the Blowup
use the codomain of `g`::

    sage: B = Blowup(f).codomain()

In this particular example we expect an additional divisor e::

    sage: B.chowring().gens()
    (e, k)
    sage: B.chowring().rels()
    [k^6, e*k^3, e^3 - 9/2*e^2*k + 15/2*e*k^2 - 4*k^3]

The tangent bundle on `B` is equally computed::

    sage: TB = B.tangent_bundle()
    sage: TB.chern_classes()
    [1, -2*e + 6*k, -15/2*e*k + 15*k^2, 9/2*e^2*k - 93/4*e*k^2 + 28*k^3, 27/4*e^2*k^2 + 27*k^4, 12*k^5]

as well as the usual topological invariants::

    sage: B.betti_numbers()
    [1, 2, 3, 3, 2, 1]
    sage: B.euler_number()
    12
    sage: B.euler_number() == TB.chern_classes()[5].integral()
    True

Finally one can answer the classical problem of finding the smooth plane conics
tangent to five given general conics. As each tangency is a degree 6 condition
on the `\mathbb{P}^5` of all (not necessarily smooth) conics, containing the
double lines, one may compute as follows::

    sage: (e, k) = B.chowring().gens()
    sage: ((6*k - 2*e)^5).integral()
    3264

There is no restriction on the number of generators of the Chow ring
of the exceptional divisor::

    sage: P2xP2 = Proj(2, 'h') * Proj(2, 'k')
    sage: P8 = Proj(8, 'l')
    sage: f = P2xP2.hom(['h+k'], P8)  # Segre map P2xP2 -> P8
    sage: g = Blowup(f)
    sage: B = g.codomain()
    sage: B.gens()
    (e1, e2, e3, l)
    sage: B.betti_numbers()
    [1, 2, 4, 7, 8, 7, 4, 2, 1]

TESTS::

    sage: P2.<h> = Proj(2, 'h')
    sage: P5.<k> = Proj(5, 'k')
    sage: f = P2.hom([2*h], P5)
    sage: g = Blowup(f)
    sage: TestSuite(g).run(skip=["_test_category","_test_pickling"])

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
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

from sage.all import QQ
from sage.rings.polynomial.term_order import TermOrder
from sage.matrix.constructor import matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.chow.finite_ring_extension import FiniteRingExtension
from sage.schemes.chow.library.proj import ProjBundle
from sage.schemes.chow.morphism import ChowSchemeMorphism, is_chowSchemeMorphism
from sage.schemes.chow.scheme import ChowScheme


class Blowup(ChowSchemeMorphism):
    r"""
    Construct the blowup of a ChowSchemeMorphism representing an embedding
    of smooth projective varieties:

    EXAMPLE::

        sage: X.<w> = Proj(1, 'w', name='X')
        sage: Y.<h> = Proj(3, 'h', name='Y')
        sage: i = X.hom([3*w], Y)
        sage: g = Blowup(i)
        sage: XX, YY = g.domain(), g.codomain()
        sage: XX.chowring().gens()  # XX is a ProjBundle hence the generator z
        (z, w)
        sage: YY.chowring().gens()  # YY is the Blowup hence the class e
        (e, h)
        sage: YY.chowring().basis()
        [h^3, h^2, e*h, h, e, 1]
        sage: YY.chowring().intersection_matrix()
        [ 0  0  0  0  0  1]
        [ 0  0  0  1  0  0]
        [ 0  0  0  0 -3  0]
        [ 0  1  0  0  0  0]
        [ 0  0 -3  0  0  0]
        [ 1  0  0  0  0  0]
    """
    def __init__(self, f, var_name='e', proj_var_name='z', verbose=False,
                 domain_name=None, codomain_name=None,
                 latex_domain_name=None, latex_codomain_name=None):
        r"""
        Construct the class `:class:Blowup`.

        INPUT:

        - ``f`` -- a ChowSchemeMorphism

        - ``var_name`` -- an optional string representing the variable name(s)
          of the exceptional divisor. Defaults to `e` or to `e_1,\dots,e_r`
          if there are several components.

        - ``proj_var_name`` -- an optional string representing the variable
          used for the ProjBundle constructed over the domain of `f`.

        OUTPUT:

        A ChowSchemeMorphism representing the Blowup.

        TESTS::

            sage: P = PointChowScheme
            sage: W = Bundle(P, 4, [1])  # Vector space of dimension 4
            sage: I = incidence_variety(W)
            sage: X = variety_of_nets_of_quadrics(W)

            sage: IE = I.sheaves["K1"] * I.sheaves["L2"].dual()
            sage: IF = I.sheaves["K2"] * I.sheaves["L2"].dual() * I.sheaves["K1"].determinant()
            sage: ie1 = IE.chern_classes()[1]
            sage: ie2 = IE.chern_classes()[2]
            sage: ie3 = IE.chern_classes()[3]
            sage: if2 = IF.chern_classes()[2]
            sage: f = I.hom([ie1, ie2, ie3, if2], X)

            sage: H = Blowup(f).codomain()
            sage: H.betti_numbers()
            [1, 2, 6, 10, 16, 19, 22, 19, 16, 10, 6, 2, 1]
            sage: c12 = H.tangent_bundle().chern_classes()[12]
            sage: c12.integral() == H.euler_number()
            True

        """

        #######################################################################
        # Validate input
        #######################################################################
        if not is_chowSchemeMorphism(f):
            m = "Morphism in Blowup expected."
            raise ValueError(m, f)
        if not isinstance(var_name, str):
            m = "String as var_name in Blowup expected."
            raise ValueError(m, var_name)
        if not isinstance(proj_var_name, str):
            m = "String as var_name in Blowup expected."
            raise ValueError(m, var_name)

        #######################################################################
        # Preparations: Get domain, codomain, chowrings and finite ring exts...
        #######################################################################

        X, Y = f.domain(), f.codomain()

        # Get the Chow rings and Hf: HY --> HX
        HX, HY = X.chowring(), Y.chowring()
        Hf = f.chowring_morphism()

        # Get the finite ring extension associated to Hf: HY --> HX
        if verbose:
            print("Computing finite ring extension for Hf: HY --> HX...")

        A_f = FiniteRingExtension(Hf, var_name=var_name)

        #######################################################################
        # Compute XX
        #######################################################################

        if verbose:
            print("Computing Exceptional locus...")
        # Get the normal bundle N_{X/Y} from the exact sequence:
        #       0 --> TX --> f^*(TY) --> N_{X/Y} --> 0
        TX, TY = X.tangent_bundle(), Y.tangent_bundle()
        NXY = f.upperstar(TY) - TX
        codim = NXY.rank()

        # Get XX : this is GProj(N_{X/Y}^*) --> X.
        XX = ProjBundle(NXY.dual(), hyperplane_class=proj_var_name,
                        name=domain_name, latex_name=latex_domain_name)
        HXX = XX.chowring()

        #######################################################################
        # Compute YY
        #######################################################################

        # Generators: Add as many new variables to Y as module generators of HX
        # viewed as HY-module under Hf.

        if verbose:
            print("Computing Blowup...")
        vnms = A_f.nvs() + Y.variable_names()    # Strings
        mds1 = tuple(d + 1 for d in A_f.mds())   # Add 1 to the module degrees
        degs = mds1 + Y.degs()
        term_order = TermOrder('wdegrevlex', mds1) + HY.term_order()
        R = PolynomialRing(QQ, len(vnms), names=vnms, order=term_order)

        s = len(A_f.nvs())
        last_nv = R.gen(s - 1)  # t_s where t_1, ..., t_s are the t_vars
        T = matrix(1, s, [R(z) for z in A_f.nvs()])

        # Build the relations:

        # 1) We start with the relations from Y
        rels = R.ideal(HY.defining_ideal())

        # 2) Add relations from Ann_HY(HX)
        new_rels = [R(z) * R(e) for e in A_f.ann().gens() for z in T.list()]
        rels = rels + R.ideal(new_rels)

        # 3) Add module relations of HX as HY module (eg presentation matrix)
        mod_matrix = T * A_f.prm().apply_map(lambda t: R(t))
        new_rels = mod_matrix.list()
        rels = rels + R.ideal(new_rels)

        # 4) Multiplicative relations from m_i * m_j
        mm = [m * n for m in A_f.mgs() for n in A_f.mgs()]
        pd = A_f.push_down(mm).apply_map(lambda t: R(t))
        mu = matrix(s, s, (T * pd).list())
        new_rels = (T.transpose() * T - last_nv * mu).list()
        rels = rels + R.ideal(new_rels)

        # 5) Multiplicative relations from the tau-powers
        # TODO : Singular code looks more complicated
        z = HXX.cover_ring()(proj_var_name)
        u = (z ** codim).reduce(HXX.defining_ideal())
        uc = [HX(str(u.coefficient({z: d}))) for d in range(codim)]
        pd = A_f.push_down(uc).apply_map(lambda t: R(t))
        M = matrix(codim, 1, [(- last_nv) ** d for d in range(codim)])
        new_rels = (mu * pd * M - T.transpose() * (- last_nv) ** codim).list()
        rels = rels + R.ideal(new_rels)

        # 6) Finally add the additive relations from HX
        # TODO : Singular code looks more complicated
        #   Consider the universal sequence on XX=P(N^*)
        #   0 --> S --> pi^*(NX^*) --> O(1) --> 0
        #   pi = Morphism(XX, X, XX.base_chow_ring_map)
        #   pi.upperstar(NXY.dual()) - XX.o(1)
        #   In Fulton Lemma 15.4 (i) : F = S^*
        S = XX.sheaves['universal_sub']

        top = S.dual().chern_classes()[codim - 1].lift()
        top = top.reduce(HXX.defining_ideal())
        topc = [HX(str(top.coefficient({z: d}))) for d in range(codim)]
        pd = A_f.push_down(topc).apply_map(lambda t: R(t))
        JX = matrix(s, 1, f.lowerstar(list(A_f.mgs()), verbose=verbose))
        JX = JX.apply_map(lambda t: R(t))
        new_rels = (mu * pd * M - JX).list()
        rels = rels + R.ideal(new_rels)

        # The point_class is that from Y
        pc = str(Y.point_class())

        YY = ChowScheme(Y.dimension(),
                        generators=list(vnms), degrees=list(degs),
                        relations=[str(r) for r in rels.gens()],
                        point_class=pc,
                        name=codomain_name, latex_name=latex_codomain_name)

        if verbose:
            print("Computing tangent bundle of blowup...")

        SY = YY.hom(Y.gens(), Y)
        YY = YY.base_change(SY)

        j_map = [(-HXX(z)) * HXX(m) for m in
                 A_f.mgs()] + Hf.morphism_from_cover().im_gens()
        j = XX.hom([str(x) for x in j_map], YY)

        # Calculate the tangent bundle of YY:
        # 0 --> T_YY --> p^* T_Y --> j_*(S^*) --> 0 Fulton Lemma 15.4
        p = YY.base_morphism()
        TY = Y.tangent_bundle()
        # The normal bundle of XX to YY is XX.o(-1)
        TYY = p.upperstar(TY) - j.lowerstar(S.dual(), normal_bundle=XX.o(-1),
                                            verbose=verbose)
        YY.sheaves['tangent'] = TYY

        ChowSchemeMorphism.__init__(self, XX.Hom(YY), j.chowring_morphism())
