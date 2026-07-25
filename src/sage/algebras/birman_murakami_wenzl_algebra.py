# sage.doctest: needs sage.combinat sage.groups sage.modules
r"""
Birman-Murakami-Wenzl Algebras

Let `R` be an integral domain with invertible elements `l` and `m`.
The  Birman-Murakami-Wenzl algebra `BMW_n = BMW_n(R; l, m)` is the
unital `R` algebra with invertible generators `G_i` and
non-invertible generators `E_i` for `1 \le i \le n-1` and relations

.. MATH::

    \begin{aligned}
    G_{i}G_{j} & = G_{j}G_{i},\mathrm {if} \left\vert i-j\right\vert \geq 2,\\
    G_{i}G_{i+1}G_{i} & = G_{i+1}G_{i}G_{i+1},\\
    E_{i}E_{i\pm 1}E_{i} & = E_{i},\\
    G_{i}+{G_{i}}^{-1} & = m(1+E_{i}),\\
    G_{i\pm 1}G_{i}E_{i\pm 1} & = E_{i}G_{i\pm 1}G_{i} = E_{i}E_{i\pm 1},\\
    G_{i\pm 1}E_{i}G_{i\pm 1} & = {G_{i}}^{-1}E_{i\pm 1}{G_{i}}^{-1},\\
    G_{i\pm 1}E_{i}E_{i\pm 1} & = {G_{i}}^{-1}E_{i\pm 1},\\
    E_{i\pm 1}E_{i}G_{i\pm 1} & = E_{i\pm 1}{G_{i}}^{-1},\\
    E_{i}G_{i\pm 1}E_{i} & = lE_{i},\\
    G_{i}E_{i} = E_{i}G_{i} & = l^{-1} E_{i}.
    \end{aligned}

There is a realizable presentation of `BMW_n` as the algebra of `n`-to-`n`-tangle
diagrams modulo regular isotopy, commonly referred to as the *Kauffman tangle
algebra* (see, for example, Section 6 in [MW2010]_).

Our implementation essentially uses this realization. It establishes a direct
connection to the *Kauffman link polynomial*
(:meth:`~sage.knots.links.Link.kauffman_polynomial`) or more precisely to its
regular isotopy version, which is implemented in the method :meth:`markov_trace`.

In this identification the relation in the fourth row corresponds to the skein
relation between these tangles and the relation in the last row to an untwist
relation. All other relations are consequences of Reidemeister moves applied
to the tangles.

Let us verify the relations above for `n = 4`::

    sage: BMW4.<G1, G2, G3, E1, E2, E3> = algebras.BirmanMurakamiWenzl(4)
    sage: G1*G3 == G3*G1, E1*E3 == E3*E1
    (True, True)
    sage: G1*G2*G1 == G2*G1*G2, G3*G2*G3 == G2*G3*G2
    (True, True)
    sage: E1*E2*E1 == E1, E2*E1*E2 == E2, E2*E3*E2 == E2, E3*E2*E3 == E3
    (True, True, True, True)
    sage: G1*G2*E1 == E2*G1*G2 == E2*E1, G2*G1*E2 == E1*G2*G1 == E1*E2
    (True, True)
    sage: G2*G3*E2 == E3*G2*G3 == E3*E2, G3*G2*E3 == E2*G3*G2 == E2*E3
    (True, True)
    sage: G1*E2*G1 == ~G2*E1*~G2, G2*E1*G2 == ~G1*E2*~G1
    (True, True)
    sage: G2*E3*G2 == ~G3*E2*~G3, G3*E2*G3 == ~G2*E3*~G2
    (True, True)
    sage: G2*E1*E2 == ~G1*E2, G1*E2*E1 == ~G2*E1, G2*E3*E2 == ~G3*E2, G3*E2*E3 == ~G2*E3
    (True, True, True, True)
    sage: E1*E2*G1 == E1*~G2, E2*E1*G2 == E2*~G1, E3*E2*G3 == E3*~G2, E2*E3*G2 == E2*~G3
    (True, True, True, True)
    sage: O = BMW4.one()
    sage: l, m = BMW4.base_ring().gens()
    sage: G1 + ~G1 == m*(O + E1), G2 + ~G2 == m*(O + E2), G3 + ~G3 == m*(O + E3)
    (True, True, True)
    sage: G1*E1 == E1*G1 == ~l*E1, G2*E2 == E2*G2 == ~l*E2,  G3*E3 == E3*G3 == ~l*E3
    (True, True, True)
    sage: E1*G2*E1 == l*E1,  E2*G1*E2 == l*E2, E3*G2*E3 == l*E3,  E2*G3*E2 == l*E2
    (True, True, True, True)

Note that there are other conventions concerning the the skein and untwist relations.
They can be choosen in the declaration::

    sage: BMW2.<g1, e1> = algebras.BirmanMurakamiWenzl(2, skein_normalization=(-1, -1, -1))
    sage: g1**2
    (-l*m)*e1 + m*g1 + o1
    sage: BMW2.<g1, e1> = algebras.BirmanMurakamiWenzl(2, skein_normalization=(-1, 1, 1))
    sage: g1**2
    (l^-1*m)*e1 + m*g1 + o1
    sage: BMW2.<g1, e1> = algebras.BirmanMurakamiWenzl(2, skein_normalization=(1, -1, 1))
    sage: g1**2
    (-l^-1*m)*e1 + m*g1 + (-1)*o1
    sage: BMW2.<g1, e1> = algebras.BirmanMurakamiWenzl(2, skein_normalization=(1, 1, -1))
    sage: g1**2
    l*m*e1 + m*g1 + (-1)*o1
    sage: BMW2.<g1, e1> = algebras.BirmanMurakamiWenzl(2, skein_normalization=(-1, -1, 1))
    sage: g1**2
    (-l^-1*m)*e1 + m*g1 + o1

Here the first two signs in the ``skein_normalization`` triple are used to change the two
``+`` signs in the skein relation to ``-``. The third sign is used to invert the untwist
factor. The convention showed in the first line is used in [MW2010]_ the one in the last
line is used in [EG2017]_. The default ``(1, 1, 1)`` corresponds the convention given on
Wikipedia.

To change the base ring variables use the keyword ``params``::

    sage: BMW2.<g1, e1> = algebras.BirmanMurakamiWenzl(2, params='a, z')
    sage: g1**2
    (a^-1*z)*e1 + z*g1 + (-1)*o1

AUTHORS:

- Sebastian Oehms Jan 2026: initial version

REFERENCES:

- :wikipedia:`Birman%E2%80%93Wenzl_algebra`
- :wikipedia:`Kauffman_polynomial`
- :wikipedia:`Reidemeister_move`
- [MW2010]_
- [EG2017]_
"""
# ###########################################################################
#       Copyright (C) 2026 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ###########################################################################
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.misc.verbose import get_verbose, verbose
from sage.modules.free_module_element import vector
from sage.monoids.tangles import KauffmanTangle, KauffmanTangles
from sage.rings.integer_ring import ZZ


class BirmanMurakamiWenzlElement(CombinatorialFreeModule.Element):
    r"""
    An element of a :class:`BirmanMurakamiWenzlAlgebra`.

    For more information see :class:`BirmanMurakamiWenzlAlgebra`.

    EXAMPLES::

        sage: BMW3 = algebras.BirmanMurakamiWenzl(3, 'G1, G2', names_idempotents='E1, E2')
        sage: BMW3.an_element()
        (l^-1*m^2)*E2 + m^2*G2*E1*G2^-1*G1^-1 + (-m)*G2*E1*G2^-1 + (-m)*G1*G2*E1
         + G2^-1*G1*G2 + (-m)*G1*G2 + m^2*G2*E1 + m^2*G2 + (-m)*o1
        sage: BMW3.<g1, g2, e1, e2> = algebras.BirmanMurakamiWenzl(3)
        sage: e1**2*~g2*g1*e2
        (l^2-l*m-l^2*m^-2+2*l*m^-1+1-l^-1*m-2*m^-2+2*l^-1*m^-1+l^-2-l^-2*m^-2)*e1*g2^-1*g1^-1
    """
    # --------------------------------------------------------------------------
    # Overloading inherited methods
    # --------------------------------------------------------------------------
    @cached_method
    def __invert__(self):
        r"""
        Return inverse of ``self`` (if possible).

        EXAMPLES::

            sage: BMW3.<g1, g2, e1, e2> = algebras.BirmanMurakamiWenzl(3)
            sage: ~g1                         # indirect doctest
            m*e1 + (-1)*g1 + m*o1
            sage: ~(g1*g2)                    # indirect doctest
            m^2*e2 + (-m)*g1*g2*e1*g2^-1 + m^2*g1*g2*e1 + (-m)*g2*e1 + g2*g1
             + (-m)*g2 + m^2*e1 + (-m)*g1 + m^2*o1
            sage: c = BMW3.base_ring().an_element()
            sage: ~(c*g1)                     # indirect doctest
            (l^-1*m)*e1 + (-l^-1)*g1 + (l^-1*m)*o1
            sage: ~(g1 + g2)                  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: only braid images can be inverted
        """
        b = self.braid_group_algebra_preimage()

        if b is None or len(b.support()) > 1:
            raise NotImplementedError('only braid images can be inverted')

        (br, coeff), = list(b._monomial_coefficients.items())
        P = self.parent()
        return ~coeff * P(~br)

    @cached_method
    def braid_group_algebra_preimage(self):
        r"""
        Return a pre image of ``self`` in the group algebra of the braid group.

        OUTPUT:

        The pre image of ``self`` as an element of the group algebra of the
        braid group.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: ele = BMW3.an_element(); ele
            (l^-1*m^2)*e1 + m^2*g1*e0*g1^-1*g0^-1 + (-m)*g1*e0*g1^-1 + (-m)*g0*g1*e0
             + g1^-1*g0*g1 + (-m)*g0*g1 + m^2*g1*e0 + m^2*g1 + (-m)*o1
            sage: b_ele = ele.braid_group_algebra_preimage(); b_ele
            (-l^-1*m^2) + (l^-1*m)*s1 + m*s1*s0 + m*s1*s0^-1*s1^-1*s0^-1
             + (-1)*s1*s0^-1*s1^-1 + m*s1*s0^-1 + (-m^2)*s0^-1 + m*s1*s0*s1^-1*s0^-1
             + (-1)*s1*s0*s1^-1 + (l^-1*m)*s1^-1 + (-1)*s0*s1*s0
            sage: ele in BMW3
            True
            sage: b_ele in BMW3
            False
            sage: b_ele in BMW3.braid_group_algebra()
            True
        """
        bmw_algebra = self.parent()
        basis = bmw_algebra.basis().keys()
        braid_group_algebra = bmw_algebra.braid_group_algebra()
        braid_group = bmw_algebra.braid_group()
        n = braid_group.strands()
        skn1, skn2, skn3 = bmw_algebra._skein_normalization
        l, m = braid_group_algebra.base_ring().gens()

        def phi(bas_ele):
            tangle = basis[bas_ele]
            w = tangle.defining_word()
            we = [i for i in w if i >= n]
            if we:
                wt = w[:w.index(we[0])]
                wb = w[w.index(we[-1]) + 1:]
                bt = braid_group_algebra(braid_group(wt))
                bb = braid_group_algebra(braid_group(wb))
                if len(we) == 1:
                    i = we[0] - n + 1
                    o = braid_group_algebra.one()
                    bep = braid_group_algebra(braid_group((i,)))
                    ben = braid_group_algebra(braid_group((-i,)))
                    # g_i + skn1*~g_i = m*(1 + skn2*e_i)
                    # => skn2*e_i = ~m*(g_i + skn1*~g_i) - 1
                    # => e_i = skn2*(~m*(g_i + skn1*~g_i) - 1)
                    return skn2*bt*(~m*bep + skn1*~m*ben - o)*bb
                we1 = (we[0],)
                we2 = tuple(we[1:])
                KT = tangle.parent()
                be1 = KT(we1).connector()[0]
                be2 = KT(we2).connector()[0]
                return bt * phi(be1) * phi(be2) * bb
            # in this case wb must be empty, too
            return braid_group_algebra(braid_group(w))

        return bmw_algebra._apply_module_morphism(self, phi,
                                                  codomain=braid_group_algebra)

    @cached_method
    def markov_trace(self):
        r"""
        Return the value of the Markov trace of ``self``.

        If ``self`` is an image of a braid then this methods yield the polynomial
        invariant under regular isotopy of the braid´s closure which is used in the
        definition of the Kauffman polynomial.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: b = BMW3(KnotInfo.K4_1.braid())
            sage: b.markov_trace()
            l^2*m^2 + l*m^3 - l^2 - l*m + 2*m^2 + l^-1*m^3 - 1 - l^-1*m + l^-2*m^2 - l^-2

        REFERENCES::

        - :wikipedia:`Kauffman_polynomial`
        """
        bmw_algebra = self.parent()
        basis = bmw_algebra.basis().keys()
        base_ring = bmw_algebra.base_ring()
        l, m = base_ring.gens()
        x = bmw_algebra._delta
        skn1, skn2, skn3 = bmw_algebra._skein_normalization

        def mtr(bas_ele):
            tangle = basis[bas_ele]
            w = tangle.writhe(closure=True)
            los = tangle.list_of_strands()
            locs = {tuple(st.closure()) for st in los}
            num_closures = len(locs)
            return l**(-skn3*w) * x**(num_closures-1)

        return bmw_algebra._apply_module_morphism(self, mtr,
                                                  codomain=base_ring)

    @cached_method
    def to_iwahori_hecke_algebra(self):
        r"""
        Return the image of ``self`` in the Iwahori-Hecke algebra of Cartan type
        ``A`` under the natural epimorphism.

        OUTPUT:

        The image of ``self`` as instance of the element class of
        :class:`~sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra`.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: b = BMW3.an_element(); b
            (l^-1*m^2)*e1 + m^2*g1*e0*g1^-1*g0^-1 + (-m)*g1*e0*g1^-1 + (-m)*g0*g1*e0
             + g1^-1*g0*g1 + (-m)*g0*g1 + m^2*g1*e0 + m^2*g1 + (-m)*o1
            sage: h = b.to_iwahori_hecke_algebra(); h
            -T[1,2,1] + (q^-2+2+q^2)*T[2] - (q^-1+q)
            sage: h.parent()
            Iwahori-Hecke algebra of type A3 in q^-1,q over Univariate Laurent
             Polynomial Ring in q over Integer Ring in the T-basis
            sage: b2 = b*b
            sage: h2 = b2.to_iwahori_hecke_algebra(); h2
            -(q^-3+2*q^-1+2*q+q^3)*T[1,2,1] + (q^-1+q)*T[1]
             + (q^-5+3*q^-3+5*q^-1+5*q+3*q^3+q^5)*T[2] - (q^-4+3*q^-2+5+3*q^2+q^4)
            sage: h2 == h*h
            True
            sage: BMW3 = algebras.BirmanMurakamiWenzl(3, skein_normalization=(-1, -1, 1))
            sage: b = BMW3.an_element()
            sage: h = b.to_iwahori_hecke_algebra(); h
            T[1,2,1] - (q^-2-2+q^2)*T[2] - (q^-1-q)
            sage: h.parent()
            Iwahori-Hecke algebra of type A3 in q^-1,-q over Univariate Laurent
             Polynomial Ring in q over Integer Ring in the T-basis
            sage: b2 = b*b
            sage: h2 = b2.to_iwahori_hecke_algebra(); h2
            -(q^-3-2*q^-1+2*q-q^3)*T[1,2,1] + (q^-1-q)*T[1]
             + (q^-5-3*q^-3+5*q^-1-5*q+3*q^3-q^5)*T[2] + (q^-4-3*q^-2+5-3*q^2+q^4)
            sage: h2 == h*h
            True
            sage: BMW3 = algebras.BirmanMurakamiWenzl(3, skein_normalization=(-1, 1, -1))
            sage: b = BMW3.an_element()
            sage: h = b.to_iwahori_hecke_algebra(); h
            T[1,2,1] - (q^-2-2+q^2)*T[2] + (q^-1-q)
            sage: h.parent()
            Iwahori-Hecke algebra of type A3 in q,-q^-1 over Univariate Laurent
             Polynomial Ring in q over Integer Ring in the T-basis
            sage: b2 = b*b
            sage: h2 = b2.to_iwahori_hecke_algebra(); h2
            (q^-3-2*q^-1+2*q-q^3)*T[1,2,1] - (q^-1-q)*T[1]
             - (q^-5-3*q^-3+5*q^-1-5*q+3*q^3-q^5)*T[2] + (q^-4-3*q^-2+5-3*q^2+q^4)
            sage: h2 == h*h
            True
        """
        from sage.algebras.iwahori_hecke_algebra import IwahoriHeckeAlgebra
        from sage.functions.generalized import sgn
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        bmw_algebra = self.parent()
        basis = bmw_algebra.basis().keys()
        skn1, skn2, skn3 = bmw_algebra._skein_normalization
        n = bmw_algebra.strands()
        S = LaurentPolynomialRing(ZZ, 'q')
        q = S.gen(0)
        FS = S.fraction_field()
        fq = FS(q)
        T = IwahoriHeckeAlgebra(['A', n], q**(-skn3), skn1*q**skn3).T()
        Ti = T.algebra_generators()
        R = bmw_algebra.base_ring()
        rho = R.hom((fq**(-skn3), fq**(-skn3) + skn1*fq**skn3))

        def phi(bas_ele):
            tangle = basis[bas_ele]
            w = tangle.defining_word()
            if not len(w):
                return T.one()
            if max(w) >= n:
                return T.zero()
            return T.prod(Ti[abs(i)]**sgn(i) for i in w)

        return T.linear_combination((phi(k), S(rho(v))) for k, v in dict(self).items())

    @cached_method
    def to_brauer_algebra(self):
        r"""
        Return the image of ``self`` in the Brauer algebra under the module
        homomorphism according to the shared basis index set.

        .. NOTE::

            This is not an algebra homomorphism.

        OUTPUT:

        The image of ``self`` as instance of the element class of the Brauer
        algebra of the :class:`BrauerAlgebra`.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: ele = BMW3.an_element(); ele
            (l^-1*m^2)*e1 + m^2*g1*e0*g1^-1*g0^-1 + (-m)*g1*e0*g1^-1 + (-m)*g0*g1*e0
             + g1^-1*g0*g1 + (-m)*g0*g1 + m^2*g1*e0 + m^2*g1 + (-m)*o1
            sage: br = ele.to_brauer_algebra(); br
            (l^-1*m^2)*B{{-3, -2}, {-1, 1}, {2, 3}} + m^2*B{{-3, -2}, {-1, 2}, {1, 3}}
             + (-m)*B{{-3, -1}, {-2, 2}, {1, 3}} + (-m)*B{{-3, 1}, {-2, -1}, {2, 3}}
             + B{{-3, 1}, {-2, 2}, {-1, 3}} + (-m)*B{{-3, 1}, {-2, 3}, {-1, 2}}
             + m^2*B{{-3, 2}, {-2, -1}, {1, 3}} + m^2*B{{-3, 2}, {-2, 3}, {-1, 1}}
             + (-m)*B{{-3, 3}, {-2, 2}, {-1, 1}}
        """
        bmw_algebra = self.parent()
        brauer_algebra = bmw_algebra.brauer_algebra()

        return bmw_algebra._apply_module_morphism(self, brauer_algebra.__call__, codomain=brauer_algebra)


class BirmanMurakamiWenzlAlgebra(CombinatorialFreeModule):
    r"""
    The Birman-Murakami-Wenzl algebra.

    INPUT:

    - ``nstrands`` -- a positive integer greater than 1, which indicates the
      number of (unclosed) strands that a monomial, considered as a tangle,
      possesses

    - ``names`` -- a string or tuple of strings containing the names of the
      braid generators, if the lenght of the tuple (or the number of
      comma-separated parts of the string) equals ``nstrands - 1``.
      Alternatively ``names`` may also contain the names of the cap-cup
      generators, if the number in sum is ``2 * nstrands - 2``.
      If this keyword is not specified, the braid generators are named
      ``g_0 ... g_{nstrands - 1}``. If the keyword is given as a single string,
      all braid generators are prefixed with that string and suffixed with
      integers from ``0`` to ``nstrands -2``

    - ``names_idempotents`` -- a string or tuple of strings containing the names
      of the cap-cup generators. The lenght of the tuple (or the number of
      comma-separated parts of the string) must equal ``nstrands - 1`` or ``1``.
      In the latter case all cap-cup generators are prefixed with the given string
      and suffixed with integers from ``0`` to ``nstrands -2``. If this keyword
      is not specified, all cap-cup generators are prefixed with ``e_``.

    - ``params`` -- a string or tuple of strings containing the names of the
      base ring generators. The lenght of the tuple (or the number of
      comma-separated parts of the string) must equal ``2``. If this keyword
      is not specified, these names are set to ``l, m``.

    - ``skein_normalization`` -- a triple of signs (specified as the integers
      ``1`` and ``-1``). This allows switching to other sign conventions in the
      skein relation and the curl-relation. By default, all three signs are
      positive, which corresponds to the convention according to Wikipedia or
      the relations shown in the module header. There you can find some examples
      for this keyword. If the three signs are called `skn1, skn2` and `skn3`
      their meaning can be derived from the following modified relations:

      .. MATH::

        \begin{aligned}
        G_{i} + skn1*{G_{i}}^{-1} & = m(1 + skn2*E_{i}),\\
        G_{i}E_{i} = E_{i}G_{i} & = l^{-skn3}E_{i},\\
        E_{i}G_{i\pm 1}E_{i} & = l^{skn3}E_{i}.
        \end{aligned}

    EXAMPLES::

        sage: BMW2 = algebras.BirmanMurakamiWenzl(2, 'G', 'E', params='Q, Z'); BMW2
        Birman-Murakami-Wenzl algebra on 2 strands
         over Multivariate Laurent Polynomial Ring in Q, Z over Integer Ring
        sage: BMW2.gens()
        (G, E)
        sage: BMW4 = algebras.BirmanMurakamiWenzl(4)
        sage: BMW4.gens()
        (g0, g1, g2, e0, e1, e2)
        sage: BMW4.base_ring()
        Multivariate Laurent Polynomial Ring in l, m over Integer Ring
        sage: BMW3 = BMW4.birman_murakami_wenzl_subalgebra()
        sage: BMW3.gens()
        (g0, g1, e0, e1)

    Element construction::

        sage: tupt = (1, -2, 3, -1, 4)
        sage: bt = BMW3(tupt); bt
        (l*m^2-l)*e1 + (-l*m)*g1*e0*g1^-1*g0^-1 + (l*m^2+m)*e0*g1^-1*g0^-1
        sage: bt4 = BMW4(bt); bt4
        (l*m^2-l)*e1 + (l*m^2+m)*e0*g1^-1*g0^-1 + (-l*m)*g1*e0*g1^-1*g0^-1
        sage: bt4.parent().strands()
        4
        sage: bt.parent().strands()
        3
        sage: tangle = BMW3.tangle_semigroup()(tupt); tangle
        g0*g1^-1*e0*g0^-1*e1
        sage: bt == BMW3(tangle)
        True
        sage: pre = bt.braid_group_algebra_preimage(); pre
        (-l*m^2+l) + (l*m-l*m^-1)*s1 + (l*m+1)*s0^-1*s1^-1*s0^-1
         + (-l)*s1*s0^-1*s1^-1*s0^-1 + l*m*s0^-1 + (-l)*s1*s0*s1^-1*s0^-1
         + (-l*m^2-m)*s1^-1*s0^-1 + (l*m+1)*s0*s1^-1*s0^-1 + (l*m-l*m^-1)*s1^-1
        sage: bt == BMW3(pre)
        True

        sage: tupb = (1, -2, 1, -2)
        sage: bb = BMW3(tupb); bb
        (l*m^3+m^4-m^2)*e1 + (-l*m^2-m^3)*g1*e0*g1^-1*g0^-1
         + (l*m^3+m^4)*e0*g1^-1*g0^-1 + m^3*g0*g1*e0*g1^-1
         + (-m^2)*g1*e0*g1^-1 + (m^3+l^-1*m^2)*e0*g1^-1
         + (-m^2)*g0*g1*e0 + (-m^2)*g0*g1 + m*g1*e0
         + (-1)*g1*g0 + m*g1 + (-m^2)*e0 + m^3*g0 + (-m^2)*o1
        sage: braid = BMW3.braid_group()(tupb); braid
        (s0*s1^-1)^2
        sage: bb == BMW3(braid)
        True

        sage: c = BMW3.cubic_hecke_algebra().an_element(); c
        (-l^-1)*c0*c1^-1 + (1+l^-1*m)*c0 + (m+l^-1)*c1 + (l*m-l^-1*m)
        sage: BMW3(c)
        (-l^-1*m^2)*e1 + (l^-1*m)*g1*e0*g1^-1*g0^-1
         + (-l^-1*m^2)*e0*g1^-1*g0^-1 + (l^-1)*g0*g1
         + (m+l^-1)*g1 + g0 + (l*m-l^-1*m)*o1

        sage: BR = BMW4.brauer_algebra()
        sage: br = BR.an_element(); br
        2*B{{-4, -3}, {-2, -1}, {1, 2}, {3, 4}}
         + 2*B{{-4, -3}, {-2, -1}, {1, 3}, {2, 4}}
         + 3*B{{-4, -3}, {-2, -1}, {1, 4}, {2, 3}}
        sage: BMW4(br)
        2*e0*e2 + 2*g1^-1*e0*e2 + 3*g2^-1*g1^-1*e0*e2
    """
    Element = BirmanMurakamiWenzlElement

    @staticmethod
    def __classcall_private__(
        cls, nstrands: int,
        names: str | tuple = 'g',
        names_idempotents: str | tuple | None = None,
        params: str | tuple = 'l, m',
        skein_normalization: tuple = (1, 1, 1)
    ):
        r"""
        Normalize input to ensure a unique representation.

        TESTS::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: TestSuite(BMW3).run()
        """
        ng_gens = nstrands - 1
        ne_gens = ng_gens
        if type(names) is str:
            nn = len(names.split(','))
        else:
            nn = len(names)

        from sage.structure.category_object import normalize_names
        if nn == 2*ng_gens:
            names = tuple(normalize_names(nn, names))
        elif nn not in (1, ng_gens):
            raise ValueError('there must be %s names' % ng_gens)
        else:
            if names_idempotents is None:
                names_idempotents = 'e'
            if type(names_idempotents) is str:
                ni = len(names_idempotents.split(','))
            else:
                ni = len(names_idempotents)
            if ni not in (1, ne_gens):
                raise ValueError('there must be %s names_idempotents' % ne_gens)
            names = tuple(list(normalize_names(ng_gens, names)) + list(normalize_names(ne_gens, names_idempotents)))

        params = tuple(normalize_names(2, params))
        return super().__classcall__(cls, nstrands, names, params=params, skein_normalization=skein_normalization)

    def __init__(
        self,
        nstrands: int,
        names: str | tuple,
        params: str | tuple = 'l, m',
        skein_normalization: tuple = (1, 1, 1),
    ):
        r"""
        Initialize ``self``.

        TESTS::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2, 'G', 'E', params='Q, X'); BMW2
            Birman-Murakami-Wenzl algebra on 2 strands
             over Multivariate Laurent Polynomial Ring in Q, X over Integer Ring
            sage: TestSuite(BMW2).run()
        """
        # ----------------------------------------------------------------------
        # Define underlying monoid of tangles
        # ----------------------------------------------------------------------
        self._nstrands = int(nstrands)
        KT = KauffmanTangles(names)
        self._tangles = KT
        from sage.groups.braid import BraidGroup
        self._braid_group = BraidGroup(self._nstrands)

        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        R = LaurentPolynomialRing(ZZ, params)
        l, m = R.gens()
        self._skein_normalization = skein_normalization
        skn1, skn2, skn3 = skein_normalization
        self._delta = ~m*(l**skn3 + skn1*l**(-skn3)) - skn2
        verbose('base_ring %s defined' % R, level=2)

        # ----------------------------------------------------------------------
        # defining associated algebras
        # ----------------------------------------------------------------------
        from sage.combinat.diagram_algebras import BrauerAlgebra
        self._brauer_algebra = BrauerAlgebra(self._nstrands, self._delta)
        from sage.algebras.group_algebra import GroupAlgebra
        self._braid_group_algebra = GroupAlgebra(self._braid_group, R=R)

        # ----------------------------------------------------------------------
        # Setup a Morton-Wasserman basis of self. This means that we use
        # the basis of the Brauer algebra viewing its elements as simple layered
        # tangles as described in [EG2017]_ and map them to elements in
        # the semigroup of self._tangles with the help of
        # :meth:`morton_wasserman_tangle` of :class:`KauffmanTangles`.
        # ----------------------------------------------------------------------
        from sage.sets.family import Family
        brauer_basis = self._brauer_algebra.basis().keys()
        basis = Family(brauer_basis, function=KT.morton_wasserman_tangle, lazy=True, name='Morton-Wasserman-Tangle')

        # ----------------------------------------------------------------------
        # defining the algebra itself
        # ----------------------------------------------------------------------
        from sage.categories.finite_dimensional_algebras_with_basis import (
            FiniteDimensionalAlgebrasWithBasis,
        )
        category = FiniteDimensionalAlgebrasWithBasis(R)

        CombinatorialFreeModule.__init__(self, R, basis, prefix='', bracket=False, category=category)
        self.print_options()['names'] = KT._mwt_names

        # ----------------------------------------------------------------------
        # init the attributes being set on demand
        # ----------------------------------------------------------------------
        self._birman_murakami_wenzl_subalgebra = None
        self._from_kauffman_tangle_cache = {}

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT: string describing ``self``

        TESTS::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: BMW3 # indirect doctest
            Birman-Murakami-Wenzl algebra on 3 strands
             over Multivariate Laurent Polynomial Ring in l, m over Integer Ring
        """
        s = 'Birman-Murakami-Wenzl algebra on %s strands over %s'
        return s % (self.strands(), self.base_ring())

    def _element_constructor_(self, x) -> BirmanMurakamiWenzlElement:
        r"""
        Extensions to the element constructor of class :class:`CombinatorialFreeModule`.

        New functionalities are:

        - constructing element from a tangle (semigroup homomorphism)
        - constructing element from a braid (group homomorphism)
        - constructing element from a tangle or braid giving in Tietze form
        - constructing element from an element of the braid group algebra
          (algebra homomorphism)
        - constructing element from an element of the cubic Hecke algebra
          (algebra homomorphism)
        - constructing element from an element of the Brauer algebra
          (module homomorphism)
        - constructing element from an element of an other Birman-Murakami-Wenzl
          algebra over an other base ring or with less strands

        INPUT:

        - ``x`` -- can be one of the following:

          * a tuple of integers interpreted as a word of a tangle or braid
          * an instance of the element class of ``self.tangle_semigroup()``
          * an instance of the element class of ``self.braid_group()``
          * an instance of the element class of ``self.braid_group_algebra()``
          * an instance of the element class of ``self.brauer_algebra()``
          * an instance of the element class of ``self.cubic_hecke_algebra()``
          * an instance of the element class of ``self`` (but possible
            to a different parent)
          * any other object which works for the element constructor
            of :class:`CombinatorialFreeModule`

        EXAMPLES::

            sage: BMW2.<g, e> = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2((1, 2, -1, 2))
            (l*m^-1-1+l^-1*m^-1)*e
        """
        braid_grp_alg = self.braid_group_algebra()
        n = self.strands()
        BA = self.brauer_algebra()
        BD = self.basis().function.domain().keys()
        T = self._tangles

        xb = x
        if isinstance(x, (tuple, list)):
            x = tuple(x)
            xb = T(x)

        from sage.groups.braid import Braid
        if isinstance(x, Braid) and x.strands() == n:
            xb = T(x)

        if isinstance(xb, T.element_class) and xb.strands() == n:
            return self._from_kauffman_tangle(xb)

        if xb in BD:
            return self.monomial(xb)

        def mwt_mon(bas_ele): return self.monomial(T.morton_wasserman_tangle(bas_ele))

        if xb in BA:
            return BA._apply_module_morphism(xb, lambda ele: mwt_mon(ele), codomain=self)

        if isinstance(xb, BirmanMurakamiWenzlElement):
            other_bmw = xb.parent()
            other_base_ring = other_bmw.base_ring()
            on = other_bmw.strands()
            if other_base_ring != self.base_ring():
                if on == n:
                    xbv = xb.to_vector()
                    img_xbv = vector([self.base_ring()(cf) for cf in xbv])
                    return self.from_vector(img_xbv)
            elif on < n:
                def fd(ele): return mwt_mon(BD([(-n, n)] + list(ele)))
                return other_bmw._apply_module_morphism(xb, fd, codomain=self)

        if xb in braid_grp_alg:
            return braid_grp_alg._apply_module_morphism(xb, self, codomain=self)

        from sage.algebras.hecke_algebras.cubic_hecke_algebra import CubicHeckeElement
        if isinstance(xb, CubicHeckeElement):
            CHA = self.cubic_hecke_algebra()
            if xb in CHA:
                def fc(ele): return self(ele.braid())
                return CHA._apply_module_morphism(xb, fc, codomain=self)

        result = CombinatorialFreeModule._element_constructor_(self, x)
        return result

    def ngens(self) -> int:
        r"""
        The number of generators of the algebra.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.ngens()
            2
        """
        return 2*(self.strands() - 1)

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.algebra_generators()
            Finite family {{{-2, 1}, {-1, 2}}: g, {{-2, -1}, {1, 2}}: e}
        """
        T = self._tangles
        d = {T(g).connector()[0]: self(T(g)) for g in T.ambient().gens()}
        from sage.sets.family import Family
        return Family(list(d), d.__getitem__)

    def gens(self) -> tuple:
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.gens()
            (g, e)
        """
        return tuple(self.algebra_generators())

    def gen(self, i) -> BirmanMurakamiWenzlElement:
        r"""
        The ``i``-th generator of the algebra.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.gen(0), BMW2.gen(1)
            (g, e)
        """
        i = int(i)
        n = 2 * self.strands() - 2
        if i < 0 or i >= n:
            raise IndexError('i must be non negative and less than %s' % n)
        return self.gens()[i]

    def g(self, i: int) -> BirmanMurakamiWenzlElement:
        r"""
        The ``i``-th g-generator of the algebra.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.g(1)
            g
        """
        i = int(i)
        n = self.strands()
        if i < 1 or i >= n:
            raise IndexError('i must be positive and less than %s' % n)
        return self.gen(i - 1)

    def e(self, i: int) -> BirmanMurakamiWenzlElement:
        r"""
        The ``i``-th e-generator of the algebra.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.e(1)
            e
        """
        i = int(i)
        n = self.strands()
        if i < 1 or i >= n:
            raise IndexError('i must be positive and less than %s' % n)
        return self.gen(i + self.strands() - 2)

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the basis element for the identity element

         EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.one_basis()
            {{-2, 2}, {-1, 1}}
        """
        BD = self.basis().function.domain().keys()
        return BD.from_involution_permutation_triple(([], [], list(range(1, self.strands() + 1))))

    def _an_element_(self) -> BirmanMurakamiWenzlElement:
        r"""
        Overwrite the original method from :mod:`~sage.combinat.free_module`
        to obtain a more interesting element for ``TestSuite``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: ele = BMW2.an_element(); ele
            (l^-1*m)*e + m*g + (-1)*o1
            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: ele = BMW3.an_element(); ele
            (l^-1*m^2)*e1 + m^2*g1*e0*g1^-1*g0^-1 + (-m)*g1*e0*g1^-1
             + (-m)*g0*g1*e0 + g1^-1*g0*g1 + (-m)*g0*g1 + m^2*g1*e0
             + m^2*g1 + (-m)*o1
            sage: BMW4 = algebras.BirmanMurakamiWenzl(4)
            sage: ele = BMW4.an_element(); ele
            m^3*e1*g2^-1*g1^-1 + (-m^2)*g1*e0*g1^-1*g2*g0^-1*g1
             + (-m^2)*e1*g2^-1 + m*g1*e0*g1^-1*g2^-1*g0^-1
             + m^2*e0*g1^-1*g2*g0^-1 + (-1)*g0*g1*g2
             + m*g2*g0 + m^3*e1 + (-m^2)*g1*e0*g1^-1*g0^-1
        """
        n = self.strands()
        gens = self.gens()
        gs = [gen for gen in gens if gens.index(gen) < n - 1]
        es = [gen for gen in gens if gens.index(gen) >= n - 1]
        if n == 2:
            g1, = gs
            e1, = es
            return g1**2
        if n == 3:
            g1, g2 = gs
            e1, e2 = es
            return g2*~g1*g2
        g1, g2, g3 = gs
        e1, e2, e3 = es
        return g1*~g2*g3

    def _from_kauffman_tangle(self, tangle: KauffmanTangle, rec_count: int = 0) -> BirmanMurakamiWenzlElement:
        r"""
        Return an element of ``self`` constructed from the given instance of
        :class:`KauffmanTangle`.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: T = BMW3.tangle_semigroup()
            sage: t = T(KnotInfo.K6_2.braid()); t
            g0^2*(g0*g1^-1)^2
            sage: BMW3._from_kauffman_tangle(t)
            (l*m^5+m^6-2*l*m^3-3*m^4+m^2)*e1 + (-l*m^4-m^5+l*m^2+2*m^3)*g1*e0*g1^-1*g0^-1
             + (l*m^5+m^6-l*m^3-m^4+l^-1*m^3-l^-1*m)*e0*g1^-1*g0^-1
             + (m^5-2*m^3)*g0*g1*e0*g1^-1 + (-m^4+m^2)*g1*e0*g1^-1
             + (m^5-m^3+2*l^-1*m^4-l^-1*m^2+2*l^-2*m^3+l^-3*m^2)*e0*g1^-1
             + (-m^4+m^2)*g0*g1*e0 + m*g1^-1*g0*g1 + (-m^4+m^2)*g0*g1
             + (m^3-m)*g1*e0 + g1*g0 + (m^3-m)*g1 + (-m^4+m^2-l^-1*m^3-l^-2*m^2)*e0
             + (m^5-2*m^3)*g0 + (-m^4+m^2)*o1
        """
        cache = self._from_kauffman_tangle_cache
        if tangle in cache:
            return cache[tangle]
        res = self._compute_from_kauffman_tangle(tangle, rec_count)
        cache[tangle] = res
        return res

    def _compute_from_kauffman_tangle(self, tangle: KauffmanTangle, rec_count: int) -> BirmanMurakamiWenzlElement:
        r"""
        Worker for :meth:`_from_kauffman_tangle`.

        This does the actual (uncached) work; go through
        :meth:`_from_kauffman_tangle` so that the result is memoized.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: T = BMW3.tangle_semigroup()
            sage: t = T(KnotInfo.K6_2.braid())
            sage: BMW3._compute_from_kauffman_tangle(t, 0) == BMW3._from_kauffman_tangle(t)
            True
        """
        base_ring = self.base_ring()
        l, m = base_ring.gens()
        x = self._delta
        skn1, skn2, skn3 = self._skein_normalization
        # only build the (expensive) verbose messages when they would be shown;
        # otherwise the eager ``%``-formatting of tangles/elements dominates
        dbg = get_verbose() >= 2

        prefix = ' '*2*rec_count + repr(rec_count) if dbg else ''
        pos = tangle.find_unlayered_crossing()
        if pos is None:
            conn, loops = tangle.connector()
            # obtain the Morton Wasserman tangle (needed because of lazy basis family)
            mwt = self.basis().keys()[conn]
            writhe = tangle.writhe()
            if dbg:
                verbose('%s Tangle %s with %s loops and writhe %s is layered and isotopic to %s' % (prefix, tangle, loops, writhe, mwt), level=2)
            return l**(-skn3*writhe)*x**loops*self(conn)

        # if there are unlayered crossings we resolve them recursively
        w = tangle.defining_word()
        P = tangle.parent()
        tangl = P(w[:pos])
        tangr = P(w[pos + 1:])

        i = w[pos]
        g = P((-i,))
        e = P((abs(i) + self.strands() - 1,))
        if dbg:
            prompt = '%s Tangle %s is not layered at position %s' % (prefix, tangle, pos)
            verbose('%s (left %s, right %s, g %s, e %s): starting recursion' % (prompt, tangl, tangr, g, e), level=2)
            verbose('%s, elem_g start recursion' % prompt, level=2)
        elem_g = self._from_kauffman_tangle(tangl * g * tangr, rec_count + 1)
        if dbg:
            verbose('%s, elem_g: %s end recursion' % (prompt, elem_g), level=2)
            verbose('%s, elem_e start recursion' % prompt, level=2)
        elem_e = self._from_kauffman_tangle(tangl * e * tangr, rec_count + 1)
        if dbg:
            verbose('%s, elem_e: %s end recursion' % (prompt, elem_e), level=2)
            verbose('%s, elem_0 start recursion' % prompt, level=2)
        elem_0 = self._from_kauffman_tangle(tangl * tangr, rec_count + 1)
        if dbg:
            verbose('%s, elem_0: %s end recursion' % (prompt, elem_0), level=2)

        # since elem_g has less unlayered crossings and  since elem_e and
        # elem_0 have less crossings at all the recursion must terminate
        # we use g_i + skn1*~g_i = m*(1 + skn2*e_i)
        if i > 0:
            res = -skn1*elem_g + m*(elem_0 + skn2*elem_e)
        else:
            res = -skn1*elem_g + skn1*m*(elem_0 + skn2*elem_e)
        if dbg:
            verbose('%s, result: %s' % (prompt, res), level=2)
        return res

    @cached_method
    def product_on_basis(self, g1, g2) -> BirmanMurakamiWenzlElement:
        r"""
        Return the product of the basis elements indexed by ``g1`` and ``g2``.

        EXAMPLES::

            sage: BMW3.<g1, g2, e1, e2> = algebras.BirmanMurakamiWenzl(3)
            sage: g = BMW3.basis().keys().keys().an_element(); g
            {{-3, 3}, {-2, -1}, {1, 2}}
            sage: gg = BMW3.product_on_basis(g, g); gg
            (l*m^-1-1+l^-1*m^-1)*e1
            sage: gg == e1**2
            True
        """
        # ----------------------------------------------------------------------
        # short way for multiplications with one
        # ----------------------------------------------------------------------
        if g1 == self.one_basis():
            return self.monomial(g2)

        if g2 == self.one_basis():
            return self.monomial(g1)

        bas = self.basis().keys()
        t1 = bas[g1]
        t2 = bas[g2]
        verbose('calculating %s * %s:' % (t1, t2), level=2)
        return self._from_kauffman_tangle(t1 * t2)

    def strands(self) -> int:
        r"""
        Return the number of (unclosed) strands a monomial considered as a
        tangle has.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.strands()
            2
        """
        return self._nstrands

    def tangle_semigroup(self):
        r"""
        Return the semigroup of tangles attached to ``self``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.tangle_semigroup()
            Semigroup of tangles with 2 (non closed) strands with generators Family (1, g, e, g^-1)
        """
        return self._tangles

    def braid_group(self):
        r"""
        Return the braid group attached to ``self``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.braid_group()
            Braid group on 2 strands
        """
        return self._braid_group

    def brauer_algebra(self):
        r"""
        Return the Brauer algebra attached to ``self`` over the
        base ring of ``self``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.brauer_algebra()
            Brauer Algebra of rank 2 with parameter l*m^-1 - 1 + l^-1*m^-1
             over Multivariate Laurent Polynomial Ring in l, m over Integer Ring
        """
        return self._brauer_algebra

    def braid_group_algebra(self):
        r"""
        Return the group algebra of braid group attached to ``self`` over the
        base ring of ``self``.

        EXAMPLES::

            sage: BMW2 = algebras.BirmanMurakamiWenzl(2)
            sage: BMW2.braid_group_algebra()
            Algebra of Braid group on 2 strands
             over Multivariate Laurent Polynomial Ring in l, m over Integer Ring
        """
        return self._braid_group_algebra

    @cached_method
    def cubic_hecke_algebra(self, extension_param: str = 'q'):
        r"""
        Return the cubic Hecke algebra attached to ``self`` over the
        base ring of ``self``.

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3)
            sage: CHA = BMW3.cubic_hecke_algebra()
            sage: cub_equ = CHA.cubic_equation(); cub_equ
            h^3 + (-m - l^-1)*h^2 + (1 + l^-1*m)*h - l^-1
            sage: roots = CHA.cubic_equation_roots(); roots
            [a, b, -b - a + m + l^-1]
            sage: r1, r2, r3 = roots
            sage: x = polygen(r1)
            sage: cub_equ_2 = (x - r1)*(x - r2)*(x - r3)
            sage: cub_equ.coefficients() == cub_equ_2.coefficients()
            True

        Compare Markov traces::

            sage: b = KnotInfo.K6_2.braid()
            sage: c = CHA(b)
            sage: m = c.formal_markov_trace(); m
            (-l^3*m-l^2*m^2-2*l*m-2*m^2-l^-1*m^3-l^-1*m-l^-2*m^2+l^-3*m+l^-4)*B[U1]
             + (l^2*m^2+l*m^3+l*m+2*m^2+l^-1*m)*B[U2] + (m^2-1+l^-1*m+l^-2)*B[K4]
            sage: g0, g1, e0, e1 = BMW3.gens()
            sage: mtU2 = e0.markov_trace()
            sage: mtK4 = BMW3(KnotInfo.K4_1.braid()).markov_trace()
            sage: mcfU1, mcfU2, mcfK4 = m.coefficients()
            sage: mcfU1 + mcfU2*mtU2 + mcfK4*mtK4 == BMW3(b).markov_trace()
            True

        Using other ``skein_normalization``::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3, skein_normalization=(-1, -1, 1))
            sage: CHA = BMW3.cubic_hecke_algebra()
            sage: CHA.cubic_equation()
            h^3 + (-m - l^-1)*h^2 + (-1 + l^-1*m)*h + l^-1
        """
        R = self.base_ring()
        l, m = R.gens()
        # g_i + skn1*~g_i = m*(1 + skn2*e_i)            (skein relation)
        # g_i^2 + skn1 = m*(g_i + skn2*~l**skn3*e_i)    (using twist relation)
        # skn2 * ~l**skn3*e_i =  ~l**skn3*~m(g_i + skn1*~g_1))  - ~l**skn3
        # g_i^3 + skn1*g_i = m*g_i^2 + ~l**skn3*(g_i^2 + skn1) - ~l**skn3*m*g_i
        # g_i^3 - (m + ~l**skn3)*g_i^2 + (~l**skn3*m + skn1)*g_i -skn1*~l**skn3 == 0
        # u = (m + ~l**skn3), v = (~l**skn3*m + skn1), w = skn1*~l**skn3
        skn1, skn2, skn3 = self._skein_normalization
        cubic_equation_parameters = (m + ~l**skn3, ~l**skn3*m + skn1, ~l**skn3*skn1)
        from sage.algebras.hecke_algebras.cubic_hecke_algebra import CubicHeckeAlgebra
        CHA = CubicHeckeAlgebra(self.strands(),
                                cubic_equation_parameters=cubic_equation_parameters,
                                warning=False)
        # induce a map from the base ring of the Markov trace module to R
        MTR = CHA._markov_trace_module().base_ring()
        u, v, w = CHA.cubic_equation_parameters(generic=True)
        MTR.create_specialization([R(u), R(v), R(w)], ~R(w))
        return CHA

    @cached_method
    def birman_murakami_wenzl_subalgebra(self, nstrands=None):
        r"""
        Return a :class:`BirmanMurakamiWenzlAlgebra` that realizes a sub-algebra
        of ``self`` on the first ``n_strands`` strands.

        INPUT:

        - ``nstrands`` -- integer at least 1 and at most :meth:`strands` giving
          the number of strands for the subgroup; the default is one strand
          less than ``self`` has

        OUTPUT: an instance of this class realizing the sub-algebra

        EXAMPLES::

            sage: BMW3 = algebras.BirmanMurakamiWenzl(3, params='a, z')
            sage: BMW3.birman_murakami_wenzl_subalgebra()
            Birman-Murakami-Wenzl algebra on 2 strands
              over Multivariate Laurent Polynomial Ring in a, z over Integer Ring
        """
        n = self.strands()
        if nstrands is None:
            nstrands = n - 1
        nstrands = ZZ(nstrands)

        if nstrands >= n or nstrands <= 0:
            raise ValueError('nstrands must be positive and less than %s' % n)

        names = [str(g) for g in self.gens()]
        names_g = tuple(g for g in names if names.index(g) < nstrands - 1)
        names_e = tuple(g for g in names if names.index(g) < 2*nstrands - 1 and names.index(g) >= n - 1)
        params = self.base_ring()._names
        skn = self._skein_normalization

        SubBMWAlg = BirmanMurakamiWenzlAlgebra(nstrands,
                                               names=names_g,
                                               names_idempotents=names_e,
                                               params=params,
                                               skein_normalization=skn)

        if nstrands == n - 1:
            self._birman_murakami_wenzl_subalgebra = SubBMWAlg
        return SubBMWAlg
