# sage.doctest: needs sage.combinat sage.groups sage.modules
r"""
Cubic Hecke Algebras

We consider the factors of the group algebra of the Artin braid groups
such that the images `s_i` of the braid generators satisfy a cubic equation:

.. MATH::

    s_i^3 = u s_i^2 - v s_i + w.

Here `u, v, w` are elements in an arbitrary integral domain and `i` is a
positive integer less than `n`, the number of the braid group's strands.
By the analogue to the *Iwahori Hecke algebras* (see
:class:`~sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra`), in which the
braid generators satisfy a quadratic relation these algebras have been called
*cubic Hecke algebras*. The relations inherited from the braid group are:

.. MATH::

    s_i s_{i+1} s_i = s_{i+1} s_i s_{i+1} \text{ for } 1 \leq i < n - 1
    \mbox{ and } s_i s_j = s_j s_i \text{ for } 1 \leq i < j - 1 < n - 1.

The algebra epimorphism from the braid group algebra over the same base ring is
realized inside the element constructor of the present class, for example in the
case of the 3 strand cubic Hecke algebra::

    sage: CHA3 = algebras.CubicHecke(3)
    sage: BG3 = CHA3.braid_group()
    sage: braid = BG3((1,2,-1,2,2,-1)); braid
    c0*c1*c0^-1*c1^2*c0^-1
    sage: braid_image = CHA3(braid); braid_image
    u*c1*c0^-1*c1 + u*v*c0*c1^-1*c0^-1 + (-u^2)*c0^-1*c1
    + ((u^2*v-v^2)/w)*c0*c1*c0^-1 + ((u^2-v)/w)*c0*c1*c0
    + ((-u^3+u*v)/w)*c0*c1 + (-u*v+w)*c1^-1

If the ring elements `u, v, w` (which will be called the *cubic equation
parameters* in the sequel) are taken to be `u = v = 0, w = 1` the cubic Hecke
algebra specializes to the group algebra of the *cubic braid group*, which is
the factor group of the Artin braid group under setting the generators order
to be three. These groups can be obtained by
:meth:`CubicHeckeAlgebra.cubic_braid_group`.

It is well known, that these algebras are free of finite rank as long as the
number of braid generators is less than six and infinite dimensional else wise.
In the former (non trivial) cases they are also known as *cyclotomic Hecke
algebras* corresponding to the complex reflection groups having Shepard-Todd
number `4`, `25` and `32`.

Since the *Broué, Malle, Rouquiere* conjecture has been proved (for references
of these cases see [Mar2012]_) there exists a finite free basis of the cubic
Hecke algebra which is in bijection to the cubic braid group and compatible
with the specialization to the cubic braid group algebra as explained above.

For the algebras corresponding to braid groups of less than five strands such
a basis has been calculated by Ivan Marin. This one is used here. In the case
of 5 strands such a basis is not available, right now. Instead the elements
of the cubic braid group class themselves are used as basis elements. This
is also the case when the cubic braid group is infinite, even though it is
not known if these elements span all of the cubic Hecke algebra.

Accordingly, be aware that the module embedding of the group algebra of the
cubicbraid groups is known to be an isomorphism of free modules only in the
cases of less than five strands.

EXAMPLES:

Consider the obstruction ``b`` of the *triple quadratic algebra* from Section 2.6
of [Mar2018]_. We verify that the third power of it is a scalar multiple
of itself (explicitly ``2*w^2`` times the *Schur element* of the three
dimensional irreducible representation)::

    sage: CHA3 = algebras.CubicHecke(3)
    sage: c1, c2 = CHA3.gens()
    sage: b = c1^2*c2 - c2*c1^2 - c1*c2^2 + c2^2*c1; b
    w*c0^-1*c1 + (-w)*c0*c1^-1 + (-w)*c1*c0^-1 + w*c1^-1*c0
    sage: b2 = b*b
    sage: b3 = b2*b
    sage: BR = CHA3.base_ring()
    sage: ER = CHA3.extension_ring()
    sage: u, v, w = BR.gens()
    sage: f =  BR(b3.coefficients()[0]/w)
    sage: try:
    ....:     sh = CHA3.schur_element(CHA3.irred_repr.W3_111)
    ....: except NotImplementedError:    # for the case GAP3 / CHEVIE not available
    ....:     sh = ER(f/(2*w^2))
    sage: ER(f/(2*w^2)) == sh
    True
    sage: b3 == f*b
    True

Defining the cubic Hecke algebra on 6 strands will need some seconds for
initializing. However, you can do calculations inside the infinite
algebra as well::

    sage: # optional - database_cubic_hecke
    sage: CHA6 = algebras.CubicHecke(6)
    sage: CHA6.inject_variables()
    Defining c0, c1, c2, c3, c4
    sage: s = c0*c1*c2*c3*c4; s
    c0*c1*c2*c3*c4
    sage: s^2
    (c0*c1*c2*c3*c4)^2
    sage: t = CHA6.an_element() * c4; t
    (-w)*c0*c1^-1*c4 + v*c0*c2^-1*c4 + u*c2*c1*c4 + ((-v*w+u)/w)*c4

REFERENCES:

- [Mar2012]_
- [Mar2018]_
- [CM2012]_

AUTHORS:

- Sebastian Oehms May 2020: initial version
"""
# ###########################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ###########################################################################
from warnings import warn

from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.misc.verbose import verbose
from sage.groups.cubic_braid import CubicBraidGroup
from sage.rings.integer_ring import ZZ
from sage.algebras.splitting_algebra import solve_with_extension
from sage.modules.free_module_element import vector
from sage.matrix.matrix_space import MatrixSpace
from sage.algebras.hecke_algebras.cubic_hecke_base_ring import CubicHeckeRingOfDefinition
from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import CubicHeckeMatrixSpace, AbsIrreducibeRep, RepresentationType


##############################################################################
#
#                  Class CubicHeckeElement (for elements)
#
##############################################################################
class CubicHeckeElement(CombinatorialFreeModule.Element):
    r"""
    An element of a :class:`CubicHeckeAlgebra`.

    For more information see :class:`CubicHeckeAlgebra`.

    EXAMPLES::

        sage: CHA3s = algebras.CubicHecke('s1, s2'); CHA3s.an_element()
        (-w)*s1*s2^-1 + v*s1 + u*s2 + ((-v*w+u)/w)
        sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
        sage: c1**3*~c2
        u*w*c1^-1*c2^-1 + (u^2-v)*c1*c2^-1 + (-u*v+w)*c2^-1
    """
    # --------------------------------------------------------------------------
    # Overloading inherited methods
    # --------------------------------------------------------------------------
    def __invert__(self):
        r"""
        Return inverse of ``self`` (if possible).

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele1 = CHA3((1,-2,1)); ele1
            c0*c1^-1*c0
            sage: ~ele1                       # indirect doctest
            c0^-1*c1*c0^-1

            sage: CHA2 = algebras.CubicHecke(2)
            sage: x = CHA2.an_element(); x
            v*c + ((-v*w+u)/w)
            sage: ~x
            Traceback (most recent call last):
            ...
            ValueError: cannot invert self (= v*c + ((-v*w+u)/w))
        """
        self_Tietze = self.Tietze()

        if self_Tietze is None:
            return super().__invert__()

        inverse_Tietze = ()
        len_self = len(self_Tietze)

        inverse_Tietze = tuple([-1 * self_Tietze[len_self - i - 1]
                                for i in range(len_self)])
        P = self.parent()
        return P(inverse_Tietze)

    def Tietze(self):
        r"""
        Return the Tietze presentation of ``self`` if ``self`` belongs to the
        basis of its parent and ``None`` otherwise.

        OUTPUT:

        A tuple representing the pre image braid of ``self`` if ``self`` is a
        monomial from the basis ``None`` else-wise

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element(); ele
            (-w)*c0*c1^-1 + v*c0 + u*c1 + ((-v*w+u)/w)
            sage: ele.Tietze() is None
            True
            sage: [CHA3(sp).Tietze() for sp in ele.support()]
            [(), (1,), (1, -2), (2,)]
        """
        vecd = self._monomial_coefficients
        if len(vecd) != 1:
            return None
        ind, coeff = next(iter(vecd.items()))
        if coeff.is_one():
            return ind.Tietze()

    def max_len(self):
        r"""
        Return the maximum of the length of Tietze expressions among the
        support of ``self``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element(); ele
            (-w)*c0*c1^-1 + v*c0 + u*c1 + ((-v*w+u)/w)
            sage: ele.max_len()
            2
        """
        return max(len(bas_ele.Tietze()) for bas_ele in self.support())

    def braid_group_algebra_pre_image(self):
        r"""
        Return a pre image of ``self`` in the group algebra of the braid group
        (with respect to the basis given by Ivan Marin).

        OUTPUT:

        The pre image of ``self`` as instance of the element class of the group
        algebra of the BraidGroup

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element(); ele
            (-w)*c0*c1^-1 + v*c0 + u*c1 + ((-v*w+u)/w)
            sage: b_ele = ele.braid_group_algebra_pre_image(); b_ele
            ((-v*w+u)/w) + v*c0 + u*c1 + (-w)*c0*c1^-1
            sage: ele in CHA3
            True
            sage: b_ele in CHA3
            False
            sage: b_ele in CHA3.braid_group_algebra()
            True
        """
        ch_algebra = self.parent()
        braid_group_algebra = ch_algebra.braid_group_algebra()
        braid_group = ch_algebra.braid_group()

        def phi(bas_ele):
            return braid_group_algebra(braid_group(bas_ele))
        return ch_algebra._apply_module_morphism(self, phi,
                                                 codomain=braid_group_algebra)

    def cubic_braid_group_algebra_pre_image(self):
        r"""
        Return a pre image of ``self`` in the group algebra of the cubic braid
        group.

        OUTPUT:

        The pre image of ``self`` as instance of the element class of the group
        algebra of the :class:`CubicBraidGroup`.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element(); ele
            (-w)*c0*c1^-1 + v*c0 + u*c1 + ((-v*w+u)/w)
            sage: cb_ele = ele.cubic_braid_group_algebra_pre_image(); cb_ele
            (-w)*c0*c1^-1 + v*c0 + u*c1 + ((-v*w+u)/w)
            sage: ele in CHA3
            True
            sage: cb_ele in CHA3
            False
            sage: cb_ele in CHA3.cubic_braid_group_algebra()
            True
        """
        ch_algebra = self.parent()
        cbraid_group_algebra = ch_algebra.cubic_braid_group_algebra()
        cbraid_group = ch_algebra.cubic_braid_group()

        def phi(bas_ele):
            return cbraid_group_algebra(cbraid_group(bas_ele))
        return ch_algebra._apply_module_morphism(self, phi,
                                                 codomain=cbraid_group_algebra)

    @cached_method
    def matrix(self, subdivide=False, representation_type=None, original=False):
        r"""
        Return certain types of matrix representations of ``self``.

        The absolutely irreducible representations of the cubic Hecke algebra
        are constructed using the ``GAP3`` interface and the ``CHEVIE`` package
        if ``GAP3`` and ``CHEVIE`` are installed on the system. Furthermore,
        the representations given on `Ivan Marin's homepage
        <http://www.lamfa.u-picardie.fr/marin/representationH4-en.html>`__
        are used:

        INPUT:

        - ``subdivide`` -- boolean (default: ``False``); this boolean is passed
          to the block_matrix function
        - ``representation_type`` -- instance of enum :class:`RepresentationType`;
          this can be obtained by the attribute :attr:`CubicHeckeAlgebra.repr_type`
          of ``self``; the following values are possible:

          - ``RegularLeft`` -- (regular left repr. from the above URL)
          - ``RegularRight`` -- (regular right repr. from the above URL)
          - ``SplitIrredChevie`` -- (split irred. repr. via CHEVIE)
          - ``SplitIrredMarin`` -- (split irred. repr. from the above URL)
          - default:  ``SplitIrredChevie`` taken if GAP3 and CHEVIE are installed
            on the system, otherwise the default will be ``SplitIrredMarin``

        - ``original`` -- boolean (default: ``False``); if set to ``True`` the base
          ring of the matrix will be the generic base_ring resp. generic extension
          ring (for the split versions) of the parent of ``self``

        OUTPUT:

        An instance of :class:`~sage.algebras.hecke_algebras.cubic_hecke_matrix_rep.CubicHeckeMatrixRep`,
        which is inherited from :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense`.
        In the case of the irreducible representations the matrix is given as a
        block matrix. Each single irreducible can be obtained as item indexed by
        the members of the enum :class:`AbsIrreducibeRep` available via
        :attr:`CubicHeckeAlgebra.irred_repr`.
        For details type: ``CubicHeckeAlgebra.irred_repr?``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: CHA3.inject_variables()
            Defining c0, c1
            sage: c0m = c0.matrix()
            sage: c0m[CHA3.irred_repr.W3_111]
            [                      -b - a + u     0    0]
            [(-2*a + u)*b - 2*a^2 + 2*u*a - v     b    0]
            [                               b     1    a]

        using the the ``representation_type`` option::

            sage: CHA3.<c0, c1> = algebras.CubicHecke(3)     #  optional gap3
            sage: chevie = CHA3.repr_type.SplitIrredChevie   #  optional gap3
            sage: c0m_ch = c0.matrix(representation_type=chevie) #  optional gap3
            sage: c0m_ch[CHA3.irred_repr.W3_011]             #  optional gap3
            [         b          0]
            [        -b -b - a + u]
            sage: c0m[CHA3.irred_repr.W3_011]
            [            b             0]
            [a^2 - u*a + v    -b - a + u]

        using the the ``original`` option::

            sage: c0mo = c0.matrix(original=True)
            sage: c0mo_ch = c0.matrix(representation_type=chevie, original=True) #  optional gap3
            sage: c0mo[CHA3.irred_repr.W3_011]
            [  b   0]
            [b*c   c]
            sage: c0mo_ch[CHA3.irred_repr.W3_011]            #  optional gap3
            [ b  0]
            [-b  c]

        specialized matrices::

            sage: t = (3,7,11)
            sage: CHA4 = algebras.CubicHecke(4, cubic_equation_roots=t)  # optional database_cubic_hecke
            sage: e = CHA4.an_element(); e                     # optional database_cubic_hecke
            -231*c0*c1^-1 + 131*c0*c2^-1 + 21*c2*c1 - 1440/11
            sage: em = e.matrix()                              # optional database_cubic_hecke
            sage: em.base_ring()                               # optional database_cubic_hecke
            Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
               over Integer Ring localized at (3, 7, 11)
            sage: em.dimensions()                              # optional database_cubic_hecke
            (108, 108)
            sage: em_irr24 = em[23]                            # optional database_cubic_hecke
            sage: em_irr24.dimensions()                        # optional database_cubic_hecke
            (9, 9)
            sage: em_irr24[3,2]                                # optional database_cubic_hecke
            -131*E3 - 393/7
            sage: emg = e.matrix(representation_type=chevie)   # optional gap3 database_cubic_hecke
            sage: emg_irr24 = emg[23]                          # optional gap3 database_cubic_hecke
            sage: emg_irr24[3,2]                               # optional gap3 database_cubic_hecke
            -131*E3 - 393/7
        """
        parent = self.parent()
        MS = CubicHeckeMatrixSpace(parent, representation_type=representation_type, subdivide=subdivide, original=original)
        return MS(self)

    def revert_garside(self):
        r"""
        Return the image of ``self`` under the Garside involution.

        .. SEEALSO::

            :meth:`CubicHeckeAlgebra.garside_involution`

        EXAMPLES::

            sage: roots = (E(3), ~E(3), 1)
            sage: CHA3.<c1, c2> = algebras.CubicHecke(3, cubic_equation_roots=roots)
            sage: e = CHA3.an_element(); e
            -c1*c2^-1
            sage: _.revert_garside()
            -c2*c1^-1
            sage: _.revert_garside()
            -c1*c2^-1
        """
        return self.parent().garside_involution(self)

    def revert_mirror(self):
        r"""
        Return the image of ``self`` under the mirror isomorphism.

        .. SEEALSO::

            :meth:`CubicHeckeAlgebra.mirror_isomorphism`

        EXAMPLES::

            sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
            sage: e = CHA3.an_element()
            sage: e.revert_mirror()
            -1/w*c0^-1*c1 + u/w*c0^-1 + v/w*c1^-1 + ((v*w-u)/w)
            sage: _.revert_mirror() == e
            True
        """
        return self.parent().mirror_isomorphism(self)

    def revert_orientation(self):
        r"""
        Return the image of ``self`` under the anti involution reverting the
        orientation of braids.

        .. SEEALSO::

            :meth:`CubicHeckeAlgebra.orientation_antiinvolution`

        EXAMPLES::

            sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
            sage: e = CHA3.an_element()
            sage: e.revert_orientation()
            (-w)*c2^-1*c1 + v*c1 + u*c2 + ((-v*w+u)/w)
            sage: _.revert_orientation() == e
            True
        """
        return self.parent().orientation_antiinvolution(self)

    def formal_markov_trace(self, extended=False, field_embedding=False):
        r"""
        Return a formal expression which can be specialized to Markov traces
        which factor through the cubic Hecke algebra.

        This covers Markov traces corresponding to the

        - HOMFLY-PT polynomial,
        - Kauffman polynomial,
        - Links-Gould polynomial.

        These expressions are elements of a sub-module of the module of linear
        forms on ``self`` the base ring of which is an extension of the
        generic base ring of ``self`` by an additional variable ``s``
        representing the writhe factor. All variables of this base ring
        extension are invertible.

        A Markov trace is a family of class functions `tr_n` on the family
        of braid groups `B_n` into some commutative ring `R` depending on
        a unit `s \in R` such that for all `b \in B_n` the following two
        conditions are satisfied (see [Kau1991]_, section 7):

        .. MATH::

            \begin{array}{lll}
            tr_{n+1}(b g_n)      &  = &  s tr_n(b), \\
            tr_{n+1}(b g^{-1}_n) &  = &  s^{-1} tr_n(b).
            \end{array}

        The unit `s` is often called the writhe factor and corresponds to the
        additional variable mentioned above.

        .. NOTE::

            Currently it is not known if all linear forms of this sub-module
            belong to a Markov trace, i.e. can be extended to the full tower
            of cubic Hecke algebras. Anyway, at least the four basis elements
            (``U1``, ``U2``, ``U3`` and ``K4``) can be reconstructed form
            the HOMFLY-PT and Kauffman polynomial.

        INPUT:

        - ``extended`` -- boolean (default: ``False``); if set to ``True`` the
          base ring of the Markov trace module is constructed as an extension
          of generic extension ring of ``self``; per default it is constructed
          upon the generic base ring
        - ``field_embedding`` -- boolean (default: ``False``); if set to ``True``
          the base ring of the module is the smallest field containing the
          generic extension ring of ``self``; ignored if ``extended=False``

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: CHA2 = algebras.CubicHecke(2)
            sage: K3_1 = KnotInfo.K3_1
            sage: b3_1 = CHA2(K3_1.braid())
            sage: mt3_1 = b3_1.formal_markov_trace(); mt3_1
            ((u^2*s^2-v*s^2+u*w)/s)*B[U1] + (-u*v+w)*B[U2]
            sage: mt3_1.parent()
            Free module generated by {U1, U2}
               over Multivariate Polynomial Ring in u, v, w, s
               over Integer Ring localized at (s, w, v, u)

            sage: f = b3_1.formal_markov_trace(extended=True); f
            (a^2*b*c*s^-1+a*b^2*c*s^-1+a*b*c^2*s^-1+a^2*s+a*b*s+b^2*s+a*c*s+b*c*s+c^2*s)*B[U1]
              + (-a^2*b-a*b^2-a^2*c+(-2)*a*b*c-b^2*c-a*c^2-b*c^2)*B[U2]
            sage: f.parent().base_ring()
            Multivariate Laurent Polynomial Ring in a, b, c, s
              over Splitting Algebra of x^2 + x + 1 with roots [e3, -e3 - 1]
              over Integer Ring

            sage: f = b3_1.formal_markov_trace(extended=True, field_embedding=True); f
            ((a^2*b*c+a*b^2*c+a*b*c^2+a^2*s^2+a*b*s^2+b^2*s^2+a*c*s^2+b*c*s^2+c^2*s^2)/s)*B[U1]
              + (-a^2*b-a*b^2-a^2*c-2*a*b*c-b^2*c-a*c^2-b*c^2)*B[U2]
            sage: f.parent().base_ring()
            Fraction Field of Multivariate Polynomial Ring in a, b, c, s
              over Cyclotomic Field of order 3 and degree 2

        Obtaining the well known link invariants from it::

            sage: MT = mt3_1.base_ring()
            sage: sup = mt3_1.support()
            sage: u, v, w, s = mt3_1.base_ring().gens()
            sage: LK3_1 = mt3_1*s**-3 # since the writhe of K3_1 is 3
            sage: f = MT.specialize_homfly()
            sage: g = sum(f(LK3_1.coefficient(b)) * b.regular_homfly_polynomial() for b in sup); g
            L^-2*M^2 - 2*L^-2 - L^-4
            sage: g == K3_1.link().homfly_polynomial()
            True

            sage: f = MT.specialize_kauffman()
            sage: g = sum(f(LK3_1.coefficient(b)) * b.regular_kauffman_polynomial() for b in sup); g
            a^-2*z^2 - 2*a^-2 + a^-3*z + a^-4*z^2 - a^-4 + a^-5*z
            sage: g == K3_1.kauffman_polynomial()
            True

            sage: f = MT.specialize_links_gould()
            sage: g = sum(f(LK3_1.coefficient(b)) * b.links_gould_polynomial() for b in sup); g
            -t0^2*t1 - t0*t1^2 + t0^2 + 2*t0*t1 + t1^2 - t0 - t1 + 1
            sage: g == K3_1.link().links_gould_polynomial()
            True
        """
        cha = self.parent()
        vs = self.to_vector()
        mtcf = cha._markov_trace_coeffs()
        M = cha._markov_trace_module(extended=extended, field_embedding=field_embedding)
        if M != mtcf[0].parent():
            if field_embedding:
                # intermediate step needed since internal coercion to the field
                # maps (u, v, w ) -> (a, b, c)
                MI = cha._markov_trace_module(extended=extended, field_embedding=False)
                RI = MI.base_ring()
                mtcf = [MI.from_vector(cf.to_vector()) for cf in mtcf]
                vs = vs.change_ring(RI)
            mtcf = [M.from_vector(cf.to_vector()) for cf in mtcf]

        R = M.base_ring()
        return M.linear_combination((mtcf[i], R(val)) for i, val in vs.iteritems())


class CubicHeckeAlgebra(CombinatorialFreeModule):
    r"""
    Return the Cubic-Hecke algebra with respect to the Artin braid group on
    `n` strands.

    This is a quotient of the group algebra of the Artin braid group, such that
    the images `s_i` (`1 \leq i < n`) of the braid generators satisfy a cubic
    equation (see :mod:`~sage.algebras.hecke_algebras.cubic_hecke_algebra`
    for more information, in a session type
    ``sage.algebras.hecke_algebras.cubic_hecke_algebra?``):

    .. MATH::

        s_i^3 = u s_i^2 - v s_i + w.

    The base ring of this algebra can be specified by giving optional keywords
    described below. If no keywords are given, the base ring will be a
    :class:`CubicHeckeRingOfDefinition`, which is constructed as the
    polynomial ring in `u, v, w` over the integers localized at `w`.
    This ring will be called the *ring of definition* or sometimes for short
    *generic base ring*. However note, that in this context the word *generic*
    should not remind in a generic point of the corresponding scheme.

    In addition to the base ring, another ring containing the roots (`a`, `b`
    and `c`) of the cubic equation will be needed to handle the split
    irreducible representations. This ring will be called the *extension ring*.
    Generically, the extension ring will be a
    :class:`~sage.algebras.hecke_algebras.cubic_hecke_base_ring.CubicHeckeExtensionRing`,
    which is constructed as the Laurent polynomial ring in `a, b` and `c` over
    the integers adjoined with a primitive third root of unity. A special form
    of this *generic extension ring* is constructed as a
    :class:`~sage.algebras.splitting_algebra.SplittingAlgebra` for the roots of
    the cubic equation and a primitive third root of unity over the ring of
    definition. This ring will be called the *default extension ring*.

    This class uses a static and a dynamic data library. The first one is defined
    as instance of :class:`~sage.databases.cubic_hecke_db.CubicHeckeDataBase`
    and contains the complete basis for the algebras with less than 5 strands
    and various types of representation matrices of the generators. These data
    have been calculated by `Ivan Marin <http://www.lamfa.u-picardie.fr/marin/anglais.html>`__
    and have been imported from his corresponding
    `web page <http://www.lamfa.u-picardie.fr/marin/representationH4-en.html>`__.

    Note that just the data for the cubic Hecke algebras on less than four
    strands is available in Sage by default. To deal with four strands and
    more you need to install the optional package
    `database_cubic_hecke <https://pypi.org/project/database-cubic-hecke/>`__
    by typing

    - ``sage -i database_cubic_hecke`` (first time installation) or
    - ``sage -f database_cubic_hecke`` (reinstallation) respective
    - ``sage -i -c database_cubic_hecke`` (for running all test in concern)
    - ``sage -f -c database_cubic_hecke``

    This will add a `Python wrapper <https://github.com/soehms/database_cubic_hecke#readme>`__
    around Ivan Marin's data to the Sage library. For more installation hints
    see the documentation of this wrapper.

    Furthermore, representation matrices can be obtained from the ``CHEVIE``
    package of ``GAP3`` via the ``GAP3`` interface if ``GAP3`` is installed
    inside Sage. For more information on how to obtain representation matrices
    to elements of this class, see the documentation of the element class
    :class:`~sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeElement`
    or its method
    :meth:`~sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeElement.matrix`:

        ``algebras.CubicHecke.Element?`` or ``algebras.CubicHecke.Element.matrix?``

    The second library is created as instance of
    :class:`~sage.databases.cubic_hecke_db.CubicHeckeFileCache` and used while
    working with the class to achieve a better performance. This file cache
    contains images of braids and representation matrices of basis elements
    from former calculations. A refresh of the file cache can be done using
    the :meth:`reset_filecache`.

    INPUT:

    - ``names`` -- string containing the names of the generators as images of
      the braid group generators
    - ``cubic_equation_parameters`` -- tuple ``(u, v, w)`` of three elements
      in an integral domain used as coefficients in the cubic equation. If this
      argument is given the base ring will be set to the common parent of
      ``u, v, w``. In addition a conversion map from the generic base ring is
      supplied. This keyword can also be used to change the variable names of
      the generic base ring (see example 3 below)
    - ``cubic_equation_roots`` -- tuple ``(a, b, c)`` of three elements in an
      integral domain which stand for the roots of the cubic equation. If this
      argument is given the extension ring will be set to the common parent of
      ``a, b, c``. In addition a conversion map from the generic extension ring
      and the generic base ring is supplied. This keyword can also be used to
      change the variable names of the generic extension ring (see example 3
      below)

    EXAMPLES:

    Cubic Hecke algebra over the ring of definition::

        sage: CHA3 = algebras.CubicHecke('s1, s2'); CHA3
        Cubic Hecke algebra on 3 strands over Multivariate Polynomial Ring
          in u, v, w
          over Integer Ring localized at (w,)
          with cubic equation: h^3 - u*h^2 + v*h - w = 0
        sage: CHA3.gens()
        (s1, s2)
        sage: GER = CHA3.extension_ring(generic=True); GER
        Multivariate Laurent Polynomial Ring in a, b, c
          over Splitting Algebra of x^2 + x + 1
            with roots [e3, -e3 - 1] over Integer Ring
        sage: ER = CHA3.extension_ring(); ER
        Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
          over Splitting Algebra of h^3 - u*h^2 + v*h - w
            with roots [a, b, -b - a + u]
          over Multivariate Polynomial Ring in u, v, w
          over Integer Ring localized at (w,)

    Element construction::

        sage: ele = CHA3.an_element(); ele
        (-w)*s1*s2^-1 + v*s1 + u*s2 + ((-v*w+u)/w)
        sage: ele2 = ele**2; ele2
        w^2*(s1^-1*s2)^2 + (-u*w^2)*s1^-1*s2*s1^-1 + (-v*w)*s2*s1^-1*s2
        + (-v*w^2)*s1^-1*s2^-1 + u*w*s1*s2*s1^-1*s2 + (-u*w)*s1^-1*s2*s1
        + (-u*v*w+2*v*w-2*u)*s1*s2^-1 + u*v*w*s2*s1^-1 + u*v*s2*s1 + v^2*w*s1^-1
        + (-u^2*w)*s1*s2*s1^-1 + ((u*v^2*w-2*v^2*w-u*w^2+2*u*v)/w)*s1
        + u*v*s1*s2 + (u^2*w+v^2*w)*s2^-1 + ((u^3*w-2*u*v*w+2*u^2)/w)*s2
        + ((-u^2*v*w^2-v^3*w^2+v^2*w^2-2*u*v*w+u^2)/w^2)
        sage: B3 = CHA3.braid_group()
        sage: braid = B3((2,-1, 2, 1)); braid
        s2*s1^-1*s2*s1
        sage: ele3 = CHA3(braid); ele3
        s1*s2*s1^-1*s2 + u*s1^-1*s2*s1 + (-v)*s1*s2^-1 + v*s2^-1*s1 + (-u)*s1*s2*s1^-1
        sage: ele3t = CHA3((2,-1, 2, 1))
        sage: ele3 == ele3t
        True
        sage: CHA4 = algebras.CubicHecke(4)     # optional database_cubic_hecke
        sage: ele4 = CHA4(ele3); ele4           # optional database_cubic_hecke
        c0*c1*c0^-1*c1 + u*c0^-1*c1*c0 + (-v)*c0*c1^-1 + v*c1^-1*c0 + (-u)*c0*c1*c0^-1

    Cubic Hecke algebra over the ring of definition using different variable
    names::

        sage: algebras.CubicHecke(3, cubic_equation_parameters='u, v, w', cubic_equation_roots='p, q, r')
        Cubic Hecke algebra on 3 strands over Multivariate Polynomial Ring
          in u, v, w
          over Integer Ring localized at (w,)
            with cubic equation: h^3 - u*h^2 + v*h - w = 0
        sage: _.extension_ring()
        Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
          over Splitting Algebra of h^3 - u*h^2 + v*h - w
            with roots [p, q, -q - p + u]
          over Multivariate Polynomial Ring in u, v, w
          over Integer Ring localized at (w,)

    Cubic Hecke algebra over a special base ring with respect to a special
    cubic equation::

        sage: algebras.CubicHecke('s1, s2', cubic_equation_parameters=(QQ(1),3,1))
        Cubic Hecke algebra on 3 strands over Rational Field
          with cubic equation: h^3 - h^2 + 3*h - 1 = 0
        sage: CHA3 = _
        sage: ER = CHA3.extension_ring(); ER
        Number Field in T with defining polynomial T^12 + 4*T^11 + 51*T^10
        + 154*T^9 + 855*T^8 + 1880*T^7 + 5805*T^6 + 8798*T^5 + 15312*T^4
        + 14212*T^3 + 13224*T^2 + 5776*T + 1444
        sage: CHA3.cubic_equation_roots()[0]
        -4321/1337904*T^11 - 4181/445968*T^10 - 4064/27873*T^9 - 51725/167238*T^8
        - 2693189/1337904*T^7 - 1272907/445968*T^6 - 704251/74328*T^5
        - 591488/83619*T^4 - 642145/83619*T^3 + 252521/111492*T^2 + 45685/5868*T
        + 55187/17604

        sage: F = GF(25,'u')
        sage: algebras.CubicHecke('s1, s2', cubic_equation_parameters=(F(1), F.gen(), F(3)))
        Cubic Hecke algebra on 3 strands over Finite Field in u of size 5^2
          with cubic equation: h^3 + 4*h^2 + u*h + 2 = 0
        sage: CHA3 = _
        sage: ER = CHA3.extension_ring(); ER
        Finite Field in S of size 5^4
        sage: CHA3.cubic_equation_roots()
        [2*S^3 + 2*S^2 + 2*S + 1, 2*S^3 + 3*S^2 + 3*S + 2, S^3 + 3]


    Cubic Hecke algebra over a special extension ring with respect to special
    roots of the cubic equation::

        sage: UCF = UniversalCyclotomicField()
        sage: e3=UCF.gen(3); e5=UCF.gen(5)
        sage: algebras.CubicHecke('s1, s2', cubic_equation_roots=(1, e5, e3))
        Cubic Hecke algebra on 3 strands over Universal Cyclotomic Field
          with cubic equation:
          h^3 + (-E(15) - E(15)^4 - E(15)^7 + E(15)^8)*h^2 + (-E(15)^2 - E(15)^8
          - E(15)^11 - E(15)^13 - E(15)^14)*h - E(15)^8 = 0

    TESTS::

        sage: CHA3 = algebras.CubicHecke(3)
        sage: TestSuite(CHA3).run()

    Note, that the ``TestSuite`` run on the cubic Hecke algebra on four strands
    would take up to half an hour if the file cache is empty. A repetition takes
    less than half a minute.
    """
    Element = CubicHeckeElement
    repr_type = RepresentationType
    irred_repr = AbsIrreducibeRep

    ############################################################################
    # private methods
    ############################################################################
    @staticmethod
    def __classcall_private__(cls, n=None, names='c', cubic_equation_parameters=None, cubic_equation_roots=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2, 'd', cubic_equation_roots=(3,5,7)); CHA2
            Cubic Hecke algebra on 2 strands
              over Integer Ring localized at (3, 5, 7)
                with cubic equation:
                h^3 - 15*h^2 + 71*h - 105 = 0
            sage: CHA2.inject_variables()
            Defining d
            sage: CHA3 = algebras.CubicHecke(3,  cubic_equation_parameters=(3,5,7)); CHA3
            Cubic Hecke algebra on 3 strands
              over Integer Ring localized at (7,)
                with cubic equation:
                h^3 - 3*h^2 + 5*h - 7 = 0
            sage: CHA3.cubic_equation_roots()
            [a, b, -b - a + 3]
        """
        # Support Freegroup('a,b') syntax
        if n is not None:
            try:
                n = ZZ(n) - 1
            except TypeError:
                names = n
                n = None

        # derive n from counting names
        if n is None:
            if type(names) is str:
                n = len(names.split(','))
            else:
                names = list(names)
                n = len(names)

        from sage.structure.category_object import normalize_names
        names = tuple(normalize_names(n, names))
        return super().__classcall__(cls, names,
                                     cubic_equation_parameters=cubic_equation_parameters,
                                     cubic_equation_roots=cubic_equation_roots)

    def __init__(self, names, cubic_equation_parameters=None, cubic_equation_roots=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: CHA2 = algebras.CubicHecke(2, 'd', cubic_equation_roots=(3,5,7))
            sage: TestSuite(CHA2).run()
            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_parameters=(3,5,7))
            sage: TestSuite(CHA2).run()
        """
        # ----------------------------------------------------------------------
        # Define underlying group
        # ----------------------------------------------------------------------
        self._cubic_braid_group = CubicBraidGroup(names)
        self._braid_group = self._cubic_braid_group.braid_group()
        n = len(self._cubic_braid_group.gens())
        self._nstrands = n + 1
        self._dim_irr_rep = sum([irr.dimension() for irr in AbsIrreducibeRep if irr.number_gens() == n])

        # ----------------------------------------------------------------------
        # preparing use of data base anf file cache
        # ----------------------------------------------------------------------
        from sage.databases.cubic_hecke_db import CubicHeckeDataBase, CubicHeckeFileCache
        self._database = CubicHeckeDataBase()
        self._filecache = CubicHeckeFileCache(self._nstrands)

        # ----------------------------------------------------------------------
        # interpretation of keywords and type verifications
        # ----------------------------------------------------------------------
        # cubic_equation_parameters
        # ----------------------------------------------------------------------
        ring_of_definition_names = ('u', 'v', 'w')
        if cubic_equation_parameters is not None:
            if isinstance(cubic_equation_parameters, str):
                # --------------------------------------------------------------
                # Input specifies names for the generic base ring
                # --------------------------------------------------------------
                ring_of_definition_names = tuple(cubic_equation_parameters.split(','))
                if len(ring_of_definition_names) != 3:
                    raise ValueError('cubic_equation_parameters must consist of exactly 3 elements')
                cubic_equation_parameters = None
            else:
                # --------------------------------------------------------------
                # Input specifies a specialized base ring
                # --------------------------------------------------------------
                if isinstance(cubic_equation_parameters, list):
                    cubic_equation_parameters = tuple(cubic_equation_parameters)
                if not isinstance(cubic_equation_parameters, tuple):
                    raise TypeError('cubic_equation_parameters must be a tuple or list')
                if len(cubic_equation_parameters) != 3:
                    raise ValueError('cubic_equation_parameters must consist of exactly 3 elements')

        # ----------------------------------------------------------------------
        # cubic_equation_roots
        # ----------------------------------------------------------------------
        generic_extension_ring_names = ('a', 'b', 'c')
        if cubic_equation_roots is not None:
            if isinstance(cubic_equation_roots, str):
                # --------------------------------------------------------------
                # Input specifies names for the generic extension ring
                # --------------------------------------------------------------
                generic_extension_ring_names = tuple(cubic_equation_roots.split(','))
                if len(generic_extension_ring_names) != 3:
                    raise ValueError('cubic_equation_roots must consist of exactly 3 elements')
                cubic_equation_roots = None
            else:
                # --------------------------------------------------------------
                # Input specifies a specialized base ring
                # --------------------------------------------------------------
                if isinstance(cubic_equation_roots, list):
                    cubic_equation_roots = tuple(cubic_equation_roots)
                if not isinstance(cubic_equation_roots, tuple):
                    raise TypeError('cubic_equation_roots must be a tuple or list')
                if len(cubic_equation_roots) != 3:
                    raise ValueError('cubic_equation_roots must consist of exactly 3 elements')

        if len(set(ring_of_definition_names + generic_extension_ring_names)) < 6:
            raise ValueError('there is an overlap of names between cubic equation '
                             'parameters (%s) and cubic equation roots (%s)'
                             % (ring_of_definition_names, generic_extension_ring_names))

        # ----------------------------------------------------------------------
        # setting the generic rings
        # ----------------------------------------------------------------------
        ring_of_definition = CubicHeckeRingOfDefinition(names=ring_of_definition_names)
        u, v, w = ring_of_definition.gens()

        generic_extension_ring = ring_of_definition.extension_ring(names=generic_extension_ring_names)
        a, b, c = generic_extension_ring.gens()

        # ----------------------------------------------------------------------
        # registering generic items as variables
        # ----------------------------------------------------------------------
        self._ring_of_definition = ring_of_definition
        self._generic_extension_ring = generic_extension_ring
        self._generic_cubic_equation_parameters = [u, v, w]
        self._generic_cubic_equation_roots = [a, b, c]

        # ----------------------------------------------------------------------
        # interpreting user given cubic equation parameters to define the
        # corresponding specialized base ring
        # ----------------------------------------------------------------------
        if cubic_equation_parameters is None and cubic_equation_roots is not None:
            pa, pb, pc = cubic_equation_roots
            cubic_equation_parameters = [pa+pb+pc, pa*pb+pb*pc+pa*pc, pa*pb*pc]
            verbose('cubic_equation_parameters %s set according to '
                    'cubic_equation_roots %s' % (cubic_equation_parameters,
                                                 cubic_equation_roots), level=2)

        if cubic_equation_parameters is not None:
            base_ring = ring_of_definition.create_specialization(cubic_equation_parameters)
            cubic_equation_parameters = [base_ring(para) for para in cubic_equation_parameters]
            verbose('base_ring %s set according to cubic_equation_parameters %s'
                    % (base_ring, cubic_equation_parameters), level=2)
        else:
            base_ring = self._ring_of_definition
            cubic_equation_parameters = self._generic_cubic_equation_parameters

        verbose('base_ring %s and cubic_equation_parameters %s defined'
                % (base_ring, cubic_equation_parameters), level=2)

        # ----------------------------------------------------------------------
        #  defining the cubic equation
        # ----------------------------------------------------------------------
        pu, pv, pw = cubic_equation_parameters
        pol_bas_ring = base_ring['h']
        cubic_equation = pol_bas_ring([-pw, pv, -pu, 1])

        verbose('cubic_equation %s defined' % cubic_equation, level=2)

        # ----------------------------------------------------------------------
        # defining cubic_equation_roots if not given using the cubic_equation
        # ----------------------------------------------------------------------
        if base_ring != ring_of_definition:
            if cubic_equation_roots is None:
                # --------------------------------------------------------------
                # No roots given
                # --------------------------------------------------------------
                ext_ring_names = list(generic_extension_ring_names)
                cubic_equation_roots = solve_with_extension(cubic_equation,
                                                            ext_ring_names,
                                                            var='S', flatten=True)

        # ----------------------------------------------------------------------
        # interpreting user given cubic equation roots to define the
        # corresponding specialized extension ring
        # ----------------------------------------------------------------------
        if cubic_equation_roots is not None:
            extension_ring = generic_extension_ring.create_specialization(cubic_equation_roots)
            cubic_equation_roots = [extension_ring(root) for root in cubic_equation_roots]
            verbose('extension_ring %s set according to cubic_equation_roots %s'
                    % (base_ring, cubic_equation_roots), level=2)

        else:
            extension_ring = generic_extension_ring.as_splitting_algebra()
            cubic_equation_roots = [extension_ring(a), extension_ring(b), extension_ring(c)]

        verbose('cubic roots %s and extension ring %s defined'
                % (cubic_equation_roots, extension_ring), level=2)
        pa, pb, pc = cubic_equation_roots

        # ----------------------------------------------------------------------
        # check keywords plausibility
        # ----------------------------------------------------------------------
        if base_ring == extension_ring:
            val_a = cubic_equation.substitute(h=pa)
            val_b = cubic_equation.substitute(h=pb)
            val_c = cubic_equation.substitute(h=pc)
            if val_a != 0 or val_b != 0 or val_c != 0:
                raise ValueError('cubic equation does not vanish on cubic equation roots')

        # ----------------------------------------------------------------------
        # defining the base ring embedding into the extension ring
        # ----------------------------------------------------------------------
        im_base_gens = [pa+pb+pc, pa*pb+pa*pc+pb*pc, pa*pb*pc]
        base_ring_embedding = extension_ring.coerce_map_from(base_ring)

        def check_base_ring_embedding(base_ring_embedding):
            if base_ring_embedding is None:
                return False
            try:
                ipu = base_ring_embedding(pu)
                ipv = base_ring_embedding(pv)
                ipw = base_ring_embedding(pw)
                if [ipu, ipv, ipw] != im_base_gens:
                    return False
            except (TypeError, ValueError):
                return False
            return True

        if check_base_ring_embedding(base_ring_embedding):
            verbose('base_ring_embedding defined via coercion', level=2)
        else:
            base_ring_embedding = extension_ring.convert_map_from(base_ring)
            if check_base_ring_embedding(base_ring_embedding):
                verbose('base_ring_embedding defined via conversion', level=2)
            else:
                try:
                    if base_ring.gens() == cubic_equation_parameters:
                        base_ring_embedding = base_ring.hom(im_base_gens, codomain=extension_ring)
                except (TypeError, ValueError):
                    base_ring_embedding = None

        if base_ring_embedding is None:
            warn('Warning: no base_ring_embedding found')

        # ----------------------------------------------------------------------
        # registering variables
        # ----------------------------------------------------------------------
        self._extension_ring = extension_ring
        self._base_ring_embedding = base_ring_embedding
        self._ring_of_definition_map = base_ring.convert_map_from(ring_of_definition)
        self._generic_extension_ring_map = extension_ring.convert_map_from(generic_extension_ring)
        self._cubic_equation_parameters = cubic_equation_parameters
        self._cubic_equation_roots = cubic_equation_roots

        # ----------------------------------------------------------------------
        # defining the associated group algebras
        # ----------------------------------------------------------------------
        from sage.algebras.group_algebra import GroupAlgebra
        self._cubic_braid_group_algebra = GroupAlgebra(self._cubic_braid_group, R=base_ring)
        self._braid_group_algebra = GroupAlgebra(self._braid_group, R=base_ring)

        # ----------------------------------------------------------------------
        # Setup of Basis
        # ----------------------------------------------------------------------
        # Fetch Ivan Marin's basis for the algebras on at most 4 strands.
        # An explicit list of basis element which represents a flat deformation
        # of the cubic braid group is only available in the cases where the
        # number of strands is less than 5. In the case of exactly 5 strands
        # it is known that such a basis exist by work of Ivan Marin in
        # [Marin2012] but it can not be calculated, right now. In the (infinite
        # dimensional) cases of more than 5 strands it is even an open problem
        # if the cubic Hecke algebra is a flat deformation of the group algebra
        # of the corresponding cubic braid group.
        #
        # But anyway, we will take the elements of the cubic braid group as a
        # basis of the cubic Hecke algebra in all cases. Beware that this might
        # not cover the whole cubic Hecke algebra if the number of strands is
        # larger than 4
        #
        # Internally the basis is implemented using two lists one of which
        # consists of fixed braid pre images of the basis elements and the other
        # (redundant) of the corresponding Tietze expressions.
        #
        # In the case of less than 5 strands these lists are directly obtained
        # from the list calculated by Ivan Marin available at
        #
        # In the other cases these lists are implemented as growing list which
        # is initialized with Marin's list and is extended on demand.
        # ----------------------------------------------------------------------
        db = self._database
        ns = min(self._nstrands, 4)
        self._basis_static = db.read(db.section.basis, nstrands=ns)

        # ----------------------------------------------------------------------
        # defining the algebra itself
        # ----------------------------------------------------------------------
        if self._cubic_braid_group.is_finite():
            from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
            category = FiniteDimensionalAlgebrasWithBasis(base_ring)
        else:
            from sage.categories.algebras_with_basis import AlgebrasWithBasis
            category = AlgebrasWithBasis(base_ring)

        CombinatorialFreeModule.__init__(self, base_ring, self._cubic_braid_group,
                                         prefix='', names=names, bracket=False,
                                         category=category)

        # ----------------------------------------------------------------------
        # init the attributes being set on demand
        # ----------------------------------------------------------------------
        self._cubic_hecke_subalgebra = None
        self._mirror_image = None
        self._is_mirror = False
        self._base_ring_mirror = None
        self._gens_reg_repres_matrix = {}

        # ----------------------------------------------------------------------
        # initializing the basis extension (in case of more than 4 strands)
        # ----------------------------------------------------------------------
        self._init_basis_extension()
        return

    ############################################################################
    # --------------------------------------------------------------------------
    # overloaded inherited methods
    # --------------------------------------------------------------------------
    ############################################################################
    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT: string describing ``self``

        TESTS::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: CHA3 # indirect doctest
            Cubic Hecke algebra on 3 strands
              over Multivariate Polynomial Ring in u, v, w
              over Integer Ring localized at (w,)
                with cubic equation: h^3 - u*h^2 + v*h - w = 0
        """
        s = 'Cubic Hecke algebra on %s strands over %s with cubic equation: %s = 0'
        return s % (self._nstrands, self.base_ring(), self.cubic_equation())

    def _element_constructor_(self, x):
        r"""
        Extensions to the element constructor of class :class:`CombinatorialFreeModule`.

        New functionalities are:

        - constructing element from a braid (group homomorphism)
        - constructing element from a braid giving in Tietze form
        - constructing element from an element of the braid group algebra
          (algebra homomorphism)
        - constructing element from an element of the cubic braid group
          algebra (module homomorphism)
        - constructing element from an element of an other cubic Hecke
          algebra over an other base ring or with less strands
        - constructing element from an element of the mirror image of
          ``self`` (see method mirror_image)

        INPUT:

        - ``x`` -- can be one of the following:

          * an instance of the element class of ``self`` (but possible
            to a different parent)
          * an instance of the element class of the braid group
          * an instance of the element class of the braid group algebra
            over the base ring of ``self``
          * an instance of the element class of the cubic braid group
          * an instance of the element class of the cubic braid group
            algebra over the base ring of ``self``
          * an instance of the element class of the mirror image of
            ``self``
          * a tuple representing a braid in Tietze form
          * any other object which works for the element constructor
            of :class:`CombinatorialFreeModule`

        EXAMPLES::

            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens()
            sage: b2 = b**2
            sage: CB2 = CubicBraidGroup(2)
            sage: cb, = CB2.gens()
            sage: cb2 = cb**2
            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2(b2)
            w*c^-1 + u*c + (-v)
            sage: CHA2(cb2)
            c^-1
            sage: CHA3  = algebras.CubicHecke(3)
            sage: B3    = CHA3.braid_group()
            sage: CB3   = CHA3.cubic_braid_group()
            sage: CB3GA = CHA3.cubic_braid_group_algebra()
            sage: braid = B3((1,2,2,-1,2,1,1,-1)); braid
            c0*c1^2*c0^-1*c1*c0
            sage: img_braid = CHA3(braid); img_braid
            u*w*(c0^-1*c1)^2 + u*v*c0*c1^-1*c0 + (-u^2*w)*c0^-1*c1*c0^-1
            + (-u*v)*c1*c0^-1*c1 + (-u*v*w+w^2)*c0^-1*c1^-1 + u^2*c0*c1*c0^-1*c1
            + (-u^2*v+u*w)*c0*c1^-1 + u^2*v*c1*c0^-1 + (u^2-v)*c1*c0
            + (-u^3)*c0*c1*c0^-1 + (u*v^2-v*w)*c1^-1
            sage: cbraid = CB3(braid); cbraid
            c0*c1^2*c0^-1*c1*c0
            sage: img_cbraid = CHA3(cbraid); img_cbraid
            c0^-1*c1^-1
            sage: img_cbraid_back = img_cbraid.cubic_braid_group_algebra_pre_image()
            sage: img_cbraid_back in CB3GA
            True
            sage: img_cbraid_back == CB3GA(cbraid)
            True
        """
        braid_grp = self.braid_group()
        braid_grp_alg = self.braid_group_algebra()
        braid_img = self._braid_image

        cbraid_grp = self.cubic_braid_group()
        cbraid_grp_alg = self.cubic_braid_group_algebra()
        cbraid_img = self._cubic_braid_image

        base_ring = self.base_ring()
        ngens = self.ngens()
        params = self.cubic_equation_parameters()

        # ----------------------------------------------------------------------
        # if x is a tuple we may interpret it as a braid in Tietze form
        # ----------------------------------------------------------------------
        xb = x
        if type(x) in (tuple, list):
            x = tuple(x)
            result = self._tietze_to_finite_sub_basis_monomial(x)
            if result is not None:
                # x represents a monomial
                verbose('end from tuple %s: %s' % (x, result), level=2)
                return result

            try:
                xb = braid_grp(x)
            except (TypeError, ValueError, NotImplementedError):
                pass

        # ----------------------------------------------------------------------
        # embedding of an element of an other cubic Hecke algebra with lower
        # number of strands but same base ring
        # ----------------------------------------------------------------------
        if isinstance(xb, CubicHeckeElement):
            other_cha = xb.parent()
            other_base_ring = other_cha.base_ring()
            other_ngens = other_cha.ngens()
            other_params = other_cha.cubic_equation_parameters()
            if other_base_ring != base_ring:
                if other_ngens == ngens:
                    xbv = xb.to_vector()
                    img_xbv = vector([self.base_ring()(cf) for cf in xbv])
                    return self.from_vector(img_xbv)
                elif other_ngens < ngens:
                    sub_alg = self.cubic_hecke_subalgebra(other_ngens+1)
                    return self(sub_alg(xb))

            elif other_ngens < ngens and other_params == params:
                cbraid_preimg = xb.cubic_braid_group_algebra_pre_image()
                other_cbga = other_cha.cubic_braid_group_algebra()

                def fc(ele):
                    return cbraid_img(cbraid_grp(ele))
                result = other_cbga._apply_module_morphism(cbraid_preimg, fc, codomain=self)
                verbose('end from smaller cubic Hecke algebra %s: %s' % (xb, result), level=2)
                return result

            elif other_cha == self._mirror_image:
                result = other_cha.mirror_isomorphism(xb)
                verbose('end from mirror image %s: %s' % (xb, result), level=2)
                return result

        # ----------------------------------------------------------------------
        # if xb is an element of the braid group or its group algebra over the
        # same base ring the algebra morphism self._braid_image is applied
        # ----------------------------------------------------------------------
        if isinstance(xb, braid_grp_alg.element_class) and xb in braid_grp_alg:

            def fb(ele):
                return braid_img(ele)
            result = braid_grp_alg._apply_module_morphism(xb, fb, codomain=self)
            verbose('end from braid_group algebra %s: %s' % (xb, result), level=2)
            return result

        from sage.groups.braid import Braid
        if isinstance(xb, Braid) and xb.strands() == self._nstrands:
            result = braid_img(xb)
            verbose('end from braid_group %s: %s' % (xb, result), level=2)
            return result

        # ----------------------------------------------------------------------
        # if xb is an element of the cubic_braid group or its group algebra over
        # the same base ring xb the module morphism self._braid_image is applied
        # ----------------------------------------------------------------------
        if isinstance(xb, cbraid_grp_alg.element_class) and xb in cbraid_grp_alg:
            result = cbraid_grp_alg._apply_module_morphism(xb, cbraid_img, codomain=self)
            verbose('end from cubic braid_group algebra %s: %s' % (xb, result), level=2)
            return result

        from sage.groups.cubic_braid import CubicBraidElement
        if isinstance(xb, CubicBraidElement) and xb.parent().strands() == self._nstrands:
            result = cbraid_img(xb)
            verbose('end from cubic braid_group %s: %s' % (xb, result), level=2)
            return result

        # ----------------------------------------------------------------------
        # doing the default construction by inheritance
        # ----------------------------------------------------------------------
        result = CombinatorialFreeModule._element_constructor_(self, x)
        verbose('end (default) %s: %s' % (xb, result), level=2)
        return result

    def get_order(self):
        r"""
        Return an ordering of the basis of ``self``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: len(CHA3.get_order())
            24
        """
        # The reason we have overridden this is that we have to care about
        # the dynamical growth of thefinite sub basis used for the
        # calculation in case of more than 4 strands.

        if self._nstrands < 5:
            return self._order

        # detect change of _order of sub algebra
        sub_alg = self.cubic_hecke_subalgebra()
        former_len_sub = len(sub_alg._order)
        sub_order = sub_alg.get_order()
        if former_len_sub == len(sub_order):
            if len(self._order) == former_len_sub + len(self._basis_extension):
                return self._order
        # _order has changed! re-calculation necessary:
        cbg = self.cubic_braid_group()
        sub_order = [cbg(cb) for cb in sub_order]
        self._order = sub_order + [cbg(tup) for tup in self._basis_extension]
        return self._order

    def _order_key(self, x):
        """
        Return a key for `x` compatible with the term order.

        INPUT:

        - ``x`` -- indices of the basis of ``self``

        EXAMPLES::

            sage: A = CombinatorialFreeModule(QQ, ['x','y','a','b'])
            sage: A.set_order(['x', 'y', 'a', 'b'])
            sage: A._order_key('x')
            0
            sage: A._order_key('y')
            1
            sage: A._order_key('a')
            2
        """
        try:
            return self._rank_basis(x)
        except AttributeError:
            from sage.combinat.ranker import rank_from_list
            self._rank_basis = rank_from_list(self._order)
            return self._rank_basis(x)

    def _dense_free_module(self, base_ring=None):
        r"""
        Return a dense free module with the same dimension as ``self``.

        This overwrites the corresponding method of :class:`CombinatorialFreeModule`.
        The only difference is, that the dimension is not the dimension of
        ``self`` but the dimension of the sub-module generated by the dynamically
        growing basis given by the :meth:`get_order`. In particular there is
        no difference if the number of strands is less than 5.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._dense_free_module()
            Ambient free module of rank 3
              over the integral domain Multivariate Polynomial Ring in u, v, w
              over Integer Ring localized at (w,)
        """
        if base_ring is None:
            base_ring = self.base_ring()
        from sage.modules.free_module import FreeModule
        return FreeModule(base_ring, len(self.get_order()))

    def ngens(self):
        r"""
        The number of generators of the algebra.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.ngens()
            1
        """
        return self._nstrands - 1

    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.algebra_generators()
            Finite family {c: c}
        """
        from sage.sets.family import Family
        return Family(self._cubic_braid_group.gens(), self.monomial)

    def gens(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.gens()
            (c,)
        """
        return tuple(self.algebra_generators())

    def gen(self, i):
        r"""
        The ``i``-th generator of the algebra.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.gen(0)
            c
        """
        return self.gens()[i]

    def one_basis(self):
        r"""
        Return the index of the basis element for the identity element
        in the cubic braid group.

         EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.one_basis()
            1
        """
        return self.cubic_braid_group().one()

    def _an_element_(self):
        r"""
        Overwrite the original method from :mod:`~sage.combinat.free_module`
        to obtain a more interesting element for ``TestSuite``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.an_element()              # indirect doctest
            v*c + ((-v*w+u)/w)
        """
        n = self.ngens() + 1
        base_ring = self.base_ring()
        u, v, w = (base_ring(para) for para in self._cubic_equation_parameters)
        const = (u*~w - v) * self.one()

        gens = self.gens()
        first_gens = [gen for gen in gens if gens.index(gen) < 3]
        if n == 2:
            c1, = first_gens
            return const + v*c1
        elif n == 3:
            c1, c2 = first_gens
            return const + v*c1 - w*c1*~c2 + u*c2
        else:
            c1, c2, c3 = first_gens
            return const + v*c1*~c3 - w*c1*~c2 + u*c3*c2

    @cached_method
    def chevie(self):
        r"""
        Return the ``GAP3``-``CHEVIE`` realization of the corresponding
        cyclotomic Hecke algebra in the finite-dimensional case.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)  # optional gap3
            sage: CHA3.chevie()                  # optional gap3
            Hecke(G4,[[a,b,c]])
        """
        from sage.combinat.root_system.reflection_group_real import is_chevie_available
        if not is_chevie_available():
            raise NotImplementedError('this functionality needs GAP3 with package CHEVIE')

        n = self._nstrands
        if n == 3:
            st_number = 4
        elif n == 4:
            st_number = 25
        elif n == 5:
            st_number = 32
        else:
            raise NotImplementedError('CHEVIE version doesn\'t exist for this cubic Hecke algebra')

        gap3_function_str = """function(st_number, na, nb,nc)
            local a, b, c,     # embedded Indeterminates
                  ReflGroup,   # Reflection group
                  HeckeAlg;    # Hecke algebra

            a := Mvp(na); b := Mvp(nb); c := Mvp(nc);
            ReflGroup := ComplexReflectionGroup(st_number);
            HeckeAlg  := Hecke(ReflGroup, [[a, b, c]] );
            return HeckeAlg;
        end;"""

        from sage.interfaces.gap3 import gap3
        gap3_function = gap3(gap3_function_str)
        na, nb, nc = ('\"%s\"' % indet for indet in self.extension_ring(generic=True).variable_names())
        return gap3_function(st_number, na, nb, nc)

    @cached_method
    def product_on_basis(self, g1, g2):
        r"""
        Return product on basis elements indexed by ``g1`` and ``g2``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: g = CHA3.basis().keys().an_element(); g
            c0*c1
            sage: CHA3.product_on_basis(g, ~g)
            1
            sage: CHA3.product_on_basis(g, g)
            w*c0^-1*c1*c0 + (-v)*c1*c0 + u*c0*c1*c0
        """
        # ----------------------------------------------------------------------
        # short way for multiplications with one
        # ----------------------------------------------------------------------
        if g1 == g1.parent().one():
            return self.monomial(g2)

        if g2 == g2.parent().one():
            return self.monomial(g1)

        # ----------------------------------------------------------------------
        # convert to monomials
        # ----------------------------------------------------------------------
        g1 = self.monomial(g1)
        g2 = self.monomial(g2)

        result = None

        g1_Tietze = g1.Tietze()
        g2_Tietze = g2.Tietze()

        verbose('Tietze established (%s, %s)' % (g1_Tietze, g2_Tietze), level=2)

        # ----------------------------------------------------------------------
        # The product is calculated from the corresponding product of the braids
        # ----------------------------------------------------------------------
        braid_group = self.braid_group()
        braid_product = braid_group(g1_Tietze+g2_Tietze)
        result = self._braid_image(braid_product)
        return result

    ############################################################################
    # --------------------------------------------------------------------------
    # local methods
    # --------------------------------------------------------------------------
    ############################################################################
    def _basis_tietze(self):
        r"""
        Return the complete finite sub basis as list of Tietze tuples.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._basis_tietze()
            [[], [1], [-1]]
        """
        if self._nstrands > 4:
            self_sub = self.cubic_hecke_subalgebra()
            result_list = self_sub._basis_tietze() + self._basis_extension
        else:
            result_list = self._basis_static
        return result_list

    def _tietze_to_finite_sub_basis_monomial(self, tietze_tup):
        r"""
        Return the monomial corresponding to a Tietze tuple
        if it is in the finite sub basis, otherwise return ``None``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._tietze_to_finite_sub_basis_monomial([-1])
            c^-1
            sage: CHA2._tietze_to_finite_sub_basis_monomial([-2]) is None
            True
        """
        tietze_list = list(tietze_tup)
        in_basis = False
        if tietze_list in self._basis_tietze():
            in_basis = True
        elif self._nstrands > 4:
            fsb_dict = self._finite_sub_basis_tuples
            if tietze_list in list(fsb_dict.values()):
                in_basis = True
            elif self.cubic_hecke_subalgebra()._tietze_to_finite_sub_basis_monomial(tietze_tup) is not None:
                in_basis = True

        if in_basis:
            # tietze_tup represents a monomial
            B = self.basis()
            cb = self.cubic_braid_group()(tietze_tup)
            return B[cb]

        return None

    @cached_method
    def _create_matrix_list_for_one(self, representation_type):
        r"""
        Return the matrix list for the given representation type
        for ``self.one()``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._create_matrix_list_for_one(CHA2.repr_type.SplitIrredMarin)
            [[1], [1], [1]]
            sage: CHA2._create_matrix_list_for_one(CHA2.repr_type.RegularLeft)
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]
        """
        n = self._nstrands
        if representation_type.is_split():
            gen_base_ring = self.extension_ring(generic=True)
            rep_ind = [rep.internal_index() for rep in AbsIrreducibeRep if rep.number_gens() == n - 1]
            rep_dim = [rep.dimension() for rep in AbsIrreducibeRep if rep.number_gens() == n - 1]
            dim_sort = [rep_dim[rep_ind.index(i)] for i in range(len(rep_ind))]
            matrix_list = [MatrixSpace(gen_base_ring, dim_sort[i]).one() for i in range(len(rep_ind))]
        else:
            gen_base_ring = self.base_ring(generic=True)
            matrix_list = [MatrixSpace(gen_base_ring, self.dimension()).one()]
        return matrix_list

    @cached_method
    def _fetch_matrix_list_from_chevie(self, number):
        r"""
        This method reads irreducible representation of the cubic Hecke algebra
        via the *GAP3* interface from *CHEVIE*.

        INPUT:

        - ``number`` -- integer; number of the representation according to
          *CHEVIE*

        OUTPUT:

        A list of representing matrices over the generic extension ring, one
        matrix for each generators.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)              # optional gap3
            sage: CHA3._fetch_matrix_list_from_chevie(5)     # optional gap3
            [
            [ a  0]  [c c]
            [-a  c], [0 a]
            ]
        """
        GER = self.extension_ring(generic=True)
        gap3_result = self.chevie().Representations(number)
        from sage.matrix.constructor import matrix
        matrix_list_gens = [matrix(GER, mat_gap) for mat_gap in gap3_result]
        for m in matrix_list_gens:
            m.set_immutable()
        return matrix_list_gens

    # --------------------------------------------------------------------------
    # Methods for test_suite
    # --------------------------------------------------------------------------
    # _test_ring_constructions
    # --------------------------------------------------------------------------
    def _test_ring_constructions(self, **options):
        r"""
        Method called by :class:`TestSuite`.

        The following is checked:

        - construction of base_ring and extension ring
        - construction of maps between generic base and extension ring and
          the user defined rings
        - application of the ring homomorphisms

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3,7,11))
            sage: CHA2._test_ring_constructions()
        """
        # ------------------------------------------------------------------------
        # testing ring constructions
        # ------------------------------------------------------------------------
        br = self.base_ring()
        er = self.extension_ring()
        A, B, C = self.cubic_equation_roots()
        a, b, c = self.cubic_equation_roots(generic=True)
        U, V, W = self.cubic_equation_parameters()
        u, v, w = self.cubic_equation_parameters(generic=True)
        eleB = U*V - W**2
        eleBgen = u*v - w**2
        eleE = A*B - C**2
        eleEgen = a*b - c**2

        mbr = self._ring_of_definition_map
        mer = self._generic_extension_ring_map
        bri = self._base_ring_embedding

        eleBgenEmb = 0
        eleEgenEmb = 0
        eleBembE = 0

        try:
            eleBgenEmb = br(eleBgen)
        except (TypeError, ValueError, NotImplementedError):
            verbose('generic base ring map not registered')
            try:
                eleBgenEmb = mbr(eleBgen)
            except (TypeError, ValueError, NotImplementedError):
                raise RuntimeError('fatal: generic base ring map %s does not work' % mbr)
        try:
            eleEgenEmb = er(eleEgen)
        except (TypeError, ValueError, NotImplementedError):
            verbose('generic extension ring map not registered')
            try:
                eleEgenEmb = mer(eleEgen)
            except (TypeError, ValueError, NotImplementedError):
                raise RuntimeError('fatal: generic extension ring map %s does not work' % mer)

        try:
            eleBembE = er(eleB)
        except (TypeError, ValueError, NotImplementedError):
            verbose('base ring embedding map not registered')
            try:
                eleBembE = bri(eleB)
            except (TypeError, ValueError, NotImplementedError):
                raise RuntimeError('fatal: base ring embedding %s does not work' % bri)

        test_eleBgenEmb = self._tester(**options)
        test_eleBgenEmb.assertEqual(eleBgenEmb, eleB)
        test_eleEgenEmb = self._tester(**options)
        test_eleEgenEmb.assertEqual(eleEgenEmb, eleE)
        test_eleBembE = self._tester(**options)
        test_eleBembE.assertEqual(eleBembE, eleB)

    # --------------------------------------------------------------------------
    # _test_matrix_constructions
    # --------------------------------------------------------------------------
    def _test_matrix_constructions(self, **options):
        r"""
        Test that the matrix constructions are valid.

        The following is checked:

        - construction of matrices of the following types:

          * ``RepresentationType.SplitIrredChevie``
          * ``RepresentationType.SplitIrredMarin``
          * ``RepresentationType.RegularLeft``

        - multiplication of matrices and compare with the matrix of the
          corresponding product of elements
        - construction of maps between generic base and extension ring and
          the user defined rings

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3,7,11))
            sage: CHA2._test_matrix_constructions()
        """
        # ----------------------------------------------------------------------
        # testing matrix constructions
        # ----------------------------------------------------------------------
        if self.ngens() > 4:
            return

        gens = self.gens()
        b1 = gens[0]
        b2 = self.an_element()
        b12 = b1*b2
        verbose('b12 %s' % b12)

        def check_matrix(representation_type):
            m1 = b1.matrix(representation_type=representation_type)
            m2 = b2.matrix(representation_type=representation_type)
            m12mult = m1*m2
            m12mat = b12.matrix(representation_type=representation_type)
            test_matrix = self._tester(**options)
            test_matrix.assertEqual(m12mult, m12mat)

        from sage.combinat.root_system.reflection_group_real import is_chevie_available

        if is_chevie_available():
            check_matrix(RepresentationType.SplitIrredChevie)
            if self.ngens() < 3:
                check_matrix(RepresentationType.SplitIrredMarin)
        elif self.ngens() < 4:
            check_matrix(RepresentationType.SplitIrredMarin)

        if self.ngens() < 3:
            check_matrix(RepresentationType.RegularLeft)
        return

    # --------------------------------------------------------------------------
    # _init_basis_extension
    # --------------------------------------------------------------------------
    def _init_basis_extension(self):
        r"""
        Return the extension of the basis for more than 4 strands hold
        in file cache.

        The basis elements from the file are added to the elements of
        the Marin basis.

        EXAMPLES::

            sage: # optional - database_cubic_hecke
            sage: CHA5 = algebras.CubicHecke(5)
            sage: fc = CHA5._filecache
            sage: be = fc.section.basis_extensions
            sage: CHA5.reset_filecache(be)
            sage: fc.read(be)
            [[4], [-4]]
            sage: ele = CHA5.an_element()
            sage: CHA5.inject_variables()
            Defining c0, c1, c2, c3
            sage: ele2 = ele * c3
            sage: bex = fc.read(be)
            sage: bex.sort(); bex
            [[-4], [1, -3, 4], [1, -2, 4], [3, 2, 4], [4]]
        """
        self._basis_extension = []
        tietze_list = self._basis_tietze()
        cbg = self._cubic_braid_group
        self._finite_sub_basis_tuples = {}

        order_list = [cbg(tup) for tup in tietze_list]
        # We avoid the call to set_order for speed.
        # This avoids hashing the keys, which can become expensive.
        self._order = order_list
        verbose('finite sub basis length: %s' % (len(order_list)), level=2)

        if self._nstrands < 5:
            return

        # ----------------------------------------------------------------------
        # loading the extension of the basis from data file
        # ----------------------------------------------------------------------
        fc = self._filecache
        former_bas_ext = fc.read(fc.section.basis_extensions)

        # ----------------------------------------------------------------------
        # pre definition of additional basis elements
        # ----------------------------------------------------------------------
        cub_braid_group = self.cubic_braid_group()
        if not former_bas_ext:
            gens = cub_braid_group.gens()
            last_gen = gens[len(gens)-1]
            self._cubic_braid_image(last_gen, check=False)
            self._cubic_braid_image(~last_gen, check=False)
            self._filecache.update_basis_extensions(self._basis_extension)
            return

        # ----------------------------------------------------------------------
        # Installing the additional basis elements from filecache via the
        # embedding of the corresponding cubic braid
        # ----------------------------------------------------------------------
        for bas_Tietze in former_bas_ext:
            cub_braid = cub_braid_group(bas_Tietze)
            self._cubic_braid_image(cub_braid, check=False)

        verbose('finite sub basis (extended) length: %s' % (len(self.get_order())), level=2)
        self._filecache.update_basis_extensions(self._basis_extension)
        return

    # --------------------------------------------------------------------------
    # _braid_image_from_filecache
    # --------------------------------------------------------------------------
    def _braid_image_from_filecache(self, braid):
        r"""
        Return the image of the given braid in ``self`` from file cache (if
        contained).

        INPUT:

        - ``braid`` -- braid as instance of the braid group of ``self``

        OUTPUT:

        Image of braid as element of ``self``. ``None`` if the product has
        not been stored.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: b1, b2 = CHA3.braid_group().gens(); br = ~b2*b1*~b2
            sage: CHA3._braid_image_from_filecache(br)
            1/w*c0*c1*c0^-1*c1 + v/w*c1^-1*c0 + ((-u)/w)*c0*c1*c0^-1
            sage: F = CHA3.base_ring().fraction_field()
            sage: par = tuple([F(p) for p in CHA3.cubic_equation_parameters()])
            sage: CHA3F = algebras.CubicHecke(3, cubic_equation_parameters=par)
            sage: CHA3F._braid_image_from_filecache(br)
            1/w*c0*c1*c0^-1*c1 + v/w*c1^-1*c0 + ((-u)/w)*c0*c1*c0^-1
            sage: section = CHA3.filecache_section().braid_images
            sage: CHA3.reset_filecache(section)
            sage: CHA3._braid_image_from_filecache(br)
        """
        base_ring = self.base_ring()
        gen_base_ring = self.base_ring(generic=True)
        result_vect = self._filecache.read_braid_image(braid.Tietze(), gen_base_ring)
        if result_vect is not None:
            if gen_base_ring != base_ring:
                base_map = self._ring_of_definition_map
                result_vect = vector(base_ring, [base_map(cf) for cf in result_vect])
            return self.from_vector(result_vect)
        return None

    # --------------------------------------------------------------------------
    # _braid_image_to_filecache
    # --------------------------------------------------------------------------
    def _braid_image_to_filecache(self, braid_tietze, braid_image_vect):
        r"""
        Write the given braid image of to file cache.

        INPUT:

        - ``braid_tietze`` -- braid in Tietze form
        - ``braid_image_vect`` -- image of the given braid in ``self`` in vector
          representation

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: br, = CHA2.braid_group().gens(); br2 = br**2
            sage: section = CHA2.filecache_section().braid_images
            sage: CHA2.is_filecache_empty(section)    # note: 2-strand images are not automatically cached in file system
            True
            sage: CHA2._braid_image_to_filecache(br2.Tietze(), CHA2(br2).to_vector())
            sage: CHA2._braid_image_from_filecache(br2)
            w*c^-1 + u*c + (-v)
            sage: CHA2.reset_filecache(CHA2.filecache_section().braid_images)
            sage: CHA2._braid_image_from_filecache(br2) == None
            True
        """
        if self.base_ring() != self.base_ring(generic=True):
            # this should not be done for specialized base rings
            return

        self._filecache.write_braid_image(braid_tietze, braid_image_vect)
        return

    # --------------------------------------------------------------------------
    # _braid_image
    # --------------------------------------------------------------------------
    @cached_method
    def _braid_image(self, braid):
        r"""
        Return the image of the given braid in ``self``.

        INPUT:

        - ``braid`` -- :class:`~sage.groups.braid.Braid` whose image
          in ``self`` should be calculated

        OUTPUT: an instance of the element class of ``self``

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: br, = CHA2.braid_group().gens(); br2 = br**2
            sage: CHA2._braid_image(br2)
            w*c^-1 + u*c + (-v)
        """
        # ----------------------------------------------------------------------
        # first use the cubic equation to express the braid as a linear
        # combination of braids having no other exponent as 1 and -1 in their
        # defining word in the braid generators
        # ----------------------------------------------------------------------
        coeffs, braids = self._reduce_all_gen_powers(braid.Tietze())

        # ----------------------------------------------------------------------
        # in the second step the images of these "reduced" braids is calculated
        # ----------------------------------------------------------------------
        result = self.zero()
        for i in range(len(coeffs)):
            braid_image = self._braid_image_from_reduced_powers(braids[i])
            result += coeffs[i]*braid_image

        return result

    # --------------------------------------------------------------------------
    # _braid_image_from_reduced_powers
    # --------------------------------------------------------------------------
    @cached_method
    def _braid_image_from_reduced_powers(self, braid_tietze):
        r"""
        Return the image of a braid in ``self`` assuming that no successive
        repetitions occur in the Tietze form of the braid.

        INPUT:

        - ``braid_tietze`` -- tuple representing the Braid whose image in
          ``self`` should be computed; it is assumed that no successive
          repetitions occur among the entries (i.e. ``(1, 1)`` is not allowed
          but ``(1, -2, 1)`` is)

        OUTPUT: the image of the braid as an element of ``self``

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: CHA3._braid_image_from_reduced_powers((1, -2, 1))
            c0*c1^-1*c0
            sage: CHA3._braid_image_from_reduced_powers((1, -2, 1, 2))
            w*c0^-1*c1*c0^-1 + (-v)*c1*c0^-1 + u*c0*c1*c0^-1
        """
        n = self.ngens()
        braid_list = list(braid_tietze)
        len_braid = len(braid_list)

        # ----------------------------------------------------------------------
        # in the case of two strands we calculate the power of the generator
        # ----------------------------------------------------------------------
        if n == 1:
            if len_braid == 0:
                return self.one()
            k = braid_tietze[0]*len_braid
            result_vect = self._reduce_gen_power(k)
            result = self.from_vector(result_vect)
            return result

        # ----------------------------------------------------------------------
        # Try to use former calculations (from dynamic library) to obtain the
        # braid image
        # ----------------------------------------------------------------------
        result, word_decomposition = self._braid_image_from_former_calculations(braid_tietze)

        if word_decomposition is None:
            return result

        # ----------------------------------------------------------------------
        # proceed the calculation by use of the regular representation matrices
        # (given by Ivan Marin) or in case of more than 4 strands by extension
        # of the finite sub basis
        # ----------------------------------------------------------------------

        if n > 3 and (n in braid_tietze or -n in braid_tietze):
            # ------------------------------------------------------------------
            # matrices for the regular representation are at the moment just
            # available in the case of less than five strands. In the higher
            # cases the basis is realized to grow up from the basis on 4 strands
            # to use the recursion, only those cubic braids are stored as new
            # basis elements if they involve the generator with largest index
            # ------------------------------------------------------------------
            return self._braid_image_by_basis_extension(braid_tietze)

        word_left, word_result, word_right = word_decomposition
        result_vect = None

        if word_left is not None:
            # ------------------------------------------------------------------
            # Operating from the left on the pre-calculated result
            # ------------------------------------------------------------------
            vect = result.to_vector()
            braid_preimage = tuple(word_result)
            result_vect = self._mult_by_regular_rep(vect, tuple(word_left), RepresentationType.RegularLeft, braid_preimage)

        if word_right is not None:
            # ------------------------------------------------------------------
            # Operating from the right on the pre-calculated result
            # ------------------------------------------------------------------
            if result_vect is not None:
                vect = result_vect
                braid_preimage = tuple(word_left + word_result)
            else:
                vect = result.to_vector()
                braid_preimage = tuple(word_result)
            result_vect = self._mult_by_regular_rep(vect, tuple(word_right), RepresentationType.RegularRight, braid_preimage)

        result = self.from_vector(result_vect)
        return result

    # --------------------------------------------------------------------------
    # _braid_image_from_former_calculations
    # --------------------------------------------------------------------------
    def _braid_image_from_former_calculations(self, braid_tietze):
        r"""
        Return the image of a braid in ``self`` as far as this can be done by
        use of former calculations and is sure not to go into an endless
        recursion, that is

          - using the cubic Hecke sub-algebra on one strand less
          - using the file cache.

        If the image can not be calculated from former registered results this
        method returns None. Therefore, it is just intended to be used as as
        step in the complete calculation.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid whose image in
          ``self`` should be computed; he generator exponents in the braid
          word are assumed to be ``1`` or ``-1``

        OUTPUT:

        A pair (image, basis_factors) where result is an element of ``self``
        representing the image of the input if calculation was possible and
        ``None`` else-wise. If ``image == None`` the output basis_factors is
        given as a list of basis element whose product equals the input.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: CHA3._braid_image_from_former_calculations((1, -2, 1))
            (c0*c1^-1*c0, None)
            sage: CHA3._braid_image_from_former_calculations((1, -2, 1, 2))
            (c0*c1^-1*c0, ([], [1, -2, 1], [2]))
        """
        braid_list = list(braid_tietze)
        len_braid = len(braid_list)
        result = None
        n = self.ngens()

        # ----------------------------------------------------------------------
        # if the braid is in the basis take its image to be the corresponding
        # monomial
        # ----------------------------------------------------------------------
        result = self._tietze_to_finite_sub_basis_monomial(braid_tietze)
        if result is not None:
            return result, None

        # ----------------------------------------------------------------------
        # if the braid lies in a sub-algebra take its image from there.
        # ----------------------------------------------------------------------
        sub_alg = self.cubic_hecke_subalgebra()
        if n not in braid_list and -n not in braid_list:
            result = self(sub_alg(braid_tietze))
            verbose('end (%s): %s in smaller algebra' % (braid_list, result), level=2)
            return result, None

        # ----------------------------------------------------------------------
        # proceed the calculation by splitting self into a product of basis
        # elements and try to simplify.
        # ----------------------------------------------------------------------
        braid_group = self.braid_group()
        braid = braid_group(braid_tietze)
        result = self._braid_image_from_filecache(braid)
        if result is not None:
            verbose('end from file cache (%s)' % (list(braid_tietze)), level=2)
            return result, None

        # ----------------------------------------------------------------------
        # If we come here len_braid must be larger than 1 (otherwise we already
        # have found in in the basis). By recursion we check if the subwords
        # with one generator removed on the left (respectively on the right)
        # side contain a subword whose image has already been calculated. We
        # choose the longest such subword as our result.
        # ----------------------------------------------------------------------
        braid_list_red_left = [braid_tietze[j] for j in range(1, len_braid)]
        braid_list_red_right = [braid_tietze[j] for j in range(len_braid - 1)]

        result_left, word_decomp_left = self._braid_image_from_former_calculations(tuple(braid_list_red_left))
        result_right, word_decomp_right = self._braid_image_from_former_calculations(tuple(braid_list_red_right))
        if word_decomp_left is None:
            return result_left, ([braid_tietze[0]], braid_list_red_left, [])

        if word_decomp_right is None:
            return result_right, ([], braid_list_red_right, [braid_tietze[len_braid - 1]])

        word_decomp_left_left, word_decomp_left_result, word_decomp_left_right = word_decomp_left
        word_decomp_right_left, word_decomp_right_result, word_decomp_right_right = word_decomp_right

        if len(word_decomp_left_result) >= len(word_decomp_right_result):
            return result_left, ([braid_tietze[0]] + word_decomp_left_left, word_decomp_left_result, word_decomp_left_right)

        return result_right, (word_decomp_right_left, word_decomp_right_result, word_decomp_right_right + [braid_tietze[len_braid - 1]])

    # --------------------------------------------------------------------------
    # _braid_image_by_basis_expansion_
    # --------------------------------------------------------------------------
    @cached_method
    def _braid_image_by_basis_extension(self, braid_tietze):
        r"""
        Return the given braid as a new basis element of ``self`` expanding the
        incomplete order (which is just a part of the whole basis) in the case
        of more than 4 strands.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid whose image in
          ``self`` should be computed; he generator exponents in the braid
          word are assumed to be ``1`` or ``-1``

        OUTPUT: an instance of the element class of ``self``

        EXAMPLES::

            sage: # optional - database_cubic_hecke
            sage: CHA5 = algebras.CubicHecke(5)
            sage: be = CHA5.filecache_section().basis_extensions
            sage: CHA5.reset_filecache(be)
            sage: CHA5._basis_extension
            [[4], [-4]]
            sage: CHA5._braid_image_by_basis_extension((4,1))
            c3*c0
            sage: CHA5._basis_extension
            [[4], [-4], [4, 1]]

        case where the braid already has a corresponding basis element::

            sage: CHA5._braid_image_by_basis_extension((1,)) # optional - database_cubic_hecke
            c0
            sage: CHA5._basis_extension           # optional - database_cubic_hecke
            [[4], [-4], [4, 1]]

        case where the braid doesn't have corresponding basis element but depends
        on them::

            sage: CHA5._braid_image_by_basis_extension((1,1)) # optional - database_cubic_hecke
            Traceback (most recent call last):
            ...
            NotImplementedError: no algorithm available to calculate braid image of (1, 1)
        """
        cubic_braid = self._cubic_braid_group(braid_tietze)
        tup = self._cubic_braid_basis_tuple(cubic_braid)
        if tup is not None:
            bgrp = self.braid_group()
            if bgrp(braid_tietze) != bgrp(tup):
                raise NotImplementedError('no algorithm available to calculate braid image of %s' % str(braid_tietze))
            B = self.basis()
            verbose('braid-image %s in Basis' % str(braid_tietze), level=2)
            return self.monomial(B[cubic_braid])

        return self._cubic_braid_append_to_basis(cubic_braid)

    # --------------------------------------------------------------------------
    # _reduce_all_gen_powers
    # --------------------------------------------------------------------------
    @cached_method
    def _reduce_all_gen_powers(self, braid_tietze):
        r"""
        Return a linear combination of braids that have no higher powers in the
        braid generators having the same image in ``self`` than the given braid.

        This linear combination is returned as a pair of lists of braids and
        corresponding coefficients.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid whose powers should
          be reduced given in Tietze form

        OUTPUT:

        A pair of two lists: ``coeffs``, ``braids``. The fist one contains the
        coefficients corresponding to the braids in Tietze form from the second
        list of braids.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: CHA3._reduce_all_gen_powers((1, 1, -2, -2))
            ([u*v/w, (-v)/w, (-v^2)/w, (-u^2)/w, u/w, u*v/w, -u, 1, v],
             [(), (2,), (-2,), (1,), (1, 2), (1, -2), (-1,), (-1, 2), (-1, -2)])
        """
        braid_list = list(braid_tietze)
        len_braid = len(braid_list)

        # ----------------------------------------------------------------------
        # find a higher power position in braid_list
        # ----------------------------------------------------------------------
        power = 0
        pos = 0
        for i in range(len_braid - 1):
            if braid_list[i] != braid_list[i + 1]:
                continue
            pos = i
            for power in range(1, len_braid - pos + 1):
                if pos+power == len_braid:
                    break
                if braid_list[pos] != braid_list[pos+power]:
                    break
            break

        if power == 0:
            verbose('end (%s) no powers' % braid_list, level=2)
            return [self.base_ring().one()], [braid_tietze]

        # ----------------------------------------------------------------------
        # eliminate this power from braid_tietze.
        # ----------------------------------------------------------------------
        val = braid_list[pos]
        if val > 0:
            gen_ind = val
            exp = power
        else:
            gen_ind = -val
            exp = -power

        braid_list_start = [braid_list[i] for i in range(pos)]
        braid_list_end = [braid_list[i] for i in range(pos+power, len_braid)]

        # ----------------------------------------------------------------------
        # merging the new reduced tuple. Note that all the new tuples are
        # smaller than the given one, which will make the recursion terminate.
        # ----------------------------------------------------------------------
        tuple_one = tuple(braid_list_start + braid_list_end)
        tuple_gen = tuple(braid_list_start + [gen_ind] + braid_list_end)
        tuple_gen_inv = tuple(braid_list_start + [-gen_ind] + braid_list_end)

        # ----------------------------------------------------------------------
        # convert them to braids (to reduce cancellation of inverses and obvious
        # braid relations) Note that this will not increase the length of the
        # word. Thus the recursion still must terminate.
        # ----------------------------------------------------------------------
        braid_group = self.braid_group()
        braid_one = braid_group(tuple_one)
        braid_gen = braid_group(tuple_gen)
        braid_gen_inv = braid_group(tuple_gen_inv)

        # ----------------------------------------------------------------------
        # eliminate all powers from braid_tietze by recursion. The recursion
        # will terminate by the length reduction of the Tietze tuple (but not
        # necessarily by the number of generators whose exponent must be reduced).
        # ----------------------------------------------------------------------
        one_coeffs, one_braids = self._reduce_all_gen_powers(braid_one.Tietze())
        gen_coeffs, gen_braids = self._reduce_all_gen_powers(braid_gen.Tietze())
        gen_inv_coeffs, gen_inv_braids = self._reduce_all_gen_powers(braid_gen_inv.Tietze())

        cf_one, cf_gen, cf_gen_inv = self._reduce_gen_power(exp)

        one_coeffs = [cf*cf_one for cf in one_coeffs]
        gen_coeffs = [cf*cf_gen for cf in gen_coeffs]
        gen_inv_coeffs = [cf*cf_gen_inv for cf in gen_inv_coeffs]

        return one_coeffs + gen_coeffs + gen_inv_coeffs, one_braids + gen_braids + gen_inv_braids

    # --------------------------------------------------------------------------
    # _reduce_gen_power
    # --------------------------------------------------------------------------
    @cached_method
    def _reduce_gen_power(self, k):
        r"""
        Return the `k`-th power on an arbitrary generator,
        for example `c_0^k`.

        INPUT:

        - ``k`` -- integer giving the power

        OUTPUT:

        A list ``[coeff_one, coeff_gen, coeff_gen_inverse]`` of the three
        coefficients of the generators power in the span of the generator.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._reduce_gen_power(5)
            (-u^3*v + 2*u*v^2 + u^2*w - 2*v*w, u^4 - 3*u^2*v
            + v^2 + 2*u*w, u^3*w - 2*u*v*w + w^2)
        """
        n = self.ngens()

        # ----------------------------------------------------------------------
        # take it from smaller sub-algebras if possible
        # ----------------------------------------------------------------------
        if n > 1:
            sub_alg = self.cubic_hecke_subalgebra()
            return sub_alg._reduce_gen_power(k)

        # ----------------------------------------------------------------------
        # calculate it in the sub-algebra on 2 strands
        # ----------------------------------------------------------------------

        if k == 0:
            result_ele = self.one()
            result = result_ele.to_vector()

        elif abs(k) == 1:
            result_ele = self._tietze_to_finite_sub_basis_monomial((k,))
            result = result_ele.to_vector()

        else:
            if k < 0:
                right_vect = self._reduce_gen_power(k + 1)
                genTietze = (-1,)
            else:
                right_vect = self._reduce_gen_power(k - 1)
                genTietze = (1,)
            result = self._mult_by_regular_rep(right_vect, genTietze, RepresentationType.RegularLeft)
        return result

    # --------------------------------------------------------------------------
    # _mult_by_regular_rep
    # --------------------------------------------------------------------------
    @cached_method
    def _mult_by_regular_rep(self, vect, gen_tuple, representation_type, braid_preimage=None):
        r"""
        Return the product of an`element of ``self`` given as a coefficient
        vector with a sequence (tuple) of generators (that is a braid word)
        using regular representation matrices.

        The multiplication will be performed form left or right according
        to the given ``representation_type``.

        INPUT:

        - ``vect`` -- element of ``self`` in vector form (obtained by
          :meth:`to_vector`)
        - ``gen_tuple`` -- list of generators (that is a braid in Tietze form)
          which operates on ``vect``
        - ``representation_type`` -- instance of :class:`RepresentationType`
          (one of `RegularLeft` or `RegularRight`)
        - ``braid_preimage`` -- (optional) a word representing a braid whose
          image is vect (if it exist). This is used to record intermediate
          results to the dynamic library

        OUTPUT:

        The coefficient vector resulting after applying the multiplication of all
        matrices corresponding to the generators from gen_tuple.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: CHA3.inject_variables()
            Defining c0, c1
            sage: repr_type = CHA3.repr_type.RegularRight
            sage: CHA3._mult_by_regular_rep(c0.to_vector(), (1, -2, -2), repr_type)
            (u*v/w, (-u^2)/w, -u, (-v)/w, (-v^2)/w, u/w, u*v/w, 1, v,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        verbose('multiply %s (pre-image %s) by %s using %s'
                % (vect, braid_preimage, gen_tuple, representation_type), level=2)
        m = len(gen_tuple)
        braid_list = None
        if braid_preimage:
            braid_list = list(braid_preimage)

        result = vect
        for i in range(m):
            verbose('multiply image of %s with position %s' % (braid_list, i), level=2)

            if representation_type == RepresentationType.RegularLeft:
                gen_ind = gen_tuple[m - i - 1]
                if braid_list:
                    braid_list = [gen_ind] + braid_list
            else:
                gen_ind = gen_tuple[i]
                if braid_list:
                    braid_list = braid_list + [gen_ind]

            if (gen_ind, representation_type) in list(self._gens_reg_repres_matrix.keys()):
                mat = self._gens_reg_repres_matrix[(gen_ind, representation_type)]
            else:
                if gen_ind > 0:
                    gen = self.gen(gen_ind - 1)
                    mat = gen.matrix(representation_type=representation_type)
                else:
                    # data of inverse of generators is stored under negative strand-index
                    gen = self.gen(-gen_ind - 1)**(-1)
                    mat = gen.matrix(representation_type=representation_type)

                self._gens_reg_repres_matrix[(gen_ind, representation_type)] = mat

            result = mat * result

            # ------------------------------------------------------------------
            # save this intermediate result to the dynamic library
            # ------------------------------------------------------------------
            if braid_list:
                verbose('save image of %s to file cache' % braid_list, level=2)
                self._braid_image_to_filecache(tuple(braid_list), result)

        verbose('multiply %s by %s using %s result %s' % (vect, gen_tuple, representation_type, result), level=2)
        return result

    # --------------------------------------------------------------------------
    # _cubic_braid_append_to_basis
    # --------------------------------------------------------------------------

    def _cubic_braid_append_to_basis(self, cubic_braid):
        r"""
        Append the given cubic braid to the finite sub basis which is used for
        calculation of products and representation matrices.

        This only makes sense if the ``cubic_braid`` is not in this finite
        sub-basis, before. This can happen if the number of strands is more
        than 4.

        INPUT:

        - ``cubic_braid`` -- :class:`~sage.groups.cubic_braid.CubicBraid`
          whose image in ``self`` should be appended

        OUTPUT: the new monomial of ``self``

        EXAMPLES::

            sage: # optional - database_cubic_hecke
            sage: CHA5 = algebras.CubicHecke(5)
            sage: be = CHA5.filecache_section().basis_extensions
            sage: CHA5.reset_filecache(be)
            sage: CHA5._basis_extension
            [[4], [-4]]
            sage: CBG = CHA5.cubic_braid_group()
            sage: CHA5._cubic_braid_append_to_basis(CBG((4,1)))
            c3*c0
            sage: CHA5._basis_extension
            [[4], [-4], [4, 1]]
        """
        cbTietze = list(cubic_braid.Tietze())
        order = self.get_order()
        next_index = len(order)
        self._basis_extension.append(cbTietze)
        try:
            self._rank_basis.update({cubic_braid: next_index})  # supporting :meth:`get_order_key`
        except AttributeError:
            pass
        order.append(cubic_braid)
        monomial = self.monomial(cubic_braid)
        self._finite_sub_basis_tuples.update({cubic_braid: cbTietze})

        verbose('registering new basis element: %s (par %s ind %s)'
                % (cubic_braid, cubic_braid.parent(), next_index), level=2)
        self._filecache.update_basis_extensions(self._basis_extension)
        return monomial

    # --------------------------------------------------------------------------
    # _cubic_braid_basis_tuple
    # --------------------------------------------------------------------------
    def _cubic_braid_basis_tuple(self, cubic_braid):
        r"""
        Return the Tietze tuple that represents the given cubic_braid in the
        basis of ``self``.

        In the case ``self`` has more than 4 strands it may happen that the
        given cubic braid is not contained in the finite subbasis, so far.
        In this case it is automatically added to it.

        INPUT:

        - ``cubic_braid`` -- :class:`~sage.groups.cubic_braid.CubicBraid`

        OUTPUT: a tuple from the basis representing the cubic braid

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CBG = CHA2.cubic_braid_group()
            sage: CHA2._cubic_braid_basis_tuple(CBG((1, 1)))
            (-1,)
        """
        tietze_list = self._basis_tietze()
        cubic_braid_tietze = cubic_braid.Tietze()
        if list(cubic_braid_tietze) in tietze_list:
            verbose('cubic_braid_tietze: %s in basis' % str(cubic_braid_tietze), level=2)
            return cubic_braid_tietze
        else:
            if cubic_braid in self._finite_sub_basis_tuples.keys():
                verbose('cubic_braid: %s in finite_sub_basis' % cubic_braid, level=2)
                return self._finite_sub_basis_tuples[cubic_braid]

        for tup in tietze_list:
            cb_tup = self.cubic_braid_group()(tup)
            if cubic_braid == cb_tup:
                self._finite_sub_basis_tuples.update({cb_tup: tup})
                verbose('cubic_braid: %s added to finite_sub_basis with tuple %s'
                        % (cubic_braid, tup), level=2)
                return tuple(tup)
        return None

    # --------------------------------------------------------------------------
    # _cubic_braid_image
    # --------------------------------------------------------------------------
    def _cubic_braid_image(self, cubic_braid, check=True):
        r"""
        Return the given cubic braid as monomial of ``self``, that is the image
        under the map onto the basis.

        If the number of strands is larger than 4, the corresponding basis
        element may not be contained in the order of ``self``. In this
        case it will be appended here.

        INPUT:

        - ``cubic_braid`` -- :class:`~sage.groups.cubic_braid.CubicBraid`
          whose image in ``self`` should be returned
        - ``check`` -- boolean (default: ``True``); check if the given cubic
          braid is already registered in the finite sub basis; if set to
          ``False`` duplicate entries can occur

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CBG = CHA2.cubic_braid_group()
            sage: CHA2._cubic_braid_image(CBG((1,1)))
            c^-1
        """
        if check:
            tup = self._cubic_braid_basis_tuple(cubic_braid)
            if tup is not None:
                return self._tietze_to_finite_sub_basis_monomial(tup)
        return self._cubic_braid_append_to_basis(cubic_braid)

    # --------------------------------------------------------------------------
    # _extend_braid_automorphism
    # --------------------------------------------------------------------------
    def _extend_braid_automorphism(self, element, braid_automorphism):
        r"""
        Return the image of element under the extension of the given braid group
        automorphism to ``self``. It is assumed that the given
        ``braid_automorphism`` factors through ``self``.

        INPUT:

        - ``element`` -- an element class of ``self``
        - ``braid_automorphism`` -- braid group automorphism factoring through
          ``self``

        OUTPUT:

        Element class of ``self`` representing the image of element
        under the extension of the given braid group automorphism.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: br, = CHA2.gens()
            sage: CHA2.mirror_isomorphism(br)   # indirect doctest
            c^-1
        """

        result = self.zero()
        for braid in element.support():
            autom_braid = braid_automorphism(braid)
            img_braid = self._braid_image_from_reduced_powers(autom_braid.Tietze())
            result += element[braid] * img_braid

        return result

    # --------------------------------------------------------------------------
    # _markov_trace_module
    # --------------------------------------------------------------------------
    def _markov_trace_module(self, extended=False, field_embedding=False):
        r"""
        Return the module that contains the formal Markov trace as elements.

        INPUT:

        - ``extended`` -- boolean (default: ``False``); if set to ``True`` the
          base ring of the module is the Markov trace version of the generic
          extension ring of ``self``.

        - ``field_embedding`` -- boolean (default: ``False``); if set to ``True`
          the base ring of the module is the smallest field containing the
          generic extension ring of ``self``. The keyword is meaningless if
          ``extended=False``.

        OUTPUT: a :class:`~sage.combinat.free_module.CombinatorialFreeModule`

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._markov_trace_module()
            Free module generated by {U1, U2}
              over Multivariate Polynomial Ring in u, v, w, s
              over Integer Ring localized at (s, w, v, u)

            sage: CHA2._markov_trace_module(extended=True)
            Free module generated by {U1, U2}
              over Multivariate Laurent Polynomial Ring in a, b, c, s
              over Splitting Algebra of x^2 + x + 1 with roots [e3, -e3 - 1]
              over Integer Ring

            sage: CHA2._markov_trace_module(extended=True, field_embedding=True)
            Free module generated by {U1, U2}
              over Fraction Field of Multivariate Polynomial Ring in a, b, c, s
              over Cyclotomic Field of order 3 and degree 2
        """
        from sage.modules.free_module import FreeModule
        from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
        basis = [b for b in MarkovTraceModuleBasis if b.strands() <= self._nstrands]
        BRM = self.base_ring(generic=True).markov_trace_version()
        if extended:
            BRM = BRM.extension_ring()
            if field_embedding:
                emb = BRM.field_embedding()
                BRM = emb.codomain()
        return FreeModule(BRM, basis)

    # --------------------------------------------------------------------------
    # _markov_trace_coeffs
    # --------------------------------------------------------------------------
    @cached_method
    def _markov_trace_coeffs(self):
        r"""
        Return a list of formal Markov traces of the basis elements of ``self``.

        OUTPUT:

        A list of elements of the Markov trace module (over the Markov trace
        version of the generic base ring). Each entry of the list corresponds
        to the according basis element of ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2._markov_trace_coeffs()
            [B[U2], s*B[U1], 1/s*B[U1]]
            sage: M = _[0].parent(); M
            Free module generated by {U1, U2}
              over Multivariate Polynomial Ring in u, v, w, s
              over Integer Ring localized at (s, w, v, u)
        """
        M = self._markov_trace_module()
        Mbas = M.basis().keys()
        db = self._database
        sec = db.section.markov_tr_cfs
        cfs = db.read(sec, variables=M.base_ring().gens(), nstrands=self._nstrands)
        d = self.dimension()
        return [sum(cfs[bas_ele][i]*M(bas_ele) for bas_ele in Mbas) for i in range(d)]

    ############################################################################
    # --------------------------------------------------------------------------
    # public methods
    # --------------------------------------------------------------------------
    ############################################################################
    def filecache_section(self):
        r"""
        Return the ``enum`` to select a section in the file cache.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: list(CHA2.filecache_section())
            [<section.matrix_representations: 'matrix_representations'>,
             <section.braid_images: 'braid_images'>,
             <section.basis_extensions: 'basis_extensions'>,
             <section.markov_trace: 'markov_trace'>]
        """
        return self._filecache.section

    def is_filecache_empty(self, section=None):
        r"""
        Return ``True`` if the file cache of the given ``section`` is empty.
        If no ``section`` is given the answer is given for the complete
        file cache.

        INPUT:

        - ``section`` -- (default: all sections) an element of enum
          :class:`~sage.databases.cubic_hecke_db.CubicHeckeFileCache.section`
          that can be selected using :meth:`filecache_section`

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.is_filecache_empty()
            False
        """
        return self._filecache.is_empty(section=section)

    def reset_filecache(self, section=None, commit=True):
        r"""
        Reset the file cache of the given ``section`` resp. the complete
        file cache if no ``section`` is given.

        INPUT:

        - ``section`` -- (default: all sections) an element of enum
          :class:`~sage.databases.cubic_hecke_db.CubicHeckeFileCache.section`
          that can be selected using :meth:`filecache_section`
        - ``commit`` -- boolean (default: ``True``); if set to ``False`` the
          reset is not written to the filesystem

        EXAMPLES::

            sage: # optional - database_cubic_hecke
            sage: CHA5 = algebras.CubicHecke(5)
            sage: be = CHA5.filecache_section().basis_extensions
            sage: CHA5.is_filecache_empty(be)
            False
            sage: CHA5.reset_filecache(be)
            sage: CHA5.is_filecache_empty(be)
            True
        """
        fc = self._filecache
        if section == fc.section.basis_extensions:
            if self._nstrands < 5:
                raise ValueError('not allowed for less than 5 strand')
        fc.reset_library(section=section)

        if section == fc.section.basis_extensions:
            self._init_basis_extension()
        if commit:
            fc.write(section=section)

    def strands(self):
        r"""
        Return the number of strands of the braid group whose group algebra
        image is ``self``.

        EXAMPLES::

            sage: CHA4 = algebras.CubicHecke(2)
            sage: CHA4.strands()
            2
        """
        return self._nstrands

    # --------------------------------------------------------------------------
    # Garside involution
    # --------------------------------------------------------------------------
    def garside_involution(self, element):
        r"""
        Return the image of the given element of ``self`` under the extension of
        the Garside involution of braids to ``self``.

        This method may be invoked by the ``revert_garside`` method of the
        element class of ``self``, alternatively.

        INPUT:

        - ``element`` -- instance of the element class of ``self``

        OUTPUT:

        Instance of the element class of ``self`` representing the image of
        ``element`` under the extension of the Garside involution to ``self``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element()
            sage: ele_gar = CHA3.garside_involution(ele); ele_gar
            (-w)*c1*c0^-1 + u*c0 + v*c1 + ((-v*w+u)/w)
            sage: ele == CHA3.garside_involution(ele_gar)
            True
        """
        braid_group = self.braid_group()
        reverse_gens = list(braid_group.gens())
        reverse_gens.reverse()
        brgrp_garside_involution = braid_group.hom(reverse_gens, check=False)
        return self._extend_braid_automorphism(element, brgrp_garside_involution)

    # --------------------------------------------------------------------------
    # orientation anti involution
    # --------------------------------------------------------------------------
    def orientation_antiinvolution(self, element):
        r"""
        Return the image of the given element of ``self`` under the extension of
        the orientation anti involution of braids to ``self``. The orientation
        anti involution of a braid is given by reversing the order of generators
        in the braid word.

        This method may be invoked by the ``revert_orientation`` method of the
        element class of ``self``, alternatively.

        INPUT:

        - ``element`` -- instance of the element class of ``self``

        OUTPUT:

        Instance of the element class of ``self`` representing the image of
        ``element`` under the extension of the orientation reversing braid
        involution to ``self``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element()
            sage: ele_ori = CHA3.orientation_antiinvolution(ele); ele_ori
            (-w)*c1^-1*c0 + v*c0 + u*c1 + ((-v*w+u)/w)
            sage: ele == CHA3.orientation_antiinvolution(ele_ori)
            True
        """
        braid_group = self.braid_group()

        def brgrp_orientation_antiinvolution(braid):
            braid_list = list(braid.Tietze())
            braid_list.reverse()
            return braid_group(tuple(braid_list))
        return self._extend_braid_automorphism(element, brgrp_orientation_antiinvolution)

    # --------------------------------------------------------------------------
    # mirror isomorphism
    # --------------------------------------------------------------------------
    def mirror_isomorphism(self, element):
        r"""
        Return the image of the given element of ``self`` under the extension
        of the mirror involution of braids to ``self``. The mirror involution
        of a braid is given by inverting all generators in the braid word. It
        does not factor through ``self`` over the base ring but it factors
        through ``self`` considered as a `\ZZ`-module relative to the mirror
        automorphism of the generic base ring. Considering ``self`` as algebra
        over its base ring this involution defines an isomorphism of ``self``
        onto a different cubic Hecke algebra with a different cubic equation.
        This is defined over a different base (and extension) ring than
        ``self``. It can be obtained by the method ``mirror_image`` or as
        parent of the output of this method.

        This method may be invoked by the ``CubicHeckeElelemnt.revert_mirror``
        method of the element class of ``self``, alternatively.

        INPUT:

        - ``element`` -- instance of the element class of ``self``

        OUTPUT:

        Instance of the element class of the mirror image of ``self``
        representing the image of element under the extension of the braid
        mirror involution to ``self``.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ele = CHA3.an_element()
            sage: ele_mirr = CHA3.mirror_isomorphism(ele); ele_mirr
            -1/w*c0^-1*c1 + u/w*c0^-1 + v/w*c1^-1 + ((v*w-u)/w)
            sage: ele_mirr2 = ele.revert_mirror()  # indirect doctest
            sage: ele_mirr == ele_mirr2
            True
            sage: par_mirr = ele_mirr.parent()
            sage: par_mirr == CHA3
            False
            sage: par_mirr == CHA3.mirror_image()
            True
            sage: ele == par_mirr.mirror_isomorphism(ele_mirr)
            True
        """
        mirror_image = self.mirror_image()
        braid_group = self.braid_group()
        mirror_involution = braid_group.hom([~g for g in braid_group.gens()], check=False)
        # Todo: have mirror_involution be a method of :class:`BraidGroup_class`
        base_ring_mirror = self._base_ring_mirror
        element_vec = vector([base_ring_mirror(cf) for cf in list(element.to_vector())])
        element_mirr = mirror_image.from_vector(element_vec)
        return mirror_image._extend_braid_automorphism(element_mirr, mirror_involution)

    # --------------------------------------------------------------------------
    # cubic_equation
    # --------------------------------------------------------------------------
    def cubic_equation(self, var='h', as_coefficients=False, generic=False):
        r"""
        Return the cubic equation attached to ``self``.

        INPUT:

        - ``var`` -- string (default: ``'h'``); setting the indeterminate of the
          equation
        - ``as_coefficients`` -- boolean (default: ``False``);  if set to ``True``
          the list of coefficients is returned
        - ``generic`` -- boolean (default: ``False``);  if set to ``True`` the
          cubic equation will be given over the generic base ring

        OUTPUT:

        A polynomial over the base ring (resp. generic base ring if ``generic``
        is set to True). In case ``as_coefficients`` is set to ``True`` a list
        of them is returned.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(E(3), ~E(3), 1))
            sage: CHA2.cubic_equation()
            h^3 - 1
            sage: CHA2.cubic_equation(generic=True)
            h^3 - u*h^2 + v*h - w
            sage: CHA2.cubic_equation(as_coefficients=True, generic=True)
            [-w, v, -u, 1]
            sage: CHA2.cubic_equation(as_coefficients=True)
            [-1, 0, 0, 1]
        """
        BaseRing = self.base_ring(generic=generic)
        if generic:
            u, v, w = BaseRing.gens()
        else:
            u, v, w = self._cubic_equation_parameters
        cf = [-w, v, -u, 1]
        if as_coefficients:
            return cf
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        P = PolynomialRing(BaseRing, var)
        return P(cf)

    # --------------------------------------------------------------------------
    # cubic_equation_roots
    # --------------------------------------------------------------------------
    def cubic_equation_roots(self, generic=False):
        r"""
        Return the roots of the underlying cubic equation.

        INPUT:

        - ``generic`` -- boolean (default: ``False``);  if set to ``True`` the
          roots are returned as elements of the generic extension ring

        OUTPUT: a triple consisting of the roots

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3, 4, 5))
            sage: CHA2.cubic_equation()
            h^3 - 12*h^2 + 47*h - 60
            sage: CHA2.cubic_equation_roots()
            [3, 4, 5]
            sage: CHA2.cubic_equation_roots(generic=True)
            [a, b, c]
        """
        if generic:
            return self._generic_cubic_equation_roots
        else:
            return self._cubic_equation_roots

    # --------------------------------------------------------------------------
    # cubic_equation_roots
    # --------------------------------------------------------------------------
    def cubic_equation_parameters(self, generic=False):
        r"""
        Return the coefficients of the underlying cubic equation.

        INPUT:

        - ``generic`` -- boolean (default: ``False``);  if set to ``True`` the
          coefficients are returned as elements of the generic base ring

        OUTPUT: a triple consisting of the coefficients

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3, 4, 5))
            sage: CHA2.cubic_equation()
            h^3 - 12*h^2 + 47*h - 60
            sage: CHA2.cubic_equation_parameters()
            [12, 47, 60]
            sage: CHA2.cubic_equation_parameters(generic=True)
            [u, v, w]
        """
        if generic:
            return self._generic_cubic_equation_parameters
        else:
            return self._cubic_equation_parameters

    # --------------------------------------------------------------------------
    # base_ring
    # --------------------------------------------------------------------------
    def base_ring(self, generic=False):
        r"""
        Return the base ring of ``self``.

        INPUT:

        - ``generic`` -- boolean (default: ``False``); if ``True`` the ring
          of definition (here often called the generic base ring) is returned

        EXAMMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3, 4, 5))
            sage: CHA2.base_ring()
            Integer Ring localized at (2, 3, 5)
            sage: CHA2.base_ring(generic=True)
            Multivariate Polynomial Ring in u, v, w
              over Integer Ring localized at (w,)
        """
        if generic:
            return self._ring_of_definition
        else:
            return super().base_ring()

    # --------------------------------------------------------------------------
    # extension_ring
    # --------------------------------------------------------------------------
    def extension_ring(self, generic=False):
        r"""
        Return the extension ring of ``self``.

        This is an extension of its base ring containing the roots
        of the cubic equation.

        INPUT:

        - ``generic`` -- boolean (default: ``False``); if ``True`` the
          extension ring of definition (here often called the generic
          extension ring) is returned

        EXAMMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3, 4, 5))
            sage: CHA2.extension_ring()
            Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
            over Integer Ring localized at (2, 3, 5)
            sage: CHA2.extension_ring(generic=True)
            Multivariate Laurent Polynomial Ring in a, b, c
            over Splitting Algebra of x^2 + x + 1
              with roots [e3, -e3 - 1] over Integer Ring
        """
        if generic:
            return self._generic_extension_ring
        else:
            return self._extension_ring

    # --------------------------------------------------------------------------
    # cyclotomic_generator
    # --------------------------------------------------------------------------
    def cyclotomic_generator(self, generic=False):
        r"""
        Return the third root of unity as element of the extension ring.

        The only thing where this is needed is in the nine dimensional
        irreducible representations of the cubic Hecke algebra on four strands
        (see the examples of :meth:`CubicHeckeElement.matrix` for instance).

        INPUT:

        - ``generic`` -- boolean (default: ``False``); if ``True`` the
          cyclotomic generator is returned as an element extension ring of
          definition

        EXAMMPLES::

            sage: CHA2 = algebras.CubicHecke(2, cubic_equation_roots=(3, 4, 5))
            sage: CHA2.cyclotomic_generator()
            E3
            sage: CHA2.cyclotomic_generator(generic=True)
            e3
        """
        e3gen = self.extension_ring(generic=True).cyclotomic_generator()
        if generic:
            return e3gen
        else:
            return self._generic_extension_ring_map(e3gen)

    # --------------------------------------------------------------------------
    # braid_group
    # --------------------------------------------------------------------------
    def braid_group(self):
        r"""
        Return the braid group attached to ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.braid_group()
            Braid group on 2 strands
        """
        return self._braid_group

    # --------------------------------------------------------------------------
    # cubic_braid_group
    # --------------------------------------------------------------------------
    def cubic_braid_group(self):
        r"""
        Return the cubic braid group attached to ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.cubic_braid_group()
            Cubic Braid group on 2 strands
        """
        return self._cubic_braid_group

    # --------------------------------------------------------------------------
    # braid_group_algebra
    # --------------------------------------------------------------------------
    def braid_group_algebra(self):
        r"""
        Return the group algebra of braid group attached to ``self`` over the
        base ring of ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.braid_group_algebra()
            Algebra of Braid group on 2 strands
             over Multivariate Polynomial Ring in u, v, w
             over Integer Ring localized at (w,)
        """
        return self._braid_group_algebra

    # --------------------------------------------------------------------------
    # cubic_braid_group_algebra
    # --------------------------------------------------------------------------
    def cubic_braid_group_algebra(self):
        r"""
        Return the group algebra of cubic braid group attached to ``self`` over
        the base ring of ``self``.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: CHA2.cubic_braid_group_algebra()
            Algebra of Cubic Braid group on 2 strands
             over Multivariate Polynomial Ring in u, v, w
             over Integer Ring localized at (w,)
        """
        return self._cubic_braid_group_algebra

    # --------------------------------------------------------------------------
    # creating a CubicHeckeAlgebra as sub-algebra of self on less strands
    # --------------------------------------------------------------------------
    def cubic_hecke_subalgebra(self, nstrands=None):
        r"""
        Return a :class:`CubicHeckeAlgebra` that realizes a sub-algebra
        of ``self`` on the first ``n_strands`` strands.

        INPUT:

        - ``nstrands`` -- integer at least 1 and at most :meth:`strands` giving
          the number of strands for the subgroup; the default is one strand
          less than ``self`` has

        OUTPUT: an instance of this class realizing the sub-algebra

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3, cubic_equation_roots=(3, 4, 5))
            sage: CHA3.cubic_hecke_subalgebra()
            Cubic Hecke algebra on 2 strands
              over Integer Ring localized at (2, 3, 5)
                with cubic equation: h^3 - 12*h^2 + 47*h - 60 = 0
        """
        if nstrands is None:
            nstrands = self._nstrands - 1
        n = self._nstrands
        nstrands = ZZ(nstrands)

        if nstrands >= n or nstrands <= 0:
            raise ValueError('nstrands must be positive and less than %s' % self._nstrands)

        names = self.variable_names()
        if nstrands == self._nstrands - 1 and self._cubic_hecke_subalgebra is not None:
            return self._cubic_hecke_subalgebra

        names_red = names[:nstrands - 1]
        if self.base_ring() == self.base_ring(generic=True):
            SubHeckeAlg = CubicHeckeAlgebra(names=names_red)
        else:
            SubHeckeAlg = CubicHeckeAlgebra(names=names_red,
                                            cubic_equation_parameters=tuple(self._cubic_equation_parameters),
                                            cubic_equation_roots=tuple(self._cubic_equation_roots))

        if nstrands == self._nstrands - 1:
            self._cubic_hecke_subalgebra = SubHeckeAlg
        return SubHeckeAlg

    # --------------------------------------------------------------------------
    # mirror image
    # --------------------------------------------------------------------------
    def mirror_image(self):
        r"""
        Return a copy of ``self`` with the mirrored cubic equation, that is:
        the cubic equation has the inverse roots to the roots with respect
        to ``self``.

        This is needed since the mirror involution of the braid group does
        not factor through ``self`` (considered as an algebra over the base
        ring, just considered as `\ZZ`-algebra). Therefore, the mirror
        involution of an element of ``self`` belongs to ``mirror_image``.

        OUTPUT:

        A cubic Hecke algebra over the same base and extension ring,
        but whose cubic equation is transformed by the mirror involution
        applied to its coefficients and roots.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: ce = CHA2.cubic_equation(); ce
            h^3 - u*h^2 + v*h - w
            sage: CHA2m = CHA2.mirror_image()
            sage: cem =  CHA2m.cubic_equation(); cem
            h^3 + ((-v)/w)*h^2 + u/w*h + (-1)/w
            sage: mi = CHA2.base_ring().mirror_involution(); mi
            Ring endomorphism of Multivariate Polynomial Ring in u, v, w
                                 over Integer Ring localized at (w,)
              Defn: u |--> v/w
                    v |--> u/w
                    w |--> 1/w
            sage: cem == cem.parent()([mi(cf) for cf in ce.coefficients()])
            True

        Note that both cubic Hecke algebras have the same ring of definition
        and identical generic cubic equation::

            sage: cemg = CHA2m.cubic_equation(generic=True)
            sage: CHA2.cubic_equation(generic=True) == cemg
            True
            sage: CHA2.cubic_equation() == cemg
            True
            sage: a, b, c = CHA2.cubic_equation_roots()
            sage: CHA2m.cubic_equation_roots(generic=True) == [a, b, c]
            True
            sage: CHA2m.cubic_equation_roots()
            [((-1)/(-w))*a^2 + (u/(-w))*a + (-v)/(-w),
             ((1/(-w))*a)*b + (1/(-w))*a^2 + ((-u)/(-w))*a,
             (((-1)/(-w))*a)*b]
            sage: ai, bi, ci = _
            sage: ai == ~a, bi == ~b, ci == ~c
            (True, True, True)
            sage: CHA2.extension_ring(generic=True).mirror_involution()
            Ring endomorphism of Multivariate Laurent Polynomial Ring in a, b, c
                                 over Splitting Algebra of x^2 + x + 1
                                   with roots [e3, -e3 - 1] over Integer Ring
              Defn: a |--> a^-1
                    b |--> b^-1
                    c |--> c^-1
                    with map of base ring

        The mirror image can not be obtained for specialized cubic Hecke
        algebras if the specialization does not factor through the mirror
        involution on the ring if definition::

            sage: CHA2s = algebras.CubicHecke(2, cubic_equation_roots=(3, 4, 5))
            sage: CHA2s
            Cubic Hecke algebra on 2 strands
              over Integer Ring localized at (2, 3, 5)
                with cubic equation: h^3 - 12*h^2 + 47*h - 60 = 0

        In the next example it is not clear what the mirror image of ``7``
        should be::

            sage: CHA2s.mirror_image()
            Traceback (most recent call last):
            ...
            RuntimeError: base ring Integer Ring localized at (2, 3, 5)
            does not factor through mirror involution
        """
        base_ring = self.base_ring()
        base_gen = self.base_ring(generic=True)

        base_gen_mirror = base_gen.mirror_involution()
        base_ring_mirror = self._base_ring_mirror
        cepg = self.cubic_equation_parameters(generic=True)
        cerg = self.cubic_equation_roots(generic=True)
        if not base_ring_mirror:
            mirr_paras_gen = [base_gen_mirror(par) for par in cepg]
            mirr_paras = [base_ring(mirr_para) for mirr_para in mirr_paras_gen]
            try:
                base_ring_mirror = base_ring.hom(mirr_paras)
            except (TypeError, ValueError, NotImplementedError):
                raise RuntimeError('base ring %s does not factor through mirror involution' % base_ring)

            # check for involution
            mirr_paras_back = [base_ring_mirror(mirr_para) for mirr_para in mirr_paras]
            if mirr_paras_back != self.cubic_equation_parameters():
                raise RuntimeError('base ring %s does not factor through mirror involution' % base_ring)
            self._base_ring_mirror = base_ring_mirror

        mirror_image = self._mirror_image
        if mirror_image is None:
            extension_ring = self.extension_ring()
            extension_gen = self.extension_ring(generic=True)
            extension_gen_mirror = extension_gen.mirror_involution()

            mirr_paras_gen = [base_gen_mirror(par) for par in cepg]
            mirr_roots_gen = [extension_gen_mirror(root) for root in cerg]

            mirr_paras = tuple([base_ring(par) for par in mirr_paras_gen])
            mirr_roots = tuple([extension_ring(root) for root in mirr_roots_gen])
            n = self._nstrands

            mirror_image = CubicHeckeAlgebra(n, cubic_equation_parameters=mirr_paras, cubic_equation_roots=mirr_roots)

            # go back by involution property
            mirror_image._mirror_image = self
            mirror_image._base_ring_mirror = base_ring_mirror
            mirror_image._ring_of_definition._mirror_ = base_gen_mirror
            mirror_image._is_mirror = True

            self._mirror_image = mirror_image

        return mirror_image

    # --------------------------------------------------------------------------
    # Schur elements
    # --------------------------------------------------------------------------
    def schur_elements(self, generic=False):
        r"""
        Return the list of Schur elements of ``self`` as elements
        of the extension ring of ``self``.

        .. NOTE::

            This method needs ``GAP3`` installed with package ``CHEVIE``.

        INPUT:

        - ``generic`` -- boolean (default: ``False``); if ``True``,
          the element is returned as element of the generic
          extension ring

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)       # optional gap3
            sage: sch_eles = CHA3.schur_elements()    # optional gap3
            sage: sch_eles[6]                         # optional gap3
            (u^3*w + v^3 - 6*u*v*w + 8*w^2)/w^2
        """
        gap3_result = self.chevie().SchurElements()
        GER = self.extension_ring(generic=True)
        generic_result = [GER(s) for s in gap3_result]
        if generic:
            return list(generic_result)
        else:
            ER = self.extension_ring()
            return [ER(s) for s in generic_result]

    # --------------------------------------------------------------------------
    # Schur element
    # --------------------------------------------------------------------------
    def schur_element(self, item, generic=False):
        r"""
        Return a single Schur element of ``self`` as elements
        of the extension ring of ``self``.

        .. NOTE::

            This method needs ``GAP3`` installed with package ``CHEVIE``.

        INPUT:

        - ``item`` -- an element of :class:`AbsIrreducibeRep` to give
          the irreducible representation of ``self`` to which the Schur
          element should be returned
        - ``generic`` -- boolean (default: ``False``); if ``True``,
          the element is returned as element of the generic
          extension ring

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)                 # optional gap3
            sage: CHA3.schur_element(CHA3.irred_repr.W3_111)    # optional gap3
            (u^3*w + v^3 - 6*u*v*w + 8*w^2)/w^2
        """
        if not isinstance(item, AbsIrreducibeRep):
            raise ValueError('item must be an instance of %s' % AbsIrreducibeRep)
        return self.schur_elements(generic=generic)[item.gap_index()]

    # --------------------------------------------------------------------------
    # characters
    # --------------------------------------------------------------------------
    def characters(self, irr=None, original=True):
        r"""
        Return the irreducible characters of ``self``.

        By default the values are given in the generic extension ring.
        Setting the keyword ``original`` to ``False`` you can obtain
        the values in the (non generic) extension ring (compare the
        same keyword for :meth:`CubicHeckeElement.matrix`).

        INPUT:

        - ``irr`` -- (optional) instance of :class:`AbsIrreducibeRep`
          selecting the irreducible representation corresponding to the
          character; if not given a list of all characters is returned
        - ``original`` -- boolean (default: ``True``); see description above

        OUTPUT:

        Function or list of Functions from the element class of ``self`` to
        the (generic or non generic) extension ring depending on the given
        keyword arguments.

        EXAMPLES::

            sage: CHA3 = algebras.CubicHecke(3)
            sage: ch = CHA3.characters()
            sage: e = CHA3.an_element()
            sage: ch[0](e)
            a^2*b + a^2*c + a^2 - b*c + b^-1*c^-1 + a^-1*c^-1 + a^-1*b^-1
            sage: _.parent()
            Multivariate Laurent Polynomial Ring in a, b, c
              over Splitting Algebra of x^2 + x + 1 with roots [e3, -e3 - 1]
              over Integer Ring
            sage: ch_w3_100 = CHA3.characters(irr=CHA3.irred_repr.W3_100)
            sage: ch_w3_100(e) == ch[0](e)
            True
            sage: ch_x = CHA3.characters(original=False)
            sage: ch_x[0](e)
            (u + v)*a + (-v*w - w^2 + u)/w
            sage: _.parent()
            Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
              over Splitting Algebra of h^3 - u*h^2 + v*h - w
                with roots [a, b, -b - a + u]
              over Multivariate Polynomial Ring in u, v, w
              over Integer Ring localized at (w,)
        """
        def char_function(ele):
            if isinstance(ele, self.element_class):
                m = ele.matrix(original=original)
                return m[irr].trace()
        if irr:
            return char_function
        irrs = [irr for irr in self.irred_repr if irr.number_gens() == self._nstrands - 1]
        return [self.characters(irrs[i], original=original) for i in range(len(irrs))]
