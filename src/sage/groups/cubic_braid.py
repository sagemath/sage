r"""
Cubic Braid Groups

This module is devoted to factor groups of the Artin braid groups, such
that the images `s_i` of the braid generators have order three:

.. MATH::

    s_i^3 = 1.

In general these groups have firstly been investigated by Coxeter, H.S.M.
in: "Factor groups of the braid groups, Proceedings of the Fourth Canadian
Mathematical Congress (Vancouver 1957), pp. 95-122".

Coxeter showed, that these groups are finite as long as the number of
strands is less than 6 and infinite else-wise. More explicitly the factor
group on three strand braids is isomorphic to `SL(2,3)`, on four strand
braids to `GU(3,2)` and on five strand braids to `Sp(4,3)  \times C_3`.
Today, these finite groups are known as irreducible complex reflection groups
enumerated in the Shephard-Todd classification as `G_{4}`, `G_{25}` and
`G_{32}`.

Coxeter realized these groups as subgroups of unitary groups with respect
to a certain Hermitian form over the complex numbers (in fact over `\QQ`
adjoined with a primitive `12`-th root of unity).

In "Einige endliche Faktorgruppen der Zopfgruppen" (Math. Z., 163 (1978),
291-302) J. Assion considered two series `S(m)` and `U(m)` of finite
factors of these groups. The additional relations on the braid group
generators `\{ s_1, \cdots, s_{m-1}\}` are

.. MATH::

    \begin{array}{lll}
    \mbox{S:} & s_3 s_1 t_2 s_1 t_2^{-1} t_3 t_2 s_1 t_2^{-1} t_3^{-1} = 1
              & \mbox{ for } m >= 5 \mbox{ in case } S(m)\\
    \mbox{U:} & t_1 t_3 = 1
              & \mbox{ for } m >= 5 \mbox{ in case } U(m)
    \end{array}

where `t_i = (s_i s_{i+1})^3`. He showed that each series of finite cubic
braid group factors must be an epimorphic image of one of his two series,
as long as the groups with less than 5 strands are the full cubic braid
groups, whereas the group on 5 strands is not. He realized the groups `S(m)`
as symplectic groups over `GF(3)` (resp. subgroups therein) and `U(m)` as
general unitary groups over `GF(4)` (resp. subgroups therein).

All the groups considered by Coxeter and Assion are considered as finitely
presented groups together with the classical realizations. It also allows
for the conversion maps between the two realizations. In addition, we can
construct other realizations and maps to matrix groups with help of the
Burau representation. In case ``gap3`` and ``CHEVIE`` are installed, the
reflection groups (via the ``gap3`` interface) are available, too. This can
be done using the methods :meth:`as_classical_group`, :meth:`as_matrix_group`,
:meth:`as_permutation_group`, and :meth:`as_reflection_group`.

TESTS::

    sage: CubicBraidGroup(4).category()
    Category of shephard groups
    sage: CubicBraidGroup(6).category()
    Category of infinite groups

REFERENCES:

- [Cox1957]_
- [Ass1978]_

AUTHORS:

- Sebastian Oehms 2019-02-16, initial version.
"""
# ****************************************************************************
#       Copyright (C) 2019 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from enum import Enum

import sage.rings.abc

from sage.categories.groups import Groups
from sage.categories.shephard_groups import ShephardGroups
from sage.groups.free_group import FreeGroup
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement
from sage.groups.braid import BraidGroup
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer

try:
    from sage.libs.gap.element import GapElement
except ImportError:
    GapElement = ()


##############################################################################
#
#   helper functions
#
##############################################################################

def _reduce_tietze(tietze_list):
    r"""
    Reduce the length of a list representing a cubic braid as much as it is
    easily possible using the second braid relation and degree reduction.

    EXAMPLES::

        sage: from sage.groups.cubic_braid import _reduce_tietze
        sage: _reduce_tietze((2, 2, -3, 5, 3, 1, 1, 5))
        [-2, -5, -1]
    """
    def eliminate_item(tietze_list):
        """
        This sub method searches for an item in the Tietze expression such
        that it together with the first entry gives a pair which can be
        replaced by the second braid relation and the generators degree
        reduction. If no such pair exists, it returns ``None``. Otherwise
        the reduced tietze list is returned.
        """
        l = len(tietze_list)
        if l < 2:
            return None
        first = tietze_list[0]
        second = None
        for i in range(1, l):
            if tietze_list[i] in (first, -first):
                if i == 1:
                    second = tietze_list[i]
                    break
                if all(abs(abs(tietze_list[j])-abs(first)) > 1 for j in range(1, i)):
                    # the entry on position i can be moved right to the first entry
                    # by the second braid relation
                    second = tietze_list[i]
                    break
        if second is None:
            return None
        middle = tietze_list[1:i]
        end = tietze_list[i+1:l]
        if first == second:
            return [-first] + middle + end
        else:
            return middle + end

    tietze_list = list(tietze_list)
    l = len(tietze_list)
    for i in range(l):
        end = tietze_list[i:l]
        tietze_list_red = eliminate_item(end)
        if tietze_list_red is not None:
            start = tietze_list[:i]
            return start + _reduce_tietze(tietze_list_red)
    return tietze_list


##############################################################################
#
#   Functions to create Instances of the CubicBraidGroup
#
##############################################################################


# ----------------------------------------------------------------------------------
# Short-Hands for Assions groups:
# ----------------------------------------------------------------------------------
def AssionGroupS(n=None, names='s'):
    r"""
    Construct cubic braid groups :class:`CubicBraidGroup` which have been
    investigated by J.Assion using the notation S(m). This function is a short hand cut
    for setting the construction arguments ``cbg_type=CubicBraidGroup.type.AssionS``
    and default ``names='s'``.

    INPUT:

    - ``n`` -- integer (optional); the number of strands
    - ``names`` -- (default: ``'s'``) string or list/tuple/iterable of strings

    .. SEEALSO::

        :class:`CubicBraidGroup`

    EXAMPLES::

        sage: S3 = AssionGroupS(3);  S3
        Assion group on 3 strands of type S
        sage: S3x = CubicBraidGroup(3, names='s', cbg_type=CubicBraidGroup.type.AssionS); S3x
        Assion group on 3 strands of type S
        sage: S3 == S3x
        True
    """
    return CubicBraidGroup(n=n, names=names, cbg_type=CubicBraidGroup.type.AssionS)


def AssionGroupU(n=None, names='u'):
    r"""
    Construct cubic braid groups as instance of :class:`CubicBraidGroup` which have been
    investigated by J.Assion using the notation U(m). This function is a short hand cut
    for setting the construction arguments ``cbg_type=CubicBraidGroup.type.AssionU``
    and default ``names='u'``.

    INPUT:

    - ``n`` -- integer (optional); the number of strands
    - ``names`` -- (default: ``'s'``) string or list/tuple/iterable of strings

    .. SEEALSO::

        :class:`CubicBraidGroup`

    EXAMPLES::

        sage: U3 = AssionGroupU(3);  U3
        Assion group on 3 strands of type U
        sage: U3x = CubicBraidGroup(3, names='u', cbg_type=CubicBraidGroup.type.AssionU); U3x
        Assion group on 3 strands of type U
        sage: U3 == U3x
        True
    """
    return CubicBraidGroup(n=n, names=names, cbg_type=CubicBraidGroup.type.AssionU)


##############################################################################
#
#                  Class CubicBraidElement (for elements)
#
##############################################################################

class CubicBraidElement(FinitelyPresentedGroupElement):
    r"""
    Elements of cubic factor groups of the braid group.

    For more information see :class:`CubicBraidGroup`.

    EXAMPLES::

        sage: C4.<c1, c2, c3> = CubicBraidGroup(4); C4
        Cubic Braid group on 4 strands
        sage: ele1 = c1*c2*c3^-1*c2^-1
        sage: ele2 = C4((1, 2, -3, -2))
        sage: ele1 == ele2
        True
    """
    def __init__(self, parent, x, check=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C6 = CubicBraidGroup(6)
            sage: C6.inject_variables()
            Defining c0, c1, c2, c3, c4
            sage: c1**2*~c2*c4*c0**2*c4   # indirect doctest
            c1^-1*c2^-1*c4^-1*c0^-1
        """
        if type(x) in (tuple, list):
            x = tuple(_reduce_tietze(tuple(x)))
        elif isinstance(x, GapElement):
            tietze_list = x.UnderlyingElement().TietzeWordAbstractWord().sage()
            tietze_red = _reduce_tietze(tietze_list)
            if tietze_red != tietze_list:
                x = tuple(tietze_red)
        super().__init__(parent, x, check=check)

    def _richcmp_(self, other, op):
        """
        Rich comparison of ``self`` with ``other``.

        Overwrite comparison since the inherited one from
        :class:`FinitelyPresentedGroupElement` (via Gap) does not terminate
        in the case of more than 5 strands (not only infinite cases).
        On less than 5 strands comparison is not assumed to be deterministic
        (see the :issue:`33498` and :gap:`section 47.3-2 of the
        Gap Reference manual <chap47>`).

        Therefore, the comparison is done via the Burau representation.

        EXAMPLES::

            sage: CBG3 = CubicBraidGroup(3)
            sage: sorted(CBG3)    # indirect doctest
            [(c0*c1^-1)^2, c0*c1^-1*c0, c0^-1*c1*c0^-1, c0^-1*c1^-1*c0,
             c1*c0^-1*c1, c0^-1*c1^-1*c0^-1, c0^-1*c1^-1, c1^-1*c0*c1^-1,
             c0*c1^-1*c0^-1, c0^-1*c1, c0^-1*c1*c0, c0*c1^-1, c1*c0^-1,
             c1^-1*c0^-1, c1^-1*c0, c1*c0, c0^-1, c0*c1*c0^-1, c0*c1*c0,
             c0, c0*c1, c1^-1, c1, 1]

            sage: C6.<c1, c2, c3, c4, c5> = CubicBraidGroup(6)
            sage: ele1 = c1*c2*c3*c4*c5
            sage: ele2 = c1*c2*c4*c5
            sage: ele1 == ele2    # indirect doctest
            False
            sage: ele1 > ele2     # indirect doctest
            True

        TESTS::

            sage: S7 = AssionGroupS(7)
            sage: all(S7(rel).is_one() for rel in S7.relations())
            True
        """
        smat = self._matrix_()
        omat = other._matrix_()
        return smat._richcmp_(omat, op)

    def __hash__(self):
        r"""
        Return a hash value for ``self``.

        EXAMPLES::

            sage: C3.<c1, c2> = CubicBraidGroup(3)
            sage: hash(~c1) == hash(c1**2)
            True
        """
        return hash(self._matrix_())

    @cached_method
    def _matrix_(self):
        r"""
        Return ``self`` as a matrix according to its parent's default matrix
        group.

        So far, this method returns the same results as :meth:`burau_matrix`
        in the default case.

        EXAMPLES::

            sage: S3.<s1, s2> = AssionGroupS(3)
            sage: matrix(S3.an_element())
            [2 1 1]
            [1 0 0]
            [0 1 0]
        """
        mat_grp = self.parent().as_matrix_group()
        return mat_grp(self).matrix()

    def braid(self):
        r"""
        Return the canonical braid preimage of ``self`` as a :class:`Braid`.

        EXAMPLES::

            sage: C3.<c1, c2> = CubicBraidGroup(3)
            sage: c1.parent()
            Cubic Braid group on 3 strands
            sage: c1.braid().parent()
            Braid group on 3 strands
        """
        braid_group = self.parent().braid_group()
        return braid_group(self)

    @cached_method
    def burau_matrix(self, root_bur=None, domain=None, characteristic=None,
                     var='t', reduced=False):
        r"""
        Return the Burau matrix of the cubic braid coset.

        This method uses the same method belonging to :class:`Braid`, but
        reduces the indeterminate to a primitive sixth (resp. twelfth in case
        ``reduced='unitary'``) root of unity.

        INPUT (all arguments are optional keywords):

        - ``root_bur`` -- six (resp. twelfth) root of unity in some field
          (default: root of unity over `\QQ`)
        - ``domain`` -- (default: cyclotomic field of order 3 and degree 2, resp.
          the domain of `root_bur` if given) base ring for the Burau matrix
        - ``characteristic`` -- integer giving the characteristic of the
          domain (default: 0 or the characteristic of ``domain`` if given)
        - ``var`` -- string used for the indeterminate name in case ``root_bur``
          must be constructed in a splitting field
        - ``reduced`` -- boolean or string (default: ``False``); for more
          information see the documentation of :meth:`burau_matrix` of
          :class:`Braid`

        OUTPUT:

        The Burau matrix of the cubic braid coset with entries in the
        domain given by the options. In case the option ``reduced='unitary'``
        is given a triple consisting of the Burau matrix, its adjoined and
        the Hermitian form is returned.

        EXAMPLES::

            sage: C3.<c1, c2> = CubicBraidGroup(3)
            sage: ele = c1*c2*c1

            sage: # needs sage.rings.number_field
            sage: BuMa = ele.burau_matrix(); BuMa
            [  -zeta3         1     zeta3]
            [  -zeta3 zeta3 + 1         0]
            [       1         0         0]
            sage: BuMa.base_ring()
            Cyclotomic Field of order 3 and degree 2
            sage: BuMa == ele.burau_matrix(characteristic=0)
            True
            sage: BuMa = ele.burau_matrix(domain=QQ); BuMa
            [-t + 1      1  t - 1]
            [-t + 1      t      0]
            [     1      0      0]
            sage: BuMa.base_ring()
            Number Field in t with defining polynomial t^2 - t + 1
            sage: BuMa = ele.burau_matrix(domain = QQ[I, sqrt(3)]); BuMa                # needs sage.symbolic
            [ 1/2*sqrt3*I + 1/2                  1 -1/2*sqrt3*I - 1/2]
            [ 1/2*sqrt3*I + 1/2 -1/2*sqrt3*I + 1/2                  0]
            [                 1                  0                  0]
            sage: BuMa.base_ring()                                                      # needs sage.symbolic
            Number Field in I with defining polynomial x^2 + 1 over its base field

            sage: BuMa = ele.burau_matrix(characteristic=7); BuMa
            [3 1 4]
            [3 5 0]
            [1 0 0]
            sage: BuMa.base_ring()
            Finite Field of size 7

            sage: # needs sage.rings.finite_rings
            sage: BuMa = ele.burau_matrix(characteristic=2); BuMa
            [t + 1     1 t + 1]
            [t + 1     t     0]
            [    1     0     0]
            sage: BuMa.base_ring()
            Finite Field in t of size 2^2
            sage: F4.<r64> = GF(4)
            sage: BuMa = ele.burau_matrix(root_bur=r64); BuMa
            [r64 + 1       1 r64 + 1]
            [r64 + 1     r64       0]
            [      1       0       0]
            sage: BuMa.base_ring()
            Finite Field in r64 of size 2^2
            sage: BuMa = ele.burau_matrix(domain=GF(5)); BuMa
            [2*t + 2       1 3*t + 3]
            [2*t + 2 3*t + 4       0]
            [      1       0       0]
            sage: BuMa.base_ring()
            Finite Field in t of size 5^2

            sage: # needs sage.rings.number_field
            sage: BuMa, BuMaAd, H = ele.burau_matrix(reduced='unitary'); BuMa
            [       0 zeta12^3]
            [zeta12^3        0]
            sage: BuMa * H * BuMaAd == H
            True
            sage: BuMa.base_ring()
            Cyclotomic Field of order 12 and degree 4
            sage: BuMa, BuMaAd, H = ele.burau_matrix(domain=QQ[I, sqrt(3)],             # needs sage.symbolic
            ....:                                    reduced='unitary'); BuMa
            [0 I]
            [I 0]
            sage: BuMa.base_ring()                                                      # needs sage.symbolic
            Number Field in I with defining polynomial x^2 + 1 over its base field
        """
        braid = self.braid()

        from sage.misc.functional import cyclotomic_polynomial
        min_pol_root_bur = cyclotomic_polynomial(6, var=var)
        unitary = False
        if isinstance(reduced, str):
            if reduced == 'unitary':
                unitary = True
                min_pol_root_bur = cyclotomic_polynomial(12, var=var)

        burau_ori = braid.burau_matrix(reduced=reduced)

        if unitary:
            burau_ori, burau_ori_adj, herm_form_ori = burau_ori

        if domain is not None:
            if isinstance(domain, sage.rings.abc.UniversalCyclotomicField):
                if root_bur is None:
                    if unitary:
                        root_bur = domain.gen(12)
                    else:
                        root_bur = domain.gen(6)

        if root_bur is None:
            def find_root(domain):
                min_pol = min_pol_root_bur.change_ring(domain)
                root_list = min_pol.roots()
                if not root_list:
                    domain = min_pol.splitting_field(min_pol_root_bur.variable_name())
                    min_pol = min_pol_root_bur.change_ring(domain)
                root_list = min_pol.roots()
                for root in root_list:
                    if root[0] == 0:
                        continue
                    root_bur = root[0]
                    if root[1] == 1:
                        break
                return root_bur

            if domain is None:
                if characteristic is None:
                    # --------------------------------------------------------------------
                    # setting the default characteristic in order to achieve the according
                    # representations being well defined
                    # --------------------------------------------------------------------
                    cbg_type = self.parent()._cbg_type
                    if cbg_type == CubicBraidGroup.type.AssionS:
                        characteristic = 3  # making Assion type S relations vanish
                    elif cbg_type == CubicBraidGroup.type.AssionU:
                        characteristic = 2  # making Assion type U relations vanish
                    else:
                        characteristic = 0
                try:
                    characteristic = Integer(characteristic)
                except ValueError:
                    raise ValueError('characteristic must be in integer')

                if not characteristic.is_zero() and not characteristic.is_prime():
                    raise ValueError('characteristic must be a prime')
                if characteristic.is_zero():
                    from sage.rings.number_field.number_field import CyclotomicField
                    if unitary:
                        domain = CyclotomicField(12)
                    else:
                        domain = CyclotomicField(3)
                else:
                    from sage.rings.finite_rings.finite_field_constructor import GF
                    domain = GF(characteristic)
                root_bur = find_root(domain)
                domain = root_bur.parent()

            else:  # domain is not None
                root_bur = find_root(domain)

        else:  # root_bur is not None
            if domain is None:
                domain = root_bur.parent()

            if 1 not in domain:
                raise ValueError('root_bur must belong to a domain containing 1')

            min_pol_root_bur = min_pol_root_bur.change_ring(domain)
            if not min_pol_root_bur(root_bur).is_zero():
                raise ValueError('root_bur must vanish on %s' % (min_pol_root_bur))

        def conv2domain(laur_pol):
            l1, l2 = laur_pol.polynomial_construction()
            p1 = l1.change_ring(domain)
            p2 = root_bur**(l2)
            res = p1(root_bur) * p2
            return res

        from sage.matrix.constructor import matrix

        d1, d2 = burau_ori.dimensions()
        burau_mat = matrix(d1, d2, lambda i, j: conv2domain(burau_ori[i, j]))

        if unitary:
            burau_mat_adj = matrix(d1, d2, lambda i, j: conv2domain(burau_ori_adj[i, j]))
            herm_form = matrix(d1, d2, lambda i, j: conv2domain(herm_form_ori[i, j]))
            return burau_mat, burau_mat_adj, herm_form

        return burau_mat


##############################################################################
#
#                  Class CubicBraidGroup
#
##############################################################################
class CubicBraidGroup(FinitelyPresentedGroup):
    r"""
    Factor groups of the Artin braid group mapping their generators to elements
    of order 3.

    These groups are implemented as a particular case of finitely presented
    groups similar to the :class:`BraidGroup_class`.

    A cubic braid group can be created by giving the number of strands, and
    the name of the generators in a similar way as it works for the
    :class:`BraidGroup_class`.

    INPUT:

    - ``names`` -- see the corresponding documentation of :class:`BraidGroup_class`
    - ``cbg_type`` -- (default: ``CubicBraidGroup.type.Coxeter``;
      see explanation below) enum type :class:`CubicBraidGroup.type`

    Setting the keyword ``cbg_type`` to one on the values
    ``CubicBraidGroup.type.AssionS`` or ``CubicBraidGroup.type.AssionU``,
    the additional relations due to Assion are added:

    .. MATH::

        \begin{array}{lll}
        \mbox{S:} & s_3 s_1 t_2 s_1 t_2^{-1} t_3 t_2 s_1 t_2^{-1} t_3^{-1} = 1
                  & \mbox{ for } m >= 5 \mbox{ in case } S(m), \\
        \mbox{U:} & t_1 t_3 = 1
                  & \mbox{ for } m >= 5 \mbox{ in case } U(m),
        \end{array}

    where `t_i = (s_i s_{i+1})^3`. If ``cbg_type == CubicBraidGroup.type.Coxeter``
    (default) only the cubic relation on the generators is active (Coxeter's
    case of investigation). Note that for `n = 2, 3, 4`, the groups do not
    differ between the three possible values of ``cbg_type`` (as finitely
    presented groups). However, the ``CubicBraidGroup.type.Coxeter``,
    ``CubicBraidGroup.type.AssionS`` and ``CubicBraidGroup.type.AssionU``
    are different, so they have different classical realizations implemented.

    .. SEEALSO::

        Instances can also be constructed more easily by using
        :func:`CubicBraidGroup`, :func:`AssionGroupS` and :func:`AssionGroupU`.

    EXAMPLES::

        sage: U3 = CubicBraidGroup(3, cbg_type=CubicBraidGroup.type.AssionU); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (c0, c1)

    Alternative possibilities defining ``U3``::

        sage: U3 = AssionGroupU(3); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (u0, u1)
        sage: U3.<u1,u2> = AssionGroupU(3); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (u1, u2)

    Alternates naming the generators::

        sage: U3 = AssionGroupU(3, 'a, b'); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (a, b)
        sage: C3 = CubicBraidGroup(3, 't'); C3
        Cubic Braid group on 3 strands
        sage: C3.gens()
        (t0, t1)
        sage: U3.is_isomorphic(C3)
        #I  Forcing finiteness test
        True
        sage: U3.as_classical_group()
        Subgroup generated by [(1,7,6)(3,19,14)(4,15,10)(5,11,18)(12,16,20),
                               (1,12,13)(2,15,19)(4,9,14)(5,18,8)(6,21,16)]
         of (The projective general unitary group of degree 3 over Finite Field of size 2)
        sage: C3.as_classical_group()
        Subgroup with 2 generators (
        [  E(3)^2        0]  [       1 -E(12)^7]
        [-E(12)^7        1], [       0   E(3)^2]
        ) of General Unitary Group of degree 2 over Universal Cyclotomic Field
        with respect to positive definite hermitian form
        [-E(12)^7 + E(12)^11                  -1]
        [                 -1 -E(12)^7 + E(12)^11]

    REFERENCES:

    - [Cox1957]_
    - [Ass1978]_
    """

    Element = CubicBraidElement

    ##############################################################################
    #                  Enum for the type of the group
    ##############################################################################

    class type(Enum):
        r"""
        Enum class to select the type of the group:

        - ``Coxeter`` -- ``'C'`` the full cubic braid group
        - ``AssionS`` -- ``'S'`` finite factor group of type S considered by Assion
        - ``AssionU`` -- ``'U'`` finite factor group of type U considered by Assion

        EXAMPLES::

            sage: S2 = CubicBraidGroup(2, cbg_type=CubicBraidGroup.type.AssionS); S2
            Assion group on 2 strands of type S
            sage: U3 = CubicBraidGroup(2, cbg_type='U')
            Traceback (most recent call last):
            ...
            TypeError: the cbg_type must be an instance of <enum 'CubicBraidGroup.type'>
        """
        Coxeter = 'C'
        AssionS = 'S'
        AssionU = 'U'

    ###########################################################################################
    # private methods
    ###########################################################################################
    @staticmethod
    def __classcall_private__(cls, n=None, names='c', cbg_type=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: C3 = CubicBraidGroup(3); C3.generators()
            (c0, c1)
            sage: CubicBraidGroup(3, 'g').generators()
            (g0, g1)
            sage: U3.<u1,u2>=CubicBraidGroup(3, cbg_type=CubicBraidGroup.type.AssionU); U3.generators()
            (u1, u2)
        """
        # this code is adapted from :func:`BraidGroup`
        # Support Freegroup('a,b') syntax
        if n is not None:
            try:
                n = Integer(n)-1
            except TypeError:
                names = n

                n = None
        # derive n from counting names
        if n is None:
            if isinstance(names, str):
                n = len(names.split(','))
            else:
                names = list(names)
                n = len(names)

        from sage.structure.category_object import normalize_names
        names = tuple(normalize_names(n, names))
        return super().__classcall__(cls, names, cbg_type=cbg_type)

    def __init__(self, names, cbg_type=None):
        r"""
        Python constructor.

        TESTS::

            sage: C3 = CubicBraidGroup(3)    # indirect doctest
            sage: TestSuite(C3).run()
            sage: C4 = CubicBraidGroup(4)    # indirect doctest
            sage: TestSuite(C4).run()        # long time
            sage: C6 = CubicBraidGroup(6)    # indirect doctest
            sage: TestSuite(C6).run()        # long time
            sage: S3 = AssionGroupS(3)       # indirect doctest
            sage: TestSuite(S3).run()
            sage: S5 = AssionGroupS(5)       # indirect doctest
            sage: TestSuite(S5).run()        # long time
            sage: U3 = AssionGroupU(3)       # indirect doctest
            sage: TestSuite(U3).run()
            sage: U4 = AssionGroupU(4)       # indirect doctest
            sage: TestSuite(U4).run()        # long time
            sage: U5 = AssionGroupU(5)       # indirect doctest
            sage: TestSuite(U5).run()        # long time
        """
        n = Integer(len(names))
        if n < 1:
            raise ValueError("the number of strands must be an integer larger than one")

        if cbg_type is None:
            cbg_type = CubicBraidGroup.type.Coxeter
        if not isinstance(cbg_type, CubicBraidGroup.type):
            raise TypeError("the cbg_type must be an instance of %s" % CubicBraidGroup.type)

        free_group = FreeGroup(names)
        self._cbg_type = cbg_type
        self._nstrands = n + 1
        self._ident = self._cbg_type.value + self._nstrands.str()
        self._braid_group = BraidGroup(names)

        # internal naming of elements for convenience
        b = [free_group([i]) for i in range(1, n+1)]
        t = [free_group([i, i+1]) ** 3 for i in range(1, n)]
        ti = [free_group([-i, -i-1]) ** 3 for i in range(1, n)]

        # first the braid relations
        rels = list(self._braid_group.relations())

        # than the cubic relations
        rels.extend(b[i]**3 for i in range(n))

        # than Assion's relation Satz 2.2 for cbg_type=CubicBraidGroup.type.AssionS
        # and Satz 2.4 for cbg_type=CubicBraidGroup.type.AssionU
        if n > 3:
            for i in range(n - 3):
                if cbg_type == CubicBraidGroup.type.AssionU:
                    rels.append((t[i]*t[i+2])**3)
                elif cbg_type == CubicBraidGroup.type.AssionS:
                    rels.append(b[i+2]*b[i]*t[i+1]*b[i]*ti[i+1]*t[i+2]*t[i+1]*b[i]*ti[i+1]*ti[i+2])

        if cbg_type != CubicBraidGroup.type.Coxeter:
            cat = Groups().Finite()
        elif self._nstrands <= 5:
            cat = ShephardGroups()
        else:
            cat = Groups().Infinite()
        FinitelyPresentedGroup.__init__(self, free_group, tuple(rels), category=cat)
        self._free_group = free_group

        # ------------------------------------------------------------------------------------------------
        # the following global pointers to classical group realizations will be set in the private method
        # _create_classical_realization
        # ------------------------------------------------------------------------------------------------
        self._classical_group = None   # This is the classical Group returned by as_classical_group
        self._classical_base_group = None   # this only differs for special cases for Assion groups from the former
        self._classical_invariant_form = None   # invariant form of the classical base group
        self._classical_embedding = None   # if self._classical_group different from self._classical_base_group
        self._centralizing_matrix = None   # for Assion groups: element in classical base group commuting with self
        self._centralizing_element = None   # image under nat. map of the former one in the proj. classical group
        return

    def _repr_(self):
        r"""
        Return a string representation.

        EXAMPLES::

            sage: CubicBraidGroup(2)
            Cubic Braid group on 2 strands
            sage: AssionGroupU(2)
            Assion group on 2 strands of type U
        """
        if self._cbg_type == CubicBraidGroup.type.Coxeter:
            return "Cubic Braid group on %s strands" % (self.strands())
        else:
            return "Assion group on %s strands of type %s" % (self.strands() ,self._cbg_type.value)

    def index_set(self):
        r"""
        Return the index set of ``self``.

        This is the set of integers `0,\dots,n-2` where `n` is
        the number of strands.

        This is only used when ``self`` is a finite reflection group.

        EXAMPLES::

            sage: CubicBraidGroup(3).index_set()
            [0, 1]
        """
        return list(range(self.strands() - 1))

    def simple_reflections(self):
        """
        Return the generators of ``self``.

        This is only used when ``self`` is a finite reflection group.

        EXAMPLES::

            sage: CubicBraidGroup(3).simple_reflections()
            (c0, c1)
        """
        return self.generators()

    def degrees(self):
        """
        Return the degrees of ``self``.

        This only makes sense when ``self`` is a finite reflection group.

        EXAMPLES::

            sage: CubicBraidGroup(4).degrees()
            (6, 9, 12)
        """
        if self._cbg_type != CubicBraidGroup.type.Coxeter:
            raise TypeError('not a finite reflection group')
        if self.strands() > 5:
            raise TypeError('not a finite reflection group')
        d_table = {1: (), 2: (3,), 3: (4, 6),
                   4: (6, 9, 12), 5: (12, 18, 24, 30)}
        return tuple(Integer(deg) for deg in d_table[self.strands()])

    def codegrees(self):
        """
        Return the codegrees of ``self``.

        This only makes sense when ``self`` is a finite reflection group.

        EXAMPLES::

            sage: CubicBraidGroup(5).codegrees()
            (0, 6, 12, 18)
        """
        if self._cbg_type != CubicBraidGroup.type.Coxeter:
            raise TypeError('not a finite reflection group')
        if self.strands() > 5:
            raise TypeError('not a finite reflection group')
        d_table = {1: (), 2: (0,), 3: (0, 2),
                   4: (0, 3, 6), 5: (0, 6, 12, 18)}
        return tuple(Integer(deg) for deg in d_table[self.strands()])

    # -------------------------------------------------------------------------------
    # Methods for test_suite
    # -------------------------------------------------------------------------------

    def _internal_test_attached_group(self, attached_group, tester):
        r"""
        Test conversion maps from ``self`` to the given attached Group,
        which must have been defined using the :meth:`as_classical_group`,
        :meth:`as_matrix_group`, :meth:`as_permutation_group` or
        :meth:`as_reflection_group`.

        INPUT:

        - ``attached_group`` -- attached group to be tested as specified above

        EXAMPLES::

            sage: CBG2 = CubicBraidGroup(2)
            sage: tester = CBG2._tester()
            sage: CBG2M = CBG2.as_matrix_group()
            sage: CBG2._internal_test_attached_group(CBG2M, tester)
        """
        elem = self.an_element()
        att_grp_elem = attached_group(elem)
        if self.is_finite() and self.strands() <= 7: # not realistic for larger number of strands
            att_grp_elem_back = self(att_grp_elem)
            tester.assertEqual(att_grp_elem_back, elem)
        return

    def _test_classical_group(self, **options):
        r"""
        Check the classical group properties.

        The following is checked:

        - Construction of classical group was faithful.
        - Coercion maps to and from classical group exist and are
          inverse to each other.

        EXAMPLES::

            sage: CBG2 = CubicBraidGroup(2)
            sage: CBG2._test_classical_group()
        """
        tester = self._tester(**options)
        classic_grp = self.as_classical_group()
        if self.is_finite():
            self._internal_test_attached_group(classic_grp, tester)
        return

    def _test_permutation_group(self, **options):
        r"""
        Check the permutation group properties.

        The following is checked:

        - Construction of permutation group was faithful.
        - Coercion maps to and from permutation group exist and are
          inverse to each other.

        EXAMPLES::

            sage: CBG2 = CubicBraidGroup(2)
            sage: CBG2._test_permutation_group()
        """
        if self.is_finite():
            tester = self._tester(**options)
            permgrp = self.as_permutation_group()
            self._internal_test_attached_group(permgrp, tester)
        return

    def _test_matrix_group(self, **options):
        r"""
        Check the matrix group properties.

        The following is checked:

        - Construction of matrix group was faithful.
        - Coercion maps to and from matrix group exist and are
          inverse to each other.

        EXAMPLES::

            sage: CBG2 = CubicBraidGroup(2)
            sage: CBG2._test_matrix_group()
        """
        tester = self._tester(**options)

        MatDEF = self.as_matrix_group()
        self._internal_test_attached_group(MatDEF, tester)

        try:
            from sage.rings.finite_rings.finite_field_constructor import GF
        except ImportError:
            return

        F3 = GF(3)
        r63 = F3(2)
        F4 = GF(4)
        r64 = F4.gen()

        if self._cbg_type != CubicBraidGroup.type.AssionU or self.strands() < 5: # not well defined else-wise
            matrix_grpF3 = self.as_matrix_group(root_bur=r63)
            self._internal_test_attached_group(matrix_grpF3, tester)

        if self._cbg_type != CubicBraidGroup.type.AssionS or self.strands() < 5: # not well defined else-wise
            matrix_grpF4 = self.as_matrix_group(root_bur=r64)
            self._internal_test_attached_group(matrix_grpF4, tester)

        if self.strands() < 5 or self._cbg_type == CubicBraidGroup.type.Coxeter:
            matrix_grpF5 = self.as_matrix_group(characteristic=5)
            self._internal_test_attached_group(matrix_grpF5, tester)

            matrix_grpF7 = self.as_matrix_group(domain=GF(7))
            self._internal_test_attached_group(matrix_grpF7, tester)
        return

    def _test_reflection_group(self, **options):
        r"""
        Check the reflection group properties.

        The following is checked:

        - Construction of reflection group was faithful.
        - Coercion maps to and from reflection group exist and are
          inverse to each other.

        EXAMPLES::

            sage: CBG2 = CubicBraidGroup(2)
            sage: CBG2._test_reflection_group()
        """
        if self._cbg_type == CubicBraidGroup.type.Coxeter and self.is_finite() and self.strands() > 2:
            from sage.combinat.root_system.reflection_group_real import is_chevie_available
            if is_chevie_available():
                tester = self._tester(**options)
                reflgrp = self.as_reflection_group()
                self._internal_test_attached_group(reflgrp, tester)
        return

    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    # local utility-methods
    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    def _create_classical_realization(self, just_embedded=False):
        r"""
        Internal method to create the classical groups attached to ``self``.

        This methods sets the following attributes of ``self``:

         - self._classical_group            This is the classical group returned by as_classical_group method.
         - self._classical_base_group       this only differs in special cases for Assion groups from the former.
         - self._classical_invariant_form   invariant form of the classical base group.
         - self._centralizing_matrix        for Assion groups: element in classical base group commuting with ``self``.
         - self._centralizing_element       image under natural map of the former one in the projective classical group.
         - self._classical_embedding        as subgroup of classical base group (if different from classical group).

        EXAMPLES::

            sage: AU2 = AssionGroupU(2)
            sage: AU2._classical_group is None
            True
            sage: AU2._classical_embedding is None
            True
            sage: AU2._classical_invariant_form is None
            True
            sage: AU2._create_classical_realization()
            sage: AU2._classical_group
            General Unitary Group of degree 1 over Finite Field in a of size 2^2
            sage: AU2._classical_embedding is AU2._classical_group
            True
            sage: AU2._classical_invariant_form
            [1]
        """

        # -------------------------------------------------------------------------------
        # Set up data of the classical Assion group (generic part)
        # -------------------------------------------------------------------------------

        def set_classical_realization(self, base_group, proj_group, centralizing_matrix, transvec_matrices):
            r"""
            Internal method to create classical group for Assion groups.

            This is a local function of :meth:`_create_classical_realization`.

            It handles the common part of symplectic and unitary version and
            creates conversion maps.

            INPUT:

            - ``base_group`` -- the symplectic or unitary groups Sp(m,3) resp. GU(m,2)
            - ``proj_group`` -- the corresponding projective group of base_group
            - ``centralizing_matrix`` -- the centralizing matrix according to Assion
            - ``transvec_matrices`` -- list of transvection matrices according to Assion

            OUTPUT: no output, but the function sets the attributes of ``self`` described above
            """
            centralizing_element = None

            # ------------------------------------------------------------------------------
            # Setting the List of Braid Images
            # ------------------------------------------------------------------------------
            im_gens = [base_group(m) for m in transvec_matrices]

            # ------------------------------------------------------------------------------
            # By the work of Assion no check on the group homomorphism is needed, at all.
            # But to take care of software bugs they are performed in cases where they are
            # not really expansive.
            # ------------------------------------------------------------------------------
            check = False
            if self.strands() < 7:
                check = True

            # ------------------------------------------------------------------------------
            # Do the projective group realization if needed
            # ------------------------------------------------------------------------------
            embedding = self._classical_embedding
            classical_group = None
            if proj_group is None:
                classical_group = base_group
                hom_to_classic = self.hom(im_gens, check=check)
                classical_group.register_conversion(hom_to_classic)
                embedding = classical_group
            else:
                if embedding is None:
                    im_gens.pop()
                    embedding = base_group.subgroup(im_gens, check=check)
                    embedding.register_conversion(self.hom(embedding.gens(), check=check))
                    hom_to_base = self.hom(im_gens, check=check, codomain=base_group)
                    base_group.register_conversion(hom_to_base)
                if not just_embedded:
                    transvec_matrices.pop()
                    nat_hom = base_group.hom(proj_group.gens(), check=check)
                    centralizing_element = nat_hom(centralizing_matrix)
                    classical_group_gens = [nat_hom(m) for m in transvec_matrices]
                    classical_group = proj_group.subgroup(classical_group_gens, canonicalize=False)
                    hom_to_classic = self.hom(classical_group.gens(), check=check)
                    classical_group.register_conversion(hom_to_classic)

            # ------------------------------------------------------------------------------
            # register constructed items
            # ------------------------------------------------------------------------------
            self._classical_group = classical_group
            self._classical_base_group = base_group
            self._classical_invariant_form = base_group.invariant_form()
            self._centralizing_matrix = centralizing_matrix
            self._centralizing_element = centralizing_element
            self._classical_embedding = embedding
            return

        # -------------------------------------------------------------------------------
        # local methods to set up the classical group (specific part)
        # -------------------------------------------------------------------------------
        # Case for symplectic groups
        # -------------------------------------------------------------------------------

        def create_sympl_realization(self, m):
            r"""
            Internal method to create classical group for symplectic
            Assion groups (`cbg_type == CubicBraidGroup.type.AssionS`).

            INPUT:

            - ``m`` -- integer; the dimension of the classical groups
              vector-space of operation

            The function calculates the centralizing matrix and the
            transvections as given by Assion and then uses
            ``set_classical_realization`` to complete the construction.
            """
            # -----------------------------------------------------------
            # getting the invariant bilinear form of the group
            # and setting constants.
            # -----------------------------------------------------------
            n = self.strands()

            from sage.groups.matrix_gps.symplectic import Sp
            base_group = Sp(m, 3)
            proj_group = None
            if m == n:
                from sage.groups.perm_gps.permgroup_named import PSp
                proj_group = PSp(m, 3)

            bform = base_group.invariant_form()
            bas = bform.column_space().basis()

            mhalf = m // 2

            # -----------------------------------------------------------
            # computing a hyperbolic decomposition basis with respect
            # to the invariant bilinear form.
            # -----------------------------------------------------------
            xbas = [bas[mhalf - i - 1] for i in range(mhalf)]
            ybas = [bas[mhalf + i] for i in range(mhalf)]

            # -----------------------------------------------------------
            # computing the List of transvection vectors according to
            # the Assion paper, page 292.
            # -----------------------------------------------------------
            transvections = [xbas[0]]                                   # t_1      = x_1
            for i in range(mhalf-1):
                transvections.append(ybas[i])                          # t_{2i}   = y_i
                transvections.append(xbas[i] + xbas[i+1])              # t_{2i+1} = x_j + x_(j+1)
            transvections.append(ybas[mhalf-1])                        # t_n      = y_m

            # -----------------------------------------------------------
            # Conversion-Map from transvection vector to transvection
            # matrix.
            # -----------------------------------------------------------
            from sage.matrix.constructor import matrix

            def transvec2mat(v, bas=bas, bform=bform, fact=1):
                t = [x + fact*(x * bform * v) * v for x in bas]
                return matrix(bform.base_ring(),  t)

            # ------------------------------------------------------------------------------
            # setting the centralizing matrix for the case of projective group realization
            # ------------------------------------------------------------------------------
            centralizing_vector = xbas[mhalf-1]
            centralizing_matrix = base_group(transvec2mat(centralizing_vector, fact=1))
            transvec_matrices = [transvec2mat(v) for v in transvections]

            set_classical_realization(self, base_group, proj_group, centralizing_matrix, transvec_matrices)
            return

        # -------------------------------------------------------------------------------
        # Case for unitary groups
        # -------------------------------------------------------------------------------

        def create_unitary_realization(self, m):
            """
            Internal method to create classical group for
            unitary Assion groups (`cbg_type == CubicBraidGroup.type.AssionU`).

            INPUT:

            - ``m`` -- integer; the dimension of the classical groups
              vector-space of operation

            The function calculates the centralizing_matrix and the
            transvections as given by Assion and then uses
            ``set_classical_realization`` to complete the construction.
            """
            # ---------------------------------------------------------------------
            # getting the invariant bilinear form of the group
            # and setting constants
            # ---------------------------------------------------------------------
            n = self.strands()

            from sage.groups.matrix_gps.unitary import GU
            base_group = GU(m, 2)
            proj_group = None
            if m == n:
                from sage.groups.perm_gps.permgroup_named import PGU
                proj_group = PGU(m, 2)

            bform = base_group.invariant_form()
            bas = bform.column_space().basis()
            F = bform.base_ring()
            a = F.gen()

            mthird = m // 3

            # -----------------------------------------------------------
            # computing a orthonormal basis with respect
            # to the invariant bilinear form.
            # -----------------------------------------------------------
            xbas = []
            for i in range(m):
                if 2*i == m-1:
                    xbas.append(bas[i])
                else:
                    xbas.append(a*bas[i] + a.frobenius()*bas[m-1 - i])

            # -----------------------------------------------------------
            # computing the List of transvection vectors according to
            # Assion paper, page 293.
            # -----------------------------------------------------------
            transvections = [xbas[0]]                                          # t_1 = x_1
            if m > 1:
                transvections.append(xbas[0]+xbas[1]+xbas[2])                 # t_2 = x_1 + x_2 + x_3
            for j in range(mthird):
                pos = 3*(j+1)-1
                transvections.append(xbas[pos-1])                             # t_{3i}   = x_{3i-1}
                if pos + 1 < m:
                    transvections.append(xbas[pos-1]+xbas[pos]+xbas[pos+1])   # t_{3i+1} = x_{3i-1} + x_{3i} + x_{3i+1}
                if pos + 3 < m:
                    transvections.append(xbas[pos+1]+xbas[pos+2]+xbas[pos+3])  # t_{3i+2} = x_{3i+1} + x_{3i+2} + x_{3i+3}

            # -----------------------------------------------------------
            # Conversion-Map from transvection vector to transvection
            # matrix.
            # -----------------------------------------------------------
            from sage.matrix.constructor import matrix

            def transvec2mat(v, bas=bas, bform=bform, fact=a):
                # note x does not change under conjugation, since it belongs to standard basis
                t = [x + fact * (x * bform * v.conjugate()) * v for x in bas]
                return matrix(F, t)

            # ------------------------------------------------------------------------------
            # setting the centralizing matrix for the case of projective group realization.
            # ------------------------------------------------------------------------------
            centralizing_vector = xbas[m-2]+xbas[m-1]
            centralizing_matrix = base_group(transvec2mat(centralizing_vector, fact=1))
            transvec_matrices = [transvec2mat(v) for v in transvections]

            set_classical_realization(self, base_group, proj_group, centralizing_matrix, transvec_matrices)
            return

        # ----------------------------------------------------------------
        # local functions declaration section finishes here
        # ----------------------------------------------------------------

        # ----------------------------------------------------------------
        # initialization of constants
        # ----------------------------------------------------------------

        n = self.strands()

        # -------------------------------------------------------------------------------
        # Setting the Classical group
        # -------------------------------------------------------------------------------
        if self._cbg_type == CubicBraidGroup.type.AssionS:
            dim_sympl_group = n-1              # S(n-1) = Sp(n-1, 3)
            if n % 2 == 0:
                dim_sympl_group = n            # S(n-1) = subgroup of PSp(n, 3)
            create_sympl_realization(self, dim_sympl_group)
        elif self._cbg_type == CubicBraidGroup.type.AssionU:
            dim_unitary_group = n-1              # U(n-1) = GU(n-1, 2)
            if n % 3 == 0:
                dim_unitary_group = n            # U(n-1) = subgroup PGU(n, 3)
            create_unitary_realization(self, dim_unitary_group)
        else:
            # -----------------------------------------------------------------------------------------------
            # connection between Coxeter realization and unitary Burau representation according to Squier:
            # -----------------------------------------------------------------------------------------------
            # Notation of Coxeter: p = 3, \theta =\pi/3 = primitive 6th root of unity
            #     i\theta = \pi/3+\pi/2 = 5\pi/6 = 5th power of primitive 12th root of unity
            #     i\theta = z12^5 = - ~z12 where z12 = UCF.gen(12)
            # Let f be the unitary Form of Coxeter and J(s) the one of Squier. Then we have
            #     J(z12) = 2 f
            # Let `buc` be the unitary Burau Matrix of Coxeter for the first braid generator and `bus`
            # the corresponding one according to Squier. Then we have:
            #     buc_[i,i-1] = buc_[i,i+1]= - i\theta, buc_[i,i] = 1 + 2*cos(\pi/6)*i\theta
            #     bus_[i,i-1] = bus_[i,i+1]=   s,       bus_[i,i] = -s^2
            # now 1 + 2*cos(\pi/6)*i\theta = 1 + sqrt(3)*(-sqrt(3)/2 + I/2) = 1- 3/2 + sqrt(3)I/2 = z12^4 = - ~z12^2
            # finally: Coxeter's Realization is the unitary Burau representation of Squier for s = ~z12
            # -----------------------------------------------------------------------------------------------
            from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

            UCF = UniversalCyclotomicField()
            z12 = UCF.gen(12)
            classical_group = self.as_matrix_group(root_bur=~z12, domain=UCF, reduced='unitary')
            self._classical_group = classical_group
            self._classical_base_group = classical_group
            self._classical_embedding = classical_group
            if self._classical_invariant_form is None:
                self._classical_invariant_form = classical_group.ambient().invariant_form()
        return

    def _element_constructor_(self, x, **kwds):
        r"""
        Return an element of ``self``.

        Extensions to the element constructor of :class:`FinitelyPresentedGroup`:
        new functionalities are:

        - constructing element from an element of the attached classical group
          (embedded and not embedded)
        - constructing element from an element of the attached permutation group
        - constructing element from an element of the attached reflection group

        INPUT:

        - ``x`` -- can be one of the following:

          * an instance of the element class of ``self`` (but possible
            to a different parent)
          * an instance of the element class of the braid group
          * a tuple representing a braid in Tietze form
          * an instance of an element class of a parent ``P`` such that there
            is a map from ``self`` to ``P`` having :meth:`lift`; for example,
            an element of an alternative realization of ``self``, such as
            the classical realization
          * any other object which works for the element constructor
            of :class:`FinitelyPresentedGroup`

        EXAMPLES::

            sage: S3 = AssionGroupS(3)
            sage: S3Cl = S3.as_classical_group()
            sage: g = mul(S3Cl.gens())
            sage: S3(g)  # indirect doctest
            s0*s1*s0^-1
        """
        if hasattr(x, 'parent'):
            parent = x.parent()
            map_to = parent.convert_map_from(self)
            if map_to is not None:
                if hasattr(map_to, 'lift'):
                    return map_to.lift(x)
        return super()._element_constructor_(x)

    #######################################################################################################################
    # ----------------------------------------------------------------------------------
    # public methods
    # ----------------------------------------------------------------------------------
    #######################################################################################################################

    # ---------------------------------------------------------------------------------------------------------------------
    # strands
    # ---------------------------------------------------------------------------------------------------------------------
    def strands(self):
        r"""
        Return the number of strands of the braid group whose image is ``self``.

        OUTPUT: :class:`Integer`

        EXAMPLES::

            sage: C4 = CubicBraidGroup(4)
            sage: C4.strands()
            4
        """
        return self._nstrands

    # ----------------------------------------------------------------------------------
    # braid_group
    # ----------------------------------------------------------------------------------
    def braid_group(self):
        r"""
        Return a :class:`BraidGroup` with identical generators, such that
        there exists an epimorphism to ``self``.

        OUTPUT: a :class:`BraidGroup` having conversion maps to and from
        ``self`` (which is just a section in the latter case)

        EXAMPLES::

            sage: U5 = AssionGroupU(5); U5
            Assion group on 5 strands of type U
            sage: B5 = U5.braid_group(); B5
            Braid group on 5 strands
            sage: b = B5([4,3,2,-4,1])
            sage: u = U5([4,3,2,-4,1])
            sage: u == b
            False
            sage: b.burau_matrix()
            [ 1 - t      t      0      0      0]
            [ 1 - t      0      t      0      0]
            [ 1 - t      0      0      0      t]
            [ 1 - t      0      0      1 -1 + t]
            [     1      0      0      0      0]
            sage: u.burau_matrix()
            [t + 1     t     0     0     0]
            [t + 1     0     t     0     0]
            [t + 1     0     0     0     t]
            [t + 1     0     0     1 t + 1]
            [    1     0     0     0     0]
            sage: bU = U5(b)
            sage: uB = B5(u)
            sage: bU == u
            True
            sage: uB == b
            True
        """
        return self._braid_group

    # ----------------------------------------------------------------------------------
    # as_matrix_group
    # ----------------------------------------------------------------------------------
    @cached_method
    def as_matrix_group(self, root_bur=None, domain=None, characteristic=None, var='t', reduced=False):
        r"""
        Create an epimorphic image of ``self`` as a matrix group by use of
        the burau representation.

        INPUT:

        - ``root_bur`` -- (default: root of unity over `\QQ`) six (resp. twelfth)
          root of unity in some field
        - ``domain`` -- (default: cyclotomic field of order 3 and degree 2, resp.
          the domain of ``root_bur`` if given) base ring for the Burau matrix
        - ``characteristic`` -- integer (optional); the characteristic of the
          domain; if none of the keywords ``root_bur``, ``domain`` and
          ``characteristic`` are given, the default characteristic is 3
          (resp. 2) if ``self`` is of ``cbg_type``
          ``CubicBraidGroup.type.AssionS`` (resp. ``CubicBraidGroup.type.AssionU``)
        - ``var`` -- string used for the indeterminate name in case ``root_bur``
          must be constructed in a splitting field
        - ``reduced`` -- boolean (default: ``False``); for more information
          see the documentation of :meth:`Braid.burau_matrix`

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: C5 = CubicBraidGroup(5)
            sage: C5Mch5 = C5.as_matrix_group(characteristic=5); C5Mch5
            Matrix group over Finite Field in t of size 5^2 with 4 generators (
            [2*t + 2 3*t + 4       0       0       0]
            [     1       0       0       0       0]
            [     0       0       1       0       0]
            [     0       0       0       1       0]
            [     0       0       0       0       1],
            <BLANKLINE>
            [     1       0       0       0       0]
            [     0 2*t + 2 3*t + 4       0       0]
            [     0       1       0       0       0]
            [     0       0       0       1       0]
            [     0       0       0       0       1],
            <BLANKLINE>
            [     1       0       0       0       0]
            [     0       1       0       0       0]
            [     0       0 2*t + 2 3*t + 4       0]
            [     0       0       1       0       0]
            [     0       0       0       0       1],
            <BLANKLINE>
            [     1       0       0       0       0]
            [     0       1       0       0       0]
            [     0       0       1       0       0]
            [     0       0       0 2*t + 2 3*t + 4]
            [     0       0       0       1       0]
            )
            sage: c = C5([3,4,-2,-3,1]); c
            c2*c3*c1^-1*c2^-1*c0
            sage: m = C5Mch5(c); m
            [2*t + 2 3*t + 4       0       0       0]
            [     0       0       0       1       0]
            [2*t + 1       0 2*t + 2     3*t 3*t + 3]
            [2*t + 2       0       0 3*t + 4       0]
            [     0       0 2*t + 2 3*t + 4       0]
            sage: m_back = C5(m)
            sage: m_back == c
            True
            sage: U5 = AssionGroupU(5); U5
            Assion group on 5 strands of type U
            sage: U5Mch3 = U5.as_matrix_group(characteristic=3)
            Traceback (most recent call last):
            ...
            ValueError: Burau representation does not factor through the relations
        """
        # ------------------------------------------------------------------
        # define matrix group by generators using the Burau representation
        # ------------------------------------------------------------------
        unitary = False
        if isinstance(reduced, str):
            if reduced == 'unitary':
                unitary = True
        gen_list = []
        for braid_gen in self.gens():
            bur_mat = braid_gen.burau_matrix(root_bur=root_bur, domain=domain,
                                             characteristic=characteristic,
                                             var=var, reduced=reduced)
            if unitary:
                bur_mat, bur_mat_ad, herm_form = bur_mat

            if domain is None:
                domain = bur_mat.base_ring()

            gen_list.append(bur_mat)

        if unitary and herm_form.is_singular():
            unitary = False  # since a degenerated hermitian form doesn't define a unitary group
            if self._classical_invariant_form is None:
                self._classical_invariant_form = herm_form

        if unitary:
            from sage.rings.finite_rings.finite_field_base import FiniteField
            from sage.groups.matrix_gps.unitary import GU
            _, d = herm_form.dimensions()
            if isinstance(domain, FiniteField):
                base_group = GU(d, domain, var=domain.gen(), invariant_form=herm_form)
            else:
                base_group = GU(d, domain, invariant_form=herm_form)

            matrix_group = base_group.subgroup(gen_list)
        else:
            from sage.groups.matrix_gps.finitely_generated import MatrixGroup
            matrix_group = MatrixGroup(gen_list, category=self.category())

        # --------------------------------------------------------------------
        # check if there is a well defined group homomorphism to matrix_group
        # Register map from ``self`` to matrix_group.
        # Since GAP' check is very expansive (on time and memory), the check is performed
        # here.
        # -------------------------------------------------------------------------------
        hom_to_mat = self.hom(matrix_group.gens(), check=False)
        if not all(hom_to_mat(rel).is_one() for rel in self.relations()):
            raise ValueError("Burau representation does not factor through the relations")
        matrix_group.register_conversion(hom_to_mat)
        return matrix_group

    # ----------------------------------------------------------------------------------
    # Although this method is available for finitely presented group
    # we use the classical group implementation (by performance reason) to get
    # the permutation_group.
    # ----------------------------------------------------------------------------------
    @cached_method
    def as_permutation_group(self, use_classical=True):
        r"""
        Return a permutation group isomorphic to ``self`` that has a
        group isomorphism from ``self`` as a conversion.

        INPUT:

        - ``use_classical`` -- boolean (default: ``True``); the permutation
          group is calculated via the attached classical matrix group as this
          results in a smaller degree; if ``False``, the permutation group will
          be calculated using ``self`` (as finitely presented group)

        EXAMPLES::

            sage: C3 = CubicBraidGroup(3)
            sage: PC3 = C3.as_permutation_group()
            sage: assert C3.is_isomorphic(PC3)  # random (with respect to the occurrence of the info message)
            #I  Forcing finiteness test
            sage: PC3.degree()
            8
            sage: c = C3([2,1-2])
            sage: C3(PC3(c)) == c
            True
        """
        if not self.is_finite():
            raise ValueError('cubic braid group is infinite!')

        if use_classical:
            CG = self.as_classical_group()
            from sage.groups.perm_gps.permgroup import PermutationGroup_generic
            if isinstance(CG, PermutationGroup_generic):
                return CG
            CGM = CG.as_matrix_group()
            PG = CGM.as_permutation_group()
            img_gens = [PG(CGM(CG(gen))) for gen in self.gens()]
        else:
            PG = super().as_permutation_group()
            img_gens = PG.gens()

        img_gens = [PG(gen) for gen in img_gens]
        hom_to_perm = self.hom(img_gens)
        PG.register_conversion(hom_to_perm)
        return PG

    # ----------------------------------------------------------------------------------
    # as_classical_group
    # ----------------------------------------------------------------------------------
    def as_classical_group(self, embedded=False):
        r"""
        Create an isomorphic image of ``self`` as a classical group according
        to the construction given by Coxeter resp. Assion.

        INPUT:

        - ``embedded`` -- boolean (default: ``False``); this boolean effects the
          cases of Assion groups when they are realized as projective groups only.
          More precisely: if ``self`` is of ``cbg_type CubicBraidGroup.type.AssionS``
          (for example) and the number of strands ``n`` is even, than its classical
          group is a subgroup of ``PSp(n,3)`` (being centralized by the element
          ``self.centralizing_element(projective=True))``. By default this group
          will be given. Setting ``embedded = True`` the classical realization
          is given as subgroup of its classical enlargement with one more strand
          (in this case as subgroup of ``Sp(n,3))``.

        OUTPUT:

        Depending on the type of ``self`` and the number of strands an
        instance of ``Sp(n-1,3)``, ``GU(n-1,2)``, subgroup of ``PSp(n,3)``,
        ``PGU(n,2)``, or a subgroup of ``GU(n-1, UCF)``
        (``cbg_type == CubicBraidGroup.type.Coxeter``) with respect to a
        certain Hermitian form attached to the Burau representation
        (used by Coxeter and Squier). Here ``UCF`` stands for the universal
        cyclotomic field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: U3 = AssionGroupU(3)
            sage: U3Cl = U3.as_classical_group(); U3Cl
            Subgroup generated by [(1,7,6)(3,19,14)(4,15,10)(5,11,18)(12,16,20),
                                   (1,12,13)(2,15,19)(4,9,14)(5,18,8)(6,21,16)] of
             (The projective general unitary group of degree 3 over Finite Field of size 2)
            sage: U3Clemb = U3.as_classical_group(embedded=True); U3Clemb
            Subgroup with 2 generators (
            [0 0 a]  [a + 1     a     a]
            [0 1 0]  [    a a + 1     a]
            [a 0 a], [    a     a a + 1]
            ) of General Unitary Group of degree 3 over Finite Field in a of size 2^2
            sage: u = U3([-2,1,-2,1]); u
            (u1^-1*u0)^2
            sage: uCl = U3Cl(u); uCl
            (1,16)(2,9)(3,10)(4,19)(6,12)(7,20)(13,21)(14,15)
            sage: uCle = U3Clemb(u); uCle
            [a + 1 a + 1     1]
            [a + 1     0     a]
            [   1     a     a]
            sage: U3(uCl) == u
            True
            sage: U3(uCle) == u
            True
            sage: U4 = AssionGroupU(4)
            sage: U4Cl = U4.as_classical_group(); U4Cl
            General Unitary Group of degree 3 over Finite Field in a of size 2^2
            sage: U3Clemb.ambient() == U4Cl
            True

            sage: # needs sage.rings.number_field
            sage: C4 = CubicBraidGroup(4)
            sage: C4Cl = C4.as_classical_group(); C4Cl
            Subgroup with 3 generators (
            [  E(3)^2        0        0]  [       1 -E(12)^7        0]
            [-E(12)^7        1        0]  [       0   E(3)^2        0]
            [       0        0        1], [       0 -E(12)^7        1],
            <BLANKLINE>
            [       1        0        0]
            [       0        1 -E(12)^7]
            [       0        0   E(3)^2]
            ) of General Unitary Group of degree 3 over Universal Cyclotomic Field
            with respect to positive definite hermitian form
            [-E(12)^7 + E(12)^11                  -1                   0]
            [                 -1 -E(12)^7 + E(12)^11                  -1]
            [                  0                  -1 -E(12)^7 + E(12)^11]
        """
        # -------------------------------------------------------------------------------
        # create the classical group if not already done
        # -------------------------------------------------------------------------------
        if self._classical_group is None:
            if embedded and self._classical_embedding is None:
                # this is separated to avoid unnecessary (for larger group exhaustive) calculations
                self._create_classical_realization(just_embedded=True)
            else:
                self._create_classical_realization()

        if embedded and self._classical_embedding is not None:
            # ----------------------------------------------------------------------------------------
            # there is a difference between self._classical_group and self._classical_embedding
            # only in the cases where self.strands() divides by 2 (AssionGroupS) resp. 3
            # (AssionGroupU). In this case the embedding is the subgroup of the classical group
            # of one strand more (self.strands() +1) generated by the first self.strands() -1
            # generators
            # ----------------------------------------------------------------------------------------
            return self._classical_embedding
        elif self._classical_group is not None:
            return self._classical_group

        raise ValueError("no classical embedding defined")

    # ----------------------------------------------------------------------------------
    # as_refection_group
    # ----------------------------------------------------------------------------------
    def as_reflection_group(self):
        r"""
        Return an isomorphic image of ``self`` as irreducible complex
        reflection group.

        This is possible only for the finite cubic braid groups of ``cbg_type``
        ``CubicBraidGroup.type.Coxeter``.

        .. NOTE::

            This method uses the sage implementation of reflection group via
            the ``gap3`` ``CHEVIE`` package. These must be installed in order
            to use this method.

        EXAMPLES::

           sage: # optional - gap3
           sage: C3.<c1,c2> = CubicBraidGroup(3)
           sage: R3 = C3.as_reflection_group(); R3
           Irreducible complex reflection group of rank 2 and type ST4
           sage: R3.cartan_matrix()
           [-2*E(3) - E(3)^2           E(3)^2]
           [        -E(3)^2 -2*E(3) - E(3)^2]
           sage: R3.simple_roots()
           Finite family {1: (0, -2*E(3) - E(3)^2), 2: (2*E(3)^2, E(3)^2)}
           sage: R3.simple_coroots()
           Finite family {1: (0, 1), 2: (1/3*E(3) - 1/3*E(3)^2, 1/3*E(3) - 1/3*E(3)^2)}

        Conversion maps::

           sage: # optional - gap3
           sage: r = R3.an_element()
           sage: cr = C3(r); cr
           c1*c2
           sage: mr = r.matrix(); mr
           [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
           [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
           sage: C3Cl = C3.as_classical_group()
           sage: C3Cl(cr)
           [ E(3)^2    -E(4)]
           [-E(12)^7        0]

        The reflection groups can also be viewed as subgroups of unitary groups
        over the universal cyclotomic field. Note that the unitary group
        corresponding to the reflection group is isomorphic but different from
        the classical group due to different hermitian forms for the unitary
        groups they live in::

           sage: # optional - gap3
           sage: C4 = CubicBraidGroup(4)
           sage: R4 = C4.as_reflection_group()
           sage: R4.invariant_form()
           [1 0 0]
           [0 1 0]
           [0 0 1]
           sage: _ == C4.classical_invariant_form()
           False
        """
        # -------------------------------------------------------------------------------
        # the reflection groups are called according to the Shephard-Todd classification:
        #    2 strands -> G(2,1,1)
        #    3 strands -> G4
        #    4 strands -> G25
        #    5 strands -> G32
        # -------------------------------------------------------------------------------
        from sage.combinat.root_system.reflection_group_real import is_chevie_available
        if not is_chevie_available():
            raise ImportError("the GAP3 package 'CHEVIE' is needed to obtain the corresponding reflection groups")

        if self._cbg_type != CubicBraidGroup.type.Coxeter or self.strands() > 5 or self.strands() < 2:
            raise ValueError("no reflection group defined")

        # -------------------------------------------------------------------------------
        # define reflection group associated to self
        # -------------------------------------------------------------------------------
        reflection_group = None

        from sage.combinat.root_system.reflection_group_real import ReflectionGroup

        if self.strands() == 2:
            reflection_group = ReflectionGroup([2, 1, 1])
        elif self.strands() == 3:
            reflection_group = ReflectionGroup(4)
        elif self.strands() == 4:
            reflection_group = ReflectionGroup(25)
        elif self.strands() == 5:
            reflection_group = ReflectionGroup(32)

        hom_to_refl = self.hom(reflection_group.gens())
        reflection_group.register_conversion(hom_to_refl)
        return reflection_group

    # ----------------------------------------------------------------------------------
    # classical invariant form returns the invariant form of the classical realization
    # ----------------------------------------------------------------------------------
    def classical_invariant_form(self):
        r"""
        Return the invariant form of the classical realization of ``self``.

        OUTPUT:

        A square matrix of dimension according to the space the classical
        realization is operating on. In the case of the full cubic braid groups
        and of the Assion groups of ``cbg_type CubicBraidGroup.type.AssionU``
        the matrix is Hermitian. In the case of the Assion groups of
        ``cbg_type CubicBraidGroup.type.AssionS`` it is alternating.
        Note that the invariant form of the full cubic braid group on more
        than 5 strands is degenerated (causing the group to be infinite).

        In the case of Assion groups having projective classical groups,
        the invariant form corresponds to the ambient group of its
        classical embedding.

        EXAMPLES::

            sage: S3 = AssionGroupS(3)
            sage: S3.classical_invariant_form()
            [0 1]
            [2 0]
            sage: S4 = AssionGroupS(4)
            sage: S4.classical_invariant_form()
            [0 0 0 1]
            [0 0 1 0]
            [0 2 0 0]
            [2 0 0 0]
            sage: S5 = AssionGroupS(5)
            sage: S4.classical_invariant_form() == S5.classical_invariant_form()
            True
            sage: U4 = AssionGroupU(4)
            sage: U4.classical_invariant_form()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: C5 = CubicBraidGroup(5)
            sage: C5.classical_invariant_form()
            [-E(12)^7 + E(12)^11                  -1                   0                   0]
            [                -1 -E(12)^7 + E(12)^11                  -1                   0]
            [                 0                  -1 -E(12)^7 + E(12)^11                  -1]
            [                 0                   0                  -1 -E(12)^7 + E(12)^11]
            sage: _.is_singular()
            False
            sage: C6 = CubicBraidGroup(6)
            sage: C6.classical_invariant_form().is_singular()
            True
        """
        # -------------------------------------------------------------------------------
        # create the classical_invariant_form if not already done
        # -------------------------------------------------------------------------------
        if self._classical_invariant_form is None:
            self._create_classical_realization()

        if self._classical_invariant_form is None:
            raise ValueError("no classical invariant form defined")

        return self._classical_invariant_form

    # ----------------------------------------------------------------------------------
    # centralizing element in the classical symplectic resp. unitary group
    # ----------------------------------------------------------------------------------
    def centralizing_element(self, embedded=False):
        r"""
        Return the centralizing element defined by the work of Assion
        (Hilfssatz 1.1.3 and 1.2.3).

        INPUT:

        - ``embedded`` -- boolean (default; ``False``); this boolean only effects
          the cases of Assion groups when they are realized as projective groups.
          More precisely: if ``self`` is of ``cbg_type CubicBraidGroup.type.AssionS``
          (for example) and the number of strands ``n`` is even, than its
          classical group is a subgroup of ``PSp(n,3)`` being centralized
          by the element return for option ``embedded=False``. Otherwise the
          image of this element inside the embedded classical group will be
          returned (see option embedded of :meth:`classical_group`).

        OUTPUT:

        Depending on the optional keyword a permutation as an element
        of ``PSp(n,3)`` (type S) or ``PGU(n,2)`` (type U) for ``n = 0 mod 2``
        (type S) resp. ``n = 0 mod 3`` (type U) is returned. Otherwise, the
        centralizing element is a matrix belonging to ``Sp(n,3)``
        resp. ``GU(n,2)``.

        EXAMPLES::

            sage: U3 = AssionGroupU(3);  U3
            Assion group on 3 strands of type U
            sage: U3Cl = U3.as_classical_group(); U3Cl
            Subgroup generated by [(1,7,6)(3,19,14)(4,15,10)(5,11,18)(12,16,20),
                                   (1,12,13)(2,15,19)(4,9,14)(5,18,8)(6,21,16)]
             of (The projective general unitary group of degree 3 over Finite Field of size 2)
            sage: c = U3.centralizing_element(); c
            (1,16)(2,9)(3,10)(4,19)(6,12)(7,20)(13,21)(14,15)
            sage: c in U3Cl
            True
            sage: P = U3Cl.ambient_group()
            sage: P.centralizer(c) == U3Cl
            True

        Embedded version::

            sage: cm = U3.centralizing_element(embedded=True); cm
            [a + 1 a + 1     1]
            [a + 1     0     a]
            [   1     a     a]
            sage: U4 = AssionGroupU(4)
            sage: U4Cl = U4.as_classical_group()
            sage: cm in U4Cl
            True
            sage: [cm * U4Cl(g) == U4Cl(g) * cm for g in U4.gens()]
            [True, True, False]
        """
        # -------------------------------------------------------------------------------
        # create the centralizing elements if not already done
        # -------------------------------------------------------------------------------
        if self._centralizing_matrix is None:
            self._create_classical_realization()

        if self._centralizing_matrix is None:
            raise ValueError("no centralizing element defined")
        else:
            if embedded or self._centralizing_element is None:
                return self._centralizing_matrix
            else:
                return self._centralizing_element

    # ----------------------------------------------------------------------------------
    # calculating the order by formula
    # ----------------------------------------------------------------------------------
    def order(self):
        r"""
        To avoid long wait-time on calculations the order will be obtained
        using the classical realization.

        OUTPUT: cardinality of the group as integer or infinity

        EXAMPLES::

            sage: S15 = AssionGroupS(15)
            sage: S15.order()
            109777561863482259035023554842176139436811616256000
            sage: C6 = CubicBraidGroup(6)
            sage: C6.order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        n = self.strands()

        if self._cbg_type == CubicBraidGroup.type.Coxeter and n > 5:
            order = infinity
        else:
            order = self.as_classical_group(embedded=True).order()

        return order

    cardinality = order

    def is_finite(self):
        r"""
        Return if ``self`` is a finite group or not.

        EXAMPLES::

            sage: CubicBraidGroup(6).is_finite()
            False
            sage: AssionGroupS(6).is_finite()
            True
        """
        return not (self._cbg_type == CubicBraidGroup.type.Coxeter and self.strands() > 5)

    # ------------------------------------------------------------------
    # creating a CubicBraidGroup as subgroup of self on less strands
    # ------------------------------------------------------------------
    def cubic_braid_subgroup(self, nstrands=None):
        r"""
        Return a cubic braid group as subgroup of ``self`` on the first
        ``nstrands`` strands.

        INPUT:

        - ``nstrands`` -- (default: ``self.strands() - 1``) integer at least 1
          and at most  ``self.strands()`` giving the number of strands of
          the subgroup

        .. WARNING::

            Since ``self`` is inherited from :class:`UniqueRepresentation`, the
            obtained instance is identical to other instances created with the
            same arguments (see example below). The ambient group corresponds
            to the last call of this method.

        EXAMPLES::

            sage: U5 = AssionGroupU(5)
            sage: U3s = U5.cubic_braid_subgroup(3)
            sage: u1, u2 = U3s.gens()
            sage: u1 in U5
            False
            sage: U5(u1) in U5.gens()
            True
            sage: U3s is AssionGroupU(3)
            True
            sage: U3s.ambient() == U5
            True
        """
        if nstrands is None:
            nstrands = self.strands() - 1

        n = self.strands()

        nstrands = Integer(nstrands)

        if nstrands >= n or nstrands <= 0:
            raise ValueError("nstrands must be positive and less than %s" % (self.strands()))

        names = self.variable_names()
        names_red = names[:nstrands - 1]
        subgrp = CubicBraidGroup(names=names_red, cbg_type=self._cbg_type)
        subgrp._ambient = self
        return subgrp
