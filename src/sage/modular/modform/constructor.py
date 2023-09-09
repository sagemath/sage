# -*- coding: utf-8 -*-
"""
Creating spaces of modular forms

EXAMPLES::

    sage: m = ModularForms(Gamma1(4),11)
    sage: m
    Modular Forms space of dimension 6 for Congruence Subgroup Gamma1(4) of weight 11 over Rational Field
    sage: m.basis()
    [
    q - 134*q^5 + O(q^6),
    q^2 + 80*q^5 + O(q^6),
    q^3 + 16*q^5 + O(q^6),
    q^4 - 4*q^5 + O(q^6),
    1 + 4092/50521*q^2 + 472384/50521*q^3 + 4194300/50521*q^4 + O(q^6),
    q + 1024*q^2 + 59048*q^3 + 1048576*q^4 + 9765626*q^5 + O(q^6)
    ]

"""

#*****************************************************************************
#       Copyright (C) 2004-2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import weakref
import re

import sage.modular.arithgroup.all as arithgroup
import sage.modular.dirichlet as dirichlet
from sage.rings.integer import Integer
from sage.rings.rational_field import Q as QQ
from sage.rings.ring import CommutativeRing

from .ambient_eps import ModularFormsAmbient_eps
from .ambient_g0 import ModularFormsAmbient_g0_Q
from .ambient_g1 import ModularFormsAmbient_g1_Q, ModularFormsAmbient_gH_Q
from . import ambient_R
from . import defaults


def canonical_parameters(group, level, weight, base_ring):
    """
    Given a group, level, weight, and base_ring as input by the user,
    return a canonicalized version of them, where level is a Sage
    integer, group really is a group, weight is a Sage integer, and
    base_ring a Sage ring. Note that we can't just get the level from
    the group, because we have the convention that the character for
    Gamma1(N) is None (which makes good sense).

    INPUT:


    -  ``group`` - int, long, Sage integer, group,
       Dirichlet character, or

    -  ``level`` - int, long, Sage integer, or group

    -  ``weight`` - coercible to Sage integer

    -  ``base_ring`` - commutative Sage ring


    OUTPUT:


    -  ``level`` - Sage integer

    -  ``group`` - congruence subgroup

    -  ``weight`` - Sage integer

    -  ``ring`` - commutative Sage ring


    EXAMPLES::

        sage: from sage.modular.modform.constructor import canonical_parameters
        sage: v = canonical_parameters(5, 5, int(7), ZZ); v
        (5, Congruence Subgroup Gamma0(5), 7, Integer Ring)
        sage: type(v[0]), type(v[1]), type(v[2]), type(v[3])
        (<class 'sage.rings.integer.Integer'>,
         <class 'sage.modular.arithgroup.congroup_gamma0.Gamma0_class_with_category'>,
         <class 'sage.rings.integer.Integer'>,
         <class 'sage.rings.integer_ring.IntegerRing_class'>)
        sage: canonical_parameters( 5, 7, 7, ZZ )
        Traceback (most recent call last):
        ...
        ValueError: group and level do not match.
    """
    weight = Integer(weight)
    if weight <= 0:
        raise NotImplementedError("weight must be at least 1")

    if isinstance(group, dirichlet.DirichletCharacter):
        if ( group.level() != Integer(level) ):
            raise ValueError("group.level() and level do not match.")
        group = group.minimize_base_ring()
        level = Integer(level)

    elif arithgroup.is_CongruenceSubgroup(group):
        if ( Integer(level) != group.level() ):
            raise ValueError("group.level() and level do not match.")
        # normalize the case of SL2Z
        if arithgroup.is_SL2Z(group) or \
           arithgroup.is_Gamma1(group) and group.level() == Integer(1):
            group = arithgroup.Gamma0(Integer(1))

    elif group is None:
        pass

    else:
        try:
            m = Integer(group)
        except TypeError:
            raise TypeError("group of unknown type.")
        level = Integer(level)
        if ( m != level ):
            raise ValueError("group and level do not match.")
        group = arithgroup.Gamma0(m)

    if not isinstance(base_ring, CommutativeRing):
        raise TypeError("base_ring (=%s) must be a commutative ring"%base_ring)

    # it is *very* important to include the level as part of the data
    # that defines the key, since Dirichlet characters of different
    # levels can compare equal, but define much different modular
    # forms spaces.
    return level, group, weight, base_ring

_cache = {}

def ModularForms_clear_cache():
    """
    Clear the cache of modular forms.

    EXAMPLES::

        sage: M = ModularForms(37,2)
        sage: sage.modular.modform.constructor._cache == {}
        False

    ::

        sage: sage.modular.modform.constructor.ModularForms_clear_cache()
        sage: sage.modular.modform.constructor._cache
        {}
    """
    global _cache
    _cache = {}


def ModularForms(group=1,
                 weight=2,
                 base_ring=None,
                 eis_only=False,
                 use_cache=True,
                 prec=defaults.DEFAULT_PRECISION):
    r"""
    Create an ambient space of modular forms.

    INPUT:

    - ``group`` - A congruence subgroup or a Dirichlet character eps.

    - ``weight`` - int, the weight, which must be an integer >= 1.

    - ``base_ring`` - the base ring (ignored if group is a Dirichlet character)

    - ``eis_only`` - if True, compute only the Eisenstein part of the space.
      Only permitted (and only useful) in weight 1, where computing dimensions
      of cusp form spaces is expensive.

    Create using the command ModularForms(group, weight, base_ring)
    where group could be either a congruence subgroup or a Dirichlet
    character.

    EXAMPLES: First we create some spaces with trivial character::

        sage: ModularForms(Gamma0(11),2).dimension()
        2
        sage: ModularForms(Gamma0(1),12).dimension()
        2

    If we give an integer N for the congruence subgroup, it defaults to
    `\Gamma_0(N)`::

        sage: ModularForms(1,12).dimension()
        2
        sage: ModularForms(11,4)
        Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(11) of weight 4 over Rational Field

    We create some spaces for `\Gamma_1(N)`.

    ::

        sage: ModularForms(Gamma1(13),2)
        Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
        sage: ModularForms(Gamma1(13),2).dimension()
        13
        sage: [ModularForms(Gamma1(7),k).dimension() for k in [2,3,4,5]]
        [5, 7, 9, 11]
        sage: ModularForms(Gamma1(5),11).dimension()
        12

    We create a space with character::

        sage: e = (DirichletGroup(13).0)^2
        sage: e.order()
        6
        sage: M = ModularForms(e, 2); M
        Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
        sage: f = M.T(2).charpoly('x'); f
        x^3 + (-2*zeta6 - 2)*x^2 - 2*zeta6*x + 14*zeta6 - 7
        sage: f.factor()
        (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)

    We can also create spaces corresponding to the groups `\Gamma_H(N)` intermediate
    between `\Gamma_0(N)` and `\Gamma_1(N)`::

        sage: G = GammaH(30, [11])
        sage: M = ModularForms(G, 2); M
        Modular Forms space of dimension 20 for Congruence Subgroup Gamma_H(30) with H generated by [11] of weight 2 over Rational Field
        sage: M.T(7).charpoly().factor()  # long time (7s on sage.math, 2011)
        (x + 4) * x^2 * (x - 6)^4 * (x + 6)^4 * (x - 8)^7 * (x^2 + 4)

    More examples of spaces with character::

        sage: e = DirichletGroup(5, RationalField()).gen(); e
        Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1

        sage: m = ModularForms(e, 2); m
        Modular Forms space of dimension 2, character [-1] and weight 2 over Rational Field
        sage: m == loads(dumps(m))
        True
        sage: m.T(2).charpoly('x')
        x^2 - 1
        sage: m = ModularForms(e, 6); m.dimension()
        4
        sage: m.T(2).charpoly('x')
        x^4 - 917*x^2 - 42284

    This came up in a subtle bug (:trac:`5923`)::

        sage: ModularForms(gp(1), gap(12))
        Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field

    This came up in another bug (related to :trac:`8630`)::

        sage: chi = DirichletGroup(109, CyclotomicField(3)).0
        sage: ModularForms(chi, 2, base_ring = CyclotomicField(15))
        Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2 over Cyclotomic Field of order 15 and degree 8

    We create some weight 1 spaces. Here modular symbol algorithms do not work.
    In some small examples we can prove using Riemann--Roch that there are no
    cusp forms anyway, so the entire space is Eisenstein::

        sage: M = ModularForms(Gamma1(11), 1); M
        Modular Forms space of dimension 5 for Congruence Subgroup Gamma1(11) of weight 1 over Rational Field
        sage: M.basis()
        [
        1 + 22*q^5 + O(q^6),
        q + 4*q^5 + O(q^6),
        q^2 - 4*q^5 + O(q^6),
        q^3 - 5*q^5 + O(q^6),
        q^4 - 3*q^5 + O(q^6)
        ]
        sage: M.cuspidal_subspace().basis()
        [
        ]
        sage: M == M.eisenstein_subspace()
        True

    When this does not work (which happens as soon as the level is more than
    about 30), we use the Hecke stability algorithm of George Schaeffer::

        sage: M = ModularForms(Gamma1(57), 1); M # long time
        Modular Forms space of dimension 38 for Congruence Subgroup Gamma1(57) of weight 1 over Rational Field
        sage: M.cuspidal_submodule().basis() # long time
        [
        q - q^4 + O(q^6),
        q^3 - q^4 + O(q^6)
        ]

    The Eisenstein subspace in weight 1 can be computed quickly, without
    triggering the expensive computation of the cuspidal part::

        sage: E = EisensteinForms(Gamma1(59), 1); E # indirect doctest
        Eisenstein subspace of dimension 29 of Modular Forms space for Congruence Subgroup Gamma1(59) of weight 1 over Rational Field
        sage: (E.0 + E.2).q_expansion(40)
        1 + q^2 + 196*q^29 - 197*q^30 - q^31 + q^33 + q^34 + q^37 + q^38 - q^39 + O(q^40)

    """
    if isinstance(group, dirichlet.DirichletCharacter):
        if base_ring is None:
            base_ring = group.minimize_base_ring().base_ring()
    if base_ring is None:
        base_ring = QQ

    if isinstance(group, dirichlet.DirichletCharacter) \
           or arithgroup.is_CongruenceSubgroup(group):
        level = group.level()
    else:
        level = group

    eis_only = bool(eis_only)
    key = canonical_parameters(group, level, weight, base_ring) + (eis_only,)

    if eis_only and weight != 1:
        raise ValueError("eis_only parameter only valid in weight 1")

    if use_cache and key in _cache:
        M = _cache[key]()
        if not (M is None):
            M.set_precision(prec)
            return M

    (level, group, weight, base_ring, eis_only) = key

    M = None
    if arithgroup.is_Gamma0(group):
        M = ModularFormsAmbient_g0_Q(group.level(), weight)
        if base_ring != QQ:
            M = ambient_R.ModularFormsAmbient_R(M, base_ring)

    elif arithgroup.is_Gamma1(group):
        M = ModularFormsAmbient_g1_Q(group.level(), weight, eis_only)
        if base_ring != QQ:
            M = ambient_R.ModularFormsAmbient_R(M, base_ring)

    elif arithgroup.is_GammaH(group):
        M = ModularFormsAmbient_gH_Q(group, weight, eis_only)
        if base_ring != QQ:
            M = ambient_R.ModularFormsAmbient_R(M, base_ring)

    elif isinstance(group, dirichlet.DirichletCharacter):
        eps = group
        if eps.base_ring().characteristic() != 0:
            # TODO -- implement this
            # Need to add a lift_to_char_0 function for characters,
            # and need to still remember eps.
            raise NotImplementedError("currently the character must be over a ring of characteristic 0.")
        eps = eps.minimize_base_ring()
        if eps.is_trivial():
            return ModularForms(eps.modulus(), weight, base_ring,
                                use_cache=use_cache,
                                prec=prec)
        M = ModularFormsAmbient_eps(eps, weight, eis_only=eis_only)
        if base_ring != eps.base_ring():
            M = M.base_extend(base_ring) # ambient_R.ModularFormsAmbient_R(M, base_ring)

    if M is None:
        raise NotImplementedError("computation of requested space of modular forms not defined or implemented")

    M.set_precision(prec)
    _cache[key] = weakref.ref(M)
    return M


def CuspForms(group=1,
              weight=2,
              base_ring=None,
              use_cache=True,
              prec=defaults.DEFAULT_PRECISION):
    """
    Create a space of cuspidal modular forms.

    See the documentation for the ModularForms command for a
    description of the input parameters.

    EXAMPLES::

        sage: CuspForms(11,2)
        Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    """
    return ModularForms(group, weight, base_ring,
                        use_cache=use_cache, prec=prec).cuspidal_submodule()


def EisensteinForms(group=1,
              weight=2,
              base_ring=None,
              use_cache=True,
              prec=defaults.DEFAULT_PRECISION):
    """
    Create a space of Eisenstein modular forms.

    See the documentation for the ModularForms command for a
    description of the input parameters.

    EXAMPLES::

        sage: EisensteinForms(11,2)
        Eisenstein subspace of dimension 1 of Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    """
    if weight == 1:
        return ModularForms(group, weight, base_ring,
                        use_cache=use_cache, eis_only=True, prec=prec).eisenstein_submodule()
    else:
        return ModularForms(group, weight, base_ring,
                        use_cache=use_cache, prec=prec).eisenstein_submodule()


def Newforms(group, weight=2, base_ring=None, names=None):
    r"""
    Returns a list of the newforms of the given weight and level (or weight,
    level and character). These are calculated as
    `\operatorname{Gal}(\overline{F} / F)`-orbits, where `F` is the given base
    field.

    INPUT:


    -  ``group`` - the congruence subgroup of the newform, or a Nebentypus
       character

    -  ``weight`` - the weight of the newform (default 2)

    -  ``base_ring`` - the base ring (defaults to `\QQ` for spaces without
       character, or the base ring of the character otherwise)

    -  ``names`` - if the newform has coefficients in a
       number field, a generator name must be specified


    EXAMPLES::

        sage: Newforms(11, 2)
        [q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)]
        sage: Newforms(65, names='a')
        [q - q^2 - 2*q^3 - q^4 - q^5 + O(q^6),
         q + a1*q^2 + (a1 + 1)*q^3 + (-2*a1 - 1)*q^4 + q^5 + O(q^6),
         q + a2*q^2 + (-a2 + 1)*q^3 + q^4 - q^5 + O(q^6)]

    A more complicated example involving both a nontrivial character, and a
    base field that is not minimal for that character::

        sage: K.<i> = QuadraticField(-1)
        sage: chi = DirichletGroup(5, K)[1]
        sage: len(Newforms(chi, 7, names='a'))
        1
        sage: x = polygen(K); L.<c> = K.extension(x^2 - 402*i)
        sage: N = Newforms(chi, 7, base_ring = L); len(N)
        2
        sage: sorted([N[0][2], N[1][2]]) == sorted([1/2*c - 5/2*i - 5/2, -1/2*c - 5/2*i - 5/2])
        True

    TESTS:

    We test that :trac:`8630` is fixed::

        sage: chi = DirichletGroup(109, CyclotomicField(3)).0
        sage: CuspForms(chi, 2, base_ring = CyclotomicField(9))
        Cuspidal subspace of dimension 8 of Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2 over Cyclotomic Field of order 9 and degree 6

    Check that :trac:`15486` is fixed (this used to take over a day)::

        sage: N = Newforms(719, names='a'); len(N)  # long time (3 s)
        3

    """
    return CuspForms(group, weight, base_ring).newforms(names)


def Newform(identifier, group=None, weight=2, base_ring=QQ, names=None):
    """
    INPUT:


    -  ``identifier`` - a canonical label, or the index of
       the specific newform desired

    -  ``group`` - the congruence subgroup of the newform

    -  ``weight`` - the weight of the newform (default 2)

    -  ``base_ring`` - the base ring

    -  ``names`` - if the newform has coefficients in a
       number field, a generator name must be specified


    EXAMPLES::

        sage: Newform('67a', names='a')
        q + 2*q^2 - 2*q^3 + 2*q^4 + 2*q^5 + O(q^6)
        sage: Newform('67b', names='a')
        q + a1*q^2 + (-a1 - 3)*q^3 + (-3*a1 - 3)*q^4 - 3*q^5 + O(q^6)
    """
    if isinstance(group, str) and names is None:
        names = group
    if isinstance(identifier, str):
        group, identifier = parse_label(identifier)
        if weight != 2:
            raise ValueError("Canonical label not implemented for higher weight forms.")
        elif base_ring != QQ:
            raise ValueError("Canonical label not implemented except for over Q.")
    elif group is None:
        raise ValueError("Must specify a group or a label.")
    return Newforms(group, weight, base_ring, names=names)[identifier]


def parse_label(s):
    """
    Given a string s corresponding to a newform label, return the
    corresponding group and index.

    EXAMPLES::

        sage: sage.modular.modform.constructor.parse_label('11a')
        (Congruence Subgroup Gamma0(11), 0)
        sage: sage.modular.modform.constructor.parse_label('11aG1')
        (Congruence Subgroup Gamma1(11), 0)
        sage: sage.modular.modform.constructor.parse_label('11wG1')
        (Congruence Subgroup Gamma1(11), 22)

    GammaH labels should also return the group and index (:trac:`20823`)::

        sage: sage.modular.modform.constructor.parse_label('389cGH[16]')
        (Congruence Subgroup Gamma_H(389) with H generated by [16], 2)
    """
    m = re.match(r'(\d+)([a-z]+)((?:G.*)?)$', s)
    if not m:
        raise ValueError("Invalid label: %s" % s)
    N, order, G = m.groups()
    N = int(N)
    index = 0
    for c in reversed(order):
        index = 26*index + ord(c)-ord('a')
    if G == '' or G == 'G0':
        G = arithgroup.Gamma0(N)
    elif G == 'G1':
        G = arithgroup.Gamma1(N)
    elif G[:2] == 'GH':
        if G[2] != '[' or G[-1] != ']':
            raise ValueError("Invalid congruence subgroup label: %s" % G)
        gens = [int(g.strip()) for g in G[3:-1].split(',')]
        return arithgroup.GammaH(N, gens), index
    else:
        raise ValueError("Invalid congruence subgroup label: %s" % G)
    return G, index
