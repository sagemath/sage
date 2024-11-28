# sage.doctest: needs sage.libs.gap
r"""
Extended Affine Weyl Groups

AUTHORS:

- Daniel Bump (2012): initial version
- Daniel Orr (2012): initial version
- Anne Schilling (2012): initial version
- Mark Shimozono (2012): initial version
- Nicolas M. Thiery (2012): initial version
- Mark Shimozono (2013): twisted affine root systems, multiple realizations, GL_n
"""
# ***************************************************************************
#       Copyright (C) 2012 Daniel Bump <bump at match.stanford.edu>,
#                     2012 Daniel Orr <danorr at live.unc.edu>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
#                     2012 Mark Shimozono <mshimo at math.vt.edu>
#                     2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#                     2013 Mark Shimozono <mshimo at math.vt.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.categories.groups import Groups
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.categories.realizations import Category_realization_of_parent
from sage.misc.bindable_class import BindableClass
from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
from sage.misc.abstract_method import abstract_method
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.groups.group_exp import GroupExp
from sage.groups.group_semidirect_product import GroupSemidirectProduct
from sage.combinat.root_system.root_system import RootSystem
from sage.rings.finite_rings.integer_mod import Mod
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ


def ExtendedAffineWeylGroup(cartan_type, general_linear=None, **print_options):
    r"""
    The extended affine Weyl group.

    INPUT:

    - ``cartan_type`` -- an affine or finite Cartan type (a finite Cartan type is an
      abbreviation for its untwisted affinization)
    - ``general_linear`` -- (default: ``None``) if ``True`` and ``cartan_type`` indicates
      untwisted type A, returns the universal central extension
    - ``print_options`` -- special instructions for printing elements (see below)

    .. RUBRIC:: Mnemonics

    - "P" -- subgroup of translations
    - "Pv" -- subgroup of translations in a dual form
    - "W0" -- classical Weyl group
    - "W" -- affine Weyl group
    - "F" -- fundamental group of length zero elements

    There are currently six realizations: "PW0", "W0P, "WF", "FW", "PvW0", and "W0Pv".

    "PW0" means the semidirect product of "P" with "W0" acting from the right.
    "W0P" is similar but with "W0" acting from the left.
    "WF" is the semidirect product of "W" with "F" acting from the right, etc.

    Recognized arguments for ``print_options`` are:

    - ``print_tuple`` -- boolean (default: ``False``); if ``True``, elements
      are printed `(a,b)`, otherwise as `a * b`
    - ``affine`` -- prefix for simple reflections in the affine Weyl group
    - ``classical`` -- prefix for simple reflections in the classical Weyl group
    - ``translation`` -- prefix for the translation elements
    - ``fundamental`` -- prefix for the elements of the fundamental group

    These options are not mutable.

    The *extended affine Weyl group* was introduced in the following references.

    REFERENCES:

    .. [Iwahori] Iwahori,
       *Generalized Tits system (Bruhat decomposition) on p-adic semisimple groups*.
       1966 Algebraic Groups and Discontinuous
       Subgroups (AMS Proc. Symp. Pure Math.., 1965) pp. 71-83 Amer. Math. Soc.,
       Providence, R.I.

    .. [Bour] Bourbaki, *Lie Groups and Lie Algebras* IV.2

    - [Ka1990]_

    .. RUBRIC:: Notation

    - `R` -- an irreducible affine root system
    - `I` -- set of nodes of the Dynkin diagram of `R`
    - `R_0` -- the classical subsystem of `R`
    - `I_0` -- set of nodes of the Dynkin diagram of `R_0`
    - `E` -- extended affine Weyl group of type `R`
    - `W` -- affine Weyl group of type `R`
    - `W_0` -- finite (classical) Weyl group (of type `R_0`)
    - `M` -- translation lattice for `W`
    - `L` -- translation lattice for `E`
    - `F` -- fundamental subgroup of `E` (the length zero elements)
    - `P` -- finite weight lattice
    - `Q` -- finite root lattice
    - `P^\vee` -- finite coweight lattice
    - `Q^\vee` -- finite coroot lattice

    .. RUBRIC:: Translation lattices

    The styles "PW0" and "W0P" use the following lattices:

    - Untwisted affine: `L = P^\vee`, `M = Q^\vee`
    - Dual of untwisted affine: `L = P`, `M = Q`
    - `BC_n` (`A_{2n}^{(2)}`): `L = M = P`
    - Dual of `BC_n` (`A_{2n}^{(2)\dagger}`): `L = M = P^\vee`

    The styles "PvW0" and "W0Pv" use the following lattices:

    - Untwisted affine: The weight lattice of the dual finite Cartan type.
    - Dual untwisted affine: The same as for "PW0" and "W0P".

    For mixed affine type (`A_{2n}^{(2)}`, aka  `\tilde{BC}_n`, and their affine duals)
    the styles "PvW0" and "W0Pv" are not implemented.

    .. RUBRIC:: Finite and affine Weyl groups `W_0` and `W`

    The finite Weyl group `W_0` is generated by the simple reflections `s_i` for `i \in I_0` where
    `s_i` is the reflection across a suitable hyperplane `H_i` through the origin in the
    real span `V` of the lattice `M`.

    `R` specifies another (affine) hyperplane `H_0`. The affine Weyl group `W` is generated by `W_0`
    and the reflection `S_0` across `H_0`.

    .. RUBRIC:: Extended affine Weyl group `E`

    The complement in `V` of the set `H` of hyperplanes obtained from the `H_i` by the action of
    `W`, has connected components called alcoves. `W` acts freely and transitively on the set
    of alcoves. After the choice of a certain alcove (the fundamental alcove),
    there is an induced bijection from `W` to the set of alcoves under which the identity
    in `W` maps to the fundamental alcove.

    Then `L` is the largest sublattice of `V`, whose translations stabilize the set of alcoves.

    There are isomorphisms

    .. MATH::

        \begin{aligned}
        W &\cong M \rtimes W_0 \cong W_0 \ltimes M \\
        E &\cong L \rtimes W_0 \cong W_0 \ltimes L
        \end{aligned}

    .. RUBRIC:: Fundamental group of affine Dynkin automorphisms

    Since `L` acts on the set of alcoves, the group `F = L/M` may be viewed as a
    subgroup of the symmetries of the fundamental alcove or equivalently the
    symmetries of the affine Dynkin diagram.
    `F` acts on the set of alcoves and hence on `W`. Conjugation by an element of `F`
    acts on `W` by permuting the indices of simple reflections.

    There are isomorphisms

    .. MATH::

        E \cong F \ltimes W \cong W \rtimes F

    An affine Dynkin node is *special* if it is conjugate to the zero node under some
    affine Dynkin automorphism.

    There is a bijection `i` `\mapsto` `\pi_i` from the set of special nodes
    to the group `F`, where `\pi_i` is the unique element of `F` that sends `0` to `i`.
    When `L=P` (resp. `L=P^\vee`) the element `\pi_i` is induced
    (under the isomorphism `F \cong L/M`) by addition of the coset of the
    `i`-th fundamental weight (resp. coweight).

    The length function of the Coxeter group `W` may be extended to `E` by
    `\ell(w \pi) = \ell(w)` where `w \in W` and `\pi\in F`.
    This is the number of hyperplanes in `H` separating the
    fundamental alcove from its image by `w \pi` (or equivalently `w`).

    It is known that if `G` is the compact Lie group of adjoint type with root
    system `R_0` then `F` is isomorphic to the fundamental group of `G`, or
    to the center of its simply-connected covering group. That is why we
    call `F` the *fundamental group*.

    In the future we may want to build an element of the group from an appropriate linear map f
    on some of the root lattice realizations for this Cartan type: W.from_endomorphism(f).

    EXAMPLES::

        sage: E = ExtendedAffineWeylGroup(["A",2,1]); E
        Extended affine Weyl group of type ['A', 2, 1]
        sage: type(E)
        <class 'sage.combinat.root_system.extended_affine_weyl_group.ExtendedAffineWeylGroup_Class_with_category'>

        sage: PW0 = E.PW0(); PW0
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
         Multiplicative form of Coweight lattice of the Root system of type ['A', 2]
         acted upon by Weyl Group of type ['A', 2]
         (as a matrix group acting on the coweight lattice)

        sage: W0P = E.W0P(); W0P
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
         Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
         acting on Multiplicative form of Coweight lattice of the Root system of type ['A', 2]

        sage: PvW0 = E.PvW0(); PvW0
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
         Multiplicative form of Weight lattice of the Root system of type ['A', 2]
         acted upon by Weyl Group of type ['A', 2]
         (as a matrix group acting on the weight lattice)

        sage: W0Pv = E.W0Pv(); W0Pv
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
         Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
         acting on Multiplicative form of Weight lattice of the Root system of type ['A', 2]

        sage: WF = E.WF(); WF
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
         Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
         acted upon by Fundamental group of type ['A', 2, 1]

        sage: FW = E.FW(); FW
        Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
         Fundamental group of type ['A', 2, 1] acting on Weyl Group of type ['A', 2, 1]
         (as a matrix group acting on the root lattice)

    When the realizations are constructed from each other as above, there are built-in coercions between them. ::

        sage: F = E.fundamental_group()
        sage: x = WF.from_reduced_word([0,1,2]) * WF(F(2)); x
        S0*S1*S2 * pi[2]
        sage: FW(x)
        pi[2] * S1*S2*S0
        sage: W0P(x)
        s1*s2*s1 * t[-2*Lambdacheck[1] - Lambdacheck[2]]
        sage: PW0(x)
        t[Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2*s1
        sage: PvW0(x)
        t[Lambda[1] + 2*Lambda[2]] * s1*s2*s1

    The translation lattice and its distinguished basis are obtained from ``E``::

        sage: L = E.lattice(); L
        Coweight lattice of the Root system of type ['A', 2]
        sage: b = E.lattice_basis(); b
        Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}

    Translation lattice elements can be coerced into any realization::

        sage: PW0(b[1]-b[2])
        t[Lambdacheck[1] - Lambdacheck[2]]
        sage: FW(b[1]-b[2])
        pi[2] * S0*S1

    The dual form of the translation lattice and its basis are similarly obtained::

        sage: Lv = E.dual_lattice(); Lv
        Weight lattice of the Root system of type ['A', 2]
        sage: bv = E.dual_lattice_basis(); bv
        Finite family {1: Lambda[1], 2: Lambda[2]}
        sage: FW(bv[1]-bv[2])
        pi[2] * S0*S1

    The abstract fundamental group is accessed from ``E``::

        sage: F = E.fundamental_group(); F
        Fundamental group of type ['A', 2, 1]

    Its elements are indexed by the set of special nodes of the affine Dynkin diagram::

        sage: E.cartan_type().special_nodes()
        (0, 1, 2)
        sage: F.special_nodes()
        (0, 1, 2)
        sage: [F(i) for i in F.special_nodes()]
        [pi[0], pi[1], pi[2]]

    There is a coercion from the fundamental group into each realization::

        sage: F(2)
        pi[2]
        sage: WF(F(2))
        pi[2]
        sage: W0P(F(2))
        s2*s1 * t[-Lambdacheck[1]]
        sage: W0Pv(F(2))
        s2*s1 * t[-Lambda[1]]

    Using ``E`` one may access the classical and affine Weyl groups and their morphisms
    into each realization::

        sage: W0 = E.classical_weyl(); W0
        Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
        sage: v = W0.from_reduced_word([1,2,1]); v
        s1*s2*s1
        sage: PW0(v)
        s1*s2*s1
        sage: WF(v)
        S1*S2*S1
        sage: W = E.affine_weyl(); W
        Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
        sage: w = W.from_reduced_word([2,1,0]); w
        S2*S1*S0
        sage: WF(w)
        S2*S1*S0
        sage: PW0(w)
        t[Lambdacheck[1] - 2*Lambdacheck[2]] * s1

    Note that for untwisted affine type, the dual form of the classical
    Weyl group is isomorphic to the usual one, but acts on a different
    lattice and is therefore different to sage::

        sage: W0v = E.dual_classical_weyl(); W0v
        Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
        sage: v = W0v.from_reduced_word([1,2])
        sage: x = PvW0(v); x
        s1*s2
        sage: y = PW0(v); y
        s1*s2
        sage: x.parent() == y.parent()
        False

    However, because there is a coercion from ``PvW0`` to ``PW0``,
    the elements ``x`` and ``y`` compare as equal::

        sage: x == y
        True

    An element can be created directly from a reduced word::

        sage: PW0.from_reduced_word([2,1,0])
        t[Lambdacheck[1] - 2*Lambdacheck[2]] * s1

    Here is a demonstration of the printing options::

        sage: E = ExtendedAffineWeylGroup(["A",2,1], affine='sx', classical='Sx',
        ....:                             translation='x', fundamental='pix')
        sage: PW0 = E.PW0()
        sage: y = PW0(E.lattice_basis()[1]); y
        x[Lambdacheck[1]]
        sage: FW = E.FW()
        sage: FW(y)
        pix[1] * sx2*sx1
        sage: PW0.an_element()
        x[2*Lambdacheck[1] + 2*Lambdacheck[2]] * Sx1*Sx2

    .. TODO::

        - Implement a "slow" action of `E` on any affine root or weight lattice realization.
        - Implement the level `m` actions of `E` and `W` on the lattices of finite type.
        - Implement the relevant methods from the usual affine Weyl group
        - Implementation by matrices: style "M".
        - Use case: implement the Hecke algebra on top of this

    The semidirect product construction in sage currently only
    admits multiplicative groups. Therefore for the styles involving "P" and "Pv", one must
    convert the additive group of translations `L` into a multiplicative group by
    applying the :class:`sage.groups.group_exp.GroupExp` functor.

    .. RUBRIC:: The general linear case

    The general linear group is not semisimple. Sage can build its extended
    affine Weyl group::

        sage: E = ExtendedAffineWeylGroup(['A',2,1], general_linear=True); E
        Extended affine Weyl group of GL(3)

    If the Cartan type is ``['A', n-1, 1]`` and the parameter ``general_linear`` is not
    True, the extended affine Weyl group that is built will be for `SL_n`, not
    `GL_n`. But if ``general_linear`` is True, let `W_a` and `W_e` be the affine and
    extended affine Weyl groups. We make the following nonstandard definition: the
    extended affine Weyl group `W_e(GL_n)` is defined by

    .. MATH::

        W_e(GL_n) = P(GL_n) \rtimes W

    where `W` is the finite Weyl group (the symmetric group `S_n`) and `P(GL_n)` is the weight lattice
    of `GL_n`, which is usually identified with the lattice `\ZZ^n` of `n`-tuples of integers::

        sage: PW0 = E.PW0(); PW0
        Extended affine Weyl group of GL(3) realized by Semidirect product of
         Multiplicative form of Ambient space of the Root system of type ['A', 2] acted upon
          by Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        sage: PW0.an_element()
        t[(2, 2, 3)] * s1*s2

    There is an isomorphism

    .. MATH::

        W_e(GL_n) = \ZZ \ltimes W_a

    where the group of integers `\ZZ` (with generator `\pi`) acts on `W_a` by

    .. MATH::

        \pi\, s_i\, \pi^{-1} = s_{i+1}

    and the indices of the simple reflections are taken modulo `n`::

        sage: FW = E.FW(); FW
        Extended affine Weyl group of GL(3) realized by
         Semidirect product of Fundamental group of GL(3) acting on
          Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
        sage: FW.an_element()
        pi[5] * S0*S1*S2

    We regard `\ZZ` as the fundamental group of affine type `GL_n`::

        sage: F = E.fundamental_group(); F
        Fundamental group of GL(3)
        sage: F.special_nodes()
        Integer Ring

        sage: x = FW.from_fundamental(F(10)); x
        pi[10]
        sage: x*x
        pi[20]
        sage: E.PvW0()(x*x)
        t[(7, 7, 6)] * s2*s1
    """
    cartan_type = CartanType(cartan_type)
    if cartan_type.is_reducible():
        raise ValueError("Extended affine Weyl groups are only implemented for irreducible affine Cartan types")
    if cartan_type.is_finite(): # a finite Cartan type is an abbreviation for its untwisted affinization
        cartan_type = cartan_type.affine()
    elif not cartan_type.is_affine():
        raise ValueError("Cartan type must be finite or affine")

    return ExtendedAffineWeylGroup_Class(cartan_type, general_linear, **print_options)


class ExtendedAffineWeylGroup_Class(UniqueRepresentation, Parent):
    r"""
    The parent-with-realization class of an extended affine Weyl group.
    """

    def __init__(self, cartan_type, general_linear, **print_options):
        r"""

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(["D",3,2])
            sage: E in Groups().Infinite()
            True
            sage: TestSuite(E).run()
        """
        if not cartan_type.is_affine():
            raise ValueError("%s is not affine" % cartan_type)

        self._cartan_type = cartan_type

        self._prefixt = "t"
        self._prefixf = "pi"
        self._prefixcl = None
        self._prefixaf = None
        self._print_tuple = False

        if general_linear is True:
            self._general_linear = True
            self._n = self._cartan_type.n + 1
        else:
            self._general_linear = False
        for option in print_options:
            if option == 'translation':
                self._prefixt = print_options['translation']
            elif option == 'fundamental':
                self._prefixf = print_options['fundamental']
            elif option == 'print_tuple':
                self._print_tuple = print_options['print_tuple']
            elif option == 'affine':
                self._prefixaf = print_options['affine']
            elif option == 'classical':
                self._prefixcl = print_options['classical']
            else:
                raise ValueError("Print option %s is unrecognized" % option)

        if self._prefixaf:
            if not self._prefixcl:
                if self._prefixaf.islower():
                    self._prefixcl = self._prefixaf.upper()
                else:
                    self._prefixcl = self._prefixaf.lower()
        elif self._prefixcl:
            if self._prefixcl.islower():
                self._prefixaf = self._prefixcl.upper()
            else:
                self._prefixaf = self._prefixcl.lower()
        else:
            self._prefixaf = "S"
            self._prefixcl = "s"

        self._ct0 = cartan_type.classical()
        self._R0 = self._ct0.root_system()
        self._I0 = self._ct0.index_set()
        self._ct0v = self._ct0.dual()
        self._R0v = self._ct0v.root_system()
        self._a0check = self._cartan_type.acheck()[self._cartan_type.special_node()]

        if self._cartan_type.is_untwisted_affine():
            self._type = 'untwisted'
        elif self._cartan_type.dual().is_untwisted_affine():
            self._type = 'dual_untwisted'
        elif self._a0check == 1:
            # if there are three root lengths with the special affine node extra short
            self._type = 'special_extra_short'
        else:
            # if there are three root lengths with the special affine node extra long
            self._type = 'special_extra_long'
        # this boolean is used to decide which translation lattice to use
        self._untwisted = (self._type in ('untwisted', 'special_extra_long'))

        # fundamental group
        self._fundamental_group = FundamentalGroupOfExtendedAffineWeylGroup(cartan_type, prefix=self._prefixf, general_linear=self._general_linear)

        # lattice data
        if self._untwisted:
            if self._general_linear:
                self._lattice = self._R0.ambient_space()
                self._simpleR0 = self._lattice.simple_roots()
            else:
                self._lattice = self._R0.coweight_lattice()
                self._basis_name = 'Lambdacheck'
                self._simpleR0 = self._R0.root_lattice().simple_roots()
            self._basis = self._lattice.fundamental_weights()
            if self._type == 'special_extra_long':
                self._special_root = self._R0.coroot_lattice().highest_root()
                # get the node adjacent to the special affine node
                # the [0] is just taking the first and only list element among the neighbors of the distinguished node
                node_adjacent_to_special = self._cartan_type.dynkin_diagram().neighbors(self._cartan_type.special_node())[0]
                self._special_translation = self._lattice.fundamental_weight(node_adjacent_to_special)
            else:
                # untwisted affine case
                self._special_root = self._R0.root_lattice().highest_root().associated_coroot()
                self._special_translation = self._special_root
            self._special_translation_covector = self._special_root.associated_coroot()
            # in the "Pv" realization for the untwisted case, the weight lattice of dual type is used for translations
            if self._general_linear:
                self._dual_lattice = self._lattice
            else:
                self._dual_lattice = self._R0v.weight_lattice()
                self._dual_basis_name = 'Lambda'
            self._dual_basis = self._dual_lattice.fundamental_weights()
        else:
            self._lattice = self._R0.weight_lattice()
            self._basis = self._lattice.fundamental_weights()
            self._basis_name = 'Lambda'
            self._simpleR0 = self._R0.coroot_lattice().simple_roots()
            if self._type == 'special_extra_short':
                self._special_root = self._R0.root_lattice().highest_root()
                node_adjacent_to_special = self._cartan_type.dynkin_diagram().neighbors(self._cartan_type.special_node())[0]
                self._special_translation = self._lattice.fundamental_weight(node_adjacent_to_special)
                self._special_translation_covector = 2*self._special_root.associated_coroot()
            else:
                # dual untwisted case
                self._special_root = self._R0.coroot_lattice().highest_root().associated_coroot()
                self._special_translation = self._special_root
                self._special_translation_covector = self._special_root.associated_coroot()

            self._dual_lattice = self._lattice
            self._dual_basis = self._basis
            self._dual_basis_name = 'Lambda'

        # classical and affine Weyl groups
        self._W0 = WeylGroup(self._lattice, prefix=self._prefixcl)
        self._W = WeylGroup(self._cartan_type.root_system().root_lattice(), prefix=self._prefixaf)
        self._special_reflection = self._W0.from_reduced_word(self._special_root.associated_reflection())

        if self._general_linear:
            self._special_root = self._special_root.to_ambient()
            self._special_translation = self._special_root
            self._special_translation_covector = self._special_root

        # "Pv" version of classical Weyl group; use same prefix as for W0
        self._W0v = WeylGroup(self._dual_lattice, prefix=self._prefixcl)

        # wrap the lattice into a multiplicative group for internal use in the semidirect product
        self._exp_lattice = GroupExp()(self._lattice)
        self._exp_dual_lattice = GroupExp()(self._dual_lattice)

        self._extended = True

        Parent.__init__(self, category=Groups().WithRealizations().Infinite())

        # create the realizations (they are cached)
        PW0 = self.PW0()
        W0P = self.W0P()
        WF = self.WF()
        FW = self.FW()
        PvW0 = self.PvW0()
        W0Pv = self.W0Pv()

        # coercions between realizations

        W0P_to_PW0 = SetMorphism(Hom(W0P, PW0, Groups()), lambda x: PW0(x.to_opposite()))
        W0P_to_PW0.register_as_coercion()

        PW0_to_W0P = SetMorphism(Hom(PW0, W0P, Groups()), lambda x: W0P(x.to_opposite()))
        PW0_to_W0P.register_as_coercion()

        FW_to_WF = SetMorphism(Hom(FW, WF, Groups()), lambda x: WF(x.to_opposite()))
        FW_to_WF.register_as_coercion()

        WF_to_FW = SetMorphism(Hom(WF, FW, Groups()), lambda x: FW(x.to_opposite()))
        WF_to_FW.register_as_coercion()

        PW0_to_WF = SetMorphism(Hom(PW0, WF, Groups()), self.PW0_to_WF_func)
        PW0_to_WF.register_as_coercion()
        WF_to_PW0 = SetMorphism(Hom(WF, PW0, Groups()), self.WF_to_PW0_func)
        WF_to_PW0.register_as_coercion()

        PvW0_to_W0Pv = SetMorphism(Hom(PvW0, W0Pv, Groups()), lambda x: W0Pv(x.to_opposite()))
        PvW0_to_W0Pv.register_as_coercion()
        W0Pv_to_PvW0 = SetMorphism(Hom(W0Pv, PvW0, Groups()), lambda x: PvW0(x.to_opposite()))
        W0Pv_to_PvW0.register_as_coercion()

        if self._general_linear:
            PW0_to_PvW0 = SetMorphism(Hom(PW0, PvW0, Groups()), lambda x: PvW0((x.cartesian_projection(0),x.cartesian_projection(1))))
            PvW0_to_PW0 = SetMorphism(Hom(PvW0, PW0, Groups()), lambda x: PW0((x.cartesian_projection(0),x.cartesian_projection(1))))
            W0P_to_W0Pv = SetMorphism(Hom(W0P, W0Pv, Groups()), lambda x: W0Pv((x.cartesian_projection(0),x.cartesian_projection(1))))
            W0Pv_to_W0P = SetMorphism(Hom(W0Pv, W0P, Groups()), lambda x: W0P((x.cartesian_projection(0),x.cartesian_projection(1))))
        elif self._untwisted:
            PW0_to_PvW0 = SetMorphism(Hom(PW0, PvW0, Groups()), lambda x: PvW0((self.exp_dual_lattice()(x.cartesian_projection(0).value.to_dual_type_cospace()),self.dual_classical_weyl().from_reduced_word(x.cartesian_projection(1).reduced_word()))))
            PvW0_to_PW0 = SetMorphism(Hom(PvW0, PW0, Groups()), lambda x: PW0((self.exp_lattice()(x.cartesian_projection(0).value.to_dual_type_cospace()),self.classical_weyl().from_reduced_word(x.cartesian_projection(1).reduced_word()))))
            W0P_to_W0Pv = SetMorphism(Hom(W0P, W0Pv, Groups()), lambda x: W0Pv((self.dual_classical_weyl().from_reduced_word(x.cartesian_projection(0).reduced_word()),self.exp_dual_lattice()(x.cartesian_projection(1).value.to_dual_type_cospace()))))
            W0Pv_to_W0P = SetMorphism(Hom(W0Pv, W0P, Groups()), lambda x: W0P((self.classical_weyl().from_reduced_word(x.cartesian_projection(0).reduced_word()),self.exp_lattice()(x.cartesian_projection(1).value.to_dual_type_cospace()))))
        else:
            PW0_to_PvW0 = SetMorphism(Hom(PW0, PvW0, Groups()), lambda x: PvW0((x.cartesian_projection(0),self.dual_classical_weyl().from_reduced_word(x.cartesian_projection(1).reduced_word()))))
            PvW0_to_PW0 = SetMorphism(Hom(PvW0, PW0, Groups()), lambda x: PW0((x.cartesian_projection(0),self.classical_weyl().from_reduced_word(x.cartesian_projection(1).reduced_word()))))
            W0P_to_W0Pv = SetMorphism(Hom(W0P, W0Pv, Groups()), lambda x: W0Pv((self.dual_classical_weyl().from_reduced_word(x.cartesian_projection(0).reduced_word()),x.cartesian_projection(1))))
            W0Pv_to_W0P = SetMorphism(Hom(W0Pv, W0P, Groups()), lambda x: W0P((self.classical_weyl().from_reduced_word(x.cartesian_projection(0).reduced_word()),x.cartesian_projection(1))))

        PW0_to_PvW0.register_as_coercion()
        PvW0_to_PW0.register_as_coercion()
        W0P_to_W0Pv.register_as_coercion()
        W0Pv_to_W0P.register_as_coercion()

        # coercions of the translation lattice into the appropriate realizations
        P_to_PW0 = SetMorphism(Hom(self.lattice(), PW0, Sets()), PW0.from_translation)
        P_to_PW0.register_as_coercion()
        P_to_W0P = SetMorphism(Hom(self.lattice(), W0P, Sets()), W0P.from_translation)
        P_to_W0P.register_as_coercion()
        Pv_to_PvW0 = SetMorphism(Hom(self.dual_lattice(), PvW0, Sets()), PvW0.from_dual_translation)
        Pv_to_PvW0.register_as_coercion()
        Pv_to_W0Pv = SetMorphism(Hom(self.dual_lattice(), W0Pv, Sets()), W0Pv.from_dual_translation)
        Pv_to_W0Pv.register_as_coercion()

        # coercions of the classical Weyl group into the appropriate realizations

        W0_to_PW0 = SetMorphism(Hom(self.classical_weyl(), PW0, Groups()), PW0.from_classical_weyl)
        W0_to_PW0.register_as_coercion()
        W0_to_W0P = SetMorphism(Hom(self.classical_weyl(), W0P, Groups()), W0P.from_classical_weyl)
        W0_to_W0P.register_as_coercion()
        W0v_to_PvW0 = SetMorphism(Hom(self.dual_classical_weyl(), PvW0, Groups()), PvW0.from_dual_classical_weyl)
        W0v_to_PvW0.register_as_coercion()
        W0v_to_W0Pv = SetMorphism(Hom(self.dual_classical_weyl(), W0Pv, Groups()), W0Pv.from_dual_classical_weyl)
        W0v_to_W0Pv.register_as_coercion()

        # coercions of the fundamental group into the appropriate realizations

        F_to_WF = SetMorphism(Hom(self.fundamental_group(), WF, Groups()), WF.from_fundamental)
        F_to_WF.register_as_coercion()
        F_to_FW = SetMorphism(Hom(self.fundamental_group(), FW, Groups()), FW.from_fundamental)
        F_to_FW.register_as_coercion()

        # coercions of the affine Weyl group into the appropriate realizations

        W_to_WF = SetMorphism(Hom(self.affine_weyl(), WF, Groups()), WF.from_affine_weyl)
        W_to_WF.register_as_coercion()

        W_to_FW = SetMorphism(Hom(self.affine_weyl(), FW, Groups()), FW.from_affine_weyl)
        W_to_FW.register_as_coercion()

    def PW0(self):
        r"""
        Realizes ``self`` in "PW0"-style.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).PW0()
            Extended affine Weyl group of type ['A', 2, 1] realized by
             Semidirect product of Multiplicative form of
              Coweight lattice of the Root system of type ['A', 2] acted upon by
               Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
        """
        return self.ExtendedAffineWeylGroupPW0()

    def W0P(self):
        r"""
        Realizes ``self`` in "W0P"-style.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).W0P()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
             acting on Multiplicative form of Coweight lattice of the Root system of type ['A', 2]
        """
        return self.ExtendedAffineWeylGroupW0P()

    def WF(self):
        r"""
        Realizes ``self`` in "WF"-style.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).WF()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
             acted upon by Fundamental group of type ['A', 2, 1]
        """
        return self.ExtendedAffineWeylGroupWF()

    def FW(self):
        r"""
        Realizes ``self`` in "FW"-style.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).FW()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Fundamental group of type ['A', 2, 1] acting on
              Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
        """
        return self.ExtendedAffineWeylGroupFW()

    def PvW0(self):
        r"""
        Realizes ``self`` in "PvW0"-style.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).PvW0()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Multiplicative form of Weight lattice of the Root system of type ['A', 2]
             acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
        """
        return self.ExtendedAffineWeylGroupPvW0()

    def W0Pv(self):
        r"""
        Realizes ``self`` in "W0Pv"-style.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).W0Pv()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
              acting on Multiplicative form of Weight lattice of the Root system of type ['A', 2]
        """
        return self.ExtendedAffineWeylGroupW0Pv()

    def cartan_type(self):
        r"""
        The Cartan type of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(["D",3,2]).cartan_type()
            ['C', 2, 1]^*
        """
        return self._cartan_type

    def _repr_(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1])
            Extended affine Weyl group of type ['A', 2, 1]
        """
        if self._general_linear:
            return "Extended affine Weyl group of GL(%s)" % self._n
        return "Extended affine Weyl group of type %s" % self.cartan_type()

    def fundamental_group(self):
        r"""
        Return the abstract fundamental group.

        EXAMPLES::

            sage: F = ExtendedAffineWeylGroup(['D',5,1]).fundamental_group(); F
            Fundamental group of type ['D', 5, 1]
            sage: [a for a in F]
            [pi[0], pi[1], pi[4], pi[5]]
        """
        return self._fundamental_group

    def lattice(self):
        r"""
        Return the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).lattice()
            Coweight lattice of the Root system of type ['A', 2]
            sage: ExtendedAffineWeylGroup(['A',5,2]).lattice()
            Weight lattice of the Root system of type ['C', 3]
            sage: ExtendedAffineWeylGroup(['A',4,2]).lattice()
            Weight lattice of the Root system of type ['C', 2]
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).lattice()
            Coweight lattice of the Root system of type ['B', 2]
            sage: ExtendedAffineWeylGroup(CartanType(['A',2,1]),
            ....:                         general_linear=True).lattice()
            Ambient space of the Root system of type ['A', 2]
        """

        return self._lattice

    def exp_lattice(self):
        r"""
        Return the multiplicative version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).exp_lattice()
            Multiplicative form of Coweight lattice of the Root system of type ['A', 2]
        """
        return self._exp_lattice

    def lattice_basis(self):
        r"""
        Return the distinguished basis of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).lattice_basis()
            Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}
            sage: ExtendedAffineWeylGroup(['A',5,2]).lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2], 3: Lambda[3]}
            sage: ExtendedAffineWeylGroup(['A',4,2]).lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2]}
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).lattice_basis()
            Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}
        """
        return self._basis

    def dual_lattice(self):
        r"""
        Return the dual version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).dual_lattice()
            Weight lattice of the Root system of type ['A', 2]
            sage: ExtendedAffineWeylGroup(['A',5,2]).dual_lattice()
            Weight lattice of the Root system of type ['C', 3]
        """
        return self._dual_lattice

    def exp_dual_lattice(self):
        r"""
        Return the multiplicative version of the dual version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).exp_dual_lattice()
            Multiplicative form of Weight lattice of the Root system of type ['A', 2]
        """
        return self._exp_dual_lattice

    def dual_lattice_basis(self):
        r"""
        Return the distinguished basis of the dual version of the translation lattice for ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).dual_lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2]}
            sage: ExtendedAffineWeylGroup(['A',5,2]).dual_lattice_basis()
            Finite family {1: Lambda[1], 2: Lambda[2], 3: Lambda[3]}
        """
        return self._dual_basis

    def classical_weyl(self):
        r"""
        Return the classical Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
            sage: ExtendedAffineWeylGroup(['A',5,2]).classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the weight lattice)
            sage: ExtendedAffineWeylGroup(['A',4,2]).classical_weyl()
            Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice)
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).classical_weyl()
            Weyl Group of type ['C', 2] (as a matrix group acting on the coweight lattice)
        """
        return self._W0

    def dual_classical_weyl(self):
        r"""
        Return the dual version of the classical Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).dual_classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
            sage: ExtendedAffineWeylGroup(['A',5,2]).dual_classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the weight lattice)
        """
        return self._W0v

    def affine_weyl(self):
        r"""
        Return the affine Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).affine_weyl()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
            sage: ExtendedAffineWeylGroup(['A',5,2]).affine_weyl()
            Weyl Group of type ['B', 3, 1]^* (as a matrix group acting on the root lattice)
            sage: ExtendedAffineWeylGroup(['A',4,2]).affine_weyl()
            Weyl Group of type ['BC', 2, 2] (as a matrix group acting on the root lattice)
            sage: ExtendedAffineWeylGroup(CartanType(['A',4,2]).dual()).affine_weyl()
            Weyl Group of type ['BC', 2, 2]^* (as a matrix group acting on the root lattice)
        """
        return self._W

    def classical_weyl_to_affine(self, w):
        r"""
        The image of `w` under the homomorphism from the classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(['A',2,1])
            sage: W0 = E.classical_weyl()
            sage: w = W0.from_reduced_word([1,2]); w
            s1*s2
            sage: v = E.classical_weyl_to_affine(w); v
            S1*S2
        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def dual_classical_weyl_to_affine(self, w):
        r"""
        The image of `w` under the homomorphism from the dual version of the classical
        Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(['A',2,1])
            sage: W0v = E.dual_classical_weyl()
            sage: w = W0v.from_reduced_word([1,2]); w
            s1*s2
            sage: v = E.dual_classical_weyl_to_affine(w); v
            S1*S2
        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def a_realization(self):
        r"""
        Return the default realization of an extended affine Weyl group.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).a_realization()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Multiplicative form of Coweight lattice of the Root system of type ['A', 2]
             acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
        """
        return self.PW0()

    def group_generators(self):
        r"""
        Return a set of generators for the default realization of ``self``.

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).group_generators()
            (t[Lambdacheck[1]], t[Lambdacheck[2]], s1, s2)
        """
        return self.a_realization().group_generators()

    @cached_method
    def PW0_to_WF_func(self, x):
        r"""
        Implement coercion from style "PW0" to "WF".

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(["A", 2, 1])
            sage: x = E.PW0().an_element(); x
            t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2
            sage: E.PW0_to_WF_func(x)
            S0*S1*S2*S0*S1*S0

        .. WARNING::

            This function cannot use coercion, because it is used to define the coercion maps.
        """
        i = x.first_descent(side='left')
        if i is None:
            t = x.to_translation_left()
            # t must be zero or a special fundamental basis element
            if self._general_linear:
                ispecial = ZZ.sum([t[j] for j in t.support()])
            elif t == self.lattice().zero():
                ispecial = 0
            else:
                supp = t.support()
                assert len(supp) == 1
                ispecial = supp[0]
            return self.WF().from_fundamental(self.fundamental_group()(ispecial))
        return self.PW0_to_WF_func(x.apply_simple_reflection(i, side='left')).apply_simple_reflection(i, side='left')

    @cached_method
    def WF_to_PW0_func(self, x):
        r"""
        Coercion from style "WF" to "PW0".

        EXAMPLES::

            sage: E = ExtendedAffineWeylGroup(["A", 2, 1])
            sage: x = E.WF().an_element(); x
            S0*S1*S2 * pi[2]
            sage: E.WF_to_PW0_func(x)
            t[Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2*s1

        .. WARNING::

            Since this is used to define some coercion maps it cannot itself use coercion.
        """
        w = x.to_affine_weyl_left()
        f = x.to_fundamental_group()
        i = w.first_descent(side='left')
        if i is not None:
            return self.WF_to_PW0_func(x.apply_simple_reflection(i, side='left')).apply_simple_reflection(i, side='left')
        # the element is in the fundamental group
        PW0 = self.PW0()
        ispecial = f.value()
        W = self.classical_weyl()
        if self._general_linear:
            r = ZZ(Mod(ispecial, self._n))
            weight = self.lattice().from_vector(vector([ZZ((ispecial-r)/self._n)]*self._n))
            if r != ZZ(0):
                weight = weight + self.lattice_basis()[r]
                wo = W.from_reduced_word(self.fundamental_group().reduced_word(r))
            else:
                wo = W.one()
        elif ispecial == 0:
            weight = self.lattice().zero()
            wo = W.one()
        else:
            weight = self.lattice_basis()[ispecial]
            wo = W.from_reduced_word(self.fundamental_group().reduced_word(ispecial))
        return PW0((weight,wo))

    class Realizations(Category_realization_of_parent):
        r"""
        The category of the realizations of an extended affine Weyl group
        """

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: R = ExtendedAffineWeylGroup(['A',2,1]).Realizations(); R
                Category of realizations of Extended affine Weyl group of type ['A', 2, 1]
                sage: R.super_categories()
                [Category of associative inverse realizations of unital magmas]
            """
            return [Groups().Realizations()]

        class ParentMethods:

            @cached_method
            def from_fundamental(self, x):
                r"""
                Return the image of `x` under the homomorphism from the fundamental group into
                ``self``.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1])
                    sage: PW0=E.PW0()
                    sage: F = E.fundamental_group()
                    sage: Is = F.special_nodes()
                    sage: [(i, PW0.from_fundamental(F(i))) for i in Is]
                    [(0, 1),
                     (1, t[Lambdacheck[1]] * s1*s2*s3),
                     (2, t[Lambdacheck[2]] * s2*s3*s1*s2),
                     (3, t[Lambdacheck[3]] * s3*s2*s1)]
                    sage: [(i, E.W0P().from_fundamental((F(i)))) for i in Is]
                    [(0, 1),
                     (1, s1*s2*s3 * t[-Lambdacheck[3]]),
                     (2, s2*s3*s1*s2 * t[-Lambdacheck[2]]),
                     (3, s3*s2*s1 * t[-Lambdacheck[1]])]
                    sage: [(i, E.WF().from_fundamental(F(i))) for i in Is]
                    [(0, 1), (1, pi[1]), (2, pi[2]), (3, pi[3])]

                .. WARNING::

                    This method must be implemented by the "WF" and "FW" realizations.
                """
                WF = self.realization_of().WF()
                return self(WF.from_fundamental(x))

            def from_translation(self, la):
                r"""
                Return the element of translation by ``la`` in ``self``.

                INPUT:

                - ``self`` -- a realization of the extended affine Weyl group
                - ``la`` -- an element of the translation lattice

                In the notation of the documentation for :meth:`ExtendedAffineWeylGroup`,
                ``la`` must be an element of "P".

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1]); PW0 = E.PW0()
                    sage: b = E.lattice_basis(); b
                    Finite family {1: Lambdacheck[1], 2: Lambdacheck[2]}
                    sage: x = PW0.from_translation(2*b[1] - b[2]); x
                    t[2*Lambdacheck[1] - Lambdacheck[2]]
                    sage: FW = E.FW()
                    sage: y = FW.from_translation(2*b[1] - b[2]); y
                    S0*S2*S0*S1
                    sage: FW(x) == y
                    True

                Since the implementation as a semidirect product requires
                wrapping the lattice group to make it multiplicative,
                we cannot declare that this map is a morphism for
                sage ``Groups()``.

                .. WARNING::

                    This method must be implemented by the "PW0" and "W0P" realizations.
                """
                PW0 = self.realization_of().PW0()
                return self(PW0.from_translation(la))

            def from_dual_translation(self, la):
                r"""
                Return the image of ``la`` under the homomorphism of the dual version of the
                translation lattice into ``self``.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1]); PvW0 = E.PvW0()
                    sage: bv = E.dual_lattice_basis(); bv
                    Finite family {1: Lambda[1], 2: Lambda[2]}
                    sage: x = PvW0.from_dual_translation(2*bv[1] - bv[2]); x
                    t[2*Lambda[1] - Lambda[2]]
                    sage: FW = E.FW()
                    sage: y = FW.from_dual_translation(2*bv[1] - bv[2]); y
                    S0*S2*S0*S1
                    sage: FW(x) == y
                    True
                """
                return self(self.realization_of().PvW0().from_dual_translation(la))

            @abstract_method
            def simple_reflections(self):
                r"""
                Return a family from the set of affine Dynkin nodes to the simple reflections
                in the realization of the extended affine Weyl group.

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).W0P().simple_reflections()
                    Finite family {0: s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]],
                                   1: s1, 2: s2, 3: s3}
                    sage: ExtendedAffineWeylGroup(['A',3,1]).WF().simple_reflections()
                    Finite family {0: S0, 1: S1, 2: S2, 3: S3}
                    sage: ExtendedAffineWeylGroup(['A',3,1],
                    ....:                         print_tuple=True).FW().simple_reflections()
                    Finite family {0: (pi[0], S0), 1: (pi[0], S1),
                                   2: (pi[0], S2), 3: (pi[0], S3)}
                    sage: ExtendedAffineWeylGroup(['A',3,1],
                    ....:                         fundamental='f',
                    ....:                         print_tuple=True).FW().simple_reflections()
                    Finite family {0: (f[0], S0), 1: (f[0], S1),
                                   2: (f[0], S2), 3: (f[0], S3)}
                    sage: ExtendedAffineWeylGroup(['A',3,1]).PvW0().simple_reflections()
                    Finite family {0: t[Lambda[1] + Lambda[3]] * s1*s2*s3*s2*s1,
                                   1: s1, 2: s2, 3: s3}
                """

            def simple_reflection(self, i):
                r"""
                Return the `i`-th simple reflection in ``self``.

                INPUT:

                - ``self`` -- a realization of the extended affine Weyl group
                - ``i`` -- an affine Dynkin node

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).PW0().simple_reflection(0)
                    t[Lambdacheck[1] + Lambdacheck[3]] * s1*s2*s3*s2*s1
                    sage: ExtendedAffineWeylGroup(['C',2,1]).WF().simple_reflection(0)
                    S0
                    sage: ExtendedAffineWeylGroup(['D',3,2]).PvW0().simple_reflection(1)
                    s1
                """
                return self.simple_reflections()[i]

            def from_classical_weyl(self, w):
                r"""
                Return the image of `w` from the finite Weyl group into ``self``.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); PW0=E.PW0()
                    sage: W0 = E.classical_weyl()
                    sage: w = W0.from_reduced_word([2,1,3])
                    sage: y = PW0.from_classical_weyl(w); y
                    s2*s3*s1
                    sage: y.parent() == PW0
                    True
                    sage: y.to_classical_weyl() == w
                    True
                    sage: W0P = E.W0P()
                    sage: z = W0P.from_classical_weyl(w); z
                    s2*s3*s1
                    sage: z.parent() == W0P
                    True
                    sage: W0P(y) == z
                    True
                    sage: FW = E.FW()
                    sage: x = FW.from_classical_weyl(w); x
                    S2*S3*S1
                    sage: x.parent() == FW
                    True
                    sage: FW(y) == x
                    True
                    sage: FW(z) == x
                    True

                .. WARNING::

                    Must be implemented in style "PW0" and "W0P".
                """
                PW0 = self.realization_of().PW0()
                return self(PW0.from_classical_weyl(w))

            def from_dual_classical_weyl(self, w):
                r"""
                Return the image of `w` from the finite Weyl group of dual form into ``self``.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); PvW0 = E.PvW0()
                    sage: W0v = E.dual_classical_weyl()
                    sage: w = W0v.from_reduced_word([2,1,3])
                    sage: y = PvW0.from_dual_classical_weyl(w); y
                    s2*s3*s1
                    sage: y.parent() == PvW0
                    True
                    sage: y.to_dual_classical_weyl() == w
                    True
                    sage: x = E.FW().from_dual_classical_weyl(w); x
                    S2*S3*S1
                    sage: PvW0(x) == y
                    True

                .. WARNING::

                    Must be implemented in style "PvW0" and "W0Pv".
                """
                return self(self.realization_of().PvW0().from_dual_classical_weyl(w))

            def from_affine_weyl(self, w):
                r"""
                Return the image of `w` under the homomorphism from the affine Weyl group
                into ``self``.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); PW0 = E.PW0()
                    sage: W = E.affine_weyl()
                    sage: w = W.from_reduced_word([2,1,3,0])
                    sage: x = PW0.from_affine_weyl(w); x
                    t[Lambdacheck[1] - 2*Lambdacheck[2] + Lambdacheck[3]] * s3*s1
                    sage: FW = E.FW()
                    sage: y = FW.from_affine_weyl(w); y
                    S2*S3*S1*S0
                    sage: FW(x) == y
                    True

                .. WARNING::

                    Must be implemented in style "WF" and "FW".
                """
                WF = self.realization_of().WF()
                return self(WF.from_affine_weyl(w))

            def from_reduced_word(self, word):
                r"""
                Convert an affine or finite reduced word into a group element.

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',2,1]).PW0().from_reduced_word([1,0,1,2])
                    t[-Lambdacheck[1] + 2*Lambdacheck[2]]
                """

                return self.from_affine_weyl(self.realization_of().affine_weyl().from_reduced_word(word))

        class ElementMethods:

            @abstract_method
            def has_descent(self, i, side='right', positive=False):
                r"""
                Return whether ``self`` * `s_i` < ``self`` where `s_i` is the `i`-th simple
                reflection in the realized group.

                INPUT:

                - ``i`` -- an affine Dynkin index

                OPTIONAL:

                - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``)
                - ``positive`` -- boolean (default: ``False``)

                If ``side='left'``, then the reflection acts
                on the left. If ``positive=True``, then the inequality is reversed.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); WF = E.WF()
                    sage: F = E.fundamental_group()
                    sage: x = WF.an_element(); x
                    S0*S1*S2*S3 * pi[3]
                    sage: I = E.cartan_type().index_set()
                    sage: [(i, x.has_descent(i)) for i in I]
                    [(0, True), (1, False), (2, False), (3, False)]
                    sage: [(i, x.has_descent(i,side='left')) for i in I]
                    [(0, True), (1, False), (2, False), (3, False)]
                    sage: [(i, x.has_descent(i,positive=True)) for i in I]
                    [(0, False), (1, True), (2, True), (3, True)]

                .. WARNING::

                    This method is abstract because it is used in the recursive coercions
                    between "PW0" and "WF" and other methods use this coercion.
                """

            def first_descent(self, side='right', positive=False, index_set=None):
                r"""
                Return the first descent of ``self``.

                INPUT:

                - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
                - ``positive`` -- boolean (default: ``False``)
                - ``index_set`` -- an optional subset of Dynkin nodes

                If ``index_set`` is not ``None``, then the descent must be in the ``index_set``.

                EXAMPLES::

                    sage: x = ExtendedAffineWeylGroup(['A',3,1]).WF().an_element(); x
                    S0*S1*S2*S3 * pi[3]
                    sage: x.first_descent()
                    0
                    sage: x.first_descent(side='left')
                    0
                    sage: x.first_descent(positive=True)
                    1
                    sage: x.first_descent(side='left',positive=True)
                    1
                """
                if index_set is None:
                    index_set = self.parent().realization_of().cartan_type().index_set()
                for i in index_set:
                    if self.has_descent(i, side=side, positive=positive):
                        return i
                return None

            def apply_simple_reflection(self, i, side='right'):
                r"""
                Apply the `i`-th simple reflection to ``self``.

                EXAMPLES::

                    sage: x = ExtendedAffineWeylGroup(['A',3,1]).WF().an_element(); x
                    S0*S1*S2*S3 * pi[3]
                    sage: x.apply_simple_reflection(1)
                    S0*S1*S2*S3*S0 * pi[3]
                    sage: x.apply_simple_reflection(0, side='left')
                    S1*S2*S3 * pi[3]
                """
                s = self.parent().simple_reflection(i)
                if side == 'right':
                    return self*s
                else:
                    return s*self

            def apply_simple_projection(self, i, side='right', length_increasing=True):
                r"""
                Return the product of ``self`` by the simple reflection `s_i` if that product is
                of greater length than ``self`` and otherwise return ``self``.

                INPUT:

                - ``self`` -- an element of the extended affine Weyl group
                - ``i`` -- a Dynkin node (index of a simple reflection `s_i`)
                - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``)
                  according to which side of ``self`` the reflection `s_i`
                  should be multiplied
                - ``length_increasing`` -- boolean (default: ``True``);
                  if ``False``, do the above with the word "greater" replaced by "less"

                EXAMPLES::

                    sage: x = ExtendedAffineWeylGroup(['A',3,1]).WF().an_element(); x
                    S0*S1*S2*S3 * pi[3]
                    sage: x.apply_simple_projection(1)
                    S0*S1*S2*S3*S0 * pi[3]
                    sage: x.apply_simple_projection(1, length_increasing=False)
                    S0*S1*S2*S3 * pi[3]
                """
                if self.has_descent(i, side=side, positive=length_increasing):
                    return self.apply_simple_reflection(i, side=side)
                return self

            def to_fundamental_group(self):
                r"""
                Return the image of ``self`` under the homomorphism to the fundamental group.

                EXAMPLES::

                    sage: PW0 = ExtendedAffineWeylGroup(['A',3,1]).PW0()
                    sage: b = PW0.realization_of().lattice_basis()
                    sage: [(x, PW0.from_translation(x).to_fundamental_group()) for x in b]
                    [(Lambdacheck[1], pi[1]), (Lambdacheck[2], pi[2]), (Lambdacheck[3], pi[3])]

                .. WARNING::

                    Must be implemented in style "WF".
                """
                WF = self.parent().realization_of().WF()
                return WF(self).to_fundamental_group()

            def to_classical_weyl(self):
                r"""
                Return the image of ``self`` under the homomorphism to the classical Weyl group.

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).WF().simple_reflection(0).to_classical_weyl()
                    s1*s2*s3*s2*s1

                .. WARNING::

                    Must be implemented in style "PW0".
                """
                PW0 = self.parent().realization_of().PW0()
                return PW0(self).to_classical_weyl()

            def to_dual_classical_weyl(self):
                r"""
                Return the image of ``self`` under the homomorphism to the dual form of the classical Weyl group.

                EXAMPLES::

                    sage: x = ExtendedAffineWeylGroup(['A',3,1]).WF().simple_reflection(0).to_dual_classical_weyl(); x
                    s1*s2*s3*s2*s1
                    sage: x.parent()
                    Weyl Group of type ['A', 3] (as a matrix group acting on the weight lattice)

                .. WARNING::

                    Must be implemented in style "PvW0".
                """
                PvW0 = self.parent().realization_of().PvW0()
                return PvW0(self).to_dual_classical_weyl()

            def to_affine_weyl_left(self):
                r"""
                Return the projection of ``self`` to the affine Weyl group on the left,
                after factorizing using the style "WF".

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); PW0 = E.PW0()
                    sage: b = E.lattice_basis()
                    sage: [(x, PW0.from_translation(x).to_affine_weyl_left()) for x in b]
                    [(Lambdacheck[1], S0*S3*S2),
                     (Lambdacheck[2], S0*S3*S1*S0),
                     (Lambdacheck[3], S0*S1*S2)]

                .. WARNING::

                    Must be implemented in style "WF".
                """
                WF = self.parent().realization_of().WF()
                return WF(self).to_affine_weyl_left()

            def to_affine_weyl_right(self):
                r"""
                Return the projection of ``self`` to the affine Weyl group on the right,
                after factorizing using the style "FW".

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); PW0 = E.PW0()
                    sage: b = E.lattice_basis()
                    sage: [(x, PW0.from_translation(x).to_affine_weyl_right()) for x in b]
                    [(Lambdacheck[1], S3*S2*S1),
                     (Lambdacheck[2], S2*S3*S1*S2),
                     (Lambdacheck[3], S1*S2*S3)]

                .. WARNING::

                    Must be implemented in style "FW".
                """
                FW = self.parent().realization_of().FW()
                return FW(self).to_affine_weyl_right()

            def to_translation_left(self):
                r"""
                Return the projection of ``self`` to the translation lattice after factorizing
                it to the left using the style "PW0".

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).PW0().simple_reflection(0).to_translation_left()
                    Lambdacheck[1] + Lambdacheck[3]

                .. WARNING::

                    Must be implemented in style "PW0".
                """
                PW0 = self.parent().realization_of().PW0()
                return PW0(self).to_translation_left()

            def to_translation_right(self):
                r"""
                Return the projection of ``self`` to the translation lattice after factorizing
                it to the right using the style "W0P".

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).PW0().simple_reflection(0).to_translation_right()
                    -Lambdacheck[1] - Lambdacheck[3]

                .. WARNING::

                    Must be implemented in style "W0P".
                """
                W0P = self.parent().realization_of().W0P()
                return W0P(self).to_translation_right()

            def to_dual_translation_left(self):
                r"""
                Return the projection of ``self`` to the dual translation lattice after factorizing
                it to the left using the style "PvW0".

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).PvW0().simple_reflection(0).to_dual_translation_left()
                    Lambda[1] + Lambda[3]

                .. WARNING::

                    Must be implemented in style "PvW0".
                """
                PvW0 = self.parent().realization_of().PvW0()
                return PvW0(self).to_dual_translation_left()

            def to_dual_translation_right(self):
                r"""
                Return the projection of ``self`` to the dual translation lattice after factorizing
                it to the right using the style "W0Pv".

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',3,1]).PW0().simple_reflection(0).to_dual_translation_right()
                    -Lambda[1] - Lambda[3]

                .. WARNING::

                    Must be implemented in style "W0Pv".
                """
                W0Pv = self.parent().realization_of().W0Pv()
                return W0Pv(self).to_dual_translation_right()

            def length(self):
                r"""
                Return the length of ``self`` in the Coxeter group sense.

                EXAMPLES::

                     sage: E = ExtendedAffineWeylGroup(['A',3,1]); PW0 = E.PW0()
                     sage: I0 = E.cartan_type().classical().index_set()
                     sage: [PW0.from_translation(E.lattice_basis()[i]).length() for i in I0]
                     [3, 4, 3]
                """
                return self.to_affine_weyl_left().length()

            def coset_representative(self, index_set, side='right'):
                r"""
                Return the minimum length representative in the coset of ``self`` with respect to
                the subgroup generated by the reflections given by ``index_set``.

                INPUT:

                - ``self`` -- an element of the extended affine Weyl group
                - ``index_set`` -- a subset of the set of Dynkin nodes
                - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``)
                  the side on which the subgroup acts

                EXAMPLES::

                     sage: E = ExtendedAffineWeylGroup(['A',3,1]); WF = E.WF()
                     sage: b = E.lattice_basis()
                     sage: I0 = E.cartan_type().classical().index_set()
                     sage: [WF.from_translation(x).coset_representative(index_set=I0) for x in b]
                     [pi[1], pi[2], pi[3]]
                """
                while True:
                    i = self.first_descent(index_set=index_set, side=side)
                    if i is None:
                        return self
                    self = self.apply_simple_reflection(i,side=side)

            def is_grassmannian(self, index_set, side='right'):
                r"""
                Return whether ``self`` is of minimum length in its coset with respect to the
                subgroup generated by the reflections of ``index_set``.

                EXAMPLES::

                     sage: E = ExtendedAffineWeylGroup(['A',3,1]); PW0 = E.PW0()
                     sage: x = PW0.from_translation(E.lattice_basis()[1]); x
                     t[Lambdacheck[1]]
                     sage: I = E.cartan_type().index_set()
                     sage: [(i, x.is_grassmannian(index_set=[i])) for i in I]
                     [(0, True), (1, False), (2, True), (3, True)]
                     sage: [(i, x.is_grassmannian(index_set=[i], side='left')) for i in I]
                     [(0, False), (1, True), (2, True), (3, True)]
                """
                return self == self.coset_representative(index_set=index_set,side=side)

            def to_affine_grassmannian(self):
                r"""
                Return the unique affine Grassmannian element in the same coset of ``self``
                with respect to the finite Weyl group acting on the right.

                EXAMPLES::

                     sage: elts = ExtendedAffineWeylGroup(['A',2,1]).PW0().some_elements()
                     sage: [(x, x.to_affine_grassmannian()) for x in elts]
                     [(t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2,
                       t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2*s1)]
                """

                return self.coset_representative(index_set=self.parent().realization_of().cartan_type().classical().index_set())

            def is_affine_grassmannian(self):
                r"""
                Return whether ``self`` is affine Grassmannian.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1]); PW0 = E.PW0()
                    sage: F = E.fundamental_group()
                    sage: [(x,PW0.from_fundamental(x).is_affine_grassmannian()) for x in F]
                    [(pi[0], True), (pi[1], True), (pi[2], True)]
                    sage: b = E.lattice_basis()
                    sage: [(-x,PW0.from_translation(-x).is_affine_grassmannian()) for x in b]
                    [(-Lambdacheck[1], True), (-Lambdacheck[2], True)]
                """

                return self == self.to_affine_grassmannian()

            def bruhat_le(self, x):
                r"""
                Return whether ``self`` <= `x` in Bruhat order.

                INPUT:

                - ``self`` -- an element of the extended affine Weyl group
                - ``x`` -- another element with the same parent as ``self``

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1], print_tuple=True); WF = E.WF()
                    sage: W = E.affine_weyl()
                    sage: v = W.from_reduced_word([2,1,0])
                    sage: w = W.from_reduced_word([2,0,1,0])
                    sage: v.bruhat_le(w)
                    True
                    sage: vx = WF.from_affine_weyl(v); vx
                    (S2*S1*S0, pi[0])
                    sage: wx = WF.from_affine_weyl(w); wx
                    (S2*S0*S1*S0, pi[0])
                    sage: vx.bruhat_le(wx)
                    True
                    sage: F = E.fundamental_group()
                    sage: f = WF.from_fundamental(F(2))
                    sage: vx.bruhat_le(wx*f)
                    False
                    sage: (vx*f).bruhat_le(wx*f)
                    True

                .. WARNING::

                    Must be implemented by "WF".
                """
                WF = self.parent().realization_of().WF()
                return WF(self).bruhat_le(WF(x))

            def is_translation(self):
                r"""
                Return whether ``self`` is a translation element or not.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1]); FW = E.FW()
                    sage: F = E.fundamental_group()
                    sage: FW.from_affine_weyl(E.affine_weyl().from_reduced_word([1,2,1,0])).is_translation()
                    True
                    sage: FW.from_translation(E.lattice_basis()[1]).is_translation()
                    True
                    sage: FW.simple_reflection(0).is_translation()
                    False
                """
                w = self.to_classical_weyl()
                return w == w.parent().one()

            def action(self, la):
                r"""
                Action of ``self`` on a lattice element ``la``.

                INPUT:

                - ``self`` -- an element of the extended affine Weyl group
                - ``la`` -- an element of the translation lattice of the extended
                  affine Weyl group, the lattice denoted by the mnemonic "P" in the
                  documentation for :meth:`ExtendedAffineWeylGroup`.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1], affine='s')
                    sage: x = E.FW().an_element(); x
                    pi[2] * s0*s1*s2
                    sage: la = E.lattice().an_element(); la
                    2*Lambdacheck[1] + 2*Lambdacheck[2]
                    sage: x.action(la)
                    5*Lambdacheck[1] - 3*Lambdacheck[2]
                    sage: E = ExtendedAffineWeylGroup(['C',2,1], affine='s')
                    sage: x = E.PW0().from_translation(E.lattice_basis()[1])
                    sage: x.action(E.lattice_basis()[2])
                    Lambdacheck[1] + Lambdacheck[2]

                .. WARNING::

                    Must be implemented by style "PW0".
                """
                PW0 = self.parent().realization_of().PW0()
                return PW0(self).action(la)

            def dual_action(self, la):
                r"""
                Action of ``self`` on a dual lattice element ``la``.

                INPUT:

                - ``self`` -- an element of the extended affine Weyl group
                - ``la`` -- an element of the dual translation lattice of the extended
                  affine Weyl group, the lattice denoted by the mnemonic "Pv" in
                  the documentation for :meth:`ExtendedAffineWeylGroup`.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1], affine='s')
                    sage: x = E.FW().an_element(); x
                    pi[2] * s0*s1*s2
                    sage: la = E.dual_lattice().an_element(); la
                    2*Lambda[1] + 2*Lambda[2]
                    sage: x.dual_action(la)
                    5*Lambda[1] - 3*Lambda[2]
                    sage: E = ExtendedAffineWeylGroup(['C',2,1], affine='s')
                    sage: x = E.PvW0().from_dual_translation(E.dual_lattice_basis()[1])
                    sage: x.dual_action(E.dual_lattice_basis()[2])
                    Lambda[1] + Lambda[2]

                .. WARNING::

                    Must be implemented by style "PvW0".
                """
                PvW0 = self.parent().realization_of().PvW0()
                return PvW0(self).dual_action(la)

            def action_on_affine_roots(self, beta):
                r"""
                Act by ``self`` on the affine root lattice element ``beta``.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',2,1])
                    sage: beta = E.cartan_type().root_system().root_lattice().an_element(); beta
                    2*alpha[0] + 2*alpha[1] + 3*alpha[2]
                    sage: x = E.FW().an_element(); x
                    pi[2] * S0*S1*S2
                    sage: x.action_on_affine_roots(beta)
                    alpha[0] + alpha[1]

                .. WARNING::

                    Must be implemented by style "FW".
                """
                E = self.parent().realization_of()
                assert beta in RootSystem(E.cartan_type()).root_lattice()
                return E.FW()(self).action_on_affine_roots(beta)

            def face_data(self, i):
                r"""
                Return a description of one of the bounding hyperplanes of the alcove of an extended affine Weyl group element.

                INPUT:

                - ``self`` -- an element of the extended affine Weyl group
                - ``i`` -- an affine Dynkin node

                OUTPUT:

                - A 2-tuple `(m,\beta)` defined as follows.

                ALGORITHM:

                Each element of the extended affine Weyl group corresponds to an alcove,
                and each alcove has a face for each affine Dynkin node. Given the data of ``self`` and `i`,
                let the extended affine Weyl group element ``self`` act on the affine simple root `\alpha_i`,
                yielding a real affine root, which can be expressed uniquely as

                .. MATH::

                    ``self`` \cdot \alpha_i = m \delta + \beta

                where `m` is an integer (the height of the `i`-th bounding hyperplane of the alcove of ``self``)
                and `\beta` is a classical root (the normal vector for the hyperplane which points towards the alcove).

                EXAMPLES::

                    sage: x = ExtendedAffineWeylGroup(['A',2,1]).PW0().an_element(); x
                    t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2
                    sage: x.face_data(0)
                    (-1, alpha[1])
                """
                Qaf = self.parent().realization_of().cartan_type().root_system().root_lattice()
                gamma = self.action_on_affine_roots(Qaf.simple_root(i))
                return gamma[0], Qaf.classical()(gamma)

            def alcove_walk_signs(self):
                r"""
                Return a signed alcove walk for ``self``.

                INPUT:

                - An element ``self`` of the extended affine Weyl group.

                OUTPUT:

                - A 3-tuple (`g`, ``rw``, ``signs``).

                ALGORITHM:

                The element ``self`` can be uniquely written  ``self`` = `g` * `w`
                where `g` has length zero and `w` is an element of the nonextended affine Weyl group.
                Let `w` have reduced word ``rw``.
                Starting with `g` and applying simple reflections from ``rw``, one obtains
                a sequence of extended affine Weyl group elements (that is, alcoves) and simple roots.
                The signs give the sequence of sides on which the alcoves lie, relative to the face
                indicated by the simple roots.

                EXAMPLES::

                    sage: E = ExtendedAffineWeylGroup(['A',3,1]); FW=E.FW()
                    sage: w = FW.from_reduced_word([0,2,1,3,0])*FW.from_fundamental(1); w
                    pi[1] * S3*S1*S2*S0*S3
                    sage: w.alcove_walk_signs()
                    (pi[1], [3, 1, 2, 0, 3], [-1, 1, -1, -1, 1])
                """
                We = self.parent()
                gw = We.realization_of().FW()(self)
                g = gw.cartesian_projection(0)
                w = gw.cartesian_projection(1)
                rw = w.reduced_word()
                u_curr = We.from_fundamental(g.value())
                signs = []
                for i in rw:
                    m, beta = u_curr.face_data(i)
                    if beta.is_positive_root():
                        signs = signs + [1]
                    else:
                        signs = signs + [-1]
                    u_curr = u_curr * We.simple_reflection(i)
                return g, rw, signs

    class ExtendedAffineWeylGroupPW0Element(GroupSemidirectProduct.Element):
        r"""
        The element class for the "PW0" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Return whether ``self`` has `i` as a descent.

            INPUT:

            - ``i`` -- an affine Dynkin node

            OPTIONAL:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: w = ExtendedAffineWeylGroup(['A',4,2]).PW0().from_reduced_word([0,1]); w
                t[Lambda[1]] * s1*s2
                sage: w.has_descent(0, side='left')
                True
            """

            E = self.parent().realization_of()
            if side == 'right':
                self = ~self
            if positive:
                return not self.has_descent(i, side='left')
            la = self.cartesian_projection(0).value
            w = self.cartesian_projection(1)
            if i == 0:
                ip = la.scalar(E._special_translation_covector) * E._a0check
                if ip > 1:
                    return True
                if ip < 1:
                    return False
                return E._special_root.weyl_action(w, inverse=True).is_positive_root()
            ip = la.scalar(E._simpleR0[i]) # test height versus simple (co)root
            if ip < 0:
                return True
            if ip > 0:
                return False
            return w.has_descent(i, side='left')

        def action(self, la):
            r"""
            Return the action of ``self`` on an element ``la`` of the translation lattice.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1]); PW0 = E.PW0()
                sage: x = PW0.an_element(); x
                t[2*Lambdacheck[1] + 2*Lambdacheck[2]] * s1*s2
                sage: la = E.lattice().an_element(); la
                2*Lambdacheck[1] + 2*Lambdacheck[2]
                sage: x.action(la)
                -2*Lambdacheck[1] + 4*Lambdacheck[2]
            """
            w = self.cartesian_projection(1)
            assert la in w.parent().domain()
            return self.cartesian_projection(0).value + w.action(la)

        def to_translation_left(self):
            r"""
            The image of ``self`` under the map that projects to the translation lattice
            factor after factoring it to the left as in style "PW0".

            EXAMPLES::

                sage: s = ExtendedAffineWeylGroup(['A',2,1]).PW0().S0(); s
                t[Lambdacheck[1] + Lambdacheck[2]] * s1*s2*s1
                sage: s.to_translation_left()
                Lambdacheck[1] + Lambdacheck[2]
            """
            return self.cartesian_projection(0).value # undo the GroupExp

        def to_classical_weyl(self):
            r"""
            Return the image of ``self`` under the homomorphism that projects to the classical
            Weyl group factor after rewriting it in either style "PW0" or "W0P".

            EXAMPLES::

                sage: s = ExtendedAffineWeylGroup(['A',2,1]).PW0().S0(); s
                t[Lambdacheck[1] + Lambdacheck[2]] * s1*s2*s1
                sage: s.to_classical_weyl()
                s1*s2*s1
            """
            return self.cartesian_projection(1)

    class ExtendedAffineWeylGroupPW0(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the translation lattice
        by the finite Weyl group.

        INPUT:

        - ``E`` -- a parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).PW0()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Multiplicative form of Coweight lattice of the Root system of type ['A', 2]
              acted upon by Weyl Group of type ['A', 2]
               (as a matrix group acting on the coweight lattice)
        """

        def __init__(self, E):
            r"""
            Create the PW0 realization of the extended affine Weyl group.

            EXAMPLES::

                sage: PW0 = ExtendedAffineWeylGroup(['D',3,2]).PW0()
                sage: TestSuite(PW0).run()
            """
            # note that we have to use the multiplicative version of the translation lattice
            # and change the twist to deal with this
            def twist(w, l):
                return E.exp_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.exp_lattice(), E.classical_weyl(), twist=twist, act_to_right=False, prefix0=E._prefixt, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "PW0"

        def _repr_(self):
            r"""
            A string representing ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',4,2]).PW0()._repr_()
                "Extended affine Weyl group of type ['BC', 2, 2] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['C', 2] acted upon by Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice)"
            """
            return self.realization_of()._repr_() + " realized by " + super()._repr_()

        def from_translation(self, la):
            r"""
            Map the translation lattice element ``la`` into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], translation='tau',
                ....:                             print_tuple=True)
                sage: la = E.lattice().an_element(); la
                2*Lambdacheck[1] + 2*Lambdacheck[2]
                sage: E.PW0().from_translation(la)
                (tau[2*Lambdacheck[1] + 2*Lambdacheck[2]], 1)
            """
            E = self.realization_of()
            return self((E.exp_lattice()(la),self.cartesian_factors()[1].one()))

        @cached_method
        def S0(self):
            r"""
            Return the affine simple reflection.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['B',2]).PW0().S0()
                t[Lambdacheck[2]] * s2*s1*s2
            """
            E = self.realization_of()
            return self((E.exp_lattice()(E.lattice()(E._special_translation)), E._special_reflection))

        @cached_method
        def simple_reflection(self, i):
            r"""
            Return the `i`-th simple reflection in ``self``.

            EXAMPLES::

               sage: E = ExtendedAffineWeylGroup("G2")
               sage: [(i, E.PW0().simple_reflection(i)) for i in E.cartan_type().index_set()]
               [(0, t[Lambdacheck[2]] * s2*s1*s2*s1*s2), (1, s1), (2, s2)]
            """
            if i == 0:
                return self.S0()
            else:
                E = self.realization_of()
                return self.from_classical_weyl(E.classical_weyl().simple_reflection(i))

        @cached_method
        def simple_reflections(self):
            r"""
            Return a family for the simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup("A3").PW0().simple_reflections()
                Finite family {0: t[Lambdacheck[1] + Lambdacheck[3]] * s1*s2*s3*s2*s1,
                               1: s1, 2: s2, 3: s3}
            """
            return Family(self.realization_of().cartan_type().index_set(), self.simple_reflection)

        def from_classical_weyl(self, w):
            r"""
            Return the image of `w` under the homomorphism of the classical Weyl group into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup("A3",print_tuple=True)
                sage: E.PW0().from_classical_weyl(E.classical_weyl().from_reduced_word([1,2]))
                (t[0], s1*s2)
            """
            return self((self.cartesian_factors()[0].one(),w))

    class ExtendedAffineWeylGroupW0PElement(GroupSemidirectProduct.Element):
        r"""
        The element class for the W0P realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Return whether ``self`` has `i` as a descent.

            INPUT:

            - ``i`` -- an index

            OPTIONAL:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: W0P = ExtendedAffineWeylGroup(['A',4,2]).W0P()
                sage: w = W0P.from_reduced_word([0,1]); w
                s1*s2 * t[Lambda[1] - Lambda[2]]
                sage: w.has_descent(0, side='left')
                True
            """
            E = self.parent().realization_of()
            if side == 'left':
                self = ~self
            if positive:
                return not self.has_descent(i, side='right')
            w = self.cartesian_projection(0)
            la = self.cartesian_projection(1).value
            if i == 0:
                ip = la.scalar(E._special_translation_covector) * E._a0check
                if ip < -1:
                    return True
                if ip > -1:
                    return False
                return E._special_root.weyl_action(w).is_positive_root()
            ip = la.scalar(E._simpleR0[i]) # test height versus simple (co)root
            if ip > 0:
                return True
            if ip < 0:
                return False
            return w.has_descent(i, side='right')

        def to_classical_weyl(self):
            r"""
            Project ``self`` into the classical Weyl group.

            EXAMPLES::

                sage: x = ExtendedAffineWeylGroup(['A',2,1]).W0P().simple_reflection(0); x
                s1*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[2]]
                sage: x.to_classical_weyl()
                s1*s2*s1
            """
            return self.cartesian_projection(0)

        def to_translation_right(self):
            r"""
            Project onto the right (translation) factor in the "W0P" style.

            EXAMPLES::

                sage: x = ExtendedAffineWeylGroup(['A',2,1]).W0P().simple_reflection(0); x
                s1*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[2]]
                sage: x.to_translation_right()
                -Lambdacheck[1] - Lambdacheck[2]
            """
            return self.cartesian_projection(1).value

    class ExtendedAffineWeylGroupW0P(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the finite Weyl group
        by the translation lattice.

        INPUT:

        - ``E`` -- a parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).W0P()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Weyl Group of type ['A', 2] (as a matrix group acting on the coweight lattice)
             acting on Multiplicative form of Coweight lattice of the Root system of type ['A', 2]
        """

        def __init__(self, E):
            r"""
            EXAMPLES::

                sage: W0P = ExtendedAffineWeylGroup(['D',3,2]).W0P()
                sage: TestSuite(W0P).run()
            """
            def twist(w, l):
                return E.exp_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.classical_weyl(), E.exp_lattice(), twist=twist, act_to_right=True, prefix1=E._prefixt, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "W0P"

        def _repr_(self):
            r"""
            A string representing ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',4,2]).W0P()._repr_()
                "Extended affine Weyl group of type ['BC', 2, 2] realized by Semidirect product of Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice) acting on Multiplicative form of Weight lattice of the Root system of type ['C', 2]"
            """
            return self.realization_of()._repr_() + " realized by " + super()._repr_()

        def S0(self):
            r"""
            Return the zero-th simple reflection in style "W0P".

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(["A",3,1]).W0P().S0()
                s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]]
            """
            E = self.realization_of()
            return self((E._special_reflection,E.exp_lattice()(E.lattice()(-E._special_translation))))

        def simple_reflection(self, i):
            r"""
            Return the `i`-th simple reflection in ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',3,1]); W0P = E.W0P()
                sage: [(i, W0P.simple_reflection(i)) for i in E.cartan_type().index_set()]
                [(0, s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]]),
                 (1, s1), (2, s2), (3, s3)]
            """
            if i == 0:
                return self.S0()
            E = self.realization_of()
            return self.from_classical_weyl(E.classical_weyl().simple_reflection(i))

        @cached_method
        def simple_reflections(self):
            r"""
            Return the family of simple reflections.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(["A",3,1]).W0P().simple_reflections()
                Finite family {0: s1*s2*s3*s2*s1 * t[-Lambdacheck[1] - Lambdacheck[3]],
                               1: s1, 2: s2, 3: s3}
            """
            return Family(self.realization_of().cartan_type().index_set(), self.simple_reflection)

        def from_classical_weyl(self, w):
            r"""
            Return the image of the classical Weyl group element `w` in ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1],print_tuple=True)
                sage: E.W0P().from_classical_weyl(E.classical_weyl().from_reduced_word([2,1]))
                (s2*s1, t[0])
            """
            return self((w,self.cartesian_factors()[1].one()))

        def from_translation(self, la):
            r"""
            Return the image of the lattice element ``la`` in ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1],print_tuple=True)
                sage: E.W0P().from_translation(E.lattice().an_element())
                (1, t[2*Lambdacheck[1] + 2*Lambdacheck[2]])
            """
            return self((self.cartesian_factors()[0].one(),self.realization_of().exp_lattice()(la)))

    class ExtendedAffineWeylGroupWFElement(GroupSemidirectProduct.Element):
        r"""
        Element class for the "WF" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Return whether ``self`` has `i` as a descent.

            INPUT:

            - ``i`` -- an affine Dynkin index

            OPTIONAL:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.WF().an_element(); x
                S0*S1*S2 * pi[2]
                sage: [(i, x.has_descent(i)) for i in E.cartan_type().index_set()]
                [(0, True), (1, False), (2, False)]
            """
            if side == 'right':
                self = ~self
            if positive:
                return not self.has_descent(i, side='left')
            return self.cartesian_projection(0).has_descent(i, side='left')

        def to_fundamental_group(self):
            r"""
            Project ``self`` to the right (fundamental group) factor in the "WF" style.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.WF().from_translation(E.lattice_basis()[1]); x
                S0*S2 * pi[1]
                sage: x.to_fundamental_group()
                pi[1]
            """
            return self.cartesian_projection(1)

        def to_affine_weyl_left(self):
            r"""
            Project ``self`` to the left (affine Weyl group) factor in the "WF" style.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.WF().from_translation(E.lattice_basis()[1]); x
                S0*S2 * pi[1]
                sage: x.to_affine_weyl_left()
                S0*S2
            """
            return self.cartesian_projection(0)

        def bruhat_le(self, x):
            r"""
            Return whether ``self`` is less than or equal to `x` in the Bruhat order.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], affine='s',
                ....:                             print_tuple=True); WF = E.WF()
                sage: r = E.affine_weyl().from_reduced_word
                sage: v = r([1,0])
                sage: w = r([1,2,0])
                sage: v.bruhat_le(w)
                True
                sage: vv = WF.from_affine_weyl(v); vv
                (s1*s0, pi[0])
                sage: ww = WF.from_affine_weyl(w); ww
                (s1*s2*s0, pi[0])
                sage: vv.bruhat_le(ww)
                True
                sage: f = E.fundamental_group()(2); f
                pi[2]
                sage: ff = WF.from_fundamental(f); ff
                (1, pi[2])
                sage: vv.bruhat_le(ww*ff)
                False
                sage: (vv*ff).bruhat_le(ww*ff)
                True
            """
            if self.cartesian_projection(1) != x.cartesian_projection(1):
                return False
            return self.cartesian_projection(0).bruhat_le(x.cartesian_projection(0))

    class ExtendedAffineWeylGroupWF(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the affine Weyl group
        by the fundamental group.

        INPUT:

        - ``E`` -- a parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).WF()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
              acted upon by Fundamental group of type ['A', 2, 1]
        """

        def __init__(self, E):
            r"""
            EXAMPLES::

                sage: WF = ExtendedAffineWeylGroup(['D',3,2]).WF()
                sage: TestSuite(WF).run()
            """
            def twist(g, w):
                return g.act_on_affine_weyl(w)

            GroupSemidirectProduct.__init__(self, E.affine_weyl(), E.fundamental_group(), twist=twist, act_to_right=False, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "WF"

        def _repr_(self):
            r"""
            A string representing ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',4,2]).WF()._repr_()
                "Extended affine Weyl group of type ['BC', 2, 2] realized by Semidirect product of Weyl Group of type ['BC', 2, 2] (as a matrix group acting on the root lattice) acted upon by Fundamental group of type ['BC', 2, 2]"
            """
            return self.realization_of()._repr_() + " realized by " + super()._repr_()

        def from_affine_weyl(self, w):
            r"""
            Return the image of the affine Weyl group element `w` in ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['C',2,1], print_tuple=True)
                sage: E.WF().from_affine_weyl(E.affine_weyl().from_reduced_word([1,2,1,0]))
                (S1*S2*S1*S0, pi[0])
            """
            return self((w,self.cartesian_factors()[1].one()))

        @cached_method
        def simple_reflections(self):
            r"""
            Return the family of simple reflections.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(["A",3,1], affine='r').WF().simple_reflections()
                Finite family {0: r0, 1: r1, 2: r2, 3: r3}
            """
            E = self.realization_of()
            W = E.affine_weyl()
            return Family(E.cartan_type().index_set(), lambda i: self.from_affine_weyl(W.simple_reflection(i)))

        @cached_method
        def from_fundamental(self, f):
            r"""
            Return the image of `f` under the homomorphism from the fundamental group into
            the right (fundamental group) factor in "WF" style.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['E',6,1], print_tuple=True); WF = E.WF()
                sage: F = E.fundamental_group()
                sage: [(x, WF.from_fundamental(x)) for x in F]
                [(pi[0], (1, pi[0])), (pi[1], (1, pi[1])), (pi[6], (1, pi[6]))]
            """
            return self((self.cartesian_factors()[0].one(),f))

    class ExtendedAffineWeylGroupFWElement(GroupSemidirectProduct.Element):
        r"""
        The element class for the "FW" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Return whether ``self`` has descent at `i`.

            INPUT:

            - ``i`` -- an affine Dynkin index

            OPTIONAL:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.FW().an_element(); x
                pi[2] * S0*S1*S2
                sage: [(i, x.has_descent(i)) for i in E.cartan_type().index_set()]
                [(0, False), (1, False), (2, True)]
            """
            if side == 'left':
                self = ~self
            if positive:
                return not self.has_descent(i, side='right')
            return self.cartesian_projection(1).has_descent(i, side='right')

        def to_fundamental_group(self):
            r"""
            Return the projection of ``self`` to the fundamental group in the "FW" style.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.FW().from_translation(E.lattice_basis()[2]); x
                pi[2] * S1*S2
                sage: x.to_fundamental_group()
                pi[2]
            """
            return self.cartesian_projection(0)

        def to_affine_weyl_right(self):
            r"""
            Project ``self`` to the right (affine Weyl group) factor in the "FW" style.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.FW().from_translation(E.lattice_basis()[1]); x
                pi[1] * S2*S1
                sage: x.to_affine_weyl_right()
                S2*S1
            """
            return self.cartesian_projection(1)

        def action_on_affine_roots(self, beta):
            r"""
            Act by ``self`` on the affine root lattice element ``beta``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], affine='s')
                sage: x = E.FW().an_element(); x
                pi[2] * s0*s1*s2
                sage: v = RootSystem(['A',2,1]).root_lattice().an_element(); v
                2*alpha[0] + 2*alpha[1] + 3*alpha[2]
                sage: x.action_on_affine_roots(v)
                alpha[0] + alpha[1]
            """
            g = self.cartesian_projection(0)
            w = self.cartesian_projection(1)
            return g.act_on_affine_lattice(w.action(beta))

    class ExtendedAffineWeylGroupFW(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the affine Weyl group
        by the fundamental group.

        INPUT:

        - ``E`` -- a parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).FW()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Fundamental group of type ['A', 2, 1] acting on Weyl Group of type ['A', 2, 1]
             (as a matrix group acting on the root lattice)
        """

        def __init__(self, E):
            r"""

            EXAMPLES::

                sage: FW = ExtendedAffineWeylGroup(['D',3,2]).FW()
                sage: TestSuite(FW).run()
            """
            def twist(g, w):
                return g.act_on_affine_weyl(w)

            GroupSemidirectProduct.__init__(self, E.fundamental_group(), E.affine_weyl(), twist=twist, act_to_right=True, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "FW"

        def _repr_(self):
            r"""
            A string representing ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',4,2]).FW()._repr_()
                "Extended affine Weyl group of type ['BC', 2, 2] realized by Semidirect product of Fundamental group of type ['BC', 2, 2] acting on Weyl Group of type ['BC', 2, 2] (as a matrix group acting on the root lattice)"
            """
            return self.realization_of()._repr_() + " realized by " + super()._repr_()

        @cached_method
        def simple_reflections(self):
            r"""
            Return the family of simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',2,1],print_tuple=True).FW().simple_reflections()
                Finite family {0: (pi[0], S0), 1: (pi[0], S1), 2: (pi[0], S2)}
            """
            E = self.realization_of()
            W = E.affine_weyl()
            return Family(E.cartan_type().index_set(), lambda i: self.from_affine_weyl(W.simple_reflection(i)))

        def from_affine_weyl(self, w):
            r"""
            Return the image of `w` under the map of the affine Weyl group into the right
            (affine Weyl group) factor in the "FW" style.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], print_tuple=True)
                sage: E.FW().from_affine_weyl(E.affine_weyl().from_reduced_word([0,2,1]))
                (pi[0], S0*S2*S1)
            """
            return self((self.cartesian_factors()[0].one(),w))

        @cached_method
        def from_fundamental(self, f):
            r"""
            Return the image of the fundamental group element `f` into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], print_tuple=True)
                sage: E.FW().from_fundamental(E.fundamental_group()(2))
                (pi[2], 1)
            """
            return self((f,self.cartesian_factors()[1].one()))

    class ExtendedAffineWeylGroupPvW0Element(GroupSemidirectProduct.Element):
        r"""
        The element class for the "PvW0" realization.
        """

        def has_descent(self, i, side='right', positive=False):
            r"""
            Return whether ``self`` has `i` as a descent.

            INPUT:

            - ``i`` -- an affine Dynkin index

            OPTIONAL:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',4,2])
                sage: w = E.PvW0().from_reduced_word([0,1]); w
                t[Lambda[1]] * s1*s2
                sage: [(i, w.has_descent(i, side='left')) for i in E.cartan_type().index_set()]
                [(0, True), (1, False), (2, False)]
            """
            return self.parent().realization_of().PW0()(self).has_descent(i, side=side, positive=positive)

        def dual_action(self, la):
            r"""
            Return the action of ``self`` on an element ``la`` of the dual version of the translation lattice.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.PvW0().an_element(); x
                t[2*Lambda[1] + 2*Lambda[2]] * s1*s2
                sage: la = E.dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: x.dual_action(la)
                -2*Lambda[1] + 4*Lambda[2]
            """
            w = self.cartesian_projection(1)
            assert la in w.parent().domain()
            return self.cartesian_projection(0).value + w.action(la)

        def to_dual_translation_left(self):
            r"""
            The image of ``self`` under the map that projects to the dual translation lattice
            factor after factoring it to the left as in style "PvW0".

            EXAMPLES::

                sage: s = ExtendedAffineWeylGroup(['A',2,1]).PvW0().simple_reflection(0); s
                t[Lambda[1] + Lambda[2]] * s1*s2*s1
                sage: s.to_dual_translation_left()
                Lambda[1] + Lambda[2]
            """
            return self.cartesian_projection(0).value # undo the GroupExp

        def to_dual_classical_weyl(self):
            r"""
            Return the image of ``self`` under the homomorphism that projects to the dual classical
            Weyl group factor after rewriting it in either style "PvW0" or "W0Pv".

            EXAMPLES::

                sage: s = ExtendedAffineWeylGroup(['A',2,1]).PvW0().simple_reflection(0); s
                t[Lambda[1] + Lambda[2]] * s1*s2*s1
                sage: s.to_dual_classical_weyl()
                s1*s2*s1
            """
            return self.cartesian_projection(1)

            def is_translation(self):
                r"""
                Return whether ``self`` is a translation element or not.

                EXAMPLES::

                    sage: PvW0 = ExtendedAffineWeylGroup(['A',2,1]).PvW0()
                    sage: t = PvW0.from_reduced_word([1,2,1,0])
                    sage: t.is_translation()
                    True
                    sage: PvW0.simple_reflection(0).is_translation()
                    False
                """
                w = self.to_dual_classical_weyl()
                return w == w.parent().one()

    class ExtendedAffineWeylGroupPvW0(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the dual form of the translation lattice
        by the finite Weyl group.

        INPUT:

        - ``E`` -- a parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).PvW0()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Multiplicative form of Weight lattice of the Root system of type ['A', 2] acted
              upon by Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
        """

        def __init__(self, E):
            r"""

            EXAMPLES::

                sage: PvW0 = ExtendedAffineWeylGroup(['D',3,2]).PvW0()
                sage: TestSuite(PvW0).run()
            """
            # note that we have to use the multiplicative version of the translation lattice
            # and change the twist to deal with this
            def twist(w, l):
                return E.exp_dual_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.exp_dual_lattice(), E.dual_classical_weyl(), twist=twist, act_to_right=False, prefix0=E._prefixt, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "PvW0"

        def _repr_(self):
            r"""
            A string representing ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',4,2]).PvW0()._repr_()
                "Extended affine Weyl group of type ['BC', 2, 2] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['C', 2] acted upon by Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice)"
            """
            return self.realization_of()._repr_() + " realized by " + super()._repr_()

        def from_dual_translation(self, la):
            r"""
            Map the dual translation lattice element ``la`` into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], translation='tau',
                ....:                             print_tuple=True)
                sage: la = E.dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: E.PvW0().from_dual_translation(la)
                (tau[2*Lambda[1] + 2*Lambda[2]], 1)
            """
            E = self.realization_of()
            return self((E.exp_dual_lattice()(la),self.cartesian_factors()[1].one()))

        @cached_method
        def simple_reflections(self):
            r"""
            Return a family for the simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',3,1]).PvW0().simple_reflections()
                Finite family {0: t[Lambda[1] + Lambda[3]] * s1*s2*s3*s2*s1,
                               1: s1, 2: s2, 3: s3}
            """
            E = self.realization_of()
            return Family(E.cartan_type().index_set(), lambda i: self(E.PW0().simple_reflection(i)))

        def from_dual_classical_weyl(self, w):
            r"""
            Return the image of `w` under the homomorphism of the dual form of the classical Weyl group into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',3,1], print_tuple=True)
                sage: E.PvW0().from_dual_classical_weyl(
                ....:     E.dual_classical_weyl().from_reduced_word([1,2]))
                (t[0], s1*s2)
            """
            return self((self.cartesian_factors()[0].one(),w))

    class ExtendedAffineWeylGroupW0PvElement(GroupSemidirectProduct.Element):
        r"""
        The element class for the "W0Pv" realization.
        """

        def dual_action(self, la):
            r"""
            Return the action of ``self`` on an element ``la`` of the dual version of the translation lattice.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1])
                sage: x = E.W0Pv().an_element(); x
                s1*s2 * t[2*Lambda[1] + 2*Lambda[2]]
                sage: la = E.dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: x.dual_action(la)
                -8*Lambda[1] + 4*Lambda[2]
            """
            w = self.cartesian_projection(0)
            assert la in w.parent().domain()
            return w.action(self.cartesian_projection(1).value + la)

        def has_descent(self, i, side='right', positive=False):
            r"""
            Return whether ``self`` has `i` as a descent.

            INPUT:

            - ``i`` -- an affine Dynkin index

            OPTIONAL:

            - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: w = ExtendedAffineWeylGroup(['A',4,2]).W0Pv().from_reduced_word([0,1]); w
                s1*s2 * t[Lambda[1] - Lambda[2]]
                sage: w.has_descent(0, side='left')
                True
            """
            return self.parent().realization_of().W0P()(self).has_descent(i, side=side, positive=positive)

        def to_dual_translation_right(self):
            r"""
            The image of ``self`` under the map that projects to the dual translation lattice
            factor after factoring it to the right as in style "W0Pv".

            EXAMPLES::

                sage: s = ExtendedAffineWeylGroup(['A',2,1]).W0Pv().simple_reflection(0); s
                s1*s2*s1 * t[-Lambda[1] - Lambda[2]]
                sage: s.to_dual_translation_right()
                -Lambda[1] - Lambda[2]
            """
            return self.cartesian_projection(1).value # undo the GroupExp

        def to_dual_classical_weyl(self):
            r"""
            Return the image of ``self`` under the homomorphism that projects to the dual classical
            Weyl group factor after rewriting it in either style "PvW0" or "W0Pv".

            EXAMPLES::

                sage: s = ExtendedAffineWeylGroup(['A',2,1]).W0Pv().simple_reflection(0); s
                s1*s2*s1 * t[-Lambda[1] - Lambda[2]]
                sage: s.to_dual_classical_weyl()
                s1*s2*s1
            """
            return self.cartesian_projection(0)

            def is_translation(self):
                r"""
                Return whether ``self`` is a translation element or not.

                EXAMPLES::

                    sage: ExtendedAffineWeylGroup(['A',2,1]).W0Pv().from_reduced_word([1,2,1,0]).is_translation()
                    True
                """
                w = self.to_dual_classical_weyl()
                return w == w.parent().one()

    class ExtendedAffineWeylGroupW0Pv(GroupSemidirectProduct, BindableClass):
        r"""
        Extended affine Weyl group, realized as the semidirect product of the finite Weyl group, acting on the
        dual form of the translation lattice.

        INPUT:

        - ``E`` -- a parent with realization in :class:`ExtendedAffineWeylGroup_Class`

        EXAMPLES::

            sage: ExtendedAffineWeylGroup(['A',2,1]).W0Pv()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of
             Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
             acting on Multiplicative form of Weight lattice of the Root system of type ['A', 2]
        """

        def __init__(self, E):
            r"""
            EXAMPLES::

                sage: W0Pv = ExtendedAffineWeylGroup(['D',3,2]).W0Pv()
                sage: TestSuite(W0Pv).run()
            """
            # note that we have to use the multiplicative version of the translation lattice
            # and change the twist to deal with this
            def twist(w, l):
                return E.exp_dual_lattice()(w.action(l.value))

            GroupSemidirectProduct.__init__(self, E.dual_classical_weyl(), E.exp_dual_lattice(), twist=twist, act_to_right=True, prefix1=E._prefixt, print_tuple=E._print_tuple, category=E.Realizations())
            self._style = "W0Pv"

        def _repr_(self):
            r"""
            A string representing ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',4,2]).W0Pv()._repr_()
                "Extended affine Weyl group of type ['BC', 2, 2] realized by Semidirect product of Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice) acting on Multiplicative form of Weight lattice of the Root system of type ['C', 2]"
            """
            return self.realization_of()._repr_() + " realized by " + super()._repr_()

        def from_dual_translation(self, la):
            r"""
            Map the dual translation lattice element ``la`` into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',2,1], translation='tau',
                ....:                             print_tuple=True)
                sage: la = E.dual_lattice().an_element(); la
                2*Lambda[1] + 2*Lambda[2]
                sage: E.W0Pv().from_dual_translation(la)
                (1, tau[2*Lambda[1] + 2*Lambda[2]])
            """
            E = self.realization_of()
            return self((self.cartesian_factors()[0].one(),E.exp_dual_lattice()(la)))

        @cached_method
        def simple_reflections(self):
            r"""
            Return a family for the simple reflections of ``self``.

            EXAMPLES::

                sage: ExtendedAffineWeylGroup(['A',3,1]).W0Pv().simple_reflections()
                Finite family {0: s1*s2*s3*s2*s1 * t[-Lambda[1] - Lambda[3]],
                               1: s1, 2: s2, 3: s3}
            """
            E = self.realization_of()
            return Family(E.cartan_type().index_set(), lambda i: self(E.PW0().simple_reflection(i)))

        def from_dual_classical_weyl(self, w):
            r"""
            Return the image of `w` under the homomorphism of the dual form of the classical Weyl group into ``self``.

            EXAMPLES::

                sage: E = ExtendedAffineWeylGroup(['A',3,1],print_tuple=True)
                sage: E.W0Pv().from_dual_classical_weyl(E.dual_classical_weyl().from_reduced_word([1,2]))
                (s1*s2, t[0])
            """
            return self((w,self.cartesian_factors()[1].one()))


ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPW0.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPW0Element
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0P.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0PElement
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupWF.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupWFElement
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupFW.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupFWElement
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPvW0.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupPvW0Element
ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0Pv.Element = ExtendedAffineWeylGroup_Class.ExtendedAffineWeylGroupW0PvElement
