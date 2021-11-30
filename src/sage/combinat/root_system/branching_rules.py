"""
Branching Rules
"""
# ****************************************************************************
#  Copyright (C) 2014 Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.combinat.root_system.weyl_characters
from sage.combinat.root_system.root_system import RootSystem
from sage.matrix.constructor import matrix
from sage.misc.flatten import flatten
from sage.structure.sage_object import SageObject
from sage.combinat.root_system.cartan_type import CartanType
from sage.modules.free_module_element import vector
from sage.rings.rational_field import QQ
from sage.misc.functional import is_even, is_odd


def branch_weyl_character(chi, R, S, rule="default"):
    r"""
    A branching rule describes the restriction of representations from
    a Lie group or algebra `G` to a subgroup `H`. See for example, R. C.
    King, Branching rules for classical Lie groups using tensor and
    spinor methods. J. Phys. A 8 (1975), 429-449, Howe, Tan and
    Willenbring, Stable branching rules for classical symmetric pairs,
    Trans. Amer. Math. Soc. 357 (2005), no. 4, 1601-1626, McKay and
    Patera, Tables of Dimensions, Indices and Branching Rules for
    Representations of Simple Lie Algebras (Marcel Dekker, 1981),
    and Fauser, Jarvis, King and Wybourne, New branching rules induced
    by plethysm. J. Phys. A 39 (2006), no. 11, 2611--2655. If `H\subset G`
    we will write `G\Rightarrow H` to denote the branching rule, which
    is a homomorphism of WeylCharacterRings.

    INPUT:

    - ``chi`` -- a character of `G`

    - ``R`` -- the Weyl Character Ring of `G`

    - ``S`` -- the Weyl Character Ring of `H`

    - ``rule`` -- an element of the ``BranchingRule`` class
      or one (most usually) a keyword such as:

      * ``"levi"``
      * ``"automorphic"``
      * ``"symmetric"``
      * ``"extended"``
      * ``"orthogonal_sum"``
      * ``"tensor"``
      * ``"triality"``
      * ``"miscellaneous"``

    The ``BranchingRule`` class is a wrapper for functions
    from the weight lattice of `G` to the weight lattice of `H`.
    An instance of this class encodes an embedding of `H` into
    `G`. The usual way to specify an embedding is to supply a
    keyword, which tells Sage to use one of the built-in rules.
    We will discuss these first.

    To explain the predefined rules, we survey the most important
    branching rules. These may be classified into several cases, and
    once this is understood, the detailed classification can be read
    off from the Dynkin diagrams. Dynkin classified the maximal
    subgroups of Lie groups in Mat. Sbornik N.S. 30(72):349-462 (1952).

    We will list give predefined rules that cover most cases where the
    branching rule is to a maximal subgroup. For convenience, we
    also give some branching rules to subgroups that are not maximal.
    For example, a Levi subgroup may or may not be maximal.

    For example, there is a "levi" branching rule defined from `SL(5)` (with
    Cartan type `A_4`) to `SL(4)` (with Cartan type `A_3`), so
    we may compute the branching rule as follows:

    EXAMPLES::

        sage: A3=WeylCharacterRing("A3",style="coroots")
        sage: A2=WeylCharacterRing("A2",style="coroots")
        sage: [A3(fw).branch(A2,rule="levi") for fw in A3.fundamental_weights()]
        [A2(0,0) + A2(1,0), A2(0,1) + A2(1,0), A2(0,0) + A2(0,1)]

    In this case the Levi branching rule is the default branching rule
    so we may omit the specification rule="levi".

    If a subgroup is not maximal, you may specify a branching rule
    by finding a chain of intermediate subgroups. For this
    purpose, branching rules may be multiplied as in the following
    example.

    EXAMPLES::

        sage: A4=WeylCharacterRing("A4",style="coroots")
        sage: A2=WeylCharacterRing("A2",style="coroots")
        sage: br=branching_rule("A4","A3")*branching_rule("A3","A2")
        sage: A4(1,0,0,0).branch(A2,rule=br)
        2*A2(0,0) + A2(1,0)

    You may try omitting the rule if it is "obvious". Default
    rules are provided for the following cases:

    .. MATH::

        \begin{aligned}
        A_{2s} & \Rightarrow B_s,
        \\ A_{2s-1} & \Rightarrow C_s,
        \\ A_{2*s-1} & \Rightarrow D_s.
        \end{aligned}

    The above default rules correspond to embedding the group
    `SO(2s+1)`, `Sp(2s)` or `SO(2s)` into the corresponding general
    or special linear group by the standard representation. Default
    rules are also specified for the following cases:

    .. MATH::

        \begin{aligned}
        B_{s+1} & \Rightarrow D_s,
        \\ D_s & \Rightarrow B_s.
        \end{aligned}

    These correspond to the embedding of `O(n)` into `O(n+1)` where
    `n = 2s` or `2s + 1`. Finally, the branching rule for the embedding of
    a Levi subgroup is also implemented as a default rule.

    EXAMPLES::

        sage: A1 = WeylCharacterRing("A1", style="coroots")
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: D4 = WeylCharacterRing("D4", style="coroots")
        sage: B3 = WeylCharacterRing("B3", style="coroots")
        sage: B4 = WeylCharacterRing("B4", style="coroots")
        sage: A6 = WeylCharacterRing("A6", style="coroots")
        sage: A7 = WeylCharacterRing("A7", style="coroots")
        sage: def try_default_rule(R,S): return [R(f).branch(S) for f in R.fundamental_weights()]
        sage: try_default_rule(A2,A1)
        [A1(0) + A1(1), A1(0) + A1(1)]
        sage: try_default_rule(D4,B3)
        [B3(0,0,0) + B3(1,0,0), B3(1,0,0) + B3(0,1,0), B3(0,0,1), B3(0,0,1)]
        sage: try_default_rule(B4,D4)
        [D4(0,0,0,0) + D4(1,0,0,0), D4(1,0,0,0) + D4(0,1,0,0),
        D4(0,1,0,0) + D4(0,0,1,1), D4(0,0,1,0) + D4(0,0,0,1)]
        sage: try_default_rule(A7,D4)
        [D4(1,0,0,0), D4(0,1,0,0), D4(0,0,1,1), D4(0,0,2,0) + D4(0,0,0,2),
        D4(0,0,1,1),
        D4(0,1,0,0),
        D4(1,0,0,0)]
        sage: try_default_rule(A6,B3)
        [B3(1,0,0), B3(0,1,0), B3(0,0,2), B3(0,0,2), B3(0,1,0), B3(1,0,0)]

    If a default rule is not known, you may cue Sage as to what the
    Lie group embedding is by supplying a rule from the list of
    predefined rules. We will treat these next.

    .. RUBRIC:: Levi Type

    These can be read off from the Dynkin diagram. If
    removing a node from the Dynkin diagram produces another Dynkin
    diagram, there is a branching rule. A Levi subgroup may
    or may not be maximal. If it is maximal, there may or may not
    be a built-in branching rule for but you may obtain the
    Levi branching rule by first branching to a suitable
    maximal subgroup. For these rules use the option ``rule="levi"``:

    .. MATH::

        \begin{aligned}
        A_r & \Rightarrow A_{r-1}
        \\ B_r & \Rightarrow A_{r-1}
        \\ B_r & \Rightarrow B_{r-1}
        \\ C_r & \Rightarrow A_{r-1}
        \\ C_r & \Rightarrow C_{r-1}
        \\ D_r & \Rightarrow A_{r-1}
        \\ D_r & \Rightarrow D_{r-1}
        \\ E_r & \Rightarrow A_{r-1} \quad r = 7,8
        \\ E_r & \Rightarrow D_{r-1} \quad r = 6,7,8
        \\ E_r & \Rightarrow E_{r-1}
        \\ F_4 & \Rightarrow B_3
        \\ F_4 & \Rightarrow C_3
        \\ G_2 & \Rightarrow A_1 \text{(short root)}
        \end{aligned}

    Not all Levi subgroups are maximal subgroups. If the Levi is not
    maximal there may or may not be a preprogrammed ``rule="levi"`` for
    it. If there is not, the branching rule may still be obtained by going
    through an intermediate subgroup that is maximal using rule="extended".
    Thus the other Levi branching rule from `G_2 \Rightarrow A_1` corresponding to the
    long root is available by first branching `G_2 \Rightarrow A_2` then `A_2 \Rightarrow A_1`.
    Similarly the branching rules to the Levi subgroup:

    .. MATH::

        E_r \Rightarrow A_{r-1} \quad r = 6,7,8

    may be obtained by first branching `E_6 \Rightarrow A_5 \times A_1`, `E_7 \Rightarrow A_7`
    or `E_8 \Rightarrow A_8`.

    EXAMPLES::

        sage: A1 = WeylCharacterRing("A1")
        sage: A2 = WeylCharacterRing("A2")
        sage: A3 = WeylCharacterRing("A3")
        sage: A4 = WeylCharacterRing("A4")
        sage: A5 = WeylCharacterRing("A5")
        sage: B2 = WeylCharacterRing("B2")
        sage: B3 = WeylCharacterRing("B3")
        sage: B4 = WeylCharacterRing("B4")
        sage: C2 = WeylCharacterRing("C2")
        sage: C3 = WeylCharacterRing("C3")
        sage: D3 = WeylCharacterRing("D3")
        sage: D4 = WeylCharacterRing("D4")
        sage: G2 = WeylCharacterRing("G2")
        sage: F4 = WeylCharacterRing("F4",style="coroots")
        sage: E6=WeylCharacterRing("E6",style="coroots")
        sage: E7=WeylCharacterRing("E7",style="coroots")
        sage: D5=WeylCharacterRing("D5",style="coroots")
        sage: D6=WeylCharacterRing("D6",style="coroots")
        sage: [B3(w).branch(A2,rule="levi") for w in B3.fundamental_weights()]
        [A2(0,0,0) + A2(1,0,0) + A2(0,0,-1),
        A2(0,0,0) + A2(1,0,0) + A2(1,1,0) + A2(1,0,-1) + A2(0,-1,-1) + A2(0,0,-1),
        A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]

    The last example must be understood as follows. The representation
    of `B_3` being branched is spin, which is not a representation of
    `SO(7)` but of its double cover `\mathrm{spin}(7)`. The group `A_2` is
    really GL(3) and the double cover of `SO(7)` induces a cover of `GL(3)`
    that is trivial over `SL(3)` but not over the center of `GL(3)`. The weight
    lattice for this `GL(3)` consists of triples `(a,b,c)` of half integers
    such that `a - b` and `b - c` are in `\ZZ`, and this is reflected in the
    last decomposition.

    ::

        sage: [C3(w).branch(A2,rule="levi") for w in C3.fundamental_weights()]
        [A2(1,0,0) + A2(0,0,-1),
        A2(1,1,0) + A2(1,0,-1) + A2(0,-1,-1),
        A2(-1,-1,-1) + A2(1,-1,-1) + A2(1,1,-1) + A2(1,1,1)]
        sage: [D4(w).branch(A3,rule="levi") for w in D4.fundamental_weights()]
        [A3(1,0,0,0) + A3(0,0,0,-1),
        A3(0,0,0,0) + A3(1,1,0,0) + A3(1,0,0,-1) + A3(0,0,-1,-1),
        A3(1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,-1/2),
        A3(-1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,1/2)]
        sage: [B3(w).branch(B2,rule="levi") for w in B3.fundamental_weights()]
        [2*B2(0,0) + B2(1,0), B2(0,0) + 2*B2(1,0) + B2(1,1), 2*B2(1/2,1/2)]
        sage: C3 = WeylCharacterRing(['C',3])
        sage: [C3(w).branch(C2,rule="levi") for w in C3.fundamental_weights()]
        [2*C2(0,0) + C2(1,0),
         C2(0,0) + 2*C2(1,0) + C2(1,1),
         C2(1,0) + 2*C2(1,1)]
        sage: [D5(w).branch(D4,rule="levi") for w in D5.fundamental_weights()]
        [2*D4(0,0,0,0) + D4(1,0,0,0),
         D4(0,0,0,0) + 2*D4(1,0,0,0) + D4(1,1,0,0),
         D4(1,0,0,0) + 2*D4(1,1,0,0) + D4(1,1,1,0),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2)]
        sage: G2(1,0,-1).branch(A1,rule="levi")
        A1(1,0) + A1(1,-1) + A1(0,-1)
        sage: E6=WeylCharacterRing("E6",style="coroots")
        sage: D5=WeylCharacterRing("D5",style="coroots")
        sage: fw = E6.fundamental_weights()
        sage: [E6(fw[i]).branch(D5,rule="levi") for i in [1,2,6]]
        [D5(0,0,0,0,0) + D5(0,0,0,0,1) + D5(1,0,0,0,0),
         D5(0,0,0,0,0) + D5(0,0,0,1,0) + D5(0,0,0,0,1) + D5(0,1,0,0,0),
         D5(0,0,0,0,0) + D5(0,0,0,1,0) + D5(1,0,0,0,0)]
         sage: E7=WeylCharacterRing("E7",style="coroots")
         sage: A3xA3xA1=WeylCharacterRing("A3xA3xA1",style="coroots")
         sage: E7(1,0,0,0,0,0,0).branch(A3xA3xA1,rule="extended") # long time (0.7s)
         A3xA3xA1(0,0,1,0,0,1,1) + A3xA3xA1(0,1,0,0,1,0,0) + A3xA3xA1(1,0,0,1,0,0,1) +
          A3xA3xA1(1,0,1,0,0,0,0) + A3xA3xA1(0,0,0,1,0,1,0) + A3xA3xA1(0,0,0,0,0,0,2)
        sage: fw = E7.fundamental_weights()
        sage: [E7(fw[i]).branch(D6,rule="levi") for i in [1,2,7]] # long time (0.3s)
        [3*D6(0,0,0,0,0,0) + 2*D6(0,0,0,0,1,0) + D6(0,1,0,0,0,0),
         3*D6(0,0,0,0,0,1) + 2*D6(1,0,0,0,0,0) + 2*D6(0,0,1,0,0,0) + D6(1,0,0,0,1,0),
         D6(0,0,0,0,0,1) + 2*D6(1,0,0,0,0,0)]
        sage: D7=WeylCharacterRing("D7",style="coroots")
        sage: E8=WeylCharacterRing("E8",style="coroots")
        sage: D7=WeylCharacterRing("D7",style="coroots")
        sage: E8(1,0,0,0,0,0,0,0).branch(D7,rule="levi") # long time (7s)
         3*D7(0,0,0,0,0,0,0) + 2*D7(0,0,0,0,0,1,0) + 2*D7(0,0,0,0,0,0,1) + 2*D7(1,0,0,0,0,0,0)
         + D7(0,1,0,0,0,0,0) + 2*D7(0,0,1,0,0,0,0) + D7(0,0,0,1,0,0,0) + D7(1,0,0,0,0,1,0) + D7(1,0,0,0,0,0,1) + D7(2,0,0,0,0,0,0)
        sage: E8(0,0,0,0,0,0,0,1).branch(D7,rule="levi") # long time (0.6s)
         D7(0,0,0,0,0,0,0) + D7(0,0,0,0,0,1,0) + D7(0,0,0,0,0,0,1) + 2*D7(1,0,0,0,0,0,0) + D7(0,1,0,0,0,0,0)
        sage: [F4(fw).branch(B3,rule="levi") for fw in F4.fundamental_weights()] # long time (1s)
         [B3(0,0,0) + 2*B3(1/2,1/2,1/2) + 2*B3(1,0,0) + B3(1,1,0),
         B3(0,0,0) + 6*B3(1/2,1/2,1/2) + 5*B3(1,0,0) + 7*B3(1,1,0) + 3*B3(1,1,1)
         + 6*B3(3/2,1/2,1/2) + 2*B3(3/2,3/2,1/2) + B3(2,0,0) + 2*B3(2,1,0) + B3(2,1,1),
         3*B3(0,0,0) + 6*B3(1/2,1/2,1/2) + 4*B3(1,0,0) + 3*B3(1,1,0) + B3(1,1,1) + 2*B3(3/2,1/2,1/2),
         3*B3(0,0,0) + 2*B3(1/2,1/2,1/2) + B3(1,0,0)]
        sage: [F4(fw).branch(C3,rule="levi") for fw in F4.fundamental_weights()] # long time (1s)
         [3*C3(0,0,0) + 2*C3(1,1,1) + C3(2,0,0),
         3*C3(0,0,0) + 6*C3(1,1,1) + 4*C3(2,0,0) + 2*C3(2,1,0) + 3*C3(2,2,0) + C3(2,2,2) + C3(3,1,0) + 2*C3(3,1,1),
         2*C3(1,0,0) + 3*C3(1,1,0) + C3(2,0,0) + 2*C3(2,1,0) + C3(2,1,1),
         2*C3(1,0,0) + C3(1,1,0)]
        sage: A1xA1 = WeylCharacterRing("A1xA1")
        sage: [A3(hwv).branch(A1xA1,rule="levi") for hwv in A3.fundamental_weights()]
        [A1xA1(1,0,0,0) + A1xA1(0,0,1,0),
        A1xA1(1,1,0,0) + A1xA1(1,0,1,0) + A1xA1(0,0,1,1),
        A1xA1(1,1,1,0) + A1xA1(1,0,1,1)]
        sage: A1xB1=WeylCharacterRing("A1xB1",style="coroots")
        sage: [B3(x).branch(A1xB1,rule="levi") for x in B3.fundamental_weights()]
        [2*A1xB1(1,0) + A1xB1(0,2),
        3*A1xB1(0,0) + 2*A1xB1(1,2) + A1xB1(2,0) + A1xB1(0,2),
        A1xB1(1,1) + 2*A1xB1(0,1)]

    .. RUBRIC:: Automorphic Type

    If the Dynkin diagram has a symmetry, then there
    is an automorphism that is a special case of a branching rule.
    There is also an exotic "triality" automorphism of `D_4` having order 3.
    Use ``rule="automorphic"`` (or for `D_4` ``rule="triality"``):

    .. MATH::

        \begin{aligned}
        A_r & \Rightarrow A_r
        \\ D_r & \Rightarrow D_r
        \\ E_6 & \Rightarrow E_6
        \end{aligned}

    EXAMPLES::

        sage: [A3(chi).branch(A3,rule="automorphic") for chi in A3.fundamental_weights()]
        [A3(0,0,0,-1), A3(0,0,-1,-1), A3(0,-1,-1,-1)]
        sage: [D4(chi).branch(D4,rule="automorphic") for chi in D4.fundamental_weights()]
        [D4(1,0,0,0), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1/2,1/2,1/2,-1/2)]

    Here is an example with `D_4` triality::

        sage: [D4(chi).branch(D4,rule="triality") for chi in D4.fundamental_weights()]
        [D4(1/2,1/2,1/2,-1/2), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1,0,0,0)]

    .. RUBRIC:: Symmetric Type

    Related to the automorphic type, when `G` admits
    an outer automorphism (usually of degree 2) we may consider
    the branching rule to the isotropy subgroup `H`. Outer
    automorphisms correspond to symmetries of the Dynkin diagram.
    For such isotropy subgroups use ``rule="symmetric"``.
    We may thus obtain the following branching rules.

    .. MATH::

        A_{2r} & \Rightarrow B_r
        \\ A_{2r-1} & \Rightarrow C_r
        \\ A_{2r-1} & \Rightarrow D_r
        \\ D_r & \Rightarrow B_{r-1}
        \\ E_6 & \Rightarrow F_4
        \\ E_6 & \Rightarrow C_4
        \\ D_4 & \Rightarrow G_2

    The last branching rule, `D_4 \Rightarrow G_2` is not to a maximal subgroup
    since `D_4 \Rightarrow B_3 \Rightarrow G_2`, but it is included for convenience.

    In some cases, two outer automorphisms that differ by an
    inner automorphism may have different fixed subgroups.
    Thus, while the Dynkin diagram of `E_6` has a single
    involutory automorphism, there are two involutions
    of the group (differing by an inner automorphism) with
    fixed subgroups `F_4` and `C_4`. Similarly
    `SL(2r)`, of Cartan type `A_{2r-1}`, has subgroups
    `SO(2r)` and `Sp(2r)`, both fixed subgroups of outer
    automorphisms that differ from each other by an inner
    automorphism.

    In many cases the Dynkin diagram of `H` can be obtained by
    folding the Dynkin diagram of `G`.

    EXAMPLES::

        sage: [w.branch(B2,rule="symmetric") for w in [A4(1,0,0,0,0),A4(1,1,0,0,0),A4(1,1,1,0,0),A4(2,0,0,0,0)]]
        [B2(1,0), B2(1,1), B2(1,1), B2(0,0) + B2(2,0)]
        sage: [A5(w).branch(C3,rule="symmetric") for w in A5.fundamental_weights()]
        [C3(1,0,0), C3(0,0,0) + C3(1,1,0), C3(1,0,0) + C3(1,1,1), C3(0,0,0) + C3(1,1,0), C3(1,0,0)]
        sage: [A5(w).branch(D3,rule="symmetric") for w in A5.fundamental_weights()]
        [D3(1,0,0), D3(1,1,0), D3(1,1,-1) + D3(1,1,1), D3(1,1,0), D3(1,0,0)]
        sage: [D4(x).branch(B3,rule="symmetric") for x in D4.fundamental_weights()]
        [B3(0,0,0) + B3(1,0,0), B3(1,0,0) + B3(1,1,0), B3(1/2,1/2,1/2), B3(1/2,1/2,1/2)]
        sage: [D4(x).branch(G2,rule="symmetric") for x in D4.fundamental_weights()]
        [G2(0,0,0) + G2(1,0,-1), 2*G2(1,0,-1) + G2(2,-1,-1), G2(0,0,0) + G2(1,0,-1), G2(0,0,0) + G2(1,0,-1)]
        sage: [E6(fw).branch(F4,rule="symmetric") for fw in E6.fundamental_weights()] # long time (4s)
        [F4(0,0,0,0) + F4(0,0,0,1),
         F4(0,0,0,1) + F4(1,0,0,0),
         F4(0,0,0,1) + F4(1,0,0,0) + F4(0,0,1,0),
         F4(1,0,0,0) + 2*F4(0,0,1,0) + F4(1,0,0,1) + F4(0,1,0,0),
         F4(0,0,0,1) + F4(1,0,0,0) + F4(0,0,1,0),
         F4(0,0,0,0) + F4(0,0,0,1)]
        sage: E6=WeylCharacterRing("E6",style="coroots")
        sage: C4=WeylCharacterRing("C4",style="coroots")
        sage: chi = E6(1,0,0,0,0,0); chi.degree()
        27
        sage: chi.branch(C4,rule="symmetric")
        C4(0,1,0,0)

    .. RUBRIC:: Extended Type

    If removing a node from the extended Dynkin diagram
    results in a Dynkin diagram, then there is a branching rule. Use
    ``rule="extended"`` for these. We will also use this classification
    for some rules that are not of this type, mainly involving type `B`,
    such as `D_6 \Rightarrow B_3 \times B_3`.

    Here is the extended Dynkin diagram for `D_6`::

            0       6
            O       O
            |       |
            |       |
        O---O---O---O---O
        1   2   3   4   6

    Removing the node 3 results in an embedding `D_3 \times D_3 \Rightarrow D_6`.
    This corresponds to the embedding `SO(6) \times SO(6) \Rightarrow SO(12)`, and
    is of extended type. On the other hand the embedding `SO(5) \times SO(7)
    \Rightarrow SO(12)` (e.g. `B_2 \times B_3 \Rightarrow D_6`) cannot be explained this way
    but for uniformity is implemented under ``rule="extended"``.

    The following rules are implemented as special cases
    of ``rule="extended"``:

    .. MATH::

        \begin{aligned}
        E_6 & \Rightarrow A_5 \times A_1, A_2 \times A_2 \times A_2
        \\ E_7 & \Rightarrow A_7, D_6 \times A_1, A_3 \times A_3 \times A_1
        \\ E_8 & \Rightarrow A_8, D_8, E_7 \times A_1, A_4 \times A_4,
        D_5 \times A_3, E_6 \times A_2
        \\ F_4 & \Rightarrow B_4, C_3 \times A_1, A_2 \times A_2, A_3 \times A_1
        \\ G_2 & \Rightarrow A_1 \times A_1
        \end{aligned}

    Note that `E_8` has only a limited number of representations of
    reasonably low degree.

    EXAMPLES::

        sage: [B3(x).branch(D3,rule="extended") for x in B3.fundamental_weights()]
        [D3(0,0,0) + D3(1,0,0),
         D3(1,0,0) + D3(1,1,0),
         D3(1/2,1/2,-1/2) + D3(1/2,1/2,1/2)]
        sage: [G2(w).branch(A2, rule="extended") for w in G2.fundamental_weights()]
        [A2(0,0,0) + A2(1/3,1/3,-2/3) + A2(2/3,-1/3,-1/3),
         A2(1/3,1/3,-2/3) + A2(2/3,-1/3,-1/3) + A2(1,0,-1)]
        sage: [F4(fw).branch(B4,rule="extended") for fw in F4.fundamental_weights()] # long time (2s)
        [B4(1/2,1/2,1/2,1/2) + B4(1,1,0,0),
         B4(1,1,0,0) + B4(1,1,1,0) + B4(3/2,1/2,1/2,1/2) + B4(3/2,3/2,1/2,1/2) + B4(2,1,1,0),
         B4(1/2,1/2,1/2,1/2) + B4(1,0,0,0) + B4(1,1,0,0) + B4(1,1,1,0) + B4(3/2,1/2,1/2,1/2),
         B4(0,0,0,0) + B4(1/2,1/2,1/2,1/2) + B4(1,0,0,0)]

        sage: E6 = WeylCharacterRing("E6", style="coroots")
        sage: A2xA2xA2 = WeylCharacterRing("A2xA2xA2",style="coroots")
        sage: A5xA1 = WeylCharacterRing("A5xA1",style="coroots")
        sage: G2 = WeylCharacterRing("G2", style="coroots")
        sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
        sage: F4 = WeylCharacterRing("F4",style="coroots")
        sage: A3xA1 = WeylCharacterRing("A3xA1", style="coroots")
        sage: A2xA2 = WeylCharacterRing("A2xA2", style="coroots")
        sage: A1xC3 = WeylCharacterRing("A1xC3",style="coroots")
        sage: E6(1,0,0,0,0,0).branch(A5xA1,rule="extended") # (0.7s)
         A5xA1(0,0,0,1,0,0) + A5xA1(1,0,0,0,0,1)
        sage: E6(1,0,0,0,0,0).branch(A2xA2xA2, rule="extended") # (0.7s)
        A2xA2xA2(0,1,1,0,0,0) + A2xA2xA2(1,0,0,0,0,1) + A2xA2xA2(0,0,0,1,1,0)
        sage: E7 = WeylCharacterRing("E7",style="coroots")
        sage: A7 = WeylCharacterRing("A7",style="coroots")
        sage: E7(1,0,0,0,0,0,0).branch(A7,rule="extended")
         A7(0,0,0,1,0,0,0) + A7(1,0,0,0,0,0,1)
        sage: D6xA1 = WeylCharacterRing("D6xA1",style="coroots")
        sage: E7(1,0,0,0,0,0,0).branch(D6xA1,rule="extended")
         D6xA1(0,0,0,0,1,0,1) + D6xA1(0,1,0,0,0,0,0) + D6xA1(0,0,0,0,0,0,2)
        sage: A5xA2 = WeylCharacterRing("A5xA2",style="coroots")
        sage: E7(1,0,0,0,0,0,0).branch(A5xA2,rule="extended")
        A5xA2(0,0,0,1,0,1,0) + A5xA2(0,1,0,0,0,0,1) + A5xA2(1,0,0,0,1,0,0) + A5xA2(0,0,0,0,0,1,1)
        sage: E8 = WeylCharacterRing("E8",style="coroots")
        sage: D8 = WeylCharacterRing("D8",style="coroots")
        sage: A8 = WeylCharacterRing("A8",style="coroots")
        sage: E8(0,0,0,0,0,0,0,1).branch(D8,rule="extended") # long time (0.56s)
         D8(0,0,0,0,0,0,1,0) + D8(0,1,0,0,0,0,0,0)
        sage: E8(0,0,0,0,0,0,0,1).branch(A8,rule="extended") # long time (0.73s)
        A8(0,0,0,0,0,1,0,0) + A8(0,0,1,0,0,0,0,0) + A8(1,0,0,0,0,0,0,1)
        sage: F4(1,0,0,0).branch(A1xC3,rule="extended") # (0.05s)
         A1xC3(1,0,0,1) + A1xC3(2,0,0,0) + A1xC3(0,2,0,0)
        sage: G2(0,1).branch(A1xA1, rule="extended")
         A1xA1(2,0) + A1xA1(3,1) + A1xA1(0,2)
        sage: F4(0,0,0,1).branch(A2xA2, rule="extended") # (0.4s)
         A2xA2(0,1,0,1) + A2xA2(1,0,1,0) + A2xA2(0,0,1,1)
        sage: F4(0,0,0,1).branch(A3xA1,rule="extended") # (0.34s)
         A3xA1(0,0,0,0) + A3xA1(0,0,1,1) + A3xA1(0,1,0,0) + A3xA1(1,0,0,1) + A3xA1(0,0,0,2)
        sage: D4=WeylCharacterRing("D4",style="coroots")
        sage: D2xD2=WeylCharacterRing("D2xD2",style="coroots") # We get D4 => A1xA1xA1xA1 by remembering that A1xA1 = D2.
        sage: [D4(fw).branch(D2xD2, rule="extended") for fw in D4.fundamental_weights()]
         [D2xD2(1,1,0,0) + D2xD2(0,0,1,1),
         D2xD2(2,0,0,0) + D2xD2(0,2,0,0) + D2xD2(1,1,1,1) + D2xD2(0,0,2,0) + D2xD2(0,0,0,2),
         D2xD2(1,0,0,1) + D2xD2(0,1,1,0),
         D2xD2(1,0,1,0) + D2xD2(0,1,0,1)]

    .. RUBRIC:: Orthogonal Sum

    Using ``rule="orthogonal_sum"``, for `n = a + b + c + \cdots`,
    you can get any branching rule

    .. MATH::

        \begin{aligned}
        SO(n) & \Rightarrow SO(a) \times SO(b) \times SO(c) \times \cdots,
        \\ Sp(2n) & \Rightarrow Sp(2a) \times Sp(2b) \times Sp(2c) x \times \cdots,
        \end{aligned}

    where `O(a)` is type `D_r` for `a = 2r` or `B_r` for `a = 2r+1`
    and `Sp(2r)` is type `C_r`. In some cases these are also of
    extended type, as in the case `D_3 \times D_3 \Rightarrow D_6` discussed above.
    But in other cases, for example `B_3 \times B_3 \Rightarrow D_7`, they are not
    of extended type.

    .. RUBRIC:: Tensor

    There are branching rules:

    .. MATH::

        \begin{aligned}
        A_{rs-1} & \Rightarrow A_{r-1} \times A_{s-1},
        \\ B_{2rs+r+s} & \Rightarrow B_r \times B_s,
        \\ D_{2rs+s} & \Rightarrow B_r \times D_s,
        \\ D_{2rs} & \Rightarrow D_r \times D_s,
        \\ D_{2rs} & \Rightarrow C_r \times C_s,
        \\ C_{2rs+s} & \Rightarrow B_r \times C_s,
        \\ C_{2rs} & \Rightarrow C_r \times D_s.
        \end{aligned}

    corresponding to the tensor product homomorphism. For type
    `A`, the homomorphism is `GL(r) \times GL(s) \Rightarrow GL(rs)`. For the
    classical types, the relevant fact is that if `V, W` are
    orthogonal or symplectic spaces, that is, spaces endowed
    with symmetric or skew-symmetric bilinear forms, then `V \otimes W`
    is also an orthogonal space (if `V` and `W` are both
    orthogonal or both symplectic) or symplectic (if one of
    `V` and `W` is orthogonal and the other symplectic).

    The corresponding branching rules are obtained using ``rule="tensor"``.

    EXAMPLES::

        sage: A5=WeylCharacterRing("A5", style="coroots")
        sage: A2xA1=WeylCharacterRing("A2xA1", style="coroots")
        sage: [A5(hwv).branch(A2xA1, rule="tensor") for hwv in A5.fundamental_weights()]
        [A2xA1(1,0,1),
        A2xA1(0,1,2) + A2xA1(2,0,0),
        A2xA1(1,1,1) + A2xA1(0,0,3),
        A2xA1(1,0,2) + A2xA1(0,2,0),
        A2xA1(0,1,1)]
        sage: B4=WeylCharacterRing("B4",style="coroots")
        sage: B1xB1=WeylCharacterRing("B1xB1",style="coroots")
        sage: [B4(f).branch(B1xB1,rule="tensor") for f in B4.fundamental_weights()]
        [B1xB1(2,2),
        B1xB1(2,0) + B1xB1(2,4) + B1xB1(4,2) + B1xB1(0,2),
        B1xB1(2,0) + B1xB1(2,2) + B1xB1(2,4) + B1xB1(4,2) + B1xB1(4,4) + B1xB1(6,0) + B1xB1(0,2) + B1xB1(0,6),
        B1xB1(1,3) + B1xB1(3,1)]
        sage: D4=WeylCharacterRing("D4",style="coroots")
        sage: C2xC1=WeylCharacterRing("C2xC1",style="coroots")
        sage: [D4(f).branch(C2xC1,rule="tensor") for f in D4.fundamental_weights()]
        [C2xC1(1,0,1),
        C2xC1(0,1,2) + C2xC1(2,0,0) + C2xC1(0,0,2),
        C2xC1(1,0,1),
        C2xC1(0,1,0) + C2xC1(0,0,2)]
        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: B1xC1=WeylCharacterRing("B1xC1",style="coroots")
        sage: [C3(f).branch(B1xC1,rule="tensor") for f in C3.fundamental_weights()]
        [B1xC1(2,1), B1xC1(2,2) + B1xC1(4,0), B1xC1(4,1) + B1xC1(0,3)]

    .. RUBRIC:: Symmetric Power

    The `k`-th symmetric and exterior power homomorphisms map

    .. MATH::

        GL(n) \Rightarrow GL\left(\binom{n+k-1}{k}\right)
        \times GL\left(\binom{n}{k}\right).

    The corresponding branching rules are not implemented but a special
    case is. The `k`-th symmetric power homomorphism `SL(2) \Rightarrow GL(k+1)`
    has its image inside of `SO(2r+1)` if `k = 2r` and inside of `Sp(2r)` if
    `k = 2r - 1`. Hence there are branching rules:

    .. MATH::

        \begin{aligned}
        B_r & \Rightarrow A_1
        \\ C_r & \Rightarrow A_1
        \end{aligned}

    and these may be obtained using the rule "symmetric_power".

    EXAMPLES::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: B3=WeylCharacterRing("B3",style="coroots")
        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: [B3(fw).branch(A1,rule="symmetric_power") for fw in B3.fundamental_weights()]
        [A1(6), A1(2) + A1(6) + A1(10), A1(0) + A1(6)]
        sage: [C3(fw).branch(A1,rule="symmetric_power") for fw in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    .. RUBRIC:: Miscellaneous

    Use ``rule="miscellaneous"`` for the following embeddings of maximal subgroups,
    all involving exceptional groups.

    .. MATH::

        \begin{aligned}
        B_3 & \Rightarrow G_2,
        \\ E_6 & \Rightarrow G_2,
        \\ E_6 & \Rightarrow A_2,
        \\ F_4 & \Rightarrow G_2 \times A_1,
        \\ E_6 & \Rightarrow G_2 \times A_2,
        \\ E_7 & \Rightarrow G_2 \times C_3,
        \\ E_7 & \Rightarrow F_4 \times A_1,
        \\ E_7 & \Rightarrow A_1 \times A_1,
        \\ E_7 & \Rightarrow G_2 \times A_1,
        \\ E_8 & \Rightarrow G_2 \times F_4.
        \\ E_8 & \Rightarrow A2 \times A_1.
        \\ E_8 & \Rightarrow B2.
        \end{aligned}

    Except for those embeddings available by ``rule="extended"``, these
    are the only embeddings of these groups as maximal subgroups.
    There may be other embeddings besides these. For example,
    there are other more obvious embeddings of `A_2` and `G_2` into `E_6`.
    However the embeddings in this table are characterized as embeddings
    as maximal subgroups. Regarding the embeddings of `A_2` and `G_2` in
    `E_6`, the embeddings in question may be characterized by the condition that the
    27-dimensional representations of `E_6` restrict irreducibly to `A_2` or
    `G_2`. Since `G_2` has a subgroup isomorphic to `A_2`, it is worth
    mentioning that the composite branching rules::

        branching_rule("E6","G2","miscellaneous")*branching_rule("G2","A2","extended")
        branching_rule("E6","A2","miscellaneous")

    are distinct.

    These embeddings are described more completely (with references
    to the literature) in the thematic tutorial at:

    https://doc.sagemath.org/html/en/thematic_tutorials/lie.html


    EXAMPLES::

        sage: G2 = WeylCharacterRing("G2")
        sage: [fw1, fw2, fw3] = B3.fundamental_weights()
        sage: B3(fw1+fw3).branch(G2, rule="miscellaneous")
        G2(1,0,-1) + G2(2,-1,-1) + G2(2,0,-2)
        sage: E6 = WeylCharacterRing("E6",style="coroots")
        sage: G2 = WeylCharacterRing("G2",style="coroots")
        sage: E6(1,0,0,0,0,0).branch(G2,"miscellaneous")
        G2(2,0)
        sage: A2=WeylCharacterRing("A2",style="coroots")
        sage: E6(1,0,0,0,0,0).branch(A2,rule="miscellaneous")
        A2(2,2)
        sage: E6(0,1,0,0,0,0).branch(A2,rule="miscellaneous")
        A2(1,1) + A2(1,4) + A2(4,1)
        sage: E6(0,0,0,0,0,2).branch(G2,"miscellaneous") # long time (0.59s)
        G2(0,0) + G2(2,0) + G2(1,1) + G2(0,2) + G2(4,0)
        sage: F4=WeylCharacterRing("F4",style="coroots")
        sage: G2xA1=WeylCharacterRing("G2xA1",style="coroots")
        sage: F4(0,0,1,0).branch(G2xA1,rule="miscellaneous")
        G2xA1(1,0,0) + G2xA1(1,0,2) + G2xA1(1,0,4) + G2xA1(1,0,6) + G2xA1(0,1,4) + G2xA1(2,0,2) + G2xA1(0,0,2) + G2xA1(0,0,6)
        sage: E6 = WeylCharacterRing("E6",style="coroots")
        sage: A2xG2 = WeylCharacterRing("A2xG2",style="coroots")
        sage: E6(1,0,0,0,0,0).branch(A2xG2,rule="miscellaneous")
        A2xG2(0,1,1,0) + A2xG2(2,0,0,0)
        sage: E7=WeylCharacterRing("E7",style="coroots")
        sage: G2xC3=WeylCharacterRing("G2xC3",style="coroots")
        sage: E7(0,1,0,0,0,0,0).branch(G2xC3,rule="miscellaneous") # long time (1.84s)
        G2xC3(1,0,1,0,0) + G2xC3(1,0,1,1,0) + G2xC3(0,1,0,0,1) + G2xC3(2,0,1,0,0) + G2xC3(0,0,1,1,0)
        sage: F4xA1=WeylCharacterRing("F4xA1",style="coroots")
        sage: E7(0,0,0,0,0,0,1).branch(F4xA1,"miscellaneous")
        F4xA1(0,0,0,1,1) + F4xA1(0,0,0,0,3)
        sage: A1xA1=WeylCharacterRing("A1xA1",style="coroots")
        sage: E7(0,0,0,0,0,0,1).branch(A1xA1,rule="miscellaneous")
        A1xA1(2,5) + A1xA1(4,1) + A1xA1(6,3)
        sage: A2=WeylCharacterRing("A2",style="coroots")
        sage: E7(0,0,0,0,0,0,1).branch(A2,rule="miscellaneous")
        A2(0,6) + A2(6,0)
        sage: G2xA1=WeylCharacterRing("G2xA1",style="coroots")
        sage: E7(1,0,0,0,0,0,0).branch(G2xA1,rule="miscellaneous")
        G2xA1(1,0,4) + G2xA1(0,1,0) + G2xA1(2,0,2) + G2xA1(0,0,2)
        sage: E8 = WeylCharacterRing("E8",style="coroots")
        sage: G2xF4 = WeylCharacterRing("G2xF4",style="coroots")
        sage: E8(0,0,0,0,0,0,0,1).branch(G2xF4,rule="miscellaneous") # long time (0.76s)
        G2xF4(1,0,0,0,0,1) + G2xF4(0,1,0,0,0,0) + G2xF4(0,0,1,0,0,0)
        sage: E8=WeylCharacterRing("E8",style="coroots")
        sage: A1xA2=WeylCharacterRing("A1xA2",style="coroots")
        sage: E8(0,0,0,0,0,0,0,1).branch(A1xA2,rule="miscellaneous") # long time (0.76s)
        A1xA2(2,0,0) + A1xA2(2,2,2) + A1xA2(4,0,3) + A1xA2(4,3,0) + A1xA2(6,1,1) + A1xA2(0,1,1)
        sage: B2=WeylCharacterRing("B2",style="coroots")
        sage: E8(0,0,0,0,0,0,0,1).branch(B2,rule="miscellaneous") # long time (0.53s)
        B2(0,2) + B2(0,6) + B2(3,2)

    .. RUBRIC:: A1 maximal subgroups of exceptional groups

    There are seven cases where the exceptional group `G_2`, `F_4`,
    `E_7` or `E_8` contains a maximal subgroup of type `A_1`.
    These are tabulated in Theorem 1 of Testerman,
    The construction of the maximal A1's in the exceptional algebraic groups,
    Proc. Amer. Math. Soc. 116 (1992), no. 3, 635-644. The
    names of these branching rules are roman numerals referring
    to the seven cases of her Theorem 1. Use these branching
    rules as in the following examples.

    EXAMPLES::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: G2=WeylCharacterRing("G2",style="coroots")
        sage: F4=WeylCharacterRing("F4",style="coroots")
        sage: E7=WeylCharacterRing("E7",style="coroots")
        sage: E8=WeylCharacterRing("E8",style="coroots")
        sage: [G2(f).branch(A1,rule="i") for f in G2.fundamental_weights()]
        [A1(6), A1(2) + A1(10)]
        sage: F4(1,0,0,0).branch(A1,rule="ii")
        A1(2) + A1(10) + A1(14) + A1(22)
        sage: E7(0,0,0,0,0,0,1).branch(A1,rule="iii")
        A1(9) + A1(17) + A1(27)
        sage: E7(0,0,0,0,0,0,1).branch(A1,rule="iv")
        A1(5) + A1(11) + A1(15) + A1(21)
        sage: E8(0,0,0,0,0,0,0,1).branch(A1,rule="v") # long time (0.6s)
        A1(2) + A1(14) + A1(22) + A1(26) + A1(34) + A1(38) + A1(46) + A1(58)
        sage: E8(0,0,0,0,0,0,0,1).branch(A1,rule="vi") # long time (0.6s)
        A1(2) + A1(10) + A1(14) + A1(18) + A1(22) + A1(26) + A1(28) + A1(34) + A1(38) + A1(46)
        sage: E8(0,0,0,0,0,0,0,1).branch(A1,rule="vii") # long time (0.6s)
        A1(2) + A1(6) + A1(10) + A1(14) + A1(16) + A1(18) + 2*A1(22) + A1(26) + A1(28) + A1(34) + A1(38)

    .. RUBRIC:: Branching Rules From Plethysms

    Nearly all branching rules `G \Rightarrow H` where `G` is of type `A`, `B`, `C`
    or `D` are covered by the preceding rules. The function
    :func:`branching_rule_from_plethysm` covers the remaining cases.

    This is a general rule that includes any branching rule
    from types `A`, `B`, `C`, or `D` as a special case. Thus it could be
    used in place of the above rules and would give the same
    results. However it is most useful when branching from `G`
    to a maximal subgroup `H` such that
    `\mathrm{rank}(H) < \mathrm{rank}(G) - 1`.

    We consider a homomorphism `H \Rightarrow G` where `G` is one of
    `SL(r+1)`, `SO(2r+1)`, `Sp(2r)` or `SO(2r)`. The function
    :func:`branching_rule_from_plethysm` produces the corresponding
    branching rule. The main ingredient is the character
    `\chi` of the representation of `H` that is the homomorphism
    to `GL(r+1)`, `GL(2r+1)` or `GL(2r)`.

    This rule is so powerful that it contains the other
    rules implemented above as special cases. First let
    us consider the symmetric fifth power representation
    of `SL(2)`.

    ::

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: chi=A1([5])
        sage: chi.degree()
         6
        sage: chi.frobenius_schur_indicator()
        -1

    This confirms that the character has degree 6 and
    is symplectic, so it corresponds to a homomorphism
    `SL(2) \Rightarrow Sp(6)`, and there is a corresponding
    branching rule `C_3 \Rightarrow A_1`.

    ::

        sage: C3 = WeylCharacterRing("C3",style="coroots")
        sage: sym5rule = branching_rule_from_plethysm(chi,"C3")
        sage: [C3(hwv).branch(A1,rule=sym5rule) for hwv in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    This is identical to the results we would obtain using
    ``rule="symmetric_power"``. The next example gives a branching
    not available by other standard rules.

    ::

        sage: G2 = WeylCharacterRing("G2",style="coroots")
        sage: D7 = WeylCharacterRing("D7",style="coroots")
        sage: ad=G2(0,1); ad.degree(); ad.frobenius_schur_indicator()
         14
         1
        sage: spin = D7(0,0,0,0,0,1,0); spin.degree()
         64
        sage: spin.branch(G2, rule=branching_rule_from_plethysm(ad, "D7"))
         G2(1,1)

    We have confirmed that the adjoint representation of `G_2`
    gives a homomorphism into `SO(14)`, and that the pullback
    of the one of the two 64 dimensional spin representations
    to `SO(14)` is an irreducible representation of `G_2`.

    We do not actually have to create the character or
    its parent WeylCharacterRing to create the
    branching rule::

        sage: b = branching_rule("C7","C3(0,0,1)","plethysm"); b
        plethysm (along C3(0,0,1)) branching rule C7 => C3

    .. RUBRIC:: Isomorphic Type

    Although not usually referred to as a branching
    rule, the effects of the accidental isomorphisms may be handled
    using ``rule="isomorphic"``:

    .. MATH::

        \begin{aligned}
        B_2 & \Rightarrow C_2
        \\ C_2 & \Rightarrow B_2
        \\ A_3 & \Rightarrow D_3
        \\ D_3 & \Rightarrow A_3
        \\ D_2 & \Rightarrow A_1 \Rightarrow A_1
        \\ B_1 & \Rightarrow A_1
        \\ C_1 & \Rightarrow A_1
        \end{aligned}

    EXAMPLES::

        sage: B2 = WeylCharacterRing("B2")
        sage: C2 = WeylCharacterRing("C2")
        sage: [B2(x).branch(C2, rule="isomorphic") for x in B2.fundamental_weights()]
        [C2(1,1), C2(1,0)]
        sage: [C2(x).branch(B2, rule="isomorphic") for x in C2.fundamental_weights()]
        [B2(1/2,1/2), B2(1,0)]
        sage: D3 = WeylCharacterRing("D3")
        sage: A3 = WeylCharacterRing("A3")
        sage: [A3(x).branch(D3,rule="isomorphic") for x in A3.fundamental_weights()]
        [D3(1/2,1/2,1/2), D3(1,0,0), D3(1/2,1/2,-1/2)]
        sage: [D3(x).branch(A3,rule="isomorphic") for x in D3.fundamental_weights()]
        [A3(1/2,1/2,-1/2,-1/2), A3(1/4,1/4,1/4,-3/4), A3(3/4,-1/4,-1/4,-1/4)]

    Here `A_3(x,y,z,w)` can be understood as a representation of `SL(4)`.
    The weights `x,y,z,w` and `x+t,y+t,z+t,w+t` represent the same
    representation of `SL(4)` - though not of `GL(4)` - since
    `A_3(x+t,y+t,z+t,w+t)` is the same as `A_3(x,y,z,w)` tensored with
    `\mathrm{det}^t`. So as a representation of `SL(4)`,
    ``A3(1/4,1/4,1/4,-3/4)`` is the same as ``A3(1,1,1,0)``. The exterior
    square representation `SL(4) \Rightarrow GL(6)` admits an invariant symmetric
    bilinear form, so is a representation `SL(4) \Rightarrow SO(6)` that lifts to
    an isomorphism `SL(4) \Rightarrow \mathrm{Spin}(6)`. Conversely, there are two
    isomorphisms `SO(6) \Rightarrow SL(4)`, of which we've selected one.

    In cases like this you might prefer ``style="coroots"``::

        sage: A3 = WeylCharacterRing("A3",style="coroots")
        sage: D3 = WeylCharacterRing("D3",style="coroots")
        sage: [D3(fw) for fw in D3.fundamental_weights()]
        [D3(1,0,0), D3(0,1,0), D3(0,0,1)]
        sage: [D3(fw).branch(A3,rule="isomorphic") for fw in D3.fundamental_weights()]
        [A3(0,1,0), A3(0,0,1), A3(1,0,0)]
        sage: D2 = WeylCharacterRing("D2", style="coroots")
        sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
        sage: [D2(fw).branch(A1xA1,rule="isomorphic") for fw in D2.fundamental_weights()]
        [A1xA1(1,0), A1xA1(0,1)]

    .. RUBRIC:: Branching From a Reducible WeylCharacterRing

    If the Cartan Type of R is reducible, we may project a character onto
    any of the components, or any combination of components. The rule to
    project on the first component is specified by the string ``"proj1"``,
    the rule to project on the second component is ``"proj2". To
    project on the first and third components, use ``"proj13"`` and so on.

    EXAMPLES::

        sage: A2xG2=WeylCharacterRing("A2xG2",style="coroots")
        sage: A2=WeylCharacterRing("A2",style="coroots")
        sage: G2=WeylCharacterRing("G2",style="coroots")
        sage: A2xG2(1,0,1,0).branch(A2,rule="proj1")
        7*A2(1,0)
        sage: A2xG2(1,0,1,0).branch(G2,rule="proj2")
        3*G2(1,0)
        sage: A2xA2xG2=WeylCharacterRing("A2xA2xG2",style="coroots")
        sage: A2xA2xG2(0,1,1,1,0,1).branch(A2xG2,rule="proj13")
        8*A2xG2(0,1,0,1)

    A more general way of specifying a branching rule from a reducible type is
    to supply a *list* of rules, one *component rule* for each
    component type in the root system. In the following example, we branch the
    fundamental representations of `D_4` down to `A_1\times A_1\times A_1
    \times A_1` through the intermediate group `D_2\times D_2`. We use
    multiplicative notation to compose the branching rules. There is no need
    to construct the intermediate WeylCharacterRing with type `D_2\times D_2`.

    EXAMPLES::

        sage: D4 = WeylCharacterRing("D4",style="coroots")
        sage: A1xA1xA1xA1 = WeylCharacterRing("A1xA1xA1xA1",style="coroots")
        sage: b = branching_rule("D2","A1xA1","isomorphic")
        sage: br = branching_rule("D4","D2xD2","extended")*branching_rule("D2xD2","A1xA1xA1xA1",[b,b])
        sage: [D4(fw).branch(A1xA1xA1xA1,rule=br) for fw in D4.fundamental_weights()]
        [A1xA1xA1xA1(1,1,0,0) + A1xA1xA1xA1(0,0,1,1),
        A1xA1xA1xA1(1,1,1,1) + A1xA1xA1xA1(2,0,0,0) + A1xA1xA1xA1(0,2,0,0) + A1xA1xA1xA1(0,0,2,0) + A1xA1xA1xA1(0,0,0,2),
        A1xA1xA1xA1(1,0,0,1) + A1xA1xA1xA1(0,1,1,0),
        A1xA1xA1xA1(1,0,1,0) + A1xA1xA1xA1(0,1,0,1)]

    In the list of rules to be supplied in branching from a reducible root
    system, we may use two key words "omit" and "identity". The term "omit"
    means that we omit one factor, projecting onto the remaining factors.
    The term "identity" is supplied when the irreducible factor Cartan Types
    of both the target and the source are the same, and the component
    branching rule is to be the identity map. For example, we have
    projection maps from `A_3\times A_2` to `A_3` and `A_2`, and
    the corresponding branching may be accomplished as follows. In
    this example the same could be accomplished using ``rule="proj2"``.

    EXAMPLES::

        sage: A3xA2=WeylCharacterRing("A3xA2",style="coroots")
        sage: A3=WeylCharacterRing("A3",style="coroots")
        sage: chi = A3xA2(0,1,0,1,0)
        sage: chi.branch(A3,rule=["identity","omit"])
        3*A3(0,1,0)
        sage: A2=WeylCharacterRing("A2",style="coroots")
        sage: chi.branch(A2,rule=["omit","identity"])
        6*A2(1,0)

    Yet another way of branching from a reducible root system with
    repeated Cartan types is to embed along the diagonal. The
    branching rule is equivalent to the tensor product, as
    the example shows::

        sage: G2=WeylCharacterRing("G2",style="coroots")
        sage: G2xG2=WeylCharacterRing("G2xG2",style="coroots")
        sage: G2=WeylCharacterRing("G2",style="coroots")
        sage: G2xG2(1,0,0,1).branch(G2,rule="diagonal")
        G2(1,0) + G2(2,0) + G2(1,1)
        sage: G2xG2(1,0,0,1).branch(G2,rule="diagonal") == G2(1,0)*G2(0,1)
        True

    .. RUBRIC:: Writing Your Own (Branching) Rules

    Suppose you want to branch from a group `G` to a subgroup `H`.
    Arrange the embedding so that a Cartan subalgebra `U` of `H` is
    contained in a Cartan subalgebra `T` of `G`. There is thus
    a mapping from the weight spaces `\mathrm{Lie}(T)^* \Rightarrow \mathrm{Lie}(U)^*`.
    Two embeddings will produce identical branching rules if they
    differ by an element of the Weyl group of `H`.

    The *rule* is this map `\mathrm{Lie}(T)^*`, which is ``G.space()``, to
    `\mathrm{Lie}(U)^*`, which is ``H.space()``,
    which you may implement as a function. As an example, let
    us consider how to implement the branching rule `A_3 \Rightarrow C_2`.
    Here `H = C_2 = Sp(4)` embedded as a subgroup in `A_3 = GL(4)`. The
    Cartan subalgebra `U` consists of diagonal matrices with
    eigenvalues `u_1, u_2, -u_2, -u_1`. The ``C2.space()`` is the
    two dimensional vector spaces consisting of the linear
    functionals `u_1` and `u_2` on `U`. On the other hand `\mathrm{Lie}(T)` is
    `\RR^4`. A convenient way to see the restriction is to
    think of it as the adjoint of the map `(u_1, u_2) \mapsto
    (u_1,u_2, -u_2, -u_1)`,
    that is, `(x_0, x_1, x_2, x_3) \Rightarrow (x_0 - x_3, x_1 - x_2)`. Hence we may
    encode the rule as follows::

       def rule(x):
           return [x[0]-x[3],x[1]-x[2]]

    or simply::

        rule = lambda x: [x[0]-x[3],x[1]-x[2]]

    We may now make and use the branching rule as follows.

    EXAMPLES::

        sage: br = BranchingRule("A3", "C2", lambda x: [x[0]-x[3],x[1]-x[2]], "homemade"); br
        homemade branching rule A3 => C2
        sage: [A3,C2]=[WeylCharacterRing(x,style="coroots") for x in ["A3","C2"]]
        sage: A3(0,1,0).branch(C2,rule=br)
        C2(0,0) + C2(0,1)
    """
    if isinstance(rule, (str, list)):
        rule = branching_rule(R._cartan_type, S._cartan_type, rule)
    if hasattr(rule, "_S"):
        if rule._S != S.cartan_type():
            raise ValueError("rule has wrong target Cartan type")
    mdict = {}
    for k in chi.weight_multiplicities():
        # TODO: Could this use the new from_vector of ambient_space ?
        if S._style == "coroots":
            if S._cartan_type.is_atomic() and S._cartan_type[0] == 'E':
                if S._cartan_type[1] == 6:
                    h = S._space(rule(list(k.to_vector()))).coerce_to_e6()
                elif S._cartan_type[1] == 7:
                    h = S._space(rule(list(k.to_vector()))).coerce_to_e7()
            else:
                h = (S._space(rule(list(k.to_vector())))).coerce_to_sl()
        else:
            h = S._space(rule(list(k.to_vector())))
        chi_mdict = chi.weight_multiplicities()
        if h in mdict:
            mdict[h] += chi_mdict[k]
        else:
            mdict[h] = chi_mdict[k]
    return S.char_from_weights(mdict)


class BranchingRule(SageObject):
    """
    A class for branching rules.
    """
    def __init__(self, R, S, f, name="default", intermediate_types=[],
                 intermediate_names=[]):
        """
        INPUT:

        - ``R, S`` -- CartanTypes
        -  ``f`` -- a function from the weight lattice of R to the weight lattice of S.
        """
        self._R = CartanType(R)
        self._S = CartanType(S)
        self._f = f
        self._intermediate_types = intermediate_types
        if intermediate_names:
            self._intermediate_names = intermediate_names
        else:
            self._intermediate_names = [name]
        self._name = name

    def _repr_(self):
        """
        EXAMPLES::

            sage: branching_rule("E6","F4","symmetric")
            symmetric branching rule E6 => F4
            sage: b=branching_rule("F4","B3",rule="levi")*branching_rule("B3","G2",rule="miscellaneous"); b
            composite branching rule F4 => (levi) B3 => (miscellaneous) G2
        """
        R_repr = self._R._repr_(compact=True)
        S_repr = self._S._repr_(compact=True)
        if self._name == "composite":
            ret = "composite branching rule %s => " % R_repr
            for i in range(len(self._intermediate_types)):
                intt = self._intermediate_types[i]
                intn = self._intermediate_names[i]
                ret += "(%s) %s => " % (intn, intt._repr_(compact=True))
            ret += "(%s) %s" % (self._intermediate_names[-1], S_repr)
            return ret

        return "%s branching rule %s => %s" % (self._name, R_repr, S_repr)

    def __call__(self, x):
        """
        EXAMPLES::

            sage: b=branching_rule("A3","C2","symmetric")
            sage: b([2,1,0,0])
            [2, 1]
        """
        try:
            return self._f(x)
        except Exception:
            return self._f(x.to_vector())

    def __eq__(self, other):
        """
        Two branching rules with the same source and target Cartan types are
        considered equal if they are the same as mappings from the weight
        lattice of the larger group to the smaller. The last example shows
        that two rules may be different by this criterion yet describe the
        same branching, if they differ by conjugation by an element of the
        Weyl group.

        EXAMPLES::

            sage: b = branching_rule("E6","F4","symmetric")*branching_rule("F4","B3","levi")*branching_rule("B3","G2","miscellaneous"); b
            composite branching rule E6 => (symmetric) F4 => (levi) B3 => (miscellaneous) G2
            sage: c = branching_rule("E6", "G2xA2", "miscellaneous")*branching_rule("G2xA2", "G2", "proj1"); c
            composite branching rule E6 => (miscellaneous) G2xA2 => (proj1) G2
            sage: b == c
            True
            sage: d = branching_rule("A3","A2","levi")*branching_rule("A2","A1","levi"); d
            composite branching rule A3 => (levi) A2 => (levi) A1
            sage: e = branching_rule("A3","D3","isomorphic")*branching_rule("D3","B2","symmetric")*branching_rule("B2","A1","levi"); e
            composite branching rule A3 => (isomorphic) D3 => (symmetric) B2 => (levi) A1
            sage: d == e
            False
            sage: b1 = BranchingRule("A2","A2",lambda x: [x[2], x[1], x[0]], "long Weyl element conjugation")
            sage: b2 = BranchingRule("A2","A2",lambda x: x, "identity map")
            sage: b1 == b2
            False
            sage: A2 = WeylCharacterRing("A2",style="coroots")
            sage: [A2(f).branch(A2,rule=b1) == A2(f).branch(A2,rule=b2) for f in A2.fundamental_weights()]
            [True, True]
        """
        if not isinstance(other, BranchingRule):
            return False

        Rspace = RootSystem(self._R).ambient_space()
        Rspace_other = RootSystem(other._R).ambient_space()
        if Rspace != Rspace_other:
            return False

        Sspace = RootSystem(self._S).ambient_space()
        Sspace_other = RootSystem(other._S).ambient_space()
        if Sspace != Sspace_other:
            return False

        for v in Rspace.fundamental_weights():
            w = list(v.to_vector())
            if Sspace(self(w)) != Sspace(other(w)):
                return False
        return True

    def __ne__(self, other):
        """
        Test inequality

        EXAMPLES::

            sage: b1 = BranchingRule("A2","A2",lambda x: [x[2], x[1], x[0]], "long Weyl element conjugation")
            sage: b2 = BranchingRule("A2","A2",lambda x: x, "identity map")
            sage: b1 != b2
            True
        """
        return not(self == other)

    def __mul__(self, other):
        """
        EXAMPLES::

            sage: E6 = WeylCharacterRing("E6",style="coroots")
            sage: A5 = WeylCharacterRing("A5",style="coroots")
            sage: br = branching_rule("E6","A5xA1",rule="extended")*branching_rule("A5xA1","A5",rule="proj1"); br
            composite branching rule E6 => (extended) A5xA1 => (proj1) A5
            sage: E6(1,0,0,0,0,0).branch(A5,rule=br)
            A5(0,0,0,1,0) + 2*A5(1,0,0,0,0)
        """
        if self._S == other._R:
            intermediates = flatten([self._intermediate_types, self._S,
                                     other._intermediate_types])
            internames = flatten([self._intermediate_names,
                                  other._intermediate_names])

            def f(x):
                return other._f(self._f(x))
            return BranchingRule(self._R, other._S, f, "composite",
                                 intermediate_types=intermediates,
                                 intermediate_names=internames)
        else:
            raise ValueError("unable to define composite: source and target don't agree")

    def Rtype(self):
        """
        In a branching rule R => S, returns the Cartan Type of the ambient group R.

        EXAMPLES::

            sage: branching_rule("A3","A2","levi").Rtype()
            ['A', 3]
        """
        return self._R

    def Stype(self):
        """
        In a branching rule R => S, returns the Cartan Type of the subgroup S.

        EXAMPLES::

            sage: branching_rule("A3","A2","levi").Stype()
            ['A', 2]
        """
        return self._S

    def describe(self, verbose=False, debug=False, no_r=False):
        """
        Describes how extended roots restrict under self.

        EXAMPLES::

            sage: branching_rule("G2","A2","extended").describe()
            <BLANKLINE>
            3
            O=<=O---O
            1   2   0
            G2~
            <BLANKLINE>
            root restrictions G2 => A2:
            <BLANKLINE>
            O---O
            1   2
            A2
            <BLANKLINE>
            0 => 2
            2 => 1
            <BLANKLINE>
            For more detailed information use verbose=True

        In this example, `0` is the affine root, that is, the negative
        of the highest root, for `"G2"`. If `i => j` is printed, this
        means that the i-th simple (or affine) root of the ambient
        group restricts to the j-th simple root of the subgroup.
        For reference the Dynkin diagrams are also printed. The
        extended Dynkin diagram of the ambient group is printed if
        the affine root restricts to a simple root. More information
        is printed if the parameter `verbose` is true.
        """
        Rspace = RootSystem(self._R).ambient_space()
        Sspace = RootSystem(self._S).ambient_space()
        if self._R.is_compound():
            raise ValueError("Cannot describe branching rule from reducible type")
        if not no_r:
            print("\n%r" % self._R.affine().dynkin_diagram())
        if self._S.is_compound():
            for j in range(len(self._S.component_types())):
                ctype = self._S.component_types()[j]
                component_rule = self*branching_rule(self._S, ctype,
                                                     "proj%s" % (j + 1))
                print("\nprojection %d on %s " % (j + 1,
                                                  ctype._repr_(compact=True)),
                      component_rule.describe(verbose=verbose, no_r=True))
            if not verbose:
                print("\nfor more detailed information use verbose=True")
        else:
            print("root restrictions %s => %s:" % (self._R._repr_(compact=True),
                                                   self._S._repr_(compact=True)))
            print("\n%r\n" % self._S.dynkin_diagram())
            for j in self._R.affine().index_set():
                if j == 0:
                    r = -Rspace.highest_root()
                else:
                    r = Rspace.simple_roots()[j]
                resr = Sspace(self(list(r.to_vector())))
                if debug:
                    print("root %d: r = %s, b(r)=%s" % (j, r, resr))
                done = False
                if resr == Sspace.zero():
                    done = True
                    print("%s => (zero)" % j)
                else:
                    for s in Sspace.roots():
                        if s == resr:
                            for i in self._S.index_set():
                                if s == Sspace.simple_root(i):
                                    print("%s => %s" % (j, i))
                                    done = True
                                    break
                            if not done:
                                done = True
                                if verbose:
                                    print("%s => root %s" % (j, s))
                if not done:
                    done = True
                    if verbose:
                        print("%s => weight %s" % (j, resr))
            if verbose:
                print("\nfundamental weight restrictions %s => %s:" % (self._R._repr_(compact=True),self._S._repr_(compact=True)))
                for j in self._R.index_set():
                    resfw = Sspace(self(list(Rspace.fundamental_weight(j).to_vector())))
                    print("%d => %s" % (j,
                                        tuple([resfw.inner_product(a)
                                               for a in Sspace.simple_coroots()])))
            if not no_r and not verbose:
                print("\nFor more detailed information use verbose=True")

    def branch(self, chi, style=None):
        """
        INPUT:

        - ``chi`` -- A character of the WeylCharacterRing with Cartan type self.Rtype().

        Returns the branched character.

        EXAMPLES::

            sage: G2=WeylCharacterRing("G2",style="coroots")
            sage: chi=G2(1,1); chi.degree()
            64
            sage: b=G2.maximal_subgroup("A2"); b
            extended branching rule G2 => A2
            sage: b.branch(chi)
            A2(0,1) + A2(1,0) + A2(0,2) + 2*A2(1,1) + A2(2,0) + A2(1,2) + A2(2,1)
            sage: A2=WeylCharacterRing("A2",style="coroots"); A2
            The Weyl Character Ring of Type A2 with Integer Ring coefficients
            sage: chi.branch(A2,rule=b)
            A2(0,1) + A2(1,0) + A2(0,2) + 2*A2(1,1) + A2(2,0) + A2(1,2) + A2(2,1)

        """
        from sage.combinat.root_system.weyl_characters import WeylCharacterRing
        if style is None:
            style = chi.parent()._style
        S = WeylCharacterRing(self.Stype(), style=style)
        return chi.branch(S, rule=self)


def branching_rule(Rtype, Stype, rule="default"):
    """
    Creates a branching rule.

    INPUT:

    - ``R`` -- the Weyl Character Ring of `G`

    - ``S`` -- the Weyl Character Ring of `H`

    - ``rule`` -- a string describing the branching rule as a map from
      the weight space of `S` to the weight space of `R`.

    If the rule parameter is omitted, in some cases, a default rule is supplied. See
    :func:`~sage.combinat.root_system.branching_rules.branch_weyl_character`.

    EXAMPLES::

       sage: rule = branching_rule(CartanType("A3"),CartanType("C2"),"symmetric")
       sage: [rule(x) for x in WeylCharacterRing("A3").fundamental_weights()]
       [[1, 0], [1, 1], [1, 0]]
    """
    if rule == "plethysm":
        try:
            S = sage.combinat.root_system.weyl_characters.WeylCharacterRing(Stype.split("(")[0], style="coroots")
            chi = S(eval("("+Stype.split("(")[1]))
        except Exception:
            S = sage.combinat.root_system.weyl_characters.WeylCharacterRing(Stype.split(".")[0], style="coroots")
            chi = eval("S." + Stype.split(".")[1])
        return branching_rule_from_plethysm(chi, Rtype)
    Rtype = CartanType(Rtype)
    Stype = CartanType(Stype)
    r = Rtype.rank()
    s = Stype.rank()
    rdim = Rtype.root_system().ambient_space().dimension()
    sdim = Stype.root_system().ambient_space().dimension()
    if Rtype.is_compound():
        Rtypes = Rtype.component_types()
        if isinstance(rule, str):
            if rule[:4] == "proj":
                name = rule
                proj = [int(j)-1 for j in rule[4:]]
                rule = []
                for j in range(len(Rtypes)):
                    if j in proj:
                        rule.append("identity")
                    else:
                        rule.append("omit")
            elif rule == "diagonal":
                if not Stype.is_compound():
                    k = len(Rtypes)
                    n = RootSystem(Stype).ambient_space().dimension()
                    return BranchingRule(Rtype, Stype, lambda x: [sum(x[i+n*j] for j in range(k)) for i in range(n)], "diagonal")
                raise ValueError("invalid Cartan types for diagonal branching rule")
            else:
                raise ValueError("Rule not found")
        else:
            name = repr(rule)
        rules = []
        stor = []
        for i in range(len(Rtypes)):
            l = rule[i]
            if l != "omit":
                if l == "identity":
                    rules.append(BranchingRule(Rtypes[i], Rtypes[i], lambda x: x, "identity"))
                else:
                    rules.append(l)
                stor.append(i)
        shifts = Rtype._shifts
        Stypes = [CartanType(ru._S) for ru in rules]
        ntypes = len(Stypes)
        if Stype.is_compound():
            def br(x):
                yl = []
                for i in range(ntypes):
                    yl.append(rules[i](x[shifts[stor[i]]:shifts[stor[i]+1]]))
                return flatten(yl)
        else:
            j = stor[0]
            rulej = rules[0]

            def br(x):
                return rulej(x[shifts[j]:shifts[j+1]])
        return BranchingRule(Rtype, Stype, br, name)
    if Stype.is_compound():
        stypes = Stype.component_types()
    if rule == "default":
        if not Rtype.is_compound():
            if Stype.is_compound() and s == r-1:
                try:
                    return branching_rule(Rtype, Stype, rule="levi")
                except Exception:
                    pass
            if Rtype[0] == "A":
                if Stype[0] == "B" and r == 2*s:
                    return branching_rule(Rtype, Stype, rule="symmetric")
                elif Stype[0] == "C" and r == 2*s-1:
                    return branching_rule(Rtype, Stype, rule="symmetric")
                elif Stype[0] == "D" and r == 2*s-1:
                    return branching_rule(Rtype, Stype, rule="symmetric")
            elif Rtype[0] == "B" and Stype[0] == "D" and r == s:
                return branching_rule(Rtype, Stype, rule="extended")
            elif Rtype[0] == "D" and Stype[0] == "B" and r == s+1:
                return branching_rule(Rtype, Stype, rule="symmetric")

            if s == r-1:
                try:
                    return branching_rule(Rtype, Stype, rule="levi")
                except Exception:
                    pass
        raise ValueError("No default rule found (you must specify the rule)")
    elif rule == "identity":
        if Rtype is not Stype:
            raise ValueError("Cartan types must match for identity rule")
        return BranchingRule(Rtype, Stype, lambda x: x, "identity")
    elif rule == "levi":
        if not s == r-1:
            raise ValueError("Incompatible ranks")
        if Rtype[0] == 'A':
            if Stype.is_compound():
                if all(ct[0] == 'A' for ct in stypes) and rdim == sdim:
                    return BranchingRule(Rtype, Stype, lambda x: x, "levi")
                else:
                    raise ValueError("Rule not found")
            elif Stype[0] == 'A':
                return BranchingRule(Rtype, Stype, lambda x: list(x)[:r], "levi")
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] in ['B', 'C', 'D']:
            if Stype.is_atomic():
                if Stype[0] == 'A':
                    return BranchingRule(Rtype, Stype, lambda x: x, "levi")
                elif Stype[0] == Rtype[0]:
                    return BranchingRule(Rtype, Stype, lambda x: list(x)[1:], "levi")
            elif stypes[-1][0] == Rtype[0] and all(t[0] == 'A' for t in stypes[:-1]):
                return BranchingRule(Rtype, Stype, lambda x: x, "levi")
            else:
                raise ValueError("Rule not found")
        elif Rtype == CartanType("E6"):
            if Stype == CartanType("D5"):
                return BranchingRule(Rtype, Stype, lambda x: [-x[4],-x[3],-x[2],-x[1],-x[0]], "levi")
            elif Stype == CartanType("A5"):  # non-maximal levi
                return branching_rule("E6","A5xA1","extended")*branching_rule("A5xA1","A5","proj1")
            elif Stype.is_compound():
                if Stype[0] == CartanType("A4") and Stype[1] == CartanType("A1"):  # non-maximal levi
                    return branching_rule("E6","A5xA1","extended")*branching_rule("A5xA1","A4xA1",[branching_rule("A5","A4","levi"),"identity"])
                if Stype[0] == CartanType("A1") and Stype[1] == CartanType("A4"):  # non-maximal levi
                    return branching_rule("E6","A1xA5","extended")*branching_rule("A1xA5","A1xA4",["identity",branching_rule("A5","A4","levi")])
                elif Stype[0] == CartanType("A2") and Stype[1] == CartanType("A2") and Stype[2] == CartanType("A1"):  # non-maximal levi
                    return branching_rule("E6","A2xA2xA2","extended")*branching_rule("A2xA2xA2","A2xA2xA2",["identity","identity",branching_rule("A2","A2","automorphic")*branching_rule("A2","A1","levi")])
                elif Stype[0] == CartanType("A2") and Stype[1] == CartanType("A1") and Stype[2] == CartanType("A2"):  # non-maximal levi
                    raise ValueError("Not implemented: use A2xA2xA1 levi or A2xA2xA2 extended rule. (Non-maximal Levi.)")
                elif Stype[0] == CartanType("A1") and Stype[1] == CartanType("A2") and Stype[2] == CartanType("A2"):  # non-maximal levi
                    raise ValueError("Not implemented: use A2xA2xA1 levi or A2xA2xA2 extended rule. (Non-maximal Levi.)")
        elif Rtype == CartanType("E7"):
            if Stype == CartanType("D6"):
                return branching_rule("E7","D6xA1","extended")*branching_rule("D6xA1","D6","proj1")  # non-maximal levi
            if Stype == CartanType("E6"):
                return BranchingRule(Rtype, Stype, lambda x: [x[0], x[1], x[2], x[3], x[4], (x[5]+x[6]-x[7])/3, (2*x[5]+5*x[6]+x[7])/6, (-2*x[5]+x[6]+5*x[7])/6], "levi")
            elif Stype == CartanType("A6"):  # non-maximal levi
                return branching_rule("E7","A7","extended")*branching_rule("A7","A7","automorphic")*branching_rule("A7","A6","levi")
            if Stype.is_compound():
                if Stype[0] == CartanType("A5") and Stype[1] == CartanType("A1"):
                    return branching_rule("E7","A5xA2","extended")*branching_rule("A5xA2","A5xA1",["identity",branching_rule("A2","A2","automorphic")*branching_rule("A2","A1","levi")])
                elif Stype[0] == CartanType("A1") and Stype[1] == CartanType("A5"):
                    raise NotImplementedError("Not implemented: use A5xA1")
        elif Rtype == CartanType("E8"):
            if Stype == CartanType("D7"):
                return BranchingRule(Rtype, Stype, lambda x: [-x[6],-x[5],-x[4],-x[3],-x[2],-x[1],-x[0]], "levi")
            elif Stype == CartanType("E7"):
                return BranchingRule(Rtype, Stype, lambda x: [x[0],x[1],x[2],x[3],x[4],x[5],(x[6]-x[7])/2,(x[7]-x[6])/2], "levi")
            elif Stype == CartanType("A7"):
                return branching_rule("E8","A8","extended")*branching_rule("A8","A7","levi")
            raise NotImplementedError("Not implemented yet: branch first using extended rule to get non-maximal levis")
        elif Rtype == CartanType("F4"):
            if Stype == CartanType("B3"):
                return BranchingRule(Rtype, Stype, lambda x: x[1:], "levi")
            elif Stype == CartanType("C3"):
                return BranchingRule(Rtype, Stype, lambda x: [x[1]-x[0],x[2]+x[3],x[2]-x[3]], "levi")
            else:
                raise NotImplementedError("Not implemented yet")
        elif Rtype == CartanType("G2") and Stype == CartanType("A1"):
            return BranchingRule(Rtype, Stype, lambda x: list(x)[1:][:2], "levi")
        else:
            raise ValueError("Rule not found")
    elif rule == "automorphic":
        if not Rtype == Stype:
            raise ValueError("Cartan types must agree for automorphic branching rule")
        elif Rtype[0] == 'A':
            def rule(x):
                y = [-i for i in x]
                y.reverse()
                return y
            return BranchingRule(Rtype, Stype, rule, "automorphic")
        elif Rtype[0] == 'D':
            def rule(x):
                x[len(x) - 1] = -x[len(x) - 1]
                return x
            return BranchingRule(Rtype, Stype, rule, "automorphic")
        elif Rtype[0] == 'E' and r == 6:
            M = matrix(QQ,[(3, 3, 3, -3, 0, 0, 0, 0),
                           (3, 3, -3, 3, 0, 0, 0, 0),
                           (3, -3, 3, 3, 0, 0, 0, 0),
                           (-3, 3, 3, 3, 0, 0, 0, 0),
                           (0, 0, 0, 0, -3, -3, -3, 3),
                           (0, 0, 0, 0, -3, 5, -1, 1),
                           (0, 0, 0, 0, -3, -1, 5, 1),
                           (0, 0, 0, 0, 3, 1, 1, 5)])/6
            return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "automorphic")
        else:
            raise ValueError("No automorphism found")
    elif rule == "triality":
        if not Rtype == Stype:
            raise ValueError("Triality is an automorphic type (for D4 only)")
        elif not Rtype[0] == 'D' and r == 4:
            raise ValueError("Triality is for D4 only")
        else:
            return BranchingRule(Rtype, Stype, lambda x: [(x[0]+x[1]+x[2]+x[3])/2,(x[0]+x[1]-x[2]-x[3])/2,(x[0]-x[1]+x[2]-x[3])/2,(-x[0]+x[1]+x[2]-x[3])/2], "triality")
    elif rule == "symmetric":
        if Rtype[0] == 'A':
            if (Stype[0] == 'C' or Stype[0] == 'D' and r == 2*s-1) or (Stype[0] == 'B' and r == 2*s):
                return BranchingRule(Rtype, Stype, lambda x: [x[i]-x[r-i] for i in range(s)], "symmetric")
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'D' and Stype[0] == 'B' and s == r-1:
            return BranchingRule(Rtype, Stype, lambda x: x[:s], "symmetric")
        elif Rtype == CartanType("D4") and Stype == CartanType("G2"):
            return BranchingRule(Rtype, Stype, lambda x: [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]], "symmetric")
        elif Rtype == CartanType("E6") and Stype == CartanType("F4"):
            return BranchingRule(Rtype, Stype, lambda x: [(x[4]-3*x[5])/2,(x[0]+x[1]+x[2]+x[3])/2,(-x[0]-x[1]+x[2]+x[3])/2,(-x[0]+x[1]-x[2]+x[3])/2], "symmetric")
        elif Rtype == CartanType("E6") and Stype == CartanType("C4"):
            def f(x):
                [x0, x1, x2, x3, x4, x5] = x[:6]
                return [(x0+x1+x2+x3+x4-3*x5)/2,
                        (-x0-x1-x2-x3+x4-3*x5)/2,
                        -x0 + x3, -x1 + x2]
            return BranchingRule(Rtype, Stype, f, "symmetric")
        else:
            raise ValueError("Rule not found")
    elif rule == "extended" or rule == "orthogonal_sum":
        if rule == "extended" and not s == r:
            raise ValueError('Ranks should be equal for rule="extended"')
        if Stype.is_compound():
            if Rtype[0] in ['B','D'] and all(t[0] in ['B','D'] for t in stypes):
                if Rtype[0] == 'D':
                    rdeg = 2*r
                else:
                    rdeg = 2*r+1
                sdeg = 0
                for t in stypes:
                    if t[0] == 'D':
                        sdeg += 2*t[1]
                    else:
                        sdeg += 2*t[1]+1
                if rdeg == sdeg:
                    return BranchingRule(Rtype, Stype, lambda x: x[:s], "orthogonal_sum")
                else:
                    raise ValueError("Rule not found")
            elif Rtype[0] == 'C':
                if all(t[0] == Rtype[0] for t in stypes):
                    return BranchingRule(Rtype, Stype, lambda x: x, "orthogonal_sum")
            if rule == "orthogonal_sum":
                raise ValueError("Rule not found")
            elif Rtype[0] == 'E':
                if r == 6:
                    if stypes == [CartanType("A5"),CartanType("A1")]:
                        M = matrix(QQ,[(-3, -3, -3, -3, -3, -5, -5, 5),
                                       (-9, 3, 3, 3, 3, 1, 1, -1),
                                       (3, -9, 3, 3, 3, 1, 1, -1),
                                       (3, 3, -9, 3, 3, 1, 1, -1),
                                       (3, 3, 3, -9, 3, 1, 1, -1),
                                       (3, 3, 3, 3, -9, 9, -3, 3),
                                       (-3, -3, -3, -3, -3, -1, 11, 1),
                                       (3, 3, 3, 3, 3, 1, 1, 11)])/12
                        return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
                    if stypes == [CartanType("A1"),CartanType("A5")]:
                        M = matrix(QQ,[(-3, -3, -3, -3, -3, -1, 11, 1),
                                       (3, 3, 3, 3, 3, 1, 1, 11),
                                       (-3, -3, -3, -3, -3, -5, -5, 5),
                                       (-9, 3, 3, 3, 3, 1, 1, -1),
                                       (3, -9, 3, 3, 3, 1, 1, -1),
                                       (3, 3, -9, 3, 3, 1, 1, -1),
                                       (3, 3, 3, -9, 3, 1, 1, -1),
                                       (3, 3, 3, 3, -9, 9, -3, 3)])/12
                        return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
                    if stypes == [CartanType("A2"),CartanType("A2"),CartanType("A2")]:
                        M = matrix(QQ,[(0, 0, -2, -2, -2, -2, -2, 2),
                                       (-3, 3, 1, 1, 1, 1, 1, -1),
                                       (3, -3, 1, 1, 1, 1, 1, -1),
                                       (0, 0, -2, -2, 4, 0, 0, 0),
                                       (0, 0, -2, 4, -2, 0, 0, 0),
                                       (0, 0, 4, -2, -2, 0, 0, 0),
                                       (0, 0, -2, -2, -2, 2, 2, -2),
                                       (3, 3, 1, 1, 1, -1, -1, 1),
                                       (-3, -3, 1, 1, 1, -1, -1, 1)])/6
                        return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
                elif r == 7:
                    if stypes == [CartanType("D6"),CartanType("A1")]:
                        return BranchingRule(Rtype, Stype, lambda x: [x[5],x[4],x[3],x[2],x[1],x[0],x[6],x[7]], "extended")
                    elif stypes == [CartanType("A1"),CartanType("D6")]:
                        return BranchingRule(Rtype, Stype, lambda x: [x[6],x[7],x[5],x[4],x[3],x[2],x[1],x[0]], "extended")
                    elif stypes == [CartanType("A5"),CartanType("A2")]:
                        M = matrix(QQ,[(5, 1, 1, 1, 1, 1, 0, 0),
                                       (-1, -5, 1, 1, 1, 1, 0, 0),
                                       (-1, 1, -5, 1, 1, 1, 0, 0),
                                       (-1, 1, 1, -5, 1, 1, 0, 0),
                                       (-1, 1, 1, 1, -5, 1, 0, 0),
                                       (-1, 1, 1, 1, 1, -5, 0, 0),
                                       (1, -1, -1, -1, -1, -1, 0, -6),
                                       (1, -1, -1, -1, -1, -1, -6, 0),
                                       (-2, 2, 2, 2, 2, 2, -3, -3)])/6
                        return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
                    elif stypes == [CartanType("A3"),CartanType("A3"),CartanType("A1")]:
                        M = matrix(QQ, [(0, 0, -1, -1, -1, -1, 2, -2),
                                        (0, 0, -1, -1, -1, -1, -2, 2),
                                        (-2, 2, 1, 1, 1, 1, 0, 0),
                                        (2, -2, 1, 1, 1, 1, 0, 0),
                                        (0, 0, -1, -1, -1, 3, 0, 0),
                                        (0, 0, -1, -1, 3, -1, 0, 0),
                                        (0, 0, -1, 3, -1, -1, 0, 0),
                                        (0, 0, 3, -1, -1, -1, 0, 0),
                                        (2, 2, 0, 0, 0, 0, -2, -2),
                                        (-2, -2, 0, 0, 0, 0, -2, -2)])/4
                        return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
                elif r == 8:
                    if stypes == [CartanType("A4"),CartanType("A4")]:
                        M = matrix(QQ,[(0, 0, 0, -4, -4, -4, -4, 4),
                                       (-5, 5, 5, 1, 1, 1, 1, -1),
                                       (5, -5, 5, 1, 1, 1, 1, -1),
                                       (5, 5, -5, 1, 1, 1, 1, -1),
                                       (-5, -5, -5, 1, 1, 1, 1, -1),
                                       (0, 0, 0, -8, 2, 2, 2, -2),
                                       (0, 0, 0, 2, -8, 2, 2, -2),
                                       (0, 0, 0, 2, 2, -8, 2, -2),
                                       (0, 0, 0, 2, 2, 2, -8, -2),
                                       (0, 0, 0, 2, 2, 2, 2, 8)])/10
                        return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
                    elif len(stypes) == 3:
                        if 5 in stypes[0][i]:  # S is A5xA2xA1
                            raise NotImplementedError("Not maximal: first branch to A7xA1")
                    elif stypes == [CartanType("D5"), CartanType("A3")]:
                        raise NotImplementedError("Not maximal: first branch to D8 then D5xD3=D5xA3")
                    elif stypes == [CartanType("A3"), CartanType("D5")]:
                        raise NotImplementedError("Not maximal: first branch to D8 then D5xD3=D5xA3")
                    elif stypes == [CartanType("E6"), CartanType("A2")]:
                        def br(x):
                            return [x[0], x[1], x[2], x[3], x[4],
                                    (x[5]+x[6]-x[7])/3,(x[5]+x[6]-x[7])/3,
                                    (-x[5]-x[6]+x[7])/3,
                                    (-x[5]-x[6]-2*x[7])/3,
                                    (-x[5]+2*x[6]+x[7])/3,
                                    (2*x[5]-x[6]+x[7])/3]
                        return BranchingRule(Rtype, Stype, br, "extended")
                    elif stypes == [CartanType("E7"), CartanType("A1")]:
                        def br(x):
                            return [x[0], x[1], x[2], x[3], x[4], x[5],
                                    (x[6]-x[7])/2, (-x[6]+x[7])/2,
                                    (-x[6]-x[7])/2, (x[6]+x[7])/2]
                        return BranchingRule(Rtype, Stype, br, "extended")
                raise ValueError("Rule not found")
            elif Rtype[0] == 'F':
                if stypes == [CartanType("C3"), CartanType("A1")]:
                    return BranchingRule(Rtype, Stype, lambda x: [x[0]-x[1],x[2]+x[3],x[2]-x[3],(-x[0]-x[1])/2,(x[0]+x[1])/2], "extended")
                elif stypes == [CartanType("A1"), CartanType("C3")]:
                    return BranchingRule(Rtype, Stype, lambda x: [(-x[0]-x[1])/2,(x[0]+x[1])/2,x[0]-x[1],x[2]+x[3],x[2]-x[3]], "extended")
                elif stypes == [CartanType("A2"), CartanType("A2")]:
                    M = matrix(QQ,[(-2, -1, -1, 0), (1, 2, -1, 0), (1, -1, 2, 0), (1, -1, -1, 3), (1, -1, -1, -3), (-2, 2, 2, 0)])/3
                elif stypes == [CartanType("A3"), CartanType("A1")]:
                    M = matrix(QQ,[(-3, -1, -1, -1), (1, 3, -1, -1), (1, -1, 3, -1), (1, -1, -1, 3), (2, -2, -2, -2), (-2, 2, 2, 2)])/4
                elif stypes == [CartanType("A1"), CartanType("A3")]:
                    M = matrix(QQ,[(2, -2, -2, -2), (-2, 2, 2, 2), (-3, -1, -1, -1), (1, 3, -1, -1), (1, -1, 3, -1), (1, -1, -1, 3)])/4
                else:
                    raise ValueError("Rule not found")
                return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
            elif Rtype[0] == 'G':
                if stypes == [CartanType("A1"), CartanType("A1")]:
                    return BranchingRule(Rtype, Stype, lambda x: [(x[1]-x[2])/2,-(x[1]-x[2])/2, x[0]/2, -x[0]/2], "extended")
            raise ValueError("Rule not found")
        else:  # irreducible Stype
            if Rtype[0] == 'B' and Stype[0] == 'D':
                return BranchingRule(Rtype, Stype, lambda x: x, "extended")
            elif Rtype == CartanType("E7"):
                if Stype == CartanType("A7"):
                    M = matrix(QQ, [(-1, -1, -1, -1, -1, -1, 2, -2),
                                    (-1, -1, -1, -1, -1, -1, -2, 2),
                                    (-3, 1, 1, 1, 1, 1, 0, 0),
                                    (1, -3, 1, 1, 1, 1, 0, 0),
                                    (1, 1, -3, 1, 1, 1, 0, 0),
                                    (1, 1, 1, -3, 1, 1, 0, 0),
                                    (1, 1, 1, 1, -3, 1, 2, 2),
                                    (1, 1, 1, 1, 1, -3, 2, 2)])/4
                    return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
            elif Rtype == CartanType("E8"):
                if Stype == CartanType("D8"):
                    return BranchingRule(Rtype, Stype, lambda x: [-x[7],x[6],x[5],x[4],x[3],x[2],x[1],x[0]], "extended")
                elif Stype == CartanType("A8"):
                    M = matrix([(-2, -2, -2, -2, -2, -2, -2, 2),
                                (-5, 1, 1, 1, 1, 1, 1, -1),
                                (1, -5, 1, 1, 1, 1, 1, -1),
                                (1, 1, -5, 1, 1, 1, 1, -1),
                                (1, 1, 1, -5, 1, 1, 1, -1),
                                (1, 1, 1, 1, -5, 1, 1, -1),
                                (1, 1, 1, 1, 1, -5, 1, -1),
                                (1, 1, 1, 1, 1, 1, -5, -1),
                                (1, 1, 1, 1, 1, 1, 1, 5)])/6
                    return BranchingRule(Rtype, Stype, lambda x: tuple(M*vector(x)), "extended")
            elif Rtype == CartanType("F4") and Stype == CartanType("B4"):
                return BranchingRule(Rtype, Stype, lambda x: [-x[0], x[1], x[2], x[3]], "extended")
            elif Rtype == CartanType("G2") and Stype == CartanType("A2"):
                return BranchingRule(Rtype, Stype, lambda x: [(-x[1]+x[2])/3, (-x[0]+x[1])/3, (x[0]-x[2])/3], "extended")
            else:
                raise ValueError("Rule not found")
    elif rule == "isomorphic":
        if r != s:
            raise ValueError("Incompatible ranks")
        if Rtype == Stype:
            return BranchingRule(Rtype, Stype, lambda x: x, "isomorphic")
        elif Rtype == CartanType("B2") and Stype == CartanType("C2"):
            def rule(x):
                [x1, x2] = x
                return [x1 + x2, x1 - x2]
            return BranchingRule(Rtype, Stype, rule, "isomorphic")
        elif Rtype == CartanType("C2") and Stype == CartanType("B2"):
            def rule(x):
                [x1, x2] = x
                return [(x1 + x2) / 2, (x1 - x2) / 2]
            return BranchingRule(Rtype, Stype, rule, "isomorphic")
        elif Rtype == CartanType("B1") and Stype == CartanType("A1"):
            return BranchingRule(Rtype, Stype, lambda x: [x[0],-x[0]], "isomorphic")
        elif Rtype == CartanType("A1") and Stype == CartanType("B1"):
            return BranchingRule(Rtype, Stype, lambda x: [(x[0]-x[1])/2], "isomorphic")
        elif Rtype == CartanType("C1") and Stype == CartanType("A1"):
            return BranchingRule(Rtype, Stype, lambda x: [x[0]/2,-x[0]/2], "isomorphic")
        elif Rtype == CartanType("A1") and Stype == CartanType("C1"):
            return BranchingRule(Rtype, Stype, lambda x: [x[0]-x[1]], "isomorphic")
        elif Rtype == CartanType("A3") and Stype == CartanType("D3"):
            def rule(x):
                [x1, x2, x3, x4] = x
                return [(x1+x2-x3-x4)/2, (x1-x2+x3-x4)/2, (x1-x2-x3+x4)/2]
            return BranchingRule(Rtype, Stype, rule, "isomorphic")
        elif Rtype == CartanType("D3") and Stype == CartanType("A3"):
            def rule(x):
                [t1, t2, t3] = x
                return [(t1+t2+t3)/2, (t1-t2-t3)/2,
                        (-t1+t2-t3)/2, (-t1-t2+t3)/2]
            return BranchingRule(Rtype, Stype, rule, "isomorphic")
        elif Rtype == CartanType("D2") and Stype == CartanType("A1xA1"):
            def rule(x):
                [t1, t2] = x
                return [(t1-t2)/2, -(t1-t2)/2, (t1+t2)/2, -(t1+t2)/2]
            return BranchingRule(Rtype, Stype, rule, "isomorphic")
        else:
            raise ValueError("Rule not found")
    elif rule == "tensor" or rule == "tensor-debug":
        if not Stype.is_compound():
            raise ValueError("Tensor product requires more than one factor")
        if len(stypes) != 2:
            raise ValueError("Not implemented")
        if Rtype[0] == 'A':
            nr = Rtype[1]+1
        elif Rtype[0] == 'B':
            nr = 2*Rtype[1]+1
        elif Rtype[0] in ['C', 'D']:
            nr = 2*Rtype[1]
        else:
            raise ValueError("Rule not found")
        [s1, s2] = [stypes[i][1] for i in range(2)]
        ns = [s1, s2]
        for i in range(2):
            if stypes[i][0] == 'A':
                ns[i] = ns[i]+1
            if stypes[i][0] == 'B':
                ns[i] = 2*ns[i]+1
            if stypes[i][0] in ['C','D']:
                ns[i] = 2*ns[i]
        if nr != ns[0]*ns[1]:
            raise ValueError("Ranks don't agree with tensor product")
        if Rtype[0] == 'A':
            if all(t[0] == 'A' for t in stypes):
                def rule(x):
                    ret = [sum(x[i*ns[1]:(i+1)*ns[1]]) for i in range(ns[0])]
                    ret.extend([sum(x[ns[1]*j+i] for j in range(ns[0])) for i in range(ns[1])])
                    return ret
                return BranchingRule(Rtype, Stype, rule, "tensor")
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'B':
            if not all(t[0] == 'B' for t in stypes):
                raise ValueError("Rule not found")
        elif Rtype[0] == 'C':
            if stypes[0][0] in ['B','D'] and stypes[1][0] == 'C':
                pass
            elif stypes[1][0] in ['B','D'] and stypes[0][0] == 'C':
                pass
            else:
                raise ValueError("Rule not found")
        elif Rtype[0] == 'D':
            if stypes[0][0] in ['B','D'] and stypes[1][0] == 'D':
                pass
            elif stypes[1][0] == 'B' and stypes[0][0] == 'D':
                pass
            elif stypes[1][0] == 'C' and stypes[0][0] == 'C':
                pass
            else:
                raise ValueError("Rule not found")
        rows = []
        for i in range(s1):
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                nextrow[s1+j] = 1
                rows.append(nextrow)
        if stypes[1][0] == 'B':
            for i in range(s1):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                rows.append(nextrow)
        for i in range(s1):
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                nextrow[s1+j] = -1
                rows.append(nextrow)
        if stypes[0][0] == 'B':
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[s1+j] = 1
                rows.append(nextrow)
        mat = matrix(rows).transpose()
        if rule == "tensor-debug":
            print(mat)
        return BranchingRule(Rtype, Stype, lambda x: tuple(mat*vector(x)), "tensor")
    elif rule == "symmetric_power":
        if Stype[0] == 'A' and s == 1:
            if Rtype[0] == 'B':
                def rule(x):
                    a = sum((r-i)*x[i] for i in range(r))
                    return [a,-a]
                return BranchingRule(Rtype, Stype, rule, "symmetric_power")
            elif Rtype[0] == 'C':
                def rule(x):
                    a = sum((2*r-2*i-1)*x[i] for i in range(r))
                    return [a/2,-a/2]
                return BranchingRule(Rtype, Stype, rule, "symmetric_power")
    elif rule == "miscellaneous":
        if Rtype[0] == 'B' and Stype[0] == 'G' and r == 3:
            return BranchingRule(Rtype, Stype, lambda x: [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]], "miscellaneous")
        elif Rtype == CartanType("E6"):
            if Stype.is_compound():
                if stypes == [CartanType("A2"),CartanType("G2")]:
                    return BranchingRule(Rtype, Stype, lambda x: [-2*x[5],x[5]+x[4],x[5]-x[4],x[2]+x[3],x[1]-x[2],-x[1]-x[3]], "miscellaneous")
                elif stypes == [CartanType("G2"),CartanType("A2")]:
                    return BranchingRule(Rtype, Stype, lambda x: [x[2]+x[3],x[1]-x[2],-x[1]-x[3],-2*x[5],x[5]+x[4],x[5]-x[4]], "miscellaneous")
            else:
                if Stype == CartanType("G2"):
                    return BranchingRule(Rtype, Stype, lambda x: [x[2]+x[3]+x[4]-3*x[5], x[1]-2*x[2]-x[3], -x[1]+x[2]-x[4]+3*x[5]],"miscellaneous")
                if Stype == CartanType("A2"):
                    return BranchingRule(Rtype, Stype, lambda x: [x[2]+x[3]+x[4]-3*x[5], x[1]-2*x[2]-x[3], -x[1]+x[2]-x[4]+3*x[5]],"miscellaneous")
        elif Rtype == CartanType("E7"):
            if Stype.is_compound():
                if stypes == [CartanType("C3"), CartanType("G2")]:
                    return BranchingRule(Rtype, Stype, lambda x: [-2*x[6],x[4]+x[5],-x[4]+x[5],x[1]+x[3],x[2]-x[3],-x[1]-x[2]], "miscellaneous")
                elif stypes == [CartanType("G2"), CartanType("C3")]:
                    return BranchingRule(Rtype, Stype, lambda x: [x[1]+x[3],x[2]-x[3],-x[1]-x[2],-2*x[6],x[4]+x[5],-x[4]+x[5]], "miscellaneous")
                elif stypes == [CartanType("F4"), CartanType("A1")]:
                    def f(x):
                        [x0, x1, x2, x3, x4, x5, x6] = x[:7]
                        return [(x4-x5)/2-x6, (x0+x1+x2+x3)/2,
                                (-x0-x1+x2+x3)/2, (-x0+x1-x2+x3)/2,
                                x5-x6, x6-x5]
                    return BranchingRule(Rtype, Stype, f, "miscellaneous")
                elif stypes == [CartanType("A1"), CartanType("F4")]:
                    def f(x):
                        [x0, x1, x2, x3, x4, x5, x6] = x[:7]
                        return [x5-x6, x6-x5, (x4-x5)/2-x6,
                                (x0+x1+x2+x3)/2,
                                (-x0-x1+x2+x3)/2,
                                (-x0+x1-x2+x3)/2]
                    return BranchingRule(Rtype, Stype, f, "miscellaneous")
                elif stypes == [CartanType("A1"), CartanType("A1")]:
                    return BranchingRule(Rtype, Stype,
                                         lambda x: [x[1]+2*x[2]-2*x[3]-x[4]-2*x[6], -x[1]-2*x[2]+2*x[3]+x[4]+2*x[6],
                                                    (x[3]+x[4]+x[5]-3*x[6]),-(x[3]+x[4]+x[5]-3*x[6])], "miscellaneous")
                elif stypes == [CartanType("G2"), CartanType("A1")]:
                    def f(x):
                        return [(x[0]-x[1]+x[2]+3*x[3]+x[4]-x[5]+2*x[6])/2,
                                (-3*x[0]-x[1]-x[2]-x[3]+x[4]+x[5]-2*x[6])/2,
                                (2*x[0]+2*x[1]-2*x[3]-2*x[4])/2,
                                (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]-4*x[6])/2,
                                -(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]-4*x[6])/2]
                    return BranchingRule(Rtype, Stype, f, "miscellaneous")
                elif stypes == [CartanType("A1"), CartanType("G2")]:
                    def f(x):
                        return [(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]-4*x[6])/2,
                                -(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]-4*x[6])/2,
                                (x[0]-x[1]+x[2]+3*x[3]+x[4]-x[5]+2*x[6])/2,
                                (-3*x[0]-x[1]-x[2]-x[3]+x[4]+x[5]-2*x[6])/2,
                                (2*x[0]+2*x[1]-2*x[3]-2*x[4])/2]
                    return BranchingRule(Rtype, Stype, f, "miscellaneous")
            elif Stype == CartanType("A2"):
                return BranchingRule(Rtype, Stype, lambda x: (x[1]+x[2]+2*x[4]-4*x[6],-2*x[1]-x[2]+x[3]-2*x[4]+2*x[5],x[1]-x[3]-2*x[5]+4*x[6]), "miscellaneous")
        elif Rtype == CartanType("E8"):
            if Stype.is_compound():
                if stypes == [CartanType("F4"),CartanType("G2")]:
                    return BranchingRule(Rtype, Stype, lambda x: [x[7], x[6], x[5], x[4], x[1]+x[3], -x[3]+x[2], -x[1]-x[2]], "miscellaneous")
                elif stypes == [CartanType("G2"),CartanType("F4")]:
                    return BranchingRule(Rtype, Stype, lambda x: [x[1]+x[3], -x[3]+x[2], -x[1]-x[2], x[7], x[6], x[5], x[4]], "miscellaneous")
                elif stypes == [CartanType("A2"), CartanType("A1")]:
                    def f(x):
                        return [(x[0]-x[1]+x[2]+x[3]+3*x[4]+x[5]-x[6]-x[7])/2,
                                (-3*x[0]-x[1]-x[2]-x[3]-x[4]+x[5]+x[6]+x[7])/2,
                                (2*x[0]+2*x[1]-2*x[4]-2*x[5])/2,
                                (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+5*x[7])/2,
                                -(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+5*x[7])/2]
                    return BranchingRule("E8","A2xA1",f,"miscellaneous")
                elif stypes == [CartanType("A1"), CartanType("A2")]:
                    def f(x):
                        return [(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+5*x[7])/2,
                                -(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+5*x[7])/2,
                                (x[0]-x[1]+x[2]+x[3]+3*x[4]+x[5]-x[6]-x[7])/2,
                                (-3*x[0]-x[1]-x[2]-x[3]-x[4]+x[5]+x[6]+x[7])/2,
                                (2*x[0]+2*x[1]-2*x[4]-2*x[5])/2]
                    return BranchingRule("E8", "A1xA2", f, "miscellaneous")
            elif Stype == CartanType("B2"):
                return BranchingRule("E8", "B2", lambda x: [-x[0] + x[2] + x[5] + 3*x[7], 2*x[0] - x[2] + x[3] + x[4] + 2*x[6] + x[7]], "miscellaneous")
        elif Rtype[0] == 'F':
            if Stype.is_compound():
                if stypes == [CartanType("A1"), CartanType("G2")]:
                    return BranchingRule("F4", "A1xG2", lambda x: [2*x[0], -2*x[0], x[1]+x[2], -x[2]+x[3], -x[1]-x[3]], "miscellaneous")
                elif stypes == [CartanType("G2"), CartanType("A1")]:
                    return BranchingRule("F4","G2xA1", lambda x: [x[1]+x[2], -x[2]+x[3], -x[1]-x[3], 2*x[0], -2*x[0]], "miscellaneous")
        raise ValueError("Rule not found")
    elif rule in ["i", "ii", "iii", "iv", "v", "vi", "vii"]:
        if Stype != CartanType("A1"):
            raise ValueError("Wrong target Cartan Type for rule %s" % rule)
        if rule == "i" and Rtype == CartanType("G2"):
            return BranchingRule(Rtype, Stype, lambda x: [(5*x[0]-x[1]-4*x[2])/3,-(5*x[0]-x[1]-4*x[2])/3], "i")
        elif rule == "ii" and Rtype == CartanType("F4"):
            return BranchingRule(Rtype, Stype, lambda x: [8*x[0]+3*x[1]+2*x[2]+x[3],-(8*x[0]+3*x[1]+2*x[2]+x[3])], "ii")
        elif rule == "iii" and Rtype == CartanType("E7"):
            return BranchingRule(Rtype, Stype,
                                 lambda x: [x[1]+2*x[2]+3*x[3]+4*x[4]+5*x[5]-17*x[6],-(x[1]+2*x[2]+3*x[3]+4*x[4]+5*x[5]-17*x[6])], "iii")
        elif rule == "iv" and Rtype == CartanType("E7"):
            return BranchingRule(Rtype, Stype,
                                 lambda x: [x[1]+x[2]+2*x[3]+3*x[4]+4*x[5]-13*x[6],-(x[1]+x[2]+2*x[3]+3*x[4]+4*x[5]-13*x[6])], "iv")
        elif rule == "v" and Rtype == CartanType("E8"):
            return BranchingRule(Rtype, Stype,
                                 lambda x: [x[1]+2*x[2]+3*x[3]+4*x[4]+5*x[5]+6*x[6]+23*x[7],-(x[1]+2*x[2]+3*x[3]+4*x[4]+5*x[5]+6*x[6]+23*x[7])], "v")
        elif rule == "vi" and Rtype == CartanType("E8"):
            return BranchingRule(Rtype, Stype,
                                 lambda x: [x[1]+x[2]+2*x[3]+3*x[4]+4*x[5]+5*x[6]+18*x[7],-(x[1]+x[2]+2*x[3]+3*x[4]+4*x[5]+5*x[6]+18*x[7])], "vi")
        elif rule == "vii" and Rtype == CartanType("E8"):
            return BranchingRule(Rtype, Stype,
                                 lambda x: [x[1]+x[2]+2*x[3]+2*x[4]+3*x[5]+4*x[6]+15*x[7],-(x[1]+x[2]+2*x[3]+2*x[4]+3*x[5]+4*x[6]+15*x[7])], "vii")
        raise ValueError("Wrong source Cartan Type for rule %s" % rule)
    raise ValueError("Rule not found")

get_branching_rule = branching_rule


def branching_rule_from_plethysm(chi, cartan_type, return_matrix=False):
    r"""
    Create the branching rule of a plethysm.

    INPUT:

    - ``chi`` -- the character of an irreducible representation `\pi` of
      a group `G`
    - ``cartan_type`` -- a classical Cartan type (`A`,`B`,`C` or `D`).

    It is assumed that the image of the irreducible representation pi
    naturally has its image in the group `G`.

    Returns a branching rule for this plethysm.

    EXAMPLES:

    The adjoint representation `SL(3) \to GL(8)` factors
    through `SO(8)`. The branching rule in question will
    describe how representations of `SO(8)` composed with
    this homomorphism decompose into irreducible characters
    of `SL(3)`::

        sage: A2 = WeylCharacterRing("A2")
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: ad = A2.adjoint_representation(); ad
        A2(1,1)
        sage: ad.degree()
        8
        sage: ad.frobenius_schur_indicator()
        1

    This confirms that `ad` has degree 8 and is orthogonal,
    hence factors through `SO(8)` which is type `D_4`::

        sage: br = branching_rule_from_plethysm(ad,"D4")
        sage: D4 = WeylCharacterRing("D4")
        sage: [D4(f).branch(A2,rule = br) for f in D4.fundamental_weights()]
        [A2(1,1), A2(0,3) + A2(1,1) + A2(3,0), A2(1,1), A2(1,1)]
    """
    ct = CartanType(cartan_type)
    if ct[0] not in ["A", "B", "C", "D"]:
        raise ValueError("not implemented for type {}".format(ct[0]))
    if ct[0] == "A":
        ret = []
        ml = chi.weight_multiplicities()
        for v in ml:
            n = ml[v]
            ret.extend(n * [v.to_vector()])
        M = matrix(ret).transpose()
        if len(M.columns()) != ct[1] + 1:
            raise ValueError("representation has wrong degree for type {}".format(ct))
        return BranchingRule(ct, chi.parent().cartan_type(), lambda x: tuple(M*vector(x)), "plethysm (along %s)" % chi)
    if ct[0] in ["B", "D"]:
        if chi.frobenius_schur_indicator() != 1:
            raise ValueError("character is not orthogonal")
    if ct[0] == "C":
        if chi.frobenius_schur_indicator() != -1:
            raise ValueError("character is not symplectic")
    if ct[0] == "B":
        if is_even(chi.degree()):
            raise ValueError("degree is not odd")
    if ct[0] in ["C", "D"]:
        if is_odd(chi.degree()):
            raise ValueError("degree is not even")
    ret = []
    ml = chi.weight_multiplicities()
    for v in ml:
        n = ml[v]
        vec = v.to_vector()
        if all(x == 0 for x in vec):
            if ct[0] == "B":
                n = (n-1)/2
            else:
                n = n/2
        elif [x for x in vec if x != 0][0] < 0:
            continue
        ret.extend(n * [vec])
    M = matrix(ret).transpose()
    if len(M.columns()) != ct.root_system().ambient_space().dimension():
        raise ValueError("representation has wrong degree for type {}".format(ct))
    if return_matrix:
        return M
    else:
        return BranchingRule(ct, chi.parent().cartan_type(), lambda x: tuple(M*vector(x)), "plethysm (along %s)" % chi)


def maximal_subgroups(ct, mode="print_rules"):
    """
    Given a classical Cartan type (of rank less than or equal to 8)
    this prints the Cartan types of maximal subgroups, with a method
    of obtaining the branching rule. The string to the right of the
    colon in the output is a command to create a branching rule.

    INPUT:

    - ``ct`` -- a classical irreducible Cartan type

    Returns a list of maximal subgroups of ct.

    EXAMPLES::

        sage: from sage.combinat.root_system.branching_rules import maximal_subgroups
        sage: maximal_subgroups("D4")
        B3:branching_rule("D4","B3","symmetric")
        A2:branching_rule("D4","A2(1,1)","plethysm")
        A1xC2:branching_rule("D4","C1xC2","tensor")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])
        A1xA1xA1xA1:branching_rule("D4","D2xD2","orthogonal_sum")*branching_rule("D2xD2","A1xA1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])

    .. SEEALSO:: :meth:`~sage.combinat.root_system.weyl_characters.WeylCharacterRing.ParentMethods.maximal_subgroups`

    """

    if CartanType(ct) == CartanType("A2"):
        rul = ["""A1:branching_rule("A2","A1","levi")"""]
    elif CartanType(ct) == CartanType("A3"):
        rul = ["""A2:branching_rule("A3","A2","levi")""",
               """A1xA1:branching_rule("A3","A1xA1","tensor")""",
               """C2:branching_rule("A3","C2","symmetric")""",
               """A1xA1:branching_rule("A3","A1xA1","levi")"""]
    elif CartanType(ct) == CartanType("A4"):
        rul = ["""A3:branching_rule("A4","A3","levi")""",
               """B2:branching_rule("A4","B2","symmetric")""",
               """A1xA2:branching_rule("A4","A1xA2","levi")"""]
    elif CartanType(ct) == CartanType("A5"):
        rul = ["""A4:branching_rule("A5","A4","levi")""",
               """A3:branching_rule("A5","D3","symmetric")*branching_rule("D3","A3","isomorphic")""",
               """A3:branching_rule("A5","A3(0,1,0)","plethysm") # alternative""",
               """C3:branching_rule("A5","C3","symmetric")""",
               """A2:branching_rule("A5","A2(2,0)","plethysm")""",
               """A1xA2:branching_rule("A5","A1xA2","tensor")""",
               """A1xA3:branching_rule("A5","A1xA3","levi")""",
               """A2xA2:branching_rule("A5","A2xA2","levi")"""]
    elif CartanType(ct) == CartanType("A6"):
        rul = ["""A5:branching_rule("A6","A5","levi")""",
               """B3:branching_rule("A6","B3","symmetric")""",
               """A1xA4:branching_rule("A6","A1xA4","levi")""",
               """A2xA3:branching_rule("A6","A2xA3","levi")"""]
    elif CartanType(ct) == CartanType("A7"):
        rul = ["""A6:branching_rule("A7","A6","levi")""",
               """C4:branching_rule("A7","C4","symmetric")""",
               """D4:branching_rule("A7","D4","symmetric")""",
               """A1xA3:branching_rule("A7","A1xA3","tensor")""",
               """A1xA5:branching_rule("A7","A1xA5","levi")""",
               """A2xA4:branching_rule("A7","A2xA4","levi")""",
               """A3xA3:branching_rule("A7","A3xA3","levi")"""]
    elif CartanType(ct) == CartanType("A8"):
        rul = ["""A7:branching_rule("A8","A7","levi")""",
               """B4:branching_rule("A8","B4","symmetric")""",
               """A2xA2:branching_rule("A8","A2xA2","tensor")""",
               """A1xA6:branching_rule("A8","A1xA6","levi")""",
               """A2xA5:branching_rule("A8","A2xA5","levi")""",
               """A3xA4:branching_rule("A8","A3xA4","levi")"""]
    elif CartanType(ct) == CartanType("B3"):
        rul = ["""G2:branching_rule("B3","G2","miscellaneous")""",
               """A3:branching_rule("B3","D3","extended")*branching_rule("D3","A3","isomorphic")""",
               """A1xA1xA1:branching_rule("B3","D2xB1","orthogonal_sum")*branching_rule("D2xB1","A1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("B1","A1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("B4"):
        rul = ["""D4:branching_rule("B4","D4","extended")""",
               """A1:branching_rule("B4","A1","symmetric_power")""",
               """A1xA1:branching_rule("B4","B1xB1","tensor")*branching_rule("B1xB1","A1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("B1","A1","isomorphic")])""",
               """A1xA1xB2:branching_rule("B4","D2xB2","extended")*branching_rule("D2xB2","A1xA1xB2",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """A1xA3:branching_rule("B4","B1xD3","extended")*branching_rule("B1xD3","A1xA3",[branching_rule("B1","A1","isomorphic"),branching_rule("D3","A3","isomorphic")])"""]
    elif CartanType(ct) == CartanType("B5"):
        rul = ["""D5:branching_rule("B5","D5","extended")""",
               """A1:branching_rule("B5","A1","symmetric_power")""",
               """A1xA2xB3:branching_rule("B5","D2xB3","extended")*branching_rule("D2xB3","A1xA2xB3",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """A1xD4:branching_rule("B5","B1xD4","orthogonal_sum")*branching_rule("B1xD4","A1xD4",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """A3xB2:branching_rule("B5","D3xB2","orthogonal_sum")*branching_rule("D3xB2","A3xB2",[branching_rule("D3","A3","isomorphic"),"identity"])"""]
    elif CartanType(ct) == CartanType("B6"):
        rul = ["""D6:branching_rule("B6","D6","extended")""",
               """A1:branching_rule("B6","A1","symmetric_power")""",
               """A1xA1xB4:branching_rule("B6","D2xB4","orthogonal_sum")*branching_rule("D2xB4","A1xA1xB4",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """A1xD5:branching_rule("B6","B1xD5","orthogonal_sum")*branching_rule("B1xD5","A1xD5",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """A3xB3:branching_rule("B6","D3xB3","orthogonal_sum")*branching_rule("D3xB3","A3xB3",[branching_rule("D3","A3","isomorphic"),"identity"])""",
               """B2xD4:branching_rule("B6","B2xD4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("B7"):
        rul = ["""D7:branching_rule("B7","D7","extended")""",
               """A3:branching_rule("B7","A3(1,0,1)","plethysm")""",
               """A1:branching_rule("B7","A1","symmetric_power")""",
               """A1xB2:branching_rule("B7","B1xB2","tensor")*branching_rule("B1xB2","A1xB2",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """A1xD6:branching_rule("B7","B1xD6","extended")*branching_rule("B1xD6","A1xD6",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """A1xA1xB5:branching_rule("B7","D2xB5","extended")*branching_rule("D2xB5","A1xA1xB5",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """B2xD5:branching_rule("B7","B2xD5","orthogonal_sum")""",
               """A3xB4:branching_rule("B7","D3xB4","orthogonal_sum")*branching_rule("D3xB4","A3xB4",[branching_rule("D3","A3","isomorphic"),"identity"])""",
               """B3xD4:branching_rule("B7","B3xD4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("B8"):
        rul = ["""D8:branching_rule("B8","D8","extended")""",
               """A1:branching_rule("B8","A1","symmetric_power")""",
               """A1xD7:branching_rule("B8","B1xD7","orthogonal_sum")*branching_rule("B1xD7","A1xD7",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """A1xA1xB6:branching_rule("B8","D2xB6","orthogonal_sum")*branching_rule("D2xB6","A1xA1xB6",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """B2xD6:branching_rule("B8","B2xD6","orthogonal_sum")""",
               """A3xB5:branching_rule("B8","D3xB5","orthogonal_sum")*branching_rule("D3xB5","A3xB5",[branching_rule("D3","A3","isomorphic"),"identity"])""",
               """B3xD5:branching_rule("B8","B3xD5","orthogonal_sum")""",
               """B4xD4:branching_rule("B8","B4xD4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("C2"):
        rul = ["""A1:branching_rule("C2","A1","symmetric_power")""",
               """A1xA1:branching_rule("C2","C1xC1","orthogonal_sum")*branching_rule("C1xC1","A1xA1",[branching_rule("C1","A1","isomorphic"),branching_rule("C1","A1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("C3"):
        rul = ["""A2:branching_rule("C3","A2","levi")""",
               """A1:branching_rule("C3","A1","symmetric_power")""",
               """A1xA1:branching_rule("C3","B1xC1","tensor")*branching_rule("B1xC1","A1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("C1","A1","isomorphic")])""",
               """A1xC2:branching_rule("C3","C1xC2","orthogonal_sum")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])"""]
    elif CartanType(ct) == CartanType("C4"):
        rul = ["""A3:branching_rule("C4","A3","levi")""",
               """A1:branching_rule("C4","A1","symmetric_power")""",
               """A1xA3:branching_rule("C4","C1xC3","orthogonal_sum")*branching_rule("C1xC3","A1xA3",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """C2xC2:branching_rule("C4","C2xC2","orthogonal_sum")""",
               """A1xA1xA1:branching_rule("C4","C1xD2","tensor")*branching_rule("C1xD2","A1xA1xA1",[branching_rule("C1","A1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("C5"):
        rul = ["""A4:branching_rule("C5","A4","levi")""",
               """A1:branching_rule("C5","A1","symmetric_power")""",
               """A1xC4:branching_rule("C5","C1xC4","orthogonal_sum")*branching_rule("C1xC4","A1xC4",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """C2xC3:branching_rule("C5","C2xC3","orthogonal_sum")""",
               """A1xB2:branching_rule("C5","C1xB2","tensor")*branching_rule("C1xB2","A1xB2",[branching_rule("C1","A1","isomorphic"),"identity"])"""]
    elif CartanType(ct) == CartanType("C6"):
        rul = ["""A5:branching_rule("C6","A5","levi")""",
               """A1:branching_rule("C6","A1","symmetric_power")""",
               """A1xA3:branching_rule("C6","C1xD3","tensor")*branching_rule("C1xD3","A1xA3",[branching_rule("C1","A1","isomorphic"),branching_rule("D3","A3","isomorphic")])""",
               """A1xC2:branching_rule("C6","B1xC2","tensor")*branching_rule("B1xC2","A1xC2",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """A1xC5:branching_rule("C6","C1xC5","orthogonal_sum")*branching_rule("C1xC5","A1xC5",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """C2xC4:branching_rule("C6","C2xC4","orthogonal_sum")""",
               """C3xC3:branching_rule("C6","C3xC3","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("C7"):
        rul = ["""A6:branching_rule("C7","A6","levi")""",
               """A1:branching_rule("C7","A1","symmetric_power")""",
               """A1xB3:branching_rule("C7","C1xB3","tensor")*branching_rule("C1xB3","A1xB3",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """A1xC6:branching_rule("C7","C1xC6","orthogonal_sum")*branching_rule("C1xC6","A1xC6",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """C2xC5:branching_rule("C7","C2xC5","orthogonal_sum")""",
               """C3xC4:branching_rule("C7","C3xC4","orthogonal_sum")""",
               """C3:branching_rule("C7","C3(0,0,1)","plethysm") # overlooked by Patera and McKay"""]
    elif CartanType(ct) == CartanType("C8"):
        rul = ["""A7:branching_rule("C8","A7","levi")""",
               """A1:branching_rule("C8","A1","symmetric_power")""",
               """C2:branching_rule("C8","C2(1,1)","plethysm")""",
               """A1xD4:branching_rule("C8","C1xD4","tensor")*branching_rule("C1xD4","A1xD4",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """A1xC7:branching_rule("C8","C1xC7","orthogonal_sum")*branching_rule("C1xC7","A1xC7",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """C2xC6:branching_rule("C8","C2xC6","orthogonal_sum")""",
               """C3xC5:branching_rule("C8","C3xC5","orthogonal_sum")""",
               """C4xC4:branching_rule("C8","C4xC4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("D4"):
        rul = ["""B3:branching_rule("D4","B3","symmetric")""",
               """A2:branching_rule("D4","A2(1,1)","plethysm")""",
               """A1xC2:branching_rule("D4","C1xC2","tensor")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """A1xA1xA1xA1:branching_rule("D4","D2xD2","orthogonal_sum")*branching_rule("D2xD2","A1xA1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("D5"):
        rul = ["""A4:branching_rule("D5","A4","levi")""",
               """B4:branching_rule("D5","B4","symmetric")""",
               """C2:branching_rule("D5","C2(2,0)","plethysm")""",
               """A1xA1xA3:branching_rule("D5","D2xD3","orthogonal_sum")*branching_rule("D2xD3","A1xA1xA3",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D3","A3","isomorphic")])""",
               """A1xA3:branching_rule("D5","B1xB3","orthogonal_sum")*branching_rule("B1xB3","A1xA3",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """B2xB2:branching_rule("D5","B2xB2","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("D6"):
        rul = ["""A5:branching_rule("D6","A5","levi")""",
               """B5:branching_rule("D6","B5","symmetric")""",
               """A1xA3:branching_rule("D6","C1xC3","tensor")*branching_rule("C1xC3","A1xA3",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """A1xA1xD4:branching_rule("D6","D2xD4","orthogonal_sum")*branching_rule("D2xD4","A1xA1xD4",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """A3xA3:branching_rule("D6","D3xD3","orthogonal_sum")*branching_rule("D3xD3","A3xA3",[branching_rule("D3","A3","isomorphic"),branching_rule("D3","A3","isomorphic")])""",
               """A1xB4:branching_rule("D6","B1xB4","orthogonal_sum")*branching_rule("B1xB4","A1xB4",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """B2xB3:branching_rule("D6","B2xB3","orthogonal_sum")""",
               """A1xA1xA1:branching_rule("D6","B1xD2","tensor")*branching_rule("B1xD2","A1xA1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("D7"):
        rul = ["""A6:branching_rule("D7","A6","levi")""",
               """B6:branching_rule("D7","B6","symmetric")""",
               """C3:branching_rule("D7","C3(0,1,0)","plethysm")""",
               """C2:branching_rule("D7","C2(0,2)","plethysm")""",
               """G2:branching_rule("D7","G2(0,1)","plethysm")""",
               """A1xA1xD5:branching_rule("D7","D2xD5","orthogonal_sum")*branching_rule("D2xD5","A1xA1xD5",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """A3xD4:branching_rule("D7","D3xD4","orthogonal_sum")*branching_rule("D3xD4","A3xD4",[branching_rule("D3","A3","isomorphic"),"identity"])""",
               """A1xB5:branching_rule("D7","B1xB5","orthogonal_sum")*branching_rule("B1xB5","A1xB5",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """B2xB4:branching_rule("D7","B2xB4","orthogonal_sum")""",
               """B3xB3:branching_rule("D7","B3xB3","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("D8"):
        rul = ["""A7:branching_rule("D8","A7","levi")""",
               """B7:branching_rule("D8","B7","symmetric")""",
               """B4:branching_rule("D8","B4(0,0,0,1)","plethysm")""",
               """A1xC4:branching_rule("D8","C1xC4","tensor")*branching_rule("C1xC4","A1xC4",[branching_rule("C1","A1","isomorphic"),"identity"])""",
               """A1xA1xD6:branching_rule("D8","D2xD6","orthogonal_sum")*branching_rule("D2xD6","A1xA1xD6",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""",
               """A3xD5:branching_rule("D8","D3xD5","orthogonal_sum")*branching_rule("D3xD5","A3xD5",[branching_rule("D3","A3","isomorphic"),"identity"])""",
               """D4xD4:branching_rule("D8","D4xD4","orthogonal_sum")""",
               """A1xB6:branching_rule("D8","B1xB6","orthogonal_sum")*branching_rule("B1xB6","A1xB6",[branching_rule("B1","A1","isomorphic"),"identity"])""",
               """B2xB5:branching_rule("D8","B2xB5","orthogonal_sum")""",
               """B3xB4:branching_rule("D8","B3xB4","orthogonal_sum")""",
               """C2xC2:branching_rule("D8","C2xC2","tensor")"""]
    elif CartanType(ct) == CartanType("G2"):
        rul = ["""A2:branching_rule("G2","A2","extended")""",
               """A1:branching_rule("G2","A1","i")""",
               """A1xA1:branching_rule("G2","A1xA1","extended")"""]
    elif CartanType(ct) == CartanType("F4"):
        rul = ["""B4:branching_rule("F4","B4","extended")""",
               """A1:branching_rule("F4","A1","ii")""",
               """A1xG2:branching_rule("F4","A1xG2","miscellaneous")""",
               """A1xC3:branching_rule("F4","A1xC3","extended")""",
               """A2xA2:branching_rule("F4","A2xA2","extended")"""]
    elif CartanType(ct) == CartanType("E6"):
        rul = ["""D5:branching_rule("E6","D5","levi")""",
               """C4:branching_rule("E6","C4","symmetric")""",
               """F4:branching_rule("E6","F4","symmetric")""",
               """A2:branching_rule("E6","A2","miscellaneous")""",
               """G2:branching_rule("E6","G2","miscellaneous")""",
               """A2xG2:branching_rule("E6","A2xG2","miscellaneous")""",
               """A1xA5:branching_rule("E6","A1xA5","extended")""",
               """A2xA2xA2:branching_rule("E6","A2xA2xA2","extended")"""]
    elif CartanType(ct) == CartanType("E7"):
        rul = ["""A7:branching_rule("E7","A7","extended")""",
               """E6:branching_rule("E7","E6","levi")""",
               """A2:branching_rule("E7","A2","miscellaneous")""",
               """A1:branching_rule("E7","A1","iii")""",
               """A1:branching_rule("E7","A1","iv")""",
               """A1xF4:branching_rule("E7","A1xF4","miscellaneous")""",
               """G2xC3:branching_rule("E7","G2xC3","miscellaneous")""",
               """A1xG2:branching_rule("E7","A1xG2","miscellaneous")""",
               """A1xA1:branching_rule("E7","A1xA1","miscellaneous")""",
               """A1xD6:branching_rule("E7","A1xD6","extended")""",
               """A5xA2:branching_rule("E7","A5xA2","extended")"""]
    elif CartanType(ct) == CartanType("E8"):
        rul = ["""A4xA4:branching_rule("E8","A4xA4","extended")""",
               """G2xF4:branching_rule("E8","G2xF4","miscellaneous")""",
               """E6xA2:branching_rule("E8","E6xA2","extended")""",
               """E7xA1:branching_rule("E8","E7xA1","extended")""",
               """D8:branching_rule("E8","D8","extended")""",
               """A8:branching_rule("E8","A8","extended")""",
               """B2:branching_rule("E8","B2","miscellaneous")""",
               """A1xA2:branching_rule("E8","A1xA2","miscellaneous")""",
               """A1:branching_rule("E8","A1","v")""",
               """A1:branching_rule("E8","A1","vi")""",
               """A1:branching_rule("E8","A1","vii")"""]
    else:
        raise ValueError("Argument must be an irreducible classical Cartan Type with rank less than or equal to 8")
    if mode == "print_rules":
        for line in rul:
            print(line)
    elif mode == "get_rule":
        d = {}
        for line in rul:
            [k, br] = line.split(":")
            br = eval(br)
            if k in d:
                if not isinstance(d[k], list):
                    d[k] = [d[k]]
                d[k].append(br)
            else:
                d[k] = br
        return d
