.. _chapter-groups:

******
Groups
******

.. index::
   pair: group; permutation

.. _section-permutation:

Permutation groups
==================

A permutation group is a subgroup of some symmetric group
:math:`S_n`. Sage has a Python class ``PermutationGroup``, so you
can work with such groups directly::

    sage: G = PermutationGroup(['(1,2,3)(4,5)'])
    sage: G
    Permutation Group with generators [(1,2,3)(4,5)]
    sage: g = G.gens()[0]; g
    (1,2,3)(4,5)
    sage: g*g
    (1,3,2)
    sage: G = PermutationGroup(['(1,2,3)'])
    sage: g = G.gens()[0]; g
    (1,2,3)
    sage: g.order()
    3

For the example of the Rubik's cube group (a permutation subgroup
of :math:`S_{48}`, where the non-center facets of the Rubik's
cube are labeled :math:`1,2,...,48` in some fixed way), you can
use the GAP-Sage interface as follows.

.. index::
   pair: group; Rubik's cube

.. skip

::

    sage: cube = "cubegp := Group(
    ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
    ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
    (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
    (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
    (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
    (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"
    sage: gap(cube)
    'permutation group with 6 generators'
    sage: gap("Size(cubegp)")
    43252003274489856000'

Another way you can choose to do this:

-  Create a file ``cubegroup.py`` containing the
   lines::

       cube = "cubegp := Group(
       ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
       ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
       (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
       (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
       (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
       (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"

   Then place the file in the subdirectory
   ``$SAGE_ROOT/local/lib/python2.4/site-packages/sage`` of your Sage home
   directory. Last, read (i.e., ``import``) it into Sage:

   .. skip

   ::

       sage: import sage.cubegroup
       sage: sage.cubegroup.cube
       'cubegp := Group(( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)
       (11,35,27,19),( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)
       ( 6,22,46,35),(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)
       ( 8,30,41,11),(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)
       ( 8,33,48,24),(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)
       ( 1,14,48,27),(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)
       (16,24,32,40) )'
       sage: gap(sage.cubegroup.cube)
       'permutation group with 6 generators'
       sage: gap("Size(cubegp)")
       '43252003274489856000'

   (You will have line wrap instead of the above carriage returns in
   your Sage output.)

-  Use the ``CubeGroup`` class::

       sage: rubik = CubeGroup()
       sage: rubik
       The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).
       sage: rubik.order()
       43252003274489856000

   (1) has implemented classical groups (such as :math:`GU(3,\GF{5})`)
   and matrix groups over a finite field with user-defined generators.

   (2) also has implemented finite and infinite (but finitely
   generated) abelian groups.

.. index::
   pair: group; conjugacy classes

.. _section-conjugacy:

Conjugacy classes
=================

You can compute conjugacy classes of a finite group using "natively"::

    sage: G = PermutationGroup(['(1,2,3)', '(1,2)(3,4)', '(1,7)'])
    sage: CG = G.conjugacy_classes_representatives()
    sage: gamma = CG[2]
    sage: CG; gamma
    [(), (4,7), (3,4,7), (2,3)(4,7), (2,3,4,7), (1,2)(3,4,7), (1,2,3,4,7)]
    (3,4,7)

You can use the Sage-GAP interface::

    sage: libgap.eval("G := Group((1,2)(3,4),(1,2,3))")
    Group([ (1,2)(3,4), (1,2,3) ])
    sage: libgap.eval("CG := ConjugacyClasses(G)")
    [ ()^G, (2,3,4)^G, (2,4,3)^G, (1,2)(3,4)^G ]
    sage: libgap.eval("gamma := CG[3]")
    (2,4,3)^G
    sage: libgap.eval("g := Representative(gamma)")
    (2,4,3)

Or, here's another (more "pythonic") way to do this type of computation::

    sage: G = libgap.eval("Group([(1,2,3), (1,2)(3,4), (1,7)])")
    sage: CG = G.ConjugacyClasses()
    sage: gamma = CG[2]
    sage: g = gamma.Representative()
    sage: CG; gamma; g
    [ ()^G, (4,7)^G, (3,4,7)^G, (2,3)(4,7)^G, (2,3,4,7)^G, (1,2)(3,4,7)^G, (1,2,3,4,7)^G ]
    (3,4,7)^G
    (3,4,7)

.. index::
   pair: group; normal subgroups

.. _section-normal:

Normal subgroups
================

If you want to find all the normal subgroups of a permutation group
:math:`G` (up to conjugacy), you can use Sage's interface to GAP::

    sage: G = AlternatingGroup( 5 )
    sage: libgap(G).NormalSubgroups()
    [ Alt( [ 1 .. 5 ] ), Group(()) ]

or

::

    sage: G = libgap.AlternatingGroup( 5 )
    sage: G.NormalSubgroups()
    [ Alt( [ 1 .. 5 ] ), Group(()) ]

Here's another way, working more directly with GAP::

    sage: libgap.eval("G := AlternatingGroup( 5 )")
    Alt( [ 1 .. 5 ] )
    sage: libgap.eval("normal := NormalSubgroups( G )")
    [ Alt( [ 1 .. 5 ] ), Group(()) ]
    sage: G = libgap.eval("DihedralGroup( 10 )")
    sage: G.NormalSubgroups().SortedList()
    [ Group([  ]), Group([ f2 ]), <pc group of size 10 with 2 generators> ]
    sage: libgap.eval("G := SymmetricGroup( 4 )")
    Sym( [ 1 .. 4 ] )
    sage: libgap.eval("normal := NormalSubgroups( G );")
    [ Sym( [ 1 .. 4 ] ), Alt( [ 1 .. 4 ] ), Group([ (1,4)(2,3),  ... ]),
          Group(()) ]

.. index::
   pair: groups; center

.. _section-center:

Centers
=======

How do you compute the center of a group in Sage?

Although Sage calls GAP to do the computation of the group center,
``center`` is "wrapped" (i.e., Sage has a class PermutationGroup with
associated class method "center"), so the user does not need to use
the ``libgap`` command. Here's an example::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.center()
    Subgroup generated by [()] of (Permutation Group with generators [(3,4), (1,2,3)(4,5)])

A similar syntax for matrix groups also works::

    sage: G = SL(2, GF(5) )
    sage: G.center()
    Subgroup with 1 generators (
    [4 0]
    [0 4]
    ) of Special Linear Group of degree 2 over Finite Field of size 5
    sage: G = PSL(2, 5 )
    sage: G.center()
    Subgroup generated by [()] of (The projective special linear group of degree 2 over Finite Field of size 5)

.. NOTE:: ``center`` can be spelled either way in GAP, not so in Sage.

The group id database
=====================

The function ``group_id`` uses the Small Groups Library of
E. A. O'Brien, B. Eick, and H. U. Besche, which is a part of GAP.

::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.order()
    120
    sage: G.group_id()
    [120, 34]

Another example of using the small groups database: ``group_id``

.. skip

::

    sage: gap_console()
    ┌───────┐   GAP 4.10.0 of 01-Nov-2018
    │  GAP  │   https://www.gap-system.org
    └───────┘   Architecture: x86_64-pc-linux-gnu-default64
    Configuration:  gmp 6.0.0, readline
    Loading the library and packages ...
    Packages:   GAPDoc 1.6.2, PrimGrp 3.3.2, SmallGrp 1.3, TransGrp 2.0.4
    Try '??help' for help. See also '?copyright', '?cite' and '?authors'
    gap> G:=Group((4,6,5)(7,8,9),(1,7,2,4,6,9,5,3));
    Group([ (4,6,5)(7,8,9), (1,7,2,4,6,9,5,3) ])
    gap> StructureDescription(G);
    "(C3 x C3) : GL(2,3)"

Construction instructions for every group of order less than 32
===============================================================

AUTHORS:

* Davis Shurbert

Every group of order less than 32 is implemented in Sage as a permutation
group. They can all be created easily. We will first show how to build direct
products and semidirect products, then give the commands necessary to build
all of these small groups.

Let ``G1``, ``G2``, ..., ``Gn`` be permutation groups already initialized in
Sage. The following command can be used to take their direct product (where,
of course, the ellipses are simply being used here as a notation, and you
actually must enter every factor in your desired product explicitly).

.. skip

::

    sage: G = direct_product_permgroups([G1, G2, ..., Gn])

The semidirect product operation can be thought of as a generalization of the
direct product operation. Given two groups, `H` and `K`, their semidirect
product, `H \ltimes_{\phi} K`, (where `\phi : H \rightarrow Aut(K)` is a
homomorphism) is a group whose underlying set is the cartersian product of
`H` and `K`, but with the operation:

.. MATH::

    (h_1, k_1) (h_2, k_2) = (h_1 h_2, k_1^{\phi(h_2)} k_2).

The output is not the group explicitly described in the definition of the
operation, but rather an isomorphic group of permutations. In the routine
below, assume ``H`` and ``K`` already have been defined and initialized in
Sage. Also, ``phi`` is a list containing two sublists that define the
underlying homomorphism by giving the images of a set of generators of ``H``.
For each semidirect product in the table below we will show you how to build
``phi``, then assume you have read this passage and understand how to go
from there.

.. skip

::

    sage: G = H.semidirect_product(K, phi)

To avoid unnecessary repetition, we will now give commands that one can use to
create the cyclic group of order `n`, `C_n`, and the dihedral group on `n`
letters, `D_n`. We will present one more example of each to ensure that the
reader understands the command, then it will be withheld.

.. skip

::

    sage: G = CyclicPermutationGroup(n)

    sage: G = DihedralGroup(n)

Note that exponential notation will be used for the direct product operation.
For example, `{C_2}^2 = C_2 \times C_2`. This table was crafted with the help
of *Group Tables*, by AD Thomas and GV Wood (1980, Shiva Publishing).


===== =============================================== =============================================================================================== ===========================
Order Group Description                                Command(s)                                                                                     GAP ID
===== =============================================== =============================================================================================== ===========================
1     The Trivial Group                               ::                                                                                              [1,1]

                                                        sage: G = SymmetricGroup(1)
2     `C_2`                                           ::                                                                                              [2,1]

                                                        sage: G = SymmetricGroup(2)
3     `C_3`                                           ::                                                                                              [3,1]

                                                        sage: G = CyclicPermutationGroup(3)
4     `C_4`                                                                                                                                           [4,1]
4     `C_2 \times C_2`                                ::                                                                                              [4,2]

                                                        sage: G = KleinFourGroup()
5     `C_5`                                                                                                                                           [5,1]
6     `C_6`                                                                                                                                           [6,2]
6     `S_3` (Symmetric Group on 3 letters)            ::                                                                                              [6,1]

                                                        sage: G = SymmetricGroup(3)
7     `C_7`                                                                                                                                           [7,1]
8     `C_8`                                                                                                                                           [8,1]
8     `C_4 \times C_2`                                                                                                                                [8,2]
8     `C_2\times C_2\times C_2`                                                                                                                       [8,5]
8     `D_4`                                           ::                                                                                              [8,3]

                                                        sage: G = DihedralGroup(4)
8     The Quaternion Group (Q)                        ::                                                                                              [8,4]

                                                        sage: G = QuaternionGroup()
9     `C_9`                                                                                                                                           [9,1]
9     `C_3 \times C_3`                                                                                                                                [9,2]
10    `C_{10}`                                                                                                                                        [10,2]
10    `D_5`                                                                                                                                           [10,1]
11    `C_{11}`                                                                                                                                        [11,1]
12    `C_{12}`                                                                                                                                        [12,2]
12    `C_6 \times C_2`                                                                                                                                [12,5]
12    `D_6`                                                                                                                                           [12,4]
12    `A_4` (Alternating Group on 4 letters)          ::                                                                                              [12,3]

                                                        sage: G = AlternatingGroup(4)
12    `Q_6` (DiCyclic group of order 12)              ::                                                                                              [12,1]

                                                        sage: G = DiCyclicGroup(3)
13    `C_{13}`                                                                                                                                        [13,1]
14    `C_{14}`                                                                                                                                        [14,2]
14    `D_{7}`                                                                                                                                         [14,1]
15    `C_{15}`                                                                                                                                        [15,1]
16    `C_{16}`                                                                                                                                        [16,1]
16    `C_8 \times C_2`                                                                                                                                [16,5]
16    `C_4 \times C_4`                                                                                                                                [16,2]
16    `C_4\times C_2\times C_2`                                                                                                                       [16,10]
16    `{C_2}^4`                                                                                                                                       [16,14]
16    `D_4 \times C_2`                                                                                                                                [16,11]
16    `Q \times C_2`                                                                                                                                  [16,12]
16    `D_8`                                                                                                                                           [16,7]
16    `Q_{8}` (Dicyclic group of order 16)            ::                                                                                              [16,9]

                                                        sage: G = DiCyclicGroup(4)
16    Semidihedral Group of order `2^4`               ::                                                                                              [16,8]

                                                        sage: G = SemidihedralGroup(4)
16    Split Metacyclic Group of order `2^4`           ::                                                                                              [16,6]

                                                        sage: G = SplitMetacyclicGroup(2,4)
16    `(C_4 \times C_2) \rtimes_{\phi} C_2`           ::                                                                                              [16,13]

                                                        sage: C2 = SymmetricGroup(2); C4 = CyclicPermutationGroup(4)
                                                        sage: A = direct_product_permgroups([C2,C4])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0],A.gens()[0]^2*A.gens()[1]])
                                                        sage: phi = [[(1,2)],[alpha]]
16    `(C_4 \times C_2) \rtimes_{\phi} C_2`           ::                                                                                              [16,3]

                                                        sage: C2 = SymmetricGroup(2); C4 = CyclicPermutationGroup(4)
                                                        sage: A = direct_product_permgroups([C2,C4])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0]^3*A.gens()[1],A.gens()[1]])
                                                        sage: phi = [[(1,2)],[alpha]]
16    `C_4 \rtimes_{\phi} C_4`                        ::                                                                                              [16,4]

                                                        sage: C4 = CyclicPermutationGroup(4)
                                                        sage: alpha = PermutationGroupMorphism(C4,C4,[C4.gen().inverse()])
                                                        sage: phi = [[(1,2,3,4)],[alpha]]
17    `C_{17}`                                                                                                                                        [17,1]
18    `C_{18}`                                                                                                                                        [18,2]
18    `C_6 \times C_3`                                                                                                                                [18,5]
18    `D_9`                                                                                                                                           [18,1]
18    `S_3 \times C_3`                                                                                                                                [18,3]
18    `Dih(C_3 \times C_3)`                           ::                                                                                              [18,4]

                                                        sage: G = GeneralDihedralGroup([3,3])
19    `C_{19}`                                                                                                                                        [19,1]
20    `C_{20}`                                                                                                                                        [20,2]
20    `C_{10} \times C_2`                                                                                                                             [20,5]
20    `D_{10}`                                                                                                                                        [20,4]
20    `Q_{10}` (Dicyclic Group of order 20)                                                                                                           [20,1]
20    `Hol(C_5)`                                      ::                                                                                              [20,3]

                                                        sage: C5 = CyclicPermutationGroup(5)
                                                        sage: G = C5.holomorph()
21    `C_{21}`                                                                                                                                        [21,2]
21    `C_7 \rtimes_{\phi} C_3`                        ::                                                                                              [21,1]

                                                        sage: C7 = CyclicPermutationGroup(7)
                                                        sage: alpha = PermutationGroupMorphism(C7,C7,[C7.gen()**4])
                                                        sage: phi = [[(1,2,3)],[alpha]]
22    `C_{22}`                                                                                                                                        [22,2]
22    `D_{11}`                                                                                                                                        [22,1]
23    `C_{23}`                                                                                                                                        [23,1]
24    `C_{24}`                                                                                                                                        [24,2]
24    `D_{12}`                                                                                                                                        [24,6]
24    `Q_{12}` (DiCyclic Group of order 24)                                                                                                           [24,4]
24    `C_{12} \times C_2`                                                                                                                             [24,9]
24    `C_6 \times C_2 \times C_2`                                                                                                                     [24,15]
24    `S_4` (Symmetric Group on 4 letters)            ::                                                                                              [24,12]

                                                        sage: G = SymmetricGroup(4)
24    `S_3 \times C_4`                                                                                                                                [24,5]
24    `S_3 \times C_2 \times C_2`                                                                                                                     [24,14]
24    `D_4 \times C_3`                                                                                                                                [24,10]
24    `Q \times C_3`                                                                                                                                  [24,11]
24    `A_4 \times C_2`                                                                                                                                [24,13]
24    `Q_6 \times C_2`                                                                                                                                [24,7]
24    `Q \rtimes_{\phi} C_3`                          ::                                                                                              [24,3]

                                                        sage: Q = QuaternionGroup()
                                                        sage: alpha = PermutationGroupMorphism(Q,Q,[Q.gens()[0]*Q.gens()[1],Q.gens()[0].inverse()])
                                                        sage: phi = [[(1,2,3)],[alpha]]
24    `C_3 \rtimes_{\phi} C_8`                        ::                                                                                              [24,1]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: alpha = PermutationGroupMorphism(C3,C3,[C3.gen().inverse()])
                                                        sage: phi = [[(1,2,3,4,5,6,7,8)],[alpha]]
24    `C_3 \rtimes_{\phi} D_4`                        ::                                                                                              [24,8]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: alpha1 = PermutationGroupMorphism(C3,C3,[C3.gen().inverse()])
                                                        sage: alpha2 = PermutationGroupMorphism(C3,C3,[C3.gen()])
                                                        sage: phi = [[(1,2,3,4),(1,3)],[alpha1,alpha2]]
25    `C_{25}`                                                                                                                                        [25,1]
25    `C_5 \times C_5`                                                                                                                                [25,2]
26    `C_{26}`                                                                                                                                        [26,2]
26    `D_{13}`                                                                                                                                        [26,1]
27    `C_{27}`                                                                                                                                        [27,1]
27    `C_9 \times C_3`                                                                                                                                [27,2]
27    `C_3 \times C_3 \times C_3`                                                                                                                     [27,5]
27    Split Metacyclic Group of order `3^3`           ::                                                                                              [27,4]

                                                        sage: G = SplitMetacyclicGroup(3,3)
27    `(C_3 \times C_3) \rtimes_{\phi} C_3`           ::                                                                                              [27,3]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: A = direct_product_permgroups([C3,C3])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0]*A.gens()[1].inverse(),A.gens()[1]])
                                                        sage: phi = [[(1,2,3)],[alpha]]
28    `C_{28}`                                                                                                                                        [28,2]
28    `C_{14} \times C_2`                                                                                                                             [28,4]
28    `D_{14}`                                                                                                                                        [28,3]
28    `Q_{14}` (DiCyclic Group of order 28)                                                                                                           [28,1]
29    `C_{29}`                                                                                                                                        [29,1]
30    `C_{30}`                                                                                                                                        [30,4]
30    `D_{15}`                                                                                                                                        [30,3]
30    `D_5 \times C_3`                                                                                                                                [30,2]
30    `D_3 \times C_5`                                                                                                                                [30,1]
31    `C_{31}`                                                                                                                                        [31,1]
===== =============================================== =============================================================================================== ===========================


Table By Kevin Halasz

Construction instructions for every finitely presented group of order 15 or less
================================================================================

Sage has the capability to easily construct every group of order 15 or less
as a finitely presented group. We will begin with some discussion on creating
finitely generated abelian groups, as well as direct and semidirect products
of finitely presented groups.

All finitely generated abelian groups can be created using the
``groups.presentation.FGAbelian(ls)`` command, where ``ls`` is a list of
non-negative integers which gets reduced to invariants defining the group
to be returned. For example, to construct
`C_4 \times C_2 \times C_2 \times C_2` we can simply use::

    sage: A = groups.presentation.FGAbelian([4,2,2,2])

The output for a given group is the same regardless of the input list of
integers.  The following example yields identical presentations for the
cyclic group of order 30.
::

    sage: A = groups.presentation.FGAbelian([2,3,5])
    sage: B = groups.presentation.FGAbelian([30])

If ``G`` and ``H`` are finitely presented groups, we can use the following
code to create the direct product of ``G`` and ``H``, `G \times H`.

.. skip

::

    sage: D = G.direct_product(H)

Suppose there exists a homomorphism `\phi` from a group `G` to the
automorphism group of a group `H`. Define the semidirect product of `G`
with `H` via `\phi`, as the Cartesian product of `G` and `H`, with the
operation `(g_1, h_1)(g_2, h_2) = (g_1 g_2, \phi_{h_1}(g_2) h_2)` where
`\phi_h = \phi(h)`. To construct this product in Sage for two finitely
presented groups, we must define `\phi` manually using a pair of lists. The
first list consists of generators of the group `G`, while the second list
consists of images of the corresponding generators in the first list. These
automorphisms are similarly defined as a pair of lists, generators in one
and images in the other. As an example, we construct the dihedral group of
order 16 as a semidirect product of cyclic groups.
::

    sage: C2 = groups.presentation.Cyclic(2)
    sage: C8 = groups.presentation.Cyclic(8)
    sage: hom = (C2.gens(), [ ([C8([1])], [C8([-1])]) ])
    sage: D = C2.semidirect_product(C8, hom)

The following table shows the groups of order 15 or less, and how to construct
them in Sage. Repeated commands have been omitted but instead are described
by the following examples.

The cyclic group of order `n` can be created with a single command:

.. skip

::

    sage: C = groups.presentation.Cyclic(n)

Similarly for the dihedral group of order `2n`:

.. skip

::

    sage: D = groups.presentation.Dihedral(n)

This table was modeled after the preceding table created by Kevin Halasz.


===== =============================================== =============================================================================================== ===========================
Order Group Description                                Command(s)                                                                                     GAP ID
===== =============================================== =============================================================================================== ===========================
1     The Trivial Group                               ::                                                                                              [1,1]

                                                        sage: G = groups.presentation.Symmetric(1)

2     `C_2`                                           ::                                                                                              [2,1]

                                                        sage: G = groups.presentation.Symmetric(2)

3     `C_3`                                           ::                                                                                              [3,1]

                                                        sage: G = groups.presentation.Cyclic(3)

4     `C_4`                                                                                                                                           [4,1]

4     `C_2 \times C_2`                                ::                                                                                              [4,2]

                                                        sage: G = groups.presentation.Klein()

5     `C_5`                                                                                                                                           [5,1]
6     `C_6`                                                                                                                                           [6,2]

6     `S_3` (Symmetric Group on 3 letters)            ::                                                                                              [6,1]

                                                        sage: G = groups.presentation.Symmetric(3)

7     `C_7`                                                                                                                                           [7,1]
8     `C_8`                                                                                                                                           [8,1]

8     `C_4 \times C_2`                                ::                                                                                              [8,2]

                                                        sage: G = groups.presentation.FGAbelian([4,2])

8     `C_2\times C_2\times C_2`                       ::                                                                                              [8,5]

                                                        sage: G = groups.presentation.FGAbelian([2,2,2])

8     `D_4`                                           ::                                                                                              [8,3]

                                                        sage: G = groups.presentation.Dihedral(4)

8     The Quaternion Group (Q)                        ::                                                                                              [8,4]

                                                        sage: G = groups.presentation.Quaternion()

9     `C_9`                                                                                                                                           [9,1]
9     `C_3 \times C_3`                                                                                                                                [9,2]
10    `C_{10}`                                                                                                                                        [10,2]
10    `D_5`                                                                                                                                           [10,1]
11    `C_{11}`                                                                                                                                        [11,1]
12    `C_{12}`                                                                                                                                        [12,2]
12    `C_6 \times C_2`                                                                                                                                [12,5]
12    `D_6`                                                                                                                                           [12,4]
12    `A_4` (Alternating Group on 4 letters)          ::                                                                                              [12,3]

                                                        sage: G = groups.presentation.Alternating(4)

12    `Q_6` (DiCyclic group of order 12)              ::                                                                                              [12,1]

                                                        sage: G = groups.presentation.DiCyclic(3)

13    `C_{13}`                                                                                                                                        [13,1]
14    `C_{14}`                                                                                                                                        [14,2]
14    `D_{7}`                                                                                                                                         [14,1]
15    `C_{15}`                                                                                                                                        [15,1]
===== =============================================== =============================================================================================== ===========================

