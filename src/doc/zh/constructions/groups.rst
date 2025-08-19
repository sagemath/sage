.. _chapter-groups:

**
群
**

.. index::
   pair: group; permutation

.. _section-permutation:

置换群
======

置换群是某一对称群 `S_n` 的子群。Sage 有一个 Python 类 ``PermutationGroup``，
因此你可以直接使用此类群::

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

对于魔方群（`S_{48}` 的置换子群，其中魔方的非中心面以某种固定方式标记为 `1,2,...,48`），
你可以按如下方式使用 GAP-Sage 接口。

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

你还可以选择另一种方式来实现：

-  创建一个包含以下内容的文件 ``cubegroup.py``::

       cube = "cubegp := Group(
       ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
       ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
       (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
       (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
       (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
       (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"

   然后将该文件放置在 Sage 目录的
   ``$SAGE_ROOT/local/lib/python2.4/site-packages/sage``
   子目录中。最后，读取（即 ``import``）该文件到 Sage 中：

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

   （在 Sage 输出中，会使用换行符来代替上面的回车符。）

-  使用 ``CubeGroup`` 类::

       sage: rubik = CubeGroup()
       sage: rubik
       The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).
       sage: rubik.order()
       43252003274489856000

   (1) Sage 实现了经典群（如 `GU(3,\GF{5})`）
   和具有用户定义生成器的有限域上的矩阵群。

   (2) Sage 还实现了有限和无限
   （但有限生成）阿贝尔群。

.. index::
   pair: group; conjugacy classes

.. _section-conjugacy:

共轭类
======

你可以“原生地”计算有限群的共轭类::

    sage: G = PermutationGroup(['(1,2,3)', '(1,2)(3,4)', '(1,7)'])
    sage: CG = G.conjugacy_classes_representatives()
    sage: gamma = CG[2]
    sage: CG; gamma
    [(), (4,7), (3,4,7), (2,3)(4,7), (2,3,4,7), (1,2)(3,4,7), (1,2,3,4,7)]
    (3,4,7)

你可以使用 Sage-GAP 接口完成这一任务::

    sage: libgap.eval("G := Group((1,2)(3,4),(1,2,3))")
    Group([ (1,2)(3,4), (1,2,3) ])
    sage: libgap.eval("CG := ConjugacyClasses(G)")
    [ ()^G, (2,3,4)^G, (2,4,3)^G, (1,2)(3,4)^G ]
    sage: libgap.eval("gamma := CG[3]")
    (2,4,3)^G
    sage: libgap.eval("g := Representative(gamma)")
    (2,4,3)

或者，这里有另一种（更符合 Python 风格的）方法来进行该计算::

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

正规子群
========

如果想要找到置换群 `G` （从共轭角度）的所有正规子群，可以使用 Sage 的 GAP 接口::

    sage: G = AlternatingGroup( 5 )
    sage: libgap(G).NormalSubgroups()
    [ Alt( [ 1 .. 5 ] ), Group(()) ]

或者

::

    sage: G = libgap.AlternatingGroup( 5 )
    sage: G.NormalSubgroups()
    [ Alt( [ 1 .. 5 ] ), Group(()) ]

这里有另一种更直接使用 GAP 的方法::

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

中心
====

如何在 Sage 中计算群的中心？

虽然 Sage 调用 GAP 来计算群的中心，
但 ``center`` 是“封装”过的方法（即 Sage 有一个类 PermutationGroup 关联 "center" 方法），
因此用户不需要使用 ``libgap`` 命令。这里有一个例子::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.center()
    Subgroup generated by [()] of (Permutation Group with generators [(3,4), (1,2,3)(4,5)])

类似的语法也适用于矩阵群::

    sage: G = SL(2, GF(5) )
    sage: G.center()
    Subgroup with 1 generators (
    [4 0]
    [0 4]
    ) of Special Linear Group of degree 2 over Finite Field of size 5
    sage: G = PSL(2, 5 )
    sage: G.center()
    Subgroup generated by [()] of (The projective special linear group of degree 2 over Finite Field of size 5)

.. NOTE:: 在 GAP 中 ``center`` 有两种拼写方式，但在 Sage 中不行。

群 id 数据库
============

函数 ``group_id`` 使用了 E. A. O'Brien、B. Eick 和 H. U. Besche 的小群库，它是 GAP 的一部分。

::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.order()
    120
    sage: G.group_id()
    [120, 34]

另一个使用小型群数据库的例子：``group_id``

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

小于 32 阶的群的构建指令
========================

作者：

* Davis Shurbert

每个小于 32 阶的群都在 Sage 中实现为置换群。这些群的构建都非常简单。
我们首先展示如何构建直积和半直积，然后给出构建这些小群所需的命令。

设 ``G1``, ``G2``, ..., ``Gn`` 是已经在 Sage 中初始化的置换群。
可以使用以下命令取它们的直积
（当然，这里省略号只是作为符号使用，实际上必须显式输入所求乘积中的每个因子）。

.. skip

::

    sage: G = direct_product_permgroups([G1, G2, ..., Gn])

半直积运算可以被视为直积运算的推广。给定两个群 `H` 和 `K`，它们的半直积 `H \ltimes_{\phi} K`
（其中 `\phi : H \rightarrow Aut(K)` 是一个同态）是一个群，其基础集合是 `H` 和 `K` 的笛卡尔积，
但具有以下运算：

.. MATH::

    (h_1, k_1) (h_2, k_2) = (h_1 h_2, k_1^{\phi(h_2)} k_2).

输出不是运算定义中明确描述的群，而是一个同构的置换群。
在下面的例程中，假设 ``H`` 和 ``K`` 已经在 Sage 中定义且初始化。
此外，``phi`` 是一个包含两个子列表的列表，通过给出 ``H`` 的生成器集合的像来定义底层同态。
对于下表中的每个半直积群，我们将展示如何构建 ``phi``，然后假设你已经阅读此段落并理解如何从那里开始。

.. skip

::

    sage: G = H.semidirect_product(K, phi)

为了避免不必要的重复，我们现在将给出创建 `n` 阶循环群 `C_n` 的命令和 `n` 个字母的二面体群 `D_n` 的命令。
我们还会为每个命令展示一个例子以确保读者理解这些命令，然后不再重复。

.. skip

::

    sage: G = CyclicPermutationGroup(n)

    sage: G = DihedralGroup(n)

请注意，直积运算中将使用指数表示法。例如 `{C_2}^2 = C_2 \times C_2`。
该表格是在 AD Thomas 和 GV Wood 的 *Group Tables* (1980, Shiva Publishing) 的帮助下制作的。


===== =============================================== =============================================================================================== ===========================
阶     群描述                                           命令                                                                                            GAP ID
===== =============================================== =============================================================================================== ===========================
1     平凡群                                           ::                                                                                              [1,1]

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
6     `S_3` （三字母对称群）                           ::                                                                                              [6,1]

                                                        sage: G = SymmetricGroup(3)
7     `C_7`                                                                                                                                           [7,1]
8     `C_8`                                                                                                                                           [8,1]
8     `C_4 \times C_2`                                                                                                                                [8,2]
8     `C_2\times C_2\times C_2`                                                                                                                       [8,5]
8     `D_4`                                           ::                                                                                              [8,3]

                                                        sage: G = DihedralGroup(4)
8     四元群 (Q)                                       ::                                                                                              [8,4]

                                                        sage: G = QuaternionGroup()
9     `C_9`                                                                                                                                           [9,1]
9     `C_3 \times C_3`                                                                                                                                [9,2]
10    `C_{10}`                                                                                                                                        [10,2]
10    `D_5`                                                                                                                                           [10,1]
11    `C_{11}`                                                                                                                                        [11,1]
12    `C_{12}`                                                                                                                                        [12,2]
12    `C_6 \times C_2`                                                                                                                                [12,5]
12    `D_6`                                                                                                                                           [12,4]
12    `A_4` （四字母交错群）                           ::                                                                                              [12,3]

                                                        sage: G = AlternatingGroup(4)
12    `Q_6` （12 阶双环群）                            ::                                                                                              [12,1]

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
16    `Q_{8}` （16 阶双环群）                          ::                                                                                              [16,9]

                                                        sage: G = DiCyclicGroup(4)
16    `2^4` 阶半二面体群                               ::                                                                                              [16,8]

                                                        sage: G = SemidihedralGroup(4)
16    `2^4` 阶分裂亚循环群                             ::                                                                                              [16,6]

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
20    `Q_{10}` （20 阶双环群）                                                                                                           [20,1]
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
24    `Q_{12}` （24 阶双环群）                                                                                                           [24,4]
24    `C_{12} \times C_2`                                                                                                                             [24,9]
24    `C_6 \times C_2 \times C_2`                                                                                                                     [24,15]
24    `S_4` （四字母对称群）                           ::                                                                                              [24,12]

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
27    `3^3` 阶分裂亚循环群                             ::                                                                                              [27,4]

                                                        sage: G = SplitMetacyclicGroup(3,3)
27    `(C_3 \times C_3) \rtimes_{\phi} C_3`           ::                                                                                              [27,3]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: A = direct_product_permgroups([C3,C3])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0]*A.gens()[1].inverse(),A.gens()[1]])
                                                        sage: phi = [[(1,2,3)],[alpha]]
28    `C_{28}`                                                                                                                                        [28,2]
28    `C_{14} \times C_2`                                                                                                                             [28,4]
28    `D_{14}`                                                                                                                                        [28,3]
28    `Q_{14}` （28 阶双环群）                                                                                                           [28,1]
29    `C_{29}`                                                                                                                                        [29,1]
30    `C_{30}`                                                                                                                                        [30,4]
30    `D_{15}`                                                                                                                                        [30,3]
30    `D_5 \times C_3`                                                                                                                                [30,2]
30    `D_3 \times C_5`                                                                                                                                [30,1]
31    `C_{31}`                                                                                                                                        [31,1]
===== =============================================== =============================================================================================== ===========================

该表由 Kevin Halasz 提供。


小于等于 15 阶的有限呈示群的构建说明
====================================

Sage 能够轻松构建阶数小于等于 15 的有限呈示群。
我们将首先探讨创建有限生成阿贝尔群，以及有限呈示群的直积和半直积。

所有有限生成阿贝尔群都可以使用 ``groups.presentation.FGAbelian(ls)`` 命令创建，
其中 ``ls`` 是一个非负整数列表，该列表被简化为定义要返回群的不变量。
例如，要构建 `C_4 \times C_2 \times C_2 \times C_2`，我们可以简单地使用::

    sage: A = groups.presentation.FGAbelian([4,2,2,2])


无论输入的整数列表如何，对于给定群的输出都是相同的。
以下示例为阶数为 30 的循环群产生相同的表示。

::

    sage: A = groups.presentation.FGAbelian([2,3,5])
    sage: B = groups.presentation.FGAbelian([30])

如果 ``G`` 和 ``H`` 是有限呈示群，我们可以使用以下代码来创建 ``G`` 和 ``H`` 的直积，`G \times H`。

.. skip

::

    sage: D = G.direct_product(H)

假设存在从群 `G` 到群 `H` 的自同构群的同态 `\phi`。
通过 `\phi` 将 `G` 与 `H` 的半直积定义为 `G` 和 `H` 的笛卡尔积，
运算为 `(g_1, h_1)(g_2, h_2) = (g_1 g_2, \phi_{h_1}(g_2) h_2)` 其中 `\phi_h = \phi(h)`。
要在 Sage 中为两个有限呈示群构造此乘积，我们必须使用一对列表手动定义 `\phi`。
第一个列表由群 `G` 的生成器组成，而第二个列表由第一个列表中相应生成器的像组成。
这些自同构同样定义为一对列表，一个列表为生成器，另一个列表为像。
作为示例，我们将阶数为 16 的二面体群构造为循环群的半直积。

::

    sage: C2 = groups.presentation.Cyclic(2)
    sage: C8 = groups.presentation.Cyclic(8)
    sage: hom = (C2.gens(), [ ([C8([1])], [C8([-1])]) ])
    sage: D = C2.semidirect_product(C8, hom)

下表显示了阶数小于等于 15 的群，以及如何在 Sage 中构造它们。重复命令已被省略，但通过以下示例进行了描述。

阶数为 `n` 的循环群可以通过单个命令创建：

.. skip

::

    sage: C = groups.presentation.Cyclic(n)

对于阶数为 `2n` 的二面体群也类似：

.. skip

::

    sage: D = groups.presentation.Dihedral(n)

该表是根据前面 Kevin Halasz 创建的表格构造的。


===== =============================================== =============================================================================================== ===========================
阶     群描述                                           命令                                                                                            GAP ID
===== =============================================== =============================================================================================== ===========================
1     平凡群                                           ::                                                                                              [1,1]

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

6     `S_3` （三字母对称群）                           ::                                                                                              [6,1]

                                                        sage: G = groups.presentation.Symmetric(3)

7     `C_7`                                                                                                                                           [7,1]
8     `C_8`                                                                                                                                           [8,1]

8     `C_4 \times C_2`                                ::                                                                                              [8,2]

                                                        sage: G = groups.presentation.FGAbelian([4,2])

8     `C_2\times C_2\times C_2`                       ::                                                                                              [8,5]

                                                        sage: G = groups.presentation.FGAbelian([2,2,2])

8     `D_4`                                           ::                                                                                              [8,3]

                                                        sage: G = groups.presentation.Dihedral(4)

8     四元群 (Q)                                       ::                                                                                              [8,4]

                                                        sage: G = groups.presentation.Quaternion()

9     `C_9`                                                                                                                                           [9,1]
9     `C_3 \times C_3`                                                                                                                                [9,2]
10    `C_{10}`                                                                                                                                        [10,2]
10    `D_5`                                                                                                                                           [10,1]
11    `C_{11}`                                                                                                                                        [11,1]
12    `C_{12}`                                                                                                                                        [12,2]
12    `C_6 \times C_2`                                                                                                                                [12,5]
12    `D_6`                                                                                                                                           [12,4]
12    `A_4` （四字母交错群）                           ::                                                                                              [12,3]

                                                        sage: G = groups.presentation.Alternating(4)

12    `Q_6` （12 阶双环群）                            ::                                                                                              [12,1]

                                                        sage: G = groups.presentation.DiCyclic(3)

13    `C_{13}`                                                                                                                                        [13,1]
14    `C_{14}`                                                                                                                                        [14,2]
14    `D_{7}`                                                                                                                                         [14,1]
15    `C_{15}`                                                                                                                                        [15,1]
===== =============================================== =============================================================================================== ===========================

