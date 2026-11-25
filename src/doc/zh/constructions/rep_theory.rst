******
表示论
******

.. index:
   pair: ordinary representation; character

.. _section-character:

普通特征标 (Ordinary characters)
================================

如何在 Sage 中计算有限群的特征标表？可以使用 Sage-GAP 接口来计算特征标表。

你可以使用 ``PermutationGroup`` 类的 ``character_table`` 方法，
或通过 GAP 命令 ``CharacterTable`` 的接口，将置换群 `G` 的特征标值表构建为 Sage 矩阵。

::

    sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
    sage: G.order()
    8
    sage: G.character_table() # random
    [ 1  1  1  1  1]
    [ 1 -1 -1  1  1]
    [ 1 -1  1 -1  1]
    [ 1  1 -1 -1  1]
    [ 2  0  0  0 -2]
    sage: CT = libgap(G).CharacterTable()
    sage: CT.Display() # random
    CT1
    <BLANKLINE>
     2  3  2  2  2  3
    <BLANKLINE>
       1a 2a 2b 4a 2c
    2P 1a 1a 1a 2c 1a
    3P 1a 2a 2b 4a 2c
    <BLANKLINE>
    X.1     1  1  1  1  1
    X.2     1 -1 -1  1  1
    X.3     1 -1  1 -1  1
    X.4     1  1 -1 -1  1
    X.5     2  .  .  . -2

下面是另一个示例:

::

    sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3)]])
    sage: G.character_table() # random
    [         1          1          1          1]
    [         1 -zeta3 - 1      zeta3          1]
    [         1      zeta3 -zeta3 - 1          1]
    [         3          0          0         -1]
    sage: G = libgap.eval("Group((1,2)(3,4),(1,2,3))"); G
    Group([ (1,2)(3,4), (1,2,3) ])
    sage: T = G.CharacterTable()
    sage: T.Display() # random
    CT2
    <BLANKLINE>
         2  2  .  .  2
         3  1  1  1  .
    <BLANKLINE>
           1a 3a 3b 2a
        2P 1a 3b 3a 1a
        3P 1a 1a 1a 2a
    <BLANKLINE>
    X.1     1  1  1  1
    X.2     1  A /A  1
    X.3     1 /A  A  1
    X.4     3  .  . -1
    <BLANKLINE>
    A = E(3)^2
      = (-1-Sqrt(-3))/2 = -1-b3

其中 `E(3)` 表示单位立方根，`ER(-3)` 表示 `-3` 的平方根，即 `i\sqrt{3}`,
而 `b3 = \frac{1}{2}(-1+i \sqrt{3})`。
请注意添加的 ``print`` Python 命令。这会令输出更美观。

.. link

::

    sage: irr = G.Irr(); sorted(irr)
    [Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, 1, 1, 1 ] ),
     Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, E(3)^2, E(3), 1 ] ),
     Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, E(3), E(3)^2, 1 ] ),
     Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 3, 0, 0, -1 ] )]
    sage: irr.Display() # random
    [ [       1,       1,       1,       1 ],
      [       1,  E(3)^2,    E(3),       1 ],
      [       1,    E(3),  E(3)^2,       1 ],
      [       3,       0,       0,      -1 ] ]
    sage: CG = G.ConjugacyClasses(); CG
    [ ()^G, (2,3,4)^G, (2,4,3)^G, (1,2)(3,4)^G ]
    sage: gamma = CG[2]; gamma
    (2,4,3)^G
    sage: g = gamma.Representative(); g
    (2,4,3)
    sage: chi = irr[1]; chi # random
    Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, E(3)^2, E(3), 1 ] )
    sage: g^chi # random
    E(3)

最后一个量是特征 ``chi`` 在群元素 ``g`` 处的值。

或者，如果你关闭 IPython 的“美观打印”，那么表格将打印得更好。

.. skip

::

    sage: %Pprint
    Pretty printing has been turned OFF
    sage: G = libgap.eval("Group((1,2)(3,4),(1,2,3))"); G
    Group([ (1,2)(3,4), (1,2,3) ])
    sage: T = G.CharacterTable(); T
    CharacterTable( Alt( [ 1 .. 4 ] ) )
    sage: T.Display()
    CT3
    <BLANKLINE>
         2  2  2  .  .
         3  1  .  1  1
    <BLANKLINE>
           1a 2a 3a 3b
        2P 1a 1a 3b 3a
        3P 1a 2a 1a 1a
    <BLANKLINE>
    X.1     1  1  1  1
    X.2     1  1  A /A
    X.3     1  1 /A  A
    X.4     3 -1  .  .
    <BLANKLINE>
    A = E(3)^2
      = (-1-Sqrt(-3))/2 = -1-b3
    sage: irr = G.Irr(); irr
    [ Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, 1, 1, 1 ] ),
      Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, 1, E(3)^2, E(3) ] ),
      Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, 1, E(3), E(3)^2 ] ),
      Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 3, -1, 0, 0 ] ) ]
    sage: irr.Display()
    [ [       1,       1,       1,       1 ],
      [       1,       1,  E(3)^2,    E(3) ],
      [       1,       1,    E(3),  E(3)^2 ],
      [       3,      -1,       0,       0 ] ]
    sage: %Pprint
    Pretty printing has been turned ON

.. index::
   pair: modular representation; character
   pair: character; Brauer

.. _section-brauer:

布劳尔特征标 (Brauer characters)
================================

GAP 中的布劳尔特征标表尚未具有“原生”接口。
要访问它们，你可以使用 ``libgap.eval`` 命令直接与 GAP 交互。

下面的示例通过使用 GAP 接口来说明语法。

::

    sage: G = libgap.eval("Group((1,2)(3,4),(1,2,3))"); G
    Group([ (1,2)(3,4), (1,2,3) ])
    sage: irr = G.IrreducibleRepresentations(GF(7)); irr   # random arch. dependent output
    [ [ (1,2)(3,4), (1,2,3) ] -> [ [ [ Z(7)^0 ] ], [ [ Z(7)^4 ] ] ],
      [ (1,2)(3,4), (1,2,3) ] -> [ [ [ Z(7)^0 ] ], [ [ Z(7)^2 ] ] ],
      [ (1,2)(3,4), (1,2,3) ] -> [ [ [ Z(7)^0 ] ], [ [ Z(7)^0 ] ] ],
      [ (1,2)(3,4), (1,2,3) ] ->
        [ [ [ Z(7)^2, Z(7)^5, Z(7) ], [ Z(7)^3, Z(7)^2, Z(7)^3 ],
            [ Z(7), Z(7)^5, Z(7)^2 ] ],
          [ [ 0*Z(7), Z(7)^0, 0*Z(7) ], [ 0*Z(7), 0*Z(7), Z(7)^0 ],
            [ Z(7)^0, 0*Z(7), 0*Z(7) ] ] ] ]
    sage: brvals = [[chi.Image(c.Representative()).BrauerCharacterValue()
    ....:            for c in G.ConjugacyClasses()] for chi in irr]
    sage: brvals         # random architecture dependent output
    [ [       1,       1,  E(3)^2,    E(3) ],
      [       1,       1,    E(3),  E(3)^2 ],
      [       1,       1,       1,       1 ],
      [       3,      -1,       0,       0 ] ]
    sage: T = G.CharacterTable()
    sage: T.Display() # random
    CT3
    <BLANKLINE>
         2  2  .  .  2
         3  1  1  1  .
    <BLANKLINE>
           1a 3a 3b 2a
        2P 1a 3b 3a 1a
        3P 1a 1a 1a 2a
    <BLANKLINE>
    X.1     1  1  1  1
    X.2     1  A /A  1
    X.3     1 /A  A  1
    X.4     3  .  . -1
    <BLANKLINE>
    A = E(3)^2
      = (-1-Sqrt(-3))/2 = -1-b3
