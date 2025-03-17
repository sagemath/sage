*********************
Representation theory
*********************

.. index:
   pair: ordinary representation; character

.. _section-character:

Ordinary characters
===================

How can you compute character tables of a finite group in Sage? The
Sage-GAP interface can be used to compute character tables.

You can construct the table of character values of a permutation
group :math:`G` as a Sage matrix, using the method
``character_table`` of the PermutationGroup class, or via the
interface to the GAP command ``CharacterTable``.

::

    sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
    sage: G.order()
    8
    sage: G.character_table()
    [ 1  1  1  1  1]
    [ 1 -1 -1  1  1]
    [ 1 -1  1 -1  1]
    [ 1  1 -1 -1  1]
    [ 2  0  0  0 -2]
    sage: CT = libgap(G).CharacterTable()
    sage: CT.Display()
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

Here is another example:

::

    sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3)]])
    sage: G.character_table()
    [         1          1          1          1]
    [         1 -zeta3 - 1      zeta3          1]
    [         1      zeta3 -zeta3 - 1          1]
    [         3          0          0         -1]
    sage: G = libgap.eval("Group((1,2)(3,4),(1,2,3))"); G
    Group([ (1,2)(3,4), (1,2,3) ])
    sage: T = G.CharacterTable()
    sage: T.Display()
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

where :math:`E(3)` denotes a cube root of unity, :math:`ER(-3)`
denotes a square root of :math:`-3`, say :math:`i\sqrt{3}`, and
:math:`b3 = \frac{1}{2}(-1+i \sqrt{3})`. Note the added ``print``
Python command. This makes the output look much nicer.

.. link

::

    sage: irr = G.Irr(); irr
    [ Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, 1, 1, 1 ] ), 
      Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, E(3)^2, E(3), 1 ] ), 
      Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, E(3), E(3)^2, 1 ] ), 
      Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 3, 0, 0, -1 ] ) ]
    sage: irr.Display()
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
    sage: chi = irr[1]; chi
    Character( CharacterTable( Alt( [ 1 .. 4 ] ) ), [ 1, E(3)^2, E(3), 1 ] )
    sage: g^chi
    E(3)

This last quantity is the value of the character ``chi`` at the group
element ``g``.

Alternatively, if you turn IPython "pretty printing" off, then the
table prints nicely.

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

Brauer characters
=================

The Brauer character tables in GAP do not yet have a "native"
interface. To access them you can directly interface with GAP using
the ``libgap.eval`` command.

The example below using the GAP interface illustrates the syntax.

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
    sage: T.Display()
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
