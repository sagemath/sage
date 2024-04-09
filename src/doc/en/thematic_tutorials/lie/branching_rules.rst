.. linkall

.. _branch_rules:

-------------------------------------
Maximal Subgroups and Branching Rules
-------------------------------------

Branching rules
---------------

If `G` is a Lie group and `H` is a subgroup, one often needs to know
how representations of `G` restrict to `H`. Irreducibles usually do
not restrict to irreducibles. In some cases the restriction is regular
and predictable, in other cases it is chaotic. In some cases it
follows a rule that can be described combinatorially, but the
combinatorial description is subtle. The description of how
irreducibles decompose into irreducibles is called a *branching rule*.

References for this topic:

- [FauserEtAl2006]_

- [King1975]_

- [HoweEtAl2005]_

- [McKayPatera1981]_

Sage can compute how a character of `G` restricts to `H`. It does
so not by memorizing a combinatorial rule, but by computing the
character and restricting the character to a maximal torus of `H`.
What Sage has memorized (in a series of built-in encoded rules)
are the various embeddings of maximal tori of maximal subgroups of `G`.
The maximal subgroups of Lie groups were determined in [Dynkin1952]_.
This approach to computing branching rules has a limitation: the
character must fit into memory and be computable by Sage's
internal code in real time.

It is sufficient to consider the case where `H` is a maximal
subgroup of `G`, since if this is known then one may branch down
successively through a series of subgroups, each maximal in its
predecessors. A problem is therefore to understand the maximal
subgroups in a Lie group, and to give branching rules for each,
and a goal of this tutorial is to explain the embeddings of
maximal subgroups.

Sage has a class ``BranchingRule`` for branching rules. The function
``branching_rule`` returns elements of this class. For example,
the natural embedding of `Sp(4)` into `SL(4)` corresponds to
the branching rule that we may create as follows::

    sage: b = branching_rule("A3","C2",rule="symmetric"); b
    symmetric branching rule A3 => C2

The name "symmetric" of this branching rule will be
explained further later, but it means that `Sp(4)` is
the fixed subgroup of an involution of `Sl(4)`.
Here ``A3`` and ``C2`` are the Cartan types of the groups
`G=SL(4)` and `H=Sp(4)`.

Now we may see how representations of `SL(4)` decompose
into irreducibles when they are restricted to `Sp(4)`::

    sage: A3 = WeylCharacterRing("A3",style="coroots")
    sage: chi = A3(1,0,1); chi.degree()
    15
    sage: C2 = WeylCharacterRing("C2",style="coroots")
    sage: chi.branch(C2,rule=b)
    C2(0,1) + C2(2,0)

Alternatively, we may pass ``chi`` to ``b`` as an
argument of its branch method, which gives the same
result::

    sage: b.branch(chi)
    C2(0,1) + C2(2,0)

It is believed that the built-in branching rules of
Sage are sufficient to handle all maximal subgroups
and this is certainly the case when the rank if
less than or equal to 8.

However, if you want to branch to a subgroup that
is not maximal you may not find a built-in 
branching rule. We may compose branching rules to build
up embeddings. For example, here are two different
embeddings of `Sp(4)` with Cartan type ``C2`` in
`Sp(8)`, with Cartan type ``C4``. One embedding
factors through `Sp(4)\times Sp(4)`, while the
other factors through `SL(4)`. To check that the embeddings
are not conjugate, we branch a (randomly chosen) representation.
Observe that we do not have to build the intermediate
:class:`Weyl character rings <WeylCharacterRing>`.

::

    sage: C4 = WeylCharacterRing("C4",style="coroots")
    sage: b1 = branching_rule("C4","A3","levi")*branching_rule("A3","C2","symmetric"); b1
    composite branching rule C4 => (levi) A3 => (symmetric) C2
    sage: b2 = branching_rule("C4","C2xC2","orthogonal_sum")*branching_rule("C2xC2","C2","proj1"); b2
    composite branching rule C4 => (orthogonal_sum) C2xC2 => (proj1) C2
    sage: C2 = WeylCharacterRing("C2",style="coroots")
    sage: C4 = WeylCharacterRing("C4",style="coroots")
    sage: [C4(2,0,0,1).branch(C2, rule=br) for br in [b1,b2]]
    [4*C2(0,0) + 7*C2(0,1) + 15*C2(2,0) + 7*C2(0,2) + 11*C2(2,1) + C2(0,3) + 6*C2(4,0) + 3*C2(2,2),
     10*C2(0,0) + 40*C2(1,0) + 50*C2(0,1) + 16*C2(2,0) + 20*C2(1,1) + 4*C2(3,0) + 5*C2(2,1)]


What's in a branching rule?
---------------------------

The essence of the branching rule is a function from the
weight lattice of `G` to the weight lattice of the subgroup `H`,
usually implemented as a function on the ambient vector
spaces. Indeed, we may conjugate the embedding so that a
Cartan subalgebra `U` of `H` is contained in a Cartan subalgebra
`T` of `G`. Since the ambient vector space of the weight
lattice of `G` is `\hbox{Lie}(T)^*`, we get map
`\hbox{Lie}(T)^*\to\hbox{Lie}(U)^*`, and this must be
implemented as a function. For speed, the function usually
just returns a list, which can be coerced into `\hbox{Lie}(U)^*`.

::

    sage: b = branching_rule("A3","C2","symmetric")
    sage: for r in RootSystem("A3").ambient_space().simple_roots():
    ....:     print("{} {}".format(r, b(r)))
    (1, -1, 0, 0) [1, -1]
    (0, 1, -1, 0) [0, 2]
    (0, 0, 1, -1) [1, -1]

We could conjugate this map by an element of the Weyl
group of `G`, and the resulting map would give the same
decomposition of representations of `G` into irreducibles
of `H`. However it is a good idea to choose the map so
that it takes dominant weights to dominant weights, and,
insofar as possible, simple roots of `G` into
simple roots of `H`. This includes sometimes the affine root `\alpha_0`
of `G`, which we recall is the negative of the highest root.

The branching rule has a ``describe()`` method that shows how
the roots (including the affine root) restrict. This is a
useful way of understanding the embedding. You might
want to try it with various branching rules of different
kinds, ``"extended"``, ``"symmetric"``, ``"levi"`` etc.

::

    sage: b.describe()
    <BLANKLINE>
    0
    O-------+
    |       |
    |       |
    O---O---O
    1   2   3
    A3~
    <BLANKLINE>
    root restrictions A3 => C2:
    <BLANKLINE>
    O=<=O
    1   2
    C2
    <BLANKLINE>
    1 => 1
    2 => 2
    3 => 1
    <BLANKLINE>
    For more detailed information use verbose=True
    
The extended Dynkin diagram of `G` and the ordinary
Dynkin diagram of `H` are shown for reference, and
``3 => 1`` means that the third simple root `\alpha_3`
of `G` restricts to the first simple root of `H`.
In this example, the affine root does not restrict to
a simple roots, so it is omitted from the list of
restrictions. If you add the parameter ``verbose=true`` you will
be shown the restriction of all simple roots and the
affine root, and also the restrictions of the fundamental weights
(in coroot notation).

Maximal subgroups
-----------------

Sage has a database of maximal subgroups for every simple Cartan
type of rank `\le 8`. You may access this with the
``maximal_subgroups`` method of the Weyl character ring::

    sage: E7 = WeylCharacterRing("E7",style="coroots")
    sage: E7.maximal_subgroups()
    A7:branching_rule("E7","A7","extended")
    E6:branching_rule("E7","E6","levi")
    A2:branching_rule("E7","A2","miscellaneous")
    A1:branching_rule("E7","A1","iii")
    A1:branching_rule("E7","A1","iv")
    A1xF4:branching_rule("E7","A1xF4","miscellaneous")
    G2xC3:branching_rule("E7","G2xC3","miscellaneous")
    A1xG2:branching_rule("E7","A1xG2","miscellaneous")
    A1xA1:branching_rule("E7","A1xA1","miscellaneous")
    A1xD6:branching_rule("E7","A1xD6","extended")
    A5xA2:branching_rule("E7","A5xA2","extended")

It should be understood that there are other ways of
embedding `A_2=\hbox{SL}(3)` into the Lie group `E_7`,
but only one way as a maximal subgroup. On the other
hand, there are but only one way to embed it as a
maximal subgroup. The embedding will be explained below.
You may obtain the branching rule as follows, and use it to determine
the decomposition of irreducible representations of `E_7`
as follows::

    sage: b = E7.maximal_subgroup("A2"); b
    miscellaneous branching rule E7 => A2
    sage: E7, A2 = [WeylCharacterRing(x,style="coroots") for x in ["E7","A2"]]
    sage: E7(1,0,0,0,0,0,0).branch(A2,rule=b)
    A2(1,1) + A2(4,4)

This gives the same branching rule as just pasting line beginning
to the right of the colon onto the command line::

    sage:branching_rule("E7","A2","miscellaneous")
    miscellaneous branching rule E7 => A2

There are two distinct embeddings of `A_1=\hbox{SL}(2)` into
`E_7` as maximal subgroups, so the ``maximal_subgroup``
method will return a list of rules::

     sage: WeylCharacterRing("E7").maximal_subgroup("A1")
     [iii branching rule E7 => A1, iv branching rule E7 => A1]

The list of maximal subgroups returned by the ``maximal_subgroups``
method for irreducible Cartan types of rank up to 8 is believed to
be complete up to outer automorphisms. You may want a list that is
complete up to inner automorphisms.  For example, `E_6` has a
nontrivial Dynkin diagram automorphism so it has an outer
automorphism that is not inner::

    sage: E6, A2xG2 = [WeylCharacterRing(x,style="coroots") for x in ["E6","A2xG2"]]
    sage: b = E6.maximal_subgroup("A2xG2"); b
    miscellaneous branching rule E6 => A2xG2
    sage: E6(1,0,0,0,0,0).branch(A2xG2,rule=b)
    A2xG2(0,1,1,0) + A2xG2(2,0,0,0)
    sage: E6(0,0,0,0,0,1).branch(A2xG2,rule=b)
    A2xG2(1,0,1,0) + A2xG2(0,2,0,0)
    
Since as we see the two 27 dimensional irreducibles (which are
interchanged by the outer automorphism) have different branching,
the `A_2\times G_2` subgroup is changed to a different one
by the outer automorphism. To obtain the second branching
rule, we compose the given one with this automorphism::

    sage: b1 = branching_rule("E6","E6","automorphic")*b; b1
    composite branching rule E6 => (automorphic) E6 => (miscellaneous) A2xG2

.. _levi_branch_rules:

Levi subgroups
--------------

A Levi subgroup may or may not be maximal. They are easily
classified. If one starts with a Dynkin diagram for `G` and removes a
single node, one obtains a smaller Dynkin diagram, which is the Dynkin
diagram of a smaller subgroup `H`.

For example, here is the A3 Dynkin diagram:

::

    sage: A3 = WeylCharacterRing("A3")
    sage: A3.dynkin_diagram()
    O---O---O
    1   2   3
    A3

We see that we may remove the node 3 and obtain `A_2`, or the node 2
and obtain `A_1 \times A_1`. These correspond to the Levi subgroups
`GL(3)` and `GL(2) \times GL(2)` of `GL(4)`.

Let us construct the irreducible
representations of `GL(4)` and branch them down to these down to
`GL(3)` and `GL(2) \times GL(2)`::

    sage: reps = [A3(v) for v in A3.fundamental_weights()]; reps
    [A3(1,0,0,0), A3(1,1,0,0), A3(1,1,1,0)]
    sage: A2 = WeylCharacterRing("A2")
    sage: A1xA1 = WeylCharacterRing("A1xA1")
    sage: [pi.branch(A2, rule="levi") for pi in reps]
    [A2(0,0,0) + A2(1,0,0), A2(1,0,0) + A2(1,1,0), A2(1,1,0) + A2(1,1,1)]
    sage: [pi.branch(A1xA1, rule="levi") for pi in reps]
    [A1xA1(1,0,0,0) + A1xA1(0,0,1,0),
     A1xA1(1,1,0,0) + A1xA1(1,0,1,0) + A1xA1(0,0,1,1),
     A1xA1(1,1,1,0) + A1xA1(1,0,1,1)]

Let us redo this calculation in coroot notation. As we have explained,
coroot notation does not distinguish between representations of
`GL(4)` that have the same restriction to `SL(4)`, so in effect we are
now working with the groups `SL(4)` and its Levi subgroups `SL(3)` and
`SL(2) \times SL(2)`, which is the derived group of its Levi subgroup::

    sage: A3 = WeylCharacterRing("A3", style="coroots")
    sage: reps = [A3(v) for v in A3.fundamental_weights()]; reps
    [A3(1,0,0), A3(0,1,0), A3(0,0,1)]
    sage: A2 = WeylCharacterRing("A2", style="coroots")
    sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
    sage: [pi.branch(A2, rule="levi") for pi in reps]
    [A2(0,0) + A2(1,0), A2(0,1) + A2(1,0), A2(0,0) + A2(0,1)]
    sage: [pi.branch(A1xA1, rule="levi") for pi in reps]
    [A1xA1(1,0) + A1xA1(0,1), 2*A1xA1(0,0) + A1xA1(1,1), A1xA1(1,0) + A1xA1(0,1)]

Now we may observe a distinction difference in branching from

.. MATH::

    GL(4) \to GL(2) \times GL(2)

versus

.. MATH::

    SL(4) \to SL(2) \times SL(2).

Consider the representation ``A3(0,1,0)``, which is the six dimensional
exterior square. In the coroot notation, the restriction contained two
copies of the trivial representation, ``2*A1xA1(0,0)``. The other way,
we had instead three distinct representations in the restriction, namely
``A1xA1(1,1,0,0)`` and ``A1xA1(0,0,1,1)``, that is,
`\det \otimes 1` and `1 \otimes \det`.

The Levi subgroup ``A1xA1`` is actually not maximal. Indeed, we may
factor the embedding:

.. MATH::

    SL(2) \times SL(2) \to Sp(4) \to SL(4).

Therefore there are branching rules ``A3 -> C2`` and ``C2 -> A2``, and
we could accomplish the branching in two steps, thus::

    sage: A3 = WeylCharacterRing("A3", style="coroots")
    sage: C2 = WeylCharacterRing("C2", style="coroots")
    sage: B2 = WeylCharacterRing("B2", style="coroots")
    sage: D2 = WeylCharacterRing("D2", style="coroots")
    sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
    sage: reps = [A3(fw) for fw in A3.fundamental_weights()]
    sage: [pi.branch(C2, rule="symmetric").branch(B2, rule="isomorphic"). \
             branch(D2, rule="extended").branch(A1xA1, rule="isomorphic") for pi in reps]
    [A1xA1(1,0) + A1xA1(0,1), 2*A1xA1(0,0) + A1xA1(1,1), A1xA1(1,0) + A1xA1(0,1)]

As you can see, we've redone the branching rather circuitously this
way, making use of the branching rules ``A3 -> C2`` and ``B2 -> D2``, and
two accidental isomorphisms ``C2 = B2`` and ``D2 = A1xA1``. It is much
easier to go in one step using ``rule="levi"``, but reassuring that we
get the same answer!

Subgroups classified by the extended Dynkin diagram
---------------------------------------------------

It is also true that if we remove one node from the extended Dynkin
diagram that we obtain the Dynkin diagram of a subgroup. For example::

    sage: G2 = WeylCharacterRing("G2", style="coroots")
    sage: G2.extended_dynkin_diagram()
      3
    O=<=O---O
    1   2   0
    G2~

Observe that by removing the 1 node that we obtain an `A_2` Dynkin
diagram. Therefore the exceptional group `G_2` contains a copy of
`SL(3)`. We branch the two representations of `G_2` corresponding to the
fundamental weights to this copy of `A_2`::

    sage: G2 = WeylCharacterRing("G2", style="coroots")
    sage: A2 = WeylCharacterRing("A2", style="coroots")
    sage: [G2(f).degree() for f in G2.fundamental_weights()]
    [7, 14]
    sage: [G2(f).branch(A2, rule="extended") for f in G2.fundamental_weights()]
    [A2(0,0) + A2(0,1) + A2(1,0), A2(0,1) + A2(1,0) + A2(1,1)]

The two representations of `G_2`, of degrees 7 and 14 respectively, are
the action on the octonions of trace zero and the adjoint
representation.

For embeddings of this type, the rank of the subgroup `H` is the same
as the rank of `G`. This is in contrast with embeddings of Levi type,
where `H` has rank one less than `G`.

Levi subgroups of `G_2`
-----------------------

The exceptional group `G_2` has two Levi subgroups of type
`A_1`. Neither is maximal, as we can see from the extended
Dynkin diagram: the subgroups `A_1\times A_1` and `A_2`
are maximal and each contains a Levi subgroup. (Actually
`A_1\times A_1` contains a conjugate of both.) Only
the Levi subgroup containing the short root is implemented
as an instance of ``rule="levi"``. To obtain the other,
use the rule::

    sage: branching_rule("G2","A2","extended")*branching_rule("A2","A1","levi")
    composite branching rule G2 => (extended) A2 => (levi) A1

which branches to the `A_1` Levi subgroup containing a long root.

Orthogonal and symplectic subgroups of orthogonal and symplectic groups
-----------------------------------------------------------------------

If `G = \hbox{SO}(n)` then `G` has a subgroup `\hbox{SO}(n-1)`. Depending on
whether `n` is even or odd, we thus have branching rules
``['D',r]`` to ``['B',r-1]`` or ``['B',r]`` to ``['D',r]``. These are
handled as follows::

     sage: branching_rule("B4","D4",rule="extended")
     extended branching rule B4 => D4
     sage: branching_rule("D4","B3",rule="symmetric")
     symmetric branching rule D4 => B3

If `G = \hbox{SO}(r+s)` then `G` has a subgroup
`\hbox{SO}(r) \times \hbox{SO}(s)`. This
lifts to an embedding of the universal covering groups

.. MATH::

    \hbox{spin}(r) \times \hbox{spin}(s) \to \hbox{spin}(r+s).

Sometimes this embedding is of extended type, and sometimes it is
not. It is of extended type unless `r` and `s` are both odd. If it is
of extended type then you may use ``rule="extended"``. In any case you
may use ``rule="orthogonal_sum"``. The name refer to the origin of the
embedding `SO(r) \times SO(s) \to SO(r+s)` from the decomposition of
the underlying quadratic space as a direct sum of two orthogonal
subspaces.

There are four cases depending on the parity of `r` and `s`. For
example, if `r = 2k` and `s = 2l` we have an embedding::

    ['D',k] x ['D',l] --> ['D',k+l]

This is of extended type. Thus consider the embedding
``D4xD3 -> D7``. Here is the extended Dynkin diagram::

      0 O           O 7
        |           |
        |           |
    O---O---O---O---O---O
    1   2   3   4   5   6

Removing the 4 vertex results in a disconnected Dynkin diagram::

      0 O           O 7
        |           |
        |           |
    O---O---O       O---O
    1   2   3       5   6

This is ``D4xD3``.  Therefore use the "extended" branching rule:

::

    sage: D7 = WeylCharacterRing("D7", style="coroots")
    sage: D4xD3 = WeylCharacterRing("D4xD3", style="coroots")
    sage: spin = D7(D7.fundamental_weights()[7]); spin
    D7(0,0,0,0,0,0,1)
    sage: spin.branch(D4xD3, rule="extended")
    D4xD3(0,0,1,0,0,1,0) + D4xD3(0,0,0,1,0,0,1)

But we could equally well use the "orthogonal_sum" rule::

    sage: spin.branch(D4xD3, rule="orthogonal_sum")
    D4xD3(0,0,1,0,0,1,0) + D4xD3(0,0,0,1,0,0,1)

Similarly we have embeddings::

    ['D',k] x ['B',l] --> ['B',k+l]

These are also of extended type. For example consider the embedding of
``D3xB2 -> B5``. Here is the ``B5`` extended Dynkin diagram::

        O 0
        |
        |
    O---O---O---O=>=O
    1   2   3   4   5

Removing the 3 node gives::

        O 0
        |
    O---O       O=>=O
    1   2       4   5

and this is the Dynkin diagram or ``D3xB2``. For such branchings we
again use either ``rule="extended"`` or ``rule="orthogonal_sum"``.

Finally, there is an embedding ::

    ['B',k] x ['B',l] --> ['D',k+l+1]

This is *not* of extended type, so you may not use ``rule="extended"``.
You *must* use ``rule="orthogonal_sum"``::

    sage: D5 = WeylCharacterRing("D5",style="coroots")
    sage: B2xB2 = WeylCharacterRing("B2xB2",style="coroots")
    sage: [D5(v).branch(B2xB2,rule="orthogonal_sum") for v in D5.fundamental_weights()]
    [B2xB2(1,0,0,0) + B2xB2(0,0,1,0),
     B2xB2(0,2,0,0) + B2xB2(1,0,1,0) + B2xB2(0,0,0,2),
     B2xB2(0,2,0,0) + B2xB2(0,2,1,0) + B2xB2(1,0,0,2) + B2xB2(0,0,0,2),
     B2xB2(0,1,0,1), B2xB2(0,1,0,1)]

Non-maximal Levi subgroups and Projection from Reducible Types
--------------------------------------------------------------

Not all Levi subgroups are maximal. Recall that the Dynkin-diagram
of a Levi subgroup `H` of `G` is obtained by removing a node
from the Dynkin diagram of `G`. Removing the same node from
the extended Dynkin diagram of `G` results in the Dynkin
diagram of a subgroup of `G` that is strictly larger than
`H`. However this subgroup may or may not be proper, so the
Levi subgroup may or may not be maximal.

If the Levi subgroup is not maximal, the branching rule
may or may not be implemented in Sage. However if it is
not implemented, it may be constructed as a composition
of two branching rules.

For example, prior to Sage-6.1 ``branching_rule("E6","A5","levi")``
returned a not-implemented error and the advice to branch to
``A5xA1``. And we can see from the extended Dynkin diagram of `E_6`
that indeed `A_5` is not a maximal subgroup, since removing node 2
from the extended Dynkin diagram (see below) gives ``A5xA1``. To
construct the branching rule to `A_5` we may proceed as follows::

    sage: b = branching_rule("E6","A5xA1","extended")*branching_rule("A5xA1","A5","proj1"); b
    composite branching rule E6 => (extended) A5xA1 => (proj1) A5
    sage: E6 = WeylCharacterRing("E6",style="coroots")
    sage: A5 = WeylCharacterRing("A5",style="coroots")
    sage: E6(0,1,0,0,0,0).branch(A5,rule=b)
    3*A5(0,0,0,0,0) + 2*A5(0,0,1,0,0) + A5(1,0,0,0,1)
    sage: b.describe()
    <BLANKLINE>
            O 0
            |
            |
            O 2
            |
            |
    O---O---O---O---O
    1   3   4   5   6
    E6~
    root restrictions E6 => A5:
    <BLANKLINE>
    O---O---O---O---O
    1   2   3   4   5
    A5
    <BLANKLINE>
    0 => (zero)
    1 => 1
    3 => 2
    4 => 3
    5 => 4
    6 => 5
    <BLANKLINE>
    For more detailed information use verbose=True

Note that it is not necessary to construct the Weyl character ring
for the intermediate group ``A5xA1``.

This last example illustrates another common problem:
how to extract one component from a reducible root system.
We used the rule ``"proj1"`` to extract the first component.
We could similarly use ``"proj2"`` to get the second, or
more generally any combination of components::

    sage: branching_rule("A2xB2xG2","A2xG2","proj13")
    proj13 branching rule A2xB2xG2 => A2xG2

Symmetric subgroups
-------------------

If `G` admits an outer automorphism (usually of order two) then we may
try to find the branching rule to the fixed subgroup `H`. It can be
arranged that this automorphism maps the maximal torus `T` to itself
and that a maximal torus `U` of `H` is contained in `T`.

Suppose that the Dynkin diagram of `G` admits an automorphism. Then
`G` itself admits an outer automorphism. The Dynkin diagram of the
group `H` of invariants may be obtained by "folding" the Dynkin
diagram of `G` along the automorphism. The exception is the branching
rule `GL(2r) \to SO(2r)`.

Here are the branching rules that can be obtained using
``rule="symmetric"``.

+------------+-------------+---------------------------+
| `G`        | `H`         | Cartan Types              |
+============+=============+===========================+
| `GL(2r)`   | `Sp(2r)`    | ``['A',2r-1] => ['C',r]`` |
+------------+-------------+---------------------------+
| `GL(2r+1)` | `SO(2r+1)`  | ``['A',2r] => ['B',r]``   |
+------------+-------------+---------------------------+
| `GL(2r)`   | `SO(2r)`    | ``['A',2r-1] => ['D',r]`` |
+------------+-------------+---------------------------+
| `SO(2r)`   | `SO(2r-1)`  | ``['D',r] => ['B',r-1]``  |
+------------+-------------+---------------------------+
| `E_6`      | `F_4`       | ``['E',6] => ['F',4]``    |
+------------+-------------+---------------------------+


Tensor products
---------------

If `G_1` and `G_2` are Lie groups, and we have representations
`\pi_1: G_1 \to GL(n)` and `\pi_2: G_2 \to GL(m)` then the tensor
product is a representation of `G_1 \times G_2`. It has its image
in `GL(nm)` but sometimes this is conjugate to a subgroup of `SO(nm)`
or `Sp(nm)`. In particular we have the following cases.

+-------------------+---------------------------+------------------------------------------+
| Group             | Subgroup                  | Cartan Types                             |
+===================+===========================+==========================================+
| `GL(rs)`          | `GL(r)\times GL(s)`       | ``['A', rs-1] => ['A',r-1] x ['A',s-1]`` |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs+2r+2s+1)` | `SO(2r+1)\times SO(2s+1)` | ``['B',2rs+r+s] => ['B',r] x ['B',s]``   |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs+2s)`      | `SO(2r+1)\times SO(2s)`   | ``['D',2rs+s] => ['B',r] x ['D',s]``     |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs)`         | `SO(2r)\times SO(2s)`     | ``['D',2rs] => ['D',r] x ['D',s]``       |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs)`         | `Sp(2r)\times Sp(2s)`     | ``['D',2rs] => ['C',r] x ['C',s]``       |
+-------------------+---------------------------+------------------------------------------+
| `Sp(4rs+2s)`      | `SO(2r+1)\times Sp(2s)`   | ``['C',2rs+s] => ['B',r] x ['C',s]``     |
+-------------------+---------------------------+------------------------------------------+
| `Sp(4rs)`         | `Sp(2r)\times SO(2s)`     | ``['C',2rs] => ['C',r] x ['D',s]``       |
+-------------------+---------------------------+------------------------------------------+

These branching rules are obtained using ``rule="tensor"``.


Symmetric powers
----------------

The `k`-th symmetric and exterior power homomorphisms map
`GL(n) \to GL \left(\binom{n+k-1}{k} \right)` and
`GL \left(\binom{n}{k} \right)`. The corresponding branching rules
are not implemented but a special case is. The `k`-th symmetric power
homomorphism `SL(2) \to GL(k+1)` has its image inside of `SO(2r+1)` if
`k = 2r` and inside of `Sp(2r)` if `k = 2r-1`. Hence there are
branching rules:

.. CODE-BLOCK:: text

    ['B',r] => A1
    ['C',r] => A1

and these may be obtained using ``rule="symmetric_power"``.


Plethysms
---------

The above branching rules are sufficient for most cases, but a few
fall between the cracks. Mostly these involve maximal subgroups of
fairly small rank.

The rule ``rule="plethysm"`` is a powerful rule that includes any
branching rule from types `A`, `B`, `C` or `D` as a special case. Thus it
could be used in place of the above rules and would give the same
results. However, it is most useful when branching from `G` to a
maximal subgroup `H` such that `rank(H) < rank(G)-1`.

We consider a homomorphism `H \to G` where `G` is one of `SL(r+1)`,
`SO(2r+1)`, `Sp(2r)` or `SO(2r)`. The function
``branching_rule_from_plethysm`` produces the corresponding branching
rule. The main ingredient is the character `\chi` of the
representation of `H` that is the homomorphism to `GL(r+1)`,
`GL(2r+1)` or `GL(2r)`.

Let us consider the symmetric fifth power representation of
`SL(2)`. This is implemented above by ``rule="symmetric_power"``, but
suppose we want to use ``rule="plethysm"``. First we construct the
homomorphism by invoking its character, to be called ``chi``::

    sage: A1 = WeylCharacterRing("A1", style="coroots")
    sage: chi = A1([5])
    sage: chi.degree()
    6
    sage: chi.frobenius_schur_indicator()
    -1

This confirms that the character has degree 6 and is symplectic, so it
corresponds to a homomorphism `SL(2) \to Sp(6)`, and there is a
corresponding branching rule ``C3 -> A1``::

    sage: A1 = WeylCharacterRing("A1", style="coroots")
    sage: C3 = WeylCharacterRing("C3", style="coroots")
    sage: chi = A1([5])
    sage: sym5rule = branching_rule_from_plethysm(chi, "C3")
    sage: [C3(hwv).branch(A1, rule=sym5rule) for hwv in C3.fundamental_weights()]
    [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

This is identical to the results we would obtain using
``rule="symmetric_power"``::

    sage: A1 = WeylCharacterRing("A1", style="coroots")
    sage: C3 = WeylCharacterRing("C3", style="coroots")
    sage: [C3(v).branch(A1, rule="symmetric_power") for v in C3.fundamental_weights()]
    [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

But the next example of plethysm gives a branching rule not available
by other methods::

    sage: G2 = WeylCharacterRing("G2", style="coroots")
    sage: D7 = WeylCharacterRing("D7", style="coroots")
    sage: ad = G2.adjoint_representation(); ad.degree()
    14
    sage: ad.frobenius_schur_indicator()
    1
    sage: for r in D7.fundamental_weights():  # long time (1.29s)
    ....:    print(D7(r).branch(G2, rule=branching_rule_from_plethysm(ad, "D7")))
    G2(0,1)
    G2(0,1) + G2(3,0)
    G2(0,0) + G2(2,0) + G2(3,0) + G2(0,2) + G2(4,0)
    G2(0,1) + G2(2,0) + G2(1,1) + G2(0,2) + G2(2,1) + G2(4,0) + G2(3,1)
    G2(1,0) + G2(0,1) + G2(1,1) + 2*G2(3,0) + 2*G2(2,1) + G2(1,2) + G2(3,1) + G2(5,0) + G2(0,3)
    G2(1,1)
    G2(1,1)

In this example, `ad` is the 14-dimensional adjoint representation of the
exceptional group `G_2`. Since the Frobenius-Schur indicator is 1, the
representation is orthogonal, and factors through `SO(14)`, that is, `D7`.

We do not actually have to create the character (or for that matter
its ambient WeylCharacterRing) in order to create the branching rule::

    sage: branching_rule("D4","A2.adjoint_representation()","plethysm")
    plethysm (along A2(1,1)) branching rule D4 => A2

The adjoint representation of any semisimple Lie group is orthogonal, so we
do not need to compute the Frobenius-Schur indicator.

Miscellaneous other subgroups
-----------------------------

Use ``rule="miscellaneous"`` for the following rules. Every maximal
subgroup `H` of an exceptional group `G` are either among these,
or the five `A_1` subgroups described in the next section,
or (if `G` and `H` have the same rank) is available using
``rule="extended"``.

    .. MATH::

        \begin{aligned}
        B_3 & \to G_2,
        \\ E_6 & \to A_2,
        \\ E_6 & \to G_2,
        \\ F_4 & \to G_2 \times A_1,
        \\ E_6 & \to G_2 \times A_2,
        \\ E_7 & \to G_2 \times C_3,
        \\ E_7 & \to F_4 \times A_1,
        \\ E_7 & \to A_1 \times A_1,
        \\ E_7 & \to G_2 \times A_1,
        \\ E_7 & \to A_2
        \\ E_8 & \to G_2 \times F_4.
        \\ E_8 & \to A_2 \times A_1.
        \\ E_8 & \to B_2
        \end{aligned}

The first rule corresponds to the embedding of `G_2` in
`\hbox{SO}(7)` in its action on the trace zero octonions.
The two branching rules from `E_6` to `G_2` or `A_2`
are described in [Testerman1989]_. We caution the reader
that Theorem G.2 of that paper, proved there in positive
characteristic is false over the complex numbers. On
the other hand, the assumption of characteristic `p`
is not important for Theorems G.1 and A.1, which
describe the torus embeddings, hence contain enough
information to compute the branching rule. There
are other ways of embedding `G_2` or `A_2` into
`E_6`.  These may embeddings be characterized by the
condition that the two 27-dimensional representations of
`E_6` restrict irreducibly to `G_2` or `A_2`.
Their images are maximal subgroups.

The remaining rules come about as follows. Let `G` be
`F_4`, `E_6`, `E_7` or `E_8`, and let `H` be `G_2`,
or else (if `G=E_7`) `F_4`. We embed `H` into `G`
in the most obvious way; that is, in the chain
of subgroups

    .. MATH::

       G_2\subset F_4\subset E_6 \subset E_7 \subset E_8

Then the centralizer of `H` is `A_1`, `A_2`, `C_3`, `F_4` (if `H=G_2`) or
`A_1` (if `G=E_7` and `H=F_4`). This gives us five of the cases.
Regarding the branching rule `E_6 \to G_2 \times A_2`, Rubenthaler
[Rubenthaler2008]_ describes the embedding and applies it in an interesting
way.

The embedding of `A_1\times A_1` into `E_7` is as
follows. Deleting the 5 node of the `E_7` Dynkin
diagram gives the Dynkin diagram of `A_4\times A_2`, so this is a Levi
subgroup. We embed `\hbox{SL}(2)` into this Levi subgroup via the
representation `[4]\otimes[2]`.  This embeds the first copy of `A_1`. The
other `A_1` is the connected centralizer. See [Seitz1991]_, particularly the
proof of (3.12).

The embedding if `G_2\times A_1` into `E_7` is as
follows. Deleting the 2 node of the `E_7` Dynkin
diagram gives the `A_6` Dynkin diagram, which is
the Levi subgroup `\hbox{SL}(7)`. We embed `G_2` into
`\hbox{SL}(7)` via the irreducible seven-dimensional representation
of `G_2`. The `A_1` is the centralizer.

The embedding if `A_2\times A_1` into `E_8` is as
follows. Deleting the 2 node of the `E_8` Dynkin
diagram gives the `A_7` Dynkin diagram, which is
the Levi subgroup `\hbox{SL}(8)`. We embed `A_2` into
`\hbox{SL}(8)` via the irreducible eight-dimensional adjoint
representation of `\hbox{SL}(2)`. The `A_1` is the centralizer.

The embedding `A_2` into `E_7` is proved in
[Seitz1991]_ (5.8). In particular, he computes the
embedding of the `\hbox{SL}(3)` torus in the
`E_7` torus, which is what is needed to implement
the branching rule. The embedding of `B_2` into
`E_8` is also constructed in [Seitz1991]_ (6.7).
The embedding of the `B_2` Cartan subalgebra,
needed to implement the branching rule, is
easily deduced from (10) on page 111.

Maximal A1 subgroups of Exceptional Groups
------------------------------------------

There are seven embeddings of `SL(2)` into an exceptional
group as a maximal subgroup: one each for `G_2` and `F_4`,
two nonconjugate embeddings for `E_7` and three for `E_8`
These are constructed in [Testerman1992]_. Create the
corresponding branching rules as follows. The names of
the rules are roman numerals referring to the seven
cases of Testerman's Theorem 1::

       sage: branching_rule("G2","A1","i")
       i branching rule G2 => A1
       sage: branching_rule("F4","A1","ii")
       ii branching rule F4 => A1
       sage: branching_rule("E7","A1","iii")
       iii branching rule E7 => A1
       sage: branching_rule("E7","A1","iv")
       iv branching rule E7 => A1
       sage: branching_rule("E8","A1","v")
       v branching rule E8 => A1
       sage: branching_rule("E8","A1","vi")
       vi branching rule E8 => A1
       sage: branching_rule("E8","A1","vii")
       vii branching rule E8 => A1

The embeddings are characterized by the root
restrictions in their branching rules: usually
a simple root of the ambient group `G` restricts
to the unique simple root of `A_1`, except for
root `\alpha_4` for rules iv, vi and vii,
and the root `\alpha_6` for root vii; this is
essentially the way Testerman characterizes
the embeddings, and this information may
be obtained from Sage by employing the 
``describe()`` method of the branching rule.
Thus::

       sage: branching_rule("E8","A1","vii").describe()
       <BLANKLINE>
               O 2
               |
               |
       O---O---O---O---O---O---O---O
       1   3   4   5   6   7   8   0
       E8~
       root restrictions E8 => A1:
       <BLANKLINE>
       O
       1
       A1
       <BLANKLINE>       
       1 => 1
       2 => 1
       3 => 1
       4 => (zero)
       5 => 1
       6 => (zero)
       7 => 1
       8 => 1
       <BLANKLINE>       
       For more detailed information use verbose=True

Writing your own branching rules
--------------------------------

Sage has many built-in branching rules. Indeed, at least
up to rank eight (including all the exceptional groups)
branching rules to all maximal subgroups are implemented
as built in rules, except for a few obtainable using
``branching_rule_from_plethysm``. This means that
all the rules in [McKayPatera1981]_ are available in Sage.

Still in this section we are including instructions for coding a rule by
hand. As we have already explained, the branching rule is a function from the
weight lattice of ``G`` to the weight lattice of ``H``, and if you supply this
you can write your own branching rules.

As an example, let us consider how to implement the branching rule
``A3 -> C2``.  Here ``H = C2 = Sp(4)`` embedded as a subgroup in
``A3 = GL(4)``. The Cartan subalgebra `\hbox{Lie}(U)` consists of
diagonal matrices with eigenvalues ``u1, u2, -u2, -u1``. Then
``C2.space()`` is the two dimensional vector spaces consisting of the
linear functionals ``u1`` and ``u2`` on ``U``. On the other hand
`Lie(T) = \mathbf{R}^4`. A convenient way to see the restriction is to
think of it as the adjoint of the map ``[u1,u2] -> [u1,u2,-u2,-u1]``,
that is, ``[x0,x1,x2,x3] -> [x0-x3,x1-x2]``. Hence we may encode the
rule:

.. CODE-BLOCK:: python

    def brule(x):
        return [x[0]-x[3], x[1]-x[2]]

or simply:

.. CODE-BLOCK:: python

    brule = lambda x: [x[0]-x[3], x[1]-x[2]]

Let us check that this agrees with the built-in rule::

    sage: A3 = WeylCharacterRing(['A', 3])
    sage: C2 = WeylCharacterRing(['C', 2])
    sage: brule = lambda x: [x[0]-x[3], x[1]-x[2]]
    sage: A3(1,1,0,0).branch(C2, rule=brule)
    C2(0,0) + C2(1,1)
    sage: A3(1,1,0,0).branch(C2, rule="symmetric")
    C2(0,0) + C2(1,1)

Although this works, it is better to make the rule
into an element of the :class:`BranchingRule` class, as follows.

::

    sage: brule = BranchingRule("A3","C2",lambda x : [x[0]-x[3], x[1]-x[2]],"custom")
    sage: A3(1,1,0,0).branch(C2, rule=brule)
    C2(0,0) + C2(1,1)

Automorphisms and triality
--------------------------

The case where `G=H` can be treated as a special case of a branching
rule. In most cases if `G` has a nontrivial outer automorphism, it
has order two, corresponding to the symmetry of the Dynkin diagram.
Such an involution exists in the cases `A_r`, `D_r`, `E_6`.

So the automorphism acts on the representations of `G`, and its
effect may be computed using the branching rule code::

    sage: A4 = WeylCharacterRing("A4",style="coroots")
    sage: A4(1,0,1,0).degree()
    45
    sage: A4(0,1,0,1).degree()
    45
    sage: A4(1,0,1,0).branch(A4,rule="automorphic")
    A4(0,1,0,1)

In the special case where ``G=D4``, the Dynkin diagram has
extra symmetries, and these correspond to outer automorphisms
of the group. These are implemented as the ``"triality"``
branching rule::

    sage: branching_rule("D4","D4","triality").describe()
    <BLANKLINE>
        O 4
        |
        |
    O---O---O
    1   |2  3
        |
        O 0
    D4~
    root restrictions D4 => D4:
    <BLANKLINE>
        O 4
        |
        |
    O---O---O
    1   2   3
    D4
    <BLANKLINE>
    1 => 3
    2 => 2
    3 => 4
    4 => 1
    <BLANKLINE>
    For more detailed information use verbose=True

Triality his is not an automorphisms of `SO(8)`, but
of its double cover `spin(8)`. Note that `spin(8)` has
three representations of degree 8, namely the standard
representation of `SO(8)` and the two
eight-dimensional spin representations. These are
permuted by triality::

    sage: D4 = WeylCharacterRing("D4",style="coroots")
    sage: D4(0,0,0,1).branch(D4,rule="triality")
    D4(1,0,0,0)
    sage: D4(0,0,0,1).branch(D4,rule="triality").branch(D4,rule="triality")
    D4(0,0,1,0)
    sage: D4(0,0,0,1).branch(D4,rule="triality").branch(D4,rule="triality").branch(D4,rule="triality")
    D4(0,0,0,1)

By contrast, ``rule="automorphic"`` simply interchanges the two
spin representations, as it always does in type `D`::

    sage: D4(0,0,0,1).branch(D4,rule="automorphic")
    D4(0,0,1,0)
    sage: D4(0,0,1,0).branch(D4,rule="automorphic")
    D4(0,0,0,1)
