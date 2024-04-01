------------------------------------------------
Weyl Groups, Coxeter Groups and the Bruhat Order
------------------------------------------------


Classical and affine Weyl groups
--------------------------------

You can create Weyl groups and affine Weyl groups for any root
system. A variety of methods are available for these. Some of these
are methods are available for general Coxeter groups.

By default, elements of the Weyl group are represented as matrices::

    sage: WeylGroup("A3").simple_reflection(1)
    [0 1 0 0]
    [1 0 0 0]
    [0 0 1 0]
    [0 0 0 1]

You may prefer a notation in which elements are written out as
products of simple reflections. In order to implement this you need to
specify a prefix, typically ``"s"``::

    sage: W = WeylGroup("A3",prefix="s")
    sage: [s1,s2,s3] = W.simple_reflections()
    sage: (s1*s2*s1).length()
    3
    sage: W.long_element()
    s1*s2*s3*s1*s2*s1
    sage: s1*s2*s3*s1*s2*s1 == s3*s2*s1*s3*s2*s3
    True

The Weyl group acts on the ambient space of the root lattice, which is
accessed by the method ``domain``. To illustrate this, recall that if `w_0` is
the long element then `\alpha \mapsto -w_0(\alpha)` is a permutation of the
simple roots. We may compute this as follows::

    sage: W = WeylGroup("E6",prefix="s")
    sage: w0 = W.long_element(); w0
    s1*s3*s4*s5*s6*s2*s4*s5*s3*s4*s1*s3*s2*s4*s5*s6*s2*s4*s5*s3*s4*s1*s3*s2*s4*s5*s3*s4*s1*s3*s2*s4*s1*s3*s2*s1
    sage: sr = W.domain().simple_roots().list(); sr
    [(1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2), (1, 1, 0, 0, 0, 0, 0, 0),
    (-1, 1, 0, 0, 0, 0, 0, 0), (0, -1, 1, 0, 0, 0, 0, 0), (0, 0, -1, 1, 0, 0, 0, 0),
    (0, 0, 0, -1, 1, 0, 0, 0)]
    sage: [-w0.action(a) for a in sr]
    [(0, 0, 0, -1, 1, 0, 0, 0), (1, 1, 0, 0, 0, 0, 0, 0), (0, 0, -1, 1, 0, 0, 0, 0),
    (0, -1, 1, 0, 0, 0, 0, 0), (-1, 1, 0, 0, 0, 0, 0, 0),
    (1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2)]

We may ask when this permutation is trivial. If it is nontrivial it
induces an automorphism of the Dynkin diagram, so it must be
nontrivial when the Dynkin diagram has no automorphism. But if there
is a nontrivial automorphism, the permutation might or might not be
trivial::

    sage: def roots_not_permuted(ct):
    ....:     W = WeylGroup(ct)
    ....:     w0 = W.long_element()
    ....:     sr = W.domain().simple_roots()
    ....:     return all(a == -w0.action(a) for a in sr)
    sage: for ct in [CartanType(['D', r]) for r in [2..8]]:
    ....:     print("{} {}".format(ct, roots_not_permuted(ct)))
    ['D', 2] True
    ['D', 3] False
    ['D', 4] True
    ['D', 5] False
    ['D', 6] True
    ['D', 7] False
    ['D', 8] True

If `\alpha` is a root let `r_\alpha` denote the reflection in the
hyperplane that is orthogonal to `\alpha`. We reserve the notation `s_\alpha`
for the simple reflections, that is, the case where `\alpha` is a simple
root. The reflections are just the conjugates of the simple reflections.

The reflections are the values in a finite family, which is a wrapper
around a python dictionary. The keys are the positive roots, so
given a positive root, you can look up the corresponding reflection.
If you want a list of all reflections, you can use the usual methods to
construct a list (e.g., using the ``list`` function) or use the method
``values`` for the family of reflections::

    sage: W = WeylGroup("B3",prefix="s")
    sage: ref = W.reflections(); ref
    Finite family {(1, -1, 0): s1, (0, 1, -1): s2, ...}
    sage: [a1,a2,a3] = W.domain().simple_roots()
    sage: a1+a2+a3
    (1, 0, 0)
    sage: ref[a1+a2+a3]
    s1*s2*s3*s2*s1
    sage: sorted(ref)
    [s1*s2*s3*s2*s1,
     s2*s3*s1*s2*s3*s1*s2,
     s3*s1*s2*s3*s1,
     s1*s2*s1,
     s1,
     s2*s3*s2,
     s3*s2*s3,
     s2,
     s3]

If instead you want a family whose keys are the reflections
and whose values are the roots, you may use the inverse family::

    sage: from pprint import pprint
    sage: W = WeylGroup("B3",prefix="s")
    sage: [s1,s2,s3] = W.simple_reflections()
    sage: altref = W.reflections().inverse_family()
    sage: altref
    Finite family {s1*s2*s3*s2*s1: (1, 0, 0),
     s2*s3*s1*s2*s3*s1*s2: (1, 1, 0),
     s3*s1*s2*s3*s1: (1, 0, 1),
     s1*s2*s1: (1, 0, -1),
     s1: (1, -1, 0),
     s2*s3*s2: (0, 1, 0),
     s3*s2*s3: (0, 1, 1),
     s2: (0, 1, -1),
     s3: (0, 0, 1)}
    sage: altref[s3*s2*s3]
    (0, 1, 1)

.. NOTE::

    The behaviour of this function was changed in :issue:`20027`.

The Weyl group is implemented as a GAP matrix group. You therefore can
display its character table as follows::

    sage: WeylGroup("D4").character_table()
    CT1
    <BLANKLINE>
          2  6  4  5  1  3  5  5  4  3  3  1  4  6
          3  1  .  .  1  .  .  .  .  .  .  1  .  1
    <BLANKLINE>
            1a 2a 2b 6a 4a 2c 2d 2e 4b 4c 3a 4d 2f
    <BLANKLINE>
    X.1      1  1  1  1  1  1  1  1  1  1  1  1  1
    X.2      1 -1  1  1 -1  1  1 -1 -1 -1  1  1  1
    X.3      2  .  2 -1  .  2  2  .  .  . -1  2  2
    X.4      3 -1 -1  .  1  3 -1 -1 -1  1  . -1  3
    X.5      3  1 -1  . -1  3 -1  1  1 -1  . -1  3
    X.6      3 -1  3  . -1 -1 -1 -1  1  1  . -1  3
    X.7      3 -1 -1  .  1 -1  3 -1  1 -1  . -1  3
    X.8      3  1  3  .  1 -1 -1  1 -1 -1  . -1  3
    X.9      3  1 -1  . -1 -1  3  1 -1  1  . -1  3
    X.10     4 -2  . -1  .  .  .  2  .  .  1  . -4
    X.11     4  2  . -1  .  .  . -2  .  .  1  . -4
    X.12     6  . -2  .  . -2 -2  .  .  .  .  2  6
    X.13     8  .  .  1  .  .  .  .  .  . -1  . -8


Affine Weyl groups
------------------

Affine Weyl groups may be created the same way. You simply begin with
an affine Cartan type::

    sage: W = WeylGroup(['A',2,1],prefix="s")
    sage: W.cardinality()
    +Infinity
    sage: [s0,s1,s2] = W.simple_reflections()
    sage: s0*s1*s2*s1*s0
    s0*s1*s2*s1*s0

The affine Weyl group differs from a classical Weyl group since it is
infinite. The associated classical Weyl group is a subgroup that may be
extracted as follows::

    sage: W = WeylGroup(['A',2,1],prefix="s")
    sage: W1 = W.classical(); W1
    Parabolic Subgroup of the Weyl Group of type ['A', 2, 1] (as a matrix group
    acting on the root space)
    sage: W1.simple_reflections()
    Finite family {1: s1, 2: s2}

Although ``W1`` in this example is isomorphic to ``WeylGroup("A2")`` it
has a different matrix realization::

    sage: for s in WeylGroup(['A',2,1]).classical().simple_reflections():
    ....:     print(s)
    ....:     print("")
    [ 1  0  0]
    [ 1 -1  1]
    [ 0  0  1]
    <BLANKLINE>
    [ 1  0  0]
    [ 0  1  0]
    [ 1  1 -1]

    sage: for s in WeylGroup(['A',2]).simple_reflections():
    ....:     print(s)
    ....:     print("")
    [0 1 0]
    [1 0 0]
    [0 0 1]
    <BLANKLINE>
    [1 0 0]
    [0 0 1]
    [0 1 0]


Bruhat order
------------

The Bruhat partial order on the Weyl group may be defined as follows.

If `u,v \in W`, find a reduced expression of `v` into a product of
simple reflections: `v = s_1 \cdots s_n`. (It is not assumed that the
`s_i` are distinct.) If omitting some of the `s_i` gives a product
that represents `u`, then `u \le v`.

The Bruhat order is implemented in Sage as a method of Coxeter groups,
and so it is available for Weyl groups, classical or affine.

If `u`, `v \in W` then ``u.bruhat_le(v)`` returns ``True`` if
`u \le v` in the Bruhat order.

If `u \le v` then the *Bruhat interval* `[u,v]` is defined to be the
set of all `t` such that `u \le t \le v`. One might try to implement
this as follows::

    sage: W = WeylGroup("A2",prefix="s")
    sage: [s1,s2] = W.simple_reflections()
    sage: def bi(u,v) : return [t for t in W if u.bruhat_le(t) and t.bruhat_le(v)]
    ...
    sage: bi(s1,s1*s2*s1)
    [s1*s2, s2*s1, s1, s1*s2*s1]

This would not be a good definition since it would fail if `W` is
affine and be inefficient of `W` is large. Sage has a Bruhat interval
method::

    sage: W = WeylGroup("A2",prefix="s")
    sage: [s1,s2] = W.simple_reflections()
    sage: W.bruhat_interval(s1,s1*s2*s1)
    [s1*s2*s1, s2*s1, s1*s2, s1]

This works even for affine Weyl groups.


The Bruhat graph
----------------

References:

- [Carrell1994]_

- [Deodhar1977]_

- [Dyer1993]_

- [BumpNakasuji2010]_

The *Bruhat graph* is a structure on the Bruhat interval. Suppose that
`u \le v`. The vertices of the graph are `x` with `u \le x \le v`.
There is a vertex connecting `x,y \in [x,y]` if `x = y \cdot r` where
`r` is a reflection. If this is true then either `x < y` or `y < x`.

If `W` is a classical Weyl group the Bruhat graph is implemented in Sage::

    sage: W = WeylGroup("A3",prefix="s")
    sage: [s1,s2,s3] = W.simple_reflections()
    sage: bg = W.bruhat_graph(s2,s2*s1*s3*s2); bg
    Digraph on 10 vertices
    sage: bg.show3d()

The Bruhat graph has interesting regularity properties that were
investigated by Carrell and Peterson. It is a regular graph if both
the Kazhdan Lusztig polynomials `P_{u,v}` and `P_{w_0v,w_0u}` are 1,
where `w_0` is the long Weyl group element. It is closely related to
the *Deodhar conjecture*, which was proved by Deodhar, Carrell and
Peterson, Dyer and Polo.

Deodhar proved that if `u < v` then the Bruhat interval `[u,v]`
contains as many elements of odd length as it does of even length. We
observe that often this can be strengthened: If there exists a
reflection `r` such that left (or right) multiplication by `r` takes
the Bruhat interval `[u,v]` to itself, then this gives an explicit
bijection between the elements of odd and even length in `[u,v]`.

Let us search for such reflections. Put the following commands in a
file and load the file:

.. CODE-BLOCK:: python

    W = WeylGroup("A3",prefix="s")
    [s1,s2,s3] = W.simple_reflections()
    ref = W.reflections().keys()

    def find_reflection(u,v):
        bi = W.bruhat_interval(u,v)
        ret = []
        for r in ref:
            if all( r*x in bi for x in bi):
                ret.append(r)
        return ret

    for v in W:
        for u in W.bruhat_interval(1,v):
            if u != v:
                print((u,v,find_reflection(u,v)))

This shows that the Bruhat interval is stabilized by a reflection for
all pairs `(u,v)` with `u < v` except the following two:
`s_3s_1,s_1s_2s_3s_2s_1` and `s_2,s_2s_3s_1s_2`. Now these are
precisely the pairs such that `u\prec v` in the notation of Kazhdan
and Lusztig, and `l(v)-l(u) > 1`. One should not rashly suppose that
this is a general characterization of the pairs `(u,v)` such that no
reflection stabilizes the Bruhat interval, for this is not true. However
it does suggest that the question is worthy of further investigation.
