=======================
The Weyl Character Ring
=======================

.. linkall

Weyl character rings
--------------------

The Weyl character ring is the representation ring of a compact Lie group. It
has a basis consisting of the irreducible representations of `G`, or
equivalently, their characters. The addition and multiplication in the Weyl
character ring correspond to direct sum and tensor product of representations.

Methods of the ambient space
----------------------------

In Sage, many useful features of the Lie group are available as
methods of the ambient space::

    sage: S = RootSystem("B2").ambient_space(); S
    Ambient space of the Root system of type ['B', 2]
    sage: S.roots()
    [(1, -1), (1, 1), (1, 0), (0, 1), (-1, 1), (-1, -1), (-1, 0), (0, -1)]
    sage: S.fundamental_weights()
    Finite family {1: (1, 0), 2: (1/2, 1/2)}
    sage: S.positive_roots()
    [(1, -1), (1, 1), (1, 0), (0, 1)]
    sage: S.weyl_group()
    Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)


Methods of the Weyl character ring
----------------------------------

If you are going to work with representations, you may want to create
a *Weyl Character ring*. Many methods of the ambient space are
available as methods of the Weyl character ring::

    sage: B3 = WeylCharacterRing("B3")
    sage: B3.fundamental_weights()
    Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1/2, 1/2, 1/2)}
    sage: B3.simple_roots()
    Finite family {1: (1, -1, 0), 2: (0, 1, -1), 3: (0, 0, 1)}
    sage: B3.dynkin_diagram()
    O---O=>=O
    1   2   3
    B3

Other useful methods of the Weyl character ring include:

- ``cartan_type``

- ``highest_root``

- ``positive_root``

- ``extended_dynkin_diagram``

- ``rank``

Some methods of the ambient space are not implemented as methods of
the Weyl character ring. However, the ambient space itself is a
method, and so you have access to its methods from the Weyl character
ring::

    sage: B3 = WeylCharacterRing("B3")
    sage: B3.space().weyl_group()
    Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space)
    sage: B3.space()
    Ambient space of the Root system of type ['B', 3]
    sage: B3.space().rho()
    (5/2, 3/2, 1/2)
    sage: B3.cartan_type()
    ['B', 3]


Coroot notation
---------------

It is useful to give the Weyl character ring a name that corresponds
to its Cartan type. This has the effect that the ring can parse its
own output::

    sage: G2 = WeylCharacterRing("G2")
    sage: [G2(fw) for fw in G2.fundamental_weights()]
    [G2(1,0,-1), G2(2,-1,-1)]
    sage: G2(1,0,-1)
    G2(1,0,-1)

Actually the prefix for the ring is configurable, so you don't really
have to call this ring ``G2``. Type ``WeylCharacterRing?`` at the
``sage:`` prompt for details.

There is one important option that you may want to know about. This
is *coroot notation*. You select this by specifying the option
``style="coroots"`` when you create the ring. With the coroot style,
the fundamental weights are represented ``(1,0,0, ...)``,
``(0,1,0,...)`` instead of as vectors in the ambient space::

    sage: B3 = WeylCharacterRing("B3", style="coroots")
    sage: [B3(fw) for fw in B3.fundamental_weights()]
    [B3(1,0,0), B3(0,1,0), B3(0,0,1)]
    sage: B3(0,0,1)
    B3(0,0,1)
    sage: B3(0,0,1).degree()
    8

The last representation is the eight dimensional spin representation
of `G = spin(7)`, the double cover of the orthogonal group `SO(7)`. In
the default notation it would be represented ``B3(1/2,1/2,1/2)``.

With the coroot notation every irreducible representation is
represented ``B3(a,b,c)`` where ``a``, ``b`` and ``c`` are nonnegative
integers. This is often convenient. For many purposes the coroot style
is preferable.

One disadvantage: in the coroot style the Lie group or Lie algebra is
treated as semisimple, so you lose the distinction between `GL(n)` and
`SL(n)`; you also some information about representations of E6 and E7
for the same reason.


Tensor products of representations
----------------------------------

The multiplication in the Weyl character ring corresponds to tensor
product. This gives us a convenient way of decomposing a tensor
product into irreducibles::

    sage: B3 = WeylCharacterRing("B3")
    sage: fw = B3.fundamental_weights()
    sage: spinweight = fw[3]; spinweight
    (1/2, 1/2, 1/2)
    sage: spin = B3(spinweight); spin
    B3(1/2,1/2,1/2)
    sage: spin.degree()
    8

The element `spin` of the WeylCharacterRing is the representation
corresponding to the third highest weight representation, the
eight-dimensional spin representation of `spin(7)`. We could
just as easily construct it with the command::

    sage: spin = B3(1/2,1/2,1/2)

We may compute its tensor product with itself, using the
multiplicative structure of the Weyl character ring::

    sage: chi = spin*spin; chi
    B3(0,0,0) + B3(1,0,0) + B3(1,1,0) + B3(1,1,1)

We have taken the eight-dimensional spin representation and tensored
with itself. We see that the tensor square splits into four
irreducibles, each with multiplicity one.

The highest weights that appear here are available (with their
coefficients) through the usual free module accessors::

    sage: from pprint import pprint
    sage: list(chi)                           # random
    [((1, 1, 1), 1), ((0, 0, 0), 1), ((1, 0, 0), 1), ((1, 1, 0), 1)]
    sage: sorted(chi, key=str)
    [((0, 0, 0), 1), ((1, 0, 0), 1), ((1, 1, 0), 1), ((1, 1, 1), 1)]
    sage: pprint(dict(chi))
    {(0, 0, 0): 1, (1, 0, 0): 1, (1, 1, 0): 1, (1, 1, 1): 1}
    sage: M = sorted(chi.monomials(), key=lambda x: tuple(x.support())); M
    [B3(0,0,0), B3(1,0,0), B3(1,1,0), B3(1,1,1)]
    sage: sorted(chi.support())
    [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)]
    sage: chi.coefficients()
    [1, 1, 1, 1]
    sage: [r.degree() for r in M]
    [1, 7, 21, 35]
    sage: sum(r.degree() for r in chi.monomials())
    64

Here we have extracted the individual representations, computed
their degrees and checked that they sum up to `64`.


Weight multiplicities
---------------------

The weights of the character are available (with their coefficients)
through the method ``weight_multiplicities``. Continuing from the
example in the last section::

    sage: pprint(chi.weight_multiplicities())
    {(0, 0, 0): 8, (-1, 0, 0): 4, (-1, -1, 0): 2, (-1, -1, -1): 1,
     (-1, -1, 1): 1, (-1, 1, 0): 2, (-1, 1, -1): 1, (-1, 1, 1): 1,
     (-1, 0, -1): 2, (-1, 0, 1): 2, (1, 0, 0): 4, (1, -1, 0): 2,
     (1, -1, -1): 1, (1, -1, 1): 1, (1, 1, 0): 2, (1, 1, -1): 1,
     (1, 1, 1): 1, (1, 0, -1): 2, (1, 0, 1): 2, (0, -1, 0): 4,
     (0, -1, -1): 2, (0, -1, 1): 2, (0, 1, 0): 4, (0, 1, -1): 2,
     (0, 1, 1): 2, (0, 0, -1): 4, (0, 0, 1): 4}

Each key of this dictionary is a weight, and its value is the
multiplicity of that weight in the character.

Example
-------

Suppose that we wish to compute the integral

.. MATH::

   \int_{U(n)} |tr(g)|^{2k}\,dg

for various `n`. Here `U(n)` is the unitary group, which is the maximal
compact subgroup of `GL(n,\mathbb{C})`. The irreducible unitary representations
of `U(n)` may be regarded as the basis elements of the WeylCharacterRing of
type `A_r`, where `r=n-1` so we might work in that ring. The trace `tr(g)` is
then just the character of the standard representation. We may realize
it in the WeylCharacterRing by taking the first fundamental weight and
coercing it into the ring. For example, if `k=5` and `n=3` so `r=2`::

    sage: A2 = WeylCharacterRing("A2")
    sage: fw = A2.fundamental_weights(); fw
    Finite family {1: (1, 0, 0), 2: (1, 1, 0)}
    sage: tr = A2(fw[1]); tr
    A2(1,0,0)

We may compute the norm square the character ``tr^5`` by decomposing it into
irreducibles, and taking the sum of the squares of their multiplicities. By
Schur orthogonality, this gives the inner product of the `tr(g)^5` with
itself, that is, the integral of `|tr(g)|^{10}`::

    sage: sum(d^2 for d in (tr^5).coefficients())
    103

So far we have been working with `n=3`. For general `n`::

    sage: def f(n,k):
    ....:     R = WeylCharacterRing(['A',n-1])
    ....:     tr = R(R.fundamental_weights()[1])
    ....:     return sum(d^2 for d in (tr^k).coefficients())
    sage: [f(n,5) for n in [2..7]]
    [42, 103, 119, 120, 120, 120]

We see that the 10-th moment of `tr(g)` is just `5!` when `n` is sufficiently
large. What if we fix `n` and vary `k`?

::

    sage: [f(2,k) for k in [1..10]]
    [1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]
    sage: [catalan_number(k) for k in [1..10]]
    [1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]


Frobenius-Schur indicator
-------------------------

The Frobeinus-Schur indicator of an irreducible representation of a
compact Lie group `G` with character `\chi` is:

.. MATH::

    \int_G\chi(g^2) \, dg

The Haar measure is normalized so that `vol(G) = 1`. The
Frobenius-Schur indicator equals `1` if the representation is real
(orthogonal), `-1` if it is quaternionic (symplectic) and `0` if it is
complex (not self-contragredient). This is a method of weight ring
elements corresponding to irreducible representations. Let us compute
the Frobenius-Schur indicators of the spin representations of some
odd spin groups::

    sage: def spinrepn(r):
    ....:     R = WeylCharacterRing(['B',r])
    ....:     return R(R.fundamental_weights()[r])
    ....:
    sage: spinrepn(3)
    B3(1/2,1/2,1/2)
    sage: for r in [1..4]:
    ....:     print("{} {}".format(r, spinrepn(r).frobenius_schur_indicator()))
    1 -1
    2 -1
    3 1
    4 1

Here we have defined a function that returns the spin representation
of the group `spin(2r+1)` with Cartan type `['B',r]`, then computed
the Frobenius-Schur indicators for a few values. From this experiment
we see that the spin representations of `spin(3)` and `spin(5)` are
symplectic, while those of `spin(7)` and `spin(9)` are orthogonal.


Symmetric and exterior powers
-----------------------------

Sage can compute symmetric and exterior powers of a representation::

    sage: B3 = WeylCharacterRing("B3",style="coroots")
    sage: spin = B3(0,0,1); spin.degree()
    8
    sage: spin.exterior_power(2)
    B3(1,0,0) + B3(0,1,0)
    sage: spin.exterior_square()
    B3(1,0,0) + B3(0,1,0)
    sage: spin.exterior_power(5)
    B3(0,0,1) + B3(1,0,1)
    sage: spin.symmetric_power(5)
    B3(0,0,1) + B3(0,0,3) + B3(0,0,5)

The `k`-th exterior square of a representation is zero if `k`
is greater than the degree of the representation. However the
`k`-th symmetric power is nonzero for all `k`.

The tensor square of any representation decomposes as the direct sum
of the symmetric and exterior squares::

    sage: C4 = WeylCharacterRing("C4",style="coroots")
    sage: chi = C4(1,0,0,0); chi.degree()
    8
    sage: chi.symmetric_square()
    C4(2,0,0,0)
    sage: chi.exterior_square()
    C4(0,0,0,0) + C4(0,1,0,0)
    sage: chi^2 == chi.symmetric_square() + chi.exterior_square()
    True

Since in this example the exterior square contains the trivial
representation we expect the Frobenius-Schur indicator to be `-1`, and
indeed it is::

    sage: chi = C4(1,0,0,0)
    sage: chi.frobenius_schur_indicator()
    -1

This is not surprising since this is the standard representation
of a symplectic group, which is symplectic *by definition*!


Weyl dimension formula
----------------------

If the representation is truly large you will not be able to construct
it in the Weyl character ring, since internally it is represented by a
dictionary of its weights. If you want to know its degree, you can
still compute that since Sage implements the Weyl dimension
formula. The degree of the representation is implemented as a method
of its highest weight vector::

    sage: L = RootSystem("E8").ambient_space()
    sage: [L.weyl_dimension(f) for f in L.fundamental_weights()]
    [3875, 147250, 6696000, 6899079264, 146325270, 2450240, 30380, 248]

It is a fact that for any compact Lie group if `\rho` is the Weyl vector
(half the sum of the positive roots) then the degree of the irreducible
representation with highest weight `\rho` equals `2^N` where `N` is the number
of positive roots. Let us check this for `E_8`. In this case `N = 120`::

    sage: L = RootSystem("E8").ambient_space()
    sage: len(L.positive_roots())
    120
    sage: 2^120
    1329227995784915872903807060280344576
    sage: L.weyl_dimension(L.rho())
    1329227995784915872903807060280344576


.. _SLvsGL:

SL versus GL
------------

Sage takes the weight space for type ``['A',r]`` to be `r+1`
dimensional. As a by-product, if you create the Weyl character ring
with the command::

    sage: A2 = WeylCharacterRing("A2")

Then you are effectively working with `GL(3)` instead of `SL(3)`. For
example, the determinant is the character ``A2(1,1,1)``. However, as
we will explain later, you can work with `SL(3)` if you like, so long
as you are willing to work with fractional weights. On the other hand
if you create the Weyl character ring with the command::

    sage: A2 = WeylCharacterRing("A2", style="coroots")

Then you are working with `SL(3)`.

There are some advantages to this arrangement:

- The group `GL(r+1)` arises frequently in practice. For example, even
  if you care mainly about semisimple groups, the group `GL(r+1)` may
  arise as a Levi subgroup.

- It avoids fractional weights. If you want to work with `SL(3)` the
  fundamental weights are ``(2/3,-1/3,-1/3)`` and
  ``(1/3,1/3,-2/3)``. If you work instead with `GL(3)`, they are
  ``(1,0,0)`` and ``(1,1,0)``. For many mathematical purposes it
  doesn't make any difference which you use. This is because the
  difference between ``(2/3,-1/3,-1/3)`` and ``(1,0,0)`` is a vector
  that is orthogonal to all the simple roots. Thus these vectors are
  interchangeable. But for convenience avoiding fractional weights is
  advantageous.

However if you want to be an `SL` purist, Sage will support you. The
weight space for `SL(3)` can be taken to be the hyperplane in
`\mathbf{Q}^3` consisting of vectors `(a,b,c)` with
`a+b+c = 0`. The fundamental weights for SL(3) are then
``(2/3,-1/3,-1/3)`` and ``(1/3,1/3,-2/3)``, though Sage will tell you
they are ``(1,0,0)`` and ``(1,1,0)``. The work-around is to filter
them through the method ``coerce_to_sl`` as follows::

    sage: A2 = WeylCharacterRing("A2")
    sage: [fw1,fw2] = [w.coerce_to_sl() for w in A2.fundamental_weights()]
    sage: [standard, contragredient] = [A2(fw1), A2(fw2)]
    sage: standard, contragredient
    (A2(2/3,-1/3,-1/3), A2(1/3,1/3,-2/3))
    sage: standard*contragredient
    A2(0,0,0) + A2(1,0,-1)

Sage is not confused by the fractional weights. Note that if you use
coroot notation, you are working with `SL` automatically::

    sage: A2 = WeylCharacterRing("A2", style="coroots")
    sage: A2(1,0).weight_multiplicities()
    {(-1/3, -1/3, 2/3): 1, (-1/3, 2/3, -1/3): 1, (2/3, -1/3, -1/3): 1}

There is no convenient way to create the determinant in the Weyl
character ring if you adopt the coroot style.

Just as we coerced the fundamental weights into the `SL` weight
lattice, you may need to coerce the Weyl vector `\rho` if you are
working with `SL`. The default value for `\rho` in type `A_r` is
`(r,r-1,\dots,0)`, but if you are an `SL` purist you want

.. MATH::

    \left(\frac{r}{2}, \frac{r}{2}-1,\dots,-\frac{r}{2}\right).

Therefore take the value of `\rho` that you get from the method of the
ambient space and coerce it into `SL`::

    sage: A2 = WeylCharacterRing("A2", style="coroots")
    sage: rho = A2.space().rho().coerce_to_sl(); rho
    (1, 0, -1)
    sage: rho == (1/2)*sum(A2.space().positive_roots())
    True

You do not need to do this for other Cartan types. If you are working
with (say) `F4` then a `\rho` is a `\rho`::

    sage: F4 = WeylCharacterRing("F4")
    sage: L = F4.space()
    sage: rho = L.rho()
    sage: rho == (1/2)*sum(L.positive_roots())
    True


Integration
-----------

Suppose that we wish to compute the integral

.. MATH::

   \int_{U(n)} |tr(g)|^{2k}\,dg

for various `n`. Here `U(n)` is the unitary group, which is the maximal
compact subgroup of `GL(n,\mathbf{C})`, and `dg` is the Haar measure on
`U(n)`, normalized so that the volume of the group is 1.

The irreducible unitary representations of `U(n)` may be regarded as the basis
elements of the WeylCharacterRing of type `A_r`, where `r=n-1` so we might
work in that ring. The trace `tr(g)` is then just the character of the
standard representation. We may realize it in the WeylCharacterRing by taking
the first fundamental weight and coercing it into the ring. For example, if
`k=5` and `n=3` so `r=2`::

    sage: A2 = WeylCharacterRing("A2")
    sage: fw = A2.fundamental_weights(); fw
    Finite family {1: (1, 0, 0), 2: (1, 1, 0)}
    sage: tr = A2(fw[1]); tr
    A2(1,0,0)

We may compute the norm square the character ``tr^5`` by decomposing it into
irreducibles, and taking the sum of the squares of their multiplicities. By
Schur orthogonality, this gives the inner product of the `tr(g)^5` with
itself, that is, the integral of `|tr(g)|^{10}`::

    sage: tr^5
    5*A2(2,2,1) + 6*A2(3,1,1) + 5*A2(3,2,0) + 4*A2(4,1,0) + A2(5,0,0)
    sage: sorted((tr^5).monomials(), key=lambda x: tuple(x.support()))
    [A2(2,2,1), A2(3,1,1), A2(3,2,0), A2(4,1,0), A2(5,0,0)]
    sage: sorted((tr^5).coefficients())
    [1, 4, 5, 5, 6]
    sage: sum(x^2 for x in (tr^5).coefficients())
    103

So far we have been working with `n=3`. For general `n`::

    sage: def f(n,k):
    ....:    R = WeylCharacterRing(['A',n-1])
    ....:    tr = R(R.fundamental_weights()[1])
    ....:    return sum(x^2 for x in (tr^k).coefficients())
    ....:
    sage: [f(n,5) for n in [2..7]]  # long time (31s on sage.math, 2012)
    [42, 103, 119, 120, 120, 120]

We see that the 10-th moment of `tr(g)` is just `5!` when `n` is sufficiently
large. What if we fix `n` and vary `k`?

::

    sage: [f(2,k) for k in [1..10]]
    [1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]
    sage: [catalan_number(k) for k in [1..10]]
    [1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]


Invariants and multiplicities
-----------------------------

Sometimes we are only interested in the multiplicity of the trivial
representation in some character. This may be found by the method
``invariant_degree``. Continuing from the preceding example,

::

    sage: A2 = WeylCharacterRing("A2",style="coroots")
    sage: ad = A2(1,1)
    sage: [ad.symmetric_power(k).invariant_degree() for k in [0..6]]
    [1, 0, 1, 1, 1, 1, 2]
    sage: [ad.exterior_power(k).invariant_degree() for k in [0..6]]
    [1, 0, 0, 1, 0, 1, 0]

If we want the multiplicity of some other representation, we may
obtain that using the method ``multiplicity``::

    sage: (ad^3).multiplicity(ad)
    8


Weight Rings
------------

You may wish to work directly with the weights of a representation.

Weyl character ring elements are represented internally by a
dictionary of their weights with multiplicities. However these are
subject to a constraint: the coefficients must be invariant under the
action of the Weyl group.

The :class:`WeightRing` is also a ring whose elements are represented
internally by a dictionary of their weights with multiplicities, but
it is not subject to this constraint of Weyl group invariance. The
weights are allowed to be fractional, that is, elements of the ambient
space. In other words, the weight ring is the group algebra over the
ambient space of the weight lattice.

To create a :class:`WeightRing` first construct the
:class:`WeylCharacterRing`, then create the
:class:`WeightRing` as follows::

    sage: A2 = WeylCharacterRing(['A',2])
    sage: a2 = WeightRing(A2)

You may coerce elements of the :class:`WeylCharacterRing` into the weight
ring. For example, if you want to see the weights of the adjoint
representation of `GL(3)`, you may use the method ``mlist``, but
another way is to coerce it into the weight ring::

    sage: from pprint import pprint
    sage: A2 = WeylCharacterRing(['A',2])
    sage: ad = A2(1,0,-1)
    sage: pprint(ad.weight_multiplicities())
    {(0, 0, 0): 2, (-1, 1, 0): 1, (-1, 0, 1): 1, (1, -1, 0): 1,
     (1, 0, -1): 1, (0, -1, 1): 1, (0, 1, -1): 1}

This command produces a dictionary of the weights that appear in
the representation, together with their multiplicities. But another
way of getting the same information, with an aim of working with it,
is to coerce it into the weight ring::

    sage: a2 = WeightRing(A2)
    sage: a2(ad)
    2*a2(0,0,0) + a2(-1,1,0) + a2(-1,0,1) + a2(1,-1,0) + a2(1,0,-1) + a2(0,-1,1) + a2(0,1,-1)

For example, the Weyl denominator formula is usually written this way:

.. MATH::

    \prod_{\alpha\in\Phi^+}\left(e^{\alpha/2}-e^{-\alpha/2}\right)
    =
    \sum_{w\in W} (-1)^{l(w)}e^{w(\rho)}.

The notation is as follows. Here if `\lambda` is a weight, or more
generally, an element of the ambient space, then `e^\lambda` means the
image of `\lambda` in the group algebra of the ambient space of the
weight lattice `\lambda`. Since this group algebra is just the weight
ring, we can interpret `e^\lambda` as its image in the weight ring.

Let us confirm the Weyl denominator formula for ``A2``::

    sage: A2 = WeylCharacterRing("A2")
    sage: a2 = WeightRing(A2)
    sage: L = A2.space()
    sage: W = L.weyl_group()
    sage: rho = L.rho().coerce_to_sl()
    sage: lhs = prod(a2(alpha/2)-a2(-alpha/2) for alpha in L.positive_roots()); lhs
    a2(-1,1,0) - a2(-1,0,1) - a2(1,-1,0) + a2(1,0,-1) + a2(0,-1,1) - a2(0,1,-1)
    sage: rhs = sum((-1)^(w.length())*a2(w.action(rho)) for w in W); rhs
    a2(-1,1,0) - a2(-1,0,1) - a2(1,-1,0) + a2(1,0,-1) + a2(0,-1,1) - a2(0,1,-1)
    sage: lhs == rhs
    True

Note that we have to be careful to use the right value of `\rho`. The
reason for this is explained in :ref:`SLvsGL`.

We have seen that elements of the :class:`WeylCharacterRing` can be coerced
into the :class:`WeightRing`. Elements of the :class:`WeightRing` can be
coerced into the :class:`WeylCharacterRing` *provided* they are invariant
under the Weyl group.

