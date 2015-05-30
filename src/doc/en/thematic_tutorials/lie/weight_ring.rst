------------
Weight Rings
------------

.. linkall

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

