Generators for Rings of Modular Forms
=====================================

Computing Generators
--------------------

For any congruence subgroup :math:`\Gamma`, the direct sum

.. math::

   M(\Gamma) =  \bigoplus_{k\geq 0} M_k(\Gamma)

is a ring, since the product of modular forms
:math:`f\in M_k(\Gamma)` and :math:`g \in M_{k'}(\Gamma)` is
an element :math:`fg \in M_{k+k'}(\Gamma)`. Sage can compute
likely generators for rings of modular forms, but currently doesn't
prove any of these results.

We verify the statement proved in Serre's "A Course in Arithmetic"
that :math:`E_4` and :math:`E_6` generate the space of level
one modular forms.

::

    sage: ModularFormsRing(SL2Z).generators(prec=4)
    [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + O(q^4)),
     (6, 1 - 504*q - 16632*q^2 - 122976*q^3 + O(q^4))]

Have you ever wondered which forms generate the ring :math:`M(\Gamma_0(2))`? It
turns out that one form of weight 2 and one form of weight 4 suffice.

::

    sage: ModularFormsRing(Gamma0(2)).generators(prec=12)
    [(2, 1 + 24*q + 24*q^2 + 96*q^3 + 24*q^4 + 144*q^5 + 96*q^6 + 192*q^7 + 24*q^8 + 312*q^9 + 144*q^10 + 288*q^11 + O(q^12)),
    (4, 1 + 240*q^2 + 2160*q^4 + 6720*q^6 + 17520*q^8 + 30240*q^10 + O(q^12))]

Here's generators for :math:`M(\Gamma_0(3))`. Notice that
elements of weight :math:`6` are now required, in addition to
weights :math:`2` and :math:`4`.

::

    sage: ModularFormsRing(Gamma0(3)).generators()
    [(2, 1 + 12*q + 36*q^2 + 12*q^3 + 84*q^4 + 72*q^5 + 36*q^6 + 96*q^7 + 180*q^8 + 12*q^9 + O(q^10)),
    (4, 1 + 240*q^3 + 2160*q^6 + 6720*q^9 + O(q^10)),
    (6, 1 - 504*q^3 - 16632*q^6 - 122976*q^9 + O(q^10))]

(*Note*: As of 2012, updates to the code mean that the output of this test is
not quite the same as it was in 2008, but of course there are multiple equally
valid answers.)

We can also handle rings of modular forms for odd congruence subgroups, but
with the usual caveat that we can't calculate forms of weight 1. So these are
elements generating the graded ring of forms of weight 0 or :math:`\ge 2`.

::

    sage: ModularFormsRing(Gamma1(3)).generators()
    [(2, 1 + 12*q + 36*q^2 + 12*q^3 + 84*q^4 + 72*q^5 + 36*q^6 + 96*q^7 + 180*q^8 + 12*q^9 + O(q^10)),
    (3, 1 + 54*q^2 + 72*q^3 + 432*q^5 + 270*q^6 + 918*q^8 + 720*q^9 + O(q^10)),
    (3, q + 3*q^2 + 9*q^3 + 13*q^4 + 24*q^5 + 27*q^6 + 50*q^7 + 51*q^8 + 81*q^9 + O(q^10)),
    (4, 1 + 240*q^3 + 2160*q^6 + 6720*q^9 + O(q^10))]


