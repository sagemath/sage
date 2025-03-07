Level One Modular Forms
=======================

Computing :math:`\Delta`
------------------------

The modular form

.. math::

   \Delta = q\prod(1-q^n)^{24} = \sum \tau(n)q^n

is perhaps the world's most famous modular form. We compute some terms
from the definition.

::

    sage: R.<q> = QQ[[]]
    sage: q * prod( 1-q^n+O(q^6) for n in (1..5) )^24
    q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 + O(q^7)

There are much better ways to compute :math:`\Delta`, which
amount to just a few polynomial multiplications over
:math:`\ZZ`.

::

    sage: D = delta_qexp(10^5)      # less than 10 seconds
    sage: D[:10]
    q - 24*q^2 + 252*q^3 - 1472*q^4 + ...
    sage: [p for p in primes(10^5) if D[p] % p == 0]
    [2, 3, 5, 7, 2411]
    sage: D[2411]
    4542041100095889012
    sage: f = eisenstein_series_qexp(12,6) - D[:6]; f
    691/65520 + 2073*q^2 + 176896*q^3 + 4197825*q^4 + 48823296*q^5 + O(q^6)
    sage: f % 691
    O(q^6)

The Victor Miller Basis
-----------------------
The Victor Miller basis for
:math:`M_k(\mathrm{SL}_2(\ZZ))` is the reduced row echelon
basis. It's a lemma that it has all integer coefficients, and a
rather nice diagonal shape.

::

    sage: victor_miller_basis(24, 6)
    [1 + 52416000*q^3 + 39007332000*q^4 + 6609020221440*q^5 + O(q^6),
     q + 195660*q^3 + 12080128*q^4 + 44656110*q^5 + O(q^6),
     q^2 - 48*q^3 + 1080*q^4 - 15040*q^5 + O(q^6)]
    sage: from sage.modular.dims import dimension_modular_forms
    sage: dimension_modular_forms(1,200)
    17
    sage: B = victor_miller_basis(200, 18) #5 seconds
    sage: B
    [1 + 79288314420681734048660707200000*q^17 + O(q^18),
     q + 2687602718106772837928968846869*q^17 + O(q^18),
     q^2 + 85789116961248834349485762560*q^17 + O(q^18),
     q^3 + 2567045661341737693075080984*q^17 + O(q^18),
     q^4 + 71629117222531314878690304*q^17 + O(q^18),
     q^5 + 1852433650992110376992650*q^17 + O(q^18),
     q^6 + 44081333191517147315712*q^17 + O(q^18),
     q^7 + 956892703246212300900*q^17 + O(q^18),
     q^8 + 18748755998771700480*q^17 + O(q^18),
     q^9 + 327218645736859401*q^17 + O(q^18),
     q^10 + 5001104379048960*q^17 + O(q^18),
     q^11 + 65427591611128*q^17 + O(q^18),
     q^12 + 709488619776*q^17 + O(q^18),
     q^13 + 6070433286*q^17 + O(q^18),
     q^14 + 37596416*q^17 + O(q^18),
     q^15 + 138420*q^17 + O(q^18),
     q^16 + 96*q^17 + O(q^18)]

Note: Craig Citro has made the above computation an order of
magnitude faster in code he has not quite got into Sage yet.

   "I'll clean those up and submit them soon, since I need them for
   something I'm working on ... I'm currently in the process of making
   spaces of modular forms of level one subclass the existing code,
   and actually take advantage of all our fast :math:`E_k` and
   :math:`\Delta` computation code, as well as cleaning things up a
   bit."
