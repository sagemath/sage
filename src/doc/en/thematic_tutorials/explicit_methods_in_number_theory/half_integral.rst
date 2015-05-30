Half Integral Weight Forms
==========================

Basmaji's Algorithm
-------------------

Basmaji
(page 55 of his Essen thesis, "Ein Algorithmus zur Berechnung von
Hecke-Operatoren und Anwendungen auf modulare Kurven",
http://wstein.org/scans/papers/basmaji/).

Let :math:`S = S_{k+1}(\varepsilon)` be the space of cusp forms of
even integer weight :math:`k+1` and character :math:`\varepsilon = \chi
\psi^{(k+1)/2}`, where :math:`\psi` is the nontrivial mod-4 Dirichlet
character. Let :math:`U` be the subspace of :math:`S \times S` of
elements :math:`(a,b)` such that :math:`\Theta_2 a = \Theta_3 b`. Then
:math:`U` is isomorphic to :math:`S_{k/2}(\chi)` via the map
:math:`(a,b) \mapsto a/\Theta_3`.

This algorithm is implemented in Sage. I'm sure it could be
implemented in a way that is much faster than the current
implementation...

::

    sage: half_integral_weight_modform_basis(DirichletGroup(16,QQ).1, 3, 10)
    []
    sage: half_integral_weight_modform_basis(DirichletGroup(16,QQ).1, 5, 10)
    [q - 2*q^3 - 2*q^5 + 4*q^7 - q^9 + O(q^10)]
    sage: half_integral_weight_modform_basis(DirichletGroup(16*7).0^2,3,30)
    [q - 2*q^2 - q^9 + 2*q^14 + 6*q^18 - 2*q^21 - 4*q^22 - q^25 + O(q^30),
     q^2 - q^14 - 3*q^18 + 2*q^22 + O(q^30),
     q^4 - q^8 - q^16 + q^28 + O(q^30), q^7 - 2*q^15 + O(q^30)]
