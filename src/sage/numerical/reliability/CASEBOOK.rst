Numerical Reliability Casebook
==============================

This document collects concrete examples of numerical computations that
appear correct but are in fact unreliable. Each case demonstrates how
numerical reliability diagnostics can explain or detect such failures.

Case 1: Multiple Roots
----------------------

Problem
^^^^^^^

Computing roots of a polynomial with a multiple root.

::

    sage: R.<x> = RR[]
    sage: p = (x - 1)^2
    sage: p.roots()
    [(1.00000000000000, 2)]

Issue
^^^^^

Although the numerical root is correct, the problem is ill-conditioned.
Small perturbations in the coefficients can significantly change the
roots.

Reliability Diagnostic
^^^^^^^^^^^^^^^^^^^^^^

::

    sage: r = p.roots()[0][0]
    sage: r.reliability()
    'low'

Explanation
^^^^^^^^^^^

Residual-based diagnostics detect that the root is unreliable because
the derivative of the polynomial vanishes at the root, indicating
ill-conditioning.

Case 2: Nearly Equal Roots
--------------------------

Problem
^^^^^^^

Roots that are very close to each other can lead to loss of numerical
significance.

Issue
^^^^^

Finite precision arithmetic may merge distinct roots or produce
misleading multiplicities.

Reliability Diagnostic
^^^^^^^^^^^^^^^^^^^^^^

Such cases are classified as unreliable due to sensitivity to small
perturbations.

