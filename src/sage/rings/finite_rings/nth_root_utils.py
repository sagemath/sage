sage.doctest: needs sage.rings.finite_rings
r"""
Utilities for :meth:`~sage.rings.finite_rings.element_base.FiniteRingElement.nth_root`.

This module previously contained helpers for intersecting a factorization of
``q-1`` with an integer exponent.  That logic now lives on the
:class:`~sage.structure.factorization.Factorization` class itself via
the :meth:`~sage.structure.factorization.Factorization.gcd_with_exponent`
method and :meth:`~sage.structure.factorization.Factorization.gcd`.

In particular, if ``F`` is a factorization of ``q-1`` then the
factorization of ``\gcd(n_{\text{exp}}, q-1)`` can be obtained with::

    sage: from sage.structure.factorization import Factorization
    sage: F = Factorization([(2, 4)])    # 16
    sage: F.gcd_with_exponent(4, q=17)
    2^2
    sage: F.gcd_with_exponent(8, q=17)
    2^3
    sage: F.gcd_with_exponent(3, q=17)
    1

This keeps the API centralized on :class:`Factorization` instead of providing
a separate helper function.
"""
