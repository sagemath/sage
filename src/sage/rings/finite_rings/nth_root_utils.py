# sage.doctest: needs sage.rings.finite_rings
r"""
Utilities for :meth:`~sage.rings.finite_rings.element_base.FiniteRingElement.nth_root`.

When ``F`` is a factorization of ``q-1`` and ``n = \gcd(n_{\text{exp}}, q-1)``,
the factorization of ``n`` is obtained with
:meth:`~sage.structure.factorization.Factorization.gcd` applied to
``n.factor()`` where ``n = \gcd(n_{\text{exp}}, q-1)`` (not
``n_{\text{exp}}.factor()``), matching the classical reduction before factoring
(see :issue:`41911`).  The implementation
in :meth:`~sage.rings.finite_rings.element_base.FiniteRingElement._nth_root_common`
also checks that the cached factorization multiplies to ``q-1``.

EXAMPLES::

    sage: from sage.structure.factorization import Factorization
    sage: from sage.rings.integer_ring import ZZ
    sage: F = Factorization([(2, 4)])    # 16 = 2^4, as for q = 17
    sage: F.gcd(ZZ(4).factor())   # gcd(4, 16) = 4 -> 2^2
    2^2
    sage: F.gcd(ZZ(8).factor())   # gcd(8, 16) = 8 -> 2^3
    2^3
    sage: F.gcd(ZZ(3).factor())   # gcd(3, 16) = 1
    1
"""
