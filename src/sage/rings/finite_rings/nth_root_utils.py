# sage.doctest: needs sage.rings.finite_rings
r"""
Utilities for :meth:`~sage.rings.finite_rings.element_base.FiniteRingElement.nth_root`.

When ``F`` is a factorization of ``q-1``, the factorization of
``\gcd(n_{\text{exp}}, q-1)`` is obtained with
:meth:`~sage.structure.factorization.Factorization.gcd` applied to
``n_{\text{exp}}.factor()`` (see :issue:`41911`).  The implementation in
:meth:`~sage.rings.finite_rings.element_base.FiniteRingElement._nth_root_common`
also checks that the cached factorization multiplies to ``q-1``.

EXAMPLES::

    sage: from sage.structure.factorization import Factorization
    sage: from sage.rings.integer_ring import ZZ
    sage: F = Factorization([(2, 4)])    # 16 = 2^4, as for q = 17
    sage: F.gcd(ZZ(4).factor())
    2^2
    sage: F.gcd(ZZ(8).factor())
    2^3
    sage: F.gcd(ZZ(3).factor())
    1
"""
