r"""
Homsets and endomorphism rings of elliptic curves

The set of homomorphisms between two elliptic curves (:class:`EllipticCurveHom`)
forms an abelian group under addition. Moreover, if the two curves are the same,
it even forms a (not always commutative) ring under composition.

This module encapsulates the set of homomorphisms between two given elliptic
curves as a Sage object.

.. NOTE::

    Currently only little nontrivial functionality is available, but this will
    hopefully change in the future.

EXAMPLES:

The only useful thing this class does at the moment is coercing integers into
the endomorphism ring as scalar multiplications::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: f = End(E)(7); f
    Scalar-multiplication endomorphism [7] of Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field
    sage: f == E.scalar_multiplication(7)
    True

::

    sage: E = EllipticCurve(GF(431^2), [0,1])
    sage: E.automorphisms()[0] == 1
    True
    sage: E.automorphisms()[1] == -1
    True
    sage: omega = E.automorphisms()[2]
    sage: omega == 1
    False
    sage: omega^3 == 1
    True
    sage: (1 + omega + omega^2) == 0
    True
    sage: (2*omega + 1)^2 == -3
    True

AUTHORS:

- Lorenz Panny (2023)
"""

# ****************************************************************************
#       Copyright (C) 2023 Lorenz Panny
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.categories.morphism import Morphism
from sage.schemes.generic.homset import SchemeHomset_generic


class EllipticCurveHomset(SchemeHomset_generic):
    r"""
    This class represents the set of all homomorphisms between two fixed
    elliptic curves.

    EXAMPLES::

        sage: E = EllipticCurve(GF(419^2), [1,0])
        sage: E.frobenius_isogeny() in End(E)
        True
        sage: phi = E.isogenies_prime_degree(7)[0]
        sage: phi in End(E)
        False
        sage: phi in Hom(E, phi.codomain())
        True

    Note that domain and codomain are *not* taken up to isomorphism::

        sage: iso = E.isomorphism_to(EllipticCurve(GF(419^2), [2,0]))
        sage: iso in End(E)
        False
    """
    def __init__(self, *args, **kwds):
        r"""
        Construct the homset for a given pair of curves.

        TESTS::

            sage: E1 = EllipticCurve(j=42)
            sage: E2 = EllipticCurve(j=43)
            sage: Hom(E1, E2)
            Group of elliptic-curve morphisms
              From: Elliptic Curve defined by y^2 = x^3 + 5901*x + 1105454 over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 1510*x - 140675 over Rational Field
        """
        super().__init__(*args, **kwds)

        if self.domain() == self.codomain():
            # set up automated coercion of integers to scalar multiplications
            from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
            class ScalarMultiplicationEmbedding(Morphism):
                def __init__(self, End):
                    assert End.domain() is End.codomain()
                    super().__init__(ZZ, End)
                def _call_(self, m):
                    return EllipticCurveHom_scalar(self.codomain().domain(), m)
            self.register_coercion(ScalarMultiplicationEmbedding(self))

    def _repr_(self):
        r"""
        Output a description of this homset, with special formatting
        for endomorphism rings.

        EXAMPLES::

            sage: E1 = EllipticCurve([1,1])
            sage: E2 = EllipticCurve([2,2])
            sage: End(E1)
            Ring of elliptic-curve endomorphisms
              From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
              To:   Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
            sage: Hom(E1, E1)
            Ring of elliptic-curve endomorphisms
              From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
              To:   Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
            sage: Hom(E1, E2)
            Group of elliptic-curve morphisms
              From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
              To:   Elliptic Curve defined by y^2 = x^3 + 2*x + 2 over Rational Field
        """
        if self.domain() == self.codomain():
            s = 'Ring of elliptic-curve endomorphisms'
        else:
            s = 'Group of elliptic-curve morphisms'
        s += f'\n  From: {self.domain()}'
        s += f'\n  To:   {self.codomain()}'
        return s

