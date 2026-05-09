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
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic


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
    def __init__(self, E1, E2, category=None):
        r"""
        Construct the homset for a given pair of elliptic curves
        defined over the same base ring.

        TESTS::

            sage: E = EllipticCurve(GF(101), [1,1])
            sage: H = End(E)
            sage: TestSuite(H).run(skip='_test_elements')

        ::

            sage: E1 = EllipticCurve(GF(101), [1,1])
            sage: E2 = EllipticCurve(GF(101), [4,9])
            sage: H = Hom(E1, E2)
            sage: TestSuite(H).run(skip='_test_elements')

        ::

            sage: E1 = EllipticCurve(j=42)
            sage: E2 = EllipticCurve(j=43)
            sage: Hom(E1, E2)
            Additive group of elliptic-curve morphisms
              From: Elliptic Curve defined by y^2 = x^3 + 5901*x + 1105454 over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 1510*x - 140675 over Rational Field
            sage: Hom(E1, E2) in CommutativeAdditiveGroups()
            True
            sage: Hom(E1, E2) in Rings()
            False
            sage: End(E1) in CommutativeRings()  # not implemented; see note in code below
            True

        ::

            sage: E0 = EllipticCurve(GF(419), [1,0])
            sage: EE0 = E0.change_ring(GF(419^2))
            sage: End(E0) in CommutativeRings()  # not implemented; see note in code below
            True
            sage: End(EE0) in CommutativeRings()
            False
            sage: End(EE0) in Rings()
            True
            sage: E1 = EllipticCurve(GF(419), [1,1])
            sage: EE1 = E1.change_ring(GF(419^2))
            sage: Hom(E0, E1) in CommutativeAdditiveGroups()
            True
            sage: Hom(EE0, EE1) in CommutativeAdditiveGroups()
            True
            sage: Hom(E0, E1) in Rings()
            False
            sage: Hom(EE0, EE1) in Rings()
            False
        """
        if not isinstance(E1, EllipticCurve_generic):
            raise ValueError('domain must be an elliptic curve')
        if not isinstance(E2, EllipticCurve_generic):
            raise ValueError('codomain must be an elliptic curve')
        base = E1.base_ring()
        if base != E2.base_ring():
            raise ValueError('domain and codomain must have the same base ring')

        super().__init__(E1, E2, category=category, base=base)

        #TODO: We should also add CommutativeRings to the category
        # of self whenever this holds true; see the method
        # EllipticCurve_field.endomorphism_ring_is_commutative().
        # Is there a way to perform this check lazily?

    def _coerce_map_from_(self, other):
        r"""
        Check if this homset has a coercion map from another
        parent ``other``.

        The only currently supported case is when this homset has
        equal domain and codomain: In this case elements from `\ZZ`
        are embedded as scalar multiplications.

        EXAMPLES::

            sage: E1 = EllipticCurve(j=42)
            sage: E2 = EllipticCurve(j=43)
            sage: Hom(E1, E1)._coerce_map_from_(ZZ)
            True
            sage: Hom(E1, E2)._coerce_map_from_(ZZ)
            False
        """
        return self.is_endomorphism_set() and other is ZZ

    def _element_constructor_(self, data):
        r"""
        Construct an element of this homset from the given ``data``.

        The only currently supported case is when this homset has
        equal domain and codomain: In this case elements from `\ZZ`
        are embedded as scalar multiplications.

        EXAMPLES::

            sage: E1 = EllipticCurve(j=42)
            sage: E2 = EllipticCurve(j=43)
            sage: Hom(E1, E1)(5)
            Scalar-multiplication endomorphism [5] of Elliptic Curve defined by y^2 = x^3 + 5901*x + 1105454 over Rational Field
            sage: Hom(E1, E2)(5)
            Traceback (most recent call last):
            ...
            ValueError: domain and codomain must be equal
        """
        m = ZZ(data)
        if not self.is_endomorphism_set():
            if m:
                raise ValueError('domain and codomain must be equal')
            from sage.schemes.elliptic_curves.hom_sum import EllipticCurveHom_sum
            return EllipticCurveHom_sum([], domain=self.domain(), codomain=self.codomain())
        from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
        return EllipticCurveHom_scalar(self.domain(), m)

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
            Additive group of elliptic-curve morphisms
              From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
              To:   Elliptic Curve defined by y^2 = x^3 + 2*x + 2 over Rational Field
        """
        if self.is_endomorphism_set():
            s = 'Ring of elliptic-curve endomorphisms'
        else:
            s = 'Additive group of elliptic-curve morphisms'
        s += f'\n  From: {self.domain()}'
        s += f'\n  To:   {self.codomain()}'
        return s

    def identity(self):
        r"""
        Return the identity morphism in this elliptic-curve homset
        as an :class:`EllipticCurveHom` object.

        EXAMPLES::

            sage: E = EllipticCurve(j=42)
            sage: End(E).identity()
            Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + 5901*x + 1105454 over Rational Field
              Via:  (u,r,s,t) = (1, 0, 0, 0)
            sage: End(E).identity() == E.scalar_multiplication(1)
            True
        """
        if not self.is_endomorphism_set():
            raise ValueError('domain and codomain must be equal')
        return self.domain().identity_morphism()

    def is_commutative(self):
        r"""
        Assuming this homset is an endomorphism ring, check whether
        it is a commutative ring.

        ALGORITHM: :meth:`EllipticCurve_field.endomorphism_ring_is_commutative`

        EXAMPLES::

            sage: End(EllipticCurve(j=123)).is_commutative()
            True
            sage: End(EllipticCurve(GF(11), [1,0])).is_commutative()
            True
            sage: End(EllipticCurve(GF(11^2), [1,0])).is_commutative()
            False
        """
        if not self.is_endomorphism_set():
            raise ValueError('commutativity does not make sense for homsets between different objects')
        return self.domain().endomorphism_ring_is_commutative()
