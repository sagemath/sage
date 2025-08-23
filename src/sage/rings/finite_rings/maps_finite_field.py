# sage.doctest: needs sage.rings.finite_rings
"""
Structure maps for finite fields

This module provides classes for isomorphisms between finite fields and vector spaces.

AUTHORS:

- Kwankyu Lee (2017-11-07): initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Kwankyu <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism


class FiniteFieldVectorSpaceIsomorphism(Morphism):
    """
    Base class of the vector space isomorphism between a finite field
    and a vector space over a subfield of the finite field.
    """
    def _repr_(self):
        """
        Return the string representation of this isomorphism
        between a finite field and a vector space.

        EXAMPLES::

            sage: E = GF(16)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi
            Isomorphism:
              From: Vector space of dimension 1 over Finite Field in z4 of size 2^4
              To:   Finite Field in z4 of size 2^4
        """
        s = "Isomorphism:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

    def is_injective(self):
        """
        EXAMPLES::

            sage: E = GF(9)
            sage: F = GF(3)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        EXAMPLES::

            sage: E = GF(9)
            sage: F = GF(3)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi.is_surjective()
            True
        """
        return True


class MorphismVectorSpaceToFiniteField(FiniteFieldVectorSpaceIsomorphism):
    """
    Isomorphisms from vector spaces to finite fields.
    """
    def __init__(self, V, K, C):
        """
        Initialize.

        INPUT:

        - ``V`` -- vector space

        - ``K`` -- finite field

        - ``C`` -- matrix

        EXAMPLES::

            sage: E = GF(16)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi
            Isomorphism:
              From: Vector space of dimension 1 over Finite Field in z4 of size 2^4
              To:   Finite Field in z4 of size 2^4
        """
        if C.is_mutable():
            C = C.__copy__()
            C.set_immutable()
        self._C = C
        FiniteFieldVectorSpaceIsomorphism.__init__(self, V, K)

    def _call_(self, v):
        r"""
        TESTS::

            sage: E = GF(64)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(F, map=True)
            sage: phi(V.zero())
            0
            sage: [phi(v) for v in V.basis()]
            [1, z6, z6^2]
        """
        E = self.codomain()  # = GF((p^n)^m)
        V = self.domain()    # = GF(p^n)^m
        m = V.dimension()
        F = V.base_ring()    # = GF(p^n)
        n = F.degree()

        if m == n == 1:
            # 1x1 matrix
            return self._C[0][0] * v[0]
        else:
            # expand v as a vector over GF(p)
            w = self._C.row_ambient_module()()
            for i in range(m):
                w[i*n:(i+1)*n] = v[i]._vector_()
            return E(w * self._C)


class MorphismFiniteFieldToVectorSpace(FiniteFieldVectorSpaceIsomorphism):
    """
    Isomorphisms from finite fields to vector spaces
    """
    def __init__(self, K, V, C):
        """
        Initialize.

        INPUT:

        - ``K`` -- finite field GF((p^m)^n)

        - ``V`` -- vector space of rank n over GF(p^m)

        - ``C`` -- matrix

        EXAMPLES::

            sage: E = GF(16)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: psi
            Isomorphism:
              From: Finite Field in z4 of size 2^4
              To:   Vector space of dimension 1 over Finite Field in z4 of size 2^4
        """
        if C.is_mutable():
            C = C.__copy__()
            C.set_immutable()
        self._C = C
        FiniteFieldVectorSpaceIsomorphism.__init__(self, K, V)

    def _call_(self, e):
        r"""
        TESTS::

            sage: E = GF(64)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(F, map=True)
            sage: psi(E.zero())
            (0, 0, 0)
            sage: psi(E.one())
            (1, 0, 0)
            sage: psi(E.gen())
            (0, 1, 0)
        """
        V = self.codomain()   # = GF(p^n)^m
        m = V.dimension()
        F = V.base_ring()     # = GF(p^n)
        n = F.degree()
        w = e._vector_() * self._C
        if F.degree() > 1:
            return V([F(w[i*n:(i+1)*n]) for i in range(m)])
        else:
            return w
