# sage_setup: distribution = sagemath-flint
r"""
Morphisms between number fields

This module provides classes to represent ring homomorphisms between number
fields (i.e. field embeddings).
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.rings.morphism import RingHomomorphism_im_gens, RingHomomorphism
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp


class NumberFieldHomomorphism_im_gens(RingHomomorphism_im_gens):
    def __invert__(self):
        r"""
        Return the inverse of an isomorphism of absolute number fields.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 5)
            sage: tau1, tau2 = K.automorphisms(); tau1, tau2
            (Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
              Defn: a |--> a,
             Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
              Defn: a |--> -a)
            sage: ~tau1
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
             Defn: a |--> a
            sage: ~tau2
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
             Defn: a |--> -a

            sage: L.<z> = CyclotomicField(5)
            sage: tau1, tau2, tau3, tau4 = L.automorphisms()
            sage: (tau1, ~tau1)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z)
            sage: (tau2, ~tau2)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^2,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^3)
            sage: (tau3, ~tau3)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^3,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^2)

             sage: M.<w> = NumberField(x^4 - 5*x + 5)
             sage: phi = M.hom([z - z^2]); phi
             Ring morphism:
               From: Number Field in w with defining polynomial x^4 - 5*x + 5
               To:   Cyclotomic Field of order 5 and degree 4
               Defn: w |--> -z^2 + z
             sage: phi^-1
             Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Number Field in w with defining polynomial x^4 - 5*x + 5
               Defn: z |--> 3/11*w^3 + 4/11*w^2 + 9/11*w - 14/11
        """
        K = self.domain()
        L = self.codomain()
        if K.degree() != L.degree():
            raise TypeError("Can only invert isomorphisms")
        V, V_into_K, _ = K.vector_space()
        _, _, L_into_W = L.vector_space()
        linear_inverse = ~V.hom([(L_into_W * self * V_into_K)(b)
                                 for b in V.basis()])
        return L.hom([(V_into_K * linear_inverse * L_into_W)(b)
                      for b in [L.gen()]])

    def preimage(self, y):
        r"""
        Compute a preimage of `y` in the domain, provided one exists.
        Raises a :exc:`ValueError` if `y` has no preimage.

        INPUT:

        - ``y`` -- an element of the codomain of ``self``

        OUTPUT:

        Returns the preimage of `y` in the domain, if one exists.
        Raises a :exc:`ValueError` if `y` has no preimage.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 7)
            sage: L.<b> = NumberField(x^4 - 7)
            sage: f = K.embeddings(L)[0]
            sage: f.preimage(3*b^2 - 12/7)
            3*a - 12/7
            sage: f.preimage(b)
            Traceback (most recent call last):
            ...
            ValueError: Element 'b' is not in the image of this homomorphism.

        ::

            sage: # needs sage.libs.linbox
            sage: F.<b> = QuadraticField(23)
            sage: G.<a> = F.extension(x^3 + 5)
            sage: f = F.embeddings(G)[0]
            sage: f.preimage(a^3 + 2*b + 3)
            2*b - 2
        """
        # Throughout this method I am using the convention that self is a homomorphism from the number field K to the number field L
        # Therefore, I use the names K and L in place of domain and codomain

        # try to get the cached transformation matrix and vector space isomorphisms if they exist
        try:
            M,LtoV,VtoK = self._transformation_data
        except Exception:
            # get the identifications of K and L with vector spaces over Q
            V,VtoL,LtoV = self.codomain().absolute_vector_space()
            V,VtoK,KtoV = self.domain().absolute_vector_space()
            # construct the transformation matrix from K to L by making the columns be the image of the basis of V_K in V_L using the homomorphism
            from sage.matrix.constructor import matrix
            from sage.rings.rational_field import QQ
            M = matrix(QQ, [LtoV(self(VtoK(e))) for e in V.basis()]).transpose()
            self._transformation_data = (M,LtoV,VtoK)

        # get the coordinate vector of y, solve the linear system, pass to domain
        yvec = LtoV(y)                  # pass from a point in L to its vector space representation
        try:
            xvec = M.solve_right(yvec)      # solve the linear system, throws an exception if there is no solution
        except ValueError:
            raise ValueError("Element '{}' is not in the image of this homomorphism.".format(y))
        return VtoK(xvec)               # pass from the vector space representation of K back to a point in K


class RelativeNumberFieldHomomorphism_from_abs(RingHomomorphism):
    r"""
    A homomorphism from a relative number field to some other ring, stored as a
    homomorphism from the corresponding absolute field.
    """

    def __init__(self, parent, abs_hom):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^3 + 2, x^2 + x + 1])
            sage: f = K.hom(-a*b - a, K); f
            Relative number field endomorphism of
             Number Field in a with defining polynomial x^3 + 2 over its base field
              Defn: a |--> (-b - 1)*a
                    b |--> b
            sage: type(f)
            <class 'sage.rings.number_field.homset.RelativeNumberFieldHomset_with_category.element_class'>
        """
        RingHomomorphism.__init__(self, parent)
        self._abs_hom = abs_hom
        K = abs_hom.domain()
        from_K, to_K = K.structure()
        self._from_K = from_K
        self._to_K = to_K

    def abs_hom(self):
        r"""
        Return the corresponding homomorphism from the absolute number field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^3 + 2, x^2 + x + 1])
            sage: K.hom(a, K).abs_hom()
            Ring morphism:
              From: Number Field in a with defining polynomial
                    x^6 - 3*x^5 + 6*x^4 - 3*x^3 - 9*x + 9
              To:   Number Field in a with defining polynomial x^3 + 2 over its base field
              Defn: a |--> a - b
        """
        return self._abs_hom

    def _repr_type(self):
        r"""
        A short string to identify the type of this homomorphism.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^3 + 2, x^2 + x + 1])
            sage: K.hom(a, K)._repr_type()
            'Relative number field'
        """
        return "Relative number field"

    @cached_method
    def im_gens(self):
        r"""
        Return the images of the generators under this map.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^3 + 2, x^2 + x + 1])
            sage: K.hom(a, K).im_gens()
            [a, b]
        """
        D = self.domain()
        C = self.codomain()
        return Sequence([self(x) for x in D.gens()], universe=C, check=False, immutable=True)

    def _richcmp_(self, other, op):
        """
        Compare.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^2 - 2, x^2 - 3])
            sage: e, u, v, w = End(K)
            sage: all([u^2 == e, u*v == w, u != e])
            True
        """
        return richcmp(self.abs_hom(), other.abs_hom(), op)

    def _repr_defn(self):
        r"""
        Return a string describing the images of the generators under this map.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^3 + 2, x^2 + x + 1])
            sage: K.hom(a, K)._repr_defn()
            'a |--> a\nb |--> b'
        """
        D = self.domain()
        ig = self.im_gens()
        return '\n'.join('%s |--> %s' % (D.gen(i), ig[i])
                         for i in range(D.ngens()))

    def _call_(self, x):
        r"""
        Evaluate this map at the element ``x``.

        This is done by first
        converting ``x`` to an element of the absolute field and then
        evaluating ``self.abs_hom()`` on it.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^3 + 2, x^2 + x + 1])
            sage: K.hom(a*b, K)(17 + 3*a + 2*b) # indirect doctest
            3*b*a + 2*b + 17
        """
        return self._abs_hom(self._to_K(x))


class CyclotomicFieldHomomorphism_im_gens(NumberFieldHomomorphism_im_gens):
    pass
