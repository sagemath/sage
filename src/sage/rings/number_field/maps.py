r"""
Structure maps for number fields

This module provides isomorphisms between relative and absolute presentations,
to and from vector spaces, name changing maps, etc.

EXAMPLES::

    sage: x = polygen(ZZ, 'x')
    sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
    sage: K = L.absolute_field('a')
    sage: from_K, to_K = K.structure()
    sage: from_K
    Isomorphism map:
      From: Number Field in a with defining polynomial
            x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1
      To:   Number Field in cuberoot2 with defining polynomial
            x^3 - 2 over its base field
    sage: to_K
    Isomorphism map:
      From: Number Field in cuberoot2 with defining polynomial
            x^3 - 2 over its base field
      To:   Number Field in a with defining polynomial
            x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.map import Map
from sage.categories.homset import Hom
from sage.categories.morphism import IdentityMorphism

import sage.rings.rational_field as rational_field

QQ = rational_field.RationalField()

IdentityMap = IdentityMorphism


class NumberFieldIsomorphism(Map):
    r"""
    A base class for various isomorphisms between number fields and
    vector spaces.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^4 + 3*x + 1)
        sage: V, fr, to = K.vector_space()
        sage: isinstance(fr, sage.rings.number_field.maps.NumberFieldIsomorphism)
        True
    """
    def _repr_type(self):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^4 + 3*x + 1)
            sage: V, fr, to = K.vector_space()
            sage: fr._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def is_injective(self):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^4 + 3*x + 1)
            sage: V, fr, to = K.vector_space()
            sage: fr.is_injective()
            True
        """
        return True

    def is_surjective(self):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^4 + 3*x + 1)
            sage: V, fr, to = K.vector_space()
            sage: fr.is_surjective()
            True
        """
        return True


class MapVectorSpaceToNumberField(NumberFieldIsomorphism):
    r"""
    The map to an absolute number field from its underlying `\QQ`-vector space.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^4 + 3*x + 1)
        sage: V, fr, to = K.vector_space()
        sage: V
        Vector space of dimension 4 over Rational Field
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 4 over Rational Field
          To:   Number Field in a with defining polynomial x^4 + 3*x + 1
        sage: to
        Isomorphism map:
          From: Number Field in a with defining polynomial x^4 + 3*x + 1
          To:   Vector space of dimension 4 over Rational Field
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapVectorSpaceToNumberField'>,
         <class 'sage.rings.number_field.maps.MapNumberFieldToVectorSpace'>)

        sage: fr.is_injective(), fr.is_surjective()
        (True, True)

        sage: fr.domain(), to.codomain()
        (Vector space of dimension 4 over Rational Field,
         Vector space of dimension 4 over Rational Field)
        sage: to.domain(), fr.codomain()
        (Number Field in a with defining polynomial x^4 + 3*x + 1,
         Number Field in a with defining polynomial x^4 + 3*x + 1)
        sage: fr * to
        Composite map:
          From: Number Field in a with defining polynomial x^4 + 3*x + 1
          To:   Number Field in a with defining polynomial x^4 + 3*x + 1
          Defn:   Isomorphism map:
                  From: Number Field in a with defining polynomial x^4 + 3*x + 1
                  To:   Vector space of dimension 4 over Rational Field
                then
                  Isomorphism map:
                  From: Vector space of dimension 4 over Rational Field
                  To:   Number Field in a with defining polynomial x^4 + 3*x + 1
        sage: to * fr
        Composite map:
          From: Vector space of dimension 4 over Rational Field
          To:   Vector space of dimension 4 over Rational Field
          Defn:   Isomorphism map:
                  From: Vector space of dimension 4 over Rational Field
                  To:   Number Field in a with defining polynomial x^4 + 3*x + 1
                then
                  Isomorphism map:
                  From: Number Field in a with defining polynomial x^4 + 3*x + 1
                  To:   Vector space of dimension 4 over Rational Field

        sage: to(a), to(a + 1)
        ((0, 1, 0, 0), (1, 1, 0, 0))
        sage: fr(to(a)), fr(V([0, 1, 2, 3]))
        (a, 3*a^3 + 2*a^2 + a)
    """

    def __init__(self, V, K):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<c> = NumberField(x^9 + 3)
            sage: V, fr, to = K.vector_space(); fr # indirect doctest
            Isomorphism map:
              From: Vector space of dimension 9 over Rational Field
              To:   Number Field in c with defining polynomial x^9 + 3
            sage: type(fr)
            <class 'sage.rings.number_field.maps.MapVectorSpaceToNumberField'>
        """
        NumberFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<c> = NumberField(x^9 + 3)
            sage: V, fr, to = K.vector_space()
            sage: list(map(fr, V.gens()))  # indirect doctest
            [1, c, c^2, c^3, c^4, c^5, c^6, c^7, c^8]
        """
        K = self.codomain()
        f = K.polynomial_ring()(v.list())
        return K._element_class(K, f)


class MapNumberFieldToVectorSpace(Map):
    r"""
    A class for the isomorphism from an absolute number field to its underlying
    `\QQ`-vector space.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: L.<a> = NumberField(x^3 - x + 1)
        sage: V, fr, to = L.vector_space()
        sage: type(to)
        <class 'sage.rings.number_field.maps.MapNumberFieldToVectorSpace'>
    """
    def __init__(self, K, V):
        r"""
        Standard initialisation function.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a> = NumberField(x^3 - x + 1)
            sage: L.vector_space()[2] # indirect doctest
            Isomorphism map:
              From: Number Field in a with defining polynomial x^3 - x + 1
              To:   Vector space of dimension 3 over Rational Field
        """
        NumberFieldIsomorphism.__init__(self, Hom(K, V))

    def _repr_type(self):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 1, x^2 - 3])
            sage: V, fr, to = L.relative_vector_space()
            sage: fr._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a> = NumberField(x^3 - x + 1)
            sage: V, _, to = L.vector_space()
            sage: v = to(a^2 - a/37 + 56); v # indirect doctest
            (56, -1/37, 1)
            sage: v.parent() is V
            True
        """
        v = x._coefficients()
        k = self.domain().degree() - len(v)
        if k > 0:
            v = v + [QQ.zero()] * k
        return self.codomain()(v)


class MapRelativeVectorSpaceToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
        sage: K = L.relativize(L.subfields(2)[0][1], 'a'); K
        Number Field in a with defining polynomial x^2 - b0*x + 1 over its base field
        sage: V, fr, to = K.relative_vector_space()
        sage: V
        Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 2
                over Number Field in b0 with defining polynomial x^2 + 1
          To:   Number Field in a
                with defining polynomial x^2 - b0*x + 1 over its base field
        sage: type(fr)
        <class 'sage.rings.number_field.maps.MapRelativeVectorSpaceToRelativeNumberField'>

        sage: a0 = K.gen(); b0 = K.base_field().gen()
        sage: fr(to(a0 + 2*b0)), fr(V([0, 1])), fr(V([b0, 2*b0]))
        (a + 2*b0, a, 2*b0*a + b0)
        sage: (fr * to)(K.gen()) == K.gen()
        True
        sage: (to * fr)(V([1, 2])) == V([1, 2])
        True
    """
    def __init__(self, V, K):
        r"""

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^2 + 1, x^2 - 2])
            sage: V, _, to = K.relative_vector_space(); to # indirect doctest
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 + 1 over its base field
              To:   Vector space of dimension 2 over Number Field in b with defining polynomial x^2 - 2
        """
        NumberFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
            sage: K = L.relativize(L.subfields(2)[0][1], 'a')
            sage: a0 = K.gen(); b0 = K.base_field().gen()
            sage: V, fr, to = K.relative_vector_space()
            sage: fr(to(a0 + 2*b0)), fr(V([0, 1])), fr(V([b0, 2*b0])) # indirect doctest
            (a + 2*b0, a, 2*b0*a + b0)
        """
        from sage.libs.pari.all import pari
        
        K = self.codomain()
        B = K.base_field().absolute_field('a')
        # Convert v to a PARI polynomial in x with coefficients that
        # are polynomials in y.
        _, to_B = B.structure()
        h = pari([to_B(a).__pari__('y') for a in v]).Polrev()
        # Rewrite the polynomial in terms of an absolute generator for
        # the relative number field.
        g = K._pari_rnfeq()._eltreltoabs(h)
        return K._element_class(K, g)


class MapRelativeNumberFieldToRelativeVectorSpace(NumberFieldIsomorphism):
    r"""
    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
        sage: V, fr, to = K.relative_vector_space()
        sage: type(to)
        <class 'sage.rings.number_field.maps.MapRelativeNumberFieldToRelativeVectorSpace'>
    """

    def __init__(self, K, V):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
            sage: K = L.relativize(L.subfields(2)[0][1], 'a')
            sage: V, fr, to = K.relative_vector_space()
            sage: to
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 - b0*x + 1 over its base field
              To:   Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
        """
        NumberFieldIsomorphism.__init__(self, Hom(K, V))

    def _call_(self, alpha):
        """
        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^5 + 2)
            sage: R.<y> = K[]
            sage: D.<x0> = K.extension(y + a + 1)
            sage: D(a)
            a
            sage: V, from_V, to_V = D.relative_vector_space()
            sage: to_V(a) # indirect doctest
            (a)
            sage: to_V(a^3) # indirect doctest
            (a^3)
            sage: to_V(x0) # indirect doctest
            (-a - 1)

            sage: K.<a> = QuadraticField(-3)
            sage: L.<b> = K.extension(x-5)
            sage: L(a)
            a
            sage: a*b
            5*a
            sage: b
            5
            sage: V, from_V, to_V = L.relative_vector_space()
            sage: to_V(a) # indirect doctest
            (a)
        """
        K = self.domain()
        # The element alpha is represented internally by an absolute
        # polynomial over QQ, and f is its PARI representation.
        f = alpha._pari_polynomial('x')
        # Convert f to a relative polynomial g; this is a polynomial
        # in x whose coefficients are polynomials in y.
        g = K._pari_rnfeq()._eltabstorel_lift(f)
        # Now g is a polynomial in the standard generator of the PARI
        # field; convert it to a polynomial in the Sage generator.
        if g.poldegree() > 0:
            beta = K._pari_relative_structure()[2]
            g = g(beta).lift()
        # Convert the coefficients to elements of the base field.
        B, from_B, _ = K.absolute_base_field()
        return self.codomain()([from_B(B(z.lift(), check=False)) for z in g.Vecrev(K.relative_degree())])


class NameChangeMap(NumberFieldIsomorphism):
    r"""
    A map between two isomorphic number fields with the same defining
    polynomial but different variable names.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 3)
        sage: L.<b> = K.change_names()
        sage: from_L, to_L = L.structure()
        sage: from_L
        Isomorphism given by variable name change map:
          From: Number Field in b with defining polynomial x^2 - 3
          To:   Number Field in a with defining polynomial x^2 - 3
        sage: to_L
        Isomorphism given by variable name change map:
          From: Number Field in a with defining polynomial x^2 - 3
          To:   Number Field in b with defining polynomial x^2 - 3
        sage: type(from_L), type(to_L)
        (<class 'sage.rings.number_field.maps.NameChangeMap'>,
         <class 'sage.rings.number_field.maps.NameChangeMap'>)
    """
    def __init__(self, K, L):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^2 - 3, x^2 + 7])
            sage: L.<c, d> = K.change_names()
            sage: L.structure()
            (Isomorphism given by variable name change map:
              From: Number Field in c with defining polynomial x^2 - 3 over its base field
              To:   Number Field in a with defining polynomial x^2 - 3 over its base field,
             Isomorphism given by variable name change map:
              From: Number Field in a with defining polynomial x^2 - 3 over its base field
              To:   Number Field in c with defining polynomial x^2 - 3 over its base field)
        """
        NumberFieldIsomorphism.__init__(self, Hom(K, L))

    def _repr_type(self):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 3)
            sage: L.<b> = K.change_names()
            sage: from_L, to_L = L.structure()
            sage: from_L._repr_type()
            'Isomorphism given by variable name change'
        """
        return "Isomorphism given by variable name change"

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a, b> = NumberField([x^2 - 3, x^2 + 7])
            sage: L.<c, d> = K.change_names()
            sage: to_K, from_K = L.structure()
            sage: from_K(a + 17*b) # indirect doctest
            c + 17*d
            sage: to_K(57*c + 19/8*d) # indirect doctest
            57*a + 19/8*b
        """
        y = x._copy_for_parent(self.codomain())
        return y


class MapRelativeToAbsoluteNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^6 + 4*x^2 + 200)
        sage: L = K.relativize(K.subfields(3)[0][1], 'b'); L
        Number Field in b with defining polynomial x^2 + a0 over its base field
        sage: fr, to = L.structure()
        sage: fr
        Relative number field morphism:
          From: Number Field in b with defining polynomial x^2 + a0 over its base field
          To:   Number Field in a with defining polynomial x^6 + 4*x^2 + 200
          Defn: b |--> a
                a0 |--> -a^2
        sage: to
        Ring morphism:
          From: Number Field in a with defining polynomial x^6 + 4*x^2 + 200
          To:   Number Field in b with defining polynomial x^2 + a0 over its base field
          Defn: a |--> b
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.homset.RelativeNumberFieldHomset_with_category.element_class'>,
         <class 'sage.rings.number_field.homset.NumberFieldHomset_with_category.element_class'>)

        sage: M.<c> = L.absolute_field(); M
        Number Field in c with defining polynomial x^6 + 4*x^2 + 200
        sage: fr, to = M.structure()
        sage: fr
        Isomorphism map:
          From: Number Field in c with defining polynomial x^6 + 4*x^2 + 200
          To:   Number Field in b with defining polynomial x^2 + a0 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in b with defining polynomial x^2 + a0 over its base field
          To:   Number Field in c with defining polynomial x^6 + 4*x^2 + 200
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapAbsoluteToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeToAbsoluteNumberField'>)
        sage: fr(M.gen()), to(fr(M.gen())) == M.gen()
        (b, True)
        sage: to(L.gen()), fr(to(L.gen())) == L.gen()
        (c, True)
        sage: (to * fr)(M.gen()) == M.gen(), (fr * to)(L.gen()) == L.gen()
        (True, True)
    """

    def __init__(self, R, A):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[1]; f
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 + 3 over its base field
              To:   Number Field in c with defining polynomial x^4 + 16*x^2 + 4
            sage: type(f)
            <class 'sage.rings.number_field.maps.MapRelativeToAbsoluteNumberField'>
        """
        NumberFieldIsomorphism.__init__(self, Hom(R, A))

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[1]
            sage: f(a + 3*b) # indirect doctest
            -c^3 - 17*c
        """
        A = self.codomain()  # absolute field
        f = x.polynomial()
        return A._element_class(A, f)


class MapAbsoluteToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    See :class:`~MapRelativeToAbsoluteNumberField` for examples.
    """
    def __init__(self, A, R):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[0] # indirect doctest
            sage: type(f)
            <class 'sage.rings.number_field.maps.MapAbsoluteToRelativeNumberField'>
        """
        NumberFieldIsomorphism.__init__(self, Hom(A, R))

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[0]
            sage: f(c + 13*c^2) # indirect doctest
            (-26*b + 1)*a - b - 104
        """
        R = self.codomain()  # relative field
        f = x.polynomial()
        return R._element_class(R, f)


class MapVectorSpaceToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    The isomorphism to a relative number field from its underlying `\QQ`-vector
    space. Compare :class:`~MapRelativeVectorSpaceToRelativeNumberField`.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
        sage: V, fr, to = L.absolute_vector_space()
        sage: type(fr)
        <class 'sage.rings.number_field.maps.MapVectorSpaceToRelativeNumberField'>
    """

    def __init__(self, V, L, from_V, from_K):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space() # indirect doctest
            sage: fr
            Isomorphism map:
              From: Vector space of dimension 4 over Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3 over its base field
        """
        self._from_V = from_V
        self._from_K = from_K
        NumberFieldIsomorphism.__init__(self, Hom(V, L))

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space()
            sage: fr(V([1,3,0,1/17])) # indirect doctest
            33/17*a - 37/17*b + 1
            sage: fr(to(a)), fr(to(b)) # indirect doctest
            (a, b)
        """
        return self._from_K(self._from_V(x))


class MapRelativeNumberFieldToVectorSpace(NumberFieldIsomorphism):
    r"""
    The isomorphism from a relative number field to its underlying `\QQ`-vector
    space. Compare :class:`~MapRelativeNumberFieldToRelativeVectorSpace`.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^8 + 100*x^6 + x^2 + 5)
        sage: L = K.relativize(K.subfields(4)[0][1], 'b'); L
        Number Field in b with defining polynomial x^2 + a0 over its base field
        sage: L_to_K, K_to_L = L.structure()

        sage: V, fr, to = L.absolute_vector_space()
        sage: V
        Vector space of dimension 8 over Rational Field
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 8 over Rational Field
          To:   Number Field in b with defining polynomial x^2 + a0 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in b with defining polynomial x^2 + a0 over its base field
          To:   Vector space of dimension 8 over Rational Field
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapVectorSpaceToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeNumberFieldToVectorSpace'>)

        sage: v = V([1, 1, 1, 1, 0, 1, 1, 1])
        sage: fr(v), to(fr(v)) == v
        ((-a0^3 + a0^2 - a0 + 1)*b - a0^3 - a0 + 1, True)
        sage: to(L.gen()), fr(to(L.gen())) == L.gen()
        ((0, 1, 0, 0, 0, 0, 0, 0), True)
    """
    def __init__(self, L, V, to_K, to_V):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space() # indirect doctest
            sage: to
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 + 3 over its base field
              To:   Vector space of dimension 4 over Rational Field
        """
        self._to_K = to_K
        self._to_V = to_V
        NumberFieldIsomorphism.__init__(self, Hom(L, V))

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space()
            sage: to(1 + 2*a + 3*b + 4*a*b) # indirect doctest
            (-15, -41/2, -2, -5/4)
            sage: to(fr(V([1,3,0,1/17]))) # indirect doctest
            (1, 3, 0, 1/17)
        """
        return self._to_V(self._to_K(x))
