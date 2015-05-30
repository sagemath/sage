r"""
Structure maps for number fields

Provides isomorphisms between relative and absolute presentations, to and from
vector spaces, name changing maps, etc.

EXAMPLES::

    sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
    sage: K = L.absolute_field('a')
    sage: from_K, to_K = K.structure()
    sage: from_K
    Isomorphism map:
      From: Number Field in a with defining polynomial x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1
      To:   Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
    sage: to_K
    Isomorphism map:
      From: Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
      To:   Number Field in a with defining polynomial x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1
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

from sage.libs.pari.all import pari


QQ = rational_field.RationalField()

IdentityMap = IdentityMorphism

class NumberFieldIsomorphism(Map):
    r"""
    A base class for various isomorphisms between number fields and
    vector spaces.

    EXAMPLE::

        sage: K.<a> = NumberField(x^4 + 3*x + 1)
        sage: V, fr, to = K.vector_space()
        sage: isinstance(fr, sage.rings.number_field.maps.NumberFieldIsomorphism)
        True
    """
    def _repr_type(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^4 + 3*x + 1)
            sage: V, fr, to = K.vector_space()
            sage: fr._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def is_injective(self):
        r"""
         EXAMPLE::

             sage: K.<a> = NumberField(x^4 + 3*x + 1)
             sage: V, fr, to = K.vector_space()
             sage: fr.is_injective()
             True
        """
        return True

    def is_surjective(self):
        r"""
         EXAMPLE::

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
        (Vector space of dimension 4 over Rational Field, Vector space of dimension 4 over Rational Field)
        sage: to.domain(), fr.codomain()
        (Number Field in a with defining polynomial x^4 + 3*x + 1, Number Field in a with defining polynomial x^4 + 3*x + 1)
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

            sage: K.<c> = NumberField(x^9 + 3)
            sage: V, fr, to = K.vector_space(); fr # indirect doctest
            Isomorphism map:
              From: Vector space of dimension 9 over Rational Field
              To:   Number Field in c with defining polynomial x^9 + 3
            sage: type(fr)
            <class 'sage.rings.number_field.maps.MapVectorSpaceToNumberField'>
        """
        self.__V = V
        self.__K = K
        self.__R = K.polynomial_ring()
        NumberFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        r"""
        EXAMPLES::

            sage: K.<c> = NumberField(x^9 + 3)
            sage: V, fr, to = K.vector_space()
            sage: map(fr, V.gens()) # indirect doctest
            [1, c, c^2, c^3, c^4, c^5, c^6, c^7, c^8]
        """
        f = self.__R(self.__V(v).list())
        return self.__K._element_class(self.__K, f)

class MapNumberFieldToVectorSpace(Map):
    r"""
    A class for the isomorphism from an absolute number field to its underlying
    `\QQ`-vector space.

    EXAMPLE::

        sage: L.<a> = NumberField(x^3 - x + 1)
        sage: V, fr, to = L.vector_space()
        sage: type(to)
        <class 'sage.rings.number_field.maps.MapNumberFieldToVectorSpace'>
    """
    def __init__(self, K, V):
        r"""
        Standard initialisation function.

        EXAMPLE::

            sage: L.<a> = NumberField(x^3 - x + 1)
            sage: L.vector_space()[2] # indirect doctest
            Isomorphism map:
              From: Number Field in a with defining polynomial x^3 - x + 1
              To:   Vector space of dimension 3 over Rational Field
        """
        self.__V = V
        self.__K = K
        self.__zero = QQ(0)
        self.__n = K.degree()
        NumberFieldIsomorphism.__init__(self, Hom(K, V))

    def _repr_type(self):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^2 + 1, x^2 - 3])
            sage: V, fr, to = L.relative_vector_space()
            sage: fr._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def _call_(self, x):
        r"""
        EXAMPLE::

            sage: L.<a> = NumberField(x^3 - x + 1)
            sage: V, _, to = L.vector_space()
            sage: v = to(a^2 - a/37 + 56); v # indirect doctest
            (56, -1/37, 1)
            sage: v.parent() is V
            True
        """
        y = self.__K(x)
        v = y._coefficients()
        w = v + [self.__zero]*(self.__n - len(v))
        return self.__V(w)

class MapRelativeVectorSpaceToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES::

        sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
        sage: K = L.relativize(L.subfields(2)[0][1], 'a'); K
        Number Field in a0 with defining polynomial x^2 - b0*x + 1 over its base field
        sage: V, fr, to = K.relative_vector_space()
        sage: V
        Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
          To:   Number Field in a0 with defining polynomial x^2 - b0*x + 1 over its base field
        sage: type(fr)
        <class 'sage.rings.number_field.maps.MapRelativeVectorSpaceToRelativeNumberField'>

        sage: a0 = K.gen(); b0 = K.base_field().gen()
        sage: fr(to(a0 + 2*b0)), fr(V([0, 1])), fr(V([b0, 2*b0]))
        (a0 + 2*b0, a0, 2*b0*a0 + b0)
        sage: (fr * to)(K.gen()) == K.gen()
        True
        sage: (to * fr)(V([1, 2])) == V([1, 2])
        True
    """
    def __init__(self, V, K):
        r"""

        EXAMPLE::

            sage: K.<a, b> = NumberField([x^2 + 1, x^2 - 2])
            sage: V, _, to = K.relative_vector_space(); to # indirect doctest
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 + 1 over its base field
              To:   Vector space of dimension 2 over Number Field in b with defining polynomial x^2 - 2
        """
        self.__V = V
        self.__K = K
        self.__R = K.polynomial_ring()
        self.__rnf = K.pari_rnf()
        self.__B = K.base_field().absolute_field('a')
        NumberFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        r"""
        EXAMPLE::

            sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
            sage: K = L.relativize(L.subfields(2)[0][1], 'a')
            sage: a0 = K.gen(); b0 = K.base_field().gen()
            sage: V, fr, to = K.relative_vector_space()
            sage: fr(to(a0 + 2*b0)), fr(V([0, 1])), fr(V([b0, 2*b0])) # indirect doctest
            (a0 + 2*b0, a0, 2*b0*a0 + b0)
        """
        # Given a relative vector space element, we have to
        # compute the corresponding number field element, in terms
        # of an absolute generator.
        w = self.__V(v).list()

        # First, construct from w a PARI polynomial in x with coefficients
        # that are polynomials in y:
        B = self.__B
        _, to_B = B.structure()
        # Apply to_B, so now each coefficient is in an absolute field,
        # and is expressed in terms of a polynomial in y, then make
        # the PARI poly in x.
        w = [to_B(a)._pari_('y') for a in w]
        h = pari(w).Polrev()

        # Next we write the poly in x over a poly in y in terms
        # of an absolute polynomial for the rnf.
        g = self.__rnf.rnfeltreltoabs(h)
        return self.__K._element_class(self.__K, g)

class MapRelativeNumberFieldToRelativeVectorSpace(NumberFieldIsomorphism):
    r"""
    EXAMPLE::

        sage: K.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
        sage: V, fr, to = K.relative_vector_space()
        sage: type(to)
        <class 'sage.rings.number_field.maps.MapRelativeNumberFieldToRelativeVectorSpace'>
    """

    def __init__(self, K, V):
        r"""
        EXAMPLE::

            sage: L.<b> = NumberField(x^4 + 3*x^2 + 1)
            sage: K = L.relativize(L.subfields(2)[0][1], 'a')
            sage: V, fr, to = K.relative_vector_space()
            sage: to
            Isomorphism map:
              From: Number Field in a0 with defining polynomial x^2 - b0*x + 1 over its base field
              To:   Vector space of dimension 2 over Number Field in b0 with defining polynomial x^2 + 1
        """
        self.__V = V
        self.__K = K
        self.__rnf = K.pari_rnf()
        self.__zero = QQ(0)
        self.__n = K.relative_degree()
        self.__x = pari("'x")
        self.__y = pari("'y")
        self.__B = K.absolute_base_field()
        NumberFieldIsomorphism.__init__(self, Hom(K, V))

    def _call_(self, alpha):
        """
        TESTS::

            sage: K.<a> = NumberField(x^5+2)
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
        # An element of a relative number field is represented
        # internally by an absolute polynomial over QQ.
        alpha = self.__K(alpha)

        # f is the absolute polynomial that defines this number field
        # element
        f = alpha.polynomial('x')
        g = self.__rnf.rnfeltabstorel(pari(f)).lift().lift()
        # Now g is a relative polynomial that defines this element.
        # This g is a polynomial in a pari variable x with
        # coefficients polynomials in a variable y.  These
        # coefficients define the coordinates of the vector we are
        # constructing.

        # The list v below has the coefficients that are the
        # components of the vector we are constructing, but each is
        # converted into polynomials in a variable x, which we will
        # use to define elements of the base field.
        (x, y) = (self.__x, self.__y)
        v = [g.polcoeff(i).subst(x,y) for i in range(self.__n)]
        B,from_B, _ = self.__B
        w = [from_B(B(z)) for z in v]

        # Now w gives the coefficients.
        return self.__V(w)

class NameChangeMap(NumberFieldIsomorphism):
    r"""
    A map between two isomorphic number fields with the same defining
    polynomial but different variable names.

    EXAMPLE::

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
        (<class 'sage.rings.number_field.maps.NameChangeMap'>, <class 'sage.rings.number_field.maps.NameChangeMap'>)
    """
    def __init__(self, K, L):
        r"""
        EXAMPLE::

            sage: K.<a, b> = NumberField([x^2 - 3, x^2 + 7])
            sage: L.<c, d> = K.change_names()
            sage: L.structure()
            (Isomorphism given by variable name change map:
              From: Number Field in c with defining polynomial x^2 - 3 over its base field
              To:   Number Field in a with defining polynomial x^2 - 3 over its base field, Isomorphism given by variable name change map:
              From: Number Field in a with defining polynomial x^2 - 3 over its base field
              To:   Number Field in c with defining polynomial x^2 - 3 over its base field)
        """
        self.__K = K
        self.__L = L
        NumberFieldIsomorphism.__init__(self, Hom(K, L))

    def _repr_type(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^2 - 3)
            sage: L.<b> = K.change_names()
            sage: from_L, to_L = L.structure()
            sage: from_L._repr_type()
            'Isomorphism given by variable name change'
        """
        return "Isomorphism given by variable name change"

    def _call_(self, x):
        r"""
        EXAMPLE::

            sage: K.<a, b> = NumberField([x^2 - 3, x^2 + 7])
            sage: L.<c, d> = K.change_names()
            sage: to_K, from_K = L.structure()
            sage: from_K(a + 17*b) # indirect doctest
            c + 17*d
            sage: to_K(57*c + 19/8*d) # indirect doctest
            57*a + 19/8*b
        """
        y = self.__K(x)
        z = y.__copy__()
        z._set_parent(self.__L)
        return z

class MapRelativeToAbsoluteNumberField(NumberFieldIsomorphism):
    r"""
    EXAMPLES::

        sage: K.<a> = NumberField(x^6 + 4*x^2 + 200)
        sage: L = K.relativize(K.subfields(3)[0][1], 'b'); L
        Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: fr, to = L.structure()
        sage: fr
        Relative number field morphism:
          From: Number Field in b0 with defining polynomial x^2 + a0 over its base field
          To:   Number Field in a with defining polynomial x^6 + 4*x^2 + 200
          Defn: b0 |--> a
                a0 |--> -a^2
        sage: to
        Ring morphism:
          From: Number Field in a with defining polynomial x^6 + 4*x^2 + 200
          To:   Number Field in b0 with defining polynomial x^2 + a0 over its base field
          Defn: a |--> b0
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.morphism.RelativeNumberFieldHomomorphism_from_abs'>,
         <class 'sage.rings.number_field.morphism.NumberFieldHomomorphism_im_gens'>)

        sage: M.<c> = L.absolute_field(); M
        Number Field in c with defining polynomial x^6 + 4*x^2 + 200
        sage: fr, to = M.structure()
        sage: fr
        Isomorphism map:
          From: Number Field in c with defining polynomial x^6 + 4*x^2 + 200
          To:   Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in b0 with defining polynomial x^2 + a0 over its base field
          To:   Number Field in c with defining polynomial x^6 + 4*x^2 + 200
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapAbsoluteToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeToAbsoluteNumberField'>)
        sage: fr(M.gen()), to(fr(M.gen())) == M.gen()
        (b0, True)
        sage: to(L.gen()), fr(to(L.gen())) == L.gen()
        (c, True)
        sage: (to * fr)(M.gen()) == M.gen(), (fr * to)(L.gen()) == L.gen()
        (True, True)
    """

    def __init__(self, R, A):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[1]; f
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 + 3 over its base field
              To:   Number Field in c with defining polynomial x^4 + 16*x^2 + 4
            sage: type(f)
            <class 'sage.rings.number_field.maps.MapRelativeToAbsoluteNumberField'>
        """
        self.__R = R          # relative field
        self.__A = A          # absolute field
        self.__poly_ring = self.__A.polynomial_ring()
        self.__zero = QQ(0)
        self.__n = A.degree()
        NumberFieldIsomorphism.__init__(self, Hom(R, A))

    def _call_(self, x):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[1]
            sage: f(a + 3*b) # indirect doctest
            -c^3 - 17*c
        """
        f = self.__R(x).polynomial()
        return self.__A._element_class(self.__A, f)

class MapAbsoluteToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    See :class:`~MapRelativeToAbsoluteNumberField` for examples.
    """
    def __init__(self, A, R):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[0] # indirect doctest
            sage: type(f)
            <class 'sage.rings.number_field.maps.MapAbsoluteToRelativeNumberField'>
        """
        self.__A = A          # absolute field
        self.__R = R          # relative field
        self.__poly_ring = self.__A.polynomial_ring()
        self.__zero = QQ(0)
        self.__n = A.degree()
        NumberFieldIsomorphism.__init__(self, Hom(A, R))

    def _call_(self, x):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: K.<c> = L.absolute_field()
            sage: f = K.structure()[0]
            sage: f(c + 13*c^2) # indirect doctest
            (-26*b + 1)*a - b - 104
        """

        f = self.__A(x).polynomial()
        return self.__R._element_class(self.__R, f)

class MapVectorSpaceToRelativeNumberField(NumberFieldIsomorphism):
    r"""
    The isomorphism to a relative number field from its underlying `\QQ`-vector
    space. Compare :class:`~MapRelativeVectorSpaceToRelativeNumberField`.

    EXAMPLES::

        sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
        sage: V, fr, to = L.absolute_vector_space()
        sage: type(fr)
        <class 'sage.rings.number_field.maps.MapVectorSpaceToRelativeNumberField'>
    """

    def __init__(self, V, L, from_V, from_K):
        r"""
        EXAMPLES::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space() # indirect doctest
            sage: fr
            Isomorphism map:
              From: Vector space of dimension 4 over Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3 over its base field
        """
        self.__V = V
        self.__L = L
        self.__from_V = from_V
        self.__from_K = from_K
        NumberFieldIsomorphism.__init__(self, Hom(V, L))

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space()
            sage: fr(V([1,3,0,1/17])) # indirect doctest
            33/17*a - 37/17*b + 1
            sage: fr(to(a)), fr(to(b)) # indirect doctest
            (a, b)
        """
        return self.__from_K(self.__from_V(x))

class MapRelativeNumberFieldToVectorSpace(NumberFieldIsomorphism):
    r"""
    The isomorphism from a relative number field to its underlying `\QQ`-vector
    space. Compare :class:`~MapRelativeNumberFieldToRelativeVectorSpace`.

    EXAMPLES::

        sage: K.<a> = NumberField(x^8 + 100*x^6 + x^2 + 5)
        sage: L = K.relativize(K.subfields(4)[0][1], 'b'); L
        Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: L_to_K, K_to_L = L.structure()

        sage: V, fr, to = L.absolute_vector_space()
        sage: V
        Vector space of dimension 8 over Rational Field
        sage: fr
        Isomorphism map:
          From: Vector space of dimension 8 over Rational Field
          To:   Number Field in b0 with defining polynomial x^2 + a0 over its base field
        sage: to
        Isomorphism map:
          From: Number Field in b0 with defining polynomial x^2 + a0 over its base field
          To:   Vector space of dimension 8 over Rational Field
        sage: type(fr), type(to)
        (<class 'sage.rings.number_field.maps.MapVectorSpaceToRelativeNumberField'>,
         <class 'sage.rings.number_field.maps.MapRelativeNumberFieldToVectorSpace'>)

        sage: v = V([1, 1, 1, 1, 0, 1, 1, 1])
        sage: fr(v), to(fr(v)) == v
        ((-a0^3 + a0^2 - a0 + 1)*b0 - a0^3 - a0 + 1, True)
        sage: to(L.gen()), fr(to(L.gen())) == L.gen()
        ((0, 1, 0, 0, 0, 0, 0, 0), True)
    """
    def __init__(self, L, V, to_K, to_V):
        r"""
        EXAMPLES::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space() # indirect doctest
            sage: to
            Isomorphism map:
              From: Number Field in a with defining polynomial x^2 + 3 over its base field
              To:   Vector space of dimension 4 over Rational Field
        """
        self.__L = L
        self.__V = V
        self.__to_K = to_K
        self.__to_V = to_V
        NumberFieldIsomorphism.__init__(self, Hom(L, V))

    def _call_(self, x):
        r"""
        EXAMPLES::

            sage: L.<a, b> = NumberField([x^2 + 3, x^2 + 5])
            sage: V, fr, to = L.absolute_vector_space()
            sage: to(1 + 2*a + 3*b + 4*a*b) # indirect doctest
            (-15, -41/2, -2, -5/4)
            sage: to(fr(V([1,3,0,1/17]))) # indirect doctest
            (1, 3, 0, 1/17)
        """
        return self.__to_V(self.__to_K(x))
