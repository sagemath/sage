# -*- coding: utf-8 -*-
r"""
Morphisms of function fields

Maps and morphisms useful for computations with function fields.

EXAMPLES::

    sage: K.<x> = FunctionField(QQ); R.<y> = K[]
    sage: K.hom(1/x)
    Function Field endomorphism of Rational function field in x over Rational Field
      Defn: x |--> 1/x
    sage: L.<y> = K.extension(y^2 - x)
    sage: K.hom(y)
    Function Field morphism:
      From: Rational function field in x over Rational Field
      To:   Function field in y defined by y^2 - x
      Defn: x |--> y
    sage: L.hom([y,x])
    Function Field endomorphism of Function field in y defined by y^2 - x
      Defn: y |--> y
            x |--> x
    sage: L.hom([x,y])
    Traceback (most recent call last):
    ...
    ValueError: invalid morphism

AUTHORS:

- William Stein (2010): initial version

- Julian Rüth (2011-09-14, 2014-06-23, 2017-08-21): refactored class hierarchy; added
  derivation classes; morphisms to/from fraction fields

"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011-2017 Julian Rüth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from sage.categories.map import Map
from sage.rings.morphism import RingHomomorphism

class FunctionFieldDerivation(Map):
    r"""
    Base class for derivations on function fields.

    A derivation on `R` is a map `R \to R` with
    `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
    D(\alpha)+\alpha D(\beta)` for all `\alpha,\beta\in R`.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: d
        Derivation map:
          From: Rational function field in x over Rational Field
          To:   Rational function field in x over Rational Field
    """
    def __init__(self, K):
        r"""
        Initialize a derivation from `K` to `K`.

        INPUT:

        - ``K`` -- function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: TestSuite(d).run(skip=['_test_category', '_test_pickling'])

        .. TODO::

            Make the caching done at the map by subclassing
            ``UniqueRepresentation``, which will then implement a
            valid equality check. Then this will pass the pickling test.
        """
        from .function_field import is_FunctionField
        if not is_FunctionField(K):
            raise ValueError("K must be a function field")
        self.__field = K
        from sage.categories.homset import Hom
        from sage.categories.sets_cat import Sets
        Map.__init__(self, Hom(K,K,Sets()))

    def _repr_type(self):
        r"""
        Return the type of this map (a derivation), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d._repr_type()
            'Derivation'

        """
        return "Derivation"

    def is_injective(self):
        r"""
        Return ``False`` since a derivation is never injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d.is_injective()
            False
        """
        return False

class FunctionFieldDerivation_rational(FunctionFieldDerivation):
    """
    Derivations on rational function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: K.derivation()
        Derivation map:
          From: Rational function field in x over Rational Field
          To:   Rational function field in x over Rational Field
    """
    def __init__(self, K, u):
        """
        Initialize a derivation of ``K`` which sends the generator of ``K`` to ``u``.

        INPUT:

        - ``K`` -- rational function field

        - ``u`` -- element of ``K``; the image of the generator of K under the
          derivation

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: TestSuite(d).run(skip=["_test_category", "_test_pickling"])

        See the comment about the test suite run in
        ``FunctionFieldDerivation.__init__``.
        """
        FunctionFieldDerivation.__init__(self, K)

        self._u = u

    def _call_(self, x):
        """
        Compute the derivation of ``x``.

        INPUT:

        - ``x`` -- element of the rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(x^3)
            3*x^2
            sage: d(1/x)
            -1/x^2
        """
        f = x.numerator()
        g = x.denominator()

        numerator = f.derivative() * g - f * g.derivative()
        if numerator.is_zero():
            return self.codomain().zero()
        else:
            return self._u * self.codomain()(numerator / g**2)

class FunctionFieldDerivation_separable(FunctionFieldDerivation):
    """
    Derivations of separable extensions.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: L.derivation()
        Derivation map:
          From: Function field in y defined by y^2 - x
          To:   Function field in y defined by y^2 - x
          Defn: y |--> 1/2/x*y
    """
    def __init__(self, L, d):
        """
        Initialize.

        INPUT:

        - ``L`` -- function field; a separable extension of the domain of ``d``

        - ``d`` -- derivation on the base function field of ``L``

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: TestSuite(d).run(skip=["_test_category", "_test_pickling"])

        See the comment about the test suite run in
        ``FunctionFieldDerivation.__init__``.
        """
        FunctionFieldDerivation.__init__(self, L)

        f = self.domain().polynomial()
        if not f.gcd(f.derivative()).is_one():
            raise ValueError("L must be a separable extension of its base field")

        x = self.domain().gen()

        self._d = d
        self._gen_image = - f.map_coefficients(lambda c: d(c))(x) / f.derivative()(x)

    def _call_(self, x):
        r"""
        Evaluate the derivation on ``x``.

        INPUT:

        - ``x`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(y)
            1/2/x*y
            sage: d(y^2)
            1
        """
        if x.is_zero():
            return self.codomain().zero()

        x = x._x
        y = self.domain().gen()

        return x.map_coefficients(self._d) + x.derivative()(y) * self._gen_image

    def _repr_defn(self):
        """
        Return the string representation of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.derivation() # indirect doctest
            Derivation map:
              From: Function field in y defined by y^2 - x
              To:   Function field in y defined by y^2 - x
              Defn: y |--> 1/2/x*y

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.derivation()
            Derivation map:
              From: Function field in z defined by z^2 - y
              To:   Function field in z defined by z^2 - y
              Defn: y |--> 1/2/x*y
                    z |--> 1/4/x*z
        """
        base = self._d._repr_defn()
        ret = '{} |--> {}'.format(self.domain().gen(), self._gen_image)
        if base:
            return base + '\n' + ret
        else:
            return ret

class FunctionFieldVectorSpaceIsomorphism(Morphism):
    """
    Base class for isomorphisms between function fields and vector spaces.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space()
        sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldVectorSpaceIsomorphism)
        True
    """
    def _repr_(self):
        """
        Return the string representation of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f
            Isomorphism:
              From: Vector space of dimension 2 over Rational function field in x over Rational Field
              To:   Function field in y defined by y^2 - x*y + 4*x^3
            sage: t
            Isomorphism:
              From: Function field in y defined by y^2 - x*y + 4*x^3
              To:   Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        s = "Isomorphism:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

    def is_injective(self):
        """
        Return ``True``, since the isomorphism is injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``True``, since the isomorphism is surjective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_surjective()
            True
        """
        return True

    def _richcmp_(self, other, op):
        r"""
        Compare this map to ``other``.

        .. NOTE::

            This implementation assumes that this isomorphism is defined by its
            domain and codomain. Isomorphisms for which this is not true must
            override this implementation.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)

            sage: K = QQbar['x'].fraction_field()
            sage: L = K.function_field()
            sage: g = K.coerce_map_from(L)

            sage: f == g
            False
            sage: f == f
            True

        """
        if type(self) != type(other):
            return NotImplemented

        from sage.structure.richcmp import richcmp
        return richcmp((self.domain(),self.codomain()), (other.domain(),other.codomain()), op)

    def __hash__(self):
        r"""
        Return a hash value of this map.

        This implementation assumes that this isomorphism is defined by its
        domain and codomain. Isomorphisms for which this is not true should
        override this implementation.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)
            sage: hash(f) == hash(f)
            True

        """
        return hash((self.domain(), self.codomain()))

class MapVectorSpaceToFunctionField(FunctionFieldVectorSpaceIsomorphism):
    """
    Isomorphism from a vector space to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); f
        Isomorphism:
          From: Vector space of dimension 2 over Rational function field in x over Rational Field
          To:   Function field in y defined by y^2 - x*y + 4*x^3
    """
    def __init__(self, V, K):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(f)
            <class 'sage.rings.function_field.maps.MapVectorSpaceToFunctionField'>
        """
        self._V = V
        self._K = K
        self._R = K.polynomial_ring()
        from sage.categories.homset import Hom
        FunctionFieldVectorSpaceIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        """
        Map ``v`` to the function field.

        INPUT:

        - ``v`` -- element of the vector space

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f(x*V.0 + (1/x^3)*V.1) # indirect doctest
            1/x^3*y + x

        TESTS:

        Test that this map is a bijection for some random inputs::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y - x)
            sage: for F in [K,L,M]:
            ....:     for base in F._intermediate_fields(K):
            ....:         V, f, t = F.vector_space(base)
            ....:         for i in range(100):
            ....:             a = F.random_element()
            ....:             assert(f(t(a)) == a)

        """
        fields = self._K._intermediate_fields(self._V.base_field())
        fields.pop()
        degrees = [k.degree() for k in fields]
        gens = [k.gen() for k in fields]

        # construct the basis composed of powers of the generators of all the
        # intermediate fields, i.e., 1, x, y, x*y, ...
        from sage.misc.misc_c import prod
        from itertools import product
        exponents = product(*[range(d) for d in degrees])
        basis = [prod(g**e for g,e in zip(gens,es)) for es in exponents]

        # multiply the entries of v with the values in basis
        coefficients = self._V(v).list()
        ret = sum([c*b for (c,b) in zip(coefficients,basis)])
        return self._K(ret)

    def domain(self):
        """
        Return the vector space which is the domain of the isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.domain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def codomain(self):
        """
        Return the function field which is the codomain of the isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.codomain()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self._K

class MapFunctionFieldToVectorSpace(FunctionFieldVectorSpaceIsomorphism):
    """
    Isomorphism from a function field to a vector space.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); t
        Isomorphism:
          From: Function field in y defined by y^2 - x*y + 4*x^3
          To:   Vector space of dimension 2 over Rational function field in x over Rational Field
    """
    def __init__(self, K, V):
        """
        Initialize.

        INPUT:

        - ``K`` -- function field

        - ``V`` -- vector space isomorphic to the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: TestSuite(t).run(skip="_test_category")
        """
        self._V = V
        self._K = K
        self._zero = K.base_ring()(0)
        self._n = K.degree()
        from sage.categories.homset import Hom
        FunctionFieldVectorSpaceIsomorphism.__init__(self, Hom(K, V))

    def _call_(self, x):
        """
        Map ``x`` to the vector space.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t(x + (1/x^3)*y) # indirect doctest
            (x, 1/x^3)

        TESTS:

        Test that this map is a bijection for some random inputs::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y - x)
            sage: for F in [K,L,M]:
            ....:     for base in F._intermediate_fields(K):
            ....:         V, f, t = F.vector_space(base)
            ....:         for i in range(100):
            ....:             a = V.random_element()
            ....:             assert(t(f(a)) == a)

        """
        ret = [x]
        fields = self._K._intermediate_fields(self._V.base_field())
        fields.pop()
        from itertools import chain
        for k in fields:
            ret = chain.from_iterable([y.list() for y in ret])
        ret = list(ret)
        assert all([t.parent() is self._V.base_field() for t in ret])
        return self._V(ret)

class FunctionFieldMorphism(RingHomomorphism):
    """
    Base class for morphisms between function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: f = K.hom(1/x); f
        Function Field endomorphism of Rational function field in x over Rational Field
          Defn: x |--> 1/x
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Rational Field
              Defn: x |--> 1/x
            sage: TestSuite(f).run(skip="_test_category")
        """
        RingHomomorphism.__init__(self, parent)

        self._im_gen = im_gen
        self._base_morphism = base_morphism

    def _repr_type(self):
        r"""
        Return the type of the morphism for the purpose of printing.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._repr_type()
            'Function Field'

        """
        return "Function Field"

    def _repr_defn(self):
        """
        Return the string containing the definition of the morphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._repr_defn()
            'y |--> 2*y'
        """
        a = '%s |--> %s'%(self.domain().variable_name(), self._im_gen)
        if self._base_morphism is not None:
            a += '\n' + self._base_morphism._repr_defn()
        return a

class FunctionFieldMorphism_polymod(FunctionFieldMorphism):
    """
    Morphism from a finite extension of a function field to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
        sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
        sage: f = L.hom(y*2); f
        Function Field endomorphism of Function field in y defined by y^3 + 6*x^3 + x
          Defn: y |--> 2*y
        sage: factor(L.polynomial())
        y^3 + 6*x^3 + x
        sage: f(y).charpoly('y')
        y^3 + 6*x^3 + x
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: TestSuite(f).run(skip="_test_category")
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)
        # Verify that the morphism is valid:
        R = self.codomain()['X']
        v = parent.domain().polynomial().list()
        if base_morphism is not None:
            v = [base_morphism(a) for a in v]
        f = R(v)
        if f(im_gen):
            raise ValueError("invalid morphism")

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f(y/x + x^2/(x+1))            # indirect doctest
            2/x*y + x^2/(x + 1)
            sage: f(y)
            2*y
        """
        v = x.list()
        if self._base_morphism is not None:
            v = [self._base_morphism(a) for a in v]
        f = v[0].parent()['X'](v)
        return f(self._im_gen)

class FunctionFieldMorphism_rational(FunctionFieldMorphism):
    """
    Morphism from a rational function field to a function field.
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> 1/x
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> 1/x
            sage: f(x+1)                          # indirect doctest
            (x + 1)/x
            sage: 1/x + 1
            (x + 1)/x

        You can specify a morphism on the base ring::

            sage: Qi = GaussianIntegers().fraction_field()
            sage: i = Qi.gen()
            sage: K.<x> = FunctionField(Qi)
            sage: phi1 = Qi.hom([CC.gen()])
            sage: phi2 = Qi.hom([-CC.gen()])
            sage: f = K.hom(CC.pi(),phi1)
            sage: f(1+i+x)
            4.14159265358979 + 1.00000000000000*I
            sage: g = K.hom(CC.pi(),phi2)
            sage: g(1+i+x)
            4.14159265358979 - 1.00000000000000*I
        """
        a = x.element()
        if self._base_morphism is None:
            return a.subs({a.parent().gen():self._im_gen})
        else:
            f = self._base_morphism
            num = a.numerator()
            den = a.denominator()
            R = self._im_gen.parent()['X']
            num = R([f(c) for c in num.list()])
            den = R([f(c) for c in den.list()])
            return num.subs(self._im_gen) / den.subs(self._im_gen)

class FunctionFieldConversionToConstantBaseField(Map):
    r"""
    Conversion map from the function field to its constant base field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: QQ.convert_map_from(K)
        Conversion map:
          From: Rational function field in x over Rational Field
          To:   Rational Field
    """
    def __init__(self, parent):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: f = QQ.convert_map_from(K)
            sage: from sage.rings.function_field.maps import FunctionFieldConversionToConstantBaseField
            sage: isinstance(f, FunctionFieldConversionToConstantBaseField)
            True
        """
        Map.__init__(self, parent)

    def _repr_type(self):
        r"""
        Return the type of this map (a conversion), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: QQ.convert_map_from(K) # indirect doctest
            Conversion map:
              From: Rational function field in x over Rational Field
              To:   Rational Field

        """
        return "Conversion"

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: QQ(K(1)) # indirect doctest
            1

        """
        return x.parent()._to_constant_base_field(x)

class FunctionFieldToFractionField(FunctionFieldVectorSpaceIsomorphism):
    r"""
    Isomorphism from rational function field to the isomorphic fraction
    field of a polynomial ring.

    EXAMPLES::

        sage: K = QQ['x'].fraction_field()
        sage: L = K.function_field()
        sage: f = K.coerce_map_from(L); f
        Isomorphism:
          From: Rational function field in x over Rational Field
          To:   Fraction Field of Univariate Polynomial Ring in x over Rational Field

    .. SEEALSO::

        :class:`FractionFieldToFunctionField`

    TESTS::

        sage: from sage.rings.function_field.maps import FunctionFieldToFractionField
        sage: isinstance(f, FunctionFieldToFractionField)
        True
        sage: TestSuite(f).run()

    """
    def _call_(self, f):
        r"""
        Return the value of this map at ``f``.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)
            sage: f(~L.gen())
            1/x

        """
        return self.codomain()(f.numerator(), f.denominator())

    def section(self):
        r"""
        Return the inverse map of this isomorphism.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)
            sage: f.section()
            Isomorphism:
                From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
                To:   Rational function field in x over Rational Field


        """
        from sage.categories.all import Hom
        parent = Hom(self.codomain(), self.domain())
        return parent.__make_element_class__(FractionFieldToFunctionField)(parent.domain(), parent.codomain())

class FractionFieldToFunctionField(FunctionFieldVectorSpaceIsomorphism):
    r"""
    Isomorphism from a fraction field of a polynomial ring to the isomorphic
    function field.

    EXAMPLES::

        sage: K = QQ['x'].fraction_field()
        sage: L = K.function_field()
        sage: f = L.coerce_map_from(K); f
        Isomorphism:
            From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
            To:   Rational function field in x over Rational Field

    .. SEEALSO::

        :class:`FunctionFieldToFractionField`

    TESTS::

        sage: from sage.rings.function_field.maps import FractionFieldToFunctionField
        sage: isinstance(f, FractionFieldToFunctionField)
        True
        sage: TestSuite(f).run()

    """
    def _call_(self, f):
        r"""
        Return the value of this morphism at ``f``.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = L.coerce_map_from(K)
            sage: f(~K.gen())
            1/x
        """
        return self.codomain()._element_constructor_(f)

    def section(self):
        r"""
        Return the inverse map of this isomorphism.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = L.coerce_map_from(K)
            sage: f.section()
            Isomorphism:
                From: Rational function field in x over Rational Field
                To:   Fraction Field of Univariate Polynomial Ring in x over Rational Field

        """
        from sage.categories.all import Hom
        parent = Hom(self.codomain(), self.domain())
        return parent.__make_element_class__(FunctionFieldToFractionField)(parent)

