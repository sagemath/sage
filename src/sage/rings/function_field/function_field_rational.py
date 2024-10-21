# sage_setup: distribution = sagemath-categories
r"""
Function Fields: rational
"""

# *****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2010      Robert Bradshaw <robertwb@math.washington.edu>
#                     2011-2018 Julian Rüth <julian.rueth@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2011      Syed Ahmad Lavasani
#                     2013-2014 Simon King
#                     2017      Dean Bisogno
#                     2017      Alyson Deines
#                     2017-2019 David Roe
#                     2017-2022 Kwankyu Lee
#                     2018      Marc Mezzarobba
#                     2018      Wilfried Luebbe
#                     2019      Brent Baccala
#                     2022      Frédéric Chapoton
#                     2022      Gonzalo Tornaría
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.arith.functions import lcm
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.structure.category_object import CategoryObject
from sage.rings.integer import Integer
from sage.categories.homset import Hom
from sage.categories.function_fields import FunctionFields

from .element import FunctionFieldElement
from .element_rational import FunctionFieldElement_rational
from .function_field import FunctionField


class RationalFunctionField(FunctionField):
    """
    Rational function field in one variable, over an arbitrary base field.

    INPUT:

    - ``constant_field`` -- arbitrary field

    - ``names`` -- string or tuple of length 1

    EXAMPLES::

        sage: K.<t> = FunctionField(GF(3)); K
        Rational function field in t over Finite Field of size 3
        sage: K.gen()
        t
        sage: 1/t + t^3 + 5
        (t^4 + 2*t + 1)/t

        sage: K.<t> = FunctionField(QQ); K
        Rational function field in t over Rational Field
        sage: K.gen()
        t
        sage: 1/t + t^3 + 5
        (t^4 + 5*t + 1)/t

    There are various ways to get at the underlying fields and rings
    associated to a rational function field::

        sage: K.<t> = FunctionField(GF(7))
        sage: K.base_field()
        Rational function field in t over Finite Field of size 7
        sage: K.field()
        Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7
        sage: K.constant_field()
        Finite Field of size 7
        sage: K.maximal_order()
        Maximal order of Rational function field in t over Finite Field of size 7

        sage: K.<t> = FunctionField(QQ)
        sage: K.base_field()
        Rational function field in t over Rational Field
        sage: K.field()
        Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: K.constant_field()
        Rational Field
        sage: K.maximal_order()
        Maximal order of Rational function field in t over Rational Field

    We define a morphism::

        sage: K.<t> = FunctionField(QQ)
        sage: L = FunctionField(QQ, 'tbar') # give variable name as second input
        sage: K.hom(L.gen())
        Function Field morphism:
          From: Rational function field in t over Rational Field
          To:   Rational function field in tbar over Rational Field
          Defn: t |--> tbar

    Here are some calculations over a number field::

        sage: R.<x> = FunctionField(QQ)
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))                                        # needs sage.rings.function_field
        sage: (y/x).divisor()                                                           # needs sage.rings.function_field
        - Place (x, y - 1)
         - Place (x, y + 1)
         + Place (x^2 + 1, y)

        sage: # needs sage.rings.number_field
        sage: A.<z> = QQ[]
        sage: NF.<i> = NumberField(z^2 + 1)
        sage: R.<x> = FunctionField(NF)
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))                                        # needs sage.rings.function_field
        sage: (x/y*x.differential()).divisor()                                          # needs sage.rings.function_field
        -2*Place (1/x, 1/x*y - 1)
         - 2*Place (1/x, 1/x*y + 1)
         + Place (x, y - 1)
         + Place (x, y + 1)
        sage: (x/y).divisor()                                                           # needs sage.rings.function_field
        - Place (x - i, y)
         + Place (x, y - 1)
         + Place (x, y + 1)
         - Place (x + i, y)
    """
    Element = FunctionFieldElement_rational

    def __init__(self, constant_field, names, category=None):
        """
        Initialize.

        EXAMPLES::

            sage: K.<t> = FunctionField(CC); K                                          # needs sage.rings.real_mpfr
            Rational function field in t over Complex Field with 53 bits of precision
            sage: TestSuite(K).run()            # long time (5s)                        # needs sage.rings.real_mpfr

            sage: FunctionField(QQ[I], 'alpha')                                         # needs sage.rings.number_field
            Rational function field in alpha over
             Number Field in I with defining polynomial x^2 + 1 with I = 1*I

        Must be over a field::

            sage: FunctionField(ZZ, 't')
            Traceback (most recent call last):
            ...
            TypeError: constant_field must be a field
        """
        if names is None:
            raise ValueError("variable name must be specified")
        elif not isinstance(names, tuple):
            names = (names, )
        if not constant_field.is_field():
            raise TypeError("constant_field must be a field")

        self._constant_field = constant_field

        FunctionField.__init__(self, self, names=names, category=FunctionFields().or_subcategory(category))

        from .place_rational import FunctionFieldPlace_rational
        self._place_class = FunctionFieldPlace_rational

        R = constant_field[names[0]]
        self._hash = hash((constant_field, names))
        self._ring = R
        self._field = R.fraction_field()

        hom = Hom(self._field, self)
        from .maps import FractionFieldToFunctionField
        self.register_coercion(hom.__make_element_class__(FractionFieldToFunctionField)(hom.domain(), hom.codomain()))

        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        from sage.categories.morphism import SetMorphism
        R.register_conversion(SetMorphism(self.Hom(R, SetsWithPartialMaps()), self._to_polynomial))

        self._gen = self(R.gen())

    def __reduce__(self):
        """
        Return the arguments which were used to create this instance. The
        rationale for this is explained in the documentation of
        :class:`UniqueRepresentation`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: clazz,args = K.__reduce__()
            sage: clazz(*args)
            Rational function field in x over Rational Field
        """
        from .constructor import FunctionField
        return FunctionField, (self._constant_field, self._names)

    def __hash__(self):
        """
        Return hash of the function field.

        The hash is formed from the constant field and the variable names.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: hash(K) == hash((K.constant_base_field(), K.variable_names()))
            True
        """
        return self._hash

    def _repr_(self):
        """
        Return string representation of the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K._repr_()
            'Rational function field in t over Rational Field'
        """
        return "Rational function field in %s over %s" % (
            self.variable_name(), self._constant_field)

    def _element_constructor_(self, x):
        r"""
        Coerce ``x`` into an element of the function field, possibly not canonically.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: a = K._element_constructor_(K.maximal_order().gen()); a
            t
            sage: a.parent()
            Rational function field in t over Rational Field

        TESTS:

        Conversion of a string::

            sage: K('t')
            t
            sage: K('1/t')
            1/t

        Conversion of a constant polynomial over the function field::

            sage: K(K.polynomial_ring().one())
            1

        Some indirect test of conversion::

            sage: S.<x, y> = K[]
            sage: I = S * [x^2 - y^2, y - t]
            sage: I.groebner_basis()                                                    # needs sage.rings.function_field
            [x^2 - t^2, y - t]
        """
        if isinstance(x, FunctionFieldElement):
            return self.element_class(self, self._field(x._x))
        try:
            x = self._field(x)
        except TypeError as Err:
            try:
                if x.parent() is self.polynomial_ring():
                    return x[0]
            except AttributeError:
                pass
            raise Err
        return self.element_class(self, x)

    def _to_constant_base_field(self, f):
        r"""
        Return ``f`` as an element of the constant base field.

        INPUT:

        - ``f`` -- element of the rational function field which is a
          constant of the underlying rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._to_constant_base_field(K(1))
            1
            sage: K._to_constant_base_field(K(x))
            Traceback (most recent call last):
            ...
            ValueError: only constants can be converted into the constant base field but x is not a constant

        TESTS:

        Verify that :issue:`21872` has been resolved::

            sage: K(1) in QQ
            True
            sage: x in QQ
            False
        """
        K = self.constant_base_field()
        if f.denominator() in K and f.numerator() in K:
            # When K is not exact, f.denominator() might not be an exact 1, so
            # we need to divide explicitly to get the correct precision
            return K(f.numerator()) / K(f.denominator())
        raise ValueError("only constants can be converted into the constant base field but %r is not a constant" % (f,))

    def _to_polynomial(self, f):
        """
        If ``f`` is integral, return it as a polynomial.

        INPUT:

        - ``f`` -- an element of this rational function field whose denominator
          is a constant

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._ring(x) # indirect doctest
            x
        """
        K = f.parent().constant_base_field()
        if f.denominator() in K:
            return f.numerator() / K(f.denominator())
        raise ValueError("only polynomials can be converted to the underlying polynomial ring")

    def _to_bivariate_polynomial(self, f):
        """
        Convert ``f`` from a univariate polynomial over the rational function
        field into a bivariate polynomial and a denominator.

        INPUT:

        - ``f`` -- univariate polynomial over the function field

        OUTPUT: bivariate polynomial, denominator

        EXAMPLES::

            sage: R.<t> = FunctionField(GF(7))
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: R._to_bivariate_polynomial(f)
            (X^7*t^2 - X^4*t^5 - X^3 + t^3, t^3)
        """
        v = f.list()
        denom = lcm([a.denominator() for a in v])
        S = denom.parent()
        x, t = S.base_ring()['%s,%s' % (f.parent().variable_name(),
                                        self.variable_name())].gens()
        phi = S.hom([t])
        return sum([phi((denom * v[i]).numerator()) * x**i for i in range(len(v))]), denom

    def _factor_univariate_polynomial(self, f, proof=None):
        """
        Factor the univariate polynomial f over the function field.

        INPUT:

        - ``f`` -- univariate polynomial over the function field

        EXAMPLES:

        We do a factorization over the function field over the rationals::

            sage: R.<t> = FunctionField(QQ)
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: f.factor()             # indirect doctest                             # needs sage.libs.singular
            (1/t) * (X - t) * (X^2 - 1/t) * (X^2 + 1/t) * (X^2 + t*X + t^2)
            sage: f.factor().prod() == f                                                # needs sage.libs.singular
            True

        We do a factorization over a finite prime field::

            sage: R.<t> = FunctionField(GF(7))
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: f.factor()                                                            # needs sage.libs.pari
            (1/t) * (X + 3*t) * (X + 5*t) * (X + 6*t) * (X^2 + 1/t) * (X^2 + 6/t)
            sage: f.factor().prod() == f                                                # needs sage.libs.pari
            True

        Factoring over a function field over a non-prime finite field::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: R.<t> = FunctionField(k)
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^3 - a*t^3)
            sage: f.factor()
            (1/t) * (X + (a + 2)*t)^3
            sage: f.factor().prod() == f
            True

        Factoring over a function field over a tower of finite fields::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(4)
            sage: R.<b> = k[]
            sage: l.<b> = k.extension(b^2 + b + a)
            sage: K.<x> = FunctionField(l)
            sage: R.<t> = K[]
            sage: F = t*x
            sage: F.factor(proof=False)
            (x) * t
        """
        old_variable_name = f.variable_name()
        # the variables of the bivariate polynomial must be distinct
        if self.variable_name() == f.variable_name():
            # replace x with xx to make the variable names distinct
            f = f.change_variable_name(old_variable_name + old_variable_name)

        F, d = self._to_bivariate_polynomial(f)
        fac = F.factor(proof=proof)
        x = f.parent().gen()
        t = f.parent().base_ring().gen()
        phi = F.parent().hom([x, t])
        v = [(phi(P), e) for P, e in fac]
        unit = phi(fac.unit()) / d
        w = []
        for a, e in v:
            c = a.leading_coefficient()
            a = a / c
            # undo any variable substitution that we introduced for the bivariate polynomial
            if old_variable_name != a.variable_name():
                a = a.change_variable_name(old_variable_name)
            unit *= (c**e)
            if a.is_unit():
                unit *= a**e
            else:
                w.append((a, e))
        from sage.structure.factorization import Factorization
        return Factorization(w, unit=unit)

    def extension(self, f, names=None):
        """
        Create an extension `L = K[y]/(f(y))` of the rational function field.

        INPUT:

        - ``f`` -- univariate polynomial over self

        - ``names`` -- string or length-1 tuple

        OUTPUT: a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)                                    # needs sage.rings.function_field
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1))                                # needs sage.rings.function_field
            Function field in y defined by y^3 + 1/t*y + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))                              # needs sage.rings.function_field
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names)

    @cached_method
    def polynomial_ring(self, var='x'):
        """
        Return a polynomial ring in one variable over the rational function field.

        INPUT:

        - ``var`` -- string; name of the variable

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.polynomial_ring()
            Univariate Polynomial Ring in x over Rational function field in x over Rational Field
            sage: K.polynomial_ring('T')
            Univariate Polynomial Ring in T over Rational function field in x over Rational Field
        """
        return self[var]

    @cached_method(key=lambda self, base, basis, map: map)
    def free_module(self, base=None, basis=None, map=True):
        """
        Return a vector space `V` and isomorphisms from the field to `V` and
        from `V` to the field.

        This function allows us to identify the elements of this field with
        elements of a one-dimensional vector space over the field itself. This
        method exists so that all function fields (rational or not) have the
        same interface.

        INPUT:

        - ``base`` -- the base field of the vector space; must be the function
          field itself (the default)

        - ``basis`` -- (ignored) a basis for the vector space

        - ``map`` -- (default: ``True``) whether to return maps to and from the
          vector space

        OUTPUT:

        - a vector space `V` over base field

        - an isomorphism from `V` to the field

        - the inverse isomorphism from the field to `V`

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.free_module()                                                                                       # needs sage.modules
            (Vector space of dimension 1 over Rational function field in x over Rational Field,
             Isomorphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field,
             Isomorphism:
              From: Rational function field in x over Rational Field
              To:   Vector space of dimension 1 over Rational function field in x over Rational Field)

        TESTS::

            sage: K.free_module()                                                                                       # needs sage.modules
            (Vector space of dimension 1 over Rational function field in x over Rational Field,
             Isomorphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field,
             Isomorphism:
              From: Rational function field in x over Rational Field
              To:   Vector space of dimension 1 over Rational function field in x over Rational Field)
        """
        if basis is not None:
            raise NotImplementedError
        from .maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        if base is None:
            base = self
        elif base is not self:
            raise ValueError("base must be the rational function field itself")
        V = base**1
        if not map:
            return V
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def random_element(self, *args, **kwds):
        """
        Create a random element of the rational function field.

        Parameters are passed to the random_element method of the
        underlying fraction field.

        EXAMPLES::

            sage: FunctionField(QQ,'alpha').random_element()   # random
            (-1/2*alpha^2 - 4)/(-12*alpha^2 + 1/2*alpha - 1/95)
        """
        return self(self._field.random_element(*args, **kwds))

    def degree(self, base=None):
        """
        Return the degree over the base field of the rational function
        field.  Since the base field is the rational function field itself, the
        degree is 1.

        INPUT:

        - ``base`` -- the base field of the vector space; must be the function
          field itself (the default)

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.degree()
            1
        """
        if base is None:
            base = self
        elif base is not self:
            raise ValueError("base must be the rational function field itself")
        from sage.rings.integer_ring import ZZ
        return ZZ(1)

    def gen(self, n=0):
        """
        Return the ``n``-th generator of the function field.  If ``n`` is not
        0, then an :class:` IndexError` is raised.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ); K.gen()
            t
            sage: K.gen().parent()
            Rational function field in t over Rational Field
            sage: K.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0:
            raise IndexError("Only one generator.")
        return self._gen

    def ngens(self):
        """
        Return the number of generators, which is 1.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.ngens()
            1
        """
        return 1

    def base_field(self):
        """
        Return the base field of the rational function field, which is just
        the function field itself.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: K.base_field()
            Rational function field in t over Finite Field of size 7
        """
        return self

    def hom(self, im_gens, base_morphism=None):
        """
        Create a homomorphism from ``self`` to another ring.

        INPUT:

        - ``im_gens`` -- exactly one element of some ring.  It must be
          invertible and transcendental over the image of
          ``base_morphism``; this is not checked.

        - ``base_morphism`` -- a homomorphism from the base field into the
          other ring; if ``None``, try to use a coercion map

        OUTPUT: a map between function fields

        EXAMPLES:

        We make a map from a rational function field to itself::

            sage: K.<x> = FunctionField(GF(7))
            sage: K.hom((x^4 + 2)/x)
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> (x^4 + 2)/x

        We construct a map from a rational function field into a
        non-rational extension field::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = K.hom(y^2 + y  + 2); f
            Function Field morphism:
              From: Rational function field in x over Finite Field of size 7
              To:   Function field in y defined by y^3 + 6*x^3 + x
              Defn: x |--> y^2 + y + 2
            sage: f(x)
            y^2 + y + 2
            sage: f(x^2)
            5*y^2 + (x^3 + 6*x + 4)*y + 2*x^3 + 5*x + 4
        """
        if isinstance(im_gens, CategoryObject):
            return self.Hom(im_gens).natural_map()
        if not isinstance(im_gens, (list, tuple)):
            im_gens = [im_gens]
        if len(im_gens) != 1:
            raise ValueError("there must be exactly one generator")
        x = im_gens[0]
        R = x.parent()
        if base_morphism is None and not R.has_coerce_map_from(self.constant_field()):
            raise ValueError("you must specify a morphism on the base field")
        from .maps import FunctionFieldMorphism_rational
        return FunctionFieldMorphism_rational(self.Hom(R), x, base_morphism)

    def field(self):
        """
        Return the underlying field, forgetting the function field
        structure.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: K.field()
            Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7

        .. SEEALSO::

            :meth:`sage.rings.fraction_field.FractionField_1poly_field.function_field`
        """
        return self._field

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        Since this is a rational function field it is of the form `K(t)`, and the
        maximal order is by definition `K[t]`, where `K` is the constant field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.maximal_order()
            Maximal order of Rational function field in t over Rational Field
            sage: K.equation_order()
            Maximal order of Rational function field in t over Rational Field
        """
        from .order_rational import FunctionFieldMaximalOrder_rational
        return FunctionFieldMaximalOrder_rational(self)

    equation_order = maximal_order

    @cached_method
    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        By definition, this is the valuation ring of the degree valuation of
        the rational function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.maximal_order_infinite()
            Maximal infinite order of Rational function field in t over Rational Field
            sage: K.equation_order_infinite()
            Maximal infinite order of Rational function field in t over Rational Field
        """
        from .order_rational import FunctionFieldMaximalOrderInfinite_rational
        return FunctionFieldMaximalOrderInfinite_rational(self)

    equation_order_infinite = maximal_order_infinite

    def constant_base_field(self):
        """
        Return the field of which the rational function field is a
        transcendental extension.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.constant_base_field()
            Rational Field
        """
        return self._constant_field

    constant_field = constant_base_field

    def different(self):
        """
        Return the different of the rational function field.

        For a rational function field, the different is simply the zero
        divisor.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.different()                                                         # needs sage.modules
            0
        """
        return self.divisor_group().zero()

    def genus(self):
        """
        Return the genus of the function field, namely 0.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.genus()
            0
        """
        return Integer(0)

    def change_variable_name(self, name):
        r"""
        Return a field isomorphic to this field with variable ``name``.

        INPUT:

        - ``name`` -- string or tuple consisting of a single string; the
          name of the new variable

        OUTPUT:

        A triple ``F,f,t`` where ``F`` is a rational function field, ``f`` is
        an isomorphism from ``F`` to this field, and ``t`` is the inverse of
        ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: L,f,t = K.change_variable_name('y')
            sage: L,f,t
            (Rational function field in y over Rational Field,
             Function Field morphism:
              From: Rational function field in y over Rational Field
              To:   Rational function field in x over Rational Field
              Defn: y |--> x,
             Function Field morphism:
              From: Rational function field in x over Rational Field
              To:   Rational function field in y over Rational Field
              Defn: x |--> y)
            sage: L.change_variable_name('x')[0] is K
            True
        """
        if isinstance(name, tuple):
            if len(name) != 1:
                raise ValueError("names must be a tuple with a single string")
            name = name[0]
        if name == self.variable_name():
            id = Hom(self, self).identity()
            return self, id, id
        else:
            from .constructor import FunctionField
            ret = FunctionField(self.constant_base_field(), name)
            return ret, ret.hom(self.gen()), self.hom(ret.gen())

    def residue_field(self, place, name=None):
        """
        Return the residue field of the place along with the maps from
        and to it.

        INPUT:

        - ``place`` -- place of the function field

        - ``name`` -- string; name of the generator of the residue field

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: p = F.places_finite(2)[0]                                             # needs sage.libs.pari
            sage: R, fr_R, to_R = F.residue_field(p)                                    # needs sage.libs.pari sage.rings.function_field
            sage: R                                                                     # needs sage.libs.pari sage.rings.function_field
            Finite Field in z2 of size 5^2
            sage: to_R(x) in R                                                          # needs sage.libs.pari sage.rings.function_field
            True
        """
        return place.residue_field(name=name)


class RationalFunctionField_char_zero(RationalFunctionField):
    """
    Rational function fields of characteristic zero.
    """
    @cached_method
    def higher_derivation(self):
        """
        Return the higher derivation for the function field.

        This is also called the Hasse-Schmidt derivation.

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: d = F.higher_derivation()                                             # needs sage.libs.singular sage.modules
            sage: [d(x^5,i) for i in range(10)]                                         # needs sage.libs.singular sage.modules
            [x^5, 5*x^4, 10*x^3, 10*x^2, 5*x, 1, 0, 0, 0, 0]
            sage: [d(x^9,i) for i in range(10)]                                         # needs sage.libs.singular sage.modules
            [x^9, 9*x^8, 36*x^7, 84*x^6, 126*x^5, 126*x^4, 84*x^3, 36*x^2, 9*x, 1]
        """
        from .derivations_polymod import FunctionFieldHigherDerivation_char_zero
        return FunctionFieldHigherDerivation_char_zero(self)


class RationalFunctionField_global(RationalFunctionField):
    """
    Rational function field over finite fields.
    """
    _differentials_space = LazyImport('sage.rings.function_field.differential', 'DifferentialsSpace_global')

    def places(self, degree=1):
        """
        Return all places of the degree.

        INPUT:

        - ``degree`` -- (default: 1) a positive integer

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.places()                                                            # needs sage.libs.pari
            [Place (1/x),
             Place (x),
             Place (x + 1),
             Place (x + 2),
             Place (x + 3),
             Place (x + 4)]
        """
        if degree == 1:
            return [self.place_infinite()] + self.places_finite(degree)
        else:
            return self.places_finite(degree)

    def places_finite(self, degree=1):
        """
        Return the finite places of the degree.

        INPUT:

        - ``degree`` -- (default: 1) a positive integer

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.places_finite()                                                     # needs sage.libs.pari
            [Place (x), Place (x + 1), Place (x + 2), Place (x + 3), Place (x + 4)]
        """
        return list(self._places_finite(degree))

    def _places_finite(self, degree=1):
        """
        Return a generator for all monic irreducible polynomials of the degree.

        INPUT:

        - ``degree`` -- (default: 1) a positive integer

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F._places_finite()
            <generator object ...>
        """
        O = self.maximal_order()
        R = O._ring
        G = R.polynomials(max_degree=degree - 1)
        lm = R.monomial(degree)
        for g in G:
            h = lm + g
            if h.is_irreducible():
                yield O.ideal(h).place()

    def place_infinite(self):
        """
        Return the unique place at infinity.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.place_infinite()
            Place (1/x)
        """
        return self.maximal_order_infinite().prime_ideal().place()

    def get_place(self, degree):
        """
        Return a place of ``degree``.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: K.get_place(1)                                                        # needs sage.libs.pari
            Place (x)
            sage: K.get_place(2)                                                        # needs sage.libs.pari
            Place (x^2 + x + 1)
            sage: K.get_place(3)                                                        # needs sage.libs.pari
            Place (x^3 + x + 1)
            sage: K.get_place(4)                                                        # needs sage.libs.pari
            Place (x^4 + x + 1)
            sage: K.get_place(5)                                                        # needs sage.libs.pari
            Place (x^5 + x^2 + 1)
        """
        for p in self._places_finite(degree):
            return p

        assert False, "there is a bug around"

    @cached_method
    def higher_derivation(self):
        """
        Return the higher derivation for the function field.

        This is also called the Hasse-Schmidt derivation.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: d = F.higher_derivation()                                             # needs sage.rings.function_field
            sage: [d(x^5,i) for i in range(10)]                                         # needs sage.rings.function_field
            [x^5, 0, 0, 0, 0, 1, 0, 0, 0, 0]
            sage: [d(x^7,i) for i in range(10)]                                         # needs sage.rings.function_field
            [x^7, 2*x^6, x^5, 0, 0, x^2, 2*x, 1, 0, 0]
        """
        from .derivations_polymod import RationalFunctionFieldHigherDerivation_global
        return RationalFunctionFieldHigherDerivation_global(self)
