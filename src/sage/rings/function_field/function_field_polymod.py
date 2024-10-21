# sage_setup: distribution = sagemath-singular
# sage.doctest: needs sage.rings.function_field
r"""
Function Fields: extension
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
from sage.rings.qqbar_decorators import handle_AA_and_QQbar
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.categories.homset import Hom
from sage.categories.function_fields import FunctionFields
from sage.categories.number_fields import NumberFields

from .element import FunctionFieldElement
from .element_polymod import FunctionFieldElement_polymod
from .function_field import FunctionField
from .function_field_rational import RationalFunctionField


class FunctionField_polymod(FunctionField):
    """
    Function fields defined by a univariate polynomial, as an extension of the
    base field.

    INPUT:

    - ``polynomial`` -- univariate polynomial over a function field

    - ``names`` -- tuple of length 1 or string; variable names

    - ``category`` -- category (default: category of function fields)

    EXAMPLES:

    We make a function field defined by a degree 5 polynomial over the
    rational function field over the rational numbers::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    We next make a function field over the above nontrivial function
    field L::

        sage: S.<z> = L[]
        sage: M.<z> = L.extension(z^2 + y*z + y); M
        Function field in z defined by z^2 + y*z + y
        sage: 1/z
        ((-x/(x^4 + 1))*y^4 + 2*x^2/(x^4 + 1))*z - 1
        sage: z * (1/z)
        1

    We drill down the tower of function fields::

        sage: M.base_field()
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
        sage: M.base_field().base_field()
        Rational function field in x over Rational Field
        sage: M.base_field().base_field().constant_field()
        Rational Field
        sage: M.constant_base_field()
        Rational Field

    .. WARNING::

        It is not checked if the polynomial used to define the function field is irreducible
        Hence it is not guaranteed that this object really is a field!
        This is illustrated below.

    ::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(x^2 - y^2)
        sage: (y - x)*(y + x)
        0
        sage: 1/(y - x)
        1
        sage: y - x == 0; y + x == 0
        False
        False
    """
    Element = FunctionFieldElement_polymod

    def __init__(self, polynomial, names, category=None):
        """
        Create a function field defined as an extension of another function
        field by adjoining a root of a univariate polynomial.

        EXAMPLES:

        We create an extension of a function field::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in y defined by y^5 + x*y - x^3 - 3*x
            sage: TestSuite(L).run(max_runs=512)   # long time (15s)

        We can set the variable name, which doesn't have to be y::

            sage: L.<w> = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in w defined by w^5 + x*w - x^3 - 3*x

        TESTS:

        Test that :issue:`17033` is fixed::

            sage: K.<t> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: M.<z> = K.extension(x^7 - x - t)
            sage: M(x)
            z
            sage: M('z')
            z
            sage: M('x')
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'x' in Fraction Field of Univariate
            Polynomial Ring in t over Rational Field
        """
        from sage.rings.polynomial.polynomial_element import Polynomial
        if polynomial.parent().ngens() > 1 or not isinstance(polynomial, Polynomial):
            raise TypeError("polynomial must be univariate a polynomial")
        if names is None:
            names = (polynomial.variable_name(), )
        elif names != polynomial.variable_name():
            polynomial = polynomial.change_variable_name(names)
        if polynomial.degree() <= 0:
            raise ValueError("polynomial must have positive degree")
        base_field = polynomial.base_ring()
        if not isinstance(base_field, FunctionField):
            raise TypeError("polynomial must be over a FunctionField")

        self._base_field = base_field
        self._polynomial = polynomial

        FunctionField.__init__(self, base_field, names=names,
                               category=FunctionFields().or_subcategory(category))

        from .place_polymod import FunctionFieldPlace_polymod
        self._place_class = FunctionFieldPlace_polymod

        self._hash = hash(polynomial)
        self._ring = self._polynomial.parent()

        self._populate_coercion_lists_(coerce_list=[base_field, self._ring])
        self._gen = self(self._ring.gen())

    def __hash__(self):
        """
        Return hash of the function field.

        The hash value is equal to the hash of the defining polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y)
            sage: hash(L) == hash(L.polynomial())
            True
        """
        return self._hash

    def _element_constructor_(self, x):
        r"""
        Make ``x`` into an element of the function field, possibly not canonically.

        INPUT:

        - ``x`` -- element

        TESTS::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L._element_constructor_(L.polynomial_ring().gen())
            y
        """
        if isinstance(x, FunctionFieldElement):
            return self.element_class(self, self._ring(x.element()))
        return self.element_class(self, self._ring(x))

    def gen(self, n=0):
        """
        Return the `n`-th generator of the function field. By default, `n` is 0; any other
        value of `n` leads to an error. The generator is the class of `y`, if we view
        the function field as being presented as `K[y]/(f(y))`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.gen()
            y
            sage: L.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        if n != 0:
            raise IndexError("there is only one generator")
        return self._gen

    def ngens(self):
        """
        Return the number of generators of the function field over its base
        field. This is by definition 1.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.ngens()
            1
        """
        return 1

    def _to_base_field(self, f):
        r"""
        Return ``f`` as an element of the :meth:`base_field`.

        INPUT:

        - ``f`` -- element of the function field which lies in the base field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L._to_base_field(L(x))
            x
            sage: L._to_base_field(y)
            Traceback (most recent call last):
            ...
            ValueError: y is not an element of the base field

        TESTS:

        Verify that :issue:`21872` has been resolved::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)

            sage: M(1) in QQ
            True
            sage: M(y) in L
            True
            sage: M(x) in K
            True
            sage: z in K
            False
        """
        K = self.base_field()
        if f.element().is_constant():
            return K(f.element())
        raise ValueError("%r is not an element of the base field" % (f,))

    def _to_constant_base_field(self, f):
        """
        Return ``f`` as an element of the :meth:`constant_base_field`.

        INPUT:

        - ``f`` -- element of the rational function field which is a
          constant

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L._to_constant_base_field(L(1))
            1
            sage: L._to_constant_base_field(y)
            Traceback (most recent call last):
            ...
            ValueError: y is not an element of the base field

        TESTS:

        Verify that :issue:`21872` has been resolved::

            sage: L(1) in QQ
            True
            sage: y in QQ
            False
        """
        return self.base_field()._to_constant_base_field(self._to_base_field(f))

    def monic_integral_model(self, names=None):
        """
        Return a function field isomorphic to this field but which is an
        extension of a rational function field with defining polynomial that is
        monic and integral over the constant base field.

        INPUT:

        - ``names`` -- string or tuple of up to two strings (default:
          ``None``); the name of the generator of the field, and the name of
          the generator of the underlying rational function field (if a tuple).
          If not given, then the names are chosen automatically.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a function field, ``f`` is an
        isomorphism from ``F`` to this field, and ``t`` is the inverse of
        ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(x^2*y^5 - 1/x); L
            Function field in y defined by x^2*y^5 - 1/x
            sage: A, from_A, to_A = L.monic_integral_model('z')
            sage: A
            Function field in z defined by z^5 - x^12
            sage: from_A
            Function Field morphism:
              From: Function field in z defined by z^5 - x^12
              To:   Function field in y defined by x^2*y^5 - 1/x
              Defn: z |--> x^3*y
                    x |--> x
            sage: to_A
            Function Field morphism:
              From: Function field in y defined by x^2*y^5 - 1/x
              To:   Function field in z defined by z^5 - x^12
              Defn: y |--> 1/x^3*z
                    x |--> x
            sage: to_A(y)
            1/x^3*z
            sage: from_A(to_A(y))
            y
            sage: from_A(to_A(1/y))
            x^3*y^4
            sage: from_A(to_A(1/y)) == 1/y
            True

        This also works for towers of function fields::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2*y - 1/x)
            sage: M.monic_integral_model()
            (Function field in z_ defined by z_^10 - x^18,
             Function Field morphism:
              From: Function field in z_ defined by z_^10 - x^18
              To:   Function field in z defined by y*z^2 - 1/x
              Defn: z_ |--> x^2*z
                    x |--> x, Function Field morphism:
              From: Function field in z defined by y*z^2 - 1/x
              To:   Function field in z_ defined by z_^10 - x^18
              Defn: z |--> 1/x^2*z_
                    y |--> 1/x^15*z_^8
                    x |--> x)

        TESTS:

        If the field is already a monic integral extension, then it is returned
        unchanged::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.monic_integral_model()
            (Function field in y defined by y^2 - x,
             Function Field endomorphism of Function field in y defined by y^2 - x
              Defn: y |--> y
                    x |--> x, Function Field endomorphism of Function field in y defined by y^2 - x
              Defn: y |--> y
                    x |--> x)

        unless ``names`` does not match with the current names::

            sage: L.monic_integral_model(names=('yy','xx'))
            (Function field in yy defined by yy^2 - xx,
             Function Field morphism:
              From: Function field in yy defined by yy^2 - xx
              To:   Function field in y defined by y^2 - x
              Defn: yy |--> y
                    xx |--> x, Function Field morphism:
              From: Function field in y defined by y^2 - x
              To:   Function field in yy defined by yy^2 - xx
              Defn: y |--> yy
                    x |--> xx)
        """
        if names:
            if not isinstance(names, tuple):
                names = (names,)
            if len(names) > 2:
                raise ValueError("names must contain at most 2 entries")

        if self.base_field() is not self.rational_function_field():
            L, from_L, to_L = self.simple_model()
            ret, ret_to_L, L_to_ret = L.monic_integral_model(names)
            from_ret = ret.hom([from_L(ret_to_L(ret.gen())),
                                from_L(ret_to_L(ret.base_field().gen()))])
            to_ret = self.hom([L_to_ret(to_L(k.gen())) for k in self._intermediate_fields(self.rational_function_field())])
            return ret, from_ret, to_ret
        else:
            if self.polynomial().is_monic() and all(c.denominator().is_one() for c in self.polynomial()):
                # self is already monic and integral
                if names is None or names == ():
                    names = (self.variable_name(),)
                return self.change_variable_name(names)
            else:
                if not names:
                    names = (self.variable_name() + "_",)
                if len(names) == 1:
                    names = (names[0], self.rational_function_field().variable_name())

                g, d = self._make_monic_integral(self.polynomial())
                K, from_K, to_K = self.base_field().change_variable_name(names[1])
                g = g.map_coefficients(to_K)
                ret = K.extension(g, names=names[0])
                from_ret = ret.hom([self.gen() * d, self.base_field().gen()])
                to_ret = self.hom([ret.gen() / d, ret.base_field().gen()])
                return ret, from_ret, to_ret

    def _make_monic_integral(self, f):
        """
        Return a monic integral polynomial `g` and an element `d` of the base
        field such that `g(y*d)=0` where `y` is a root of `f`.

        INPUT:

        - ``f`` -- polynomial

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(x^2*y^5 - 1/x)
            sage: g, d = L._make_monic_integral(L.polynomial()); g,d
            (y^5 - x^12, x^3)
            sage: (y*d).is_integral()
            True
            sage: g.is_monic()
            True
            sage: g(y*d)
            0
        """
        R = f.base_ring()
        if not isinstance(R, RationalFunctionField):
            raise NotImplementedError

        # make f monic
        n = f.degree()
        c = f.leading_coefficient()
        if c != 1:
            f = f / c

        # find lcm of denominators
        # would be good to replace this by minimal...
        d = lcm([b.denominator() for b in f.list() if b])
        if d != 1:
            x = f.parent().gen()
            g = (d**n) * f(x/d)
        else:
            g = f
        return g, d

    def constant_field(self):
        """
        Return the algebraic closure of the constant field of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^5 - x)                                          # needs sage.rings.finite_rings
            sage: L.constant_field()                                                    # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def constant_base_field(self):
        """
        Return the base constant field of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
            Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
            sage: L.constant_base_field()
            Rational Field
            sage: S.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.constant_base_field()
            Rational Field
        """
        return self.base_field().constant_base_field()

    @cached_method(key=lambda self, base: self.base_field() if base is None else base)
    def degree(self, base=None):
        """
        Return the degree of the function field over the function field ``base``.

        INPUT:

        - ``base`` -- a function field (default: ``None``), a function field
          from which this field has been constructed as a finite extension

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
            Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
            sage: L.degree()
            5
            sage: L.degree(L)
            1

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.degree(L)
            2
            sage: M.degree(K)
            10

        TESTS::

            sage: L.degree(M)
            Traceback (most recent call last):
            ...
            ValueError: base must be the rational function field itself
        """
        if base is None:
            base = self.base_field()
        if base is self:
            from sage.rings.integer_ring import ZZ
            return ZZ(1)
        return self._polynomial.degree() * self.base_field().degree(base)

    def _repr_(self):
        """
        Return the string representation of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L._repr_()
            'Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x'
        """
        return "Function field in %s defined by %s" % (self.variable_name(), self._polynomial)

    def base_field(self):
        """
        Return the base field of the function field. This function field is
        presented as `L = K[y]/(f(y))`, and the base field is by definition the
        field `K`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.base_field()
            Rational function field in x over Rational Field
        """
        return self._base_field

    def random_element(self, *args, **kwds):
        """
        Create a random element of the function field. Parameters are passed
        onto the random_element method of the base_field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x^2 + x))
            sage: L.random_element() # random
            ((x^2 - x + 2/3)/(x^2 + 1/3*x - 1))*y^2 + ((-1/4*x^2 + 1/2*x - 1)/(-5/2*x + 2/3))*y
            + (-1/2*x^2 - 4)/(-12*x^2 + 1/2*x - 1/95)
        """
        return self(self._ring.random_element(degree=self.degree(), *args, **kwds))

    def polynomial(self):
        """
        Return the univariate polynomial that defines the function field, that
        is, the polynomial `f(y)` so that the function field is of the form
        `K[y]/(f(y))`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.polynomial()
            y^5 - 2*x*y + (-x^4 - 1)/x
        """
        return self._polynomial

    def is_separable(self, base=None):
        r"""
        Return whether this is a separable extension of ``base``.

        INPUT:

        - ``base`` -- a function field from which this field has been created
          as an extension or ``None`` (default: ``None``); if ``None``, then
          return whether this is a separable extension over its base field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.is_separable()
            False
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y)
            sage: M.is_separable()
            True
            sage: M.is_separable(K)
            False

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(5))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.is_separable()
            True

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(5))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - 1)
            sage: L.is_separable()
            False
        """
        if base is None:
            base = self.base_field()
        for k in self._intermediate_fields(base)[:-1]:
            f = k.polynomial()
            g = f.derivative()
            if f.gcd(g).degree() != 0:
                return False
        return True

    def polynomial_ring(self):
        """
        Return the polynomial ring used to represent elements of the
        function field.  If we view the function field as being presented
        as `K[y]/(f(y))`, then this function returns the ring `K[y]`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.polynomial_ring()
            Univariate Polynomial Ring in y over Rational function field in x over Rational Field
        """
        return self._ring

    @cached_method(key=lambda self, base, basis, map: (self.base_field() if base is None else base, basis, map))
    def free_module(self, base=None, basis=None, map=True):
        """
        Return a vector space and isomorphisms from the field to and from the
        vector space.

        This function allows us to identify the elements of this field with
        elements of a vector space over the base field, which is useful for
        representation and arithmetic with orders, ideals, etc.

        INPUT:

        - ``base`` -- a function field (default: ``None``), the returned vector
          space is over this subfield `R`, which defaults to the base field of this
          function field.

        - ``basis`` -- a basis for this field over the base

        - ``maps`` -- boolean (default: ``True``); whether to return
          `R`-linear maps to and from `V`

        OUTPUT:

        - a vector space over the base function field

        - an isomorphism from the vector space to the field (if requested)

        - an isomorphism from the field to the vector space (if requested)

        EXAMPLES:

        We define a function field::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
            Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

        We get the vector spaces, and maps back and forth::

            sage: # needs sage.modules
            sage: V, from_V, to_V = L.free_module()
            sage: V
            Vector space of dimension 5 over Rational function field in x over Rational Field
            sage: from_V
            Isomorphism:
              From: Vector space of dimension 5 over Rational function field in x over Rational Field
              To:   Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
            sage: to_V
            Isomorphism:
              From: Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
              To:   Vector space of dimension 5 over Rational function field in x over Rational Field

        We convert an element of the vector space back to the function field::

            sage: from_V(V.1)                                                           # needs sage.modules
            y

        We define an interesting element of the function field::

            sage: a = 1/L.0; a                                                          # needs sage.modules
            (x/(x^4 + 1))*y^4 - 2*x^2/(x^4 + 1)

        We convert it to the vector space, and get a vector over the base field::

            sage: to_V(a)                                                               # needs sage.modules
            (-2*x^2/(x^4 + 1), 0, 0, 0, x/(x^4 + 1))

        We convert to and back, and get the same element::

            sage: from_V(to_V(a)) == a                                                  # needs sage.modules
            True

        In the other direction::

            sage: v = x*V.0 + (1/x)*V.1                                                 # needs sage.modules
            sage: to_V(from_V(v)) == v                                                  # needs sage.modules
            True

        And we show how it works over an extension of an extension field::

            sage: R2.<z> = L[]; M.<z> = L.extension(z^2 - y)
            sage: M.free_module()                                                       # needs sage.modules
            (Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x,
             Isomorphism:
              From: Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
              To:   Function field in z defined by z^2 - y,
             Isomorphism:
              From: Function field in z defined by z^2 - y
              To:   Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x)

        We can also get the vector space of ``M`` over ``K``::

            sage: M.free_module(K)                                                      # needs sage.modules
            (Vector space of dimension 10 over Rational function field in x over Rational Field,
             Isomorphism:
              From: Vector space of dimension 10 over Rational function field in x over Rational Field
              To:   Function field in z defined by z^2 - y,
             Isomorphism:
              From: Function field in z defined by z^2 - y
              To:   Vector space of dimension 10 over Rational function field in x over Rational Field)
        """
        if basis is not None:
            raise NotImplementedError
        from .maps import MapVectorSpaceToFunctionField, MapFunctionFieldToVectorSpace
        if base is None:
            base = self.base_field()
        degree = self.degree(base)
        V = base**degree
        if not map:
            return V
        from_V = MapVectorSpaceToFunctionField(V, self)
        to_V = MapFunctionFieldToVectorSpace(self, V)
        return (V, from_V, to_V)

    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.maximal_order()
            Maximal order of Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
        """
        from .order_polymod import FunctionFieldMaximalOrder_polymod
        return FunctionFieldMaximalOrder_polymod(self)

    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]                             # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)                        # needs sage.rings.finite_rings
            sage: F.maximal_order_infinite()                                            # needs sage.rings.finite_rings
            Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings
            sage: L.maximal_order_infinite()                                            # needs sage.rings.finite_rings
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .order_polymod import FunctionFieldMaximalOrderInfinite_polymod
        return FunctionFieldMaximalOrderInfinite_polymod(self)

    def different(self):
        """
        Return the different of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)                        # needs sage.rings.finite_rings
            sage: F.different()                                                         # needs sage.rings.finite_rings
            2*Place (x, (1/(x^3 + x^2 + x))*y^2)
             + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)
        """
        O = self.maximal_order()
        Oinf = self.maximal_order_infinite()
        return O.different().divisor() + Oinf.different().divisor()

    def equation_order(self):
        """
        Return the equation order of the function field.

        If we view the function field as being presented as `K[y]/(f(y))`, then
        the order generated by the class of `y` is returned.  If `f`
        is not monic, then :meth:`_make_monic_integral` is called, and instead
        we get the order generated by some integral multiple of a root of `f`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: O = L.equation_order()
            sage: O.basis()
            (1, x*y, x^2*y^2, x^3*y^3, x^4*y^4)

        We try an example, in which the defining polynomial is not
        monic and is not integral::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(x^2*y^5 - 1/x); L
            Function field in y defined by x^2*y^5 - 1/x
            sage: O = L.equation_order()
            sage: O.basis()
            (1, x^3*y, x^6*y^2, x^9*y^3, x^12*y^4)
        """
        d = self._make_monic_integral(self.polynomial())[1]
        return self.order(d*self.gen(), check=False)

    def hom(self, im_gens, base_morphism=None):
        """
        Create a homomorphism from the function field to another function field.

        INPUT:

        - ``im_gens`` -- list of images of the generators of the function field
          and of successive base rings

        - ``base_morphism`` -- homomorphism of the base ring, after the
          ``im_gens`` are used.  Thus if ``im_gens`` has length 2, then
          ``base_morphism`` should be a morphism from the base ring of the base
          ring of the function field.

        EXAMPLES:

        We create a rational function field, and a quadratic extension of it::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)

        We make the field automorphism that sends y to -y::

            sage: f = L.hom(-y); f
            Function Field endomorphism of Function field in y defined by y^2 - x^3 - 1
              Defn: y |--> -y

        Evaluation works::

            sage: f(y*x - 1/x)
            -x*y - 1/x

        We try to define an invalid morphism::

            sage: f = L.hom(y + 1)
            Traceback (most recent call last):
            ...
            ValueError: invalid morphism

        We make a morphism of the base rational function field::

            sage: phi = K.hom(x + 1); phi
            Function Field endomorphism of Rational function field in x over Rational Field
              Defn: x |--> x + 1
            sage: phi(x^3 - 3)
            x^3 + 3*x^2 + 3*x - 2
            sage: (x+1)^3 - 3
            x^3 + 3*x^2 + 3*x - 2

        We make a morphism by specifying where the generators and the
        base generators go::

            sage: L.hom([-y, x])
            Function Field endomorphism of Function field in y defined by y^2 - x^3 - 1
              Defn: y |--> -y
                    x |--> x

        You can also specify a morphism on the base::

            sage: R1.<q> = K[]
            sage: L1.<q> = K.extension(q^2 - (x+1)^3 - 1)
            sage: L.hom(q, base_morphism=phi)
            Function Field morphism:
              From: Function field in y defined by y^2 - x^3 - 1
              To:   Function field in q defined by q^2 - x^3 - 3*x^2 - 3*x - 2
              Defn: y |--> q
                    x |--> x + 1

        We make another extension of a rational function field::

            sage: K2.<t> = FunctionField(QQ); R2.<w> = K2[]
            sage: L2.<w> = K2.extension((4*w)^2 - (t+1)^3 - 1)

        We define a morphism, by giving the images of generators::

            sage: f = L.hom([4*w, t + 1]); f
            Function Field morphism:
              From: Function field in y defined by y^2 - x^3 - 1
              To:   Function field in w defined by 16*w^2 - t^3 - 3*t^2 - 3*t - 2
              Defn: y |--> 4*w
                    x |--> t + 1

        Evaluation works, as expected::

            sage: f(y+x)
            4*w + t + 1
            sage: f(x*y + x/(x^2+1))
            (4*t + 4)*w + (t + 1)/(t^2 + 2*t + 2)

        We make another extension of a rational function field::

            sage: K3.<yy> = FunctionField(QQ); R3.<xx> = K3[]
            sage: L3.<xx> = K3.extension(yy^2 - xx^3 - 1)

        This is the function field L with the generators exchanged. We define a morphism to L::

            sage: g = L3.hom([x,y]); g
            Function Field morphism:
              From: Function field in xx defined by -xx^3 + yy^2 - 1
              To:   Function field in y defined by y^2 - x^3 - 1
              Defn: xx |--> x
                    yy |--> y
        """
        if not isinstance(im_gens, (list, tuple)):
            im_gens = [im_gens]
        if len(im_gens) == 0:
            raise ValueError("no images specified")

        if len(im_gens) > 1:
            base_morphism = self.base_field().hom(im_gens[1:], base_morphism)

        # the codomain of this morphism is the field containing all the im_gens
        codomain = im_gens[0].parent()
        if base_morphism is not None:
            from sage.categories.pushout import pushout
            codomain = pushout(codomain, base_morphism.codomain())

        from .maps import FunctionFieldMorphism_polymod
        return FunctionFieldMorphism_polymod(self.Hom(codomain), im_gens[0], base_morphism)

    @cached_method
    def genus(self):
        """
        Return the genus of the function field.

        For now, the genus is computed using Singular.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: L.genus()
            3
        """
        # Unfortunately Singular can not compute the genus with the
        # polynomial_ring()._singular_ object because genus method only accepts
        # a ring of transcendental degree 2 over a prime field not a ring of
        # transcendental degree 1 over a rational function field of one variable

        if (isinstance(self._base_field, RationalFunctionField) and
                self._base_field.constant_field().is_prime_field()):
            from sage.interfaces.singular import singular

            # making the auxiliary ring which only has polynomials
            # with integral coefficients.
            tmpAuxRing = PolynomialRing(self._base_field.constant_field(),
                            str(self._base_field.gen()) + ',' + str(self._ring.gen()))
            intMinPoly, d = self._make_monic_integral(self._polynomial)
            curveIdeal = tmpAuxRing.ideal(intMinPoly)

            singular.lib('normal.lib')  # loading genus method in Singular
            return int(curveIdeal._singular_().genus())

        else:
            raise NotImplementedError("computation of genus over non-prime "
                                      "constant fields not implemented yet")

    def _simple_model(self, name='v'):
        r"""
        Return a finite extension `N/K(x)` isomorphic to the tower of
        extensions `M/L/K(x)` with `K` perfect.

        Helper method for :meth:`simple_model`.

        INPUT:

        - ``name`` -- string; the name of the generator of `N`

        ALGORITHM:

        Since `K` is perfect, the extension `M/K(x)` is simple, i.e., generated
        by a single element [BM1940]_. Therefore, there are only finitely many
        intermediate fields (Exercise 3.6.7 in [Bo2009]_).
        Let `a` be a generator of `M/L` and let `b` be a generator of `L/K(x)`.
        For some `i` the field `N_i=K(x)(a+x^ib)` is isomorphic to `M` and so
        it is enough to test for all terms of the form `a+x^ib` whether they
        generate a field of the right degree.
        Indeed, suppose for contradiction that for all `i` we had `N_i\neq M`.
        Then `N_i=N_j` for some `i,j`.  Thus `(a+x^ib)-(a+x^jb)=b(x^i-x^j)\in
        N_j` and so `b\in N_j`.  Similarly,
        `a+x^ib-x^{i-j}(a+x^jb)=a(1+x^{i-j})\in N_j` and so `a\in N_j`.
        Therefore, `N_j=M`.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M._simple_model()
            (Function field in v defined by v^4 - x,
             Function Field morphism:
              From: Function field in v defined by v^4 - x
              To:   Function field in z defined by z^2 - y
              Defn: v |--> z,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in v defined by v^4 - x
              Defn: z |--> v
                    y |--> v^2)

        Check that this also works for inseparable extensions::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M._simple_model()
            (Function field in v defined by v^4 + x,
             Function Field morphism:
              From: Function field in v defined by v^4 + x
              To:   Function field in z defined by z^2 + y
              Defn: v |--> z,
             Function Field morphism:
              From: Function field in z defined by z^2 + y
              To:   Function field in v defined by v^4 + x
              Defn: z |--> v
                    y |--> v^2)

        An example where the generator of the last extension does not generate
        the extension of the rational function field::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - 1)
            sage: M._simple_model()
            (Function field in v defined by v^6 + x*v^4 + x^2*v^2 + x^3 + 1,
             Function Field morphism:
               From: Function field in v defined by v^6 + x*v^4 + x^2*v^2 + x^3 + 1
               To:   Function field in z defined by z^3 + 1
               Defn: v |--> z + y,
             Function Field morphism:
               From: Function field in z defined by z^3 + 1
               To:   Function field in v defined by v^6 + x*v^4 + x^2*v^2 + x^3 + 1
               Defn: z |--> v^4 + x^2
                     y |--> v^4 + v + x^2)
        """
        M = self
        L = M.base_field()
        K = L.base_field()

        assert (isinstance(K, RationalFunctionField))
        assert (K is not L)
        assert (L is not M)

        if not K.constant_field().is_perfect():
            raise NotImplementedError("simple_model() only implemented over perfect constant fields")

        x = K.gen()
        b = L.gen()
        a = M.gen()

        # using a+x^i*b tends to lead to huge powers of x in the minimal
        # polynomial of the resulting field; it is better to try terms of
        # the form a+i*b first (but in characteristic p>0 there are only
        # finitely many of these)
        # We systematically try elements of the form a+b*factor*x^exponent
        factor = self.constant_base_field().zero()
        exponent = 0
        while True:
            v = M(a+b*factor*x**exponent)
            minpoly = v.matrix(K).minpoly()
            if minpoly.degree() == M.degree()*L.degree():
                break
            factor += 1
            if factor == 0:
                factor = self.constant_base_field().one()
                exponent += 1

        N = K.extension(minpoly, names=(name,))

        # the morphism N -> M, v |-> v
        N_to_M = N.hom(v)

        # the morphism M -> N, b |-> M_b, a |-> M_a
        V, V_to_M, M_to_V = M.free_module(K)
        V, V_to_N, N_to_V = N.free_module(K)
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(V.base_field(), V.dimension())
        # the power basis of v over K
        B = [M_to_V(v**i) for i in range(V.dimension())]
        B = MS(B)
        M_b = V_to_N(B.solve_left(M_to_V(b)))
        M_a = V_to_N(B.solve_left(M_to_V(a)))
        M_to_N = M.hom([M_a, M_b])

        return N, N_to_M, M_to_N

    @cached_method
    def simple_model(self, name=None):
        """
        Return a function field isomorphic to this field which is a simple
        extension of a rational function field.

        INPUT:

        - ``name`` -- string (default: ``None``); the name of generator of
          the simple extension. If ``None``, then the name of the generator
          will be the same as the name of the generator of this function field.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a field isomorphic to this field,
        ``f`` is an isomorphism from ``F`` to this function field and ``t`` is
        the inverse of ``f``.

        EXAMPLES:

        A tower of four function fields::

            sage: K.<x> = FunctionField(QQ); R.<z> = K[]
            sage: L.<z> = K.extension(z^2 - x); R.<u> = L[]
            sage: M.<u> = L.extension(u^2 - z); R.<v> = M[]
            sage: N.<v> = M.extension(v^2 - u)

        The fields N and M as simple extensions of K::

            sage: N.simple_model()
            (Function field in v defined by v^8 - x,
             Function Field morphism:
              From: Function field in v defined by v^8 - x
              To:   Function field in v defined by v^2 - u
              Defn: v |--> v,
             Function Field morphism:
              From: Function field in v defined by v^2 - u
              To:   Function field in v defined by v^8 - x
              Defn: v |--> v
                    u |--> v^2
                    z |--> v^4
                    x |--> x)
            sage: M.simple_model()
            (Function field in u defined by u^4 - x,
             Function Field morphism:
              From: Function field in u defined by u^4 - x
              To:   Function field in u defined by u^2 - z
              Defn: u |--> u,
             Function Field morphism:
              From: Function field in u defined by u^2 - z
              To:   Function field in u defined by u^4 - x
              Defn: u |--> u
                    z |--> u^2
                    x |--> x)

        An optional parameter ``name`` can be used to set the name of the
        generator of the simple extension::

            sage: M.simple_model(name='t')
            (Function field in t defined by t^4 - x, Function Field morphism:
              From: Function field in t defined by t^4 - x
              To:   Function field in u defined by u^2 - z
              Defn: t |--> u, Function Field morphism:
              From: Function field in u defined by u^2 - z
              To:   Function field in t defined by t^4 - x
              Defn: u |--> t
                    z |--> t^2
                    x |--> x)

        An example with higher degrees::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3)); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - x); R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - x)
            sage: M.simple_model()
            (Function field in z defined by z^15 + x*z^12 + x^2*z^9 + 2*x^3*z^6 + 2*x^4*z^3 + 2*x^5 + 2*x^3,
             Function Field morphism:
               From: Function field in z defined by z^15 + x*z^12 + x^2*z^9 + 2*x^3*z^6 + 2*x^4*z^3 + 2*x^5 + 2*x^3
               To:   Function field in z defined by z^3 + 2*x
               Defn: z |--> z + y,
             Function Field morphism:
               From: Function field in z defined by z^3 + 2*x
               To:   Function field in z defined by z^15 + x*z^12 + x^2*z^9 + 2*x^3*z^6 + 2*x^4*z^3 + 2*x^5 + 2*x^3
               Defn: z |--> 2/x*z^6 + 2*z^3 + z + 2*x
                     y |--> 1/x*z^6 + z^3 + x
                     x |--> x)

        This also works for inseparable extensions::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x); R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.simple_model()
            (Function field in z defined by z^4 + x,
             Function Field morphism:
               From: Function field in z defined by z^4 + x
               To:   Function field in z defined by z^2 + y
               Defn: z |--> z,
             Function Field morphism:
               From: Function field in z defined by z^2 + y
               To:   Function field in z defined by z^4 + x
               Defn: z |--> z
                     y |--> z^2
                     x |--> x)
        """
        if name is None:
            name = self.variable_name()

        if isinstance(self.base_field(), RationalFunctionField):
            # the extension is simple already
            if name == self.variable_name():
                id = Hom(self, self).identity()
                return self, id, id
            else:
                ret = self.base_field().extension(self.polynomial(), names=(name,))
                f = ret.hom(self.gen())
                t = self.hom(ret.gen())
                return ret, f, t
        else:
            # recursively collapse the tower of fields
            base = self.base_field()
            base_, from_base_, to_base_ = base.simple_model()
            self_ = base_.extension(self.polynomial().map_coefficients(to_base_), names=(name,))
            gens_in_base_ = [to_base_(k.gen())
                             for k in base._intermediate_fields(base.rational_function_field())]
            to_self_ = self.hom([self_.gen()] + gens_in_base_)
            from_self_ = self_.hom([self.gen(), from_base_(base_.gen())])

            # now collapse self_/base_/K(x)
            ret, ret_to_self_, self__to_ret = self_._simple_model(name)
            ret_to_self = ret.hom(from_self_(ret_to_self_(ret.gen())))
            gens_in_ret = [self__to_ret(to_self_(k.gen()))
                           for k in self._intermediate_fields(self.rational_function_field())]
            self_to_ret = self.hom(gens_in_ret)
            return ret, ret_to_self, self_to_ret

    @cached_method
    def primitive_element(self):
        r"""
        Return a primitive element over the underlying rational function field.

        If this is a finite extension of a rational function field `K(x)` with
        `K` perfect, then this is a simple extension of `K(x)`, i.e., there is
        a primitive element `y` which generates this field over `K(x)`. This
        method returns such an element `y`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: R.<z> = L[]
            sage: N.<u> = L.extension(z^2 - x - 1)
            sage: N.primitive_element()
            u + y
            sage: M.primitive_element()
            z
            sage: L.primitive_element()
            y

        This also works for inseparable extensions::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x)
            sage: R.<Z> = L[]
            sage: M.<z> = L.extension(Z^2 - y)
            sage: M.primitive_element()
            z
        """
        N, f, t = self.simple_model()
        return f(N.gen())

    @cached_method
    def separable_model(self, names=None):
        r"""
        Return a function field isomorphic to this field which is a separable
        extension of a rational function field.

        INPUT:

        - ``names`` -- tuple of two strings or ``None`` (default: ``None``);
          the second entry will be used as the variable name of the rational
          function field, the first entry will be used as the variable name of
          its separable extension. If ``None``, then the variable names will be
          chosen automatically.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a function field, ``f`` is an
        isomorphism from ``F`` to this function field, and ``t`` is the inverse
        of ``f``.

        ALGORITHM:

        Suppose that the constant base field is perfect. If this is a monic
        integral inseparable extension of a rational function field, then the
        defining polynomial is separable if we swap the variables (Proposition
        4.8 in Chapter VIII of [Lan2002]_.)
        The algorithm reduces to this case with :meth:`monic_integral_model`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3)
            sage: L.separable_model(('t','w'))
            (Function field in t defined by t^3 + w^2,
             Function Field morphism:
               From: Function field in t defined by t^3 + w^2
               To:   Function field in y defined by y^2 + x^3
               Defn: t |--> x
                     w |--> y,
             Function Field morphism:
               From: Function field in y defined by y^2 + x^3
               To:   Function field in t defined by t^3 + w^2
               Defn: y |--> w
                     x |--> t)

        This also works for non-integral polynomials::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2/x - x^2)
            sage: L.separable_model()
            (Function field in y_ defined by y_^3 + x_^2,
             Function Field morphism:
               From: Function field in y_ defined by y_^3 + x_^2
               To:   Function field in y defined by 1/x*y^2 + x^2
               Defn: y_ |--> x
                     x_ |--> y,
             Function Field morphism:
               From: Function field in y defined by 1/x*y^2 + x^2
               To:   Function field in y_ defined by y_^3 + x_^2
               Defn: y |--> x_
                     x |--> y_)

        If the base field is not perfect this is only implemented in trivial cases::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = FunctionField(GF(2))
            sage: k.is_perfect()
            False
            sage: K.<x> = FunctionField(k)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - t)
            sage: L.separable_model()
            (Function field in y defined by y^3 + t,
             Function Field endomorphism of Function field in y defined by y^3 + t
               Defn: y |--> y
                     x |--> x,
             Function Field endomorphism of Function field in y defined by y^3 + t
               Defn: y |--> y
                     x |--> x)

        Some other cases for which a separable model could be constructed are
        not supported yet::

            sage: R.<y> = K[]                                                           # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(y^2 - t)                                          # needs sage.rings.finite_rings
            sage: L.separable_model()                                                   # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NotImplementedError: constructing a separable model is only implemented
            for function fields over a perfect constant base field

        TESTS:

        Check that this also works in characteristic zero::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3)
            sage: L.separable_model()
            (Function field in y defined by y^2 - x^3,
             Function Field endomorphism of Function field in y defined by y^2 - x^3
               Defn: y |--> y
                     x |--> x,
             Function Field endomorphism of Function field in y defined by y^2 - x^3
               Defn: y |--> y
                     x |--> x)

        Check that this works for towers of inseparable extensions::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.separable_model()
            (Function field in z_ defined by z_ + x_^4,
             Function Field morphism:
               From: Function field in z_ defined by z_ + x_^4
               To:   Function field in z defined by z^2 + y
               Defn: z_ |--> x
                     x_ |--> z,
             Function Field morphism:
               From: Function field in z defined by z^2 + y
               To:   Function field in z_ defined by z_ + x_^4
               Defn: z |--> x_
                     y |--> x_^2
                     x |--> x_^4)

        Check that this also works if only the first extension is inseparable::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y)
            sage: M.separable_model()
            (Function field in z_ defined by z_ + x_^6,
             Function Field morphism:
               From: Function field in z_ defined by z_ + x_^6
               To:   Function field in z defined by z^3 + y
               Defn: z_ |--> x
                     x_ |--> z,
             Function Field morphism:
               From: Function field in z defined by z^3 + y
               To:   Function field in z_ defined by z_ + x_^6
               Defn: z |--> x_
                     y |--> x_^3
                     x |--> x_^6)
        """
        if names is None:
            pass
        elif not isinstance(names, tuple):
            raise TypeError("names must be a tuple consisting of two strings")
        elif len(names) != 2:
            raise ValueError("must provide exactly two variable names")

        if self.base_ring() is not self.rational_function_field():
            L, from_L, to_L = self.simple_model()
            K, from_K, to_K = L.separable_model(names=names)
            f = K.hom([from_L(from_K(K.gen())), from_L(from_K(K.base_field().gen()))])
            t = self.hom([to_K(to_L(k.gen())) for k in self._intermediate_fields(self.rational_function_field())])
            return K, f, t

        if self.polynomial().gcd(self.polynomial().derivative()).is_one():
            # the model is already separable
            if names is None:
                names = self.variable_name(), self.base_field().variable_name()
            return self.change_variable_name(names)

        if not self.constant_base_field().is_perfect():
            raise NotImplementedError("constructing a separable model is only implemented for function fields over a perfect constant base field")

        if names is None:
            names = (self.variable_name()+"_", self.rational_function_field().variable_name()+"_")

        L, from_L, to_L = self.monic_integral_model()

        if L.polynomial().gcd(L.polynomial().derivative()).is_one():
            # L is separable
            ret, ret_to_L, L_to_ret = L.change_variable_name(names)
            f = ret.hom([from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))])
            t = self.hom([L_to_ret(to_L(self.gen())), L_to_ret(to_L(self.base_field().gen()))])
            return ret, f, t
        else:
            # otherwise, the polynomial of L must be separable in the other variable
            from .constructor import FunctionField
            K = FunctionField(self.constant_base_field(), names=(names[1],))
            # construct a field isomorphic to L on top of K

            # turn the minpoly of K into a bivariate polynomial
            if names[0] == names[1]:
                raise ValueError("names of generators must be distinct")
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.constant_base_field(), names=names)
            S = R.remove_var(names[1])
            f = R(L.polynomial().change_variable_name(names[1]).map_coefficients(
                lambda c: c.numerator().change_variable_name(names[0]), S))
            f = f.polynomial(R.gen(0)).change_ring(K)
            f /= f.leading_coefficient()
            # f must be separable in the other variable (otherwise it would factor)
            assert f.gcd(f.derivative()).is_one()

            ret = K.extension(f, names=(names[0],))
            # isomorphisms between L and ret are given by swapping generators
            ret_to_L = ret.hom([L(L.base_field().gen()), L.gen()])
            L_to_ret = L.hom([ret(K.gen()), ret.gen()])
            # compose with from_L and to_L to get the desired isomorphisms between self and ret
            f = ret.hom([from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))])
            t = self.hom([L_to_ret(to_L(self.gen())), L_to_ret(to_L(self.base_field().gen()))])
            return ret, f, t

    def change_variable_name(self, name):
        r"""
        Return a field isomorphic to this field with variable(s) ``name``.

        INPUT:

        - ``name`` -- string or tuple consisting of a strings; the names of
          the new variables starting with a generator of this field and going
          down to the rational function field

        OUTPUT:

        A triple ``F,f,t`` where ``F`` is a function field, ``f`` is an
        isomorphism from ``F`` to this field, and ``t`` is the inverse of
        ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)

            sage: M.change_variable_name('zz')
            (Function field in zz defined by zz^2 - y,
             Function Field morphism:
              From: Function field in zz defined by zz^2 - y
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    y |--> y
                    x |--> x,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in zz defined by zz^2 - y
              Defn: z |--> zz
                    y |--> y
                    x |--> x)
            sage: M.change_variable_name(('zz','yy'))
            (Function field in zz defined by zz^2 - yy,
             Function Field morphism:
              From: Function field in zz defined by zz^2 - yy
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    yy |--> y
                    x |--> x,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in zz defined by zz^2 - yy
              Defn: z |--> zz
                    y |--> yy
                    x |--> x)
            sage: M.change_variable_name(('zz','yy','xx'))
            (Function field in zz defined by zz^2 - yy,
             Function Field morphism:
              From: Function field in zz defined by zz^2 - yy
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    yy |--> y
                    xx |--> x,
             Function Field morphism:
              From: Function field in z defined by z^2 - y
              To:   Function field in zz defined by zz^2 - yy
              Defn: z |--> zz
                    y |--> yy
                    x |--> xx)
        """
        if not isinstance(name, tuple):
            name = (name,)
        if len(name) == 0:
            raise ValueError("name must contain at least one string")
        elif len(name) == 1:
            base = self.base_field()
            from_base = to_base = Hom(base, base).identity()
        else:
            base, from_base, to_base = self.base_field().change_variable_name(name[1:])

        ret = base.extension(self.polynomial().map_coefficients(to_base), names=(name[0],))
        f = ret.hom([k.gen() for k in self._intermediate_fields(self.rational_function_field())])
        t = self.hom([k.gen() for k in ret._intermediate_fields(ret.rational_function_field())])
        return ret, f, t


class FunctionField_simple(FunctionField_polymod):
    """
    Function fields defined by irreducible and separable polynomials
    over rational function fields.
    """
    @cached_method
    def _inversion_isomorphism(self):
        r"""
        Return an inverted function field isomorphic to ``self`` and isomorphisms
        between them.

        An *inverted* function field `M` is an extension of the base rational
        function field `k(x)` of ``self``, and isomorphic to ``self`` by an
        isomorphism sending `x` to `1/x`, which we call an *inversion*
        isomorphism.  Also the defining polynomial of the function field `M` is
        required to be monic and integral.

        The inversion isomorphism is for internal use to reposition infinite
        places to finite places.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]                             # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)                        # needs sage.rings.finite_rings
            sage: F._inversion_isomorphism()                                            # needs sage.rings.finite_rings
            (Function field in s defined by s^3 + x^16 + x^14 + x^12, Composite map:
               From: Function field in s defined by s^3 + x^16 + x^14 + x^12
               To:   Function field in y defined by y^3 + x^6 + x^4 + x^2
               Defn:   Function Field morphism:
                       From: Function field in s defined by s^3 + x^16 + x^14 + x^12
                       To:   Function field in T defined by T^3 + (x^4 + x^2 + 1)/x^6
                       Defn: s |--> x^6*T
                             x |--> x
                     then
                       Function Field morphism:
                       From: Function field in T defined by T^3 + (x^4 + x^2 + 1)/x^6
                       To:   Function field in y defined by y^3 + x^6 + x^4 + x^2
                       Defn: T |--> y
                             x |--> 1/x,
             Composite map:
               From: Function field in y defined by y^3 + x^6 + x^4 + x^2
               To:   Function field in s defined by s^3 + x^16 + x^14 + x^12
               Defn:   Function Field morphism:
                       From: Function field in y defined by y^3 + x^6 + x^4 + x^2
                       To:   Function field in T defined by T^3 + (x^4 + x^2 + 1)/x^6
                       Defn: y |--> T
                             x |--> 1/x
                     then
                       Function Field morphism:
                       From: Function field in T defined by T^3 + (x^4 + x^2 + 1)/x^6
                       To:   Function field in s defined by s^3 + x^16 + x^14 + x^12
                       Defn: T |--> 1/x^6*s
                             x |--> x)
        """
        K = self.base_field()
        R = PolynomialRing(K, 'T')
        x = K.gen()
        xinv = 1/x

        h = K.hom(xinv)
        F_poly = R([h(c) for c in self.polynomial().list()])
        F = K.extension(F_poly)

        self2F = self.hom([F.gen(), xinv])
        F2self = F.hom([self.gen(), xinv])

        M, M2F, F2M = F.monic_integral_model('s')

        return M, F2self*M2F, F2M*self2F

    def places_above(self, p):
        """
        Return places lying above ``p``.

        INPUT:

        - ``p`` -- place of the base rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)                        # needs sage.rings.finite_rings
            sage: all(q.place_below() == p                                              # needs sage.rings.finite_rings
            ....:     for p in K.places() for q in F.places_above(p))
            True

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = K.maximal_order()
            sage: pls = [O.ideal(x - c).place() for c in [-2, -1, 0, 1, 2]]
            sage: all(q.place_below() == p
            ....:     for p in pls for q in F.places_above(p))
            True

            sage: # needs sage.rings.number_field
            sage: K.<x> = FunctionField(QQbar); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = K.maximal_order()
            sage: pls = [O.ideal(x - QQbar(sqrt(c))).place()
            ....:        for c in [-2, -1, 0, 1, 2]]
            sage: all(q.place_below() == p      # long time (4s)
            ....:     for p in pls for q in F.places_above(p))
            True
        """
        R = self.base_field()

        if p not in R.place_set():
            raise TypeError("not a place of the base rational function field")

        if p.is_infinite_place():
            dec = self.maximal_order_infinite().decomposition()
        else:
            dec = self.maximal_order().decomposition(p.prime_ideal())

        return tuple([q.place() for q, deg, exp in dec])

    def constant_field(self):
        """
        Return the algebraic closure of the base constant field in the function
        field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3)); _.<y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))                        # needs sage.rings.finite_rings
            sage: L.constant_field()                                                    # needs sage.rings.finite_rings
            Finite Field of size 3
        """
        return self.exact_constant_field()[0]

    def exact_constant_field(self, name='t'):
        """
        Return the exact constant field and its embedding into the function field.

        INPUT:

        - ``name`` -- name (default: `t`) of the generator of the exact constant field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3)); _.<Y> = K[]
            sage: f = Y^2 - x*Y + x^2 + 1 # irreducible but not absolutely irreducible
            sage: L.<y> = K.extension(f)
            sage: L.genus()
            0
            sage: L.exact_constant_field()
            (Finite Field in t of size 3^2, Ring morphism:
               From: Finite Field in t of size 3^2
               To:   Function field in y defined by y^2 + 2*x*y + x^2 + 1
               Defn: t |--> y + x)
            sage: (y+x).divisor()
            0
        """
        # A basis of the full constant field is obtained from
        # computing a Riemann-Roch basis of zero divisor.
        basis = self.divisor_group().zero().basis_function_space()

        dim = len(basis)

        for e in basis:
            _min_poly = e.minimal_polynomial(name)
            if _min_poly.degree() == dim:
                break
        k = self.constant_base_field()
        R = k[name]
        min_poly = R([k(c) for c in _min_poly.list()])

        k_ext = k.extension(min_poly, name)

        if k_ext.is_prime_field():
            # The cover of the quotient ring k_ext is the integer ring
            # whose generator is 1. This is different from the generator
            # of k_ext.
            embedding = k_ext.hom([self(1)], self)
        else:
            embedding = k_ext.hom([e], self)

        return k_ext, embedding

    def genus(self):
        """
        Return the genus of the function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(16)
            sage: K.<x> = FunctionField(F); K
            Rational function field in x over Finite Field in a of size 2^4
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L.genus()
            6

            sage: # needs sage.rings.number_field
            sage: R.<T> = QQ[]
            sage: N.<a> = NumberField(T^2 + 1)
            sage: K.<x> = FunctionField(N); K
            Rational function field in x over Number Field in a with defining polynomial T^2 + 1
            sage: K.genus()
            0
            sage: S.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^2 - x^3 + x)
            sage: L.genus()
            1

        The genus is computed by the Hurwitz genus formula.
        """
        k, _ = self.exact_constant_field()
        if k in NumberFields():
            k_degree = k.relative_degree()
        else:
            k_degree = k.degree()
        different_degree = self.different().degree()  # must be even
        return Integer(different_degree // 2 - self.degree() / k_degree) + 1

    def residue_field(self, place, name=None):
        """
        Return the residue field associated with the place along with the maps
        from and to the residue field.

        INPUT:

        - ``place`` -- place of the function field

        - ``name`` -- string; name of the generator of the residue field

        The domain of the map to the residue field is the discrete valuation
        ring associated with the place.

        The discrete valuation ring is defined as the ring of all elements of
        the function field with nonnegative valuation at the place. The maximal
        ideal is the set of elements of positive valuation.  The residue field
        is then the quotient of the discrete valuation ring by its maximal
        ideal.

        If an element not in the valuation ring is applied to the map, an
        exception :exc:`TypeError` is raised.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R, fr_R, to_R = L.residue_field(p)
            sage: R
            Finite Field of size 2
            sage: f = 1 + y
            sage: f.valuation(p)
            -1
            sage: to_R(f)
            Traceback (most recent call last):
            ...
            TypeError: ...
            sage: (1+1/f).valuation(p)
            0
            sage: to_R(1 + 1/f)
            1
            sage: [fr_R(e) for e in R]
            [0, 1]
        """
        return place.residue_field(name=name)

    def places_infinite(self, degree=1):
        """
        Return a list of the infinite places with ``degree``.

        INPUT:

        - ``degree`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L.places_infinite(1)
            [Place (1/x, 1/x^4*y^3)]
        """
        return list(self._places_infinite(degree))

    def _places_infinite(self, degree):
        """
        Return a generator of *infinite* places with ``degree``.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L._places_infinite(1)
            <generator object ...>
        """
        Oinf = self.maximal_order_infinite()
        for prime, _, _ in Oinf.decomposition():
            place = prime.place()
            if place.degree() == degree:
                yield place


class FunctionField_char_zero(FunctionField_simple):
    """
    Function fields of characteristic zero.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
        sage: L
        Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
        sage: L.characteristic()
        0
    """
    @cached_method
    def higher_derivation(self):
        """
        Return the higher derivation (also called the Hasse-Schmidt derivation)
        for the function field.

        The higher derivation of the function field is uniquely determined with
        respect to the separating element `x` of the base rational function
        field `k(x)`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.higher_derivation()                                                 # needs sage.modules
            Higher derivation map:
              From: Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
              To:   Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
        """
        from .derivations_polymod import FunctionFieldHigherDerivation_char_zero
        return FunctionFieldHigherDerivation_char_zero(self)


class FunctionField_global(FunctionField_simple):
    """
    Global function fields.

    INPUT:

    - ``polynomial`` -- monic irreducible and separable polynomial

    - ``names`` -- name of the generator of the function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                            # needs sage.rings.finite_rings
        sage: L                                                                         # needs sage.rings.finite_rings
        Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)

    The defining equation needs not be monic::

        sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension((1 - x)*Y^7 - x^3)                                    # needs sage.rings.finite_rings
        sage: L.gaps()                          # long time (6s)                        # needs sage.rings.finite_rings
        [1, 2, 3]

    or may define a trivial extension::

        sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(Y-1)                                                  # needs sage.rings.finite_rings
        sage: L.genus()                                                                 # needs sage.rings.finite_rings
        0
    """
    _differentials_space = LazyImport('sage.rings.function_field.differential', 'DifferentialsSpace_global')

    def __init__(self, polynomial, names):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                        # needs sage.rings.finite_rings
            sage: TestSuite(L).run()            # long time (7s)                        # needs sage.rings.finite_rings
        """
        FunctionField_polymod.__init__(self, polynomial, names)

    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18)
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y)
        """
        from .order_polymod import FunctionFieldMaximalOrder_global
        return FunctionFieldMaximalOrder_global(self)

    @cached_method
    def higher_derivation(self):
        """
        Return the higher derivation (also called the Hasse-Schmidt derivation)
        for the function field.

        The higher derivation of the function field is uniquely determined with
        respect to the separating element `x` of the base rational function
        field `k(x)`.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                        # needs sage.rings.finite_rings
            sage: L.higher_derivation()                                                 # needs sage.modules sage.rings.finite_rings
            Higher derivation map:
              From: Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
              To:   Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        from .derivations_polymod import FunctionFieldHigherDerivation_global
        return FunctionFieldHigherDerivation_global(self)

    def get_place(self, degree):
        """
        Return a place of ``degree``.

        INPUT:

        - ``degree`` -- positive integer

        OUTPUT: a place of ``degree`` if any exists; otherwise ``None``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<Y> = PolynomialRing(K)
            sage: L.<y> = K.extension(Y^4 + Y - x^5)
            sage: L.get_place(1)
            Place (x, y)
            sage: L.get_place(2)
            Place (x, y^2 + y + 1)
            sage: L.get_place(3)
            Place (x^3 + x^2 + 1, y + x^2 + x)
            sage: L.get_place(4)
            Place (x + 1, x^5 + 1)
            sage: L.get_place(5)
            Place (x^5 + x^3 + x^2 + x + 1, y + x^4 + 1)
            sage: L.get_place(6)
            Place (x^3 + x^2 + 1, y^2 + y + x^2)
            sage: L.get_place(7)
            Place (x^7 + x + 1, y + x^6 + x^5 + x^4 + x^3 + x)
            sage: L.get_place(8)
        """
        for p in self._places_finite(degree):
            return p

        for p in self._places_infinite(degree):
            return p

        return None

    def places(self, degree=1):
        """
        Return a list of the places with ``degree``.

        INPUT:

        - ``degree`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L.places(1)
            [Place (1/x, 1/x^4*y^3), Place (x, y), Place (x, y + 1)]
        """
        return self.places_infinite(degree) + self.places_finite(degree)

    def places_finite(self, degree=1):
        """
        Return a list of the finite places with ``degree``.

        INPUT:

        - ``degree`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L.places_finite(1)
            [Place (x, y), Place (x, y + 1)]
        """
        return list(self._places_finite(degree))

    def _places_finite(self, degree):
        """
        Return a generator of finite places with ``degree``.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4 + t - x^5)
            sage: L._places_finite(1)
            <generator object ...>
        """
        O = self.maximal_order()
        K = self.base_field()

        degree = Integer(degree)

        for d in degree.divisors():
            for p in K._places_finite(degree=d):
                for prime, _, _ in O.decomposition(p.prime_ideal()):
                    place = prime.place()
                    if place.degree() == degree:
                        yield place

    def gaps(self):
        """
        Return the gaps of the function field.

        These are the gaps at the ordinary places, that is, places which are
        not Weierstrass places.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)                                # needs sage.rings.finite_rings
            sage: L.gaps()                                                              # needs sage.modules sage.rings.finite_rings
            [1, 2, 3]
        """
        return self._weierstrass_places()[1]

    def weierstrass_places(self):
        """
        Return all Weierstrass places of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)                                # needs sage.rings.finite_rings
            sage: L.weierstrass_places()                                                # needs sage.modules sage.rings.finite_rings
            [Place (1/x, 1/x^3*y^2 + 1/x),
             Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
             Place (x, y),
             Place (x + 1, (x^3 + 1)*y + x + 1),
             Place (x^3 + x + 1, y + 1),
             Place (x^3 + x + 1, y + x^2),
             Place (x^3 + x + 1, y + x^2 + 1),
             Place (x^3 + x^2 + 1, y + x),
             Place (x^3 + x^2 + 1, y + x^2 + 1),
             Place (x^3 + x^2 + 1, y + x^2 + x + 1)]
        """
        return self._weierstrass_places()[0].support()

    @cached_method
    def _weierstrass_places(self):
        """
        Return the Weierstrass places together with the gap sequence for
        ordinary places.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)                                # needs sage.rings.finite_rings
            sage: len(L.weierstrass_places())  # indirect doctest                       # needs sage.modules sage.rings.finite_rings
            10

        This method implements Algorithm 30 in [Hes2002b]_.
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        W = self(self.base_field().gen()).differential().divisor()
        basis = W._basis()

        if not basis:
            return [], []
        d = len(basis)

        der = self.higher_derivation()
        M = matrix([basis])
        e = 1
        gaps = [1]
        while M.nrows() < d:
            row = vector([der._derive(basis[i], e) for i in range(d)])
            if row not in M.row_space():
                M = matrix(M.rows() + [row])
                gaps.append(e + 1)
            e += 1

        # This is faster than M.determinant(). Note that Mx
        # is a matrix over univariate polynomial ring.
        Mx = matrix(M.nrows(), [c._x for c in M.list()])
        detM = self(Mx.determinant() % self._polynomial)

        R = detM.divisor() + sum(gaps)*W  # ramification divisor

        return R, gaps

    @cached_method
    def L_polynomial(self, name='t'):
        """
        Return the L-polynomial of the function field.

        INPUT:

        - ``name`` -- (default: ``t``) name of the variable of the polynomial

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings
            sage: F.L_polynomial()                                                      # needs sage.rings.finite_rings
            2*t^2 + t + 1
        """
        from sage.rings.integer_ring import ZZ
        q = self.constant_field().order()
        g = self.genus()

        B = [len(self.places(i+1)) for i in range(g)]
        N = [sum(d * B[d-1] for d in ZZ(i+1).divisors()) for i in range(g)]
        S = [N[i] - q**(i+1) - 1 for i in range(g)]

        a = [1]
        for i in range(1, g+1):
            a.append(sum(S[j] * a[i-j-1] for j in range(i)) / i)
        for j in range(1, g+1):
            a.append(q**j * a[g-j])

        return ZZ[name](a)

    def number_of_rational_places(self, r=1):
        """
        Return the number of rational places of the function field whose
        constant field extended by degree ``r``.

        INPUT:

        - ``r`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: F.number_of_rational_places()
            4
            sage: [F.number_of_rational_places(r) for r in [1..10]]
            [4, 8, 4, 16, 44, 56, 116, 288, 508, 968]
        """
        from sage.rings.integer_ring import IntegerRing

        q = self.constant_field().order()
        L = self.L_polynomial()
        Lp = L.derivative()

        R = IntegerRing()[[L.parent().gen()]]  # power series ring

        f = R(Lp / L, prec=r)
        n = f[r-1] + q**r + 1

        return n


@handle_AA_and_QQbar
def _singular_normal(ideal):
    r"""
    Compute the normalization of the affine algebra defined by ``ideal`` using
    Singular.

    The affine algebra is the quotient algebra of a multivariate polynomial
    ring `R` by the ideal. The normalization is by definition the integral
    closure of the algebra in its total ring of fractions.

    INPUT:

    - ``ideal`` -- a radical ideal in a multivariate polynomial ring

    OUTPUT:

    a list of lists, one list for each ideal in the equidimensional
    decomposition of the ``ideal``, each list giving a set of generators of the
    normalization of each ideal as an R-module by dividing all elements of the
    list by the final element. Thus the list ``[x, y]`` means that `\{x/y, 1\}`
    is the set of generators of the normalization of `R/(x,y)`.

    ALGORITHM:

    Singular's implementation of the normalization algorithm described in G.-M.
    Greuel, S. Laplagne, F. Seelisch: Normalization of Rings (2009).

    EXAMPLES::

        sage: from sage.rings.function_field.function_field_polymod import _singular_normal
        sage: R.<x,y> = QQ[]

        sage: f = (x^2 - y^3) * x
        sage: _singular_normal(ideal(f))
        [[x, y], [1]]

        sage: f = y^2 - x
        sage: _singular_normal(ideal(f))
        [[1]]
    """
    from sage.libs.singular.function import (singular_function,
                                             lib as singular_lib,
                                             get_printlevel, set_printlevel)
    singular_lib('normal.lib')
    normal = singular_function('normal')

    # verbose unless printlevel is -1.
    saved_printlevel = get_printlevel()
    set_printlevel(-1)
    nor = normal(ideal)
    set_printlevel(saved_printlevel)

    return nor[1]


class FunctionField_integral(FunctionField_simple):
    """
    Integral function fields.

    A function field is integral if it is defined by an irreducible separable
    polynomial, which is integral over the maximal order of the base rational
    function field.
    """
    def _maximal_order_basis(self):
        """
        Return a basis of the maximal order of the function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18)
            sage: F._maximal_order_basis()
            [1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y]

        The basis of the maximal order *always* starts with 1. This is assumed
        in some algorithms.
        """
        from sage.matrix.constructor import matrix
        from .hermite_form_polynomial import reversed_hermite_form

        k = self.constant_base_field()
        K = self.base_field()  # rational function field
        n = self.degree()

        # Construct the defining polynomial of the function field as a
        # two-variate polynomial g in the ring k[y,x] where k is the constant
        # base field.
        S, (y, x) = PolynomialRing(k, names='y,x', order='lex').objgens()
        v = self.polynomial().list()
        g = sum([v[i].numerator().subs(x) * y**i for i in range(len(v))])

        if self.is_global():
            from sage.libs.singular.function import singular_function, lib
            from sage.env import SAGE_EXTCODE
            lib(SAGE_EXTCODE + '/singular/function_field/core.lib')
            normalize = singular_function('core_normalize')

            # Singular "normalP" algorithm assumes affine domain over
            # a prime field. So we construct gflat lifting g as in
            # k_prime[yy,xx,zz]/(k_poly) where k = k_prime[zz]/(k_poly)
            R = PolynomialRing(k.prime_subfield(), names='yy,xx,zz')
            gflat = R.zero()
            for m in g.monomials():
                c = g.monomial_coefficient(m).polynomial('zz')
                gflat += R(c) * R(m)  # R(m) is a monomial in yy and xx

            k_poly = R(k.polynomial('zz'))

            # invoke Singular
            pols_in_R = normalize(R.ideal([k_poly, gflat]))

            # reconstruct polynomials in S
            h = R.hom([y, x, k.gen()], S)
            pols_in_S = [h(f) for f in pols_in_R]
        else:
            # Call Singular. Singular's "normal" function returns a basis
            # of the integral closure of k(x,y)/(g) as a k[x,y]-module.
            pols_in_S = _singular_normal(S.ideal(g))[0]

        # reconstruct the polynomials in the function field
        x = K.gen()
        y = self.gen()
        pols = []
        for f in pols_in_S:
            p = f.polynomial(S.gen(0))
            s = 0
            for i in range(p.degree()+1):
                s += p[i].subs(x) * y**i
            pols.append(s)

        # Now if pols = [g1,g2,...gn,g0], then the g1/g0,g2/g0,...,gn/g0,
        # and g0/g0=1 are the module generators of the integral closure
        # of the equation order Sb = k[xb,yb] in its fraction field,
        # that is, the function field. The integral closure of k[x]
        # is then obtained by multiplying these generators with powers of y
        # as the equation order itself is an integral extension of k[x].
        d = ~ pols[-1]
        _basis = []
        for f in pols:
            b = d * f
            for i in range(n):
                _basis.append(b)
                b *= y

        # Finally we reduce _basis to get a basis over k[x]. This is done of
        # course by Hermite normal form computation. Here we apply a trick to
        # get a basis that starts with 1 and is ordered in increasing
        # y-degrees. The trick is to use the reversed Hermite normal form.
        # Note that it is important that the overall denominator l lies in k[x].
        V, fr_V, to_V = self.free_module()
        basis_V = [to_V(bvec) for bvec in _basis]
        l = lcm([vvec.denominator() for vvec in basis_V])

        _mat = matrix([[coeff.numerator() for coeff in l*v] for v in basis_V])
        reversed_hermite_form(_mat)

        basis = [fr_V(v) / l for v in _mat if not v.is_zero()]
        return basis

    @cached_method
    def equation_order(self):
        """
        Return the equation order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)               # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)                            # needs sage.rings.finite_rings
            sage: F.equation_order()                                                    # needs sage.rings.finite_rings
            Order in Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(QQ); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: F.equation_order()
            Order in Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2
        """
        from .order_basis import FunctionFieldOrder_basis
        a = self.gen()
        basis = [a**i for i in range(self.degree())]
        return FunctionFieldOrder_basis(tuple(basis))

    @cached_method
    def primitive_integal_element_infinite(self):
        """
        Return a primitive integral element over the base maximal infinite order.

        This element is integral over the maximal infinite order of the base
        rational function field and the function field is a simple extension by
        this element over the base order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: b = F.primitive_integal_element_infinite(); b
            1/x^2*y
            sage: b.minimal_polynomial('t')
            t^3 + (x^4 + x^2 + 1)/x^4
        """
        f = self.polynomial()
        n = f.degree()
        y = self.gen()
        x = self.base_field().gen()

        cf = max([(f[i].numerator().degree()/(n-i)).ceil() for i in range(n)
                  if f[i] != 0])
        return y*x**(-cf)

    @cached_method
    def equation_order_infinite(self):
        """
        Return the infinite equation order of the function field.

        This is by definition `o[b]` where `b` is the primitive integral
        element from :meth:`primitive_integal_element_infinite()` and `o` is
        the maximal infinite order of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)               # needs sage.rings.finite_rings
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)                            # needs sage.rings.finite_rings
            sage: F.equation_order_infinite()                                           # needs sage.rings.finite_rings
            Infinite order in Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(QQ); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: F.equation_order_infinite()
            Infinite order in Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2
        """
        from .order_basis import FunctionFieldOrderInfinite_basis
        b = self.primitive_integal_element_infinite()
        basis = [b**i for i in range(self.degree())]
        return FunctionFieldOrderInfinite_basis(tuple(basis))


class FunctionField_char_zero_integral(FunctionField_char_zero, FunctionField_integral):
    """
    Function fields of characteristic zero, defined by an irreducible and
    separable polynomial, integral over the maximal order of the base rational
    function field with a finite constant field.
    """
    pass


class FunctionField_global_integral(FunctionField_global, FunctionField_integral):
    """
    Global function fields, defined by an irreducible and separable polynomial,
    integral over the maximal order of the base rational function field with a
    finite constant field.
    """
    pass
