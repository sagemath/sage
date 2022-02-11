# -*- coding: utf-8 -*-
r"""
Function Fields

A function field (of one variable) is a finitely generated field extension of
transcendence degree one. In Sage, a function field can be a rational function
field or a finite extension of a function field.

EXAMPLES:

We create a rational function field::

    sage: K.<x> = FunctionField(GF(5^2,'a')); K
    Rational function field in x over Finite Field in a of size 5^2
    sage: K.genus()
    0
    sage: f = (x^2 + x + 1) / (x^3 + 1)
    sage: f
    (x^2 + x + 1)/(x^3 + 1)
    sage: f^3
    (x^6 + 3*x^5 + x^4 + 2*x^3 + x^2 + 3*x + 1)/(x^9 + 3*x^6 + 3*x^3 + 1)

Then we create an extension of the rational function field, and do some
simple arithmetic in it::

    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x)); L
    Function field in y defined by y^3 + 3*x*y + (4*x^4 + 4)/x
    sage: y^2
    y^2
    sage: y^3
    2*x*y + (x^4 + 1)/x
    sage: a = 1/y; a
    (x/(x^4 + 1))*y^2 + 3*x^2/(x^4 + 1)
    sage: a * y
    1

We next make an extension of the above function field, illustrating
that arithmetic with a tower of three fields is fully supported::

    sage: S.<t> = L[]
    sage: M.<t> = L.extension(t^2 - x*y)
    sage: M
    Function field in t defined by t^2 + 4*x*y
    sage: t^2
    x*y
    sage: 1/t
    ((1/(x^4 + 1))*y^2 + 3*x/(x^4 + 1))*t
    sage: M.base_field()
    Function field in y defined by y^3 + 3*x*y + (4*x^4 + 4)/x
    sage: M.base_field().base_field()
    Rational function field in x over Finite Field in a of size 5^2

It is also possible to construct function fields over an imperfect base field::

    sage: N.<u> = FunctionField(K)

and inseparable extension function fields::

    sage: J.<x> = FunctionField(GF(5)); J
    Rational function field in x over Finite Field of size 5
    sage: T.<v> = J[]
    sage: O.<v> = J.extension(v^5 - x); O
    Function field in v defined by v^5 + 4*x

Function fields over the rational field are supported::

    sage: F.<x> = FunctionField(QQ)
    sage: R.<Y> = F[]
    sage: L.<y> = F.extension(Y^2 - x^8 - 1)
    sage: O = L.maximal_order()
    sage: I = O.ideal(x, y - 1)
    sage: P = I.place()
    sage: D = P.divisor()
    sage: D.basis_function_space()
    [1]
    sage: (2*D).basis_function_space()
    [1]
    sage: (3*D).basis_function_space()
    [1]
    sage: (4*D).basis_function_space()
    [1, 1/x^4*y + 1/x^4]

    sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
    sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
    sage: O = F.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    2*Place (x, y, (1/(x^3 + x^2 + x))*y^2)
     + 2*Place (x^2 + x + 1, y, (1/(x^3 + x^2 + x))*y^2)

    sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: O = L.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    - Place (x, x*y)
     + Place (x^2 + 1, x*y)

Function fields over the algebraic field are supported::

    sage: K.<x> = FunctionField(QQbar); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: O = L.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    Place (x - I, x*y)
     - Place (x, x*y)
     + Place (x + I, x*y)
    sage: pl = I.divisor().support()[0]
    sage: m = L.completion(pl, prec=5)
    sage: m(x)
    I + s + O(s^5)
    sage: m(y)
    -2*s + (-4 - I)*s^2 + (-15 - 4*I)*s^3 + (-75 - 23*I)*s^4 + (-413 - 154*I)*s^5 + O(s^6)
    sage: m(y)^2 + m(y) + m(x) + 1/m(x)
    O(s^5)

TESTS::

    sage: TestSuite(J).run()
    sage: TestSuite(K).run(max_runs=1024) # long time (5s)
    sage: TestSuite(L).run(max_runs=64)   # long time (10s)
    sage: TestSuite(M).run(max_runs=32)   # long time (30s)
    sage: TestSuite(N).run(max_runs=64, skip = '_test_derivation')  # long time (8s)
    sage: TestSuite(O).run(max_runs=128, skip = '_test_derivation') # long time (8s)

    sage: TestSuite(R).run()
    sage: TestSuite(S).run() # long time (3s)

Global function fields
----------------------

A global function field in Sage is an extension field of a rational function field
over a *finite* constant field by an irreducible separable polynomial over the
rational function field.

EXAMPLES:

A fundamental computation for a global or any function field is to get a basis
of its maximal order and maximal infinite order, and then do arithmetic with
ideals of those maximal orders::

    sage: K.<x> = FunctionField(GF(3)); _.<t> = K[]
    sage: L.<y> = K.extension(t^4 + t - x^5)
    sage: O = L.maximal_order()
    sage: O.basis()
    (1, y, 1/x*y^2 + 1/x*y, 1/x^3*y^3 + 2/x^3*y^2 + 1/x^3*y)
    sage: I = O.ideal(x,y); I
    Ideal (x, y) of Maximal order of Function field in y defined by y^4 + y + 2*x^5
    sage: J = I^-1
    sage: J.basis_matrix()
    [  1   0   0   0]
    [1/x 1/x   0   0]
    [  0   0   1   0]
    [  0   0   0   1]
    sage: L.maximal_order_infinite().basis()
    (1, 1/x^2*y, 1/x^3*y^2, 1/x^4*y^3)

As an example of the most sophisticated computations that Sage can do with a
global function field, we compute all the Weierstrass places of the Klein
quartic over `\GF{2}` and gap numbers for ordinary places::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
    sage: L.genus()
    3
    sage: L.weierstrass_places()
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
    sage: L.gaps()
    [1, 2, 3]

The gap numbers for Weierstrass places are of course not ordinary::

    sage: p1,p2,p3 = L.weierstrass_places()[:3]
    sage: p1.gaps()
    [1, 2, 4]
    sage: p2.gaps()
    [1, 2, 4]
    sage: p3.gaps()
    [1, 2, 4]

AUTHORS:

- William Stein (2010): initial version

- Robert Bradshaw (2010-05-30): added is_finite()

- Julian Rüth (2011-06-08, 2011-09-14, 2014-06-23, 2014-06-24, 2016-11-13):
  fixed hom(), extension(); use @cached_method; added derivation(); added
  support for relative vector spaces; fixed conversion to base fields

- Maarten Derickx (2011-09-11): added doctests

- Syed Ahmad Lavasani (2011-12-16): added genus(), is_RationalFunctionField()

- Simon King (2014-10-29): Use the same generator names for a function field
  extension and the underlying polynomial ring.

- Kwankyu Lee (2017-04-30): added global function fields

- Brent Baccala (2019-12-20): added function fields over number fields and QQbar

"""
# ****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2011-2018 Julian Rüth <julian.rueth@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method

from sage.interfaces.all import singular

from sage.arith.all import lcm

from sage.rings.ring import Field
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar_decorators import handle_AA_and_QQbar

from sage.modules.free_module_element import vector

from sage.categories.homset import Hom
from sage.categories.function_fields import FunctionFields
from sage.structure.category_object import CategoryObject

from .differential import DifferentialsSpace, DifferentialsSpace_global

from .element import (
    FunctionFieldElement,
    FunctionFieldElement_rational,
    FunctionFieldElement_polymod)

def is_FunctionField(x):
    """
    Return ``True`` if ``x`` is a function field.

    EXAMPLES::

        sage: from sage.rings.function_field.function_field import is_FunctionField
        sage: is_FunctionField(QQ)
        False
        sage: is_FunctionField(FunctionField(QQ, 't'))
        True
    """
    if isinstance(x, FunctionField):
        return True
    return x in FunctionFields()


class FunctionField(Field):
    """
    Abstract base class for all function fields.

    INPUT:

    - ``base_field`` -- field; the base of this function field

    - ``names`` -- string that gives the name of the generator

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: K
        Rational function field in x over Rational Field
    """
    _differentials_space = DifferentialsSpace

    def __init__(self, base_field, names, category=FunctionFields()):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: TestSuite(K).run()
        """
        Field.__init__(self, base_field, names=names, category=category)

        # allow conversion into the constant base field
        from .maps import FunctionFieldConversionToConstantBaseField
        to_constant_base_field = FunctionFieldConversionToConstantBaseField(Hom(self, self.constant_base_field()))
        # the conversion map must not keep the field alive if that is the only reference to it
        to_constant_base_field._make_weak_references()
        self.constant_base_field().register_conversion(to_constant_base_field)

    def is_perfect(self):
        r"""
        Return whether the field is perfect, i.e., its characteristic `p` is zero
        or every element has a `p`-th root.

        EXAMPLES::

            sage: FunctionField(QQ, 'x').is_perfect()
            True
            sage: FunctionField(GF(2), 'x').is_perfect()
            False
        """
        return self.characteristic() == 0

    def some_elements(self):
        """
        Return some elements in this function field.

        EXAMPLES::

           sage: K.<x> = FunctionField(QQ)
           sage: K.some_elements()
           [1,
            x,
            2*x,
            x/(x^2 + 2*x + 1),
            1/x^2,
            x/(x^2 - 1),
            x/(x^2 + 1),
            1/2*x/(x^2 + 1),
            0,
            1/x,
            ...]

        ::

           sage: R.<y> = K[]
           sage: L.<y> = K.extension(y^2 - x)
           sage: L.some_elements()
           [1,
            y,
            1/x*y,
            ((x + 1)/(x^2 - 2*x + 1))*y - 2*x/(x^2 - 2*x + 1),
            1/x,
            (1/(x - 1))*y,
            (1/(x + 1))*y,
            (1/2/(x + 1))*y,
            0,
            ...]
        """
        elements = []

        polynomials = [self(f) for f in self._ring.some_elements()]

        for numerator in polynomials:
            for denominator in polynomials:
                if denominator:
                    some_element = numerator/denominator
                    if some_element not in elements:
                        elements.append(some_element)

        return elements

    def characteristic(self):
        """
        Return the characteristic of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.characteristic()
            0
            sage: K.<x> = FunctionField(QQbar)
            sage: K.characteristic()
            0
            sage: K.<x> = FunctionField(GF(7))
            sage: K.characteristic()
            7
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.characteristic()
            7
        """
        return self.constant_base_field().characteristic()

    def is_finite(self):
        """
        Return whether the function field is finite, which is false.

        EXAMPLES::

            sage: R.<t> = FunctionField(QQ)
            sage: R.is_finite()
            False
            sage: R.<t> = FunctionField(GF(7))
            sage: R.is_finite()
            False
        """
        return False

    def is_global(self):
        """
        Return whether the function field is global, that is, whether
        the constant field is finite.

        EXAMPLES::

            sage: R.<t> = FunctionField(QQ)
            sage: R.is_global()
            False
            sage: R.<t> = FunctionField(QQbar)
            sage: R.is_global()
            False
            sage: R.<t> = FunctionField(GF(7))
            sage: R.is_global()
            True
        """
        return self.constant_base_field().is_finite()

    def extension(self, f, names=None):
        """
        Create an extension `K(y)` of this function field `K` extended with
        a root `y` of the univariate polynomial `f` over `K`.

        INPUT:

        - ``f`` -- univariate polynomial over `K`

        - ``names`` -- string or tuple of length 1 that names the variable `y`

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1), 'z')
            Function field in z defined by z^3 + 1/t*z + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names)

    def order_with_basis(self, basis, check=True):
        """
        Return the order with given basis over the maximal order of
        the base field.

        INPUT:

        - ``basis`` -- list of elements of this function field

        - ``check`` -- boolean (default: ``True``); if ``True``, check that the
          basis is really linearly independent and that the module it spans is
          closed under multiplication, and contains the identity element.

        OUTPUT:

        - an order in the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order_with_basis([1, y, y^2]); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, y, y^2)

        Note that 1 does not need to be an element of the basis, as long it is in the module spanned by it::

            sage: O = L.order_with_basis([1+y, y, y^2]); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (y + 1, y, y^2)

        The following error is raised when the module spanned by the basis is not closed under multiplication::

            sage: O = L.order_with_basis([1, x^2 + x*y, (2/3)*y^2]); O
            Traceback (most recent call last):
            ...
            ValueError: the module generated by basis (1, x*y + x^2, 2/3*y^2) must be closed under multiplication

        and this happens when the identity is not in the module spanned by the basis::

            sage: O = L.order_with_basis([x, x^2 + x*y, (2/3)*y^2])
            Traceback (most recent call last):
            ...
            ValueError: the identity element must be in the module spanned by basis (x, x*y + x^2, 2/3*y^2)
        """
        from .order import FunctionFieldOrder_basis
        return FunctionFieldOrder_basis(tuple([self(a) for a in basis]), check=check)

    def order(self, x, check=True):
        """
        Return the order generated by ``x`` over the base maximal order.

        INPUT:

        - ``x`` -- element or list of elements of the function field

        - ``check`` -- boolean; if ``True``, check that ``x`` really generates an order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order(y); O
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, y, y^2)

            sage: Z = K.order(x); Z
            Order in Rational function field in x over Rational Field
            sage: Z.basis()
            (1,)

        Orders with multiple generators are not yet supported::

            sage: Z = K.order([x,x^2]); Z
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(x, (list, tuple)):
            x = [x]
        if len(x) == 1:
            g = x[0]
            basis = [self(1)]
            for i in range(self.degree()-1):
                basis.append(basis[-1]*g)
        else:
            raise NotImplementedError
        return self.order_with_basis(basis, check=check)

    def order_infinite_with_basis(self, basis, check=True):
        """
        Return the order with given basis over the maximal infinite order of
        the base field.

        INPUT:

        - ``basis`` -- list of elements of the function field

        - ``check`` -- boolean (default: ``True``); if ``True``, check that the basis
          is really linearly independent and that the module it spans is closed
          under multiplication, and contains the identity element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order_infinite_with_basis([1, 1/x*y, 1/x^2*y^2]); O
            Infinite order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, 1/x*y, 1/x^2*y^2)

        Note that 1 does not need to be an element of the basis, as long it is
        in the module spanned by it::

            sage: O = L.order_infinite_with_basis([1+1/x*y,1/x*y, 1/x^2*y^2]); O
            Infinite order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1/x*y + 1, 1/x*y, 1/x^2*y^2)

        The following error is raised when the module spanned by the basis is
        not closed under multiplication::

            sage: O = L.order_infinite_with_basis([1,y, 1/x^2*y^2]); O
            Traceback (most recent call last):
            ...
            ValueError: the module generated by basis (1, y, 1/x^2*y^2) must be closed under multiplication

        and this happens when the identity is not in the module spanned by the
        basis::

            sage: O = L.order_infinite_with_basis([1/x,1/x*y, 1/x^2*y^2])
            Traceback (most recent call last):
            ...
            ValueError: the identity element must be in the module spanned by basis (1/x, 1/x*y, 1/x^2*y^2)
        """
        from .order import FunctionFieldOrderInfinite_basis
        return FunctionFieldOrderInfinite_basis(tuple([self(g) for g in basis]), check=check)

    def order_infinite(self, x, check=True):
        """
        Return the order generated by ``x`` over the maximal infinite order.

        INPUT:

        - ``x`` -- element or a list of elements of the function field

        - ``check`` -- boolean; if ``True``, check that ``x`` really generates an order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: L.order_infinite(y)  # todo: not implemented

            sage: Z = K.order(x); Z
            Order in Rational function field in x over Rational Field
            sage: Z.basis()
            (1,)

        Orders with multiple generators, not yet supported::

            sage: Z = K.order_infinite([x,x^2]); Z
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(x, (list, tuple)):
            x = [x]
        if len(x) == 1:
            g = x[0]
            basis = [self(1)]
            for i in range(self.degree()-1):
                basis.append(basis[-1]*g)
        else:
            raise NotImplementedError
        return self.order_infinite_with_basis(tuple(basis), check=check)

    def _coerce_map_from_(self, source):
        """
        Return ``True`` if there is a coerce map from ``R`` to the function field.

        INPUT:

        - ``source`` -- ring

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: L.equation_order()
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: L._coerce_map_from_(L.equation_order())
            Conversion map:
              From: Order in Function field in y defined by y^3 + x^3 + 4*x + 1
              To:   Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: L._coerce_map_from_(GF(7))

            sage: K.<x> = FunctionField(QQ)
            sage: L.<x> = FunctionField(GaussianIntegers().fraction_field())
            sage: L.has_coerce_map_from(K)
            True

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 1)
            sage: K.<x> = FunctionField(GaussianIntegers().fraction_field())
            sage: R.<y> = K[]
            sage: M.<y> = K.extension(y^3 + 1)
            sage: M.has_coerce_map_from(L) # not tested (the constant field including into a function field is not yet known to be injective)
            True

            sage: K.<x> = FunctionField(QQ)
            sage: R.<I> = K[]
            sage: L.<I> = K.extension(I^2 + 1)
            sage: M.<x> = FunctionField(GaussianIntegers().fraction_field())
            sage: M.has_coerce_map_from(L)
            True

        Check that :trac:`31072` is fixed::

            sage: L.<t> = FunctionField(QQ)
            sage: L(Sequence([1, 2]))
            2*t + 1
        """
        from .order import FunctionFieldOrder_base
        if isinstance(source, FunctionFieldOrder_base):
            K = source.fraction_field()
            if K is self:
                return self._generic_coerce_map(source)
            source_to_K = K.coerce_map_from(source)
            K_to_self = self.coerce_map_from(K)
            if source_to_K and K_to_self:
                return K_to_self * source_to_K
        if isinstance(source, CategoryObject) and source in FunctionFields():
            if source.base_field() is source:
                if self.base_field() is self:
                    # source and self are rational function fields
                    if source.variable_name() == self.variable_name():
                        # ... in the same variable
                        base_coercion = self.constant_field().coerce_map_from(source.constant_field())
                        if base_coercion is not None:
                            return source.hom([self.gen()], base_morphism=base_coercion)
            else:
                # source is an extensions of rational function fields
                base_coercion = self.coerce_map_from(source.base_field())
                if base_coercion is not None and base_coercion.is_injective():
                    # the base field of source coerces into the base field of self
                    self_polynomial = source.polynomial().map_coefficients(base_coercion)
                    # try to find a root of the defining polynomial in self
                    if self_polynomial(self.gen()) == 0:
                        # The defining polynomial of source has a root in self,
                        # therefore there is a map. To be sure that it is
                        # canonical, we require a root of the defining polynomial
                        # of self to be a root of the defining polynomial of
                        # source (and that the variables are named equally):
                        if source.variable_name() == self.variable_name():
                            return source.hom([self.gen()], base_morphism=base_coercion)

                    try:
                        sourcegen_in_self = self(source.variable_name())
                    except TypeError:
                        pass
                    else:
                        if self_polynomial(sourcegen_in_self) == 0:
                            # The defining polynomial of source has a root in self,
                            # therefore there is a map. To be sure that it is
                            # canonical, we require the names of the roots to match
                            return source.hom([sourcegen_in_self], base_morphism=base_coercion)

    def _test_derivation(self, **options):
        """
        Test the correctness of the derivations of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: TestSuite(K).run() # indirect doctest
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        K = self.constant_base_field().some_elements()
        d = self.derivation()
        from itertools import product
        # Non-zero
        tester.assertFalse(d.is_zero())
        # Well-defined
        if hasattr(self, "polynomial"):
            f = self.polynomial()
            tester.assertEqual(0, d(f))
        # Leibniz's law
        for x,y in tester.some_elements(product(S, S)):
            tester.assertEqual(d(x*y), x*d(y) + d(x)*y)
        # Linearity
        for x,y in tester.some_elements(product(S, S)):
            tester.assertEqual(d(x+y), d(x) + d(y))
        for c,x in tester.some_elements(product(K, S)):
            tester.assertEqual(d(c*x), c*d(x))
        # Constants map to zero
        for c in tester.some_elements(K):
            tester.assertEqual(d(c), 0)

    def _convert_map_from_(self, R):
        """
        Return a conversion from ``R`` to this function field if one exists.

        INPUT:

        - ``R`` -- ring

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: K(L(x)) # indirect doctest
            x
        """
        if isinstance(R, FunctionField_polymod):
            base_conversion = self.convert_map_from(R.base_field())
            if base_conversion is not None:
                from sage.categories.morphism import SetMorphism
                return base_conversion * SetMorphism(R.Hom(R.base_field()), R._to_base_field)

    def _intermediate_fields(self, base):
        """
        Return the fields which lie in between base and the function field in the
        tower of function fields.

        INPUT:

        - ``base`` -- function field, either this field or a field from which
          this field has been created as an extension

        OUTPUT:

        - a list of fields; the first entry is this field, the last entry is ``base``

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._intermediate_fields(K)
            [Rational function field in x over Rational Field]

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: L._intermediate_fields(K)
            [Function field in y defined by y^2 - x, Rational function field in x over Rational Field]

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M._intermediate_fields(L)
            [Function field in z defined by z^2 - y, Function field in y defined by y^2 - x]
            sage: M._intermediate_fields(K)
            [Function field in z defined by z^2 - y, Function field in y defined by y^2 - x,
            Rational function field in x over Rational Field]

        TESTS::

            sage: K._intermediate_fields(M)
            Traceback (most recent call last):
            ...
            ValueError: field has not been constructed as a finite extension of base
            sage: K._intermediate_fields(QQ)
            Traceback (most recent call last):
            ...
            TypeError: base must be a function field
        """
        if not is_FunctionField(base):
            raise TypeError("base must be a function field")

        ret = [self]
        while ret[-1] is not base:
            ret.append(ret[-1].base_field())
            if ret[-1] is ret[-2]:
                raise ValueError("field has not been constructed as a finite extension of base")
        return ret

    def rational_function_field(self):
        r"""
        Return the rational function field from which this field has been
        created as an extension.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.rational_function_field()
            Rational function field in x over Rational Field

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: L.rational_function_field()
            Rational function field in x over Rational Field

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M.rational_function_field()
            Rational function field in x over Rational Field
        """
        return self if isinstance(self, RationalFunctionField) else self.base_field().rational_function_field()

    def valuation(self, prime):
        r"""
        Return the discrete valuation on this function field defined by
        ``prime``.

        INPUT:

        - ``prime`` -- a place of the function field, a valuation on a subring,
          or a valuation on another function field together with information
          for isomorphisms to and from that function field

        EXAMPLES:

        We create valuations that correspond to finite rational places of a
        function field::

            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(1); v
            (x - 1)-adic valuation
            sage: v(x)
            0
            sage: v(x - 1)
            1

        A place can also be specified with an irreducible polynomial::

            sage: v = K.valuation(x - 1); v
            (x - 1)-adic valuation

        Similarly, for a finite non-rational place::

            sage: v = K.valuation(x^2 + 1); v
            (x^2 + 1)-adic valuation
            sage: v(x^2 + 1)
            1
            sage: v(x)
            0

        Or for the infinite place::

            sage: v = K.valuation(1/x); v
            Valuation at the infinite place
            sage: v(x)
            -1

        Instead of specifying a generator of a place, we can define a valuation on a
        rational function field by giving a discrete valuation on the underlying
        polynomial ring::

            sage: R.<x> = QQ[]
            sage: w = valuations.GaussValuation(R, valuations.TrivialValuation(QQ)).augmentation(x - 1, 1)
            sage: v = K.valuation(w); v
            (x - 1)-adic valuation

        Note that this allows us to specify valuations which do not correspond to a
        place of the function field::

            sage: w = valuations.GaussValuation(R, QQ.valuation(2))
            sage: v = K.valuation(w); v
            2-adic valuation

        The same is possible for valuations with `v(1/x) > 0` by passing in an
        extra pair of parameters, an isomorphism between this function field and an
        isomorphic function field. That way you can, for example, indicate that the
        valuation is to be understood as a valuation on `K[1/x]`, i.e., after
        applying the substitution `x \mapsto 1/x` (here, the inverse map is also `x
        \mapsto 1/x`)::

            sage: w = valuations.GaussValuation(R, QQ.valuation(2)).augmentation(x, 1)
            sage: w = K.valuation(w)
            sage: v = K.valuation((w, K.hom([~K.gen()]), K.hom([~K.gen()]))); v
            Valuation on rational function field induced by [ Gauss valuation induced by 2-adic valuation, v(x) = 1 ] (in Rational function field in x over Rational Field after x |--> 1/x)

        Note that classical valuations at finite places or the infinite place are
        always normalized such that the uniformizing element has valuation 1::

            sage: K.<t> = FunctionField(GF(3))
            sage: M.<x> = FunctionField(K)
            sage: v = M.valuation(x^3 - t)
            sage: v(x^3 - t)
            1

        However, if such a valuation comes out of a base change of the ground
        field, this is not the case anymore. In the example below, the unique
        extension of ``v`` to ``L`` still has valuation 1 on `x^3 - t` but it has
        valuation ``1/3`` on its uniformizing element  `x - w`::

            sage: R.<w> = K[]
            sage: L.<w> = K.extension(w^3 - t)
            sage: N.<x> = FunctionField(L)
            sage: w = v.extension(N) # missing factorization, :trac:`16572`
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: w(x^3 - t) # not tested
            1
            sage: w(x - w) # not tested
            1/3

        There are several ways to create valuations on extensions of rational
        function fields::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x); L
            Function field in y defined by y^2 - x

        A place that has a unique extension can just be defined downstairs::

            sage: v = L.valuation(x); v
            (x)-adic valuation

        """
        from sage.rings.function_field.function_field_valuation import FunctionFieldValuation
        return FunctionFieldValuation(self, prime)

    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.space_of_differentials()
            Space of differentials of Rational function field in t over Rational Field

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.space_of_differentials()
            Space of differentials of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        return self._differentials_space(self)

    def divisor_group(self):
        """
        Return the group of divisors attached to the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.divisor_group()
            Divisor group of Rational function field in t over Rational Field

            sage: _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (t^3 - 1)/(t^3 - 2))
            sage: L.divisor_group()
            Divisor group of Function field in y defined by y^3 + (-t^3 + 1)/(t^3 - 2)

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.divisor_group()
            Divisor group of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        from .divisor import DivisorGroup
        return DivisorGroup(self)

    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: K.place_set()
            Set of places of Rational function field in t over Finite Field of size 7

            sage: K.<t> = FunctionField(QQ)
            sage: K.place_set()
            Set of places of Rational function field in t over Rational Field

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.place_set()
            Set of places of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def completion(self, place, name=None, prec=None, gen_name=None):
        """
        Return the completion of the function field at the place.

        INPUT:

        - ``place`` -- place

        - ``name`` -- string; name of the series variable

        - ``prec`` -- positive integer; default precision

        - ``gen_name`` -- string; name of the generator of the residue field;
          used only when the place is non-rational

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p); m
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(x,10)
            s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + s^9 + s^10 + O(s^12)
            sage: m(y,10)
            s^-1 + 1 + s^3 + s^5 + s^7 + O(s^9)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p); m
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(x,10)
            s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + s^9 + s^10 + O(s^12)
            sage: m(y,10)
            s^-1 + 1 + s^3 + s^5 + s^7 + O(s^9)

            sage: K.<x> = FunctionField(GF(2))
            sage: p = K.places_finite()[0]; p
            Place (x)
            sage: m = K.completion(p); m
            Completion map:
              From: Rational function field in x over Finite Field of size 2
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(1/(x+1))
            1 + s + s^2 + s^3 + s^4 + s^5 + s^6 + s^7 + s^8 + s^9 + s^10 + s^11 + s^12
            + s^13 + s^14 + s^15 + s^16 + s^17 + s^18 + s^19 + O(s^20)

            sage: p = K.place_infinite(); p
            Place (1/x)
            sage: m = K.completion(p); m
            Completion map:
              From: Rational function field in x over Finite Field of size 2
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(x)
            s^-1 + O(s^19)

            sage: m = K.completion(p, prec=infinity); m
            Completion map:
              From: Rational function field in x over Finite Field of size 2
              To:   Lazy Laurent Series Ring in s over Finite Field of size 2
            sage: f = m(x); f
            s^-1 + ...
            sage: f.coefficient(100)
            0

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x)
            sage: O = L.maximal_order()
            sage: decomp = O.decomposition(K.maximal_order().ideal(x - 1))
            sage: pls = (decomp[0][0].place(), decomp[1][0].place())
            sage: m = L.completion(pls[0]); m
            Completion map:
              From: Function field in y defined by y^2 - x
              To:   Laurent Series Ring in s over Rational Field
            sage: xe = m(x)
            sage: ye = m(y)
            sage: ye^2 - xe == 0
            True

            sage: decomp2 = O.decomposition(K.maximal_order().ideal(x^2 + 1))
            sage: pls2 = decomp2[0][0].place()
            sage: m = L.completion(pls2); m
            Completion map:
              From: Function field in y defined by y^2 - x
              To:   Laurent Series Ring in s over Number Field in a with defining polynomial x^4 + 2*x^2 + 4*x + 2
            sage: xe = m(x)
            sage: ye = m(y)
            sage: ye^2 - xe == 0
            True
        """
        from .maps import FunctionFieldCompletion
        return FunctionFieldCompletion(self, place, name=name, prec=prec, gen_name=gen_name)


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
            sage: TestSuite(L).run()  # long time

        We can set the variable name, which doesn't have to be y::

            sage: L.<w> = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in w defined by w^5 + x*w - x^3 - 3*x

        TESTS:

        Test that :trac:`17033` is fixed::

            sage: K.<t> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: M.<z> = K.extension(x^7-x-t)
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
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        if polynomial.parent().ngens()>1 or not is_Polynomial(polynomial):
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

        from .place import FunctionFieldPlace_polymod
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

        - ``f`` -- element of the function field which lies in the base
          field.

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

        Verify that :trac:`21872` has been resolved::

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
        raise ValueError("%r is not an element of the base field"%(f,))

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

        Verify that :trac:`21872` has been resolved::

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

        - ``names`` -- a string or a tuple of up to two strings (default:
          ``None``), the name of the generator of the field, and the name of
          the generator of the underlying rational function field (if a tuple);
          if not given, then the names are chosen automatically.

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
            (Function field in z_ defined by z_^10 - x^18, Function Field morphism:
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
            sage: L.<y> = K.extension(y^2-x)
            sage: L.monic_integral_model()
            (Function field in y defined by y^2 - x, Function Field endomorphism of Function field in y defined by y^2 - x
              Defn: y |--> y
                    x |--> x, Function Field endomorphism of Function field in y defined by y^2 - x
              Defn: y |--> y
                    x |--> x)

        unless ``names`` does not match with the current names::

            sage: L.monic_integral_model(names=('yy','xx'))
            (Function field in yy defined by yy^2 - xx, Function Field morphism:
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
            L,from_L,to_L = self.simple_model()
            ret,ret_to_L,L_to_ret = L.monic_integral_model(names)
            from_ret = ret.hom( [from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))] )
            to_ret = self.hom( [L_to_ret(to_L(k.gen())) for k in self._intermediate_fields(self.rational_function_field())] )
            return ret, from_ret, to_ret
        else:
            if self.polynomial().is_monic() and all(c.denominator().is_one() for c in self.polynomial()):
                # self is already monic and integral
                if names is None or names == ():
                    names = (self.variable_name(),)
                return self.change_variable_name(names)
            else:
                if not names:
                    names = (self.variable_name()+"_",)
                if len(names) == 1:
                    names = (names[0], self.rational_function_field().variable_name())

                g, d = self._make_monic_integral(self.polynomial())
                K,from_K,to_K = self.base_field().change_variable_name(names[1])
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

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^5 - x)
            sage: L.constant_field()
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
          from which this field has been constructed as a finite extension.

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
        return "Function field in %s defined by %s"%(self.variable_name(), self._polynomial)

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

            sage: K.<x> = FunctionField(GF(5))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.is_separable()
            True

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

        - ``basis`` -- a basis for this field over the base.

        - ``maps`` -- boolean (default ``True``), whether to return
          `R`-linear maps to and from `V`.

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

            sage: from_V(V.1)
            y

        We define an interesting element of the function field::

            sage: a = 1/L.0; a
            (x/(x^4 + 1))*y^4 - 2*x^2/(x^4 + 1)

        We convert it to the vector space, and get a vector over the base field::

            sage: to_V(a)
            (-2*x^2/(x^4 + 1), 0, 0, 0, x/(x^4 + 1))

        We convert to and back, and get the same element::

            sage: from_V(to_V(a)) == a
            True

        In the other direction::

            sage: v = x*V.0 + (1/x)*V.1
            sage: to_V(from_V(v)) == v
            True

        And we show how it works over an extension of an extension field::

            sage: R2.<z> = L[]; M.<z> = L.extension(z^2 -y)
            sage: M.free_module()
            (Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x, Isomorphism:
              From: Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
              To:   Function field in z defined by z^2 - y, Isomorphism:
              From: Function field in z defined by z^2 - y
              To:   Vector space of dimension 2 over Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x)

        We can also get the vector space of ``M`` over ``K``::

            sage: M.free_module(K)
            (Vector space of dimension 10 over Rational function field in x over Rational Field, Isomorphism:
              From: Vector space of dimension 10 over Rational function field in x over Rational Field
              To:   Function field in z defined by z^2 - y, Isomorphism:
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
        to_V   = MapFunctionFieldToVectorSpace(self, V)
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
        from .order import FunctionFieldMaximalOrder_polymod
        return FunctionFieldMaximalOrder_polymod(self)

    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: F.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .order import FunctionFieldMaximalOrderInfinite_polymod
        return FunctionFieldMaximalOrderInfinite_polymod(self)

    def different(self):
        """
        Return the different of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: F.different()
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
          and of successive base rings.

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

            sage: f = L.hom(y+1)
            Traceback (most recent call last):
            ...
            ValueError: invalid morphism

        We make a morphism of the base rational function field::

            sage: phi = K.hom(x+1); phi
            Function Field endomorphism of Rational function field in x over Rational Field
              Defn: x |--> x + 1
            sage: phi(x^3 - 3)
            x^3 + 3*x^2 + 3*x - 2
            sage: (x+1)^3-3
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

            sage: f = L.hom([4*w, t+1]); f
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
        if not isinstance(im_gens, (list,tuple)):
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

            # making the auxiliary ring which only has polynomials
            # with integral coefficients.
            tmpAuxRing = PolynomialRing(self._base_field.constant_field(),
                            str(self._base_field.gen())+','+str(self._ring.gen()))
            intMinPoly, d = self._make_monic_integral(self._polynomial)
            curveIdeal = tmpAuxRing.ideal(intMinPoly)

            singular.lib('normal.lib') #loading genus method in Singular
            return int(curveIdeal._singular_().genus())

        else:
            raise NotImplementedError("computation of genus over non-prime "
                                      "constant fields not implemented yet")

    @cached_method
    def derivation(self):
        r"""
        Return a generator of the space of derivations over the constant base
        ring of this function field `K`.

        A derivation on `K` is map `K\to K` with
        `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
        D(\alpha)+\alpha D(\beta)` for all `\alpha,\beta\in K`.

        If the base field `k` of `K` is perfect, then the derivations on `K`
        form a one-dimensional `K`-vector space. (More generally, this is true
        if `K` is separable over `k`, or in other words if the corresponding
        curve is geometrically reduced over `k`; this is automatically the case
        if `k` is perfect.) We apply the techniques from [GT1996]_ to find a
        generator of this one-dimensional vector space, which is then returned
        by the algorithm.

        ALGORITHM:

        If `K` is a separable extension of another function field `F` between
        `K` and `k`, then Proposition 11 of [GT1996]_ describes how to compute
        the unique extension of a derivation on `F` to `K`; we then apply this
        algorithm to the generator of the space of derivations on `F`, which we
        may calculate inductively.
        If `K` is not given as a separable extension of another function field,
        then we find a field isomorphic to `K` that is a separable extension of
        a rational function field over `k` by using :meth:`separable_model`.
        This part of the algorithm uses the assumption that `k` is perfect.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation(); d
            Derivation map:
                From: Function field in y defined by y^2 + 2*x
                To:   Function field in y defined by y^2 + 2*x
                Defn: y |--> 2/x*y
                      x |--> 1
            sage: d(x)
            1
            sage: d(x^3)
            0
            sage: d(x*y)
            0
            sage: d(y)
            2/x*y

        This also works for inseparable extensions::

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - x)
            sage: d = L.derivation(); d
            Derivation map:
                From: Function field in y defined by y^3 + 2*x
                To:   Function field in y defined by y^3 + 2*x
                Defn: y |--> 1
                      x |--> 0
            sage: d(x^2)
            0
            sage: d(y^2)
            2*y
            sage: d(x*y)
            x

        """
        from .maps import FunctionFieldDerivation_separable, FunctionFieldDerivation_inseparable
        if self.is_separable():
            return FunctionFieldDerivation_separable(self, self.base_ring().derivation())
        else:
            return FunctionFieldDerivation_inseparable(self)

    def _simple_model(self, name='v'):
        r"""
        Return a finite extension `N/K(x)` isomorphic to the tower of
        extensions `M/L/K(x)` with `K` perfect.

        Helper method for :meth:`simple_model`.

        INPUT:

        - ``name`` -- a string, the name of the generator of `N`

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
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
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

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
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

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3-1)
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

        assert(isinstance(K, RationalFunctionField))
        assert(K is not L)
        assert(L is not M)

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
        M_to_N = M.hom([M_a,M_b])

        return N, N_to_M, M_to_N

    @cached_method
    def simple_model(self, name=None):
        """
        Return a function field isomorphic to this field which is a simple
        extension of a rational function field.

        INPUT:

        - ``name`` -- a string (default: ``None``), the name of generator of
          the simple extension. If ``None``, then the name of the generator
          will be the same as the name of the generator of this function field.

        OUTPUT:

        A triple ``(F,f,t)`` where ``F`` is a field isomorphic to this field,
        ``f`` is an isomorphism from ``F`` to this function field and ``t`` is
        the inverse of ``f``.

        EXAMPLES:

        A tower of four function fields::

            sage: K.<x> = FunctionField(QQ); R.<z> = K[]
            sage: L.<z> = K.extension(z^2-x); R.<u> = L[]
            sage: M.<u> = L.extension(u^2-z); R.<v> = M[]
            sage: N.<v> = M.extension(v^2-u)

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

            sage: K.<x> = FunctionField(GF(3)); R.<y> = K[]
            sage: L.<y> = K.extension(y^5-x); R.<z> = L[]
            sage: M.<z> = L.extension(z^3-x)
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

            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2-x); R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: M.simple_model()
            (Function field in z defined by z^4 + x, Function Field morphism:
               From: Function field in z defined by z^4 + x
               To:   Function field in z defined by z^2 + y
               Defn: z |--> z, Function Field morphism:
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
                id = Hom(self,self).identity()
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
            to_self_ = self.hom([self_.gen()]+gens_in_base_)
            from_self_ = self_.hom([self.gen(),from_base_(base_.gen())])

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
            sage: L.<y> = K.extension(y^2-x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2-y)
            sage: R.<z> = L[]
            sage: N.<u> = L.extension(z^2-x-1)
            sage: N.primitive_element()
            u + y
            sage: M.primitive_element()
            z
            sage: L.primitive_element()
            y

        This also works for inseparable extensions::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<Y> = K[]
            sage: L.<y> = K.extension(Y^2-x)
            sage: R.<Z> = L[]
            sage: M.<z> = L.extension(Z^2-y)
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

        - ``names`` -- a tuple of two strings or ``None`` (default: ``None``);
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

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - t)
            sage: L.separable_model()
            Traceback (most recent call last):
            ...
            NotImplementedError: constructing a separable model is only implemented for function fields over a perfect constant base field

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

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y)
            sage: M.separable_model()
            (Function field in z_ defined by z_ + x_^6, Function Field morphism:
               From: Function field in z_ defined by z_ + x_^6
               To:   Function field in z defined by z^3 + y
               Defn: z_ |--> x
                     x_ |--> z, Function Field morphism:
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
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(self.constant_base_field(), names=names)
            S = R.remove_var(names[1])
            f = R( L.polynomial().change_variable_name(names[1]).map_coefficients(
                     lambda c:c.numerator().change_variable_name(names[0]), S))
            f = f.polynomial(R.gen(0)).change_ring(K)
            f /= f.leading_coefficient()
            # f must be separable in the other variable (otherwise it would factor)
            assert f.gcd(f.derivative()).is_one()

            ret = K.extension(f, names=(names[0],))
            # isomorphisms between L and ret are given by swapping generators
            ret_to_L = ret.hom( [L(L.base_field().gen()), L.gen()] )
            L_to_ret = L.hom( [ret(K.gen()), ret.gen()] )
            # compose with from_L and to_L to get the desired isomorphisms between self and ret
            f = ret.hom( [from_L(ret_to_L(ret.gen())), from_L(ret_to_L(ret.base_field().gen()))] )
            t = self.hom( [L_to_ret(to_L(self.gen())), L_to_ret(to_L(self.base_field().gen()))] )
            return ret, f, t

    def change_variable_name(self, name):
        r"""
        Return a field isomorphic to this field with variable(s) ``name``.

        INPUT:

        - ``name`` -- a string or a tuple consisting of a strings, the names of
          the new variables starting with a generator of this field and going
          down to the rational function field.

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
            (Function field in zz defined by zz^2 - yy, Function Field morphism:
              From: Function field in zz defined by zz^2 - yy
              To:   Function field in z defined by z^2 - y
              Defn: zz |--> z
                    yy |--> y
                    x |--> x, Function Field morphism:
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
            from_base = to_base = Hom(base,base).identity()
        else:
            base, from_base, to_base = self.base_field().change_variable_name(name[1:])

        ret = base.extension(self.polynomial().map_coefficients(to_base), names=(name[0],))
        f = ret.hom( [k.gen() for k in self._intermediate_fields(self.rational_function_field())] )
        t = self.hom( [k.gen() for k in ret._intermediate_fields(ret.rational_function_field())] )
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

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: F._inversion_isomorphism()
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
                             x |--> 1/x, Composite map:
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
        R = PolynomialRing(K,'T')
        x = K.gen()
        xinv = 1/x

        h = K.hom(xinv)
        F_poly = R([h(c) for c in self.polynomial().list()])
        F = K.extension(F_poly)

        self2F = self.hom([F.gen(),xinv])
        F2self = F.hom([self.gen(),xinv])

        M, M2F, F2M = F.monic_integral_model('s')

        return M, F2self*M2F, F2M*self2F

    def places_above(self, p):
        """
        Return places lying above ``p``.

        INPUT:

        - ``p`` -- place of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: all(q.place_below() == p for p in K.places() for q in F.places_above(p))
            True

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = K.maximal_order()
            sage: pls = [O.ideal(x-c).place() for c in [-2, -1, 0, 1, 2]]
            sage: all(q.place_below() == p for p in pls for q in F.places_above(p))
            True

            sage: K.<x> = FunctionField(QQbar); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = K.maximal_order()
            sage: pls = [O.ideal(x-QQbar(sqrt(c))).place() for c in [-2, -1, 0, 1, 2]]
            sage: all(q.place_below() == p for p in pls for q in F.places_above(p))
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

            sage: K.<x> = FunctionField(GF(3)); _.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: L.constant_field()
            Finite Field of size 3
        """
        return self.exact_constant_field()[0]

    def exact_constant_field(self, name='t'):
        """
        Return the exact constant field and its embedding into the function field.

        INPUT:

        - ``name`` -- name (default: `t`) of the generator of the exact constant field

        EXAMPLES::

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

            sage: F.<a> = GF(16)
            sage: K.<x> = FunctionField(F); K
            Rational function field in x over Finite Field in a of size 2^4
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4+t-x^5)
            sage: L.genus()
            6

        The genus is computed by the Hurwitz genus formula.
        """
        k, _ = self.exact_constant_field()
        different_degree = self.different().degree() # must be even
        return different_degree // 2 - self.degree() / k.degree() + 1

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
        exception ``TypeError`` is raised.

        EXAMPLES::

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
            sage: L.higher_derivation()
            Higher derivation map:
              From: Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
              To:   Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
        """
        from .maps import FunctionFieldHigherDerivation_char_zero
        return FunctionFieldHigherDerivation_char_zero(self)


class FunctionField_global(FunctionField_simple):
    """
    Global function fields.

    INPUT:

    - ``polynomial`` -- monic irreducible and separable polynomial

    - ``names`` -- name of the generator of the function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
        sage: L
        Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)

    The defining equation needs not be monic::

        sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
        sage: L.<y> = K.extension((1 - x)*Y^7 - x^3)
        sage: L.gaps()
        [1, 2, 3]

    or may define a trivial extension::

        sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y-1)
        sage: L.genus()
        0
    """
    _differentials_space = DifferentialsSpace_global

    def __init__(self, polynomial, names):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: TestSuite(L).run()
        """
        FunctionField_polymod.__init__(self, polynomial, names)

    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2));
            sage: R.<t> = PolynomialRing(K);
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18);
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y)
        """
        from .order import FunctionFieldMaximalOrder_global
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

            sage: K.<x>=FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.higher_derivation()
            Higher derivation map:
              From: Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
              To:   Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        from .maps import FunctionFieldHigherDerivation_global
        return FunctionFieldHigherDerivation_global(self)

    def get_place(self, degree):
        """
        Return a place of ``degree``.

        INPUT:

        - ``degree`` -- a positive integer

        OUTPUT: a place of ``degree`` if any exists; otherwise ``None``

        EXAMPLES::

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

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4+t-x^5)
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

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4+t-x^5)
            sage: L._places_finite(1)
            <generator object ...>
        """
        O = self.maximal_order()
        K = self.base_field()

        from sage.rings.integer import Integer
        degree = Integer(degree)

        for d in degree.divisors():
            for p in K._places_finite(degree=d):
                for prime,_,_ in O.decomposition(p.prime_ideal()):
                    place = prime.place()
                    if place.degree() == degree:
                        yield place

    def places_infinite(self, degree=1):
        """
        Return a list of the infinite places with ``degree``.

        INPUT:

        - ``degree`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4+t-x^5)
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

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: R.<t> = PolynomialRing(K)
            sage: L.<y> = K.extension(t^4+t-x^5)
            sage: L._places_infinite(1)
            <generator object ...>
        """
        Oinf = self.maximal_order_infinite()
        for prime,_,_ in Oinf.decomposition():
            place = prime.place()
            if place.degree() == degree:
                yield place

    def gaps(self):
        """
        Return the gaps of the function field.

        These are the gaps at the ordinary places, that is, places which are
        not Weierstrass places.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: L.gaps()
            [1, 2, 3]
        """
        return self._weierstrass_places()[1]

    def weierstrass_places(self):
        """
        Return all Weierstrass places of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: L.weierstrass_places()
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

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: len(L.weierstrass_places())  # indirect doctest
            10

        This method implements Algorithm 30 in [Hes2002b]_.
        """
        from sage.matrix.constructor import matrix

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

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: F.L_polynomial()
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

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: F.number_of_rational_places()
            4
            sage: [F.number_of_rational_places(r) for r in [1..10]]
            [4, 8, 4, 16, 44, 56, 116, 288, 508, 968]
        """
        from sage.rings.all import IntegerRing

        q = self.constant_field().order()
        L = self.L_polynomial()
        Lp = L.derivative()

        R = IntegerRing()[[L.parent().gen()]] # power series ring

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

        sage: from sage.rings.function_field.function_field import _singular_normal
        sage: R.<x,y> = QQ[]

        sage: f = (x^2-y^3) * x
        sage: _singular_normal(ideal(f))
        [[x, y], [1]]

        sage: f = (y^2-x)
        sage: _singular_normal(ideal(f))
        [[1]]
    """
    from sage.libs.singular.function import singular_function, lib
    lib('normal.lib')
    normal = singular_function('normal')
    execute = singular_function('execute')

    try:
        get_printlevel = singular_function('get_printlevel')
    except NameError:
        execute('proc get_printlevel {return (printlevel);}')
        get_printlevel = singular_function('get_printlevel')

    # It's fairly verbose unless printlevel is -1.
    saved_printlevel = get_printlevel()
    execute('printlevel=-1')
    nor = normal(ideal)
    execute('printlevel={}'.format(saved_printlevel))

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
        K = self.base_field() # rational function field
        n = self.degree()

        # Construct the defining polynomial of the function field as a
        # two-variate polynomial g in the ring k[y,x] where k is the constant
        # base field.
        S,(y,x) = PolynomialRing(k, names='y,x', order='lex').objgens()
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
                gflat += R(c) * R(m) # R(m) is a monomial in yy and xx

            k_poly = R(k.polynomial('zz'))

            # invoke Singular
            pols_in_R = normalize(R.ideal([k_poly, gflat]))

            # reconstruct polynomials in S
            h = R.hom([y,x,k.gen()],S)
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

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.equation_order()
            Order in Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(QQ); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.equation_order()
            Order in Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2
        """
        from .order import FunctionFieldOrder_basis
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

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
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

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.equation_order_infinite()
            Infinite order in Function field in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(QQ); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.equation_order_infinite()
            Infinite order in Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2
        """
        from .order import FunctionFieldOrderInfinite_basis
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
        sage: F.<y> = R.extension(y^2 - (x^2+1))
        sage: (y/x).divisor()
        - Place (x, y - 1)
         - Place (x, y + 1)
         + Place (x^2 + 1, y)

        sage: A.<z> = QQ[]
        sage: NF.<i> = NumberField(z^2+1)
        sage: R.<x> = FunctionField(NF)
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))

        sage: (x/y*x.differential()).divisor()
        -2*Place (1/x, 1/x*y - 1)
         - 2*Place (1/x, 1/x*y + 1)
         + Place (x, y - 1)
         + Place (x, y + 1)

        sage: (x/y).divisor()
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

            sage: K.<t> = FunctionField(CC); K
            Rational function field in t over Complex Field with 53 bits of precision
            sage: TestSuite(K).run()

            sage: FunctionField(QQ[I], 'alpha')
            Rational function field in alpha over Number Field in I with defining polynomial x^2 + 1 with I = 1*I

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

        from .place import FunctionFieldPlace_rational
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
        return "Rational function field in %s over %s"%(
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
            sage: I = S*[x^2 - y^2, y-t]
            sage: I.groebner_basis()
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
          constant of the underlying rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._to_constant_base_field(K(1))
            1
            sage: K._to_constant_base_field(K(x))
            Traceback (most recent call last):
            ...
            ValueError: only constants can be converted into the constant base field but x is not a constant

        TESTS:

        Verify that :trac:`21872` has been resolved::

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
        raise ValueError("only constants can be converted into the constant base field but %r is not a constant"%(f,))

    def _to_polynomial(self, f):
        """
        If ``f`` is integral, return it as a polynomial.

        INPUT:

        - ``f`` -- an element of this rational function field whose denominator is a constant.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._ring(x) # indirect doctest
            x
        """
        K = f.parent().constant_base_field()
        if f.denominator() in K:
            return f.numerator()/K(f.denominator())
        raise ValueError("only polynomials can be converted to the underlying polynomial ring")

    def _to_bivariate_polynomial(self, f):
        """
        Convert ``f`` from a univariate polynomial over the rational function
        field into a bivariate polynomial and a denominator.

        INPUT:

        - ``f`` -- univariate polynomial over the function field

        OUTPUT:

        - bivariate polynomial, denominator

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
        x,t = S.base_ring()['%s,%s'%(f.parent().variable_name(),self.variable_name())].gens()
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
            sage: f.factor()             # indirect doctest
            (1/t) * (X - t) * (X^2 - 1/t) * (X^2 + 1/t) * (X^2 + t*X + t^2)
            sage: f.factor().prod() == f
            True

        We do a factorization over a finite prime field::

            sage: R.<t> = FunctionField(GF(7))
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^4 - 1/t^2)*(X^3 - t^3)
            sage: f.factor()
            (1/t) * (X + 3*t) * (X + 5*t) * (X + 6*t) * (X^2 + 1/t) * (X^2 + 6/t)
            sage: f.factor().prod() == f
            True

        Factoring over a function field over a non-prime finite field::

            sage: k.<a> = GF(9)
            sage: R.<t> = FunctionField(k)
            sage: S.<X> = R[]
            sage: f = (1/t)*(X^3 - a*t^3)
            sage: f.factor()
            (1/t) * (X + (a + 2)*t)^3
            sage: f.factor().prod() == f
            True

        Factoring over a function field over a tower of finite fields::

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
        v = [(phi(P),e) for P, e in fac]
        unit = phi(fac.unit())/d
        w = []
        for a, e in v:
            c = a.leading_coefficient()
            a = a/c
            # undo any variable substitution that we introduced for the bivariate polynomial
            if old_variable_name != a.variable_name():
                a = a.change_variable_name(old_variable_name)
            unit *= (c**e)
            if a.is_unit():
                unit *= a**e
            else:
                w.append((a,e))
        from sage.structure.factorization import Factorization
        return Factorization(w, unit=unit)

    def extension(self, f, names=None):
        """
        Create an extension `L = K[y]/(f(y))` of the rational function field.

        INPUT:

        - ``f`` -- univariate polynomial over self

        - ``names`` -- string or length-1 tuple

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by y^3 + 1/t*y + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))
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

        - ``map`` -- (default ``True``), whether to return maps to and from the vector space

        OUTPUT:

        - a vector space `V` over base field

        - an isomorphism from `V` to the field

        - the inverse isomorphism from the field to `V`

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.free_module()
            (Vector space of dimension 1 over Rational function field in x over Rational Field, Isomorphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field, Isomorphism:
              From: Rational function field in x over Rational Field
              To:   Vector space of dimension 1 over Rational function field in x over Rational Field)

        TESTS::

            sage: K.free_module()
            (Vector space of dimension 1 over Rational function field in x over Rational Field, Isomorphism:
              From: Vector space of dimension 1 over Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field, Isomorphism:
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
        to_V   = MapFunctionFieldToVectorSpace(self, V)
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
        0, then an IndexError is raised.

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
          other ring.  If ``None``, try to use a coercion map.

        OUTPUT:

        - a map between function fields

        EXAMPLES:

        We make a map from a rational function field to itself::

            sage: K.<x> = FunctionField(GF(7))
            sage: K.hom( (x^4 + 2)/x)
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> (x^4 + 2)/x

        We construct a map from a rational function field into a
        non-rational extension field::

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
        if not isinstance(im_gens, (list,tuple)):
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
        from .order import FunctionFieldMaximalOrder_rational
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
        from .order import FunctionFieldMaximalOrderInfinite_rational
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
            sage: K.different()
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
        return 0

    def change_variable_name(self, name):
        r"""
        Return a field isomorphic to this field with variable ``name``.

        INPUT:

        - ``name`` -- a string or a tuple consisting of a single string, the
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
            id = Hom(self,self).identity()
            return self,id,id
        else:
            from .constructor import FunctionField
            ret = FunctionField(self.constant_base_field(), name)
            return ret, ret.hom(self.gen()), self.hom(ret.gen())

    @cached_method
    def derivation(self):
        r"""
        Return a derivation of the rational function field over the constant
        base field.

        The derivation maps the generator of the rational function field to 1.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: m = K.derivation(); m
            Derivation map:
              From: Rational function field in x over Finite Field of size 3
              To:   Rational function field in x over Finite Field of size 3
              Defn: x |--> 1
            sage: m(x)
            1

        TESTS::

            sage: L.<y> = FunctionField(K)
            sage: L.derivation()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented for non-perfect base fields
        """
        from .maps import FunctionFieldDerivation_rational
        if not self.constant_base_field().is_perfect():
            raise NotImplementedError("not implemented for non-perfect base fields")
        return FunctionFieldDerivation_rational(self, self.one())

    def residue_field(self, place, name=None):
        """
        Return the residue field of the place along with the maps from
        and to it.

        INPUT:

        - ``place`` -- place of the function field

        - ``name`` -- string; name of the generator of the residue field

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: p = F.places_finite(2)[0]
            sage: R, fr_R, to_R = F.residue_field(p)
            sage: R
            Finite Field in z2 of size 5^2
            sage: to_R(x) in R
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
            sage: d = F.higher_derivation()
            sage: [d(x^5,i) for i in range(10)]
            [x^5, 5*x^4, 10*x^3, 10*x^2, 5*x, 1, 0, 0, 0, 0]
            sage: [d(x^9,i) for i in range(10)]
            [x^9, 9*x^8, 36*x^7, 84*x^6, 126*x^5, 126*x^4, 84*x^3, 36*x^2, 9*x, 1]
        """
        from .maps import FunctionFieldHigherDerivation_char_zero
        return FunctionFieldHigherDerivation_char_zero(self)


class RationalFunctionField_global(RationalFunctionField):
    """
    Rational function field over finite fields.
    """
    _differentials_space = DifferentialsSpace_global

    def places(self, degree=1):
        """
        Return all places of the degree.

        INPUT:

        - ``degree`` -- (default: 1) a positive integer

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(5))
            sage: F.places()
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
            sage: F.places_finite()
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

        - ``degree`` -- a positive integer

        EXAMPLES::

            sage: F.<a> = GF(2)
            sage: K.<x> = FunctionField(F)
            sage: K.get_place(1)
            Place (x)
            sage: K.get_place(2)
            Place (x^2 + x + 1)
            sage: K.get_place(3)
            Place (x^3 + x + 1)
            sage: K.get_place(4)
            Place (x^4 + x + 1)
            sage: K.get_place(5)
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
            sage: d = F.higher_derivation()
            sage: [d(x^5,i) for i in range(10)]
            [x^5, 0, 0, 0, 0, 1, 0, 0, 0, 0]
            sage: [d(x^7,i) for i in range(10)]
            [x^7, 2*x^6, x^5, 0, 0, x^2, 2*x, 1, 0, 0]
        """
        from .maps import RationalFunctionFieldHigherDerivation_global
        return RationalFunctionFieldHigherDerivation_global(self)
