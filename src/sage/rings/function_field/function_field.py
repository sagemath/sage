r"""
Function Fields

A function field (of one variable) is a finitely generated field extension of
transcendence degree one. In Sage, a function field can be a rational function
field or a finite extension of a function field.

EXAMPLES:

We create a rational function field::

    sage: # needs sage.rings.finite_rings
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

    sage: # needs sage.rings.finite_rings sage.rings.function_field
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

    sage: # needs sage.rings.finite_rings sage.rings.function_field
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

    sage: N.<u> = FunctionField(K)                                                      # needs sage.rings.finite_rings

and inseparable extension function fields::

    sage: J.<x> = FunctionField(GF(5)); J
    Rational function field in x over Finite Field of size 5
    sage: T.<v> = J[]
    sage: O.<v> = J.extension(v^5 - x); O                                               # needs sage.rings.function_field
    Function field in v defined by v^5 + 4*x

Function fields over the rational field are supported::

    sage: # needs sage.rings.function_field
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

    sage: # needs sage.rings.function_field
    sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
    sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
    sage: O = F.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    2*Place (x, y, (1/(x^3 + x^2 + x))*y^2)
     + 2*Place (x^2 + x + 1, y, (1/(x^3 + x^2 + x))*y^2)

    sage: # needs sage.rings.function_field
    sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: O = L.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    - Place (x, x*y)
     + Place (x^2 + 1, x*y)

Function fields over the algebraic field are supported::

    sage: # needs sage.rings.function_field sage.rings.number_field
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
    sage: m(y)                                  # long time (4s)
    -2*s + (-4 - I)*s^2 + (-15 - 4*I)*s^3 + (-75 - 23*I)*s^4 + (-413 - 154*I)*s^5 + O(s^6)
    sage: m(y)^2 + m(y) + m(x) + 1/m(x)         # long time (8s)
    O(s^5)

TESTS::

    sage: TestSuite(J).run()
    sage: TestSuite(K).run(max_runs=256)        # long time (10s)                       # needs sage.rings.function_field sage.rings.number_field
    sage: TestSuite(L).run(max_runs=8)          # long time (25s)                       # needs sage.rings.function_field sage.rings.number_field

    sage: # needs sage.rings.finite_rings sage.rings.function_field
    sage: TestSuite(M).run(max_runs=8)                                  # long time (35s)
    sage: TestSuite(N).run(max_runs=8, skip='_test_derivation')         # long time (15s)
    sage: TestSuite(O).run()
    sage: TestSuite(R).run()
    sage: TestSuite(S).run()                                            # long time (4s)

Global function fields
----------------------

A global function field in Sage is an extension field of a rational function field
over a *finite* constant field by an irreducible separable polynomial over the
rational function field.

EXAMPLES:

A fundamental computation for a global or any function field is to get a basis
of its maximal order and maximal infinite order, and then do arithmetic with
ideals of those maximal orders::

    sage: # needs sage.rings.function_field
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

    sage: # needs sage.rings.function_field
    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
    sage: L.genus()
    3
    sage: L.weierstrass_places()                                                        # needs sage.modules
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
    sage: L.gaps()                                                                      # needs sage.modules
    [1, 2, 3]

The gap numbers for Weierstrass places are of course not ordinary::

    sage: # needs sage.rings.function_field
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

- Sebastian A. Spindler (2024-03-06): implemented Hilbert symbols for global
  function fields
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.rings.ring import Field
from sage.categories.homset import Hom
from sage.categories.function_fields import FunctionFields
from sage.structure.category_object import CategoryObject


def is_FunctionField(x):
    """
    Return ``True`` if ``x`` is a function field.

    EXAMPLES::

        sage: from sage.rings.function_field.function_field import is_FunctionField
        sage: is_FunctionField(QQ)
        doctest:warning...
        DeprecationWarning: The function is_FunctionField is deprecated; use '... in FunctionFields()' instead.
        See https://github.com/sagemath/sage/issues/38289 for details.
        False
        sage: is_FunctionField(FunctionField(QQ, 't'))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38289,
                "The function is_FunctionField is deprecated; "
                "use '... in FunctionFields()' instead.")
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
    _differentials_space = LazyImport('sage.rings.function_field.differential', 'DifferentialsSpace')

    def __init__(self, base_field, names, category=FunctionFields()):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: TestSuite(K).run()               # long time (3s)
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
           sage: L.<y> = K.extension(y^2 - x)                                           # needs sage.rings.function_field
           sage: L.some_elements()                                                      # needs sage.rings.function_field
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
            sage: K.<x> = FunctionField(QQbar)                                          # needs sage.rings.number_field
            sage: K.characteristic()
            0
            sage: K.<x> = FunctionField(GF(7))
            sage: K.characteristic()
            7
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)                                          # needs sage.rings.function_field
            sage: L.characteristic()                                                    # needs sage.rings.function_field
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
            sage: R.<t> = FunctionField(QQbar)                                          # needs sage.rings.number_field
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

        OUTPUT: a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)                                    # needs sage.rings.function_field
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1), 'z')                           # needs sage.rings.function_field
            Function field in z defined by z^3 + 1/t*z + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))                              # needs sage.rings.function_field
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

        OUTPUT: an order in the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)                              # needs sage.rings.function_field
            sage: O = L.order_with_basis([1, y, y^2]); O                                # needs sage.rings.function_field
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()                                                             # needs sage.rings.function_field
            (1, y, y^2)

        Note that 1 does not need to be an element of the basis, as long it is in the module spanned by it::

            sage: O = L.order_with_basis([1+y, y, y^2]); O                              # needs sage.rings.function_field
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()                                                             # needs sage.rings.function_field
            (y + 1, y, y^2)

        The following error is raised when the module spanned by the basis is not closed under multiplication::

            sage: O = L.order_with_basis([1, x^2 + x*y, (2/3)*y^2]); O                  # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: the module generated by basis (1, x*y + x^2, 2/3*y^2) must be closed under multiplication

        and this happens when the identity is not in the module spanned by the basis::

            sage: O = L.order_with_basis([x, x^2 + x*y, (2/3)*y^2])                     # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: the identity element must be in the module spanned by basis (x, x*y + x^2, 2/3*y^2)
        """
        from .order_basis import FunctionFieldOrder_basis
        return FunctionFieldOrder_basis(tuple([self(a) for a in basis]), check=check)

    def order(self, x, check=True):
        """
        Return the order generated by ``x`` over the base maximal order.

        INPUT:

        - ``x`` -- element or list of elements of the function field

        - ``check`` -- boolean; if ``True``, check that ``x`` really generates an order

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order(y); O                                                     # needs sage.modules
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()                                                             # needs sage.modules
            (1, y, y^2)

            sage: Z = K.order(x); Z                                                     # needs sage.rings.function_field
            Order in Rational function field in x over Rational Field
            sage: Z.basis()                                                             # needs sage.rings.function_field
            (1,)

        Orders with multiple generators are not yet supported::

            sage: Z = K.order([x, x^2]); Z                                              # needs sage.rings.function_field
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

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: O = L.order_infinite_with_basis([1, 1/x*y, 1/x^2*y^2]); O
            Infinite order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()
            (1, 1/x*y, 1/x^2*y^2)

        Note that 1 does not need to be an element of the basis, as long it is
        in the module spanned by it::

            sage: O = L.order_infinite_with_basis([1+1/x*y,1/x*y, 1/x^2*y^2]); O        # needs sage.rings.function_field
            Infinite order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: O.basis()                                                             # needs sage.rings.function_field
            (1/x*y + 1, 1/x*y, 1/x^2*y^2)

        The following error is raised when the module spanned by the basis is
        not closed under multiplication::

            sage: O = L.order_infinite_with_basis([1,y, 1/x^2*y^2]); O                  # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: the module generated by basis (1, y, 1/x^2*y^2) must be closed under multiplication

        and this happens when the identity is not in the module spanned by the
        basis::

            sage: O = L.order_infinite_with_basis([1/x,1/x*y, 1/x^2*y^2])               # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: the identity element must be in the module spanned by basis (1/x, 1/x*y, 1/x^2*y^2)
        """
        from .order_basis import FunctionFieldOrderInfinite_basis
        return FunctionFieldOrderInfinite_basis(tuple([self(g) for g in basis]), check=check)

    def order_infinite(self, x, check=True):
        """
        Return the order generated by ``x`` over the maximal infinite order.

        INPUT:

        - ``x`` -- element or a list of elements of the function field

        - ``check`` -- boolean; if ``True``, check that ``x`` really generates an order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)                              # needs sage.rings.function_field
            sage: L.order_infinite(y)           # not implemented                       # needs sage.rings.function_field

            sage: Z = K.order(x); Z                                                     # needs sage.modules
            Order in Rational function field in x over Rational Field
            sage: Z.basis()                                                             # needs sage.modules
            (1,)

        Orders with multiple generators, not yet supported::

            sage: Z = K.order_infinite([x, x^2]); Z
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

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)
            sage: L.equation_order()
            Order in Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: L._coerce_map_from_(L.equation_order())
            Conversion map:
              From: Order in Function field in y defined by y^3 + x^3 + 4*x + 1
              To:   Function field in y defined by y^3 + x^3 + 4*x + 1
            sage: L._coerce_map_from_(GF(7))

            sage: K.<x> = FunctionField(QQ)
            sage: L.<x> = FunctionField(GaussianIntegers().fraction_field())            # needs sage.rings.number_field
            sage: L.has_coerce_map_from(K)                                              # needs sage.rings.number_field
            True

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 1)                                          # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GaussianIntegers().fraction_field())            # needs sage.rings.number_field
            sage: R.<y> = K[]
            sage: M.<y> = K.extension(y^3 + 1)                                          # needs sage.rings.function_field
            sage: M.has_coerce_map_from(L)      # not tested (the constant field including into a function field is not yet known to be injective), needs sage.rings.function_field
            True

            sage: K.<x> = FunctionField(QQ)
            sage: R.<I> = K[]
            sage: L.<I> = K.extension(I^2 + 1)                                          # needs sage.rings.function_field
            sage: M.<x> = FunctionField(GaussianIntegers().fraction_field())            # needs sage.rings.number_field
            sage: M.has_coerce_map_from(L)                                              # needs sage.rings.function_field sage.rings.number_field
            True

        Check that :issue:`31072` is fixed::

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
            sage: TestSuite(K).run()    # indirect doctest, long time (3s)
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        K = self.constant_base_field().some_elements()
        try:
            d = self.derivation()
        except ImportError:
            return
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
            sage: L.<y> = K.extension(y^3 + x^3 + 4*x + 1)                              # needs sage.rings.function_field
            sage: K(L(x))  # indirect doctest                                           # needs sage.rings.function_field
            x
        """
        try:
            from .function_field_polymod import FunctionField_polymod
        except ImportError:
            FunctionField_polymod = ()
        if isinstance(R, FunctionField_polymod):
            base_conversion = self.convert_map_from(R.base_field())
            if base_conversion is not None:
                from sage.categories.morphism import SetMorphism
                return base_conversion * SetMorphism(R.Hom(R.base_field()), R._to_base_field)

    def _intermediate_fields(self, base):
        """
        Return the fields which lie in between base and the function field in
        the tower of function fields.

        INPUT:

        - ``base`` -- function field, either this field or a field from which
          this field has been created as an extension

        OUTPUT: list of fields; the first entry is this field, the last entry is ``base``

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K._intermediate_fields(K)
            [Rational function field in x over Rational Field]

            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)                                          # needs sage.rings.function_field
            sage: L._intermediate_fields(K)                                             # needs sage.rings.function_field
            [Function field in y defined by y^2 - x,
             Rational function field in x over Rational Field]

            sage: # needs sage.rings.function_field
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M._intermediate_fields(L)
            [Function field in z defined by z^2 - y,
             Function field in y defined by y^2 - x]
            sage: M._intermediate_fields(K)
            [Function field in z defined by z^2 - y,
             Function field in y defined by y^2 - x,
            Rational function field in x over Rational Field]

        TESTS::

            sage: K._intermediate_fields(M)                                             # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: field has not been constructed as a finite extension of base
            sage: K._intermediate_fields(QQ)
            Traceback (most recent call last):
            ...
            TypeError: base must be a function field
        """
        if base not in FunctionFields():
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
            sage: L.<y> = K.extension(y^2 - x)                                          # needs sage.rings.function_field
            sage: L.rational_function_field()                                           # needs sage.rings.function_field
            Rational function field in x over Rational Field

            sage: R.<z> = L[]                                                           # needs sage.rings.function_field
            sage: M.<z> = L.extension(z^2 - y)                                          # needs sage.rings.function_field
            sage: M.rational_function_field()                                           # needs sage.rings.function_field
            Rational function field in x over Rational Field
        """
        from .function_field_rational import RationalFunctionField

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
            sage: v = K.valuation(1); v                                                 # needs sage.rings.function_field
            (x - 1)-adic valuation
            sage: v(x)                                                                  # needs sage.rings.function_field
            0
            sage: v(x - 1)                                                              # needs sage.rings.function_field
            1

        A place can also be specified with an irreducible polynomial::

            sage: v = K.valuation(x - 1); v                                             # needs sage.rings.function_field
            (x - 1)-adic valuation

        Similarly, for a finite non-rational place::

            sage: v = K.valuation(x^2 + 1); v                                           # needs sage.rings.function_field
            (x^2 + 1)-adic valuation
            sage: v(x^2 + 1)                                                            # needs sage.rings.function_field
            1
            sage: v(x)                                                                  # needs sage.rings.function_field
            0

        Or for the infinite place::

            sage: v = K.valuation(1/x); v                                               # needs sage.rings.function_field
            Valuation at the infinite place
            sage: v(x)                                                                  # needs sage.rings.function_field
            -1

        Instead of specifying a generator of a place, we can define a valuation on a
        rational function field by giving a discrete valuation on the underlying
        polynomial ring::

            sage: # needs sage.rings.function_field
            sage: R.<x> = QQ[]
            sage: u = valuations.GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: w = u.augmentation(x - 1, 1)
            sage: v = K.valuation(w); v
            (x - 1)-adic valuation

        Note that this allows us to specify valuations which do not correspond to a
        place of the function field::

            sage: w = valuations.GaussValuation(R, QQ.valuation(2))                     # needs sage.rings.function_field
            sage: v = K.valuation(w); v                                                 # needs sage.rings.function_field
            2-adic valuation

        The same is possible for valuations with `v(1/x) > 0` by passing in an
        extra pair of parameters, an isomorphism between this function field and an
        isomorphic function field. That way you can, for example, indicate that the
        valuation is to be understood as a valuation on `K[1/x]`, i.e., after
        applying the substitution `x \mapsto 1/x` (here, the inverse map is also `x
        \mapsto 1/x`)::

            sage: # needs sage.rings.function_field
            sage: w = valuations.GaussValuation(R, QQ.valuation(2)).augmentation(x, 1)
            sage: w = K.valuation(w)
            sage: v = K.valuation((w, K.hom([~K.gen()]), K.hom([~K.gen()]))); v
            Valuation on rational function field
            induced by [ Gauss valuation induced by 2-adic valuation, v(x) = 1 ]
            (in Rational function field in x over Rational Field after x |--> 1/x)

        Note that classical valuations at finite places or the infinite place are
        always normalized such that the uniformizing element has valuation 1::

            sage: # needs sage.rings.function_field
            sage: K.<t> = FunctionField(GF(3))
            sage: M.<x> = FunctionField(K)
            sage: v = M.valuation(x^3 - t)
            sage: v(x^3 - t)
            1

        However, if such a valuation comes out of a base change of the ground
        field, this is not the case anymore. In the example below, the unique
        extension of ``v`` to ``L`` still has valuation 1 on `x^3 - t` but it has
        valuation ``1/3`` on its uniformizing element  `x - w`::

            sage: # needs sage.rings.function_field
            sage: R.<w> = K[]
            sage: L.<w> = K.extension(w^3 - t)
            sage: N.<x> = FunctionField(L)
            sage: w = v.extension(N)  # missing factorization, :issue:`16572`
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: w(x^3 - t)                    # not tested
            1
            sage: w(x - w)                      # not tested
            1/3

        There are several ways to create valuations on extensions of rational
        function fields::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x); L                                       # needs sage.rings.function_field
            Function field in y defined by y^2 - x

        A place that has a unique extension can just be defined downstairs::

            sage: v = L.valuation(x); v                                                 # needs sage.rings.function_field
            (x)-adic valuation
        """
        from sage.rings.function_field.valuation import FunctionFieldValuation
        return FunctionFieldValuation(self, prime)

    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.space_of_differentials()                                            # needs sage.modules
            Space of differentials of Rational function field in t over Rational Field

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                        # needs sage.rings.function_field
            sage: L.space_of_differentials()                                            # needs sage.rings.function_field
            Space of differentials of Function field in y
             defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
        """
        return self._differentials_space(self)

    def space_of_holomorphic_differentials(self):
        """
        Return the space of holomorphic differentials of this function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.space_of_holomorphic_differentials()                                # needs sage.libs.pari sage.modules
            (Vector space of dimension 0 over Rational Field,
             Linear map:
               From: Vector space of dimension 0 over Rational Field
               To:   Space of differentials of Rational function field in t over Rational Field,
             Section of linear map:
               From: Space of differentials of Rational function field in t over Rational Field
               To:   Vector space of dimension 0 over Rational Field)

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                        # needs sage.rings.function_field
            sage: L.space_of_holomorphic_differentials()                                # needs sage.rings.function_field
            (Vector space of dimension 4 over Finite Field of size 5,
             Linear map:
               From: Vector space of dimension 4 over Finite Field of size 5
               To:   Space of differentials of Function field in y
                      defined by y^3 + (4*x^3 + 1)/(x^3 + 3),
             Section of linear map:
               From: Space of differentials of Function field in y
                      defined by y^3 + (4*x^3 + 1)/(x^3 + 3)
               To:   Vector space of dimension 4 over Finite Field of size 5)
        """
        return self.divisor_group().zero().differential_space()

    space_of_differentials_of_first_kind = space_of_holomorphic_differentials

    def basis_of_holomorphic_differentials(self):
        """
        Return a basis of the space of holomorphic differentials of this function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.basis_of_holomorphic_differentials()                                # needs sage.libs.pari sage.modules
            []

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                        # needs sage.rings.function_field
            sage: L.basis_of_holomorphic_differentials()                                # needs sage.rings.function_field
            [((x/(x^3 + 4))*y) d(x),
             ((1/(x^3 + 4))*y) d(x),
             ((x/(x^3 + 4))*y^2) d(x),
             ((1/(x^3 + 4))*y^2) d(x)]
        """
        return self.divisor_group().zero().basis_differential_space()

    basis_of_differentials_of_first_kind = basis_of_holomorphic_differentials

    def divisor_group(self):
        """
        Return the group of divisors attached to the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: K.divisor_group()                                                     # needs sage.modules
            Divisor group of Rational function field in t over Rational Field

            sage: _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (t^3 - 1)/(t^3 - 2))                        # needs sage.rings.function_field
            sage: L.divisor_group()                                                     # needs sage.rings.function_field
            Divisor group of Function field in y defined by y^3 + (-t^3 + 1)/(t^3 - 2)

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))                        # needs sage.rings.function_field
            sage: L.divisor_group()                                                     # needs sage.rings.function_field
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
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.function_field
            sage: L.place_set()                                                         # needs sage.rings.function_field
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

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p); m
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(x, 10)
            s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + s^9 + s^10 + O(s^12)
            sage: m(y, 10)
            s^-1 + 1 + s^3 + s^5 + s^7 + O(s^9)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p); m
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(x, 10)
            s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + s^9 + s^10 + O(s^12)
            sage: m(y, 10)
            s^-1 + 1 + s^3 + s^5 + s^7 + O(s^9)

            sage: K.<x> = FunctionField(GF(2))
            sage: p = K.places_finite()[0]; p                                           # needs sage.libs.pari
            Place (x)
            sage: m = K.completion(p); m                                                # needs sage.rings.function_field
            Completion map:
              From: Rational function field in x over Finite Field of size 2
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(1/(x+1))                                                            # needs sage.rings.function_field
            1 + s + s^2 + s^3 + s^4 + s^5 + s^6 + s^7 + s^8 + s^9 + s^10 + s^11 + s^12
            + s^13 + s^14 + s^15 + s^16 + s^17 + s^18 + s^19 + O(s^20)

            sage: p = K.place_infinite(); p
            Place (1/x)
            sage: m = K.completion(p); m                                                # needs sage.rings.function_field
            Completion map:
              From: Rational function field in x over Finite Field of size 2
              To:   Laurent Series Ring in s over Finite Field of size 2
            sage: m(x)                                                                  # needs sage.rings.function_field
            s^-1 + O(s^19)

            sage: m = K.completion(p, prec=infinity); m                                 # needs sage.rings.function_field
            Completion map:
              From: Rational function field in x over Finite Field of size 2
              To:   Lazy Laurent Series Ring in s over Finite Field of size 2
            sage: f = m(x); f                                                           # needs sage.rings.function_field
            s^-1 + ...
            sage: f.coefficient(100)                                                    # needs sage.rings.function_field
            0

            sage: # needs sage.rings.function_field
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

            sage: # needs sage.rings.function_field
            sage: decomp2 = O.decomposition(K.maximal_order().ideal(x^2 + 1))
            sage: pls2 = decomp2[0][0].place()
            sage: m = L.completion(pls2); m
            Completion map:
              From: Function field in y defined by y^2 - x
              To:   Laurent Series Ring in s over
                     Number Field in a with defining polynomial x^4 + 2*x^2 + 4*x + 2
            sage: xe = m(x)
            sage: ye = m(y)
            sage: ye^2 - xe == 0
            True
        """
        from .maps import FunctionFieldCompletion
        return FunctionFieldCompletion(self, place, name=name, prec=prec, gen_name=gen_name)

    def hilbert_symbol(self, a, b, P):
        r"""
        Return the Hilbert symbol `(a,b)_{F_P}` for the local field `F_P`.

        The local field `F_P` is the completion of this function field `F`
        at the place `P`.

        INPUT:

        - ``a``, ``b`` -- elements of this function field

        - ``P`` -- a place of this function field

        The Hilbert symbol `(a,b)_{F_P}` is `0` if `a` or `b` is zero.
        Otherwise it takes the value `1` if  the quaternion algebra
        defined by `(a,b)` over `F_P` is split, and `-1` if said
        algebra is a division ring.

        ALGORITHM:

        For the valuation `\nu = \nu_P` of `F`, we compute the valuations
        `\nu(a)` and `\nu(b)` as well as elements `a_0` and `b_0` of the
        residue field such that for a uniformizer `\pi` at `P`,
        `a\pi^{-\nu(a))}` respectively `b\pi^{-\nu(b)}` has the residue class
        `a_0` respectively `b_0` modulo `\pi`. Then the Hilbert symbol is
        computed by formula 12.4.10 in [Voi2021]_.

        Currently only implemented for global function fields.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(17))
            sage: P = K.places()[0]; P
            Place (1/x)
            sage: a = (5*x + 6)/(x + 15)
            sage: b = 7/x
            sage: K.hilbert_symbol(a, b, P)
            -1

            sage: Q = K.places()[7]; Q
            Place (x + 6)
            sage: c = 15*x + 12
            sage: d = 16/(x + 13)
            sage: K.hilbert_symbol(c, d, Q)
            1

        Check that the Hilbert symbol is symmetric and bimultiplicative::

            sage: K.<x> = FunctionField(GF(5)); R.<T> = PolynomialRing(K)
            sage: f = ((x^2 + 2*x + 2)*T^5 + (4*x^2 + 2*x + 3)*T^4 + 3*T^3 + 4*T^2
            ....:     + (2/(x^2 + 4*x + 1))*T + 3*x^2 + 2*x + 4)
            sage: L.<y> = K.extension(f)
            sage: a = L.random_element()
            sage: b = L.random_element()
            sage: c = L.random_element()
            sage: P = L.places_above(K.places()[0])[1]
            sage: Q = L.places_above(K.places()[1])[0]

            sage: hP_a_c = L.hilbert_symbol(a, c, P)
            sage: hP_a_c == L.hilbert_symbol(c, a, P)
            True
            sage: L.hilbert_symbol(a, b, P) * hP_a_c == L.hilbert_symbol(a, b*c, P)
            True
            sage: hP_a_c * L.hilbert_symbol(b, c, P) == L.hilbert_symbol(a*b, c, P)
            True

            sage: hQ_a_c = L.hilbert_symbol(a, c, Q)
            sage: hQ_a_c == L.hilbert_symbol(c, a, Q)
            True
            sage: L.hilbert_symbol(a, b, Q) * hQ_a_c == L.hilbert_symbol(a, b*c, Q)
            True
            sage: hQ_a_c * L.hilbert_symbol(b, c, Q) == L.hilbert_symbol(a*b, c, Q)
            True
        """
        if not self.is_global():
            raise NotImplementedError('only supported for global function fields')

        if self.characteristic() == 2:
            raise ValueError('Hilbert symbol is only defined for'
                            ' odd characteristic function fields')

        if not (a in self and b in self):
            raise ValueError('a and b must be elements of the function field')

        if a.is_zero() or b.is_zero():
            return 0

        # Compute the completion map to precision 1 for computation of the
        # valuations v(a), v(b) as well as the elements a0, b0
        try:
            sigma = self.completion(P, prec=1, gen_name='i')
        except AttributeError:
            raise ValueError('P must be a place of the function field F')

        # Apply the completion map to a to get v(a) and a0
        ser_a = sigma(a)
        v_a = ser_a.valuation()
        a0 = ser_a.coefficients()[0]

        # Apply the completion map to b to get v(b) and b0
        ser_b = sigma(b)
        v_b = ser_b.valuation()
        b0 = ser_b.coefficients()[0]

        # Get the residue field of the completion together with the necessary exponent
        k = sigma.codomain().base_ring()
        e = (k.order() - 1) // 2

        # Use Euler's criterion to compute the powers of Legendre symbols
        a_rd_pw = a0**(v_b * e)
        b_rd_pw = b0**(v_a * e)

        # Finally, put the result together and transform it into the correct output
        res = k(-1)**(v_a * v_b * e) * a_rd_pw * b_rd_pw

        from sage.rings.integer import Integer
        return Integer(1) if res.is_one() else Integer(-1)

    def extension_constant_field(self, k):
        """
        Return the constant field extension with constant field `k`.

        INPUT:

        - ``k`` -- an extension field of the constant field of this function field

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: E = F.extension_constant_field(GF(2^4))
            sage: E
            Function field in y defined by y^2 + y + (x^2 + 1)/x over its base
            sage: E.constant_base_field()
            Finite Field in z4 of size 2^4
        """
        from .extensions import ConstantFieldExtension
        return ConstantFieldExtension(self, k)

    @cached_method
    def jacobian(self, model=None, base_div=None, **kwds):
        """
        Return the Jacobian of the function field.

        INPUT:

        - ``model`` -- (default: ``'hess'``) model to use for arithmetic

        - ``base_div`` -- an effective divisor

        The degree of the base divisor should satisfy certain degree condition
        corresponding to the model used. The following table lists these
        conditions. Let `g` be the genus of the function field.

        - ``hess``: ideal-based arithmetic; requires base divisor of degree `g`

        - ``km_large``: Khuri-Makdisi's large model; requires base divisor of
          degree at least `2g + 1`

        - ``km_medium``: Khuri-Makdisi's medium model; requires base divisor of
          degree at least `2g + 1`

        - ``km_small``: Khuri-Makdisi's small model requires base divisor of
          degree at least `g + 1`

        We assume the function field has a rational place. If a base divisor is
        not given, one is constructed using an arbitrary rational place.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(5), 2)
            sage: C = Curve(y^2*(x^3 - 1) - (x^3 - 2))
            sage: F = C.function_field()
            sage: F.jacobian()
            Jacobian of Function field in y defined by (x^3 + 4)*y^2 + 4*x^3 + 2 (Hess model)

        TESTS:

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve(y^2 - x^3 - 1, A).projective_closure()
            sage: C.jacobian(model='hess')
            Traceback (most recent call last):
            ...
            ValueError: failed to obtain a rational place; provide a base divisor
        """
        from .place import FunctionFieldPlace

        if model is None:
            model = 'hess'

        if base_div is None:
            try:
                base_place = self.get_place(1)
            except AttributeError:
                raise ValueError('failed to obtain a rational place; provide a base divisor')
            if base_place is None:
                raise ValueError('the function field has no rational place')
            # appropriate base divisor is constructed below.
        else:
            if isinstance(base_div, FunctionFieldPlace):
                base_div = base_div.divisor()

        g = self.genus()
        curve = kwds.get('curve')

        if model.startswith('km'):
            from .jacobian_khuri_makdisi import Jacobian
            if model == 'km' or model.endswith('large'):
                if base_div is None:
                    base_div = (2*g + 1) * base_place
                if not base_div.degree() >= 2*g + 1:
                    raise ValueError("Khuri-Makdisi large model requires base divisor of degree "
                                     "at least 2*g + 1 for genus g")
                return Jacobian(self, base_div, model='large', curve=curve)
            elif model.endswith('medium'):
                if base_div is None:
                    base_div = (2*g + 1) * base_place
                if not base_div.degree() >= 2*g + 1:
                    raise ValueError("Khuri-Makdisi medium model requires base divisor of degree "
                                     "at least 2*g + 1 for genus g")
                return Jacobian(self, base_div, model='medium', curve=curve)
            elif model.endswith('small'):
                if base_div is None:
                    base_div = (g + 1) * base_place
                if not base_div.degree() >= g + 1:
                    raise ValueError("Khuri-Makdisi small model requires base divisor of degree "
                                     "at least g + 1 for genus g")
                return Jacobian(self, base_div, model='small', curve=curve)
        elif model == 'hess':
            from .jacobian_hess import Jacobian
            if base_div is None:
                base_div = g * base_place
            if base_div.degree() != g:
                raise ValueError("Hess model requires base divisor of degree g for genus g")
            return Jacobian(self, base_div, curve=curve)

        raise ValueError("unknown model")
