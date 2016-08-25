r"""
Orders in Function Fields

AUTHORS:

- William Stein (2010): initial version

- Maarten Derickx (2011-09-14): fixed ideal_with_gens_over_base() for rational function fields

- Julian Rueth (2011-09-14): added check in _element_constructor_

EXAMPLES:

Maximal orders in rational function fields::

    sage: K.<x> = FunctionField(QQ)
    sage: O = K.maximal_order()
    sage: I = O.ideal(1/x); I
    Ideal (1/x) of Maximal order in Rational function field in x over Rational Field
    sage: 1/x in O
    False

Equation orders in extensions of rational function fields::

    sage: K.<x> = FunctionField(GF(3)); R.<y> = K[]
    sage: L.<y> = K.extension(y^3-y-x)
    sage: O = L.equation_order()
    sage: 1/y in O
    False
    sage: x/y in O
    True
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#       Copyright (C) 2011 Julian Rueth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import IntegralDomain, PrincipalIdealDomain
from sage.rings.ideal import is_Ideal

class FunctionFieldOrder(IntegralDomain):
    """
    Base class for orders in function fields.
    """
    def __init__(self, fraction_field):
        """
        INPUT:

            - ``fraction_field`` -- the function field in which this is an order.

        EXAMPLES::

            sage: R = FunctionField(QQ,'y').maximal_order()
            sage: isinstance(R, sage.rings.function_field.function_field_order.FunctionFieldOrder)
            True
        """
        IntegralDomain.__init__(self, self)
        self._fraction_field = fraction_field

    def _repr_(self):
        """
        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()._repr_()
            'Maximal order in Rational function field in y over Rational Field'
        """
        return "Order in %s"%self.fraction_field()

    def is_finite(self):
        """
        Returns False since orders are never finite.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_finite()
            False
        """
        return False

    def is_field(self, proof=True):
        """
        Returns False since orders are never fields.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_field()
            False
        """
        return False

    def is_noetherian(self):
        """
        Returns True since orders in function fields are noetherian.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_noetherian()
            True
        """
        return True

    def fraction_field(self):
        """
        Returns the function field in which this is an order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().fraction_field()
            Rational function field in y over Rational Field
        """
        return self._fraction_field

    function_field = fraction_field

    def ideal_with_gens_over_base(self, gens):
        """
        Returns the fractional ideal with basis ``gens`` over the
        maximal order of the base field. That this is really an ideal
        is not checked.

        INPUT:

            - ``gens`` -- list of elements that are a basis for the
              ideal over the maximal order of the base field

        EXAMPLES:

        We construct an ideal in a rational function field::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal_with_gens_over_base([y]); I
            Ideal (y) of Maximal order in Rational function field in y over Rational Field
            sage: I*I
            Ideal (y^2) of Maximal order in Rational function field in y over Rational Field

        We construct some ideals in a nontrivial function field::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order(); O
            Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I = O.ideal_with_gens_over_base([1, y]);  I
            Ideal (1, y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I.module()
            Free module of degree 2 and rank 2 over Maximal order in Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [1 0]
            [0 1]

        There is no check if the resulting object is really an ideal::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([y]); I
            Ideal (y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y in I
            True
            sage: y^2 in I
            False
        """
        from .function_field_ideal import ideal_with_gens_over_base
        return ideal_with_gens_over_base(self, [self(a) for a in gens])

    def ideal(self, *gens):
        """
        Returns the fractional ideal generated by the elements in ``gens``.

        INPUT:

            - ``gens`` -- a list of generators or an ideal in a ring which
                          coerces to this order.

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O.ideal(y)
            Ideal (y) of Maximal order in Rational function field in y over Rational Field
            sage: O.ideal([y,1/y]) == O.ideal(y,1/y) # multiple generators may be given as a list
            True

        A fractional ideal of a nontrivial extension::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2-4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.equation_order()
            sage: S.ideal(1/y)
            Ideal (1, (6/(x^3 + 1))*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 = S.ideal(x^2-4); I2
            Ideal (x^2 + 3, (x^2 + 3)*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 == S.ideal(I)
            True
        """
        if len(gens) == 1:
            gens = gens[0]
            if not isinstance(gens, (list, tuple)):
                if is_Ideal(gens):
                    gens = gens.gens()
                else:
                    gens = [gens]
        from .function_field_ideal import ideal_with_gens
        return ideal_with_gens(self, gens)

class FunctionFieldOrder_basis(FunctionFieldOrder):
    """
    An order given by a basis over the maximal order of the base
    field.
    """
    def __init__(self, basis, check=True):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order(); O
            Order in Function field in y defined by y^4 + x*y + 4*x + 1
            sage: type(O)
            <class 'sage.rings.function_field.function_field_order.FunctionFieldOrder_basis_with_category'>

        The basis only defines an order if the module it generates is closed under multiplication
         and contains the identity element (only checked when ``check`` is True)::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x));
            sage: y.is_integral()
            False
            sage: L.order(y)
            Traceback (most recent call last):
            ...
            ValueError: The module generated by basis [1, y, y^2, y^3, y^4] must be closed under multiplication

        The basis also has to be linearly independent and of the same rank as the degree of the function field of its elements (only checked when ``check`` is True)::

            sage: L.order(L(x))
            Traceback (most recent call last):
            ...
            ValueError: Basis [1, x, x^2, x^3, x^4] is not linearly independent
            sage: sage.rings.function_field.function_field_order.FunctionFieldOrder_basis([y,y,y^3,y^4,y^5])
            Traceback (most recent call last):
            ...
            ValueError: Basis [y, y, y^3, y^4, 2*x*y + (x^4 + 1)/x] is not linearly independent
        """
        if len(basis) == 0:
            raise ValueError("basis must have positive length")

        fraction_field = basis[0].parent()
        if len(basis) != fraction_field.degree():
            raise ValueError("length of basis must equal degree of field")

        FunctionFieldOrder.__init__(self, fraction_field)

        self._basis = tuple(basis)
        V, fr, to = fraction_field.vector_space()
        R = fraction_field.base_field().maximal_order()
        self._module = V.span([to(b) for b in basis], base_ring=R)
        self._ring = fraction_field.polynomial_ring()
        self._populate_coercion_lists_(coerce_list=[self._ring])
        if check:
            if self._module.rank() != fraction_field.degree():
                raise ValueError("Basis %s is not linearly independent"%(basis))
            if not to(fraction_field(1)) in self._module:
                raise ValueError("The identity element must be in the module spanned by basis %s"%(basis))
            if not all(to(a*b) in self._module for a in basis for b in basis):
                raise ValueError("The module generated by basis %s must be closed under multiplication"%(basis))
        IntegralDomain.__init__(self, self, names = fraction_field.variable_names(), normalize = False)

    def _element_constructor_(self, f, check=True):
        """
        Make ``f`` into an element of this order.

        INPUT:

        - ``f`` -- the element
        - ``check`` -- check if the element is in the order

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.maximal_order()._element_constructor_(x)
            x
        """
        fraction_field=self.fraction_field()

        if f.parent() is fraction_field:
            f = f.element()
        f = self._ring(f)
        if check:
            V, fr, to = fraction_field.vector_space()
            f_vector = to(fraction_field(f))
            if not f_vector in self._module:
                raise TypeError("%r is not an element of %r"%(f_vector,self))
        return fraction_field._element_class(self, f)

    def fraction_field(self):
        """
        Returns the function field in which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.fraction_field()
            Function field in y defined by y^4 + x*y + 4*x + 1
        """
        return self._fraction_field

    def polynomial(self):
        """
        Returns the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._fraction_field.polynomial()

    def basis(self):
        """
        Returns a basis of self over the maximal order of the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)
        """
        return self._basis

    def free_module(self):
        """
        Returns the free module formed by the basis over the maximal order of the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.free_module()
            Free module of degree 4 and rank 4 over Maximal order in Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        return self._module

class FunctionFieldOrder_rational(PrincipalIdealDomain, FunctionFieldOrder):
    """
    The maximal order in a rational function field.
    """
    def __init__(self, function_field):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19)); K
            Rational function field in t over Finite Field of size 19
            sage: R = K.maximal_order(); R
            Maximal order in Rational function field in t over Finite Field of size 19
            sage: type(R)
            <class 'sage.rings.function_field.function_field_order.FunctionFieldOrder_rational_with_category'>
        """
        FunctionFieldOrder.__init__(self, function_field)
        PrincipalIdealDomain.__init__(self, self, names = function_field.variable_names(), normalize = False)
        self._ring = function_field._ring
        self._populate_coercion_lists_(coerce_list=[self._ring])
        self._gen = self(self._ring.gen())
        self._basis = (self(1),)

    def basis(self):
        """
        Returns the basis (=1) for this order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order()
            sage: O.basis()
            (1,)
            sage: parent(O.basis()[0])
            Maximal order in Rational function field in t over Finite Field of size 19
        """
        return self._basis

    def ideal(self, *gens):
        """
        Returns the fractional ideal generated by ``gens``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O.ideal(x)
            Ideal (x) of Maximal order in Rational function field in x over Rational Field
            sage: O.ideal([x,1/x]) == O.ideal(x,1/x) # multiple generators may be given as a list
            True
            sage: O.ideal(x^3+1,x^3+6)
            Ideal (1) of Maximal order in Rational function field in x over Rational Field
            sage: I = O.ideal((x^2+1)*(x^3+1),(x^3+6)*(x^2+1)); I
            Ideal (x^2 + 1) of Maximal order in Rational function field in x over Rational Field
            sage: O.ideal(I)
            Ideal (x^2 + 1) of Maximal order in Rational function field in x over Rational Field
        """
        if len(gens) == 1:
            gens = gens[0]
            if not isinstance(gens, (list, tuple)):
                if is_Ideal(gens):
                    gens = gens.gens()
                else:
                    gens = (gens,)
        from .function_field_ideal import ideal_with_gens
        return ideal_with_gens(self, gens)

    def _repr_(self):
        """
        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()._repr_()
            'Maximal order in Rational function field in y over Rational Field'
        """
        return "Maximal order in %s"%self.fraction_field()

    def gen(self, n=0):
        """
        Returns the ``n``-th generator of self. Since there is only one generator ``n`` must be 0.

        EXAMPLES::

            sage: O = FunctionField(QQ,'y').maximal_order()
            sage: O.gen()
            y
            sage: O.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0: raise IndexError("Only one generator.")
        return self._gen

    def ngens(self):
        """
        Returns 1, the number of generators of self.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().ngens()
            1
        """
        return 1

    def _element_constructor_(self, f):
        """
        Make ``f`` into an element of this order.

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O._element_constructor_(y)
            y
            sage: O._element_constructor_(1/y)
            Traceback (most recent call last):
            ...
            TypeError: 1/y is not an element of Maximal order in Rational function field in y over Rational Field
        """
        if f.parent() is self.fraction_field():
            if not f.denominator() in self.fraction_field().constant_base_field():
                raise TypeError("%r is not an element of %r"%(f,self))
            f = f.element()
        from .function_field_element import FunctionFieldElement_rational
        return FunctionFieldElement_rational(self, self._ring(f))
