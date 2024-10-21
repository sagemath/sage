# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules
"""
Differentials of function fields

Sage provides arithmetic with differentials of function fields.

EXAMPLES:

The module of differentials on a function field forms an one-dimensional vector space over
the function field::

    sage: # needs sage.rings.finite_rings sage.rings.function_field
    sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
    sage: f = x + y
    sage: g = 1 / y
    sage: df = f.differential()
    sage: dg = g.differential()
    sage: dfdg = f.derivative() / g.derivative()
    sage: df == dfdg * dg
    True
    sage: df
    (x*y^2 + 1/x*y + 1) d(x)
    sage: df.parent()
    Space of differentials of Function field in y defined by y^3 + x^3*y + x

We can compute a canonical divisor::

    sage: # needs sage.rings.finite_rings sage.rings.function_field
    sage: k = df.divisor()
    sage: k.degree()
    4
    sage: k.degree() == 2 * L.genus() - 2
    True

Exact differentials vanish and logarithmic differentials are stable under the
Cartier operation::

    sage: # needs sage.rings.finite_rings sage.rings.function_field
    sage: df.cartier()
    0
    sage: w = 1/f * df
    sage: w.cartier() == w
    True

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version
"""

# ****************************************************************************
#       Copyright (C) 2016-2019 Kwankyu Lee <ekwankyu@gmail.com>
#                     2019      Brent Baccala
#                     2019      Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.latex import latex

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp
from sage.sets.family import Family

from sage.categories.modules import Modules
from sage.categories.morphism import Morphism


class FunctionFieldDifferential(ModuleElement):
    """
    Base class for differentials on function fields.

    INPUT:

    - ``f`` -- element of the function field

    - ``t`` -- element of the function field; if `t` is not specified, the generator
      of the base differential is assumed

    EXAMPLES::

        sage: F.<x> = FunctionField(QQ)
        sage: f = x/(x^2 + x + 1)
        sage: f.differential()
        ((-x^2 + 1)/(x^4 + 2*x^3 + 3*x^2 + 2*x + 1)) d(x)

    ::

        sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                      # needs sage.rings.function_field
        sage: L(x).differential()                                                       # needs sage.rings.function_field
        d(x)
        sage: y.differential()                                                          # needs sage.rings.function_field
        ((21/4*x/(x^7 + 27/4))*y^2 + ((3/2*x^7 + 9/4)/(x^8 + 27/4*x))*y + 7/2*x^4/(x^7 + 27/4)) d(x)
    """
    def __init__(self, parent, f, t=None):
        """
        Initialize the differential `fdt`.

        TESTS::

            sage: F.<x> = FunctionField(GF(7))
            sage: f = x/(x^2 + x + 1)
            sage: w = f.differential()
            sage: TestSuite(w).run()
        """
        ModuleElement.__init__(self, parent)

        if t is not None:
            f *= parent._derivation(t) * parent._gen_derivative_inv

        self._f = f

    def _repr_(self):
        """
        Return the string representation of the differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)                                      # needs sage.rings.finite_rings sage.rings.function_field
            sage: y.differential()                                                      # needs sage.rings.finite_rings sage.rings.function_field
            (x*y^2 + 1/x*y) d(x)

            sage: F.<x> = FunctionField(QQ)
            sage: f = 1/x
            sage: f.differential()
            (-1/x^2) d(x)
        """
        if self._f.is_zero(): # zero differential
            return '0'

        r = 'd({})'.format(self.parent()._gen_base_differential)

        if self._f.is_one():
            return r

        return '({})'.format(self._f) + ' ' + r

    def _latex_(self):
        r"""
        Return a latex representation of the differential.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: w = y.differential()                                                  # needs sage.rings.function_field
            sage: latex(w)                                                              # needs sage.rings.function_field
            \left( x y^{2} + \frac{1}{x} y \right)\, dx
        """
        if self._f.is_zero(): # zero differential
            return '0'

        r = 'd{}'.format(self.parent()._gen_base_differential)

        if self._f.is_one():
            return r

        return '\\left(' + latex(self._f) + '\\right)\\,' + r

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: {x.differential(): 1}
            {d(x): 1}
            sage: {y.differential(): 1}                                                 # needs sage.rings.function_field
            {(x*y^2 + 1/x*y) d(x): 1}
        """
        return hash((self.parent(), self._f))

    def _richcmp_(self, other, op):
        """
        Compare the differential and the other differential with respect to the
        comparison operator.

        INPUT:

        - ``other`` -- differential

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = L(x).differential()
            sage: w3 = (x*y).differential()
            sage: w1 < w2
            False
            sage: w2 < w1
            True
            sage: w3 == x * w1 + y * w2
            True

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = ((x^2+x+1)^10).differential()
            sage: w2 = (x^2+x+1).differential()
            sage: w1 < w2
            False
            sage: w1 > w2
            True
            sage: w1 == 10*(x^2+x+1)^9 * w2
            True
        """
        return richcmp(self._f, other._f, op)

    def _add_(self, other):
        """
        Return the sum of the differential and the other differential.

        INPUT:

        - ``other`` -- differential

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = (1/y).differential()
            sage: w1 + w2
            (((x^3 + 1)/x^2)*y^2 + 1/x*y) d(x)

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = (1/x).differential()
            sage: w2 = (x^2+1).differential()
            sage: w1 + w2
            ((2*x^3 - 1)/x^2) d(x)
        """
        W = self.parent()
        return W.element_class(W, self._f + other._f)

    def _div_(self, other):
        """
        Return the quotient of ``self`` and ``other``.

        INPUT:

        - ``other`` -- differential

        OUTPUT: an element of the function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = (1/y).differential()
            sage: w1 / w2
            y^2

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = (1/x).differential()
            sage: w2 = (x^2+1).differential()
            sage: w1 / w2
            -1/2/x^3
        """
        if other._f.is_zero():
            raise ZeroDivisionError("division by zero differential")

        return self._f / other._f

    def _neg_(self):
        """
        Return the negation of the differential.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = (-y).differential()
            sage: -w1 == w2
            True

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = (1/x).differential()
            sage: w2 = (-1/x).differential()
            sage: -w1 == w2
            True
        """
        W = self.parent()
        return W.element_class(W, -self._f)

    def _rmul_(self, f):
        """
        Return the differential multiplied by the element of the function
        field.

        INPUT:

        - ``f`` -- element of the function field

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: w1 = (1/y).differential()
            sage: w2 = (-1/y^2) * y.differential()
            sage: w1 == w2
            True

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = (x^2*(x^2+x+1)).differential()
            sage: w2 = (x^2).differential()
            sage: w3 = (x^2+x+1).differential()
            sage: w1 == (x^2) * w3 + (x^2+x+1) * w2
            True
        """
        W = self.parent()
        return W.element_class(W, f * self._f)

    def _acted_upon_(self, f, self_on_left):
        """
        Define multiplication of ``self`` with an element ``f`` of a function
        field, by coercing ``self`` to the space of differentials of the
        function field.

        INPUT:

        - ``f`` -- an element of a function field

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(31)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x); _.<Z> = L[]
            sage: M.<z> = L.extension(Z^2 - y)
            sage: z.differential()
            (8/x*z) d(x)
            sage: 1/(2*z) * y.differential()
            (8/x*z) d(x)

            sage: # needs sage.rings.function_field
            sage: z * x.differential()
            (z) d(x)
            sage: z * (y^2).differential()
            (z) d(x)
            sage: z * (z^4).differential()
            (z) d(x)

        ::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.finite_rings sage.rings.function_field
            sage: y * x.differential()                                                  # needs sage.rings.finite_rings sage.rings.function_field
            (y) d(x)
        """
        F = f.parent()
        if self.parent().function_field() is not F:
            phi = F.space_of_differentials().coerce_map_from(self.parent())
            if phi is not None:
                return phi(self)._rmul_(f)
        raise TypeError

    def divisor(self):
        """
        Return the divisor of the differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: w = (1/y) * y.differential()                                          # needs sage.rings.function_field
            sage: w.divisor()                                                           # needs sage.rings.function_field
            - Place (1/x, 1/x^3*y^2 + 1/x)
             - Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             - Place (x, y)
             + Place (x + 2, y + 3)
             + Place (x^6 + 3*x^5 + 4*x^4 + 2*x^3 + x^2 + 3*x + 4, y + x^5)

        ::

            sage: F.<x> = FunctionField(QQ)
            sage: w = (1/x).differential()
            sage: w.divisor()                                                           # needs sage.libs.pari
            -2*Place (x)
        """
        F = self.parent().function_field()
        x = F.base_field().gen()
        return self._f.divisor() + (-2)*F(x).divisor_of_poles() + F.different()

    def valuation(self, place):
        """
        Return the valuation of the differential at the place.

        INPUT:

        - ``place`` -- a place of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: w = (1/y) * y.differential()                                          # needs sage.rings.function_field
            sage: [w.valuation(p) for p in L.places()]                                  # needs sage.rings.function_field
            [-1, -1, -1, 0, 1, 0]
        """
        F = self.parent().function_field()
        x = F.base_field().gen()
        return (self._f.valuation(place) + 2*min(F(x).valuation(place), 0)
                + F.different().valuation(place))

    def residue(self, place):
        """
        Return the residue of the differential at the place.

        INPUT:

        - ``place`` -- a place of the function field

        OUTPUT: an element of the residue field of the place

        EXAMPLES:

        We verify the residue theorem in a rational function field::

            sage: # needs sage.rings.finite_rings
            sage: F.<x> = FunctionField(GF(4))
            sage: f = 0
            sage: while f == 0:
            ....:     f = F.random_element()
            sage: w = 1/f * f.differential()
            sage: d = f.divisor()
            sage: s = d.support()
            sage: sum([w.residue(p).trace() for p in s])                                # needs sage.rings.function_field
            0

        and in an extension field::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(7)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: f = 0
            sage: while f == 0:
            ....:     f = L.random_element()
            sage: w = 1/f * f.differential()
            sage: d = f.divisor()
            sage: s = d.support()
            sage: sum([w.residue(p).trace() for p in s])
            0

        and also in a function field of characteristic zero::

            sage: # needs sage.rings.function_field
            sage: R.<x> = FunctionField(QQ)
            sage: L.<Y> = R[]
            sage: F.<y> = R.extension(Y^2 - x^4 - 4*x^3 - 2*x^2 - 1)
            sage: a = 6*x^2 + 5*x + 7
            sage: b = 2*x^6 + 8*x^5 + 3*x^4 - 4*x^3 - 1
            sage: w = y*a/b*x.differential()
            sage: d = w.divisor()
            sage: sum([QQ(w.residue(p)) for p in d.support()])
            0
        """
        R,fr_R,to_R = place._residue_field()

        # Step 1: compute f such that fds equals this differential.
        s = place.local_uniformizer()
        dxds = ~(s.derivative())
        g = self._f * dxds

        # Step 2: compute c that is the coefficient of s^-1 in
        # the power series expansion of f
        r = g.valuation(place)
        if r >= 0:
            return R.zero()
        else:
            g_shifted = g * s**(-r)
            c = g_shifted.higher_derivative(-r-1, s)
            return to_R(c)

    def monomial_coefficients(self, copy=True):
        """
        Return a dictionary whose keys are indices of basis elements in the
        support of ``self`` and whose values are the corresponding coefficients.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: d = y.differential()                                                  # needs sage.rings.function_field
            sage: d                                                                     # needs sage.rings.function_field
            ((4*x/(x^7 + 3))*y^2 + ((4*x^7 + 1)/(x^8 + 3*x))*y + x^4/(x^7 + 3)) d(x)
            sage: d.monomial_coefficients()                                             # needs sage.rings.function_field
            {0: (4*x/(x^7 + 3))*y^2 + ((4*x^7 + 1)/(x^8 + 3*x))*y + x^4/(x^7 + 3)}
        """
        return {0: self._f}


class FunctionFieldDifferential_global(FunctionFieldDifferential):
    """
    Differentials on global function fields.

    EXAMPLES::

        sage: F.<x> = FunctionField(GF(7))
        sage: f = x/(x^2 + x + 1)
        sage: f.differential()
        ((6*x^2 + 1)/(x^4 + 2*x^3 + 3*x^2 + 2*x + 1)) d(x)

    ::

        sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                      # needs sage.rings.finite_rings sage.rings.function_field
        sage: y.differential()                                                          # needs sage.rings.finite_rings sage.rings.function_field
        (x*y^2 + 1/x*y) d(x)
    """
    def cartier(self):
        r"""
        Return the image of the differential by the Cartier operator.

        The Cartier operator operates on differentials. Let `x` be a separating
        element of the function field.  If a differential `\omega` is written
        in prime-power representation
        `\omega=(f_0^p+f_1^px+\dots+f_{p-1}^px^{p-1})dx`, then the Cartier
        operator maps `\omega` to `f_{p-1}dx`. It is known that this definition
        does not depend on the choice of `x`.

        The Cartier operator has interesting properties. Notably, the set of
        exact differentials is precisely the kernel of the Cartier operator and
        logarithmic differentials are stable under the Cartier operation.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: f = x/y
            sage: w = 1/f*f.differential()
            sage: w.cartier() == w
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<x> = FunctionField(GF(4))
            sage: f = x/(x^2 + x + 1)
            sage: w = 1/f*f.differential()
            sage: w.cartier() == w                                                      # needs sage.rings.function_field
            True
        """
        W = self.parent()
        F = W.function_field()
        der = F.higher_derivation()
        power_repr = der._prime_power_representation(self._f)
        return W.element_class(W, power_repr[-1])


class DifferentialsSpace(UniqueRepresentation, Parent):
    """
    Space of differentials of a function field.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                      # needs sage.rings.finite_rings sage.rings.function_field
        sage: L.space_of_differentials()                                                # needs sage.rings.finite_rings sage.rings.function_field
        Space of differentials of Function field in y defined by y^3 + x^3*y + x

    The space of differentials is a one-dimensional module over the function
    field. So a base differential is chosen to represent differentials.
    Usually the generator of the base rational function field is a separating
    element and used to generate the base differential. Otherwise a separating
    element is automatically found and used to generate the base differential
    relative to which other differentials are denoted::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(GF(5))
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - 1/x)
        sage: L(x).differential()
        0
        sage: y.differential()
        d(y)
        sage: (y^2).differential()
        (2*y) d(y)
    """
    Element = FunctionFieldDifferential

    def __init__(self, field, category=None):
        """
        Initialize the space of differentials of the function field.

        TESTS::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: W = L.space_of_differentials()
            sage: TestSuite(W).run()
        """
        Parent.__init__(self, base=field, category=Modules(field).FiniteDimensional().WithBasis().or_subcategory(category))

        # Starting from the base rational function field, find the first
        # generator x of an intermediate function field that doesn't map to zero
        # by the derivation map and use dx as our base differential.

        der = field.derivation()
        for F in reversed(field._intermediate_fields(field.rational_function_field())):
            if der(F.gen()) != 0:
                break

        self._derivation = der
        self._gen_base_differential = F.gen()
        self._gen_derivative_inv = ~der(F.gen())  # used for fast computation

    def _repr_(self):
        """
        Return the string representation of the space of differentials.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: w = y.differential()
            sage: w.parent()
            Space of differentials of Function field in y defined by y^3 + x^3*y + x
        """
        return "Space of differentials of {}".format(self.base())

    def _element_constructor_(self, f):
        """
        Construct differential `df` in the space from `f`.

        INPUT:

        - ``f`` -- element of the function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S(y)
            (x*y^2 + 1/x*y) d(x)
            sage: S(y) in S
            True
            sage: S(1)
            0
        """
        if f in self.base():
            return self.element_class(self, self.base().one(), f)

        raise ValueError

    def _coerce_map_from_(self, S):
        """
        Define coercions.

        We can coerce from any DifferentialsSpace whose underlying field
        can be coerced into our underlying field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)                                            # needs sage.rings.function_field
            sage: L.space_of_differentials().coerce_map_from(K.space_of_differentials())            # needs sage.rings.function_field
            Inclusion morphism:
              From: Space of differentials of Rational function field in x over Rational Field
              To:   Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3
        """
        if isinstance(S, DifferentialsSpace):
            if self.function_field().has_coerce_map_from(S.function_field()):
                return DifferentialsSpaceInclusion(S, self)

    def function_field(self):
        """
        Return the function field to which the space of differentials
        is attached.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: S = L.space_of_differentials()
            sage: S.function_field()
            Function field in y defined by y^3 + x^3*y + x
        """
        return self.base()

    def _an_element_(self):
        """
        Return a differential.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S.an_element()  # random
            (x*y^2 + 1/x*y) d(x)
        """
        F = self.base()
        return self.element_class(self, F.one(), F.an_element())

    def basis(self):
        """
        Return a basis.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: S = L.space_of_differentials()
            sage: S.basis()
            Family (d(x),)
        """
        return Family([self.element_class(self, self.base().one())])


class DifferentialsSpace_global(DifferentialsSpace):
    """
    Space of differentials of a global function field.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                                 # needs sage.rings.finite_rings
        sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                      # needs sage.rings.finite_rings sage.rings.function_field
        sage: L.space_of_differentials()                                                # needs sage.rings.finite_rings sage.rings.function_field
        Space of differentials of Function field in y defined by y^3 + x^3*y + x
    """
    Element = FunctionFieldDifferential_global


class DifferentialsSpaceInclusion(Morphism):
    """
    Inclusion morphisms for extensions of function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)                                    # needs sage.rings.function_field
        sage: OK = K.space_of_differentials()
        sage: OL = L.space_of_differentials()                                           # needs sage.rings.function_field
        sage: OL.coerce_map_from(OK)                                                    # needs sage.rings.function_field
        Inclusion morphism:
          From: Space of differentials of Rational function field in x over Rational Field
          To:   Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3
    """

    def _repr_(self):
        """
        Return the string representation of this morphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)                                # needs sage.rings.function_field
            sage: OK = K.space_of_differentials()
            sage: OL = L.space_of_differentials()                                       # needs sage.rings.function_field
            sage: OL.coerce_map_from(OK)                                                # needs sage.rings.function_field
            Inclusion morphism:
              From: Space of differentials of Rational function field in x over Rational Field
              To:   Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3
        """
        s = "Inclusion morphism:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

    def is_injective(self):
        """
        Return ``True``, since the inclusion morphism is injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)                                # needs sage.rings.function_field
            sage: OK = K.space_of_differentials()
            sage: OL = L.space_of_differentials()                                       # needs sage.rings.function_field
            sage: OL.coerce_map_from(OK).is_injective()                                 # needs sage.rings.function_field
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``True`` if the inclusion morphism is surjective.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: OK = K.space_of_differentials()
            sage: OL = L.space_of_differentials()
            sage: OL.coerce_map_from(OK).is_surjective()
            False
            sage: S.<z> = L[]
            sage: M.<z> = L.extension(z - 1)
            sage: OM = M.space_of_differentials()
            sage: OM.coerce_map_from(OL).is_surjective()
            True
        """
        K = self.domain().function_field()
        L = self.codomain().function_field()
        return L.degree(K) == 1

    def _call_(self, v):
        """
        Map ``v`` to a differential in the codomain.

        INPUT:

        - ``v`` -- a differential in the domain

        EXAMPLES::

            sage: # needs sage.rings.function_field sage.rings.number_field
            sage: K.<x> = FunctionField(QQbar); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x*Y + 4*x^3)
            sage: OK = K.space_of_differentials()
            sage: OL = L.space_of_differentials()
            sage: mor = OL.coerce_map_from(OK)
            sage: mor(x.differential()).parent()
            Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3
        """
        domain = self.domain()
        F = self.codomain().function_field()
        return F(v._f)*F(domain._gen_base_differential).differential()
