r"""
Elements of function fields

Sage provides arithmetic with elements of function fields.

EXAMPLES:

Arithmetic with rational functions::

    sage: K.<t> = FunctionField(QQ)
    sage: f = t - 1
    sage: g = t^2 - 3
    sage: h = f^2/g^3
    sage: h.valuation(t-1)
    2
    sage: h.valuation(t)
    0
    sage: h.valuation(t^2 - 3)
    -3

Derivatives of elements in separable extensions::

    sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                                     # needs sage.rings.finite_rings
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                        # needs sage.rings.finite_rings sage.rings.function_field
    sage: (y^3 + x).derivative()                                                        # needs sage.rings.finite_rings sage.rings.function_field
    ((x^2 + 1)/x^2)*y + (x^4 + x^3 + 1)/x^3

The divisor of an element of a global function field::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                        # needs sage.rings.function_field
    sage: y.divisor()                                                                   # needs sage.rings.function_field
    - Place (1/x, 1/x*y)
     - Place (x, x*y)
     + 2*Place (x + 1, x*y)

AUTHORS:

- William Stein: initial version

- Robert Bradshaw (2010-05-27): cythonize function field elements

- Julian Rueth (2011-06-28, 2020-09-01): treat zero correctly; implement nth_root/is_nth_power

- Maarten Derickx (2011-09-11): added doctests, fixed pickling

- Kwankyu Lee (2017-04-30): added elements for global function fields

- Vincent Macri (2024-09-03): added subs method
"""
# *****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2010      Robert Bradshaw <robertwb@math.washington.edu>
#                     2011-2020 Julian Rueth <julian.rueth@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2015      Nils Bruin
#                     2016      Frédéric Chapoton
#                     2017-2019 Kwankyu Lee
#                     2018-2020 Travis Scrimshaw
#                     2019      Brent Baccala
#                     2021      Saher Amasha
#                     2024      Vincent Macri
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.function_fields import FunctionFields
from sage.misc.cachefunc import cached_method
from sage.structure.element cimport FieldElement


def is_FunctionFieldElement(x):
    """
    Return ``True`` if ``x`` is any type of function field element.

    EXAMPLES::

        sage: t = FunctionField(QQ,'t').gen()
        sage: sage.rings.function_field.element.is_FunctionFieldElement(t)
        doctest:warning...
        DeprecationWarning: The function is_FunctionFieldElement is deprecated;
        use '....parent() in FunctionFields()' instead.
        See https://github.com/sagemath/sage/issues/38289 for details.
        True
        sage: sage.rings.function_field.element.is_FunctionFieldElement(0)
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38289,
                       "The function is_FunctionFieldElement is deprecated; "
                       "use '....parent() in FunctionFields()' instead.")
    if isinstance(x, FunctionFieldElement):
        return True
    from sage.rings.function_field.function_field import FunctionField
    if isinstance(x.parent(), FunctionField):
        return True
    return x.parent() in FunctionFields()


def make_FunctionFieldElement(parent, element_class, representing_element):
    """
    Used for unpickling FunctionFieldElement objects (and subclasses).

    EXAMPLES::

        sage: from sage.rings.function_field.element import make_FunctionFieldElement
        sage: K.<x> = FunctionField(QQ)
        sage: make_FunctionFieldElement(K, K.element_class, (x+1)/x)
        (x + 1)/x
    """
    return element_class(parent, representing_element, reduce=False)


cdef class FunctionFieldElement(FieldElement):
    """
    Abstract base class for function field elements.

    EXAMPLES::

        sage: t = FunctionField(QQ,'t').gen()
        sage: isinstance(t, sage.rings.function_field.element.FunctionFieldElement)
        True
    """
    def __reduce__(self):
        """
        EXAMPLES::

            sage: K = FunctionField(QQ,'x')
            sage: f = K.random_element()
            sage: loads(f.dumps()) == f
            True
        """
        return (make_FunctionFieldElement,
                (self._parent, type(self), self._x))

    cdef FunctionFieldElement _new_c(self):
        cdef type t = type(self)
        cdef FunctionFieldElement x = <FunctionFieldElement>t.__new__(t)
        x._parent = self._parent
        return x

    def __pari__(self):
        r"""
        Coerce the element to PARI.

        PARI does not know about general function field elements, so this
        raises an Exception.

        TESTS:

        Check that :issue:`16369` has been resolved::

            sage: K.<a> = FunctionField(QQ)
            sage: R.<b> = K[]
            sage: L.<b> = K.extension(b^2 - a)                                          # needs sage.rings.function_field
            sage: b.__pari__()                                                          # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            NotImplementedError: PARI does not support general function field elements.
        """
        raise NotImplementedError("PARI does not support general function field elements.")

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: latex((t+1)/t)
            \frac{t + 1}{t}
            sage: latex((t+1)/t^67)
            \frac{t + 1}{t^{67}}
            sage: latex((t+1/2)/t^67)
            \frac{t + \frac{1}{2}}{t^{67}}
        """
        return self._x._latex_()

    def subs(self, in_dict=None, **kwds):
        r"""
        Substitute the given generators with given values while not touching
        other generators.

        INPUT:

        - ``in_dict`` -- (optional) dictionary of inputs

        - ``**kwds`` -- named parameters

        OUTPUT: new object if substitution is possible, otherwise ``self``

        EXAMPLES:

        Basic substitution::

            sage: K = GF(7)
            sage: Kx.<x> = FunctionField(K)
            sage: y = polygen(Kx)
            sage: f = x^6 + 3; f
            x^6 + 3

        We also substitute the generators in any base fields::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: S.<t> = L[]
            sage: M.<t> = L.extension(t^2 - x*y)
            sage: f = 7 * t + 3*x*y
            sage: f.subs(t=9)
            3*x*y + 63
            sage: f.subs(x=2, y=4)
            7*t + 24
            sage: f.subs(t=1, x=2, y=3)
            25

        Because of the possibility of extension fields, a generator to
        substitute must be specified::

            sage: K.<x> = FunctionField(QQ)
            sage: f = x
            sage: f.subs(2)
            Traceback (most recent call last):
            ...
            TypeError: in_dict must be a dict
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: f = x + y
            sage: f.subs(0)
            Traceback (most recent call last):
            ...
            TypeError: in_dict must be a dict

        We can also substitute using dictionary syntax::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: S.<t> = L[]
            sage: M.<t> = L.extension(t^2 - x*y)
            sage: f = x + y + t
            sage: f.subs({x: 1, y: 3, t: 4})
            8
            sage: f.subs({x: 1, t: 4})
            y + 5

        TESTS:

        Check that we correctly handle extension fields::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: S.<t> = L[]
            sage: M.<t> = L.extension(t^2 - x*y)
            sage: f = t + x*y
            sage: f.subs(x=1, y=3, t=5)
            8
            sage: f_sub = f.subs(x=1); f_sub
            t + y
            sage: f_sub.parent() == f.parent()
            True
            sage: f.subs(y=2)
            t + 2*x
            sage: f_sub = f.subs(x=1, y=1, t=1); f_sub
            2
            sage: f_sub.parent() == M
            True

        Test that substitution works for rational functions::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^4 - 3)
            sage: f = x / y
            sage: f.subs(x=2) == 2 / y
            True
            sage: f.subs(y=3)
            9*x
            sage: f.subs(t=-1) is f
            True
            sage: f.subs({x: 2, y: 4})
            128/3

        Make sure that we return the same object when there is no
        substitution::

            sage: K = GF(7)
            sage: Kx.<x> = FunctionField(K)
            sage: y = polygen(Kx)
            sage: f = x^6 + 3
            sage: g = f.subs(z=2)
            sage: g == f
            True
            sage: g is f
            True

        Same purpose as above but over an extension field over the rationals::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: S.<t> = L[]
            sage: M.<t> = L.extension(t^2 - x*y)
            sage: f = t + x*y
            sage: f.subs() is f
            True
            sage: f.subs(w=7) is f
            True
            sage: f.subs(w=7) is f.subs(w=7)
            True
            sage: f.subs(y=y) is f
            True
            sage: f.subs({y: y}) is f
            True
            sage: f.subs(x=x, y=y, t=t) is f
            True

        Test proper handling of not making substitutions::

            sage: K.<x> = FunctionField(QQ)
            sage: f = x
            sage: f.subs() is f
            True
            sage: f.subs(dict()) is f
            True
            sage: f.subs(w=0) is f
            True
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x))
            sage: f = 3*y
            sage: f.subs(x=0)
            3*y
            sage: f = 3*y
            sage: f.subs(x=0, y=y)
            3*y

        Test error handling for wrong argument type::

            sage: K.<x> = FunctionField(QQ)
            sage: f = x
            sage: f.subs(0)
            Traceback (most recent call last):
            ...
            TypeError: in_dict must be a dict

        Test error handling for dictionary with keys that don't match
        generators::

            sage: K.<x> = FunctionField(QQ)
            sage: f = x
            sage: f.subs({1: 1})
            Traceback (most recent call last):
            ...
            TypeError: key does not match any field generators

        Test error handling with ambiguously named generators::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = K[]
            sage: L.<x> = K.extension(x^3 - x)
            sage: str(L.gen()) == str(K.gen())
            True
            sage: f = K.gen() - L.gen()
            sage: f.subs(x=2)
            Traceback (most recent call last):
            ...
            TypeError: multiple generators have the same name, making substitution ambiguous. Rename generators or pass substitution values in using dictionary format
            sage: f.subs({K.gen(): 1})
            -x + 1
            sage: f.subs({L.gen(): 2})
            x - 2
            sage: f.subs({K.gen(): 1, L.gen(): 2})
            -1
            sage: f.subs({K.gen(): 2, L.gen(): 1})
            1
        """
        # Helper method to do the recursion through base fields.

        def sub_recurse(self, sub_dict):
            ff = self.parent()
            if ff.base_field() == ff:
                return ff(self._x.subs({ff.gen(): sub_dict[ff.gen()]}))
            total = ff.zero()
            for i, v in enumerate(list(self._x)):
                total += sub_recurse(v, sub_dict) * sub_dict[ff.gen()]**i
            return ff(total)

        if in_dict is None and kwds is None:
            return self

        if in_dict is not None and not isinstance(in_dict, dict):
            raise TypeError('in_dict must be a dict')

        field_tower = [self.parent()]
        ff = self.parent()

        while ff.base_field() != ff:
            ff = ff.base_field()
            field_tower.append(ff)
        sub_dict = {f.gen(): f.gen() for f in field_tower}

        made_substitution = False
        if in_dict is not None:
            for k, v in in_dict.items():
                if k in sub_dict:
                    sub_dict[k] = v
                    if v != k:
                        made_substitution = True
                else:
                    raise TypeError('key does not match any field generators')
        else:
            used_kwds = {k: False for k in kwds}
            for g in sub_dict:
                strg = str(g)
                if strg not in kwds:
                    continue
                v = kwds[strg]
                        sub_dict[g] = v
                        if used_kwds[k]:
                            raise TypeError('multiple generators have the '
                                            'same name, making substitution '
                                            'ambiguous. Rename generators '
                                            'or pass substitution values in '
                                            'using dictionary format')
                        used_kwds[k] = True
                        if g != v:
                            made_substitution = True

        if made_substitution:
            return sub_recurse(self, sub_dict)
        return self

    @cached_method
    def matrix(self, base=None):
        r"""
        Return the matrix of multiplication by this element, interpreting this
        element as an element of a vector space over ``base``.

        INPUT:

        - ``base`` -- a function field (default: ``None``); if ``None``, then
          the matrix is formed over the base field of this function field

        EXAMPLES:

        A rational function field::

            sage: K.<t> = FunctionField(QQ)
            sage: t.matrix()                                                            # needs sage.modules
            [t]
            sage: (1/(t+1)).matrix()                                                    # needs sage.modules
            [1/(t + 1)]

        Now an example in a nontrivial extension of a rational function field::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y.matrix()
            [     0      1]
            [-4*x^3      x]
            sage: y.matrix().charpoly('Z')
            Z^2 - x*Z + 4*x^3

        An example in a relative extension, where neither function
        field is rational::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: M.<T> = L[]
            sage: Z.<alpha> = L.extension(T^3 - y^2*T + x)
            sage: alpha.matrix()
            [          0           1           0]
            [          0           0           1]
            [         -x x*y - 4*x^3           0]
            sage: alpha.matrix(K)
            [           0            0            1            0            0            0]
            [           0            0            0            1            0            0]
            [           0            0            0            0            1            0]
            [           0            0            0            0            0            1]
            [          -x            0       -4*x^3            x            0            0]
            [           0           -x       -4*x^4 -4*x^3 + x^2            0            0]
            sage: alpha.matrix(Z)
            [alpha]

        We show that this matrix does indeed work as expected when making a
        vector space from a function field::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: V, from_V, to_V = L.vector_space()
            sage: y5 = to_V(y^5); y5
            ((x^4 + 1)/x, 2*x, 0, 0, 0)
            sage: y4y = to_V(y^4) * y.matrix(); y4y
            ((x^4 + 1)/x, 2*x, 0, 0, 0)
            sage: y5 == y4y
            True
        """
        # multiply each element of the vector space isomorphic to the parent
        # with this element; make matrix whose rows are the coefficients of the
        # result, and transpose
        V, f, t = self.parent().vector_space(base)
        rows = [t(self*f(b)) for b in V.basis()]
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(V.base_field(), V.dimension())
        ret = MS(rows)
        ret.transpose()
        ret.set_immutable()
        return ret

    def trace(self):
        """
        Return the trace of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)                                # needs sage.rings.function_field
            sage: y.trace()                                                             # needs sage.rings.function_field
            x
        """
        return self.matrix().trace()

    def norm(self):
        """
        Return the norm of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)                                # needs sage.rings.function_field
            sage: y.norm()                                                              # needs sage.rings.function_field
            4*x^3

        The norm is relative::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3); R.<z> = L[]                   # needs sage.rings.function_field
            sage: M.<z> = L.extension(z^3 - y^2*z + x)                                  # needs sage.rings.function_field
            sage: z.norm()                                                              # needs sage.rings.function_field
            -x
            sage: z.norm().parent()                                                     # needs sage.rings.function_field
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self.matrix().determinant()

    def degree(self):
        """
        Return the max degree between the denominator and numerator.

        EXAMPLES::

            sage: FF.<t> = FunctionField(QQ)
            sage: f = (t^2 + 3) / (t^3 - 1/3); f
            (t^2 + 3)/(t^3 - 1/3)
            sage: f.degree()
            3

            sage: FF.<t> = FunctionField(QQ)
            sage: f = (t+8); f
            t + 8
            sage: f.degree()
            1

        TESTS::

            sage: FF.<t> = FunctionField(QQ)
            sage: f = FF(0); f
            0
            sage: f.degree()
            0
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.degree()
            2
            sage: f = (t+1); f
            t + 1
            sage: f.degree()
            1
        """
        return max(self._x.denominator().degree(), self._x.numerator().degree())

    def characteristic_polynomial(self, *args, **kwds):
        """
        Return the characteristic polynomial of the element. Give an optional
        input string to name the variable in the characteristic polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: x.characteristic_polynomial('W')                                      # needs sage.modules
            W - x

            sage: # needs sage.rings.function_field
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3); R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y^2*z + x)
            sage: y.characteristic_polynomial('W')
            W^2 - x*W + 4*x^3
            sage: z.characteristic_polynomial('W')
            W^3 + (-x*y + 4*x^3)*W + x
        """
        return self.matrix().characteristic_polynomial(*args, **kwds)

    charpoly = characteristic_polynomial

    def minimal_polynomial(self, *args, **kwds):
        """
        Return the minimal polynomial of the element. Give an optional input
        string to name the variable in the characteristic polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: x.minimal_polynomial('W')                                             # needs sage.modules
            W - x

            sage: # needs sage.rings.function_field
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3); R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y^2*z + x)
            sage: y.minimal_polynomial('W')
            W^2 - x*W + 4*x^3
            sage: z.minimal_polynomial('W')
            W^3 + (-x*y + 4*x^3)*W + x
        """
        return self.matrix().minimal_polynomial(*args, **kwds)

    minpoly = minimal_polynomial

    def is_integral(self):
        r"""
        Determine if the element is integral over the maximal order of the base field.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y.is_integral()
            True
            sage: (y/x).is_integral()
            True
            sage: (y/x)^2 - (y/x) + 4*x
            0
            sage: (y/x^2).is_integral()
            False
            sage: (y/x).minimal_polynomial('W')
            W^2 - W + 4*x
        """
        R = self.parent().base_field().maximal_order()
        return all(a in R for a in self.minimal_polynomial())

    def differential(self):
        """
        Return the differential `dx` where `x` is the element.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = 1 / t
            sage: f.differential()                                                      # needs sage.modules
            (-1/t^2) d(t)

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x +1/x)                                 # needs sage.rings.finite_rings sage.rings.function_field
            sage: (y^3 + x).differential()                                              # needs sage.rings.finite_rings sage.rings.function_field
            (((x^2 + 1)/x^2)*y + (x^4 + x^3 + 1)/x^3) d(x)

        TESTS:

        Verify that :issue:`27712` is resolved::

            sage: K.<x> = FunctionField(GF(31))
            sage: x.differential()                                                      # needs sage.modules
            d(x)

            sage: # needs sage.rings.function_field
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: y.differential()
            (16/x*y) d(x)
            sage: z.differential()
            (8/x*z) d(x)
        """
        F = self.parent()
        W = F.space_of_differentials()
        return W.element_class(W, F.one(), self)

    def derivative(self):
        """
        Return the derivative of the element.

        The derivative is with respect to the generator of the base rational
        function field, over which the function field is a separable extension.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t + 1) / (t^2 - 1/3)
            sage: f.derivative()                                                        # needs sage.modules
            (-t^2 - 2*t - 1/3)/(t^4 - 2/3*t^2 + 1/9)

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings sage.rings.function_field
            sage: (y^3 + x).derivative()                                                # needs sage.rings.finite_rings sage.rings.function_field
            ((x^2 + 1)/x^2)*y + (x^4 + x^3 + 1)/x^3
        """
        D = self.parent().derivation()
        return D(self)

    def higher_derivative(self, i, separating_element=None):
        """
        Return the `i`-th derivative of the element with respect to the
        separating element.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``separating_element`` -- a separating element of the function field;
          the default is the generator of the rational function field

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(2))
            sage: f = t^2
            sage: f.higher_derivative(2)                                                # needs sage.rings.function_field
            1

        ::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings sage.rings.function_field
            sage: (y^3 + x).higher_derivative(2)                                        # needs sage.rings.finite_rings sage.rings.function_field
            1/x^3*y + (x^6 + x^4 + x^3 + x^2 + x + 1)/x^5
        """
        D = self.parent().higher_derivation()
        return D(self, i, separating_element)

    @cached_method
    def divisor(self):
        """
        Return the divisor of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.divisor()                                                           # needs sage.libs.pari sage.modules
            3*Place (1/x)
             - Place (x)
             - Place (x^2 + x + 1)

        ::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.function_field
            sage: y.divisor()                                                           # needs sage.rings.function_field
            - Place (1/x, 1/x*y)
             - Place (x, x*y)
             + 2*Place (x + 1, x*y)
        """
        if self.is_zero():
            raise ValueError("divisor not defined for zero")

        F = self.parent()
        I = F.maximal_order().ideal(self)
        J = F.maximal_order_infinite().ideal(self)
        return I.divisor() + J.divisor()

    def divisor_of_zeros(self):
        """
        Return the divisor of zeros for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.divisor_of_zeros()                                                  # needs sage.libs.pari sage.modules
            3*Place (1/x)

        ::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings sage.rings.function_field
            sage: (x/y).divisor_of_zeros()                                              # needs sage.rings.finite_rings sage.rings.function_field
            3*Place (x, x*y)
        """
        if self.is_zero():
            raise ValueError("divisor of zeros not defined for zero")

        F = self.parent()
        I = F.maximal_order().ideal(self)
        J = F.maximal_order_infinite().ideal(self)
        return I.divisor_of_zeros() + J.divisor_of_zeros()

    def divisor_of_poles(self):
        """
        Return the divisor of poles for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.divisor_of_poles()                                                  # needs sage.libs.pari sage.modules
            Place (x)
             + Place (x^2 + x + 1)

        ::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings sage.rings.function_field
            sage: (x/y).divisor_of_poles()                                              # needs sage.rings.finite_rings sage.rings.function_field
            Place (1/x, 1/x*y) + 2*Place (x + 1, x*y)
        """
        if self.is_zero():
            raise ValueError("divisor of poles not defined for zero")

        F = self.parent()
        I = F.maximal_order().ideal(self)
        J = F.maximal_order_infinite().ideal(self)
        return I.divisor_of_poles() + J.divisor_of_poles()

    def zeros(self):
        """
        Return the list of the zeros of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.zeros()                                                             # needs sage.libs.pari sage.modules
            [Place (1/x)]

        ::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings sage.rings.function_field
            sage: (x/y).zeros()                                                         # needs sage.rings.finite_rings sage.rings.function_field
            [Place (x, x*y)]
        """
        return self.divisor_of_zeros().support()

    def poles(self):
        """
        Return the list of the poles of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.poles()                                                             # needs sage.libs.pari sage.modules
            [Place (x), Place (x^2 + x + 1)]

        ::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]                             # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.finite_rings sage.rings.function_field
            sage: (x/y).poles()                                                         # needs sage.rings.finite_rings sage.rings.function_field
            [Place (1/x, 1/x*y), Place (x + 1, x*y)]
        """
        return self.divisor_of_poles().support()

    def valuation(self, place):
        """
        Return the valuation of the element at the place.

        INPUT:

        - ``place`` -- a place of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.function_field
            sage: p = L.places_infinite()[0]                                            # needs sage.rings.function_field
            sage: y.valuation(p)                                                        # needs sage.rings.function_field
            -1

        ::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x - 1).place()
            sage: y.valuation(p)
            0
        """
        prime = place.prime_ideal()
        ideal = prime.ring().ideal(self)
        return prime.valuation(ideal)

    def evaluate(self, place):
        """
        Return the value of the element at the place.

        INPUT:

        - ``place`` -- a function field place

        OUTPUT:

        If the element is in the valuation ring at the place, then an element
        in the residue field at the place is returned. Otherwise, a
        :exc:`ValueError` is raised.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(5))
            sage: p = K.place_infinite()
            sage: f = 1/t^2 + 3
            sage: f.evaluate(p)
            3

        ::

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p, = L.places_infinite()
            sage: p, = L.places_infinite()
            sage: (y + x).evaluate(p)
            Traceback (most recent call last):
            ...
            ValueError: has a pole at the place
            sage: (y/x + 1).evaluate(p)
            1
        """
        R, _, to_R = place._residue_field()

        v = self.valuation(place)
        if v > 0:
            return R.zero()
        if v == 0:
            return to_R(self)
        # v < 0
        raise ValueError('has a pole at the place')

    cpdef bint is_nth_power(self, n) noexcept:
        r"""
        Return whether this element is an ``n``-th power in the rational
        function field.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        Returns ``True`` if there is an element `a` in the function field such
        that this element equals `a^n`.

        .. SEEALSO::

            :meth:`nth_root`

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: f = (x+1)/(x-1)
            sage: f.is_nth_power(2)
            False
        """
        raise NotImplementedError("is_nth_power() not implemented for generic elements")

    cpdef FunctionFieldElement nth_root(self, n):
        """
        Return an ``n``-th root of this element in the function field.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        Returns an element ``a`` in the function field such that this element
        equals `a^n`. Raises an error if no such element exists.

        .. SEEALSO::

            :meth:`is_nth_power`

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)                                          # needs sage.rings.function_field
            sage: L(y^27).nth_root(27)                                                  # needs sage.rings.function_field
            y
        """
        raise NotImplementedError("nth_root() not implemented for generic elements")
