# sage.doctest: optional - sage.rings.function_field
r"""
Elements of function fields: extension
"""

#*****************************************************************************
#       Copyright (C) 2023 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp cimport richcmp, richcmp_not_equal
from sage.structure.element cimport FieldElement, RingElement, ModuleElement, Element

from sage.rings.function_field.element cimport FunctionFieldElement


cdef class FunctionFieldElement_polymod(FunctionFieldElement):
    """
    Elements of a finite extension of a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: x*y + 1/x^3
        x*y + 1/x^3
    """
    def __init__(self, parent, x, reduce=True):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: TestSuite(x*y + 1/x^3).run()
        """
        FieldElement.__init__(self, parent)
        if reduce:
            self._x = x % self._parent.polynomial()
        else:
            self._x = x

    def element(self):
        """
        Return the underlying polynomial that represents the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<T> = K[]
            sage: L.<y> = K.extension(T^2 - x*T + 4*x^3)
            sage: f = y/x^2 + x/(x^2+1); f
            1/x^2*y + x/(x^2 + 1)
            sage: f.element()
            1/x^2*y + x/(x^2 + 1)
        """
        return self._x

    def _repr_(self):
        """
        Return the string representation of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y._repr_()
            'y'
        """
        return self._x._repr(name=self.parent().variable_name())

    def __bool__(self):
        """
        Return True if the element is not zero.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: bool(y)
            True
            sage: bool(L(0))
            False
            sage: bool(L.coerce(L.polynomial()))
            False
        """
        return not not self._x

    def __hash__(self):
        """
        Return the hash of the element.

        TESTS::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: len({hash(y^i+x^j) for i in [-2..2] for j in [-2..2]}) >= 24
            True
        """
        return hash(self._x)

    cpdef _richcmp_(self, other, int op):
        """
        Do rich comparison with the other element with respect to ``op``

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: L(0) == 0
            True
            sage: y != L(2)
            True
        """
        cdef FunctionFieldElement left = <FunctionFieldElement>self
        cdef FunctionFieldElement right = <FunctionFieldElement>other
        return richcmp(left._x, right._x, op)

    cpdef _add_(self, right):
        """
        Add the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*y + x/(1+x^3))  +  (3*y + 5*x*y)         # indirect doctest
            (5*x + 5)*y + x/(x^3 + 1)
            sage: (y^2 - x*y + 4*x^3)==0                      # indirect doctest
            True
            sage: -y + y
            0
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x + (<FunctionFieldElement>right)._x
        return res

    cpdef _sub_(self, right):
        """
        Subtract the other element from the element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*y + x/(1+x^3))  -  (3*y + 5*x*y)         # indirect doctest
            (-5*x - 1)*y + x/(x^3 + 1)
            sage: y - y
            0
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x - (<FunctionFieldElement>right)._x
        return res

    cpdef _mul_(self, right):
        """
        Multiply the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y  *  (3*y + 5*x*y)                          # indirect doctest
            (5*x^2 + 3*x)*y - 20*x^4 - 12*x^3
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = (self._x * (<FunctionFieldElement>right)._x) % self._parent.polynomial()
        return res

    cpdef _div_(self, right):
        """
        Divide the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*y + x/(1+x^3))  /  (2*y + x/(1+x^3))       # indirect doctest
            1
            sage: 1 / (y^2 - x*y + 4*x^3)                       # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert 0
        """
        return self * ~right

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*y + 1/x); a                           # indirect doctest
            (-1/8*x^2/(x^5 + 1/8*x^2 + 1/16))*y + (1/8*x^3 + 1/16*x)/(x^5 + 1/8*x^2 + 1/16)
            sage: a*(2*y + 1/x)
            1
        """
        if self.is_zero():
            raise ZeroDivisionError("Cannot invert 0")
        P = self._parent
        return P(self._x.xgcd(P._polynomial)[1])

    cpdef list list(self):
        """
        Return the list of the coefficients representing the element.

        If the function field is `K[y]/(f(y))`, then return the coefficients of
        the reduced presentation of the element as a polynomial in `K[y]`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*y + 1/x); a
            (-1/8*x^2/(x^5 + 1/8*x^2 + 1/16))*y + (1/8*x^3 + 1/16*x)/(x^5 + 1/8*x^2 + 1/16)
            sage: a.list()
            [(1/8*x^3 + 1/16*x)/(x^5 + 1/8*x^2 + 1/16), -1/8*x^2/(x^5 + 1/8*x^2 + 1/16)]
            sage: (x*y).list()
            [0, x]
        """
        return self._x.padded_list(self._parent.degree())

    cpdef FunctionFieldElement nth_root(self, n):
        r"""
        Return an ``n``-th root of this element in the function field.

        INPUT:

        - ``n`` -- an integer

        OUTPUT:

        Returns an element ``a`` in the function field such that this element
        equals `a^n`. Raises an error if no such element exists.

        ALGORITHM:

        If ``n`` is a power of the characteristic of the field and the constant
        base field is perfect, then this uses the algorithm described in
        Proposition 12 of [GiTr1996]_.

        .. SEEALSO::

            :meth:`is_nth_power`

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))                                          # optional - sage.rings.finite_rings
            sage: R.<y> = K[]                                                           # optional - sage.rings.finite_rings
            sage: L.<y> = K.extension(y^2 - x)                                          # optional - sage.rings.finite_rings
            sage: L(y^3).nth_root(3)                                                    # optional - sage.rings.finite_rings
            y
            sage: L(y^9).nth_root(-9)                                                   # optional - sage.rings.finite_rings
            1/x*y

        This also works for inseparable extensions::

            sage: K.<x> = FunctionField(GF(3))                                          # optional - sage.rings.finite_rings
            sage: R.<y> = K[]                                                           # optional - sage.rings.finite_rings
            sage: L.<y> = K.extension(y^3 - x^2)                                        # optional - sage.rings.finite_rings
            sage: L(x).nth_root(3)^3                                                    # optional - sage.rings.finite_rings
            x
            sage: L(x^9).nth_root(-27)^-27                                              # optional - sage.rings.finite_rings
            x^9

        """
        if n == 1:
            return self
        if n < 0:
            return (~self).nth_root(-n)
        if n == 0:
            if not self.is_one():
                raise ValueError("element is not a 0-th power")
            return self

        # reduce to the separable case
        poly = self._parent._polynomial
        if not poly.gcd(poly.derivative()).is_one():
            L, from_L, to_L = self._parent.separable_model(('t', 'w'))
            return from_L(to_L(self).nth_root(n))

        constant_base_field = self._parent.constant_base_field()
        p = constant_base_field.characteristic()
        if p.divides(n) and constant_base_field.is_perfect():
            return self._pth_root().nth_root(n//p)

        raise NotImplementedError("nth_root() not implemented for this n")

    cpdef bint is_nth_power(self, n):
        r"""
        Return whether this element is an ``n``-th power in the function field.

        INPUT:

        - ``n`` -- an integer

        ALGORITHM:

        If ``n`` is a power of the characteristic of the field and the constant
        base field is perfect, then this uses the algorithm described in
        Proposition 12 of [GiTr1996]_.

        .. SEEALSO::

            :meth:`nth_root`

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))                                          # optional - sage.rings.finite_rings
            sage: R.<y> = K[]                                                           # optional - sage.rings.finite_rings
            sage: L.<y> = K.extension(y^2 - x)                                          # optional - sage.rings.finite_rings
            sage: y.is_nth_power(2)                                                     # optional - sage.rings.finite_rings
            False
            sage: L(x).is_nth_power(2)                                                  # optional - sage.rings.finite_rings
            True

        """
        if n == 0:
            return self.is_one()
        if n == 1:
            return True
        if n < 0:
            return self.is_unit() and (~self).is_nth_power(-n)

        # reduce to the separable case
        poly = self._parent._polynomial
        if not poly.gcd(poly.derivative()).is_one():
            L, from_L, to_L = self._parent.separable_model(('t', 'w'))
            return to_L(self).is_nth_power(n)

        constant_base_field = self._parent.constant_base_field()
        p = constant_base_field.characteristic()
        if p.divides(n) and constant_base_field.is_perfect():
            return self._parent.derivation()(self).is_zero() and self._pth_root().is_nth_power(n//p)

        raise NotImplementedError("is_nth_power() not implemented for this n")

    cdef FunctionFieldElement _pth_root(self):
        r"""
        Helper method for :meth:`nth_root` and :meth:`is_nth_power` which
        computes a `p`-th root if the characteristic is `p` and the constant
        base field is perfect.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))                                          # optional - sage.rings.finite_rings
            sage: R.<y> = K[]                                                           # optional - sage.rings.finite_rings
            sage: L.<y> = K.extension(y^2 - x)                                          # optional - sage.rings.finite_rings
            sage: (y^3).nth_root(3)  # indirect doctest                                 # optional - sage.rings.finite_rings
            y
        """
        cdef Py_ssize_t deg = self._parent.degree()
        if deg == 1:
            return self._parent(self._x[0].nth_root(self._parent.characteristic()))

        from .function_field_rational import RationalFunctionField
        if not isinstance(self.base_ring(), RationalFunctionField):
            raise NotImplementedError("only implemented for simple extensions of function fields")
        # compute a representation of the generator y of the field in terms of powers of y^p
        cdef Py_ssize_t i
        cdef list v = []
        char = self._parent.characteristic()
        cdef FunctionFieldElement_polymod yp = self._parent.gen() ** char
        val = self._parent.one()._x
        poly = self._parent.polynomial()
        for i in range(deg):
            v += val.padded_list(deg)
            val = (val * yp._x) % poly
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(self._parent._base, deg)
        M = MS(v)
        y = self._parent._base.polynomial_ring()(M.solve_left(MS.column_space()([0,1]+[0]*(deg-2))).list())

        f = self._x(y).map_coefficients(lambda c: c.nth_root(char))
        return self._parent(f)
