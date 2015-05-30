"""
Elements of Finite Algebras
"""

#*****************************************************************************
#  Copyright (C) 2011 Johan Bosman <johan.g.bosman@gmail.com>
#  Copyright (C) 2011, 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#  Copyright (C) 2011 Michiel Kosters <kosters@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re

from sage.misc.lazy_attribute import lazy_attribute
from sage.matrix.constructor import Matrix
from sage.matrix.matrix import is_Matrix
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.structure.element import AlgebraElement, is_Vector, parent


class FiniteDimensionalAlgebraElement(AlgebraElement):
    r"""
    Create an element of a :class:`FiniteDimensionalAlgebra` using a multiplication table.

    INPUT:

    - ``A`` -- a :class:`FiniteDimensionalAlgebra` which will be the parent

    - ``elt`` -- vector, matrix or element of the base field
      (default: ``None``)

    - ``check`` -- boolean (default: ``True``); if ``False`` and ``elt`` is a
      matrix, assume that it is known to be the matrix of an element

    If ``elt`` is a vector, it is interpreted as a vector of
    coordinates with respect to the given basis of ``A``.  If ``elt`` is
    a matrix, it is interpreted as a multiplication matrix with
    respect to this basis.

    EXAMPLES::

        sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
        sage: A(17)
        2*e0
        sage: A([1,1])
        e0 + e1
    """
    def __init__(self, A, elt=None, check=True):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: A(QQ(4))
            Traceback (most recent call last):
            ...
            TypeError: elt should be a vector, a matrix, or an element of the base field

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [-1,0]])])
            sage: elt = B(Matrix([[1,1], [-1,1]])); elt
            e0 + e1
            sage: TestSuite(elt).run()
            sage: B(Matrix([[0,1], [1,0]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an element of the algebra
        """
        AlgebraElement.__init__(self, A)
        k = A.base_ring()
        n = A.degree()
        if elt is None:
            self._vector = vector(k, n)
            self._matrix = Matrix(k, n)
        else:
            if isinstance(elt, int):
                elt = Integer(elt)
            elif isinstance(elt, list):
                elt = vector(elt)
            if A == elt.parent():
                self._vector = elt._vector.base_extend(k)
                self._matrix = elt._matrix.base_extend(k)
            elif k.has_coerce_map_from(elt.parent()):
                e = k(elt)
                if e == 0:
                    self._vector = vector(k, n)
                    self._matrix = Matrix(k, n)
                elif A.is_unitary():
                    self._vector = A._one * e
                    self._matrix = Matrix.identity(k, n) * e
                else:
                    raise TypeError("algebra is not unitary")
            elif is_Vector(elt):
                self._vector = elt.base_extend(k)
                self._matrix = Matrix(k, sum([elt[i] * A.table()[i] for i in xrange(n)]))
            elif is_Matrix(elt):
                if not A.is_unitary():
                    raise TypeError("algebra is not unitary")
                self._vector = A._one * elt
                if not check or sum([self._vector[i]*A.table()[i] for i in xrange(n)]) == elt:
                    self._matrix = elt
                else:
                    raise ValueError("matrix does not define an element of the algebra")
            else:
                raise TypeError("elt should be a vector, a matrix, " +
                                "or an element of the base field")

    def vector(self):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(5).vector()
            (5, 0, 5)
        """
        return self._vector

    def matrix(self):
        """
        Return the matrix for multiplication by ``self`` from the right.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(5).matrix()
            [5 0 0]
            [0 5 0]
            [0 0 5]
        """
        return self._matrix

    def left_matrix(self):
        """
        Return the matrix for multiplication by ``self`` from the left.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: C([1,2,0]).left_matrix()
            [1 0 0]
            [0 1 0]
            [0 2 0]

        """
        A = self.parent()
        if A.is_commutative():
            return self._matrix
        return sum([self.vector()[i] * A.left_table()[i] for
                    i in xrange(A.degree())])

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: A(1)
            e0
        """
        s = " "
        A = self.parent()
        m = A.degree()
        coeffs = list(self.vector())
        atomic = A.base_ring()._repr_option('element_is_atomic')
        non_zero = False
        for n in xrange(m):
            x = coeffs[n]
            if x:
                if non_zero:
                    s += " + "
                non_zero = True
                x = y = repr(x)
                if y.find('-') == 0:
                    y = y[1:]
                if not atomic and (y.find("+") != -1 or y.find("-") != -1):
                    x = "({})".format(x)
                var = "*{}".format(A._names[n])
                s += "{}{}".format(x, var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*',' ', s)
        s = re.sub(r' -1(\.0+)?\*',' -', s)
        if s == " ":
            return "0"
        return s[1:]

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: latex(A(1))  # indirect doctest
            \left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)
        """
        from sage.misc.latex import latex
        return latex(self.matrix())

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: A(2) == 2
            True
            sage: A(2) == 3
            False
            sage: A(2) == GF(5)(2)
            False
        """
        A = self.parent()
        return (A.has_coerce_map_from(parent(other)) and
                self.vector() == A(other).vector())

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(1) != 0
            True
        """
        return not self.__eq__(other)

    def __gt__(self, other):
        """
        Raise a ``TypeError`` as there is no (natural) ordering defined on a
        finite-dimensional algebra::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: A(1) > 0
            Traceback (most recent call last):
            ...
            TypeError: there is no ordering defined on a finite-dimensional algebra
            sage: A(1) < 0
            Traceback (most recent call last):
            ...
            TypeError: there is no ordering defined on a finite-dimensional algebra
            sage: A(1) >= 0
            Traceback (most recent call last):
            ...
            TypeError: there is no ordering defined on a finite-dimensional algebra
            sage: A(1) <= 0
            Traceback (most recent call last):
            ...
            TypeError: there is no ordering defined on a finite-dimensional algebra
        """
        raise TypeError("there is no ordering defined on a finite-dimensional algebra")

    __lt__ = __ge__ = __le__ = __gt__

    def _add_(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: A.basis()[0] + A.basis()[1]
            e0 + e1
        """
        return self.__class__(self.parent(), self._vector + other._vector)

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: A.basis()[0] - A.basis()[1]
            e0 + 2*e1
        """
        return self.__class__(self.parent(), self._vector - other._vector)

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: C.basis()[1] * C.basis()[2]
            e1
        """
        return self.__class__(self.parent(), self._vector * other._matrix)

    def _lmul_(self, other):
        """
        TESTS::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: c = C.random_element()
            sage: c * 2 == c + c
            True
        """
        if not self.parent().base_ring().has_coerce_map_from(other.parent()):
            raise TypeError("unsupported operand parent(s) for '*': '{}' and '{}'"
                            .format(self.parent(), other.parent()))
        return self.__class__(self.parent(), self._vector * other)

    def _rmul_(self, other):
        """
        TESTS::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: c = C.random_element()
            sage: 2 * c == c + c
            True
        """
        if not self.parent().base_ring().has_coerce_map_from(other.parent()):
            raise TypeError("unsupported operand parent(s) for '*': '{}' and '{}'"
                            .format(self.parent(), other.parent()))
        return FiniteDimensionalAlgebraElement(self.parent(), other * self._vector) # Note the different order

    def __pow__(self, n):
        """
        Return ``self`` raised to the power ``n``.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: b = B(vector(QQ, [2,3,4]))
            sage: b^6
            64*e0 + 576*e1 + 4096*e2
        """
        A = self.parent()
        if not (A._assume_associative or A.is_associative()):
            raise TypeError("algebra is not associative")
        if n > 0:
            return self.__class__(A, self.vector() * self._matrix.__pow__(n - 1))
        if not A.is_unitary():
            raise TypeError("algebra is not unitary")
        if n == 0:
            return A.one()
        a = self.inverse()
        return self.__class__(A, a.vector() * a.matrix().__pow__(-n - 1))

    def is_invertible(self):
        """
        Return ``True`` if ``self`` has a two-sided multiplicative
        inverse.

        .. NOTE::

            If an element of a unitary finite-dimensional algebra over a field
            admits a left inverse, then this is the unique left
            inverse, and it is also a right inverse.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [-1,0]])])
            sage: C([1,2]).is_invertible()
            True
            sage: C(0).is_invertible()
            False
        """
        return self._inverse is not None

    @lazy_attribute
    def _inverse(self):
        """
        The two-sided inverse of ``self``, if it exists; otherwise this
        is ``None``.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [-1,0]])])
            sage: C([1,2])._inverse
            1/5*e0 - 2/5*e1
            sage: C(0)._inverse is None
            True
        """
        A = self.parent()
        if not A.is_unitary():
            return None

        try:
            a = self.matrix().inverse()
            y = FiniteDimensionalAlgebraElement(A, a, check=True)
            y._inverse = self
            return y
        except (ZeroDivisionError, ValueError):
            return None

    def inverse(self):
        """
        Return the two-sided multiplicative inverse of ``self``, if it
        exists.

        .. NOTE::

            If an element of a unitary finite-dimensional algebra over a field
            admits a left inverse, then this is the unique left
            inverse, and it is also a right inverse.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [-1,0]])])
            sage: C([1,2]).inverse()
            1/5*e0 - 2/5*e1
        """
        A = self.parent()
        if not A.is_unitary():
            raise TypeError("algebra is not unitary")

        if self._inverse is None:
            raise ZeroDivisionError("element is not invertible")
        return self._inverse

    def is_zerodivisor(self):
        """
        Return ``True`` if ``self`` is a left or right zero-divisor.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: C([1,0]).is_zerodivisor()
            False
            sage: C([0,1]).is_zerodivisor()
            True
        """
        return self.matrix().det() == 0 or self.left_matrix().det() == 0

    def is_nilpotent(self):
        """
        Return ``True`` if ``self`` is nilpotent.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [0,0]])])
            sage: C([1,0]).is_nilpotent()
            False
            sage: C([0,1]).is_nilpotent()
            True

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([0])])
            sage: A([1]).is_nilpotent()
            True
        """
        A = self.parent()
        if not (A._assume_associative or A.is_associative()):
            raise TypeError("algebra is not associative")
        return self.matrix().__pow__(A.degree()) == 0

    def minimal_polynomial(self):
        """
        Return the minimal polynomial of ``self``.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(0).minimal_polynomial()
            x
            sage: b = B.random_element()
            sage: f = b.minimal_polynomial(); f  # random
            x^3 + 1/2*x^2 - 7/16*x + 1/16
            sage: f(b) == 0
            True
        """
        A = self.parent()
        if not A.is_unitary():
            raise TypeError("algebra is not unitary")
        if not (A._assume_associative or A.is_associative()):
            raise TypeError("algebra is not associative")
        return self.matrix().minimal_polynomial()

    def characteristic_polynomial(self):
        """
        Return the characteristic polynomial of ``self``.

        .. NOTE::

            This function just returns the characteristic polynomial
            of the matrix of right multiplication by ``self``.  This
            may not be a very meaningful invariant if the algebra is
            not unitary and associative.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(0).characteristic_polynomial()
            x^3
            sage: b = B.random_element()
            sage: f = b.characteristic_polynomial(); f  # random
            x^3 - 8*x^2 + 16*x
            sage: f(b) == 0
            True
        """
        return self.matrix().characteristic_polynomial()

