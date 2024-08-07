"""
Elements of Finite Algebras
"""
# ****************************************************************************
#  Copyright (C) 2011 Johan Bosman <johan.g.bosman@gmail.com>
#  Copyright (C) 2011, 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#  Copyright (C) 2011 Michiel Kosters <kosters@gmail.com>
#  Copyright (C) 2017 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import re

from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer import Integer

from cpython.object cimport PyObject_RichCompare as richcmp

cpdef FiniteDimensionalAlgebraElement unpickle_FiniteDimensionalAlgebraElement(A, vec, mat):
    """
    Helper for unpickling of finite dimensional algebra elements.

    TESTS::

        sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
        ....:                                   Matrix([[1,1,0], [0,1,1], [0,1,1]]),
        ....:                                   Matrix([[0,0,1], [0,1,0], [1,0,0]])])
        sage: x = B([1,2,3])
        sage: loads(dumps(x)) == x      # indirect doctest
        True
    """
    cdef FiniteDimensionalAlgebraElement x = A.element_class.__new__(A.element_class)
    AlgebraElement.__init__(x, A)
    x._vector = vec
    x.__matrix = mat
    return x

cdef class FiniteDimensionalAlgebraElement(AlgebraElement):
    r"""
    Create an element of a :class:`FiniteDimensionalAlgebra` using a multiplication table.

    INPUT:

    - ``A`` -- a :class:`FiniteDimensionalAlgebra` which will be the parent

    - ``elt`` -- vector, matrix or element of the base field
      (default: ``None``)

    - ``check`` -- boolean (default: ``True``); if ``False`` and ``elt`` is a
      matrix, assume that it is known to be the matrix of an element

    If ``elt`` is a vector or a matrix consisting of a single row, it is
    interpreted as a vector of coordinates with respect to the given basis
    of ``A``.  If ``elt`` is a square matrix, it is interpreted as a
    multiplication matrix with respect to this basis.

    EXAMPLES::

        sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
        ....:                                      Matrix([[0,1], [0,0]])])
        sage: A(17)
        2*e0
        sage: A([1,1])
        e0 + e1
    """
    def __init__(self, A, elt=None, check=True):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
            ....:                                      Matrix([[0,1], [0,0]])])
            sage: A(QQ(4))
            Traceback (most recent call last):
            ...
            TypeError: elt should be a vector, a matrix, or an element of the base field

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [-1,0]])])
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
            self._vector = MatrixSpace(k, 1, n)()
            self.__matrix = MatrixSpace(k, n)()
        else:
            if isinstance(elt, int):
                elt = Integer(elt)
            elif isinstance(elt, list):
                elt = MatrixSpace(k, 1, n)(elt)
            if A == elt.parent():
                mat = (<FiniteDimensionalAlgebraElement> elt).__matrix
                if mat is None:
                    self.__matrix = None
                else:
                    self.__matrix = mat.base_extend(k)
                self._vector = elt._vector.base_extend(k)
            elif k.has_coerce_map_from(elt.parent()):
                e = k(elt)
                if e == 0:
                    self._vector = MatrixSpace(k, 1, n)()
                    self.__matrix = MatrixSpace(k, n)()
                elif A.is_unitary():
                    self._vector = A._one * e
                    self.__matrix = MatrixSpace(k, n)(1) * e
                else:
                    raise TypeError("algebra is not unitary")
            elif isinstance(elt, Vector):
                self._vector = MatrixSpace(k, 1, n)(list(elt))
            elif isinstance(elt, Matrix):
                if elt.ncols() != n:
                    raise ValueError("matrix does not define an element of the algebra")
                if elt.nrows() == 1:
                    self._vector = elt.__copy__()
                else:
                    if not A.is_unitary():
                        raise TypeError("algebra is not unitary")
                    self._vector = A._one * elt
                    if check and self._matrix != elt:
                        raise ValueError("matrix does not define an element of the algebra")
            else:
                raise TypeError("elt should be a vector, a matrix, " +
                                "or an element of the base field")

    def __reduce__(self):
        """
        TESTS::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[1,1,0], [0,1,1], [0,1,1]]),
            ....:                                   Matrix([[0,0,1], [0,1,0], [1,0,0]])])
            sage: x = B([1,2,3])
            sage: loads(dumps(x)) == x      # indirect doctest
            True
            sage: loads(dumps(x)) is x
            False
        """
        return unpickle_FiniteDimensionalAlgebraElement, (self._parent, self._vector, self.__matrix)

    def __setstate__(self, state):
        """
        This method serves at unpickling old pickles.

        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: x = A.element_class.__new__(A.element_class)
            sage: x.__setstate__((A, {'_vector':vector([1,1,1]), '_matrix':matrix(QQ,3,[1,1,0, 0,1,0, 0,0,1])}))
            sage: x
            e0 + e1 + e2
            sage: x.matrix()
            [1 1 0]
            [0 1 0]
            [0 0 1]

        Note that in old pickles, the vector actually is a vector. However,
        it is converted into a single-row matrix, in the new implementation::

            sage: x.vector()
            (1, 1, 1)
        """
        self._parent, D = state
        v = D.pop('_vector')
        if isinstance(v, Vector):
            self._vector = MatrixSpace(self._parent.base_ring(), 1, len(v))(list(v))
        else:
            self._vector = v
        self.__matrix = D.pop('_matrix', None)
        try:
            self.__dict__ = D
        except AttributeError:
            pass

    @property
    def _matrix(self):
        """
        TESTS::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[1,1,0], [0,1,1], [0,1,1]]),
            ....:                                   Matrix([[0,0,1], [0,1,0], [1,0,0]])])
            sage: x = B([1,2,3])
            sage: x._matrix
            [3 2 3]
            [0 6 2]
            [3 2 2]
        """
        cdef Py_ssize_t i
        cdef tuple table
        if self.__matrix is None:
            A = self.parent()
            table = <tuple> A.table()
            ret = sum(self._vector[0, i] * table[i] for i in range(A.degree()))
            self.__matrix = MatrixSpace(A.base_ring(), A.degree())(ret)
        return self.__matrix

    def vector(self):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(5).vector()
            (5, 0, 5)
        """
        # By :issue:`23707`, ``self._vector`` now is a single row matrix,
        # not a vector, which results in a speed-up.
        # For backwards compatibility, this method still returns a vector.
        return self._vector[0]

    def matrix(self):
        """
        Return the matrix for multiplication by ``self`` from the right.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(5).matrix()
            [5 0 0]
            [0 5 0]
            [0 0 5]
        """
        return self._matrix

    def monomial_coefficients(self, copy=True):
        """
        Return a dictionary whose keys are indices of basis elements in
        the support of ``self`` and whose values are the corresponding
        coefficients.

        INPUT:

        - ``copy`` -- ignored

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [-1,0]])])
            sage: elt = B(Matrix([[1,1], [-1,1]]))
            sage: elt.monomial_coefficients()
            {0: 1, 1: 1}
        """
        cdef Py_ssize_t i
        return {i: self._vector[0, i] for i in range(self._vector.ncols())}

    def left_matrix(self):
        """
        Return the matrix for multiplication by ``self`` from the left.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: C([1,2,0]).left_matrix()
            [1 0 0]
            [0 1 0]
            [0 2 0]
        """
        A = self.parent()
        if A.is_commutative():
            return self._matrix
        return sum([self._vector[0, i] * A.left_table()[i] for
                    i in range(A.degree())])

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
            ....:                                      Matrix([[0,1], [0,0]])])
            sage: A(1)
            e0
        """
        s = " "
        A = self.parent()
        m = A.degree()
        coeffs = self._vector.list()
        atomic = A.base_ring()._repr_option('element_is_atomic')
        non_zero = False
        for n in range(m):
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
        s = re.sub(r' 1(\.0+)?\*', ' ', s)
        s = re.sub(r' -1(\.0+)?\*', ' -', s)
        if s == " ":
            return "0"
        return s[1:]

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
            ....:                                      Matrix([[0,1], [0,0]])])
            sage: latex(A(1))  # indirect doctest
            \left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)
        """
        from sage.misc.latex import latex
        return latex(self.matrix())

    def __getitem__(self, m):
        """
        Return the `m`-th coefficient of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: A([2,1/4,3])[2]
            3
        """
        return self._vector[0, m]

    def __len__(self):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: len(A([2,1/4,3]))
            3
        """
        return self._vector.ncols()

    # (Rich) comparison
    cpdef _richcmp_(self, right, int op):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
            ....:                                      Matrix([[0,1], [0,0]])])
            sage: A(2) == 2
            True
            sage: A(2) == 3
            False
            sage: A(2) == GF(5)(2)
            False

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(1) != 0
            True

        By :issue:`23707`, an ordering is defined on finite-dimensional algebras, corresponding
        to the ordering of the defining vectors; this may be handy if the vector space basis of
        the algebra corresponds to the standard monomials of the relation ideal, when
        the algebra is considered as a quotient of a path algebra. ::

            sage: A(1) > 0
            True
            sage: A(1) < 0
            False
            sage: A(1) >= 0
            True
            sage: A(1) <= 0
            False
        """
        return richcmp(self._vector, <FiniteDimensionalAlgebraElement>right._vector, op)

    cpdef _add_(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
            ....:                                      Matrix([[0,1], [0,0]])])
            sage: A.basis()[0] + A.basis()[1]
            e0 + e1
        """
        return self._parent.element_class(self._parent, self._vector + <FiniteDimensionalAlgebraElement>other._vector)

    cpdef _sub_(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1,0], [0,1]]),
            ....:                                      Matrix([[0,1], [0,0]])])
            sage: A.basis()[0] - A.basis()[1]
            e0 + 2*e1
        """
        return self._parent.element_class(self._parent, self._vector - <FiniteDimensionalAlgebraElement>other._vector)

    cpdef _mul_(self, other):
        """
        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: C.basis()[1] * C.basis()[2]
            e1
        """
        return self._parent.element_class(self._parent, self._vector * <FiniteDimensionalAlgebraElement>(other)._matrix)

    cpdef _lmul_(self, Element other):
        """
        TESTS::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: c = C.random_element()
            sage: c * 2 == c + c
            True
        """
        if not self._parent.base_ring().has_coerce_map_from(other.parent()):
            raise TypeError("unsupported operand parent(s) for *: '{}' and '{}'"
                            .format(self.parent(), other.parent()))
        return self._parent.element_class(self._parent, self._vector * other)

    cpdef _rmul_(self, Element other):
        """
        TESTS::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: c = C.random_element()
            sage: 2 * c == c + c
            True
        """
        if not self._parent.base_ring().has_coerce_map_from(other.parent()):
            raise TypeError("unsupported operand parent(s) for *: '{}' and '{}'"
                            .format(self.parent(), other.parent()))
        # Note the different order below
        return self._parent.element_class(self._parent, other * self._vector)

    def __pow__(self, n, m):
        """
        Return ``self`` raised to the power ``n``.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: b = B([2,3,4])
            sage: b^6
            64*e0 + 576*e1 + 4096*e2
        """
        A = self.parent()
        if not (A._assume_associative or A.is_associative()):
            raise TypeError("algebra is not associative")
        if n > 0:
            return A.element_class(A, self._vector * self._matrix ** (n - 1))
        if not A.is_unitary():
            raise TypeError("algebra is not unitary")
        if n == 0:
            return A.one()
        cdef FiniteDimensionalAlgebraElement a = <FiniteDimensionalAlgebraElement>(~self)
        return A.element_class(A, a._vector * a.__matrix ** (-n - 1))

    def __invert__(self):
        """
        TESTS::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [-1,0]])])
            sage: x = C([1,2])
            sage: y = ~x; y                 # indirect doctest
            1/5*e0 - 2/5*e1
            sage: x*y
            e0
            sage: C.one()
            e0
        """
        return self.inverse()

    def is_invertible(self):
        """
        Return ``True`` if ``self`` has a two-sided multiplicative
        inverse.

        This assumes that the algebra to which ``self`` belongs is
        associative.

        .. NOTE::

            If an element of a unitary finite-dimensional algebra over a field
            admits a left inverse, then this is the unique left
            inverse, and it is also a right inverse.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [-1,0]])])
            sage: C([1,2]).is_invertible()
            True
            sage: C(0).is_invertible()
            False
        """
        return self._inverse is not None

    @property
    def _inverse(self):
        """
        The two-sided inverse of ``self``, if it exists; otherwise this
        is ``None``.

        This assumes that the algebra to which ``self`` belongs is
        associative.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [-1,0]])])
            sage: C([1,2])._inverse
            1/5*e0 - 2/5*e1
            sage: C(0)._inverse is None
            True
        """
        cdef FiniteDimensionalAlgebraElement y
        if self.__inverse is None:
            A = self.parent()
            if not A.is_unitary():
                self.__inverse = False
            try:
                a = self._matrix.inverse()
                y = FiniteDimensionalAlgebraElement(A, a, check=True)
                y.__inverse = self
                self.__inverse = y
            except (ZeroDivisionError, ValueError):
                self.__inverse = False
        if self.__inverse is False:
            return None
        return self.__inverse

    def inverse(self):
        """
        Return the two-sided multiplicative inverse of ``self``, if it
        exists.

        This assumes that the algebra to which ``self`` belongs is
        associative.

        .. NOTE::

            If an element of a finite-dimensional unitary associative
            algebra over a field admits a left inverse, then this is the
            unique left inverse, and it is also a right inverse.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [-1,0]])])
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

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [0,0]])])
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

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]),
            ....:                                   Matrix([[0,1], [0,0]])])
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
        return self.matrix() ** A.degree() == 0

    def minimal_polynomial(self):
        """
        Return the minimal polynomial of ``self``.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(0).minimal_polynomial()                                             # needs sage.libs.pari
            x
            sage: b = B.random_element()
            sage: f = b.minimal_polynomial(); f  # random                               # needs sage.libs.pari
            x^3 + 1/2*x^2 - 7/16*x + 1/16
            sage: f(b) == 0                                                             # needs sage.libs.pari
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

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]),
            ....:                                   Matrix([[0,1,0], [0,0,0], [0,0,0]]),
            ....:                                   Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(0).characteristic_polynomial()                                      # needs sage.libs.pari
            x^3
            sage: b = B.random_element()
            sage: f = b.characteristic_polynomial(); f  # random                        # needs sage.libs.pari
            x^3 - 8*x^2 + 16*x
            sage: f(b) == 0                                                             # needs sage.libs.pari
            True
        """
        return self.matrix().characteristic_polynomial()
