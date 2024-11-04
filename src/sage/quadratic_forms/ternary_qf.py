r"""
Ternary quadratic form with integer coefficients

AUTHOR:

- Gustavo Rama

Based in code of Gonzalo Tornaria

The form `a\cdot x^2 + b\cdot y^2 + c\cdot z^2 + r\cdot yz + s\cdot xz + t\cdot xy` is stored as a tuple ``(a, b, c, r, s, t)`` of integers.
"""

# ****************************************************************************
#       Copyright (C) 2012 Gustavo Rama
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import gcd, kronecker as kronecker_symbol
from sage.categories.rings import Rings
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.prandom import randint
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.quadratic_forms.ternary import (_basic_lemma,
                                          _find_a_ternary_qf_by_level_disc,
                                          _find_all_ternary_qf_by_level_disc,
                                          _find_p_neighbor_from_vec,
                                          _find_zeros_mod_p_2,
                                          _find_zeros_mod_p_odd,
                                          _reduced_ternary_form_eisenstein_with_matrix,
                                          _reduced_ternary_form_eisenstein_without_matrix)

from sage.rings.finite_rings.integer_mod import mod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygens
from sage.structure.element import Vector, Matrix
from sage.structure.sage_object import SageObject


class TernaryQF(SageObject):
    r"""
    The ``TernaryQF`` class represents a quadratic form in 3 variables with coefficients in `\ZZ`.

    INPUT:

    - ``v`` -- list or tuple of 6 entries:  ``[a,b,c,r,s,t]``

    OUTPUT: the ternary quadratic form `a\cdot x^2 + b\cdot y^2 + c\cdot z^2 + r\cdot y\cdot z + s\cdot x\cdot z + t\cdot x\cdot y`

    EXAMPLES::

        sage: Q = TernaryQF([1, 2, 3, 4, 5, 6]); Q
        Ternary quadratic form with integer coefficients:
        [1 2 3]
        [4 5 6]
        sage: A = matrix(ZZ, 3, [1, -7, 1, 0, -2, 1, 0, -1, 0])
        sage: Q(A)
        Ternary quadratic form with integer coefficients:
        [1 187 9]
        [-85 8 -31]
        sage: TestSuite(TernaryQF).run()
    """

    __slots__ = ['_a', '_b', '_c', '_r', '_s', '_t', '_automorphisms', '_number_of_automorphisms']

    possible_automorphisms = None

    def __init__(self, v):
        r"""
        Create the ternary quadratic form `a\cdot x^2 + b\cdot y^2 + c\cdot z^2 + r\cdot y\cdot z + s\cdot x\cdot z + t\cdot x\cdot y` from the
        tuple ``v=[a,b,c,r,s,t]`` over `\ZZ`.

        INPUT:

        - ``v`` -- 6-tuple of integers

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 4, 5, 6]); Q
            Ternary quadratic form with integer coefficients:
            [1 2 3]
            [4 5 6]
        """

        if len(v) != 6:
            # Check we have six coefficients
            raise ValueError("Ternary quadratic form must be given by a list of six coefficients")
        self._a, self._b, self._c, self._r, self._s, self._t = (ZZ(x) for x in v)
        self._automorphisms = None
        self._number_of_automorphisms = None

    def coefficients(self):
        r"""
        Return the list of coefficients of the ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 4, 5, 6]); Q
            Ternary quadratic form with integer coefficients:
            [1 2 3]
            [4 5 6]
            sage: Q.coefficients()
            (1, 2, 3, 4, 5, 6)
        """
        return self._a, self._b, self._c, self._r, self._s, self._t

    def __hash__(self):
        """
        Return a hash for ``self``.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 4, 5, 6])
            sage: Q.__hash__()
            5881802312257552497  # 64-bit
            1770036893           # 32-bit
        """

        return hash(self.coefficients())

    def coefficient(self, n):
        r"""
        Return the `n`-th coefficient of the ternary quadratic form.

        INPUT:

        - ``n`` -- integer with `0 \leq n \leq 5`

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 4, 5, 6]); Q
            Ternary quadratic form with integer coefficients:
            [1 2 3]
            [4 5 6]
            sage: Q.coefficient(2)
            3
            sage: Q.coefficient(5)
            6
        """
        return self.coefficients()[n]

    def polynomial(self, names='x,y,z'):
        r"""
        Return the polynomial associated to the ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 0, 2, -3, -1]); Q
            Ternary quadratic form with integer coefficients:
            [1 1 0]
            [2 -3 -1]
            sage: p = Q.polynomial(); p
            x^2 - x*y + y^2 - 3*x*z + 2*y*z
            sage: p.parent()
            Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        x, y, z = polygens(ZZ, names)
        return self._a * x**2 + self._b * y**2 + self._c * z**2 + self._t * x*y + self._s * x*z + self._r * y*z

    def _repr_(self):
        r"""
        Display the quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 0, 2, -3, -1])
            sage: print(Q._repr_())
            Ternary quadratic form with integer coefficients:
            [1 1 0]
            [2 -3 -1]
            sage: Q = TernaryQF([0, 0, 0, 0, 0, 0]); Q
            Ternary quadratic form with integer coefficients:
            [0 0 0]
            [0 0 0]
        """
        rep = 'Ternary quadratic form with integer coefficients:\n'
        rep += '[' + str(self._a) + ' ' + str(self._b) + ' ' + str(self._c) + ']\n'
        rep += '[' + str(self._r) + ' ' + str(self._s) + ' ' + str(self._t) + ']'
        return rep

    def __call__(self, v):
        r"""
        Evaluate this ternary quadratic form `Q` on a vector of 3 elements,
        or matrix of elements in Z, with 3 rows.

        OUTPUT:

        If a vector is given, then the output will be an integer `Q(v)`,
        but if a matrix is given, the output will be a ternary quadratic form
        if the matrix has 3 columns, or a quadratic form if not.
        The quadratic form in matrix notation will be:

        .. MATH::

            Q' = v^t\cdot Q\cdot v.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 1, -1, -2, -3])
            sage: Q((1, 1, 1))
            -3
            sage: M = matrix(ZZ, 3, 2, [358, 6, 2, 0, 0, 4])
            sage: Q(M)
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 126020 1388 ]
            [ * 4 ]
            sage: M = matrix(ZZ, 3, 3, [1, 3, 0, -1, 4, 2, 1, -1, -1])
            sage: M
            [ 1  3  0]
            [-1  4  2]
            [ 1 -1 -1]
            sage: Q(M)
            Ternary quadratic form with integer coefficients:
            [5 0 7]
            [12 -13 -16]
        """
        if isinstance(v, Matrix):
            # Check that v has 3 rows
            if v.nrows() != 3:
                raise TypeError("the matrix must have 3 rows")
            # Check if v has 3 cols
            if v.ncols() == 3:
                M = v.transpose() * self.matrix() * v
                return TernaryQF([M[0, 0]//2, M[1, 1]//2, M[2, 2]//2,
                                  M[1, 2], M[0, 2], M[0, 1]])
            else:
                return QuadraticForm(ZZ, v.transpose() * self.matrix() * v)
        elif isinstance(v, (Vector, list, tuple)):
            # Check that v has length 3
            if len(v) != 3:
                raise TypeError("your vector needs to have length 3")
            v0, v1, v2 = v
            a, b, c, r, s, t = self.coefficients()
            return a*v0**2 + b*v1**2 + c*v2**2 + r*v1*v2 + s*v0*v2 + t*v0*v1
        else:
            raise TypeError("presently we can only evaluate a quadratic form on a list, tuple, vector or matrix")

    def quadratic_form(self):
        r"""
        Return a :class:`QuadraticForm` with the same coefficients as ``self`` over `\ZZ`.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 1, 1, 1])
            sage: QF1 = Q.quadratic_form(); QF1
            Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 1 1 1 ]
            [ * 2 1 ]
            [ * * 3 ]
            sage: QF2 = QuadraticForm(ZZ, 3, [1, 1, 1, 2, 1, 3])
            sage: bool(QF1 == QF2)
            True
        """
        return QuadraticForm(ZZ, 3, [self._a, self._t, self._s, self._b, self._r, self._c])

    def matrix(self):
        r"""
        Return the Hessian matrix associated to the ternary quadratic form.
        That is, if `Q` is a ternary quadratic form, `Q(x,y,z) = a\cdot x^2 + b\cdot y^2 + c\cdot z^2 + r\cdot y\cdot z + s\cdot x\cdot z + t\cdot x\cdot y`,
        then the Hessian matrix associated to `Q` is
        ::

            [2\cdot a t s]
            [t 2\cdot b r]
            [s r 2\cdot c]

        EXAMPLES::

            sage: Q = TernaryQF([1,1,2,0,-1,4]); Q
            Ternary quadratic form with integer coefficients:
            [1 1 2]
            [0 -1 4]
            sage: M = Q.matrix(); M
            [ 2  4 -1]
            [ 4  2  0]
            [-1  0  4]
            sage: v = vector((1, 2, 3))
            sage: Q(v)
            28
            sage: (v*M*v.column())[0]//2
            28
        """
        M = matrix(ZZ, 3, [2*self._a, self._t, self._s, self._t, 2*self._b, self._r, self._s, self._r, 2*self._c])
        return M

    def disc(self):
        r"""
        Return the discriminant of the ternary quadratic form, this is the determinant of the matrix divided by 2.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 2, 0, -1, 4])
            sage: Q.disc()
            -25
            sage: Q.matrix().det()
            -50
        """
        return 4*self._a*self._b*self._c + self._r*self._s*self._t - self._a*self._r**2 - self._b*self._s**2 - self._c*self._t**2

    def is_definite(self) -> bool:
        """
        Determine if the ternary quadratic form is definite.

        EXAMPLES::

            sage: Q = TernaryQF([10, 10, 1, -1, 2, 3])
            sage: Q.is_definite()
            True
            sage: (-Q).is_definite()
            True
            sage: Q = TernaryQF([1, 1, 2, -3, 0, -1])
            sage: Q.is_definite()
            False
        """
        d1 = self._a
        if d1 == 0:
            return False
        d2 = 4*self._a*self._b-self._t**2
        if d2 == 0:
            return False
        d3 = self.disc()
        if d3 == 0:
            return False
        if d1 > 0:
            if d2 > 0:
                if d3 > 0:
                    return True
                else:
                    return False
            else:
                return False
        else:
            if d2 > 0:
                if d3 < 0:
                    return True
                else:
                    return False
            else:
                return False

    def is_positive_definite(self) -> bool:
        """
        Determine if the ternary quadratic form is positive definite.

        EXAMPLES::

            sage: Q = TernaryQF([10, 10, 1, -1, 2, 3])
            sage: Q.is_positive_definite()
            True
            sage: (-Q).is_positive_definite()
            False
            sage: Q = TernaryQF([1, 1, 0, 0, 0, 0])
            sage: Q.is_positive_definite()
            False
            sage: Q = TernaryQF([1, 1, 1, -1, -2, -3])
            sage: Q((1,1,1))
            -3
            sage: Q.is_positive_definite()
            False
        """
        d1 = self._a
        if d1 == 0:
            return False
        d2 = 4*self._a*self._b-self._t**2
        if d2 == 0:
            return False
        d3 = self.disc()
        if d3 == 0:
            return False
        if d1 > 0:
            if d2 > 0:
                if d3 > 0:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def is_negative_definite(self) -> bool:
        """
        Determine if the ternary quadratic form is negative definite.

        EXAMPLES::

            sage: Q = TernaryQF([-8, -9, -10, 1, 9, -3])
            sage: Q.is_negative_definite()
            True
            sage: Q = TernaryQF([-4, -1, 6, -5, 1, -5])
            sage: Q((0, 0, 1))
            6
            sage: Q.is_negative_definite()
            False
        """
        d1 = self._a
        if d1 == 0:
            return False
        d2 = 4*self._a*self._b-self._t**2
        if d2 == 0:
            return False
        d3 = self.disc()
        if d3 == 0:
            return False
        if d1 < 0:
            if d2 > 0:
                if d3 < 0:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def __neg__(self):
        """
        Return the ternary quadratic form with coefficients negatives of ``self``.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 2, -2, 0, -1]); Q
            Ternary quadratic form with integer coefficients:
            [1 1 2]
            [-2 0 -1]
            sage: -Q
            Ternary quadratic form with integer coefficients:
            [-1 -1 -2]
            [2 0 1]
            sage: Q = TernaryQF([0, 0, 0, 0, 0, 0])
            sage: Q == -Q
            True
        """
        return TernaryQF([-a for a in self.coefficients()])

    def is_primitive(self) -> bool:
        """
        Determine if the ternary quadratic form is primitive.

        This means that the greatest common divisor of the coefficients
        of the form is 1.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 4, 5, 6])
            sage: Q.is_primitive()
            True
            sage: Q.content()
            1
            sage: Q = TernaryQF([10, 10, 10, 5, 5, 5])
            sage: Q.content()
            5
            sage: Q.is_primitive()
            False
        """
        return self.content() == 1

    def primitive(self):
        """
        Return the primitive version of the ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([2, 2, 2, 1, 1, 1])
            sage: Q.is_primitive()
            True
            sage: Q.primitive()
            Ternary quadratic form with integer coefficients:
            [2 2 2]
            [1 1 1]
            sage: Q.primitive() == Q
            True
            sage: Q = TernaryQF([10, 10, 10, 5, 5, 5])
            sage: Q.primitive()
            Ternary quadratic form with integer coefficients:
            [2 2 2]
            [1 1 1]
        """
        l = self.coefficients()
        g = gcd(l)
        return TernaryQF([a//g for a in l])

    def scale_by_factor(self, k):
        """
        Scale the values of the ternary quadratic form by the number ``k``.

        OUTPUT:

        If ``k`` times the content of the ternary quadratic form is an integer, return a ternary quadratic form;
        otherwise, return a quadratic form of dimension 3.

        EXAMPLES::

            sage: Q = TernaryQF([2, 2, 4, 0, -2, 8])
            sage: Q
            Ternary quadratic form with integer coefficients:
            [2 2 4]
            [0 -2 8]
            sage: Q.scale_by_factor(5)
            Ternary quadratic form with integer coefficients:
            [10 10 20]
            [0 -10 40]
            sage: Q.scale_by_factor(1/2)
            Ternary quadratic form with integer coefficients:
            [1 1 2]
            [0 -1 4]
            sage: Q.scale_by_factor(1/3)
            Quadratic form in 3 variables over Rational Field with coefficients:
            [ 2/3 8/3 -2/3 ]
            [ * 2/3 0 ]
            [ * * 4/3 ]
        """
        if k * self.content() in ZZ:
            return TernaryQF([ZZ(k*self._a), ZZ(k*self._b), ZZ(k*self._c),
                              ZZ(k*self._r), ZZ(k*self._s), ZZ(k*self._t)])

        R = k.parent()
        if R not in Rings():
            raise TypeError(f"{k} does not belong to a ring")

        return QuadraticForm(R, 3, [k * self._a, k * self._t, k * self._s,
                                    k * self._b, k * self._r, k * self._c])

    def reciprocal(self):
        """
        Return the reciprocal quadratic form associated to the given form.

        This is defined as the multiple of the primitive adjoint with the same
        content as the given form.

        EXAMPLES::

            sage: Q = TernaryQF([2, 2, 14, 0, 0, 0])
            sage: Q.reciprocal()
            Ternary quadratic form with integer coefficients:
            [14 14 2]
            [0 0 0]
            sage: Q.content()
            2
            sage: Q.reciprocal().content()
            2
            sage: Q.adjoint().content()
            16
        """
        return self.adjoint().primitive().scale_by_factor(self.content())

    def reciprocal_reduced(self):
        """
        Return the reduced form of the reciprocal form of the given ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 3, 0, -1, 0])
            sage: Qrr = Q.reciprocal_reduced(); Qrr
            Ternary quadratic form with integer coefficients:
            [4 11 12]
            [0 -4 0]
            sage: Q.is_eisenstein_reduced()
            True
            sage: Qr = Q.reciprocal()
            sage: Qr.reduced_form_eisenstein(matrix=False) == Qrr
            True
        """
        return self.reciprocal().reduced_form_eisenstein(matrix=False)

    def divisor(self):
        """
        Return the content of the adjoint form associated to the given form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 17, 0, 0, 0])
            sage: Q.divisor()
            4
        """
        A11 = 4*self._b*self._c - self._r**2
        A22 = 4*self._a*self._c - self._s**2
        A33 = 4*self._a*self._b - self._t**2
        A23 = self._s*self._t - 2*self._a*self._r
        A13 = self._r*self._t - 2*self._b*self._s
        A12 = self._r*self._s - 2*self._c*self._t
        m = gcd([A11, A22, A33, 2*A12, 2*A13, 2*A23])
        return m

    def __eq__(self, right) -> bool:
        """
        Determine if two ternary quadratic forms are equal.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 3, 1, 2, 3])
            sage: Q == Q
            True
            sage: Q1 = TernaryQF([1, 2, 3, 1, 2, 2])
            sage: Q == Q1
            False
        """
        if not isinstance(right, TernaryQF):
            return False
        return self.coefficients() == right.coefficients()

    def adjoint(self):
        """
        Return the adjoint form associated to the given ternary quadratic
        form.

        That is, the Hessian matrix of the adjoint form is twice the
        classical adjoint matrix of the Hessian matrix of the given form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 17, 0, 0, 1])
            sage: Q.adjoint()
            Ternary quadratic form with integer coefficients:
            [68 68 3]
            [0 0 -68]
            sage: Q.adjoint().matrix() == 2*Q.matrix().adjoint_classical()
            True
        """
        A11 = 4*self._b*self._c - self._r**2
        A22 = 4*self._a*self._c - self._s**2
        A33 = 4*self._a*self._b - self._t**2
        A23 = self._s*self._t - 2*self._a*self._r
        A13 = self._r*self._t - 2*self._b*self._s
        A12 = self._r*self._s - 2*self._c*self._t
        return TernaryQF([A11, A22, A33, 2*A23, 2*A13, 2*A12])

    def content(self):
        """
        Return the greatest common divisor of the coefficients of the given ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 2, 0, 0, 0])
            sage: Q.content()
            1
            sage: Q = TernaryQF([2, 4, 6, 0, 0, 0])
            sage: Q.content()
            2
            sage: Q.scale_by_factor(100).content()
            200
        """
        return gcd(self.coefficients())

    def omega(self):
        """
        Return the content of the adjoint of the primitive associated
        ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([4, 11, 12, 0, -4, 0])
            sage: Q.omega()
            176
            sage: Q.primitive().adjoint().content()
            176
        """
        return self.primitive().adjoint().content()

    def delta(self):
        """
        Return the omega of the adjoint of the given ternary quadratic form,
        which is the same as the omega of the reciprocal form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 2, -1, 0, -1])
            sage: Q.delta()
            208
            sage: Q.adjoint().omega()
            208
            sage: Q = TernaryQF([1, -1, 1, 0, 0, 0])
            sage: Q.delta()
            4
            sage: Q.omega()
            4
        """
        return self.adjoint().omega()

    def level(self):
        """
        Return the level of the ternary quadratic form, which is 4 times the discriminant divided by the divisor.

        EXAMPLES::

            sage: Q = TernaryQF([1, 2, 2, -1, 0, -1])
            sage: Q.level()
            52
            sage: 4*Q.disc()/Q.divisor()
            52
        """
        return 4 * self.disc() // self.divisor()

    def is_eisenstein_reduced(self) -> bool:
        r"""
        Determine if the ternary quadratic form is Eisenstein reduced.

        That is, if we have a ternary quadratic form:
        ::

        [a b c]
        [r s t]

        then

        1. `a \leq b \leq c`;
        2. `r`, `s`, and `t` are all positive or all nonpositive;
        3. `a \geq |t|`; `a \geq |s|`; `b \geq |r|`;
        4. `a+b+r+s+t \geq 0`;
        5. `a=t` implies `s \leq 2\cdot r`; `a=s` implies `t \leq 2\cdot r`; `b=r` implies `t \leq 2\cdot s`;
        6. `a=-t` implies `s=0`; `a=-s` implies `t=0`; `b=-r` implies `t=0`;
        7. `a+b+r+s+t = 0` implies `2\cdot a+2\cdot s+t \leq 0`;
        8. `a=b` implies `|r| \leq |s|`; `b=c` implies `|s| \leq |t|`.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 1, 0, 0, 0])
            sage: Q.is_eisenstein_reduced()
            True
            sage: Q = TernaryQF([34, 14, 44, 12, 25, -22])
            sage: Q.is_eisenstein_reduced()
            False
        """
        a, b, c, r, s, t = [self._a, self._b, self._c,
                            self._r, self._s, self._t]

        # cond 2
        if not (r > 0 and t > 0 and s > 0):
            if not (r <= 0 and s <= 0 and t <= 0):
                return False

        # cond 1 & 4
        if not (a <= b <= c and 0 <= a+b+r+s+t):
            return False

        # cond 3
        if not (a >= abs(s) and a >= abs(t) and b >= abs(r)):
            return False

        # cond 8
        if a == b and abs(r) > abs(s):
            return False
        if b == c and abs(s) > abs(t):
            return False
        if a+b+r+s+t == 0 and 2*a+2*s+t > 0:
            return False

        # cond 6
        # r, s, t <= 0
        if r <= 0:
            if a == -t and s != 0:
                return False
            if a == -s and t != 0:
                return False
            if b == -r and t != 0:
                return False

        # cond 7
        # r, s, t > 0
        if a == t and s > 2*r:
            return False
        if a == s and t > 2*r:
            return False
        if b == r and t > 2*s:
            return False

        return True

    def reduced_form_eisenstein(self, matrix=True):
        r"""
        Return the Eisenstein reduced form equivalent to the given positive ternary quadratic form,
        which is unique.

        EXAMPLES::

            sage: Q = TernaryQF([293, 315, 756, 908, 929, 522])
            sage: Qr, m = Q.reduced_form_eisenstein()
            sage: Qr
            Ternary quadratic form with integer coefficients:
            [1 2 2]
            [-1 0 -1]
            sage: Qr.is_eisenstein_reduced()
            True
            sage: m
            [ -54  137  -38]
            [ -23   58  -16]
            [  47 -119   33]
            sage: m.det()
            1
            sage: Q(m) == Qr
            True
            sage: Q = TernaryQF([12,36,3,14,-7,-19])
            sage: Q.reduced_form_eisenstein(matrix = False)
            Ternary quadratic form with integer coefficients:
            [3 8 20]
            [3 2 1]
        """
        if matrix:
            v, M = _reduced_ternary_form_eisenstein_with_matrix(self._a, self._b, self._c, self._r, self._s, self._t)
            return TernaryQF(v), M

        v = _reduced_ternary_form_eisenstein_without_matrix(self._a, self._b, self._c, self._r, self._s, self._t)
        return TernaryQF(v)

    def pseudorandom_primitive_zero_mod_p(self, p):
        """
        Return a tuple of the form `v = (a, b, 1)` such that is a zero of the given ternary quadratic
        positive definite form modulo an odd prime `p`, where `p` doesn't divides the discriminant of the form.

        EXAMPLES::

             sage: Q = TernaryQF([1, 1, 11, 0, -1, 0])
             sage: Q.disc()
             43
             sage: Q.pseudorandom_primitive_zero_mod_p(3)  # random
             (1, 2, 1)
             sage: Q((1, 2, 1))
             15
             sage: v = Q.pseudorandom_primitive_zero_mod_p(1009)                        # needs sage.libs.pari
             sage: Q(v) % 1009                                                          # needs sage.libs.pari
             0
             sage: v[2]                                                                 # needs sage.libs.pari
             1
        """
        a, b, c, r, s, t = self.coefficients()
        while True:

            r1 = randint(0, p-1)
            r2 = randint(0, p-1)
            alpha = (b*r1**2+t*r1+a) % p
            if alpha != 0:

                beta = (2*b*r1*r2+t*r2+r*r1+s) % p
                gamma = (b*r2**2+r*r2+c) % p
                disc = beta**2-4*alpha*gamma
                if mod(disc, p).is_square():

                    z = (-beta+mod(disc, p).sqrt().lift())*(2*alpha).inverse_mod(p)
                    # return vector((z,r1*z+r2,1))%p
                    return z % p, (r1*z+r2) % p, 1

    def find_zeros_mod_p(self, p):
        """
        Find the zeros of the given ternary quadratic positive definite form modulo a prime `p`, where `p` doesn't divide the discriminant of the form.

        EXAMPLES::

            sage: Q = TernaryQF([4, 7, 8, -4, -1, -3])
            sage: Q.is_positive_definite()
            True
            sage: Q.disc().factor()
            3 * 13 * 19
            sage: Q.find_zeros_mod_p(2)
            [(1, 0, 0), (1, 1, 0), (0, 0, 1)]
            sage: zeros_17 = Q.find_zeros_mod_p(17)                                     # needs sage.libs.pari
            sage: len(zeros_17)                                                         # needs sage.libs.pari
            18
            sage: [Q(v)%17 for v in zeros_17]                                           # needs sage.libs.pari
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        if p == 2:
            return _find_zeros_mod_p_2(self._a, self._b, self._c,
                                       self._r, self._s, self._t)

        v = self.pseudorandom_primitive_zero_mod_p(p)
        a, b, c, r, s, t = self.coefficients()
        return _find_zeros_mod_p_odd(a, b, c, r, s, t, p, v)

    def find_p_neighbor_from_vec(self, p, v, mat=False):
        r"""
        Finds the reduced equivalent of the `p`-neighbor of this ternary quadratic form associated to a given
        vector `v` satisfying:

        1. `Q(v) = 0`  (mod `p`)

        2. `v` is a non-singular point of the conic `Q(v) = 0` (mod `p`).

        REFERENCES:

        Gonzalo Tornaria's Thesis, Thrm 3.5, p34.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: Q = TernaryQF([1, 3, 3, -2, 0, -1]); Q
            Ternary quadratic form with integer coefficients:
            [1 3 3]
            [-2 0 -1]
            sage: Q.disc()
            29
            sage: v = (9, 7, 1)
            sage: v in Q.find_zeros_mod_p(11)
            True
            sage: Q11, M = Q.find_p_neighbor_from_vec(11, v, mat=True)
            sage: Q11
            Ternary quadratic form with integer coefficients:
            [1 2 4]
            [-1 -1 0]
            sage: M
            [    -1  -5/11   7/11]
            [     0 -10/11   3/11]
            [     0  -3/11  13/11]
            sage: Q(M) == Q11
            True

        Test that it works with (0, 0, 1)::

            sage: Q.find_p_neighbor_from_vec(3, (0,0,1))                                # needs sage.libs.pari
            Ternary quadratic form with integer coefficients:
            [1 3 3]
            [-2 0 -1]
        """
        if mat:
            q, M = _find_p_neighbor_from_vec(self._a, self._b, self._c, self._r, self._s, self._t, p, v, mat)
            M = matrix(3, M)
            return TernaryQF(q), M*M.det()
        else:
            return TernaryQF(_find_p_neighbor_from_vec(self._a, self._b, self._c, self._r, self._s, self._t, p, v, mat))

    def find_p_neighbors(self, p, mat=False):
        """
        Find a list with all the reduced equivalent of the `p`-neighbors of this ternary quadratic form, given by the zeros mod `p` of the form.
        See :meth:`find_p_neighbor_from_vec` for more information.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: Q0 = TernaryQF([1, 3, 3, -2, 0, -1]); Q0
            Ternary quadratic form with integer coefficients:
            [1 3 3]
            [-2 0 -1]
            sage: neig = Q0.find_p_neighbors(5)
            sage: len(neig)
            6
            sage: Q1 = TernaryQF([1, 1, 10, 1, 1, 1])
            sage: Q2 = TernaryQF([1, 2, 4, -1, -1, 0])
            sage: neig.count(Q0)
            2
            sage: neig.count(Q1)
            1
            sage: neig.count(Q2)
            3
        """

        z = self.find_zeros_mod_p(p)
        return [self.find_p_neighbor_from_vec(p, v, mat) for v in z]

    def basic_lemma(self, p):
        """
        Find a number represented by ``self`` and coprime to the prime `p`.

        EXAMPLES::

            sage: Q = TernaryQF([3, 3, 3, -2, 0, -1])
            sage: Q.basic_lemma(3)
            4
        """

        return _basic_lemma(self._a, self._b, self._c, self._r, self._s, self._t, p)

    def xi(self, p):
        """
        Return the value of the genus characters Xi_p... which may be
        missing one character. We allow `-1` as a prime.

        REFERENCES:

        Dickson's "Studies in the Theory of Numbers"

        EXAMPLES::

            sage: Q1 = TernaryQF([26, 42, 53, -36, -17, -3])
            sage: Q2 = Q1.find_p_neighbors(2)[1]
            sage: Q1.omega()
            3
            sage: Q1.xi(3), Q2.xi(3)
            (-1, -1)
        """
        if p == 4:
            p = -1
        if p == 8:
            p = 2

        if self.omega() % p != 0:
            raise ValueError("not a valid character")

        if p == -1 and self.omega() % 2**4 != 0:
            raise ValueError("not a valid character")

        if p == 2 and self.omega() % 2**5 != 0:
            raise ValueError("not a valid character")

        if (p == -1) or (p == 2):
            return kronecker_symbol(p, self.basic_lemma(2))

        return kronecker_symbol(self.basic_lemma(p), p)

    def xi_rec(self, p):
        """
        Return Xi(p) for the reciprocal form.

        EXAMPLES::

            sage: Q1 = TernaryQF([1, 1, 7, 0, 0, 0])
            sage: Q2 = Q1.find_p_neighbors(3)[0]
            sage: Q1.delta()
            28
            sage: Q1.xi_rec(7), Q2.xi_rec(7)
            (1, 1)
        """
        return self.reciprocal().xi(p)

    def symmetry(self, v):
        """
        Return `A`, the automorphism of the ternary quadratic form such that:

        - `Av = -v`,
        - `Au = 0`, if `u` is orthogonal to `v`,

        where `v` is a given vector.

        EXAMPLES::

            sage: Q = TernaryQF([4, 5, 8, 5, 2, 2])
            sage: v = vector((1,1,1))
            sage: M = Q.symmetry(v)
            sage: M
            [  7/13 -17/26 -23/26]
            [ -6/13   9/26 -23/26]
            [ -6/13 -17/26   3/26]
            sage: M.det()
            -1
            sage: M*v
            (-1, -1, -1)
            sage: v1 = vector((23, 0, -12))
            sage: v2 = vector((0, 23, -17))
            sage: v1*Q.matrix()*v
            0
            sage: v2*Q.matrix()*v
            0
            sage: M*v1 == v1
            True
            sage: M*v2 == v2
            True
        """

        return identity_matrix(3) - v.column()*matrix(v)*self.matrix()/self(v)

    def automorphism_symmetries(self, A):
        """
        Given the automorphism `A`, if `A` is the identity, return the empty list.
        Otherwise, return a list of two vectors `v_1`, `v_2` such that the product of
        the symmetries of the ternary quadratic form given by the two vectors is `A`.

        EXAMPLES::

            sage: Q = TernaryQF([9, 12, 30, -26, -28, 20])
            sage: A = matrix(ZZ, 3, [9, 10, -10, -6, -7, 6, 2, 2, -3])
            sage: Q(A) == Q
            True
            sage: v1, v2 = Q.automorphism_symmetries(A)
            sage: v1, v2
            ((8, -6, 2), (1, -5/4, -1/4))
            sage: A1 = Q.symmetry(v1)
            sage: A1
            [    9     9   -13]
            [   -6 -23/4  39/4]
            [    2   9/4  -9/4]
            sage: A2 = Q.symmetry(v2)
            sage: A2
            [    1     1     3]
            [    0  -1/4 -15/4]
            [    0  -1/4   1/4]
            sage: A1*A2 == A
            True
            sage: Q.automorphism_symmetries(identity_matrix(ZZ,3))
            []
        """

        if A == identity_matrix(3):
            return []
        else:
            bs = (A - 1).columns()
            for b1 in bs:
                if b1 != 0:
                    break
            A1 = self.symmetry(b1)*A
            bs = (A1 - 1).columns()
            for b2 in bs:
                if b2 != 0:
                    break
            return [b1, b2]

    def automorphism_spin_norm(self, A):
        """
        Return the spin norm of the automorphism `A`.

        EXAMPLES::

            sage: Q = TernaryQF([9, 12, 30, -26, -28, 20])
            sage: A = matrix(ZZ, 3, [9, 10, -10, -6, -7, 6, 2, 2, -3])
            sage: A.det()
            1
            sage: Q(A) == Q
            True
            sage: Q.automorphism_spin_norm(A)
            7
        """
        if A == identity_matrix(ZZ, 3):
            return 1
        bs = self.automorphism_symmetries(A)
        s = self(bs[0]) * self(bs[1])
        return s.squarefree_part()

    def _border(self, n):
        """
        Auxiliary function to find the automorphisms of a positive definite ternary quadratic form.
        It return a boolean whether the n-condition is true. If Q = TernaryQF([a,b,c,r,s,t]), the conditions are:

        1.  a = t, s = 2r.
        2.  a = s, t = 2r.
        3.  b = r, t = 2s.
        4.  a = -t.
        5.  a = -s.
        6.  b = -r.
        7.  a + b + r + s + t = 0, 2a + 2s + t = 0.
        8.  a = b, r = s.
        9.  b = c, s = t.
        10. r = s, r = 0.
        11. r = t, r = 0.
        12. s = t, s = 0.
        13. r = s, s = t, t = a.
        14. a = s, a = t.
        15. a = b, a + b + r + s + t = 0.
        16. a = b, b = c, a + b + r + s + t = 0.

        EXAMPLES::

            sage: Q01 = TernaryQF([5, 5, 9, 2, 4, 5])
            sage: Q01._border(1)
            True
            sage: Q02 = TernaryQF([6, 7, 8, 2, 6, 4])
            sage: Q02._border(2)
            True
            sage: Q03 = TernaryQF([6, 9, 9, 9, 3, 6])
            sage: Q03._border(3)
            True
            sage: Q04 = TernaryQF([1, 2, 3, -1, 0, -1])
            sage: Q04._border(4)
            True
            sage: Q05 = TernaryQF([2, 3, 5, -1, -2, 0])
            sage: Q05._border(5)
            True
            sage: Q06 = TernaryQF([1, 5, 7, -5, 0, 0])
            sage: Q06._border(6)
            True
            sage: Q07 = TernaryQF([1, 1, 7, -1, -1, 0])
            sage: Q07._border(7)
            True
            sage: Q08 = TernaryQF([2, 2, 5, -1, -1, -1])
            sage: Q08._border(8)
            True
            sage: Q09 = TernaryQF([3, 8, 8, 6, 2, 2])
            sage: Q09._border(9)
            True
            sage: Q10 = TernaryQF([1, 3, 4, 0, 0, 0])
            sage: Q10._border(10)
            True
            sage: Q11 = TernaryQF([3, 5, 8, 0, -1, 0])
            sage: Q11._border(11)
            True
            sage: Q12 = TernaryQF([2, 6, 7, -5, 0, 0])
            sage: Q12._border(12)
            True
            sage: Q13 = TernaryQF([1, 1, 2, 1, 1, 1])
            sage: Q13._border(13)
            True
            sage: Q14 = TernaryQF([1, 3, 4, 3, 1, 1])
            sage: Q14._border(14)
            True
            sage: Q15 = TernaryQF([3, 3, 6, -3, -3, 0])
            sage: Q15._border(15)
            True
            sage: Q16 = TernaryQF([4, 4, 4, -2, -3, -3])
            sage: Q16._border(16)
            True
        """

        a, b, c, r, s, t = self.coefficients()
        if n == 1:
            return (a == t) and (s == 2*r)
        elif n == 2:
            return (a == s) and (t == 2*r)
        elif n == 3:
            return (b == r) and (t == 2*s)
        elif n == 4:
            return (a == -t)
        elif n == 5:
            return (a == -s)
        elif n == 6:
            return (b == -r)
        elif n == 7:
            return (a + b + r + s + t == 0) and (2*a + 2*s + t == 0)
        elif n == 8:
            return (a == b) and (r == s)
        elif n == 9:
            return (b == c) and (s == t)
        elif n == 10:
            return (r == s) and (r == 0)
        elif n == 11:
            return (r == t) and (r == 0)
        elif n == 12:
            return (s == t) and (s == 0)
        elif n == 13:
            return (r == s) and (s == t) and (t == a)
        elif n == 14:
            return (a == s) and (a == t)
        elif n == 15:
            return (a == b) and (a + b + r + s + t == 0)
        elif n == 16:
            return (a == b) and (b == c) and (a + b + r + s + t == 0)

    def _borders(self):
        """
        Return the borders that the ternary quadratic form meet.

        See: TernaryQF._border

        EXAMPLES::

            sage: Q01 = TernaryQF([5, 5, 9, 2, 4, 5])
            sage: Q01._borders()
            (1,)
            sage: Q02 = TernaryQF([6, 7, 8, 2, 6, 4])
            sage: Q02._borders()
            (2,)
            sage: Q03 = TernaryQF([6, 9, 9, 9, 3, 6])
            sage: Q03._borders()
            (3,)
            sage: Q04 = TernaryQF([1, 2, 3, -1, 0, -1])
            sage: Q04._borders()
            (4,)
            sage: Q05 = TernaryQF([2, 3, 5, -1, -2, 0])
            sage: Q05._borders()
            (5,)
            sage: Q06 = TernaryQF([1, 5, 7, -5, 0, 0])
            sage: Q06._borders()
            (6, 12)
            sage: Q07 = TernaryQF([1, 1, 7, -1, -1, 0])
            sage: Q07._borders()
            (5, 6, 7, 8, 15)
            sage: Q08 = TernaryQF([2, 2, 5, -1, -1, -1])
            sage: Q08._borders()
            (8,)
            sage: Q09 = TernaryQF([3, 8, 8, 6, 2, 2])
            sage: Q09._borders()
            (9,)
            sage: Q10 = TernaryQF([1, 3, 4, 0, 0, 0])
            sage: Q10._borders()
            (10, 11, 12)
            sage: Q11 = TernaryQF([3, 5, 8, 0, -1, 0])
            sage: Q11._borders()
            (11,)
            sage: Q12 = TernaryQF([2, 6, 7, -5, 0, 0])
            sage: Q12._borders()
            (12,)
            sage: Q13 = TernaryQF([1, 1, 2, 1, 1, 1])
            sage: Q13._borders()
            (8, 13, 14)
            sage: Q14 = TernaryQF([1, 3, 4, 3, 1, 1])
            sage: Q14._borders()
            (14,)
            sage: Q15 = TernaryQF([3, 3, 6, -3, -3, 0])
            sage: Q15._borders()
            (5, 6, 7, 8, 15)
            sage: Q16 = TernaryQF([4, 4, 4, -2, -3, -3])
            sage: Q16._borders()
            (9, 15, 16)
        """
        return tuple(n for n in range(1, 17) if self._border(n))

    def _automorphisms_reduced_fast(self):
        """
        Return the coefficients of the matrices of the automorphisms of the reduced ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 7, 0, 0, 0])
            sage: Q.is_eisenstein_reduced()
            True
            sage: auts = Q._automorphisms_reduced_fast()
            sage: len(auts)
            8
            sage: A = matrix(3, auts[randint(0,7)])
            sage: Q(A) == Q
            True
            sage: Q = TernaryQF([3, 4, 5, 3, 3, 2])
            sage: Q._automorphisms_reduced_fast()
            [(1, 0, 0, 0, 1, 0, 0, 0, 1)]
        """

        if self._border(1):
            if self._border(2):
                if self._border(14):
                    if self._border(9):
                        # borders 1, 2, 9, 14
                        return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                                (-1, -1, -1, 0, 0, 1, 0, 1, 0),
                                (-1, -1, 0, 0, 1, 0, 0, 0, -1),
                                (-1, 0, -1, 0, -1, 0, 0, 0, 1),
                                (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                                (1, 0, 1, 0, 0, -1, 0, 1, 0),
                                (1, 1, 0, 0, 0, 1, 0, -1, 0),
                                (1, 1, 1, 0, -1, 0, 0, 0, -1)]
                    else:
                        # borders 1, 2, 14
                        return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                                (-1, -1, 0, 0, 1, 0, 0, 0, -1),
                                (-1, 0, -1, 0, -1, 0, 0, 0, 1),
                                (1, 1, 1, 0, -1, 0, 0, 0, -1)]
            else:
                # borders 1
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, -1, 0, 0, 1, 0, 0, 0, -1)]

        if self._border(2):
            # borders 2
            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                    (-1, 0, -1, 0, -1, 0, 0, 0, 1)]

        if self._border(3):
            # borders 3
            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                    (-1, 0, 0, 0, -1, -1, 0, 0, 1)]

        if self._border(4):
            if self._border(10):
                if self._border(8):
                    # borders 4, 8, 10
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, -1, 1, 0, 0, 0, -1),
                            (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                            (-1, 1, 0, -1, 0, 0, 0, 0, 1),
                            (-1, 1, 0, 0, 1, 0, 0, 0, -1),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, -1, 0, 1, -1, 0, 0, 0, 1),
                            (0, 1, 0, -1, 1, 0, 0, 0, 1),
                            (0, 1, 0, 1, 0, 0, 0, 0, -1),
                            (1, -1, 0, 0, -1, 0, 0, 0, -1),
                            (1, -1, 0, 1, 0, 0, 0, 0, 1),
                            (1, 0, 0, 1, -1, 0, 0, 0, -1)]
                else:
                    # borders 4, 10
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                            (-1, 1, 0, 0, 1, 0, 0, 0, -1),
                            (1, -1, 0, 0, -1, 0, 0, 0, -1)]
            else:
                # borders 4
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (1, -1, 0, 0, -1, 0, 0, 0, -1)]

        if self._border(5):
            if self._border(6):
                if self._border(7):
                    if self._border(8):
                        if self._border(15):
                            # borders 5, 6, 7, 8, 15
                            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                                    (-1, 0, 0, 0, 1, -1, 0, 0, -1),
                                    (-1, 0, 1, 0, -1, 1, 0, 0, 1),
                                    (0, -1, 0, -1, 0, 0, 0, 0, -1),
                                    (0, -1, 1, 1, 0, 0, 0, 0, 1),
                                    (0, 1, -1, 1, 0, -1, 0, 0, -1),
                                    (0, 1, 0, -1, 0, 1, 0, 0, 1),
                                    (1, 0, -1, 0, -1, 0, 0, 0, -1)]
                    else:
                        # borders 5, 6, 7
                        return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                                (-1, 0, 0, 0, 1, -1, 0, 0, -1),
                                (-1, 0, 1, 0, -1, 1, 0, 0, 1),
                                (1, 0, -1, 0, -1, 0, 0, 0, -1)]
            elif self._border(11):
                # borders 5, 11
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, 1, 0, 0, 0, -1),
                        (-1, 0, 1, 0, -1, 0, 0, 0, 1),
                        (1, 0, -1, 0, -1, 0, 0, 0, -1)]
            else:
                # borders 5
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (1, 0, -1, 0, -1, 0, 0, 0, -1)]

        if self._border(6):
            if self._border(12):
                if self._border(9):
                    # borders 6, 9, 12
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 0, -1, 1),
                            (-1, 0, 0, 0, -1, 1, 0, 0, 1),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (-1, 0, 0, 0, 0, 1, 0, 1, 0),
                            (-1, 0, 0, 0, 1, -1, 0, 0, -1),
                            (-1, 0, 0, 0, 1, 0, 0, 1, -1),
                            (1, 0, 0, 0, -1, 0, 0, 0, -1),
                            (1, 0, 0, 0, -1, 1, 0, -1, 0),
                            (1, 0, 0, 0, 0, -1, 0, 1, -1),
                            (1, 0, 0, 0, 0, 1, 0, -1, 1),
                            (1, 0, 0, 0, 1, -1, 0, 1, 0)]
                else:
                    # borders 6, 12
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 1, 0, 0, 1),
                            (-1, 0, 0, 0, 1, -1, 0, 0, -1),
                            (1, 0, 0, 0, -1, 0, 0, 0, -1)]
            else:
                # borders 6
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, 1, -1, 0, 0, -1)]

        if self._border(7):
            if self._border(8) and self._border(15):
                if self._border(16):
                    if self._border(9):
                        # borders 7, 8, 9, 15, 16
                        return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                                (-1, 0, 0, -1, 0, 1, -1, 1, 0),
                                (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                                (-1, 0, 1, -1, 1, 0, -1, 0, 0),
                                (-1, 0, 1, 0, -1, 1, 0, 0, 1),
                                (-1, 1, 0, -1, 0, 0, -1, 0, 1),
                                (-1, 1, 0, 0, 1, 0, 0, 1, -1),
                                (0, -1, 0, -1, 0, 0, 0, 0, -1),
                                (0, -1, 0, 1, -1, 0, 0, -1, 1),
                                (0, -1, 1, 0, -1, 0, 1, -1, 0),
                                (0, -1, 1, 0, 0, 1, -1, 0, 1),
                                (0, 0, -1, 0, -1, 0, -1, 0, 0),
                                (0, 0, -1, 0, 1, -1, 1, 0, -1),
                                (0, 0, 1, -1, 0, 1, 0, -1, 1),
                                (0, 0, 1, 1, 0, 0, 0, 1, 0),
                                (0, 1, -1, -1, 1, 0, 0, 1, 0),
                                (0, 1, -1, 1, 0, -1, 0, 0, -1),
                                (0, 1, 0, 0, 0, 1, 1, 0, 0),
                                (0, 1, 0, 0, 1, -1, -1, 1, 0),
                                (1, -1, 0, 0, -1, 1, 0, -1, 0),
                                (1, -1, 0, 1, 0, -1, 1, 0, 0),
                                (1, 0, -1, 0, 0, -1, 0, 1, -1),
                                (1, 0, -1, 1, 0, 0, 1, -1, 0),
                                (1, 0, 0, 1, -1, 0, 1, 0, -1)]
                    else:
                        # borders 7, 8, 15, 16
                        return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                                (-1, 0, 0, -1, 0, 1, -1, 1, 0),
                                (-1, 0, 1, 0, -1, 1, 0, 0, 1),
                                (0, -1, 0, -1, 0, 0, 0, 0, -1),
                                (0, -1, 1, 0, -1, 0, 1, -1, 0),
                                (0, 1, -1, 1, 0, -1, 0, 0, -1),
                                (0, 1, 0, 0, 1, -1, -1, 1, 0),
                                (1, 0, -1, 1, 0, 0, 1, -1, 0)]
                else:
                    # borders 7, 8, 15
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 1, 0, -1, 1, 0, 0, 1),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, 1, -1, 1, 0, -1, 0, 0, -1)]
            elif self._border(9):
                # borders 7, 9
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                        (-1, 0, 1, 0, -1, 1, 0, 0, 1),
                        (-1, 1, 0, 0, 1, 0, 0, 1, -1),
                        (1, -1, 0, 0, -1, 1, 0, -1, 0),
                        (1, 0, -1, 0, 0, -1, 0, 1, -1)]
            else:
                # borders 7
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 1, 0, -1, 1, 0, 0, 1)]

        if self._border(8):
            if self._border(9):
                if self._border(10) and self._border(11) and self._border(12):
                    # borders 8, 9, 10, 11, 12
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (-1, 0, 0, 0, 0, 1, 0, 1, 0),
                            (-1, 0, 0, 0, 1, 0, 0, 0, -1),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, -1, 0, 0, 0, -1, 1, 0, 0),
                            (0, -1, 0, 0, 0, 1, -1, 0, 0),
                            (0, -1, 0, 1, 0, 0, 0, 0, 1),
                            (0, 0, -1, -1, 0, 0, 0, 1, 0),
                            (0, 0, -1, 0, -1, 0, -1, 0, 0),
                            (0, 0, -1, 0, 1, 0, 1, 0, 0),
                            (0, 0, -1, 1, 0, 0, 0, -1, 0),
                            (0, 0, 1, -1, 0, 0, 0, -1, 0),
                            (0, 0, 1, 0, -1, 0, 1, 0, 0),
                            (0, 0, 1, 0, 1, 0, -1, 0, 0),
                            (0, 0, 1, 1, 0, 0, 0, 1, 0),
                            (0, 1, 0, -1, 0, 0, 0, 0, 1),
                            (0, 1, 0, 0, 0, -1, -1, 0, 0),
                            (0, 1, 0, 0, 0, 1, 1, 0, 0),
                            (0, 1, 0, 1, 0, 0, 0, 0, -1),
                            (1, 0, 0, 0, -1, 0, 0, 0, -1),
                            (1, 0, 0, 0, 0, -1, 0, 1, 0),
                            (1, 0, 0, 0, 0, 1, 0, -1, 0)]
                elif self._border(13) and self._border(14):
                    # borders 8, 9, 13, 14
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, -1, -1, 0, 0, 1, 0, 1, 0),
                            (-1, -1, -1, 0, 1, 0, 1, 0, 0),
                            (-1, -1, -1, 1, 0, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 1, 1, 1),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (-1, 0, 0, 1, 1, 1, 0, 0, -1),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, -1, 0, 0, 0, -1, 1, 1, 1),
                            (0, -1, 0, 1, 1, 1, -1, 0, 0),
                            (0, 0, -1, -1, 0, 0, 1, 1, 1),
                            (0, 0, -1, 0, -1, 0, -1, 0, 0),
                            (0, 0, -1, 1, 1, 1, 0, -1, 0),
                            (0, 0, 1, -1, -1, -1, 1, 0, 0),
                            (0, 0, 1, 0, 1, 0, -1, -1, -1),
                            (0, 0, 1, 1, 0, 0, 0, 1, 0),
                            (0, 1, 0, -1, -1, -1, 0, 0, 1),
                            (0, 1, 0, 0, 0, 1, 1, 0, 0),
                            (0, 1, 0, 1, 0, 0, -1, -1, -1),
                            (1, 0, 0, -1, -1, -1, 0, 1, 0),
                            (1, 0, 0, 0, 0, 1, -1, -1, -1),
                            (1, 1, 1, -1, 0, 0, 0, -1, 0),
                            (1, 1, 1, 0, -1, 0, 0, 0, -1),
                            (1, 1, 1, 0, 0, -1, -1, 0, 0)]
                else:
                    # borders 8, 9
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, 0, -1, 0, -1, 0, -1, 0, 0),
                            (0, 0, 1, 1, 0, 0, 0, 1, 0),
                            (0, 1, 0, 0, 0, 1, 1, 0, 0)]
            elif self._border(10):
                if self._border(11) and self._border(12):
                    # borders 8, 10, 11, 12
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, 1, 0, 0, 0, -1),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, -1, 0, 1, 0, 0, 0, 0, 1),
                            (0, 1, 0, -1, 0, 0, 0, 0, 1),
                            (0, 1, 0, 1, 0, 0, 0, 0, -1),
                            (1, 0, 0, 0, -1, 0, 0, 0, -1)]
                else:
                    # borders 8, 10
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                            (0, -1, 0, -1, 0, 0, 0, 0, -1),
                            (0, 1, 0, 1, 0, 0, 0, 0, -1)]
            elif self._border(14):
                # borders 8, 13, 14
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, -1, -1, 1, 0, 0, 0, 0, 1),
                        (-1, 0, 0, 1, 1, 1, 0, 0, -1),
                        (0, -1, 0, -1, 0, 0, 0, 0, -1),
                        (0, 1, 0, -1, -1, -1, 0, 0, 1),
                        (1, 1, 1, 0, -1, 0, 0, 0, -1)]
            else:
                # borders 8
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (0, -1, 0, -1, 0, 0, 0, 0, -1)]

        if self._border(9):
            if self._border(12):
                if self._border(10) and self._border(11):
                    # borders 9, 10, 11, 12
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (-1, 0, 0, 0, 0, 1, 0, 1, 0),
                            (-1, 0, 0, 0, 1, 0, 0, 0, -1),
                            (1, 0, 0, 0, -1, 0, 0, 0, -1),
                            (1, 0, 0, 0, 0, -1, 0, 1, 0),
                            (1, 0, 0, 0, 0, 1, 0, -1, 0)]
                else:
                    # borders 9, 12
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (-1, 0, 0, 0, 0, 1, 0, 1, 0),
                            (1, 0, 0, 0, -1, 0, 0, 0, -1)]
            elif self._border(14):
                if self._border(13):
                    # borders 9, 13, 14
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, -1, -1, 0, 0, 1, 0, 1, 0),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (1, 1, 1, 0, -1, 0, 0, 0, -1)]
                else:
                    # borders 9, 14
                    return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                            (-1, -1, -1, 0, 0, 1, 0, 1, 0),
                            (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                            (1, 1, 1, 0, -1, 0, 0, 0, -1)]
            elif self._border(15):
                # borders 9, 15, 16
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, -1, 0, 1, -1, 1, 0),
                        (-1, 0, 0, 0, 0, -1, 0, -1, 0),
                        (0, -1, 1, 0, -1, 0, 1, -1, 0),
                        (0, -1, 1, 0, 0, 1, -1, 0, 1),
                        (0, 1, -1, -1, 1, 0, 0, 1, 0),
                        (0, 1, -1, 1, 0, -1, 0, 0, -1),
                        (1, 0, 0, 1, -1, 0, 1, 0, -1)]
            else:
                # borders 9
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, 0, -1, 0, -1, 0)]

        if self._border(10):
            if self._border(11) and self._border(12):
                # borders 10, 11, 12
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, -1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, 1, 0, 0, 0, -1),
                        (1, 0, 0, 0, -1, 0, 0, 0, -1)]
            else:
                # borders 10
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, 0, -1, 0, 0, 0, 1)]

        if self._border(11):
            # borders 11
            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                    (-1, 0, 0, 0, 1, 0, 0, 0, -1)]

        if self._border(12):
            # border 12
            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                    (1, 0, 0, 0, -1, 0, 0, 0, -1)]

        if self._border(13) and self._border(14):
            # border 13, 14
            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                    (1, 1, 1, 0, -1, 0, 0, 0, -1)]

        if self._border(14):
            # border 14
            return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                    (1, 1, 1, 0, -1, 0, 0, 0, -1)]

        if self._border(15):
            if self._border(16):
                # borders 15, 16
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (-1, 0, 0, -1, 0, 1, -1, 1, 0),
                        (0, -1, 1, 0, -1, 0, 1, -1, 0),
                        (0, 1, -1, 1, 0, -1, 0, 0, -1)]
            else:
                # borders 15
                return [(1, 0, 0, 0, 1, 0, 0, 0, 1),
                        (0, 1, -1, 1, 0, -1, 0, 0, -1)]

        return [(1, 0, 0, 0, 1, 0, 0, 0, 1)]

    def _automorphisms_reduced_slow(self):
        """
        Return the automorphisms of the reduced ternary quadratic form.
        It searches over all 3x3 matrices with coefficients -1, 0, 1,
        determinant 1 and finite order, because Eisenstein reduced forms
        are Minkowski reduced. See Cassels.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 7, 0, 0, 0])
            sage: Q.is_eisenstein_reduced()
            True

            sage: # long time
            sage: auts = Q._automorphisms_reduced_slow()  # 3s on sage.math, 2014
            sage: len(auts)
            8
            sage: A = auts[randint(0,7)]
            sage: Q(A) == Q
            True
            sage: Q = TernaryQF([3, 4, 5, 3, 3, 2])
            sage: Q._automorphisms_reduced_slow()
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]
        """
        from itertools import product
        if TernaryQF.possible_automorphisms is None:
            auts = (matrix(ZZ, 3, 3, m) for m in product([-1, 0, 1], repeat=9))
            auts = [m for m in auts if m.det() == 1]
            auts = [m for m in auts if m**2 in auts]
            auts = [m for m in auts if m**2 in auts]
            auts = [m for m in auts if m**2 in auts]
            TernaryQF.possible_automorphisms = auts
        return [m for m in TernaryQF.possible_automorphisms if self(m) == self]

    def automorphisms(self, slow=True):
        """
        Return a list with the automorphisms of the definite ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 7, 0, 0, 0])
            sage: auts = Q.automorphisms(); auts
            [
            [-1  0  0]  [-1  0  0]  [ 0 -1  0]  [ 0 -1  0]  [ 0  1  0]  [ 0  1  0]
            [ 0 -1  0]  [ 0  1  0]  [-1  0  0]  [ 1  0  0]  [-1  0  0]  [ 1  0  0]
            [ 0  0  1], [ 0  0 -1], [ 0  0 -1], [ 0  0  1], [ 0  0  1], [ 0  0 -1],
            [ 1  0  0]  [1 0 0]
            [ 0 -1  0]  [0 1 0]
            [ 0  0 -1], [0 0 1]
            ]
            sage: all(Q == Q(A) for A in auts)
            True
            sage: Q = TernaryQF([3, 4, 5, 3, 3, 2])
            sage: Q.automorphisms(slow=False)
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]
            sage: Q = TernaryQF([4, 2, 4, 3, -4, -5])
            sage: auts = Q.automorphisms(slow=False)
            sage: auts
            [
            [1 0 0]  [ 2 -1 -1]
            [0 1 0]  [ 3 -2 -1]
            [0 0 1], [ 0  0 -1]
            ]
            sage: A = auts[1]
            sage: Q(A) == Q
            True
            sage: Qr, M_red = Q.reduced_form_eisenstein()
            sage: Qr
            Ternary quadratic form with integer coefficients:
            [1 2 3]
            [-1 0 -1]
            sage: Q(A*M_red) == Qr
            True
        """
        if not self.is_definite():
            raise ValueError("only implemented for definite forms")

        if self._automorphisms is not None:
            return self._automorphisms

        if self.is_positive_definite():
            if self.is_eisenstein_reduced():
                if slow:
                    self._automorphisms = self._automorphisms_reduced_slow()
                else:
                    auts = self._automorphisms_reduced_fast()
                    self._automorphisms = [matrix(ZZ, 3, A) for A in auts]
            else:
                Qr, M = self.reduced_form_eisenstein()
                auts = Qr.automorphisms(slow)
                M_inv = M.inverse()
                self._automorphisms = [M*m*M_inv for m in auts]
        else:
            self._automorphisms = (-self).automorphisms()
        return self._automorphisms

    def _number_of_automorphisms_reduced(self):
        """
        Return the number of automorphisms of the reduced definite ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 7, 0, 0, 0])
            sage: Q._number_of_automorphisms_reduced()
            8
            sage: len(Q.automorphisms(slow=False))
            8
            sage: Q = TernaryQF([3, 4, 5, 3, 3, 2])
            sage: Q._number_of_automorphisms_reduced()
            1
        """
        if self._border(1):
            if self._border(2):
                if self._border(14):
                    if self._border(9):
                        # borders 1, 2, 9, 14
                        return 8
                    else:
                        # borders 1, 2, 14
                        return 4
            else:
                # borders 1
                return 2

        if self._border(2):
            # borders 2
            return 2

        if self._border(3):
            # borders 3
            return 2

        if self._border(4):
            if self._border(10):
                if self._border(8):
                    # borders 4, 8, 10
                    return 12
                else:
                    # borders 4, 10
                    return 4
            else:
                # borders 4
                return 2

        if self._border(5):
            if self._border(6):
                if self._border(7):
                    if self._border(8):
                        if self._border(15):
                            # borders 5, 6, 7, 8, 15
                            return 8
                    else:
                        # borders 5, 6, 7
                        return 4
            elif self._border(11):
                # borders 5, 11
                return 4
            else:
                # borders 5
                return 2

        if self._border(6):
            if self._border(12):
                if self._border(9):
                    # borders 6, 9, 12
                    return 12
                else:
                    # borders 6, 12
                    return 4
            else:
                # borders 6
                return 2

        if self._border(7):
            if self._border(8) and self._border(15):
                if self._border(16):
                    if self._border(9):
                        # borders 7, 8, 9, 15, 16
                        return 24
                    else:
                        # borders 7, 8, 15, 16
                        return 8
                else:
                    # borders 7, 8, 15
                    return 4
            elif self._border(9):
                # borders 7, 9
                return 6
            else:
                # borders 7
                return 2

        if self._border(8):
            if self._border(9):
                if self._border(10) and self._border(11) and self._border(12):
                    # borders 8, 9, 10, 11, 12
                    return 24
                elif self._border(13) and self._border(14):
                    # borders 8, 9, 13, 14
                    return 24
                else:
                    # borders 8, 9
                    return 6
            elif self._border(10):
                if self._border(11) and self._border(12):
                    # borders 8, 10, 11, 12
                    return 8
                else:
                    # borders 8, 10
                    return 4
            elif self._border(14):
                # borders 8, 13, 14
                return 6
            else:
                # borders 8
                return 2

        if self._border(9):
            if self._border(12):
                if self._border(10) and self._border(11):
                    # borders 9, 10, 11, 12
                    return 8
                else:
                    # borders 9, 12
                    return 4
            elif self._border(14):
                if self._border(13):
                    # borders 9, 13, 14
                    return 4
                else:
                    # borders 9, 14
                    return 4
            elif self._border(15):
                # borders 9, 15, 16
                return 8
            else:
                # borders 9
                return 2

        if self._border(10):
            if self._border(11) and self._border(12):
                # borders 10, 11, 12
                return 4
            else:
                # borders 10
                return 2

        if self._border(11):
            # borders 11
            return 2

        if self._border(12):
            # border 12
            return 2

        if self._border(13) and self._border(14):
            # border 13, 14
            return 2

        if self._border(14):
            # border 14
            return 2

        if self._border(15):
            if self._border(16):
                # borders 15, 16
                return 4
            else:
                # borders 15
                return 2

        return 1

    def number_of_automorphisms(self, slow=True):
        """
        Return the number of automorphisms of the definite ternary quadratic form.

        EXAMPLES::

            sage: Q = TernaryQF([1, 1, 7, 0, 0, 0])
            sage: A = matrix(ZZ, 3, [0, 1, 0, -1, 5, 0, -8, -1, 1])
            sage: A.det()
            1
            sage: Q1 = Q(A); Q1
            Ternary quadratic form with integer coefficients:
            [449 33 7]
            [-14 -112 102]
            sage: Q1.number_of_automorphisms()
            8
            sage: Q = TernaryQF([-19, -7, -6, -12, 20, 23])
            sage: Q.is_negative_definite()
            True
            sage: Q.number_of_automorphisms(slow=False)
            24
        """
        if not self.is_definite():
            raise ValueError("only implemented for definite forms")

        if self._number_of_automorphisms is not None:
            return self._number_of_automorphisms

        if slow:
            self._number_of_automorphisms = len(self.automorphisms())
        else:
            if self.is_negative_definite():
                self._number_of_automorphisms = (-self).reduced_form_eisenstein(False)._number_of_automorphisms_reduced()
            else:
                self._number_of_automorphisms = self.reduced_form_eisenstein(False)._number_of_automorphisms_reduced()

        return self._number_of_automorphisms


def find_all_ternary_qf_by_level_disc(N, d):
    """
    Find the coefficients of all the reduced ternary quadratic forms given its discriminant `d` and level `N`.

    If `N|4d` and `d|N^2`, then it may be some forms with that discriminant and level.

    EXAMPLES::

        sage: find_all_ternary_qf_by_level_disc(44, 11)
        [Ternary quadratic form with integer coefficients:
        [1 1 3]
        [0 -1 0], Ternary quadratic form with integer coefficients:
        [1 1 4]
        [1 1 1]]
        sage: find_all_ternary_qf_by_level_disc(44, 11^2 * 16)
        [Ternary quadratic form with integer coefficients:
        [3 15 15]
        [-14 -2 -2], Ternary quadratic form with integer coefficients:
        [4 11 12]
        [0 -4 0]]
        sage: Q = TernaryQF([1, 1, 3, 0, -1, 0])
        sage: Q.is_eisenstein_reduced()
        True
        sage: Q.reciprocal_reduced()
        Ternary quadratic form with integer coefficients:
        [4 11 12]
        [0 -4 0]
        sage: find_all_ternary_qf_by_level_disc(44, 22)
        []
        sage: find_all_ternary_qf_by_level_disc(44, 33)
        Traceback (most recent call last):
        ...
        ValueError: There are no ternary forms of this level and discriminant
    """
    return [TernaryQF(_) for _ in _find_all_ternary_qf_by_level_disc(N, d)]


def find_a_ternary_qf_by_level_disc(N, d):
    """
    Find a reduced ternary quadratic form given its discriminant `d` and level `N`.
    If `N|4d` and `d|N^2`, then it may be a form with that discriminant and level.

    EXAMPLES::

        sage: Q1 = find_a_ternary_qf_by_level_disc(44, 11); Q1
        Ternary quadratic form with integer coefficients:
        [1 1 3]
        [0 -1 0]
        sage: Q2 = find_a_ternary_qf_by_level_disc(44, 11^2 * 16)
        sage: Q2
        Ternary quadratic form with integer coefficients:
        [3 15 15]
        [-14 -2 -2]
        sage: Q1.is_eisenstein_reduced()
        True
        sage: Q1.level()
        44
        sage: Q1.disc()
        11
        sage: find_a_ternary_qf_by_level_disc(44, 22)
        sage: find_a_ternary_qf_by_level_disc(44, 33)
        Traceback (most recent call last):
        ...
        ValueError: There are no ternary forms of this level and discriminant
    """
    q = _find_a_ternary_qf_by_level_disc(N, d)
    if q is not None:
        return TernaryQF(q)
