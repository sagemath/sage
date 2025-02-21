# distutils: language = c++
# distutils: libraries = gmp m NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
"""
Elements of Quaternion Algebras

Sage allows for computation with elements of quaternion algebras over
a nearly arbitrary base field of characteristic not 2.  Sage also has
very highly optimized implementation of arithmetic in rational
quaternion algebras and quaternion algebras over number fields.

TESTS:

Check that :issue:`20829` is fixed::

    sage: D.<i,j,k> = QuaternionAlgebra(QQ,-1,-3)
    sage: hash(i)
    184301497
"""

# ****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#       Copyright (C) 2009 Jonathan Bober <jwbober@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element cimport AlgebraElement, Element
from sage.structure.richcmp cimport rich_to_bool, rich_to_bool_sgn, richcmp_item
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.number_field.number_field_element cimport NumberFieldElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix


from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.ntl.convert cimport mpz_to_ZZ, ZZ_to_mpz
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_poly_sage cimport *
from sage.libs.flint.ntl_interface cimport *

# variables for holding temporary values computed in
# QuaternionAlgebraElement_rational_field._mul_()
cdef mpz_t T1, T2, t3, t4, t5, t6, t7, t8, s1, s2, U1, U2
cdef fmpz_poly_t fT1, fT2, ft3, ft4, ft5, ft6, ft7, ft8, fs1, fs2, fU1, fU2

cdef _clear_globals():
    """
    Clear all global variables allocated for optimization of
    quaternion algebra arithmetic.

    Do *not* call this yourself.
    """
    mpz_clear(T1)
    mpz_clear(T2)
    mpz_clear(t3)
    mpz_clear(t4)
    mpz_clear(t5)
    mpz_clear(t6)
    mpz_clear(t7)
    mpz_clear(t8)

    mpz_clear(s1)
    mpz_clear(s2)

    mpz_clear(U1)
    mpz_clear(U2)

    fmpz_poly_clear(fT1)
    fmpz_poly_clear(fT2)
    fmpz_poly_clear(ft3)
    fmpz_poly_clear(ft4)
    fmpz_poly_clear(ft5)
    fmpz_poly_clear(ft6)
    fmpz_poly_clear(ft7)
    fmpz_poly_clear(ft8)

    fmpz_poly_clear(fs1)
    fmpz_poly_clear(fs2)

    fmpz_poly_clear(fU1)
    fmpz_poly_clear(fU2)

cdef _init_globals():
    """
    Initialize all global variables allocated for optimization of
    quaternion algebra arithmetic, and register a hook to eventually
    clear them.

    Do *not* call this yourself.
    """
    # over QQ
    mpz_init(T1)
    mpz_init(T2)
    mpz_init(t3)
    mpz_init(t4)
    mpz_init(t5)
    mpz_init(t6)
    mpz_init(t7)
    mpz_init(t8)

    mpz_init(s1)
    mpz_init(s2)

    mpz_init(U1)
    mpz_init(U2)

    # Number fields
    fmpz_poly_init(fT1)
    fmpz_poly_init(fT2)
    fmpz_poly_init(ft3)
    fmpz_poly_init(ft4)
    fmpz_poly_init(ft5)
    fmpz_poly_init(ft6)
    fmpz_poly_init(ft7)
    fmpz_poly_init(ft8)

    fmpz_poly_init(fs1)
    fmpz_poly_init(fs2)

    fmpz_poly_init(fU1)
    fmpz_poly_init(fU2)

    # ...and clear them when sage terminates.
    import atexit
    atexit.register(_clear_globals)

# Initialize module-scope global C variables.
_init_globals()

cdef to_quaternion(R, x):
    """
    Internal function used implicitly by quaternion algebra creation.

    INPUT:

    - ``R`` -- callable
    - ``x`` -- element or 4-tuple

    Given a callable R and an x that defines a quaternion, which can be a
    4-tuple, list of length 4, or something that coerces to R, return
    4-tuple of elements of R.

    EXAMPLES::

        sage: Q.<i,j,kkkk> = QuaternionAlgebra(QQ,-7, 13)
        sage: kkkk._repr_()   # implicit doctest
        'kkkk'
    """
    if isinstance(x, (list, tuple)):
        return R(x[0]), R(x[1]), R(x[2]), R(x[3])
    else:
        return R(x), R(0), R(0), R(0)

cdef inline print_coeff(y, i, bint atomic):
    """
    Internal function used implicitly by all quaternion algebra printing.

    INPUT:

    - ``y`` -- coefficient
    - ``i`` -- string (name of a generator)
    - ``atomic`` -- boolean int; whether or not elements of base ring
      print atomically

    EXAMPLES::

        sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-7, 13)
        sage: i._repr_()   # implicit doctest
        'i'
    """
    if not y:
        return ''
    if y == 1:
        return i
    elif y == -1:
        return "-%s" % i
    y = str(y)
    if not atomic and ('+' in y or '-' in y):
        return '(%s)*%s' % (y, i)
    else:
        return '%s*%s' % (y, i)


cdef class QuaternionAlgebraElement_abstract(AlgebraElement):
    def __hash__(self):
        r"""
        TESTS::

            sage: from itertools import product
            sage: for K in [QQ, QuadraticField(2), AA, Frac(QQ['x'])]:
            ....:     Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            ....:     assert hash(Q.one()) == hash(K.one())
            ....:     assert hash(Q(2)) == hash(K(2))
            ....:     elts = []
            ....:     for (x,y,z,w) in product([K(0), K(1), K(2), K(-1)], repeat=4):
            ....:         elts.append(x + y*i + z*j + w*k)
            ....:     assert len(set(map(hash, elts))) == len(elts)
        """
        cdef long h
        h = hash(self[0])
        x = self[1]
        if x:
            h = ((h+14152L)*13023L) ^ hash(x)
        x = self[2]
        if x:
            h = ((h+33325L)*31321L) ^ hash(x)
        x = self[3]
        if x:
            h = ((h+34125L)*51125L) ^ hash(x)
        return h

    cpdef bint is_constant(self) noexcept:
        """
        Return ``True`` if this quaternion is constant, i.e., has no `i`, `j`,
        or `k` term.

        OUTPUT: boolean

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(-1,-2)
            sage: A(1).is_constant()
            True
            sage: A(1+i).is_constant()
            False
            sage: A(i).is_constant()
            False
        """
        return not (self[1] or self[2] or self[3])

    def __int__(self):
        """
        Try to coerce this quaternion to a Python int.

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: int(A(-3))
            -3
            sage: int(A(-3/2))
            -1
            sage: int(-3 + i)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if self.is_constant():
            return int(self[0])
        raise TypeError

    def __float__(self):
        """
        Try to coerce this quaternion to a Python float.

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: float(A(-3/2))
            -1.5
            sage: float(A(-3))
            -3.0
            sage: float(-3 + i)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if self.is_constant():
            return float(self[0])
        raise TypeError

    def _integer_(self, ZZ=None):
        """
        Try to coerce this quaternion to an Integer.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(-1,-2)
            sage: Integer(A(-3))               # indirect doctest
            -3
            sage: type(Integer(A(-3)))
            <class 'sage.rings.integer.Integer'>
            sage: Integer(A(-3/2))
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: Integer(-3 + i)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if self.is_constant():
            return Integer(self[0])
        raise TypeError

    def _rational_(self):
        """
        Try to coerce this quaternion to a Rational.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x']),-5,-2)
            sage: Rational(Q(2/3))                  # indirect doctest
            2/3
            sage: Rational(2/3 + i)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if self.is_constant():
            return Rational(self[0])
        raise TypeError

    def __bool__(self):
        """
        Return ``True`` if this quaternion is nonzero.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x']),-5,-2)
            sage: bool(i+j)
            True
            sage: bool(Q(0))
            False
        """
        return self[0] or self[1] or self[2] or self[3]

    cdef _do_print(self, x, y, z, w):
        """
        Used internally by the print function.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(-17,-19)
            sage: str(i+j+k-3/4)                            # indirect doctest
            '-3/4 + i + j + k'
        """
        cdef bint atomic = self._parent._base._repr_option('element_is_atomic')
        v = []
        i, j, k = self._parent.variable_names()
        if x:
            v.append(str(x))
        c = print_coeff(y, i, atomic)
        if c:
            v.append(c)
        c = print_coeff(z, j, atomic)
        if c:
            v.append(c)
        c = print_coeff(w, k, atomic)
        if c:
            v.append(c)
        if not v:
            return '0'
        return ' + '.join(v).replace('+ -', '- ')

    def _repr_(self):
        """
        Return string representation of this quaternion.

        EXAMPLES::

            sage: R.<x> = Frac(QQ['x'])
            sage: Q.<i,j,k> = QuaternionAlgebra(R,-5*x,-2)
            sage: a = x + i*x^3 + j*x^2 + k*x
            sage: a._repr_()
            'x + x^3*i + x^2*j + x*k'
            sage: a = x+2/3 + i*x^3 + j*(x^2-5/2) + k*x
            sage: a._repr_()
            'x + 2/3 + x^3*i + (x^2 - 5/2)*j + x*k'
            sage: type(a)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: Q(0)._repr_()
            '0'
        """
        return self._do_print(self[0], self[1], self[2], self[3])

    cpdef _richcmp_(self, right, int op):
        """
        Comparing elements.

        TESTS::

            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x']),-5,-2)
            sage: a = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: a == a
            True
            sage: a == a + 1
            False
            sage: a < a + 1
            True

            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x); a = x + 2*x*i + 3*j + (x-2)*k
            sage: (1/a)*a == 1
            True
            sage: (1/a)*a
            1
        """
        cdef int i
        for i in range(4):
            res = richcmp_item(self[i], right[i], op)
            if res is not NotImplemented:
                return res
        return rich_to_bool(op, 0)

    cpdef conjugate(self):
        """
        Return the conjugate of the quaternion: if `\\theta = x + yi + zj + wk`,
        return `x - yi - zj - wk`; that is, return theta.reduced_trace() - theta.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: a = 3*i - j + 2
            sage: type(a)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: a.conjugate()
            2 - 3*i + j

        The "universal" test::

            sage: K.<x,y,z,w,a,b> = QQ[]
            sage: Q.<i,j,k> = QuaternionAlgebra(a,b)
            sage: theta = x+y*i+z*j+w*k
            sage: theta.conjugate()
            x + (-y)*i + (-z)*j + (-w)*k
        """
        return self.__class__(self._parent, (self[0], -self[1], -self[2], -self[3]), check=False)

    cpdef reduced_trace(self):
        """
        Return the reduced trace of self: if `\\theta = x + yi + zj +
        wk`, then `\\theta` has reduced trace `2x`.

        EXAMPLES::

            sage: K.<x,y,z,w,a,b> = QQ[]
            sage: Q.<i,j,k> = QuaternionAlgebra(a,b)
            sage: theta = x+y*i+z*j+w*k
            sage: theta.reduced_trace()
            2*x
        """
        return 2*self[0]

    cpdef reduced_norm(self):
        """
        Return the reduced norm of self: if `\\theta = x + yi + zj +
        wk`, then `\\theta` has reduced norm `x^2 - ay^2 - bz^2 +
        abw^2`.

        EXAMPLES::

            sage: K.<x,y,z,w,a,b> = QQ[]
            sage: Q.<i,j,k> = QuaternionAlgebra(a,b)
            sage: theta = x+y*i+z*j+w*k
            sage: theta.reduced_norm()
            w^2*a*b - y^2*a - z^2*b + x^2
        """
        a, b = self._parent._a, self._parent._b
        x, y, z, w = self[0], self[1], self[2], self[3]
        return w*w*a*b - y*y*a - z*z*b + x*x

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-7,-13)
            sage: theta = 1/3 - 2/3*i + 4/19*j - 17/3*k
            sage: (1/theta) * theta
            1
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: 1/Q(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        Note that the quaternion algebra need not be a division
        algebra, in which case we can get a ZeroDivisionException::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,4,9)
            sage: theta = 2-i
            sage: theta.reduced_norm()
            0
            sage: 1/theta
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        The ``universal`` test:

            sage: K.<x,y,z,w,a,b> = QQ[]
            sage: Q.<i,j,k> = QuaternionAlgebra(a,b)
            sage: theta = x+y*i+z*j+w*k
            sage: 1/theta == theta.conjugate()/theta.reduced_norm()
            True
        """
        return ~self.reduced_norm() * self.conjugate()

    cpdef _rmul_(self, Element left):
        """
        Return left*self, where left is in the base ring.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2-x-1); Q.<i,j,k> = QuaternionAlgebra(K,-1,-1); z=2*i+3*j+4/3*k+5/8
            sage: a*z
            5/8*a + 2*a*i + 3*a*j + 4/3*a*k
            sage: type(z)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
        """
        return self.__class__(self._parent, (left*self[0], left*self[1], left*self[2], left*self[3]), check=False)

    cpdef _lmul_(self, Element right):
        """
        Return self*right, where right is in the base ring.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2-x-1); Q.<i,j,k> = QuaternionAlgebra(K,-1,-1); z=2*i+3*j+4/3*k+5/8
            sage: z*a
            5/8*a + 2*a*i + 3*a*j + 4/3*a*k
            sage: type(z)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
        """
        return self.__class__(self._parent, (self[0]*right, self[1]*right, self[2]*right, self[3]*right), check=False)

    cpdef _div_(self, right):
        """
        Return quotient of ``self`` by ``right``.

        EXAMPLES::

            sage: K.<x> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(x, 2*x)
            sage: theta = x + 2*x*i + 3*j + (x-2)*k
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: theta._div_(theta)
            1
            sage: theta._div_(theta) == 1
            True
        """
        return self * ~right

    def reduced_characteristic_polynomial(self, var='x'):
        """
        Return the reduced characteristic polynomial of this
        quaternion algebra element, which is `X^2 - tX + n`, where `t`
        is the reduced trace and `n` is the reduced norm.

        INPUT:

        - ``var`` -- string (default: ``'x'``); indeterminate of characteristic
          polynomial

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: i.reduced_characteristic_polynomial()
            x^2 + 1
            sage: j.reduced_characteristic_polynomial()
            x^2 + 2
            sage: (i+j).reduced_characteristic_polynomial()
            x^2 + 3
            sage: (2+j+k).reduced_trace()
            4
            sage: (2+j+k).reduced_characteristic_polynomial('T')
            T^2 - 4*T + 8
        """
        R = PolynomialRing(self.base_ring(), var)
        return R([self.reduced_norm(), -self.reduced_trace(), 1])

    def matrix(self, action='right'):
        """
        Return the matrix of right or left multiplication of ``self`` on
        the basis for the ambient quaternion algebra.

        In particular, if action is ``'right'`` (the default), returns the
        matrix of the mapping sending ``x`` to ``x*self``.

        INPUT:

        - ``action`` -- (default: ``'right'``) ``'right'`` or ``'left'``

        OUTPUT: a matrix

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(-3,-19)
            sage: a = 2/3 -1/2*i + 3/5*j - 4/3*k
            sage: a.matrix()
            [  2/3  -1/2   3/5  -4/3]
            [  3/2   2/3     4   3/5]
            [-57/5 -76/3   2/3   1/2]
            [   76 -57/5  -3/2   2/3]
            sage: a.matrix() == a.matrix(action='right')
            True
            sage: a.matrix(action='left')
            [  2/3  -1/2   3/5  -4/3]
            [  3/2   2/3    -4  -3/5]
            [-57/5  76/3   2/3  -1/2]
            [   76  57/5   3/2   2/3]
            sage: (i*a,j*a,k*a)
            (3/2 + 2/3*i + 4*j + 3/5*k, -57/5 - 76/3*i + 2/3*j + 1/2*k, 76 - 57/5*i - 3/2*j + 2/3*k)
            sage: a.matrix(action='foo')
            Traceback (most recent call last):
            ...
            ValueError: action must be either 'left' or 'right'

        We test over a more generic base field::

            sage: K.<x> = QQ['x']
            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K),-5,-2)
            sage: a = 1/2*x^2 + 2/3*x*i - 3/4*j + 5/7*k
            sage: type(a)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.matrix()
            [1/2*x^2   2/3*x    -3/4     5/7]
            [-10/3*x 1/2*x^2   -25/7    -3/4]
            [    3/2    10/7 1/2*x^2  -2/3*x]
            [  -50/7     3/2  10/3*x 1/2*x^2]
        """
        if action == 'right':
            v = [(a*self).coefficient_tuple() for a in self._parent.basis()]
        elif action == 'left':
            v = [(self*a).coefficient_tuple() for a in self._parent.basis()]
        else:
            raise ValueError("action must be either 'left' or 'right'")
        return matrix(self.base_ring(), 4, v)

    def coefficient_tuple(self):
        """
        Return 4-tuple of coefficients of this quaternion.

        EXAMPLES::

            sage: K.<x> = QQ['x']
            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K),-5,-2)
            sage: a = 1/2*x^2 + 2/3*x*i - 3/4*j + 5/7*k
            sage: type(a)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: a.coefficient_tuple()
            (1/2*x^2, 2/3*x, -3/4, 5/7)
        """
        return (self[0], self[1], self[2], self[3])

    def pair(self, right):
        """
        Return the result of pairing ``self`` and ``right``, which should both
        be elements of a quaternion algebra.  The pairing is
        ``(x,y) = (x.conjugate()*y).reduced_trace()``.

        INPUT:

        - ``right`` -- quaternion

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: (1+i+j-2*k).pair(2/3+5*i-3*j+k)
            -26/3
            sage: x = 1+i+j-2*k; y = 2/3+5*i-3*j+k
            sage: x.pair(y)
            -26/3
            sage: y.pair(x)
            -26/3
            sage: (x.conjugate()*y).reduced_trace()
            -26/3
        """
        return (self.conjugate() * right).reduced_trace()

    def _im_gens_(self, codomain, im_gens, base_map):
        r"""
        Return the image of this quaternion under the morphism
        defined by ``im_gens`` in ``codomain``, where elements
        of the base ring are mapped by ``base_map``.

        EXAMPLES::

            sage: Quat.<i,j,k> = QuaternionAlgebra(-1, -19)
            sage: b = 5 + 6*i + 7*j + 8*k
            sage: ims = [~b * g * b for g in Quat.gens()]
            sage: Quat(1)._im_gens_(Quat, ims, None) == 1
            True
            sage: i._im_gens_(Quat, ims, None) == ims[0]
            True
            sage: j._im_gens_(Quat, ims, None) == ims[1]
            True
            sage: k._im_gens_(Quat, ims, None) == ims[2]
            True
        """
        if base_map is None:
            base_map = lambda v: v
        return sum(base_map(c)*g for c,g in zip(self, [1] + list(im_gens)))


cdef class QuaternionAlgebraElement_generic(QuaternionAlgebraElement_abstract):
    """
    TESTS:

    Test operations on quaternions over a base ring that is not a field::

        sage: A.<t> = LaurentPolynomialRing(GF(3))
        sage: B = QuaternionAlgebra(A, -1, t)
        sage: i, j, k = B.gens()
        sage: i*j
        k
        sage: (j + k).reduced_norm()
        t

    Inverting an element is currently only possible if its reduced
    norm is a unit::

        sage: ~k
        (t^-1)*k
        sage: ~(i + j)
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for *: 'Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 3' and 'Quaternion Algebra (2, t) with base ring Univariate Laurent Polynomial Ring in t over Finite Field of size 3'

    We test pickling::

        sage: R.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(R,-5*x,-2)
        sage: theta = x + i*x^3 + j*x^2 + k*x
        sage: theta == loads(dumps(theta))
        True
    """
    def __init__(self, parent, v, bint check=True):
        """
        Create a quaternion over a general base ring.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic(Q, (x,1,-7,2/3*x^3))
            x + i + (-7)*j + 2/3*x^3*k
        """
        self._parent = parent
        if check:
            self.x, self.y, self.z, self.w = to_quaternion(parent._base, v)
        else:
            self.x, self.y, self.z, self.w = v

    def __getitem__(self, int i):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(QQ['x']),-5,-2)
            sage: theta = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: list(theta)
            [1/2, 2/3, -3/4, 5/7]
        """
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        elif i == 3:
            return self.w
        else:
            raise IndexError("quaternion element index out of range")

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: theta = 1/x + x*i - (x+1)*j + 2/(3*x^3+5)*k
            sage: loads(dumps(theta)) == theta
            True
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
        """
        return (unpickle_QuaternionAlgebraElement_generic_v0,
                (self._parent, (self.x, self.y, self.z, self.w)))

    cpdef _add_(self, _right):
        """
        Return the sum of ``self`` and ``_right``.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: (x+i+j+x^3*k) + (x-i-j+ (2/3*x^3+x)*k)              # indirect doctest
            2*x + (5/3*x^3 + x)*k
            sage: type(i)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
        """
        cdef QuaternionAlgebraElement_generic right = _right
        # TODO -- make this, etc. use __new__
        return QuaternionAlgebraElement_generic(self._parent, (self.x + right.x, self.y + right.y, self.z + right.z, self.w + right.w), check=False)

    cpdef _sub_(self, _right):
        """
        Return the difference of ``self`` and ``_right``.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: type(i)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: (x+i+j+x^3*k)._sub_(x-i-j+ (2/3*x^3+x)*k)
            2*i + 2*j + (1/3*x^3 - x)*k
        """
        cdef QuaternionAlgebraElement_generic right = _right
        return QuaternionAlgebraElement_generic(self._parent, (self.x - right.x, self.y - right.y, self.z - right.z, self.w - right.w), check=False)

    cpdef _mul_(self, _right):
        """
        Return the product of ``self`` and ``_right``.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: type(i)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: (x+i+j+x^3*k)._mul_(x-i-j+ (2/3*x^3+x)*k)
            -20/3*x^6 - 10*x^4 + x^2 + 7 + (10/3*x^3 + 2*x)*i + (-25/3*x^3 - 5*x)*j + (5/3*x^4 + x^2)*k
        """
        cdef QuaternionAlgebraElement_generic right = _right

        a = self._parent._a
        b = self._parent._b

        x1, y1, z1, w1 = self.x, self.y, self.z, self.w
        x2, y2, z2, w2 = right.x, right.y, right.z, right.w

        # x = x1*x2 + y1*y2*a + z1*z2*b - w1*w2*a*b
        # y = x1*y2 + y1*x2 - z1*w2*b + w1*z2*b
        # z = x1*z2 + y1*w2 + z1*x2 - w1*y2*a
        # w = x1*w2 + y1*z2 - z1*y2 + w1*x2

        t1 = x1 * x2
        t2 = y1 * y2
        t3 = z1 * z2
        t4 = w1 * w2
        t5 = x2 * z1
        t6 = y2 * w1
        t7 = x1 * z2
        t8 = y1 * w2

        x = t1 + a * t2 + b * (t3 - a*t4)
        y = (x1 + y1)*(x2 + y2) - t1 - t2 + b*((z1 + w1)*(z2 - w2) - t3 + t4)
        z = t5 - a*t6 + t7 + a*t8
        w = (x2 - y2)*(z1 + w1) - t5 + t6 + (x1 + y1)*(z2 + w2) - t7 - t8

        return QuaternionAlgebraElement_generic(self._parent, (x, y, z, w), check=False)

    def _repr_(self):
        """
        Print representation of ``self``.

        EXAMPLES::

            sage: K.<x> = Frac(QQ['x']); Q.<i,j,k> = QuaternionAlgebra(K,-5,-2)
            sage: theta = 1/x + x*i - (x+1)*j + 2/(3*x^3+5)*k
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: theta._repr_()
            '1/x + x*i + (-x - 1)*j + (2/3/(x^3 + 5/3))*k'
        """
        return self._do_print(self.x, self.y, self.z, self.w)


cdef class QuaternionAlgebraElement_rational_field(QuaternionAlgebraElement_abstract):
    """
    TESTS:

    We test pickling::

        sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
        sage: i + j + k == loads(dumps(i+j+k))
        True
    """
    # Implementation Notes:
    #
    # A Quaternion algebra element (call it a) over Q are implemented as a 4-tuple of
    # integers x, y, z, w and a denominator d, all of type mpz_t, such that
    #
    #           theta = (1/d) * (x + y * i + z * j + w * k)
    #
    # (although different letters may be specified instead of i, j, and k, if desired).
    #
    # Inside the element we also store mpz_t integers a and b, where
    #
    #       i^2 = a   and   j^2 = b

    def __cinit__(self):
        """
        Initialize C variables.

        EXAMPLES::

            sage: QuaternionAlgebra(QQ,-5,-2)([1/2,-1/3,2/3,4/5])  # implicit doctest
            1/2 - 1/3*i + 2/3*j + 4/5*k
        """
        mpz_init(self.x)
        mpz_init(self.y)
        mpz_init(self.z)
        mpz_init(self.w)
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.d)

    def __dealloc__(self):
        mpz_clear(self.x)
        mpz_clear(self.y)
        mpz_clear(self.z)
        mpz_clear(self.w)
        mpz_clear(self.a)
        mpz_clear(self.b)
        mpz_clear(self.d)

    def _rational_(self):
        """
        Try to coerce this quaternion to a Rational.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(-1,-2)
            sage: Rational(A(-2/3))                             # indirect doctest
            -2/3
            sage: Rational(i)
            Traceback (most recent call last):
            ...
            TypeError
        """
        cdef Rational x = Rational()
        if self.is_constant():
            mpq_set_num(x.value, self.x)
            mpq_set_den(x.value, self.d)
            mpq_canonicalize(x.value)
            return x
        raise TypeError

    cpdef bint is_constant(self) noexcept:
        """
        Return ``True`` if this quaternion is constant, i.e., has no `i`, `j`,
        or `k` term.

        OUTPUT: boolean

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: A(1/3).is_constant()
            True
            sage: A(-1).is_constant()
            True
            sage: (1+i).is_constant()
            False
            sage: j.is_constant()
            False
        """
        return not (mpz_sgn(self.y) or mpz_sgn(self.z) or mpz_sgn(self.w))

    def __bool__(self):
        """
        Return ``True`` if this quaternion is nonzero.

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: bool(1+j+k)
            True
            sage: bool(A(0))
            False
        """
        return bool(mpz_sgn(self.x) or mpz_sgn(self.y) or mpz_sgn(self.z) or mpz_sgn(self.w))

    cpdef _richcmp_(self, _right, int op):
        """
        Compare two quaternions.

        The comparison is fairly arbitrary
        -- first the denominators are compared and if equal then each
        of the other coefficients are compared.

        TESTS::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: i < j
            False
            sage: -i < j
            True
            sage: i == i
            True
            sage: Q.one() != -Q.one()
            True
        """
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef int i
        i = mpz_cmp(self.d, right.d)
        if i:
            return rich_to_bool_sgn(op, i)
        i = mpz_cmp(self.x, right.x)
        if i:
            return rich_to_bool_sgn(op, i)
        i = mpz_cmp(self.y, right.y)
        if i:
            return rich_to_bool_sgn(op, i)
        i = mpz_cmp(self.z, right.z)
        if i:
            return rich_to_bool_sgn(op, i)
        i = mpz_cmp(self.w, right.w)
        if i:
            return rich_to_bool_sgn(op, i)
        return rich_to_bool(op, 0)

    def __init__(self, parent, v, bint check=True):
        """
        Setup element data from parent and coordinates.

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-4,-5)
            sage: A(2/3)
            2/3
            sage: type(A(2/3))
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>

            sage: A([-1/2,-10/3,-2/3,-4/5])     # implicit doctest
            -1/2 - 10/3*i - 2/3*j - 4/5*k
            sage: A(vector([1,2/3,3/4,4/5]))
            1 + 2/3*i + 3/4*j + 4/5*k

        ::

            sage: QA = QuaternionAlgebra(QQ, -1, -1)
            sage: foo = QA(3.0); foo
            3
            sage: parent(foo)
            Quaternion Algebra (-1, -1) with base ring Rational Field
            sage: foo[0]
            3
            sage: parent(foo[0])
            Rational Field
        """
        self._parent = parent

        # cache a and b
        mpz_set(self.a, (<Integer>parent._a).value)
        mpz_set(self.b, (<Integer>parent._b).value)

        cdef Rational x, y, z, w
        cdef mpz_t lcm
        cdef mpq_t lcm_rat
        if not isinstance(v, (list, tuple)):
            try:
                x = Rational(v)
                mpz_set(self.x, mpq_numref(x.value))
                mpz_set_si(self.y, 0)
                mpz_set_si(self.z, 0)
                mpz_set_si(self.w, 0)
                mpz_set(self.d, mpq_denref(x.value))
                return
            except TypeError:
                pass
        if check:
            v = tuple(v)
            # Now v is definitely a list or tuple, and we convert each
            # entry to a rational, then clear denominators, etc.
            x = Rational(v[0])
            y = Rational(v[1])
            z = Rational(v[2])
            w = Rational(v[3])
        else:
            x, y, z, w = v
        mpz_init(lcm)
        mpz_lcm(lcm, mpq_denref(x.value), mpq_denref(y.value))
        mpz_lcm(lcm, lcm, mpq_denref(z.value))
        mpz_lcm(lcm, lcm, mpq_denref(w.value))
        if mpz_cmp_si(lcm, 1):
            mpz_init_set(mpq_numref(lcm_rat), lcm)
            mpz_init_set_si(mpq_denref(lcm_rat), 1)
            mpq_mul(x.value, x.value, lcm_rat)
            mpq_mul(y.value, y.value, lcm_rat)
            mpq_mul(z.value, z.value, lcm_rat)
            mpq_mul(w.value, w.value, lcm_rat)
            mpq_clear(lcm_rat)
        mpz_set(self.x, mpq_numref(x.value))
        mpz_set(self.y, mpq_numref(y.value))
        mpz_set(self.z, mpq_numref(z.value))
        mpz_set(self.w, mpq_numref(w.value))
        mpz_set(self.d, lcm)
        mpz_clear(lcm)

    def __getitem__(self, int i):
        """
        TESTS::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: theta = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: list(theta)
            [1/2, 2/3, -3/4, 5/7]
        """
        cdef Rational r = Rational()
        if i == 0:
            mpq_set_num(r.value, self.x)
        elif i == 1:
            mpq_set_num(r.value, self.y)
        elif i == 2:
            mpq_set_num(r.value, self.z)
        elif i == 3:
            mpq_set_num(r.value, self.w)
        else:
            raise IndexError("quaternion element index out of range")
        mpq_set_den(r.value, self.d)
        mpq_canonicalize(r.value)
        return r

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: K.<x> = QQ[]
            sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K),-5,-19)
            sage: theta = 1/2 + 2/3*i - 3/4*j + 5/7*k
            sage: type(theta)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: loads(dumps(theta)) == theta
            True
        """
        return (unpickle_QuaternionAlgebraElement_rational_field_v0,
                (self._parent, (self[0], self[1], self[2], self[3])))

    cpdef _add_(self, _right):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: type(i)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: (2/3 + 3/4*i + 5/6*j + 7/8*k)._add_(-2/3 - 3/4*i + 5/6*j + 7/8*k)
            5/3*j + 7/4*k
        """

        #   Given two quaternion algebra elements
        #       theta = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #       nu = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #   we compute their sum as
        #
        #   (theta + nu) = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #   with    d3 = d1 * d2
        #           x3 = d1 * x2 + d2 * x1
        #           y3 = d1 * y2 + d2 * y1
        #           z3 = d1 * z2 + d2 * z1
        #           w3 = d1 * w2 + d2 * w1
        #
        #   and then we reduce the sum by dividing everything
        #   by the gcd of d3, x3, y3, z3, and w3

        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_mul(U1, self.x, right.d)        # U1 = x1 * d2
        mpz_mul(U2, right.x, self.d)        # U2 = x2 * d1
        mpz_add(result.x, U1, U2)           # x3 = x1 * d2 + x2 * d1

        mpz_mul(U1, self.y, right.d)        # U1 = y1 * d2
        mpz_mul(U2, right.y, self.d)        # U2 = y2 * d1
        mpz_add(result.y, U1, U2)           # x3 = y1 * d2 + y2 * d1

        mpz_mul(U1, self.z, right.d)        # U1 = z1 * d2
        mpz_mul(U2, right.z, self.d)        # U2 = z2 * d1
        mpz_add(result.z, U1, U2)           # z3 = z1 * d2 + z2 * d1

        mpz_mul(U1, self.w, right.d)        # U1 = w1 * d2
        mpz_mul(U2, right.w, self.d)        # U2 = w2 * d1
        mpz_add(result.w, U1, U2)           # w3 = w1 * d2 + w2 * d1

        mpz_mul(result.d, self.d, right.d)  # d3 = d1 * d2

        result.canonicalize()

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        return result

    cpdef _sub_(self, _right):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: type(i)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: (2/3 + 3/4*i + 5/6*j + 7/8*k)._sub_(-2/3 - 3/4*i + 5/6*j + 7/8*k)
            4/3 + 3/2*i
        """
        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        # Implementation Note: To obtain _sub_, we simply replace every occurrence of
        # "add" in _add_ with "sub"; that is, we s/add/sub to get _sub_

        mpz_mul(U1, self.x, right.d)
        mpz_mul(U2, right.x, self.d)
        mpz_sub(result.x, U1, U2)

        mpz_mul(U1, self.y, right.d)
        mpz_mul(U2, right.y, self.d)
        mpz_sub(result.y, U1, U2)

        mpz_mul(U1, self.z, right.d)
        mpz_mul(U2, right.z, self.d)
        mpz_sub(result.z, U1, U2)

        mpz_mul(U1, self.w, right.d)
        mpz_mul(U2, right.w, self.d)
        mpz_sub(result.w, U1, U2)

        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        return result

    cpdef _mul_(self, _right):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: type(i)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: (2/3 + 3/4*i + 5/6*j + 7/8*k)._mul_(-2/3 - 3/4*i + 5/6*j + 7/8*k)
            9331/576 - i - 63/16*j + 5/4*k
        """

        # We use the following formula for multiplication:
        #
        #    Given two quaternion algebra elements
        #
        #        theta = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #        nu = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #    we compute their product as
        #
        #    theta*nu = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #    where
        #       x3 = t1 + a * t2 + b * (t3 - a*t4)
        #       y3 = s1*(x2 + y2) - t1 - t2 + b*( s2*(z2 - w2) - t3 + t4)
        #       z3 = t5 - a*t6 + t7 + a*t8
        #       w3 = (x2 - y2)*s2 - t5 + t6 + s1*(z2 + w2) - t7 - t8
        #
        #       and where
        #           t1 = x1 * x2
        #           t2 = y1 * y2
        #           t3 = z1 * z2
        #           t4 = w1 * w2
        #           t5 = x2 * z1
        #           t6 = y2 * w1
        #           t7 = x1 * z2
        #           t8 = y1 * w2
        #
        #           s1 = x1 + y1
        #           s2 = z1 + w1
        #
        # This takes more integer addition operations but fewer integer multiplication
        # operations than the "straightforward" multiplication method.
        #
        # There might be a way to optimize this formula further.

        cdef QuaternionAlgebraElement_rational_field right = _right
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)

        mpz_mul(T1, self.x, right.x)    # t1 = x1 * x2
        mpz_mul(T2, self.y, right.y)    # t2 = y1 * y2
        mpz_mul(t3, self.z, right.z)    # t3 = z1 * z2
        mpz_mul(t4, self.w, right.w)    # t4 = w1 * w2
        mpz_mul(t5, right.x, self.z)    # t5 = x2 * z1
        mpz_mul(t6, right.y, self.w)    # t6 = y2 * w1
        mpz_mul(t7, self.x, right.z)    # t7 = x1 * z2
        mpz_mul(t8, self.y, right.w)    # t8 = y1 * w2

        mpz_add(s1, self.x, self.y)     # s1 = x1 + y1
        mpz_add(s2, self.z, self.w)     # s2 = z1 + w1

        # ------------------

        mpz_mul(U1, self.a, t4)
        mpz_sub(U1, t3, U1)
        mpz_mul(U1, U1, self.b)
        mpz_mul(U2, self.a, T2)
        mpz_add(result.x, T1, U2)
        mpz_add(result.x, result.x, U1)

        # ------------------

        mpz_sub(U1, right.z, right.w)
        mpz_mul(U1, U1, s2)
        mpz_sub(U1, U1, t3)
        mpz_add(U1, U1, t4)
        mpz_mul(U1, U1, self.b)
        mpz_sub(U1, U1, T2)
        mpz_sub(U1, U1, T1)
        mpz_add(U2, right.x, right.y)
        mpz_mul(U2, s1, U2)
        mpz_add(result.y, U1, U2)

        # ------------------

        mpz_mul(U1, self.a, t8)
        mpz_add(U1, U1, t7)
        mpz_mul(U2, self.a, t6)
        mpz_sub(U1, U1, U2)
        mpz_add(result.z, U1, t5)

        # ------------------

        mpz_add(U1, right.z, right.w)
        mpz_mul(U1, U1, s1)
        mpz_sub(U1, U1, t7)
        mpz_sub(U1, U1, t8)
        mpz_add(U1, U1, t6)
        mpz_sub(U1, U1, t5)
        mpz_sub(U2, right.x, right.y)
        mpz_mul(U2, U2, s2)
        mpz_add(result.w, U1, U2)

        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        return result

    cpdef reduced_norm(self):
        """
        Return the reduced norm of ``self``.

        Given a quaternion `x+yi+zj+wk`, this is `x^2 - ay^2 - bz^2 + abw^2`.

        EXAMPLES::

            sage: K.<i,j,k> = QuaternionAlgebra(QQ, -5, -2)
            sage: i.reduced_norm()
            5
            sage: j.reduced_norm()
            2
            sage: a = 1/3 + 1/5*i + 1/7*j + k
            sage: a.reduced_norm()
            22826/2205
        """
        mpz_mul(U1, self.x, self.x)         # U1 = x*x
        mpz_mul(U2, self.b, self.z)         # U2 = b*x
        mpz_mul(U2, U2, self.z)             # U2 = b*z*z
        mpz_sub(U2, U1, U2)                 # U2 = -b*z*z + x*x

        mpz_mul(U1, self.y, self.a)         # U1 = a*y
        mpz_mul(U1, U1, self.y)             # U1 = a*y*y
        mpz_sub(U2, U2, U1)                 # U2 = -a*y*y - b*z*z + x*x

        mpz_mul(U1, self.w, self.w)         # U1 = w*w
        mpz_mul(U1, U1, self.a)             # U1 = w*w*a
        mpz_mul(U1, U1, self.b)             # U1 = w*w*a*b

        mpz_add(U1, U1, U2)                 # U1 = w*w*a*b - a*y*y - b*z*z + x*x

        mpz_mul(U2, self.d, self.d)

        cdef Rational result = Rational.__new__(Rational)
        mpq_set_num(result.value, U1)
        mpq_set_den(result.value, U2)
        mpq_canonicalize(result.value)

        return result

    cpdef conjugate(self):
        """
        Return the conjugate of this quaternion.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: a = 3*i - j + 2
            sage: type(a)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: a.conjugate()
            2 - 3*i + j
            sage: b = 1 + 1/3*i + 1/5*j - 1/7*k
            sage: b.conjugate()
            1 - 1/3*i - 1/5*j + 1/7*k
        """

        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)
        mpz_set(result.d, self.d)

        mpz_set(result.x, self.x)
        mpz_mul_si(result.y, self.y, -1)
        mpz_mul_si(result.z, self.z, -1)
        mpz_mul_si(result.w, self.w, -1)

        return result

    cpdef reduced_trace(self):
        """
        Return the reduced trace of ``self``.

        This is `2x` if ``self`` is `x+iy+zj+wk`.

        EXAMPLES::

            sage: K.<i,j,k> = QuaternionAlgebra(QQ, -5, -2)
            sage: i.reduced_trace()
            0
            sage: j.reduced_trace()
            0
            sage: a = 1/3 + 1/5*i + 1/7*j + k
            sage: a.reduced_trace()
            2/3
        """
        mpz_mul_si(U1, self.x, 2)
        cdef Rational result = Rational.__new__(Rational)
        mpq_set_num(result.value, U1)
        mpq_set_den(result.value, self.d)
        mpq_canonicalize(result.value)
        return result

    cdef inline canonicalize(self):
        """
        Put the representation of this quaternion element into
        smallest form. For `a = (1/d)(x + yi + zj + wk)` we
        divide `a`, `x`, `y`, `z`, and `w` by the gcd of all of them.

        TESTS::

            sage: K.<i,j,k> = QuaternionAlgebra(QQ, -10, -7)
            sage: (1/4 + 1/2 * i + 1/7 * j + 1/28 * k)*14*i     # implicit doctest
            -70 + 7/2*i + 5*j - 2*k
        """

        # Note: this function changes the module-level global variable
        # U1, so it isn't always safe to use this in the middle of
        # another function. Normally this function is called
        # at the end of an arithmetic routine, so this is fine.

        # Implementation-wise, we compute the GCD's one at a time,
        # and quit if it ever becomes one

        mpz_gcd(U1, self.d, self.x)
        if mpz_cmp_ui(U1, 1) != 0:
            mpz_gcd(U1, U1, self.y)
            if mpz_cmp_ui(U1, 1) != 0:
                mpz_gcd(U1, U1, self.z)
                if mpz_cmp_ui(U1, 1) != 0:
                    mpz_gcd(U1, U1, self.w)
                    if mpz_cmp_ui(U1, 1) != 0:
                        # at this point U1 actually contains the gcd of all the terms, and we divide
                        mpz_divexact(self.d, self.d, U1)
                        mpz_divexact(self.x, self.x, U1)
                        mpz_divexact(self.y, self.y, U1)
                        mpz_divexact(self.z, self.z, U1)
                        mpz_divexact(self.w, self.w, U1)

    def denominator(self):
        """
        Return the lowest common multiple of the denominators of the coefficients
        of i, j and k for this quaternion.

        EXAMPLES::

            sage: A = QuaternionAlgebra(QQ, -1, -1)
            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1, -1)
            sage: a = (1/2) + (1/5)*i + (5/12)*j + (1/13)*k
            sage: a
            1/2 + 1/5*i + 5/12*j + 1/13*k
            sage: a.denominator()
            780
            sage: lcm([2, 5, 12, 13])
            780
            sage: (a * a).denominator()
            608400
            sage: (a + a).denominator()
            390
        """
        cdef Integer d = Integer()
        mpz_set(d.value, self.d)
        return d

    def denominator_and_integer_coefficient_tuple(self):
        """
        Return 5-tuple d, x, y, z, w, where this rational quaternion
        is equal to `(x + yi + zj + wk)/d` and x, y, z, w do not share
        a common factor with d.

        OUTPUT: 5-tuple of Integers

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: (2 + 3*i + 4/3*j - 5*k).denominator_and_integer_coefficient_tuple()
            (3, 6, 9, 4, -15)
        """
        cdef Integer d = Integer()
        cdef Integer y = Integer()
        cdef Integer x = Integer()
        cdef Integer z = Integer()
        cdef Integer w = Integer()

        mpz_set(d.value, self.d)
        mpz_set(x.value, self.x)
        mpz_set(y.value, self.y)
        mpz_set(z.value, self.z)
        mpz_set(w.value, self.w)

        return (d, x, y, z, w)

    def integer_coefficient_tuple(self):
        """
        Return the integer part of this quaternion, ignoring the common denominator.

        OUTPUT: 4-tuple of Integers

        EXAMPLES::

            sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
            sage: (2 + 3*i + 4/3*j - 5*k).integer_coefficient_tuple()
            (6, 9, 4, -15)
        """
        cdef Integer y = Integer()
        cdef Integer x = Integer()
        cdef Integer z = Integer()
        cdef Integer w = Integer()

        mpz_set(x.value, self.x)
        mpz_set(y.value, self.y)
        mpz_set(z.value, self.z)
        mpz_set(w.value, self.w)

        return (x, y, z, w)

    def coefficient_tuple(self):
        """
        Return 4-tuple of rational numbers which are the coefficients of this quaternion.

        EXAMPLES::

            sage: A.<i,j,k> = QuaternionAlgebra(-1,-2)
            sage: (2/3 + 3/5*i + 4/3*j - 5/7*k).coefficient_tuple()
            (2/3, 3/5, 4/3, -5/7)
        """
        cdef Rational x = Rational()
        cdef Rational y = Rational()
        cdef Rational z = Rational()
        cdef Rational w = Rational()

        mpq_set_num(x.value, self.x)
        mpq_set_num(y.value, self.y)
        mpq_set_num(z.value, self.z)
        mpq_set_num(w.value, self.w)

        mpq_set_den(x.value, self.d)
        mpq_set_den(y.value, self.d)
        mpq_set_den(z.value, self.d)
        mpq_set_den(w.value, self.d)

        mpq_canonicalize(x.value)
        mpq_canonicalize(y.value)
        mpq_canonicalize(z.value)
        mpq_canonicalize(w.value)
        return (x, y, z, w)

    def _multiply_by_integer(self, Integer n):
        """
        Return the product of ``self`` times the integer `n`.

        EXAMPLES::

            sage: A = QuaternionAlgebra(7)
            sage: a = A.random_element()
            sage: 5*a == a._multiply_by_integer(5)
            True
        """
        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)

        if mpz_divisible_p(self.d, n.value) != 0:
            mpz_divexact(result.d, self.d, n.value)
            mpz_set(result.x, self.x)
            mpz_set(result.y, self.y)
            mpz_set(result.z, self.z)
            mpz_set(result.w, self.w)
            return result
        if mpz_divisible_p(n.value, self.d):
            mpz_divexact(T1, n.value, self.d)
        else:
            mpz_set(T1, n.value)

        mpz_set(result.d, self.d)
        mpz_mul(result.x, self.x, T1)
        mpz_mul(result.y, self.y, T1)
        mpz_mul(result.z, self.z, T1)
        mpz_mul(result.w, self.w, T1)

        return result

    def _divide_by_integer(self, Integer n):
        """
        Return the quotient of ``self`` by the integer `n`.

        EXAMPLES::

            sage: A = QuaternionAlgebra(7)
            sage: a = A.random_element()
            sage: a/5 == a._divide_by_integer(5)
            True
            sage: a._divide_by_integer(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        if mpz_sgn(n.value) == 0:
            raise ZeroDivisionError

        cdef QuaternionAlgebraElement_rational_field result = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        result._parent = self._parent

        mpz_set(result.a, self.a)
        mpz_set(result.b, self.b)

        mpz_mul(result.d, self.d, n.value)
        mpz_set(result.x, self.x)
        mpz_set(result.y, self.y)
        mpz_set(result.z, self.z)
        mpz_set(result.w, self.w)

        result.canonicalize()

        return result


cdef class QuaternionAlgebraElement_number_field(QuaternionAlgebraElement_abstract):
    def __cinit__(self):
        """
        Allocate memory for this quaternion over a number field.

        EXAMPLES::

            sage: K.<a> = QQ[2^(1/5)]; Q.<i,j,k> = QuaternionAlgebra(K,-a,a*17/3)       # needs sage.symbolic
            sage: Q([a,-2/3,a^2-1/2,a*2])           # implicit doctest                  # needs sage.symbolic
            a + (-2/3)*i + (a^2 - 1/2)*j + 2*a*k
        """
        fmpz_poly_init(self.x)
        fmpz_poly_init(self.y)
        fmpz_poly_init(self.z)
        fmpz_poly_init(self.w)
        fmpz_poly_init(self.a)
        fmpz_poly_init(self.b)
        fmpz_poly_init(self.modulus)
        mpz_init(self.d)

    def __dealloc__(self):
        """
        Free memory used by this quaternion over a number field.
        """
        fmpz_poly_clear(self.x)
        fmpz_poly_clear(self.y)
        fmpz_poly_clear(self.z)
        fmpz_poly_clear(self.w)
        fmpz_poly_clear(self.a)
        fmpz_poly_clear(self.b)
        fmpz_poly_clear(self.modulus)
        mpz_clear(self.d)

    def __init__(self, parent, v, bint check=True):
        """
        EXAMPLES::

            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K,-a,a+1)          # needs sage.symbolic
            sage: Q([a,-2/3,a^2-1/2,a*2])           # implicit doctest                  # needs sage.symbolic
            a + (-2/3)*i + (a^2 - 1/2)*j + 2*a*k
        """
        self._parent = parent
        if check:
            x, y, z, w = to_quaternion(parent._base, v)
        else:
            x, y, z, w = v
        cdef NumberFieldElement a = <NumberFieldElement>(parent._base(parent._a))
        cdef NumberFieldElement b = <NumberFieldElement>(parent._base(parent._b))
        fmpz_poly_set_ZZX(self.x, (<NumberFieldElement>x)._numerator)
        fmpz_poly_set_ZZX(self.y, (<NumberFieldElement>y)._numerator)
        fmpz_poly_set_ZZX(self.z, (<NumberFieldElement>z)._numerator)
        fmpz_poly_set_ZZX(self.w, (<NumberFieldElement>w)._numerator)

        ZZ_to_mpz(T1, &(<NumberFieldElement>x)._denominator)
        ZZ_to_mpz(T2, &(<NumberFieldElement>y)._denominator)
        ZZ_to_mpz(t3, &(<NumberFieldElement>z)._denominator)
        ZZ_to_mpz(t4, &(<NumberFieldElement>w)._denominator)

        mpz_lcm(self.d, T1, T2)
        mpz_lcm(self.d, self.d, t3)
        mpz_lcm(self.d, self.d, t4)

        mpz_divexact(T1, self.d, T1)
        mpz_divexact(T2, self.d, T2)
        mpz_divexact(t3, self.d, t3)
        mpz_divexact(t4, self.d, t4)

        fmpz_poly_scalar_mul_mpz(self.x, self.x, T1)
        fmpz_poly_scalar_mul_mpz(self.y, self.y, T2)
        fmpz_poly_scalar_mul_mpz(self.z, self.z, t3)
        fmpz_poly_scalar_mul_mpz(self.w, self.w, t4)

        fmpz_poly_set_ZZX(self.a, a._numerator)     # we will assume that the denominator of a and b are 1
        fmpz_poly_set_ZZX(self.b, b._numerator)

        fmpz_poly_set_ZZX(self.modulus, (<NumberFieldElement>x)._fld_numerator.x)  # and same for the modulus

    def __getitem__(self, int i):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K,-a,a+1)
            sage: Q([a,-2/3,a^2-1/2,a*2])
            a + (-2/3)*i + (a^2 - 1/2)*j + 2*a*k
            sage: x = Q([a,-2/3,a^2-1/2,a*2])
            sage: type(x)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: x[0]
            a
            sage: x[1]
            -2/3
            sage: x[2]
            a^2 - 1/2
            sage: x[3]
            2*a
            sage: list(x)
            [a, -2/3, a^2 - 1/2, 2*a]
        """
        # general number -- this code assumes that the number field is not quadratic!!

        cdef NumberFieldElement el = <NumberFieldElement>(self._parent.base_ring().an_element())
        cdef NumberFieldElement item = el._new()

        if i == 0:
            fmpz_poly_get_ZZX(item._numerator, self.x)
        elif i == 1:
            fmpz_poly_get_ZZX(item._numerator, self.y)
        elif i == 2:
            fmpz_poly_get_ZZX(item._numerator, self.z)
        elif i == 3:
            fmpz_poly_get_ZZX(item._numerator, self.w)
        else:
            raise IndexError("quaternion element index out of range")

        mpz_to_ZZ(&item._denominator, self.d)

        return item

    def __reduce__(self):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = (i+j+k+a)^2; z
            a^2 + 4*a - 3 + 2*a*i + 2*a*j + 2*a*k
            sage: f, t = z.__reduce__()
            sage: f(*t)
            a^2 + 4*a - 3 + 2*a*i + 2*a*j + 2*a*k
            sage: loads(dumps(z)) == z
            True
        """
        return (unpickle_QuaternionAlgebraElement_number_field_v0,
                (self._parent, (self[0], self[1], self[2], self[3])))

    cpdef _add_(self, _right):
        """
        Add ``self`` and ``_right``:

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = a + i + (2/3)*a^3*j + (1+a)*k
            sage: w = a - i - (2/3)*a^3*j + (1/3+a)*k
            sage: type(z)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: z._add_(w)
            2*a + (2*a + 4/3)*k

        Check that the fix in :issue:`17099` is correct::

            sage: x = polygen(QQ, 'x')
            sage: K = NumberField(x**3 + x - 1, 'a')
            sage: D.<i,j,k> = QuaternionAlgebra(K, -1, -3)
            sage: j/3 + (2*j)/3 == j
            True
        """

        #   Given two quaternion algebra elements
        #       a = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #       b = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #   we compute their sum as
        #
        #   (a + b) = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #   with    d3 = d1 * d2
        #           x3 = d1 * x2 + d2 * x1
        #           y3 = d1 * y2 + d2 * y1
        #           z3 = d1 * z2 + d2 * z1
        #           w3 = d1 * w2 + d2 * w1
        #
        #   and then we reduce the sum by calling canonicalize().

        # Note: We are assuming in this routine that the modulus is monic. This shouldn't
        # currently be an issue because it is impossible to create a number field with
        # a modulus that is not monic.

        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> QuaternionAlgebraElement_number_field.__new__(QuaternionAlgebraElement_number_field)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)
        result._parent = self._parent

        fmpz_poly_scalar_mul_mpz(fU1, self.x, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.x, self.d)
        fmpz_poly_add(result.x, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.y, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.y, self.d)
        fmpz_poly_add(result.y, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.w, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.w, self.d)
        fmpz_poly_add(result.w, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.z, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.z, self.d)
        fmpz_poly_add(result.z, fU1, fU2)

        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        return result

    cpdef _sub_(self, _right):
        """
        Subtract ``_right`` from ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = a + i + (2/3)*a^3*j + (1+a)*k
            sage: w = a - i - (2/3)*a^3*j + (1/3+a)*k
            sage: type(z)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: z._sub_(w)
            2*i + 8/3*j + 2/3*k
        """
        # Implementation Note: To obtain _sub_, we simply replace every occurrence of
        # "add" in _add_ with "sub"; that is, we s/add/sub to get _sub_

        # Note: We are assuming in this routine that the modulus is monic. This shouldn't
        # currently be an issue because it is impossible to create a number field with
        # a modulus that is not monic.
        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> QuaternionAlgebraElement_number_field.__new__(QuaternionAlgebraElement_number_field)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)
        result._parent = self._parent

        fmpz_poly_scalar_mul_mpz(fU1, self.x, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.x, self.d)
        fmpz_poly_sub(result.x, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.y, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.y, self.d)
        fmpz_poly_sub(result.y, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.w, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.w, self.d)
        fmpz_poly_sub(result.w, fU1, fU2)

        fmpz_poly_scalar_mul_mpz(fU1, self.z, right.d)
        fmpz_poly_scalar_mul_mpz(fU2, right.z, self.d)
        fmpz_poly_sub(result.z, fU1, fU2)

        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        return result

    cpdef _mul_(self, _right):
        """
        Multiply ``self`` and ``_right``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a)
            sage: z = a + i + (2/3)*a^3*j + (1+a)*k
            sage: w = a - i - (2/3)*a^3*j + (1/3+a)*k
            sage: type(z)
            <class 'sage.algebras.quatalg.quaternion_algebra_element.QuaternionAlgebraElement_number_field'>
            sage: z._mul_(w)
            5*a^2 - 7/9*a + 9 + (-8/3*a^2 - 16/9*a)*i + (-6*a - 4)*j + (2*a^2 + 4/3*a)*k
        """
        # We use the following formula for multiplication:
        #
        #    Given two quaternion algebra elements
        #
        #        a = (1/d1)*(x1 + y1 * i + z1 * j + w1 * k)
        #        b = (1/d2)*(x2 + y2 * i + z2 * j + w2 * k)
        #
        #    we compute their product as
        #
        #    ab = (1/d3)*(x3 + y3 * i + z3 * j + w3 * k)
        #
        #    where
        #       x3 = t1 + a * t2 + b * (t3 - a*t4)
        #       y3 = s1*(x2 + y2) - t1 - t2 + b*( s2*(z2 - w2) - t3 + t4)
        #       z3 = t5 - a*t6 + t7 + a*t8
        #       w3 = (x2 - y2)*s2 - t5 + t6 + s1*(z2 + w2) - t7 - t8
        #
        #       and where
        #           t1 = x1 * x2
        #           t2 = y1 * y2
        #           t3 = z1 * z2
        #           t4 = w1 * w2
        #           t5 = x2 * z1
        #           t6 = y2 * w1
        #           t7 = x1 * z2
        #           t8 = y1 * w2
        #
        #           s1 = x1 + y1
        #           s2 = z1 + w1
        #
        # This takes more polynomial addition operations but fewer polynomial multiplication
        # operations than the "straightforward" multiplication method.
        #
        # There might be a way to optimize this formula further.

        cdef QuaternionAlgebraElement_number_field right = _right
        cdef QuaternionAlgebraElement_number_field result = <QuaternionAlgebraElement_number_field> QuaternionAlgebraElement_number_field.__new__(QuaternionAlgebraElement_number_field)

        mpz_set_si(result.d, 1)

        fmpz_poly_set(result.a, self.a)
        fmpz_poly_set(result.b, self.b)
        fmpz_poly_set(result.modulus, self.modulus)
        result._parent = self._parent

        fmpz_poly_mul(fT1, self.x, right.x)    # t1 = x1 * x2
        fmpz_poly_mul(fT2, self.y, right.y)    # t2 = y1 * y2
        fmpz_poly_mul(ft3, self.z, right.z)    # t3 = x1 * x2
        fmpz_poly_mul(ft4, self.w, right.w)    # t4 = w1 * w2
        fmpz_poly_mul(ft5, right.x, self.z)    # t5 = x2 * z1
        fmpz_poly_mul(ft6, right.y, self.w)    # t6 = y2 * w1
        fmpz_poly_mul(ft7, self.x, right.z)    # t7 = x1 * z2
        fmpz_poly_mul(ft8, self.y, right.w)    # t8 = y1 * w2

        fmpz_poly_add(fs1, self.x, self.y)     # s1 = x1 + y1
        fmpz_poly_add(fs2, self.z, self.w)     # s2 = z1 + w

        # ------------------

        fmpz_poly_mul(fU1, self.a, ft4)
        fmpz_poly_sub(fU1, ft3, fU1)
        fmpz_poly_mul(fU1, fU1, self.b)
        fmpz_poly_mul(fU2, self.a, fT2)
        fmpz_poly_add(result.x, fT1, fU2)
        fmpz_poly_add(result.x, result.x, fU1)

        # ------------------

        fmpz_poly_sub(fU1, right.z, right.w)
        fmpz_poly_mul(fU1, fU1, fs2)
        fmpz_poly_sub(fU1, fU1, ft3)
        fmpz_poly_add(fU1, fU1, ft4)
        fmpz_poly_mul(fU1, fU1, self.b)
        fmpz_poly_sub(fU1, fU1, fT2)
        fmpz_poly_sub(fU1, fU1, fT1)
        fmpz_poly_add(fU2, right.x, right.y)
        fmpz_poly_mul(fU2, fs1, fU2)
        fmpz_poly_add(result.y, fU1, fU2)

        # ------------------

        fmpz_poly_mul(fU1, self.a, ft8)
        fmpz_poly_add(fU1, fU1, ft7)
        fmpz_poly_mul(fU2, self.a, ft6)
        fmpz_poly_sub(fU1, fU1, fU2)
        fmpz_poly_add(result.z, fU1, ft5)

        # ------------------

        fmpz_poly_add(fU1, right.z, right.w)
        fmpz_poly_mul(fU1, fU1, fs1)
        fmpz_poly_sub(fU1, fU1, ft7)
        fmpz_poly_sub(fU1, fU1, ft8)
        fmpz_poly_add(fU1, fU1, ft6)
        fmpz_poly_sub(fU1, fU1, ft5)
        fmpz_poly_sub(fU2, right.x, right.y)
        fmpz_poly_mul(fU2, fU2, fs2)
        fmpz_poly_add(result.w, fU1, fU2)

        # At this point we have essentially computed the product, but we still
        # need to reduce modulo the modulus, which is what the following 12 lines do.
        #
        # When this was written, the version of flint in Sage had problems with
        # fpmz_poly_divrem(). This should be fixed in the newest version of
        # flint, which also should have some new functions which should do
        # this faster (Bill Hart sent an email to Bober and William about this).
        #
        # This should be fixed in the near future. (I don't know how much
        # faster it will be when it is updated, but the following code is
        # currently quite a bottleneck.

        fmpz_poly_div(fT1, result.x, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.x, result.x, fT1)

        fmpz_poly_div(fT1, result.y, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.y, result.y, fT1)

        fmpz_poly_div(fT1, result.z, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.z, result.z, fT1)

        fmpz_poly_div(fT1, result.w, result.modulus)
        fmpz_poly_mul(fT1, fT1, result.modulus)
        fmpz_poly_sub(result.w, result.w, fT1)

        mpz_mul(result.d, self.d, right.d)

        result.canonicalize()

        return result

    cdef inline canonicalize(self):
        """
        Put the representation of this quaternion element into
        smallest form. For a = `(1/d)(x + yi + zj + wk)` we
        divide `a`, `x`, `y`, `z`, and `w` by the gcd of all of them.

        TESTS::

            sage: # needs sage.symbolic
            sage: F = QQ[3^(1/3)]
            sage: a = F.gen()
            sage: K.<i,j,k> = QuaternionAlgebra(F, -10 + a, -7 - a)
            sage: ((1/4 + 1/2 * i + a^3/7 * j + a/28 * k)*14*i)^3   # implicit doctest
            34503/2*a^2 + 132195/2*a + 791399/4 + (203/8*a^2 - 10591*a + 169225/4)*i
             + (-84695/4*a^2 + 483413/8*a + 18591/4)*j + (-87/2*a^2 + 18156*a - 72525)*k
        """

        # Note: this function changes the module-level global variables
        # U1 and U2, so it isn't always safe to use this in the middle of
        # another function. Normally this function is called
        # at the end of an arithmetic routine, so this is fine.

        # Implementation-wise, we compute the GCD's one at a time,
        # and quit if it ever becomes one

        cdef fmpz_t content
        fmpz_init(content)
        fmpz_poly_content(content, self.x)
        fmpz_get_mpz(U1, content)
        mpz_gcd(U1, self.d, U1)
        if mpz_cmp_ui(U1, 1) != 0:
            fmpz_poly_content(content, self.y)
            fmpz_get_mpz(U2, content)
            mpz_gcd(U1, U1, U2)
            if mpz_cmp_ui(U1, 1) != 0:
                fmpz_poly_content(content, self.z)
                fmpz_get_mpz(U2, content)
                mpz_gcd(U1, U1, U2)
                if mpz_cmp_ui(U1, 1) != 0:
                    fmpz_poly_content(content, self.w)
                    fmpz_get_mpz(U2, content)
                    mpz_gcd(U1, U1, U2)
                    if mpz_cmp_ui(U1, 1) != 0:
                        fmpz_poly_scalar_divexact_mpz(self.x, self.x, U1)
                        fmpz_poly_scalar_divexact_mpz(self.y, self.y, U1)
                        fmpz_poly_scalar_divexact_mpz(self.z, self.z, U1)
                        fmpz_poly_scalar_divexact_mpz(self.w, self.w, U1)
                        mpz_divexact(self.d, self.d, U1)

        fmpz_clear(content)


#######################################################################
# Versioned unpickle functions
#######################################################################

def unpickle_QuaternionAlgebraElement_generic_v0(*args):
    """
    EXAMPLES::

        sage: K.<X> = QQ[]
        sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K), -5,-19); z = 2/3 + i*X - X^2*j + X^3*k
        sage: f, t = z.__reduce__()
        sage: sage.algebras.quatalg.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_generic_v0(*t)
        2/3 + X*i + (-X^2)*j + X^3*k
        sage: sage.algebras.quatalg.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_generic_v0(*t) == z
        True
    """
    return QuaternionAlgebraElement_generic(*args, check=False)


def unpickle_QuaternionAlgebraElement_rational_field_v0(*args):
    """
    EXAMPLES::

        sage: Q.<i,j,k> = QuaternionAlgebra(-5,-19); a = 2/3 + i*5/7 - j*2/5 +19/2
        sage: f, t = a.__reduce__()
        sage: sage.algebras.quatalg.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_rational_field_v0(*t)
        61/6 + 5/7*i - 2/5*j
    """
    return QuaternionAlgebraElement_rational_field(*args, check=False)


def unpickle_QuaternionAlgebraElement_number_field_v0(*args):
    """
    EXAMPLES::

        sage: # needs sage.symbolic
        sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a); z = i + j
        sage: f, t = z.__reduce__()
        sage: sage.algebras.quatalg.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_number_field_v0(*t)
        i + j
        sage: sage.algebras.quatalg.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_number_field_v0(*t) == z
        True
    """
    return QuaternionAlgebraElement_number_field(*args, check=False)
