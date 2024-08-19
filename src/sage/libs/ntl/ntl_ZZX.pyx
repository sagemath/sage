# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

from cysignals.signals cimport sig_on, sig_off

from sage.cpython.string cimport char_to_str
from sage.ext.cplusplus cimport ccreadstr

include "decl.pxi"
include 'misc.pxi'

from cpython.object cimport Py_EQ, Py_NE
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ import unpickle_class_value

from sage.rings.integer_ring import IntegerRing
from sage.arith.power cimport generic_power_pos

ZZ = IntegerRing()

cdef inline ntl_ZZ make_ZZ(ZZ_c* x):
    """ These make_XXXX functions are deprecated and should be phased out."""
    cdef ntl_ZZ y
    y = ntl_ZZ()
    y.x = x[0]
    del x
    return y

# You must do sig_on() before calling this function
cdef inline ntl_ZZ make_ZZ_sig_off(ZZ_c* x):
    cdef ntl_ZZ y = make_ZZ(x)
    sig_off()
    return y

cdef inline ntl_ZZX make_ZZX(ZZX_c* x):
    """ These make_XXXX functions are deprecated and should be phased out."""
    cdef ntl_ZZX y
    y = ntl_ZZX()
    y.x = x[0]
    del x
    return y

# You must do sig_on() before calling this function
cdef inline ntl_ZZX make_ZZX_sig_off(ZZX_c* x):
    cdef ntl_ZZX y = make_ZZX(x)
    sig_off()
    return y

from sage.structure.proof.proof import get_flag
cdef proof_flag(t):
    return get_flag(t, "polynomial")

##############################################################################
#
# ZZX: polynomials over the integers
#
##############################################################################


cdef class ntl_ZZX():
    r"""
    The class \class{ZZX} implements polynomials in `\Z[X]`, i.e.,
    univariate polynomials with integer coefficients.

    Polynomial multiplication is very fast, and is implemented using
    one of 4 different algorithms:
    \begin{enumerate}
    \item\hspace{1em} Classical
    \item\hspace{1em} Karatsuba
    \item\hspace{1em} Schoenhage and Strassen --- performs an FFT by working
          modulo a "Fermat number" of appropriate size...
          good for polynomials with huge coefficients
          and moderate degree
    \item\hspace{1em} CRT/FFT --- performs an FFT by working modulo several
         small primes.  This is good for polynomials with moderate
         coefficients and huge degree.
    \end{enumerate}

    The choice of algorithm is somewhat heuristic, and may not always be
    perfect.

    Many thanks to Juergen Gerhard {\tt
    <jngerhar@plato.uni-paderborn.de>} for pointing out the deficiency
    in the NTL-1.0 ZZX arithmetic, and for contributing the
    Schoenhage/Strassen code.

    Extensive use is made of modular algorithms to enhance performance
    (e.g., the GCD algorithm and many others).
    """

    # See ntl_ZZX.pxd for definition of data members
    def __init__(self, v=None):
        """
        EXAMPLES::

            sage: f = ntl.ZZX([1,2,5,-9])
            sage: f
            [1 2 5 -9]
            sage: g = ntl.ZZX([0,0,0]); g
            []
            sage: g[10]=5
            sage: g
            [0 0 0 0 0 0 0 0 0 0 5]
            sage: g[10]
            5
        """
        cdef ntl_ZZ cc
        cdef Py_ssize_t i

        if v is None:
            return
        elif isinstance(v, (list, tuple)):
            for i from 0 <= i < len(v):
                x = v[i]
                if not isinstance(x, ntl_ZZ):
                    cc = ntl_ZZ(x)
                else:
                    cc = x
                ZZX_SetCoeff(self.x, i, cc.x)
        else:
            ccreadstr(self.x, str(v))

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.ntl.ntl_ZZX import ntl_ZZX
            sage: f = ntl_ZZX([1,2,0,4])
            sage: loads(dumps(f)) == f
            True
        """
        return unpickle_class_value, (ntl_ZZX, self.list())

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ntl.ZZX([1,3,0,5]).__repr__()
            '[1 3 0 5]'
        """
        cdef char * val
        val = ZZX_repr(&self.x)
        result = char_to_str(val)
        cpp_delete_array(val)
        return result

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: x = ntl.ZZX([1,32])
            sage: y = copy(x)
            sage: y == x
            True
            sage: y is x
            False
        """
        return make_ZZX(ZZX_copy(&self.x))

    def __setitem__(self, long i, a):
        """
        EXAMPLES::

            sage: n=ntl.ZZX([1,2,3])
            sage: n
            [1 2 3]
            sage: n[1] = 4
            sage: n
            [1 4 3]
        """
        if i < 0:
            raise IndexError("index (i=%s) must be >= 0" % i)
        cdef ntl_ZZ cc
        if isinstance(a, ntl_ZZ):
            cc = a
        else:
            cc = ntl_ZZ(a)
        ZZX_SetCoeff(self.x, i, cc.x)

    cdef void setitem_from_int(ntl_ZZX self, long i, int value) noexcept:
        r"""
        Set `i`-th coefficient to value.

        AUTHOR: David Harvey (2006-08-05)
        """
        ZZX_setitem_from_int(&self.x, i, value)

    def setitem_from_int_doctest(self, i, value):
        r"""
        This method exists solely for automated testing of setitem_from_int().

        sage: x = ntl.ZZX([2, 3, 4])
        sage: x.setitem_from_int_doctest(5, 42)
        sage: x
         [2 3 4 0 0 42]
        """
        self.setitem_from_int(int(i), int(value))

    def __getitem__(self, long i):
        r"""
        Retrieves coefficient number i as an NTL ZZ.

        sage: x = ntl.ZZX([129381729371289371237128318293718237, 2, -3, 0, 4])
        sage: x[0]
         129381729371289371237128318293718237
        sage: type(x[0])
         <class 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
        sage: x[1]
         2
        sage: x[2]
         -3
        sage: x[3]
         0
        sage: x[4]
         4
        sage: x[5]
         0
        """
        cdef ntl_ZZ r = ntl_ZZ()
        sig_on()
        r.x = ZZX_coeff(self.x, <long>i)
        sig_off()
        return r

    cdef int getitem_as_int(ntl_ZZX self, long i) noexcept:
        r"""
        Return `i`-th coefficient as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        return ZZX_getitem_as_int(&self.x, i)

    def getitem_as_int_doctest(self, i):
        r"""
        This method exists solely for automated testing of getitem_as_int().

        sage: x = ntl.ZZX([2, 3, 5, -7, 11])
        sage: i = x.getitem_as_int_doctest(3)
        sage: i
         -7
        sage: type(i)
         <... 'int'>
        sage: x.getitem_as_int_doctest(15)
         0
        """
        return self.getitem_as_int(i)

    def list(self):
        r"""
        Retrieves coefficients as a list of ntl.ZZ Integers.

        EXAMPLES::

            sage: x = ntl.ZZX([129381729371289371237128318293718237, 2, -3, 0, 4])
            sage: L = x.list(); L
            [129381729371289371237128318293718237, 2, -3, 0, 4]
            sage: type(L[0])
            <class 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
            sage: x = ntl.ZZX()
            sage: L = x.list(); L
            []
        """
        cdef int i
        return [self[i] for i from 0 <= i <= ZZX_deg(self.x)]

    def __add__(ntl_ZZX self, ntl_ZZX other):
        """
        EXAMPLES::

            sage: ntl.ZZX(list(range(5))) + ntl.ZZX(list(range(6)))
            [0 2 4 6 8 5]
        """
        cdef ntl_ZZX r = ntl_ZZX.__new__(ntl_ZZX)
        if not isinstance(self, ntl_ZZX):
            self = ntl_ZZX(self)
        if not isinstance(other, ntl_ZZX):
            other = ntl_ZZX(other)
        ZZX_add(r.x, (<ntl_ZZX>self).x, (<ntl_ZZX>other).x)
        return r

    def __sub__(ntl_ZZX self, ntl_ZZX other):
        """
        EXAMPLES::

            sage: ntl.ZZX(list(range(5))) - ntl.ZZX(list(range(6)))
            [0 0 0 0 0 -5]
        """
        cdef ntl_ZZX r = ntl_ZZX.__new__(ntl_ZZX)
        if not isinstance(self, ntl_ZZX):
            self = ntl_ZZX(self)
        if not isinstance(other, ntl_ZZX):
            other = ntl_ZZX(other)
        ZZX_sub(r.x, (<ntl_ZZX>self).x, (<ntl_ZZX>other).x)
        return r

    def __mul__(ntl_ZZX self, ntl_ZZX other):
        """
        EXAMPLES::

            sage: ntl.ZZX(list(range(5))) * ntl.ZZX(list(range(6)))
            [0 0 1 4 10 20 30 34 31 20]
        """
        cdef ntl_ZZX r = ntl_ZZX.__new__(ntl_ZZX)
        if not isinstance(self, ntl_ZZX):
            self = ntl_ZZX(self)
        if not isinstance(other, ntl_ZZX):
            other = ntl_ZZX(other)
        sig_on()
        ZZX_mul(r.x, (<ntl_ZZX>self).x, (<ntl_ZZX>other).x)
        sig_off()
        return r

    def __truediv__(ntl_ZZX self, ntl_ZZX other):
        """
        Compute quotient ``self / other``, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])**2
            sage: g = ntl.ZZX([4,5])
            sage: f/g
            [4 13 22 15]
            sage: ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])
            [4 13 22 15]

            sage: f = ntl.ZZX(list(range(10))); g = ntl.ZZX([-1,0,1])
            sage: f/g
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0 1 2 3 4 5 6 7 8 9]) is not divisible by other (=[-1 0 1])
        """
        sig_on()
        cdef int divisible
        cdef ZZX_c* q
        q = ZZX_div(&self.x, &other.x, &divisible)
        if not divisible:
            del q
            sig_off()
            raise ArithmeticError("self (=%s) is not divisible by other (=%s)"%(self, other))
        result = make_ZZX_sig_off(q)
        return result

    def __mod__(ntl_ZZX self, ntl_ZZX other):
        """
        Given polynomials a, b in ZZ[X], there exist polynomials q, r
        in QQ[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns q if q lies in ZZ[X], and otherwise raises an
        Exception.

        EXAMPLES::

            sage: f = ntl.ZZX([2,4,6]); g = ntl.ZZX([2])
            sage: f % g   # 0
            []

            sage: f = ntl.ZZX(list(range(10))); g = ntl.ZZX([-1,0,1])
            sage: f % g
            [20 25]
        """
        cdef ntl_ZZX r = ntl_ZZX.__new__(ntl_ZZX)
        if not isinstance(self, ntl_ZZX):
            self = ntl_ZZX(self)
        if not isinstance(other, ntl_ZZX):
            other = ntl_ZZX(other)
        sig_on()
        ZZX_rem(r.x, (<ntl_ZZX>self).x, (<ntl_ZZX>other).x)
        sig_off()
        return r

    def quo_rem(self, ntl_ZZX other):
        """
        Return the unique integral q and r such that ``self = q*other +
        r``, if they exist.  Otherwise raises an Exception.

        EXAMPLES::

           sage: f = ntl.ZZX(list(range(10))); g = ntl.ZZX([-1,0,1])
           sage: q, r = f.quo_rem(g)
           sage: q, r
           ([20 24 18 21 14 16 8 9], [20 25])
           sage: q*g + r == f
           True
        """
        cdef ZZX_c *r
        cdef ZZX_c *q
        sig_on()
        try:
            ZZX_quo_rem(&self.x, &other.x, &r, &q)
            return (make_ZZX(q), make_ZZX(r))
        finally:
            sig_off()

    def square(self):
        """
        Return f*f.

        EXAMPLES::

            sage: f = ntl.ZZX([-1,0,1])
            sage: f*f
            [1 0 -2 0 1]
        """
        sig_on()
        return make_ZZX_sig_off(ZZX_square(&self.x))

    def __pow__(ntl_ZZX self, long n, ignored):
        """
        Return the `n`-th nonnegative power of ``self``.

        EXAMPLES::

            sage: g = ntl.ZZX([-1,0,1])
            sage: g ^ 10
            [1 0 -10 0 45 0 -120 0 210 0 -252 0 210 0 -120 0 45 0 -10 0 1]
            sage: g ^ 0
            [1]
            sage: g ^ 1
            [-1 0 1]
            sage: g ^ (-1)
            Traceback (most recent call last):
            ...
            ArithmeticError
        """
        if n == 0:
            from copy import copy
            return copy(one_ZZX)
        if n < 0:
            raise ArithmeticError
        return generic_power_pos(self, <unsigned long>n)

    def __richcmp__(ntl_ZZX self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3])
            sage: g = ntl.ZZX([1,2,3,0])
            sage: f == g
            True
            sage: g = ntl.ZZX([0,1,2,3])
            sage: f == g
            False
            sage: f == ntl.ZZ(0)
            False
        """
        if op != Py_EQ and op != Py_NE:
            raise TypeError("polynomials are not ordered")

        cdef ntl_ZZX b
        try:
            b = <ntl_ZZX?>other
        except TypeError:
            b = ntl_ZZX(other)

        return (op == Py_EQ) == (self.x == b.x)

    def is_zero(self):
        """
        Return ``True`` exactly if this polynomial is 0.

        EXAMPLES::

            sage: f = ntl.ZZX([0,0,0,0])
            sage: f.is_zero()
            True
            sage: f = ntl.ZZX([0,0,1])
            sage: f
            [0 0 1]
            sage: f.is_zero()
            False
        """
        return bool(ZZX_IsZero(self.x))

    def is_one(self):
        """
        Return ``True`` exactly if this polynomial is 1.

        EXAMPLES::

            sage: f = ntl.ZZX([1,1])
            sage: f.is_one()
            False
            sage: f = ntl.ZZX([1])
            sage: f.is_one()
            True
        """
        return bool(ZZX_IsOne(self.x))

    def is_monic(self):
        """
        Return ``True`` exactly if this polynomial is monic.

        EXAMPLES::

            sage: f = ntl.ZZX([2,0,0,1])
            sage: f.is_monic()
            True
            sage: g = f.reverse()
            sage: g.is_monic()
            False
            sage: g
            [1 0 0 2]
        """
        if ZZX_IsZero(self.x):
            return False
        cdef ZZ_c lc
        lc = ZZX_LeadCoeff(self.x)
        return <bint>ZZ_IsOne(lc)

        # return bool(ZZX_is_monic(&self.x))

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: f = ntl.ZZX([2,0,0,1])
            sage: -f
            [-2 0 0 -1]
        """
        return make_ZZX(ZZX_neg(&self.x))

    def left_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the left n positions.

        EXAMPLES::

            sage: f = ntl.ZZX([2,0,0,1])
            sage: f
            [2 0 0 1]
            sage: f.left_shift(2)
            [0 0 2 0 0 1]
            sage: f.left_shift(5)
            [0 0 0 0 0 2 0 0 1]

        A negative left shift is a right shift.
            sage: f.left_shift(-2)
            [0 1]
        """
        return make_ZZX(ZZX_left_shift(&self.x, n))

    def right_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the right n positions.

        EXAMPLES::

            sage: f = ntl.ZZX([2,0,0,1])
            sage: f
            [2 0 0 1]
            sage: f.right_shift(2)
            [0 1]
            sage: f.right_shift(5)
            []
            sage: f.right_shift(-2)
            [0 0 2 0 0 1]
        """
        return make_ZZX(ZZX_right_shift(&self.x, n))

    def content(self):
        """
        Return the content of f, which has sign the same as the
        leading coefficient of f.  Also, our convention is that the
        content of 0 is 0.

        EXAMPLES::

            sage: f = ntl.ZZX([2,0,0,2])
            sage: f.content()
            2
            sage: f = ntl.ZZX([2,0,0,-2])
            sage: f.content()
            -2
            sage: f = ntl.ZZX([6,12,3,9])
            sage: f.content()
            3
            sage: f = ntl.ZZX([])
            sage: f.content()
            0
        """
        cdef ntl_ZZ r = ntl_ZZ.__new__(ntl_ZZ)
        ZZX_content(r.x, self.x)
        return r

    def primitive_part(self):
        """
        Return the primitive part of f.  Our convention is that the leading
        coefficient of the primitive part is nonnegative, and the primitive
        part of 0 is 0.

        EXAMPLES::

            sage: f = ntl.ZZX([6,12,3,9])
            sage: f.primitive_part()
            [2 4 1 3]
            sage: f
            [6 12 3 9]
            sage: f = ntl.ZZX([6,12,3,-9])
            sage: f
            [6 12 3 -9]
            sage: f.primitive_part()
            [-2 -4 -1 3]
            sage: f = ntl.ZZX()
            sage: f.primitive_part()
            []
        """
        return make_ZZX(ZZX_primitive_part(&self.x))

    def pseudo_quo_rem(self, ntl_ZZX other):
        r"""
        Perform pseudo-division: computes q and r with deg(r) <
        deg(b), and \code{LeadCoeff(b)\^(deg(a)-deg(b)+1) a = b q + r}.
        Only the classical algorithm is used.

        EXAMPLES::

            sage: f = ntl.ZZX([0,1])
            sage: g = ntl.ZZX([1,2,3])
            sage: g.pseudo_quo_rem(f)
            ([2 3], [1])
            sage: f = ntl.ZZX([1,1])
            sage: g.pseudo_quo_rem(f)
            ([-1 3], [2])
        """
        cdef ZZX_c *r
        cdef ZZX_c *q
        sig_on()
        try:
            ZZX_pseudo_quo_rem(&self.x, &other.x, &r, &q)
            return (make_ZZX(q), make_ZZX(r))
        finally:
            sig_off()

    def gcd(self, ntl_ZZX other):
        """
        Return the gcd d = gcd(a, b), where by convention the leading coefficient
        of d is >= 0.  We use a multi-modular algorithm.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])**2
            sage: g = ntl.ZZX([1,1,1])**3 * ntl.ZZX([1,2,3])
            sage: f.gcd(g)
            [1 2 3]
            sage: g.gcd(f)
            [1 2 3]
        """
        sig_on()
        return make_ZZX_sig_off(ZZX_gcd(&self.x, &other.x))

    def lcm(self, ntl_ZZX other):
        """
        Return the least common multiple of ``self`` and ``other``.

        EXAMPLES::

            sage: x1 = ntl.ZZX([-1,0,0,1])
            sage: x2 = ntl.ZZX([-1,0,0,0,0,0,1])
            sage: x1.lcm(x2)
            [-1 0 0 0 0 0 1]
        """
        g = self.gcd(other)
        return (self*other).quo_rem(g)[0]

    def xgcd(self, ntl_ZZX other, proof=None):
        r"""
        If ``self`` and ``other`` are coprime over the rationals,
        return ``r, s, t`` such that ``r = s*self + t*other``.
        Otherwise return 0.

        This is \emph{not} the same as the \sage function on
        polynomials over the integers, since here the return value r
        is always an integer.

        Here r is the resultant of a and b; if r != 0, then this
        function computes s and t such that: a*s + b*t = r; otherwise
        s and t are both 0.  If proof = False (*not* the default),
        then resultant computation may use a randomized strategy that
        errors with probability no more than `2^{-80}`.  The default is
        default is proof=None, see proof.polynomial or sage.structure.proof,
        but the global default is True), then this function may use a
        randomized strategy that errors with probability no more than
        `2^{-80}`.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])**2
            sage: g = ntl.ZZX([1,1,1])**3 * ntl.ZZX([1,2,3])
            sage: f.xgcd(g)   # nothing since they are not coprime
            (0, [], [])

        In this example the input quadratic polynomials have a common root modulo 13::

            sage: f = ntl.ZZX([5,0,1])
            sage: g = ntl.ZZX([18,0,1])
            sage: f.xgcd(g)
            (169, [-13], [13])
        """
        proof = proof_flag(proof)

        cdef ZZX_c *s
        cdef ZZX_c *t
        cdef ZZ_c *r
        sig_on()
        try:
            ZZX_xgcd(&self.x, &other.x, &r, &s, &t, proof)
            return (make_ZZ(r), make_ZZX(s), make_ZZX(t))
        finally:
            sig_off()

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES::

            sage: f = ntl.ZZX([5,0,1])
            sage: f.degree()
            2
            sage: f = ntl.ZZX(list(range(100)))
            sage: f.degree()
            99
            sage: f = ntl.ZZX()
            sage: f.degree()
            -1
            sage: f = ntl.ZZX([1])
            sage: f.degree()
            0
        """
        return ZZX_deg(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES::

            sage: f = ntl.ZZX([3,6,9])
            sage: f.leading_coefficient()
            9
            sage: f = ntl.ZZX()
            sage: f.leading_coefficient()
            0
        """
        cdef ntl_ZZ r = ntl_ZZ.__new__(ntl_ZZ)
        r.x = ZZX_LeadCoeff(self.x)
        return r

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES::

            sage: f = ntl.ZZX([3,6,9])
            sage: f.constant_term()
            3
            sage: f = ntl.ZZX()
            sage: f.constant_term()
            0
        """
        cdef ntl_ZZ r = ntl_ZZ.__new__(ntl_ZZ)
        r.x = ZZX_ConstTerm(self.x)
        return r

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES::

            sage: f = ntl.ZZX()
            sage: f.set_x()
            sage: f
            [0 1]
            sage: g = ntl.ZZX([0,1])
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory::

            sage: f is g
            False
        """
        ZZX_set_x(&self.x)

    def is_x(self):
        """
        ``True`` if this is the polynomial "x".

        EXAMPLES::

            sage: f = ntl.ZZX()
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = ntl.ZZX([0,1])
            sage: f.is_x()
            True
            sage: f = ntl.ZZX([1])
            sage: f.is_x()
            False
        """
        return bool(ZZX_is_x(&self.x))

    def derivative(self):
        """
        Return the derivative of this polynomial.

        EXAMPLES::

            sage: f = ntl.ZZX([1,7,0,13])
            sage: f.derivative()
            [7 0 39]
        """
        return make_ZZX(ZZX_derivative(&self.x))

    def reverse(self, hi=None):
        """
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.  If hi is set then this function behaves
        as if this polynomial has degree hi.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: f.reverse()
            [5 4 3 2 1]
            sage: f.reverse(hi=10)
            [0 0 0 0 0 0 5 4 3 2 1]
            sage: f.reverse(hi=2)
            [3 2 1]
            sage: f.reverse(hi=-2)
            []
        """
        if not (hi is None):
            return make_ZZX(ZZX_reverse_hi(&self.x, int(hi)))
        else:
            return make_ZZX(ZZX_reverse(&self.x))

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: f.truncate(3)
            [1 2 3]
            sage: f.truncate(8)
            [1 2 3 4 5]
            sage: f.truncate(1)
            [1]
            sage: f.truncate(0)
            []
            sage: f.truncate(-1)
            []
            sage: f.truncate(-5)
            []
        """
        if m <= 0:
            from copy import copy
            return copy(zero_ZZX)
        sig_on()
        return make_ZZX_sig_off(ZZX_truncate(&self.x, m))

    def multiply_and_truncate(self, ntl_ZZX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: g = ntl.ZZX([10])
            sage: f.multiply_and_truncate(g, 2)
            [10 20]
            sage: g.multiply_and_truncate(f, 2)
            [10 20]
        """
        if m <= 0:
            from copy import copy
            return copy(zero_ZZX)
        return make_ZZX(ZZX_multiply_and_truncate(&self.x, &other.x, m))

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: f.square_and_truncate(4)
            [1 4 10 20]
            sage: (f*f).truncate(4)
            [1 4 10 20]
        """
        if m < 0:
            from copy import copy
            return copy(zero_ZZX)
        return make_ZZX(ZZX_square_and_truncate(&self.x, m))

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of ``self`` modulo `x^m`.
        The constant term of ``self`` must be 1 or -1.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3,4,5,6,7])
            sage: f.invert_and_truncate(20)
            [1 -2 1 0 0 0 0 8 -23 22 -7 0 0 0 64 -240 337 -210 49]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -512 1344 -1176 343]
        """
        if m < 0:
            raise ArithmeticError("m (=%s) must be positive" % m)
        n = self.constant_term()
        if n != ntl_ZZ(1) and n != ntl_ZZ(-1):
            raise ArithmeticError("The constant term of self must be 1 or -1.")
        sig_on()
        return make_ZZX_sig_off(ZZX_invert_and_truncate(&self.x, m))

    def multiply_mod(self, ntl_ZZX other, ntl_ZZX modulus):
        """
        Return ``self*other % modulus``.  The modulus must be monic with
        deg(modulus) > 0, and ``self`` and ``other`` must have smaller degree.

        EXAMPLES::

            sage: modulus = ntl.ZZX([1,2,0,1])    # must be monic
            sage: g = ntl.ZZX([-1,0,1])
            sage: h = ntl.ZZX([3,7,13])
            sage: h.multiply_mod(g, modulus)
            [-10 -34 -36]
        """
        sig_on()
        return make_ZZX_sig_off(ZZX_multiply_mod(&self.x, &other.x, &modulus.x))

    def trace_mod(self, ntl_ZZX modulus):
        """
        Return the trace of this polynomial modulus the modulus.
        The modulus must be monic, and of positive degree bigger
        than the degree of ``self``.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,0,3])
            sage: mod = ntl.ZZX([5,3,-1,1,1])
            sage: f.trace_mod(mod)
            -37
        """
        sig_on()
        return make_ZZ_sig_off(ZZX_trace_mod(&self.x, &modulus.x))

    def trace_list(self):
        """
        Return the list of traces of the powers `x^i` of the
        monomial x modulo this polynomial for i = 0, ..., deg(f)-1.
        This polynomial must be monic.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,0,3,0,1])
            sage: f.trace_list()
            [5, 0, -6, 0, 10]

        The input polynomial must be monic or a :exc:`ValueError` is raised::

            sage: f = ntl.ZZX([1,2,0,3,0,2])
            sage: f.trace_list()
            Traceback (most recent call last):
            ...
            ValueError: polynomial must be monic.
        """
        if not self.is_monic():
            raise ValueError("polynomial must be monic.")
        sig_on()
        cdef char* t
        t = ZZX_trace_list(&self.x)
        r = eval(char_to_str(t).replace(' ', ','))
        string_delete(t)
        return r

    def resultant(self, ntl_ZZX other, proof=None):
        """
        Return the resultant of ``self`` and ``other``.  If proof = False (the
        default is proof=None, see proof.polynomial or sage.structure.proof,
        but the global default is True), then this function may use a
        randomized strategy that errors with probability no more than
        `2^{-80}`.

        EXAMPLES::

            sage: f = ntl.ZZX([17,0,1,1])
            sage: g = ntl.ZZX([34,-17,18,2])
            sage: f.resultant(g)
            1345873
            sage: f.resultant(g, proof=False)
            1345873
        """
        proof = proof_flag(proof)
        # NOTE: Within a factor of 2 in speed compared to MAGMA.
        sig_on()
        return make_ZZ_sig_off(ZZX_resultant(&self.x, &other.x, proof))

    def norm_mod(self, ntl_ZZX modulus, proof=None):
        """
        Return the norm of this polynomial modulo the modulus.

        The modulus must be monic, and of positive degree strictly
        greater than the degree of ``self``.  If proof=False (the default
        is proof=None, see proof.polynomial or sage.structure.proof,
        but the global default is proof=True) then it may use a
        randomized strategy that errors with probability no more than
        `2^{-80}`.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,0,3])
            sage: mod = ntl.ZZX([-5,2,0,0,1])
            sage: f.norm_mod(mod)
            -8846

        The norm is the constant term of the characteristic polynomial::

            sage: f.charpoly_mod(mod)
            [-8846 -594 -60 14 1]
        """
        proof = proof_flag(proof)
        sig_on()
        return make_ZZ_sig_off(ZZX_norm_mod(&self.x, &modulus.x, proof))

    def discriminant(self, proof=None):
        r"""
        Return the discriminant of self, which is by definition
        $$
                (-1)^{m(m-1)/2} {\mbox{\tt resultant}}(a, a')/lc(a),
        $$
        where m = deg(a), and lc(a) is the leading coefficient of a.
        If proof is False (the default is proof=None, see
        proof.polynomial or sage.structure.proof, but the global
        default is proof=True), then this function may use a
        randomized strategy that errors with probability no more than
        `2^{-80}`.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,0,3])
            sage: f.discriminant()
            -339
            sage: f.discriminant(proof=False)
            -339
        """
        proof = proof_flag(proof)
        sig_on()
        return make_ZZ_sig_off(ZZX_discriminant(&self.x, proof))

    #def __call__(self, ntl_ZZ a):
    #    sig_on()
    #    return make_ZZ_sig_off(ZZX_polyeval(&self.x, a.x))

    def charpoly_mod(self, ntl_ZZX modulus, proof=None):
        """
        Return the characteristic polynomial of this polynomial modulo
        the modulus.  The modulus must be monic of degree bigger than
        ``self``. If proof=False (the default is proof=None, see
        proof.polynomial or sage.structure.proof, but the global
        default is proof=True), then this function may use a
        randomized strategy that errors with probability no more than
        `2^{-80}`.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,0,3])
            sage: mod = ntl.ZZX([-5,2,0,0,1])
            sage: f.charpoly_mod(mod)
            [-8846 -594 -60 14 1]
        """
        proof = proof_flag(proof)
        sig_on()
        return make_ZZX_sig_off(ZZX_charpoly_mod(&self.x, &modulus.x, proof))

    def minpoly_mod_noproof(self, ntl_ZZX modulus):
        """
        Return the minimal polynomial of this polynomial modulo the
        modulus.  The modulus must be monic of degree bigger than
        ``self``.  In all cases, this function may use a randomized
        strategy that errors with probability no more than `2^{-80}`.

        EXAMPLES::

            sage: f = ntl.ZZX([0,0,1])
            sage: g = f*f
            sage: f.charpoly_mod(g)
            [0 0 0 0 1]

        However, since `f^2 = 0` modulo `g`, its minimal polynomial
        is of degree `2`::

            sage: f.minpoly_mod_noproof(g)
            [0 0 1]
        """
        sig_on()
        return make_ZZX_sig_off(ZZX_minpoly_mod(&self.x, &modulus.x))

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3])
            sage: f
            [1 2 3]
            sage: f.clear()
            sage: f
            []
        """
        ZZX_clear(&self.x)

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES::

            sage: f = ntl.ZZX([1,2,3])
            sage: f.preallocate_space(20)
            sage: f
            [1 2 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1 2 3 0 0 0 0 0 0 0 5]
        """
        sig_on()
        ZZX_preallocate_space(&self.x, n)
        sig_off()

    def squarefree_decomposition(self):
        """
        Return the square-free decomposition of ``self`` (a partial
        factorization into square-free, relatively prime polynomials)
        as a list of 2-tuples, where the first element in each tuple
        is a factor, and the second is its exponent.
        Assumes that ``self`` is primitive.

        EXAMPLES::

            sage: f = ntl.ZZX([0, 1, 2, 1])
            sage: f.squarefree_decomposition()
            [([0 1], 1), ([1 1], 2)]
        """
        cdef ZZX_c** v
        cdef long* e
        cdef long i, n
        sig_on()
        ZZX_squarefree_decomposition(&v, &e, &n, &self.x)
        sig_off()
        F = []
        for i from 0 <= i < n:
            F.append((make_ZZX(v[i]), e[i]))
        sig_free(v)
        sig_free(e)
        return F


one_ZZX = ntl_ZZX([1])
zero_ZZX = ntl_ZZX()
