"""
Univariate Polynomial Base Class

AUTHORS:

-  William Stein: first version.

-  Martin Albrecht: Added singular coercion.

-  Robert Bradshaw: Move Polynomial_generic_dense to Cython.

-  Miguel Marco: Implemented resultant in the case where PARI fails.

-  Simon King: Use a faster way of conversion from the base ring.

-  Julian Rueth (2012-05-25): Fixed is_squarefree() for imperfect fields.

TESTS::

    sage: R.<x> = ZZ[]
    sage: f = x^5 + 2*x^2 + (-1)
    sage: f == loads(dumps(f))
    True
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef is_FractionField, is_RealField, is_ComplexField
cdef coerce_binop, generic_power, parent
cdef ZZ, QQ, RR, CC, RDF, CDF

import operator, copy, re

import sage.rings.rational
import sage.rings.integer
import polynomial_ring
import sage.rings.arith
#import sage.rings.ring_element
import sage.rings.integer_ring
import sage.rings.rational_field
import sage.rings.finite_rings.integer_mod_ring
import sage.rings.complex_field
import sage.rings.fraction_field_element
import sage.rings.infinity as infinity
#import sage.misc.misc as misc
from sage.misc.sage_eval import sage_eval
from sage.misc.latex import latex
from sage.structure.factorization import Factorization
from sage.structure.element import coerce_binop

from sage.interfaces.all import singular as singular_default, is_SingularElement
from sage.libs.all import pari, pari_gen, PariError

from sage.rings.real_mpfr import RealField, is_RealField, RR

from sage.rings.complex_field import is_ComplexField, ComplexField
CC = ComplexField()

from sage.rings.real_double import is_RealDoubleField, RDF
from sage.rings.complex_double import is_ComplexDoubleField, CDF
from sage.rings.real_mpfi import is_RealIntervalField

from sage.structure.element import RingElement, generic_power, parent
from sage.structure.element cimport Element, RingElement, ModuleElement, MonoidElement

from sage.rings.rational_field import QQ, is_RationalField
from sage.rings.integer_ring import ZZ, is_IntegerRing
from sage.rings.fraction_field import is_FractionField
from sage.rings.padics.generic_nodes import is_pAdicRing, is_pAdicField

from sage.rings.integral_domain import is_IntegralDomain
from sage.structure.parent_gens cimport ParentWithGens

from sage.misc.derivative import multi_derivative

from sage.rings.arith import sort_complex_numbers_for_display

import polynomial_fateman

from sage.rings.integer cimport Integer
from sage.rings.ideal import is_Ideal
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

from sage.categories.map cimport Map
from sage.categories.morphism cimport Morphism

cpdef is_Polynomial(f):
    """
    Return True if f is of type univariate polynomial.

    INPUT:


    -  ``f`` - an object


    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_element import is_Polynomial
        sage: R.<x> = ZZ[]
        sage: is_Polynomial(x^3 + x + 1)
        True
        sage: S.<y> = R[]
        sage: f = y^3 + x*y -3*x; f
        y^3 + x*y - 3*x
        sage: is_Polynomial(f)
        True

    However this function does not return True for genuine multivariate
    polynomial type objects or symbolic polynomials, since those are
    not of the same data type as univariate polynomials::

        sage: R.<x,y> = QQ[]
        sage: f = y^3 + x*y -3*x; f
        y^3 + x*y - 3*x
        sage: is_Polynomial(f)
        False
        sage: var('x,y')
        (x, y)
        sage: f = y^3 + x*y -3*x; f
        y^3 + x*y - 3*x
        sage: is_Polynomial(f)
        False
    """
    return PY_TYPE_CHECK(f, Polynomial)

from polynomial_compiled cimport CompiledPolynomialFunction

from polydict import ETuple

cdef object is_AlgebraicRealField
cdef object is_AlgebraicField
cdef object is_AlgebraicField_common
cdef object NumberField_quadratic
cdef object is_ComplexIntervalField

cdef void late_import():
    # A hack to avoid circular imports.
    global is_AlgebraicRealField
    global is_AlgebraicField
    global is_AlgebraicField_common
    global NumberField_quadratic
    global is_ComplexIntervalField

    if is_AlgebraicRealField is not None:
        return

    import sage.rings.qqbar
    is_AlgebraicRealField = sage.rings.qqbar.is_AlgebraicRealField
    is_AlgebraicField = sage.rings.qqbar.is_AlgebraicField
    is_AlgebraicField_common = sage.rings.qqbar.is_AlgebraicField_common
    import sage.rings.number_field.number_field
    NumberField_quadratic = sage.rings.number_field.number_field.NumberField_quadratic
    import sage.rings.complex_interval_field
    is_ComplexIntervalField = sage.rings.complex_interval_field.is_ComplexIntervalField

cdef class Polynomial(CommutativeAlgebraElement):
    """
    A polynomial.

    EXAMPLE::

        sage: R.<y> = QQ['y']
        sage: S.<x> = R['x']
        sage: f = x*y; f
        y*x
        sage: type(f)
        <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        sage: p = (y+1)^10; p(1)
        1024
    """

    def __init__(self, parent, is_gen = False, construct=False):
        """
        The following examples illustrate creation of elements of
        polynomial rings, and some basic arithmetic.

        First we make a polynomial over the integers and do some
        arithmetic::

            sage: R.<x> = ZZ[]
            sage: f = x^5 + 2*x^2 + (-1); f
            x^5 + 2*x^2 - 1
            sage: f^2
            x^10 + 4*x^7 - 2*x^5 + 4*x^4 - 4*x^2 + 1

        Next we do arithmetic in a sparse polynomial ring over the
        integers::

            sage: R.<x> = ZZ[ ]; R
            Univariate Polynomial Ring in x over Integer Ring
            sage: S.<Z> = R[ ]; S
            Univariate Polynomial Ring in Z over Univariate Polynomial Ring in x over Integer Ring
            sage: f = Z^3 + (x^2-2*x+1)*Z - 3; f
            Z^3 + (x^2 - 2*x + 1)*Z - 3
            sage: f*f
            Z^6 + (2*x^2 - 4*x + 2)*Z^4 - 6*Z^3 + (x^4 - 4*x^3 + 6*x^2 - 4*x + 1)*Z^2 + (-6*x^2 + 12*x - 6)*Z + 9
            sage: f^3 == f*f*f
            True
        """
        CommutativeAlgebraElement.__init__(self, parent)
        self._is_gen = is_gen

    def _dict_to_list(self, x, zero):
          if len(x) == 0:
              return []
          n = max(x.keys())
          if PY_TYPE_CHECK(n, tuple): # a mpoly dict
              n = n[0]
              v = [zero] * (n+1)
              for i, z in x.iteritems():
                  v[i[0]] = z
          else:
              v = [zero] * (n+1)
              for i, z in x.iteritems():
                  v[i] = z
          return v

    cpdef ModuleElement _add_(self, ModuleElement right):
        cdef Py_ssize_t i, min
        x = self.list()
        y = right.list()

        if len(x) > len(y):
            min = len(y)
            high = x[min:]
        elif len(x) < len(y):
            min = len(x)
            high = y[min:]
        else:
            min = len(x)
            high = []

        low = [x[i] + y[i] for i from 0 <= i < min]
        return self._parent(low + high)

    cpdef ModuleElement _neg_(self):
        return self._parent([-x for x in self.list()])

    def plot(self, xmin=None, xmax=None, *args, **kwds):
        """
        Return a plot of this polynomial.

        INPUT:


        -  ``xmin`` - float

        -  ``xmax`` - float

        -  ``*args, **kwds`` - passed to either plot or
           point


        OUTPUT: returns a graphic object.

        EXAMPLES::

            sage: x = polygen(GF(389))
            sage: plot(x^2 + 1, rgbcolor=(0,0,1))
            sage: x = polygen(QQ)
            sage: plot(x^2 + 1, rgbcolor=(1,0,0))
        """
        R = self.base_ring()
        from sage.plot.all import plot, point, line
        if R.characteristic() == 0:
            if xmin is None and xmax is None:
                (xmin, xmax) = (-1,1)
            elif xmin is None or xmax is None:
                raise AttributeError("must give both plot endpoints")
            return plot(self.__call__, (xmin, xmax), *args, **kwds)
        else:
            if R.is_finite():
                v = list(R)
                v.sort()
                w = dict([(v[i],i) for i in range(len(v))])
                z = [(i, w[self(v[i])]) for i in range(len(v))]
                return point(z, *args, **kwds)
        raise NotImplementedError("plotting of polynomials over %s not implemented"%R)

    cpdef ModuleElement _lmul_(self, RingElement left):
        """
        Multiply self on the left by a scalar.

        EXAMPLE::

            sage: R.<x> = ZZ[]
            sage: f = (x^3 + x + 5)
            sage: f._lmul_(7)
            7*x^3 + 7*x + 35
            sage: 7*f
            7*x^3 + 7*x + 35
        """
        # todo -- should multiply individual coefficients??
        #         that could be in derived class.
        #         Note that we are guaranteed that right is in the base ring, so this could be fast.
        if left == 0:
            return self.parent().zero_element()
        return self.parent()(left) * self

    cpdef ModuleElement _rmul_(self, RingElement right):
        """
        Multiply self on the right by a scalar.

        EXAMPLE::

            sage: R.<x> = ZZ[]
            sage: f = (x^3 + x + 5)
            sage: f._rmul_(7)
            7*x^3 + 7*x + 35
            sage: f*7
            7*x^3 + 7*x + 35
        """
        # todo -- Should multiply individual coefficients??
        #         that could be in derived class.
        #         Note that we are guaranteed that right is in the base ring, so this could be fast.
        if right == 0:
            return self.parent().zero_element()
        return self * self.parent()(right)

    def subs(self, *x, **kwds):
        r"""
        Identical to self(\*x).

        See the docstring for ``self.__call__``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + x - 3
            sage: f.subs(x=5)
            127
            sage: f.subs(5)
            127
            sage: f.subs({x:2})
            7
            sage: f.subs({})
            x^3 + x - 3
            sage: f.subs({'x':2})
            Traceback (most recent call last):
            ...
            TypeError: keys do not match self's parent
        """
        if len(x) == 1 and isinstance(x[0], dict):
            g = self.parent().gen()
            if x[0].has_key(g):
                return self(x[0][g])
            elif len(x[0]) > 0:
                raise TypeError("keys do not match self's parent")
            return self
        return self.__call__(*x, **kwds)
    substitute = subs

    def __call__(self, *x, **kwds):
        """
        Evaluate polynomial at x=a.

        INPUT:


        -  ``x``:

           - ring element a; need not be in the
             coefficient ring of the polynomial.

           - a dictionary for kwds:value pairs. If the variable name
             of the polynomial is a kwds it is substituted in;
             otherwise this polynomial is returned unchanged.


        OUTPUT: the value of f at a.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x/2 - 5
            sage: f(3)
            -7/2
            sage: R.<x> = ZZ[]
            sage: f = (x-1)^5
            sage: f(2/3)
            -1/243

        We evaluate a polynomial over a quaternion algebra::

            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: R.<w> = PolynomialRing(A,sparse=True)
            sage: f = i*j*w^5 - 13*i*w^2 + (i+j)*w + i
            sage: f(i+j+1)
            24 + 26*i - 10*j - 25*k
            sage: w = i+j+1; i*j*w^5 - 13*i*w^2 + (i+j)*w + i
            24 + 26*i - 10*j - 25*k

        The parent ring of the answer always "starts" with the parent of
        the object at which we are evaluating. Thus, e.g., if we input a
        matrix, we are guaranteed to get a matrix out, though the base ring
        of that matrix may change depending on the base of the polynomial
        ring. ::

            sage: R.<x> = QQ[]
            sage: f = R(2/3)
            sage: a = matrix(ZZ,2)
            sage: b = f(a); b
            [2/3   0]
            [  0 2/3]
            sage: b.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: f = R(1)
            sage: b = f(a); b
            [1 0]
            [0 1]
            sage: b.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        ::

            sage: R.<w> = GF(17)[]
            sage: f = w^3 + 3*w +2
            sage: f(5)
            6
            sage: f(w=5)
            6
            sage: f(x=10)   # x isn't mentioned
            w^3 + 3*w + 2

        Nested polynomial ring elements can be called like multi-variate
        polynomials. ::

            sage: R.<x> = QQ[]; S.<y> = R[]
            sage: f = x+y*x+y^2
            sage: f.parent()
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: f(2)
            3*x + 4
            sage: f(2,4)
            16
            sage: f(y=2,x=4)
            16
            sage: f(2,x=4)
            16
            sage: f(2,x=4,z=5)
            16
            sage: f(2,4, z=10)
            16
            sage: f(y=x)
            2*x^2 + x
            sage: f(x=y)
            2*y^2 + y

        Polynomial ring elements can also, like multi-variate
        polynomials, be called with an argument that is a list or
        tuple containing the values to be substituted, though it is
        perhaps more natural to just unpack the list::

            sage: f([2]) # calling with a list
            3*x + 4
            sage: f((2,)) # calling with a tuple
            3*x + 4
            sage: f(*[2]) # unpacking the list to call normally
            3*x + 4

        The following results in an element of the symbolic ring. ::

            sage: f(x=sqrt(2))
            (y + sqrt(2))*y + sqrt(2)

        ::

            sage: R.<t> = PowerSeriesRing(QQ, 't'); S.<x> = R[]
            sage: f = 1 + x*t^2 + 3*x*t^4
            sage: f(2)
            1 + 2*t^2 + 6*t^4
            sage: f(2, 1/2)
            15/8

        Some special cases are optimized. ::

            sage: R.<x> = PolynomialRing(QQ, sparse=True)
            sage: f = x^3-2*x
            sage: f(x) is f
            True
            sage: f(1/x)
            (-2*x^2 + 1)/x^3

            sage: f = x^100 + 3
            sage: f(0)
            3
            sage: parent(f(0))
            Rational Field
            sage: parent(f(Qp(5)(0)))
            5-adic Field with capped relative precision 20

        TESTS:

        The following shows that :trac:`2360` is indeed fixed. ::

            sage: R.<x,y> = ZZ[]
            sage: P.<a> = ZZ[]
            sage: e = [x^2,y^3]
            sage: f = 6*a^4
            sage: f(x)
            6*x^4
            sage: f(e)
            Traceback (most recent call last):
            ...
            TypeError: Wrong number of arguments
            sage: f(x)
            6*x^4

        The following shows that :trac:`9006` is also fixed. ::

            sage: f = ZZ['x'](1000000 * [1])
            sage: f(1)
            1000000

        The following test came up in :trac:`9051`::

            sage: Cif = ComplexIntervalField(64)
            sage: R.<x> = Cif[]
            sage: f = 2*x-1
            sage: jj = Cif(RIF(0,2))
            sage: f(jj).center(), f(jj).diameter()
            (1.00000000000000000, 4.00000000000000000)

        The following failed before the patch to :trac:`3979`

        ::

            sage: R.<x> = ZZ[]
            sage: S.<y> = R[]
            sage: g = x*y + 1
            sage: g(x=3)
            3*y + 1

        AUTHORS:

        -  David Joyner (2005-04-10)

        -  William Stein (2006-01-22): change so parent is determined by the
           arithmetic

        -  William Stein (2007-03-24): fix parent being determined in the
           constant case!

        -  Robert Bradshaw (2007-04-09): add support for nested calling

        -  Tom Boothby (2007-05-01): evaluation done by
           CompiledPolynomialFunction

        -  William Stein (2007-06-03): add support for keyword arguments.

        -  Francis Clarke (2012-08-26): fix keyword substitution in the
           leading coefficient.
        """
        cdef long i, d = self.degree()

        if kwds:
            P = self.parent()
            name = P.variable_name()
            if kwds.has_key(name):
                if len(x) > 0:
                    raise ValueError("must not specify both a keyword and positional argument")
                a = self(kwds[name])
                del kwds[name]
                try:
                    return a(**kwds)
                except TypeError:
                    return a
            elif len(x) > 0:
                a = self(*x)
                try:
                    return a(**kwds)
                except TypeError:
                    return a
            else:
                try:
                    result = self[d](**kwds)
                except TypeError:
                    result = self[d]
                a = P.gen()
                i = d - 1
                while i >= 0:
                    try:
                        s = self[i](**kwds)
                    except TypeError:
                        s = self[i]
                    result = result * a + s
                    i -= 1
                return result

        if len(x) == 0:
            return self

        if isinstance(x[0], (list, tuple)):
            x = x[0]
        a = x[0]

        result = self[d]
        if len(x) > 1:
            other_args = x[1:]
            try: #if hasattr(result, '__call__'):
                result = result(other_args)
            except TypeError: #else:
                raise TypeError("Wrong number of arguments")

        if d == -1:
            try:
                return a.parent().zero_element()
            except AttributeError:
                pass
            try: # for things like the interface to the PARI C library
                return a.parent()(0)
            except AttributeError:
                return result

        if d == 0:
            try:
                return a.parent().one_element() * result
            except AttributeError:
                pass
            try: # for things like the interface to the PARI C library
                return a.parent()(1) * result
            except AttributeError:
                return result

        if parent(a) is self._parent and a.is_gen():
            return self
        elif is_FractionField(parent(a)):
            try:
                a_inverse = ~a
                if a_inverse.denominator() == 1:
                    a_inverse = a_inverse.numerator()
                    return self.reverse()(a_inverse) / a_inverse**self.degree()
            except AttributeError:
                pass

        i = d - 1
        if len(x) > 1:
            while i >= 0:
                result = result * a + self[i](other_args)
                i -= 1
        elif not a and PY_TYPE_CHECK(a, Element) and a.parent().is_exact(): #without the exactness test we can run into trouble with interval fields.
            c = self[0]
            return c + c*a
        elif (d < 4 or d > 50000) and self._compiled is None:
            while i >= 0:
                result = result * a + self[i]
                i -= 1
        else:
            if self._compiled is None:
                self._compiled = CompiledPolynomialFunction(self.list())
            result = self._compiled.eval(a)
        return result

    def _compile(self):
        # For testing
        self._compiled = CompiledPolynomialFunction(self.list())
        return self._compiled

    def _fast_float_(self, *vars):
        """
        Returns a quickly-evaluating function on floats.

        EXAMPLE::

            sage: R.<t> = QQ[]
            sage: f = t^3-t
            sage: ff = f._fast_float_()
            sage: ff(10)
            990.0

        Horner's method is used::

            sage: f = (t+10)^3; f
            t^3 + 30*t^2 + 300*t + 1000
            sage: list(f._fast_float_())
            ['load 0', 'push 30.0', 'add', 'load 0', 'mul', 'push 300.0', 'add', 'load 0', 'mul', 'push 1000.0', 'add']

        TESTS::

            sage: f = t + 2 - t
            sage: ff = f._fast_float_()
            sage: ff(3)
            2.0
            sage: list(f._fast_float_())
            ['push 2.0']

            sage: f = t - t
            sage: ff = f._fast_float_()
            sage: ff(3)
            0.0
            sage: list(f._fast_float_())
            ['push 0.0']

        """
        from sage.ext.fast_eval import fast_float_arg, fast_float_constant
        var = (<ParentWithGens>self._parent)._names[0]
        if len(vars) == 0:
            x = fast_float_arg(0)
        elif var in vars:
            x = fast_float_arg(list(vars).index(var))
        else:
            raise ValueError("free variable: %s" % var)
        cdef int i, d = self.degree()
        expr = x
        coeff = self[d]
        if d <= 0:
            return fast_float_constant(coeff)
        if coeff != 1:
            expr *= fast_float_constant(coeff)
        for i from d > i >= 0:
            coeff = self[i]
            if coeff:
                expr += fast_float_constant(coeff)
            if i > 0:
                expr *= x
        return expr

    def _fast_callable_(self, etb):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['t'])
            sage: R.<t> = QQ[]
            sage: v = R.random_element(6); v
            -t^6 - 12*t^5 + 1/2*t^4 - 1/95*t^3 - 1/2*t^2 - 4
            sage: v._fast_callable_(etb)
            add(mul(mul(add(mul(add(mul(add(mul(add(mul(v_0, -1), -12), v_0), 1/2), v_0), -1/95), v_0), -1/2), v_0), v_0), -4)

        TESTS::

            sage: R(2)._fast_callable_(etb)
            2
            sage: R(0)._fast_callable_(etb)
            0
            sage: fast_callable(R(2))(3)
            2
        """
        x = etb.var(self.variable_name())
        expr = x
        cdef int i, d = self.degree()
        coeff = self[d]
        # We handle polynomial rings like QQ['x']['y']; that gives us some
        # slowdown.  Optimize away some of that:
        if len(etb._vars) == 1:
            # OK, we're in the (very common) univariate case.
            coeff_maker = etb.constant
        else:
            # There may be variables in our coefficients...
            coeff_maker = etb.make
        if d <= 0:
            return coeff_maker(coeff)
        if coeff != 1:
            expr *= coeff_maker(coeff)
        for i from d > i >= 0:
            coeff = self[i]
            if coeff:
                expr += coeff_maker(coeff)
            if i > 0:
                expr *= x
        return expr

    cdef int _cmp_c_impl(self, Element other) except -2:
        """
        Compare the two polynomials self and other.

        We order polynomials first by degree, then in dictionary order
        starting with the coefficient of largest degree.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: 3*x^3  + 5 > 10*x^2 + 19
            True
            sage: x^2 - 2*x - 1 < x^2 - 1
            True
            sage: x^2 - 2*x - 1 > x^2 - 1
            False
            sage: R(-1) < 0
            False
            sage: x^3 - 3 > 393939393
            True
        """
        d1 = self.degree(); d2 = other.degree()
        c = cmp(d1, d2)
        if c: return c
        for i in reversed(xrange(d1+1)):
            c = cmp(self[i], other[i])
            if c: return c
        return 0

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: P = PolynomialRing(ZZ,'x')(0)
            sage: bool(P)
            False
            sage: P = PolynomialRing(ZZ, 'x')([1,2,3])
            sage: bool(P)
            True
        """
        return self.degree() >= 0

    def __getitem__(self, n):
        raise NotImplementedError

    def __iter__(self):
        return iter(self.list())

    # you may have to replicate this boilerplate code in derived classes if you override
    # __richcmp__.  The python documentation at  http://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    cdef long _hash_c(self):
        """
        This hash incorporates the variable name in an effort to respect
        the obvious inclusions into multi-variable polynomial rings.

        The tuple algorithm is borrowed from
        http://effbot.org/zone/python-hash.htm.

        EXAMPLES::

            sage: R.<x>=ZZ[]
            sage: hash(R(1))==hash(1)  # respect inclusions of the integers
            True
            sage: hash(R.0)==hash(FractionField(R).0)  # respect inclusions into the fraction field
            True
            sage: R.<x>=QQ[]
            sage: hash(R(1/2))==hash(1/2)  # respect inclusions of the rationals
            True
            sage: hash(R.0)==hash(FractionField(R).0)  # respect inclusions into the fraction field
            True
            sage: R.<x>=IntegerModRing(11)[]
            sage: hash(R.0)==hash(FractionField(R).0)  # respect inclusions into the fraction field
            True
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash
        cdef int i
        for i from 0<= i <= self.degree():
            if i == 1:
                # we delay the hashing until now to not waste it one a constant poly
                var_name_hash = hash((<ParentWithGens>self._parent)._names[0])
            #  I'm assuming (incorrectly) that hashes of zero indicate that the element is 0.
            # This assumption is not true, but I think it is true enough for the purposes and it
            # it allows us to write fast code that omits terms with 0 coefficients.  This is
            # important if we want to maintain the '==' relationship with sparse polys.
            c_hash = hash(self[i])
            if c_hash != 0:
                if i == 0:
                    result += c_hash
                else:
                    # Hash (self[i], generator, i) as a tuple according to the algorithm.
                    result_mon = c_hash
                    result_mon = (1000003 * result_mon) ^ var_name_hash
                    result_mon = (1000003 * result_mon) ^ i
                    result += result_mon
        if result == -1:
            return -2
        return result

    def __float__(self):
         if self.degree() > 0:
             raise TypeError("cannot coerce nonconstant polynomial to float")
         return float(self[0])

    def __int__(self):
        if self.degree() > 0:
            raise TypeError("cannot coerce nonconstant polynomial to int")
        return int(self[0])

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: H = Hom(R, QQ); H
            Set of Homomorphisms from Univariate Polynomial Ring in x over Integer Ring to Rational Field
            sage: f = H([5]); f
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn: x |--> 5
            sage: f(x)
            5
            sage: f(x^2 + 3)
            28
        """
        a = im_gens[0]
        P = a.parent()
        d = self.degree()
        result = P._coerce_(self[d])
        i = d - 1
        while i >= 0:
            result = result * a + P._coerce_(self[i])
            i -= 1
        return result

    def _integer_(self, ZZ):
        r"""
        EXAMPLES::

            sage: k = GF(47)
            sage: R.<x> = PolynomialRing(k)
            sage: ZZ(R(45))
            45
            sage: ZZ(3*x + 45)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial
        """
        if self.degree() > 0:
            raise TypeError("cannot coerce nonconstant polynomial")
        return ZZ(self[0])

    def _rational_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: QQ(R(45/4))
            45/4
            sage: QQ(3*x + 45)
            Traceback (most recent call last):
            ...
            TypeError: not a constant polynomial
        """
        if self.degree() > 0:
            raise TypeError("not a constant polynomial")
        return sage.rings.rational.Rational(self[0])

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + x
            sage: g = f._symbolic_(SR); g
            (x^2 + 1)*x
            sage: g(x=2)
            10

            sage: g = SR(f)
            sage: g(x=2)
            10

        The polynomial does not have to be over a field of
        characteristic 0::

            sage: R.<w> = GF(7)[]
            sage: f = SR(2*w^3 + 1); f
            2*w^3 + 1
            sage: f.variables()
            (w,)
        """
        d = dict([(repr(g), R.var(g)) for g in self.parent().gens()])
        return self.subs(**d)

    def __invert__(self):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x - 90283
            sage: f.__invert__()
            1/(x - 90283)
            sage: ~f
            1/(x - 90283)
        """
        return self.parent().one_element()/self

    def inverse_of_unit(self):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x - 90283
            sage: f.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ValueError: self is not a unit.
            sage: f = R(-90283); g = f.inverse_of_unit(); g
            -1/90283
            sage: parent(g)
            Univariate Polynomial Ring in x over Rational Field
        """
        if self.degree() > 0:
            if not self.is_unit():
                raise ValueError("self is not a unit.")
            else:
                raise NotImplementedError("polynomial inversion over non-integral domains not implemented")
        return self.parent()(~(self[0]))

    def inverse_mod(a, m):
        """
        Inverts the polynomial a with respect to m, or raises a ValueError
        if no such inverse exists. The parameter m may be either a single
        polynomial or an ideal (for consistency with inverse_mod in other
        rings).

        EXAMPLES::

            sage: S.<t> = QQ[]
            sage: f = inverse_mod(t^2 + 1, t^3 + 1); f
            -1/2*t^2 - 1/2*t + 1/2
            sage: f * (t^2 + 1) % (t^3 + 1)
            1
            sage: f = t.inverse_mod((t+1)^7); f
            -t^6 - 7*t^5 - 21*t^4 - 35*t^3 - 35*t^2 - 21*t - 7
            sage: (f * t) + (t+1)^7
            1
            sage: t.inverse_mod(S.ideal((t + 1)^7)) == f
            True

        This also works over inexact rings, but note that due to rounding
        error the product may not always exactly equal the constant
        polynomial 1 and have extra terms with coefficients close to zero. ::

            sage: R.<x> = RDF[]
            sage: epsilon = RDF(1).ulp()*50   # Allow an error of up to 50 ulp
            sage: f = inverse_mod(x^2 + 1, x^5 + x + 1); f
            0.4*x^4 - 0.2*x^3 - 0.4*x^2 + 0.2*x + 0.8
            sage: poly = f * (x^2 + 1) % (x^5 + x + 1)
            sage: # Remove noisy zero terms:
            sage: parent(poly)([ 0.0 if abs(c)<=epsilon else c for c in poly.coeffs() ])
            1.0
            sage: f = inverse_mod(x^3 - x + 1, x - 2); f
            0.142857142857
            sage: f * (x^3 - x + 1) % (x - 2)
            1.0
            sage: g = 5*x^3+x-7; m = x^4-12*x+13; f = inverse_mod(g, m); f
            -0.0319636125...*x^3 - 0.0383269759...*x^2 - 0.0463050900...*x + 0.346479687...
            sage: poly = f*g % m
            sage: # Remove noisy zero terms:
            sage: parent(poly)([ 0.0 if abs(c)<=epsilon else c for c in poly.coeffs() ])
            1.0

        ALGORITHM: Solve the system as + mt = 1, returning s as the inverse
        of a mod m.

        Uses the Euclidean algorithm for exact rings, and solves a linear
        system for the coefficients of s and t for inexact rings (as the
        Euclidean algorithm may not converge in that case).

        AUTHORS:

        - Robert Bradshaw (2007-05-31)
        """
        from sage.rings.ideal import is_Ideal
        if is_Ideal(m):
            v = m.gens_reduced()
            if len(v) > 1:
                raise NotImplementedError("Don't know how to invert modulo non-principal ideal %s" % m)
            m = v[0]
        if m.degree() == 1 and m[1].is_unit():
            # a(x) mod (x-r) = a(r)
            r = -m[0]
            if not m[1].is_one():
                r *= m.base_ring()(~m[1])
            u = a(r)
            if u.is_unit():
                return a.parent()(~u)
        if a.parent().is_exact():
            # use xgcd
            g, s, _ = a.xgcd(m)
            if g == 1:
                return s
            elif g.is_unit():
                return g.inverse_of_unit() * s
            else:
                raise ValueError("Impossible inverse modulo")
        else:
            # xgcd may not converge for inexact rings.
            # Instead solve for the coefficients of
            # s (degree n-1) and t (degree n-2) in
            #               as + mt = 1
            # as a linear system.
            from sage.matrix.constructor import matrix
            from sage.modules.free_module_element import vector
            a %= m
            n = m.degree()
            R = a.parent().base_ring()
            M = matrix(R, 2*n-1)
            # a_i s_j x^{i+j} terms
            for i in range(n):
                for j in range(n):
                    M[i+j, j] = a[i]
            # m_i t_j x^{i+j} terms
            for i in range(n+1):
                for j in range(n-1):
                    M[i+j, j+n] = m[i]
            v = vector(R, [R.one_element()] + [R.zero_element()]*(2*n-2)) # the constant polynomial 1
            if M.is_invertible():
                x = M.solve_right(v) # there has to be a better way to solve
                return a.parent()(list(x)[0:n])
            else:
                raise ValueError("Impossible inverse modulo")

    def __long__(self):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x - 902384
            sage: long(f)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
            sage: long(R(939392920202))
            939392920202L
        """
        if self.degree() > 0:
            raise TypeError("cannot coerce nonconstant polynomial to long")
        return long(self[0])

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x - 4)*(x^2 - 8*x + 16)
            x^3 - 12*x^2 + 48*x - 64
            sage: C.<t> = PowerSeriesRing(ZZ)
            sage: D.<s> = PolynomialRing(C)
            sage: z = (1 + O(t)) + t*s^2
            sage: z*z
            t^2*s^4 + (2*t + O(t^2))*s^2 + 1 + O(t)

            ## More examples from trac 2943, added by Kiran S. Kedlaya 2 Dec 09
            sage: C.<t> = PowerSeriesRing(Integers())
            sage: D.<s> = PolynomialRing(C)
            sage: z = 1 + (t + O(t^2))*s + (t^2 + O(t^3))*s^2
            sage: z*z
            (t^4 + O(t^5))*s^4 + (2*t^3 + O(t^4))*s^3 + (3*t^2 + O(t^3))*s^2 + (2*t + O(t^2))*s + 1
        """
        if right == 0 or self == 0:
            return self._parent.zero_element()

        if self._parent.is_exact():
            return self._mul_karatsuba(right)
        else:
            return self._mul_generic(right)

    def square(self):
        """
        Returns the square of this polynomial.

        TODO:

        - This is just a placeholder; for now it just uses ordinary
          multiplication. But generally speaking, squaring is faster than
          ordinary multiplication, and it's frequently used, so subclasses
          may choose to provide a specialised squaring routine.

        - Perhaps this even belongs at a lower level? ring_element or
          something?

        AUTHORS:

        - David Harvey (2006-09-09)

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + 1
            sage: f.square()
            x^6 + 2*x^3 + 1
            sage: f*f
            x^6 + 2*x^3 + 1
        """
        return self * self

    def squarefree_decomposition(self):
        """
        Return the square-free decomposition of self.  This is a
        partial factorization of self into square-free, coprime
        polynomials.

        ALGORITHM: In characteristic 0, we use Yun's algorithm,
        which works for arbitrary rings of characteristic 0.
        If the characteristic is a prime number `p > 0`, we use
        [Coh]_, algorithm 3.4.2.  This is basically Yun's algorithm
        with special treatment for powers divisible by `p`.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: p = 37 * (x-1)^3 * (x-2)^3 * (x-1/3)^7 * (x-3/7)
            sage: p.squarefree_decomposition()
            (37*x - 111/7) * (x^2 - 3*x + 2)^3 * (x - 1/3)^7
            sage: p = 37 * (x-2/3)^2
            sage: p.squarefree_decomposition()
            (37) * (x - 2/3)^2
            sage: x = polygen(GF(3))
            sage: x.squarefree_decomposition()
            x
            sage: f = QQbar['x'](1)
            sage: f.squarefree_decomposition()
            1

        TESTS::

            sage: K.<a> = GF(3^2)
            sage: R.<x> = K[]
            sage: f = x^243+2*x^81+x^9+1
            sage: f.squarefree_decomposition()
            (x^27 + 2*x^9 + x + 1)^9
            sage: f = x^243+a*x^27+1
            sage: f.squarefree_decomposition()
            (x^9 + (2*a + 1)*x + 1)^27

            sage: for K in [GF(2^18,'a'), GF(3^2,'a'), GF(47^3,'a')]:
            ....:     R.<x> = K[]
            ....:     if K.characteristic() < 5: m = 4
            ....:     else: m = 1
            ....:     for _ in range(m):
            ....:         f = (R.random_element(4)^3*R.random_element(m)^(m+1))(x^6)
            ....:         F = f.squarefree_decomposition()
            ....:         assert F.prod() == f
            ....:         for i in range(len(F)):
            ....:             assert gcd(F[i][0], F[i][0].derivative()) == 1
            ....:             for j in range(len(F)):
            ....:                 if i == j: continue
            ....:                 assert gcd(F[i][0], F[j][0]) == 1
            ....:

        REFERENCES:

        .. [Coh] H. Cohen, A Course in Computational Algebraic Number
           Theory.  Springer-Verlag, 1993.

        """
        if not self.base_ring().is_unique_factorization_domain():
            raise NotImplementedError, "Squarefree decomposition not implemented for " + str(self.parent())

        if self.degree() == 0:
            return Factorization([], unit=self[0])

        p = self.base_ring().characteristic()
        factors = []
        if p == 0:
            f = [self]
            cur = self
            while cur.degree() > 0:
                cur = cur.gcd(cur.derivative())
                f.append(cur)

            g = []
            for i in range(len(f) - 1):
                g.append(f[i] // f[i+1])

            a = []
            for i in range(len(g) - 1):
                a.append(g[i] // g[i+1])
            a.append(g[-1])

            unit = f[-1]
            for i in range(len(a)):
                if a[i].degree() > 0:
                    factors.append((a[i], i+1))
                else:
                    unit = unit * a[i].constant_coefficient() ** (i + 1)
        else:
            # Beware that `p`-th roots might not exist.
            unit = self.leading_coefficient()
            T0 = self.monic()
            e = 1
            if T0.degree() > 0:
                der = T0.derivative()
                while der.is_zero():
                    T0 = T0.parent()([T0[p*i].pth_root() for i in range(T0.degree()//p + 1)])
                    if T0 == 1:
                        raise RuntimeError
                    der = T0.derivative()
                    e = e*p
                T = T0.gcd(der)
                V = T0 // T
                k = 0
                while T0.degree() > 0:
                    k += 1
                    if p.divides(k):
                        T = T // V
                        k += 1
                    W = V.gcd(T)
                    if W.degree() < V.degree():
                        factors.append((V // W, e*k))
                        V = W
                        T = T // V
                        if V.degree() == 0:
                            if T.degree() == 0:
                                break
                            # T is of the form sum_{i=0}^n t_i X^{pi}
                            T0 = T0.parent()([T[p*i].pth_root() for i in range(T.degree()//p + 1)])
                            der = T0.derivative()
                            e = p*e
                            while der.is_zero():
                                T0 = T0.parent()([T0[p*i].pth_root() for i in range(T0.degree()//p + 1)])
                                der = T0.derivative()
                                e = p*e
                            T = T0.gcd(der)
                            V = T0 // T
                            k = 0
                    else:
                        T = T//V

        return Factorization(factors, unit=unit, sort=False)

    def is_square(self, root=False):
        """
        Returns whether or not polynomial is square. If the optional
        argument root is True, then also returns the square root (or None,
        if the polynomial is not square).

        INPUT:


        -  ``root`` - whether or not to also return a square
           root (default: False)


        OUTPUT:


        -  ``bool`` - whether or not a square

        -  ``object`` - (optional) an actual square root if
           found, and None otherwise.


        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 12*(x+1)^2 * (x+3)^2
            sage: S.<y> = PolynomialRing(RR)
            sage: g = 12*(y+1)^2 * (y+3)^2
            sage: f.is_square()
            False
            sage: f.is_square(root=True)
            (False, None)
            sage: g.is_square()
            True
            sage: h = f/3; h
            4*x^4 + 32*x^3 + 88*x^2 + 96*x + 36
            sage: h.is_square(root=True)
            (True, 2*x^2 + 8*x + 6)

        TESTS:

        Make sure ticket #9093 is fixed::

            sage: R(1).is_square()
            True
            sage: R(4/9).is_square()
            True
            sage: R(-1/3).is_square()
            False
            sage: R(0).is_square()
            True

        """
        if self.degree() < 0:
            if root:
                return True, self
            else:
                return True
        from sage.rings.arith import gcd
        R = self.base_ring()
        P = self.parent()
        try:
            f = self.squarefree_decomposition()
        except NotImplementedError:
            f = self.factor()
        u = R(f.unit())
        # are all powers in the factorization even?
        if (gcd([a[1] for a in f])%2 == 0
            # is the unit a square?
            and u.is_square()):
            # build square root
            g = P(u.sqrt())
            for a in f:
                g *= a[0]**(a[1]/2)
            if root:
                return True, g
            else:
                return True
        else:
            if root:
                return False, None
            else:
                return False

    def any_root(self, ring=None, degree=None, assume_squarefree=False):
        """
        Return a root of this polynomial in the given ring.

        INPUT:

        - ``ring`` -- The ring in which a root is sought.  By default
          this is the coefficient ring.

        - ``degree`` (None or nonzero integer) -- Used for polynomials
          over finite fields.  Returns a root of degree
          ``abs(degree)`` over the ground field.  If negative, also
          assumes that all factors of this polynomial are of degree
          ``abs(degree)``.  If None, returns a root of minimal degree
          contained within the given ring.

        - ``assume_squarefree`` (bool) -- Used for polynomials over
          finite fields.  If True, this polynomial is assumed to be
          squarefree.

        EXAMPLES::

            sage: R.<x> = GF(11)[]
            sage: f = 7*x^7 + 8*x^6 + 4*x^5 + x^4 + 6*x^3 + 10*x^2 + 8*x + 5
            sage: f.any_root()
            2
            sage: f.factor()
            (7) * (x + 9) * (x^6 + 10*x^4 + 6*x^3 + 5*x^2 + 2*x + 2)
            sage: f = x^6 + 10*x^4 + 6*x^3 + 5*x^2 + 2*x + 2
            sage: f.any_root(GF(11^6, 'a'))
            a^5 + a^4 + 7*a^3 + 2*a^2 + 10*a
            sage: sorted(f.roots(GF(11^6, 'a')))
            [(10*a^5 + 2*a^4 + 8*a^3 + 9*a^2 + a, 1), (10*a^5 + 3*a^4 + 8*a^3 + a^2 + 10*a + 4, 1), (2*a^5 + 8*a^4 + 3*a^3 + 6*a + 2, 1), (9*a^5 + 5*a^4 + 10*a^3 + 8*a^2 + 3*a + 1, 1), (a^5 + 3*a^4 + 8*a^3 + 2*a^2 + 3*a + 4, 1), (a^5 + a^4 + 7*a^3 + 2*a^2 + 10*a, 1)]
            sage: f.any_root(GF(11^6, 'a'))
            a^5 + a^4 + 7*a^3 + 2*a^2 + 10*a

            sage: g = (x-1)*(x^2 + 3*x + 9) * (x^5 + 5*x^4 + 8*x^3 + 5*x^2 + 3*x + 5)
            sage: g.any_root(ring=GF(11^10, 'b'), degree=1)
            1
            sage: g.any_root(ring=GF(11^10, 'b'), degree=2)
            6*b^9 + 7*b^7 + 7*b^6 + 3*b^5 + b^2 + b + 3
            sage: g.any_root(ring=GF(11^10, 'b'), degree=5)
            7*b^9 + 2*b^8 + 6*b^7 + 3*b^6 + 2*b^5 + 7*b^3 + 7*b^2 + 2

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: K.<a> = GF(5^12)
            sage: for _ in range(40):
            ....:     f = R.random_element(degree=4)
            ....:     assert f(f.any_root(K)) == 0

        """
        if self.base_ring().is_finite() and self.base_ring().is_field():
            if self.degree() < 0:
                return ring(0)
            if self.degree() == 0:
                raise ValueError, "no roots A %s"%self
            if not assume_squarefree:
                SFD = self.squarefree_decomposition()
                SFD.sort()
                for f, e in SFD:
                    try:
                        return f.any_root(ring, degree, True)
                    except ValueError:
                        pass
            if self.degree() == 1 and (degree is None or degree == 1):
                if ring is None:
                    return -self[0]/self[1]
                else:
                    return ring(-self[0]/self[1])
            q = self.base_ring().order()
            if ring is None:
                allowed_deg_mult = Integer(1)
            else:
                if not self.base_ring().is_prime_field():
                    raise NotImplementedError
                if ring.characteristic() != self.base_ring().characteristic():
                    raise ValueError, "ring must be an extension of the base ring"
                if not (ring.is_field() and ring.is_finite()):
                    raise NotImplementedError
                allowed_deg_mult = Integer(ring.factored_order()[0][1]) # generally it will be the quotient of this by the degree of the base ring.
            if degree is None:
                x = self.parent().gen()
                if allowed_deg_mult == 1:
                    self = self.gcd(x**q-x)
                    degree = -1
                    if self.degree() == 0:
                        raise ValueError, "no roots B %s"%self
                else:
                    xq = x
                    d = Integer(0)
                    while True:
                        # one pass for roots that actually lie within ring.
                        e = self.degree()
                        if 2*d+2 > e:
                            # this polynomial has no factors dividing allowed_deg_mult
                            if allowed_deg_mult % e == 0:
                                degree = -e
                            break
                        while d < allowed_deg_mult:
                            d = d+1
                            xq = xq**q % self
                            if d.divides(allowed_deg_mult):
                                break
                        A = self.gcd(xq-x)
                        if A != 1:
                            self = A
                            degree = -d
                            break
                        if d == allowed_deg_mult:
                            break
                    if degree is None:
                        if allowed_deg_mult == 1:
                            raise ValueError, "no roots C %s"%self
                        xq = x
                        d = Integer(0)
                        while True:
                            # now another for roots that will lie in an extension.
                            e = self.degree()
                            if 2*d+2 > e:
                                # this polynomial is irreducible.
                                degree = -e
                                break
                            while True:
                                # we waste a little effort here in computing the xq again.
                                d = d+1
                                xq = xq**q % self
                                if allowed_deg_mult.divides(d):
                                    break
                            A = self.gcd(xq-x)
                            if A != 1:
                                self = A
                                degree = -d
                                break
            if degree == 0:
                raise ValueError, "degree should be nonzero"
            R = self.parent()
            x = R.gen()
            if degree > 0:
                xq = x
                d = 0
                while True:
                    e = self.degree()
                    if 2*d > e:
                        if degree != e:
                            raise ValueError, "no roots D %s"%self
                        break
                    d = d+1
                    xq = xq**q % self
                    if d == degree:
                        break
                    A = self.gcd(xq-x)
                    if A != 1:
                        self = self // A
                if d == degree:
                    self = self.gcd(xq-x)
                    if self.degree() == 0:
                        raise ValueError, "no roots E %s"%self
            else:
                degree = -degree
            if ring is None:
                if degree == 1:
                    ring = self.base_ring()
                else:
                    ring = self.base_ring().extension(degree) # this won't work yet.
            # now self has only roots of degree ``degree``.
            # for now, we only implement the Cantor-Zassenhaus split
            k = self.degree() // degree
            if k == 1:
                try:
                    return self.roots(ring, multiplicities=False)[0] # is there something better to do here?
                except IndexError:
                    raise ValueError, "no roots F %s"%self
            if q % 2 == 0:
                T = x
                while True:
                    C = T
                    for i in range(degree-1):
                        C = T + C**q % self
                    h = self.gcd(C)
                    hd = h.degree()
                    if hd == 0 or hd == self.degree():
                        T = T * x**q
                    else:
                        if 2*hd <= self.degree():
                            return h.any_root(ring, -degree, True)
                        else:
                            return (self//h).any_root(ring, -degree, True)
            else:
                from sage.rings.arith import power_mod
                while True:
                    T = R.random_element(2*degree-1)
                    if T == 0:
                        continue
                    T = T.monic()
                    h = self.gcd(power_mod(T, Integer((q**degree-1)/2), self)-1)
                    hd = h.degree()
                    if hd != 0 and hd != self.degree():
                        if 2*hd <= self.degree():
                            return h.any_root(ring, -degree, True)
                        else:
                            return (self//h).any_root(ring, -degree, True)
        else:
            return self.roots(ring=ring, multiplicities=False)[0]


    def __div__(self, right):
        """
        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = (x^3 + 5)/3; f
            1/3*x^3 + 5/3
            sage: f.parent()
            Univariate Polynomial Ring in x over Rational Field

        If we do the same over `\ZZ` the result is in the
        polynomial ring over `\QQ`.

        ::

            sage: x  = ZZ['x'].0
            sage: f = (x^3 + 5)/3; f
            1/3*x^3 + 5/3
            sage: f.parent()
            Univariate Polynomial Ring in x over Rational Field

        Divides can make elements of the fraction field::

            sage: R.<x> = QQ['x']
            sage: f = x^3 + 5
            sage: g = R(3)
            sage: h = f/g; h
            1/3*x^3 + 5/3
            sage: h.parent()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field

        This is another example over a non-prime finite field (submitted by
        a student of Jon Hanke). It illustrates cancellation between the
        numerator and denominator over a non-prime finite field.

        ::

            sage: R.<x> = PolynomialRing(GF(5^2, 'a'), 'x')
            sage: f = x^3 + 4*x
            sage: f / (x - 1)
            x^2 + x

        Be careful about coercions (this used to be broken)::

            sage: R.<x> = ZZ['x']
            sage: f = x / Mod(2,5); f
            3*x
            sage: f.parent()
            Univariate Polynomial Ring in x over Ring of integers modulo 5

        TESTS:

        Check that :trac:`12217` is fixed::

            sage: P.<x> = GF(5)[]
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

            sage: P.<x> = GF(25, 'a')[]
            sage: x/5
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in Finite Field in a of size 5^2

        """
        try:
            if not isinstance(right, Element) or right.parent() != self.parent():
                R = self.parent().base_ring()
                x = R._coerce_(right)
                return self * ~x
        except (TypeError, ValueError):
            pass
        return RingElement.__div__(self, right)


    def __pow__(self, right, modulus):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x - 1
            sage: f._pow(3)
            x^3 - 3*x^2 + 3*x - 1
            sage: f^3
            x^3 - 3*x^2 + 3*x - 1

            sage: R = PolynomialRing(GF(2), x)
            sage: f = R(x^9 + x^7 + x^6 + x^5 + x^4 + x^2 + x)
            sage: h = R(x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 1)
            sage: pow(f, 2, h)
            x^9 + x^8 + x^7 + x^5 + x^3

        TESTS::

            sage: x^(1/2)
            Traceback (most recent call last):
            ...
            TypeError: non-integral exponents not supported

        ::

            sage: x^x
            Traceback (most recent call last):
            ...
            TypeError: non-integral exponents not supported
        """
        if not PY_TYPE_CHECK_EXACT(right, Integer) or \
                PY_TYPE_CHECK_EXACT(right, int):
                    try:
                        right = Integer(right)
                    except TypeError:
                        raise TypeError("non-integral exponents not supported")

        if self.degree() <= 0:
            r = self.parent()(self[0]**right)
        elif right < 0:
            r = (~self)**(-right)
        elif (<Polynomial>self) == self.parent().gen():   # special case x**n should be faster!
            P = self.parent()
            R = P.base_ring()
            if P.is_sparse():
                v = {right:R.one_element()}
            else:
                v = [R.zero_element()]*right + [R.one_element()]
            r = self.parent()(v, check=False)
        else:
            r = generic_power(self, right)
        if modulus:
            return r % modulus
        else:
            return r

    def _pow(self, right):
        # TODO: fit __pow__ into the arithmetic structure
        if self.degree() <= 0:
            return self.parent()(self[0]**right)
        if right < 0:
            return (~self)**(-right)
        if (<Polynomial>self) == self.parent().gen():   # special case x**n should be faster!
            v = [0]*right + [1]
            return self.parent()(v, check=True)
        return generic_power(self, right)

    def _repr(self, name=None):
        """
        Return the string representation of this polynomial.

        INPUT:

        - ``name`` - None or a string; used for printing the variable.

        EXAMPLES::

            sage: S.<t> = QQ[]
            sage: R.<x> = S[]
            sage: f = (1 - t^3)*x^3 - t^2*x^2 - x + 1
            sage: f._repr()
            '(-t^3 + 1)*x^3 - t^2*x^2 - x + 1'
            sage: f._repr('z')
            '(-t^3 + 1)*z^3 - t^2*z^2 - z + 1'
            sage: P.<y> = RR[]
            sage: y, -y
            (y, -y)
        """
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = self.list()
        for n in reversed(xrange(m)):
            x = coeffs[n]
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*',' ', s)
        s = re.sub(r' -1(\.0+)?\*',' -', s)
        if s == " ":
            return "0"
        return s[1:]

    def _repr_(self):
        r"""
        Return string representation of this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ, implementation="FLINT")
            sage: f = x^3+2/3*x^2 - 5/3
            sage: f._repr_()
            'x^3 + 2/3*x^2 - 5/3'
            sage: f.rename('vaughn')
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support renaming: x^3 + 2/3*x^2 - 5/3
        """
        return self._repr()

    def _latex_(self, name=None):
        r"""
        Return the latex representation of this polynomial.

        EXAMPLES:

        A fairly simple example over `\QQ`.

        ::

            sage: C3.<omega> = CyclotomicField(3)
            sage: R.<X> = C3[]
            sage: f = X^3 - omega*X
            sage: latex(f)
            X^{3} - \omega X
            sage: R.<x> = RDF[]
            sage: latex(x+2)
            x + 2.0

        The following illustrates the fix of trac #2586::

            sage: latex(ZZ['alpha']['b']([0, ZZ['alpha'].0]))
            \alpha b

        The following illustrates a (non-intentional) superfluity of parentheses

            sage: K.<I>=QuadraticField(-1)
            sage: R.<x>=K[]
            sage: latex(I*x^2-I*x)
            \left(\sqrt{-1}\right) x^{2} + \left(-\sqrt{-1}\right) x
        """
        s = " "
        coeffs = self.list()
        m = len(coeffs)
        if name is None:
            name = self.parent().latex_variable_names()[0]
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        for n in reversed(xrange(m)):
            x = coeffs[n]
            x = y = latex(x)
            if x != '0':
                if n != m-1:
                    s += " + "
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "|%s^{%s}"%(name,n)
                elif n==1:
                    var = "|%s"%name
                else:
                    var = ""
                s += "%s %s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(" 1(\.0+)? \|"," ", s)
        s = re.sub(" -1(\.0+)? \|", " -", s)
        s = s.replace("|","")
        if s==" ":
            return "0"
        return s[1:].lstrip().rstrip()

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: K.<x> = ZZ[]
            sage: sage_input(K(0), verify=True)
            # Verified
            ZZ['x'](0)
            sage: sage_input(K(-54321), preparse=False, verify=True)
            # Verified
            ZZ['x'](-54321)
            sage: sage_input(x, verify=True)
            # Verified
            R.<x> = ZZ[]
            x
            sage: sage_input(x, preparse=False)
            R = ZZ['x']
            x = R.gen()
            x
            sage: sage_input((3*x-2)^3, verify=True)
            # Verified
            R.<x> = ZZ[]
            27*x^3 - 54*x^2 + 36*x - 8
            sage: L.<y> = K[]
            sage: sage_input(L(0), verify=True)
            # Verified
            ZZ['x']['y'](0)
            sage: sage_input((x+y+1)^2, verify=True)
            # Verified
            R1.<x> = ZZ[]
            R2.<y> = R1[]
            y^2 + (2*x + 2)*y + (x^2 + 2*x + 1)
            sage: sage_input(RR(pi) * polygen(RR), verify=True)
            # Verified
            R.<x> = RR[]
            3.1415926535897931*x
            sage: sage_input(polygen(GF(7)) + 12, verify=True)
            # Verified
            R.<x> = GF(7)[]
            x + 5
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: K(0)._sage_input_(SageInputBuilder(), True)
            {atomic:0}
            sage: (x^2 - 1)._sage_input_(SageInputBuilder(), False)
            {binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:1}}
        """
        if self.degree() > 0:
            gen = sib.gen(self.parent())
            coeffs = self.list()
            terms = []
            for i in range(len(coeffs)-1, -1, -1):
                if i > 0:
                    if i > 1:
                        gen_pow = gen**sib.int(i)
                    else:
                        gen_pow = gen
                    terms.append(sib.prod((sib(coeffs[i], True), gen_pow), simplify=True))
                else:
                    terms.append(sib(coeffs[i], True))
            return sib.sum(terms, simplify=True)
        elif coerced:
            return sib(self.constant_coefficient(), True)
        else:
            return sib(self.parent())(sib(self.constant_coefficient(), True))

    def __setitem__(self, n, value):
        """
        Set the n-th coefficient of this polynomial. This always raises an
        IndexError, since in Sage polynomials are immutable.

        INPUT:


        -  ``n`` - an integer

        -  ``value`` - value to set the n-th coefficient to


        OUTPUT: an IndexError is always raised.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^3 + x + 1
            sage: f[2] = 3
            Traceback (most recent call last):
            ...
            IndexError: polynomials are immutable
        """
        raise IndexError("polynomials are immutable")


    def __floordiv__(self,right):
        """
        Quotient of division of self by other. This is denoted //.

        If self = quotient \* right + remainder, this function returns
        quotient.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^3 + x + 1
            sage: g = f*(x^2-2) + x
            sage: g.__floordiv__(f)
            x^2 - 2
            sage: g//f
            x^2 - 2
        """
        Q, _ = self.quo_rem(right)
        return Q

    def __mod__(self, other):
        """
        Remainder of division of self by other.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x % (x+1)
            -1
            sage: (x^3 + x - 1) % (x^2 - 1)
            2*x - 1
        """
        _, R = self.quo_rem(other)
        return R

    def _is_atomic(self):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<y> = R[]
            sage: S(x+2)
            x + 2
            sage: S(x+2)._is_atomic()
            False
            sage: S(x)._is_atomic()
            True
        """
        return (self.degree() == self.valuation() and
                self.leading_coefficient()._is_atomic())

    def _mul_generic(self, right):
        if self is right:
            return self._square_generic()
        x = self.list()
        y = right.list()
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t d1 = len(x)-1, d2 = len(y)-1
        if d1 == -1:
            return self
        elif d2 == -1:
            return right
        elif d1 == 0:
            c = x[0]
            return self._parent([c*a for a in y])
        elif d2 == 0:
            c = y[0]
            return self._parent([a*c for a in x])
        coeffs = []
        for k from 0 <= k <= d1+d2:
            start = 0 if k <= d2 else k-d2 # max(0, k-d2)
            end =   k if k <= d1 else d1    # min(k, d1)
            sum = x[start] * y[k-start]
            for i from start < i <= end:
                sum += x[i] * y[k-i]
            coeffs.append(sum)
        return self._parent(coeffs)

    def _square_generic(self):
        x = self.list()
        cdef Py_ssize_t i, j
        cdef Py_ssize_t d = len(x)-1
        zero = self._parent.base_ring().zero_element()
        two = self._parent.base_ring()(2)
        coeffs = [zero] * (2 * d + 1)
        for i from 0 <= i <= d:
            coeffs[2*i] = x[i] * x[i]
            for j from 0 <= j < i:
                coeffs[i+j] += two * x[i] * x[j]
        return self._parent(coeffs)

    def _mul_fateman(self, right):
        r"""
        Returns the product of two polynomials using Kronecker's trick to
        do the multiplication. This could be used over a generic base
        ring.

        .. note::

           -  Since this is implemented in interpreted Python, it could be
              hugely sped up by reimplementing it in Pyrex.

           -  Over the reals there is precision loss, at least in the current
              implementation.


        INPUT:

        -  ``self`` - Polynomial

        -  ``right`` - Polynomial (over same base ring as
           self)


        OUTPUT: Polynomial - The product self\*right.

        ALGORITHM: Based on a paper by R. Fateman

        http://www.cs.berkeley.edu/~fateman/papers/polysbyGMP.pdf

        The idea is to encode dense univariate polynomials as big integers,
        instead of sequences of coefficients. The paper argues that because
        integer multiplication is so cheap, that encoding 2 polynomials to
        big numbers and then decoding the result might be faster than
        popular multiplication algorithms. This seems true when the degree
        is larger than 200.

        EXAMPLES::

            sage: S.<y> = PolynomialRing(RR)
            sage: f = y^10 - 1.393493*y + 0.3
            sage: f._mul_karatsuba(f)
            y^20 - 2.78698600000000*y^11 + 0.600000000000000*y^10 + 1.11022302462516e-16*y^8 - 1.11022302462516e-16*y^6 - 1.11022302462516e-16*y^3 + 1.94182274104900*y^2 - 0.836095800000000*y + 0.0900000000000000
            sage: f._mul_fateman(f)
            y^20 - 2.78698600000000*y^11 + 0.600000000000000*y^10 + 1.94182274104900*y^2 - 0.836095800000000*y + 0.0900000000000000

        Advantages:


        -  Faster than Karatsuba over `\QQ` and
           `\ZZ` (but much slower still than calling NTL's
           optimized C++ implementation, which is the default over
           `\ZZ`)

        -  Potentially less complicated.


        Drawbacks:


        -  Slower over R when the degree of both of polynomials is less
           than 250 (roughly).

        -  Over R, results may not be as accurate as the Karatsuba case.
           This is because we represent coefficients of polynomials over R as
           fractions, then convert them back to floating-point numbers.


        AUTHORS:

        - Didier Deshommes (2006-05-25)
        """
        return self.parent()(polynomial_fateman._mul_fateman_mul(self,right))

    def _mul_karatsuba(self, right):
        r"""
        Returns the product of two polynomials using the Karatsuba divide
        and conquer multiplication algorithm. This is only used over a
        generic base ring. (Special libraries like NTL are used, e.g., for
        the integers and rationals, which are much faster.)

        INPUT:

          - ``self`` - Polynomial
          - ``right`` - Polynomial (over same base ring as
            self)

        OUTPUT: Polynomial - The product self\*right.

        ALGORITHM: The basic idea is to use that

        .. math::

                            (aX + b) (cX + d) = acX^2 + ((a+b)(c+d)-ac-bd)X + bd


        where ac=a\*c and bd=b\*d, which requires three multiplications
        instead of the naive four. (In my examples, strangely just doing
        the above with four multiplications does tend to speed things up
        noticeably.) Given f and g of arbitrary degree bigger than one, let
        e be min(deg(f),deg(g))/2. Write

        .. math::

                               f = a X^e + b   \text{ and }   g = c X^e + d


        and use the identity

        .. math::

                              (aX^e + b) (cX^e + d) = ac X^{2e} +((a+b)(c+d) - ac - bd)X^e + bd


        to recursively compute `fg`.

        TIMINGS: On a Pentium M 1.8Ghz laptop::

          f=R.random(1000,bound=100)
          g=R.random(1000,bound=100)
          time h=f._mul_karatsuba(g)
          Time: 0.42 seconds

        The naive multiplication algorithm takes 14.58 seconds. In
        contrast, MAGMA does this sort of product almost instantly, and can
        easily deal with degree 5000. Basically MAGMA is 100 times faster
        at polynomial multiplication.

        Over Z using NTL, multiplying two polynomials constructed using
        R.random(10000,bound=100) takes 0.10 seconds. Using MAGMA V2.11-10
        the same takes 0.14 seconds. So in this case NTL is somewhat faster
        than MAGMA.

        Over Q using PARI, multiplying two polynomials constructed using
        R.random(10000,bound=100) takes 1.23 seconds. Not good! TODO: use
        NTL polynomials over Z with a denominator instead of PARI.

        .. note::

           - Karatsuba multiplication of polynomials is also
             implemented in PARI in src/basemath/polarit3.c

           - The MAGMA documentation appears to give no information
             about how polynomial multiplication is implemented.
        """
        return self._parent(do_karatsuba(self.list(), right.list()))

    def base_ring(self):
        """
        Return the base ring of the parent of self.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x.base_ring()
            Integer Ring
            sage: (2*x+3).base_ring()
            Integer Ring
        """
        return self.parent().base_ring()

    def base_extend(self, R):
        """
        Return a copy of this polynomial but with coefficients in R, if
        there is a natural map from coefficient ring of self to R.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 - 17*x + 3
            sage: f.base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no such base extension
            sage: f.change_ring(GF(7))
            x^3 + 4*x + 3
        """
        S = self.parent().base_extend(R)
        return S(self)

    def change_variable_name(self, var):
        """
        Return a new polynomial over the same base ring but in a different
        variable.

        EXAMPLES::

            sage: x = polygen(QQ,'x')
            sage: f = -2/7*x^3 + (2/3)*x - 19/993; f
            -2/7*x^3 + 2/3*x - 19/993
            sage: f.change_variable_name('theta')
            -2/7*theta^3 + 2/3*theta - 19/993
        """
        R = self.parent().base_ring()[var]
        return R(self.list())

    def change_ring(self, R):
        """
        Return a copy of this polynomial but with coefficients in R, if at
        all possible.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: f = K.defining_polynomial()
            sage: f.change_ring(GF(7))
            x^2 + x + 1
        """
        S = self.parent().change_ring(R)
        return S(self)

    def _mpoly_dict_recursive(self, variables=None, base_ring=None):
        """
        Return a dict of coefficient entries suitable for construction of a
        MPolynomial_polydict with the given variables.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0)._mpoly_dict_recursive()
            {}
            sage: f = 7*x^5 + x^2 - 2*x - 3
            sage: f._mpoly_dict_recursive()
            {(0,): -3, (1,): -2, (5,): 7, (2,): 1}
        """
        if not self:
            return {}

        var = self.parent().variable_name()
        if variables is None:
            variables = self.parent().variable_names_recursive()
        if not var in variables:
            x = base_ring(self) if base_ring else self
            const_ix = ETuple((0,)*len(variables))
            return { const_ix: x }

        prev_variables = variables[:list(variables).index(var)]
        const_ix = ETuple((0,)*len(prev_variables))
        mpolys = None

        if len(prev_variables) > 0:
            try:
                mpolys = [a._mpoly_dict_recursive(prev_variables, base_ring) for a in self]
            except AttributeError, msg:
                pass

        if mpolys is None:
            if base_ring is not None and base_ring is not self.base_ring():
                mpolys = [{const_ix:base_ring(a)} if a else {} for a in self]
            else:
                mpolys = [{const_ix:a} if a else {} for a in self]

        D = {}
        leftovers = (0,) * (len(variables) - len(prev_variables) - 1)
        for k in range(len(mpolys)):
            for i,a in mpolys[k].iteritems():
                j = ETuple((k,) + leftovers)
                D[i + j] = a

        return D

    def __copy__(self):
        """
        Return a "copy" of self. This is just self, since in Sage
        polynomials are immutable this just returns self again.

        EXAMPLES:

        We create the polynomial `f=x+3`, then note that
        the copy is just the same polynomial again, which is fine since
        polynomials are immutable.

        ::

            sage: x = ZZ['x'].0
            sage: f = x + 3
            sage: g = copy(f)
            sage: g is f
            True
        """
        return self

    def degree(self, gen=None):
        """
        Return the degree of this polynomial. The zero polynomial has
        degree -1.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: f = x^93 + 2*x + 1
            sage: f.degree()
            93
            sage: x = PolynomialRing(QQ, 'x', sparse=True).0
            sage: f = x^100000
            sage: f.degree()
            100000

        ::

            sage: x = QQ['x'].0
            sage: f = 2006*x^2006 - x^2 + 3
            sage: f.degree()
            2006
            sage: f = 0*x
            sage: f.degree()
            -1
            sage: f = x + 33
            sage: f.degree()
            1

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        """
        raise NotImplementedError

    def denominator(self):
        """
        Return a denominator of self.

        First, the lcm of the denominators of the entries of self
        is computed and returned. If this computation fails, the
        unit of the parent of self is returned.

        Note that some subclasses may implement their own
        denominator function. For example, see
        :class:`sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint`

        .. warning::

           This is not the denominator of the rational function
           defined by self, which would always be 1 since self is a
           polynomial.

        EXAMPLES:

        First we compute the denominator of a polynomial with
        integer coefficients, which is of course 1.

        ::

            sage: R.<x> = ZZ[]
            sage: f = x^3 + 17*x + 1
            sage: f.denominator()
            1

        Next we compute the denominator of a polynomial with rational
        coefficients.

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = (1/17)*x^19 - (2/3)*x + 1/3; f
            1/17*x^19 - 2/3*x + 1/3
            sage: f.denominator()
            51

        Finally, we try to compute the denominator of a polynomial with
        coefficients in the real numbers, which is a ring whose elements do
        not have a denominator method.

        ::

            sage: R.<x> = RR[]
            sage: f = x + RR('0.3'); f
            x + 0.300000000000000
            sage: f.denominator()
            1.00000000000000

        Check that the denominator is an element over the base whenever the base
        has no denominator function. This closes #9063.

        ::

            sage: R.<a> = GF(5)[]
            sage: x = R(0)
            sage: x.denominator()
            1
            sage: type(x.denominator())
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: isinstance(x.numerator() / x.denominator(), Polynomial)
            True
            sage: isinstance(x.numerator() / R(1), Polynomial)
            False
        """

        if self.degree() == -1:
            return self.base_ring().one_element()
        #This code was in the original algorithm, but seems irrelevant
        #R = self.base_ring()
        x = self.list()
        try:
            d = x[0].denominator()
            for y in x:
                d = d.lcm(y.denominator())
            #If we return self.parent(d) instead, automatic test
            #start to fail, ex. "sage/rings/polynomial/complex_roots.py"
            return d
        except(AttributeError):
            return self.base_ring().one_element()

    def numerator(self):
        """
        Return a numerator of self computed as self * self.denominator()

        Note that some subclases may implement its own numerator
        function. For example, see
        :class:`sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint`

        .. warning::

          This is not the numerator of the rational function
          defined by self, which would always be self since self is a
          polynomial.

        EXAMPLES:

        First we compute the numerator of a polynomial with
        integer coefficients, which is of course self.

        ::

            sage: R.<x> = ZZ[]
            sage: f = x^3 + 17*x + 1
            sage: f.numerator()
            x^3 + 17*x + 1
            sage: f == f.numerator()
            True

        Next we compute the numerator of a polynomial with rational
        coefficients.

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = (1/17)*x^19 - (2/3)*x + 1/3; f
            1/17*x^19 - 2/3*x + 1/3
            sage: f.numerator()
            3*x^19 - 34*x + 17
            sage: f == f.numerator()
            False

        We try to compute the denominator of a polynomial with
        coefficients in the real numbers, which is a ring whose elements do
        not have a denominator method.

        ::

            sage: R.<x> = RR[]
            sage: f = x + RR('0.3'); f
            x + 0.300000000000000
            sage: f.numerator()
            x + 0.300000000000000

        We check that the computation the numerator and denominator
        are valid

        ::

            sage: K=NumberField(symbolic_expression('x^3+2'),'a')['s,t']['x']
            sage: f=K.random_element()
            sage: f.numerator() / f.denominator() == f
            True
            sage: R=RR['x']
            sage: f=R.random_element()
            sage: f.numerator() / f.denominator() == f
            True
        """
        return self * self.denominator()

    def derivative(self, *args):
        r"""
        The formal derivative of this polynomial, with respect to variables
        supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        .. seealso::

           :meth:`_derivative`

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: g = -x^4 + x^2/2 - x
            sage: g.derivative()
            -4*x^3 + x - 1
            sage: g.derivative(x)
            -4*x^3 + x - 1
            sage: g.derivative(x, x)
            -12*x^2 + 1
            sage: g.derivative(x, 2)
            -12*x^2 + 1

        ::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = PolynomialRing(R)
            sage: f = t^3*x^2 + t^4*x^3
            sage: f.derivative()
            3*t^4*x^2 + 2*t^3*x
            sage: f.derivative(x)
            3*t^4*x^2 + 2*t^3*x
            sage: f.derivative(t)
            4*t^3*x^3 + 3*t^2*x^2
        """
        return multi_derivative(self, args)

    # add .diff(), .differentiate() as aliases for .derivative()
    diff = differentiate = derivative

    def _derivative(self, var=None):
        r"""
        Return the formal derivative of this polynomial with respect to the
        variable var.

        If var is the generator of this polynomial ring (or the default
        value None), this is the usual formal derivative.

        Otherwise, _derivative(var) is called recursively for each of the
        coefficients of this polynomial.

        .. seealso::

           :meth:`derivative`

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0)._derivative()
            0
            sage: parent(R(0)._derivative())
            Univariate Polynomial Ring in x over Integer Ring

        ::

            sage: f = 7*x^5 + x^2 - 2*x - 3
            sage: f._derivative()
            35*x^4 + 2*x - 2
            sage: f._derivative(None)
            35*x^4 + 2*x - 2
            sage: f._derivative(x)
            35*x^4 + 2*x - 2

        In the following example, it doesn't recognise 2\*x as the
        generator, so it tries to differentiate each of the coefficients
        with respect to 2\*x, which doesn't work because the integer
        coefficients don't have a _derivative() method::

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.integer.Integer' object has no attribute '_derivative'

        Examples illustrating recursive behaviour::

            sage: R.<x> = ZZ[]
            sage: S.<y> = PolynomialRing(R)
            sage: f = x^3 + y^3
            sage: f._derivative()
            3*y^2
            sage: f._derivative(y)
            3*y^2
            sage: f._derivative(x)
            3*x^2

        ::

            sage: R = ZZ['x']
            sage: S = R.fraction_field(); x = S.gen()
            sage: R(1).derivative(R(x))
            0
        """
        if var is not None and var != self._parent.gen():
            # call _derivative() recursively on coefficients
            return self._parent([coeff._derivative(var) for coeff in self.list()])

        # compute formal derivative with respect to generator
        if self.is_zero():
            return self
        cdef Py_ssize_t n, degree = self.degree()
        if degree == 0:
            return self.parent().zero_element()
        coeffs = self.list()
        return self._parent([n*coeffs[n] for n from 1 <= n <= degree])

    def integral(self,var=None):
        """
        Return the integral of this polynomial.

        By default, the integration variable is the variable of the
        polynomial.

        Otherwise, the integration variable is the optional parameter ``var``

        .. NOTE::

            The integral is always chosen so the constant term is 0.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0).integral()
            0
            sage: f = R(2).integral(); f
            2*x

        Note that the integral lives over the fraction field of the
        scalar coefficients::

            sage: f.parent()
            Univariate Polynomial Ring in x over Rational Field
            sage: R(0).integral().parent()
            Univariate Polynomial Ring in x over Rational Field

            sage: f = x^3 + x - 2
            sage: g = f.integral(); g
            1/4*x^4 + 1/2*x^2 - 2*x
            sage: g.parent()
            Univariate Polynomial Ring in x over Rational Field

        This shows that the issue at :trac:`7711` is resolved::

            sage: P.<x,z> = PolynomialRing(GF(2147483647))
            sage: Q.<y> = PolynomialRing(P)
            sage: p=x+y+z
            sage: p.integral()
            -1073741823*y^2 + (x + z)*y

            sage: P.<x,z> = PolynomialRing(GF(next_prime(2147483647)))
            sage: Q.<y> = PolynomialRing(P)
            sage: p=x+y+z
            sage: p.integral()
            1073741830*y^2 + (x + z)*y

        A truly convoluted example::

            sage: A.<a1, a2> = PolynomialRing(ZZ)
            sage: B.<b> = PolynomialRing(A)
            sage: C.<c> = PowerSeriesRing(B)
            sage: R.<x> = PolynomialRing(C)
            sage: f = a2*x^2 + c*x - a1*b
            sage: f.parent()
            Univariate Polynomial Ring in x over Power Series Ring in c
            over Univariate Polynomial Ring in b over Multivariate Polynomial
            Ring in a1, a2 over Integer Ring
            sage: f.integral()
            1/3*a2*x^3 + 1/2*c*x^2 - a1*b*x
            sage: f.integral().parent()
            Univariate Polynomial Ring in x over Power Series Ring in c
            over Univariate Polynomial Ring in b over Multivariate Polynomial
            Ring in a1, a2 over Rational Field
            sage: g = 3*a2*x^2 + 2*c*x - a1*b
            sage: g.integral()
            a2*x^3 + c*x^2 - a1*b*x
            sage: g.integral().parent()
            Univariate Polynomial Ring in x over Power Series Ring in c
            over Univariate Polynomial Ring in b over Multivariate Polynomial
            Ring in a1, a2 over Rational Field

        Integration with respect to a variable in the base ring::

            sage: R.<x> = QQ[]
            sage: t = PolynomialRing(R,'t').gen()
            sage: f = x*t +5*t^2
            sage: f.integral(x)
            5*x*t^2 + 1/2*x^2*t
        """
        if var is not None and var != self._parent.gen():
            # call integral() recursively on coefficients
            return self._parent([coeff.integral(var) for coeff in self.list()])
        cdef Py_ssize_t n, degree = self.degree()
        R = self.parent()
        Q = (self.constant_coefficient()/1).parent()
        coeffs = self.list()
        v = [0] + [coeffs[n]/(n+1) for n from 0 <= n <= degree]
        S = R.change_ring(Q)
        return S(v)

    def dict(self):
        """
        Return a sparse dictionary representation of this univariate
        polynomial.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + -1/7*x + 13
            sage: f.dict()
            {0: 13, 1: -1/7, 3: 1}
        """
        X = {}
        Y = self.list()
        for i in xrange(len(Y)):
            c = Y[i]
            if c:
                X[i] = c
        return X

    def factor(self, proof=True):
        r"""
        Return the factorization of self over the base ring of this
        polynomial. Factoring polynomials over
        `\ZZ/n\ZZ` for `n` composite is at
        the moment not implemented.

        INPUT: a polynomial

        OUTPUT:


        -  ``Factorization`` - the factorization of self, which
           is a product of a unit with a product of powers of irreducible
           factors.


        Over a field the irreducible factors are all monic.

        EXAMPLES:

        We factor some polynomials over `\QQ`::

            sage: x = QQ['x'].0
            sage: f = (x^3 - 1)^2
            sage: f.factor()
            (x - 1)^2 * (x^2 + x + 1)^2

        Notice that over the field `\QQ` the irreducible factors are
        monic::

            sage: f = 10*x^5 - 1
            sage: f.factor()
            (10) * (x^5 - 1/10)
            sage: f = 10*x^5 - 10
            sage: f.factor()
            (10) * (x - 1) * (x^4 + x^3 + x^2 + x + 1)

        Over `\ZZ` the irreducible factors need not be monic::

            sage: x = ZZ['x'].0
            sage: f = 10*x^5 - 1
            sage: f.factor()
            10*x^5 - 1

        We factor a non-monic polynomial over the finite field
        `F_{25}`::

            sage: k.<a> = GF(25)
            sage: R.<x> = k[]
            sage: f = 2*x^10 + 2*x + 2*a
            sage: F = f.factor(); F
            (2) * (x + a + 2) * (x^2 + 3*x + 4*a + 4) * (x^2 + (a + 1)*x + a + 2) * (x^5 + (3*a + 4)*x^4 + (3*a + 3)*x^3 + 2*a*x^2 + (3*a + 1)*x + 3*a + 1)

        Notice that the unit factor is included when we multiply `F`
        back out::

            sage: expand(F)
            2*x^10 + 2*x + 2*a

        A new ring.  In the example below, we add the special method
        _factor_univariate_polynomial to the base ring, and observe
        that this method is called instead to factor univariate
        polynomials over this ring.  This facility can be used to
        easily extend polynomial factorization to work over new rings
        you introduce. ::

             sage: R.<x> = QQ[]
             sage: (x^2 + 1).factor()
             x^2 + 1
             sage: QQ._factor_univariate_polynomial = lambda f, proof: f.change_ring(CDF).factor()
             sage: fz = (x^2 + 1).factor(); fz # random order of factors, with noise
             (x - ... + I) * (x - I)
             sage: # Change noisy zero term which affects the order of factors:
             sage: Factorization(
             ...       [ ( parent(f)(
             ...               [ (0.0,im) if abs(re)<=1e-16 else (re,im)
             ...                 for re,im in f
             ...               ]),
             ...           e )
             ...         for f,e in fz
             ...       ])
             (x - I) * (x + I)
             sage: # Clean up:
             sage: del QQ._factor_univariate_polynomial

        Arbitrary precision real and complex factorization::

            sage: R.<x> = RealField(100)[]
            sage: F = factor(x^2-3); F
            (x - 1.7320508075688772935274463415) * (x + 1.7320508075688772935274463415)
            sage: expand(F)
            x^2 - 3.0000000000000000000000000000
            sage: factor(x^2 + 1)
            x^2 + 1.0000000000000000000000000000
            sage: C = ComplexField(100)
            sage: R.<x> = C[]
            sage: F = factor(x^2+3); F
            (x - 1.7320508075688772935274463415*I) * (x + 1.7320508075688772935274463415*I)
            sage: expand(F)
            x^2 + 3.0000000000000000000000000000
            sage: factor(x^2+1)
            (x - I) * (x + I)
            sage: f = C.0 * (x^2 + 1) ; f
            I*x^2 + I
            sage: F = factor(f); F
            (1.0000000000000000000000000000*I) * (x - I) * (x + I)
            sage: expand(F)
            I*x^2 + I

        Over a number field::

            sage: K.<z> = CyclotomicField(15)
            sage: x = polygen(K)
            sage: f = (x^3 + z*x + 1)^3*(x - z); f.factor()
            (x - z) * (x^3 + z*x + 1)^3
            sage: cyclotomic_polynomial(12).change_ring(K).factor()
            (x^2 - z^5 - 1) * (x^2 + z^5)
            sage: f = (x^3 + z*x + 1)^3*(x/(z+2) - 1/3); f.factor()
            (-1/331*z^7 + 3/331*z^6 - 6/331*z^5 + 11/331*z^4 - 21/331*z^3 + 41/331*z^2 - 82/331*z + 165/331) * (x - 1/3*z - 2/3) * (x^3 + z*x + 1)^3

        Over a relative number field::

            sage: x = polygen(QQ)
            sage: K.<z> = CyclotomicField(3)
            sage: L.<a> = K.extension(x^3 - 2)
            sage: t = polygen(L, 't')
            sage: f = (t^3 + t + a)*(t^5 + t + z); f
            t^8 + t^6 + a*t^5 + t^4 + z*t^3 + t^2 + (a + z)*t + z*a
            sage: f.factor()
            (t^3 + t + a) * (t^5 + t + z)

        Over the real double field::

            sage: x = polygen(RDF)
            sage: f = (x-1)^3
            sage: f.factor() # random output (unfortunately)
            (1.0*x - 1.00000859959) * (1.0*x^2 - 1.99999140041*x + 0.999991400484)
            sage: (-2*x^2 - 1).factor()
            (-2.0) * (x^2 + 0.5)
            sage: (-2*x^2 - 1).factor().expand()
            -2.0*x^2 - 1.0

        Note that this factorization suffers from the roots function::

            sage: f.roots() # random output (unfortunately)
            [1.00000859959, 0.999995700205 + 7.44736245561e-06*I, 0.999995700205 - 7.44736245561e-06*I]

        Over the complex double field. Because this is approximate, all
        factors will occur with multiplicity 1::

            sage: x = CDF['x'].0; i = CDF.0
            sage: f = (x^2 + 2*i)^3
            sage: f.factor()    # random low order bits
            (1.0*x + -0.999994409957 + 1.00001040378*I) * (1.0*x + -0.999993785062 + 0.999989956987*I) * (1.0*x + -1.00001180498 + 0.999999639235*I) * (1.0*x + 0.999995530902 - 0.999987780431*I) * (1.0*x + 1.00001281704 - 1.00000223945*I) * (1.0*x + 0.999991652054 - 1.00000998012*I)
            sage: f(-f.factor()[0][0][0])   # random low order bits
            -2.38358052913e-14 - 2.57571741713e-14*I

        Factoring polynomials over `\ZZ/n\ZZ` for
        composite `n` is not implemented::

            sage: R.<x> = PolynomialRing(Integers(35))
            sage: f = (x^2+2*x+2)*(x^2+3*x+9)
            sage: f.factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: factorization of polynomials over rings with composite characteristic is not implemented

        Factoring polynomials over QQbar and AA is implemented (see #8544)::

            sage: x = polygen(QQbar)
            sage: (x^8-1).factor()
            (x - 1) * (x - 0.7071067811865475? - 0.7071067811865475?*I) * (x - 0.7071067811865475? + 0.7071067811865475?*I) * (x - I) * (x + I) * (x + 0.7071067811865475? - 0.7071067811865475?*I) * (x + 0.7071067811865475? + 0.7071067811865475?*I) * (x + 1)
            sage: (12*x^2-4).factor()
            (12) * (x - 0.5773502691896258?) * (x + 0.5773502691896258?)
            sage: x.parent()(-1).factor()
            -1
            sage: EllipticCurve('11a1').change_ring(QQbar).division_polynomial(5).factor()
            (5) * (x - 16) * (x - 5) * (x - 1.959674775249769?) * (x - 1.427050983124843? - 3.665468789467727?*I) * (x - 1.427050983124843? + 3.665468789467727?*I) * (x + 0.9549150281252629? - 0.8652998037182486?*I) * (x + 0.9549150281252629? + 0.8652998037182486?*I) * (x + 1.927050983124843? - 1.677599044300515?*I) * (x + 1.927050983124843? + 1.677599044300515?*I) * (x + 2.959674775249769?) * (x + 6.545084971874737? - 7.106423590645660?*I) * (x + 6.545084971874737? + 7.106423590645660?*I)

        ::

            sage: x = polygen(AA)
            sage: (x^8+1).factor()
            (x^2 - 1.847759065022574?*x + 1.000000000000000?) * (x^2 - 0.7653668647301795?*x + 1.000000000000000?) * (x^2 + 0.7653668647301795?*x + 1.000000000000000?) * (x^2 + 1.847759065022574?*x + 1.000000000000000?)
            sage: x.parent()(3).factor()
            3
            sage: (12*x^2-4).factor()
            (12) * (x - 0.5773502691896258?) * (x + 0.5773502691896258?)
            sage: (12*x^2+4).factor()
            (12) * (x^2 + 0.3333333333333334?)
            sage: EllipticCurve('11a1').change_ring(AA).division_polynomial(5).factor()
            (5) * (x - 16.00000000000000?) * (x - 5.000000000000000?) * (x - 1.959674775249769?) * (x + 2.959674775249769?) * (x^2 - 2.854101966249685?*x + 15.47213595499958?) * (x^2 + 1.909830056250526?*x + 1.660606461254312?) * (x^2 + 3.854101966249685?*x + 6.527864045000421?) * (x^2 + 13.09016994374948?*x + 93.33939353874569?)

        TESTS:

        This came up in ticket #7088::

            sage: R.<x>=PolynomialRing(ZZ)
            sage: f = 12*x^10 + x^9 + 432*x^3 + 9011
            sage: g = 13*x^11 + 89*x^3 + 1
            sage: F = f^2 * g^3
            sage: F = f^2 * g^3; F.factor()
            (12*x^10 + x^9 + 432*x^3 + 9011)^2 * (13*x^11 + 89*x^3 + 1)^3
            sage: F = f^2 * g^3 * 7; F.factor()
            7 * (12*x^10 + x^9 + 432*x^3 + 9011)^2 * (13*x^11 + 89*x^3 + 1)^3

        This example came up in ticket #7097::

            sage: x = polygen(QQ)
            sage: f = 8*x^9 + 42*x^6 + 6*x^3 - 1
            sage: g = x^24 - 12*x^23 + 72*x^22 - 286*x^21 + 849*x^20 - 2022*x^19 + 4034*x^18 - 6894*x^17 + 10182*x^16 - 13048*x^15 + 14532*x^14 - 13974*x^13 + 11365*x^12 - 7578*x^11 + 4038*x^10 - 1766*x^9 + 762*x^8 - 408*x^7 + 236*x^6 - 126*x^5 + 69*x^4 - 38*x^3 + 18*x^2 - 6*x + 1
            sage: assert g.is_irreducible()
            sage: K.<a> = NumberField(g)
            sage: len(f.roots(K))
            9
            sage: f.factor()
            (8) * (x^3 + 1/4) * (x^6 + 5*x^3 - 1/2)
            sage: f.change_ring(K).factor()
            (8) * (x - 3260097/3158212*a^22 + 35861067/3158212*a^21 - 197810817/3158212*a^20 + 722970825/3158212*a^19 - 1980508347/3158212*a^18 + 4374189477/3158212*a^17 - 4059860553/1579106*a^16 + 6442403031/1579106*a^15 - 17542341771/3158212*a^14 + 20537782665/3158212*a^13 - 20658463789/3158212*a^12 + 17502836649/3158212*a^11 - 11908953451/3158212*a^10 + 6086953981/3158212*a^9 - 559822335/789553*a^8 + 194545353/789553*a^7 - 505969453/3158212*a^6 + 338959407/3158212*a^5 - 155204647/3158212*a^4 + 79628015/3158212*a^3 - 57339525/3158212*a^2 + 26692783/3158212*a - 1636338/789553) * ...
            sage: f = QQbar['x'](1)
            sage: f.factor()
            1

        Factorization also works even if the variable of the finite
        field is nefariously labeled "x"::

            sage: x = GF(3^2, 'a')['x'].0
            sage: f = x^10 +7*x -13
            sage: G = f.factor(); G
            (x + a) * (x + 2*a + 1) * (x^4 + (a + 2)*x^3 + (2*a + 2)*x + 2) * (x^4 + 2*a*x^3 + (a + 1)*x + 2)
            sage: prod(G) == f
            True

        ::

            sage: f.parent().base_ring()._assign_names(['a'])
            sage: f.factor()
            (x + a) * (x + 2*a + 1) * (x^4 + (a + 2)*x^3 + (2*a + 2)*x + 2) * (x^4 + 2*a*x^3 + (a + 1)*x + 2)

        ::

            sage: k = GF(9,'x')    # purposely calling it x to test robustness
            sage: x = PolynomialRing(k,'x0').gen()
            sage: f = x^3 + x + 1
            sage: f.factor()
            (x0 + 2) * (x0 + x) * (x0 + 2*x + 1)
            sage: f = 0*x
            sage: f.factor()
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined

        ::

            sage: f = x^0
            sage: f.factor()
            1

        Over a complicated number field::

            sage: x = polygen(QQ, 'x')
            sage: f = x^6 + 10/7*x^5 - 867/49*x^4 - 76/245*x^3 + 3148/35*x^2 - 25944/245*x + 48771/1225
            sage: K.<a> = NumberField(f)
            sage: S.<T> = K[]
            sage: ff = S(f); ff
            T^6 + 10/7*T^5 - 867/49*T^4 - 76/245*T^3 + 3148/35*T^2 - 25944/245*T + 48771/1225
            sage: F = ff.factor()
            sage: len(F)
            4
            sage: F[:2]
            [(T - a, 1), (T - 40085763200/924556084127*a^5 - 145475769880/924556084127*a^4 + 527617096480/924556084127*a^3 + 1289745809920/924556084127*a^2 - 3227142391585/924556084127*a - 401502691578/924556084127, 1)]
            sage: expand(F)
            T^6 + 10/7*T^5 - 867/49*T^4 - 76/245*T^3 + 3148/35*T^2 - 25944/245*T + 48771/1225

        ::

            sage: f = x^2 - 1/3 ; K.<a> = NumberField(f) ; A.<T> = K[] ; g = A(x^2-1)
            sage: g.factor()
            (T - 1) * (T + 1)

        ::

            sage: h = A(3*x^2-1) ; h.factor()
            (3) * (T - a) * (T + a)

        ::

            sage: h = A(x^2-1/3) ; h.factor()
            (T - a) * (T + a)

        Test that ticket #10279 is fixed::

            sage: R.<t> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(t^4 - t^2 + 1)
            sage: pol = t^3 + (-4*a^3 + 2*a)*t^2 - 11/3*a^2*t + 2/3*a^3 - 4/3*a
            sage: pol.factor()
            (t - 2*a^3 + a) * (t - 4/3*a^3 + 2/3*a) * (t - 2/3*a^3 + 1/3*a)

        Test that this factorization really uses ``nffactor()`` internally::

            sage: pari.default("debug", 3)
            sage: F = pol.factor()
            <BLANKLINE>
            Entering nffactor:
            ...
            sage: pari.default("debug", 0)

        Test that ticket #10369 is fixed::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
            sage: R.<t> = PolynomialRing(K)

            sage: pol = (-1/7*a^5 - 1/7*a^4 - 1/7*a^3 - 1/7*a^2 - 2/7*a - 1/7)*t^10 + (4/7*a^5 - 2/7*a^4 - 2/7*a^3 - 2/7*a^2 - 2/7*a - 6/7)*t^9 + (90/49*a^5 + 152/49*a^4 + 18/49*a^3 + 24/49*a^2 + 30/49*a + 36/49)*t^8 + (-10/49*a^5 + 10/7*a^4 + 198/49*a^3 - 102/49*a^2 - 60/49*a - 26/49)*t^7 + (40/49*a^5 + 45/49*a^4 + 60/49*a^3 + 277/49*a^2 - 204/49*a - 78/49)*t^6 + (90/49*a^5 + 110/49*a^4 + 2*a^3 + 80/49*a^2 + 46/7*a - 30/7)*t^5 + (30/7*a^5 + 260/49*a^4 + 250/49*a^3 + 232/49*a^2 + 32/7*a + 8)*t^4 + (-184/49*a^5 - 58/49*a^4 - 52/49*a^3 - 66/49*a^2 - 72/49*a - 72/49)*t^3 + (18/49*a^5 - 32/49*a^4 + 10/49*a^3 + 4/49*a^2)*t^2 + (2/49*a^4 - 4/49*a^3 + 2/49*a^2)*t
            sage: pol.factor()
            (-1/7*a^5 - 1/7*a^4 - 1/7*a^3 - 1/7*a^2 - 2/7*a - 1/7) * t * (t - a^5 - a^4 - a^3 - a^2 - a - 1)^4 * (t^5 + (-12/7*a^5 - 10/7*a^4 - 8/7*a^3 - 6/7*a^2 - 4/7*a - 2/7)*t^4 + (12/7*a^5 - 8/7*a^3 + 16/7*a^2 + 2/7*a + 20/7)*t^3 + (-20/7*a^5 - 20/7*a^3 - 20/7*a^2 + 4/7*a - 2)*t^2 + (12/7*a^5 + 12/7*a^3 + 2/7*a + 16/7)*t - 4/7*a^5 - 4/7*a^3 - 4/7*a - 2/7)

            sage: pol = (1/7*a^2 - 1/7*a)*t^10 + (4/7*a - 6/7)*t^9 + (102/49*a^5 + 99/49*a^4 + 96/49*a^3 + 93/49*a^2 + 90/49*a + 150/49)*t^8 + (-160/49*a^5 - 36/49*a^4 - 48/49*a^3 - 8/7*a^2 - 60/49*a - 60/49)*t^7 + (30/49*a^5 - 55/49*a^4 + 20/49*a^3 + 5/49*a^2)*t^6 + (6/49*a^4 - 12/49*a^3 + 6/49*a^2)*t^5
            sage: pol.factor()
            (1/7*a^2 - 1/7*a) * t^5 * (t^5 + (-40/7*a^5 - 38/7*a^4 - 36/7*a^3 - 34/7*a^2 - 32/7*a - 30/7)*t^4 + (60/7*a^5 - 30/7*a^4 - 18/7*a^3 - 9/7*a^2 - 3/7*a)*t^3 + (60/7*a^4 - 40/7*a^3 - 16/7*a^2 - 4/7*a)*t^2 + (30/7*a^3 - 25/7*a^2 - 5/7*a)*t + 6/7*a^2 - 6/7*a)

            sage: pol = x^10 + (4/7*a - 6/7)*x^9 + (9/49*a^2 - 3/7*a + 15/49)*x^8 + (8/343*a^3 - 32/343*a^2 + 40/343*a - 20/343)*x^7 + (5/2401*a^4 - 20/2401*a^3 + 40/2401*a^2 - 5/343*a + 15/2401)*x^6 + (-6/16807*a^4 + 12/16807*a^3 - 18/16807*a^2 + 12/16807*a - 6/16807)*x^5
            sage: pol.factor()
            x^5 * (x^5 + (4/7*a - 6/7)*x^4 + (9/49*a^2 - 3/7*a + 15/49)*x^3 + (8/343*a^3 - 32/343*a^2 + 40/343*a - 20/343)*x^2 + (5/2401*a^4 - 20/2401*a^3 + 40/2401*a^2 - 5/343*a + 15/2401)*x - 6/16807*a^4 + 12/16807*a^3 - 18/16807*a^2 + 12/16807*a - 6/16807)

        Factoring over a number field over which we cannot factor the
        discriminant by trial division::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^16 - x - 6)
            sage: R.<x> = PolynomialRing(K)
            sage: f = (x+a)^50 - (a-1)^50
            sage: len(factor(f))
            6
            sage: pari(K.discriminant()).factor(limit=0)
            [-1, 1; 3, 15; 23, 1; 887, 1; 12583, 1; 2354691439917211, 1]
            sage: factor(K.discriminant())
            -1 * 3^15 * 23 * 887 * 12583 * 6335047 * 371692813

        Factoring over a number field over which we cannot factor the
        discriminant and over which `nffactor()` fails::

            sage: p = next_prime(10^50); q = next_prime(10^51); n = p*q;
            sage: K.<a> = QuadraticField(p*q)
            sage: R.<x> = PolynomialRing(K)
            sage: K.pari_polynomial('a').nffactor("x^2+1")
            Traceback (most recent call last):
            ...
            PariError: precision too low in floorr (precision loss in truncation)
            sage: factor(x^2 + 1)
            x^2 + 1
            sage: factor( (x - a) * (x + 2*a) )
            (x - a) * (x + 2*a)

        A test where nffactor used to fail without a nf structure::

            sage: x = polygen(QQ)
            sage: K = NumberField([x^2-1099511627777, x^3-3],'a')
            sage: x = polygen(K)
            sage: f = x^3 - 3
            sage: factor(f)
            (x - a1) * (x^2 + a1*x + a1^2)
        """
        # PERFORMANCE NOTE:
        #     In many tests with SMALL degree PARI is substantially
        #     better than NTL.  (And magma is better yet.)  And the
        #     timing difference has nothing to do with moving Python
        #     data to NTL and back.
        #     For large degree ( > 1500) in the one test I tried, NTL was
        #     *much* better than MAGMA, and far better than PARI.  So probably
        #     NTL's implementation is asymptotically better.  I could use
        #     PARI for smaller degree over other rings besides Z, and use
        #     NTL in general.
        # A remark from Bill Hart (2007-09-25) about the above observation:
        ## NTL uses the Berlekamp-Zassenhaus method with van Hoeij's improvements.
        ## But so does Magma since about Jul 2001.
        ##
        ## But here's the kicker. PARI also uses this algorithm. Even Maple uses
        ## it!
        ##
        ## NTL's LLL algorithms are extremely well developed (van Hoeij uses
        ## LLL). There is also a possible speed difference in whether one uses
        ## quadratic convergence or not in the Hensel lift. But the right choice
        ## is not always what one thinks.
        ##
        ## But more than likely NTL is just better for large problems because
        ## Victor Shoup was very careful with the choice of strategies and
        ## parameters he used. Paul Zimmerman supplied him with a pile of
        ## polynomials to factor for comparison purposes and these seem to have
        ## been used to tune the algorithm for a wide range of inputs, including
        ## cases that van Hoeij's algorithm doesn't usually like.
        ##
        ## If you have a bound on the coefficients of the factors, one can surely
        ## do better than a generic implementation, but probably not much better
        ## if there are many factors.
        ##

        ## HUGE TODO, refactor the code below here such that this method will
        ## have as only the following code
        ##
        ## R = self.parent().base_ring()
        ## return R._factor_univariate_polynomial(self)
        ##
        ## in this way we can move the specific logic of factoring to the
        ## self.parent().base_ring() and get rid of all the ugly
        ## is_SomeType(R) checks and get way nicer structured code
        ## 200 lines of spagetti code is just way to much!

        if self.degree() < 0:
            raise ValueError("factorization of 0 not defined")
        if self.degree() == 0:
            return Factorization([], unit=self[0])

        R = self.parent().base_ring()
        if hasattr(R, '_factor_univariate_polynomial'):
            return R._factor_univariate_polynomial(self,proof=proof)

        G = None
        ch = R.characteristic()
        if not (ch == 0 or sage.rings.arith.is_prime(ch)):
            raise NotImplementedError("factorization of polynomials over rings with composite characteristic is not implemented")

        from sage.rings.number_field.all import is_NumberField, \
             is_RelativeNumberField, NumberField
        from sage.rings.finite_rings.constructor import is_FiniteField
        from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
        from sage.rings.integer_ring import is_IntegerRing
        from sage.rings.rational_field import is_RationalField

        n = None

        if is_IntegerModRing(R) or is_IntegerRing(R) or is_RationalField(R):

            try:
                G = list(self._pari_with_name().factor())
            except PariError:
                raise NotImplementedError

        elif is_RelativeNumberField(R):

            M = R.absolute_field('a')
            from_M, to_M = M.structure()
            g = M['x']([to_M(x) for x in self.list()])
            F = g.factor()
            S = self.parent()
            v = [(S([from_M(x) for x in f.list()]), e) for f, e in F]
            return Factorization(v, from_M(F.unit()))

        elif is_FiniteField(R):
            v = [x._pari_("a") for x in self.list()]
            f = pari(v).Polrev()
            G = list(f.factor())

        elif is_AlgebraicField(R):
            Rx = self.parent()
            return Factorization([(Rx([-r,1]),e) for r,e in self.roots()],
                                 unit=self.leading_coefficient())

        elif is_AlgebraicRealField(R):
            from sage.rings.qqbar import QQbar
            Rx = self.parent()
            rr = self.roots()
            cr = [(r,e) for r,e in self.roots(QQbar) if r.imag()>0]
            return Factorization(
                [(Rx([-r,1]),e) for r,e in rr] +
                [(Rx([r.norm(),-2*r.real(),1]),e) for r,e in cr],
                unit=self.leading_coefficient())


        elif is_NumberField(R):

            if R.degree() == 1:
                factors = self.change_ring(QQ).factor()
                return Factorization([(self._parent(p), e) for p, e in factors], R(factors.unit()))

            # Convert the polynomial we want to factor to PARI
            f = self._pari_with_name()
            try:
                try:
                    # Try to compute the PARI nf structure with important=False.
                    # This will raise RuntimeError if the computation is too
                    # difficult.  It will raise TypeError if the defining
                    # polynomial is not integral.
                    Rpari = R.pari_nf(important=False)
                except RuntimeError:
                    # Cannot easily compute the nf structure, use the defining
                    # polynomial instead.
                    Rpari = R.pari_polynomial("y")
                # nffactor() can fail with PariError "precision too low"
                G = list(Rpari.nffactor(f))
            except (PariError, TypeError):
                # Use factornf() which only needs the defining polynomial,
                # which does not require an integral polynomial and which
                # has no problems with floating-point precision.
                G = list(f.factornf(R.pari_polynomial("y")))
            # PARI's nffactor() ignores the unit, _factor_pari_helper()
            # adds back the unit of the factorization.
            return self._factor_pari_helper(G)

        elif is_RealField(R):
            n = pari.set_real_precision(int(3.5*R.prec()) + 1)
            G = list(self._pari_with_name().factor())

        elif sage.rings.complex_double.is_ComplexDoubleField(R):
            unit = self.leading_coefficient()
            f = (~unit)*self
            roots = f.roots(multiplicities=False)
            assert len(roots) == self.degree() # all roots appear with multiplicity one
            x = self.parent().gen()
            v = [(x - a, 1) for a in roots]
            return Factorization(v, unit)

        elif sage.rings.real_double.is_RealDoubleField(R):
            roots = self.roots(sage.rings.complex_double.CDF, multiplicities=False)
            assert len(roots) == self.degree() # all roots appear with multiplicity one
            G = [[],[]]
            real_roots = []
            non_real_roots = []
            for r in roots:
                if r.imag().is_zero():
                    for i in xrange(len(real_roots)):
                        if real_roots[i][0] == r:
                            real_roots[i][1] += 1
                            r = None
                            break
                    if r is not None:
                        real_roots.append([r,1])
                else:
                    for i in xrange(len(non_real_roots)):
                        if non_real_roots[i][0] == r or non_real_roots[i][0] == r.conj():
                            non_real_roots[i][1] += 1
                            r = None
                            break
                    if r is not None:
                        non_real_roots.append([r,1])
            x = self.parent().objgen()[1]
            for r in real_roots:
                G[0].append( x - r[0].real() )
                G[1].append( r[1] )
            for z in non_real_roots:
                a = ( z[0] + z[0].conj() ).real()
                b = ( z[0]*(z[0].conj()) ).real()
                G[0].append( x**2 - a*x + b )
                assert z[1] % 2 == 0, "Bug in root finding code over RDF"
                G[1].append( z[1] // 2 )

        elif sage.rings.complex_field.is_ComplexField(R):
            # This is a hack to make the polynomial have complex coefficients, since
            # otherwise PARI will factor over RR.
            n = pari.set_real_precision(int(3.5*R.prec()) + 1)
            if self.leading_coefficient() != R.gen():
                G = list((pari(R.gen())*self._pari_with_name()).factor())
            else:
                G = self._pari_with_name().factor()

        #elif padic_field.is_pAdicField(R):
        #    G = list(self._pari_with_name().factorpadic(R.prime(), R.prec()))

        if G is None:
            raise NotImplementedError

        return self._factor_pari_helper(G, n)

    def _factor_pari_helper(self, G, n=None, unit=None):
        """
        Fix up and normalize a factorization that came from PARI.

        TESTS::

            sage: R.<x>=PolynomialRing(ZZ)
            sage: f = (2*x + 1) * (3*x^2 - 5)^2
            sage: f._factor_pari_helper(pari(f).factor())
            (2*x + 1) * (3*x^2 - 5)^2
            sage: f._factor_pari_helper(pari(f).factor(), unit=11)
            11 * (2*x + 1) * (3*x^2 - 5)^2
            sage: (8*f)._factor_pari_helper(pari(f).factor())
            8 * (2*x + 1) * (3*x^2 - 5)^2
            sage: (8*f)._factor_pari_helper(pari(f).factor(), unit=11)
            88 * (2*x + 1) * (3*x^2 - 5)^2
            sage: QQ['x'](f)._factor_pari_helper(pari(f).factor())
            (18) * (x + 1/2) * (x^2 - 5/3)^2
            sage: QQ['x'](f)._factor_pari_helper(pari(f).factor(), unit=11)
            (198) * (x + 1/2) * (x^2 - 5/3)^2

            sage: f = prod((k^2*x^k + k)^(k-1) for k in primes(10))
            sage: F = f._factor_pari_helper(pari(f).factor()); F
            1323551250 * (2*x^2 + 1) * (3*x^3 + 1)^2 * (5*x^5 + 1)^4 * (7*x^7 + 1)^6
            sage: F.prod() == f
            True
            sage: QQ['x'](f)._factor_pari_helper(pari(f).factor())
            (1751787911376562500) * (x^2 + 1/2) * (x^3 + 1/3)^2 * (x^5 + 1/5)^4 * (x^7 + 1/7)^6

            sage: g = GF(19)['x'](f)
            sage: G = g._factor_pari_helper(pari(g).factor()); G
            (4) * (x + 3) * (x + 16)^5 * (x + 11)^6 * (x^2 + 7*x + 9)^4 * (x^2 + 15*x + 9)^4 * (x^3 + 13)^2 * (x^6 + 8*x^5 + 7*x^4 + 18*x^3 + 11*x^2 + 12*x + 1)^6
            sage: G.prod() == g
            True
        """
        pols = G[0]
        exps = G[1]
        R = self.parent()
        F = [(R(f), int(e)) for f, e in zip(pols, exps)]

        if unit is None:
            unit = self.leading_coefficient()
        else:
            unit *= self.leading_coefficient()

        if R.base_ring().is_field():
            # When the base ring is a field we normalize
            # the irreducible factors so they have leading
            # coefficient 1.
            for i, (f, e) in enumerate(F):
                if not f.is_monic():
                    F[i] = (f.monic(), e)

        else:
            # Otherwise we have to adjust for
            # the content ignored by PARI.
            content_fix = R.base_ring().one_element()
            for f, e in F:
                if not f.is_monic():
                    content_fix *= f.leading_coefficient()**e
            unit //= content_fix
            if not unit.is_unit():
                F.append((R(unit), ZZ(1)))
                unit = R.base_ring().one_element()

        if not n is None:
            pari.set_real_precision(n)  # restore precision
        return Factorization(F, unit)

    @coerce_binop
    def lcm(self, other):
        """
        Let f and g be two polynomials. Then this function returns the
        monic least common multiple of f and g.
        """
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q.leading_coefficient())*q

    def _lcm(self, other):
        """
        Let f and g be two polynomials. Then this function returns the
        monic least common multiple of f and g.
        """
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q.leading_coefficient())*q  # make monic  (~ is inverse in python)

    def is_primitive(self, n=None, n_prime_divs=None):
        """
        Returns ``True`` if the polynomial is primitive.  The semantics of
        "primitive" depend on the polynomial coefficients.

        - (field theory) A polynomial of degree `m` over a finite field
          `\GF{q}` is primitive if it is irreducible and its root in
          `\GF{q^m}` generates the multiplicative group `\GF{q^m}^*`.

        - (ring theory) A polynomial over a ring is primitive if its
          coefficients generate the unit ideal.

        Calling `is_primitive` on a polynomial over an infinite field will
        raise an error.

        The additional inputs to this function are to speed up computation for
        field semantics (see note).

        INPUTS:

          - ``n`` (default: ``None``) - if provided, should equal
            `q-1` where ``self.parent()`` is the field with `q`
            elements;  otherwise it will be computed.

          - ``n_prime_divs`` (default: ``None``) - if provided, should
            be a list of the prime divisors of ``n``; otherwise it
            will be computed.

        .. note::

          Computation of the prime divisors of ``n`` can dominate the running
          time of this method, so performing this computation externally
          (e.g. ``pdivs=n.prime_divisors()``) is a good idea for repeated calls
          to is_primitive for polynomials of the same degree.

          Results may be incorrect if the wrong ``n`` and/or factorization are
          provided.

        EXAMPLES::

          Field semantics examples.

          ::

            sage: R.<x> = GF(2)['x']
            sage: f = x^4+x^3+x^2+x+1
            sage: f.is_irreducible(), f.is_primitive()
            (True, False)
            sage: f = x^3+x+1
            sage: f.is_irreducible(), f.is_primitive()
            (True, True)
            sage: R.<x> = GF(3)[]
            sage: f = x^3-x+1
            sage: f.is_irreducible(), f.is_primitive()
            (True, True)
            sage: f = x^2+1
            sage: f.is_irreducible(), f.is_primitive()
            (True, False)
            sage: R.<x> = GF(5)[]
            sage: f = x^2+x+1
            sage: f.is_primitive()
            False
            sage: f = x^2-x+2
            sage: f.is_primitive()
            True
            sage: x=polygen(QQ); f=x^2+1
            sage: f.is_primitive()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_primitive() not defined for polynomials over infinite fields.

          Ring semantics examples.

          ::

            sage: x=polygen(ZZ)
            sage: f = 5*x^2+2
            sage: f.is_primitive()
            True
            sage: f = 5*x^2+5
            sage: f.is_primitive()
            False

            sage: K=NumberField(x^2+5,'a')
            sage: R=K.ring_of_integers()
            sage: a=R.gen(1)
            sage: a^2
            -5
            sage: f=a*x+2
            sage: f.is_primitive()
            True
            sage: f=(1+a)*x+2
            sage: f.is_primitive()
            False

            sage: x=polygen(Integers(10));
            sage: f=5*x^2+2
            sage: #f.is_primitive()  #BUG:: elsewhere in Sage, should return True
            sage: f=4*x^2+2
            sage: #f.is_primitive()  #BUG:: elsewhere in Sage, should return False

        TESTS::

            sage: R.<x> = GF(2)['x']
            sage: f = x^4+x^3+x^2+x+1
            sage: f.is_primitive(15)
            False
            sage: f.is_primitive(15, [3,5])
            False
            sage: f.is_primitive(n_prime_divs=[3,5])
            False
            sage: f = x^3+x+1
            sage: f.is_primitive(7, [7])
            True
            sage: R.<x> = GF(3)[]
            sage: f = x^3-x+1
            sage: f.is_primitive(26, [2,13])
            True
            sage: f = x^2+1
            sage: f.is_primitive(8, [2])
            False
            sage: R.<x> = GF(5)[]
            sage: f = x^2+x+1
            sage: f.is_primitive(24, [2,3])
            False
            sage: f = x^2-x+2
            sage: f.is_primitive(24, [2,3])
            True
            sage: x=polygen(Integers(103)); f=x^2+1
            sage: f.is_primitive()
            False
        """
        R = self.base_ring()
        if R.is_field():
            if not R.is_finite():
                raise NotImplementedError("is_primitive() not defined for polynomials over infinite fields.")

            if not self.is_irreducible():
                return False
            if n is None:
                q = self.base_ring().order()
                n = q ** self.degree() - 1
            y = self.parent().quo(self).gen()
            from sage.groups.generic import order_from_multiple
            return n == order_from_multiple(y, n, n_prime_divs, operation="*")
        else:
            return R.ideal(self.coefficients())==R.ideal(1)

    def is_constant(self):
        """
        Return True if this is a constant polynomial.

        OUTPUT:


        -  ``bool`` - True if and only if this polynomial is
           constant


        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x.is_constant()
            False
            sage: R(2).is_constant()
            True
            sage: R(0).is_constant()
            True
        """
        return self.degree() <= 0

    def is_monomial(self):
        """
        Returns True if self is a monomial, i.e., a power of the generator.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.is_monomial()
            True
            sage: (x+1).is_monomial()
            False
            sage: (x^2).is_monomial()
            True
            sage: R(1).is_monomial()
            True

        The coefficient must be 1::

            sage: (2*x^5).is_monomial()
            False

        To allow a non-1 leading coefficient, use is_term()::

            sage: (2*x^5).is_term()
            True

        .. warning::

           The definition of is_monomial in Sage up to 4.7.1 was the
           same as is_term, i.e., it allowed a coefficient not equal
           to 1.
        """
        return len(self.exponents()) == 1 and self.leading_coefficient() == 1

    def is_term(self):
        """
        Return True if self is an element of the base ring times a
        power of the generator.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.is_term()
            True
            sage: R(1).is_term()
            True
            sage: (3*x^5).is_term()
            True
            sage: (1+3*x^5).is_term()
            False

        To require that the coefficient is 1, use is_monomial() instead::

            sage: (3*x^5).is_monomial()
            False
        """
        return len(self.exponents()) == 1

    def root_field(self, names, check_irreducible=True):
        """
        Return the field generated by the roots of the irreducible
        polynomial self. The output is either a number field, relative
        number field, a quotient of a polynomial ring over a field, or the
        fraction field of the base ring.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('a')
            Number Field in a with defining polynomial x^3 + x + 17

        ::

            sage: R.<x> = QQ['x']
            sage: f = x - 3
            sage: f.root_field('b')
            Rational Field

        ::

            sage: R.<x> = ZZ['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('b')
            Number Field in b with defining polynomial x^3 + x + 17

        ::

            sage: y = QQ['x'].0
            sage: L.<a> = NumberField(y^3-2)
            sage: R.<x> = L['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('c')
            Number Field in c with defining polynomial x^3 + x + 17 over its base field

        ::

            sage: R.<x> = PolynomialRing(GF(9,'a'))
            sage: f = x^3 + x^2 + 8
            sage: K.<alpha> = f.root_field(); K
            Univariate Quotient Polynomial Ring in alpha over Finite Field in a of size 3^2 with modulus x^3 + x^2 + 2
            sage: alpha^2 + 1
            alpha^2 + 1
            sage: alpha^3 + alpha^2
            1

        ::

            sage: R.<x> = QQ[]
            sage: f = x^2
            sage: K.<alpha> = f.root_field()
            Traceback (most recent call last):
            ...
            ValueError: polynomial must be irreducible

        TESTS::

            sage: (PolynomialRing(Integers(31),name='x').0+5).root_field('a')
            Ring of integers modulo 31
        """
        from sage.rings.number_field.number_field import is_NumberField, NumberField

        R = self.base_ring()
        if not R.is_integral_domain():
            raise ValueError("the base ring must be a domain")

        if check_irreducible and not self.is_irreducible():
            raise ValueError("polynomial must be irreducible")

        if self.degree() <= 1:
            return R.fraction_field()

        if sage.rings.integer_ring.is_IntegerRing(R):
            return NumberField(self, names)

        if sage.rings.rational_field.is_RationalField(R) or is_NumberField(R):
            return NumberField(self, names)

        return R.fraction_field()[self.parent().variable_name()].quotient(self, names)

    def sylvester_matrix(self, right, variable = None):
        """
        Returns the Sylvester matrix of self and right.

        Note that the Sylvester matrix is not defined if one of the polynomials
        is zero.

        INPUT:

        - right: a polynomial in the same ring as self.
        - variable: optional, included for compatibility with the multivariate
          case only. The variable of the polynomials.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = (6*x + 47)*(7*x^2 - 2*x + 38)
            sage: g = (6*x + 47)*(3*x^3 + 2*x + 1)
            sage: M = f.sylvester_matrix(g)
            sage: M
            [  42  317  134 1786    0    0    0]
            [   0   42  317  134 1786    0    0]
            [   0    0   42  317  134 1786    0]
            [   0    0    0   42  317  134 1786]
            [  18  141   12  100   47    0    0]
            [   0   18  141   12  100   47    0]
            [   0    0   18  141   12  100   47]

        If the polynomials share a non-constant common factor then the
        determinant of the Sylvester matrix will be zero::

            sage: M.determinant()
            0

        If self and right are polynomials of positive degree, the determinant
        of the Sylvester matrix is the resultant of the polynomials.::

            sage: h1 = R.random_element()
            sage: h2 = R.random_element()
            sage: M1 = h1.sylvester_matrix(h2)
            sage: M1.determinant() == h1.resultant(h2)
            True

        The rank of the Sylvester matrix is related to the degree of the
        gcd of self and right::

            sage: f.gcd(g).degree() == f.degree() + g.degree() - M.rank()
            True
            sage: h1.gcd(h2).degree() == h1.degree() + h2.degree() - M1.rank()
            True

        TESTS:

        The variable is optional, but must be the same in both rings::

            sage: K.<x> = QQ['x']
            sage: f = x+1
            sage: g = QQ['y']([1, 0, 1])
            sage: f.sylvester_matrix(f, x)
            [1 1]
            [1 1]
            sage: f.sylvester_matrix(g, x)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

        Polynomials must be defined over compatible base rings::

            sage: f = QQ['x']([1, 0, 1])
            sage: g = ZZ['x']([1, 0, 1])
            sage: h = GF(25, 'a')['x']([1, 0, 1])
            sage: f.sylvester_matrix(g)
            [1 0 1 0]
            [0 1 0 1]
            [1 0 1 0]
            [0 1 0 1]
            sage: g.sylvester_matrix(h)
            [1 0 1 0]
            [0 1 0 1]
            [1 0 1 0]
            [0 1 0 1]
            sage: f.sylvester_matrix(h)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in x over Finite Field in a of size 5^2'

        We can compute the sylvester matrix of a univariate and multivariate
        polynomial::

            sage: K.<x,y> = QQ['x,y']
            sage: g = K.random_element()
            sage: f.sylvester_matrix(g) == K(f).sylvester_matrix(g,x)
            True

        Corner cases::

            sage: K.<x>=QQ[]
            sage: f = x^2+1
            sage: g = K(0)
            sage: f.sylvester_matrix(g)
            Traceback (most recent call last):
            ...
            ValueError: The Sylvester matrix is not defined for zero polynomials
            sage: g.sylvester_matrix(f)
            Traceback (most recent call last):
            ...
            ValueError: The Sylvester matrix is not defined for zero polynomials
            sage: g.sylvester_matrix(g)
            Traceback (most recent call last):
            ...
            ValueError: The Sylvester matrix is not defined for zero polynomials
            sage: K(3).sylvester_matrix(x^2)
            [3 0]
            [0 3]
            sage: K(3).sylvester_matrix(K(4))
            []
        """

        # This code is almost exactly the same as that of
        # sylvester_matrix() in multi_polynomial.pyx.

        if self.parent() != right.parent():
            coercion_model = sage.structure.element.get_coercion_model()
            a, b = coercion_model.canonical_coercion(self,right)
            variable = a.parent()(self.variables()[0])
            #We add the variable to cover the case that right is a multivariate
            #polynomial
            return a.sylvester_matrix(b, variable)

        if variable:
            if variable.parent() != self.parent():
                variable = self.parent()(variable)

        from sage.matrix.constructor import matrix

        # The dimension of the sage matrix is self.degree() + right.degree()

        if self.is_zero() or right.is_zero():
            raise ValueError("The Sylvester matrix is not defined for zero polynomials")

        m = self.degree()
        n = right.degree()

        M = matrix(self.base_ring(), n + m, n + m)

        r = 0
        offset = 0
        for _ in range(n):
            for c in range(m, -1, -1):
                M[r, m - c + offset] = self[c]
            offset += 1
            r += 1

        offset = 0
        for _ in range(m):
            for c in range(n, -1, -1):
                M[r, n - c + offset] = right[c]
            offset += 1
            r += 1

        return M

    cpdef constant_coefficient(self):
        """
        Return the constant coefficient of this polynomial.

        OUTPUT: element of base ring

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = -2*x^3 + 2*x - 1/3
            sage: f.constant_coefficient()
            -1/3
        """
        return self[0]

    cpdef Polynomial _new_constant_poly(self, a, Parent P):
        """
        Create a new constant polynomial from a in P, which MUST be an
        element of the base ring of P (this is not checked).

        EXAMPLE::

            sage: R.<w> = PolynomialRing(GF(9,'a'), sparse=True)
            sage: a = w._new_constant_poly(0, R); a
            0
            sage: a.coeffs()
            []

        """
        if a:
            return self.__class__(P,[a], check=False) #P._element_constructor(a, check=False)
        return self.__class__(P,[], check=False)

    def is_monic(self):
        """
        Returns True if this polynomial is monic. The zero polynomial is by
        definition not monic.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = x + 33
            sage: f.is_monic()
            True
            sage: f = 0*x
            sage: f.is_monic()
            False
            sage: f = 3*x^3 + x^4 + x^2
            sage: f.is_monic()
            True
            sage: f = 2*x^2 + x^3 + 56*x^5
            sage: f.is_monic()
            False

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        """
        return not self.is_zero() and self[self.degree()] == 1

    def is_unit(self):
        r"""
        Return True if this polynomial is a unit.

        EXAMPLES::

            sage: a = Integers(90384098234^3)
            sage: b = a(2*191*236607587)
            sage: b.is_nilpotent()
            True
            sage: R.<x> = a[]
            sage: f = 3 + b*x + b^2*x^2
            sage: f.is_unit()
            True
            sage: f = 3 + b*x + b^2*x^2 + 17*x^3
            sage: f.is_unit()
            False

        EXERCISE (Atiyah-McDonald, Ch 1): Let `A[x]` be a
        polynomial ring in one variable. Then
        `f=\sum a_i x^i \in A[x]` is a unit if and only if
        `a_0` is a unit and `a_1,\ldots, a_n` are
        nilpotent.
        """
        if self.degree() > 0:
            for i in range(1,self.degree()+1):
                if not self[i].is_nilpotent():
                    return False
        return self[0].is_unit()

    def is_nilpotent(self):
        r"""
        Return True if this polynomial is nilpotent.

        EXAMPLES::

            sage: R = Integers(12)
            sage: S.<x> = R[]
            sage: f = 5 + 6*x
            sage: f.is_nilpotent()
            False
            sage: f = 6 + 6*x^2
            sage: f.is_nilpotent()
            True
            sage: f^2
            0

        EXERCISE (Atiyah-McDonald, Ch 1): Let `A[x]` be a
        polynomial ring in one variable. Then
        `f=\sum a_i x^i \in A[x]` is nilpotent if and only if
        every `a_i` is nilpotent.
        """
        for i in range(self.degree()+1):
            if not self[i].is_nilpotent():
                return False
        return True

    def is_gen(self):
        r"""
        Return True if this polynomial is the distinguished generator of
        the parent polynomial ring.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R(1).is_gen()
            False
            sage: R(x).is_gen()
            True

        Important - this function doesn't return True if self equals the
        generator; it returns True if self *is* the generator.

        ::

            sage: f = R([0,1]); f
            x
            sage: f.is_gen()
            False
            sage: f is x
            False
            sage: f == x
            True
        """
        return bool(self._is_gen)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        OUTPUT: element of the base ring

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (-2/5)*x^3 + 2*x - 1/3
            sage: f.leading_coefficient()
            -2/5
        """
        return self[self.degree()]

    def monic(self):
        """
        Return this polynomial divided by its leading coefficient. Does not
        change this polynomial.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = 2*x^2 + x^3 + 56*x^5
            sage: f.monic()
            x^5 + 1/56*x^3 + 1/28*x^2
            sage: f = (1/4)*x^2 + 3*x + 1
            sage: f.monic()
            x^2 + 12*x + 4

        The following happens because `f = 0` cannot be made into a
        monic polynomial

        ::

            sage: f = 0*x
            sage: f.monic()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        Notice that the monic version of a polynomial over the integers is
        defined over the rationals.

        ::

            sage: x = ZZ['x'].0
            sage: f = 3*x^19 + x^2 - 37
            sage: g = f.monic(); g
            x^19 + 1/3*x^2 - 37/3
            sage: g.parent()
            Univariate Polynomial Ring in x over Rational Field

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        """
        if self.is_monic():
            return self
        a = ~self.leading_coefficient()
        R = self.parent()
        if a.parent() != R.base_ring():
            S = R.base_extend(a.parent())
            return a*S(self)
        else:
            return a*self

    def coefficients(self):
        """
        Return the coefficients of the monomials appearing in self.

        EXAMPLES::

            sage: _.<x> = PolynomialRing(ZZ)
            sage: f = x^4+2*x^2+1
            sage: f.coefficients()
            [1, 2, 1]
        """
        zero = self.parent().base_ring().zero_element()
        return [c for c in self.list() if c != zero]

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES::

            sage: _.<x> = PolynomialRing(ZZ)
            sage: f = x^4+2*x^2+1
            sage: f.exponents()
            [0, 2, 4]
        """
        zero = self.parent().base_ring().zero_element()
        l = self.list()
        return [i for i in range(len(l)) if l[i] != zero]

    def list(self):
        """
        Return a new copy of the list of the underlying elements of self.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (-2/5)*x^3 + 2*x - 1/3
            sage: v = f.list(); v
            [-1/3, 2, 0, -2/5]

        Note that v is a list, it is mutable, and each call to the list
        method returns a new list::

            sage: type(v)
            <type 'list'>
            sage: v[0] = 5
            sage: f.list()
            [-1/3, 2, 0, -2/5]

        Here is an example with a generic polynomial ring::

            sage: R.<x> = QQ[]
            sage: S.<y> = R[]
            sage: f = y^3 + x*y -3*x; f
            y^3 + x*y - 3*x
            sage: type(f)
            <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
            sage: v = f.list(); v
            [-3*x, x, 0, 1]
            sage: v[0] = 10
            sage: f.list()
            [-3*x, x, 0, 1]
        """
        raise NotImplementedError

    def prec(self):
        """
        Return the precision of this polynomial. This is always infinity,
        since polynomials are of infinite precision by definition (there is
        no big-oh).

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: (x^5 + x + 1).prec()
            +Infinity
            sage: x.prec()
            +Infinity
        """
        return infinity.infinity

    def padded_list(self, n=None):
        """
        Return list of coefficients of self up to (but not including)
        `q^n`.

        Includes 0's in the list on the right so that the list has length
        `n`.

        INPUT:


        -  ``n`` - (default: None); if given, an integer that
           is at least 0


        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = 1 + x^3 + 23*x^5
            sage: f.padded_list()
            [1, 0, 0, 1, 0, 23]
            sage: f.padded_list(10)
            [1, 0, 0, 1, 0, 23, 0, 0, 0, 0]
            sage: len(f.padded_list(10))
            10
            sage: f.padded_list(3)
            [1, 0, 0]
            sage: f.padded_list(0)
            []
            sage: f.padded_list(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be at least 0
        """
        v = self.list()
        if n is None:
            return v
        if n < 0:
            raise ValueError("n must be at least 0")
        if len(v) < n:
            z = self._parent.base_ring().zero_element()
            return v + [z]*(n - len(v))
        else:
            return v[:int(n)]

    def coeffs(self):
        r"""
        Returns ``self.list()``.

        (It is potentially slightly faster to use
        ``self.list()`` directly.)

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = 10*x^3 + 5*x + 2/17
            sage: f.coeffs()
            [2/17, 5, 0, 10]
        """
        return self.list()

    def newton_raphson(self, n, x0):
        """
        Return a list of n iterative approximations to a root of this
        polynomial, computed using the Newton-Raphson method.

        The Newton-Raphson method is an iterative root-finding algorithm.
        For f(x) a polynomial, as is the case here, this is essentially the
        same as Horner's method.

        INPUT:


        -  ``n`` - an integer (=the number of iterations),

        -  ``x0`` - an initial guess x0.


        OUTPUT: A list of numbers hopefully approximating a root of
        f(x)=0.

        If one of the iterates is a critical point of f then a
        ZeroDivisionError exception is raised.

        EXAMPLES::

            sage: x = PolynomialRing(RealField(), 'x').gen()
            sage: f = x^2 - 2
            sage: f.newton_raphson(4, 1)
            [1.50000000000000, 1.41666666666667, 1.41421568627451, 1.41421356237469]

        AUTHORS:

        - David Joyner and William Stein (2005-11-28)
        """
        n = sage.rings.integer.Integer(n)
        df = self.derivative()
        K = self.parent().base_ring()
        a = K(x0)
        L = []
        for i in range(n):
            a -= self(a) / df(a)
            L.append(a)
        return L

    def polynomial(self, var):
        r"""
        Let var be one of the variables of the parent of self. This returns
        self viewed as a univariate polynomial in var over the polynomial
        ring generated by all the other variables of the parent.

        For univariate polynomials, if var is the generator of the parent
        ring, we return this polynomial, otherwise raise an error.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (x+1).polynomial(x)
            x + 1

        TESTS::

            sage: x.polynomial(1)
            Traceback (most recent call last):
            ...
            ValueError: given variable is not the generator of parent.
        """
        if self._parent.ngens() == 1:
            if self._parent.gen() == var:
                return self
            raise ValueError("given variable is not the generator of parent.")
        raise NotImplementedError

    def newton_slopes(self, p):
        """
        Return the `p`-adic slopes of the Newton polygon of self,
        when this makes sense.

        OUTPUT: list of rational numbers

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = x^3 + 2
            sage: f.newton_slopes(2)
            [1/3, 1/3, 1/3]

        ALGORITHM: Uses PARI.
        """
        f = self._pari_()
        v = list(f.newtonpoly(p))
        return [sage.rings.rational.Rational(x) for x in v]

    #####################################################################
    # Conversions to other systems
    #####################################################################
    def _pari_(self):
        r"""
        Return polynomial as a PARI object.

        Sage does not handle PARI's variable ordering requirements
        gracefully at this time. In practice, this means that the variable
        ``x`` needs to be the topmost variable, as in the
        example.

        EXAMPLES::

            sage: f = QQ['x']([0,1,2/3,3])
            sage: pari(f)
            3*x^3 + 2/3*x^2 + x

        ::

            sage: S.<a> = QQ['a']
            sage: R.<x> = S['x']
            sage: f = R([0, a]) + R([0, 0, 2/3])
            sage: pari(f)
            2/3*x^2 + a*x

        Polynomials over a number field work, provided that the variable is
        called 'x'::

            sage: x = polygen(QQ)
            sage: K.<b> = NumberField(x^2 + x + 1)
            sage: R.<x> = PolynomialRing(K)
            sage: pol = (b + x)^3; pol
            x^3 + 3*b*x^2 + (-3*b - 3)*x + 1
            sage: pari(pol)
            Mod(1, y^2 + y + 1)*x^3 + Mod(3*y, y^2 + y + 1)*x^2 + Mod(-3*y - 3, y^2 + y + 1)*x + Mod(1, y^2 + y + 1)

        TESTS:

        Unfortunately, variable names matter::

            sage: R.<x, y> = QQ[]
            sage: S.<a> = R[]
            sage: f = x^2 + a; g = y^3 + a
            sage: pari(f)
            Traceback (most recent call last):
            ...
            PariError: variable must have higher priority in gtopoly

        Stacked polynomial rings, first with a univariate ring on the
        bottom::

            sage: S.<a> = QQ['a']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 2*x)
            x^2 + 2*x
            sage: pari(a*x + 2*x^3)
            2*x^3 + a*x

        Stacked polynomial rings, second with a multivariate ring on the
        bottom::

            sage: S.<a, b> = ZZ['a', 'b']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 2*x)
            x^2 + 2*x
            sage: pari(a*x + 2*b*x^3)
            2*b*x^3 + a*x

        Stacked polynomial rings with exotic base rings::

            sage: S.<a, b> = GF(7)['a', 'b']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 9*x)
            x^2 + 2*x
            sage: pari(a*x + 9*b*x^3)
            2*b*x^3 + a*x

        ::

            sage: S.<a> = Integers(8)['a']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 2*x)
            Mod(1, 8)*x^2 + Mod(2, 8)*x
            sage: pari(a*x + 10*x^3)
            Mod(2, 8)*x^3 + Mod(1, 8)*a*x
        """
        return self._pari_with_name(self.parent().variable_name())

    def _pari_or_constant(self, name=None):
        r"""
        Convert ``self`` to PARI.  This behaves identical to :meth:`_pari_`
        or :meth:`_pari_with_name` except for constant polynomials:
        then the constant is returned instead of a constant polynomial.

        INPUT:

        - ``name`` -- (default: None) Variable name.  If not given, use
          ``self.parent().variable_name()``.  This argument is irrelevant
          for constant polynomials.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: pol = 2*x^2 + 7*x - 5
            sage: pol._pari_or_constant()
            2*x^2 + 7*x - 5
            sage: pol._pari_or_constant('a')
            2*a^2 + 7*a - 5
            sage: pol = R(7)
            sage: pol._pari_or_constant()
            7
            sage: pol._pari_or_constant().type()
            't_INT'
            sage: pol._pari_().type()
            't_POL'
            sage: PolynomialRing(IntegerModRing(101), 't')()._pari_or_constant()
            Mod(0, 101)
        """
        if self.is_constant():
            return self[0]._pari_()
        if name is None:
            name = self.parent().variable_name()
        return self._pari_with_name(name)

    def _pari_with_name(self, name='x'):
        r"""
        Return polynomial as a PARI object with topmost variable
        ``name``.  By default, use 'x' for the variable name.

        For internal use only.

        EXAMPLES:

            sage: R.<a> = PolynomialRing(ZZ)
            sage: (2*a^2 + a)._pari_with_name()
            2*x^2 + x
            sage: (2*a^2 + a)._pari_with_name('y')
            2*y^2 + y
        """
        vals = [x._pari_() for x in self.list()]
        return pari(vals).Polrev(name)

    def _pari_init_(self):
        return repr(self._pari_())

    def _magma_init_(self, magma):
        """
        Return a string that evaluates in Magma to this polynomial.

        EXAMPLES::

            sage: magma = Magma()  # new session
            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: f._magma_init_(magma)        # optional - magma
            '_sage_[...]![5,-17,0,1]'
            sage: g = magma(f); g              # optional - magma
            y^3 - 17*y + 5

        Note that in Magma there is only one polynomial ring over each
        base, so if we make the polynomial ring over ZZ with variable
        `z`, then this changes the variable name of the polynomial
        we already defined::

            sage: R.<z> = ZZ[]
            sage: magma(R)                     # optional - magma
            Univariate Polynomial Ring in z over Integer Ring
            sage: g                            # optional - magma
            z^3 - 17*z + 5

        In Sage the variable name does not change::

            sage: f
            y^3 - 17*y + 5

        A more complicated nested example::

            sage: k.<a> = GF(9); R.<s,t> = k[]; S.<W> = R[]
            sage: magma(a*W^20 + s*t/a)        # optional - magma
            a*W^20 + a^7*s*t
        """
        # Get a reference to Magma version of parent.
        R = magma(self.parent())
        # Get list of coefficients.
        v = ','.join([a._magma_init_(magma) for a in self.list()])
        return '%s![%s]'%(R.name(), v)

    def _gap_init_(self):
        return repr(self)

    def _gap_(self, G):
        """
        EXAMPLES::

            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: g = gap(f); g
            y^3-17*y+5
            sage: f._gap_init_()
            'y^3 - 17*y + 5'
            sage: R.<z> = ZZ[]
            sage: gap(R)
            PolynomialRing( Integers, ["z"] )
            sage: g
            y^3-17*y+5
            sage: gap(z^2 + z)
            z^2+z

        We coerce a polynomial with coefficients in a finite field::

            sage: R.<y> = GF(7)[]
            sage: f = y^3 - 17*y + 5
            sage: g = gap(f); g
            y^3+Z(7)^4*y+Z(7)^5
            sage: g.Factors()
            [ y+Z(7)^0, y+Z(7)^0, y+Z(7)^5 ]
            sage: f.factor()
            (y + 5) * (y + 1)^2
        """
        if G is None:
            import sage.interfaces.gap
            G = sage.interfaces.gap.gap
        self.parent()._gap_(G)
        return G(self._gap_init_())

    ######################################################################

    @coerce_binop
    def resultant(self, other):
        r"""
        Returns the resultant of self and other.

        INPUT:


        -  ``other`` - a polynomial


        OUTPUT: an element of the base ring of the polynomial ring

        .. note::

           Implemented using PARI's ``polresultant`` function.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            -8
            sage: r.parent() is QQ
            True

        We can also compute resultants over univariate and multivariate
        polynomial rings, provided that PARI's variable ordering
        requirements are respected. Usually, your resultants will work if
        you always ask for them in the variable ``x``::

            sage: R.<a> = QQ[]
            sage: S.<x> = R[]
            sage: f = x^2 + a; g = x^3 + a
            sage: r = f.resultant(g); r
            a^3 + a^2
            sage: r.parent() is R
            True

        ::

            sage: R.<a, b> = QQ[]
            sage: S.<x> = R[]
            sage: f = x^2 + a; g = x^3 + b
            sage: r = f.resultant(g); r
            a^3 + b^2
            sage: r.parent() is R
            True

        Unfortunately Sage does not handle PARI's variable ordering
        requirements gracefully, so the following has to be done
        through Singular::

            sage: R.<x, y> = QQ[]
            sage: S.<a> = R[]
            sage: f = x^2 + a; g = y^3 + a
            sage: h = f.resultant(g); h
            y^3 - x^2
            sage: h.parent() is R
            True

        Check that :trac:`13672` is fixed::

            sage: R.<t> = GF(2)[]
            sage: S.<x> = R[]
            sage: f = (t^2 + t)*x + t^2 + t
            sage: g = (t + 1)*x + t^2
            sage: f.resultant(g)
            t^4 + t
        """
        variable = self.variable_name()
        if variable != 'x' and self.parent()._mpoly_base_ring() != self.parent().base_ring():
            bigring = sage.rings.polynomial.multi_polynomial.PolynomialRing(self.parent()._mpoly_base_ring(),list(self.parent().variable_names_recursive()))
            newself = bigring(self)
            newother = bigring(other)
            return self.parent().base_ring()(newself.resultant(newother,bigring(self.parent().gen())))
        # Single-variable polynomial or main variable is "x": we can use PARI to compute the resultant
        res = self._pari_with_name().polresultant(other._pari_with_name())
        return self.parent().base_ring()(res)

    def discriminant(self):
        r"""
        Returns the discriminant of self.

        The discriminant is

        .. math::

            R_n := a_n^{2 n-2} \prod_{1<i<j<n} (r_i-r_j)^2,


        where `n` is the degree of self, `a_n` is the
        leading coefficient of self and the roots of self are
        `r_1, \ldots, r_n`.

        OUTPUT: An element of the base ring of the polynomial ring.

        .. note::

           Uses the identity `R_n(f) := (-1)^(n (n-1)/2) R(f,
           f') a_n^(n-k-2)`, where `n` is the degree of self,
           `a_n` is the leading coefficient of self, `f'`
           is the derivative of `f`, and `k` is the degree
           of `f'`. Calls :meth:`.resultant`.

        EXAMPLES:

        In the case of elliptic curves in special form, the discriminant is
        easy to calculate::

            sage: R.<x> = QQ[]
            sage: f = x^3 + x + 1
            sage: d = f.discriminant(); d
            -31
            sage: d.parent() is QQ
            True
            sage: EllipticCurve([1, 1]).discriminant()/16
            -31

        ::

            sage: R.<x> = QQ[]
            sage: f = 2*x^3 + x + 1
            sage: d = f.discriminant(); d
            -116

        We can also compute discriminants over univariate and multivariate
        polynomial rings, provided that PARI's variable ordering
        requirements are respected. Usually, your discriminants will work
        if you always ask for them in the variable ``x``::

            sage: R.<a> = QQ[]
            sage: S.<x> = R[]
            sage: f = a*x + x + a + 1
            sage: d = f.discriminant(); d
            1
            sage: d.parent() is R
            True

        ::

            sage: R.<a, b> = QQ[]
            sage: S.<x> = R[]
            sage: f = x^2 + a + b
            sage: d = f.discriminant(); d
            -4*a - 4*b
            sage: d.parent() is R
            True

        Unfortunately Sage does not handle PARI's variable ordering
        requirements gracefully, so the following has to be done
        through Singular::

            sage: R.<x, y> = QQ[]
            sage: S.<a> = R[]
            sage: f = x^2 + a
            sage: f.discriminant()
            1

        Check that :trac:`13672` is fixed::

            sage: R.<t> = GF(5)[]
            sage: S.<x> = R[]
            sage: f = x^10 + 2*x^6 + 2*x^5 + x + 2
            sage: (f-t).discriminant()
            4*t^5

        The following examples show that :trac:`11782` has been fixed::

            sage: ZZ.quo(81)[x](3*x^2 + 3*x + 3).discriminant()
            54
            sage: ZZ.quo(9)[x](2*x^3 + x^2 + x).discriminant()
            2

        This was fixed by :trac:`15422`::

            sage: R.<s> = PolynomialRing(Qp(2))
            sage: (s^2).discriminant()
            0
        """
        if self.is_zero():
            return self.parent().zero_element()
        n = self.degree()
        d = self.derivative()
        k = d.degree()

        r = n % 4
        u = -1 # (-1)**(n*(n-1)/2)
        if r == 0 or r == 1:
            u = 1
        try:
            an = self[n]**(n - k - 2)
        except ZeroDivisionError:
            assert(n-k-2 == -1)
            # Rather than dividing the resultant by the leading coefficient,
            # we alter the Sylvester matrix (see #11782).
            mat = self.sylvester_matrix(d)
            mat[0, 0] = self.base_ring()(1)
            mat[n - 1, 0] = self.base_ring()(n)
            return u * mat.determinant()
        else:
            return self.base_ring()(u * self.resultant(d) * an)

    def reverse(self, degree=None):
        """
        Return polynomial but with the coefficients reversed.

        If an optional degree argument is given the coefficient list will be
        truncated or zero padded as necessary and the reverse polynomial will
        have the specified degree.

        EXAMPLES::

            sage: R.<x> = ZZ[]; S.<y> = R[]
            sage: f = y^3 + x*y -3*x; f
            y^3 + x*y - 3*x
            sage: f.reverse()
            -3*x*y^3 + x*y^2 + 1
            sage: f.reverse(degree=2)
            -3*x*y^2 + x*y
            sage: f.reverse(degree=5)
            -3*x*y^5 + x*y^4 + y^2

        TESTS::

            sage: f.reverse(degree=1.5r)
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be a non-negative integer, got 1.5
        """
        v = list(self.list())

        cdef unsigned long d
        if degree:
            d = degree
            if d != degree:
                raise ValueError("degree argument must be a non-negative integer, got %s"%(degree))
            if len(v) < degree+1:
                v.reverse()
                v = [0]*(degree+1-len(v)) + v
            elif len(v) > degree+1:
                v = v[:degree+1]
                v.reverse()
            else: # len(v) == degree + 1
                v.reverse()
        else:
            v.reverse()

        return self.parent()(v)

    def roots(self, ring=None, multiplicities=True, algorithm=None):
        """
        Return the roots of this polynomial (by default, in the base ring
        of this polynomial).

        INPUT:


        -  ``ring`` - the ring to find roots in

        -  ``multiplicities`` - bool (default: True) if True
           return list of pairs (r, n), where r is the root and n is the
           multiplicity. If False, just return the unique roots, with no
           information about multiplicities.

        -  ``algorithm`` - the root-finding algorithm to use.
           We attempt to select a reasonable algorithm by default, but this
           lets the caller override our choice.


        By default, this finds all the roots that lie in the base ring of
        the polynomial. However, the ring parameter can be used to specify
        a ring to look for roots in.

        If the polynomial and the output ring are both exact (integers,
        rationals, finite fields, etc.), then the output should always be
        correct (or raise an exception, if that case is not yet handled).

        If the output ring is approximate (floating-point real or complex
        numbers), then the answer will be estimated numerically, using
        floating-point arithmetic of at least the precision of the output
        ring. If the polynomial is ill-conditioned, meaning that a small
        change in the coefficients of the polynomial will lead to a
        relatively large change in the location of the roots, this may give
        poor results. Distinct roots may be returned as multiple roots,
        multiple roots may be returned as distinct roots, real roots may be
        lost entirely (because the numerical estimate thinks they are
        complex roots). Note that polynomials with multiple roots are
        always ill-conditioned; there's a footnote at the end of the
        docstring about this.

        If the output ring is a RealIntervalField or ComplexIntervalField
        of a given precision, then the answer will always be correct (or an
        exception will be raised, if a case is not implemented). Each root
        will be contained in one of the returned intervals, and the
        intervals will be disjoint. (The returned intervals may be of
        higher precision than the specified output ring.)

        At the end of this docstring (after the examples) is a description
        of all the cases implemented in this function, and the algorithms
        used. That section also describes the possibilities for
        "algorithm=", for the cases where multiple algorithms exist.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = x^3 - 1
            sage: f.roots()
            [(1, 1)]
            sage: f.roots(ring=CC)   # note -- low order bits slightly different on ppc.
            [(1.00000000000000, 1), (-0.500000000000000 - 0.86602540378443...*I, 1), (-0.500000000000000 + 0.86602540378443...*I, 1)]
            sage: f = (x^3 - 1)^2
            sage: f.roots()
            [(1, 2)]

        ::

            sage: f = -19*x + 884736
            sage: f.roots()
            [(884736/19, 1)]
            sage: (f^20).roots()
            [(884736/19, 20)]

        ::

            sage: K.<z> = CyclotomicField(3)
            sage: f = K.defining_polynomial()
            sage: f.roots(ring=GF(7))
            [(4, 1), (2, 1)]
            sage: g = f.change_ring(GF(7))
            sage: g.roots()
            [(4, 1), (2, 1)]
            sage: g.roots(multiplicities=False)
            [4, 2]

        A new ring.  In the example below, we add the special method
        _roots_univariate_polynomial to the base ring, and observe
        that this method is called instead to find roots of
        polynomials over this ring.  This facility can be used to
        easily extend root finding to work over new rings you
        introduce::

             sage: R.<x> = QQ[]
             sage: (x^2 + 1).roots()
             []
             sage: g = lambda f, *args, **kwds: f.change_ring(CDF).roots()
             sage: QQ._roots_univariate_polynomial = g
             sage: (x^2 + 1).roots()
             [(... - 1.0*I, 1), (1.0*I, 1)]
             sage: del QQ._roots_univariate_polynomial

        An example over RR, which illustrates that only the roots in RR are
        returned::

            sage: x = RR['x'].0
            sage: f = x^3 -2
            sage: f.roots()
            [(1.25992104989487, 1)]
            sage: f.factor()
            (x - 1.25992104989487) * (x^2 + 1.25992104989487*x + 1.58740105196820)
            sage: x = RealField(100)['x'].0
            sage: f = x^3 -2
            sage: f.roots()
            [(1.2599210498948731647672106073, 1)]

        ::

            sage: x = CC['x'].0
            sage: f = x^3 -2
            sage: f.roots()
            [(1.25992104989487, 1), (-0.62996052494743... - 1.09112363597172*I, 1), (-0.62996052494743... + 1.09112363597172*I, 1)]
            sage: f.roots(algorithm='pari')
            [(1.25992104989487, 1), (-0.629960524947437 - 1.09112363597172*I, 1), (-0.629960524947437 + 1.09112363597172*I, 1)]

        Another example showing that only roots in the base ring are
        returned::

            sage: x = polygen(ZZ)
            sage: f = (2*x-3) * (x-1) * (x+1)
            sage: f.roots()
            [(1, 1), (-1, 1)]
            sage: f.roots(ring=QQ)
            [(3/2, 1), (1, 1), (-1, 1)]

        An example involving large numbers::

            sage: x = RR['x'].0
            sage: f = x^2 - 1e100
            sage: f.roots()
            [(-1.00000000000000e50, 1), (1.00000000000000e50, 1)]
            sage: f = x^10 - 2*(5*x-1)^2
            sage: f.roots(multiplicities=False)
            [-1.6772670339941..., 0.19995479628..., 0.20004530611..., 1.5763035161844...]

        ::

            sage: x = CC['x'].0
            sage: i = CC.0
            sage: f = (x - 1)*(x - i)
            sage: f.roots(multiplicities=False) #random - this example is numerically rather unstable
            [2.22044604925031e-16 + 1.00000000000000*I, 1.00000000000000 + 8.32667268468867e-17*I]
            sage: g=(x-1.33+1.33*i)*(x-2.66-2.66*i)
            sage: g.roots(multiplicities=False)
            [1.33000000000000 - 1.33000000000000*I, 2.66000000000000 + 2.66000000000000*I]

        A purely symbolic roots example::

            sage: X = var('X')
            sage: f = expand((X-1)*(X-I)^3*(X^2 - sqrt(2))); f
            X^6 - (3*I + 1)*X^5 - sqrt(2)*X^4 + (3*I - 3)*X^4 + (3*I + 1)*sqrt(2)*X^3 + (I + 3)*X^3 - (3*I - 3)*sqrt(2)*X^2 - I*X^2 - (I + 3)*sqrt(2)*X + I*sqrt(2)
            sage: print f.roots()
            [(I, 3), (-2^(1/4), 1), (2^(1/4), 1), (1, 1)]

        A couple of examples where the base ring doesn't have a
        factorization algorithm (yet). Note that this is currently done via
        naive enumeration, so could be very slow::

            sage: R = Integers(6)
            sage: S.<x> = R['x']
            sage: p = x^2-1
            sage: p.roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)
            sage: p.roots(multiplicities=False)
            [1, 5]
            sage: R = Integers(9)
            sage: A = PolynomialRing(R, 'y')
            sage: y = A.gen()
            sage: f = 10*y^2 - y^3 - 9
            sage: f.roots(multiplicities=False)
            [0, 1, 3, 6]

        An example over the complex double field (where root finding is
        fast, thanks to NumPy)::

            sage: R.<x> = CDF[]
            sage: f = R.cyclotomic_polynomial(5); f
            x^4 + x^3 + x^2 + x + 1.0
            sage: f.roots(multiplicities=False)   # slightly random
            [0.309016994375 + 0.951056516295*I, 0.309016994375 - 0.951056516295*I, -0.809016994375 + 0.587785252292*I, -0.809016994375 - 0.587785252292*I]
            sage: [z^5 for z in f.roots(multiplicities=False)]     # slightly random
            [1.0 - 2.44929359829e-16*I, 1.0 + 2.44929359829e-16*I, 1.0 - 4.89858719659e-16*I, 1.0 + 4.89858719659e-16*I]
            sage: f = CDF['x']([1,2,3,4]); f
            4.0*x^3 + 3.0*x^2 + 2.0*x + 1.0
            sage: r = f.roots(multiplicities=False)
            sage: [f(a) for a in r]    # slightly random
            [2.55351295664e-15, -4.4408920985e-16 - 2.08166817117e-16*I, -4.4408920985e-16 + 2.08166817117e-16*I]

        Another example over RDF::

            sage: x = RDF['x'].0
            sage: ((x^3 -1)).roots()
            [(1.0, 1)]
            sage: ((x^3 -1)).roots(multiplicities=False)
            [1.0]

        More examples involving the complex double field::

            sage: x = CDF['x'].0
            sage: i = CDF.0
            sage: f = x^3 + 2*i; f
            x^3 + 2.0*I
            sage: f.roots()
            [(-1.09112363597 - 0.629960524947*I, 1),
             (... + 1.25992104989*I, 1),
             (1.09112363597 - 0.629960524947*I, 1)]
            sage: f.roots(multiplicities=False)
            [-1.09112363597 - 0.629960524947*I, ... + 1.25992104989*I, 1.09112363597 - 0.629960524947*I]
            sage: [f(z) for z in f.roots(multiplicities=False)] # random, too close to 0
            [-2.56337823492e-15 - 6.66133814775e-15*I,
             3.96533069372e-16 + 1.99840144433e-15*I,
             4.19092485179e-17 - 8.881784197e-16*I]
            sage: f = i*x^3 + 2; f
            I*x^3 + 2.0
            sage: f.roots()
            [(-1.09112363597 + 0.629960524947*I, 1),
             (... - 1.25992104989*I, 1),
             (1.09112363597 + 0.629960524947*I, 1)]
            sage: f(f.roots()[0][0]) # random, too close to 0
            -2.56337823492e-15 - 6.66133814775e-15*I

        Examples using real root isolation::

            sage: x = polygen(ZZ)
            sage: f = x^2 - x - 1
            sage: f.roots()
            []
            sage: f.roots(ring=RIF)
            [(-0.6180339887498948482045868343657?, 1), (1.6180339887498948482045868343657?, 1)]
            sage: f.roots(ring=RIF, multiplicities=False)
            [-0.6180339887498948482045868343657?, 1.6180339887498948482045868343657?]
            sage: f.roots(ring=RealIntervalField(150))
            [(-0.6180339887498948482045868343656381177203091798057628621354486227?, 1), (1.618033988749894848204586834365638117720309179805762862135448623?, 1)]
            sage: f.roots(ring=AA)
            [(-0.618033988749895?, 1), (1.618033988749895?, 1)]
            sage: f = f^2 * (x - 1)
            sage: f.roots(ring=RIF)
            [(-0.6180339887498948482045868343657?, 2), (1.0000000000000000000000000000000?, 1), (1.6180339887498948482045868343657?, 2)]
            sage: f.roots(ring=RIF, multiplicities=False)
            [-0.6180339887498948482045868343657?, 1.0000000000000000000000000000000?, 1.6180339887498948482045868343657?]

        Examples using complex root isolation::

            sage: x = polygen(ZZ)
            sage: p = x^5 - x - 1
            sage: p.roots()
            []
            sage: p.roots(ring=CIF)
            [(1.167303978261419?, 1), (-0.764884433600585? - 0.352471546031727?*I, 1), (-0.764884433600585? + 0.352471546031727?*I, 1), (0.181232444469876? - 1.083954101317711?*I, 1), (0.181232444469876? + 1.083954101317711?*I, 1)]
            sage: p.roots(ring=ComplexIntervalField(200))
            [(1.167303978261418684256045899854842180720560371525489039140082?, 1), (-0.76488443360058472602982318770854173032899665194736756700778? - 0.35247154603172624931794709140258105439420648082424733283770?*I, 1), (-0.76488443360058472602982318770854173032899665194736756700778? + 0.35247154603172624931794709140258105439420648082424733283770?*I, 1), (0.18123244446987538390180023778112063996871646618462304743774? - 1.08395410131771066843034449298076657427364024315511565430114?*I, 1), (0.18123244446987538390180023778112063996871646618462304743774? + 1.08395410131771066843034449298076657427364024315511565430114?*I, 1)]
            sage: rts = p.roots(ring=QQbar); rts
            [(1.167303978261419?, 1), (-0.7648844336005847? - 0.3524715460317263?*I, 1), (-0.7648844336005847? + 0.3524715460317263?*I, 1), (0.1812324444698754? - 1.083954101317711?*I, 1), (0.1812324444698754? + 1.083954101317711?*I, 1)]
            sage: p.roots(ring=AA)
            [(1.167303978261419?, 1)]
            sage: p = (x - rts[4][0])^2 * (3*x^2 + x + 1)
            sage: p.roots(ring=QQbar)
            [(-0.1666666666666667? - 0.552770798392567?*I, 1), (-0.1666666666666667? + 0.552770798392567?*I, 1), (0.1812324444698754? + 1.083954101317711?*I, 2)]
            sage: p.roots(ring=CIF)
            [(-0.1666666666666667? - 0.552770798392567?*I, 1), (-0.1666666666666667? + 0.552770798392567?*I, 1), (0.1812324444698754? + 1.083954101317711?*I, 2)]

        Note that coefficients in a number field with defining polynomial
        `x^2 + 1` are considered to be Gaussian rationals (with the
        generator mapping to +I), if you ask for complex roots.

        ::

            sage: K.<im> = NumberField(x^2 + 1)
            sage: y = polygen(K)
            sage: p = y^4 - 2 - im
            sage: p.roots(ring=CC)
            [(-1.2146389322441... - 0.14142505258239...*I, 1), (-0.14142505258239... + 1.2146389322441...*I, 1), (0.14142505258239... - 1.2146389322441...*I, 1), (1.2146389322441... + 0.14142505258239...*I, 1)]
            sage: p = p^2 * (y^2 - 2)
            sage: p.roots(ring=CIF)
            [(-1.414213562373095?, 1), (1.414213562373095?, 1), (-1.214638932244183? - 0.141425052582394?*I, 2), (-0.141425052582394? + 1.214638932244183?*I, 2), (0.141425052582394? - 1.214638932244183?*I, 2), (1.214638932244183? + 0.141425052582394?*I, 2)]

        Note that one should not use NumPy when wanting high precision
        output as it does not support any of the high precision types::

            sage: R.<x> = RealField(200)[]
            sage: f = x^2 - R(pi)
            sage: f.roots()
            [(-1.7724538509055160272981674833411451827975494561223871282138, 1), (1.7724538509055160272981674833411451827975494561223871282138, 1)]
            sage: f.roots(algorithm='numpy')
            doctest... UserWarning: NumPy does not support arbitrary precision arithmetic.  The roots found will likely have less precision than you expect.
            [(-1.77245385090551..., 1), (1.77245385090551..., 1)]

        We can also find roots over number fields:

            sage: K.<z> = CyclotomicField(15)
            sage: R.<x> = PolynomialRing(K)
            sage: (x^2 + x + 1).roots()
            [(z^5, 1), (-z^5 - 1, 1)]

        There are many combinations of floating-point input and output
        types that work. (Note that some of them are quite pointless like using
        ``algorithm='numpy'`` with high-precision types.)

        ::

            sage: rflds = (RR, RDF, RealField(100))
            sage: cflds = (CC, CDF, ComplexField(100))
            sage: def cross(a, b):
            ...       return list(cartesian_product_iterator([a, b]))
            sage: flds = cross(rflds, rflds) + cross(rflds, cflds) + cross(cflds, cflds)
            sage: for (fld_in, fld_out) in flds:
            ...       x = polygen(fld_in)
            ...       f = x^3 - fld_in(2)
            ...       x2 = polygen(fld_out)
            ...       f2 = x2^3 - fld_out(2)
            ...       for algo in (None, 'pari', 'numpy'):
            ...           rts = f.roots(ring=fld_out, multiplicities=False)
            ...           if fld_in == fld_out and algo is None:
            ...               print fld_in, rts
            ...           for rt in rts:
            ...               assert(abs(f2(rt)) <= 1e-10)
            ...               assert(rt.parent() == fld_out)
            Real Field with 53 bits of precision [1.25992104989487]
            Real Double Field [1.25992104989]
            Real Field with 100 bits of precision [1.2599210498948731647672106073]
            Complex Field with 53 bits of precision [1.25992104989487, -0.62996052494743... - 1.09112363597172*I, -0.62996052494743... + 1.09112363597172*I]
            Complex Double Field [1.25992104989, -0.62996052494... - 1.09112363597*I, -0.62996052494... + 1.09112363597*I]
            Complex Field with 100 bits of precision [1.2599210498948731647672106073, -0.62996052494743658238360530364 - 1.0911236359717214035600726142*I, -0.62996052494743658238360530364 + 1.0911236359717214035600726142*I]

        Note that we can find the roots of a polynomial with algebraic
        coefficients::

            sage: rt2 = sqrt(AA(2))
            sage: rt3 = sqrt(AA(3))
            sage: x = polygen(AA)
            sage: f = (x - rt2) * (x - rt3); f
                x^2 - 3.146264369941973?*x + 2.449489742783178?
            sage: rts = f.roots(); rts
            [(1.414213562373095?, 1), (1.732050807568878?, 1)]
            sage: rts[0][0] == rt2
            True
            sage: f.roots(ring=RealIntervalField(150))
            [(1.414213562373095048801688724209698078569671875376948073176679738?, 1), (1.732050807568877293527446341505872366942805253810380628055806980?, 1)]

        We can handle polynomials with huge coefficients.

        This number doesn't even fit in an IEEE double-precision float, but
        RR and CC allow a much larger range of floating-point numbers::

            sage: bigc = 2^1500
            sage: CDF(bigc)
            +infinity
            sage: CC(bigc)
            3.50746621104340e451

        Polynomials using such large coefficients can't be handled by
        numpy, but pari can deal with them::

            sage: x = polygen(QQ)
            sage: p = x + bigc
            sage: p.roots(ring=RR, algorithm='numpy')
            Traceback (most recent call last):
            ...
            LinAlgError: Array must not contain infs or NaNs
            sage: p.roots(ring=RR, algorithm='pari')
            [(-3.50746621104340e451, 1)]
            sage: p.roots(ring=AA)
            [(-3.5074662110434039?e451, 1)]
            sage: p.roots(ring=QQbar)
            [(-3.5074662110434039?e451, 1)]
            sage: p = bigc*x + 1
            sage: p.roots(ring=RR)
            [(0.000000000000000, 1)]
            sage: p.roots(ring=AA)
            [(-2.8510609648967059?e-452, 1)]
            sage: p.roots(ring=QQbar)
            [(-2.8510609648967059?e-452, 1)]
            sage: p = x^2 - bigc
            sage: p.roots(ring=RR)
            [(-5.92238652153286e225, 1), (5.92238652153286e225, 1)]
            sage: p.roots(ring=QQbar)
            [(-5.9223865215328558?e225, 1), (5.9223865215328558?e225, 1)]

        Algorithms used:

        For brevity, we will use RR to mean any RealField of any precision;
        similarly for RIF, CC, and CIF. Since Sage has no specific
        implementation of Gaussian rationals (or of number fields with
        embedding, at all), when we refer to Gaussian rationals below we
        will accept any number field with defining polynomial
        `x^2+1`, mapping the field generator to +I.

        We call the base ring of the polynomial K, and the ring given by
        the ring= argument L. (If ring= is not specified, then L is the
        same as K.)

        If K and L are floating-point (RDF, CDF, RR, or CC), then a
        floating-point root-finder is used. If L is RDF or CDF then we
        default to using NumPy's roots(); otherwise, we use PARI's
        polroots(). This choice can be overridden with
        algorithm='pari' or algorithm='numpy'. If the algorithm is
        unspecified and NumPy's roots() algorithm fails, then we fall
        back to pari (numpy will fail if some coefficient is infinite,
        for instance).

        If L is AA or RIF, and K is ZZ, QQ, or AA, then the root isolation
        algorithm sage.rings.polynomial.real_roots.real_roots() is used.
        (You can call real_roots() directly to get more control than this
        method gives.)

        If L is QQbar or CIF, and K is ZZ, QQ, AA, QQbar, or the Gaussian
        rationals, then the root isolation algorithm
        sage.rings.polynomial.complex_roots.complex_roots() is used. (You
        can call complex_roots() directly to get more control than this
        method gives.)

        If L is AA and K is QQbar or the Gaussian rationals, then
        complex_roots() is used (as above) to find roots in QQbar, then
        these roots are filtered to select only the real roots.

        If L is floating-point and K is not, then we attempt to change the
        polynomial ring to L (using .change_ring()) (or, if L is complex
        and K is not, to the corresponding real field). Then we use either
        PARI or numpy as specified above.

        For all other cases where K is different than L, we just use
        .change_ring(L) and proceed as below.

        The next method, which is used if K is an integral domain, is to
        attempt to factor the polynomial. If this succeeds, then for every
        degree-one factor a\*x+b, we add -b/a as a root (as long as this
        quotient is actually in the desired ring).

        If factoring over K is not implemented (or K is not an integral
        domain), and K is finite, then we find the roots by enumerating all
        elements of K and checking whether the polynomial evaluates to zero
        at that value.

        .. note::

           We mentioned above that polynomials with multiple roots are
           always ill-conditioned; if your input is given to n bits of
           precision, you should not expect more than n/k good bits
           for a k-fold root. (You can get solutions that make the
           polynomial evaluate to a number very close to zero;
           basically the problem is that with a multiple root, there
           are many such numbers, and it's difficult to choose between
           them.)

           To see why this is true, consider the naive floating-point
           error analysis model where you just pretend that all
           floating-point numbers are somewhat imprecise - a little
           'fuzzy', if you will.  Then the graph of a floating-point
           polynomial will be a fuzzy line.  Consider the graph of
           `(x-1)^3`; this will be a fuzzy line with a
           horizontal tangent at `x=1`, `y=0`. If the
           fuzziness extends up and down by about j, then it will
           extend left and right by about cube_root(j).

        TESTS::

            sage: K.<zeta> = CyclotomicField(2)
            sage: R.<x> = K[]
            sage: factor(x^3-1)
            (x - 1) * (x^2 + x + 1)

        This shows that the issue from :trac:`6237` is fixed::

            sage: R.<u> = QQ[]
            sage: g = -27*u^14 - 32*u^9
            sage: g.roots(CDF, multiplicities=False)
            [-1.03456371594, 0.0, -0.31969776999 - 0.983928563571*I, -0.31969776999 + 0.983928563571*I, 0.836979627962 - 0.608101294789*I, 0.836979627962 + 0.608101294789*I]
            sage: g.roots(CDF)
            [(-1.03456371594, 1), (0.0, 9), (-0.31969776999 - 0.983928563571*I, 1), (-0.31969776999 + 0.983928563571*I, 1), (0.836979627962 - 0.608101294789*I, 1), (0.836979627962 + 0.608101294789*I, 1)]

        This shows that the issue at :trac:`2418` is fixed::

            sage: x = polygen(QQ)
            sage: p = (x^50/2^100 + x^10 + x + 1).change_ring(ComplexField(106))
            sage: rts = (p/2^100).roots(multiplicities=False)
            sage: eps = 2^(-50)   # we test the roots numerically
            sage: [abs(p(rt)) < eps for rt in rts] == [True]*50
            True

        This shows that the issue at :trac:`10901` is fixed::

            sage: a = var('a'); R.<x> = SR[]
            sage: f = x - a
            sage: f.roots(RR)
            Traceback (most recent call last):
            ...
            TypeError: Cannot evaluate symbolic expression to a numeric value.
            sage: f.roots(CC)
            Traceback (most recent call last):
            ...
            TypeError: Cannot evaluate symbolic expression to a numeric value.

        We can find roots of polynomials defined over `\ZZ` or `\QQ`
        over the `p`-adics, see :trac:`15422`::

            sage: R.<x> = ZZ[]
            sage: pol = (x - 1)^2
            sage: pol.roots(Qp(3,5))
            [(1 + O(3^5), 2)]

        This doesn't work if we first change coefficients to `\QQ_p`::

            sage: pol.change_ring(Qp(3,5)).roots()
            Traceback (most recent call last):
            ...
            PrecisionError: p-adic factorization not well-defined since the discriminant is zero up to the requestion p-adic precision

            sage: (pol - 3^6).roots(Qp(3,5))
            [(1 + 2*3^3 + 2*3^4 + O(3^5), 1), (1 + 3^3 + O(3^5), 1)]
            sage: r = pol.roots(Zp(3,5), multiplicities=False); r
            [1 + O(3^5)]
            sage: parent(r[0])
            3-adic Ring with capped relative precision 5
        """
        K = self.parent().base_ring()
        if hasattr(K, '_roots_univariate_polynomial'):
            return K._roots_univariate_polynomial(self, ring=ring, multiplicities=multiplicities, algorithm=algorithm)

        L = K if ring is None else ring

        late_import()

        input_fp = (is_RealField(K)
                    or is_ComplexField(K)
                    or is_RealDoubleField(K)
                    or is_ComplexDoubleField(K))
        output_fp = (is_RealField(L)
                     or is_ComplexField(L)
                     or is_RealDoubleField(L)
                     or is_ComplexDoubleField(L))
        input_complex = (is_ComplexField(K)
                         or is_ComplexDoubleField(K))
        output_complex = (is_ComplexField(L)
                          or is_ComplexDoubleField(L))
        input_gaussian = (isinstance(K, NumberField_quadratic)
                          and list(K.polynomial()) == [1, 0, 1])

        if input_fp and output_fp:
            # allow for possibly using a fast but less reliable
            # floating point algorithm from numpy
            low_prec = is_RealDoubleField(K) or is_ComplexDoubleField(K)
            if algorithm is None:
                if low_prec:
                    algorithm = 'either'
                else:
                    algorithm = 'pari'

            if algorithm != 'numpy' and algorithm != 'either' and algorithm != 'pari':
                raise ValueError("Unknown algorithm '%s'" % algorithm)

            # We should support GSL, too.  We could also support PARI's
            # old Newton-iteration algorithm.

            input_arbprec = (is_RealField(K) or
                             is_ComplexField(K))

            if algorithm == 'numpy' or algorithm == 'either':
                if K.prec() > 53 and L.prec() > 53:
                    from warnings import warn
                    warn('NumPy does not support arbitrary precision arithmetic.  ' +
                         'The roots found will likely have less precision than ' +
                         'you expect.')

                import numpy
                from numpy.linalg.linalg import LinAlgError
                numpy_dtype = ('complex' if input_complex else 'double')
                ty = (complex if input_complex else float)
                coeffs = self.list()
                numpy_array = numpy.array([ty(c) for c in reversed(coeffs)], dtype=numpy_dtype)
                try:
                    ext_rts1 = numpy.roots(numpy_array)
                    rts = []
                    for rt in ext_rts1:
                        rts.append(CDF(rt))
                    rts.sort()
                    ext_rts = rts
                except (ValueError, LinAlgError):
                    if algorithm == 'either':
                        algorithm = 'pari'
                    else:
                        raise

            if algorithm == 'pari':
                if not input_arbprec:
                    self = self.change_ring(CC if input_complex else RR)
                ext_rts = pari(self.monic()).polroots(precision = L.prec())

            if output_complex:
                rts = sort_complex_numbers_for_display([L(root) for root in ext_rts])
            else:
                rts = sorted([L(root.real()) for root in ext_rts if root.imag() == 0])

            rts_mult = []
            j = 0
            while j < len(rts):
                rt = rts[j]
                mult = rts.count(rt)
                rts_mult.append((rt, mult))
                j += mult

            if multiplicities:
                return rts_mult
            else:
                return [rt for (rt, mult) in rts_mult]

        if L != K or is_AlgebraicField_common(L):
            # So far, the only "special" implementations are for real
            # and complex root isolation and for p-adic factorization
            if (is_IntegerRing(K) or is_RationalField(K)
                or is_AlgebraicRealField(K)) and \
                (is_AlgebraicRealField(L) or is_RealIntervalField(L)):

                from sage.rings.polynomial.real_roots import real_roots

                if is_AlgebraicRealField(L):
                    rts = real_roots(self, retval='algebraic_real')
                else:
                    diam = ~(ZZ(1) << L.prec())
                    rts1 = real_roots(self, retval='interval', max_diameter=diam)

                    # We (essentially) promise in the docstring above
                    # that returned intervals will be at least the precision
                    # of the given ring.  But real_roots() does not guarantee
                    # this; for instance, if it returns exactly zero,
                    # it may return this with a low-precision
                    # RealIntervalFieldElement.

                    rts = []
                    for (rt, mult) in rts1:
                        if rt.prec() < L.prec():
                            rt = L(rt)
                        rts.append((rt, mult))

                if multiplicities:
                    return rts
                else:
                    return [rt for (rt, mult) in rts]

            if (is_IntegerRing(K) or is_RationalField(K)
                or is_AlgebraicField_common(K) or input_gaussian) and \
                (is_ComplexIntervalField(L) or is_AlgebraicField_common(L)):

                from sage.rings.polynomial.complex_roots import complex_roots

                if is_ComplexIntervalField(L):
                    rts = complex_roots(self, min_prec=L.prec())
                elif is_AlgebraicField(L):
                    rts = complex_roots(self, retval='algebraic')
                else:
                    rts = complex_roots(self, retval='algebraic_real')

                if multiplicities:
                    return rts
                else:
                    return [rt for (rt, mult) in rts]

            if output_fp and output_complex and not input_gaussian:
                # If we want the complex roots, and the input is not
                # floating point, we convert to a real polynomial
                # (except when the input coefficients are Gaussian rationals).
                if is_ComplexDoubleField(L):
                    real_field = RDF
                else:
                    real_field = RealField(L.prec())

                return self.change_ring(real_field).roots(ring=L, multiplicities=multiplicities, algorithm=algorithm)
            elif is_pAdicRing(L) or is_pAdicField(L):
                p = L.prime()
                n = L.precision_cap()
                try:
                    F = self.factor_padic(p, n)
                except AttributeError:
                    pass
                else:
                    return self.change_ring(L)._roots_from_factorization(F, multiplicities)

            return self.change_ring(L).roots(multiplicities=multiplicities, algorithm=algorithm)

        try:
            if K.is_integral_domain():
                if not K.is_field():
                    try:
                        # get rid of the content of self since we don't need it
                        # and we really don't want to factor it if it's a huge
                        # integer
                        c = self.content()
                        self = self//c
                    except AttributeError:
                        pass
                return self._roots_from_factorization(self.factor(), multiplicities)
            else:
                raise NotImplementedError
        except NotImplementedError:
            if K.is_finite():
                if multiplicities:
                    raise NotImplementedError("root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)")
                else:
                    return [a for a in K if not self(a)]

            raise NotImplementedError("root finding for this polynomial not implemented")

    def _roots_from_factorization(self, F, multiplicities):
        """
        Given a factorization ``F`` of the polynomial ``self``, return
        the roots of ``self``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: pol = 20*x^3 - 50*x^2 + 20*x
            sage: F = pol.factor(); F
            2 * 5 * (x - 2) * x * (2*x - 1)
            sage: pol._roots_from_factorization(F, multiplicities=True)
            [(2, 1), (0, 1)]
            sage: pol.change_ring(QQ)._roots_from_factorization(F, multiplicities=False)
            [2, 0, 1/2]
        """
        seq = []
        K = self.parent().base_ring()
        for fac in F:
            g = fac[0]
            if g.degree() == 1:
                rt = -g[0]/g[1]
                # We need to check that this root is actually in K;
                # otherwise we'd return roots in the fraction field of K.
                if rt in K:
                    rt = K(rt)
                    if multiplicities:
                        seq.append((rt,fac[1]))
                    else:
                        seq.append(rt)
        return seq

    def real_roots(self):
        """
        Return the real roots of this polynomial, without multiplicities.

        Calls self.roots(ring=RR), unless this is a polynomial with
        floating-point real coefficients, in which case it calls
        self.roots().

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: (x^2 - x - 1).real_roots()
            [-0.618033988749895, 1.61803398874989]

        TESTS::

            sage: x = polygen(RealField(100))
            sage: (x^2 - x - 1).real_roots()[0].parent()
                Real Field with 100 bits of precision
            sage: x = polygen(RDF)
            sage: (x^2 - x - 1).real_roots()[0].parent()
            Real Double Field

            sage: x=polygen(ZZ,'x'); v=(x^2-x-1).real_roots()
            sage: v[0].parent() is RR
            True
        """
        K = self.base_ring()
        if is_RealField(K) or is_RealDoubleField(K):
            return self.roots(multiplicities=False)

        return self.roots(ring=RR, multiplicities=False)

    def complex_roots(self):
        """
        Return the complex roots of this polynomial, without
        multiplicities.

        Calls self.roots(ring=CC), unless this is a polynomial with
        floating-point coefficients, in which case it is uses the
        appropriate precision from the input coefficients.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: (x^3 - 1).complex_roots()   # note: low order bits slightly different on ppc.
            [1.00000000000000, -0.500000000000000 - 0.86602540378443...*I, -0.500000000000000 + 0.86602540378443...*I]

        TESTS::

            sage: x = polygen(RR)
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Field with 53 bits of precision
            sage: x = polygen(RDF)
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Double Field
            sage: x = polygen(RealField(200))
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Field with 200 bits of precision
            sage: x = polygen(CDF)
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Double Field
            sage: x = polygen(ComplexField(200))
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Field with 200 bits of precision
            sage: x=polygen(ZZ,'x'); v=(x^2-x-1).complex_roots()
            sage: v[0].parent() is CC
            True
        """
        K = self.base_ring()
        if is_RealField(K):
            return self.roots(ring=ComplexField(K.prec()), multiplicities=False)
        if is_RealDoubleField(K):
            return self.roots(ring=CDF, multiplicities=False)
        if is_ComplexField(K) or is_ComplexDoubleField(K):
            return self.roots(multiplicities=False)

        return self.roots(ring=CC, multiplicities=False)

    def variable_name(self):
        """
        Return name of variable used in this polynomial as a string.

        OUTPUT: string

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = t^3 + 3/2*t + 5
            sage: f.variable_name()
            't'
        """
        return self.parent().variable_name()

    def variables(self):
        """
        Returns the tuple of variables occurring in this polynomial.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.variables()
            (x,)

        A constant polynomial has no variables.

        ::

            sage: R(2).variables()
            ()
        """
        if self.is_constant():
            return ()
        else:
            return self._parent.gens()

    def args(self):
        """
        Returns the generator of this polynomial ring, which is the (only)
        argument used when calling self.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.args()
            (x,)

        A constant polynomial has no variables, but still takes a single
        argument.

        ::

            sage: R(2).args()
            (x,)
        """
        return self._parent.gens()

    def valuation(self, p=None):
        r"""
        If `f = a_r x^r + a_{r+1}x^{r+1} + \cdots`, with
        `a_r` nonzero, then the valuation of `f` is
        `r`. The valuation of the zero polynomial is
        `\infty`.

        If a prime (or non-prime) `p` is given, then the valuation
        is the largest power of `p` which divides self.

        The valuation at `\infty` is -self.degree().

        EXAMPLES::

            sage: P,x=PolynomialRing(ZZ,'x').objgen()
            sage: (x^2+x).valuation()
            1
            sage: (x^2+x).valuation(x+1)
            1
            sage: (x^2+1).valuation()
            0
            sage: (x^3+1).valuation(infinity)
            -3
            sage: P(0).valuation()
            +Infinity
        """
        cdef int k

        if not self:
            return infinity.infinity

        if p is infinity.infinity:
            return -self.degree()

        if p is None:
            for k from 0 <= k <= self.degree():
                if self[k]:
                    return ZZ(k)
        if isinstance(p, Polynomial):
            p = self.parent().coerce(p)
        elif is_Ideal(p) and p.ring() is self.parent(): # eventually need to handle fractional ideals in the fraction field
            if self.parent().base_ring().is_field(): # common case
                p = p.gen()
            else:
                raise NotImplementedError
        else:
            from sage.rings.fraction_field import is_FractionField
            if is_FractionField(p.parent()) and self.parent().has_coerce_map_from(p.parent().ring()):
                p = self.parent().coerce(p.parent().ring()(p)) # here we require that p be integral.
            else:
                raise TypeError("The polynomial, p, must have the same parent as self.")

        if p.degree() == 0:
            raise ArithmeticError("The polynomial, p, must have positive degree.")
        k = 0
        while self % p == 0:
            k = k + 1
            self = self.__floordiv__(p)
        return sage.rings.integer.Integer(k)

    def ord(self, p=None):
        r"""
        This is the same as the valuation of self at p. See the
        documentation for ``self.valuation``.

        EXAMPLES::

            sage: P,x=PolynomialRing(ZZ,'x').objgen()
            sage: (x^2+x).ord(x+1)
            1
        """
        return self.valuation(p)

    def add_bigoh(self, prec):
        r"""
        Returns the power series of precision at most prec got by adding
        `O(q^\text{prec})` to self, where q is its variable.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = 1 + 4*x + x^3
            sage: f.add_bigoh(7)
            1 + 4*x + x^3 + O(x^7)
            sage: f.add_bigoh(2)
            1 + 4*x + O(x^2)
            sage: f.add_bigoh(2).parent()
            Power Series Ring in x over Integer Ring
        """
        return self.parent().completion(self.parent().gen())(self).add_bigoh(prec)

    def name(self):
        """
        Note: This function is deprecated. It will be removed in a future
        release of Sage. Please use the .variable_name() function
        instead.

        Return the string variable name of the indeterminate of this
        polynomial.

        EXAMPLES::

            sage: R.<theta> = ZZ[];
            sage: f = (2-theta)^3; f
            -theta^3 + 6*theta^2 - 12*theta + 8
            sage: f.name()
            doctest:...: DeprecationWarning: This function is deprecated. It will be removed in a future release of Sage. Please use the .variable_name() function instead.
            See http://trac.sagemath.org/4522 for details.
            'theta'
        """
        from sage.misc.superseded import deprecation
        deprecation(4522, "This function is deprecated. It will be removed in a future release of Sage. Please use the .variable_name() function instead.")
        return self.parent().variable_name()

    def is_irreducible(self):
        """
        Return True precisely if this polynomial is irreducible over its
        base ring.

        Testing irreducibility over `\ZZ/n\ZZ` for composite `n` is not
        implemented.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x^3 + 1).is_irreducible()
            False
            sage: (x^2 - 1).is_irreducible()
            False
            sage: (x^3 + 2).is_irreducible()
            True
            sage: R(0).is_irreducible()
            False

        See :trac:`5140`,

        ::

            sage: R(1).is_irreducible()
            False
            sage: R(4).is_irreducible()
            False
            sage: R(5).is_irreducible()
            True

        The base ring does matter:  for example, 2x is irreducible as a
        polynomial in QQ[x], but not in ZZ[x],

        ::

            sage: R.<x> = ZZ[]
            sage: R(2*x).is_irreducible()
            False
            sage: R.<x> = QQ[]
            sage: R(2*x).is_irreducible()
            True

        TESTS::

            sage: F.<t> = NumberField(x^2-5)
            sage: Fx.<xF> = PolynomialRing(F)
            sage: f = Fx([2*t - 5, 5*t - 10, 3*t - 6, -t, -t + 2, 1])
            sage: f.is_irreducible()
            False
            sage: f = Fx([2*t - 3, 5*t - 10, 3*t - 6, -t, -t + 2, 1])
            sage: f.is_irreducible()
            True
        """
        if self.is_zero():
            return False
        if self.is_unit():
            return False
        if self.degree() == 0:
            return self.base_ring()(self).is_irreducible()

        F = self.factor()
        if len(F) > 1 or F[0][1] > 1:
            return False
        return True

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power `x^n`. If
        `n` is negative, terms below `x^n` will be
        discarded. Does not change this polynomial (since polynomials are
        immutable).

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: p = x^2 + 2*x + 4
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(-5)
             0
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        One can also use the infix shift operator::

            sage: f = x^3 + x
            sage: f >> 2
            x
            sage: f << 2
            x^5 + x^3

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True

        AUTHORS:

        - David Harvey (2006-08-06)

        - Robert Bradshaw (2007-04-18): Added support for infix
          operator.
        """
        if n == 0 or self.degree() < 0:
            return self   # safe because immutable.
        if n > 0:
            output = [self.base_ring().zero_element()] * n
            output.extend(self.coeffs())
            return self._parent(output, check=False)
        if n < 0:
            if n > self.degree():
                return self._parent([])
            else:
                return self._parent(self.coeffs()[-int(n):], check=False)

    def __lshift__(self, k):
        return self.shift(k)

    def __rshift__(self, k):
        return self.shift(-k)

    cpdef Polynomial truncate(self, long n):
        r"""
        Returns the polynomial of degree ` < n` which is equivalent
        to self modulo `x^n`.

        EXAMPLES::

            sage: R.<x> = ZZ[]; S.<y> = PolynomialRing(R, sparse=True)
            sage: f = y^3 + x*y -3*x; f
            y^3 + x*y - 3*x
            sage: f.truncate(2)
            x*y - 3*x
            sage: f.truncate(1)
            -3*x
            sage: f.truncate(0)
            0
        """
        # We must not have check=False, since 0 must not have __coeffs = [0].
        return <Polynomial>self._parent(self[:n])#, check=False)

    cdef _inplace_truncate(self, long prec):
        return self.truncate(prec)

    def is_squarefree(self):
        """
        Return False if this polynomial is not square-free, i.e., if there is a
        non-unit `g` in the polynomial ring such that `g^2` divides ``self``.

        .. WARNING::

            This method is not consistent with
            :meth:`.squarefree_decomposition` since the latter does not factor
            the content of a polynomial. See the examples below.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (x-1)*(x-2)*(x^2-5)*(x^17-3); f
            x^21 - 3*x^20 - 3*x^19 + 15*x^18 - 10*x^17 - 3*x^4 + 9*x^3 + 9*x^2 - 45*x + 30
            sage: f.is_squarefree()
            True
            sage: (f*(x^2-5)).is_squarefree()
            False

        A generic implementation is available for polynomials defined over
        principal ideal domains of characteristic 0; the algorithm relies on
        gcd computations::

            sage: R.<x> = ZZ[]
            sage: (2*x).is_squarefree()
            True
            sage: (4*x).is_squarefree()
            False
            sage: (2*x^2).is_squarefree()
            False
            sage: R(0).is_squarefree()
            False

            sage: S.<y> = QQ[]
            sage: R.<x> = S[]
            sage: (2*x*y).is_squarefree() # R does not provide a gcd implementation
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense' object has no attribute 'gcd'
            sage: (2*x*y^2).is_squarefree()
            False

        Over principal ideal domains of positive characteristic, we compute the
        square-free decomposition or a full factorization depending on which is
        available::

            sage: K.<t> = FunctionField(GF(3))
            sage: R.<x> = K[]
            sage: (x^3-x).is_squarefree()
            True
            sage: (x^3-1).is_squarefree()
            False
            sage: (x^3+t).is_squarefree()
            True
            sage: (x^3+t^3).is_squarefree()
            False

        In the following example, `t^2` is a unit in the base field::

            sage: R(t^2).is_squarefree()
            True

        This method is not consistent with :meth:`.squarefree_decomposition`::

            sage: R.<x> = ZZ[]
            sage: f = 4 * x
            sage: f.is_squarefree()
            False
            sage: f.squarefree_decomposition()
            (4) * x

        If you want this method equally not to consider the content, you can
        remove it as in the following example::

            sage: c = f.content()
            sage: (f/c).is_squarefree()
            True

        """
        if self.parent().base_ring() not in sage.categories.principal_ideal_domains.PrincipalIdealDomains():
            raise NotImplementedError("is_squarefree() is only implemented for polynomials over principal ideal domains")

        # a square-free polynomial has a square-free content
        if not self.parent().base_ring().is_field():
            content = self.content()
            if content not in self.parent().base_ring():
                content = content.gen()
            if not content.is_squarefree():
                return False

        # separable polynomials are square-free
        if self.derivative().gcd(self).is_constant():
            return True

        # for characteristic zero rings, square-free polynomials have to be separable
        if self.parent().characteristic().is_zero():
            return False

        # over rings of positive characteristic, we rely on the square-free decomposition if available
        try:
            F = self.squarefree_decomposition()
        # We catch:
        # - NotImplementedError in case squarefree decomposition is not implemented
        # - AttributeError in case p-th roots are not (or do not exist)
        except (NotImplementedError, AttributeError):
            F = self.factor()
        return all([e<=1 for (f,e) in F])

    def radical(self):
        """
        Returns the radical of self; over a field, this is the product of
        the distinct irreducible factors of self. (This is also sometimes
        called the "square-free part" of self, but that term is ambiguous;
        it is sometimes used to mean the quotient of self by its maximal
        square factor.)

        EXAMPLES::

            sage: P.<x> = ZZ[]
            sage: t = (x^2-x+1)^3 * (3*x-1)^2
            sage: t.radical()
            3*x^3 - 4*x^2 + 4*x - 1
            sage: radical(12 * x^5)
            6*x

        If self has a factor of multiplicity divisible by the characteristic (see Trac 8736)::

            sage: P.<x> = GF(2)[]
            sage: (x^3 + x^2).radical()
            x^2 + x
        """
        P = self.parent()
        R = P.base_ring()
        p = R.characteristic()
        if p == 0 or p > self.degree():
            if R.is_field():
                return self // self.gcd(self.derivative())
            else:
                # Be careful with the content: return the
                # radical of the content times the radical of
                # (self/content)
                content = self.content()
                self_1 = (self//content)
                return (self_1 // self_1.gcd(self_1.derivative())) * content.radical()
        else:  # The above method is not always correct (see Trac 8736)
            return self.factor().radical_value()

    def content(self):
        """
        Return the content of ``self``, which is the ideal generated by the coefficients of ``self``.

        EXAMPLES::

            sage: R.<x> = IntegerModRing(4)[]
            sage: f = x^4 + 3*x^2 + 2
            sage: f.content()
            Ideal (2, 3, 1) of Ring of integers modulo 4
        """
        return self.base_ring().ideal(self.coefficients())

    def norm(self, p):
        r"""
        Return the `p`-norm of this polynomial.

        DEFINITION: For integer `p`, the `p`-norm of a
        polynomial is the `p`\th root of the sum of the
        `p`\th powers of the absolute values of the coefficients of
        the polynomial.

        INPUT:


        -  ``p`` - (positive integer or +infinity) the degree
           of the norm


        EXAMPLES::

            sage: R.<x> = RR[]
            sage: f = x^6 + x^2 + -x^4 - 2*x^3
            sage: f.norm(2)
            2.64575131106459
            sage: (sqrt(1^2 + 1^2 + (-1)^2 + (-2)^2)).n()
            2.64575131106459

        ::

            sage: f.norm(1)
            5.00000000000000
            sage: f.norm(infinity)
            2.00000000000000

        ::

            sage: f.norm(-1)
            Traceback (most recent call last):
            ...
            ValueError: The degree of the norm must be positive

        TESTS::

            sage: R.<x> = RR[]
            sage: f = x^6 + x^2 + -x^4 -x^3
            sage: f.norm(int(2))
            2.00000000000000

        AUTHORS:

        - Didier Deshommes

        - William Stein: fix bugs, add definition, etc.
        """
        if p <= 0 :
            raise ValueError("The degree of the norm must be positive")

        coeffs = self.coeffs()
        if p == infinity.infinity:
            return RR(max([abs(i) for i in coeffs]))

        p = sage.rings.integer.Integer(p)  # because we'll do 1/p below.

        if p == 1:
            return RR(sum([abs(i) for i in coeffs]))

        return RR(sum([abs(i)**p for i in coeffs]))**(1/p)

    def hamming_weight(self):
        """
        Returns the number of non-zero coefficients of self.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^3 - x
            sage: f.hamming_weight()
            2
            sage: R(0).hamming_weight()
            0
            sage: f = (x+1)^100
            sage: f.hamming_weight()
            101
            sage: S = GF(5)['y']
            sage: S(f).hamming_weight()
            5
            sage: cyclotomic_polynomial(105).hamming_weight()
            33
        """
        cdef long w = 0
        for a in self.coeffs():
            if a:
                w += 1
        return w

    def map_coefficients(self, f, new_base_ring = None):
        """
        Returns the polynomial obtained by applying ``f`` to the non-zero
        coefficients of self.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        polynomial will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as self. Set
        ``new_base_ring`` to override this behaviour.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of self.

        - ``new_base_ring`` (optional) -- if given, the resulting polynomial
          will be defined over this ring.

        EXAMPLES::

            sage: R.<x> = SR[]
            sage: f = (1+I)*x^2 + 3*x - I
            sage: f.map_coefficients(lambda z: z.conjugate())
            (-I + 1)*x^2 + 3*x + I
            sage: R.<x> = ZZ[]
            sage: f = x^2 + 2
            sage: f.map_coefficients(lambda a: a + 42)
            43*x^2 + 44
            sage: R.<x> = PolynomialRing(SR, sparse=True)
            sage: f = (1+I)*x^(2^32) - I
            sage: f.map_coefficients(lambda z: z.conjugate())
            (-I + 1)*x^4294967296 + I
            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: f = x^(2^32) + 2
            sage: f.map_coefficients(lambda a: a + 42)
            43*x^4294967296 + 44

        Examples with different base ring::

            sage: R.<x> = ZZ[]
            sage: k = GF(2)
            sage: residue = lambda x: k(x)
            sage: f = 4*x^2+x+3
            sage: g = f.map_coefficients(residue); g
            x + 1
            sage: g.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: g = f.map_coefficients(residue, new_base_ring = k); g
            x + 1
            sage: g.parent()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)
            sage: residue = k.coerce_map_from(ZZ)
            sage: g = f.map_coefficients(residue); g
            x + 1
            sage: g.parent()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)
        """
        R = self.parent()
        if new_base_ring is not None:
            R = R.change_ring(new_base_ring)
        elif isinstance(f, Map):
            R = R.change_ring(f.codomain())
        return R(dict([(k,f(v)) for (k,v) in self.dict().items()]))

    def is_cyclotomic(self):
        r"""
        Return ``True`` if ``self`` is a cyclotomic polynomial.

        A cyclotomic polynomial is a monic, irreducible polynomial such that
        all roots are roots of unity.

        ALGORITHM:

        The first cyclotomic polynomial ``x-1`` is treated apart,
        otherwise the first algorithm of [BD89]_ is used.

        EXAMPLES:

        Quick tests::

            sage: P.<x> = ZZ[x]
            sage: (x - 1).is_cyclotomic()
            True
            sage: (x + 1).is_cyclotomic()
            True
            sage: (x^2 - 1).is_cyclotomic()
            False

        Test first 100 cyclotomic polynomials::

            sage: all(cyclotomic_polynomial(i).is_cyclotomic() for i in xrange(1,101))
            True

        Some more tests::

            sage: (x^16 + x^14 - x^10 + x^8 - x^6 + x^2 + 1).is_cyclotomic()
            False
            sage: (x^16 + x^14 - x^10 - x^8 - x^6 + x^2 + 1).is_cyclotomic()
            True
            sage: y = polygen(QQ)
            sage: (y/2 - 1/2).is_cyclotomic()
            False
            sage: (2*(y/2 - 1/2)).is_cyclotomic()
            True

        Test using other rings::

            sage: z = polygen(GF(5))
            sage: (z - 1).is_cyclotomic()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented in non-zero characteristic

        REFERENCES:

        .. [BD89] R. J. Bradford and J. H. Davenport, Effective tests
           for cyclotomic polynomials, Symbolic and Algebraic Computation (1989)
           pp. 244 -- 251, :doi:`10.1007/3-540-51084-2_22`
        """
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError("not implemented in non-zero characteristic")
        if self.base_ring() != ZZ:
            try:
                f = self.change_ring(ZZ)
            except TypeError:
                return False
            return f.is_cyclotomic()

        if not self.is_irreducible():
            return False

        if not self.is_monic():
            return False

        gen = self.parent().gen()

        if (self == gen - 1):  # the first cyc. pol. is treated apart
            return True

        coefs = self.coeffs()
        if coefs[0] != 1:
            return False

        # construct the odd and even part of self
        po_odd = sum(coefs[i]*(gen**((i-1)/2)) for i in xrange(1,len(coefs),2))
        po_even = sum(coefs[i]*(gen**(i/2)) for i in xrange(0,len(coefs),2))

        # f1 = graeffe(self)
        f1  = po_even**2 - gen*(po_odd**2)

        # first case
        if f1 == self:
            return True

        # second case
        selfminus = self(-gen)
        if f1 == selfminus:
            if selfminus.leading_coefficient() < 0 and ((-1)*selfminus).is_cyclotomic():
                return True
            elif selfminus.is_cyclotomic():
                return True

        # third case, we need to factorise
        ff1 = f1.factor()
        if len(ff1) == 1 and ff1[0][1] == 2 and ff1[0][0].is_cyclotomic():
            return True

        # otherwise not cyclotomic
        return False

# ----------------- inner functions -------------
# Cython can't handle function definitions inside other function

cdef _karatsuba_sum(v,w):
    if len(v)>=len(w):
        x = list(v)
        y = w
    else:
        x = list(w)
        y = v
    for i in range(len(y)):
        x[i] = x[i] + y[i]
    return x

cdef _karatsuba_dif(v,w):
    if len(v)>=len(w):
        x = list(v)
        y = w
    else:
        x = list(w)
        y = v
    for i in range(len(y)):
        x[i] -= y[i]
    return x

cdef do_karatsuba(left, right):
    if len(left) == 0 or len(right) == 0:
        return []
    if len(left) == 1:
        c = left[0]
        return [c*a for a in right]
    if len(right) == 1:
        c = right[0]
        return [c*a for a in left]
    if len(left) == 2 and len(right) == 2:
        b = left[0]
        a = left[1]
        d = right[0]
        c = right[1]
        ac = a*c
        bd = b*d
        return [bd,(a+b)*(c+d)-ac-bd,ac]
    e = min(len(left), len(right))/2
    assert e>=1, "bug in karatsuba"
    a, b = left[e:], left[:e]
    c, d = right[e:], right[:e]
    ac = do_karatsuba(a,c)
    bd = do_karatsuba(b,d)
    zeros = [0] * e
    t2 = zeros + zeros + ac
    t1 = zeros + _karatsuba_dif(do_karatsuba(_karatsuba_sum(a,b),_karatsuba_sum(c,d)),_karatsuba_sum(ac,bd))
    t0 = bd
    return _karatsuba_sum(t0,_karatsuba_sum(t1,t2))

cpdef Polynomial_generic_dense _new_constant_dense_poly(list coeffs, Parent P, sample):
    cdef Polynomial_generic_dense f = <Polynomial_generic_dense>PY_NEW_SAME_TYPE(sample)
    f._parent = P
    f.__coeffs = coeffs
    return f

cdef class Polynomial_generic_dense(Polynomial):
    """
    A generic dense polynomial.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(PolynomialRing(QQ,'y'))
        sage: f = x^3 - x + 17
        sage: type(f)
        <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, int check=1, is_gen=False, int construct=0, **kwds):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = []
            return

        R = parent.base_ring()
        if PY_TYPE_CHECK(x, list):
            if check:
                self.__coeffs = [R(t) for t in x]
                self.__normalize()
            else:
                self.__coeffs = x
            return

        if sage.rings.fraction_field_element.is_FractionFieldElement(x):
            if x.denominator() != 1:
                raise TypeError("denominator must be 1")
            else:
                x = x.numerator()

        if PY_TYPE_CHECK(x, Polynomial):
            if (<Element>x)._parent is self._parent:
                x = list(x.list())
            elif R.has_coerce_map_from((<Element>x)._parent):# is R or (<Element>x)._parent == R:
                try:
                    if x.is_zero():
                        self.__coeffs = []
                        return
                except (AttributeError, TypeError):
                    pass
                x = [x]
            else:
                self.__coeffs = [R(a, **kwds) for a in x.list()]
                if check:
                    self.__normalize()
                return

        elif PY_TYPE_CHECK(x, int) and x == 0:
            self.__coeffs = []
            return

        elif isinstance(x, dict):
            x = self._dict_to_list(x, R.zero_element())

        elif isinstance(x, pari_gen):
            x = [R(w, **kwds) for w in x.list()]
            check = 0
        elif not isinstance(x, list):
            # We trust that the element constructors do not send x=0
#            if x:
            x = [x]   # constant polynomials
#            else:
#                x = []    # zero polynomial
        if check:
            self.__coeffs = [R(z, **kwds) for z in x]
            self.__normalize()
        else:
            self.__coeffs = x

    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P):
        cdef Polynomial_generic_dense f = <Polynomial_generic_dense>PY_NEW_SAME_TYPE(self)
        f._parent = P
        f.__coeffs = coeffs
        return f

    cpdef Polynomial _new_constant_poly(self, a, Parent P):
        """
        Create a new constant polynomial in P with value a.

        ASSUMPTION:

        The given value **must** be an element of the base ring. That
        assumption is not verified.

        EXAMPLES::

            sage: S.<y> = QQ[]
            sage: R.<x> = S[]
            sage: x._new_constant_poly(y+1, R)
            y + 1
            sage: parent(x._new_constant_poly(y+1, R))
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in y over Rational Field
        """
        if a:
            return self._new_c([a],P)
        else:
            return self._new_c([],P)

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: R.<x> = QQ['a,b']['x']
            sage: f = x^3-x
            sage: loads(dumps(f)) == f
            True

        Make sure we're testing the right method.
            sage: type(f)
            <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        """
        return make_generic_polynomial, (self._parent, self.__coeffs)

    def __nonzero__(self):
        return len(self.__coeffs) > 0

    cdef void __normalize(self):
        x = self.__coeffs
        cdef Py_ssize_t n = len(x) - 1
        while n >= 0 and not x[n]:
            del x[n]
            n -= 1

    # you may have to replicate this boilerplate code in derived classes if you override
    # __richcmp__.  The python documentation at  http://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: R.<x> = RDF[]
            sage: f = (1+2*x)^5; f
            32.0*x^5 + 80.0*x^4 + 80.0*x^3 + 40.0*x^2 + 10.0*x + 1.0
            sage: f[-1]
            0.0
            sage: f[2]
            40.0
            sage: f[6]
            0.0
            sage: f[:3]
            40.0*x^2 + 10.0*x + 1.0
            sage: f[2:5]
            80.0*x^4 + 80.0*x^3 + 40.0*x^2
            sage: f[2:]
            32.0*x^5 + 80.0*x^4 + 80.0*x^3 + 40.0*x^2
        """
        if isinstance(n, slice):
            start, stop = n.start, n.stop
            if start <= 0:
                start = 0
                zeros = []
            elif start > 0:
                zeros = [self._parent.base_ring().zero_element()] * start
            if stop is None:
                stop = len(self.__coeffs)
            return self._parent(zeros + self.__coeffs[start:stop])
        else:
            if n < 0 or n >= len(self.__coeffs):
                return self.base_ring().zero_element()
            return self.__coeffs[n]

    def _unsafe_mutate(self, n, value):
        """
        Never use this unless you really know what you are doing.

        .. warning::

           This could easily introduce subtle bugs, since Sage assumes
           everywhere that polynomials are immutable. It's OK to use
           this if you really know what you're doing.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = (1+2*x)^2; f
            4*x^2 + 4*x + 1
            sage: f._unsafe_mutate(1, -5)
            sage: f
            4*x^2 - 5*x + 1
        """
        n = int(n)
        value = self.base_ring()(value)
        if n >= 0 and n < len(self.__coeffs):
            self.__coeffs[n] = value
            if n == len(self.__coeffs) and value == 0:
                self.__normalize()
        elif n < 0:
            raise IndexError("polynomial coefficient index must be nonnegative")
        elif value != 0:
            zero = self.base_ring().zero_element()
            for _ in xrange(len(self.__coeffs), n):
                self.__coeffs.append(zero)
            self.__coeffs.append(value)

    def __floordiv__(self, right):
        """
        Return the quotient upon division (no remainder).

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (1+2*x)^3 + 3*x; f
            8*x^3 + 12*x^2 + 9*x + 1
            sage: g = f // (1+2*x); g
            4*x^2 + 4*x + 5/2
            sage: f - g * (1+2*x)
            -3/2
            sage: f.quo_rem(1+2*x)
            (4*x^2 + 4*x + 5/2, -3/2)
        """
        P = (<Element>self)._parent
        if right.parent() == P:
            return Polynomial.__floordiv__(self, right)
        d = P.base_ring()(right)
        cdef Polynomial_generic_dense res = self._new_c([c // d for c in self.__coeffs], P)
        res.__normalize()
        return res

    cpdef ModuleElement _add_(self, ModuleElement right):
        cdef Polynomial_generic_dense res
        cdef Py_ssize_t check=0, i, min
        x = (<Polynomial_generic_dense>self).__coeffs
        y = (<Polynomial_generic_dense>right).__coeffs
        if len(x) > len(y):
            min = len(y)
            high = x[min:]
        elif len(x) < len(y):
            min = len(x)
            high = y[min:]
        else:
            min = len(x)
        cdef list low = [x[i] + y[i] for i from 0 <= i < min]
        if len(x) == len(y):
            res = self._new_c(low, self._parent)
            res.__normalize()
            return res
        else:
            return self._new_c(low + high, self._parent)

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        cdef Py_ssize_t check=0, i, min
        x = (<Polynomial_generic_dense>self).__coeffs
        y = (<Polynomial_generic_dense>right).__coeffs
        if len(x) >= len(y):
            for i from 0 <= i < len(y):
                x[i] += y[i]
        else:
            for i from 0 <= i < len(x):
                x[i] += y[i]
            x += y[len(x):]
        if len(x) == len(y):
            self.__normalize()
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        cdef Polynomial_generic_dense res
        cdef Py_ssize_t check=0, i, min
        x = (<Polynomial_generic_dense>self).__coeffs
        y = (<Polynomial_generic_dense>right).__coeffs
        if len(x) > len(y):
            min = len(y)
            high = x[min:]
        elif len(x) < len(y):
            min = len(x)
            high = [-y[i] for i from min <= i < len(y)]
        else:
            min = len(x)
        low = [x[i] - y[i] for i from 0 <= i < min]
        if len(x) == len(y):
            res = self._new_c(low, self._parent)
            res.__normalize()
            return res
        else:
            return self._new_c(low + high, self._parent)

    cpdef ModuleElement _isub_(self, ModuleElement right):
        cdef Py_ssize_t check=0, i, min
        x = (<Polynomial_generic_dense>self).__coeffs
        y = (<Polynomial_generic_dense>right).__coeffs
        if len(x) >= len(y):
            for i from 0 <= i < len(y):
                x[i] -= y[i]
        else:
            for i from 0 <= i < len(x):
                x[i] -= y[i]
            x += [-c for c in y[len(x):]]
        if len(x) == len(y):
            self.__normalize()
        return self

    cpdef ModuleElement _rmul_(self, RingElement c):
        if len(self.__coeffs) == 0:
            return self
        if c._parent is not (<Element>self.__coeffs[0])._parent:
            c = (<Element>self.__coeffs[0])._parent._coerce_c(c)
        v = [c * a for a in self.__coeffs]
        cdef Polynomial_generic_dense res = self._new_c(v, self._parent)
        #if not v[len(v)-1]:
        # "normalize" checks this anyway...
        res.__normalize()
        return res

    cpdef ModuleElement _lmul_(self, RingElement c):
        if len(self.__coeffs) == 0:
            return self
        if c._parent is not (<Element>self.__coeffs[0])._parent:
            c = (<Element>self.__coeffs[0])._parent._coerce_c(c)
        v = [a * c for a in self.__coeffs]
        cdef Polynomial_generic_dense res = self._new_c(v, self._parent)
        #if not v[len(v)-1]:
        # "normalize" checks this anyway...
        res.__normalize()
        return res

    cpdef ModuleElement _ilmul_(self, RingElement c):
        if len(self.__coeffs) == 0:
            return self
        if c._parent is not (<Element>self.__coeffs[0])._parent:
            c = (<Element>self.__coeffs[0])._parent._coerce_c(c)
        cdef Py_ssize_t i, deg = len(self.__coeffs)
        for i from 0 <= i < deg:
            self.__coeffs[i] *= c
        if not self.__coeffs[deg-1]:
            self.__normalize()
        return self

    cpdef constant_coefficient(self):
        """
        Return the constant coefficient of this polynomial.

        OUTPUT:
            element of base ring

        EXAMPLES:
            sage: R.<t> = QQ[]
            sage: S.<x> = R[]
            sage: f = x*t + x + t
            sage: f.constant_coefficient()
            t
        """
        if len(self.__coeffs) == 0:
            return self.base_ring().zero_element()
        else:
            return self.__coeffs[0]

    def list(self, copy=True):
        """
        Return a new copy of the list of the underlying elements of self.

        EXAMPLES::

            sage: R.<x> = GF(17)[]
            sage: f = (1+2*x)^3 + 3*x; f
            8*x^3 + 12*x^2 + 9*x + 1
            sage: f.list()
            [1, 9, 12, 8]
        """
        if copy:
            return list(self.__coeffs)
        else:
            return self.__coeffs

    def degree(self, gen=None):
        """
        EXAMPLES::

            sage: R.<x> = RDF[]
            sage: f = (1+2*x^7)^5
            sage: f.degree()
            35
        """
        return len(self.__coeffs) - 1

    def shift(self, Py_ssize_t n):
        r"""
        Returns this polynomial multiplied by the power `x^n`. If
        `n` is negative, terms below `x^n` will be
        discarded. Does not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(PolynomialRing(QQ,'y'), 'x')
            sage: p = x^2 + 2*x + 4
            sage: type(p)
            <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True

        AUTHORS:

        - David Harvey (2006-08-06)
        """
        if n == 0 or self.degree() < 0:
            return self
        if n > 0:
            output = [self.base_ring().zero_element()] * n
            output.extend(self.__coeffs)
            return self._new_c(output, self._parent)
        if n < 0:
            if n > len(self.__coeffs) - 1:
                return self._parent([])
            else:
                return self._new_c(self.__coeffs[-int(n):], self._parent)

    @coerce_binop
    def quo_rem(self, other):
        """
        Returns the quotient and remainder of the Euclidean division of
        ``self`` and ``other``.

        Raises ZerodivisionError if ``other`` is zero. Raises ArithmeticError if ``other`` has
        a nonunit leading coefficient.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: R.<y> = P[]
            sage: f = R.random_element(10)
            sage: g = y^5+R.random_element(4)
            sage: q,r = f.quo_rem(g)
            sage: f == q*g + r
            True
            sage: g = x*y^5
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ArithmeticError: Nonunit leading coefficient
            sage: g = 0
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by zero polynomial
        """
        if other.is_zero():
            raise ZeroDivisionError("Division by zero polynomial")
        if not other.leading_coefficient().is_unit():
            raise ArithmeticError("Nonunit leading coefficient")
        if self.is_zero():
            return self, self

        R = self.parent().base_ring()
        x = (<Polynomial_generic_dense>self).__coeffs[:] # make a copy
        y = (<Polynomial_generic_dense>other).__coeffs
        m = len(x)  # deg(self)=m-1
        n = len(y)  # deg(other)=n-1
        if m < n:
            return self.parent()(0), self

        quo = list()
        for k from m-n >= k >= 0:
            q = x[n+k-1]/y[n-1]
            x[n+k-1] = R.zero_element()
            for j from n+k-2 >= j >= k:
                x[j] -= q * y[j-k]
            quo.insert(0,q)

        return self._new_c(quo,self._parent), self._new_c(x,self._parent)._inplace_truncate(n-1)

    cpdef Polynomial truncate(self, long n):
        r"""
        Returns the polynomial of degree ` < n` which is equivalent
        to self modulo `x^n`.

        EXAMPLES::

            sage: S.<q> = QQ['t']['q']
            sage: f = (1+q^10+q^11+q^12).truncate(11); f
            q^10 + 1
            sage: f = (1+q^10+q^100).truncate(50); f
            q^10 + 1
            sage: f.degree()
            10
            sage: f = (1+q^10+q^100).truncate(500); f
            q^100 + q^10 + 1

        TESTS:

        Make sure we're not actually testing a specialized
        implementation.

        ::

            sage: type(f)
            <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        """
        l = len(self.__coeffs)
        if n > l:
            n = l
        while n > 0 and not self.__coeffs[n-1]:
            n -= 1
        return self._new_c(self.__coeffs[:n], self._parent)

    cdef _inplace_truncate(self, long n):
        if n < len(self.__coeffs):
            while n > 0 and not self.__coeffs[n-1]:
                n -= 1
        self.__coeffs = self.__coeffs[:n]
        return self

def make_generic_polynomial(parent, coeffs):
    return parent(coeffs)

cdef class ConstantPolynomialSection(Map):
    """
    This class is used for conversion from a polynomial ring to its base ring.

    Since trac ticket #9944, it calls the constant_coefficient method,
    which can be optimized for a particular polynomial type.

    EXAMPLES::

        sage: P0.<y_1> = GF(3)[]
        sage: P1.<y_2,y_1,y_0> = GF(3)[]
        sage: P0(-y_1)    # indirect doctest
        2*y_1
        sage: phi = GF(3).convert_map_from(P0); phi
        Generic map:
          From: Univariate Polynomial Ring in y_1 over Finite Field of size 3
          To:   Finite Field of size 3
        sage: type(phi)
        <type 'sage.rings.polynomial.polynomial_element.ConstantPolynomialSection'>
        sage: phi(P0.one_element())
        1
        sage: phi(y_1)
        Traceback (most recent call last):
        ...
        TypeError: not a constant polynomial

    """
    cpdef Element _call_(self, x):
        """
        TESTS:
            sage: from sage.rings.polynomial.polynomial_element import ConstantPolynomialSection
            sage: R.<x> = QQ[]
            sage: m = ConstantPolynomialSection(R, QQ); m
            Generic map:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
            sage: m(x-x+1/2) # implicit
            1/2
            sage: m(x-x)
            0
            sage: m(x)
            Traceback (most recent call last):
            ...
            TypeError: not a constant polynomial
        """
        if x.degree() <= 0:
            try:
                return <Element>(x.constant_coefficient())
            except AttributeError:
                return <Element>((<Polynomial>x).constant_coefficient())
        else:
            raise TypeError("not a constant polynomial")

cdef class PolynomialBaseringInjection(Morphism):
    """
    This class is used for conversion from a ring to a polynomial
    over that ring.

    It calls the _new_constant_poly method on the generator,
    which should be optimized for a particular polynomial type.

    Technically, it should be a method of the polynomial ring, but
    few polynomial rings are cython classes, and so, as a method
    of a cython polynomial class, it is faster.

    EXAMPLES:

    We demonstrate that most polynomial ring classes use
    polynomial base injection maps for coercion. They are
    supposed to be the fastest maps for that purpose. See
    trac ticket #9944::

        sage: R.<x> = Qp(3)[]
        sage: R.coerce_map_from(R.base_ring())
        Polynomial base injection morphism:
          From: 3-adic Field with capped relative precision 20
          To:   Univariate Polynomial Ring in x over 3-adic Field with capped relative precision 20
        sage: R.<x,y> = Qp(3)[]
        sage: R.coerce_map_from(R.base_ring())
        Polynomial base injection morphism:
          From: 3-adic Field with capped relative precision 20
          To:   Multivariate Polynomial Ring in x, y over 3-adic Field with capped relative precision 20
        sage: R.<x,y> = QQ[]
        sage: R.coerce_map_from(R.base_ring())
        Polynomial base injection morphism:
          From: Rational Field
          To:   Multivariate Polynomial Ring in x, y over Rational Field
        sage: R.<x> = QQ[]
        sage: R.coerce_map_from(R.base_ring())
        Polynomial base injection morphism:
          From: Rational Field
          To:   Univariate Polynomial Ring in x over Rational Field

    By trac ticket #9944, there are now only very few exceptions::

        sage: PolynomialRing(QQ,names=[]).coerce_map_from(QQ)
        Generic morphism:
          From: Rational Field
          To:   Multivariate Polynomial Ring in no variables over Rational Field

    """

    cdef RingElement _an_element
    cdef object _new_constant_poly_

    def __init__(self, domain, codomain):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
            sage: PolynomialBaseringInjection(QQ, QQ['x'])
            Polynomial base injection morphism:
              From: Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field
            sage: PolynomialBaseringInjection(ZZ, QQ['x'])
            Traceback (most recent call last):
            ...
            AssertionError: domain must be basering

        ::

            sage: R.<t> = Qp(2)[]
            sage: f = R.convert_map_from(R.base_ring())    # indirect doctest
            sage: f(Qp(2).one()*3)
            (1 + 2 + O(2^20))
            sage: (Qp(2).one()*3)*t
            (1 + 2 + O(2^20))*t

        """
        assert codomain.base_ring() is domain, "domain must be basering"
        Morphism.__init__(self, domain, codomain)
        self._an_element = codomain.gen()
        self._repr_type_str = "Polynomial base injection"
        self._new_constant_poly_ = self._an_element._new_constant_poly

    cpdef Element _call_(self, x):
        """
        TESTS:
            sage: from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
            sage: m = PolynomialBaseringInjection(ZZ, ZZ['x']); m
            Polynomial base injection morphism:
              From: Integer Ring
              To:   Univariate Polynomial Ring in x over Integer Ring
            sage: m(2) # indirect doctest
            2
            sage: parent(m(2))
            Univariate Polynomial Ring in x over Integer Ring
        """
        return self._new_constant_poly_(x, self._codomain)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        TESTS:
            sage: from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
            sage: m = PolynomialBaseringInjection(Qp(5), Qp(5)['x'])
            sage: m(1 + O(5^11), absprec = 5)   # indirect doctest
            (1 + O(5^11))
        """
        try:
            return self._codomain._element_constructor_(x, *args, **kwds)
        except AttributeError:
            # if there is no element constructor,
            # there is a custom call method.
            return self._codomain(x, *args, **kwds)

    def section(self):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
            sage: m = PolynomialBaseringInjection(RDF, RDF['x'])
            sage: m.section()
            Generic map:
              From: Univariate Polynomial Ring in x over Real Double Field
              To:   Real Double Field
            sage: type(m.section())
            <type 'sage.rings.polynomial.polynomial_element.ConstantPolynomialSection'>
        """
        return ConstantPolynomialSection(self._codomain, self._domain)
