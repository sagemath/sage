# sage_setup: distribution = sagemath-categories
"""
Univariate polynomial base class

TESTS::

    sage: R.<x> = ZZ[]
    sage: f = x^5 + 2*x^2 + (-1)
    sage: f == loads(dumps(f))
    True

    sage: PolynomialRing(ZZ,'x').objgen()
    (Univariate Polynomial Ring in x over Integer Ring, x)

AUTHORS:

- William Stein: first version

- Martin Albrecht: added singular coercion

- Robert Bradshaw: moved Polynomial_generic_dense to Cython

- Miguel Marco: implemented resultant in the case where PARI fails

- Simon King: used a faster way of conversion from the base ring

- Kwankyu Lee (2013-06-02): enhanced :meth:`quo_rem`

- Julian Rueth (2012-05-25,2014-05-09): fixed is_squarefree() for imperfect
  fields, fixed division without remainder over QQbar; added ``_cache_key``
  for polynomials with unhashable coefficients

- Simon King (2013-10): implemented copying of :class:`PolynomialBaseringInjection`

- Bruno Grenet (2014-07-13): enhanced :meth:`quo_rem`

- Kiran Kedlaya (2016-03): added root counting

- Edgar Costa (2017-07): added rational reconstruction

- Kiran Kedlaya (2017-09): added reciprocal transform, trace polynomial

- David Zureick-Brown (2017-09): added is_weil_polynomial

- Sebastian Oehms (2018-10): made :meth:`roots` and  :meth:`factor` work over more
  cases of proper integral domains (see :issue:`26421`)
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef ZZ, QQ, RR, CC, RDF, CDF

cimport cython
from cpython.number cimport PyNumber_Check

import operator
import copy
import re
from io import StringIO

from sage.cpython.wrapperdescr cimport wrapperdescr_fastcall
import sage.rings.rational
import sage.rings.integer
import sage.rings.integer_ring
import sage.rings.rational_field
import sage.rings.fraction_field_element
import sage.rings.infinity as infinity
from sage.misc.latex import latex
from sage.arith.power cimport generic_power
from sage.arith.misc import crt
from sage.arith.long cimport pyobject_to_long
from sage.misc.mrange import cartesian_product_iterator
from sage.structure.factorization import Factorization
from sage.structure.richcmp cimport (richcmp, richcmp_item,
        rich_to_bool, rich_to_bool_sgn)

try:
    from sage.libs.pari.all import pari, pari_gen, PariError
except ImportError:
    pari_gen = ()
    pari = None

    class PariError(Exception):
        pass

cimport sage.rings.abc

from sage.rings.real_double import RDF
import sage.rings.abc

import sage.interfaces.abc

from sage.structure.coerce cimport coercion_model
from sage.structure.element import coerce_binop
from sage.structure.element cimport (parent, have_same_parent,
        Element, RingElement)

from sage.rings.rational_field import QQ, RationalField
from sage.rings.integer_ring import ZZ, IntegerRing_class
from sage.rings.integer cimport Integer, smallInteger
from sage.libs.gmp.mpz cimport *
from sage.rings.fraction_field import FractionField
from sage.rings.padics.padic_generic import pAdicGeneric

from sage.structure.category_object cimport normalize_names

from sage.misc.derivative import multi_derivative

from sage.arith.misc import sort_complex_numbers_for_display, power_mod, is_prime
from sage.arith.functions import lcm

from sage.rings.polynomial import polynomial_fateman

from sage.rings.ideal import Ideal_generic
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement
from sage.misc.cachefunc import cached_function

from sage.categories.map cimport Map
from sage.categories.morphism cimport Morphism

from sage.misc.superseded import deprecation_cython as deprecation, deprecated_function_alias
from sage.misc.cachefunc import cached_method

cpdef is_Polynomial(f):
    """
    Return ``True`` if ``f`` is of type univariate polynomial.

    This function is deprecated.

    INPUT:

    - ``f`` -- an object

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_element import is_Polynomial
        sage: R.<x> = ZZ[]
        sage: is_Polynomial(x^3 + x + 1)
        doctest:...: DeprecationWarning: the function is_Polynomial is deprecated;
        use isinstance(x, sage.rings.polynomial.polynomial_element.Polynomial) instead
        See https://github.com/sagemath/sage/issues/32709 for details.
        True
        sage: S.<y> = R[]
        sage: f = y^3 + x*y - 3*x; f
        y^3 + x*y - 3*x
        sage: is_Polynomial(f)
        True

    However this function does not return ``True`` for genuine multivariate
    polynomial type objects or symbolic polynomials, since those are
    not of the same data type as univariate polynomials::

        sage: R.<x,y> = QQ[]
        sage: f = y^3 + x*y - 3*x; f
        y^3 + x*y - 3*x
        sage: is_Polynomial(f)
        False

        sage: # needs sage.symbolic
        sage: var('x,y')
        (x, y)
        sage: f = y^3 + x*y - 3*x; f
        y^3 + x*y - 3*x
        sage: is_Polynomial(f)
        False
    """
    deprecation(32709, "the function is_Polynomial is deprecated; use isinstance(x, sage.rings.polynomial.polynomial_element.Polynomial) instead")

    return isinstance(f, Polynomial)

from sage.rings.polynomial.polynomial_compiled cimport CompiledPolynomialFunction

from sage.rings.polynomial.polydict cimport ETuple


cdef class Polynomial(CommutativePolynomial):
    """
    A polynomial.

    EXAMPLES::

        sage: R.<y> = QQ['y']
        sage: S.<x> = R['x']
        sage: S
        Univariate Polynomial Ring in x over Univariate Polynomial Ring in y
        over Rational Field
        sage: f = x*y; f
        y*x
        sage: type(f)
        <class 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        sage: p = (y+1)^10; p(1)
        1024

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _lmul_
    .. automethod:: _rmul_
    .. automethod:: _mul_
    .. automethod:: _mul_trunc_
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

    cdef Polynomial _new_generic(self, list coeffs):
        r"""
        Quickly construct a new polynomial of the same type as ``self``,
        bypassing the parent's element constructor.

        The new polynomial takes ownership of the coefficient list
        given on input.
        """
        cdef Py_ssize_t n = len(coeffs) - 1
        while n >= 0 and not coeffs[n]:
            del coeffs[n]
            n -= 1
        return type(self)(self._parent, coeffs, check=False)

    cpdef _add_(self, right):
        r"""
        Add two polynomials.

        EXAMPLES::

            sage: R = ZZ['x']
            sage: p = R([1,2,3,4])
            sage: q = R([4,-3,2,-1])
            sage: p + q    # indirect doctest
            3*x^3 + 5*x^2 - x + 5
        """
        cdef Py_ssize_t i, min
        cdef list x = self.list(copy=False)
        cdef list y = right.list(copy=False)

        if len(x) > len(y):
            min = len(y)
            high = x[min:]
        elif len(x) < len(y):
            min = len(x)
            high = y[min:]
        else:
            min = len(x)
            high = []

        low = [x[i] + y[i] for i in range(min)]
        return self._new_generic(low + high)

    cpdef _neg_(self):
        return self._new_generic([-x for x in self.list(copy=False)])

    cpdef bint is_zero(self) except -1:
        r"""
        Test whether this polynomial is zero.

        EXAMPLES::

            sage: R = GF(2)['x']['y']
            sage: R([0,1]).is_zero()
            False
            sage: R([0]).is_zero()
            True
            sage: R([-1]).is_zero()
            False
        """
        return self.degree() < 0

    cpdef bint is_one(self) except -1:
        r"""
        Test whether this polynomial is 1.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (x - 3).is_one()
            False
            sage: R(1).is_one()
            True

            sage: R2.<y> = R[]
            sage: R2(x).is_one()
            False
            sage: R2(1).is_one()
            True
            sage: R2(-1).is_one()
            False
        """
        return self.degree() == 0 and self.get_unsafe(0).is_one()

    def plot(self, xmin=None, xmax=None, *args, **kwds):
        """
        Return a plot of this polynomial.

        INPUT:

        - ``xmin`` -- float

        - ``xmax`` -- float

        - ``*args, **kwds`` -- passed to either plot or point

        OUTPUT: a graphic object

        EXAMPLES::

            sage: x = polygen(GF(389))
            sage: plot(x^2 + 1, rgbcolor=(0,0,1))                                       # needs sage.plot
            Graphics object consisting of 1 graphics primitive
            sage: x = polygen(QQ)
            sage: plot(x^2 + 1, rgbcolor=(1,0,0))                                       # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        R = self.base_ring()
        from sage.plot.all import plot, point
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

    cpdef _lmul_(self, Element left):
        """
        Multiply ``self`` on the left by a scalar.

        EXAMPLES::

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
        if not left:
            return self._parent.zero()
        return self._parent(left) * self

    cpdef _rmul_(self, Element right):
        """
        Multiply ``self`` on the right by a scalar.

        EXAMPLES::

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
        if not right:
            return self._parent.zero()
        return self * self._parent(right)

    def subs(self, in_dict=None, *args, **kwds):
        r"""
        Substitute the variable in ``self``.

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
        if not in_dict:
            return self(*args, **kwds)

        if isinstance(in_dict, dict):
            if len(in_dict) > 1:
                raise TypeError("only the generator can be substituted, use __call__ instead")
            k, v = next(iter(in_dict.items()))
            if not k or k == self._parent.gen():
                return self(v, *args, **kwds)
            raise TypeError("keys do not match self's parent")

        return self(in_dict, *args, **kwds)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __call__(self, *args, **kwds):
        """
        Evaluate this polynomial.

        INPUT:

        - ``*args`` -- ring elements, need not be in the coefficient ring of
            the polynomial. The **first** positional argument is substituted
            for the polynomial's indeterminate. Remaining arguments, if any,
            are used **from left to right** to evaluate the coefficients.
        - ``**kwds`` -- variable name-value pairs

        OUTPUT: the value of the polynomial at the point specified by the arguments

        ALGORITHM:

        By default, use Horner's method or create a
        :class:`~sage.rings.polynomial.polynomial_compiled.CompiledPolynomialFunction`
        depending on the polynomial's degree.

        Element classes may define a method called ``_evaluate_polynomial``
        to provide an alternative evaluation algorithm for a given argument
        type. Note that ``_evaluate_polynomial`` may not always be used:
        for instance, subclasses dedicated to specific coefficient rings
        typically do not call it.

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

            sage: # needs sage.combinat sage.modules
            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1, -1)
            sage: R.<w> = PolynomialRing(A, sparse=True)
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

            sage: # needs sage.combinat sage.modules
            sage: R.<x> = QQ[]
            sage: f = R(2/3)
            sage: a = matrix(ZZ, 2)
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
            sage: f(x=10)  # x isn't mentioned
            w^3 + 3*w + 2

        Nested polynomial ring elements can be called like multivariate
        polynomials. Note the order of the arguments::

            sage: R.<x> = QQ[]; S.<y> = R[]
            sage: f = x + y*x + y^2
            sage: f.parent()
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: f(2)
            3*x + 4
            sage: f(2,4)
            16
            sage: f(y=2, x=4)
            16
            sage: f(2, x=4)
            16
            sage: f(2, x=4, z=5)
            16
            sage: f(2, 4, z=10)
            16
            sage: f(y=x)
            2*x^2 + x
            sage: f(x=y)
            2*y^2 + y

        Also observe that ``f(y0, x0)`` means ``f(x=x0)(y=y0)``, not
        ``f(y=y0)(x=x0)``. The two expressions may take different values::

            sage: f(y, x)
            y^2 + x*y + x
            sage: f(y)(x)
            2*x^2 + x

        Polynomial ring elements can also, like multivariate
        polynomials, be called with an argument that is a list or
        tuple containing the values to be substituted, though it is
        perhaps more natural to just unpack the list::

            sage: f([2])   # calling with a list
            3*x + 4
            sage: f((2,))  # calling with a tuple
            3*x + 4
            sage: f(*[2])  # unpacking the list to call normally
            3*x + 4

        The following results in an element of the symbolic ring. ::

            sage: f(x=sqrt(2))                                                          # needs sage.symbolic
            y^2 + sqrt(2)*y + sqrt(2)

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
            sage: parent(f(Qp(5)(0)))                                                   # needs sage.rings.padics
            5-adic Field with capped relative precision 20

        TESTS:

        One test for a simple evaluation::

            sage: x, y = polygens(ZZ, 'x,y')
            sage: t = polygen(x.parent(), 't')
            sage: F = x*y*t
            sage: F(y=1)
            x*t

        The following shows that :issue:`2360` is indeed fixed. ::

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

        The following shows that :issue:`9006` is also fixed. ::

            sage: f = ZZ['x'](1000000 * [1])
            sage: f(1)
            1000000

        The following test came up in :issue:`9051`::

            sage: # needs sage.rings.complex_interval_field
            sage: Cif = ComplexIntervalField(64)
            sage: R.<x> = Cif[]
            sage: f = 2*x-1
            sage: jj = Cif(RIF(0,2))
            sage: f(jj).center(), f(jj).diameter()
            (1.00000000000000000, 4.00000000000000000)

        The following failed before the patch to :issue:`3979`

        ::

            sage: R.<x> = ZZ[]
            sage: S.<y> = R[]
            sage: g = x*y + 1
            sage: g(x=3)
            3*y + 1

        ::

            sage: Pol_x.<x> = QQ[]
            sage: Pol_xy.<y> = Pol_x[]
            sage: pol = 1000*x^2*y^2 + 100*y + 10*x + 1

            sage: pol(y, 0)
            100*y + 1

            sage: pol(~y, 0)
            (y + 100)/y

            sage: pol(y=x, x=1)
            1000*x^2 + 100*x + 11

            sage: zero = Pol_xy(0)
            sage: zero(1).parent()
            Univariate Polynomial Ring in x over Rational Field

            sage: zero = QQ['x'](0)
            sage: a = matrix(ZZ, [[1]])                                                 # needs sage.modules
            sage: zero(a).parent()                                                      # needs sage.modules
            Full MatrixSpace of 1 by 1 dense matrices over Rational Field

            sage: pol(y, x).parent() is pol(x, y).parent() is pol(y, y).parent() is Pol_xy
            True

            sage: pol(x, x).parent()
            Univariate Polynomial Ring in x over Rational Field

            sage: one = Pol_xy(1)
            sage: one(1, 1.).parent()                                                   # needs sage.rings.real_mpfr
            Real Field with 53 bits of precision

            sage: zero = GF(2)['x'](0)
            sage: zero(1.).parent()                                                     # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Finite Field of size 2' and 'Real Field with 53 bits of precision'

            sage: pol(x, y, x=1)
            Traceback (most recent call last):
            ...
            TypeError: Wrong number of arguments

        Check that :issue:`22317` is fixed::

            sage: R = ZZ['x']['y']['z']
            sage: d = R.gens_dict_recursive()
            sage: p = d['x'] * d['z']
            sage: p(x=QQ(0))
            0

        Check that :issue:`32513` is fixed::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: x.substitute([])
            x
            sage: Polynomial.__call__(x, [])
            x

        These were drastically slower prior to :issue:`33165`::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(31337)[]
            sage: f = R(list(range(100, 201)))
            sage: g = R(list(range(1, 1001)))
            sage: S.<y> = R.quotient(f)
            sage: g(y)
            22537*y^99 + 4686*y^98 + 13285*y^97 + 4216*y^96 + ... + 6389*y^3 + 30062*y^2 + 13755*y + 11875

        ::

            sage: # needs sage.rings.finite_rings
            sage: T.<z> = GF(31337)[]
            sage: g(z)
            1000*z^999 + 999*z^998 + 998*z^997 + 997*z^996 + ... + 5*z^4 + 4*z^3 + 3*z^2 + 2*z + 1
            sage: g(z^2)
            1000*z^1998 + 999*z^1996 + 998*z^1994 + 997*z^1992 + ... + 5*z^8 + 4*z^6 + 3*z^4 + 2*z^2 + 1
            sage: g(T([0, 1]))
            1000*z^999 + 999*z^998 + 998*z^997 + 997*z^996 + ... + 5*z^4 + 4*z^3 + 3*z^2 + 2*z + 1
            sage: g(T.zero())
            1
            sage: g(T(2))
            23069

        ::

            sage: # needs sage.rings.finite_rings
            sage: U.<u,v> = GF(31337)[]
            sage: g(u)
            1000*u^999 + 999*u^998 + 998*u^997 + 997*u^996 + ... + 5*u^4 + 4*u^3 + 3*u^2 + 2*u + 1
            sage: g(u*v^2)
            1000*u^999*v^1998 + 999*u^998*v^1996 + 998*u^997*v^1994 + ... + 4*u^3*v^6 + 3*u^2*v^4 + 2*u*v^2 + 1
            sage: g(U.zero())
            1
            sage: g(U(2))
            -8268

        Sparse tests for :issue:`33165`::

            sage: R.<x> = PolynomialRing(QQ, sparse=True)
            sage: f = x^1000000 + 1
            sage: S.<y> = PolynomialRing(QQ, sparse=True)
            sage: f(y)
            y^1000000 + 1
            sage: f(y^100)
            y^100000000 + 1

            sage: U.<u,v> = PolynomialRing(QQ, sparse=True)
            sage: g = x^10000 + 1
            sage: g(u)
            u^10000 + 1
            sage: g(u*v^2)
            u^10000*v^20000 + 1

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
        cdef long i, j, d, deg
        cdef Polynomial pol = self
        cdef ETuple etup
        cdef list cs
        cdef dict coeff_sparse, coeff_dict

        deg = self.degree()
        if deg < 0:
            top = self._parent._base.one()
            cst = self._parent._base.zero()
        else:
            top = self.get_unsafe(deg)
            cst = self.get_unsafe(0)

        a = args[0] if len(args) == 1 else None
        if kwds or not (isinstance(a, Element) or PyNumber_Check(a)):
            # slow path

            # Isolate the variable we are interested in, check remaining
            # arguments

            a = kwds.pop(self.variable_name(), None)
            if args:
                if a is not None:
                    raise TypeError("unsupported mix of keyword and positional arguments")
                if isinstance(args[0], (list, tuple)):
                    if len(args) > 1:
                        raise TypeError("invalid arguments")
                    args = args[0]
                if args:
                    a, args = args[0], args[1:]
            if a is None:
                a = self._parent.gen()

            eval_coeffs = False
            if args or kwds:
                try:
                    # Note that we may be calling a different implementation that
                    # is more permissive about its arguments than we are.
                    top = top(*args, **kwds)
                except TypeError:
                    if args: # bwd compat: nonsense *keyword* arguments are okay
                        raise TypeError("Wrong number of arguments")
                else:
                    eval_coeffs = True

            # Evaluate the coefficients, then fall through to evaluate the
            # resulting univariate polynomial

            if eval_coeffs:
                new_base = parent(top)
                # tentative common parent of the evaluated coefficients
                pol = pol.map_coefficients(lambda c: c(*args, **kwds),
                                           new_base_ring=new_base)
                cst = cst(*args, **kwds)

        R = parent(a)

        if R is not pol._parent._base:
            # If a is a monomial in an isomorphic polynomial ring (or quotient
            # of an isomorphic polynomial ring), the code below wastes a lot of
            # time on conversions and unnecessary intermediate reductions.
            # We can do dramatically better by simply doing the conversion once
            # and returning the result. See #33165.
            if isinstance(a, Polynomial) and a.base_ring() is pol._parent._base:
                if (<Polynomial> a).is_gen():
                    return R(pol)
                d = (<Polynomial> a).degree()
                if d < 0:  # f(0)
                    return R(cst)
                if d == 0:  # f(const)
                    return R(pol((<Polynomial> a).constant_coefficient()))
                if (<Polynomial> a).is_monomial():
                    if d == 1:  # f(x); is_gen() does not catch all such cases
                        return R(pol)
                    if pol._parent.is_sparse():
                        coeff_sparse = {}
                        coeff_dict = <dict> pol.monomial_coefficients()
                        for i in coeff_dict:
                            coeff_sparse[i * d] = coeff_dict[i]
                        return R(coeff_sparse, check=False)
                    # Dense version
                    cs = [pol._parent._base.zero()] * (pol.degree() * d + 1)
                    j = 0
                    for i in range(pol.degree()+1):
                        cs[j] = pol.get_unsafe(i)
                        j += d
                    return R(cs, check=False)

            elif isinstance(a, MPolynomial) and a.base_ring() is pol.base_ring() and a.is_monomial():
                if (<MPolynomial> a).is_constant():
                    if (<MPolynomial> a).is_zero():  # f(0) is just the constant term of f
                        return R(cst)
                    # f(const), so do the evaluation in the base ring
                    return R(pol((<MPolynomial> a).constant_coefficient()))
                if a.is_monomial():
                    etup = ETuple((<MPolynomial> a).degrees())
                    if pol._parent.is_sparse():
                        coeff_dict = <dict> pol.monomial_coefficients()
                        return R({etup.emul(key): coeff_dict[key] for key in coeff_dict})
                    # Dense version
                    return R({etup.emul(i): c for i, c in enumerate(pol)})

            elif isinstance(a, PolynomialQuotientRingElement) and a.lift().is_gen() and R.polynomial_ring() is pol._parent:
                return R(pol)

            # Coerce a once and for all to a parent containing the coefficients.
            # This can save lots of coercions when the common parent is the
            # polynomial's base ring (e.g., for evaluations at integers).
            cst, aa = coercion_model.canonical_coercion(cst, a)
            # Use fast multiplication actions like matrix × scalar.
            # If there is no action, replace a by an element of the
            # target parent.
            S = parent(aa)
            if coercion_model.get_action(S, R) is None:
                a = aa
                R = S

        d = pol.degree()

        if d <= 0 or (isinstance(a, Element) and R.is_exact() and a.is_zero()):
            return cst # with the right parent thanks to the above coercion
        elif pol._parent is R and a.is_gen():
            return pol
        elif hasattr(a, "_evaluate_polynomial"):
            try:
                return a._evaluate_polynomial(pol)
            except NotImplementedError:
                pass

        if pol._compiled is None:
            if d < 4 or d > 50000:
                result = pol.get_unsafe(d)
                for i in range(d - 1, -1, -1):
                    result = result * a + pol.get_unsafe(i)
                return result
            pol._compiled = CompiledPolynomialFunction(pol.list())
        return pol._compiled.eval(a)

    def compose_trunc(self, Polynomial other, long n):
        r"""
        Return the composition of ``self`` and ``other``, truncated to `O(x^n)`.

        This method currently works for some specific coefficient rings only.

        EXAMPLES::

            sage: Pol.<x> = CBF[]                                                       # needs sage.libs.flint
            sage: (1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120).compose_trunc(1 + x, 2)    # needs sage.libs.flint
            ([2.708333333333333 +/- ...e-16])*x + [2.71666666666667 +/- ...e-15]

            sage: Pol.<x> = QQ['y'][]
            sage: (1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120).compose_trunc(1 + x, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: truncated composition is not implemented
            for this subclass of polynomials
        """
        raise NotImplementedError("truncated composition is not implemented "
                                  "for this subclass of polynomials")

    def _compile(self):
        # For testing
        self._compiled = CompiledPolynomialFunction(self.list())
        return self._compiled

    def _get_compiled(self):
        # For testing
        return self._compiled

    def _fast_callable_(self, etb):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['t'])
            sage: R.<t> = QQ[]
            sage: v = -1/4*t^6 + 1/2*t^5 - t^4 - 12*t^3 + 1/2*t^2 - 1/95*t - 1/2
            sage: v._fast_callable_(etb)
            add(mul(add(mul(add(mul(add(mul(add(mul(add(mul(v_0, -1/4), 1/2), v_0), -1), v_0), -12), v_0), 1/2), v_0), -1/95), v_0), -1/2)

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
        # We handle polynomial rings like QQ['x']['y']; that gives us some
        # slowdown.  Optimize away some of that:
        if len(etb._vars) == 1:
            # OK, we're in the (very common) univariate case.
            coeff_maker = etb.constant
        else:
            # There may be variables in our coefficients...
            coeff_maker = etb.make
        if d == -1:
            return coeff_maker(0)
        coeff = self.get_unsafe(d)
        if d == 0:
            return coeff_maker(coeff)
        if coeff != 1:
            expr *= coeff_maker(coeff)
        for i from d > i >= 0:
            coeff = self.get_unsafe(i)
            if coeff:
                expr += coeff_maker(coeff)
            if i > 0:
                expr *= x
        return expr

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare the two polynomials ``self`` and ``other``.

        We order polynomials first by degree (but treating 0 as having
        degree 0), then in dictionary order starting with the
        coefficient of largest degree.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: 3*x^3  + 5 > 10*x^2 + 19
            True
            sage: x^2 - 2*x - 1 < x^2 - 1
            True
            sage: x^2 - 2*x - 1 > x^2 - 1
            False
            sage: x^3 - 3 > 393939393
            True

        Test comparison with zero (:issue:`18633`)::

            sage: 0 < R(1)
            True
            sage: R(-1) < 0
            True
            sage: -x < 0
            False
            sage: R(0) == R(0)
            True

        TESTS:

        Test that comparisons are consistent when using interval
        coefficients::

            sage: # needs sage.rings.real_interval_field
            sage: R.<x> = RIF[]
            sage: a = RIF(0,1) * x
            sage: b = RIF(1,2) * x
            sage: a == a
            False
            sage: a != a
            False
            sage: a == b
            False
            sage: a < b
            False
            sage: a > b
            False
            sage: a <= b
            True
            sage: a >= b
            False
            sage: a != b
            False

            sage: R.<x> = RBF[]                                                         # needs sage.libs.flint
            sage: pol = RBF(1.0, 0.1)                                                   # needs sage.libs.flint
            sage: pol == pol                                                            # needs sage.libs.flint
            False
        """
        cdef Polynomial pol = <Polynomial?>other

        cdef Py_ssize_t d1 = self.degree()
        cdef Py_ssize_t d2 = pol.degree()

        # Special case constant polynomials
        if d1 == -1: # self is the 0 polynomial
            if d2 == -1:
                return rich_to_bool(op, 0) # both polynomials are 0
            elif d2 == 0:
                return richcmp(self._parent._base.zero(), pol.get_unsafe(0), op)
            return rich_to_bool_sgn(op, -1) # we have d2 > 0
        elif d1 == 0: # self is a nonzero constant
            if d2 == -1:
                return richcmp(self.get_unsafe(0), pol._parent._base.zero(), op)
            elif d2 == 0:
                return richcmp(self.get_unsafe(0), pol.get_unsafe(0), op)
            return rich_to_bool_sgn(op, -1) # we have d2 > d1 == 0

        # For different degrees, compare the degree
        if d1 != d2:
            return rich_to_bool_sgn(op, d1 - d2)

        cdef Py_ssize_t i
        for i in reversed(range(d1+1)):
            x = self.get_unsafe(i)
            y = pol.get_unsafe(i)
            res = richcmp_item(x, y, op)
            if res is not NotImplemented:
                return res
        return rich_to_bool(op, 0)

    def __bool__(self):
        """
        EXAMPLES::

            sage: P = PolynomialRing(ZZ, 'x')(0)
            sage: bool(P)
            False
            sage: P = PolynomialRing(ZZ, 'x')([1,2,3])
            sage: bool(P)
            True
        """
        return self.degree() >= 0

    def __getitem__(self, n):
        r"""
        Return the `n`-th coefficient of ``self``.

        .. WARNING::

            If `P` is a polynomial of degree `d`, then ``P[i]``
            returns `0` when `i < 0` or `i > d`.  This behaviour
            intentionally differs from that of lists: if `L` is a list
            of length `n`, then Python defines ``L[-i] = L[n - i]``
            for `0 < i \le n``.  The definition used here is more
            meaningful for polynomials, since it can be extended
            immediately to Laurent series, for example.

        EXAMPLES:

        We illustrate the difference between polynomials and lists
        when negative indices are involved::

            sage: R.<x> = QQ[]
            sage: f = x + 2
            sage: f[-1]
            0
            sage: list(f)[-1]
            1

        Slices can be used to truncate polynomials::

            sage: pol = R(range(8)); pol
            7*x^7 + 6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x
            sage: pol[:6]
            5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x

        Any other kind of slicing is an error, see :issue:`18940`::

            sage: f[1:3]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a start is not defined

            sage: f[1:3:2]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a step is not defined
        """
        cdef Py_ssize_t d = self.degree() + 1
        if isinstance(n, slice):
            start, stop, step = n.start, n.stop, n.step
            if step is not None:
                raise IndexError("polynomial slicing with a step is not defined")
            if start is not None:
                raise IndexError("polynomial slicing with a start is not defined")
            if stop is None or stop > d:
                stop = d
            values = [self.get_unsafe(i) for i in range(stop)]
            return self._new_generic(values)

        return self.get_coeff_c(pyobject_to_long(n))

    cdef get_coeff_c(self, Py_ssize_t i):
        """
        Return the `i`-th coefficient of ``self``.
        """
        cdef Py_ssize_t d = self.degree()
        if 0 <= i <= d:
            return self.get_unsafe(i)
        else:
            return self._parent._base.zero()

    cdef get_unsafe(self, Py_ssize_t i):
        """
        Return the `i`-th coefficient of ``self``.

        Used as building block for a generic :meth:`__getitem__`. Should be
        overridden by Cython subclasses. The default implementation makes it
        possible to implement concrete subclasses in Python.
        """
        return self[i]

    def __iter__(self):
        """
        EXAMPLES::

            sage: P = PolynomialRing(ZZ, 'x')([1,2,3])
            sage: [y for y in iter(P)]
            [1, 2, 3]
        """
        return iter(self.list(copy=False))

    def _cache_key(self):
        """
        Return a hashable key which identifies this element.

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: K.<u> = Qq(4)
            sage: R.<x> = K[]
            sage: f = x
            sage: hash(f)
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'
            sage: f._cache_key()
            (Univariate Polynomial Ring in x over 2-adic Unramified Extension Field in u defined by x^2 + x + 1,
             0,
             1 + O(2^20))
            sage: @cached_function
            ....: def foo(t): return t
            ....:
            sage: foo(x)
            (1 + O(2^20))*x
        """
        return (self._parent,) + tuple(self)

    def __hash__(self):
        return self._hash_c()

    cdef long _hash_c(self) except -1:
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

        TESTS:

        Verify that :issue:`16251` has been resolved, i.e., polynomials with
        unhashable coefficients are unhashable::

            sage: K.<a> = Qq(9)                                                         # needs sage.rings.padics
            sage: R.<t> = K[]                                                           # needs sage.rings.padics
            sage: hash(t)                                                               # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash
        cdef int i
        for i from 0<= i <= self.degree():
            if i == 1:
                # we delay the hashing until now to not waste it on a constant poly
                var_name_hash = hash(self._parent._names[0])
            # I'm assuming (incorrectly) that hashes of zero indicate that the element is 0.
            # This assumption is not true, but I think it is true enough for the purposes and it
            # it allows us to write fast code that omits terms with 0 coefficients.  This is
            # important if we want to maintain the '==' relationship with sparse polys.
            c_hash = hash(self.get_unsafe(i))
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

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this element under the morphism defined by
        ``im_gens`` in ``codomain``, where elements of the
        base ring are mapped by ``base_map``.

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
            sage: f(x^2 + 3)  # indirect doctest
            28

            sage: # needs sage.rings.number_field
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: S.<y> = K[]
            sage: phi = S.hom([y^2], base_map=cc)
            sage: phi(i*y)
            -i*y^2
        """
        a = im_gens[0]
        d = self.degree()
        if d == -1:
            return codomain.zero() # Special case: 0 should always coerce to 0
        if base_map is None:
            base_map = codomain.coerce_map_from(self.base_ring())
        result = base_map(self.get_unsafe(d))
        i = d - 1
        while i >= 0:
            result = result * a + base_map(self.get_unsafe(i))
            i -= 1
        return result

    def _scalar_conversion(self, R):
        r"""
        Generic conversion of constant polynomial to the ring ``R``.

        Ideally such conversions should go through
        :class:`ConstantPolynomialSection`. However, it does not work when
        there are additional conversions involved.

        EXAMPLES::

            sage: a = QQ['x'](1/5)
            sage: QQ(a)
            1/5
            sage: AA(a)                                                                 # needs sage.rings.number_field
            1/5
            sage: QQbar(a)                                                              # needs sage.rings.number_field
            1/5
            sage: RDF(a)
            0.2
            sage: CDF(a)                                                                # needs sage.rings.complex_double
            0.2
            sage: RR(a)                                                                 # needs sage.rings.real_mpfr
            0.200000000000000
            sage: CC(a)                                                                 # needs sage.rings.real_mpfr
            0.200000000000000
            sage: RBF(a)                                                                # needs sage.libs.flint
            [0.2000000000000000 +/- 4.45e-17]
            sage: CBF(a)                                                                # needs sage.libs.flint
            [0.2000000000000000 +/- 4.45e-17]
            sage: RIF(a)                                                                # needs sage.rings.real_interval_field
            0.2000000000000000?
            sage: CIF(a)                                                                # needs sage.rings.complex_interval_field
            0.2000000000000000?
            sage: float(a)
            0.2
            sage: complex(a)
            (0.2+0j)

            sage: # needs sage.rings.number_field
            sage: b = AA['x'](AA(2/3).sqrt())
            sage: AA(b)
            0.8164965809277260?
            sage: RR(b)                                                                 # needs sage.rings.real_mpfr
            0.816496580927726
            sage: RBF(b)                                                                # needs sage.libs.flint
            [0.816496580927726 +/- 2.44e-16]
            sage: RIF(b)                                                                # needs sage.rings.real_interval_field
            0.8164965809277260?
            sage: float(b)
            0.816496580927726

            sage: # needs sage.rings.number_field
            sage: c = QQbar['x'](QQbar(-2/5).sqrt())
            sage: QQbar(c)
            0.6324555320336758?*I
            sage: CDF(c)
            0.6324555320336758*I
            sage: CC(c)                                                                 # needs sage.rings.real_mpfr
            0.632455532033676*I
            sage: CBF(c)  # abs tol 1e-16                                               # needs sage.libs.flint
            [0.6324555320336759 +/- 3.38e-17]*I
            sage: CIF(c)                                                                # needs sage.rings.complex_interval_field
            0.6324555320336758?*I
            sage: complex(c)
            0.6324555320336758j

            sage: K.<x> = Frac(RR['x'])
            sage: ZZ(2*x/x)              # indirect doctest
            2
            sage: ZZ(x)
            Traceback (most recent call last):
            ...
            TypeError: cannot convert nonconstant polynomial

            sage: x = polygen(QQ)
            sage: A.<u> = NumberField(x^3 - 2)                                          # needs sage.rings.number_field
            sage: A(A['x'](u))                                                          # needs sage.rings.number_field
            u
        """
        if self.degree() > 0:
            raise TypeError("cannot convert nonconstant polynomial")
        return R(self.get_coeff_c(0))

    _real_double_ = _scalar_conversion
    _complex_double_ = _scalar_conversion
    _mpfr_ = _scalar_conversion
    _complex_mpfr_ = _scalar_conversion
    _real_mpfi_ = _scalar_conversion
    _complex_mpfi_ = _scalar_conversion
    _arb_ = _scalar_conversion
    _acb_ = _scalar_conversion
    _integer_ = _scalar_conversion
    _algebraic_ = _scalar_conversion

    def __int__(self):
        """
        EXAMPLES::

            sage: P = PolynomialRing(ZZ, 'x')([3])
            sage: int(P)
            3
        """
        return self._scalar_conversion(int)

    def __float__(self):
        """
        EXAMPLES::

            sage: P = PolynomialRing(ZZ, 'x')([1])
            sage: float(P)
            1.0
        """
        return self._scalar_conversion(float)

    def __complex__(self):
        r"""
        EXAMPLES::

            sage: p = PolynomialRing(QQbar, 'x')(1+I)                                   # needs sage.rings.number_field
            sage: complex(p)                                                            # needs sage.rings.number_field
            (1+1j)
        """
        return self._scalar_conversion(complex)

    def _rational_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: QQ(R(45/4))
            45/4
            sage: QQ(3*x + 45)
            Traceback (most recent call last):
            ...
            TypeError: 3*x + 45 is not a constant polynomial
        """
        return self._scalar_conversion(sage.rings.rational.Rational)

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R.<x> = QQ[]
            sage: f = x^3 + x
            sage: g = f._symbolic_(SR); g
            x^3 + x
            sage: g(x=2)
            10
            sage: g = SR(f)
            sage: g(x=2)
            10

        The polynomial has to be over a field of characteristic 0 (see
        :issue:`24072`)::

            sage: R.<w> = GF(7)[]
            sage: f = SR(2*w^3 + 1); f                                                  # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations
        """
        d = dict([(repr(g), R.var(g)) for g in self._parent.gens()])
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
        return self.parent().one() / self

    def inverse_of_unit(self):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x - 90283
            sage: f.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: x - 90283 is not a unit
            in Univariate Polynomial Ring in x over Rational Field
            sage: f = R(-90283); g = f.inverse_of_unit(); g
            -1/90283
            sage: parent(g)
            Univariate Polynomial Ring in x over Rational Field

        TESTS::

            sage: Integers(1)['x'](0).inverse_of_unit()
            0
        """
        d = self.degree()
        if d > 0:
            if not self.is_unit():
                raise ArithmeticError(f"{self} is not a unit in {self.parent()}")
            else:
                raise NotImplementedError("polynomial inversion over non-integral domains not implemented")
        elif d == -1:
            cst = self._parent._base.zero()
        else:
            cst = self.get_unsafe(0)
        inv = cst.inverse_of_unit()
        return self._parent([inv])

    def inverse_mod(a, m):
        """
        Invert the polynomial ``a`` with respect to ``m``, or raise a
        :exc:`ValueError` if no such inverse exists.

        The parameter ``m`` may be either a single polynomial or an ideal
        (for consistency with :meth:`inverse_mod` in other rings).

        .. SEEALSO::

            If you are only interested in the inverse modulo a monomial `x^k`
            then you might use the specialized method
            :meth:`inverse_series_trunc` which is much faster.

        EXAMPLES::

            sage: S.<t> = QQ[]
            sage: f = inverse_mod(t^2 + 1, t^3 + 1); f
            -1/2*t^2 - 1/2*t + 1/2
            sage: f * (t^2 + 1) % (t^3 + 1)
            1
            sage: f = t.inverse_mod((t + 1)^7); f
            -t^6 - 7*t^5 - 21*t^4 - 35*t^3 - 35*t^2 - 21*t - 7
            sage: (f * t) + (t + 1)^7
            1
            sage: t.inverse_mod(S.ideal((t + 1)^7)) == f
            True

        This also works over inexact rings, but note that due to rounding
        error the product may not always exactly equal the constant
        polynomial 1 and have extra terms with coefficients close to zero. ::

            sage: # needs scipy sage.modules
            sage: R.<x> = RDF[]
            sage: epsilon = RDF(1).ulp()*50   # Allow an error of up to 50 ulp
            sage: f = inverse_mod(x^2 + 1, x^5 + x + 1); f  # abs tol 1e-14
            0.4*x^4 - 0.2*x^3 - 0.4*x^2 + 0.2*x + 0.8
            sage: poly = f * (x^2 + 1) % (x^5 + x + 1)
            sage: # Remove noisy zero terms:
            sage: parent(poly)([0.0 if abs(c) <= epsilon else c
            ....:               for c in poly.coefficients(sparse=False)])
            1.0
            sage: f = inverse_mod(x^3 - x + 1, x - 2); f
            0.14285714285714285
            sage: f * (x^3 - x + 1) % (x - 2)
            1.0
            sage: g = 5*x^3 + x - 7; m = x^4 - 12*x + 13; f = inverse_mod(g, m); f
            -0.0319636125...*x^3 - 0.0383269759...*x^2 - 0.0463050900...*x + 0.346479687...
            sage: poly = f*g % m
            sage: # Remove noisy zero terms:
            sage: parent(poly)([0.0 if abs(c) <= epsilon else c  # abs tol 1e-14
            ....:               for c in poly.coefficients(sparse=False)])
            1.0000000000000004

        ALGORITHM: Solve the system `as + mt = 1`, returning `s` as the inverse
        of `a` mod `m`.

        Uses the Euclidean algorithm for exact rings, and solves a linear
        system for the coefficients of `s` and `t` for inexact rings (as the
        Euclidean algorithm may not converge in that case).

        AUTHORS:

        - Robert Bradshaw (2007-05-31)
        """
        from sage.rings.ideal import Ideal_generic
        if isinstance(m, Ideal_generic):
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
            v = vector(R, [R.one()] + [R.zero()]*(2*n-2)) # the constant polynomial 1
            if M.is_invertible():
                x = M.solve_right(v) # there has to be a better way to solve
                return a.parent()(list(x)[0:n])
            else:
                raise ValueError("Impossible inverse modulo")

    cpdef Polynomial inverse_series_trunc(self, long prec):
        r"""
        Return a polynomial approximation of precision ``prec`` of the inverse
        series of this polynomial.

        .. SEEALSO::

            The method :meth:`inverse_mod` allows more generally to invert this
            polynomial with respect to any ideal.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: s = (1 + x).inverse_series_trunc(5)
            sage: s
            x^4 - x^3 + x^2 - x + 1
            sage: s * (1 + x)
            x^5 + 1

        Note that the constant coefficient needs to be a unit::

            sage: ZZx.<x> = ZZ[]
            sage: ZZxy.<y> = ZZx[]
            sage: (1+x + y**2).inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ValueError: constant term x + 1 is not a unit
            sage: (1+x + y**2).change_ring(ZZx.fraction_field()).inverse_series_trunc(4)
            (-1/(x^2 + 2*x + 1))*y^2 + 1/(x + 1)

        The method works over any polynomial ring::

            sage: R = Zmod(4)
            sage: Rx.<x> = R[]
            sage: Rxy.<y> = Rx[]

            sage: p = 1 + (1+2*x)*y + x**2*y**4
            sage: q = p.inverse_series_trunc(10)
            sage: (p*q).truncate(11)
            (2*x^4 + 3*x^2 + 3)*y^10 + 1

        Even noncommutative ones::

            sage: # needs sage.modules
            sage: M = MatrixSpace(ZZ, 2)
            sage: x = polygen(M)
            sage: p = M([1,2,3,4])*x^3 + M([-1,0,0,1])*x^2 + M([1,3,-1,0])*x + M.one()
            sage: q = p.inverse_series_trunc(5)
            sage: (p*q).truncate(5) == M.one()
            True
            sage: q = p.inverse_series_trunc(13)
            sage: (p*q).truncate(13) == M.one()
            True

        TESTS::

            sage: x = polygen(ZZ['a','b'])
            sage: (x + 1).inverse_series_trunc(0)
            Traceback (most recent call last):
            ...
            ValueError: the precision must be positive, got 0

        AUTHORS:

        - David Harvey (2006-09-09): Newton's method implementation for power
          series

        - Vincent Delecroix (2014-2015): move the implementation directly in
          polynomial
        """
        if prec <= 0:
            raise ValueError("the precision must be positive, got {}".format(prec))

        const_term = self.get_coeff_c(0)
        if not const_term.is_unit():
            raise ValueError("constant term {} is not a unit".format(const_term))

        R = self._parent
        A = R.base_ring()
        try:
            first_coeff = const_term.inverse_of_unit()
        except AttributeError:
            first_coeff = A(~const_term)

        current = R(first_coeff)
        for next_prec in sage.misc.misc.newton_method_sizes(prec)[1:]:
            z = current._mul_trunc_(self, next_prec)._mul_trunc_(current, next_prec)
            current = current + current - z
        return current

    def revert_series(self, n):
        r"""
        Return a polynomial ``f`` such that
        ``f(self(x)) = self(f(x)) = x mod x^n``.

        Currently, this is only implemented over some coefficient rings.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: (x + x^3/6 + x^5/120).revert_series(6)
            3/40*x^5 - 1/6*x^3 + x
            sage: Pol.<x> = CBF[]                                                       # needs sage.libs.flint
            sage: (x + x^3/6 + x^5/120).revert_series(6)                                # needs sage.libs.flint
            ([0.075000000000000 +/- ...e-17])*x^5 + ([-0.166666666666667 +/- ...e-16])*x^3 + x

            sage: # needs sage.symbolic
            sage: Pol.<x> = SR[]
            sage: x.revert_series(6)
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for certain base rings
        """
        raise NotImplementedError("only implemented for certain base rings")

    cpdef _mul_(self, right):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x - 4) * (x^2 - 8*x + 16)
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

        TESTS::

            sage: # needs sage.modules
            sage: Pol.<x> = MatrixSpace(ZZ, 2)[]
            sage: a = matrix([[1,0], [0,0]])
            sage: b = matrix([[1,2], [3,4]])
            sage: list((a*x)*(b*x + 1))
            [
            [0 0]  [1 0]  [1 2]
            [0 0], [0 0], [0 0]
            ]
            sage: list((b*x + 1)*(a*x))
            [
            [0 0]  [1 0]  [1 0]
            [0 0], [0 0], [3 0]
            ]
            sage: list((a*x + 1)*(b*x))
            [
            [0 0]  [1 2]  [1 2]
            [0 0], [3 4], [0 0]
            ]
        """
        if not self or not right:
            return self._parent.zero()

        cdef Polynomial _right = right
        if self.is_term():
            return _right._mul_term(self, term_on_right=False)
        elif _right.is_term():
            return self._mul_term(_right, term_on_right=True)

        elif self._parent.is_exact():
            return self._mul_karatsuba(right)
        else:
            return self._mul_generic(right)

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n):
        r"""
        Return the truncated multiplication of two polynomials up to ``n``.

        This is the default implementation that does the multiplication and then
        truncate! There are custom implementations in several subclasses:

        - :meth:`on dense polynomial over integers (via FLINT) <sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint._mul_trunc_>`

        - :meth:`on dense polynomial over Z/nZ (via FLINT)
          <sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint._mul_trunc_>`

        - :meth:`on dense rational polynomial (via FLINT)
          <sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint._mul_trunc_>`

        - :meth:`on dense polynomial on Z/nZ (via NTL)
          <sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_zz._mul_trunc_>`

        EXAMPLES::

            sage: R = QQ['x']['y']
            sage: y = R.gen()
            sage: x = R.base_ring().gen()
            sage: p1 = 1 - x*y + 2*y**3
            sage: p2 = -1/3 + y**5
            sage: p1._mul_trunc_(p2, 5)
            -2/3*y^3 + 1/3*x*y - 1/3

        .. TODO::

            implement a generic truncated Karatsuba and use it here.
        """
        cdef Polynomial pol
        cdef list x, y
        if not self or not right:
            return self._parent.zero()
        elif n < self._parent._Karatsuba_threshold:
            x = self.list(copy=False)
            y = right.list(copy=False)
            return self._new_generic(do_schoolbook_product(x, y, n))
        else:
            pol = self.truncate(n) * right.truncate(n)
            return pol._inplace_truncate(n)

    def multiplication_trunc(self, other, n):
        r"""
        Truncated multiplication.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x^10 + 5*x^5 + x^2 - 3).multiplication_trunc(x^7 - 3*x^3 + 1, 11)
            x^10 + x^9 - 15*x^8 - 3*x^7 + 2*x^5 + 9*x^3 + x^2 - 3

        Check that coercion is working::

            sage: R2 = QQ['x']
            sage: x2 = R2.gen()
            sage: p1 = (x^3 + 1).multiplication_trunc(x2^3 - 2, 5); p1
            -x^3 - 2
            sage: p2 = (x2^3 + 1).multiplication_trunc(x^3 - 2, 5); p2
            -x^3 - 2
            sage: parent(p1) == parent(p2) == R2
            True
        """
        if not have_same_parent(self, other):
            self, other = coercion_model.canonical_coercion(self, other)
        return self._mul_trunc_(other, pyobject_to_long(n))

    def square(self):
        """
        Return the square of this polynomial.

        .. TODO::

            - This is just a placeholder; for now it just uses ordinary
              multiplication. But generally speaking, squaring is faster than
              ordinary multiplication, and it's frequently used, so subclasses
              may choose to provide a specialised squaring routine.

            - Perhaps this even belongs at a lower level? RingElement or
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
        Return the square-free decomposition of this polynomial.  This is a
        partial factorization into square-free, coprime polynomials.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: p = 37 * (x - 1)^3 * (x - 2)^3 * (x - 1/3)^7 * (x - 3/7)
            sage: p.squarefree_decomposition()
            (37*x - 111/7) * (x^2 - 3*x + 2)^3 * (x - 1/3)^7
            sage: p = 37 * (x - 2/3)^2
            sage: p.squarefree_decomposition()
            (37) * (x - 2/3)^2
            sage: x = polygen(GF(3))
            sage: x.squarefree_decomposition()
            x
            sage: f = QQbar['x'](1)                                                     # needs sage.rings.number_field
            sage: f.squarefree_decomposition()                                          # needs sage.rings.number_field
            1
        """
        if self.degree() < 0:
            raise ValueError("square-free decomposition not defined for zero polynomial")
        if hasattr(self.base_ring(),'_squarefree_decomposition_univariate_polynomial'):
            return self.base_ring()._squarefree_decomposition_univariate_polynomial(self)
        raise NotImplementedError("square-free decomposition not implemented for this polynomial")

    def is_square(self, root=False):
        """
        Return whether or not polynomial is square.

        If the optional
        argument ``root`` is set to ``True``, then also returns the square root
        (or ``None``, if the polynomial is not square).

        INPUT:

        - ``root`` -- whether or not to also return a square
          root (default: ``False``)

        OUTPUT:

        - boolean; whether or not a square

        - (optional) an actual square root if found, and ``None`` otherwise

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: (x^2 + 2*x + 1).is_square()
            True
            sage: (x^4 + 2*x^3 - x^2 - 2*x + 1).is_square(root=True)
            (True, x^2 + x - 1)

            sage: f = 12 * (x + 1)^2 * (x + 3)^2
            sage: f.is_square()
            False
            sage: f.is_square(root=True)
            (False, None)

            sage: h = f/3; h
            4*x^4 + 32*x^3 + 88*x^2 + 96*x + 36
            sage: h.is_square(root=True)
            (True, 2*x^2 + 8*x + 6)

            sage: S.<y> = PolynomialRing(RR)
            sage: g = 12 * (y + 1)^2 * (y + 3)^2

            sage: g.is_square()
            True

        TESTS:

        Make sure :issue:`9093` is fixed::

            sage: R(1).is_square()
            True
            sage: R(4/9).is_square(root=True)
            (True, 2/3)
            sage: R(-1/3).is_square()
            False
            sage: R(0).is_square()
            True

        Make sure :issue:`35860` is fixed::

            sage: S.<x> = PolynomialRing(ZZ)
            sage: is_square(S(1), True)[1].parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        if self.is_zero():
            return (True, self) if root else True

        try:
            f = self.squarefree_decomposition()
        except NotImplementedError:
            f = self.factor()

        u = self._parent.base_ring()(f.unit())

        if all(a[1] % 2 == 0 for a in f) and u.is_square():
            g = self._parent(u.sqrt())
            for a in f:
                g *= a[0] ** (a[1] // 2)
            return (True, g) if root else True
        else:
            return (False, None) if root else False

    def _distinct_degree_factorisation_squarefree(self):
        """
        Helper function for any_irreducible_factor which computes
        the distinct degree factorisation of ``self``.

        Creates an iterator for all valid degrees `d`, and returns
        tuples of the form `(a_d, d)` for a polynomial `a_d` the
        product of irreducible polynomials of degree `d`.

        Assumes that this polynomial is squarefree.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: f = (x + 162) * (x^3 + 7*x + 161) * (x^7 + 9*x + 161)
            sage: list(f._distinct_degree_factorisation_squarefree())
            [(x + 162, 1), (x^3 + 7*x + 161, 3), (x^7 + 9*x + 161, 7)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<z8> = GF(2^8)
            sage: R.<x> = K[]
            sage: u1 = x + z8^6 + z8^3 + z8
            sage: u2 = x^2 + (z8^5 + z8^3 + z8^2 + 1)*x + z8^4 + z8^3 + z8 + 1
            sage: u3 = x^3 + z8^7 + z8^4 + z8^3 + z8^2 + z8
            sage: f = u1 * u2 * u3
            sage: list(f._distinct_degree_factorisation_squarefree())
            [(x + z8^6 + z8^3 + z8, 1),
            (x^2 + (z8^5 + z8^3 + z8^2 + 1)*x + z8^4 + z8^3 + z8 + 1, 2),
            (x^3 + z8^7 + z8^4 + z8^3 + z8^2 + z8, 3)]
        """
        q = self.base_ring().order() # p^k
        R, x = self.parent().objgen()

        # Initialise values
        v = self
        w = x
        d = 0
        e = v.degree()

        # Iterate over all possible degrees with degree
        while 2*(d + 1) <= e:
            d = d + 1
            w = pow(w, q, v)

            ad = v.gcd(w - x)
            if not ad.is_one():
                yield (ad, d)
                v = v // ad

            e = v.degree()

        # Last case, v itself might be irreducible
        if e > 0:
            yield (v, e)
        return

    def _cantor_zassenhaus_split_to_irreducible(self, degree):
        """
        Helper function for :meth:`any_irreducible_factor` which computes
        a factor from a polynomial of the form `self = prod g_i(x)`
        with all `g_i(x)` irreducible of the same degree. Uses the
        Cantor-Zassenhaus splitting method.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: f = (x^2 + 36*x + 34) * (x^2 + 36*x + 132) * (x^2 + 112*x + 48)
            sage: f._cantor_zassenhaus_split_to_irreducible(2)   # random
            x^2 + 112*x + 48
            sage: f._cantor_zassenhaus_split_to_irreducible(2)   # random
            x^2 + 36*x + 34
            sage: f._cantor_zassenhaus_split_to_irreducible(2)   # random
            x^2 + 36*x + 132
            sage: factor = f._cantor_zassenhaus_split_to_irreducible(2)
            sage: factor.is_irreducible()
            True
            sage: f % factor == 0
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = F[]
            sage: f = (x + z4^2 + z4) * (x + z4^2 + z4 + 1)
            sage: f._cantor_zassenhaus_split_to_irreducible(1)   # random
            x + z4^2 + z4
            sage: f._cantor_zassenhaus_split_to_irreducible(1)   # random
            x + z4^2 + z4
            sage: f._cantor_zassenhaus_split_to_irreducible(1)   # random
            x + z4^2 + z4 + 1
            sage: factor = f._cantor_zassenhaus_split_to_irreducible(1)
            sage: factor.is_irreducible()
            True
            sage: f % factor == 0
            True

        TESTS:

        Ensure that :issue:`37445` is fixed::

            sage: R.<x> = GF(13)[]
            sage: def irr(d, R): return f.monic() if (f := R.random_element(d)).is_irreducible() else irr(d, R)
            sage: f = prod(irr(6, R) for _ in range(10))
            sage: irr = f._cantor_zassenhaus_split_to_irreducible(6)
            sage: assert irr.degree() == 6
            sage: assert f % irr == 0
            sage: assert irr.is_irreducible()
        """
        R = self.parent()
        q = self.base_ring().order()

        # Polynomial is already irreducible by the assumptions of the function
        if self.degree() == degree:
            return self

        # We expect to succeed with greater than 1/2 probability,
        # so if we try 100 times and fail, there's a bug somewhere.
        for _ in range(100):
            # Sample a polynomial "uniformly" from R
            # TODO: once #37118 has been merged, this can be made cleaner,
            #       as we will actually have access to uniform sampling.
            #       At the moment, we make an ugly call for polynomials of
            #       degree exactly 2*degree + 1, instead of polynomials of degree
            #       less than `self.degree()` as the original desc. of C-Z.
            #       Additionally, we can use the `monic` flag:
            #       R.random_element(degree=(self.degree() - 1), monic=True)
            T = R.random_element(2*degree + 1).monic()

            # Need to handle odd and even characteristic separately
            if q % 2:
                h = self.gcd(pow(T, (q**degree-1)//2, self) - 1)
            else:
                # Compute the trace of T with field of order 2^k
                # sum T^(2^i) for i in range (degree * k)
                # We use repeated squaring to avoid redundent multiplications
                C, TT = T, T
                for _ in range(degree * self.base_ring().degree() - 1):
                    TT = TT * TT % self
                    C += TT
                h = self.gcd(C)

            hd = h.degree()

            # If we found a factor of desired degree, return it
            if hd == degree:
                return h

            # Else check if we have a non-trivial factor and keep going
            if 0 < hd < self.degree():
                if 2*hd > self.degree():
                    h = self // h
                return h._cantor_zassenhaus_split_to_irreducible(degree)

        # If you are reaching this error, chances are there's a bug in the code.
        raise AssertionError(f"no splitting of degree {degree} found for {self}")

    def _any_irreducible_factor_squarefree(self, degree=None, ext_degree=None):
        """
        Helper function for :meth:`any_irreducible_factor` which computes
        an irreducible factor from ``self``, assuming the input is
        squarefree.

        Does this by first computing the distinct degree factorisations
        of ``self`` and then finds a factor with Cantor-Zassenhaus
        splitting.

        If ``degree`` is not ``None``, then only irreducible factors of degree
        ``degree`` are searched for, otherwise the smallest degree factor
        is found.

        If ``ext_degree`` is not ``None``, then only irreducible factors whose
        degree divides ``ext_degree`` are returned.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: f = (x^2 + 36*x + 34) * (x^2 + 36*x + 132) * (x^2 + 112*x + 48)
            sage: f = (x^2 + 50*x + 28) * (x^5 + 75*x^4 + 87*x^3 + 96*x^2 + 98*x + 40) * (x^7 + 20*x^6 + 16*x^5 + 36*x^4 + 59*x^3 + 3*x^2 + 46*x + 84)
            sage: f._any_irreducible_factor_squarefree()
            x^2 + 50*x + 28
            sage: f._any_irreducible_factor_squarefree(degree=2)
            x^2 + 50*x + 28
            sage: f._any_irreducible_factor_squarefree(degree=5)
            x^5 + 75*x^4 + 87*x^3 + 96*x^2 + 98*x + 40
            sage: f._any_irreducible_factor_squarefree(degree=7)
            x^7 + 20*x^6 + 16*x^5 + 36*x^4 + 59*x^3 + 3*x^2 + 46*x + 84
            sage: g = (x^2 + 50*x + 28)
            sage: g._any_irreducible_factor_squarefree(degree=1)
            Traceback (most recent call last):
            ...
            ValueError: no irreducible factor of degree 1 could be computed from x^2 + 50*x + 28

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = F[]
            sage: f = (x + z4^3 + z4^2 + z4) * (x^2 + x + z4^3 + 1) * (x^3 + (z4^3 + z4)*x^2 + z4^2 + 1)
            sage: f._any_irreducible_factor_squarefree()
            x + z4^3 + z4^2 + z4
            sage: f._any_irreducible_factor_squarefree(degree=1)
            x + z4^3 + z4^2 + z4
            sage: f._any_irreducible_factor_squarefree(degree=2)
            x^2 + x + z4^3 + 1
            sage: f._any_irreducible_factor_squarefree(degree=3)
            x^3 + (z4^3 + z4)*x^2 + z4^2 + 1
        """
        # If the degree is not None we only want to check a single polynomial
        if degree is not None:
            for (poly, d) in self._distinct_degree_factorisation_squarefree():
                if d == degree:
                    return poly._cantor_zassenhaus_split_to_irreducible(degree)
                # Stop iterating early if the degree is too large
                elif d > degree:
                    raise ValueError(f"no irreducible factor of degree {degree} could be computed from {self}")
            raise ValueError(f"no irreducible factor of degree {degree} could be computed from {self}")

        # Otherwise we use the smallest possible d value
        for (poly, d) in self._distinct_degree_factorisation_squarefree():
            if ext_degree is None:
                return poly._cantor_zassenhaus_split_to_irreducible(d)
            elif ZZ(d).divides(ext_degree):
                return poly._cantor_zassenhaus_split_to_irreducible(d)
            if d > ext_degree:
                raise ValueError(f"no irreducible factor of degree {degree} dividing {ext_degree} could be computed from {self}")
        raise AssertionError(f"no irreducible factor could be computed from {self}")

    def any_irreducible_factor(self, degree=None, assume_squarefree=False, assume_equal_deg=False, ext_degree=None):
        """
        Return an irreducible factor of this polynomial.

        INPUT:

        - ``degree`` -- ``None`` or positive integer (default: ``None``).
          Used for polynomials over finite fields. If ``None``, returns
          the the first factor found (usually the smallest). Otherwise,
          attempts to return an irreducible factor of ``self`` of chosen
          degree ``degree``.

        - ``assume_squarefree`` -- boolean (default: ``False``);
          Used for polynomials over finite fields.  If ``True``,
          this polynomial is assumed to be squarefree.

        - ``assume_equal_deg`` -- boolean (default: ``False``);
          Used for polynomials over finite fields.  If ``True``,
          this polynomial is assumed to be the product of irreducible
          polynomials of degree equal to ``degree``.

        - ``ext_degree`` -- positive integer or ``None`` (default);
          used for polynomials over finite fields. If not ``None`` only returns
          irreducible factors of ``self`` whose degree divides ``ext_degree``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: f = (x + 40)^3 * (x^5 + 76*x^4 + 93*x^3 + 112*x^2 + 22*x + 27)^2 * (x^6 + 50*x^5 + 143*x^4 + 162*x^2 + 109*x + 140)
            sage: f.any_irreducible_factor()
            x + 40
            sage: f.any_irreducible_factor(degree=5)
            x^5 + 76*x^4 + 93*x^3 + 112*x^2 + 22*x + 27

        When the polynomial is known to be squarefree we can optimise the call
        by setting ``assume_squarefree`` to be ``True``::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: g = (x - 1) * (x^3 + 7*x + 161)
            sage: g.any_irreducible_factor(assume_squarefree=True)
            x + 162
            sage: g.any_irreducible_factor(degree=3, assume_squarefree=True)
            x^3 + 7*x + 161

        If we ask for an irreducible factor which does not exist, the function
        will throw a :exc:`ValueError`::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: g = (x - 1) * (x^3 + 7*x + 161)
            sage: g.any_irreducible_factor(degree=2, assume_squarefree=True)
            Traceback (most recent call last):
            ...
            ValueError: no irreducible factor of degree 2 could be computed from x^4 + 162*x^3 + 7*x^2 + 154*x + 2

        If we assume that the polynomial is product of irreducible polynomials of the
        same degree, we must also supply the degree::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(163)
            sage: R.<x> = F[]
            sage: h = (x + 57) * (x + 98) * (x + 117) * (x + 145)
            sage: h.any_irreducible_factor(degree=1, assume_equal_deg=True)   # random
            x + 98
            sage: h.any_irreducible_factor(assume_equal_deg=True)
            Traceback (most recent call last):
            ...
            ValueError: degree must be known if distinct degree factorisation is assumed

        Also works for extension fields and even characteristic::

            sage: F.<z4> = GF(2^4)
            sage: R.<x> = F[]
            sage: f = (x + z4^3 + z4^2)^4 * (x^2 + z4*x + z4) * (x^2 + (z4^3 + z4^2 + z4)*x + z4^2 + z4 + 1)
            sage: f.any_irreducible_factor()
            x + z4^3 + z4^2
            sage: f.any_irreducible_factor(degree=2)  # random
            x^2 + (z4^3 + z4^2 + z4)*x + z4^2 + z4 + 1

        We can also use this function for polynomials which are not defined over finite
        fields, but this simply falls back to a slow method of factorisation::

            sage: R.<x> = ZZ[]
            sage: f = 3*x^4 + 2*x^3
            sage: f.any_irreducible_factor()
            3*x + 2
        """
        # Make sure the user inputted something reasonable for degree
        if degree is not None:
            degree = ZZ(degree)
            if degree < 1:
                raise ValueError(f"{degree = } must be positive")

        if assume_equal_deg and degree is None:
            raise ValueError("degree must be known if distinct degree factorisation is assumed")

        # When not working over a finite field, do the simple thing of factoring.
        # If degree is none, we return the first factor, otherwise we iterate
        # through and look for a polynomial with the desired degree.
        from sage.categories.finite_fields import FiniteFields
        if self.base_ring() not in FiniteFields():
            try:
                factorisation = self.factor()
            except (NotImplementedError, ValueError):
                raise ValueError(f"cannot factor {self} over the base ring {self.base_ring()}")
            if degree is None:
                return factorisation[0][0]
            for (poly, e) in factorisation:
                if poly.degree() == degree:
                    return poly
            raise ValueError(f"polynomial {self} has no irreducible factor of degree {degree}")

        # For finite fields, we find irreducible factors in the following three steps:
        #
        # 1. Compute the squarefree decomposition of the polynomial ``self``
        # 2. For each squarefree polynomial find the distinct degree `d`
        #    factorisation, F, which is the product of degree `d` polynomials
        #    dividing the squarefree polynomial
        # 3. Using Cantor-Zassenhaus splitting with degree `d` to find a
        #    single linear factor and return the root.
        #
        # When degree is None, we check all degrees smaller than the degree of the
        # squarefree polynomial, otherwise we work with only a single degree set by
        # the user.

        # Initial checks for bad input
        if self.degree() <= 0:
            raise ValueError(f"there are no irreducible factors of {self}")

        # If we know the polynomial is square-free, we can start here
        if assume_squarefree:
            if assume_equal_deg:
                return self._cantor_zassenhaus_split_to_irreducible(degree)
            return self._any_irreducible_factor_squarefree(degree, ext_degree)

        # Otherwise we compute the squarefree decomposition and check each
        # polynomial for a root. If no poly has a root, we raise an error.
        SFD = self.squarefree_decomposition()
        SFD.sort()
        for poly, _ in SFD:
            try:
                return poly._any_irreducible_factor_squarefree(degree, ext_degree)
            except ValueError:
                pass

        # If degree has been set, there could just be no factor of the desired degree
        if degree:
            raise ValueError(f"polynomial {self} has no irreducible factor of degree {degree}")
        # If ext_degree has been set, then there may be no irreducible factor of degree dividing ext_degree
        if ext_degree:
            raise ValueError(f"polynomial {self} has no irreducible factor of degree dividing {ext_degree}")
        # But if any degree is allowed then there should certainly be a factor if self has degree > 0
        raise AssertionError(f"no irreducible factor was computed for {self}. Bug.")

    def any_root(self, ring=None, degree=None, assume_squarefree=False, assume_equal_deg=False):
        """
        Return a root of this polynomial in the given ring.

        INPUT:

        - ``ring`` -- the ring in which a root is sought; by default
          this is the coefficient ring

        - ``degree`` -- ``None`` or nonzero integer; used for polynomials
          over finite fields.  Return a root of degree
          ``abs(degree)`` over the ground field.  If negative, also
          assumes that all factors of this polynomial are of degree
          ``abs(degree)``.  If ``None``, returns a root of minimal degree
          contained within the given ring.

        - ``assume_squarefree`` -- boolean; used for polynomials over
          finite fields.  If ``True``, this polynomial is assumed to be
          squarefree.

        - ``assume_equal_deg`` -- boolean; used for polynomials over
          finite fields.  If ``True``, all factors of this polynomial
          are assumed to have degree ``degree``. Note that ``degree``
          must be set.

        .. WARNING::

            Negative degree input will be deprecated. Instead use
            ``assume_equal_deg``.

        .. NOTE::

            For finite fields, ``any_root()`` is non-deterministic when
            finding linear roots of a polynomial over the base ring.
            However, if ``degree`` is greater than one, or ``ring`` is an
            extension of the base ring, then the root computed is found
            by attempting to return a root after factorisation. Roots found
            in this way are deterministic. This may change in the future.
            For all other rings or fields, roots are found by first
            fully-factoring ``self`` and the output is deterministic.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(11)[]
            sage: f = 7*x^7 + 8*x^6 + 4*x^5 + x^4 + 6*x^3 + 10*x^2 + 8*x + 5
            sage: f.any_root()
            2
            sage: f.factor()
            (7) * (x + 9) * (x^6 + 10*x^4 + 6*x^3 + 5*x^2 + 2*x + 2)
            sage: f = x^6 + 10*x^4 + 6*x^3 + 5*x^2 + 2*x + 2
            sage: root = f.any_root(GF(11^6, 'a'))
            sage: roots = sorted(f.roots(GF(11^6, 'a'), multiplicities=False))
            sage: roots
            [10*a^5 + 2*a^4 + 8*a^3 + 9*a^2 + a,
            a^5 + a^4 + 7*a^3 + 2*a^2 + 10*a,
            9*a^5 + 5*a^4 + 10*a^3 + 8*a^2 + 3*a + 1,
            2*a^5 + 8*a^4 + 3*a^3 + 6*a + 2,
            a^5 + 3*a^4 + 8*a^3 + 2*a^2 + 3*a + 4,
            10*a^5 + 3*a^4 + 8*a^3 + a^2 + 10*a + 4]
            sage: root in roots
            True

            sage: # needs sage.rings.finite_rings
            sage: g = (x-1) * (x^2 + 3*x + 9) * (x^5 + 5*x^4 + 8*x^3 + 5*x^2 + 3*x + 5)
            sage: g.any_root(ring=GF(11^10, 'b'), degree=1)
            1
            sage: root = g.any_root(ring=GF(11^10, 'b'), degree=2)
            sage: roots = (x^2 + 3*x + 9).roots(ring=GF(11^10, 'b'), multiplicities=False)
            sage: root in roots
            True
            sage: root = g.any_root(ring=GF(11^10, 'b'), degree=5)
            sage: roots = (x^5 + 5*x^4 + 8*x^3 + 5*x^2 + 3*x + 5).roots(ring=GF(11^10, 'b'), multiplicities=False)
            sage: root in roots
            True

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: K.<a> = GF(5^12)                                                      # needs sage.rings.finite_rings
            sage: for _ in range(40):                                                   # needs sage.rings.finite_rings
            ....:     f = R.random_element(degree=4)
            ....:     assert f(f.any_root(K)) == 0

        Check that our Cantor-Zassenhaus implementation does not loop
        over finite fields of even characteristic (see :issue:`16162`)::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(2**8)
            sage: x = polygen(K)
            sage: r = (x**2+x+1).any_root()  # used to loop
            sage: r**2 + r
            1
            sage: (x**2+a+1).any_root()
            a^7 + a^2

        Also check that such computations can be interrupted::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(2^8)
            sage: x = polygen(K)
            sage: pol = x^1000000 + x + a
            sage: alarm(0.5); pol.any_root()
            Traceback (most recent call last):
            ...
            AlarmInterrupt

        Check root computation over large finite fields::

            sage: # needs sage.libs.m4ri sage.rings.finite_rings
            sage: K.<a> = GF(2**50)
            sage: x = polygen(K)
            sage: (x**10+x+a).any_root()
            a^49 + a^47 + a^44 + a^42 + a^41 + a^39 + a^38 + a^37 + a^36
            + a^34 + a^33 + a^29 + a^27 + a^26 + a^25 + a^23 + a^18
            + a^13 + a^7 + a^5 + a^4 + a^3 + a^2 + a
            sage: K.<a> = GF(2**150)
            sage: x = polygen(K)
            sage: (x**10+x+a).any_root()
            a^149 + a^148 + a^146 + a^144 + a^143 + a^140 + a^138 + a^136 + a^134
            + a^132 + a^131 + a^130 + a^129 + a^127 + a^123 + a^120 + a^118 + a^114
            + a^113 + a^112 + a^111 + a^108 + a^104 + a^103 + a^102 + a^99 + a^98
            + a^94 + a^91 + a^90 + a^88 + a^79 + a^78 + a^75 + a^73 + a^72 + a^67
            + a^65 + a^64 + a^63 + a^62 + a^61 + a^59 + a^57 + a^52 + a^50 + a^48
            + a^47 + a^46 + a^45 + a^43 + a^41 + a^39 + a^37 + a^34 + a^31 + a^29
            + a^27 + a^25 + a^23 + a^22 + a^20 + a^18 + a^16 + a^14 + a^11 + a^10
            + a^8 + a^6 + a^5 + a^4 + a + 1

        Check that :issue:`21998` has been resolved::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(2^4)
            sage: R.<x> = K[]
            sage: f = x^2 + x + a^2 + a
            sage: r = f.any_root()
            sage: r^2 + r
            a^2 + a

        Check for :issue:`37034`::

            sage: R.<x> = Zmod(55)[]
            sage: (x^2 + 1).any_root()
            Traceback (most recent call last):
            ...
            ValueError: polynomial x^2 + 1 has no roots
        """
        # When not working over a finite field, do the simple thing of factoring for
        # roots and picking the first root. If none are available, raise an error.
        from sage.categories.finite_fields import FiniteFields
        if self.base_ring() not in FiniteFields():
            if ring not in FiniteFields():
                rs = self.roots(ring=ring, multiplicities=False)
                if rs:
                    return rs[0]
                raise ValueError(f"polynomial {self} has no roots")

        # Ensure that a provided ring is appropriate for the function. From the
        # above we know it is either None or a finite field. When it's a finite
        # field we ensure there's a coercion from the base ring to ring.
        if ring is not None:
            if ring.coerce_map_from(self.base_ring()) is None:
                raise ValueError(f"no coercion map can be computed from {self.base_ring()} to {ring}")

        # When the degree is none, we only look for a linear factor
        if degree is None:
            # When ring is None, we attempt to find a linear factor of self
            if ring is None:
                try:
                    f = self.any_irreducible_factor(degree=1, assume_squarefree=assume_squarefree)
                except ValueError:
                    raise ValueError(f"no root of polynomial {self} can be computed")
                return - f[0] / f[1]

            # When we have a ring, then we can find an irreducible factor of degree `d` providing
            # that d divides the degree of the extension from the base ring to the given ring
            allowed_extension_degree = ring.degree() // self.base_ring().degree()
            try:
                f = self.any_irreducible_factor(assume_squarefree=assume_squarefree, ext_degree=allowed_extension_degree)
            except ValueError:
                raise ValueError(f"no root of polynomial {self} can be computed over the ring {ring}")
            # When d != 1 we then find the smallest extension
            # TODO: What we should do here is compute some minimal
            #       extension F_ext = self.base_ring().extension(d, names='a') and find a
            #       root here and then coerce this root into the parent ring. This means we
            #       would work with the smallest possible extension.
            #       However, if we have some element of GF(p^k) and we try and coerce this to
            #       some element GF(p^(k*n)) this can fail, even though mathematically it
            #       should be fine.
            # TODO: Additionally, if the above was solved, it would be faster to extend the base
            #       ring with the irreducible factor however, if the base ring is an extension
            #       then the type of self.base_ring().extension(f) is a Univariate Quotient Polynomial Ring
            #       and not a finite field.

            # When f has degree one we simply return the roots
            # TODO: should we write something fast for degree two using
            #       the quadratic formula?
            if f.degree().is_one():
                root = - f[0] / f[1]
                return ring(root)

            # TODO: The proper thing to do here would be to call
            #       return f.change_ring(ring).any_root()
            #       but as we cannot work in the minimal extension (see above) working
            #       in the extension for f.change_ring(ring).any_root() is almost always
            #       much much much slower than using the call for roots() which uses
            #       C library bindings for all finite fields.
            #       Until the coercion system for finite fields works better,
            #       this will be the most performant
            return f.roots(ring, multiplicities=False)[0]

        # The old version of `any_root()` allowed degree < 0 to indicate that the input polynomial
        # had a distinct degree factorisation, we pass this to any_irreducible_factor as a bool and
        # ensure that the degree is positive.
        degree = ZZ(degree)
        if degree < 0:
            deprecation(37170, "negative ``degree`` will be disallowed. Instead use the bool `assume_equal_deg`.")
            degree = -degree
            assume_equal_deg = True

        # If a certain degree is requested, then we find an irreducible factor of degree `degree`
        # use this to compute a field extension and return the generator as root of this polynomial
        # if the degree and a ring is given however, instead compute a degree `degree` factor in the
        # base ring and then find a factor from this in the supplied ring.
        try:
            f = self.any_irreducible_factor(degree=degree,
                                            assume_squarefree=assume_squarefree,
                                            assume_equal_deg=assume_equal_deg)
        except ValueError:
            raise ValueError(f"no irreducible factor of degree {degree} can be computed from {self}")

        # For the case when the degree is one, just return the root
        if degree.is_one():
            root = - f[0] / f[1]
            if ring is None:
                return root
            return ring(root)

        # We are now in the case where the irreducible polynomial we have found is
        # of degree > 1. The old version of this function simply computed the roots
        # of this by calling f.roots(ring)... I don't really understand why
        # though, as we can simply ask for f.any_root() for this polynomial over the
        # new ring?
        if ring is None:
            # TODO: a faster option would be to create an extension with `f`
            #       as F_ext = self.base_ring().extension(f, names='a')
            #       however this returns a quotient ring rather than a
            #       FiniteField type if the base field is a non-prime field,
            #       so this slower option is chosen to ensure the root is
            #       over explicitly a FiniteField type.
            ring = self.base_ring().extension(degree, names='a')

        # TODO: The proper thing to do here would be to call
        #       return f.change_ring(ring).any_root()
        #       but as we cannot work in the minimal extension (see above) working
        #       in the extension for f.change_ring(ring).any_root() is almost always
        #       much much much slower than using the call for roots() which uses
        #       C library bindings for all finite fields.
        #       Until the coercion system for finite fields works better,
        #       this will be the most performant
        return f.roots(ring, multiplicities=False)[0]

    def __truediv__(left, right):
        r"""
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

            sage: R.<x> = PolynomialRing(GF(5^2, 'a'), 'x')                             # needs sage.rings.finite_rings
            sage: f = x^3 + 4*x
            sage: f / (x - 1)                                                           # needs sage.rings.finite_rings
            x^2 + x

        Be careful about coercions (this used to be broken)::

            sage: R.<x> = ZZ['x']
            sage: f = x / Mod(2,5); f
            3*x
            sage: f.parent()
            Univariate Polynomial Ring in x over Ring of integers modulo 5

        TESTS:

        Check that :issue:`12217` is fixed::

            sage: P.<x> = GF(5)[]
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(0, 5) does not exist

            sage: P.<x> = GF(25, 'a')[]                                                 # needs sage.rings.finite_rings
            sage: x/5                                                                   # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field

        Check that :issue:`23611` is fixed::

            sage: int(1) / x
            1/x
        """
        # Same parents => bypass coercion
        if have_same_parent(left, right):
            return (<Element>left)._div_(right)

        # Try division of polynomial by a scalar
        if isinstance(left, Polynomial):
            R = (<Polynomial>left)._parent._base
            try:
                x = R.coerce(right)
                return left * ~x
            except TypeError:
                pass

        # Delegate to coercion model. The line below is basically
        # RingElement.__truediv__(left, right), except that it also
        # works if left is not of type RingElement.
        return wrapperdescr_fastcall(RingElement.__truediv__,
                left, (right,), <object>NULL)

    def __pow__(left, right, modulus):
        """
        EXAMPLES::

            sage: x = polygen(QQ['u']['v'])
            sage: f = x - 1
            sage: f._pow(3)
            x^3 - 3*x^2 + 3*x - 1
            sage: f^3
            x^3 - 3*x^2 + 3*x - 1

            sage: R = PolynomialRing(GF(2), 'x')
            sage: f = R(x^9 + x^7 + x^6 + x^5 + x^4 + x^2 + x)
            sage: h = R(x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 1)
            sage: pow(f, 2, h)
            x^9 + x^8 + x^7 + x^5 + x^3

        TESTS::

            sage: x = polygen(QQ['u']['v'])
            sage: x^(1/2)
            Traceback (most recent call last):
            ...
            TypeError: non-integral exponents not supported

        ::

            sage: int(1)^x
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for ** or pow(): 'int' and 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'

        ::

            sage: x^x
            Traceback (most recent call last):
            ...
            TypeError: non-integral exponents not supported

        ::

            sage: # needs sage.rings.finite_rings
            sage: k = GF(5)
            sage: D.<x> = k[]
            sage: l.<x> = k.extension(x^2 + 2)
            sage: R.<t> = l[]
            sage: f = t^4 + (2*x - 1)*t^3 + (2*x + 1)*t^2 + 3
            sage: h = t^4 - x*t^3 + (3*x + 1)*t^2 + 2*t + 2*x - 1
            sage: pow(f, 2, h)
            3*t^3 + (2*x + 3)*t^2 + (2*x + 2)*t + 2*x + 2
            sage: pow(f, 10**7, h)
            4*x*t^3 + 2*x*t^2 + 4*x*t + 4

        Check that :issue:`18457` is fixed::

            sage: R.<x> = PolynomialRing(GF(5), sparse=True)
            sage: (1+x)^(5^10)
            x^9765625 + 1
            sage: S.<t> = GF(3)[]
            sage: R1.<x> = PolynomialRing(S, sparse=True)
            sage: (1+x+t)^(3^10)
            x^59049 + t^59049 + 1
            sage: R2.<x> = PolynomialRing(S, sparse=False)
            sage: (1+x+t)^(3^10)
            x^59049 + t^59049 + 1

        Check that the algorithm used is indeed correct::

            sage: # needs sage.rings.finite_rings
            sage: from sage.arith.power import generic_power
            sage: R1 = PolynomialRing(GF(8,'a'), 'x')
            sage: R2 = PolynomialRing(GF(9,'b'), 'x', sparse=True)
            sage: R3 = PolynomialRing(R2, 'y')
            sage: R4 = PolynomialRing(R1, 'y', sparse=True)
            sage: for d in range(20,40):        # long time
            ....:     for R in [R1, R2, R3, R3]:
            ....:         a = R.random_element()
            ....:         assert a^d == generic_power(a, d)

        Test the powering modulo ``x^n`` (calling :meth:`power_trunc`)::

            sage: R.<x> = GF(3)[]
            sage: pow(x + 1, 51, x^7)
            x^6 + 2*x^3 + 1

            sage: S.<y> = QQ[]
            sage: R.<x> = S[]
            sage: pow(y*x+1, 51, x^7)
            18009460*y^6*x^6 + 2349060*y^5*x^5 + ... + 51*y*x + 1

        Check that fallback method is used when it is not possible to compute
        the characteristic of the base ring (:issue:`24308`)::

            sage: # needs sage.libs.singular
            sage: kk.<a,b> = GF(2)[]
            sage: k.<y,w> = kk.quo(a^2 + a + 1)
            sage: K.<T> = k[]
            sage: (T*y)^21
            T^21
        """
        if not isinstance(left, Polynomial):
            return NotImplemented
        cdef Polynomial self = <Polynomial>left

        if type(right) is not Integer:
            try:
                right = Integer(right)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        d = self.degree()
        if d == -1:
            return self._parent(self._parent._base.zero() ** right)
        elif d == 0:
            return self._parent(self.get_unsafe(0) ** right)
        if right < 0:
            return (~self)**(-right)
        if modulus:
            if right > 0 and \
               parent(modulus) == self._parent and \
               modulus.number_of_terms() == 1 and \
               modulus.leading_coefficient().is_one():
                return self.power_trunc(right, modulus.degree())
            return power_mod(self, right, modulus)
        if self.is_gen():   # special case x**n should be faster!
            P = self._parent
            R = P._base
            if P.is_sparse():
                v = {right: R.one()}
            else:
                v = [R.zero()]*right + [R.one()]
            return self.parent()(v, check=False)
        if right > 20: # no gain below
            try:
                p = self.parent().characteristic()
            except (AttributeError, NotImplementedError):
                # some quotients do not implement characteristic
                # see github issue 24308
                p = -1
            if 0 < p <= right and (self.base_ring() in sage.categories.integral_domains.IntegralDomains() or p.is_prime()):
                x = self.parent().gen()
                ret = self.parent().one()
                e = 1
                q = right
                sparse = self.parent().is_sparse()
                if sparse:
                    d = self.monomial_coefficients()
                else:
                    c = self.list(copy=False)
                while q > 0:
                    q, r = q.quo_rem(p)
                    if r != 0:
                        if sparse:
                            tmp = self.parent()({e*k: d[k]**e for k in d})
                        else:
                            tmp = [0] * (e * len(c) - e + 1)
                            for i in range(len(c)):
                                tmp[e*i] = c[i]**e
                            tmp = self.parent()(tmp)
                        ret *= generic_power(tmp, r)
                    e *= p
                return ret

        return generic_power(self, right)

    def power_trunc(self, n, prec):
        r"""
        Truncated ``n``-th power of this polynomial up to precision ``prec``.

        INPUT:

        - ``n`` -- nonnegative integer; power to be taken

        - ``prec`` -- integer; the precision

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (3*x^2 - 2*x + 1).power_trunc(5, 8)
            -1800*x^7 + 1590*x^6 - 1052*x^5 + 530*x^4 - 200*x^3 + 55*x^2 - 10*x + 1
            sage: ((3*x^2 - 2*x + 1)^5).truncate(8)
            -1800*x^7 + 1590*x^6 - 1052*x^5 + 530*x^4 - 200*x^3 + 55*x^2 - 10*x + 1

            sage: S.<y> = R[]
            sage: (x + y).power_trunc(5,5)
            5*x*y^4 + 10*x^2*y^3 + 10*x^3*y^2 + 5*x^4*y + x^5
            sage: ((x + y)^5).truncate(5)
            5*x*y^4 + 10*x^2*y^3 + 10*x^3*y^2 + 5*x^4*y + x^5

            sage: R.<x> = GF(3)[]
            sage: p = x^2 - x + 1
            sage: q = p.power_trunc(80, 20); q
            x^19 + x^18 + ... + 2*x^4 + 2*x^3 + x + 1
            sage: (p^80).truncate(20) == q
            True

            sage: R.<x> = GF(7)[]
            sage: p = (x^2 + x + 1).power_trunc(2^100, 100); p
            2*x^99 + x^98 + x^95 + 2*x^94 + ... + 3*x^2 + 2*x + 1

            sage: for i in range(100):
            ....:    q1 = (x^2 + x + 1).power_trunc(2^100 + i, 100)
            ....:    q2 = p * (x^2 + x + 1).power_trunc(i, 100)
            ....:    q2 = q2.truncate(100)
            ....:    assert q1 == q2, "i = {}".format(i)

        TESTS::

            sage: x = polygen(QQ)
            sage: (3*x - 5).power_trunc(2^200, 0)
            0
            sage: x.power_trunc(-1, 10)
            Traceback (most recent call last):
            ...
            ValueError: n must be a nonnegative integer
            sage: R.<y> = QQ['x']
            sage: y.power_trunc(2**32 - 1, 2)
            0
            sage: y.power_trunc(2**64 - 1, 2)
            0
        """
        cdef Integer ZZn = ZZ(n)
        if mpz_fits_ulong_p(ZZn.value):
            return self._power_trunc(mpz_get_ui(ZZn.value), prec)
        return generic_power_trunc(self, ZZn, pyobject_to_long(prec))

    cpdef Polynomial _power_trunc(self, unsigned long n, long prec):
        r"""
        Truncated ``n``-th power of this polynomial up to precision ``prec``.

        This method is overridden for certain subclasses when a library function
        is available.

        INPUT:

        - ``n`` -- nonnegative integer; power to be taken

        - ``prec`` -- integer; the precision

        TESTS::

            sage: R.<x> = QQ['y'][]
            sage: for p in [R.one(), x, x+1, x-1, x^2 - 1]:
            ....:     for n in range(0, 20):
            ....:         for prec in [1, 2, 3, 10]:
            ....:             assert p._power_trunc(n, prec) == (p**n).truncate(prec)
        """
        return generic_power_trunc(self, Integer(n), prec)

    def _pow(self, right):
        # TODO: fit __pow__ into the arithmetic structure
        d = self.degree()
        if d == -1:
            return self._parent(self._parent._base.zero() ** right)
        elif d == 0:
            return self._parent(self.get_unsafe(0) ** right)
        if right < 0:
            return (~self)**(-right)
        if (<Polynomial>self) == self._parent.gen():   # special case x**n should be faster!
            v = [0]*right + [1]
            return self._parent(v, check=True)
        return generic_power(self, right)

    def _repr(self, name=None):
        """
        Return the string representation of this polynomial.

        INPUT:

        - ``name`` -- ``None`` or a string; used for printing the variable

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

        TESTS:

        We verify that :issue:`23020` has been resolved. (There are no
        elements in the Sage library yet that do not implement
        ``__bool__``, so we have to create one artificially.)::

            sage: class PatchedAlgebraicNumber(sage.rings.qqbar.AlgebraicNumber):       # needs sage.rings.number_field
            ....:     def __bool__(self): raise NotImplementedError()
            sage: R.<x> = QQbar[]                                                       # needs sage.rings.number_field
            sage: R([PatchedAlgebraicNumber(0), 1])                                     # needs sage.rings.number_field
            x + 0
        """
        if name is None:
            name = self._parent._names[0]
        # p-adic polynomials like (1 + O(p^20))*x have _is_gen set to True but
        # want their coefficient printed with its O() term
        if self._is_gen and not isinstance(self._parent._base, pAdicGeneric):
            return name
        sbuf = StringIO()
        sbuf.write(" ")
        m = self.degree() + 1
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        coeffs = self.list(copy=False)
        for n in reversed(range(m)):
            x = coeffs[n]
            is_nonzero = False
            try:
                is_nonzero = bool(x)
            except NotImplementedError:
                # for some elements it is not possible/feasible to determine
                # whether they are zero or not; we just print them anyway in
                # such cases
                is_nonzero = True
            if is_nonzero:
                if n != m-1:
                    sbuf.write(" + ")
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
                sbuf.write(x)
                sbuf.write(var)
        s = sbuf.getvalue()
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

            sage: R.<x> = PolynomialRing(QQ)
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

            sage: # needs sage.rings.number_field
            sage: C3.<omega> = CyclotomicField(3)
            sage: R.<X> = C3[]
            sage: f = X^3 - omega*X
            sage: latex(f)
            X^{3} - \omega X

            sage: R.<x> = RDF[]
            sage: latex(x+2)
            x + 2.0

        The following illustrates the fix of :issue:`2586`::

            sage: latex(ZZ['alpha']['b']([0, ZZ['alpha'].0]))
            \alpha b

        The following illustrates a (non-intentional) superfluity of parentheses

            sage: K.<I> = QuadraticField(-1)                                            # needs sage.rings.number_field
            sage: R.<x> = K[]                                                           # needs sage.rings.number_field
            sage: latex(I*x^2 - I*x)                                                    # needs sage.rings.number_field
            \left(\sqrt{-1}\right) x^{2} + \left(-\sqrt{-1}\right) x
        """
        s = " "
        coeffs = self.list(copy=False)
        m = len(coeffs)
        if name is None:
            name = self._parent.latex_variable_names()[0]
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        for n in reversed(range(m)):
            x = coeffs[n]
            x = y = latex(x)
            if x != '0':
                if n != m-1:
                    s += " + "
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "\\left(%s\\right)" % x
                if n > 1:
                    var = "|%s^{%s}" % (name, n)
                elif n==1:
                    var = "|%s" % name
                else:
                    var = ""
                s += "%s %s" % (x, var)
        s = s.replace(" + -", " - ")
        s = re.sub(r" 1(\.0+)? \|", " ", s)
        s = re.sub(r" -1(\.0+)? \|", " -", s)
        s = s.replace("|", "")
        if s == " ":
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
            sage: sage_input(RR(pi) * polygen(RR), verify=True)                         # needs sage.symbolic
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
            gen = sib.gen(self._parent)
            coeffs = self.list(copy=False)
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
            return sib(self._parent)(sib(self.constant_coefficient(), True))

    def __setitem__(self, n, value):
        """
        Set the `n`-th coefficient of this polynomial. This always raises an
        :exc:`IndexError`, since in Sage polynomials are immutable.

        INPUT:

        - ``n`` -- integer

        - ``value`` -- value to set the `n`-th coefficient to

        OUTPUT: an :exc:`IndexError` is always raised

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^3 + x + 1
            sage: f[2] = 3
            Traceback (most recent call last):
            ...
            IndexError: polynomials are immutable
        """
        raise IndexError("polynomials are immutable")

    cpdef _floordiv_(self, right):
        r"""
        Quotient of division of ``self`` by ``other``. This is denoted //.

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
        if (<Polynomial> right).is_one():
            # quite typical when removing gcds...
            return self
        Q, _ = self.quo_rem(right)
        return Q

    def __mod__(self, other):
        """
        Remainder of division of ``self`` by ``other``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x % (x+1)
            -1
            sage: (x^3 + x - 1) % (x^2 - 1)
            2*x - 1
        """
        _, R = self.quo_rem(other)
        return R

    def mod(self, other):
        """
        Remainder of division of ``self`` by ``other``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x % (x+1)
            -1
            sage: (x^3 + x - 1) % (x^2 - 1)
            2*x - 1

        TESTS:

        Check the problem reported at :issue:`12529` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: gens = 'y a0 a1 a2 b0 b1 b2 c1 c2 d0 d1 d2 d3 d4 d5 d6 d7'.split()
            sage: R = PolynomialRing(GF(8), 17, gens)
            sage: R.inject_variables(verbose=False)
            sage: A, B, C = a0 + a1*y + a2*y^2, b0 + b1*y + b2*y^2, c1*y + c2*y^2
            sage: D = d0 + d1*y + d2*y^2 + d3*y^3 + d4*y^4 + d5*y^5 + d6*y^6 + d7*y^7
            sage: F = D.subs({y: B})
            sage: G = A.subs({y: F}) + C
            sage: g = G.mod(y^8 + y)                                                    # needs sage.libs.singular
            sage: g.degree(y)                                                           # needs sage.libs.singular
            7
        """
        return self % other

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

    cpdef _mul_generic(self, right):
        """
        Compute the product of ``self`` and ``right`` using the classical quadratic
        algorithm. This method is the default for inexact rings.

        For two polynomials of degree n and m this method needs
        (m+1)*(n+1) products and n*m additions

        EXAMPLES::

            sage: K.<x> = QQ[]
            sage: f = 1+3*x+4*x^2+x^3
            sage: g = x^2+3*x^5
            sage: f._mul_generic(g)
            3*x^8 + 12*x^7 + 9*x^6 + 4*x^5 + 4*x^4 + 3*x^3 + x^2

        Show the product in the symbolic ring::

            sage: L = SR['x']                                                           # needs sage.symbolic
            sage: var('a0,a1,b0,b1')                                                    # needs sage.symbolic
            (a0, a1, b0, b1)
            sage: L([a0, a1])._mul_generic(L([b0, b1]))                                 # needs sage.symbolic
            a1*b1*x^2 + (a1*b0 + a0*b1)*x + a0*b0

        A non-commutative example::

            sage: # needs sage.combinat sage.modules
            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: R.<w> = PolynomialRing(A)
            sage: f = i*w + j
            sage: g = k*w + 1
            sage: f._mul_generic(g)
            -j*w^2 + 2*i*w + j
            sage: g._mul_generic(f)
            j*w^2 + j


        TESTS::

            sage: K.<x> = QQ[]
            sage: f = K(0)
            sage: g = K.random_element(10)
            sage: f._mul_generic(g)
            0
            sage: g._mul_generic(f)
            0
            sage: f._mul_generic(K(0))
            0
            sage: g._mul_generic(g) - g._mul_karatsuba(g)
            0
            sage: h = K(QQ.random_element(100,100))
            sage: f._mul_generic(h)
            0
            sage: K([h*c for c in g.list()]) - g._mul_generic(h)
            0
            sage: g._mul_generic(h) - K([h*c for c in g.list()])
            0
        """
        if self is right:
            return self._square_generic()
        cdef list x = self.list(copy=False)
        cdef list y = right.list(copy=False)
        return self._new_generic(do_schoolbook_product(x, y, -1))

    cdef _square_generic(self):
        cdef list x = self.list(copy=False)
        cdef Py_ssize_t i, j
        cdef Py_ssize_t d = len(x)-1
        zero = self._parent.base_ring().zero()
        two = self._parent.base_ring()(2)
        cdef list coeffs = [zero] * (2 * d + 1)
        for i from 0 <= i <= d:
            coeffs[2*i] = x[i] * x[i]
            for j from 0 <= j < i:
                coeffs[i+j] += two * x[i] * x[j]
        return self._new_generic(coeffs)

    def _mul_fateman(self, right):
        r"""
        Return the product of two polynomials using Kronecker's trick to
        do the multiplication. This could be used over a generic base
        ring.

        .. NOTE::

           -  Since this is implemented in interpreted Python, it could be
              hugely sped up by reimplementing it in Pyrex.

           -  Over the reals there is precision loss, at least in the current
              implementation.


        INPUT:

        - ``self`` -- polynomial

        - ``right`` -- polynomial (over same base ring as ``self``)

        OUTPUT: polynomial; the product ``self * right``

        ALGORITHM: Based on a paper by R. Fateman

        http://www.cs.berkeley.edu/~fateman/papers/polysbyGMP.pdf

        The idea is to encode dense univariate polynomials as big integers,
        instead of sequences of coefficients. The paper argues that because
        integer multiplication is so cheap, that encoding 2 polynomials to
        big numbers and then decoding the result might be faster than
        popular multiplication algorithms. This seems true when the degree
        is larger than 200.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: S.<y> = PolynomialRing(RR)
            sage: f = y^10 - 1.393493*y + 0.3
            sage: f._mul_karatsuba(f,0)
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

        - Slower over R when the degree of both of polynomials is less
          than 250 (roughly).

        - Over R, results may not be as accurate as the Karatsuba case.
          This is because we represent coefficients of polynomials over R as
          fractions, then convert them back to floating-point numbers.

        AUTHORS:

        - Didier Deshommes (2006-05-25)
        """
        return self._parent(polynomial_fateman._mul_fateman_mul(self, right))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.overflowcheck(False)
    def _mul_karatsuba(self, right, K_threshold = None):
        r"""
        Compute the product of two polynomials using the Karatsuba divide
        and conquer multiplication algorithm. This is only used over a
        generic base ring. (Special libraries like Flint are used, e.g., for
        the integers and rationals, which are much faster.)

        INPUT:

        - ``self`` -- polynomial
        - ``right`` -- polynomial (over same base ring as ``self``)
        - ``K_threshold`` -- (optional) Integer. A threshold to fall back to
          schoolbook algorithm. In the recursion, if one of the polynomials
          is of degree less that K_threshold then the classic quadratic
          polynomial is used.

        OUTPUT: the product ``self * right``

        ALGORITHM:

        The basic idea is to use that

        .. MATH::

            (aX + b) (cX + d) = acX^2 + ((a+b)(c+d)-ac-bd)X + bd


        where ac=a\*c and bd=b\*d, which requires three multiplications
        instead of the naive four. Given f and g of arbitrary degree bigger
        than one, let e be min(deg(f),deg(g))/2. Write

        .. MATH::

            f = a X^e + b   \text{ and }   g = c X^e + d


        and use the identity

        .. MATH::

            (aX^e + b) (cX^e + d) = ac X^{2e} +((a+b)(c+d) - ac - bd)X^e + bd


        to recursively compute `fg`.

        If ``self`` is a polynomial of degree n and `right` is a polynomial of
        degree m with n < m, then we interpret `right` as

        .. MATH::

            g0 + g1 * x^n + g2 * x^{2n} + ... + gq * x^{nq}

        where `gi` are polynomials of degree <= n. We then compute each product
        ``gi*right`` with Karatsuba multiplication and reconstruct ``self*right``
        from the partial products.

        The theoretical complexity for multiplying two polynomials of the same
        degree n is O(n^log(3,2)). Through testing of polynomials of degree up
        to 5000 we get that the number of operations for two polynomials of
        degree up to n-1 is bounded by:

        7.53*n**1.59 additions and 1.46*n**1.59 products on the base ring.

        For polynomials of degree m-1 and n-1 with m<n the number of operations
        is bounded by:

        8.11*m**0.59*n additions and 1.56*m**0.59*n products.

        (The bound might be worse for larger degrees.)

        EXAMPLES::

            sage: K.<x> = QQ[]
            sage: f = 1+3*x+4*x^2+x^3
            sage: g = x^2+3*x^5
            sage: f._mul_karatsuba(g,0)
            3*x^8 + 12*x^7 + 9*x^6 + 4*x^5 + 4*x^4 + 3*x^3 + x^2
            sage: f._mul_karatsuba(g,2)
            3*x^8 + 12*x^7 + 9*x^6 + 4*x^5 + 4*x^4 + 3*x^3 + x^2

        Show the product in the symbolic ring::

            sage: # needs sage.symbolic
            sage: L = SR['x']
            sage: var('a0,a1,b0,b1')
            (a0, a1, b0, b1)
            sage: L([a0, a1])._mul_karatsuba(L([b0, b1]), 0)
            a1*b1*x^2 + ((a0 + a1)*(b0 + b1) - a0*b0 - a1*b1)*x + a0*b0
            sage: L([a0, a1])._mul_karatsuba(L([b0, b1]), 2)
            a1*b1*x^2 + (a1*b0 + a0*b1)*x + a0*b0

        A noncommutative example::

            sage: # needs sage.combinat sage.modules
            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: R.<w> = PolynomialRing(A)
            sage: f = i*w + j
            sage: g = k*w + 1
            sage: f._mul_karatsuba(g,0)
            -j*w^2 + 2*i*w + j
            sage: g._mul_karatsuba(f,0)
            j*w^2 + j

        TESTS::

            sage: K.<x> = QQ[]
            sage: f = K(0)
            sage: g = K.random_element(10)
            sage: f._mul_karatsuba(g,0)
            0
            sage: g._mul_karatsuba(f,0)
            0
            sage: f._mul_karatsuba(K(0),0)
            0
            sage: g._mul_generic(g) - g._mul_karatsuba(g,0)
            0
            sage: h = K(QQ.random_element(100,100))
            sage: f._mul_karatsuba(h)
            0
            sage: K([h*c for c in g.list()]) - g._mul_generic(h)
            0
            sage: g._mul_karatsuba(h) - K([h*c for c in g.list()])
            0

        Random tests for noncommutative rings::

            sage: # needs sage.combinat sage.modules
            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: R.<w> = PolynomialRing(A)
            sage: f = R.random_element(randint(10,100))
            sage: g = R.random_element(randint(10,100))
            sage: f._mul_generic(g) == f._mul_karatsuba(g,0)
            True
            sage: f._mul_generic(g) == f._mul_karatsuba(g,16)
            True
            sage: g = R.random_element(0)
            sage: f._mul_karatsuba(g,0) == f._mul_generic(g)
            True
            sage: g._mul_karatsuba(f,0) == g._mul_generic(f)
            True

        Polynomials over matrices::

            sage: # needs sage.modules
            sage: K = PolynomialRing(MatrixSpace(QQ, 2), 'x')
            sage: f = K.random_element(randint(5, 10))
            sage: g = K.random_element(randint(5, 10))
            sage: h1 = f._mul_generic(g)
            sage: h2 = f._mul_karatsuba(g,randint(0, 10))
            sage: h1 == h2
            True
        """
        if self.is_zero():
            return self
        elif right.is_zero():
            return right
        cdef list f = self.list(copy=False)
        cdef list g = right.list(copy=False)
        n = len(f)
        m = len(g)
        if n == 1:
            c = f[0]
            return self._new_generic([c*a for a in g])
        if m == 1:
            c = g[0]
            return self._new_generic([a*c for a in f])
        if K_threshold is None:
            K_threshold = self._parent._Karatsuba_threshold
        if n <= K_threshold or m <= K_threshold:
            return self._new_generic(do_schoolbook_product(f, g, -1))
        if n == m:
            return self._new_generic(do_karatsuba(f,g, K_threshold, 0, 0, n))
        return self._new_generic(do_karatsuba_different_size(f,g, K_threshold))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Polynomial _mul_term(self, Polynomial term, bint term_on_right):
        """
        Return the product ``self * term``, where ``term`` is a polynomial
        with a single term.
        """
        cdef Py_ssize_t d = term.degree()
        cdef Py_ssize_t i
        cdef list coeffs, output
        c = term.get_unsafe(d)
        if term_on_right:
            coeffs = [a * c for a in self.list(copy=False)]
        else:
            coeffs = [c * a for a in self.list(copy=False)]
        if d == 0:
            return self._new_generic(coeffs)
        # shift by d
        output = [self.base_ring().zero()] * d
        output.extend(coeffs)
        return self._new_generic(output)

    def base_ring(self):
        """
        Return the base ring of the parent of ``self``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x.base_ring()
            Integer Ring
            sage: (2*x + 3).base_ring()
            Integer Ring
        """
        return self._parent.base_ring()

    cpdef base_extend(self, R):
        """
        Return a copy of this polynomial but with coefficients in ``R``, if
        there is a natural map from the coefficient ring of ``self`` to ``R``.

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
        S = self._parent.base_extend(R)
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
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        R = PolynomialRing(self._parent.base_ring(), names=var)
        return R(self)

    def change_ring(self, R):
        """
        Return a copy of this polynomial but with coefficients in ``R``, if at
        all possible.

        INPUT:

        - ``R`` -- a ring or morphism

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)                                            # needs sage.rings.number_field
            sage: f = K.defining_polynomial()                                           # needs sage.rings.number_field
            sage: f.change_ring(GF(7))                                                  # needs sage.rings.finite_rings sage.rings.number_field
            x^2 + x + 1

        ::

            sage: # needs sage.rings.number_field
            sage: K.<z> = CyclotomicField(3)
            sage: R.<x> = K[]
            sage: f = x^2 + z
            sage: f.change_ring(K.embeddings(CC)[1])                                    # needs sage.rings.real_mpfr
            x^2 - 0.500000000000000 - 0.866025403784438*I

        ::

            sage: R.<x> = QQ[]
            sage: f = x^2 + 1
            sage: f.change_ring(QQ.embeddings(CC)[0])                                   # needs sage.rings.real_mpfr
            x^2 + 1.00000000000000

        TESTS:

        Check that :issue:`25022` is fixed::

            sage: K.<x> = ZZ[]
            sage: x.change_ring(SR) == SR['x'].gen()                                    # needs sage.symbolic
            True
            sage: x.change_ring(ZZ['x']) == ZZ['x']['x'].gen()
            True

        Check that :issue:`28541` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(7^2)
            sage: S.<x> = F[]
            sage: P = x^2 + a*x + a^2
            sage: P.change_ring(F.frobenius_endomorphism())
            x^2 + (6*a + 1)*x + 6*a + 5
        """
        if isinstance(R, Map):
            return self.map_coefficients(R)
        else:
            return self._parent.change_ring(R)(self.list(copy=False))

    cpdef dict _mpoly_dict_recursive(self, tuple variables=None, base_ring=None):
        """
        Return a dict of coefficient entries suitable for construction of a
        MPolynomial_polydict with the given variables.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0)._mpoly_dict_recursive()
            {}
            sage: f = 7*x^5 + x^2 - 2*x - 3
            sage: f._mpoly_dict_recursive()
            {(0,): -3, (1,): -2, (2,): 1, (5,): 7}
        """
        if not self:
            return {}

        cdef ETuple const_ix
        cdef list mpolys
        cdef Py_ssize_t k

        var = self._parent.variable_name()
        if variables is None:
            variables = self._parent.variable_names_recursive()
        if var not in variables:
            x = base_ring(self) if base_ring else self
            const_ix = ETuple((0,)*len(variables))
            return { const_ix: x }

        cdef tuple prev_variables = variables[:variables.index(var)]
        const_ix = ETuple((0,)*len(prev_variables))
        mpolys = None

        if prev_variables:
            try:
                mpolys = [a._mpoly_dict_recursive(prev_variables, base_ring)
                          for a in self]
            except AttributeError:
                pass

        if mpolys is None:
            if base_ring is not None and base_ring is not self.base_ring():
                mpolys = [{const_ix: base_ring(a)} if a else {} for a in self]
            else:
                mpolys = [{const_ix: a} if a else {} for a in self]

        cdef dict D = {}
        cdef tuple leftovers = (0,) * (len(variables) - len(prev_variables) - 1)
        for k in range(len(mpolys)):
            for i, a in mpolys[k].items():
                j = ETuple((k,) + leftovers)
                D[i + j] = a

        return D

    def __copy__(self):
        """
        Return a "copy" of ``self``.

        This is just ``self``, since in Sage polynomials are immutable.

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
        degree `-1`.

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

    def euclidean_degree(self):
        r"""
        Return the degree of this element as an element of a Euclidean domain.

        If this polynomial is defined over a field, this is simply its :meth:`degree`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.euclidean_degree()
            1
            sage: R.<x> = ZZ[]
            sage: x.euclidean_degree()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        from sage.categories.fields import Fields
        if self.base_ring() in Fields():
            return self.degree()
        raise NotImplementedError

    def denominator(self):
        """
        Return a denominator of ``self``.

        First, the lcm of the denominators of the entries of ``self``
        is computed and returned. If this computation fails, the
        unit of the parent of ``self`` is returned.

        Note that some subclasses may implement their own
        :meth:`denominator` method. For example, see
        :class:`sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint`

        .. warning::

           This is not the denominator of the rational function
           defined by ``self``, which would always be 1 since ``self`` is a
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
        not have a :meth:`denominator` method.

        ::

            sage: # needs sage.rings.real_mpfr
            sage: R.<x> = RR[]
            sage: f = x + RR('0.3'); f
            x + 0.300000000000000
            sage: f.denominator()
            1.00000000000000

        Check that the denominator is an element over the base whenever the base
        has no :meth:`denominator` method. This closes :issue:`9063`. ::

            sage: R.<a> = GF(5)[]
            sage: x = R(0)
            sage: x.denominator()
            1
            sage: type(x.denominator())
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: isinstance(x.numerator() / x.denominator(), Polynomial)
            True
            sage: isinstance(x.numerator() / R(1), Polynomial)
            False

        TESTS:

        Check that :issue:`18518` is fixed::

            sage: R.<x> = PolynomialRing(QQ, sparse=True)
            sage: p = x^(2^100) - 1/2
            sage: p.denominator()
            2
        """

        if self.degree() == -1:
            return self.base_ring().one()
        x = self.coefficients()
        try:
            d = x[0].denominator()
            for y in x:
                d = d.lcm(y.denominator())
            return d
        except(AttributeError):
            return self.base_ring().one()

    def numerator(self):
        """
        Return a numerator of ``self``, computed as ``self * self.denominator()``.

        Note that some subclasses may implement its own numerator
        function. For example, see
        :class:`sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint`

        .. warning::

          This is not the numerator of the rational function
          defined by ``self``, which would always be ``self`` since ``self`` is a
          polynomial.

        EXAMPLES:

        First we compute the numerator of a polynomial with
        integer coefficients, which is of course ``self``.

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

            sage: # needs sage.rings.real_mpfr
            sage: R.<x> = RR[]
            sage: f = x + RR('0.3'); f
            x + 0.300000000000000
            sage: f.numerator()
            x + 0.300000000000000

        We check that the computation of the numerator and denominator
        are valid.

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: K = NumberField(symbolic_expression('x^3+2'), 'a')['s,t']['x']
            sage: f = K.random_element()
            sage: f.numerator() / f.denominator() == f
            True

            sage: R = RR['x']
            sage: f = R.random_element()
            sage: f.numerator() / f.denominator() == f
            True
        """
        return self * self.denominator()

    def derivative(self, *args):
        r"""
        The formal derivative of this polynomial, with respect to variables
        supplied in ``args``.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global :func:`derivative` function for more
        details.

        .. SEEALSO::

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

        .. SEEALSO::

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

        It is not possible to differentiate with respect to 2\*x for instance::

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

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

        Check that :issue:`28147` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(65537)[]
            sage: p = x^4 - 17*x^3 + 2*x^2 - x + 7
            sage: p.derivative()
            4*x^3 + 65486*x^2 + 4*x + 65536
            sage: R.<x> = GF(19^2)[]
            sage: p = x^4 - 17*x^3 + 2*x^2 - x + 7
            sage: p.derivative()
            4*x^3 + 6*x^2 + 4*x + 18

            sage: R.<x> = GF(2)[]
            sage: p = x^4 + x^2 + x
            sage: p.derivative()
            1

            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f._derivative()
            4*x^3 + 76
            sage: f._derivative(None)
            4*x^3 + 76

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: # needs sage.symbolic
            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y


        Check that :issue:`26844` is fixed by :issue:`28147`::

            sage: A = PolynomialRing(GF(3), name='t')
            sage: K = A.fraction_field()
            sage: t = K.gen()
            sage: t.derivative(t)
            1

        Check that :issue:`28187` is fixed::

            sage: R.<x> = GF(65537)[]                                                   # needs sage.rings.finite_rings
            sage: x._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: # needs sage.symbolic
            sage: y = var('y')
            sage: R.gen()._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var != self._parent.gen():
            try:
                # call _derivative() recursively on coefficients
                return self._parent([coeff._derivative(var) for coeff in self.list(copy=False)])
            except AttributeError:
                raise ValueError(f'cannot differentiate with respect to {var}')

        # compute formal derivative with respect to generator
        if self.is_zero():
            return self
        cdef Py_ssize_t n, degree = self.degree()
        if degree == 0:
            return self._parent.zero()
        coeffs = self.list(copy=False)
        return self._new_generic([n*coeffs[n] for n from 1 <= n <= degree])

    def gradient(self):
        """
        Return a list of the partial derivative of ``self``
        with respect to the variable of this univariate polynomial.

        There is only one partial derivative.

        EXAMPLES::

           sage: P.<x> = QQ[]
           sage: f = x^2 + (2/3)*x + 1
           sage: f.gradient()
           [2*x + 2/3]
           sage: f = P(1)
           sage: f.gradient()
           [0]
        """
        return [self.diff()]

    def integral(self,var=None):
        """
        Return the integral of this polynomial.

        By default, the integration variable is the variable of the
        polynomial.

        Otherwise, the integration variable is the optional parameter ``var``

        .. NOTE::

            The integral is always chosen so that the constant term is 0.

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

        This shows that the issue at :issue:`7711` is resolved::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,z> = PolynomialRing(GF(2147483647))
            sage: Q.<y> = PolynomialRing(P)
            sage: p = x + y + z
            sage: p.integral()
            -1073741823*y^2 + (x + z)*y

            sage: # needs sage.rings.finite_rings
            sage: P.<x,z> = PolynomialRing(GF(next_prime(2147483647)))
            sage: Q.<y> = PolynomialRing(P)
            sage: p = x + y + z
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
            sage: f = x*t + 5*t^2
            sage: f.integral(x)
            5*x*t^2 + 1/2*x^2*t

        TESTS:

        Check that :issue:`18600` is fixed::

            sage: Sx.<x> = ZZ[]
            sage: Sxy.<y> = Sx[]
            sage: Sxyz.<z> = Sxy[]
            sage: p = 1 + x*y + x*z + y*z^2
            sage: q = p.integral()
            sage: q
            1/3*y*z^3 + 1/2*x*z^2 + (x*y + 1)*z
            sage: q.parent()
            Univariate Polynomial Ring in z over Univariate Polynomial Ring in y
            over Univariate Polynomial Ring in x over Rational Field
            sage: q.derivative() == p
            True
            sage: p.integral(y)
            1/2*y^2*z^2 + x*y*z + 1/2*x*y^2 + y
            sage: p.integral(y).derivative(y) == p
            True
            sage: p.integral(x).derivative(x) == p
            True

        Check that it works with non-integral domains (:issue:`18600`)::

            sage: x = polygen(Zmod(4))
            sage: p = x**4 + 1
            sage: p.integral()
            x^5 + x
            sage: p.integral().derivative() == p
            True
        """
        R = self._parent

        # TODO:
        # calling the coercion model bin_op is much more accurate than using the
        # true division (which is bypassed by polynomials). But it does not work
        # in all cases!!
        cm = coercion_model
        try:
            S = cm.bin_op(R.one(), ZZ.one(), operator.truediv).parent()
            Q = S.base_ring()
        except TypeError:
            Q = (R.base_ring().one() / ZZ.one()).parent()
            S = R.change_ring(Q)
        if var is not None and var != R.gen():
            # call integral() recursively on coefficients
            return S([coeff.integral(var) for coeff in self])
        cdef Py_ssize_t n
        zero = Q.zero()
        p = [zero] + [cm.bin_op(Q(self.get_unsafe(n)), n + 1, operator.truediv)
                      if self.get_unsafe(n) else zero for n in range(self.degree() + 1)]
        return S(p)

    def monomial_coefficients(self):
        """
        Return a sparse dictionary representation of this univariate
        polynomial.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + -1/7*x + 13
            sage: f.monomial_coefficients()
            {0: 13, 1: -1/7, 3: 1}

        ``dict`` is an alias::

            sage: f.dict()
            {0: 13, 1: -1/7, 3: 1}
        """
        cdef dict X = {}
        cdef list Y = self.list(copy=False)
        cdef Py_ssize_t i
        for i in range(len(Y)):
            c = Y[i]
            if c:
                X[i] = c
        return X

    dict = monomial_coefficients

    def factor(self, **kwargs):
        r"""
        Return the factorization of ``self`` over its base ring.

        INPUT:

        - ``kwargs`` -- any keyword arguments are passed to the method
          :meth:`_factor_univariate_polynomial` of the base ring if it
          defines such a method.

        OUTPUT:

        A factorization of ``self`` over its parent into a unit and
        irreducible factors.  If the parent is a polynomial ring
        over a field, these factors are monic.

        EXAMPLES:

        Factorization is implemented over various rings. Over `\QQ`::

            sage: x = QQ['x'].0
            sage: f = (x^3 - 1)^2
            sage: f.factor()                                                            # needs sage.libs.pari
            (x - 1)^2 * (x^2 + x + 1)^2

        Since `\QQ` is a field, the irreducible factors are monic::

            sage: f = 10*x^5 - 1
            sage: f.factor()                                                            # needs sage.libs.pari
            (10) * (x^5 - 1/10)
            sage: f = 10*x^5 - 10
            sage: f.factor()                                                            # needs sage.libs.pari
            (10) * (x - 1) * (x^4 + x^3 + x^2 + x + 1)

        Over `\ZZ` the irreducible factors need not be monic::

            sage: x = ZZ['x'].0
            sage: f = 10*x^5 - 1
            sage: f.factor()                                                            # needs sage.libs.pari
            10*x^5 - 1

        We factor a non-monic polynomial over a finite field of 25
        elements::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(25)
            sage: R.<x> = k[]
            sage: f = 2*x^10 + 2*x + 2*a
            sage: F = f.factor(); F
            (2) * (x + a + 2) * (x^2 + 3*x + 4*a + 4) * (x^2 + (a + 1)*x + a + 2)
            * (x^5 + (3*a + 4)*x^4 + (3*a + 3)*x^3 + 2*a*x^2 + (3*a + 1)*x + 3*a + 1)

        Notice that the unit factor is included when we multiply `F`
        back out::

            sage: expand(F)                                                             # needs sage.rings.finite_rings sage.symbolic
            2*x^10 + 2*x + 2*a

        A new ring.  In the example below, we set the special method
        ``_factor_univariate_polynomial()`` in the base ring which is
        called to factor univariate polynomials.  This facility can be
        used to easily extend polynomial factorization to work over
        new rings you introduce::

             sage: # needs sage.libs.ntl
             sage: R.<x> = PolynomialRing(IntegerModRing(4), implementation="NTL")
             sage: (x^2).factor()
             Traceback (most recent call last):
             ...
             NotImplementedError: factorization of polynomials over rings with
             composite characteristic is not implemented
             sage: def my_factor(f):
             ....:     return f.change_ring(ZZ).factor()
             sage: R.base_ring()._factor_univariate_polynomial = my_factor
             sage: (x^2).factor()                                                       # needs sage.libs.pari
             x^2
             sage: del R.base_ring()._factor_univariate_polynomial  # clean up

        Arbitrary precision real and complex factorization::

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: R.<x> = RealField(100)[]
            sage: F = factor(x^2 - 3); F
            (x - 1.7320508075688772935274463415) * (x + 1.7320508075688772935274463415)
            sage: expand(F)
            x^2 - 3.0000000000000000000000000000
            sage: factor(x^2 + 1)
            x^2 + 1.0000000000000000000000000000

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: R.<x> = ComplexField(100)[]
            sage: F = factor(x^2 + 3); F
            (x - 1.7320508075688772935274463415*I) * (x + 1.7320508075688772935274463415*I)
            sage: expand(F)
            x^2 + 3.0000000000000000000000000000
            sage: factor(x^2 + 1)
            (x - I) * (x + I)
            sage: f = R(I) * (x^2 + 1) ; f
            I*x^2 + I
            sage: F = factor(f); F
            (1.0000000000000000000000000000*I) * (x - I) * (x + I)
            sage: expand(F)
            I*x^2 + I

        Over a number field::

            sage: # needs sage.rings.number_field
            sage: K.<z> = CyclotomicField(15)
            sage: x = polygen(K)
            sage: ((x^3 + z*x + 1)^3 * (x - z)).factor()
            (x - z) * (x^3 + z*x + 1)^3
            sage: cyclotomic_polynomial(12).change_ring(K).factor()
            (x^2 - z^5 - 1) * (x^2 + z^5)
            sage: ((x^3 + z*x + 1)^3 * (x/(z+2) - 1/3)).factor()
            (-1/331*z^7 + 3/331*z^6 - 6/331*z^5 + 11/331*z^4
                - 21/331*z^3 + 41/331*z^2 - 82/331*z + 165/331)
            * (x - 1/3*z - 2/3) * (x^3 + z*x + 1)^3

        Over a relative number field::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<z> = CyclotomicField(3)
            sage: L.<a> = K.extension(x^3 - 2)
            sage: t = polygen(L, 't')
            sage: f = (t^3 + t + a) * (t^5 + t + z); f
            t^8 + t^6 + a*t^5 + t^4 + z*t^3 + t^2 + (a + z)*t + z*a
            sage: f.factor()
            (t^3 + t + a) * (t^5 + t + z)

        Over the real double field::

            sage: # needs numpy
            sage: R.<x> = RDF[]
            sage: (-2*x^2 - 1).factor()
            (-2.0) * (x^2 + 0.5000000000000001)
            sage: (-2*x^2 - 1).factor().expand()
            -2.0*x^2 - 1.0000000000000002
            sage: f = (x - 1)^3
            sage: f.factor()  # abs tol 2e-5
            (x - 1.0000065719436413) * (x^2 - 1.9999934280563585*x + 0.9999934280995487)

        The above output is incorrect because it relies on the
        :meth:`.roots` method, which does not detect that all the roots
        are real::

            sage: f.roots()  # abs tol 2e-5                                             # needs numpy
            [(1.0000065719436413, 1)]

        Over the complex double field the factors are approximate and
        therefore occur with multiplicity 1::

            sage: # needs numpy sage.rings.complex_double
            sage: R.<x> = CDF[]
            sage: f = (x^2 + 2*R(I))^3
            sage: F = f.factor()
            sage: F  # abs tol 3e-5
            (x - 1.0000138879287663 + 1.0000013435286879*I)
            * (x - 0.9999942196864997 + 0.9999873009803959*I)
            * (x - 0.9999918923847313 + 1.0000113554909125*I)
            * (x + 0.9999908759550227 - 1.0000069659624138*I)
            * (x + 0.9999985293216753 - 0.9999886153831807*I)
            * (x + 1.0000105947233 - 1.0000044186544053*I)
            sage: [f(t[0][0]).abs() for t in F]  # abs tol 1e-13
            [1.979365054e-14, 1.97936298566e-14, 1.97936990747e-14,
             3.6812407475e-14, 3.65211563729e-14, 3.65220890052e-14]

        Factoring polynomials over `\ZZ/n\ZZ` for
        composite `n` is not implemented::

            sage: R.<x> = PolynomialRing(Integers(35))
            sage: f = (x^2 + 2*x + 2) * (x^2 + 3*x + 9)
            sage: f.factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: factorization of polynomials over
            rings with composite characteristic is not implemented

        Factoring polynomials over the algebraic numbers (see
        :issue:`8544`)::

            sage: R.<x> = QQbar[]                                                       # needs sage.rings.number_field
            sage: (x^8 - 1).factor()                                                    # needs sage.rings.number_field
            (x - 1) * (x - 0.7071067811865475? - 0.7071067811865475?*I)
            * (x - 0.7071067811865475? + 0.7071067811865475?*I) * (x - I) * (x + I)
            * (x + 0.7071067811865475? - 0.7071067811865475?*I)
            * (x + 0.7071067811865475? + 0.7071067811865475?*I) * (x + 1)

        Factoring polynomials over the algebraic reals (see
        :issue:`8544`)::

            sage: R.<x> = AA[]                                                          # needs sage.rings.number_field
            sage: (x^8 + 1).factor()                                                    # needs sage.rings.number_field
            (x^2 - 1.847759065022574?*x + 1.000000000000000?)
            * (x^2 - 0.7653668647301795?*x + 1.000000000000000?)
            * (x^2 + 0.7653668647301795?*x + 1.000000000000000?)
            * (x^2 + 1.847759065022574?*x + 1.000000000000000?)

        TESTS:

        This came up in :issue:`7088`::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 12*x^10 + x^9 + 432*x^3 + 9011
            sage: g = 13*x^11 + 89*x^3 + 1
            sage: F = f^2 * g^3
            sage: F = f^2 * g^3; F.factor()                                             # needs sage.libs.pari
            (12*x^10 + x^9 + 432*x^3 + 9011)^2 * (13*x^11 + 89*x^3 + 1)^3
            sage: F = f^2 * g^3 * 7; F.factor()                                         # needs sage.libs.pari
            7 * (12*x^10 + x^9 + 432*x^3 + 9011)^2 * (13*x^11 + 89*x^3 + 1)^3

        This example came up in :issue:`7097`::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: f = 8*x^9 + 42*x^6 + 6*x^3 - 1
            sage: g = (x^24 - 12*x^23 + 72*x^22 - 286*x^21 + 849*x^20 - 2022*x^19 + 4034*x^18
            ....:       - 6894*x^17 + 10182*x^16 - 13048*x^15 + 14532*x^14 - 13974*x^13
            ....:       + 11365*x^12 - 7578*x^11 + 4038*x^10 - 1766*x^9 + 762*x^8 - 408*x^7
            ....:       + 236*x^6 - 126*x^5 + 69*x^4 - 38*x^3 + 18*x^2 - 6*x + 1)
            sage: assert g.is_irreducible()                                             # needs sage.libs.pari
            sage: K.<a> = NumberField(g)
            sage: len(f.roots(K))
            9
            sage: f.factor()                                                            # needs sage.libs.pari
            (8) * (x^3 + 1/4) * (x^6 + 5*x^3 - 1/2)
            sage: f.change_ring(K).factor()
            (8) * (x - 3260097/3158212*a^22 + 35861067/3158212*a^21 - 197810817/3158212*a^20
                     + 722970825/3158212*a^19 - 1980508347/3158212*a^18 + 4374189477/3158212*a^17
                     - 4059860553/1579106*a^16 + 6442403031/1579106*a^15 - 17542341771/3158212*a^14
                     + 20537782665/3158212*a^13 - 20658463789/3158212*a^12 + 17502836649/3158212*a^11
                     - 11908953451/3158212*a^10 + 6086953981/3158212*a^9 - 559822335/789553*a^8
                     + 194545353/789553*a^7 - 505969453/3158212*a^6 + 338959407/3158212*a^5
                     - 155204647/3158212*a^4 + 79628015/3158212*a^3 - 57339525/3158212*a^2
                     + 26692783/3158212*a - 1636338/789553) * ...
            sage: f = QQbar['x'](1)
            sage: f.factor()
            1

        Factorization also works even if the variable of the finite
        field is nefariously labeled `x`::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(3^2, 'x')[]
            sage: f = x^10 + 7*x - 13
            sage: G = f.factor(); G
            (x + x) * (x + 2*x + 1) * (x^4 + (x + 2)*x^3 + (2*x + 2)*x + 2)
            * (x^4 + 2*x*x^3 + (x + 1)*x + 2)
            sage: prod(G) == f
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x0> = GF(9,'x')[]  # purposely calling it x to test robustness
            sage: f = x0^3 + x0 + 1
            sage: f.factor()
            (x0 + 2) * (x0 + x) * (x0 + 2*x + 1)
            sage: f = 0*x0
            sage: f.factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: factorization of 0 is not defined

        ::

            sage: f = x0^0                                                              # needs sage.rings.finite_rings
            sage: f.factor()                                                            # needs sage.rings.finite_rings
            1

        Over a complicated number field::

            sage: # needs sage.rings.number_field
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
            [(T - a, 1),
             (T - 40085763200/924556084127*a^5 - 145475769880/924556084127*a^4
               + 527617096480/924556084127*a^3 + 1289745809920/924556084127*a^2
               - 3227142391585/924556084127*a - 401502691578/924556084127, 1)]
            sage: expand(F)
            T^6 + 10/7*T^5 - 867/49*T^4 - 76/245*T^3 + 3148/35*T^2 - 25944/245*T + 48771/1225

        ::

            sage: # needs sage.rings.number_field
            sage: f = x^2 - 1/3
            sage: K.<a> = NumberField(f)
            sage: A.<T> = K[]
            sage: A(x^2 - 1).factor()
            (T - 1) * (T + 1)
            sage: A(3*x^2 - 1).factor()
            (3) * (T - a) * (T + a)
            sage: A(x^2 - 1/3).factor()
            (T - a) * (T + a)

        Test that :issue:`10279` is fixed::

            sage: # needs sage.rings.number_field
            sage: R.<t> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(t^4 - t^2 + 1)
            sage: pol = t^3 + (-4*a^3 + 2*a)*t^2 - 11/3*a^2*t + 2/3*a^3 - 4/3*a
            sage: pol.factor()
            (t - 2*a^3 + a) * (t - 4/3*a^3 + 2/3*a) * (t - 2/3*a^3 + 1/3*a)

        Test that this factorization really uses ``nffactor()`` internally::

            sage: # needs sage.libs.pari sage.rings.number_field
            sage: pari.default("debug", 3)
            sage: F = pol.factor()
            <BLANKLINE>
            Entering nffactor:
            ...
            sage: pari.default("debug", 0)

        Test that :issue:`10369` is fixed::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
            sage: R.<t> = PolynomialRing(K)
            sage: pol = ((-1/7*a^5 - 1/7*a^4 - 1/7*a^3 - 1/7*a^2 - 2/7*a - 1/7)*t^10
            ....:        + (4/7*a^5 - 2/7*a^4 - 2/7*a^3 - 2/7*a^2 - 2/7*a - 6/7)*t^9
            ....:        + (90/49*a^5 + 152/49*a^4 + 18/49*a^3 + 24/49*a^2 + 30/49*a + 36/49)*t^8
            ....:        + (-10/49*a^5 + 10/7*a^4 + 198/49*a^3 - 102/49*a^2 - 60/49*a - 26/49)*t^7
            ....:        + (40/49*a^5 + 45/49*a^4 + 60/49*a^3 + 277/49*a^2 - 204/49*a - 78/49)*t^6
            ....:        + (90/49*a^5 + 110/49*a^4 + 2*a^3 + 80/49*a^2 + 46/7*a - 30/7)*t^5
            ....:        + (30/7*a^5 + 260/49*a^4 + 250/49*a^3 + 232/49*a^2 + 32/7*a + 8)*t^4
            ....:        + (-184/49*a^5 - 58/49*a^4 - 52/49*a^3 - 66/49*a^2 - 72/49*a - 72/49)*t^3
            ....:        + (18/49*a^5 - 32/49*a^4 + 10/49*a^3 + 4/49*a^2)*t^2
            ....:        + (2/49*a^4 - 4/49*a^3 + 2/49*a^2)*t)
            sage: pol.factor()
            (-1/7*a^5 - 1/7*a^4 - 1/7*a^3 - 1/7*a^2 - 2/7*a - 1/7) * t
             * (t - a^5 - a^4 - a^3 - a^2 - a - 1)^4
             * (t^5 + (-12/7*a^5 - 10/7*a^4 - 8/7*a^3 - 6/7*a^2 - 4/7*a - 2/7)*t^4
                    + (12/7*a^5 - 8/7*a^3 + 16/7*a^2 + 2/7*a + 20/7)*t^3
                    + (-20/7*a^5 - 20/7*a^3 - 20/7*a^2 + 4/7*a - 2)*t^2
                    + (12/7*a^5 + 12/7*a^3 + 2/7*a + 16/7)*t
                    - 4/7*a^5 - 4/7*a^3 - 4/7*a - 2/7)
            sage: pol = ((1/7*a^2 - 1/7*a)*t^10 + (4/7*a - 6/7)*t^9
            ....:        + (102/49*a^5 + 99/49*a^4 + 96/49*a^3 + 93/49*a^2 + 90/49*a + 150/49)*t^8
            ....:        + (-160/49*a^5 - 36/49*a^4 - 48/49*a^3 - 8/7*a^2 - 60/49*a - 60/49)*t^7
            ....:        + (30/49*a^5 - 55/49*a^4 + 20/49*a^3 + 5/49*a^2)*t^6
            ....:        + (6/49*a^4 - 12/49*a^3 + 6/49*a^2)*t^5)
            sage: pol.factor()
            (1/7*a^2 - 1/7*a) * t^5
             * (t^5 + (-40/7*a^5 - 38/7*a^4 - 36/7*a^3 - 34/7*a^2 - 32/7*a - 30/7)*t^4
                    + (60/7*a^5 - 30/7*a^4 - 18/7*a^3 - 9/7*a^2 - 3/7*a)*t^3
                    + (60/7*a^4 - 40/7*a^3 - 16/7*a^2 - 4/7*a)*t^2
                    + (30/7*a^3 - 25/7*a^2 - 5/7*a)*t + 6/7*a^2 - 6/7*a)
            sage: pol = (x^10 + (4/7*a - 6/7)*x^9 + (9/49*a^2 - 3/7*a + 15/49)*x^8
            ....:         + (8/343*a^3 - 32/343*a^2 + 40/343*a - 20/343)*x^7
            ....:         + (5/2401*a^4 - 20/2401*a^3 + 40/2401*a^2 - 5/343*a + 15/2401)*x^6
            ....:         + (-6/16807*a^4 + 12/16807*a^3 - 18/16807*a^2 + 12/16807*a - 6/16807)*x^5)
            sage: pol.factor()
            x^5 * (x^5 + (4/7*a - 6/7)*x^4 + (9/49*a^2 - 3/7*a + 15/49)*x^3
                    + (8/343*a^3 - 32/343*a^2 + 40/343*a - 20/343)*x^2
                    + (5/2401*a^4 - 20/2401*a^3 + 40/2401*a^2 - 5/343*a + 15/2401)*x
                    - 6/16807*a^4 + 12/16807*a^3 - 18/16807*a^2 + 12/16807*a - 6/16807)

        Factoring over a number field over which we cannot factor the
        discriminant by trial division::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^16 - x - 6)
            sage: R.<x> = PolynomialRing(K)
            sage: f = (x+a)^50 - (a-1)^50
            sage: len(factor(f))
            6
            sage: pari(K.discriminant()).factor(limit=10^6)
            [-1, 1; 3, 15; 23, 1; 887, 1; 12583, 1; 2354691439917211, 1]
            sage: factor(K.discriminant())
            -1 * 3^15 * 23 * 887 * 12583 * 6335047 * 371692813

        Factoring over a number field over which we cannot factor the
        discriminant and over which `nffactor()` fails::

            sage: # needs sage.libs.pari sage.rings.number_field
            sage: p = next_prime(10^50); q = next_prime(10^51); n = p*q
            sage: K.<a> = QuadraticField(p*q)
            sage: R.<x> = PolynomialRing(K)
            sage: K.pari_polynomial('a').nffactor("x^2+1")
            Mat([x^2 + 1, 1])
            sage: factor(x^2 + 1)
            x^2 + 1
            sage: factor((x - a) * (x + 2*a))
            (x - a) * (x + 2*a)

        A test where nffactor used to fail without a nf structure::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K = NumberField([x^2 - 1099511627777, x^3 - 3], 'a')
            sage: x = polygen(K)
            sage: f = x^3 - 3
            sage: factor(f)
            (x - a1) * (x^2 + a1*x + a1^2)

        We check that :issue:`7554` is fixed::

            sage: L.<q> = LaurentPolynomialRing(QQ)
            sage: F = L.fraction_field()
            sage: R.<x> = PolynomialRing(F)
            sage: factor(x)                                                             # needs sage.libs.pari
            x
            sage: factor(x^2 - q^2)
            (x - q) * (x + q)
            sage: factor(x^2 - q^-2)
            (x - 1/q) * (x + 1/q)

            sage: P.<a,b,c> = PolynomialRing(ZZ)
            sage: R.<x> = PolynomialRing(FractionField(P))
            sage: p = (x - a) * (b*x + c) * (a*b*x + a*c) / (a + 2)
            sage: factor(p)                                                             # needs sage.libs.singular
            (a/(a + 2)) * (x - a) * (b*x + c)^2

        Check that :issue:`24973` is fixed::

            sage: x1 = ZZ['x'].gen()
            sage: x2 = ZZ['x']['x'].gen()
            sage: (x1 - x2).factor()                                                    # needs sage.libs.singular
            -x + x

        Check that :issue:`26421' is fixed::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: P.<x> = R[]
            sage: p = x^4 + (-5 - 2*t)*x^3 + (-2 + 10*t)*x^2 + (10 + 4*t)*x - 20*t
            sage: p.factor()                                                            # needs sage.libs.singular
            (x - 5) * (x - 2*t) * (x^2 - 2)

        Check that :issue:`29266` is fixed:

            sage: f = t*x + t
            sage: f.is_irreducible()                                                    # needs sage.libs.singular
            True
            sage: f = 2*x + 4
            sage: f.is_irreducible()                                                    # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            NotImplementedError
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
        ## 200 lines of spaghetti code is just way to much!

        if self.degree() < 0:
            raise ArithmeticError("factorization of {!r} is not defined".format(self))
        if self.degree() == 0:
            return Factorization([], unit=self.get_unsafe(0))

        # Use multivariate implementations for polynomials over polynomial rings
        flatten = self._parent.flattening_morphism()
        tgt = flatten.codomain()
        if tgt is not self._parent and tgt._has_singular:
            try:
                F = flatten(self).factor(**kwargs)
                unflatten = flatten.section()
                return Factorization(((unflatten(f),m) for (f,m) in F), unit = F.unit())
            except NotImplementedError:
                pass

        R = self._parent.base_ring()
        if hasattr(R, '_factor_univariate_polynomial'):
            return R._factor_univariate_polynomial(self, **kwargs)

        G = None
        ch = R.characteristic()
        if not (ch == 0 or is_prime(ch)):
            raise NotImplementedError("factorization of polynomials over rings with composite characteristic is not implemented")

        from sage.rings.finite_rings.finite_field_base import FiniteField

        n = None

        if isinstance(R, sage.rings.abc.IntegerModRing) or isinstance(R, IntegerRing_class):
            try:
                G = list(self._pari_with_name().factor())
            except PariError:
                raise NotImplementedError

        elif isinstance(R, FiniteField):
            v = [x.__pari__("a") for x in self.list()]
            f = pari(v).Polrev()
            G = list(f.factor())

        if G is None:
            # See if we can do this as a singular polynomial as a fallback
            # This was copied from the general multivariate implementation
            try:
                if R.is_finite():
                    if R.characteristic() > 1<<29:
                        raise NotImplementedError("Factorization of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

                P = self._parent
                P._singular_().set_ring()
                S = self._singular_().factorize()
                factors = S[1]
                exponents = S[2]
                v = sorted([( P(factors[i+1]),
                              sage.rings.integer.Integer(exponents[i+1]))
                            for i in range(len(factors))])
                unit = P.one()
                for i in range(len(v)):
                    if v[i][0].is_unit():
                        unit = unit * v[i][0]
                        del v[i]
                        break
                F = Factorization(v, unit=unit)
                F.sort()
                return F
            except (TypeError, AttributeError):
                if R.is_integral_domain() and not R.is_field():
                    try:
                        F = R.fraction_field()
                        PF = F[self.variable_name()]
                        pol_frac = PF(self)
                        pol_frac_fact = pol_frac.factor(**kwargs)
                        if R(pol_frac_fact.unit()).is_unit():
                            # Note: :meth:`base_change` may convert the unit to a non unit
                            return pol_frac_fact.base_change(self.parent())
                    except (TypeError, AttributeError, NotImplementedError):
                        raise NotImplementedError

                raise NotImplementedError

        return self._factor_pari_helper(G, n)

    def _factor_pari_helper(self, G, n=None, unit=None):
        """
        Fix up and normalize a factorization that came from PARI.

        TESTS::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ)
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

            sage: # needs sage.libs.pari
            sage: f = prod((k^2*x^k + k)^(k-1) for k in primes(10))
            sage: F = f._factor_pari_helper(pari(f).factor()); F
            1323551250 * (2*x^2 + 1) * (3*x^3 + 1)^2 * (5*x^5 + 1)^4 * (7*x^7 + 1)^6
            sage: F.prod() == f
            True
            sage: QQ['x'](f)._factor_pari_helper(pari(f).factor())
            (1751787911376562500) * (x^2 + 1/2) * (x^3 + 1/3)^2 * (x^5 + 1/5)^4 * (x^7 + 1/7)^6

            sage: # needs sage.libs.pari
            sage: g = GF(19)['x'](f)
            sage: G = g._factor_pari_helper(pari(g).factor()); G
            (4) * (x + 3) * (x + 16)^5 * (x + 11)^6 * (x^2 + 7*x + 9)^4
            * (x^2 + 15*x + 9)^4 * (x^3 + 13)^2
            * (x^6 + 8*x^5 + 7*x^4 + 18*x^3 + 11*x^2 + 12*x + 1)^6
            sage: G.prod() == g
            True
        """
        pols, exps = G
        R = self._parent
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
            content_fix = R.base_ring().one()
            for f, e in F:
                if not f.is_monic():
                    content_fix *= f.leading_coefficient()**e
            unit //= content_fix
            if not unit.is_unit():
                F.append((R(unit), ZZ(1)))
                unit = R.base_ring().one()

        if n is not None:
            pari.set_real_precision(n)  # restore precision
        return Factorization(F, unit)

    def splitting_field(self, names=None, map=False, **kwds):
        """
        Compute the absolute splitting field of a given polynomial.

        INPUT:

        - ``names`` -- (default: ``None``) a variable name for the splitting field

        - ``map`` -- boolean (default: ``False``); also return an embedding of
          ``self`` into the resulting field

        - ``kwds`` -- additional keywords depending on the type.
          Currently, only number fields are implemented. See
          :func:`sage.rings.number_field.splitting_field.splitting_field`
          for the documentation of these keywords.

        OUTPUT:

        If ``map`` is ``False``, the splitting field as an absolute field.
        If ``map`` is ``True``, a tuple ``(K, phi)`` where ``phi`` is an
        embedding of the base field of ``self`` in ``K``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: K.<a> = (x^3 + 2).splitting_field(); K                                # needs sage.rings.number_field
            Number Field in a with defining polynomial
             x^6 + 3*x^5 + 6*x^4 + 11*x^3 + 12*x^2 - 3*x + 1
            sage: K.<a> = (x^3 - 3*x + 1).splitting_field(); K                          # needs sage.rings.number_field
            Number Field in a with defining polynomial x^3 - 3*x + 1

        Relative situation::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 + 2)
            sage: S.<t> = PolynomialRing(K)
            sage: L.<b> = (t^2 - a).splitting_field()
            sage: L
            Number Field in b with defining polynomial t^6 + 2

        With ``map=True``, we also get the embedding of the base field
        into the splitting field::

            sage: L.<b>, phi = (t^2 - a).splitting_field(map=True)                      # needs sage.rings.number_field
            sage: phi                                                                   # needs sage.rings.number_field
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Number Field in b with defining polynomial t^6 + 2
              Defn: a |--> b^2

        An example over a finite field::

            sage: P.<x> = PolynomialRing(GF(7))
            sage: t = x^2 + 1
            sage: t.splitting_field('b')                                                # needs sage.rings.finite_rings
            Finite Field in b of size 7^2

            sage: P.<x> = PolynomialRing(GF(7^3, 'a'))                                  # needs sage.rings.finite_rings
            sage: t = x^2 + 1
            sage: t.splitting_field('b', map=True)                                      # needs sage.rings.finite_rings
            (Finite Field in b of size 7^6,
             Ring morphism:
               From: Finite Field in a of size 7^3
               To:   Finite Field in b of size 7^6
               Defn: a |--> 2*b^4 + 6*b^3 + 2*b^2 + 3*b + 2)

        If the extension is trivial and the generators have the same
        name, the map will be the identity::

            sage: t = 24*x^13 + 2*x^12 + 14
            sage: t.splitting_field('a', map=True)                                      # needs sage.rings.finite_rings
            (Finite Field in a of size 7^3,
             Identity endomorphism of Finite Field in a of size 7^3)

            sage: t = x^56 - 14*x^3
            sage: t.splitting_field('b', map=True)                                      # needs sage.rings.finite_rings
            (Finite Field in b of size 7^3,
             Ring morphism:
             From: Finite Field in a of size 7^3
               To:   Finite Field in b of size 7^3
               Defn: a |--> b)

        .. SEEALSO::

            :func:`sage.rings.number_field.splitting_field.splitting_field` for more examples over number fields

        TESTS::

            sage: K.<a,b> = x.splitting_field()
            Traceback (most recent call last):
            ...
            IndexError: the number of names must equal the number of generators
            sage: polygen(RR).splitting_field('x')
            Traceback (most recent call last):
            ...
            NotImplementedError: splitting_field() is only implemented over number fields and finite fields

            sage: # needs sage.rings.finite_rings
            sage: P.<x> = PolynomialRing(GF(11^5, 'a'))
            sage: t = x^2 + 1
            sage: t.splitting_field('b')
            Finite Field in b of size 11^10
            sage: t = 24*x^13 + 2*x^12 + 14
            sage: t.splitting_field('b')
            Finite Field in b of size 11^30
            sage: t = x^56 - 14*x^3
            sage: t.splitting_field('b')
            Finite Field in b of size 11^130

            sage: # needs sage.rings.finite_rings
            sage: P.<x> = PolynomialRing(GF(19^6, 'a'))
            sage: t = -x^6 + x^2 + 1
            sage: t.splitting_field('b')
            Finite Field in b of size 19^6
            sage: t = 24*x^13 + 2*x^12 + 14
            sage: t.splitting_field('b')
            Finite Field in b of size 19^18
            sage: t = x^56 - 14*x^3
            sage: t.splitting_field('b')
            Finite Field in b of size 19^156

            sage: # needs sage.rings.finite_rings
            sage: P.<x> = PolynomialRing(GF(83^6, 'a'))
            sage: t = 2*x^14 - 5 + 6*x
            sage: t.splitting_field('b')
            Finite Field in b of size 83^84
            sage: t = 24*x^13 + 2*x^12 + 14
            sage: t.splitting_field('b')
            Finite Field in b of size 83^78
            sage: t = x^56 - 14*x^3
            sage: t.splitting_field('b')
            Finite Field in b of size 83^12

            sage: # needs sage.rings.finite_rings
            sage: P.<x> = PolynomialRing(GF(401^13, 'a'))
            sage: t = 2*x^14 - 5 + 6*x
            sage: t.splitting_field('b')                                                # long time (16s)
            Finite Field in b of size 401^104
            sage: t = 24*x^13 + 2*x^12 + 14
            sage: t.splitting_field('b')                                                # long time (39s)
            Finite Field in b of size 401^156
            sage: t = x^56 - 14*x^3
            sage: t.splitting_field('b')                                                # long time (2s)
            Finite Field in b of size 401^52

            sage: R.<x> = QQ[]
            sage: f = x^2 - 2
            sage: f.splitting_field()
            Traceback (most recent call last):
            ...
            TypeError: You must specify the name of the generator.
        """
        if names is None:
            raise TypeError("You must specify the name of the generator.")
        name = normalize_names(1, names)[0]

        from sage.rings.number_field.number_field_base import NumberField
        from sage.rings.finite_rings.finite_field_base import FiniteField

        f = self.monic()            # Given polynomial, made monic
        F = f.parent().base_ring()  # Base field
        if not F.is_field():
            F = F.fraction_field()
            f = self.change_ring(F)

        if isinstance(F, NumberField):
            from sage.rings.number_field.splitting_field import splitting_field
            return splitting_field(f, name, map, **kwds)
        elif isinstance(F, FiniteField):
            degree = lcm([f.degree() for f, _ in self.factor()])
            return F.extension(degree, name, map=map, **kwds)

        raise NotImplementedError("splitting_field() is only implemented over number fields and finite fields")

    def pseudo_quo_rem(self,other):
        r"""
        Compute the pseudo-division of two polynomials.

        INPUT:

        - ``other`` -- a nonzero polynomial

        OUTPUT:

        `Q` and `R` such that `l^{m-n+1} \mathrm{self} = Q \cdot\mathrm{other} + R`
        where `m` is the degree of this polynomial, `n` is the degree of
        ``other``, `l` is the leading coefficient of ``other``. The result is
        such that `\deg(R) < \deg(\mathrm{other})`.

        ALGORITHM:

        Algorithm 3.1.2 in [Coh1993]_.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^4 + 6*x^3 + x^2 - x + 2
            sage: q = 2*x^2 - 3*x - 1
            sage: quo, rem = p.pseudo_quo_rem(q); quo, rem
            (4*x^2 + 30*x + 51, 175*x + 67)
            sage: 2^(4-2+1)*p == quo*q + rem
            True

            sage: S.<T> = R[]
            sage: p = (-3*x^2 - x)*T^3 - 3*x*T^2 + (x^2 - x)*T + 2*x^2 + 3*x - 2
            sage: q = (-x^2 - 4*x - 5)*T^2 + (6*x^2 + x + 1)*T + 2*x^2 - x
            sage: quo, rem = p.pseudo_quo_rem(q); quo, rem
            ((3*x^4 + 13*x^3 + 19*x^2 + 5*x)*T + 18*x^4 + 12*x^3 + 16*x^2 + 16*x,
             (-113*x^6 - 106*x^5 - 133*x^4 - 101*x^3 - 42*x^2 - 41*x)*T
                  - 34*x^6 + 13*x^5 + 54*x^4 + 126*x^3 + 134*x^2 - 5*x - 50)
            sage: (-x^2 - 4*x - 5)^(3-2+1) * p == quo*q + rem
            True
        """
        if other.is_zero():
            raise ZeroDivisionError("Pseudo-division by zero is not possible")

        # if other is a constant, then R = 0 and Q = self * other^(deg(self))
        if other in self._parent.base_ring():
            return (self * other**(self.degree()), self._parent.zero())

        R = self
        B = other
        Q = self._parent.zero()
        e = self.degree() - other.degree() + 1
        d = B.leading_coefficient()

        while not R.degree() < B.degree():
            c = R.leading_coefficient()
            diffdeg = R.degree() - B.degree()
            Q = d*Q + self._parent(c).shift(diffdeg)
            R = d*R - c*B.shift(diffdeg)
            e -= 1

        q = d**e
        return (q*Q,q*R)

    @coerce_binop
    def gcd(self, other):
        """
        Return a greatest common divisor of this polynomial and ``other``.

        INPUT:

        - ``other`` -- a polynomial in the same ring as this polynomial

        OUTPUT:

        A greatest common divisor as a polynomial in the same ring as
        this polynomial. If the base ring is a field, the return value
        is a monic polynomial.

        .. NOTE::

            The actual algorithm for computing greatest common divisors depends
            on the base ring underlying the polynomial ring. If the base ring
            defines a method :meth:`_gcd_univariate_polynomial`, then this method
            will be called (see examples below).

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (2*x^2).gcd(2*x)
            x
            sage: R.zero().gcd(0)
            0
            sage: (2*x).gcd(0)
            x

        One can easily add gcd functionality to new rings by providing a method
        ``_gcd_univariate_polynomial``::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: O = ZZ[-sqrt(5)]
            sage: R.<x> = O[]
            sage: a = O.1
            sage: p = x + a
            sage: q = x^2 - 5
            sage: p.gcd(q)
            Traceback (most recent call last):
            ...
            NotImplementedError: Order of conductor 2 generated by a in Number
            Field in a with defining polynomial x^2 - 5 with a = -2.236067977499790?
            does not provide a gcd implementation for univariate polynomials
            sage: S.<x> = O.number_field()[]
            sage: O._gcd_univariate_polynomial = lambda f, g: R(S(f).gcd(S(g)))
            sage: p.gcd(q)
            x + a
            sage: del O._gcd_univariate_polynomial

        Use multivariate implementation for polynomials over polynomials rings::

            sage: R.<x> = ZZ[]
            sage: S.<y> = R[]
            sage: T.<z> = S[]
            sage: r = 2*x*y + z
            sage: p = r * (3*x*y*z - 1)
            sage: q = r * (x + y + z - 2)
            sage: p.gcd(q)                                                              # needs sage.libs.singular
            z + 2*x*y

            sage: R.<x> = QQ[]
            sage: S.<y> = R[]
            sage: r = 2*x*y + 1
            sage: p = r * (x - 1/2 * y)
            sage: q = r * (x*y^2 - x + 1/3)
            sage: p.gcd(q)                                                              # needs sage.libs.singular
            2*x*y + 1

        TESTS::

            sage: Pol = QQ['x','y']['x']
            sage: Pol.one().gcd(1)
            1
        """
        cdef Polynomial _other = <Polynomial> other
        if _other.is_one():
            return other
        elif self.is_one():
            return self
        flatten = self._parent.flattening_morphism()
        tgt = flatten.codomain()
        if tgt.ngens() > 1 and tgt._has_singular:
            g = flatten(self).gcd(flatten(other))
            return flatten.section()(g)
        try:
            doit = self._parent._base._gcd_univariate_polynomial
        except AttributeError:
            raise NotImplementedError("%s does not provide a gcd implementation for univariate polynomials"%self._parent._base)
        else:
            return doit(self, other)

    @coerce_binop
    def lcm(self, other):
        """
        Let `f` and `g` be two polynomials. Then this function returns the
        monic least common multiple of `f` and `g`.

        TESTS:

        Check that :issue:`32033` has been fixed::

            sage: R.<t> = GF(3)[]
            sage: lcm(R(0), R(0))
            0

        ::

            sage: R.<t> = RR[]
            sage: lcm(R(0), R(0))
            0
        """
        if self.is_zero() or other.is_zero():
            P = self.parent()
            return P.zero()
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q.leading_coefficient())*q

    def _lcm(self, other):
        """
        Let `f` and `g` be two polynomials. Then this function returns the
        monic least common multiple of `f` and `g`.
        """
        if self.is_zero() or other.is_zero():
            P = self.parent()
            return P.zero()
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q.leading_coefficient())*q  # make monic  (~ is inverse in python)

    def is_primitive(self, n=None, n_prime_divs=None):
        r"""
        Return ``True`` if the polynomial is primitive.

        The semantics of "primitive" depend on the polynomial coefficients.

        - (field theory) A polynomial of degree `m` over a finite field
          `\GF{q}` is primitive if it is irreducible and its root in
          `\GF{q^m}` generates the multiplicative group `\GF{q^m}^*`.

        - (ring theory) A polynomial over a ring is primitive if its
          coefficients generate the unit ideal.

        Calling :meth:`is_primitive` on a polynomial over an infinite field will
        raise an error.

        The additional inputs to this function are to speed up computation for
        field semantics (see note).

        INPUT:

        - ``n`` -- (default: ``None``) if provided, should equal
          `q-1` where ``self.parent()`` is the field with `q`
          elements; otherwise it will be computed

        - ``n_prime_divs`` -- (default: ``None``) if provided, should be a list
          of the prime divisors of `n`; otherwise it will be computed

        .. NOTE::

          Computation of the prime divisors of `n` can dominate the running
          time of this method, so performing this computation externally
          (e.g., ``pdivs = n.prime_divisors()``) is a good idea for repeated calls
          to :meth:`is_primitive` for polynomials of the same degree.

          Results may be incorrect if the wrong `n` and/or factorization are
          provided.

        EXAMPLES:

        Field semantics examples.

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(2)['x']
            sage: f = x^4 + x^3 + x^2 + x + 1
            sage: f.is_irreducible(), f.is_primitive()
            (True, False)
            sage: f = x^3 + x + 1
            sage: f.is_irreducible(), f.is_primitive()
            (True, True)
            sage: R.<x> = GF(3)[]
            sage: f = x^3 - x + 1
            sage: f.is_irreducible(), f.is_primitive()
            (True, True)
            sage: f = x^2 + 1
            sage: f.is_irreducible(), f.is_primitive()
            (True, False)
            sage: R.<x> = GF(5)[]
            sage: f = x^2 + x + 1
            sage: f.is_primitive()
            False
            sage: f = x^2 - x + 2
            sage: f.is_primitive()
            True
            sage: x = polygen(QQ); f = x^2 + 1
            sage: f.is_primitive()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_primitive() not defined for polynomials over infinite fields.

        Ring semantics examples.

        ::

            sage: x = polygen(ZZ)
            sage: f = 5*x^2 + 2
            sage: f.is_primitive()
            True
            sage: f = 5*x^2 + 5
            sage: f.is_primitive()
            False

            sage: # needs sage.rings.number_field
            sage: K = NumberField(x^2 + 5, 'a')
            sage: R = K.ring_of_integers()
            sage: a = R.gen(1)
            sage: a^2
            -5
            sage: f = a*x + 2
            sage: f.is_primitive()
            True
            sage: f = (1+a)*x + 2
            sage: f.is_primitive()
            False

            sage: x = polygen(Integers(10))
            sage: f = 5*x^2 + 2
            sage: #f.is_primitive()  #BUG:: elsewhere in Sage, should return True
            sage: f = 4*x^2 + 2
            sage: #f.is_primitive()  #BUG:: elsewhere in Sage, should return False

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(2)['x']
            sage: f = x^4 + x^3 + x^2 + x + 1
            sage: f.is_primitive(15)
            False
            sage: f.is_primitive(15, [3,5])
            False
            sage: f.is_primitive(n_prime_divs=[3,5])
            False
            sage: f = x^3 + x + 1
            sage: f.is_primitive(7, [7])
            True
            sage: R.<x> = GF(3)[]
            sage: f = x^3 - x + 1
            sage: f.is_primitive(26, [2,13])
            True
            sage: f = x^2 + 1
            sage: f.is_primitive(8, [2])
            False
            sage: R.<x> = GF(5)[]
            sage: f = x^2 + x + 1
            sage: f.is_primitive(24, [2,3])
            False
            sage: f = x^2 - x + 2
            sage: f.is_primitive(24, [2,3])
            True
            sage: x = polygen(Integers(103)); f = x^2 + 1
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
            y = self._parent.quo(self).gen()
            from sage.groups.generic import order_from_multiple
            return n == order_from_multiple(y, n, n_prime_divs, operation='*')
        elif isinstance(R, sage.rings.abc.Order):
            K = R.number_field()
            return K.fractional_ideal(self.coefficients()) == K.fractional_ideal(1)
        else:
            return R.ideal(self.coefficients()) == R.ideal(1)

    def is_constant(self):
        """
        Return ``True`` if this is a constant polynomial.

        OUTPUT: boolean; ``True`` if and only if this polynomial is constant

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
        Return ``True`` if ``self`` is a monomial, i.e., a power of the generator.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.is_monomial()
            True
            sage: (x + 1).is_monomial()
            False
            sage: (x^2).is_monomial()
            True
            sage: R(1).is_monomial()
            True

        The coefficient must be 1::

            sage: (2*x^5).is_monomial()
            False

        To allow a non-1 leading coefficient, use :meth:`is_term`::

            sage: (2*x^5).is_term()
            True

        .. warning::

           The definition of :meth:`is_monomial` in Sage up to 4.7.1 was the
           same as :meth:`is_term`, i.e., it allowed a coefficient not equal
           to 1.
        """
        return len(self.exponents()) == 1 and self.leading_coefficient() == 1

    cpdef bint is_term(self) except -1:
        """
        Return ``True`` if this polynomial is a nonzero element of the
        base ring times a power of the variable.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: x.is_term()
            True
            sage: R(0).is_term()
            False
            sage: R(1).is_term()
            True
            sage: (3*x^5).is_term()
            True
            sage: (1 + 3*x^5).is_term()
            False

        To require that the coefficient is 1, use :meth:`is_monomial`
        instead::

            sage: (3*x^5).is_monomial()
            False
        """
        return len(self.exponents()) == 1

    def root_field(self, names, check_irreducible=True):
        """
        Return the field generated by the roots of the irreducible
        polynomial ``self``. The output is either a number field, relative
        number field, a quotient of a polynomial ring over a field, or the
        fraction field of the base ring.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('a')                                                     # needs sage.rings.number_field
            Number Field in a with defining polynomial x^3 + x + 17

        ::

            sage: R.<x> = QQ['x']
            sage: f = x - 3
            sage: f.root_field('b')                                                     # needs sage.rings.number_field
            Rational Field

        ::

            sage: R.<x> = ZZ['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('b')                                                     # needs sage.rings.number_field
            Number Field in b with defining polynomial x^3 + x + 17

        ::

            sage: # needs sage.rings.number_field
            sage: y = QQ['x'].0
            sage: L.<a> = NumberField(y^3 - 2)
            sage: R.<x> = L['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('c')
            Number Field in c with defining polynomial x^3 + x + 17 over its base field

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = PolynomialRing(GF(9, 'a'))
            sage: f = x^3 + x^2 + 8
            sage: K.<alpha> = f.root_field(); K
            Univariate Quotient Polynomial Ring in alpha
             over Finite Field in a of size 3^2 with modulus x^3 + x^2 + 2
            sage: alpha^2 + 1
            alpha^2 + 1
            sage: alpha^3 + alpha^2
            1

        ::

            sage: R.<x> = QQ[]
            sage: f = x^2
            sage: K.<alpha> = f.root_field()                                            # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: polynomial must be irreducible

        TESTS::

            sage: (PolynomialRing(Integers(31), name='x').0 + 5).root_field('a')        # needs sage.rings.finite_rings
            Ring of integers modulo 31
        """
        R = self.base_ring()
        if not R.is_integral_domain():
            raise ValueError("the base ring must be a domain")

        if check_irreducible and not self.is_irreducible():
            raise ValueError("polynomial must be irreducible")

        if self.degree() <= 1:
            return R.fraction_field()

        if isinstance(R, IntegerRing_class):
            from sage.rings.number_field.number_field import NumberField
            return NumberField(self, names)

        from sage.rings.number_field.number_field_base import NumberField as NumberField_base
        if isinstance(R, NumberField_base):
            from sage.rings.number_field.number_field import NumberField
            return NumberField(self, names)

        return R.fraction_field()[self._parent.variable_name()].quotient(self, names)

    def sylvester_matrix(self, right, variable = None):
        """
        Return the Sylvester matrix of ``self`` and ``right``.

        Note that the Sylvester matrix is not defined if one of the polynomials
        is zero.

        INPUT:

        - ``right`` -- a polynomial in the same ring as ``self``
        - ``variable`` -- (optional) included for compatibility with the
          multivariate case only; the variable of the polynomials

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = (6*x + 47) * (7*x^2 - 2*x + 38)
            sage: g = (6*x + 47) * (3*x^3 + 2*x + 1)
            sage: M = f.sylvester_matrix(g); M                                          # needs sage.modules
            [  42  317  134 1786    0    0    0]
            [   0   42  317  134 1786    0    0]
            [   0    0   42  317  134 1786    0]
            [   0    0    0   42  317  134 1786]
            [  18  141   12  100   47    0    0]
            [   0   18  141   12  100   47    0]
            [   0    0   18  141   12  100   47]

        If the polynomials share a non-constant common factor then the
        determinant of the Sylvester matrix will be zero::

            sage: M.determinant()                                                       # needs sage.modules
            0

        If ``self`` and ``right`` are polynomials of positive degree, the determinant
        of the Sylvester matrix is the resultant of the polynomials.::

            sage: h1 = R._random_nonzero_element()
            sage: h2 = R._random_nonzero_element()
            sage: M1 = h1.sylvester_matrix(h2)                                          # needs sage.modules
            sage: M1.determinant() == h1.resultant(h2)                                  # needs sage.libs.pari sage.modules
            True

        The rank of the Sylvester matrix is related to the degree of the
        gcd of ``self`` and ``right``::

            sage: f.gcd(g).degree() == f.degree() + g.degree() - M.rank()               # needs sage.modules
            True
            sage: h1.gcd(h2).degree() == h1.degree() + h2.degree() - M1.rank()          # needs sage.modules
            True

        TESTS:

        The variable is optional, but must be the same in both rings::

            sage: K.<x> = QQ['x']
            sage: f = x + 1
            sage: g = QQ['y']([1, 0, 1])
            sage: f.sylvester_matrix(f, x)                                              # needs sage.modules
            [1 1]
            [1 1]
            sage: f.sylvester_matrix(g, x)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
             'Univariate Polynomial Ring in x over Rational Field' and
             'Univariate Polynomial Ring in y over Rational Field'

        Polynomials must be defined over compatible base rings::

            sage: f = QQ['x']([1, 0, 1])
            sage: g = ZZ['x']([1, 0, 1])
            sage: h = GF(25, 'a')['x']([1, 0, 1])                                       # needs sage.rings.finite_rings
            sage: f.sylvester_matrix(g)                                                 # needs sage.modules
            [1 0 1 0]
            [0 1 0 1]
            [1 0 1 0]
            [0 1 0 1]
            sage: g.sylvester_matrix(h)                                                 # needs sage.modules sage.rings.finite_rings
            [1 0 1 0]
            [0 1 0 1]
            [1 0 1 0]
            [0 1 0 1]
            sage: f.sylvester_matrix(h)                                                 # needs sage.modules sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
             'Univariate Polynomial Ring in x over Rational Field' and
             'Univariate Polynomial Ring in x over Finite Field in a of size 5^2'

        We can compute the sylvester matrix of a univariate and multivariate
        polynomial::

            sage: K.<x,y> = QQ['x,y']
            sage: g = K.random_element()
            sage: f.sylvester_matrix(g) == K(f).sylvester_matrix(g, x)                  # needs sage.modules
            True

        Corner cases::

            sage: # needs sage.modules
            sage: K.<x> = QQ[]
            sage: f = x^2 + 1
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

        if self._parent != right.parent():
            a, b = coercion_model.canonical_coercion(self, right)
            variable = a.parent()(self.variables()[0])
            #We add the variable to cover the case that right is a multivariate
            #polynomial
            return a.sylvester_matrix(b, variable)

        if variable:
            if variable.parent() != self._parent:
                variable = self._parent(variable)

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
                M[r, m - c + offset] = self.get_unsafe(c)
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
        if self.is_zero():
            return self._parent._base.zero()
        # self.degree() >= 0
        return self.get_unsafe(0)

    cpdef Polynomial _new_constant_poly(self, a, Parent P):
        """
        Create a new constant polynomial from a in P, which MUST be an
        element of the base ring of P (this is not checked).

        EXAMPLES::

            sage: R.<w> = PolynomialRing(GF(9, 'a'), sparse=True)                       # needs sage.rings.finite_rings
            sage: a = w._new_constant_poly(0, R); a                                     # needs sage.rings.finite_rings
            0
            sage: a.coefficients()                                                      # needs sage.rings.finite_rings
            []
        """
        t = type(self)
        return t(P, [a] if a else [], check=False)

    def is_monic(self):
        """
        Return ``True`` if this polynomial is monic. The zero polynomial is by
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
        return not self.is_zero() and self[self.degree()] == self.base_ring().one()

    def is_unit(self):
        r"""
        Return ``True`` if this polynomial is a unit.

        EXAMPLES::

            sage: a = Integers(90384098234^3)
            sage: b = a(2*191*236607587)
            sage: b.is_nilpotent()
            True

            sage: # needs sage.libs.pari
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

        TESTS:

        Check that :issue:`18600` is fixed::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: c = x^2^100 + 1
            sage: c.is_unit()
            False
        """
        d = self.degree()
        if d > 0:
            try:
                if self._parent.base_ring().is_integral_domain():
                    return False
            except NotImplementedError:
                pass
            for c in self.coefficients()[1:]:
                if not c.is_nilpotent():
                    return False
        if d == -1:
            return self._parent._base.zero().is_unit()
        return self.get_unsafe(0).is_unit()

    def is_nilpotent(self):
        r"""
        Return ``True`` if this polynomial is nilpotent.

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

        TESTS:

        Check that :issue:`18600` is fixed::

            sage: R.<x> = PolynomialRing(Zmod(4), sparse=True)
            sage: (2*x^2^100 + 2).is_nilpotent()
            True
        """
        for c in self.coefficients():
            if not c.is_nilpotent():
                return False
        return True

    def is_gen(self):
        r"""
        Return ``True`` if this polynomial is the distinguished generator of
        the parent polynomial ring.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R(1).is_gen()
            False
            sage: R(x).is_gen()
            True

        Important - this function doesn't return ``True`` if ``self`` equals the
        generator; it returns ``True`` if ``self`` *is* the generator.

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

    def lc(self):
        """
        Return the leading coefficient of this polynomial.

        OUTPUT: element of the base ring

        This method is the same as :meth:`leading_coefficient`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (-2/5)*x^3 + 2*x - 1/3
            sage: f.lc()
            -2/5
        """
        return self[self.degree()]

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

    def lm(self):
        """
        Return the leading monomial of this polynomial.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (-2/5)*x^3 + 2*x - 1/3
            sage: f.lm()
            x^3
            sage: R(5).lm()
            1
            sage: R(0).lm()
            0
            sage: R(0).lm().parent() is R
            True
        """
        if self.degree() < 0:
            return self
        output = [self.base_ring().zero()] * self.degree() + [self.base_ring().one()]
        return self._new_generic(output)

    def lt(self):
        """
        Return the leading term of this polynomial.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (-2/5)*x^3 + 2*x - 1/3
            sage: f.lt()
            -2/5*x^3
            sage: R(5).lt()
            5
            sage: R(0).lt()
            0
            sage: R(0).lt().parent() is R
            True
        """
        return self.lc() * self.lm()

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
        R = self._parent
        if a.parent() != R.base_ring():
            S = R.base_extend(a.parent())
            return a*S(self)
        else:
            return a*self

    def coefficients(self, sparse=True):
        r"""
        Return the coefficients of the monomials appearing in ``self``.

        If ``sparse=True`` (the default), it returns only the nonzero coefficients.
        Otherwise, it returns the same value as ``self.list()``.
        (In this case, it may be slightly faster to invoke ``self.list()`` directly.)
        In either case, the coefficients are ordered by increasing degree.

        EXAMPLES::

            sage: _.<x> = PolynomialRing(ZZ)
            sage: f = 3*x^4 + 2*x^2 + 1
            sage: f.coefficients()
            [1, 2, 3]
            sage: f.coefficients(sparse=False)
            [1, 0, 2, 0, 3]
        """
        zero = self._parent.base_ring().zero()
        if sparse:
            return [c for c in self.list() if c != zero]
        else:
            return self.list()

    def global_height(self, prec=None):
        """
        Return the (projective) global height of the polynomial.

        This returns the absolute logarithmic height of the coefficients
        thought of as a projective point.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default :class:`RealField` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 3*x^3 + 2*x^2 + x
            sage: exp(f.global_height())                                                # needs sage.symbolic
            3.00000000000000

        Scaling should not change the result::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 1/25*x^2 + 25/3*x + 1
            sage: f.global_height()                                                     # needs sage.symbolic
            6.43775164973640
            sage: g = 100 * f
            sage: g.global_height()                                                     # needs sage.symbolic
            6.43775164973640

        ::

            sage: R.<x> = PolynomialRing(QQbar)                                         # needs sage.rings.number_field
            sage: f = QQbar(i)*x^2 + 3*x                                                # needs sage.rings.number_field
            sage: f.global_height()                                                     # needs sage.rings.number_field
            1.09861228866811

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<k> = NumberField(x^2 + 5)
            sage: T.<t> = PolynomialRing(K)
            sage: f = 1/1331 * t^2 + 5 * t + 7
            sage: f.global_height()
            9.13959596745043

        ::

            sage: R.<x> = QQ[]
            sage: f = 1/123*x^2 + 12
            sage: f.global_height(prec=2)                                               # needs sage.symbolic
            8.0

        ::

            sage: R.<x> = QQ[]
            sage: f = 0*x
            sage: f.global_height()                                                     # needs sage.rings.real_mpfr
            0.000000000000000
        """
        if prec is None:
            prec = 53

        if self.is_zero():
            from sage.rings.real_mpfr import RealField
            return RealField(prec).zero()

        from sage.categories.number_fields import NumberFields

        K = self.base_ring()
        if K in NumberFields() or isinstance(K, sage.rings.abc.Order) or isinstance(K, IntegerRing_class):
            from sage.schemes.projective.projective_space import ProjectiveSpace
            P = ProjectiveSpace(K, self.number_of_terms() - 1)
            return P.point(self.coefficients()).global_height(prec=prec)
        elif isinstance(K, sage.rings.abc.AlgebraicField):
            from sage.rings.qqbar import number_field_elements_from_algebraics
            from sage.schemes.projective.projective_space import ProjectiveSpace

            K_pre, P, phi = number_field_elements_from_algebraics(self.coefficients())
            Pr = ProjectiveSpace(K_pre, len(P) - 1)
            return Pr.point(P).global_height(prec=prec)
        raise TypeError("Must be over a Numberfield or a Numberfield Order.")

    def local_height(self, v, prec=None):
        """
        Return the maximum of the local height of the coefficients of
        this polynomial.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring

        - ``prec`` -- desired floating point precision (default:
          default :class:`RealField` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 1/1331*x^2 + 1/4000*x
            sage: f.local_height(1331)                                                  # needs sage.rings.real_mpfr
            7.19368581839511

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<k> = NumberField(x^2 - 5)
            sage: T.<t> = K[]
            sage: I = K.ideal(3)
            sage: f = 1/3*t^2 + 3
            sage: f.local_height(I)
            1.09861228866811

        ::

            sage: R.<x> = QQ[]
            sage: f = 1/2*x^2 + 2
            sage: f.local_height(2, prec=2)                                             # needs sage.rings.real_mpfr
            0.75
        """
        from sage.categories.number_fields import NumberFields

        if prec is None:
            prec = 53

        K = FractionField(self.base_ring())
        if K not in NumberFields():
            raise TypeError("must be over a Numberfield or a Numberfield order")

        return max([K(c).local_height(v, prec=prec) for c in self.coefficients()])

    def local_height_arch(self, i, prec=None):
        """
        Return the maximum of the local height at the ``i``-th infinite place
        of the coefficients of this polynomial.

        INPUT:

        - ``i`` -- integer

        - ``prec`` -- desired floating point precision (default:
          default :class:`RealField` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 210*x^2
            sage: f.local_height_arch(0)                                                # needs sage.rings.real_mpfr
            5.34710753071747

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<k> = NumberField(x^2 - 5)
            sage: T.<t> = K[]
            sage: f = 1/2*t^2 + 3
            sage: f.local_height_arch(1, prec=52)
            1.09861228866811

        ::

            sage: R.<x> = QQ[]
            sage: f = 1/2*x^2 + 3
            sage: f.local_height_arch(0, prec=2)                                        # needs sage.rings.real_mpfr
            1.0
        """
        from sage.categories.number_fields import NumberFields

        if prec is None:
            prec = 53

        K = FractionField(self.base_ring())
        if K not in NumberFields():
            return TypeError("must be over a Numberfield or a Numberfield Order")

        if K == QQ:
            return max([K(c).local_height_arch(prec=prec) for c in self.coefficients()])
        return max([K(c).local_height_arch(i, prec=prec) for c in self.coefficients()])

    def exponents(self):
        """
        Return the exponents of the monomials appearing in ``self``.

        EXAMPLES::

            sage: _.<x> = PolynomialRing(ZZ)
            sage: f = x^4 + 2*x^2 + 1
            sage: f.exponents()
            [0, 2, 4]

        TESTS::

            sage: a = RIF['x'](1/3)                                                     # needs sage.rings.real_interval_field
            sage: (a - a).exponents()                                                   # needs sage.rings.real_interval_field
            [0]
        """
        cdef Py_ssize_t i
        return [i for i, c in enumerate(self.list(copy=False)) if c]

    cpdef list list(self, bint copy=True):
        """
        Return a new copy of the list of the underlying elements of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (-2/5)*x^3 + 2*x - 1/3
            sage: v = f.list(); v
            [-1/3, 2, 0, -2/5]

        Note that ``v`` is a list, it is mutable, and each call to the :meth:`list`
        method returns a new list::

            sage: type(v)
            <... 'list'>
            sage: v[0] = 5
            sage: f.list()
            [-1/3, 2, 0, -2/5]

        Here is an example with a generic polynomial ring::

            sage: R.<x> = QQ[]
            sage: S.<y> = R[]
            sage: f = y^3 + x*y - 3*x; f
            y^3 + x*y - 3*x
            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
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
        Return list of coefficients of ``self`` up to (but not including) `q^n`.

        Includes 0s in the list on the right so that the list has length `n`.

        INPUT:

        - ``n`` -- (default: ``None``) if given, an integer that is at least 0

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

        TESTS:

        Check that :issue:`18600` is fixed::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^2^100 + x^8 - 1).padded_list(10)
            [-1, 0, 0, 0, 0, 0, 0, 0, 1, 0]
        """
        if n is None:
            return self.list()
        if n < 0:
            raise ValueError("n must be at least 0")
        if self.degree() < n:
            v = self.list()
            z = self._parent.base_ring().zero()
            return v + [z]*(n - len(v))
        else:
            return self[:int(n)].padded_list(n)

    def monomial_coefficient(self, m):
        """
        Return the coefficient in the base ring of the monomial ``m`` in
        ``self``, where ``m`` must have the same parent as ``self``.

        INPUT:

        - ``m`` -- a monomial

        OUTPUT: coefficient in base ring

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: f = 2 * x
            sage: c = f.monomial_coefficient(x); c
            2
            sage: c.parent()
            Rational Field

            sage: f = x^9 - 1/2*x^2 + 7*x + 5/11
            sage: f.monomial_coefficient(x^9)
            1
            sage: f.monomial_coefficient(x^2)
            -1/2
            sage: f.monomial_coefficient(x)
            7
            sage: f.monomial_coefficient(x^0)
            5/11
            sage: f.monomial_coefficient(x^3)
            0
        """
        if not m.parent() is self._parent:
            raise TypeError("monomial must have same parent as self.")

        d = m.degree()
        coeffs = self.list()
        if 0 <= d < len(coeffs):
            return coeffs[d]
        else:
            return self._parent.base_ring().zero()

    def monomials(self):
        """
        Return the list of the monomials in ``self`` in a decreasing order of their degrees.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: f = x^2 + (2/3)*x + 1
            sage: f.monomials()
            [x^2, x, 1]
            sage: f = P(3/2)
            sage: f.monomials()
            [1]
            sage: f = P(0)
            sage: f.monomials()
            []
            sage: f = x
            sage: f.monomials()
            [x]
            sage: f = - 1/2*x^2 + x^9 + 7*x + 5/11
            sage: f.monomials()
            [x^9, x^2, x, 1]

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<rho> = NumberField(x**2 + 1)
            sage: R.<y> = QQ[]
            sage: p = rho * y
            sage: p.monomials()
            [y]
        """
        if self.is_zero():
            return []
        v = self._parent.gen()
        zero = self._parent.base_ring().zero()
        coeffs = self.list()
        return [v**i for i in range(self.degree(), -1, -1) if coeffs[i] != zero]

    def newton_raphson(self, n, x0):
        """
        Return a list of n iterative approximations to a root of this
        polynomial, computed using the Newton-Raphson method.

        The Newton-Raphson method is an iterative root-finding algorithm.
        For `f(x)` a polynomial, as is the case here, this is essentially the
        same as Horner's method.

        INPUT:

        - ``n`` -- integer; the number of iterations

        - ``x0`` -- an initial guess `x_0`

        OUTPUT: list of numbers hopefully approximating a root of `f(x)=0`

        If one of the iterates is a critical point of `f`, a
        :exc:`ZeroDivisionError` exception is raised.

        EXAMPLES::

            sage: x = PolynomialRing(RealField(), 'x').gen()                            # needs sage.rings.real_mpfr
            sage: f = x^2 - 2                                                           # needs sage.rings.real_mpfr
            sage: f.newton_raphson(4, 1)                                                # needs sage.rings.real_mpfr
            [1.50000000000000, 1.41666666666667, 1.41421568627451, 1.41421356237469]

        AUTHORS:

        - David Joyner and William Stein (2005-11-28)
        """
        n = sage.rings.integer.Integer(n)
        df = self.derivative()
        K = self._parent.base_ring()
        a = K(x0)
        L = []
        for i in range(n):
            a -= self(a) / df(a)
            L.append(a)
        return L

    def polynomial(self, var):
        r"""
        Let ``var`` be one of the variables of the parent of ``self``. This returns
        ``self`` viewed as a univariate polynomial in ``var`` over the polynomial
        ring generated by all the other variables of the parent.

        For univariate polynomials, if ``var`` is the generator of the parent
        ring, we return this polynomial, otherwise raise an error.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (x + 1).polynomial(x)
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

    def newton_slopes(self, p, lengths=False):
        """
        Return the `p`-adic slopes of the Newton polygon of ``self``,
        when this makes sense.

        OUTPUT:

        If ``lengths`` is ``False``, a list of rational numbers. If ``lengths`` is
        ``True``, a list of couples `(s,l)` where `s` is the slope and `l` the
        length of the corresponding segment in the Newton polygon.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: f = x^3 + 2
            sage: f.newton_slopes(2)                                                    # needs sage.libs.pari
            [1/3, 1/3, 1/3]
            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^5 + 6*x^2 + 4
            sage: p.newton_slopes(2)                                                    # needs sage.libs.pari
            [1/2, 1/2, 1/3, 1/3, 1/3]
            sage: p.newton_slopes(2, lengths=True)
            [(1/2, 2), (1/3, 3)]
            sage: (x^2^100 + 27).newton_slopes(3, lengths=True)
            [(3/1267650600228229401496703205376, 1267650600228229401496703205376)]

        ALGORITHM: Uses PARI if ``lengths`` is ``False``.
        """
        if not lengths:
            f = self.__pari__()
            v = list(f.newtonpoly(p))
            return [sage.rings.rational.Rational(x) for x in v]

        e = self.exponents()
        c = self.coefficients()
        if len(e) == 0: return []
        if len(e) == 1:
            if e[0] == 0: return []
            else:         return [(infinity.infinity, e[0])]

        if e[0] == 0: slopes = []
        else:         slopes = [(infinity.infinity, e[0])]

        points = [(e[0], c[0].valuation(p)), (e[1], c[1].valuation(p))]
        slopes.append((-(c[1].valuation(p)-c[0].valuation(p))/(e[1] - e[0]), e[1]-e[0]))
        for i in range(2, len(e)):
            v = c[i].valuation(p)
            s = -(v-points[-1][1])/(e[i]-points[-1][0])
            while slopes and s >= slopes[-1][0]:
                slopes = slopes[:-1]
                points = points[:-1]
                s = -(v-points[-1][1])/(e[i]-points[-1][0])
            slopes.append((s,e[i]-points[-1][0]))
            points.append((e[i],v))

        return slopes

    def dispersion_set(self, other=None):
        r"""
        Compute the dispersion set of two polynomials.

        The dispersion set of `f` and `g` is the set of nonnegative integers
        `n` such that `f(x + n)` and `g(x)` have a nonconstant common factor.

        When ``other`` is ``None``, compute the auto-dispersion set of
        ``self``, i.e., its dispersion set with itself.

        ALGORITHM:

        See Section 4 of Man & Wright [MW1994]_.

        .. SEEALSO:: :meth:`dispersion`

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x.dispersion_set(x + 1)                                               # needs sage.libs.pari
            [1]
            sage: (x + 1).dispersion_set(x)                                             # needs sage.libs.pari
            []

            sage: pol = x^3 + x - 7
            sage: (pol*pol(x+3)^2).dispersion_set()                                     # needs sage.libs.pari
            [0, 3]
        """
        other = self if other is None else self._parent.coerce(other)
        x = self._parent.gen()
        dispersions = set()
        for p, _ in self.factor():
            # need both due to the semantics of is_primitive() over fields
            assert p.is_monic() or p.is_primitive()
            for q, _ in other.factor():
                m, n = p.degree(), q.degree()
                assert q.is_monic() or q.is_primitive()
                if m != n or p[n] != q[n]:
                    continue
                alpha = (q[n-1] - p[n-1])/(n*p[n])
                if alpha.is_integer(): # ZZ() might work for non-integers...
                    alpha = ZZ(alpha)
                else:
                    continue
                if alpha < 0 or alpha in dispersions:
                    continue
                if n >= 1 and p(x + alpha) != q:
                    continue
                dispersions.add(alpha)
        return list(dispersions)

    def dispersion(self, other=None):
        r"""
        Compute the dispersion of a pair of polynomials.

        The dispersion of `f` and `g` is the largest nonnegative integer `n`
        such that `f(x + n)` and `g(x)` have a nonconstant common factor.

        When ``other`` is ``None``, compute the auto-dispersion of ``self``,
        i.e., its dispersion with itself.

        .. SEEALSO:: :meth:`dispersion_set`

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x.dispersion(x + 1)                                                   # needs sage.libs.pari
            1
            sage: (x + 1).dispersion(x)                                                 # needs sage.libs.pari
            -Infinity

            sage: # needs sage.libs.pari sage.rings.number_field sage.symbolic
            sage: Pol.<x> = QQbar[]
            sage: pol = Pol([sqrt(5), 1, 3/2])
            sage: pol.dispersion()
            0
            sage: (pol*pol(x+3)).dispersion()
            3
        """
        dispersions = self.dispersion_set(other)
        return max(dispersions) if len(dispersions) > 0 else infinity.minus_infinity

    #####################################################################
    # Conversions to other systems
    #####################################################################
    def __pari__(self):
        r"""
        Return polynomial as a PARI object.

        Sage does not handle PARI's variable ordering requirements
        gracefully at this time. In practice, this means that the variable
        ``x`` needs to be the topmost variable, as in the
        example.

        EXAMPLES::

            sage: f = QQ['x']([0,1,2/3,3])
            sage: pari(f)                                                               # needs sage.libs.pari
            3*x^3 + 2/3*x^2 + x

        ::

            sage: S.<a> = QQ['a']
            sage: R.<x> = S['x']
            sage: f = R([0, a]) + R([0, 0, 2/3])
            sage: pari(f)                                                               # needs sage.libs.pari
            2/3*x^2 + a*x

        Polynomials over a number field work, provided that the variable is
        called 'x'::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<b> = NumberField(x^2 + x + 1)
            sage: R.<x> = PolynomialRing(K)
            sage: pol = (b + x)^3; pol
            x^3 + 3*b*x^2 + (-3*b - 3)*x + 1
            sage: pari(pol)                                                             # needs sage.libs.pari
            Mod(1, y^2 + y + 1)*x^3 + Mod(3*y, y^2 + y + 1)*x^2
             + Mod(-3*y - 3, y^2 + y + 1)*x + Mod(1, y^2 + y + 1)

        TESTS:

        Unfortunately, variable names matter::

            sage: R.<x, y> = QQ[]
            sage: S.<a> = R[]
            sage: f = x^2 + a; g = y^3 + a
            sage: pari(f)                                                               # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            PariError: incorrect priority in gtopoly: variable x <= a

        Stacked polynomial rings, first with a univariate ring on the
        bottom::

            sage: S.<a> = QQ['a']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 2*x)                                                       # needs sage.libs.pari
            x^2 + 2*x
            sage: pari(a*x + 2*x^3)                                                     # needs sage.libs.pari
            2*x^3 + a*x

        Stacked polynomial rings, second with a multivariate ring on the
        bottom::

            sage: S.<a, b> = ZZ['a', 'b']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 2*x)                                                       # needs sage.libs.pari
            x^2 + 2*x
            sage: pari(a*x + 2*b*x^3)                                                   # needs sage.libs.pari
            2*b*x^3 + a*x

        Stacked polynomial rings with exotic base rings::

            sage: S.<a, b> = GF(7)['a', 'b']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 9*x)                                                       # needs sage.libs.pari
            x^2 + 2*x
            sage: pari(a*x + 9*b*x^3)                                                   # needs sage.libs.pari
            2*b*x^3 + a*x

        ::

            sage: S.<a> = Integers(8)['a']
            sage: R.<x> = S['x']
            sage: pari(x^2 + 2*x)                                                       # needs sage.libs.pari
            Mod(1, 8)*x^2 + Mod(2, 8)*x
            sage: pari(a*x + 10*x^3)                                                    # needs sage.libs.pari
            Mod(2, 8)*x^3 + Mod(1, 8)*a*x
        """
        return self._pari_with_name(self._parent.variable_name())

    def _pari_or_constant(self, name=None):
        r"""
        Convert ``self`` to PARI.  This behaves identical to :meth:`__pari__`
        or :meth:`_pari_with_name` except for constant polynomials:
        then the constant is returned instead of a constant polynomial.

        INPUT:

        - ``name`` -- (default: ``None``) variable name.  If not given, use
          ``self.parent().variable_name()``.  This argument is irrelevant
          for constant polynomials.

        EXAMPLES::

            sage: # needs sage.libs.pari
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
            sage: pol.__pari__().type()
            't_POL'
            sage: PolynomialRing(IntegerModRing(101), 't')()._pari_or_constant()
            Mod(0, 101)
        """
        if self.is_constant():
            return self.get_coeff_c(0).__pari__()
        if name is None:
            name = self._parent.variable_name()
        return self._pari_with_name(name)

    def _pari_with_name(self, name='x'):
        r"""
        Return polynomial as a PARI object with topmost variable
        ``name``.  By default, use 'x' for the variable name.

        For internal use only.

        EXAMPLES::

            sage: R.<a> = PolynomialRing(ZZ)
            sage: (2*a^2 + a)._pari_with_name()                                         # needs sage.libs.pari
            2*x^2 + x
            sage: (2*a^2 + a)._pari_with_name('y')                                      # needs sage.libs.pari
            2*y^2 + y
        """
        vals = [x.__pari__() for x in self.list()]
        return pari(vals).Polrev(name)

    def _pari_init_(self):
        return repr(self.__pari__())

    def _magma_init_(self, magma):
        """
        Return a string that evaluates in Magma to this polynomial.

        EXAMPLES::

            sage: # optional - magma
            sage: magma = Magma()  # new session
            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: f._magma_init_(magma)
            '_sage_[...]![5,-17,0,1]'
            sage: g = magma(f); g
            y^3 - 17*y + 5

        Note that in Magma there is only one polynomial ring over each
        base, so if we make the polynomial ring over ZZ with variable
        `z`, then this changes the variable name of the polynomial
        we already defined::

            sage: R.<z> = ZZ[]
            sage: magma(R)                      # optional - magma
            Univariate Polynomial Ring in z over Integer Ring
            sage: g                             # optional - magma
            z^3 - 17*z + 5

        In Sage the variable name does not change::

            sage: f                             # optional - magma
            y^3 - 17*y + 5

        A more complicated nested example::

            sage: # optional - magma, needs sage.rings.finite_rings
            sage: k.<a> = GF(9); R.<s,t> = k[]; S.<W> = R[]
            sage: magma(a*W^20 + s*t/a)
            a*W^20 + a^7*s*t
        """
        # Get a reference to Magma version of parent.
        R = magma(self._parent)
        # Get list of coefficients.
        v = ','.join(a._magma_init_(magma) for a in self.list())
        return '%s![%s]' % (R.name(), v)

    def _gap_(self, gap):
        """
        Return this polynomial in GAP.

        INPUT:

        - ``gap`` -- a GAP or libgap instance

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: g = gap(f); g   # indirect doctest
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
            sage: libgap(z^2 + z)
            z^2+z

        Coefficients in a finite field::

            sage: # needs sage.libs.gap
            sage: R.<y> = GF(7)[]
            sage: f = y^3 - 17*y + 5
            sage: g = gap(f); g
            y^3+Z(7)^4*y+Z(7)^5
            sage: h = libgap(f); h
            y^3+Z(7)^4*y+Z(7)^5
            sage: g.Factors()
            [ y+Z(7)^0, y+Z(7)^0, y+Z(7)^5 ]
            sage: h.Factors()
            [ y+Z(7)^0, y+Z(7)^0, y+Z(7)^5 ]
            sage: f.factor()
            (y + 5) * (y + 1)^2
        """
        R = gap(self._parent)
        var = list(R.IndeterminatesOfPolynomialRing())[0]
        return self(var)

    def _libgap_(self):
        r"""
        TESTS::

            sage: R.<x> = ZZ[]
            sage: libgap(-x^3 + 3*x)   # indirect doctest                               # needs sage.libs.gap
            -x^3+3*x
            sage: libgap(R.zero())     # indirect doctest                               # needs sage.libs.gap
            0
        """
        from sage.libs.gap.libgap import libgap
        return self._gap_(libgap)

    def _giac_init_(self):
        r"""
        Return a Giac string representation of this polynomial.

        TESTS::

            sage: # needs sage.libs.giac
            sage: R.<x> = GF(101)['e,i'][]
            sage: f = R('e*i') * x + x^2
            sage: f._giac_init_()
            '((1)*1)*sageVARx^2+((1)*sageVARe*sageVARi)*sageVARx'
            sage: giac(f)
            sageVARx^2+sageVARe*sageVARi*sageVARx
            sage: giac(R.zero())
            0
        """
        g = 'sageVAR' + self.variable_name()
        s = '+'.join('(%s)*%s' % (self.monomial_coefficient(m)._giac_init_(),
                                  m._repr(g)) for m in self.monomials())
        return s if s else '0'

    ######################################################################

    @coerce_binop
    def resultant(self, other):
        r"""
        Return the resultant of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

        ALGORITHM:

        Uses PARI's function :pari:`polresultant`.  For base rings that
        are not supported by PARI, the resultant is computed as the
        determinant of the Sylvester matrix.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r                                                 # needs sage.libs.pari
            -8
            sage: r.parent() is QQ                                                      # needs sage.libs.pari
            True

        We can compute resultants over univariate and multivariate
        polynomial rings::

            sage: R.<a> = QQ[]
            sage: S.<x> = R[]
            sage: f = x^2 + a; g = x^3 + a
            sage: r = f.resultant(g); r                                                 # needs sage.libs.pari
            a^3 + a^2
            sage: r.parent() is R                                                       # needs sage.libs.pari
            True

        ::

            sage: R.<a, b> = QQ[]
            sage: S.<x> = R[]
            sage: f = x^2 + a; g = x^3 + b
            sage: r = f.resultant(g); r                                                 # needs sage.libs.pari
            a^3 + b^2
            sage: r.parent() is R                                                       # needs sage.libs.pari
            True

        TESTS::

            sage: R.<x, y> = QQ[]
            sage: S.<a> = R[]
            sage: f = x^2 + a; g = y^3 + a
            sage: h = f.resultant(g); h                                                 # needs sage.libs.pari sage.modules
            y^3 - x^2
            sage: h.parent() is R                                                       # needs sage.libs.pari sage.modules
            True

        Check that :issue:`13672` is fixed::

            sage: R.<t> = GF(2)[]
            sage: S.<x> = R[]
            sage: f = (t^2 + t)*x + t^2 + t
            sage: g = (t + 1)*x + t^2
            sage: f.resultant(g)                                                        # needs sage.libs.pari
            t^4 + t

        Check that :issue:`15061` is fixed::

            sage: R.<T> = PowerSeriesRing(QQ)
            sage: F = R([1,1],2)
            sage: RP.<x> = PolynomialRing(R)
            sage: P = x^2 - F
            sage: P.resultant(P.derivative())                                           # needs sage.libs.pari
            -4 - 4*T + O(T^2)

        Check that :issue:`16360` is fixed::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: y.resultant(y + x)                                                    # needs sage.libs.pari sage.modules
            x

            sage: # needs sage.libs.singular
            sage: K.<a> = FunctionField(QQ)
            sage: R.<b> = K[]
            sage: L.<b> = K.extension(b^2 - a)
            sage: R.<x> = L[]
            sage: f = x^2 - a
            sage: g = x - b
            sage: f.resultant(g)                                                        # needs sage.libs.pari
            0

        Check that :issue:`17817` is fixed::

            sage: A.<a,b,c> = Frac(PolynomialRing(QQ,'a,b,c'))
            sage: B.<d,e,f> = PolynomialRing(A,'d,e,f')
            sage: R.<x> = PolynomialRing(B,'x')
            sage: S.<y> = PolynomialRing(R,'y')
            sage: p = ((1/b^2*d^2+1/a)*x*y^2+a*b/c*y+e+x^2)
            sage: q = -4*c^2*y^3+1
            sage: p.resultant(q)                                                        # needs sage.libs.pari sage.modules
            (16*c^4)*x^6 + (48*c^4)*e*x^4 + (1/(b^6)*d^6 + 3/(a*b^4)*d^4
            + (-12*a^3*b*c + 3)/(a^2*b^2)*d^2 + (-12*a^3*b*c + 1)/(a^3))*x^3
            + (48*c^4)*e^2*x^2 + ((-12*a*c)/b*d^2*e + (-12*b*c)*e)*x + (16*c^4)*e^3 + (4*a^3*b^3)/c


        Test for :issue:`10978`::

            sage: # needs sage.libs.pari sage.rings.complex_double sage.symbolic
            sage: R.<x> = PolynomialRing(CDF)
            sage: f = R(1 - I*x + (0.5)*x^2 + (1.7)*x^3)
            sage: g = f.derivative()
            sage: f.resultant(g)                                                        # needs sage.modules
            133.92599999999996 + 37.56999999999999*I
        """
        variable = self.variable_name()
        try:
            res = self.__pari__().polresultant(other, variable)
            return self._parent.base_ring()(res)
        except (TypeError, ValueError, PariError, NotImplementedError):
            return self.sylvester_matrix(other).det()

    @coerce_binop
    def subresultants(self, other):
        r"""
        Return the nonzero subresultant polynomials of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: list of polynomials in the same ring as ``self``

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^8 + x^6 - 3*x^4 - 3*x^3 + 8*x^2 + 2*x - 5
            sage: g = 3*x^6 + 5*x^4 - 4*x^2 - 9*x + 21
            sage: f.subresultants(g)
            [260708,
             9326*x - 12300,
             169*x^2 + 325*x - 637,
             65*x^2 + 125*x - 245,
             25*x^4 - 5*x^2 + 15,
             15*x^4 - 3*x^2 + 9]

        ALGORITHM:

        We use the schoolbook algorithm with Lazard's optimization described in [Duc1998]_

        REFERENCES:

        :wikipedia:`Polynomial_greatest_common_divisor#Subresultants`
        """
        P, Q = self, other
        if P.degree() < Q.degree():
            P, Q = Q, P
        S = []
        s = Q.leading_coefficient()**(P.degree()-Q.degree())
        A = Q
        B = P.pseudo_quo_rem(-Q)[1]
        ring = self.parent()
        while True:
            d = A.degree()
            e = B.degree()
            if B.is_zero():
                return S
            S = [ring(B)] + S
            delta = d - e
            if delta > 1:
                if len(S) > 1:
                    n = S[1].degree() - S[0].degree() - 1
                    if n == 0:
                        C = S[0]
                    else:
                        x = S[0].leading_coefficient()
                        y = S[1].leading_coefficient()
                        a = 1 << (int(n).bit_length()-1)
                        c = x
                        n = n - a
                        while a > 1:
                            a /= 2
                            c = c**2 / y
                            if n >= a:
                                c = c * x / y
                                n = n - a
                        C = c * S[0] / y
                else:
                    C = B.leading_coefficient()**(delta-1) * B / s**(delta-1)
                S = [ring(C)] + S
            else:
                C = B
            if e == 0:
                return S
            B = A.pseudo_quo_rem(-B)[1] / (s**delta * A.leading_coefficient())
            A = C
            s = A.leading_coefficient()

    def composed_op(p1, p2, op, algorithm=None, monic=False):
        r"""
        Return the composed sum, difference, product or quotient of this
        polynomial with another one.

        In the case of two monic polynomials `p_1` and `p_2` over an integral
        domain, the composed sum, difference, etc. are given by

        .. MATH::

            \prod_{p_1(a)=p_2(b)=0}(x - (a \ast b)), \qquad
            \ast ∈ \{ +, -, ×, / \}

        where the roots `a` and `b` are to be considered in the algebraic
        closure of the fraction field of the coefficients and counted with
        multiplicities. If the polynomials are not monic this quantity is
        multiplied by `\alpha_1^{\deg(p_2)} \alpha_2^{\deg(p_1)}` where
        `\alpha_1` and `\alpha_2` are the leading coefficients of `p_1` and
        `p_2` respectively.

        INPUT:

        - ``p2`` -- univariate polynomial belonging to the same polynomial ring
          as this polynomial

        - ``op`` -- ``operator.OP`` where ``OP=add`` or ``sub`` or ``mul`` or
          ``truediv``

        - ``algorithm`` -- can be ``'resultant'`` or ``'BFSS'``;
          by default the former is used when the polynomials have few nonzero
          coefficients and small degrees or if the base ring is not `\ZZ` or
          `\QQ`. Otherwise the latter is used.

        - ``monic`` -- whether to return a monic polynomial. If ``True`` the
          coefficients of the result belong to the fraction field of the
          coefficients.

        ALGORITHM:

        The computation is straightforward using resultants. Indeed for the
        composed sum it would be `Res_y(p_1(x-y), p_2(y))`. However, the method
        from [BFSS2006]_ using series expansions is asymptotically much faster.

        Note that the algorithm ``BFSS`` with polynomials with coefficients in
        `\ZZ` needs to perform operations over `\QQ`.

        .. TODO::

            - The [BFSS2006]_ algorithm has been implemented here only in the case of
              polynomials over rationals. For other rings of zero characteristic
              (or if the characteristic is larger than the product of the degrees),
              one needs to implement a generic method ``_exp_series``. In the
              general case of nonzero characteristic there is an alternative
              algorithm in the same paper.

            - The Newton series computation can be done much more efficiently!
              See [BFSS2006]_.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: p1 = x^2 - 1
            sage: p2 = x^4 - 1
            sage: p1.composed_op(p2, operator.add)                                      # needs sage.libs.singular
            x^8 - 4*x^6 + 4*x^4 - 16*x^2
            sage: p1.composed_op(p2, operator.mul)                                      # needs sage.libs.singular
            x^8 - 2*x^4 + 1
            sage: p1.composed_op(p2, operator.truediv)                                  # needs sage.libs.singular
            x^8 - 2*x^4 + 1

        This function works over any field. However for base rings other than
        `\ZZ` and `\QQ` only the resultant algorithm is available::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQbar)
            sage: p1 = x**2 - AA(2).sqrt()
            sage: p2 = x**3 - AA(3).sqrt()
            sage: r1 = p1.roots(multiplicities=False)
            sage: r2 = p2.roots(multiplicities=False)
            sage: p = p1.composed_op(p2, operator.add); p
            x^6 - 4.242640687119285?*x^4 - 3.464101615137755?*x^3 + 6*x^2
             - 14.69693845669907?*x + 0.1715728752538099?
            sage: all(p(x+y).is_zero() for x in r1 for y in r2)
            True

            sage: x = polygen(GF(2))
            sage: p1 = x**2 + x - 1
            sage: p2 = x**3 + x - 1
            sage: p_add = p1.composed_op(p2, operator.add); p_add                       # needs sage.libs.singular
            x^6 + x^5 + x^3 + x^2 + 1
            sage: p_mul = p1.composed_op(p2, operator.mul); p_mul                       # needs sage.libs.singular
            x^6 + x^4 + x^2 + x + 1
            sage: p_div = p1.composed_op(p2, operator.truediv); p_div                   # needs sage.libs.singular
            x^6 + x^5 + x^4 + x^2 + 1

            sage: # needs sage.rings.finite_rings
            sage: K = GF(2**6, 'a')
            sage: r1 = p1.roots(K, multiplicities=False)
            sage: r2 = p2.roots(K, multiplicities=False)
            sage: all(p_add(x1+x2).is_zero() for x1 in r1 for x2 in r2)                 # needs sage.libs.singular
            True
            sage: all(p_mul(x1*x2).is_zero() for x1 in r1 for x2 in r2)                 # needs sage.libs.singular
            True
            sage: all(p_div(x1/x2).is_zero() for x1 in r1 for x2 in r2)                 # needs sage.libs.singular
            True

        TESTS:

        ::

            sage: y = polygen(ZZ)
            sage: for p1 in [2*y^3 - y + 3, -y^5 - 2, 4*y - 3]:                         # needs sage.libs.singular
            ....:   for p2 in [5*y^2 - 7, -3*y - 1]:
            ....:     for monic in [True,False]:
            ....:       for op in [operator.add, operator.sub, operator.mul, operator.truediv]:
            ....:         pr = p1.composed_op(p2, op, "resultant", monic=monic)
            ....:         pb = p1.composed_op(p2, op, "BFSS", monic=monic)
            ....:         assert ((pr == pb) or ((not monic) and pr == -pb) and (parent(pr) is parent(pb)))
        """
        cdef long j
        cdef long prec

        try:
            if op is operator.div:
                op = operator.truediv
        except AttributeError:
            pass

        if op not in (operator.add, operator.sub, operator.mul, operator.truediv):
            raise ValueError("op must be operator.OP where OP=add, sub, mul or truediv")

        if not isinstance(p2, Polynomial):
            raise TypeError("p2 must be a polynomial")
        p1, p2 = coercion_model.canonical_coercion(p1, p2)
        K = p1.parent()
        assert isinstance(p1.parent(), PolynomialRing_general)
        S = K.base_ring()
        Sf = S.fraction_field()

        cdef long d1 = p1.degree()
        cdef long d2 = p2.degree()
        if d1 <= 0 or d2 <= 0:
            raise ValueError('the polynomials must have positive degree')

        if op is operator.truediv and p2.valuation() > 0:
            raise ZeroDivisionError('p2 must have zero valuation')
        if algorithm is None:
            # choose the algorithm observing that the "resultant" one
            # is fast when there are few terms and the degrees are not high
            N = 7
            if Sf is not QQ or (d1 <= N and d2 <= N):
                algorithm = "resultant"
            else:
                c = d1*sum(bool(p1[i]) for i in range(d1 + 1))*\
                    d2*sum(bool(p2[i]) for i in range(d2 + 1))
                if c <= N**4:
                    algorithm = "resultant"
                else:
                    algorithm = "BFSS"

        if algorithm == "resultant":
            R = S['x', 'y']
            x = R.gen(0)
            y = R.gen(1)
            if op is operator.add:
                lp = p1(x - y)
            elif op is operator.sub:
                lp = p1(x + y)
            elif op is operator.mul:
                lp = p1(x).homogenize(y)
            else:
                lp = p1(x * y)
            q = p2(y).resultant(lp, y).univariate_polynomial(K)
            return q.monic() if monic else q

        elif algorithm == "BFSS":
            if Sf is not QQ:
                raise ValueError("BFSS algorithm is available only for the base ring ZZ or QQ")
            if op is operator.sub:
                p2 = p2(-K.gen())
            elif op is operator.truediv:
                p2 = p2.reverse()
            # the computation below needs must be done in the fraction field
            # even though the result would have the same ring
            if Sf is not S:
                K = K.change_ring(Sf)
                p1 = p1.change_ring(Sf)
                p2 = p2.change_ring(Sf)
            prec = d1*d2 + 1
            np1 = p1.reverse().inverse_series_trunc(prec)
            np1 = np1._mul_trunc_(p1.derivative().reverse(), prec)
            np2 = p2.reverse().inverse_series_trunc(prec)
            np2 = np2._mul_trunc_(p2.derivative().reverse(), prec)
            if op in (operator.add, operator.sub):
                # compute np1e and np2e, the Hadamard products of respectively
                # np1 and np2 with the exponential series. That is
                #  a0 + a1 x + a2 x^2 + ...
                #  ->
                #  a0 + a1/1! x + a2/2! x^2 + ...
                fj = Sf.one()
                a1, a2 = [np1[0]], [np2[0]]
                for j in range(1, prec):
                    fj = fj*j
                    a1.append(np1[j] / fj)
                    a2.append(np2[j] / fj)
                np1e = K(a1)
                np2e = K(a2)

                # recover the polynomial from its Newton series
                np3e = np1e*np2e
                fj = -Sf.one()
                a3 = [Sf.zero()]
                for j in range(1, prec):
                    a3.append(np3e[j] * fj)
                    fj = fj*j
                np = K(a3)
                q = np
            else:
                np = K([-np1[j]*np2[j] for j in range(1, prec)])
                q = np.integral()

            q = q._exp_series(prec).reverse()
            q = q.shift(prec - q.degree() - 1)
            if monic:
                return q
            else:
                return (p1.leading_coefficient()**p2.degree() *
                        p2.leading_coefficient()**p1.degree() * q).change_ring(S)

        else:
            raise ValueError('algorithm must be "resultant" or "BFSS"')

    def compose_power(self, k, algorithm=None, monic=False):
        r"""
        Return the `k`-th iterate of the composed product of this
        polynomial with itself.

        INPUT:

        - ``k`` -- nonnegative integer

        - ``algorithm`` -- ``None`` (default), ``'resultant'``, or ``'BFSS'``;
          see :meth:`.composed_op`

        - ``monic`` -- boolean (default: ``False``); see :meth:`.composed_op`

        OUTPUT:

        The polynomial of degree `d^k` where `d` is the degree, whose
        roots are all `k`-fold products of roots of this polynomial.
        That is, `f*f*\dots*f` where this is `f` and
        `f*f=` ``f.composed_op(f, operator.mul)``.

        EXAMPLES::

            sage: R.<a,b,c> = ZZ[]
            sage: x = polygen(R)
            sage: f = (x - a) * (x - b) * (x - c)
            sage: f.compose_power(2).factor()                                           # needs sage.libs.singular sage.modules
            (x - c^2) * (x - b^2) * (x - a^2) * (x - b*c)^2 * (x - a*c)^2 * (x - a*b)^2

            sage: # needs sage.libs.singular sage.modules
            sage: x = polygen(QQ)
            sage: f = x^2 - 2*x + 2
            sage: f2 = f.compose_power(2); f2
            x^4 - 4*x^3 + 8*x^2 - 16*x + 16
            sage: f2 == f.composed_op(f, operator.mul)
            True
            sage: f3 = f.compose_power(3); f3
            x^8 - 8*x^7 + 32*x^6 - 64*x^5 + 128*x^4 - 512*x^3 + 2048*x^2 - 4096*x + 4096
            sage: f3 == f2.composed_op(f, operator.mul)
            True
            sage: f4 = f.compose_power(4)
            sage: f4 == f3.composed_op(f, operator.mul)
            True
        """
        try:
            k = ZZ(k)
        except ValueError("Cannot iterate {} times".format(k)):
            return self
        if k < 0:
            raise ValueError("Cannot iterate a negative number {} of times".format(k))
        if k == 0:
            return self.variables()[0] - 1
        if k == 1:
            return self
        if k == 2:
            return self.composed_op(self, operator.mul,
                                    algorithm=algorithm, monic=monic)
        k2, k1 = k.quo_rem(2)
        # recurse to get the k/2 -iterate where k=2*k2+k1:
        R = self.compose_power(k2, algorithm=algorithm, monic=monic)
        # square:
        R = R.composed_op(R, operator.mul, algorithm=algorithm, monic=monic)
        # one more factor if k odd:
        if k1:
            R = R.composed_op(self, operator.mul)
        return R

    def adams_operator_on_roots(self, n, monic=False):
        r"""
        Return the polynomial whose roots are the `n`-th powers
        of the roots of ``self``.

        INPUT:

        - ``n`` -- integer

        - ``monic`` -- boolean (default: ``False``)
          if set to ``True``, force the output to be monic

        EXAMPLES::

            sage: # needs sage.libs.pari sage.libs.singular
            sage: f = cyclotomic_polynomial(30)
            sage: f.adams_operator_on_roots(7) == f
            True
            sage: f.adams_operator_on_roots(6) == cyclotomic_polynomial(5)**2
            True
            sage: f.adams_operator_on_roots(10) == cyclotomic_polynomial(3)**4
            True
            sage: f.adams_operator_on_roots(15) == cyclotomic_polynomial(2)**8
            True
            sage: f.adams_operator_on_roots(30) == cyclotomic_polynomial(1)**8
            True

            sage: x = polygen(QQ)
            sage: f = x^2 - 2*x + 2
            sage: f.adams_operator_on_roots(10)                                                  # needs sage.libs.singular
            x^2 + 1024

        When ``self`` is monic, the output will have leading coefficient
        `\pm1` depending on the degree, but we can force it to be
        monic::

            sage: R.<a,b,c> = ZZ[]
            sage: x = polygen(R)
            sage: f = (x - a) * (x - b) * (x - c)
            sage: f.adams_operator_on_roots(3).factor()                                          # needs sage.libs.singular
            (-1) * (x - c^3) * (x - b^3) * (x - a^3)
            sage: f.adams_operator_on_roots(3, monic=True).factor()                              # needs sage.libs.singular
            (x - c^3) * (x - b^3) * (x - a^3)
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        u, v = PolynomialRing(self._parent.base_ring(), ['u', 'v']).gens()
        R = (u - v**n).resultant(self(v), v)
        R = R([self.variables()[0], 0])
        if monic:
            R = R.monic()
        return R

    adams_operator = deprecated_function_alias(36396, adams_operator_on_roots)

    def symmetric_power(self, k, monic=False):
        r"""
        Return the polynomial whose roots are products of `k`-th distinct
        roots of this.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^4 - x + 2
            sage: [f.symmetric_power(k) for k in range(5)]                              # needs sage.libs.singular
            [x - 1, x^4 - x + 2, x^6 - 2*x^4 - x^3 - 4*x^2 + 8, x^4 - x^3 + 8, x - 2]

            sage: f = x^5 - 2*x + 2
            sage: [f.symmetric_power(k) for k in range(6)]                              # needs sage.libs.singular
            [x - 1,
             x^5 - 2*x + 2,
             x^10 + 2*x^8 - 4*x^6 - 8*x^5 - 8*x^4 - 8*x^3 + 16,
             x^10 + 4*x^7 - 8*x^6 + 16*x^5 - 16*x^4 + 32*x^2 + 64,
             x^5 + 2*x^4 - 16,
             x + 2]

            sage: R.<a,b,c,d> = ZZ[]
            sage: x = polygen(R)
            sage: f = (x - a) * (x - b) * (x - c) * (x - d)
            sage: [f.symmetric_power(k).factor() for k in range(5)]                     # needs sage.libs.singular
            [x - 1,
             (-x + d) * (-x + c) * (-x + b) * (-x + a),
             (x - c*d) * (x - b*d) * (x - a*d) * (x - b*c) * (x - a*c) * (x - a*b),
             (x - b*c*d) * (x - a*c*d) * (x - a*b*d) * (x - a*b*c),
             x - a*b*c*d]
        """
        try:
            k = ZZ(k)
        except (ValueError, TypeError):
            raise ValueError(f"Cannot compute k'th symmetric power for k={k}")
        n = self.degree()
        if k < 0 or k > n:
            raise ValueError(f"Cannot compute k'th symmetric power for k={k}")
        x = self.variables()[0]
        if k == 0:
            return x - 1
        if k == 1:
            if monic:
                return self.monic()
            return self
        c = (-1)**n * self(0)
        if k == n:
            return x - c
        if k > n - k:  # use (n-k)'th symmetric power
            g = self.symmetric_power(n - k, monic=monic)
            from sage.arith.misc import binomial
            g = ((-x)**binomial(n,k) * g(c/x) / c**binomial(n-1,k)).numerator()
            if monic:
                g = g.monic()
            return g

        def star(g, h):
            return g.composed_op(h, operator.mul, monic=True)

        def rpow(g, n):
            return g.adams_operator_on_roots(n, monic=True)
        if k == 2:
            g = (star(self, self) // rpow(self, 2)).nth_root(2)
            if monic:
                g = g.monic()
            return g
        if k == 3:
            g = star(self.symmetric_power(2, monic=monic), self) * rpow(self, 3)
            h = star(rpow(self, 2), self)
            g = (g // h).nth_root(3)
            if monic:
                g = g.monic()
            return g

        fkn = fkd = self._parent.one()
        for j in range(1, k + 1):
            g = star(rpow(self, j), self.symmetric_power(k - j))
            if j % 2:
                fkn *= g
            else:
                fkd *= g

        fk = fkn // fkd
        assert fk * fkd == fkn
        g = fk.nth_root(k)
        if monic:
            g = g.monic()
        return g

    def discriminant(self):
        r"""
        Return the discriminant of ``self``.

        The discriminant is

        .. MATH::

            R_n := a_n^{2 n-2} \prod_{1<i<j<n} (r_i-r_j)^2,

        where `n` is the degree of ``self``, `a_n` is the
        leading coefficient of ``self``, and the roots of ``self`` are
        `r_1, \ldots, r_n`.

        OUTPUT: an element of the base ring of the polynomial ring

        ALGORITHM:

        Uses the identity `R_n(f) := (-1)^{n (n-1)/2} R(f, f')
        a_n^{n-k-2}`, where `n` is the degree of ``self``, `a_n` is the
        leading coefficient of ``self``, `f'` is the derivative of `f`,
        and `k` is the degree of `f'`. Calls :meth:`.resultant`.

        EXAMPLES:

        In the case of elliptic curves in special form, the discriminant is
        easy to calculate::

            sage: R.<x> = QQ[]
            sage: f = x^3 + x + 1
            sage: d = f.discriminant(); d                                               # needs sage.libs.pari
            -31
            sage: d.parent() is QQ                                                      # needs sage.libs.pari
            True
            sage: EllipticCurve([1, 1]).discriminant()/16                               # needs sage.libs.pari sage.schemes
            -31

        ::

            sage: R.<x> = QQ[]
            sage: f = 2*x^3 + x + 1
            sage: d = f.discriminant(); d                                               # needs sage.libs.pari
            -116

        We can compute discriminants over univariate and multivariate
        polynomial rings::

            sage: R.<a> = QQ[]
            sage: S.<x> = R[]
            sage: f = a*x + x + a + 1
            sage: d = f.discriminant(); d                                               # needs sage.libs.pari
            1
            sage: d.parent() is R                                                       # needs sage.libs.pari
            True

        ::

            sage: R.<a, b> = QQ[]
            sage: S.<x> = R[]
            sage: f = x^2 + a + b
            sage: d = f.discriminant(); d                                               # needs sage.libs.pari
            -4*a - 4*b
            sage: d.parent() is R                                                       # needs sage.libs.pari
            True

        TESTS::

            sage: R.<x, y> = QQ[]
            sage: S.<a> = R[]
            sage: f = x^2 + a
            sage: f.discriminant()                                                      # needs sage.libs.pari
            1

        Check that :issue:`13672` is fixed::

            sage: R.<t> = GF(5)[]
            sage: S.<x> = R[]
            sage: f = x^10 + 2*x^6 + 2*x^5 + x + 2
            sage: (f - t).discriminant()                                                # needs sage.rings.finite_rings
            4*t^5

        The following examples show that :issue:`11782` has been fixed::

            sage: # needs sage.libs.pari sage.symbolic
            sage: var('x')
            x
            sage: ZZ.quo(81)['x'](3*x^2 + 3*x + 3).discriminant()
            54
            sage: ZZ.quo(9)['x'](2*x^3 + x^2 + x).discriminant()
            2

        This was fixed by :issue:`15422`::

            sage: R.<s> = PolynomialRing(Qp(2))                                         # needs sage.rings.padics
            sage: (s^2).discriminant()                                                  # needs sage.rings.padics
            0

        This was fixed by :issue:`16014`::

            sage: # needs sage.modules
            sage: PR.<b,t1,t2,x1,y1,x2,y2> = QQ[]
            sage: PRmu.<mu> = PR[]
            sage: E1 = diagonal_matrix(PR, [1, b^2, -b^2])
            sage: M = matrix(PR, [[1,-t1,x1-t1*y1], [t1,1,y1+t1*x1], [0,0,1]])
            sage: E1 = M.transpose()*E1*M
            sage: E2 = E1.subs(t1=t2, x1=x2, y1=y2)
            sage: det(mu*E1 + E2).discriminant().degrees()                              # needs sage.libs.pari
            (24, 12, 12, 8, 8, 8, 8)

        This addresses an issue raised by :issue:`15061`::

            sage: R.<T> = PowerSeriesRing(QQ)
            sage: F = R([1,1],2)
            sage: RP.<x> = PolynomialRing(R)
            sage: P = x^2 - F
            sage: P.discriminant()                                                      # needs sage.libs.pari
            4 + 4*T + O(T^2)
        """
        # Late import to avoid cyclic dependencies:
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        try:
            from sage.rings.lazy_series_ring import LazyPowerSeriesRing
        except ImportError:
            LazyPowerSeriesRing = ()
        if self.is_zero():
            return self  # return 0
        n = self.degree()
        base_ring = self._parent.base_ring()
        if isinstance(base_ring, (MPolynomialRing_base, PowerSeriesRing_generic, LazyPowerSeriesRing)):
            # It is often cheaper to compute discriminant of simple
            # multivariate polynomial and substitute the real
            # coefficients into that result (see #16014).
            return universal_discriminant(n)(list(self))
        d = self.derivative()
        k = d.degree()

        r = n % 4
        u = -1 # (-1)**(n*(n-1)/2)
        if r == 0 or r == 1:
            u = 1
        try:
            an = self.get_coeff_c(n)**(n - k - 2)
        except ZeroDivisionError:
            assert(n-k-2 == -1)
            # Rather than dividing the resultant by the leading coefficient,
            # we alter the Sylvester matrix (see #11782).
            mat = self.sylvester_matrix(d)
            mat[0, 0] = base_ring.one()
            mat[n - 1, 0] = base_ring(n)
            return u * mat.determinant()
        else:
            return base_ring(u * self.resultant(d) * an)

    def reverse(self, degree=None):
        """
        Return polynomial but with the coefficients reversed.

        If an optional ``degree`` argument is given, the coefficient list will be
        truncated or zero padded as necessary before reversing it. Assuming
        that the constant coefficient of ``self`` is nonzero, the reverse
        polynomial will have the specified degree.

        EXAMPLES::

            sage: R.<x> = ZZ[]; S.<y> = R[]
            sage: f = y^3 + x*y - 3*x; f
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
            ValueError: degree argument must be a nonnegative integer, got 1.5

            sage: f.reverse(0)
            -3*x
            sage: f
            y^3 + x*y - 3*x
        """
        v = self.list()

        cdef unsigned long d
        if degree is not None:
            if degree < 0:
                raise ValueError("degree argument must be a nonnegative integer, got %s" % (degree))
            d = degree
            if d != degree:
                raise ValueError("degree argument must be a nonnegative integer, got %s" % (degree))
            if len(v) < degree+1:
                v.reverse()
                v = [self.base_ring().zero()]*(degree+1-len(v)) + v
            elif len(v) > degree+1:
                v = v[:degree+1]
                v.reverse()
            else: # len(v) == degree + 1
                v.reverse()
        else:
            v.reverse()

        return self._new_generic(v)

    def roots(self, ring=None, multiplicities=True, algorithm=None, **kwds):
        r"""
        Return the roots of this polynomial (by default, in the base ring
        of this polynomial).

        INPUT:

        - ``ring`` -- the ring to find roots in

        - ``multiplicities`` -- boolean (default: ``True``); if ``True``
          return list of pairs `(r, n)`, where `r` is the root and `n` is the
          multiplicity. If ``False``, just return the unique roots, with no
          information about multiplicities.

        - ``algorithm`` -- the root-finding algorithm to use.
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

        If the output ring is a :class:`RealIntervalField` or :class:`ComplexIntervalField`
        of a given precision, then the answer will always be correct (or an
        exception will be raised, if a case is not implemented). Each root
        will be contained in one of the returned intervals, and the
        intervals will be disjoint. (The returned intervals may be of
        higher precision than the specified output ring.)

        At the end of this docstring (after the examples) is a description
        of all the cases implemented in this function, and the algorithms
        used. That section also describes the possibilities for the
        ``algorithm`` keyword, for the cases where multiple algorithms exist.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: x = QQ['x'].0
            sage: f = x^3 - 1
            sage: f.roots()
            [(1, 1)]
            sage: f.roots(ring=CC)   # ... - low order bits slightly different on ppc   # needs sage.rings.real_mpfr
            [(1.00000000000000, 1),
             (-0.500000000000000 - 0.86602540378443...*I, 1),
             (-0.500000000000000 + 0.86602540378443...*I, 1)]
            sage: f = (x^3 - 1)^2
            sage: f.roots()
            [(1, 2)]
            sage: f = -19*x + 884736
            sage: f.roots()
            [(884736/19, 1)]
            sage: (f^20).roots()
            [(884736/19, 20)]

        ::

            sage: # needs sage.rings.number_field
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
        :meth:`_roots_univariate_polynomial` to the base ring, and observe
        that this method is called instead to find roots of
        polynomials over this ring.  This facility can be used to
        easily extend root finding to work over new rings you
        introduce::

             sage: R.<x> = QQ[]
             sage: (x^2 + 1).roots()                                                    # needs sage.libs.pari
             []
             sage: def my_roots(f, *args, **kwds):
             ....:     return f.change_ring(CDF).roots()
             sage: QQ._roots_univariate_polynomial = my_roots
             sage: (x^2 + 1).roots()  # abs tol 1e-14                                   # needs numpy
             [(2.7755575615628914e-17 - 1.0*I, 1), (0.9999999999999997*I, 1)]
             sage: del QQ._roots_univariate_polynomial

        An example over RR, which illustrates that only the roots in RR are
        returned::

            sage: # needs numpy sage.rings.real_mpfr
            sage: x = RR['x'].0
            sage: f = x^3 - 2
            sage: f.roots()
            [(1.25992104989487, 1)]
            sage: f.factor()
            (x - 1.25992104989487) * (x^2 + 1.25992104989487*x + 1.58740105196820)
            sage: x = RealField(100)['x'].0
            sage: f = x^3 - 2
            sage: f.roots()
            [(1.2599210498948731647672106073, 1)]

        ::

            sage: # needs sage.rings.real_mpfr
            sage: x = CC['x'].0
            sage: f = x^3 - 2
            sage: f.roots()                                                             # needs numpy
            [(1.25992104989487, 1),
             (-0.62996052494743... - 1.09112363597172*I, 1),
             (-0.62996052494743... + 1.09112363597172*I, 1)]
            sage: f.roots(algorithm='pari')                                             # needs sage.libs.pari
            [(1.25992104989487, 1),
             (-0.629960524947437 - 1.09112363597172*I, 1),
             (-0.629960524947437 + 1.09112363597172*I, 1)]

        Another example showing that only roots in the base ring are
        returned::

            sage: x = polygen(ZZ)
            sage: f = (2*x - 3) * (x - 1) * (x + 1)
            sage: f.roots()                                                             # needs sage.libs.pari
            [(1, 1), (-1, 1)]
            sage: f.roots(ring=QQ)                                                      # needs sage.libs.pari
            [(3/2, 1), (1, 1), (-1, 1)]

        An example where we compute the roots lying in a subring of the
        base ring::

            sage: Pols.<n> = QQ[]
            sage: pol = (n - 1/2)^2 * (n - 1)^2 * (n - 2)
            sage: pol.roots(ZZ)                                                         # needs sage.libs.pari
            [(2, 1), (1, 2)]

        An example involving large numbers::

            sage: # needs numpy sage.rings.real_mpfr
            sage: x = RR['x'].0
            sage: f = x^2 - 1e100
            sage: f.roots()
            [(-1.00000000000000e50, 1), (1.00000000000000e50, 1)]
            sage: f = x^10 - 2 * (5*x - 1)^2
            sage: f.roots(multiplicities=False)
            [-1.6772670339941..., 0.19995479628..., 0.20004530611..., 1.5763035161844...]

        ::

            sage: # needs numpy sage.rings.real_mpfr
            sage: x = CC['x'].0
            sage: i = CC.0
            sage: f = (x - 1) * (x - i)
            sage: f.roots(multiplicities=False)
            [1.00000000000000, 1.00000000000000*I]
            sage: g = (x - 1.33 + 1.33*i) * (x - 2.66 - 2.66*i)
            sage: g.roots(multiplicities=False)
            [1.33000000000000 - 1.33000000000000*I, 2.66000000000000 + 2.66000000000000*I]

        Describing roots using radical expressions::

            sage: x = QQ['x'].0
            sage: f = x^2 + 2
            sage: f.roots(SR)                                                           # needs sage.symbolic
            [(-I*sqrt(2), 1), (I*sqrt(2), 1)]
            sage: f.roots(SR, multiplicities=False)                                     # needs sage.symbolic
            [-I*sqrt(2), I*sqrt(2)]

        The roots of some polynomials cannot be described using radical
        expressions::

            sage: (x^5 - x + 1).roots(SR)                                               # needs sage.symbolic
            []

        For some other polynomials, no roots can be found at the moment
        due to the way roots are computed. :issue:`17516` addresses
        these defects. Until that gets implemented, one such example
        is the following::

            sage: f = x^6 - 300*x^5 + 30361*x^4 - 1061610*x^3 + 1141893*x^2 - 915320*x + 101724
            sage: f.roots()                                                             # needs sage.libs.pari
            []

        A purely symbolic roots example::

            sage: # needs sage.symbolic
            sage: X = var('X')
            sage: f = expand((X - 1) * (X - I)^3 * (X^2 - sqrt(2))); f
            X^6 - (3*I + 1)*X^5 - sqrt(2)*X^4 + (3*I - 3)*X^4 + (3*I + 1)*sqrt(2)*X^3
            + (I + 3)*X^3 - (3*I - 3)*sqrt(2)*X^2 - I*X^2 - (I + 3)*sqrt(2)*X + I*sqrt(2)
            sage: f.roots()
            [(I, 3), (-2^(1/4), 1), (2^(1/4), 1), (1, 1)]

        The same operation, performed over a polynomial ring
        with symbolic coefficients::

            sage: # needs sage.symbolic
            sage: X = SR['X'].0
            sage: f = (X - 1) * (X - I)^3 * (X^2 - sqrt(2)); f
            X^6 + (-3*I - 1)*X^5 + (-sqrt(2) + 3*I - 3)*X^4 + ((3*I + 1)*sqrt(2) + I + 3)*X^3
            + (-(3*I - 3)*sqrt(2) - I)*X^2 + (-(I + 3)*sqrt(2))*X + I*sqrt(2)
            sage: f.roots()
            [(I, 3), (-2^(1/4), 1), (2^(1/4), 1), (1, 1)]
            sage: f.roots(multiplicities=False)
            [I, -2^(1/4), 2^(1/4), 1]

        A couple of examples where the base ring does not have a
        factorization algorithm (yet). Note that this is currently done via
        a rather naive enumeration, so could be very slow::

            sage: R = Integers(6)
            sage: S.<x> = R['x']
            sage: p = x^2 - 1
            sage: p.roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial
            not implemented (try the multiplicities=False option)
            sage: p.roots(multiplicities=False)                                         # needs sage.libs.pari
            [5, 1]
            sage: R = Integers(9)
            sage: A = PolynomialRing(R, 'y')
            sage: y = A.gen()
            sage: f = 10*y^2 - y^3 - 9
            sage: f.roots(multiplicities=False)                                         # needs sage.libs.pari
            [1, 0, 3, 6]

        An example over the complex double field (where root finding is
        fast, thanks to NumPy)::

            sage: # needs numpy sage.rings.complex_double
            sage: R.<x> = CDF[]
            sage: f = R.cyclotomic_polynomial(5); f
            x^4 + x^3 + x^2 + x + 1.0
            sage: f.roots(multiplicities=False)  # abs tol 1e-9
            [-0.8090169943749469 - 0.5877852522924724*I, -0.8090169943749473 + 0.5877852522924724*I,
             0.30901699437494773 - 0.951056516295154*I, 0.30901699437494756 + 0.9510565162951525*I]
            sage: [z^5 for z in f.roots(multiplicities=False)]  # abs tol 2e-14
            [0.9999999999999957 - 1.2864981197413038e-15*I, 0.9999999999999976 + 3.062854959141552e-15*I,
             1.0000000000000024 + 1.1331077795295987e-15*I, 0.9999999999999953 - 2.0212861992297117e-15*I]
            sage: f = CDF['x']([1,2,3,4]); f
            4.0*x^3 + 3.0*x^2 + 2.0*x + 1.0
            sage: r = f.roots(multiplicities=False)
            sage: [f(a).abs() for a in r]  # abs tol 1e-14
            [2.574630599127759e-15, 1.457101633618084e-15, 1.1443916996305594e-15]

        Another example over RDF::

            sage: x = RDF['x'].0
            sage: ((x^3 - 1)).roots()  # abs tol 4e-16                                  # needs numpy
            [(1.0000000000000002, 1)]
            sage: ((x^3 - 1)).roots(multiplicities=False)  # abs tol 4e-16              # needs numpy
            [1.0000000000000002]

        More examples involving the complex double field::

            sage: # needs numpy sage.rings.complex_double sage.rings.real_mpfr
            sage: x = CDF['x'].0
            sage: i = CDF.0
            sage: f = x^3 + 2*i; f
            x^3 + 2.0*I
            sage: f.roots()
            [(-1.09112363597172... - 0.62996052494743...*I, 1),
             (...1.25992104989487...*I, 1),
             (1.09112363597172... - 0.62996052494743...*I, 1)]
            sage: f.roots(multiplicities=False)
            [-1.09112363597172... - 0.62996052494743...*I, ...1.25992104989487...*I,
              1.09112363597172... - 0.62996052494743...*I]
            sage: [abs(f(z)) for z in f.roots(multiplicities=False)]  # abs tol 1e-14
            [8.95090418262362e-16, 8.728374398092689e-16, 1.0235750533041806e-15]
            sage: f = i*x^3 + 2; f
            I*x^3 + 2.0
            sage: f.roots()
            [(-1.09112363597172... + 0.62996052494743...*I, 1),
             (...1.25992104989487...*I, 1),
             (1.09112363597172... + 0.62996052494743...*I, 1)]
            sage: abs(f(f.roots()[0][0]))  # abs tol 1e-13
            1.1102230246251565e-16

        Examples using real root isolation::

            sage: x = polygen(ZZ)
            sage: f = x^2 - x - 1
            sage: f.roots()                                                             # needs sage.libs.pari
            []
            sage: f.roots(ring=AA)                                                      # needs sage.rings.number_field
            [(-0.618033988749895?, 1), (1.618033988749895?, 1)]

            sage: # needs sage.rings.real_interval_field
            sage: f.roots(ring=RIF)
            [(-0.6180339887498948482045868343657?, 1), (1.6180339887498948482045868343657?, 1)]
            sage: f.roots(ring=RIF, multiplicities=False)
            [-0.6180339887498948482045868343657?, 1.6180339887498948482045868343657?]
            sage: f.roots(ring=RealIntervalField(150))
            [(-0.6180339887498948482045868343656381177203091798057628621354486227?, 1),
             (1.618033988749894848204586834365638117720309179805762862135448623?, 1)]
            sage: f = f^2 * (x - 1)
            sage: f.roots(ring=RIF)
            [(-0.6180339887498948482045868343657?, 2),
             (1.0000000000000000000000000000000?, 1),
             (1.6180339887498948482045868343657?, 2)]
            sage: f.roots(ring=RIF, multiplicities=False)
            [-0.6180339887498948482045868343657?,
             1.0000000000000000000000000000000?,
             1.6180339887498948482045868343657?]

        Examples using complex root isolation::

            sage: x = polygen(ZZ)
            sage: p = x^5 - x - 1
            sage: p.roots()                                                             # needs sage.libs.pari
            []
            sage: p.roots(ring=CIF)                                                     # needs sage.rings.complex_interval_field
            [(1.167303978261419?, 1),
             (-0.764884433600585? - 0.352471546031727?*I, 1),
             (-0.764884433600585? + 0.352471546031727?*I, 1),
             (0.181232444469876? - 1.083954101317711?*I, 1),
             (0.181232444469876? + 1.083954101317711?*I, 1)]
            sage: p.roots(ring=ComplexIntervalField(200))                               # needs sage.rings.complex_interval_field
            [(1.167303978261418684256045899854842180720560371525489039140082?, 1),
             (-0.76488443360058472602982318770854173032899665194736756700778? - 0.35247154603172624931794709140258105439420648082424733283770?*I, 1),
             (-0.76488443360058472602982318770854173032899665194736756700778? + 0.35247154603172624931794709140258105439420648082424733283770?*I, 1),
             (0.18123244446987538390180023778112063996871646618462304743774? - 1.08395410131771066843034449298076657427364024315511565430114?*I, 1),
             (0.18123244446987538390180023778112063996871646618462304743774? + 1.08395410131771066843034449298076657427364024315511565430114?*I, 1)]
            sage: rts = p.roots(ring=QQbar); rts                                        # needs sage.rings.number_field
            [(1.167303978261419?, 1),
             (-0.7648844336005847? - 0.3524715460317263?*I, 1),
             (-0.7648844336005847? + 0.3524715460317263?*I, 1),
             (0.1812324444698754? - 1.083954101317711?*I, 1),
             (0.1812324444698754? + 1.083954101317711?*I, 1)]
            sage: p.roots(ring=AA)                                                      # needs sage.rings.number_field
            [(1.167303978261419?, 1)]
            sage: p = (x - rts[4][0])^2 * (3*x^2 + x + 1)                               # needs sage.rings.number_field
            sage: p.roots(ring=QQbar)                                                   # needs sage.rings.number_field
            [(-0.1666666666666667? - 0.552770798392567?*I, 1),
             (-0.1666666666666667? + 0.552770798392567?*I, 1),
             (0.1812324444698754? + 1.083954101317711?*I, 2)]
            sage: p.roots(ring=CIF)                                                     # needs sage.rings.complex_interval_field
            [(-0.1666666666666667? - 0.552770798392567?*I, 1),
             (-0.1666666666666667? + 0.552770798392567?*I, 1),
             (0.1812324444698754? + 1.083954101317711?*I, 2)]

        In some cases, it is possible to isolate the roots of polynomials over
        complex ball fields::

            sage: # needs sage.libs.flint
            sage: Pol.<x> = CBF[]
            sage: (x^2 + 2).roots(multiplicities=False)
            [[+/- ...e-19] + [-1.414213562373095 +/- ...e-17]*I,
            [+/- ...e-19] + [1.414213562373095 +/- ...e-17]*I]
            sage: (x^3 - 1/2).roots(RBF, multiplicities=False)
            [[0.7937005259840997 +/- ...e-17]]
            sage: ((x - 1)^2).roots(multiplicities=False, proof=False)
            doctest:...
            UserWarning: roots may have been lost...
            [[1.00000000000 +/- ...e-12] + [+/- ...e-11]*I,
             [1.0000000000 +/- ...e-12] + [+/- ...e-12]*I]
            sage: ((x - 1)^2).roots(multiplicities=False, proof=False, warn=False)
            [[1.00000000000 +/- ...e-12] + [+/- ...e-11]*I,
             [1.0000000000 +/- ...e-12] + [+/- ...e-12]*I]

        Note that coefficients in a number field with defining polynomial
        `x^2 + 1` are considered to be Gaussian rationals (with the
        generator mapping to `+I`), if you ask for complex roots.

        ::

            sage: # needs sage.rings.number_field
            sage: K.<im> = QuadraticField(-1)
            sage: y = polygen(K)
            sage: p = y^4 - 2 - im
            sage: p.roots(ring=CC)
            [(-1.2146389322441... - 0.14142505258239...*I, 1),
             (-0.14142505258239... + 1.2146389322441...*I, 1),
             (0.14142505258239... - 1.2146389322441...*I, 1),
             (1.2146389322441... + 0.14142505258239...*I, 1)]
            sage: p = p^2 * (y^2 - 2)
            sage: p.roots(ring=CIF)
            [(-1.414213562373095?, 1), (1.414213562373095?, 1),
             (-1.214638932244183? - 0.141425052582394?*I, 2),
             (-0.141425052582394? + 1.214638932244183?*I, 2),
             (0.141425052582394? - 1.214638932244183?*I, 2),
             (1.214638932244183? + 0.141425052582394?*I, 2)]

        Note that one should not use NumPy when wanting high precision
        output as it does not support any of the high precision types::

            sage: # needs numpy sage.rings.real_mpfr sage.symbolic
            sage: R.<x> = RealField(200)[]
            sage: f = x^2 - R(pi)
            sage: f.roots()
            [(-1.7724538509055160272981674833411451827975494561223871282138, 1),
             (1.7724538509055160272981674833411451827975494561223871282138, 1)]
            sage: f.roots(algorithm='numpy')
            doctest... UserWarning: NumPy does not support arbitrary precision arithmetic.
            The roots found will likely have less precision than you expect.
            [(-1.77245385090551..., 1), (1.77245385090551..., 1)]

        We can also find roots over number fields::

            sage: K.<z> = CyclotomicField(15)                                           # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(K)                                             # needs sage.rings.number_field
            sage: (x^2 + x + 1).roots()                                                 # needs sage.rings.number_field
            [(z^5, 1), (-z^5 - 1, 1)]

        There are many combinations of floating-point input and output
        types that work. (Note that some of them are quite pointless like using
        ``algorithm='numpy'`` with high-precision types.)

        ::

            sage: # needs numpy sage.rings.complex_double sage.rings.real_mpfr
            sage: rflds = (RR, RDF, RealField(100))
            sage: cflds = (CC, CDF, ComplexField(100))
            sage: def cross(a, b):
            ....:     return list(cartesian_product_iterator([a, b]))
            sage: flds = cross(rflds, rflds) + cross(rflds, cflds) + cross(cflds, cflds)
            sage: for (fld_in, fld_out) in flds:
            ....:     x = polygen(fld_in)
            ....:     f = x^3 - fld_in(2)
            ....:     x2 = polygen(fld_out)
            ....:     f2 = x2^3 - fld_out(2)
            ....:     for algo in (None, 'pari', 'numpy'):
            ....:         rts = f.roots(ring=fld_out, multiplicities=False)
            ....:         rts = sorted(rts, key=lambda x: x.imag())
            ....:         if fld_in == fld_out and algo is None:
            ....:             print("{} {}".format(fld_in, rts))
            ....:         for rt in rts:
            ....:             assert(abs(f2(rt)) <= 1e-10)
            ....:             assert(rt.parent() == fld_out)
            Real Field with 53 bits of precision [1.25992104989487]
            Real Double Field [1.25992104989...]
            Real Field with 100 bits of precision [1.2599210498948731647672106073]
            Complex Field with 53 bits of precision [-0.62996052494743... - 1.09112363597172*I, 1.25992104989487, -0.62996052494743... + 1.09112363597172*I]
            Complex Double Field [-0.629960524947... - 1.0911236359717...*I, 1.25992104989487..., -0.629960524947... + 1.0911236359717...*I]
            Complex Field with 100 bits of precision [-0.62996052494743658238360530364 - 1.0911236359717214035600726142*I, 1.2599210498948731647672106073, -0.62996052494743658238360530364 + 1.0911236359717214035600726142*I]

        Note that we can find the roots of a polynomial with algebraic
        coefficients::

            sage: # needs sage.rings.number_field
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
            [(1.414213562373095048801688724209698078569671875376948073176679738?, 1),
             (1.732050807568877293527446341505872366942805253810380628055806980?, 1)]

        We can handle polynomials with huge coefficients.

        This number doesn't even fit in an IEEE double-precision float, but
        ``RR`` and ``CC`` allow a much larger range of floating-point numbers::

            sage: bigc = 2^1500
            sage: CDF(bigc)                                                             # needs sage.rings.complex_double
            +infinity
            sage: CC(bigc)                                                              # needs sage.rings.real_mpfr
            3.50746621104340e451

        Polynomials using such large coefficients can't be handled by
        numpy, but pari can deal with them::

            sage: x = polygen(QQ)
            sage: p = x + bigc
            sage: p.roots(ring=RR, algorithm='numpy')                                   # needs numpy sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            LinAlgError: Array must not contain infs or NaNs
            sage: p.roots(ring=RR, algorithm='pari')                                    # needs sage.libs.pari sage.rings.real_mpfr
            [(-3.50746621104340e451, 1)]
            sage: p.roots(ring=AA)                                                      # needs sage.rings.number_field
            [(-3.5074662110434039?e451, 1)]
            sage: p.roots(ring=QQbar)                                                   # needs sage.rings.number_field
            [(-3.5074662110434039?e451, 1)]
            sage: p = bigc*x + 1
            sage: p.roots(ring=RR)                                                      # needs numpy
            [(-2.85106096489671e-452, 1)]
            sage: p.roots(ring=AA)                                                      # needs sage.rings.number_field
            [(-2.8510609648967059?e-452, 1)]
            sage: p.roots(ring=QQbar)                                                   # needs sage.rings.number_field
            [(-2.8510609648967059?e-452, 1)]
            sage: p = x^2 - bigc
            sage: p.roots(ring=RR)                                                      # needs numpy
            [(-5.92238652153286e225, 1), (5.92238652153286e225, 1)]
            sage: p.roots(ring=QQbar)                                                   # needs sage.rings.number_field
            [(-5.9223865215328558?e225, 1), (5.9223865215328558?e225, 1)]

        Check that :issue:`30522` is fixed::

            sage: PolynomialRing(SR, names="x")("x^2").roots()                          # needs sage.symbolic
            [(0, 2)]

        Check that :issue:`30523` is fixed::

            sage: PolynomialRing(SR, names="x")("x^2 + q").roots()                      # needs sage.symbolic
            [(-sqrt(-q), 1), (sqrt(-q), 1)]

        ALGORITHM:

        For brevity, we will use ``RR`` to mean any :class:`RealField` of any precision;
        similarly for ``RIF``, ``CC``, and ``CIF``. Since Sage has no specific
        implementation of Gaussian rationals (or of number fields with
        embedding, at all), when we refer to Gaussian rationals below we
        will accept any number field with defining polynomial
        `x^2+1`, mapping the field generator to +I.

        We call the base ring of the polynomial `K`, and the ring given by
        the ``ring`` argument `L`. (If ``ring`` is not specified, then `L` is the
        same as `K`.)

        If `K` and `L` are floating-point (``RDF``, ``CDF``, ``RR``, or ``CC``), then a
        floating-point root-finder is used. If `L` is ``RDF`` or ``CDF``, then we
        default to using NumPy's :func:`roots`; otherwise, we use PARI's
        function :pari:`polroots`. This choice can be overridden with
        ``algorithm='pari'`` or ``algorithm='numpy'``. If the algorithm is
        unspecified and NumPy's :func:`roots` algorithm fails, then we fall
        back to PARI (NumPy will fail if some coefficient is infinite,
        for instance).

        If `L` is ``SR`` (or one of its subrings), then the roots will be radical
        expressions, computed as the solutions of a symbolic polynomial
        expression. At the moment this delegates to
        :meth:`sage.symbolic.expression.Expression.solve`
        which in turn uses Maxima to find radical solutions.
        Some solutions may be lost in this approach.
        Once :issue:`17516` gets implemented, all possible radical
        solutions should become available.

        If `L` is ``AA`` or ``RIF``, and `K` is ``ZZ``, ``QQ``, or ``AA``, then the root isolation
        algorithm :func:`sage.rings.polynomial.real_roots.real_roots` is used.
        (You can call :func:`real_roots` directly to get more control than this
        method gives.)

        If `L` is ``QQbar`` or ``CIF``, and `K` is ``ZZ``, ``QQ``, ``AA``, ``QQbar``, or the Gaussian
        rationals, then the root isolation algorithm
        :func:`sage.rings.polynomial.complex_roots.complex_roots` is used. (You
        can call :func:`complex_roots` directly to get more control than this
        method gives.)

        If `L` is ``AA`` and `K` is ``QQbar`` or the Gaussian rationals, then
        :func:`complex_roots` is used (as above) to find roots in ``QQbar``, then
        these roots are filtered to select only the real roots.

        If `L` is floating-point and `K` is not, then we attempt to change the
        polynomial ring to `L` (using :meth:`change_ring`) (or, if `L` is complex
        and `K` is not, to the corresponding real field). Then we use either
        PARI or NumPy as specified above.

        For all other cases where `K` is different from `L`, we attempt to use
        ``.change_ring(L)``. When that fails but `L` is a subring of `K`, we also
        attempt to compute the roots over `K` and filter the ones belonging
        to `L`.

        The next method, which is used if `K` is an integral domain, is to
        attempt to factor the polynomial. If this succeeds, then for every
        degree-one factor `ax+b`, we add `-b/a` as a root (as long as this
        quotient is actually in the desired ring).

        If factoring over `K` is not implemented (or `K` is not an integral
        domain), and `K` is finite, then we find the roots by enumerating all
        elements of `K` and checking whether the polynomial evaluates to zero
        at that value.

        .. NOTE::

           We mentioned above that polynomials with multiple roots are
           always ill-conditioned; if your input is given to `n` bits of
           precision, you should not expect more than `n/k` good bits
           for a `k`-fold root. (You can get solutions that make the
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

            sage: # needs sage.rings.number_field
            sage: K.<zeta> = CyclotomicField(2)
            sage: R.<x> = K[]
            sage: factor(x^3 - 1)
            (x - 1) * (x^2 + x + 1)

        This shows that the issue from :issue:`6237` is fixed::

            sage: R.<u> = QQ[]
            sage: g = -27*u^14 - 32*u^9
            sage: g.roots(CDF, multiplicities=False)  # abs tol 2e-15                   # needs numpy
            [-1.0345637159435719, 0.0, -0.3196977699902601 - 0.9839285635706636*I, -0.3196977699902601 + 0.9839285635706636*I,
             0.8369796279620465 - 0.6081012947885318*I, 0.8369796279620465 + 0.6081012947885318*I]
            sage: g.roots(CDF)  # abs tol 2e-15                                         # needs numpy
            [(-1.0345637159435719, 1), (0.0, 9), (-0.3196977699902601 - 0.9839285635706636*I, 1), (-0.3196977699902601 + 0.9839285635706636*I, 1),
             (0.8369796279620465 - 0.6081012947885318*I, 1), (0.8369796279620465 + 0.6081012947885318*I, 1)]

        This shows that the issue at :issue:`2418` is fixed::

            sage: x = polygen(QQ)
            sage: p = (x^50/2^100 + x^10 + x + 1).change_ring(ComplexField(106))        # needs sage.rings.real_mpfr
            sage: rts = (p/2^100).roots(multiplicities=False)                           # needs sage.libs.pari
            sage: eps = 2^(-50)   # we test the roots numerically
            sage: [abs(p(rt)) < eps for rt in rts] == [True]*50                         # needs sage.rings.number_field
            True

        This shows that the issue at :issue:`10901` is fixed::

            sage: # needs sage.symbolic
            sage: a = var('a'); R.<x> = SR[]
            sage: f = x - a
            sage: f.roots(RR)
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression to a numeric value
            sage: f.roots(CC)
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression to a numeric value

        We can find roots of polynomials defined over `\ZZ` or `\QQ`
        over the `p`-adics, see :issue:`15422`::

            sage: R.<x> = ZZ[]
            sage: pol = (x - 1)^2
            sage: pol.roots(Qp(3, 5))                                                   # needs sage.rings.padics
            [(1 + O(3^5), 2)]

        We lose precision if we first change coefficients to `\QQ_p`::

            sage: # needs sage.rings.padics
            sage: pol.change_ring(Qp(3, 5)).roots()
            [(1 + O(3^3), 2)]
            sage: (pol - 3^6).roots(Qp(3, 5))
            [(1 + 2*3^3 + 2*3^4 + O(3^5), 1), (1 + 3^3 + O(3^5), 1)]
            sage: r = pol.roots(Zp(3, 5), multiplicities=False); r
            [1 + O(3^5)]
            sage: parent(r[0])
            3-adic Ring with capped relative precision 5

        Spurious crash with pari-2.5.5, see :issue:`16165`::

            sage: f = (1+x+x^2)^3
            sage: f.roots(ring=CC)                                                      # needs sage.libs.pari sage.rings.real_mpfr
            [(-0.500000000000000 - 0.866025403784439*I, 3),
             (-0.500000000000000 + 0.866025403784439*I, 3)]

        Test a crash reported at :issue:`19649`::

            sage: polRing.<x> = PolynomialRing(ZZ)
            sage: j = (x+1)^2 * (x-1)^7 * (x^2-x+1)^5
            sage: j.roots(CC)                                                           # needs sage.libs.pari sage.rings.real_mpfr
            [(-1.00000000000000, 2),
             (1.00000000000000, 7),
             (0.500000000000000 - 0.866025403784439*I, 5),
             (0.500000000000000 + 0.866025403784439*I, 5)]

        Test that some large finite rings can be handled (:issue:`13825`)::

            sage: R.<x> = IntegerModRing(20000009)[]                                    # needs sage.libs.pari
            sage: eq = x^6+x-17
            sage: eq.roots(multiplicities=False)                                        # needs sage.libs.pari
            [3109038, 17207405]

        Test that roots in fixed modulus `p`-adic fields work (:issue:`17598`)::

            sage: len(cyclotomic_polynomial(3).roots(ZpFM(739, 566)))                   # needs sage.rings.padics
            2

        Check that :issue:`26421` is fixed::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: P.<x> = R[]
            sage: p = x^4 + (-5 - 2*t)*x^3 + (-2 + 10*t)*x^2 + (10 + 4*t)*x - 20*t
            sage: p.roots()                                                             # needs sage.libs.singular
            [(5, 1), (2*t, 1)]

        Check that :issue:`31040` is fixed::

            sage: R.<x> = QQ[]
            sage: K.<a> = Qq(3).extension(x^2 + 1)                                      # needs sage.rings.padics
            sage: (x^2 + 1).roots(K)                                                    # needs sage.rings.padics
            [(a + O(3^20), 1),
             (2*a + 2*a*3 + 2*a*3^2 + 2*a*3^3 + 2*a*3^4 + 2*a*3^5 + 2*a*3^6 + 2*a*3^7 + 2*a*3^8 + 2*a*3^9 + 2*a*3^10
                + 2*a*3^11 + 2*a*3^12 + 2*a*3^13 + 2*a*3^14 + 2*a*3^15 + 2*a*3^16 + 2*a*3^17 + 2*a*3^18 + 2*a*3^19 + O(3^20),
              1)]

        Check that :issue:`31710` is fixed::

            sage: CBF['x'].zero().roots(multiplicities=False)                           # needs sage.libs.flint
            Traceback (most recent call last):
            ...
            ArithmeticError: taking the roots of the zero polynomial

        Check that :issue:`33979` is fixed::

            sage: # needs sage.libs.pari
            sage: n = randint(2, 10^6)
            sage: K = Integers(n)
            sage: R.<x> = PolynomialRing(K)
            sage: a = randint(0, n - 1)
            sage: b = randint(0, n - 1)
            sage: f = (x - a) * (x - b)
            sage: all(r.parent() is K for r in f.roots(multiplicities=False))           # needs sage.rings.finite_rings
            True
        """
        from sage.rings.finite_rings.finite_field_constructor import GF
        K = self._parent.base_ring()

        # If the base ring has a method _roots_univariate_polynomial,
        # try to use it. An exception is raised if the method does not
        # handle the current parameters
        if hasattr(K, '_roots_univariate_polynomial'):
            try:
                return K._roots_univariate_polynomial(self, ring=ring, multiplicities=multiplicities, algorithm=algorithm, **kwds)
            except NotImplementedError:
                # This does not handle something, so keep calm and continue on
                pass

        if kwds:
            raise TypeError("roots() got unexpected keyword argument(s): {}".format(kwds.keys()))

        L = K if ring is None else ring

        input_fp = isinstance(K, (sage.rings.abc.RealField,
                                  sage.rings.abc.ComplexField,
                                  sage.rings.abc.RealDoubleField,
                                  sage.rings.abc.ComplexDoubleField))
        output_fp = isinstance(L, (sage.rings.abc.RealField,
                                  sage.rings.abc.ComplexField,
                                  sage.rings.abc.RealDoubleField,
                                  sage.rings.abc.ComplexDoubleField))
        input_complex = isinstance(K, (sage.rings.abc.ComplexField, sage.rings.abc.ComplexDoubleField))
        output_complex = isinstance(L, (sage.rings.abc.ComplexField, sage.rings.abc.ComplexDoubleField))
        input_gaussian = (isinstance(K, sage.rings.abc.NumberField_quadratic)
                          and list(K.polynomial()) == [1, 0, 1])

        if input_fp and output_fp:
            # allow for possibly using a fast but less reliable
            # floating point algorithm from numpy
            low_prec = isinstance(K, (sage.rings.abc.RealDoubleField, sage.rings.abc.ComplexDoubleField))
            if algorithm is None:
                if low_prec:
                    algorithm = 'either'
                else:
                    algorithm = 'pari'

            if algorithm != 'numpy' and algorithm != 'either' and algorithm != 'pari':
                raise ValueError("Unknown algorithm '%s'" % algorithm)

            # We should support GSL, too.  We could also support PARI's
            # old Newton-iteration algorithm.

            input_arbprec = isinstance(K, (sage.rings.abc.RealField, sage.rings.abc.ComplexField))

            if algorithm == 'numpy' or algorithm == 'either':
                if K.prec() > 53 and L.prec() > 53:
                    from warnings import warn
                    warn('NumPy does not support arbitrary precision arithmetic.  ' +
                         'The roots found will likely have less precision than ' +
                         'you expect.')

                import numpy
                from numpy.linalg import LinAlgError
                from sage.rings.complex_double import CDF

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
                    from sage.rings.real_mpfr import RR
                    from sage.rings.cc import CC
                    self = self.change_ring(CC if input_complex else RR)
                ext_rts = self.__pari__().polroots(precision=L.prec())

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

        if isinstance(L, sage.rings.abc.SymbolicRing):
            if self.degree() == 2:
                from sage.misc.functional import sqrt
                from sage.symbolic.constants import I
                coeffs = self.list()
                D = coeffs[1]*coeffs[1] - 4*coeffs[0]*coeffs[2]
                l = None
                if D > 0:
                    l = [((-coeffs[1]-sqrt(D))/2/coeffs[2], 1),
                         ((-coeffs[1]+sqrt(D))/2/coeffs[2], 1)]
                elif D < 0:
                    l = [((-coeffs[1]-I*sqrt(-D))/2/coeffs[2], 1),
                         ((-coeffs[1]+I*sqrt(-D))/2/coeffs[2], 1)]
                elif D == 0:
                    l = [(-coeffs[1]/2/coeffs[2], 2)]
                if l:
                    if multiplicities:
                        return l
                    else:
                        return [val for val,m in l]
            from sage.symbolic.ring import SR
            vname = 'do_not_use_this_name_in_a_polynomial_coefficient'
            var = SR(vname)
            expr = self(var)
            rts = expr.solve(var,
                             explicit_solutions=True,
                             multiplicities=multiplicities)
            if multiplicities:
                return [(rt.rhs(), mult) for rt, mult in zip(*rts)]
            else:
                return [rt.rhs() for rt in rts]

        if L != K or isinstance(L, sage.rings.abc.AlgebraicField_common):
            # So far, the only "special" implementations are for real
            # and complex root isolation and for p-adic factorization
            if (isinstance(K, IntegerRing_class) or isinstance(K, RationalField)
                or isinstance(K, sage.rings.abc.AlgebraicRealField)) and \
                isinstance(L, (sage.rings.abc.AlgebraicRealField, sage.rings.abc.RealIntervalField)):

                from sage.rings.polynomial.real_roots import real_roots

                if isinstance(L, sage.rings.abc.AlgebraicRealField):
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

            if (isinstance(K, IntegerRing_class) or isinstance(K, RationalField)
                or isinstance(K, sage.rings.abc.AlgebraicField_common) or input_gaussian) and \
                isinstance(L, (sage.rings.abc.ComplexIntervalField, sage.rings.abc.AlgebraicField_common)):

                from sage.rings.polynomial.complex_roots import complex_roots

                if isinstance(L, sage.rings.abc.ComplexIntervalField):
                    rts = complex_roots(self, min_prec=L.prec())
                elif isinstance(L, sage.rings.abc.AlgebraicField):
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
                if isinstance(L, sage.rings.abc.ComplexDoubleField):
                    real_field = RDF
                else:
                    from sage.rings.real_mpfr import RealField
                    real_field = RealField(L.prec())

                return self.change_ring(real_field).roots(ring=L, multiplicities=multiplicities, algorithm=algorithm)
            elif isinstance(L, (sage.rings.abc.pAdicRing, sage.rings.abc.pAdicField)) and L.absolute_degree() == 1:
                p = L.prime()
                n = L.precision_cap()
                try:
                    F = self.factor_padic(p, n)
                except AttributeError:
                    pass
                else:
                    return self.change_ring(L)._roots_from_factorization(F, multiplicities)

            try:
                self_L = self.change_ring(L)
            except (TypeError, ValueError):
                if L.is_exact() and L.is_subring(K):
                    return self._roots_in_subring(L, multiplicities, algorithm)
                else:
                    raise

            return self_L.roots(multiplicities=multiplicities, algorithm=algorithm)

        try:
            if K.is_integral_domain():
                if not K.is_field():
                    try:
                        # get rid of the content of self since we don't need it
                        # and we really don't want to factor it if it's a huge
                        # integer
                        c = self.content_ideal().gen()
                        self = self//c
                    except (AttributeError, NotImplementedError, TypeError):
                        pass
                return self._roots_from_factorization(self.factor(), multiplicities)
            else:
                raise NotImplementedError
        except NotImplementedError:
            if K.is_finite():
                if multiplicities:
                    raise NotImplementedError("root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)")
                elif isinstance(K, sage.rings.abc.IntegerModRing):
                    # handling via the chinese remainders theorem
                    N = K.cardinality()
                    primes = N.prime_divisors()
                    residue_roots = []
                    for p in primes:
                        local_self = self.change_ring(GF(p))
                        local_roots = local_self.roots(multiplicities=False)
                        residue_roots.append([a.lift() for a in local_roots])
                    result = []
                    P = ZZ.prod(primes)
                    for res in cartesian_product_iterator(residue_roots):
                        lifted = crt(list(res), primes)
                        candidate = lifted
                        for k in range(N // P):
                            if not self(candidate):
                                result.append(K(candidate))
                            candidate += P
                    return result
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
            sage: F = pol.factor(); F                                                   # needs sage.libs.pari
            2 * 5 * (x - 2) * x * (2*x - 1)
            sage: pol._roots_from_factorization(F, multiplicities=True)                 # needs sage.libs.pari
            [(2, 1), (0, 1)]
            sage: pol.change_ring(QQ)._roots_from_factorization(F, multiplicities=False)            # needs sage.libs.pari
            [2, 0, 1/2]
        """
        seq = []
        K = self._parent.base_ring()
        for fac in F:
            g = fac[0]
            if g.degree() == 1:
                try:
                    # We need to check that this root is actually in K;
                    # otherwise we'd return roots in the fraction field of K.
                    rt = K(-g[0] / g[1])
                except (ValueError, ArithmeticError, TypeError):
                    pass
                else:
                    if multiplicities:
                        seq.append((rt,fac[1]))
                    else:
                        seq.append(rt)
        return seq

    def _roots_in_subring(self, L, multiplicities, algorithm):
        r"""
        Helper method for :meth:`roots` that filters the roots belonging to
        a subring of the base ring.

        TESTS::

            sage: Pols.<n> = QQ[]
            sage: pol = (n - 1/2)^2 * (n - 1)^2 * (n - 2)
            sage: rts = pol.roots(ZZ, multiplicities=False); rts                        # needs sage.libs.pari
            [2, 1]
            sage: rts[0].parent()                                                       # needs sage.libs.pari
            Integer Ring

            sage: Pols_x.<x> = QQ[]
            sage: Pols_xy.<y> = Pols_x[]
            sage: ((y - 1)*(y - x))._roots_in_subring(QQ, True, None)                   # needs sage.libs.singular
            [(1, 1)]
        """
        K = self._parent.base_ring()
        sec = K.coerce_map_from(L).section()
        if sec is None:
            raise NotImplementedError("embedding of {} into {} has no implemented section", L, K)
        rts_K = self.roots(multiplicities=multiplicities, algorithm=algorithm)
        if not multiplicities:
            rts_K = [(rt, None) for rt in rts_K]
        rts_L = []
        for xK, mult in rts_K:
            try:
                xL = sec(xK)
            except (TypeError, ValueError):
                continue
            rts_L.append((xL, mult))
        if multiplicities:
            return rts_L
        else:
            return [rt for (rt, _) in rts_L]

    def real_roots(self):
        """
        Return the real roots of this polynomial, without multiplicities.

        Calls ``self.roots(ring=RR)``, unless this is a polynomial with
        floating-point real coefficients, in which case it calls
        ``self.roots()``.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: (x^2 - x - 1).real_roots()                                            # needs sage.libs.pari sage.rings.real_mpfr
            [-0.618033988749895, 1.61803398874989]

        TESTS::

            sage: x = polygen(RealField(100))                                           # needs sage.libs.pari sage.rings.real_mpfr
            sage: (x^2 - x - 1).real_roots()[0].parent()                                # needs sage.libs.pari sage.rings.real_mpfr
                Real Field with 100 bits of precision
            sage: x = polygen(RDF)
            sage: (x^2 - x - 1).real_roots()[0].parent()                                # needs numpy
            Real Double Field

            sage: x = polygen(ZZ,'x'); v = (x^2 - x - 1).real_roots()                   # needs sage.libs.pari sage.rings.real_mpfr
            sage: v[0].parent() is RR                                                   # needs sage.libs.pari sage.rings.real_mpfr
            True
        """
        K = self.base_ring()
        if isinstance(K, (sage.rings.abc.RealField, sage.rings.abc.RealDoubleField)):
            return self.roots(multiplicities=False)

        from sage.rings.real_mpfr import RR

        return self.roots(ring=RR, multiplicities=False)

    def complex_roots(self):
        """
        Return the complex roots of this polynomial, without
        multiplicities.

        Calls ``self.roots(ring=CC)``, unless this is a polynomial with
        floating-point coefficients, in which case it is uses the
        appropriate precision from the input coefficients.

        EXAMPLES::

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: x = polygen(ZZ)
            sage: (x^3 - 1).complex_roots()   # note: low order bits slightly different on ppc.
            [1.00000000000000,
             -0.500000000000000 - 0.86602540378443...*I,
             -0.500000000000000 + 0.86602540378443...*I]

        TESTS::

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: x = polygen(RR)
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Field with 53 bits of precision

            sage: x = polygen(RDF)
            sage: (x^3 - 1).complex_roots()[0].parent()                                 # needs numpy
            Complex Double Field

            sage: x = polygen(RealField(200))                                           # needs sage.rings.real_mpfr
            sage: (x^3 - 1).complex_roots()[0].parent()                                 # needs sage.rings.real_mpfr
            Complex Field with 200 bits of precision

            sage: x = polygen(CDF)                                                      # needs sage.rings.complex_double
            sage: (x^3 - 1).complex_roots()[0].parent()                                 # needs numpy sage.rings.complex_double
            Complex Double Field


            sage: # needs sage.rings.real_mpfr
            sage: x = polygen(ComplexField(200))
            sage: (x^3 - 1).complex_roots()[0].parent()
            Complex Field with 200 bits of precision
            sage: x = polygen(ZZ,'x'); v=(x^2-x-1).complex_roots()
            sage: v[0].parent() is CC
            True
        """
        K = self.base_ring()
        if isinstance(K, sage.rings.abc.RealField):
            from sage.rings.complex_mpfr import ComplexField
            return self.roots(ring=ComplexField(K.prec()), multiplicities=False)
        if isinstance(K, sage.rings.abc.RealDoubleField):
            from sage.rings.complex_double import CDF
            return self.roots(ring=CDF, multiplicities=False)
        if isinstance(K, (sage.rings.abc.ComplexField, sage.rings.abc.ComplexDoubleField)):
            return self.roots(multiplicities=False)
        from sage.rings.cc import CC
        return self.roots(ring=CC, multiplicities=False)

    def number_of_roots_in_interval(self, a=None, b=None):
        r"""
        Return the number of roots of this polynomial in the interval
        `[a,b]`, counted without multiplicity. The endpoints `a`, `b` default to
        ``-Infinity``, ``Infinity`` (which are also valid input values).

        Calls the PARI routine :pari:`polsturm`.

        Note that as of version 2.8, PARI includes the left endpoint
        of the interval (and no longer uses Sturm's algorithm on exact
        inputs). :pari:`polsturm` requires a polynomial with real
        coefficients; in case PARI returns an error, we try again
        after taking the GCD of ``self`` with its complex conjugate.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ)
            sage: pol = (x - 1)^2 * (x - 2)^2 * (x - 3)
            sage: pol.number_of_roots_in_interval(1, 2)
            2
            sage: pol.number_of_roots_in_interval(1.01, 2)
            1
            sage: pol.number_of_roots_in_interval(None, 2)
            2
            sage: pol.number_of_roots_in_interval(1, Infinity)
            3
            sage: pol.number_of_roots_in_interval()
            3
            sage: pol = (x - 1) * (x - 2) * (x - 3)

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: pol2 = pol.change_ring(CC)
            sage: pol2.number_of_roots_in_interval()
            3
            sage: R.<x> = PolynomialRing(CC)
            sage: pol = (x - 1) * (x - CC(I))
            sage: pol.number_of_roots_in_interval(0, 2)
            1

        TESTS::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ)
            sage: pol = (x - 1)^2 * (x - 2)^2 * (x - 3)
            sage: pol.number_of_roots_in_interval(1, 2)
            2
            sage: pol = chebyshev_T(5,x)
            sage: pol.number_of_roots_in_interval(-1, 2)
            5
            sage: pol.number_of_roots_in_interval(0, 2)
            3
        """
        pol = self // self.gcd(self.derivative())  # squarefree part
        if a is None:
            a1 = -infinity.infinity
        else:
            a1 = a
        if b is None:
            b1 = infinity.infinity
        else:
            b1 = b
        try:
            return pari(pol).polsturm([a1, b1])
        except PariError:
            # Take GCD with the conjugate, to extract the maximum factor
            # with real coefficients.
            pol2 = pol.gcd(pol.map_coefficients(lambda z: z.conjugate()))
            return pari(pol2).polsturm([a1, b1])

    def number_of_real_roots(self):
        r"""
        Return the number of real roots of this polynomial, counted
        without multiplicity.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ)
            sage: pol = (x - 1)^2 * (x - 2)^2 * (x - 3)
            sage: pol.number_of_real_roots()
            3
            sage: pol = (x - 1) * (x - 2) * (x - 3)

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: pol2 = pol.change_ring(CC)
            sage: pol2.number_of_real_roots()
            3
            sage: R.<x> = PolynomialRing(CC)
            sage: pol = (x - 1) * (x - CC(I))
            sage: pol.number_of_real_roots()
            1
        """
        return self.number_of_roots_in_interval()

    def all_roots_in_interval(self, a=None, b=None):
        r"""
        Return ``True`` if the roots of this polynomial are all real and
        contained in the given interval.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ)
            sage: pol = (x - 1)^2 * (x - 2)^2 * (x - 3)
            sage: pol.all_roots_in_interval(1, 3)
            True
            sage: pol.all_roots_in_interval(1.01, 3)
            False
            sage: pol = chebyshev_T(5, x)
            sage: pol.all_roots_in_interval(-1, 1)
            True
            sage: pol = chebyshev_T(5, x/2)
            sage: pol.all_roots_in_interval(-1, 1)
            False
            sage: pol.all_roots_in_interval()
            True
        """
        pol = self // self.gcd(self.derivative())
        return pol.number_of_roots_in_interval(a, b) == pol.degree()

    def is_real_rooted(self):
        r"""
        Return ``True`` if the roots of this polynomial are all real.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ)
            sage: pol = chebyshev_T(5, x)
            sage: pol.is_real_rooted()
            True
            sage: pol = x^2 + 1
            sage: pol.is_real_rooted()
            False
        """
        return self.all_roots_in_interval()

    def reciprocal_transform(self, R=1, q=1):
        r"""
        Transform a general polynomial into a self-reciprocal polynomial.

        The input `Q` and output `P` satisfy the relation

        .. MATH::

            P(x) = Q(x + q/x) x^{\deg(Q)} R(x).

        In this relation, `Q` has all roots in the real interval
        `[-2\sqrt{q}, 2\sqrt{q}]` if and only if `P` has all roots on the
        circle `|x| = \sqrt{q}` and `R` divides `x^2-q`.

        .. SEEALSO::

            The inverse operation is :meth:`trace_polynomial`.

        INPUT:

        - ``R`` -- polynomial
        - ``q`` -- scalar (default: `1`)

        EXAMPLES::

            sage: pol.<x> = PolynomialRing(Rationals())
            sage: u = x^2 + x - 1
            sage: u.reciprocal_transform()
            x^4 + x^3 + x^2 + x + 1
            sage: u.reciprocal_transform(R=x-1)
            x^5 - 1
            sage: u.reciprocal_transform(q=3)
            x^4 + x^3 + 5*x^2 + 3*x + 9
        """
        S = self.parent()
        x = S.gen()
        return S(x**(self.degree()) * self(x + q/x)) * R

    def trace_polynomial(self):
        r"""
        Compute the trace polynomial and cofactor.

        The input `P` and output `Q` satisfy the relation

        .. MATH::

            P(x) = Q(x + q/x) x^{\deg(Q)} R(x).

        In this relation, `Q` has all roots in the real interval
        `[-2\sqrt{q}, 2\sqrt{q}]` if and only if `P` has all roots on the
        circle `|x| = \sqrt{q}` and `R` divides `x^2-q`.  We thus require
        that the base ring of this polynomial have a coercion to the real
        numbers.

        .. SEEALSO::

            The inverse operation is :meth:`reciprocal_transform`.

        OUTPUT:

        - ``Q`` -- trace polynomial
        - ``R`` -- cofactor
        - ``q`` -- scaling factor

        EXAMPLES::

            sage: pol.<x> = PolynomialRing(Rationals())
            sage: u = x^5 - 1; u.trace_polynomial()
            (x^2 + x - 1, x - 1, 1)
            sage: u = x^4 + x^3 + 5*x^2 + 3*x + 9
            sage: u.trace_polynomial()
            (x^2 + x - 1, 1, 3)

        We check that this function works for rings
        that have a coercion to the reals::

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(x^2 - 2, embedding=1.4)
            sage: u = x^4 + a*x^3 + 3*x^2 + 2*a*x + 4
            sage: u.trace_polynomial()
            (x^2 + a*x - 1, 1, 2)
            sage: (u*(x^2-2)).trace_polynomial()
            (x^2 + a*x - 1, x^2 - 2, 2)
            sage: (u*(x^2-2)^2).trace_polynomial()
            (x^4 + a*x^3 - 9*x^2 - 8*a*x + 8, 1, 2)
            sage: (u*(x^2-2)^3).trace_polynomial()
            (x^4 + a*x^3 - 9*x^2 - 8*a*x + 8, x^2 - 2, 2)
            sage: u = x^4 + a*x^3 + 3*x^2 + 4*a*x + 16
            sage: u.trace_polynomial()
            (x^2 + a*x - 5, 1, 4)
            sage: (u*(x-2)).trace_polynomial()
            (x^2 + a*x - 5, x - 2, 4)
            sage: (u*(x+2)).trace_polynomial()
            (x^2 + a*x - 5, x + 2, 4)

        TESTS:

        Check that :issue:`28395` is fixed::

            sage: P.<t> = QQ[]
            sage: u = t^4 + 3*t^2 + 1
            sage: u.trace_polynomial()
            (t^2 + 1, 1, 1)
        """
        S = self.parent()
        A = S.base_ring()
        x = S.gen()
        if self.get_coeff_c(0) == 0:
            raise ValueError("Polynomial not self-reciprocal")
        d = self.degree()
        # previous if statement implies d >= 0
        val = self.get_unsafe(0) / self.get_unsafe(d)
        sg = val.sign()
        try:
            q = A(abs(val)**(2/d))
        except (TypeError, ValueError):
            raise ValueError("Polynomial not self-reciprocal")
        for i in range(d//2 + 1):
            if self.get_unsafe(d - i) != sg * self.get_unsafe(i) / q**(d/2-i):
                raise ValueError("Polynomial not self-reciprocal")
        Q = self
        if sg == -1 and Q.degree() % 2 == 0:
            cofactor = x**2 - q
        elif sg == -1:
            cofactor = x - q.sqrt()
        elif Q.degree() % 2:
            cofactor = x + q.sqrt()
        else:
            cofactor = S.one()
        Q //= cofactor
        coeffs = []
        m = Q.degree() // 2
        for i in reversed(range(m + 1)):
            coeffs.insert(0, Q[2*i]) # Note: degree of Q may be less than 2*i
            Q = (Q % (x**2 + q)**i) // x
        return S(coeffs), cofactor, q

    def is_weil_polynomial(self, return_q=False):
        r"""
        Return ``True`` if this is a Weil polynomial.

        This polynomial must have rational or integer coefficients.

        INPUT:

        - ``self`` -- polynomial with rational or integer coefficients

        - ``return_q`` -- (default: ``False``) if ``True``, return a second
          value `q` which is the prime power with respect to which this is
          `q`-Weil, or 0 if there is no such value

        EXAMPLES::

            sage: polRing.<x> = PolynomialRing(Rationals())
            sage: P0 = x^4 + 5*x^3 + 15*x^2 + 25*x + 25
            sage: P1 = x^4 + 25*x^3 + 15*x^2 + 5*x + 25
            sage: P2 = x^4 + 5*x^3 + 25*x^2 + 25*x + 25
            sage: P0.is_weil_polynomial(return_q=True)                                  # needs sage.libs.pari
            (True, 5)
            sage: P0.is_weil_polynomial(return_q=False)                                 # needs sage.libs.pari
            True
            sage: P1.is_weil_polynomial(return_q=True)
            (False, 0)
            sage: P1.is_weil_polynomial(return_q=False)
            False
            sage: P2.is_weil_polynomial()                                               # needs sage.libs.pari
            False

        .. SEEALSO::

            Polynomial rings have a method :meth:`weil_polynomials` to compute sets of Weil
            polynomials. This computation uses the iterator
            :class:`sage.rings.polynomial.weil.weil_polynomials.WeilPolynomials`.

        TESTS:

        Check that :issue:`28395` is fixed::

            sage: P.<t> = QQ[]
            sage: u = t^10 + 4*t^9 + 8*t^8 + 18*t^7 + 81*t^6 + 272*t^5 + 567*t^4 + 882*t^3 + 2744*t^2 + 9604*t + 16807
            sage: u.is_weil_polynomial()                                                # needs sage.libs.pari
            True

        AUTHORS:

        David Zureick-Brown (2017-10-01)
        """
        from sage.rings.rational_field import QQ
        if not QQ.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError
        polRing = self.parent()
        x = polRing.gen()
        # the following is the polynomial whose roots are the squares
        # of the roots of self.
        Q = polRing(list(self(x)*self(-x))[::2])
        try:
            Q, _, q = Q.trace_polynomial()
        except ValueError:
            b = False
        else:
            b = Q.all_roots_in_interval(-2*q.sqrt(), 2*q.sqrt())
        if return_q:
            return (b, self.base_ring()(q.sqrt())) if b else (b, 0)
        else:
            return b

    def is_lorentzian(self, explain=False):
        r"""
        Return ``True`` if this is a Lorentzian polynomial.

        A univariate real polynomial is Lorentzian if and only if it is a
        monomial with positive coefficient, or zero.  The definition is more
        involved for multivariate real polynomials.

        INPUT:

        - ``explain`` -- boolean (default: ``False``); if ``True``
          return a tuple whose first element is the boolean result of the test,
          and the second element is a string describing the reason the test failed,
          or ``None`` if the test succeeded

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: p1 = x^2
            sage: p1.is_lorentzian()
            True
            sage: p2 = 1 + x^2
            sage: p2.is_lorentzian()
            False
            sage: p3 = P.zero()
            sage: p3.is_lorentzian()
            True
            sage: p4 = -2*x^3
            sage: p4.is_lorentzian()
            False

        It is an error to check if a polynomial is Lorentzian if its base ring
        is not a subring of the real numbers, as the notion is not defined in
        this case::

            sage: # needs sage.rings.real_mpfr
            sage: Q.<y> = CC[]
            sage: q = y^2
            sage: q.is_lorentzian()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_lorentzian only implemented for real polynomials

        The method can give a reason for a polynomial failing to be Lorentzian::

            sage: p = x^2 + 2*x
            sage: p.is_lorentzian(explain=True)
            (False, 'inhomogeneous')

        REFERENCES:

        For full definitions and related discussion, see [BrHu2019]_ and
        [HMMS2019]_.
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        R = PolynomialRing(self.base_ring(), 1, [self.variable_name()])
        return R(self).is_lorentzian(explain=explain)

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
        return self._parent.variable_name()

    @coerce_binop
    def xgcd(self, other):
        r"""
        Return an extended gcd of this polynomial and ``other``.

        INPUT:

        - ``other`` -- a polynomial in the same ring as this polynomial

        OUTPUT:

        A tuple ``(r, s, t)`` where ``r`` is a greatest common divisor
        of this polynomial and ``other``, and ``s`` and ``t`` are such
        that ``r = s*self + t*other`` holds.

        .. NOTE::

            The actual algorithm for computing the extended gcd depends on the
            base ring underlying the polynomial ring. If the base ring defines
            a method :meth:`_xgcd_univariate_polynomial`, then this method will be
            called (see examples below).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQbar[]
            sage: (2*x^2).gcd(2*x)
            x
            sage: R.zero().gcd(0)
            0
            sage: (2*x).gcd(0)
            x

        One can easily add xgcd functionality to new rings by providing a
        method :meth:`_xgcd_univariate_polynomial`::

            sage: R.<x> = QQ[]
            sage: S.<y> = R[]
            sage: h1 = y*x
            sage: h2 = y^2*x^2
            sage: h1.xgcd(h2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Univariate Polynomial Ring in x over Rational Field
            does not provide an xgcd implementation for univariate polynomials
            sage: T.<x,y> = QQ[]
            sage: def poor_xgcd(f, g):
            ....:     ret = S(T(f).gcd(g))
            ....:     if ret == f: return ret, S.one(), S.zero()
            ....:     if ret == g: return ret, S.zero(), S.one()
            ....:     raise NotImplementedError
            sage: R._xgcd_univariate_polynomial = poor_xgcd
            sage: h1.xgcd(h2)
            (x*y, 1, 0)
            sage: del R._xgcd_univariate_polynomial
        """
        if hasattr(self.base_ring(), '_xgcd_univariate_polynomial'):
            return self.base_ring()._xgcd_univariate_polynomial(self, other)
        else:
            raise NotImplementedError("%s does not provide an xgcd implementation for univariate polynomials"%self.base_ring())

    def rational_reconstruction(self, m, n_deg=None, d_deg=None):
        r"""
        Return a tuple of two polynomials ``(n, d)``
        where ``self * d`` is congruent to ``n`` modulo ``m`` and
        ``n.degree() <= n_deg`` and ``d.degree() <= d_deg``.

        INPUT:

        - ``m`` -- a univariate polynomial

        - ``n_deg`` -- (optional) an integer; the default is `\lfloor (\deg(m) - 1)/2 \rfloor`

        - ``d_deg`` -- (optional) an integer; the default is `\lfloor (\deg(m) - 1)/2 \rfloor`

        ALGORITHM:

        The algorithm is based on the extended Euclidean algorithm for the polynomial greatest common divisor.

        EXAMPLES:

        Over `\QQ[z]`::

            sage: z  = PolynomialRing(QQ, 'z').gen()
            sage: p = -z**16 - z**15 - z**14 + z**13 + z**12 + z**11 - z**5 - z**4 - z**3 + z**2 + z + 1
            sage: m = z**21
            sage: n, d = p.rational_reconstruction(m); n, d
            (z^4 + 2*z^3 + 3*z^2 + 2*z + 1,
             z^10 + z^9 + z^8 + z^7 + z^6 + z^5 + z^4 + z^3 + z^2 + z + 1)
            sage: ((p*d - n) % m).is_zero()
            True

        Over `\ZZ[z]`::

            sage: z  = PolynomialRing(ZZ, 'z').gen()
            sage: p = -z**16 - z**15 - z**14 + z**13 + z**12 + z**11 - z**5 - z**4 - z**3 + z**2 + z + 1
            sage: m = z**21
            sage: n, d = p.rational_reconstruction(m); n, d
            (z^4 + 2*z^3 + 3*z^2 + 2*z + 1,
             z^10 + z^9 + z^8 + z^7 + z^6 + z^5 + z^4 + z^3 + z^2 + z + 1)
            sage: ((p*d - n) % m).is_zero()
            True

        Over an integral domain, ``d`` might not be monic::

            sage: P = PolynomialRing(ZZ, 'x')
            sage: x = P.gen()
            sage: p = 7*x^5 - 10*x^4 + 16*x^3 - 32*x^2 + 128*x + 256
            sage: m = x^5
            sage: n, d = p.rational_reconstruction(m, 3, 2); n, d
            (-32*x^3 + 384*x^2 + 2304*x + 2048, 5*x + 8)
            sage: ((p*d - n) % m).is_zero()
            True
            sage: n, d = p.rational_reconstruction(m, 4, 0); n, d
            (-10*x^4 + 16*x^3 - 32*x^2 + 128*x + 256, 1)
            sage: ((p*d - n) % m).is_zero()
            True

        Over `\QQ(t)[z]`::

            sage: P = PolynomialRing(QQ, 't')
            sage: t = P.gen()
            sage: Pz = PolynomialRing(P.fraction_field(), 'z')
            sage: z = Pz.gen()
            sage: # p = (1 + t^2*z + z^4) / (1 - t*z)
            sage: p = (1 + t^2*z + z^4)*(1 - t*z).inverse_mod(z^9)
            sage: m = z^9
            sage: n, d = p.rational_reconstruction(m); n, d
            (-1/t*z^4 - t*z - 1/t, z - 1/t)
            sage: ((p*d - n) % m).is_zero()
            True
            sage: w = PowerSeriesRing(P.fraction_field(), 'w').gen()
            sage: n = -10*t^2*z^4 + (-t^2 + t - 1)*z^3 + (-t - 8)*z^2 + z + 2*t^2 - t
            sage: d = z^4 + (2*t + 4)*z^3 + (-t + 5)*z^2 + (t^2 + 2)*z + t^2 + 2*t + 1
            sage: prec = 9
            sage: x = n.subs(z=w)/d.subs(z=w) + O(w^prec)

            sage: # needs sage.libs.flint (otherwise timeout)
            sage: nc, dc = Pz(x.list()).rational_reconstruction(z^prec)
            sage: (nc, dc) == (n, d)
            True

        Over `\QQ[t][z]`::

            sage: P = PolynomialRing(QQ, 't')
            sage: t = P.gen()
            sage: z = PolynomialRing(P, 'z').gen()
            sage: # p = (1 + t^2*z + z^4) / (1 - t*z) mod z^9
            sage: p = (1 + t^2*z + z^4) * sum((t*z)**i for i in range(9))
            sage: m = z^9
            sage: n, d = p.rational_reconstruction(m,); n, d
            (-z^4 - t^2*z - 1, t*z - 1)
            sage: ((p*d - n) % m).is_zero()
            True

        Over `\QQ_5`::

            sage: # needs sage.rings.padics
            sage: x = PolynomialRing(Qp(5), 'x').gen()
            sage: p = 4*x^5 + 3*x^4 + 2*x^3 + 2*x^2 + 4*x + 2
            sage: m = x^6
            sage: n, d = p.rational_reconstruction(m, 3, 2)
            sage: ((p*d - n) % m).is_zero()
            True

        Can also be used to obtain known Padé approximations::

            sage: z = PowerSeriesRing(QQ, 'z').gen()
            sage: P = PolynomialRing(QQ, 'x')
            sage: x = P.gen()
            sage: p = P(z.exp().list())
            sage: m = x^5
            sage: n, d = p.rational_reconstruction(m, 4, 0); n, d
            (1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1, 1)
            sage: ((p*d - n) % m).is_zero()
            True
            sage: m = x^3
            sage: n, d = p.rational_reconstruction(m, 1, 1); n, d
            (-x - 2, x - 2)
            sage: ((p*d - n) % m).is_zero()
            True
            sage: p = P(log(1-z).list())
            sage: m = x^9
            sage: n, d = p.rational_reconstruction(m, 4, 4); n, d
            (25/6*x^4 - 130/3*x^3 + 105*x^2 - 70*x, x^4 - 20*x^3 + 90*x^2 - 140*x + 70)
            sage: ((p*d - n) % m).is_zero()
            True
            sage: p = P(sqrt(1+z).list())
            sage: m = x^6
            sage: n, d = p.rational_reconstruction(m, 3, 2); n, d
            (1/6*x^3 + 3*x^2 + 8*x + 16/3, x^2 + 16/3*x + 16/3)
            sage: ((p*d - n) % m).is_zero()
            True
            sage: p = P((2*z).exp().list())
            sage: m = x^7
            sage: n, d = p.rational_reconstruction(m, 3, 3); n, d
            (-x^3 - 6*x^2 - 15*x - 15, x^3 - 6*x^2 + 15*x - 15)
            sage: ((p*d - n) % m).is_zero()
            True

        Over `\RR[z]`::

            sage: # needs sage.rings.real_mpfr
            sage: z = PowerSeriesRing(RR, 'z').gen()
            sage: P = PolynomialRing(RR, 'x')
            sage: x = P.gen()
            sage: p = P((2*z).exp().list())
            sage: m = x^7
            sage: n, d = p.rational_reconstruction(m, 3, 3); n, d  # absolute tolerance 1e-10
            (-x^3 - 6.0*x^2 - 15.0*x - 15.0, x^3 - 6.0*x^2 + 15.0*x - 15.0)

        .. SEEALSO::

            * :mod:`sage.matrix.berlekamp_massey`,
            * :meth:`sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint.rational_reconstruction`
        """
        P = self.parent()
        if not P.base_ring().is_field():
            if not P.base_ring().is_integral_domain():
                raise NotImplementedError("rational_reconstruction() "
                        "is only implemented when the base ring is a field "
                        "or a integral domain, "
                        "a workaround is to do a multimodular approach")
            Pf = P.base_extend(P.base_ring().fraction_field())
            sF = Pf(self)
            mF = Pf(m)
            n, d = sF.rational_reconstruction(mF, n_deg, d_deg)
            l = lcm([n.denominator(), d.denominator()])
            n *= l
            d *= l
            return P(n), P(d)

        # n and d are unique if m.degree() > (n.degree() + d.degree())
        if n_deg is None:
            n_deg = (m.degree() - 1) // 2
        if d_deg is None:
            d_deg = (m.degree() - 1) // 2

        if n_deg < 0 or d_deg < 0:
            raise ValueError("the degree bounds "
                    "n_deg and d_deg should be positive")

        #XGCD until degree the degree of t1 surpasses the degree of n
        s0 = P(0)
        t0 = P(1)
        s1 = P(m)
        t1 = self.mod(s1)

        while n_deg < t1.degree() and t1 != 0:
            q, r1 = s1.quo_rem(t1)
            r0 = s0 - q*t0
            s0 = t0
            s1 = t1
            t0 = r0
            t1 = r1

        assert t0 != 0
        if d_deg < t0.degree():
            raise ValueError("could not complete rational reconstruction")

        # make the denominator monic
        c = t0.leading_coefficient()
        t0 = t0.monic()
        t1 = t1 / c
        return t1, t0

    rational_reconstruct = deprecated_function_alias(12696, rational_reconstruction)

    def variables(self):
        """
        Return the tuple of variables occurring in this polynomial.

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
        Return the generator of this polynomial ring, which is the (only)
        argument used when calling ``self``.

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
        is the largest power of `p` which divides ``self``.

        The valuation at `\infty` is ``-self.degree()``.

        EXAMPLES::

            sage: P.<x> = ZZ[]
            sage: (x^2 + x).valuation()
            1
            sage: (x^2 + x).valuation(x + 1)
            1
            sage: (x^2 + 1).valuation()
            0
            sage: (x^3 + 1).valuation(infinity)
            -3
            sage: P(0).valuation()
            +Infinity
        """
        cdef int k
        cdef Polynomial _p

        if not self:
            return infinity.infinity

        if p is None:
            for k from 0 <= k <= self.degree():
                if self.get_unsafe(k):
                    return ZZ(k)
        if isinstance(p, Polynomial):
            _p = self._parent.coerce(p)
        elif p is infinity.infinity:
            return -self.degree()
        elif isinstance(p, Ideal_generic) and p.ring() is self._parent: # eventually need to handle fractional ideals in the fraction field
            if self._parent.base_ring().is_field(): # common case
                _p = p.gen()
            else:
                raise NotImplementedError
        else:
            from sage.rings.fraction_field import FractionField_generic
            if isinstance(p.parent(), FractionField_generic) and self._parent.has_coerce_map_from(p.parent().ring()):
                _p = self._parent.coerce(p.parent().ring()(p)) # here we require that p be integral.
            else:
                raise TypeError("The polynomial, p, must have the same parent as self.")

        if _p.degree() == 0:
            raise ArithmeticError("The polynomial, p, must be non-constant.")
        k = -1
        cdef Polynomial rem = self._parent.zero()
        while rem.is_zero():
            self, rem = self.quo_rem(_p)
            k += 1
        return sage.rings.integer.Integer(k)

    def ord(self, p=None):
        r"""
        This is the same as the valuation of ``self`` at `p`. See the
        documentation for :meth:`valuation`.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x^2 + x).ord(x + 1)
            1
        """
        return self.valuation(p)

    def add_bigoh(self, prec):
        r"""
        Return the power series of precision at most ``prec`` got by adding
        `O(q^\text{prec})` to self, where `q` is its variable.

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
        return self._parent.completion(self._parent.gen())(self).add_bigoh(prec)

    @cached_method
    def is_irreducible(self):
        r"""
        Return whether this polynomial is irreducible.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x^3 + 1).is_irreducible()                                            # needs sage.libs.pari
            False
            sage: (x^2 - 1).is_irreducible()                                            # needs sage.libs.pari
            False
            sage: (x^3 + 2).is_irreducible()                                            # needs sage.libs.pari
            True
            sage: R(0).is_irreducible()
            False

        The base ring does matter:  for example, `2x` is irreducible as a
        polynomial in `\QQ[x]`, but not in `\ZZ[x]`::

            sage: R.<x> = ZZ[]
            sage: R(2*x).is_irreducible()                                               # needs sage.libs.pari
            False
            sage: R.<x> = QQ[]
            sage: R(2*x).is_irreducible()                                               # needs sage.libs.pari
            True

        TESTS::

            sage: # needs sage.rings.number_field
            sage: F.<t> = NumberField(x^2 - 5)
            sage: Fx.<xF> = PolynomialRing(F)
            sage: f = Fx([2*t - 5, 5*t - 10, 3*t - 6, -t, -t + 2, 1])
            sage: f.is_irreducible()
            False
            sage: f = Fx([2*t - 3, 5*t - 10, 3*t - 6, -t, -t + 2, 1])
            sage: f.is_irreducible()
            True

        If the base ring implements `_is_irreducible_univariate_polynomial`,
        then this method gets used instead of the generic algorithm which just
        factors the input::

            sage: R.<x> = QQbar[]                                                       # needs sage.rings.number_field
            sage: hasattr(QQbar, "_is_irreducible_univariate_polynomial")               # needs sage.rings.number_field
            True
            sage: (x^2 + 1).is_irreducible()                                            # needs sage.rings.number_field
            False

        Constants can be irreducible if they are not units::

            sage: R.<x> = ZZ[]
            sage: R(1).is_irreducible()
            False
            sage: R(4).is_irreducible()
            False
            sage: R(5).is_irreducible()
            True

        Check that caching works::

            sage: R.<x> = ZZ[]
            sage: x.is_irreducible()                                                    # needs sage.libs.pari
            True
            sage: x.is_irreducible.cache                                                # needs sage.libs.pari
            True
        """
        if self.is_zero():
            return False
        if self.is_unit():
            return False
        if self.degree() == 0:
            return self.base_ring()(self).is_irreducible()

        B = self.parent().base_ring()
        if hasattr(B, '_is_irreducible_univariate_polynomial'):
            return B._is_irreducible_univariate_polynomial(self)

        F = self.factor()
        if len(F) > 1 or F[0][1] > 1:
            return False
        return True

    def shift(self, n):
        r"""
        Return this polynomial multiplied by the power `x^n`. If
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
            output = [self.base_ring().zero()] * n
            output.extend(self.coefficients(sparse=False))
            return self._new_generic(output)
        if n < 0:
            if n > self.degree():
                return self._new_generic([])
            else:
                return self._new_generic(self.coefficients(sparse=False)[-int(n):])

    def __lshift__(self, k):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x + 2
            sage: f << 3
            x^4 + 2*x^3
        """
        return self.shift(k)

    def __rshift__(self, k):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^4 + 2*x^3
            sage: f >> 3
            x + 2
        """
        return self.shift(-k)

    cpdef Polynomial truncate(self, long n):
        r"""
        Return the polynomial of degree `< n` which is equivalent
        to ``self`` modulo `x^n`.

        EXAMPLES::

            sage: R.<x> = ZZ[]; S.<y> = PolynomialRing(R, sparse=True)
            sage: f = y^3 + x*y - 3*x; f
            y^3 + x*y - 3*x
            sage: f.truncate(2)
            x*y - 3*x
            sage: f.truncate(1)
            -3*x
            sage: f.truncate(0)
            0
        """
        # __getitem__ already returns a polynomial!!
        # We must not have check=False, since 0 must not have _coeffs = [0].
        return <Polynomial>self._parent(self[:n])  # , check=False)

    cdef _inplace_truncate(self, long prec):
        return self.truncate(prec)

    @cached_method
    def is_squarefree(self):
        """
        Return ``False`` if this polynomial is not square-free, i.e., if there is a
        non-unit `g` in the polynomial ring such that `g^2` divides ``self``.

        .. WARNING::

            This method is not consistent with
            :meth:`.squarefree_decomposition` since the latter does not factor
            the content of a polynomial. See the examples below.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = (x-1) * (x-2) * (x^2-5) * (x^17-3); f
            x^21 - 3*x^20 - 3*x^19 + 15*x^18 - 10*x^17 - 3*x^4 + 9*x^3 + 9*x^2 - 45*x + 30
            sage: f.is_squarefree()
            True
            sage: (f * (x^2-5)).is_squarefree()
            False

        A generic implementation is available, which relies on gcd
        computations::

            sage: # needs sage.libs.pari
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
            sage: (2*x*y).is_squarefree()
            True
            sage: (2*x*y^2).is_squarefree()
            False

        In positive characteristic, we compute the square-free
        decomposition or a full factorization, depending on which is
        available::

            sage: K.<t> = FunctionField(GF(3))
            sage: R.<x> = K[]
            sage: (x^3 - x).is_squarefree()
            True
            sage: (x^3 - 1).is_squarefree()                                             # needs sage.libs.pari
            False
            sage: (x^3 + t).is_squarefree()                                             # needs sage.libs.pari
            True
            sage: (x^3 + t^3).is_squarefree()                                           # needs sage.libs.pari
            False

        In the following example, `t^2` is a unit in the base field::

            sage: R(t^2).is_squarefree()
            True

        This method is not consistent with :meth:`.squarefree_decomposition`::

            sage: R.<x> = ZZ[]
            sage: f = 4 * x
            sage: f.is_squarefree()                                                     # needs sage.libs.pari
            False
            sage: f.squarefree_decomposition()                                          # needs sage.libs.pari
            (4) * x

        If you want this method equally not to consider the content, you can
        remove it as in the following example::

            sage: c = f.content()
            sage: (f/c).is_squarefree()                                                 # needs sage.libs.pari
            True

        If the base ring is not an integral domain, the question is not
        mathematically well-defined::

            sage: R.<x> = IntegerModRing(9)[]
            sage: pol = (x + 3) * (x + 6); pol
            x^2
            sage: pol.is_squarefree()
            Traceback (most recent call last):
            ...
            TypeError: is_squarefree() is not defined for
            polynomials over Ring of integers modulo 9

        TESTS:

        Check that the results are cached::

            sage: R.<x> = ZZ[]
            sage: f = x^2
            sage: f.is_squarefree()                                                     # needs sage.libs.pari
            False
            sage: f.is_squarefree.cache
            False

        If the base ring implements `_is_squarefree_univariate_polynomial`,
        then this method gets used instead of the generic algorithm in
        :meth:`_is_squarefree_generic`::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQbar[]
            sage: (x^2).is_squarefree()
            False
            sage: hasattr(QQbar, '_is_squarefree_univariate_polynomial')
            False
            sage: QQbar._is_squarefree_univariate_polynomial = lambda self: True
            sage: (x^2).is_squarefree()
            True
            sage: del(QQbar._is_squarefree_univariate_polynomial)
        """
        B = self._parent.base_ring()
        if B not in sage.categories.integral_domains.IntegralDomains():
            raise TypeError("is_squarefree() is not defined for polynomials over {}".format(B))

        B = self.parent().base_ring()
        if hasattr(B, '_is_squarefree_univariate_polynomial'):
            return B._is_squarefree_univariate_polynomial(self)

        return self._is_squarefree_generic()

    def _is_squarefree_generic(self):
        r"""
        Return False if this polynomial is not square-free, i.e., if there is a
        non-unit `g` in the polynomial ring such that `g^2` divides ``self``.

        EXAMPLES::

            sage: R.<x> = QQbar[]                                                       # needs sage.rings.number_field
            sage: (x^2*(x + 1)).is_squarefree()  # indirect doctest                     # needs sage.rings.number_field
            False
            sage: (x*(x+1)).is_squarefree()  # indirect doctest                         # needs sage.rings.number_field
            True
        """
        B = self.parent().base_ring()

        # a square-free polynomial has a square-free content
        if not B.is_field():
            content = self.content_ideal().gen()
            if not content.is_squarefree():
                return False

        # separable polynomials are square-free
        if self.derivative().gcd(self).is_constant():
            return True

        # for characteristic zero rings, square-free polynomials have to be separable
        if B.characteristic().is_zero():
            return False

        # over rings of positive characteristic, we rely on the square-free decomposition if available
        try:
            F = self.squarefree_decomposition()
        # We catch:
        # - NotImplementedError in case squarefree decomposition is not implemented
        # - AttributeError in case p-th roots are not (or do not exist)
        except (NotImplementedError, AttributeError):
            F = self.factor()
        return all(e <= 1 for (f, e) in F)

    def radical(self):
        """
        Return the radical of ``self``.

        Over a field, this is the product of
        the distinct irreducible factors of ``self``. (This is also sometimes
        called the "square-free part" of ``self``, but that term is ambiguous;
        it is sometimes used to mean the quotient of ``self`` by its maximal
        square factor.)

        EXAMPLES::

            sage: P.<x> = ZZ[]
            sage: t = (x^2-x+1)^3 * (3*x-1)^2
            sage: t.radical()
            3*x^3 - 4*x^2 + 4*x - 1
            sage: radical(12 * x^5)
            6*x

        If ``self`` has a factor of multiplicity divisible by the
        characteristic (see :issue:`8736`)::

            sage: P.<x> = GF(2)[]
            sage: (x^3 + x^2).radical()                                                 # needs sage.rings.finite_rings
            x^2 + x
        """
        P = self._parent
        R = P.base_ring()
        p = R.characteristic()
        if p == 0 or p > self.degree():
            if R.is_field():
                return self // self.gcd(self.derivative())
            else:
                # Be careful with the content: return the
                # radical of the content times the radical of
                # (self/content)
                content = self.content_ideal().gen()
                self_1 = (self//content)
                return (self_1 // self_1.gcd(self_1.derivative())) * content.radical()
        else:  # The above method is not always correct (see Issue 8736)
            return self.factor().radical_value()

    def content_ideal(self):
        """
        Return the content ideal of this polynomial, defined as the ideal
        generated by its coefficients.

        EXAMPLES::

            sage: R.<x> = IntegerModRing(4)[]
            sage: f = x^4 + 3*x^2 + 2
            sage: f.content_ideal()
            Ideal (2, 3, 1) of Ring of integers modulo 4

        When the base ring is a gcd ring, the content as a ring element is
        the generator of the content ideal::

            sage: R.<x> = ZZ[]
            sage: f = 2*x^3 - 4*x^2 + 6*x - 10
            sage: f.content_ideal().gen()
            2
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

        - ``p`` -- positive integer or +infinity; the degree of the norm

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: R.<x> = RR[]
            sage: f = x^6 + x^2 + -x^4 - 2*x^3
            sage: f.norm(2)
            2.64575131106459
            sage: (sqrt(1^2 + 1^2 + (-1)^2 + (-2)^2)).n()                               # needs sage.symbolic
            2.64575131106459

        ::

            sage: f.norm(1)                                                             # needs sage.rings.real_mpfr
            5.00000000000000
            sage: f.norm(infinity)                                                      # needs sage.rings.real_mpfr
            2.00000000000000

        ::

            sage: f.norm(-1)                                                            # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            ValueError: The degree of the norm must be positive

        TESTS::

            sage: # needs sage.rings.real_mpfr
            sage: R.<x> = RR[]
            sage: f = x^6 + x^2 + -x^4 -x^3
            sage: f.norm(int(2))
            2.00000000000000

        Check that :issue:`18600` is fixed::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^2^100 + 1).norm(1)                                                 # needs sage.rings.real_mpfr
            2.00000000000000

        AUTHORS:

        - Didier Deshommes
        - William Stein: fix bugs, add definition, etc.
        """
        from sage.rings.real_mpfr import RR

        if p <= 0:
            raise ValueError("The degree of the norm must be positive")

        coeffs = self.coefficients()
        if p == infinity.infinity:
            return RR(max([abs(i) for i in coeffs]))

        p = sage.rings.integer.Integer(p)  # because we'll do 1/p below.

        if p == 1:
            return RR(sum([abs(i) for i in coeffs]))

        return RR(sum([abs(i)**p for i in coeffs]))**(1/p)

    cpdef long number_of_terms(self) noexcept:
        """
        Return the number of nonzero coefficients of ``self``.

        Also called weight, Hamming weight or sparsity.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^3 - x
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x + 1)^100
            sage: f.number_of_terms()
            101
            sage: S = GF(5)['y']
            sage: S(f).number_of_terms()
            5
            sage: cyclotomic_polynomial(105).number_of_terms()
            33

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        cdef long w = 0
        for a in self.coefficients(sparse=False):
            if a:
                w += 1
        return w

    # alias hamming_weight for number_of_terms:
    hamming_weight = number_of_terms

    def map_coefficients(self, f, new_base_ring=None):
        """
        Return the polynomial obtained by applying ``f`` to the nonzero
        coefficients of ``self``.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        polynomial will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as ``self``. Set
        ``new_base_ring`` to override this behaviour.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of ``self``

        - ``new_base_ring`` -- (optional) if given, the resulting polynomial
          will be defined over this ring

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^2 + 2
            sage: f.map_coefficients(lambda a: a + 42)
            43*x^2 + 44
            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: f = x^(2^32) + 2
            sage: f.map_coefficients(lambda a: a + 42)
            43*x^4294967296 + 44

            sage: # needs sage.symbolic
            sage: R.<x> = SR[]
            sage: f = (1+I)*x^2 + 3*x - I
            sage: f.map_coefficients(lambda z: z.conjugate())
            (-I + 1)*x^2 + 3*x + I
            sage: R.<x> = PolynomialRing(SR, sparse=True)
            sage: f = (1+I)*x^(2^32) - I
            sage: f.map_coefficients(lambda z: z.conjugate())
            (-I + 1)*x^4294967296 + I

        Examples with different base ring::

            sage: R.<x> = ZZ[]
            sage: k = GF(2)
            sage: residue = lambda x: k(x)
            sage: f = 4*x^2 + x + 3
            sage: g = f.map_coefficients(residue); g
            x + 1
            sage: g.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: g = f.map_coefficients(residue, new_base_ring=k); g
            x + 1
            sage: g.parent()                                                            # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
            sage: residue = k.coerce_map_from(ZZ)
            sage: g = f.map_coefficients(residue); g
            x + 1
            sage: g.parent()                                                            # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
        """
        R = self._parent
        if new_base_ring is not None:
            R = R.change_ring(new_base_ring)
        elif isinstance(f, Map):
            R = R.change_ring(f.codomain())
        return R({k: f(v) for k, v in self.monomial_coefficients().items()})

    def is_cyclotomic(self, certificate=False, algorithm='pari'):
        r"""
        Test if this polynomial is a cyclotomic polynomial.

        A *cyclotomic polynomial* is a monic, irreducible polynomial such that
        all roots are roots of unity.

        By default the answer is a boolean. But if ``certificate`` is ``True``,
        the result is a nonnegative integer: it is ``0`` if ``self`` is not
        cyclotomic, and a positive integer ``n`` if ``self`` is the `n`-th
        cyclotomic polynomial.

        .. SEEALSO::

            :meth:`is_cyclotomic_product`
            :meth:`cyclotomic_part`
            :meth:`has_cyclotomic_factor`

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); only works with
          ``algorithm`` set to ``'pari'``

        - ``algorithm`` -- either ``'pari'`` (default) or ``'sage'``

        ALGORITHM:

        The native algorithm implemented in Sage uses the first
        algorithm of [BD1989]_. The algorithm in PARI (using
        :pari:`poliscyclo`) is more subtle since it does compute the
        inverse of the Euler `\phi` function to determine the `n` such
        that the polynomial is the `n`-th cyclotomic polynomial.

        EXAMPLES:

        Quick tests::

            sage: # needs sage.libs.pari
            sage: P.<x> = ZZ['x']
            sage: (x - 1).is_cyclotomic()
            True
            sage: (x + 1).is_cyclotomic()
            True
            sage: (x^2 - 1).is_cyclotomic()
            False
            sage: (x^2 + x + 1).is_cyclotomic(certificate=True)
            3
            sage: (x^2 + 2*x + 1).is_cyclotomic(certificate=True)
            0

        Test first 100 cyclotomic polynomials::

            sage: all(cyclotomic_polynomial(i).is_cyclotomic() for i in range(1, 101))  # needs sage.libs.pari
            True

        Some more tests::

            sage: # needs sage.libs.pari
            sage: f = x^16 + x^14 - x^10 + x^8 - x^6 + x^2 + 1
            sage: f.is_cyclotomic(algorithm='pari')
            False
            sage: f.is_cyclotomic(algorithm='sage')
            False
            sage: g = x^16 + x^14 - x^10 - x^8 - x^6 + x^2 + 1
            sage: g.is_cyclotomic(algorithm='pari')
            True
            sage: g.is_cyclotomic(algorithm='sage')
            True

            sage: y = polygen(QQ)
            sage: (y/2 - 1/2).is_cyclotomic()
            False
            sage: (2*(y/2 - 1/2)).is_cyclotomic()                                       # needs sage.libs.pari
            True

        Invalid arguments::

            sage: (x - 3).is_cyclotomic(algorithm='sage', certificate=True)             # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: no implementation of the certificate within Sage

        Test using other rings::

            sage: z = polygen(GF(5))
            sage: (z - 1).is_cyclotomic()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented in nonzero characteristic

        TESTS::

            sage: R = ZZ['x']
            sage: for _ in range(20):                                                   # needs sage.libs.pari
            ....:     p = R.random_element(degree=randint(10,20))
            ....:     ans_pari = p.is_cyclotomic(algorithm='pari')
            ....:     ans_sage = p.is_cyclotomic(algorithm='sage')
            ....:     assert ans_pari == ans_sage, "problem with p={}".format(p)
            sage: for d in range(2, 20):                                                # needs sage.libs.pari
            ....:     p = cyclotomic_polynomial(d)
            ....:     assert p.is_cyclotomic(algorithm='pari'), "pari problem with p={}".format(p)
            ....:     assert p.is_cyclotomic(algorithm='sage'), "sage problem with p={}".format(p)

        Test the output type when ``certificate=True``::

            sage: type((x^2 - 2).is_cyclotomic(certificate=True))                       # needs sage.libs.pari
            <class 'sage.rings.integer.Integer'>
            sage: type((x - 1).is_cyclotomic(certificate=True))                         # needs sage.libs.pari
            <class 'sage.rings.integer.Integer'>

        Check that the arguments are forwarded when the input is not a
        polynomial with coefficients in `\ZZ`::

            sage: x = polygen(QQ)
            sage: (x - 1).is_cyclotomic(certificate=True)                               # needs sage.libs.pari
            1
        """
        S = self.base_ring()
        if S.characteristic() != 0:
            raise NotImplementedError("not implemented in nonzero characteristic")
        if S != ZZ:
            try:
                f = self.change_ring(ZZ)
            except TypeError:
                return False
            return f.is_cyclotomic(certificate=certificate, algorithm=algorithm)

        if algorithm == "pari":
            ans = self.__pari__().poliscyclo()
            return Integer(ans) if certificate else bool(ans)

        elif algorithm != "sage":
            raise ValueError("algorithm must be either 'pari' or 'sage'")

        elif certificate:
            raise ValueError("no implementation of the certificate within Sage")

        if not self.is_monic():
            return False

        P = self._parent
        gen = P.gen()

        if self == gen - 1:  # the first cyc. pol. is treated apart
            return True

        if self.constant_coefficient() != 1:
            return False

        if not self.is_irreducible():
            return False

        coefs = self.coefficients(sparse=False)

        # compute the graeffe transform of self
        po_odd = P(coefs[1::2])
        po_even = P(coefs[0::2])
        f1 = po_even * po_even - gen * (po_odd * po_odd)

        # first case
        if f1 == self:
            return True

        # second case
        selfminus = self(-gen)
        if f1 == selfminus:
            if selfminus.leading_coefficient() < 0 and (-selfminus).is_cyclotomic(algorithm='sage'):
                return True
            elif selfminus.is_cyclotomic(algorithm='sage'):
                return True

        # third case, we need to take a square root
        ans, ff1 = f1.is_square(True)
        return ans and ff1.is_cyclotomic(algorithm='sage')

    def is_cyclotomic_product(self):
        r"""
        Test whether this polynomial is a product of cyclotomic polynomials.

        This method simply calls the function :pari:`poliscycloprod`
        from the Pari library.

        .. SEEALSO::

            :meth:`is_cyclotomic`
            :meth:`cyclotomic_part`
            :meth:`has_cyclotomic_factor`

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: (x^5 - 1).is_cyclotomic_product()                                     # needs sage.libs.pari
            True
            sage: (x^5 + x^4 - x^2 + 1).is_cyclotomic_product()                         # needs sage.libs.pari
            False

            sage: p = prod(cyclotomic_polynomial(i) for i in [2, 5, 7, 12])
            sage: p.is_cyclotomic_product()                                             # needs sage.libs.pari
            True

            sage: (x^5 - 1/3).is_cyclotomic_product()
            False

            sage: x = polygen(Zmod(5))
            sage: (x - 1).is_cyclotomic_product()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented in nonzero characteristic
        """
        if self.base_ring().characteristic() != 0:
            raise NotImplementedError("not implemented in nonzero characteristic")
        if self.base_ring() != ZZ:
            try:
                f = self.change_ring(ZZ)
            except TypeError:
                return False
            return f.is_cyclotomic_product()

        return bool(self.__pari__().poliscycloprod())

    def cyclotomic_part(self):
        r"""
        Return the product of the irreducible factors of this polynomial
        which are cyclotomic polynomials.

        The algorithm assumes that the polynomial has rational coefficients.

        .. SEEALSO::

            :meth:`is_cyclotomic`
            :meth:`is_cyclotomic_product`
            :meth:`has_cyclotomic_factor`

        EXAMPLES::

            sage: P.<x> = PolynomialRing(Integers())
            sage: pol = 2*(x^4 + 1)
            sage: pol.cyclotomic_part()
            x^4 + 1
            sage: pol = x^4 + 2
            sage: pol.cyclotomic_part()
            1
            sage: pol = (x^4 + 1)^2 * (x^4 + 2)
            sage: pol.cyclotomic_part()
            x^8 + 2*x^4 + 1

            sage: P.<x> = PolynomialRing(QQ)
            sage: pol = (x^4 + 1)^2 * (x^4 + 2)
            sage: pol.cyclotomic_part()
            x^8 + 2*x^4 + 1

            sage: pol = (x - 1) * x * (x + 2)
            sage: pol.cyclotomic_part()
            x - 1

        TESTS::

            sage: P.<x> = PolynomialRing(RR)
            sage: pol = (x^4 + 1)^2 * (x^4 + 2)
            sage: pol.cyclotomic_part()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented for inexact base rings

            sage: x = polygen(Zmod(5))
            sage: (x-1).cyclotomic_part()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented in nonzero characteristic
        """
        S = self.base_ring()
        if S.characteristic():
            raise NotImplementedError("not implemented in nonzero characteristic")
        if not S.is_exact():
            raise NotImplementedError("not implemented for inexact base rings")
        R = self._parent
        x = R.gen()
        # Extract Phi_n when n is odd.
        t1 = self
        while True:
            t2 = t1.gcd(t1(x**2))
            if t1.degree() == t2.degree():
                break
            t1 = t2
        ans = t1
        ans //= x**ans.valuation()
        # Extract Phi_n when v_2(n) = 1, 2, ...
        t0 = self // t1
        i = 0
        while t0.degree():
            t1 = t0
            while True:
                t2 = t1.gcd(t1(-x**2))
                if t1.degree() == t2.degree():
                    break
                t1 = t2
            ans *= t1(x**(2**i))
            t0 = t0 // t1
            t1 = t0.gcd(t0(-x))
            t0 = R(list(t1)[::2])
            i += 1
        return ans // ans.leading_coefficient()

    def has_cyclotomic_factor(self):
        r"""
        Return ``True`` if the given polynomial has a nontrivial cyclotomic factor.

        The algorithm assumes that the polynomial has rational coefficients.

        If the polynomial is known to be irreducible, it may be slightly more
        efficient to call :meth:`is_cyclotomic` instead.

        .. SEEALSO::

            :meth:`is_cyclotomic`
            :meth:`is_cyclotomic_product`
            :meth:`cyclotomic_part`

        EXAMPLES::

            sage: pol.<x> = PolynomialRing(Rationals())
            sage: u = x^5 - 1; u.has_cyclotomic_factor()
            True
            sage: u = x^5 - 2; u.has_cyclotomic_factor()
            False
            sage: u = pol(cyclotomic_polynomial(7)) * pol.random_element()  # random
            sage: u.has_cyclotomic_factor()                                 # random
            True
        """
        if not QQ.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError("coefficients not rational")
        polRing = self.parent()
        x = polRing.gen()

        pol1 = self
        # First, while pol1 has a nontrivial even factor, replace
        # that factor with the polynomials whose roots are the squares of
        # the roots of that factor. This replaces any roots of unity of order
        # divisible by 4 with roots of unity of order not divisible by 4.

        pol2 = pol1.gcd(pol1(-x))
        while not pol2.is_constant():
            pol1 = (pol1 // pol2) * polRing(pol2.list()[::2])
            pol2 = pol1.gcd(pol1(-x))

        # Next, replace pol1 with the polynomial whose roots are the
        # squares of pol1. This replaces any roots of unity of even order
        # with roots of unity of odd order.
        pol1 = polRing((pol1*pol1(-x)).list()[::2])

        # Finally, find the largest factor of pol1 whose roots are
        # stable under squaring. This factor is constant if and only if
        # the original polynomial has no cyclotomic factor.
        while True:
            if pol1.is_constant():
                return False
            pol2 = pol1.gcd(polRing((pol1*pol1(-x)).list()[::2]))
            if pol1.degree() == pol2.degree():
                return True
            pol1 = pol2

    def homogenize(self, var='h'):
        r"""
        Return the homogenization of this polynomial.

        The polynomial itself is returned if it is homogeneous already. Otherwise,
        its monomials are multiplied with the smallest powers of ``var`` such
        that they all have the same total degree.

        INPUT:

        - ``var`` -- a variable in the polynomial ring (as a string, an element
          of the ring, or ``0``) or a name for a new variable (default:
          ``'h'``)

        OUTPUT:

        If ``var`` specifies the variable in the polynomial ring, then a
        homogeneous element in that ring is returned. Otherwise, a homogeneous
        element is returned in a polynomial ring with an extra last variable
        ``var``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^2 + 1
            sage: f.homogenize()
            x^2 + h^2

        The parameter ``var`` can be used to specify the name of the variable::

            sage: g = f.homogenize('z'); g
            x^2 + z^2
            sage: g.parent()
            Multivariate Polynomial Ring in x, z over Rational Field

        However, if the polynomial is homogeneous already, then that parameter
        is ignored and no extra variable is added to the polynomial ring::

            sage: f = x^2
            sage: g = f.homogenize('z'); g
            x^2
            sage: g.parent()
            Univariate Polynomial Ring in x over Rational Field

        For compatibility with the multivariate case, if ``var`` specifies the
        variable of the polynomial ring, then the monomials are multiplied with
        the smallest powers of ``var`` such that the result is homogeneous; in
        other words, we end up with a monomial whose leading coefficient is the
        sum of the coefficients of the polynomial::

            sage: f = x^2 + x + 1
            sage: f.homogenize('x')
            3*x^2

        In positive characteristic, the degree can drop in this case::

            sage: R.<x> = GF(2)[]
            sage: f = x + 1
            sage: f.homogenize(x)
            0

        For compatibility with the multivariate case, the parameter ``var`` can
        also be 0 to specify the variable in the polynomial ring::

            sage: R.<x> = QQ[]
            sage: f = x^2 + x + 1
            sage: f.homogenize(0)
            3*x^2
        """
        if self.is_homogeneous():
            return self

        x, = self.variables()

        if isinstance(var, int) or isinstance(var, Integer):
            if var:
                raise TypeError("Variable index %d must be < 1." % var)
            else:
                return sum(self.coefficients())*x**self.degree()

        x_name = self.variable_name()
        var = str(var)

        if var == x_name:
            return sum(self.coefficients())*x**self.degree()

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        P = PolynomialRing(self.base_ring(), [x_name, var])
        return P(self)._homogenize(1)

    def is_homogeneous(self):
        r"""
        Return ``True`` if this polynomial is homogeneous.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: x.is_homogeneous()
            True
            sage: P(0).is_homogeneous()
            True
            sage: (x + 1).is_homogeneous()
            False
        """
        return len(self.exponents()) < 2

    def nth_root(self, n):
        r"""
        Return a `n`-th root of this polynomial.

        This is computed using Newton method in the ring of power
        series. This method works only when the base ring is an
        integral domain. Moreover, for polynomial whose coefficient of
        lower degree is different from 1, the elements of the base
        ring should have a method :meth:`nth_root` implemented.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: a = 27 * (x+3)**6 * (x+5)**3
            sage: a.nth_root(3)
            3*x^3 + 33*x^2 + 117*x + 135

            sage: b = 25 * (x^2 + x + 1)
            sage: b.nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power
            sage: R(0).nth_root(3)
            0
            sage: R.<x> = QQ[]
            sage: a = 1/4 * (x/7 + 3/2)^2 * (x/2 + 5/3)^4
            sage: a.nth_root(2)
            1/56*x^3 + 103/336*x^2 + 365/252*x + 25/12

            sage: # needs sage.rings.number_field
            sage: K.<sqrt2> = QuadraticField(2)
            sage: R.<x> = K[]
            sage: a = (x + sqrt2)^3 * ((1+sqrt2)*x - 1/sqrt2)^6
            sage: b = a.nth_root(3); b
            (2*sqrt2 + 3)*x^3 + (2*sqrt2 + 2)*x^2 + (-2*sqrt2 - 3/2)*x + 1/2*sqrt2
            sage: b^3 == a
            True

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQbar[]
            sage: p = x**3 + QQbar(2).sqrt() * x - QQbar(3).sqrt()
            sage: r = (p**5).nth_root(5)
            sage: r * p[0] == p * r[0]
            True
            sage: p = (x+1)^20 + x^20
            sage: p.nth_root(20)
            Traceback (most recent call last):
            ...
            ValueError: not a 20th power

            sage: # needs sage.rings.finite_rings
            sage: z = GF(4).gen()
            sage: R.<x> = GF(4)[]
            sage: p = z*x**4 + 2*x - 1
            sage: r = (p**15).nth_root(15)
            sage: r * p[0] == p * r[0]
            True
            sage: ((x+1)**2).nth_root(2)
            x + 1
            sage: ((x+1)**4).nth_root(4)
            x + 1
            sage: ((x+1)**12).nth_root(12)
            x + 1
            sage: (x^4 + x^3 + 1).nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power
            sage: p = (x+1)^17 + x^17
            sage: r = p.nth_root(17)
            Traceback (most recent call last):
            ...
            ValueError: not a 17th power

            sage: R1.<x> = QQ[]
            sage: R2.<y> = R1[]
            sage: R3.<z> = R2[]
            sage: (((y**2+x)*z^2 + x*y*z + 2*x)**3).nth_root(3)
            (y^2 + x)*z^2 + x*y*z + 2*x
            sage: ((x+y+z)**5).nth_root(5)
            z + y + x

        Here we consider a base ring without ``nth_root`` method. The third
        example with a non-trivial coefficient of lowest degree raises an error::

            sage: # needs sage.libs.pari
            sage: R.<x> = QQ[]
            sage: R2 = R.quotient(x**2 + 1)
            sage: x = R2.gen()
            sage: R3.<y> = R2[]
            sage: (y**2 - 2*y + 1).nth_root(2)
            -y + 1
            sage: (y**3).nth_root(3)
            y
            sage: (y**2 + x).nth_root(2)
            Traceback (most recent call last):
            ...
            AttributeError: ... has no attribute 'nth_root'...

        TESTS::

            sage: R.<x> = ZZ[]
            sage: (x^12).nth_root(6)
            x^2
            sage: ((3*x)^15).nth_root(5)
            27*x^3
            sage: parent(R.one().nth_root(3))
            Univariate Polynomial Ring in x over Integer Ring
            sage: p = (x+1)**20 + x^20
            sage: p.nth_root(20)
            Traceback (most recent call last):
            ...
            ValueError: not a 20th power
            sage: (x^3 - 1).nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power
            sage: (x^3 - 1).nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power

            sage: Zmod(4)['x'].one().nth_root(4)
            Traceback (most recent call last):
            ...
            ValueError: n-th root of polynomials over rings with zero divisors
            not implemented

        Some random tests::

            sage: for R in [QQ['x'], GF(4)['x']]:                                       # needs sage.modules sage.rings.finite_rings
            ....:     for _ in range(30):
            ....:         p = R.random_element(degree=randint(10,20))
            ....:         n = ZZ.random_element(2,20)
            ....:         r = (p**n).nth_root(n)
            ....:         assert r.parent() is R, "R={}\nn={}\np={}".format(R,n,p)
            ....:         pl = p.leading_coefficient()
            ....:         rl = r.leading_coefficient()
            ....:         assert p == r * pl/rl, "R={}\np={}\nr={}".format(R,p,r)
        """
        R = self.base_ring()
        S = self.parent()
        if R not in sage.categories.integral_domains.IntegralDomains():
            raise ValueError("n-th root of polynomials over rings with zero divisors not implemented")

        if n <= 0:
            raise ValueError("n (={}) must be positive".format(n))
        elif n == 1 or self.is_zero() or self.is_one():
            return self
        elif self.degree() % n:
            raise ValueError("not a %s power"%Integer(n).ordinal_str())
        elif self.get_unsafe(0).is_zero(): # We know that self is not 0, so it must have degree >= 0
            # p = x^k q
            # p^(1/n) = x^(k/n) q^(1/n)
            i = self.valuation()
            if i%n:
                raise ValueError("not a %s power"%Integer(n).ordinal_str())
            return (self >> i).nth_root(n) << (i // n)

        if self.get_unsafe(0).is_one():
            start = S.one()
        else:
            start = S(self.get_unsafe(0).nth_root(n))

        cdef Polynomial p, q
        p = self.change_ring(R.fraction_field())
        q = p._nth_root_series(n, self.degree() // n + 1, start)

        # (possible) TODO: below we check that the result is the
        # n-th root. But in ZZ[x] that can be anticipated: we can
        # detect that a given polynomial is not a n-th root inside
        # the iteration of Newton method by looking at denominators.
        if q**n == p:
            return S(q)
        else:
            raise ValueError("not a %s power"%Integer(n).ordinal_str())

    def _nth_root_series(self, long n, long prec, start=None):
        r"""
        Return the first ``prec`` coefficients of the ``n``-th root series of this polynomial.

        The method might fail if the exponent ``n`` or the coefficient of
        lowest degree is not invertible in the base ring. In both cases an
        :exc:`ArithmeticError` is raised.

        INPUT:

        - ``n`` -- positive integer; the exponent of the root

        - ``prec`` -- positive integer; the precision of the result

        - ``start`` -- (optional) the first term of the result. This
          is only considered when the valuation is zero, i.e. when the
          polynomial has a nonzero constant term.

        ALGORITHM:

        Let us denote by `a` the polynomial from which we wish to extract
        a `n`-th root. The algorithm uses the Newton method for the fixed
        point of `F(x) = x^{-n} - a^{-1}`. The advantage of this approach
        compared to the more naive `x^n - a` is that it does require only
        one polynomial inversion instead of one per iteration of the Newton
        method.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (1 + x)._nth_root_series(2, 5)
            -5/128*x^4 + 1/16*x^3 - 1/8*x^2 + 1/2*x + 1
            sage: R.zero()._nth_root_series(3, 5)
            0
            sage: R.one()._nth_root_series(3, 5)
            1

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQbar[]
            sage: p = 2 + 3*x^2
            sage: q = p._nth_root_series(3, 20)
            sage: (q**3).truncate(20)
            3*x^2 + 2

        The exponent must be invertible in the base ring::

            sage: R.<x> = ZZ[x]
            sage: (1 + x)._nth_root_series(2, 5)
            Traceback (most recent call last):
            ...
            ArithmeticError: exponent not invertible in base ring

        Though, the base ring needs not be a field::

            sage: Ru.<u> = QQ[]
            sage: Rux.<x> = Ru[]
            sage: (4 + u*x)._nth_root_series(2,5)
            -5/16384*u^4*x^4 + 1/512*u^3*x^3 - 1/64*u^2*x^2 + 1/4*u*x + 2
            sage: ((4 + u*x)._nth_root_series(2,5)**2).truncate(5)
            u*x + 4

        Finite characteristic::

            sage: R.<x> = GF(2)[]
            sage: (1 + x)._nth_root_series(3, 10)
            x^9 + x^8 + x^3 + x^2 + x + 1
            sage: (1 + x^2)._nth_root_series(2, 10)
            x + 1
            sage: (1 + x)._nth_root_series(2, 10)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power

        TESTS::

            sage: QQ['x'].zero()._nth_root_series(3, 5).parent()
            Univariate Polynomial Ring in x over Rational Field
            sage: QQ['x'].one()._nth_root_series(3, 5).parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        cdef Integer c, cc, e, m, mp1
        cdef Polynomial p, q

        R = self.base_ring()
        S = self.parent()

        m = ZZ.coerce(n)
        if m <= 0:
            raise ValueError("n (={}) must be positive".format(m))
        elif m.is_one() or self.is_zero() or self.is_one():
            return self
        elif self.get_unsafe(0).is_zero(): # we know that self is not zero, so the degree >= 0
            # p = x^i q
            # p^(1/m) = x^(i/m) q^(1/m)
            i = self.valuation()
            if i % m:
                raise ValueError("not a %s power"%m.ordinal_str())
            return (self >> i)._nth_root_series(m, prec - i // m) << (i // m)
        else:
            c = R.characteristic()
            if c and not n % c:
                # characteristic divides n
                e = m.valuation(c)
                cc = c**e
                ans = {}
                for i in range(self.degree()+1):
                    if self.get_unsafe(i):
                        if i % cc:
                            raise ValueError("not a %s power"%m.ordinal_str())
                        ans[i//cc] = self.get_unsafe(i).nth_root(cc)
                p = self._parent(ans)
                m = m // cc
                if m.is_one():
                    return p
            else:
                p = self

            # beginning of Newton method
            # (we can safely assume that the valuation is 0)
            S = p.parent()

            if start is not None:
                a = R(start)
            elif p[0].is_one():
                a = R.one()
            else:
                a = p[0].nth_root(m)

            try:
                q = S(a.inverse_of_unit())
            except ArithmeticError:
                raise ArithmeticError("constant coefficient not invertible in base ring")

            try:
                mi = R(m).inverse_of_unit()
            except ArithmeticError:
                raise ArithmeticError("exponent not invertible in base ring")

            from sage.misc.misc import newton_method_sizes
            mp1 = m + 1
            for i in newton_method_sizes(prec):
                q = mi * (mp1 * q - p._mul_trunc_(q._power_trunc(mp1, i), i))
            return q.inverse_series_trunc(prec)

    @coerce_binop
    def divides(self, p):
        r"""
        Return ``True`` if this polynomial divides `p`.

        This method is only implemented for polynomials over an integral domain.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (2*x + 1).divides(4*x**2 - 1)
            True
            sage: (2*x + 1).divides(4*x**2 + 1)
            False
            sage: (2*x + 1).divides(R(0))
            True
            sage: R(0).divides(2*x + 1)
            False
            sage: R(0).divides(R(0))
            True
            sage: S.<y> = R[]
            sage: p = x * y**2 + (2*x + 1) * y + x + 1
            sage: q = (x + 1) * y + (3*x + 2)
            sage: q.divides(p)
            False
            sage: q.divides(p * q)
            True
            sage: R.<x> = Zmod(6)[]
            sage: p = 4*x + 3
            sage: q = 5*x**2 + x + 2
            sage: q.divides(p)
            False
            sage: p.divides(q)
            False

        TESTS::

            sage: R.<x> = PolynomialRing(ZZ, implementation="NTL")                      # needs sage.libs.ntl
            sage: (2*x + 1).divides(4*x**2 + 1)                                         # needs sage.libs.ntl
            False

            sage: # needs sage.rings.finite_rings
            sage: K.<z> = GF(4)
            sage: R.<x> = K[]
            sage: S.<y> = R[]
            sage: p = ((3*z + 2)*x + 2*z - 1) * y + 2*x + z
            sage: q = y^2 + z*y*x + 2*y + z
            sage: p.divides(q), p.divides(p*q)
            (False, True)

            sage: R.<x,y> = GF(2)[]
            sage: S.<z> = R[]
            sage: p = (x+y+1) * z + x*y
            sage: q = (y^2-x^2) * z^2 + z + x-y
            sage: p.divides(q), p.divides(p*q)                                          # needs sage.libs.singular
            (False, True)
            sage: R.<x> = Zmod(6)[]
            sage: p = 4*x + 3
            sage: q = 2*x**2 + x + 2
            sage: p.divides(q)
            Traceback (most recent call last):
            ...
            NotImplementedError: divisibility test only implemented for polynomials over an integral domain unless obvious non divisibility of leading terms
        """
        if p.is_zero():
            return True          # everything divides 0
        if self.is_zero():
            return False      # 0 only divides 0
        try:
            if self.is_unit():
                return True   # units divide everything
        except NotImplementedError:
            if self.is_one():
                return True    # if is_unit is not implemented

        if self.degree() > p.degree():
            return False

        if not self.leading_coefficient().divides(p.leading_coefficient()):
            return False

        if not self.base_ring().is_integral_domain():
            raise NotImplementedError("divisibility test only implemented for polynomials over an integral domain unless obvious non divisibility of leading terms")

        try:
            return (p % self).is_zero()      # if quo_rem is defined
        except ArithmeticError:
            return False

    def specialization(self, D=None, phi=None):
        r"""
        Specialization of this polynomial.

        Given a family of polynomials defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional)

        - ``phi`` -- :class:`SpecializationMorphism` (optional)

        OUTPUT: a new polynomial

        EXAMPLES::

            sage: R.<c> = PolynomialRing(ZZ)
            sage: S.<z> = PolynomialRing(R)
            sage: F = c*z^2 + c^2
            sage: F.specialization({c:2})
            2*z^2 + 4

        ::

            sage: A.<c> = QQ[]
            sage: R.<x> = Frac(A)[]
            sage: X = (1 + x/c).specialization({c:20})
            sage: X
            1/20*x + 1
            sage: X.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if D is None:
            if phi is None:
                raise ValueError("either the dictionary or the specialization must be provided")
        else:
            from sage.rings.polynomial.flatten import SpecializationMorphism
            phi = SpecializationMorphism(self._parent,D)
        return phi(self)

    def _log_series(self, long n):
        r"""
        Return the power series expansion of logarithm of this polynomial,
        truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]                                                       # needs sage.libs.flint
            sage: (1 + x)._log_series(3)                                                # needs sage.libs.flint
            -0.5000000000000000*x^2 + x
        """
        raise NotImplementedError

    def _exp_series(self, long n):
        r"""
        Return the power series expansion of exponential of this polynomial,
        truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]                                                       # needs sage.libs.flint
            sage: x._exp_series(3)                                                      # needs sage.libs.flint
            0.5000000000000000*x^2 + x + 1.000000000000000
        """
        raise NotImplementedError

    def _atan_series(self, long n):
        r"""
        Return the power series expansion of arctangent of this polynomial,
        truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._atan_series(4)
            -1/3*x^3 + x
        """
        raise NotImplementedError

    def _atanh_series(self, long n):
        r"""
        Return the power series expansion of hyperbolic arctangent of this
        polynomial, truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._atanh_series(4)
            1/3*x^3 + x
        """
        raise NotImplementedError

    def _asin_series(self, long n):
        r"""
        Return the power series expansion of arcsine of this polynomial,
        truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._asin_series(4)
            1/6*x^3 + x
        """
        raise NotImplementedError

    def _asinh_series(self, long n):
        r"""
        Return the power series expansion of hyperbolic arcsine of this
        polynomial, truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._asinh_series(4)
            -1/6*x^3 + x
        """
        raise NotImplementedError

    def _tan_series(self, long n):
        r"""
        Return the power series expansion of tangent of this polynomial,
        truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._tan_series(4)
            1/3*x^3 + x
        """
        raise NotImplementedError

    def _sin_series(self, long n):
        r"""
        Return the power series expansion of sine of this polynomial, truncated
        to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._sin_series(4)
            -1/6*x^3 + x
        """
        raise NotImplementedError

    def _cos_series(self, long n):
        r"""
        Return the power series expansion of cosine of this polynomial,
        truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._cos_series(4)
            -1/2*x^2 + 1
        """
        raise NotImplementedError

    def _sinh_series(self, long n):
        r"""
        Return the power series expansion of hyperbolic sine of this
        polynomial, truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._sinh_series(4)
            1/6*x^3 + x
        """
        raise NotImplementedError

    def _cosh_series(self, long n):
        r"""
        Return the power series expansion of hyperbolic cosine of this
        polynomial, truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._cosh_series(4)
            1/2*x^2 + 1
        """
        raise NotImplementedError

    def _tanh_series(self, long n):
        r"""
        Return the power series expansion of hyperbolic tangent of this
        polynomial, truncated to `O(x^n)`.

        EXAMPLES::

            sage: Pol.<x> = QQ[]
            sage: x._tanh_series(4)
            -1/3*x^3 + x
        """
        raise NotImplementedError


# ----------------- inner functions -------------
# Cython can't handle function definitions inside other function

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck(False)
cdef list do_schoolbook_product(list x, list y, Py_ssize_t deg):
    """
    Compute the truncated multiplication of two polynomials represented by
    lists, using the schoolbook algorithm.

    This is the core of _mul_generic and the code that is used by
    _mul_karatsuba below a threshold.

    INPUT:

    - ``x``, ``y`` -- lists of coefficients
    - ``deg`` -- degree at which the output should be truncated,
      negative values mean not to truncate at all

    TESTS:

    Doctested indirectly in _mul_generic and _mul_karatsuba. For the doctest we
    use a ring such that default multiplication calls external libraries::

        sage: K = ZZ['x']
        sage: f = K.random_element(8)
        sage: g = K.random_element(8)
        sage: f*g - f._mul_generic(g)
        0
    """
    cdef Py_ssize_t i, k, start, end
    cdef Py_ssize_t d1 = len(x)-1, d2 = len(y)-1
    if deg < 0 or deg > d1 + d2 + 1:
        deg = d1 + d2 + 1
    if d1 == -1:
        return x
    elif d2 == -1:
        return y
    elif d1 == 0:
        c = x[0]
        return [c*a for a in y[:deg]] # beware of noncommutative rings
    elif d2 == 0:
        c = y[0]
        return [a*c for a in x[:deg]] # beware of noncommutative rings
    coeffs = [None]*deg
    for k in range(deg):
        start = 0 if k <= d2 else k-d2  # max(0, k-d2)
        end = k if k <= d1 else d1    # min(k, d1)
        sum = x[start] * y[k-start]
        for i from start < i <= end:
            sum = sum + x[i] * y[k-i]
        coeffs[k] = sum
    return coeffs


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck(False)
cdef list do_karatsuba_different_size(list left, list right, Py_ssize_t K_threshold):
    """
    Multiply two polynomials of different degrees by splitting the one of
    largest degree in chunks that are multiplied with the other using the
    Karatsuba algorithm, as explained in _mul_karatsuba.

    INPUT:

    - ``left`` -- list representing a polynomial
    - ``right`` -- list representing a polynomial
    - ``K_threshold`` -- an Integer, a threshold to pass to the classical
      quadratic algorithm. During Karatsuba recursion, if one of the lists
      has length <= K_threshold the classical product is used instead.

    TESTS:

    This method is indirectly doctested in _mul_karatsuba.

    Here, we use Fibonacci numbers that need deepest recursion in this method.

        sage: K = ZZ['x']
        sage: f = K.random_element(21)
        sage: g = K.random_element(34)
        sage: f*g - f._mul_karatsuba(g,0)
        0
    """
    cdef Py_ssize_t n = len(left), m = len(right)
    cdef Py_ssize_t r, q, i, j, mi
    if n == 0 or m == 0:
        return []
    if n == 1:
        c = left[0]
        return [c*a for a in right]
    if m == 1:
        c = right[0]
        return [a*c for a in left] # beware of noncommutative rings
    if n <= K_threshold or m <= K_threshold:
        return do_schoolbook_product(left, right, -1)
    if n == m:
        return do_karatsuba(left, right, K_threshold, 0, 0, n)
    if n > m:
        # left is the bigger list
        # n is the bigger number
        q = n // m
        r = n % m
        output = do_karatsuba(left, right, K_threshold, 0, 0, m)
        for i from 1 <= i < q:
            mi = m*i
            carry = do_karatsuba(left, right, K_threshold, mi, 0, m)
            for j from 0 <= j < m-1:
                output[mi+j] = output[mi+j] + carry[j]
            output.extend(carry[m-1:])
        if r:
            mi = m*q
            carry = do_karatsuba_different_size(left[mi:], right, K_threshold)
            for j from 0 <= j < m-1:
                output[mi+j] = output[mi+j] + carry[j]
            output.extend(carry[m-1:])
        return output
    else:
        # n < m, I need to repeat the code due to the case
        # of noncommutative rings.
        q = m // n
        r = m % n
        output = do_karatsuba(left, right, K_threshold, 0, 0, n)
        for i from 1 <= i < q:
            mi = n*i
            carry = do_karatsuba(left, right, K_threshold, 0, mi, n)
            for j from 0 <= j < n-1:
                output[mi+j] = output[mi+j] + carry[j]
            output.extend(carry[n-1:])
        if r:
            mi = n*q
            carry = do_karatsuba_different_size(left, right[mi:], K_threshold)
            for j from 0 <= j < n-1:
                output[mi+j] = output[mi+j] + carry[j]
            output.extend(carry[n-1:])
        return output


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck(False)
cdef list do_karatsuba(list left, list right, Py_ssize_t K_threshold,Py_ssize_t start_l, Py_ssize_t start_r,Py_ssize_t num_elts):
    """
    Core routine for Karatsuba multiplication. This function works for two
    polynomials of the same degree.

    INPUT:

    - ``left`` -- list containing a slice representing a polynomial
    - ``right`` -- list containing the slice representing a polynomial with the
      same length as ``left``
    - ``K_threshold`` -- integer; for lists of ``length <= K_threshold``, the
      quadratic polynomial multiplication is used
    - ``start_l`` -- the index of ``left`` where the actual polynomial starts
    - ``start_r`` -- the index of ``right`` where the actual polynomial starts
    - ``num_elts`` -- the length of the polynomials

    Thus, the actual polynomials we want to multiply are represented by the
    slices: ``left[ start_l: start_l+num_elts ]``, ``right[ right_l: right_l+num_elts ]``.
    We use this representation in order to avoid creating slices of lists and
    create smaller lists.

    OUTPUT: list representing the product of the polynomials

    Doctested indirectly in _mul_karatsuba

    TESTS::

        sage: K.<x> = ZZ[]
        sage: f = K.random_element(50) + x^51
        sage: g = K.random_element(50) + x^51
        sage: f*g - f._mul_karatsuba(g,0)
        0

    Notes on the local variables:

    - ac will always be a list of length lenac
    - bd will always be a list of length lenbd
    - a_m_b and c_m_d are lists of length ne, we only make necessary additions
    - tt1 has length lenac
    """
    cdef Py_ssize_t e, ne, lenac, lenbd, start_le, start_re, i
    if num_elts == 0:
        return []
    if num_elts == 1:
        return [left[start_l]*right[start_r]]
    if num_elts <= K_threshold:
        # Special case of degree 2, no loop, no function call
        if num_elts == 2:
            b = left[start_l]
            a = left[start_l+1]
            d = right[start_r]
            c = right[start_r+1]
            return [b*d, a*d+b*c, a*c]
        return do_schoolbook_product(left[start_l:start_l+num_elts],
                right[start_r:start_r+num_elts], -1)
    if num_elts == 2:
        # beware of noncommutative rings
        b = left[start_l]
        a = left[start_l+1]
        d = right[start_r]
        c = right[start_r+1]
        ac = a*c
        bd = b*d
        return [bd, (a+b)*(c+d)-ac-bd, ac]
    e = num_elts//2
    ne = num_elts-e
    lenac = 2*ne-1
    lenbd = 2*e-1
    start_le = start_l+e
    start_re = start_r+e
    ac = do_karatsuba(left, right, K_threshold, start_le, start_re, ne)
    bd = do_karatsuba(left, right, K_threshold, start_l,  start_r,  e)
    a_m_b = left[start_le:start_le+ne]
    c_m_d = right[start_re:start_re+ne]
    for i from 0 <= i < e:
        a_m_b[i] = a_m_b[i] + left[start_l+i]
        c_m_d[i] = c_m_d[i] + right[start_r+i]
    tt1 = do_karatsuba(a_m_b, c_m_d, K_threshold, 0, 0, ne)
    # bd might be shorter than ac, we divide the operations in two loops
    for i from 0 <= i < lenbd:
        tt1[i] = tt1[i] - (ac[i]+bd[i])
    for i from lenbd <= i < lenac:
        tt1[i] = tt1[i] - ac[i]
    # Reconstruct the product from the lists bd, tt1, ac.
    for i from 0 <= i < e-1:
        bd[e+i] = bd[e+i] + tt1[i]
    bd.append(tt1[e-1])
    for i from 0 <= i < lenac -e:
        ac[i] = ac[i] + tt1[e+i]
    return bd + ac


cdef class Polynomial_generic_dense(Polynomial):
    """
    A generic dense polynomial.

    EXAMPLES::

        sage: f = QQ['x']['y'].random_element()
        sage: loads(f.dumps()) == f
        True

    TESTS::

        sage: # needs sage.rings.real_mpfr
        sage: from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense
        sage: isinstance(f, Polynomial_generic_dense)
        True
        sage: f = CC['x'].random_element()
        sage: isinstance(f, Polynomial_generic_dense)
        True

        sage: R.<x> = QQ[]
        sage: S = R['y']
        sage: S((x^2, 2, 1 + x))
        (x + 1)*y^2 + 2*y + x^2
    """
    def __init__(self, parent, x=None, int check=1, is_gen=False, int construct=0, **kwds):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self._coeffs = []
            return

        R = parent.base_ring()
        if isinstance(x, (list, tuple)):
            if check:
                self._coeffs = [R(t) for t in x]
                self._normalize()
            else:
                self._coeffs = x
            return

        if isinstance(x, sage.rings.fraction_field_element.FractionFieldElement):
            if x.denominator() != 1:
                raise TypeError("denominator must be 1")
            else:
                x = x.numerator()

        if isinstance(x, Polynomial):
            if (<Element>x)._parent is self._parent:
                x = x.list(copy=True)
            elif R.has_coerce_map_from((<Element>x)._parent):# is R or (<Element>x)._parent == R:
                try:
                    if x.is_zero():
                        self._coeffs = []
                        return
                except (AttributeError, TypeError):
                    pass
                x = [x]
            else:
                self._coeffs = [R(a, **kwds) for a in x.list(copy=False)]
                if check:
                    self._normalize()
                return

        elif isinstance(x, int) and x == 0:
            self._coeffs = []
            return

        elif isinstance(x, dict):
            x = _dict_to_list(x, R.zero())

        elif isinstance(x, pari_gen):
            x = [R(w, **kwds) for w in x.list()]
            check = 0
        elif not isinstance(x, (list, tuple)):
            # We trust that the element constructors do not send x=0
#            if x:
            x = [x]   # constant polynomials
#            else:
#                x = []    # zero polynomial
        if check:
            self._coeffs = [R(z, **kwds) for z in x]
            self._normalize()
        else:
            self._coeffs = x

    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P):
        cdef type t = type(self)
        cdef Polynomial_generic_dense f = <Polynomial_generic_dense>t.__new__(t)
        f._parent = P
        f._coeffs = coeffs
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

        Make sure we're testing the right method::

            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        """
        return make_generic_polynomial, (self._parent, self._coeffs)

    def __bool__(self):
        return bool(self._coeffs)

    cpdef bint is_term(self) except -1:
        """
        Return ``True`` if this polynomial is a nonzero element of the
        base ring times a power of the variable.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R.<x> = SR[]
            sage: R(0).is_term()
            False
            sage: R(1).is_term()
            True
            sage: (3*x^5).is_term()
            True
            sage: (1 + 3*x^5).is_term()
            False
        """
        if not self._coeffs:
            return False

        for c in self._coeffs[:-1]:
            if c:
                return False
        return True

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Polynomial _mul_term(self, Polynomial term, bint term_on_right):
        """
        Return the product ``self * term``, where ``term`` is a polynomial
        with a single term.
        """
        cdef Py_ssize_t d = len( (<Polynomial_generic_dense> term)._coeffs ) - 1
        cdef Py_ssize_t i
        cdef list x = self._coeffs
        cdef Py_ssize_t ell = len(x)
        c = term.get_unsafe(d)
        cdef list v = [self.base_ring().zero()] * (d + ell)
        if term_on_right:
            for i in range(ell):
                v[i+d] = x[i] * c
        else:
            for i in range(ell):
                v[i+d] = c * x[i]
        cdef Polynomial_generic_dense res = self._new_c(v, self._parent)
        #if not v[len(v)-1]:
        # "normalize" checks this anyway...
        res._normalize()
        return res

    cdef int _normalize(self) except -1:
        """
        TESTS:

        Check that exceptions are propagated correctly (:issue:`18274`)::

            sage: class BrokenRational(Rational):
            ....:     def __bool__(self):
            ....:         raise NotImplementedError("cannot check whether number is nonzero")
            ....:
            sage: z = BrokenRational()
            sage: R.<x> = QQ[]
            sage: from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense
            sage: Polynomial_generic_dense(R, [z])
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot check whether number is nonzero
        """
        cdef list x = self._coeffs
        cdef Py_ssize_t n = len(x) - 1
        while n >= 0 and not x[n]:
            del x[n]
            n -= 1

    def __hash__(self):
        return self._hash_c()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef get_unsafe(self, Py_ssize_t n):
        """
        Return the `n`-th coefficient of ``self``.

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
        """
        return self._coeffs[n]

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
        if n >= 0 and n < len(self._coeffs):
            self._coeffs[n] = value
            if n == len(self._coeffs) and value == 0:
                self._normalize()
        elif n < 0:
            raise IndexError("polynomial coefficient index must be nonnegative")
        elif value != 0:
            zero = self.base_ring().zero()
            for _ in range(len(self._coeffs), n):
                self._coeffs.append(zero)
            self._coeffs.append(value)

    def __floordiv__(self, right):
        """
        Return the quotient upon division (no remainder).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQbar[]
            sage: f = (1+2*x)^3 + 3*x; f
            8*x^3 + 12*x^2 + 9*x + 1
            sage: g = f // (1+2*x); g
            4*x^2 + 4*x + 5/2
            sage: f - g * (1+2*x)
            -3/2
            sage: f.quo_rem(1+2*x)
            (4*x^2 + 4*x + 5/2, -3/2)

        TESTS:

        Check that :issue:`13048` and :issue:`2034` are fixed::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQbar[]
            sage: x // x
            1
            sage: x // 1
            x
            sage: x // int(1)
            x
            sage: x //= int(1); x
            x
            sage: int(1) // x  # check that this doesn't segfault
            Traceback (most recent call last):
            ...
            AttributeError: type object 'int' has no attribute 'base_ring'...
        """
        if have_same_parent(self, right):
            return (<Polynomial_generic_dense>self)._floordiv_(<Polynomial_generic_dense>right)
        P = parent(self)
        d = P.base_ring()(right)
        cdef Polynomial_generic_dense res = (<Polynomial_generic_dense>self)._new_c([c // d for c in (<Polynomial_generic_dense>self)._coeffs], P)
        res._normalize()
        return res

    cpdef _add_(self, right):
        r"""
        Add two polynomials.

        EXAMPLES::

            sage: R.<y> = QQ[]
            sage: S.<x> = R[]
            sage: S([0,1,y,2*y]) + S([1,-2*y,3])   # indirect doctest
            2*y*x^3 + (y + 3)*x^2 + (-2*y + 1)*x + 1
        """
        cdef Polynomial_generic_dense res
        cdef Py_ssize_t check=0, i, min
        x = (<Polynomial_generic_dense>self)._coeffs
        y = (<Polynomial_generic_dense>right)._coeffs
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
            res._normalize()
            return res
        else:
            return self._new_c(low + high, self._parent)

    cpdef _sub_(self, right):
        cdef Polynomial_generic_dense res
        cdef Py_ssize_t check=0, i, min
        x = (<Polynomial_generic_dense>self)._coeffs
        y = (<Polynomial_generic_dense>right)._coeffs
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
            res._normalize()
            return res
        else:
            return self._new_c(low + high, self._parent)

    cpdef _rmul_(self, Element c):
        if not self._coeffs:
            return self
        if c._parent is not (<Element>self._coeffs[0])._parent:
            c = (<Element>self._coeffs[0])._parent.coerce(c)
        v = [c * a for a in self._coeffs]
        cdef Polynomial_generic_dense res = self._new_c(v, self._parent)
        #if not v[len(v)-1]:
        # "normalize" checks this anyway...
        res._normalize()
        return res

    cpdef _lmul_(self, Element c):
        if not self._coeffs:
            return self
        if c._parent is not (<Element>self._coeffs[0])._parent:
            c = (<Element>self._coeffs[0])._parent.coerce(c)
        v = [a * c for a in self._coeffs]
        cdef Polynomial_generic_dense res = self._new_c(v, self._parent)
        #if not v[len(v)-1]:
        # "normalize" checks this anyway...
        res._normalize()
        return res

    cpdef constant_coefficient(self):
        """
        Return the constant coefficient of this polynomial.

        OUTPUT: element of base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: S.<x> = R[]
            sage: f = x*t + x + t
            sage: f.constant_coefficient()
            t
        """
        if not self._coeffs:
            return self.base_ring().zero()
        else:
            return self._coeffs[0]

    cpdef list list(self, bint copy=True):
        """
        Return a new copy of the list of the underlying elements of ``self``.

        EXAMPLES::

            sage: R.<x> = GF(17)[]
            sage: f = (1+2*x)^3 + 3*x; f
            8*x^3 + 12*x^2 + 9*x + 1
            sage: f.list()
            [1, 9, 12, 8]
        """
        if copy:
            return list(self._coeffs)
        else:
            return self._coeffs

    def degree(self, gen=None):
        """
        EXAMPLES::

            sage: R.<x> = RDF[]
            sage: f = (1+2*x^7)^5
            sage: f.degree()
            35

        TESTS:

        Check that :issue:`12552` is fixed::

            sage: type(f.degree())
            <class 'sage.rings.integer.Integer'>
        """
        return smallInteger(len(self._coeffs) - 1)

    def shift(self, Py_ssize_t n):
        r"""
        Return this polynomial multiplied by the power `x^n`.

        If `n` is negative, terms below `x^n` will be discarded. Does
        not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(PolynomialRing(QQ,'y'), 'x')
            sage: p = x^2 + 2*x + 4
            sage: type(p)
            <class 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
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
            output = [self.base_ring().zero()] * n
            output.extend(self._coeffs)
            return self._new_c(output, self._parent)
        if n < 0:
            if n > len(self._coeffs) - 1:
                return self._parent([])
            else:
                return self._new_c(self._coeffs[-int(n):], self._parent)

    @coerce_binop
    def quo_rem(self, other):
        """
        Return the quotient and remainder of the Euclidean division of
        ``self`` and ``other``.

        Raises a :exc:`ZeroDivisionError` if ``other`` is zero. Raises an
        :exc:`ArithmeticError` if the division is not exact.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: R.<y> = P[]
            sage: f = y^10 + R.random_element(9)
            sage: g = y^5 + R.random_element(4)
            sage: q, r = f.quo_rem(g)
            sage: f == q*g + r
            True
            sage: g = x*y^5
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ArithmeticError: division non exact (consider coercing
            to polynomials over the fraction field)
            sage: g = 0
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero polynomial

        Polynomials over noncommutative rings are also allowed
        (after :issue:`34733`)::

            sage: # needs sage.combinat sage.modules
            sage: HH = QuaternionAlgebra(QQ, -1, -1)
            sage: P.<x> = HH[]
            sage: f = P.random_element(5)
            sage: g = P.random_element((0, 5))
            sage: q, r = f.quo_rem(g)
            sage: f == q*g + r
            True

        TESTS:

        The following shows that :issue:`16649` is indeed fixed. ::

            sage: P.<x> = QQ[]
            sage: R.<y> = P[]
            sage: f = (2*x^3+1)*y^2 + (x^2-x+3)*y + (3*x+2)
            sage: g = (-1/13*x^2 - x)*y^2 + (-x^2 + 3*x - 155/4)*y - x - 1
            sage: h = f * g
            sage: h.quo_rem(f)
            ((-1/13*x^2 - x)*y^2 + (-x^2 + 3*x - 155/4)*y - x - 1, 0)
            sage: h += (2/3*x^2-3*x+1)*y + 7/17*x+6/5
            sage: q, r = h.quo_rem(f)
            sage: h == q*f + r and r.degree() < f.degree()
            True

        :issue:`26907`::

            sage: P.<x> = ZZ[]
            sage: R.<y> = P[]
            sage: a = 3*y + 1
            sage: a//a
            1
        """
        if other.is_zero():
            raise ZeroDivisionError("division by zero polynomial")
        if self.is_zero():
            return self, self

        R = self._parent.base_ring()
        cdef list x = list((<Polynomial_generic_dense>self)._coeffs) # make a copy
        cdef list y = (<Polynomial_generic_dense>other)._coeffs
        cdef Py_ssize_t m = len(x)  # deg(self)=m-1
        cdef Py_ssize_t n = len(y)  # deg(other)=n-1
        if m < n:
            return self._parent.zero(), self

        cdef list quo = []
        cdef bint convert = False
        cdef Py_ssize_t j, k
        try:
            inv = y[n-1].inverse_of_unit()
        except (ArithmeticError, ValueError):
            inv = ~(y[n-1])
            convert = True
        if convert:
            for k from m-n >= k >= 0:
                q = x[n+k-1] * inv
                try:
                    q = R(q)
                except TypeError:
                    raise ArithmeticError("division non exact (consider coercing to polynomials over the fraction field)")
                for j from n+k-2 >= j >= k:
                    x[j] -= q * y[j-k]
                quo.append(q)
        else:
            for k from m-n >= k >= 0:
                q = x[n+k-1] * inv
                for j from n+k-2 >= j >= k:
                    x[j] -= q * y[j-k]
                quo.append(q)
        quo.reverse()

        return self._new_c(quo,self._parent), self._new_c(x,self._parent)._inplace_truncate(n-1)

    cpdef Polynomial truncate(self, long n):
        r"""
        Return the polynomial of degree `< n` which is equivalent
        to ``self`` modulo `x^n`.

        EXAMPLES::

            sage: S.<q> = QQ['t']['q']
            sage: f = (1 + q^10 + q^11 + q^12).truncate(11); f
            q^10 + 1
            sage: f = (1 + q^10 + q^100).truncate(50); f
            q^10 + 1
            sage: f.degree()
            10
            sage: f = (1 + q^10 + q^100).truncate(500); f
            q^100 + q^10 + 1

        TESTS:

        Make sure we're not actually testing a specialized
        implementation.

        ::

            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
        """
        l = len(self._coeffs)
        if n > l:
            n = l
        while n > 0 and not self._coeffs[n-1]:
            n -= 1
        return self._new_c(self._coeffs[:n], self._parent)

    cdef _inplace_truncate(self, long n):
        if n < len(self._coeffs):
            while n > 0 and not self._coeffs[n-1]:
                n -= 1
        self._coeffs = self._coeffs[:n]
        return self


def make_generic_polynomial(parent, coeffs):
    return parent(coeffs)


@cached_function
def universal_discriminant(n):
    r"""
    Return the discriminant of the 'universal' univariate polynomial
    `a_n x^n + \cdots + a_1 x + a_0` in `\ZZ[a_0, \ldots, a_n][x]`.

    INPUT:

    - ``n`` -- degree of the polynomial

    OUTPUT:

    The discriminant as a polynomial in `n + 1` variables over `\ZZ`.
    The result will be cached, so subsequent computations of
    discriminants of the same degree will be faster.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.rings.polynomial.polynomial_element import universal_discriminant
        sage: universal_discriminant(1)
        1
        sage: universal_discriminant(2)
        a1^2 - 4*a0*a2
        sage: universal_discriminant(3)
        a1^2*a2^2 - 4*a0*a2^3 - 4*a1^3*a3 + 18*a0*a1*a2*a3 - 27*a0^2*a3^2
        sage: universal_discriminant(4).degrees()
        (3, 4, 4, 4, 3)

    .. SEEALSO::
        :meth:`Polynomial.discriminant`
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    pr1 = PolynomialRing(ZZ, n + 1, 'a')
    pr2 = PolynomialRing(pr1, 'x')
    p = pr2(list(pr1.gens()))
    return (1 - (n&2))*p.resultant(p.derivative())//pr1.gen(n)


cpdef Polynomial generic_power_trunc(Polynomial p, Integer n, long prec):
    r"""
    Generic truncated power algorithm.

    INPUT:

    - ``p`` -- a polynomial

    - ``n`` -- integer (of type :class:`sage.rings.integer.Integer`)

    - ``prec`` -- a precision (should fit into a C ``long``)

    TESTS:

    Comparison with flint for polynomials over integers and finite field::

        sage: from sage.rings.polynomial.polynomial_element import generic_power_trunc

        sage: for S in [ZZ, GF(3)]:             # known bug, not tested (see :issue:`32075`)
        ....:     R = PolynomialRing(S, 'x')
        ....:     for _ in range(100):
        ....:         p = R.random_element()
        ....:         n = ZZ.random_element(0, 100)
        ....:         prec = ZZ.random_element(0, 100)
        ....:         assert p.power_trunc(n, prec) == generic_power_trunc(p, n, prec), "p = {} n = {} prec = {}".format(p, n, prec)
    """
    if mpz_sgn(n.value) < 0:
        raise ValueError("n must be a nonnegative integer")
    elif prec <= 0:
        return p._parent.zero()

    if mpz_cmp_ui(n.value, 4) < 0:
        # These cases will probably be called often
        # and don't benefit from the code below
        if mpz_cmp_ui(n.value, 0) == 0:
            return p.parent().one()
        elif mpz_cmp_ui(n.value, 1) == 0:
            return p.truncate(prec)
        elif mpz_cmp_ui(n.value, 2) == 0:
            return p._mul_trunc_(p, prec)
        elif mpz_cmp_ui(n.value, 3) == 0:
            return p._mul_trunc_(p, prec)._mul_trunc_(p, prec)

    # check for idempotence, and store the result otherwise
    cdef Polynomial a = p.truncate(prec)
    cdef Polynomial aa = a._mul_trunc_(a, prec)
    if aa == a:
        return a

    # since we've computed a^2, let's start squaring there
    # so, let's keep the least-significant bit around, just
    # in case.
    cdef int mul_to_do = mpz_tstbit(n.value, 0)
    cdef mp_bitcnt_t i = 1
    cdef mp_bitcnt_t size = mpz_sizeinbase(n.value, 2)

    # One multiplication can be saved by starting with
    # the second-smallest power needed rather than with 1
    # we've already squared a, so let's start there.
    cdef Polynomial apow = aa
    while not mpz_tstbit(n.value, i):
        apow = apow._mul_trunc_(apow, prec)
        i += 1
    cdef Polynomial power = apow
    i += 1

    # now multiply that least-significant bit in...
    if mul_to_do:
        power = power._mul_trunc_(a, prec)

    # and this is straight from the book.
    while i < size:
        apow = apow._mul_trunc_(apow, prec)
        if mpz_tstbit(n.value, i):
            power = power._mul_trunc_(apow, prec)
        i += 1

    return power

cpdef list _dict_to_list(dict x, zero):
    """
    Convert a dict to a list.

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_element import _dict_to_list
        sage: _dict_to_list({3:-1, 0:5}, 0)
        [5, 0, 0, -1]
    """
    if not x:
        return []
    n = max(x.keys())
    cdef list v
    if isinstance(n, tuple): # a mpoly dict
        n = n[0]
        v = [zero] * (n+1)
        for i, z in x.items():
            v[i[0]] = z
    else:
        v = [zero] * (n+1)
        for i, z in x.items():
            v[i] = z
    return v

cdef class Polynomial_generic_dense_inexact(Polynomial_generic_dense):
    """
    A dense polynomial over an inexact ring.

    AUTHOR:

    - Xavier Caruso (2013-03)
    """
    cdef int _normalize(self) except -1:
        r"""
        TESTS::

        Coefficients indistinguishable from 0 are not removed.

            sage: R = Zp(5)                                                             # needs sage.rings.padics
            sage: S.<x> = R[]                                                           # needs sage.rings.padics
            sage: S([1, R(0, 20)])                                                      # needs sage.rings.padics
            O(5^20)*x + 1 + O(5^20)
        """
        cdef list x = self._coeffs
        cdef Py_ssize_t n = len(x) - 1
        cdef RingElement c
        while n >= 0:
            c = x[n]
            if c.is_zero() and c.precision_absolute() is infinity.Infinity:
                del x[n]
                n -= 1
            else:
                break

    def degree(self, secure=False):
        r"""
        INPUT:

        - ``secure`` -- boolean (default: ``False``)

        OUTPUT: the degree of ``self``

        If ``secure`` is ``True`` and the degree of this polynomial
        is not determined (because the leading coefficient is
        indistinguishable from 0), an error is raised

        If ``secure`` is ``False``, the returned value is the largest
        `n` so that the coefficient of `x^n` does not compare equal
        to `0`.

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: K = Qp(3, 10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.degree()
            1
            sage: (f - T).degree()
            0
            sage: (f - T).degree(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: the leading coefficient is indistinguishable from 0

            sage: # needs sage.rings.padics
            sage: x = O(3^5)
            sage: li = [3^i * x for i in range(0,5)]; li
            [O(3^5), O(3^6), O(3^7), O(3^8), O(3^9)]
            sage: f = R(li); f
            O(3^9)*T^4 + O(3^8)*T^3 + O(3^7)*T^2 + O(3^6)*T + O(3^5)
            sage: f.degree()
            -1
            sage: f.degree(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: the leading coefficient is indistinguishable from 0

        TESTS:

        Ensure that :issue:`37621` is fixed::

            sage: k.<x> = QQ[]
            sage: K = Qp(11,5)
            sage: L.<a> = K.extension(x^20-11)
            sage: R.<x> = L[]
            sage: f = R.random_element()
            sage: type(f.degree())
            <class 'sage.rings.integer.Integer'>

        AUTHOR:

        - Xavier Caruso (2013-03)
        """
        coeffs = self._coeffs
        d = len(coeffs) - 1
        while d >= 0:
            c = coeffs[d]
            if c.is_zero():
                if secure:
                    from sage.rings.padics.precision_error import PrecisionError
                    raise PrecisionError("the leading coefficient is indistinguishable from 0")
                else:
                    d -= 1
            else:
                break
        return Integer(d)

    def prec_degree(self):
        r"""
        Return the largest `n` so that precision information is
        stored about the coefficient of `x^n`.

        Always greater than or equal to degree.

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: K = Qp(3, 10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.degree()
            1
            sage: f.prec_degree()
            1

            sage: g = f - T; g                                                          # needs sage.rings.padics
            O(3^10)*T + 2 + O(3^10)
            sage: g.degree()                                                            # needs sage.rings.padics
            0
            sage: g.prec_degree()                                                       # needs sage.rings.padics
            1

        AUTHOR:

        - Xavier Caruso (2013-03)
        """
        return len(self._coeffs) - 1


cdef class ConstantPolynomialSection(Map):
    """
    This class is used for conversion from a polynomial ring to its base ring.

    Since :issue:`9944`, it calls the :meth:`constant_coefficient` method,
    which can be optimized for a particular polynomial type.

    EXAMPLES::

        sage: P0.<y_1> = GF(3)[]
        sage: P1.<y_2,y_1,y_0> = GF(3)[]
        sage: P0(-y_1)
        2*y_1

        sage: phi = GF(3).convert_map_from(P0); phi
        Generic map:
          From: Univariate Polynomial Ring in y_1 over Finite Field of size 3
          To:   Finite Field of size 3
        sage: type(phi)
        <class 'sage.rings.polynomial.polynomial_element.ConstantPolynomialSection'>
        sage: phi(P0.one())
        1
        sage: phi(y_1)
        Traceback (most recent call last):
        ...
        TypeError: y_1 is not a constant polynomial
    """
    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_element import ConstantPolynomialSection
            sage: R.<x> = QQ[]
            sage: m = ConstantPolynomialSection(R, QQ); m
            Generic map:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
            sage: m(x-x+1/2)  # implicit
            1/2
            sage: m(x-x)
            0
            sage: m(x)
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial
        """
        if x.degree() <= 0:
            try:
                return <Element>(x.constant_coefficient())
            except AttributeError:
                return <Element>((<Polynomial>x).constant_coefficient())
        else:
            raise TypeError(f"{x} is not a constant polynomial")

cdef class PolynomialBaseringInjection(Morphism):
    """
    This class is used for conversion from a ring to a polynomial
    over that ring.

    It calls the :meth:`_new_constant_poly` method on the generator,
    which should be optimized for a particular polynomial type.

    Technically, it should be a method of the polynomial ring, but
    few polynomial rings are Cython classes, and so, as a method
    of a Cython polynomial class, it is faster.

    EXAMPLES:

    We demonstrate that most polynomial ring classes use
    polynomial base injection maps for coercion. They are
    supposed to be the fastest maps for that purpose. See
    :issue:`9944`. ::

        sage: # needs sage.rings.padics
        sage: R.<x> = Qp(3)[]
        sage: R.coerce_map_from(R.base_ring())
        Polynomial base injection morphism:
          From: 3-adic Field with capped relative precision 20
          To:   Univariate Polynomial Ring in x over
                3-adic Field with capped relative precision 20
        sage: R.<x,y> = Qp(3)[]
        sage: R.coerce_map_from(R.base_ring())
        Polynomial base injection morphism:
          From: 3-adic Field with capped relative precision 20
          To:   Multivariate Polynomial Ring in x, y over
                3-adic Field with capped relative precision 20

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

    By :issue:`9944`, there are now only very few exceptions::

        sage: PolynomialRing(QQ,names=[]).coerce_map_from(QQ)
        Call morphism:
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

            sage: # needs sage.rings.padics
            sage: R.<t> = Qp(2)[]
            sage: f = R.convert_map_from(R.base_ring())    # indirect doctest
            sage: f(Qp(2).one()*3)
            1 + 2 + O(2^20)
            sage: (Qp(2).one()*3)*t
            (1 + 2 + O(2^20))*t
        """
        assert codomain.base_ring() is domain, "domain must be basering"
        Morphism.__init__(self, domain, codomain)
        self._an_element = codomain.gen()
        self._repr_type_str = "Polynomial base injection"
        self._new_constant_poly_ = self._an_element._new_constant_poly

    cdef dict _extra_slots(self):
        """
        EXAMPLES::

            sage: phi = QQ['x'].coerce_map_from(QQ)   # indirect doctest
            sage: phi
            Polynomial base injection morphism:
              From: Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field
            sage: phi(3/1)
            3
        """
        slots = Morphism._extra_slots(self)
        slots.update(
                _an_element=self._an_element,
                _new_constant_poly_=self._new_constant_poly_)
        return slots

    cdef _update_slots(self, dict _slots):
        """
        EXAMPLES::

            sage: phi = QQ['x'].coerce_map_from(QQ)  # indirect doctest
            sage: phi
            Polynomial base injection morphism:
              From: Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field
            sage: phi(3/1)
            3
        """
        Morphism._update_slots(self, _slots)
        self._an_element = _slots['_an_element']
        self._new_constant_poly_ = _slots['_new_constant_poly_']

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
            sage: m = PolynomialBaseringInjection(ZZ, ZZ['x']); m
            Polynomial base injection morphism:
              From: Integer Ring
              To:   Univariate Polynomial Ring in x over Integer Ring
            sage: m(2)  # indirect doctest
            2
            sage: parent(m(2))
            Univariate Polynomial Ring in x over Integer Ring
        """
        return self._new_constant_poly_(x, self._codomain)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection
            sage: m = PolynomialBaseringInjection(Qp(5), Qp(5)['x'])                    # needs sage.rings.padics
            sage: m(1 + O(5^11), absprec=5)   # indirect doctest                        # needs sage.rings.padics
            1 + O(5^11)
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
            <class 'sage.rings.polynomial.polynomial_element.ConstantPolynomialSection'>
        """
        return ConstantPolynomialSection(self._codomain, self.domain())

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: S.<y> = R[]
            sage: S.coerce_map_from(R).is_injective()
            True

        Check that :issue:`23203` has been resolved::

            sage: R.is_subring(S)  # indirect doctest
            True
        """
        return True

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R.coerce_map_from(ZZ).is_surjective()
            False
        """
        return False

cpdef bint polynomial_is_variable(x) noexcept:
    r"""
    Test whether the given polynomial is a variable of its parent ring.

    Implemented for instances of :class:`Polynomial` and :class:`MPolynomial`.

    .. SEEALSO::

        - :meth:`sage.rings.polynomial.polynomial_element.Polynomial.is_gen`
        - :meth:`sage.rings.polynomial.multi_polynomial.MPolynomial.is_generator`

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_element import polynomial_is_variable
        sage: R.<x> = QQ[]
        sage: polynomial_is_variable(x)
        True
        sage: polynomial_is_variable(R([0,1]))
        True
        sage: polynomial_is_variable(x^2)
        False
        sage: polynomial_is_variable(R(42))
        False

    ::

        sage: R.<y,z> = QQ[]
        sage: polynomial_is_variable(y)
        True
        sage: polynomial_is_variable(z)
        True
        sage: polynomial_is_variable(y^2)
        False
        sage: polynomial_is_variable(y+z)
        False
        sage: polynomial_is_variable(R(42))
        False

    ::

        sage: polynomial_is_variable(42)
        False
    """
    if isinstance(x, Polynomial):
        return (x.is_gen()
                or (x.degree() == 1 and x[0].is_zero() and x[1].is_one()))
    if isinstance(x, MPolynomial):
        return x.is_generator()
    return False
