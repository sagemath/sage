"""
Finite field elements implemented via PARI's FFELT type

AUTHORS:

- Peter Bruin (June 2013): initial version, based on
  element_ext_pari.py by William Stein et al. and
  element_ntl_gf2e.pyx by Martin Albrecht.
"""

#*****************************************************************************
#      Copyright (C) 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "sage/ext/stdsage.pxi"
include "sage/libs/pari/pari_err.pxi"

from element_base cimport FinitePolyExtElement
from integer_mod import IntegerMod_abstract

import sage.libs.pari
import sage.rings.integer
from sage.interfaces.gap import is_GapElement
from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport PariInstance
from sage.modules.free_module_element import FreeModuleElement
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial
from sage.rings.rational import Rational
from sage.structure.element cimport Element, ModuleElement, RingElement

cdef long mpz_t_offset = sage.rings.integer.mpz_t_offset_python

cdef PariInstance pari = sage.libs.pari.pari_instance.pari

cdef extern from "sage/libs/pari/misc.h":
    int gcmp_sage(GEN x, GEN y)


cdef GEN _INT_to_FFELT(GEN g, GEN x) except NULL:
    """
    Convert the t_INT `x` to an element of the field of definition of
    the t_FFELT `g`.

    This function must be called within pari_catch_sig_on() ... pari_catch_sig_off().
    """
    cdef GEN f, p = gel(g, 4), result
    cdef long t

    x = modii(x, p)
    if gequal0(x):
        return FF_zero(g)
    elif gequal1(x):
        return FF_1(g)
    else:
        # In characteristic 2, we have already dealt with the
        # two possible values of x, so we may assume that the
        # characteristic is > 2.
        t = g[1]  # codeword: t_FF_FpXQ, t_FF_Flxq, t_FF_F2xq
        if t == t_FF_FpXQ:
            f = cgetg(3, t_POL)
            set_gel(f, 1, gmael(g, 2, 1))
            set_gel(f, 2, x)
        elif t == t_FF_Flxq:
            f = cgetg(3, t_VECSMALL)
            set_gel(f, 1, gmael(g, 2, 1))
            f[2] = itos(x)
        else:
            pari_catch_sig_off()
            raise TypeError("unknown PARI finite field type")
        result = cgetg(5, t_FFELT)
        result[1] = t
        set_gel(result, 2, f)
        set_gel(result, 3, gel(g, 3))  # modulus
        set_gel(result, 4, p)
        return result

cdef class FiniteFieldElement_pari_ffelt(FinitePolyExtElement):
    """
    An element of a finite field.

    EXAMPLE::

        sage: K = FiniteField(10007^10, 'a', impl='pari_ffelt')
        sage: a = K.gen(); a
        a
        sage: type(a)
        <type 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt'>

    TESTS::

        sage: n = 63
        sage: m = 3;
        sage: K.<a> = GF(2^n, impl='pari_ffelt')
        sage: f = conway_polynomial(2, n)
        sage: f(a) == 0
        True
        sage: e = (2^n - 1) / (2^m - 1)
        sage: conway_polynomial(2, m)(a^e) == 0
        True

        sage: K.<a> = FiniteField(2^16, impl='pari_ffelt')
        sage: K(0).is_zero()
        True
        sage: (a - a).is_zero()
        True
        sage: a - a
        0
        sage: a == a
        True
        sage: a - a == 0
        True
        sage: a - a == K(0)
        True
        sage: TestSuite(a).run()

    Test creating elements from basic Python types::

        sage: K.<a> = FiniteField(7^20, impl='pari_ffelt')
        sage: K(int(8))
        1
        sage: K(long(-2^300))
        6
    """

    def __init__(FiniteFieldElement_pari_ffelt self, object parent, object x):
        """
        Initialise ``self`` with the given ``parent`` and value
        converted from ``x``.

        This is called when constructing elements from Python.

        TEST::

            sage: from sage.rings.finite_rings.element_pari_ffelt import FiniteFieldElement_pari_ffelt
            sage: K = FiniteField(101^2, 'a', impl='pari_ffelt')
            sage: x = FiniteFieldElement_pari_ffelt(K, 'a + 1')
            sage: x
            a + 1
        """
        # FinitePolyExtElement.__init__(self, parent)
        self._parent = parent
        self.construct_from(x)

    # The Cython constructor __cinit__ is not necessary: according to
    # the Cython documentation, C attributes are initialised to 0.
    # def __cinit__(FiniteFieldElement_pari_ffelt self):
    #     self.block = NULL

    def __dealloc__(FiniteFieldElement_pari_ffelt self):
        """
        Cython deconstructor.
        """
        if self.block:
            sage_free(self.block)

    cdef FiniteFieldElement_pari_ffelt _new(FiniteFieldElement_pari_ffelt self):
        """
        Create an empty element with the same parent as ``self``.
        """
        cdef FiniteFieldElement_pari_ffelt x
        x = FiniteFieldElement_pari_ffelt.__new__(FiniteFieldElement_pari_ffelt)
        x._parent = self._parent
        return x

    cdef void construct(FiniteFieldElement_pari_ffelt self, GEN g):
        """
        Initialise ``self`` to the FFELT ``g``, reset the PARI stack,
        and call pari_catch_sig_off().

        This should be called exactly once on every instance.
        """
        self.val = pari.deepcopy_to_python_heap(g, <pari_sp*>&self.block)
        pari.clear_stack()

    cdef void construct_from(FiniteFieldElement_pari_ffelt self, object x) except *:
        """
        Initialise ``self`` to an FFELT constructed from the Sage
        object `x`.
        """
        cdef GEN f, g, result, x_GEN
        cdef long i, n, t
        cdef Integer xi

        if isinstance(x, FiniteFieldElement_pari_ffelt):
            if self._parent is (<FiniteFieldElement_pari_ffelt>x)._parent:
                pari_catch_sig_on()
                self.construct((<FiniteFieldElement_pari_ffelt>x).val)
            else:
                # This is where we *would* do coercion from one finite field to another...
                raise TypeError("no coercion defined")

        elif isinstance(x, Integer):
            g = (<pari_gen>self._parent._gen_pari).g
            pari_catch_sig_on()
            x_GEN = pari._new_GEN_from_mpz_t(<void *>x + mpz_t_offset)
            self.construct(_INT_to_FFELT(g, x_GEN))

        elif isinstance(x, int) or isinstance(x, long):
            g = (<pari_gen>self._parent._gen_pari).g
            x = pari(x)
            pari_catch_sig_on()
            x_GEN = (<pari_gen>x).g
            self.construct(_INT_to_FFELT(g, x_GEN))

        elif isinstance(x, IntegerMod_abstract):
            if self._parent.characteristic().divides(x.modulus()):
                g = (<pari_gen>self._parent._gen_pari).g
                x = Integer(x)
                pari_catch_sig_on()
                x_GEN = pari._new_GEN_from_mpz_t(<void *>x + mpz_t_offset)
                self.construct(_INT_to_FFELT(g, x_GEN))
            else:
                raise TypeError("no coercion defined")

        elif x is None:
            g = (<pari_gen>self._parent._gen_pari).g
            pari_catch_sig_on()
            self.construct(FF_zero(g))

        elif isinstance(x, pari_gen):
            g = (<pari_gen>self._parent._gen_pari).g
            x_GEN = (<pari_gen>x).g

            pari_catch_sig_on()
            if gequal0(x_GEN):
                self.construct(FF_zero(g))
                return
            elif gequal1(x_GEN):
                self.construct(FF_1(g))
                return

            t = typ(x_GEN)
            if t == t_FFELT and FF_samefield(x_GEN, g):
                self.construct(x_GEN)
            elif t == t_INT:
                self.construct(_INT_to_FFELT(g, x_GEN))
            elif t == t_INTMOD and gequal0(modii(gel(x_GEN, 1), FF_p_i(g))):
                self.construct(_INT_to_FFELT(g, gel(x_GEN, 2)))
            elif t == t_FRAC and not gequal0(modii(gel(x_GEN, 2), FF_p_i(g))):
                self.construct(FF_div(_INT_to_FFELT(g, gel(x_GEN, 1)),
                                      _INT_to_FFELT(g, gel(x_GEN, 2))))
            else:
                pari_catch_sig_off()
                raise TypeError("no coercion defined")

        elif (isinstance(x, FreeModuleElement)
              and x.parent() is self._parent.vector_space()):
            g = (<pari_gen>self._parent._gen_pari).g
            t = g[1]  # codeword: t_FF_FpXQ, t_FF_Flxq, t_FF_F2xq
            n = len(x)
            while n > 0 and x[n - 1] == 0:
                n -= 1
            pari_catch_sig_on()
            if n == 0:
                self.construct(FF_zero(g))
                return
            if t == t_FF_FpXQ:
                f = cgetg(n + 2, t_POL)
                set_gel(f, 1, gmael(g, 2, 1))
                for i in xrange(n):
                    xi = Integer(x[i])
                    set_gel(f, i + 2, pari._new_GEN_from_mpz_t(<void *>xi + mpz_t_offset))
            elif t == t_FF_Flxq or t == t_FF_F2xq:
                f = cgetg(n + 2, t_VECSMALL)
                set_gel(f, 1, gmael(g, 2, 1))
                for i in xrange(n):
                    f[i + 2] = x[i]
                if t == t_FF_F2xq:
                    f = Flx_to_F2x(f)
            else:
                pari_catch_sig_off()
                raise TypeError("unknown PARI finite field type")
            result = cgetg(5, t_FFELT)
            result[1] = t
            set_gel(result, 2, f)
            set_gel(result, 3, gel(g, 3))  # modulus
            set_gel(result, 4, gel(g, 4))  # p
            self.construct(result)

        elif isinstance(x, Rational):
            self.construct_from(x % self._parent.characteristic())

        elif isinstance(x, Polynomial):
            if x.base_ring() is not self._parent.base_ring():
                x = x.change_ring(self._parent.base_ring())
            self.construct_from(x.substitute(self._parent.gen()))

        elif isinstance(x, MPolynomial) and x.is_constant():
            self.construct_from(x.constant_coefficient())

        elif isinstance(x, list):
            if len(x) == self._parent.degree():
                self.construct_from(self._parent.vector_space()(x))
            else:
                Fp = self._parent.base_ring()
                self.construct_from(self._parent.polynomial_ring()([Fp(y) for y in x]))

        elif isinstance(x, str):
            self.construct_from(self._parent.polynomial_ring()(x))

        elif is_GapElement(x):
            from sage.interfaces.gap import gfq_gap_to_sage
            try:
                self.construct_from(gfq_gap_to_sage(x, self._parent))
            except (ValueError, IndexError, TypeError):
                raise TypeError("no coercion defined")

        else:
            raise TypeError("no coercion defined")

    def _repr_(FiniteFieldElement_pari_ffelt self):
        """
        Return the string representation of ``self``.

        EXAMPLE::

            sage: k.<c> = GF(3^17, impl='pari_ffelt')
            sage: c^20  # indirect doctest
            c^4 + 2*c^3
        """
        pari_catch_sig_on()
        return str(pari.new_gen(self.val))

    def __hash__(FiniteFieldElement_pari_ffelt self):
        """
        Return the hash of ``self``.  This is by definition equal to
        the hash of ``self.polynomial()``.

        EXAMPLE::

            sage: k.<a> = GF(3^15, impl='pari_ffelt')
            sage: R = GF(3)['a']; aa = R.gen()
            sage: hash(a^2 + 1) == hash(aa^2 + 1)
            True
        """
        return hash(self.polynomial())

    def __reduce__(FiniteFieldElement_pari_ffelt self):
        """
        For pickling.

        TEST::

            sage: K.<a> = FiniteField(10007^10, impl='pari_ffelt')
            sage: loads(a.dumps()) == a
            True
        """
        return unpickle_FiniteFieldElement_pari_ffelt, (self._parent, str(self))

    def __copy__(FiniteFieldElement_pari_ffelt self):
        """
        Return a copy of ``self``.

        TESTS::

            sage: k.<a> = FiniteField(3^3, impl='pari_ffelt')
            sage: a
            a
            sage: b = copy(a); b
            a
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(self.val)
        return x

    cdef int _cmp_c_impl(FiniteFieldElement_pari_ffelt self, Element other) except -2:
        """
        Comparison of finite field elements.

        TESTS::

            sage: k.<a> = FiniteField(3^3, impl='pari_ffelt')
            sage: a == 1
            False
            sage: a^0 == 1
            True
            sage: a == a
            True
            sage: a < a^2
            True
            sage: a > a^2
            False
        """
        return gcmp_sage(self.val, (<FiniteFieldElement_pari_ffelt>other).val)

    def __richcmp__(FiniteFieldElement_pari_ffelt left, object right, int op):
        """
        Rich comparison of finite field elements.

        EXAMPLE::

            sage: k.<a> = GF(2^20, impl='pari_ffelt')
            sage: e = k.random_element()
            sage: f = loads(dumps(e))
            sage: e is f
            False
            sage: e == f
            True
            sage: e != (e + 1)
            True

        .. NOTE::

            Finite fields are unordered.  However, for the purpose of
            this function, we adopt the lexicographic ordering on the
            representing polynomials.

        EXAMPLE::

            sage: K.<a> = GF(2^100, impl='pari_ffelt')
            sage: a < a^2
            True
            sage: a > a^2
            False
            sage: a+1 > a^2
            False
            sage: a+1 < a^2
            True
            sage: a+1 < a
            False
            sage: a+1 == a
            False
            sage: a == a
            True
        """
        return (<Element>left)._richcmp(right, op)

    cpdef ModuleElement _add_(FiniteFieldElement_pari_ffelt self, ModuleElement right):
        """
        Addition.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: a + a^2 # indirect doctest
            a^2 + a
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_add((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    cpdef ModuleElement _sub_(FiniteFieldElement_pari_ffelt self, ModuleElement right):
        """
        Subtraction.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: a - a # indirect doctest
            0
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_sub((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    cpdef RingElement _mul_(FiniteFieldElement_pari_ffelt self, RingElement right):
        """
        Multiplication.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: (a^12 + 1)*(a^15 - 1) # indirect doctest
            a^15 + 2*a^12 + a^11 + 2*a^10 + 2
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_mul((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    cpdef RingElement _div_(FiniteFieldElement_pari_ffelt self, RingElement right):
        """
        Division.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: (a - 1) / (a + 1) # indirect doctest
            2*a^16 + a^15 + 2*a^14 + a^13 + 2*a^12 + a^11 + 2*a^10 + a^9 + 2*a^8 + a^7 + 2*a^6 + a^5 + 2*a^4 + a^3 + 2*a^2 + a + 1
        """
        if FF_equal0((<FiniteFieldElement_pari_ffelt>right).val):
            raise ZeroDivisionError
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_div((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    def is_zero(FiniteFieldElement_pari_ffelt self):
        """
        Return ``True`` if ``self`` equals 0.

        EXAMPLE::

            sage: F.<a> = FiniteField(5^3, impl='pari_ffelt')
            sage: a.is_zero()
            False
            sage: (a - a).is_zero()
            True
        """
        return bool(FF_equal0(self.val))

    def is_one(FiniteFieldElement_pari_ffelt self):
        """
        Return ``True`` if ``self`` equals 1.

        EXAMPLE::

            sage: F.<a> = FiniteField(5^3, impl='pari_ffelt')
            sage: a.is_one()
            False
            sage: (a/a).is_one()
            True
        """
        return bool(FF_equal1(self.val))

    def is_unit(FiniteFieldElement_pari_ffelt self):
        """
        Return ``True`` if ``self`` is non-zero.

        EXAMPLE::

            sage: F.<a> = FiniteField(5^3, impl='pari_ffelt')
            sage: a.is_unit()
            True
        """
        return not bool(FF_equal0(self.val))

    __nonzero__ = is_unit

    def __pos__(FiniteFieldElement_pari_ffelt self):
        """
        Unitary positive operator...

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: +a
            a
        """
        return self

    def __neg__(FiniteFieldElement_pari_ffelt self):
        """
        Negation.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: -a
            2*a
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_neg_i((<FiniteFieldElement_pari_ffelt>self).val))
        return x

    def __invert__(FiniteFieldElement_pari_ffelt self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLE::

            sage: k.<a> = FiniteField(3^2, impl='pari_ffelt')
            sage: ~a
            a + 2
            sage: (a+1)*a
            2*a + 1
            sage: ~((2*a)/a)
            2
        """
        if FF_equal0(self.val):
            raise ZeroDivisionError
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_inv((<FiniteFieldElement_pari_ffelt>self).val))
        return x

    def __pow__(FiniteFieldElement_pari_ffelt self, object exp, object other):
        """
        Exponentiation.

        TESTS::

            sage: K.<a> = GF(5^10, impl='pari_ffelt')
            sage: n = (2*a)/a
            sage: n^-15
            2

        Large exponents are not a problem::

            sage: e = 3^10000
            sage: a^e
            2*a^9 + a^5 + 4*a^4 + 4*a^3 + a^2 + 3*a
            sage: a^(e % (5^10 - 1))
            2*a^9 + a^5 + 4*a^4 + 4*a^3 + a^2 + 3*a
        """
        if exp == 0:
            return self._parent.one_element()
        if exp < 0 and FF_equal0(self.val):
            raise ZeroDivisionError
        exp = pari(exp)
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        pari_catch_sig_on()
        x.construct(FF_pow(self.val, (<pari_gen>exp).g))
        return x

    def polynomial(FiniteFieldElement_pari_ffelt self):
        """
        Return the unique representative of ``self`` as a polynomial
        over the prime field whose degree is less than the degree of
        the finite field over its prime field.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^2, impl='pari_ffelt')
            sage: pol = a.polynomial()
            sage: pol
            a
            sage: parent(pol)
            Univariate Polynomial Ring in a over Finite Field of size 3

        ::

            sage: k = FiniteField(3^4, 'alpha', impl='pari_ffelt')
            sage: a = k.gen()
            sage: a.polynomial()
            alpha
            sage: (a**2 + 1).polynomial()
            alpha^2 + 1
            sage: (a**2 + 1).polynomial().parent()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        pari_catch_sig_on()
        return self._parent.polynomial_ring()(pari.new_gen(FF_to_FpXQ_i(self.val)))

    def charpoly(FiniteFieldElement_pari_ffelt self, object var='x'):
        """
        Return the characteristic polynomial of ``self``.

        INPUT:

        - ``var`` -- string (default: 'x'): variable name to use.

        EXAMPLE::

            sage: R.<x> = PolynomialRing(FiniteField(3))
            sage: F.<a> = FiniteField(3^2, modulus=x^2 + 1)
            sage: a.charpoly('y')
            y^2 + 1
        """
        pari_catch_sig_on()
        return self._parent.polynomial_ring(var)(pari.new_gen(FF_charpoly(self.val)))

    def is_square(FiniteFieldElement_pari_ffelt self):
        """
        Return ``True`` if and only if ``self`` is a square in the
        finite field.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^2, impl='pari_ffelt')
            sage: a.is_square()
            False
            sage: (a**2).is_square()
            True

            sage: k.<a> = FiniteField(2^2, impl='pari_ffelt')
            sage: (a**2).is_square()
            True

            sage: k.<a> = FiniteField(17^5, impl='pari_ffelt')
            sage: (a**2).is_square()
            True
            sage: a.is_square()
            False
            sage: k(0).is_square()
            True
        """
        cdef long i
        pari_catch_sig_on()
        i = FF_issquare(self.val)
        pari_catch_sig_off()
        return bool(i)

    def sqrt(FiniteFieldElement_pari_ffelt self, extend=False, all=False):
        """
        Return a square root of ``self``, if it exists.

        INPUT:

        - ``extend`` -- bool (default: ``False``)

           .. WARNING::

               This option is not implemented.

        - ``all`` - bool (default: ``False``)

        OUTPUT:

        A square root of ``self``, if it exists.  If ``all`` is
        ``True``, a list containing all square roots of ``self``
        (of length zero, one or two) is returned instead.

        If ``extend`` is ``True``, a square root is chosen in an
        extension field if necessary.  If ``extend`` is ``False``, a
        ValueError is raised if the element is not a square in the
        base field.

        .. WARNING::

           The ``extend`` option is not implemented (yet).

        EXAMPLES::

            sage: F = FiniteField(7^2, 'a', impl='pari_ffelt')
            sage: F(2).sqrt()
            4
            sage: F(3).sqrt()
            5*a + 1
            sage: F(3).sqrt()**2
            3
            sage: F(4).sqrt(all=True)
            [2, 5]

            sage: K = FiniteField(7^3, 'alpha', impl='pari_ffelt')
            sage: K(3).sqrt()
            Traceback (most recent call last):
            ...
            ValueError: element is not a square
            sage: K(3).sqrt(all=True)
            []

            sage: K.<a> = GF(3^17, impl='pari_ffelt')
            sage: (a^3 - a - 1).sqrt()
            a^16 + 2*a^15 + a^13 + 2*a^12 + a^10 + 2*a^9 + 2*a^8 + a^7 + a^6 + 2*a^5 + a^4 + 2*a^2 + 2*a + 2
        """
        if extend:
            raise NotImplementedError
        cdef GEN s
        cdef FiniteFieldElement_pari_ffelt x, mx
        pari_catch_sig_on()
        if FF_issquareall(self.val, &s):
            x = self._new()
            x.construct(s)
            if not all:
                return x
            elif gequal0(x.val) or self._parent.characteristic() == 2:
                return [x]
            else:
                pari_catch_sig_on()
                mx = self._new()
                mx.construct(FF_neg_i(x.val))
                return [x, mx]
        else:
            pari_catch_sig_off()
            if all:
                return []
            else:
                raise ValueError("element is not a square")

    def log(FiniteFieldElement_pari_ffelt self, object base):
        """
        Return a discrete logarithm of ``self`` with respect to the
        given base.

        INPUT:

        - ``base`` -- non-zero field element

        OUTPUT:

        An integer `x` such that ``self`` equals ``base`` raised to
        the power `x`.  If no such `x` exists, a ``ValueError`` is
        raised.

        EXAMPLES::

            sage: F.<g> = FiniteField(2^10, impl='pari_ffelt')
            sage: b = g; a = g^37
            sage: a.log(b)
            37
            sage: b^37; a
            g^8 + g^7 + g^4 + g + 1
            g^8 + g^7 + g^4 + g + 1

        ::

            sage: F.<a> = FiniteField(5^2, impl='pari_ffelt')
            sage: F(-1).log(F(2))
            2
            sage: F(1).log(a)
            0

        Some cases where the logarithm is not defined or does not exist::

            sage: F.<a> = GF(3^10, impl='pari_ffelt')
            sage: a.log(-1)
            Traceback (most recent call last):
            ...
            ArithmeticError: element a does not lie in group generated by 2
            sage: a.log(0)
            Traceback (most recent call last):
            ...
            ArithmeticError: discrete logarithm with base 0 is not defined
            sage: F(0).log(1)
            Traceback (most recent call last):
            ...
            ArithmeticError: discrete logarithm of 0 is not defined
        """
        base = self._parent(base)
        if self.is_zero():
            raise ArithmeticError("discrete logarithm of 0 is not defined")
        if base.is_zero():
            raise ArithmeticError("discrete logarithm with base 0 is not defined")

        # Compute the orders of self and base to check whether self
        # actually lies in the cyclic group generated by base. PARI
        # requires that this is the case.
        # We also have to specify the order of the base anyway
        # because PARI assumes by default that this element generates
        # the multiplicative group.
        cdef GEN x, base_order, self_order
        pari_catch_sig_on()
        base_order = FF_order((<FiniteFieldElement_pari_ffelt>base).val, NULL)
        self_order = FF_order(self.val, NULL)
        if not dvdii(base_order, self_order):
            # self_order does not divide base_order
            pari.clear_stack()
            raise ArithmeticError("element %s does not lie in group generated by %s"%(self, base))
        x = FF_log(self.val, (<FiniteFieldElement_pari_ffelt>base).val, base_order)
        return Integer(pari.new_gen(x))

    def multiplicative_order(FiniteFieldElement_pari_ffelt self):
        """
        Returns the order of ``self`` in the multiplicative group.

        EXAMPLE::

            sage: a = FiniteField(5^3, 'a', impl='pari_ffelt').0
            sage: a.multiplicative_order()
            124
            sage: a**124
            1
        """
        if self.is_zero():
            raise ArithmeticError("Multiplicative order of 0 not defined.")
        cdef GEN order
        pari_catch_sig_on()
        order = FF_order(self.val, NULL)
        return Integer(pari.new_gen(order))

    def lift(FiniteFieldElement_pari_ffelt self):
        """
        If ``self`` is an element of the prime field, return a lift of
        this element to an integer.

        EXAMPLE::

            sage: k = FiniteField(next_prime(10^10)^2, 'u', impl='pari_ffelt')
            sage: a = k(17)/k(19)
            sage: b = a.lift(); b
            7894736858
            sage: b.parent()
            Integer Ring
        """
        if FF_equal0(self.val):
            return Integer(0)
        f = self.polynomial()
        if f.degree() == 0:
            return f.constant_coefficient().lift()
        else:
            raise ValueError("element is not in the prime field")

    def _integer_(self, ZZ=None):
        """
        Lift to a Sage integer, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: b = k(2)
            sage: b._integer_()
            2
            sage: a._integer_()
            Traceback (most recent call last):
            ...
            ValueError: element is not in the prime field
        """
        return self.lift()

    def __int__(self):
        """
        Lift to a python int, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: b = k(2)
            sage: int(b)
            2
            sage: int(a)
            Traceback (most recent call last):
            ...
            ValueError: element is not in the prime field
        """
        return int(self.lift())

    def __long__(self):
        """
        Lift to a python long, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: b = k(2)
            sage: long(b)
            2L
        """
        return long(self.lift())

    def __float__(self):
        """
        Lift to a python float, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: b = k(2)
            sage: float(b)
            2.0
        """
        return float(self.lift())

    def _pari_(self, var=None):
        """
        Return a PARI object representing ``self``.

        INPUT:

        - var -- ignored

        EXAMPLE::

            sage: k.<a> = FiniteField(3^3, impl='pari_ffelt')
            sage: b = a**2 + 2*a + 1
            sage: b._pari_()
            a^2 + 2*a + 1
        """
        pari_catch_sig_on()
        return pari.new_gen(self.val)

    def _pari_init_(self):
        """
        Return a string representing ``self`` in PARI.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: a._pari_init_()
            'subst(a+3*a,a,ffgen(Mod(1, 3)*x^17 + Mod(2, 3)*x + Mod(1, 3),a))'
            sage: k(1)._pari_init_()
            'subst(1+3*a,a,ffgen(Mod(1, 3)*x^17 + Mod(2, 3)*x + Mod(1, 3),a))'

        This is used for conversion to GP. The element is displayed
        as "a" but has correct arithmetic::

            sage: gp(a)
            a
            sage: gp(a).type()
            t_FFELT
            sage: gp(a)^100
            2*a^16 + 2*a^15 + a^4 + a + 1
            sage: gp(a^100)
            2*a^16 + 2*a^15 + a^4 + a + 1
            sage: gp(k(0))
            0
            sage: gp(k(0)).type()
            t_FFELT
        """
        ffgen = "ffgen(%s,a)" % self._parent.modulus()._pari_init_()
        # Add this "zero" to ensure that the polynomial is not constant
        zero = "%s*a" % self._parent.characteristic()
        return "subst(%s+%s,a,%s)" % (self, zero, ffgen)

    def _magma_init_(self, magma):
        """
        Return a string representing ``self`` in Magma.

        EXAMPLE::

            sage: GF(7)(3)._magma_init_(magma)            # optional - magma
            'GF(7)!3'
        """
        k = self._parent
        km = magma(k)
        return str(self).replace(k.variable_name(), km.gen(1).name())

    def _gap_init_(self):
        r"""
        Return the a string representing ``self`` in GAP.

        .. NOTE::

           The order of the parent field must be `\leq 65536`.  This
           function can be slow since elements of non-prime finite
           fields are represented in GAP as powers of a generator for
           the multiplicative group, so a discrete logarithm must be
           computed.

        EXAMPLE::

            sage: F = FiniteField(2^3, 'a', impl='pari_ffelt')
            sage: a = F.multiplicative_generator()
            sage: gap(a) # indirect doctest
            Z(2^3)
            sage: b = F.multiplicative_generator()
            sage: a = b^3
            sage: gap(a)
            Z(2^3)^3
            sage: gap(a^3)
            Z(2^3)^2

        You can specify the instance of the Gap interpreter that is used::

            sage: F = FiniteField(next_prime(200)^2, 'a', impl='pari_ffelt')
            sage: a = F.multiplicative_generator ()
            sage: a._gap_ (gap)
            Z(211^2)
            sage: (a^20)._gap_(gap)
            Z(211^2)^20

        Gap only supports relatively small finite fields::

            sage: F = FiniteField(next_prime(1000)^2, 'a', impl='pari_ffelt')
            sage: a = F.multiplicative_generator ()
            sage: gap._coerce_(a)
            Traceback (most recent call last):
            ...
            TypeError: order must be at most 65536
        """
        F = self._parent
        if F.order() > 65536:
            raise TypeError("order must be at most 65536")

        if self == 0:
            return '0*Z(%s)'%F.order()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = self.log(g)
        return 'Z(%s)^%s'%(F.order(), n)


def unpickle_FiniteFieldElement_pari_ffelt(parent, elem):
    """
    EXAMPLE::

        sage: k.<a> = GF(2^20, impl='pari_ffelt')
        sage: e = k.random_element()
        sage: f = loads(dumps(e)) # indirect doctest
        sage: e == f
        True
    """
    return parent(elem)
