"""
Finite field elements implemented via PARI's FFELT type

AUTHORS:

- Peter Bruin (June 2013): initial version, based on
  element_ext_pari.py by William Stein et al. and
  element_ntl_gf2e.pyx by Martin Albrecht.
"""
# ****************************************************************************
#      Copyright (C) 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from cypari2.paridecl cimport *
from cypari2.paripriv cimport *
from sage.libs.pari.convert_gmp cimport _new_GEN_from_mpz_t
from cypari2.stack cimport new_gen, new_gen_noclear, clear_stack
from cypari2.gen cimport Gen as pari_gen, objtogen

from .element_base cimport FinitePolyExtElement
from .integer_mod import IntegerMod_abstract

import sage.rings.integer
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial
from sage.rings.rational import Rational
from sage.structure.richcmp cimport rich_to_bool

try:
    from sage.modules.free_module_element import FreeModuleElement
except ImportError:
    FreeModuleElement = ()

from sage.interfaces.abc import GapElement

cdef GEN _INT_to_FFELT(GEN g, GEN x) except NULL:
    """
    Convert the t_INT `x` to an element of the field of definition of
    the t_FFELT `g`.

    This function must be called within ``sig_on()`` ... ``sig_off()``.

    TESTS:

    Converting large integers to finite field elements does not lead
    to overflow errors (see :trac:`16807`)::

        sage: p = previous_prime(2^64)
        sage: F.<x> = GF(p^2)
        sage: x * 2^63
        9223372036854775808*x

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
            f[2] = itou(x)
        else:
            sig_off()
            raise TypeError("unknown PARI finite field type")
        result = cgetg(5, t_FFELT)
        result[1] = t
        set_gel(result, 2, f)
        set_gel(result, 3, gel(g, 3))  # modulus
        set_gel(result, 4, p)
        return result

cdef class FiniteFieldElement_pari_ffelt(FinitePolyExtElement):
    """
    An element of a finite field implemented using PARI.

    EXAMPLES::

        sage: K = FiniteField(10007^10, 'a', impl='pari_ffelt')
        sage: a = K.gen(); a
        a
        sage: type(a)
        <class 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt'>

    TESTS::

        sage: n = 63
        sage: m = 3
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

    ::

        sage: k = FiniteField(3^4, 'a', impl='pari_ffelt')
        sage: b = k(5) # indirect doctest
        sage: b.parent()
        Finite Field in a of size 3^4
        sage: a = k.gen()
        sage: k(a + 2)
        a + 2

    Univariate polynomials coerce into finite fields by evaluating
    the polynomial at the field's generator::

        sage: R.<x> = QQ[]
        sage: k.<a> = FiniteField(5^2, 'a', impl='pari_ffelt')
        sage: k(R(2/3))
        4
        sage: k(x^2)
        a + 3

        sage: R.<x> = GF(5)[]
        sage: k(x^3-2*x+1)
        2*a + 4

        sage: x = polygen(QQ)
        sage: k(x^25)
        a

        sage: Q.<q> = FiniteField(5^7, 'q', impl='pari_ffelt')
        sage: L = GF(5)
        sage: LL.<xx> = L[]
        sage: Q(xx^2 + 2*xx + 4)
        q^2 + 2*q + 4

        sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
        sage: k.polynomial()
        t^11 + 2*t^2 + 1
        sage: P = k.polynomial_ring()
        sage: k(P.0^11)
        t^2 + 2

    An element can be specified by its vector of coordinates with
    respect to the basis consisting of powers of the generator:

        sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
        sage: V = k.vector_space(map=False)
        sage: V
        Vector space of dimension 11 over Finite Field of size 3
        sage: v = V([0,1,2,0,1,2,0,1,2,0,1])
        sage: k(v)
        t^10 + 2*t^8 + t^7 + 2*t^5 + t^4 + 2*t^2 + t

    Multivariate polynomials only coerce if constant::

        sage: k = FiniteField(5^2, 'a', impl='pari_ffelt')
        sage: R = k['x,y,z']; R
        Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 5^2
        sage: k(R(2))
        2
        sage: R = QQ['x,y,z']
        sage: k(R(1/5))
        Traceback (most recent call last):
        ...
        ZeroDivisionError: inverse of Mod(0, 5) does not exist

    Gap elements can also be coerced into finite fields::

        sage: F = FiniteField(2^3, 'a', impl='pari_ffelt')
        sage: a = F.multiplicative_generator(); a
        a
        sage: b = gap(a^3); b
        Z(2^3)^3
        sage: F(b)
        a + 1
        sage: a^3
        a + 1

        sage: a = GF(13)(gap('0*Z(13)')); a
        0
        sage: a.parent()
        Finite Field of size 13

        sage: F = FiniteField(2^4, 'a', impl='pari_ffelt')
        sage: F(gap('Z(16)^3'))
        a^3
        sage: F(gap('Z(16)^2'))
        a^2

    You can also call a finite extension field with a string
    to produce an element of that field, like this::

        sage: k = GF(2^8, 'a')
        sage: k('a^200')
        a^4 + a^3 + a^2

    This is especially useful for conversion from Singular etc.

    TESTS::

        sage: k = FiniteField(3^2, 'a', impl='pari_ffelt')
        sage: a = k(11); a
        2
        sage: a.parent()
        Finite Field in a of size 3^2
        sage: V = k.vector_space(map=False); v = V((1,2))
        sage: k(v)
        2*a + 1

    We create elements using a list and verify that :trac:`10486` has
    been fixed::

        sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
        sage: x = k([1,0,2,1]); x
        t^3 + 2*t^2 + 1
        sage: x + x + x
        0
        sage: pari(x)
        t^3 + 2*t^2 + 1

    If the list is longer than the degree, we just get the result
    modulo the modulus::

        sage: from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt
        sage: R.<a> = PolynomialRing(GF(5))
        sage: k = FiniteField_pari_ffelt(5, a^2 - 2, 't')
        sage: x = k([0,0,0,1]); x
        2*t
        sage: pari(x)
        2*t

    When initializing from a list, the elements are first coerced
    to the prime field (:trac:`11685`)::

        sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
        sage: k([ 0, 1/2 ])
        2*t
        sage: k([ 0, 1/2, 0, 0, 0, 0, 0, 0, 0, -1, 0 ])
        2*t^9 + 2*t
        sage: k([ k(0), k(1) ])
        t
        sage: k([ GF(3)(2), GF(3^5,'u')(1) ])
        t + 2
        sage: R.<x> = PolynomialRing(k)
        sage: k([ x/x ])
        1
        sage: k([ R(-1), x/x ])
        t + 2
        sage: k([ R(-1), R(0), 0 ])
        2

    Check that zeros are created correctly (:trac:`11685`)::

        sage: K = FiniteField(3^11, 't', impl='pari_ffelt'); a = K.0
        sage: v = 0; pari(K(v))
        0
        sage: v = Mod(0,3); pari(K(v))
        0
        sage: v = pari(0); pari(K(v))
        0
        sage: v = pari("Mod(0,3)"); pari(K(v))
        0
        sage: v = []; pari(K(v))
        0
        sage: v = [0]; pari(K(v))
        0
        sage: v = [0,0]; pari(K(v))
        0
        sage: v = pari("Pol(0)"); pari(K(v))
        0
        sage: v = pari("Mod(0, %s)"%K.modulus()); pari(K(v))
        0
        sage: v = pari("Mod(Pol(0), %s)"%K.modulus()); pari(K(v))
        0
        sage: v = K(1) - K(1); pari(K(v))
        0
        sage: v = K([1]) - K([1]); pari(K(v))
        0
        sage: v = a - a; pari(K(v))
        0
        sage: v = K(1)*0; pari(K(v))
        0
        sage: v = K([1])*K([0]); pari(K(v))
        0
        sage: v = a*0; pari(K(v))
        0
    """

    def __init__(self, parent, x):
        """
        Initialise ``self`` with the given ``parent`` and value
        converted from ``x``.

        This is called when constructing elements from Python.

        TESTS::

            sage: from sage.rings.finite_rings.element_pari_ffelt import FiniteFieldElement_pari_ffelt
            sage: K = FiniteField(101^2, 'a', impl='pari_ffelt')
            sage: x = FiniteFieldElement_pari_ffelt(K, 'a + 1')
            sage: x
            a + 1
        """
        # FinitePolyExtElement.__init__(self, parent)
        self._parent = parent
        self.construct_from(x)

    def __dealloc__(self):
        """
        Cython destructor.
        """
        if self.val is not NULL:
            gunclone_deep(self.val)

    cdef FiniteFieldElement_pari_ffelt _new(self):
        """
        Create an empty element with the same parent as ``self``.
        """
        cdef FiniteFieldElement_pari_ffelt x
        x = FiniteFieldElement_pari_ffelt.__new__(FiniteFieldElement_pari_ffelt)
        x._parent = self._parent
        return x

    cdef void construct(self, GEN g):
        """
        Initialise ``self`` to the FFELT ``g``, reset the PARI stack,
        and call sig_off().

        This should be called exactly once on every instance.
        """
        self.val = gcloneref(g)
        clear_stack()

    cdef int construct_from(self, x) except -1:
        """
        Initialise ``self`` to an FFELT constructed from the Sage
        object `x`.

        TESTS:

        Conversion of elements of the underlying vector space works in
        large characteristic (see :trac:`21186`)::

            sage: p = 13189065031705623239
            sage: Fq = FiniteField(p^3, "a")
            sage: Fq_X = PolynomialRing(Fq, "x")
            sage: pol = Fq_X("x^9 + 13189065031705622723*x^7 + 13189065031705622723*x^6 + 9288*x^5 + 18576*x^4 + 13189065031705590731*x^3 + 13189065031705497851*x^2 + 13189065031705497851*x + 13189065031705581443")
            sage: R = [r[0] for r in pol.roots()]
            sage: prod(Fq_X.gen() - r for r in R) == pol
            True

        """
        cdef GEN f, g, result, x_GEN
        cdef long i, n, t
        cdef Integer xi

        if isinstance(x, FiniteFieldElement_pari_ffelt):
            if self._parent is (<FiniteFieldElement_pari_ffelt>x)._parent:
                sig_on()
                self.construct((<FiniteFieldElement_pari_ffelt>x).val)
            else:
                raise TypeError("no coercion defined")

        elif isinstance(x, Integer):
            g = (<pari_gen>self._parent._gen_pari).g
            sig_on()
            x_GEN = _new_GEN_from_mpz_t((<Integer>x).value)
            self.construct(_INT_to_FFELT(g, x_GEN))

        elif isinstance(x, int):
            g = (<pari_gen>self._parent._gen_pari).g
            x = objtogen(x)
            sig_on()
            x_GEN = (<pari_gen>x).g
            self.construct(_INT_to_FFELT(g, x_GEN))

        elif isinstance(x, IntegerMod_abstract):
            if self._parent.characteristic().divides(x.modulus()):
                g = (<pari_gen>self._parent._gen_pari).g
                sig_on()
                x_GEN = _new_GEN_from_mpz_t(Integer(x).value)
                self.construct(_INT_to_FFELT(g, x_GEN))
            else:
                raise TypeError("no coercion defined")

        elif x is None:
            g = (<pari_gen>self._parent._gen_pari).g
            sig_on()
            self.construct(FF_zero(g))

        elif isinstance(x, pari_gen):
            g = (<pari_gen>self._parent._gen_pari).g
            x_GEN = (<pari_gen>x).g

            sig_on()
            if gequal0(x_GEN):
                self.construct(FF_zero(g))
                return 0
            elif gequal1(x_GEN):
                self.construct(FF_1(g))
                return 0

            t = typ(x_GEN)
            if t == t_FFELT:
                if FF_samefield(x_GEN, g):
                    self.construct(x_GEN)
                    return 0
            elif t == t_INT:
                self.construct(_INT_to_FFELT(g, x_GEN))
                return 0
            elif t == t_INTMOD:
                if gequal0(modii(gel(x_GEN, 1), FF_p_i(g))):
                    self.construct(_INT_to_FFELT(g, gel(x_GEN, 2)))
                    return 0
            elif t == t_FRAC:
                if not gequal0(modii(gel(x_GEN, 2), FF_p_i(g))):
                    elt = FF_div(_INT_to_FFELT(g, gel(x_GEN, 1)),
                                 _INT_to_FFELT(g, gel(x_GEN, 2)))
                    self.construct(elt)
                    return 0
            sig_off()
            raise TypeError(f"unable to convert PARI {x.type()} to finite field element")

        elif (isinstance(x, FreeModuleElement)
              and x.parent() is self._parent.vector_space(map=False)):
            g = (<pari_gen>self._parent._gen_pari).g
            t = g[1]  # codeword: t_FF_FpXQ, t_FF_Flxq, t_FF_F2xq
            n = len(x)
            while n > 0 and x[n - 1] == 0:
                n -= 1
            sig_on()
            if n == 0:
                self.construct(FF_zero(g))
                return 0
            if t == t_FF_FpXQ:
                f = cgetg(n + 2, t_POL)
                set_gel(f, 1, gmael(g, 2, 1))
                for i in range(n):
                    xi = Integer(x[i])
                    set_gel(f, i + 2, _new_GEN_from_mpz_t(xi.value))
            elif t == t_FF_Flxq or t == t_FF_F2xq:
                f = cgetg(n + 2, t_VECSMALL)
                set_gel(f, 1, gmael(g, 2, 1))
                for i in range(n):
                    set_uel(f, i + 2, x[i])
                if t == t_FF_F2xq:
                    f = Flx_to_F2x(f)
            else:
                sig_off()
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
            n = len(x)
            if n == 0:
                self.construct_from(None)
            elif n == 1:
                Fp = self._parent.base_ring()
                self.construct_from(Fp(x[0]))
            elif n == self._parent.degree():
                self.construct_from(self._parent.vector_space(map=False)(x))
            else:
                Fp = self._parent.base_ring()
                self.construct_from(self._parent.polynomial_ring()([Fp(y) for y in x]))

        elif isinstance(x, str):
            self.construct_from(self._parent.polynomial_ring()(x))

        elif isinstance(x, GapElement):
            try:
                from sage.libs.gap.libgap import libgap
                self.construct_from(libgap(x).sage(ring=self._parent))
            except (ValueError, IndexError, TypeError):
                raise TypeError("no coercion defined")

        elif isinstance(x, sage.libs.gap.element.GapElement_FiniteField):
            try:
                self.construct_from(x.sage(ring=self._parent))
            except (ValueError, IndexError, TypeError):
                raise TypeError("no coercion defined")

        else:
            raise TypeError("no coercion defined")

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: k.<c> = GF(3^17, impl='pari_ffelt')
            sage: c^20  # indirect doctest
            c^4 + 2*c^3
        """
        return str(new_gen_noclear(self.val))

    def __hash__(self):
        """
        Return the hash of ``self``.  This is by definition equal to
        the hash of ``self.polynomial()``.

        EXAMPLES::

            sage: k.<a> = GF(3^15, impl='pari_ffelt')
            sage: R = GF(3)['a']; aa = R.gen()
            sage: hash(a^2 + 1) == hash(aa^2 + 1)
            True
        """
        return hash(self.polynomial())

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: K.<a> = FiniteField(10007^10, impl='pari_ffelt')
            sage: loads(a.dumps()) == a
            True
        """
        return unpickle_FiniteFieldElement_pari_ffelt, (self._parent, str(self))

    def __copy__(self):
        """
        TESTS::

            sage: k.<a> = FiniteField(3^3, impl='pari_ffelt')
            sage: a
            a
            sage: b = copy(a); b
            a
            sage: a is b
            True
        """
        # immutable
        return self

    def __deepcopy__(self, memo):
        """
        TESTS::

            sage: k.<a> = FiniteField(3^3, impl='pari_ffelt')
            sage: a
            a
            sage: b = deepcopy(a); b
            a
            sage: a is b
            True
        """
        # immutable
        return self

    cpdef _richcmp_(self, other, int op):
        """
        Comparison of finite field elements.

        .. NOTE::

            Finite fields are unordered.  However, for the purpose of
            this function, we adopt the lexicographic ordering on the
            representing polynomials.

        EXAMPLES::

            sage: k.<a> = GF(2^20, impl='pari_ffelt')
            sage: e = k.random_element()
            sage: f = loads(dumps(e))
            sage: e is f
            False
            sage: e == f
            True
            sage: e != (e + 1)
            True

        ::

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
        cdef int r
        sig_on()
        r = cmp_universal(self.val, (<FiniteFieldElement_pari_ffelt>other).val)
        sig_off()
        return rich_to_bool(op, r)

    cpdef _add_(self, right):
        """
        Addition.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: a + a^2 # indirect doctest
            a^2 + a
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(FF_add((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    cpdef _sub_(self, right):
        """
        Subtraction.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: a - a # indirect doctest
            0
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(FF_sub((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    cpdef _mul_(self, right):
        """
        Multiplication.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: (a^12 + 1)*(a^15 - 1) # indirect doctest
            a^15 + 2*a^12 + a^11 + 2*a^10 + 2
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(FF_mul((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    cpdef _div_(self, right):
        """
        Division.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: (a - 1) / (a + 1) # indirect doctest
            2*a^16 + a^15 + 2*a^14 + a^13 + 2*a^12 + a^11 + 2*a^10 + a^9 + 2*a^8 + a^7 + 2*a^6 + a^5 + 2*a^4 + a^3 + 2*a^2 + a + 1
        """
        if FF_equal0((<FiniteFieldElement_pari_ffelt>right).val):
            raise ZeroDivisionError
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(FF_div((<FiniteFieldElement_pari_ffelt>self).val,
                           (<FiniteFieldElement_pari_ffelt>right).val))
        return x

    def is_zero(self):
        """
        Return ``True`` if ``self`` equals 0.

        EXAMPLES::

            sage: F.<a> = FiniteField(5^3, impl='pari_ffelt')
            sage: a.is_zero()
            False
            sage: (a - a).is_zero()
            True
        """
        return bool(FF_equal0(self.val))

    def is_one(self):
        """
        Return ``True`` if ``self`` equals 1.

        EXAMPLES::

            sage: F.<a> = FiniteField(5^3, impl='pari_ffelt')
            sage: a.is_one()
            False
            sage: (a/a).is_one()
            True
        """
        return bool(FF_equal1(self.val))

    def is_unit(self):
        """
        Return ``True`` if ``self`` is non-zero.

        EXAMPLES::

            sage: F.<a> = FiniteField(5^3, impl='pari_ffelt')
            sage: a.is_unit()
            True
        """
        return not bool(FF_equal0(self.val))

    __bool__ = is_unit

    def __pos__(self):
        """
        Unitary positive operator...

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: +a
            a
        """
        return self

    def __neg__(self):
        """
        Negation.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: -a
            2*a
        """
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(FF_neg_i((<FiniteFieldElement_pari_ffelt>self).val))
        return x

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

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
        sig_on()
        x.construct(FF_inv((<FiniteFieldElement_pari_ffelt>self).val))
        return x

    def __pow__(FiniteFieldElement_pari_ffelt self, exp, other):
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

        The exponent is converted to an integer (see :trac:`16540`)::

            sage: q = 11^23
            sage: F.<a> = FiniteField(q)
            sage: a^Mod(1, q - 1)
            a

        .. WARNING::

            For efficiency reasons, we do not verify that the
            exponentiation is well defined before converting the
            exponent to an integer.  This means that ``a^Mod(1, n)``
            returns `a` even if `n` is not a multiple of the
            multiplicative order of `a`.

        """
        if exp == 0:
            return self._parent.one()
        if exp < 0 and FF_equal0(self.val):
            raise ZeroDivisionError
        exp = Integer(exp).__pari__()
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(FF_pow(self.val, (<pari_gen>exp).g))
        return x

    def pth_power(FiniteFieldElement_pari_ffelt self, int k=1):
        r"""
        Return the `(p^k)^{th}` power of ``self``, where `p` is the
        characteristic of the field.

        INPUT:

        - ``k`` -- integer (default: 1); must fit in a C ``int``

        Note that if `k` is negative, then this computes the appropriate root.

        TESTS::

            sage: F.<a> = GF(13^64, impl='pari_ffelt'); F
            Finite Field in a of size 13^64
            sage: x = F.random_element()
            sage: x.pth_power(0) == x
            True
            sage: x.pth_power(1) == x**13
            True
            sage: x.pth_power(2) == x**(13**2)
            True
            sage: x.pth_power(-1)**13 == x
            True

            sage: F.<a> = GF(127^16, impl='pari_ffelt'); F
            Finite Field in a of size 127^16
            sage: x = F.random_element()
            sage: x.pth_power(0) == x
            True
            sage: x.pth_power(1) == x**127
            True
            sage: x.pth_power(2) == x**(127**2)
            True
            sage: x.pth_power(-1)**127 == x
            True
        """
        cdef int n = int(self._parent.degree())
        if k % n == 0:
            return self
        cdef Integer p = self._parent.characteristic()
        if k == 1 and (p < 100 or p.bit_length()**2 < n):
            # For extremely small primes or very large extension degrees,
            # exponentiation is faster.
            return self**p
        # Otherwise use PARI field morphism (evaluation of a Fp polynomial
        # at the image of the generator).
        f = self._parent._pari_frobenius(k)
        cdef FiniteFieldElement_pari_ffelt x = self._new()
        sig_on()
        x.construct(ffmap((<pari_gen>f).g, self.val))
        return x

    frobenius = pth_power

    def polynomial(self, name=None):
        """
        Return the unique representative of ``self`` as a polynomial
        over the prime field whose degree is less than the degree of
        the finite field over its prime field.

        INPUT:

        - ``name`` -- (optional) variable name

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
            sage: (a**2 + 1).polynomial('beta')
            beta^2 + 1
            sage: (a**2 + 1).polynomial().parent()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
            sage: (a**2 + 1).polynomial('beta').parent()
            Univariate Polynomial Ring in beta over Finite Field of size 3
        """
        sig_on()
        pol = new_gen(FF_to_FpXQ(self.val))
        return self._parent.polynomial_ring(name)(pol)

    def minpoly(self, var='x'):
        """
        Return the minimal polynomial of ``self``.

        INPUT:

        - ``var`` -- string (default: 'x'): variable name to use.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(FiniteField(3))
            sage: F.<a> = FiniteField(3^2, modulus=x^2 + 1, impl='pari_ffelt')
            sage: a.minpoly('y')
            y^2 + 1
        """
        sig_on()
        pol = new_gen(FF_minpoly(self.val))
        return self._parent.polynomial_ring(var)(pol)

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of ``self``.

        INPUT:

        - ``var`` -- string (default: 'x'): variable name to use.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(FiniteField(3))
            sage: F.<a> = FiniteField(3^2, modulus=x^2 + 1, impl='pari_ffelt')
            sage: a.charpoly('y')
            y^2 + 1
        """
        sig_on()
        pol = new_gen(FF_charpoly(self.val))
        return self._parent.polynomial_ring(var)(pol)

    def is_square(self):
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
        sig_on()
        i = FF_issquare(self.val)
        sig_off()
        return bool(i)

    def sqrt(self, extend=False, all=False):
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
            sage: F(3).sqrt() in (2*F.gen() + 6, 5*F.gen() + 1)
            True
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
        sig_on()
        if FF_issquareall(self.val, &s):
            x = self._new()
            x.construct(s)
            if not all:
                return x
            elif gequal0(x.val) or self._parent.characteristic() == 2:
                return [x]
            else:
                sig_on()
                mx = self._new()
                mx.construct(FF_neg_i(x.val))
                return [x, mx]
        else:
            sig_off()
            if all:
                return []
            else:
                raise ValueError("element is not a square")

    def log(self, base):
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
        sig_on()
        base_order = FF_order((<FiniteFieldElement_pari_ffelt>base).val, NULL)
        self_order = FF_order(self.val, NULL)
        if not dvdii(base_order, self_order):
            # self_order does not divide base_order
            clear_stack()
            raise ArithmeticError("element %s does not lie in group generated by %s"%(self, base))
        x = FF_log(self.val, (<FiniteFieldElement_pari_ffelt>base).val, base_order)
        return Integer(new_gen(x))

    def multiplicative_order(self):
        """
        Returns the order of ``self`` in the multiplicative group.

        EXAMPLES::

            sage: a = FiniteField(5^3, 'a', impl='pari_ffelt').0
            sage: a.multiplicative_order()
            124
            sage: a**124
            1
        """
        if self.is_zero():
            raise ArithmeticError("Multiplicative order of 0 not defined.")
        cdef GEN order
        sig_on()
        order = FF_order(self.val, NULL)
        return Integer(new_gen(order))

    def lift(self):
        """
        If ``self`` is an element of the prime field, return a lift of
        this element to an integer.

        EXAMPLES::

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

        EXAMPLES::

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

        EXAMPLES::

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

    def __float__(self):
        """
        Lift to a python float, if possible.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_ffelt')
            sage: b = k(2)
            sage: float(b)
            2.0
        """
        return float(self.lift())

    def __pari__(self, var=None):
        """
        Return a PARI object representing ``self``.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^3, impl='pari_ffelt')
            sage: b = a**2 + 2*a + 1
            sage: b.__pari__()
            a^2 + 2*a + 1
        """
        return new_gen_noclear(self.val)

    def _pari_init_(self):
        """
        Return a string representing ``self`` in PARI.

        EXAMPLES::

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

        EXAMPLES::

            sage: GF(7)(3)._magma_init_(magma)            # optional - magma
            'GF(7)!3'
        """
        k = self._parent
        km = magma(k)
        return str(self).replace(k.variable_name(), km.gen(1).name())

    def _gap_init_(self):
        r"""
        Return the string representing ``self`` in GAP.

        .. NOTE::

            The order of the parent field must be `\leq 65536`.  This
            function can be slow since elements of non-prime finite
            fields are represented in GAP as powers of a generator for
            the multiplicative group, so a discrete logarithm must be
            computed.

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: F = FiniteField(2^3, 'aa', impl='pari_ffelt')
            sage: aa = F.multiplicative_generator()
            sage: gap(aa)  # indirect doctest
            Z(2^3)
            sage: b = F.multiplicative_generator()
            sage: a = b^3
            sage: gap(a)
            Z(2^3)^3
            sage: gap(a^3)
            Z(2^3)^2
            sage: F(gap('Z(8)^3'))
            aa + 1
            sage: F(libgap.Z(8)^3)
            aa + 1

        You can specify the instance of the Gap interpreter that is used::

            sage: # needs sage.libs.gap
            sage: F = FiniteField(next_prime(200)^2, 'a', impl='pari_ffelt')
            sage: a = F.multiplicative_generator()
            sage: a._gap_ (gap)
            Z(211^2)
            sage: (a^20)._gap_(gap)
            Z(211^2)^20

        Gap only supports relatively small finite fields::

            sage: # needs sage.libs.gap
            sage: F = FiniteField(next_prime(1000)^2, 'a', impl='pari_ffelt')
            sage: a = F.multiplicative_generator()
            sage: a._gap_init_()
            Traceback (most recent call last):
            ...
            TypeError: order must be at most 65536
            sage: gap.coerce(a)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field in a of size 1009^2 to Gap
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
    EXAMPLES::

        sage: k.<a> = GF(2^20, impl='pari_ffelt')
        sage: e = k.random_element()
        sage: f = loads(dumps(e)) # indirect doctest
        sage: e == f
        True
    """
    return parent(elem)
