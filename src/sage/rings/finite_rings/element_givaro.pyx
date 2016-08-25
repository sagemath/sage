r"""
Givaro Field Elements

Sage includes the Givaro finite field library, for highly optimized
arithmetic in finite fields.

.. NOTE::

    The arithmetic is performed by the Givaro C++ library which uses
    Zech logs internally to represent finite field elements. This
    implementation is the default finite extension field implementation
    in Sage for the cardinality less than `2^{16}`, as it is a lot
    faster than the PARI implementation. Some functionality in this
    class however is implemented using PARI.

EXAMPLES::

    sage: k = GF(5); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>
    sage: k = GF(5^2,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
    sage: k = GF(2^16,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>
    sage: k = GF(3^16,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>

    sage: n = previous_prime_power(2^16 - 1)
    sage: while is_prime(n):
    ...    n = previous_prime_power(n)
    sage: factor(n)
    251^2
    sage: k = GF(n,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>

AUTHORS:

- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-06-05)
- William Stein (2006-12-07): editing, lots of docs, etc.
- Robert Bradshaw (2007-05-23): is_square/sqrt, pow.

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "cysignals/signals.pxi"
include "sage/libs/ntl/decl.pxi"
from sage.libs.pari.paridecl cimport *

from sage.misc.randstate cimport randstate, current_randstate
from sage.rings.finite_rings.finite_field_base cimport FiniteField
from sage.rings.ring cimport Ring
from element_ext_pari import FiniteField_ext_pariElement
from element_pari_ffelt import FiniteFieldElement_pari_ffelt
from sage.structure.sage_object cimport SageObject
import operator
import sage.arith.all
import finite_field_constructor as finite_field

import sage.interfaces.gap
from sage.libs.pari.all import pari
from sage.libs.pari.gen cimport gen

from sage.structure.parent cimport Parent

from sage.misc.superseded import deprecated_function_alias

cdef object is_IntegerMod
cdef object Integer
cdef object Rational
cdef object ConwayPolynomials
cdef object conway_polynomial
cdef object MPolynomial
cdef object Polynomial
cdef object FreeModuleElement

cdef void late_import():
    """
    Late import of modules
    """
    global is_IntegerMod, \
           Integer, \
           Rational, \
           ConwayPolynomials, \
           conway_polynomial, \
           MPolynomial, \
           Polynomial, \
           FreeModuleElement

    if is_IntegerMod is not None:
        return

    import sage.rings.finite_rings.integer_mod
    is_IntegerMod = sage.rings.finite_rings.integer_mod.is_IntegerMod

    import sage.rings.integer
    Integer = sage.rings.integer.Integer

    import sage.rings.rational
    Rational = sage.rings.rational.Rational

    import sage.databases.conway
    ConwayPolynomials = sage.databases.conway.ConwayPolynomials

    import sage.rings.finite_rings.finite_field_constructor
    conway_polynomial = sage.rings.finite_rings.conway_polynomials.conway_polynomial

    import sage.rings.polynomial.multi_polynomial_element
    MPolynomial = sage.rings.polynomial.multi_polynomial_element.MPolynomial

    import sage.rings.polynomial.polynomial_element
    Polynomial = sage.rings.polynomial.polynomial_element.Polynomial

    import sage.modules.free_module_element
    FreeModuleElement = sage.modules.free_module_element.FreeModuleElement

cdef class Cache_givaro(SageObject):
    def __init__(self, parent, unsigned int p, unsigned int k, modulus, repr="poly", cache=False):
        """
        Finite Field.

        These are implemented using Zech logs and the
        cardinality must be less than `2^{16}`. By default conway polynomials
        are used as minimal polynomial.

        INPUT:

        - ``q`` -- `p^n` (must be prime power)

        - ``name`` -- variable used for poly_repr (default: ``'a'``)

        - ``modulus`` -- a polynomial to use as modulus.

        - ``repr``  -- (default: 'poly') controls the way elements are printed
          to the user:

          - 'log': repr is :meth:`~FiniteField_givaroElement.log_repr()`
          - 'int': repr is :meth:`~FiniteField_givaroElement.int_repr()`
          - 'poly': repr is :meth:`~FiniteField_givaroElement.poly_repr()`

        - ``cache`` -- (default: ``False``) if ``True`` a cache of all
          elements of this field is created. Thus, arithmetic does not
          create new elements which speeds calculations up. Also, if many
          elements are needed during a calculation this cache reduces the
          memory requirement as at most :meth:`order()` elements are created.

        OUTPUT:

        Givaro finite field with characteristic `p` and cardinality `p^n`.

        EXAMPLES:

        By default Conway polynomials are used::

            sage: k.<a> = GF(2**8)
            sage: -a ^ k.degree()
            a^4 + a^3 + a^2 + 1
            sage: f = k.modulus(); f
            x^8 + x^4 + x^3 + x^2 + 1

        You may enforce a modulus::

            sage: P.<x> = PolynomialRing(GF(2))
            sage: f = x^8 + x^4 + x^3 + x + 1 # Rijndael polynomial
            sage: k.<a> = GF(2^8, modulus=f)
            sage: k.modulus()
            x^8 + x^4 + x^3 + x + 1
            sage: a^(2^8)
            a

        You may enforce a random modulus::

            sage: k = GF(3**5, 'a', modulus='random')
            sage: k.modulus() # random polynomial
            x^5 + 2*x^4 + 2*x^3 + x^2 + 2

        For binary fields, you may ask for a minimal weight polynomial::

            sage: k = GF(2**10, 'a', modulus='minimal_weight')
            sage: k.modulus()
            x^10 + x^3 + 1
        """
        # we are calling late_import here because this constructor is
        # called at least once before any arithmetic is performed.
        late_import()

        cdef intvec cPoly

        self.parent = <Parent?> parent

        if repr=='poly':
            self.repr = 0
        elif repr=='log':
            self.repr = 1
        elif repr=='int':
            self.repr = 2
        else:
            raise RuntimeError

        if k == 1:
            sig_on()
            self.objectptr = gfq_factorypk(p, k)
        else:
            # Givaro does not support this when k == 1
            for coeff in modulus:
                cPoly.push_back(<int>coeff)
            sig_on()
            self.objectptr = gfq_factorypkp(p, k, cPoly)

        self._zero_element = make_FiniteField_givaroElement(self, self.objectptr.zero)
        self._one_element = make_FiniteField_givaroElement(self, self.objectptr.one)
        sig_off()

        parent._zero_element = self._zero_element
        parent._one_element = self._one_element
        if cache:
            self._array = self.gen_array()
            self._has_array = True

    cdef gen_array(self):
        """
        Generates an array/list/tuple containing all elements of ``self``
        indexed by their power with respect to the internal generator.
        """
        cdef int i

        array = list()
        for i from 0 <= i < self.order_c():
            array.append(make_FiniteField_givaroElement(self,i) )
        return tuple(array)

    def __dealloc__(self):
        """
        Free the memory occupied by this Givaro finite field.
        """
        delete(self.objectptr)

    cpdef int characteristic(self):
        """
        Return the characteristic of this field.

        EXAMPLES::

            sage: p = GF(19^3,'a')._cache.characteristic(); p
            19
        """
        return self.objectptr.characteristic()

    def order(self):
        """
        Returns the order of this field.

        EXAMPLES::

            sage: K.<a> = GF(9)
            sage: K._cache.order()
            9
        """
        return Integer(self.order_c())

    cpdef int order_c(self):
        """
        Returns the order of this field.

        EXAMPLES::

            sage: K.<a> = GF(9)
            sage: K._cache.order_c()
            9
        """
        return self.objectptr.cardinality()

    cpdef int exponent(self):
        r"""
        Returns the degree of this field over `\GF{p}`.

        EXAMPLES::

            sage: K.<a> = GF(9); K._cache.exponent()
            2
        """
        return self.objectptr.exponent()

    def random_element(self, *args, **kwds):
        """
        Return a random element of ``self``.

        EXAMPLES::

            sage: k = GF(23**3, 'a')
            sage: e = k._cache.random_element(); e
            2*a^2 + 14*a + 21
            sage: type(e)
            <type 'sage.rings.finite_rings.element_givaro.FiniteField_givaroElement'>

            sage: P.<x> = PowerSeriesRing(GF(3^3, 'a'))
            sage: P.random_element(5)
            2*a + 2 + (a^2 + a + 2)*x + (2*a + 1)*x^2 + (2*a^2 + a)*x^3 + 2*a^2*x^4 + O(x^5)
        """
        cdef int seed = current_randstate().c_random()
        cdef int res
        cdef GivRandom generator = GivRandomSeeded(seed)
        res = self.objectptr.random(generator,res)
        return make_FiniteField_givaroElement(self,res)

    cpdef FiniteField_givaroElement element_from_data(self, e):
        """
        Coerces several data types to ``self``.

        INPUT:

        - ``e`` -- data to coerce in.

        EXAMPLES::

            sage: k = GF(3^8, 'a')
            sage: type(k)
            <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
            sage: e = k.vector_space().gen(1); e
            (0, 1, 0, 0, 0, 0, 0, 0)
            sage: k(e) #indirect doctest
            a

        TESTS:

        Check coercion of large integers::

            sage: k(-5^13)
            1
            sage: k(2^31)
            2
            sage: k(int(10^19))
            1
            sage: k(2^63)
            2
            sage: k(2^100)
            1
            sage: k(int(2^100))
            1
            sage: k(long(2^100))
            1
            sage: k(-2^100)
            2

        For more examples, see
        ``finite_field_givaro.FiniteField_givaro._element_constructor_``
        """
        cdef int res
        cdef int g
        cdef int x
        cdef int e_int

        cdef FiniteField_givaroElement to_add
        ########

        if isinstance(e, FiniteField_givaroElement):
            if e.parent() is self.parent:
                return e
            if e.parent() == self.parent:
                return make_FiniteField_givaroElement(self,(<FiniteField_givaroElement>e).element)
            if e.parent() is self.parent.prime_subfield() or e.parent() == self.parent.prime_subfield():
                res = self.int_to_log(int(e))
            else:
                raise TypeError("unable to coerce from a finite field other than the prime subfield")

        elif isinstance(e, int) or \
             isinstance(e, Integer) or \
             isinstance(e, long) or is_IntegerMod(e):
            try:
                e_int = e
            except OverflowError:
                e_int = e % self.characteristic()
            res = self.objectptr.initi(res, e_int)

        elif e is None:
            e_int = 0
            res = self.objectptr.initi(res,e_int)

        elif isinstance(e, float):
            res = self.objectptr.initd(res,e)

        elif isinstance(e, str):
            return self.parent(eval(e.replace("^","**"),self.parent.gens_dict()))

        elif isinstance(e, FreeModuleElement):
            if self.parent.vector_space() != e.parent():
                raise TypeError("e.parent must match self.vector_space")
            ret = self._zero_element
            for i in range(len(e)):
                e_entry = e[i] % self.characteristic()
                res = self.objectptr.initi(res, int(e_entry))
                to_add = make_FiniteField_givaroElement(self, res)
                ret = ret + to_add*self.parent.gen()**i
            return ret

        elif isinstance(e, MPolynomial):
            if e.is_constant():
                return self.parent(e.constant_coefficient())
            else:
                raise TypeError("no coercion defined")

        elif isinstance(e, Polynomial):
            if e.is_constant():
                return self.parent(e.constant_coefficient())
            else:
                return e.change_ring(self.parent)(self.parent.gen())

        elif isinstance(e, Rational):
            num = e.numer()
            den = e.denom()
            return self.parent(num)/self.parent(den)

        elif isinstance(e, gen):
            pass # handle this in next if clause

        elif isinstance(e, FiniteFieldElement_pari_ffelt) or isinstance(e, FiniteField_ext_pariElement):
            # Reduce to pari
            e = e._pari_()

        elif sage.interfaces.gap.is_GapElement(e):
            from sage.interfaces.gap import gfq_gap_to_sage
            return gfq_gap_to_sage(e, self.parent)

        elif isinstance(e, list):
            if len(e) > self.exponent():
                # could reduce here...
                raise ValueError("list is too long")
            ret = self._zero_element
            for i in range(len(e)):
                e_entry = e[i] % self.characteristic()
                res = self.objectptr.initi(res, int(e_entry))
                to_add = make_FiniteField_givaroElement(self, res)
                ret = ret + to_add*self.parent.gen()**i
            return ret

        else:
            raise TypeError("unable to coerce %r" % type(e))

        cdef GEN t
        cdef long c
        if isinstance(e, gen):
            sig_on()
            t = (<gen>e).g
            if typ(t) == t_FFELT:
                t = FF_to_FpXQ(t)
            else:
                t = liftall_shallow(t)

            if typ(t) == t_INT:
                res = self.int_to_log(itos(t))
                sig_off()
            elif typ(t) == t_POL:
                res = self._zero_element

                g = self.objectptr.indeterminate()
                x = self.objectptr.one

                for i from 0 <= i <= degpol(t):
                    c = gtolong(gel(t, i+2))
                    res = self.objectptr.axpyin(res, self.int_to_log(c), x)
                    x = self.objectptr.mul(x,x,g)
                sig_off()
            else:
                raise TypeError("bad PARI type %r" % e.type())

        return make_FiniteField_givaroElement(self,res)

    cpdef FiniteField_givaroElement gen(self):
        """
        Returns a generator of the field.

        EXAMPLES::

            sage: K.<a> = GF(625)
            sage: K._cache.gen()
            a
        """
        cdef int g
        if self.objectptr.exponent() == 1:
            self.objectptr.initi(g, -self.parent.modulus()[0])
        else:
            g = self.objectptr.indeterminate()
        return make_FiniteField_givaroElement(self, g)

    cpdef int log_to_int(self, int n) except -1:
        r"""
        Given an integer `n` this method returns `i` where `i`
        satisfies `g^n = i` where `g` is the generator of ``self``; the
        result is interpreted as an integer.

        INPUT:

        - ``n`` -- log representation of a finite field element

        OUTPUT:

        integer representation of a finite field element.

        EXAMPLES::

            sage: k = GF(2**8, 'a')
            sage: k._cache.log_to_int(4)
            16
            sage: k._cache.log_to_int(20)
            180
        """
        if n < 0:
            raise IndexError("Cannot serve negative exponent %d"%n)
        elif n >= self.order_c():
            raise IndexError("n=%d must be < self.order()"%n)

        cdef int r
        sig_on()
        self.objectptr.convert(r, n)
        sig_off()
        return r

    cpdef int int_to_log(self, int n) except -1:
        r"""
        Given an integer `n` this method returns `i` where `i` satisfies
        `g^i = n \mod p` where `g` is the generator and `p` is the
        characteristic of ``self``.

        INPUT:

        - ``n`` -- integer representation of an finite field element

        OUTPUT:

        log representation of ``n``

        EXAMPLES::

            sage: k = GF(7**3, 'a')
            sage: k._cache.int_to_log(4)
            228
            sage: k._cache.int_to_log(3)
            57
            sage: k.gen()^57
            3
        """
        cdef int r
        sig_on()
        self.objectptr.initi(r,n)
        sig_off()
        return r

    def fetch_int(self, int n):
        r"""
        Given an integer ``n`` return a finite field element in ``self``
        which equals ``n`` under the condition that :meth:`gen()` is set to
        :meth:`characteristic()`.

        EXAMPLES::

            sage: k.<a> = GF(2^8)
            sage: k._cache.fetch_int(8)
            a^3
            sage: e = k._cache.fetch_int(151); e
            a^7 + a^4 + a^2 + a + 1
            sage: 2^7 + 2^4 + 2^2 + 2 + 1
            151
        """
        cdef GivaroGfq *k = self.objectptr
        cdef int ret = k.zero
        cdef int a = k.indeterminate()
        cdef int at = k.one
        cdef unsigned int ch = k.characteristic()
        cdef int _n, t, i

        if n<0 or n>k.cardinality():
            raise TypeError("n must be between 0 and self.order()")

        _n = n

        for i from 0 <= i < k.exponent():
            t = k.initi(t, _n%ch)
            ret = k.axpy(ret, t, at, ret)
            at = k.mul(at,at,a)
            _n = _n/ch
        return make_FiniteField_givaroElement(self, ret)

    def _element_repr(self, FiniteField_givaroElement e):
        """
        Wrapper for log, int, and poly representations.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._cache._element_repr(a^20)
            '2*a^3 + 2*a^2 + 2'

            sage: k = FiniteField(3^4,'a', impl='givaro', repr='int')
            sage: a = k.gen()
            sage: k._cache._element_repr(a^20)
            '74'

            sage: k = FiniteField(3^4,'a', impl='givaro', repr='log')
            sage: a = k.gen()
            sage: k._cache._element_repr(a^20)
            '20'
        """
        if self.repr==0:
            return self._element_poly_repr(e)
        elif self.repr==1:
            return self._element_log_repr(e)
        else:
            return self._element_int_repr(e)

    def _element_log_repr(self, FiniteField_givaroElement e):
        """
        Return ``str(i)`` where ``self` is ``gen^i`` with ``gen``
        being the *internal* multiplicative generator of this finite
        field.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._cache._element_log_repr(a^20)
            '20'
            sage: k._cache._element_log_repr(a)
            '1'
        """
        return str(int(e.element))

    def _element_int_repr(self, FiniteField_givaroElement e):
        r"""
        Return integer representation of ``e``.

        Elements of this field are represented as ints in as follows:
        for `e \in \GF{p}[x]` with `e = a_0 + a_1x + a_2x^2 + \cdots`, `e` is
        represented as: `n = a_0 + a_1  p + a_2  p^2 + \cdots`.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._cache._element_int_repr(a^20)
            '74'
        """
        return str(e.integer_representation())

    def _element_poly_repr(self, FiniteField_givaroElement e, varname = None):
        """
        Return a polynomial expression in the generator of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: k._cache._element_poly_repr(a^20)
            '2*a^3 + 2*a^2 + 2'
        """
        if varname is None:
            variable = self.parent.variable_name()
        else:
            variable = varname

        quo = self.log_to_int(e.element)
        b   = int(self.characteristic())

        ret = ""
        for i in range(self.exponent()):
            coeff = quo%b
            if coeff != 0:
                if i>0:
                    if coeff==1:
                        coeff=""
                    else:
                        coeff=str(coeff)+"*"
                    if i>1:
                        ret = coeff + variable + "^" + str(i) + " + " + ret
                    else:
                        ret = coeff + variable + " + " + ret
                else:
                    ret = str(coeff) + " + " + ret
            quo = quo/b
        if ret == '':
            return "0"
        return ret[:-3]

    def a_times_b_plus_c(self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return ``a*b + c``. This is faster than multiplying ``a`` and ``b``
        first and adding ``c`` to the result.

        INPUT:

        - ``a,b,c`` -- :class:`FiniteField_givaroElement`

        EXAMPLES::

            sage: k.<a> = GF(2**8)
            sage: k._cache.a_times_b_plus_c(a,a,k(1))
            a^2 + 1
        """
        cdef int r

        r = self.objectptr.axpy(r, a.element, b.element, c.element)
        return make_FiniteField_givaroElement(self,r)

    def a_times_b_minus_c(self,FiniteField_givaroElement a, FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return ``a*b - c``.

        INPUT:

        - ``a,b,c`` -- :class:`FiniteField_givaroElement`

        EXAMPLES::

            sage: k.<a> = GF(3**3)
            sage: k._cache.a_times_b_minus_c(a,a,k(1))
            a^2 + 2
        """
        cdef int r

        r = self.objectptr.axmy(r, a.element, b.element, c.element, )
        return make_FiniteField_givaroElement(self,r)

    def c_minus_a_times_b(self,FiniteField_givaroElement a,
                          FiniteField_givaroElement b, FiniteField_givaroElement c):
        """
        Return ``c - a*b``.

        INPUT:

        - ``a,b,c`` -- :class:`FiniteField_givaroElement`

        EXAMPLES::

            sage: k.<a> = GF(3**3)
            sage: k._cache.c_minus_a_times_b(a,a,k(1))
            2*a^2 + 1
        """
        cdef int r

        r = self.objectptr.maxpy(r , a.element, b.element, c.element, )
        return make_FiniteField_givaroElement(self,r)

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: k.<a> = GF(3^8)
            sage: TestSuite(a).run()
        """
        p, k = self.order().factor()[0]
        if self.repr == 0:
            rep = 'poly'
        elif self.repr == 1:
            rep = 'log'
        elif self.repr == 2:
            rep = 'int'
        return unpickle_Cache_givaro, (self.parent, p, k, self.parent.polynomial(), rep, self._has_array)

    cdef FiniteField_givaroElement _new_c(self, int value):
        return make_FiniteField_givaroElement(self, value)


def unpickle_Cache_givaro(parent, p, k, modulus, rep, cache):
    """
    EXAMPLES::

       sage: k = GF(3**7, 'a')
       sage: loads(dumps(k)) == k # indirect doctest
       True
    """
    return Cache_givaro(parent, p, k, modulus, rep, cache)

cdef class FiniteField_givaro_iterator:
    """
    Iterator over :class:`FiniteField_givaro` elements.  We iterate
    multiplicatively, as powers of a fixed internal generator.

    EXAMPLES::

        sage: for x in GF(2^2,'a'): print(x)
        0
        a
        a + 1
        1
    """

    def __init__(self, Cache_givaro cache):
        """
        EXAMPLE::

            sage: k.<a> = GF(3^4)
            sage: i = iter(k) # indirect doctest
            sage: i
            Iterator over Finite Field in a of size 3^4
        """
        self._cache = cache
        self.iterator = -1

    def __next__(self):
        """
        EXAMPLE::

            sage: k.<a> = GF(3^4)
            sage: i = iter(k) # indirect doctest
            sage: next(i)
            0
            sage: next(i)
            a
        """

        self.iterator=self.iterator+1

        if self.iterator==self._cache.order_c():
            self.iterator = -1
            raise StopIteration

        return make_FiniteField_givaroElement(self._cache,self.iterator)

    def __repr__(self):
        """
        EXAMPLE::

            sage: k.<a> = GF(3^4)
            sage: i = iter(k)
            sage: i # indirect doctest
            Iterator over Finite Field in a of size 3^4
        """
        return "Iterator over %s"%self._cache.parent

    def __iter__(self):
        """
        EXAMPLE::

            sage: K.<a> = GF(4)
            sage: K.list() # indirect doctest
            [0, a, a + 1, 1]
        """
        return self

cdef class FiniteField_givaroElement(FinitePolyExtElement):
    """
    An element of a (Givaro) finite field.
    """

    def __init__(FiniteField_givaroElement self, parent ):
        """
        Initializes an element in parent. It's much better to use
        parent(<value>) or any specialized method of parent
        like gen() instead. In general do not call this
        constructor directly.

        Alternatively you may provide a value which is directly
        assigned to this element. So the value must represent the
        log_g of the value you wish to assign.

        INPUT:

        - ``parent`` -- base field

        OUTPUT:

        A finite field element.

        EXAMPLES::

            sage: k.<a> = GF(5^2)
            sage: from sage.rings.finite_rings.element_givaro import FiniteField_givaroElement
            sage: FiniteField_givaroElement(k)
            0

        """
        FinitePolyExtElement.__init__(self, parent)
        self._cache = parent._cache
        self.element = 0

    cdef FiniteField_givaroElement _new_c(self, int value):
        return make_FiniteField_givaroElement(self._cache, value)

    def __dealloc__(FiniteField_givaroElement self):
        pass

    def _repr_(FiniteField_givaroElement self):
        """
        EXAMPLE::

            sage: k.<FOOBAR> = GF(3^4)
            sage: FOOBAR #indirect doctest
            FOOBAR

            sage: k.<FOOBAR> = GF(3^4,repr='log')
            sage: FOOBAR
            1

            sage: k.<FOOBAR> = GF(3^4,repr='int')
            sage: FOOBAR
            3
        """
        return self._cache._element_repr(self)

    def _element(self):
        """
        Returns the int interally representing this element.

        EXAMPLES::

            sage: k.<a> = GF(3^4)
            sage: (a^2 + 1)._element()
            58
        """
        return self.element

    def __nonzero__(FiniteField_givaroElement self):
        r"""
        Return ``True`` if ``self != k(0)``.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: a.is_zero()
            False
            sage: k(0).is_zero()
            True
        """
        return not self._cache.objectptr.isZero(self.element)

    def is_one(FiniteField_givaroElement self):
        r"""
        Return ``True`` if ``self == k(1)``.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: a.is_one()
            False
            sage: k(1).is_one()
            True
        """
        return self._cache.objectptr.isOne(self.element)

    def is_unit(FiniteField_givaroElement self):
        """
        Return ``True`` if self is nonzero, so it is a unit as an element of
        the finite field.

        EXAMPLES::

            sage: k.<a> = GF(3^4); k
            Finite Field in a of size 3^4
            sage: a.is_unit()
            True
            sage: k(0).is_unit()
            False
        """
        return not (<Cache_givaro>self._cache).objectptr.isZero(self.element)
        # **WARNING** Givaro seems to define unit to mean in the prime field,
        # which is totally wrong!  It's a confusion with the underlying polynomial
        # representation maybe??  That's why the following is commented out.
        # return (<FiniteField_givaro>self._parent).objectptr.isunit(self.element)


    def is_square(FiniteField_givaroElement self):
        """
        Return ``True`` if ``self`` is a square in ``self.parent()``

        ALGORITHM:

        Elements are stored as powers of generators, so we simply check
        to see if it is an even power of a generator.

        EXAMPLES::

            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: a.is_square()
            False
            sage: v = set([x^2 for x in k])
            sage: [x.is_square() for x in v]
            [True, True, True, True, True]
            sage: [x.is_square() for x in k if not x in v]
            [False, False, False, False]

        TESTS::

            sage: K = GF(27, 'a')
            sage: set([a*a for a in K]) == set([a for a in K if a.is_square()])
            True
            sage: K = GF(25, 'a')
            sage: set([a*a for a in K]) == set([a for a in K if a.is_square()])
            True
            sage: K = GF(16, 'a')
            sage: set([a*a for a in K]) == set([a for a in K if a.is_square()])
            True
        """
        cdef Cache_givaro cache = <Cache_givaro>self._cache
        if cache.objectptr.characteristic() == 2:
            return True
        elif self.element == cache.objectptr.one:
            return True
        else:
            return self.element % 2 == 0

    def sqrt(FiniteField_givaroElement self, extend=False, all=False):
        """
        Return a square root of this finite field element in its
        parent, if there is one.  Otherwise, raise a ``ValueError``.

        INPUT:

        - ``extend`` -- bool (default: ``True``); if ``True``, return a
          square root in an extension ring, if necessary. Otherwise,
          raise a ``ValueError`` if the root is not in the base ring.

          .. WARNING::

              this option is not implemented!

        - ``all`` -- bool (default: ``False``); if ``True``, return all square
          roots of ``self``, instead of just one.

        .. WARNING::

            The ``extend`` option is not implemented (yet).

        ALGORITHM:

        ``self`` is stored as `a^k` for some generator `a`.
        Return `a^{k/2}` for even `k`.

        EXAMPLES::

            sage: k.<a> = GF(7^2)
            sage: k(2).sqrt()
            3
            sage: k(3).sqrt()
            2*a + 6
            sage: k(3).sqrt()**2
            3
            sage: k(4).sqrt()
            2
            sage: k.<a> = GF(7^3)
            sage: k(3).sqrt()
            Traceback (most recent call last):
            ...
            ValueError: must be a perfect square.

        TESTS::

            sage: K = GF(49, 'a')
            sage: all([a.sqrt()*a.sqrt() == a for a in K if a.is_square()])
            True
            sage: K = GF(27, 'a')
            sage: all([a.sqrt()*a.sqrt() == a for a in K if a.is_square()])
            True
            sage: K = GF(8, 'a')
            sage: all([a.sqrt()*a.sqrt() == a for a in K if a.is_square()])
            True
            sage: K.<a>=FiniteField(9)
            sage: a.sqrt(extend = False, all = True)
            []

        """
        if all:
            if self.is_square():
                a = self.sqrt()
                return [a, -a] if -a != a else [a]
            return []
        cdef Cache_givaro cache = <Cache_givaro>self._cache
        if self.element == cache.objectptr.one:
            return make_FiniteField_givaroElement(cache, cache.objectptr.one)
        elif self.element % 2 == 0:
            return make_FiniteField_givaroElement(cache, self.element/2)
        elif cache.objectptr.characteristic() == 2:
            return make_FiniteField_givaroElement(cache, (cache.objectptr.cardinality() - 1 + self.element)/2)
        elif extend:
            raise NotImplementedError # TODO: fix this once we have nested embeddings of finite fields
        else:
            raise ValueError("must be a perfect square.")

    cpdef _add_(self, right):
        """
        Add two elements.

        EXAMPLES::

            sage: k.<b> = GF(9**2)
            sage: b^10 + 2*b # indirect doctest
            2*b^3 + 2*b^2 + 2*b + 1
        """
        cdef int r
        r = self._cache.objectptr.add(r, self.element ,
                                              (<FiniteField_givaroElement>right).element )
        return make_FiniteField_givaroElement(self._cache,r)

    cpdef _mul_(self, right):
        """
        Multiply two elements.

        EXAMPLES::

            sage: k.<c> = GF(7**4)
            sage: 3*c # indirect doctest
            3*c
            sage: c*c
            c^2
        """
        cdef int r
        r = self._cache.objectptr.mul(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(self._cache,r)


    cpdef _div_(self, right):
        """
        Divide two elements

        EXAMPLES::

            sage: k.<g> = GF(2**8)
            sage: g/g # indirect doctest
            1

            sage: k(1) / k(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field.
        """
        cdef int r
        if (<FiniteField_givaroElement>right).element == 0:
            raise ZeroDivisionError('division by zero in finite field.')
        r = self._cache.objectptr.div(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(self._cache,r)

    cpdef _sub_(self, right):
        """
        Subtract two elements.

        EXAMPLES::

            sage: k.<a> = GF(3**4)
            sage: k(3) - k(1) # indirect doctest
            2
            sage: 2*a - a^2
            2*a^2 + 2*a
        """
        cdef int r
        r = self._cache.objectptr.sub(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(self._cache,r)

    def __neg__(FiniteField_givaroElement self):
        """
        Negative of an element.

        EXAMPLES::

            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: -a
            2*a
        """
        cdef int r

        r = self._cache.objectptr.neg(r, self.element)
        return make_FiniteField_givaroElement(self._cache,r)

    def __invert__(FiniteField_givaroElement self):
        """
        Return the multiplicative inverse of an element.

        EXAMPLES::

            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: ~a
            a + 2
            sage: ~a*a
            1

        TEST:

        Check that trying to invert zero raises an error
        (see :trac:`12217`)::

            sage: F = GF(25, 'a')
            sage: z = F(0)
            sage: ~z
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in Finite Field in a of size 5^2

        """
        cdef int r
        if (<FiniteField_givaroElement>self).element == 0:
            raise ZeroDivisionError('division by zero in %s'
                                    % self.parent())
        self._cache.objectptr.inv(r, self.element)
        return make_FiniteField_givaroElement(self._cache,r)

    def __pow__(FiniteField_givaroElement self, exp, other):
        """
        EXAMPLES::

            sage: K.<a> = GF(3^3, 'a')
            sage: a^3 == a*a*a
            True
            sage: b = a+1
            sage: b^5 == b^2 * b^3
            True
            sage: b^(-1) == 1/b
            True
            sage: b = K(-1)
            sage: b^2 == 1
            True

        TESTS:

        The following checks that :trac:`7923` is resolved::

            sage: K.<a> = GF(3^10)
            sage: b = a^9 + a^7 + 2*a^6 + a^4 + a^3 + 2*a^2 + a + 2
            sage: b^(71*7381) == (b^71)^7381
            True

        We define ``0^0`` to be unity, :trac:`13897`::

            sage: K.<a> = GF(3^10)
            sage: K(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: K.<a> = GF(3^10)
            sage: type(K(0)^0) == type(K(0))
            True

        ALGORITHM:

        Givaro objects are stored as integers `i` such that ``self`` `= a^i`,
        where `a` is a generator of `K` (though not necessarily the one
        returned by ``K.gens()``).  Now it is trivial to compute
        `(a^i)^e = a^{i \cdot e}`, and reducing the exponent
        mod the multiplicative order of `K`.

        AUTHOR:

        - Robert Bradshaw
        """
        if not isinstance(exp, (int, Integer)):
            _exp = Integer(exp)
            if _exp != exp:
                raise ValueError("exponent must be an integer")
            exp = _exp

        cdef Cache_givaro cache = self._cache

        if (cache.objectptr).isOne(self.element):
            return self

        elif exp == 0:
            return make_FiniteField_givaroElement(cache, cache.objectptr.one)

        elif (cache.objectptr).isZero(self.element):
            if exp < 0:
                raise ZeroDivisionError
            return make_FiniteField_givaroElement(cache, cache.objectptr.zero)

        cdef int order = (cache.order_c()-1)
        cdef int r = exp % order

        if r == 0:
            return make_FiniteField_givaroElement(cache, cache.objectptr.one)

        cdef unsigned int r_unsigned
        if r < 0:
            r_unsigned = <unsigned int> r + order
        else:
            r_unsigned = <unsigned int>r
        cdef unsigned int elt_unsigned = <unsigned int>self.element
        cdef unsigned int order_unsigned = <unsigned int>order
        r = <int>(r_unsigned * elt_unsigned) % order_unsigned
        if r == 0:
            return make_FiniteField_givaroElement(cache, cache.objectptr.one)
        return make_FiniteField_givaroElement(cache, r)

    cpdef int _cmp_(left, right) except -2:
        """
        Comparison of finite field elements is correct or equality
        tests and somewhat random for ``<`` and ``>`` type of
        comparisons. This implementation performs these tests by
        comparing the underlying int representations.

        EXAMPLES::

            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: a == k('a')
            True
            sage: a == a + 1
            False

        Even though inequality tests do return answers, they
        really make no sense as finite fields are unordered. Thus,
        you cannot rely on the result as it is implementation specific.

        ::

            sage: a < a^2
            True
            sage: a^2 < a
            False
        """
        cdef Cache_givaro cache = (<FiniteField_givaroElement>left)._cache

        return cmp(cache.log_to_int(left.element), cache.log_to_int((<FiniteField_givaroElement>right).element) )

    def __int__(FiniteField_givaroElement self):
        """
        Return the int representation of ``self``.  When ``self`` is in the
        prime subfield, the integer returned is equal to ``self``, otherwise
        an error is raised.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: int(k(4))
            4
            sage: int(b)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce element to an integer.

        """
        cdef int self_int = self._cache.log_to_int(self.element)
        if self_int%self._cache.characteristic() != self_int:
            raise TypeError("Cannot coerce element to an integer.")
        return self_int

    def integer_representation(FiniteField_givaroElement self):
        """
        Return the integer representation of ``self``.  When ``self`` is in the
        prime subfield, the integer returned is equal to ``self``.

        Elements of this field are represented as integers as follows:
        given the element `e \in \GF{p}[x]` with
        `e = a_0 + a_1 x + a_2 x^2 + \cdots`, the integer representation
        is `a_0 + a_1 p + a_2 p^2 + \cdots`.

        OUTPUT: A Python ``int``.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: k(4).integer_representation()
            4
            sage: b.integer_representation()
            5
            sage: type(b.integer_representation())
            <type 'int'>
        """
        return self._cache.log_to_int(self.element)

    def _integer_(FiniteField_givaroElement self, ZZ=None):
        """
        Convert ``self`` to an integer if it is in the prime subfield.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: k(4)._integer_()
            4
            sage: ZZ(b)
            Traceback (most recent call last):
            ...
            TypeError: not in prime subfield
        """
        cdef int a = self._cache.log_to_int(self.element)
        if a < self._cache.objectptr.characteristic():
            return Integer(a)
        raise TypeError("not in prime subfield")

    def _log_to_int(FiniteField_givaroElement self):
        r"""
        Return the int representation of ``self``, as a Sage integer.

        Elements of this field are represented as ints as follows:
        given the element `e \in \GF{p}[x]` with
        `e = a_0 + a_1x + a_2x^2 + \cdots`, the int representation is
        `a_0 + a_1 p + a_2 p^2 + \cdots`.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: k(4)._log_to_int()
            4
            sage: b._log_to_int()
            5
            sage: type(b._log_to_int())
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(self._cache.log_to_int(self.element))

    log_to_int = deprecated_function_alias(11295, _log_to_int)

    def log(FiniteField_givaroElement self, base):
        """
        Return the log to the base `b` of ``self``, i.e., an integer `n`
        such that `b^n =` ``self``.

        .. WARNING::

            TODO -- This is currently implemented by solving the discrete
            log problem -- which shouldn't be needed because of how finite field
            elements are represented.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: a = b^7
            sage: a.log(b)
            7
        """
        b = self.parent()(base)
        return sage.groups.generic.discrete_log(self, b)

    def _int_repr(FiniteField_givaroElement self):
        r"""
        Return the string representation of ``self`` as an int (as returned
        by :meth:`log_to_int`).

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+1)._int_repr()
            '6'
        """
        return self._cache._element_int_repr(self)

    int_repr = deprecated_function_alias(11295, _int_repr)

    def _log_repr(FiniteField_givaroElement self):
        r"""
        Return the log representation of ``self`` as a string.  See the
        documentation of the ``_element_log_repr`` function of the
        parent field.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+2)._log_repr()
            '15'
        """
        return self._cache._element_log_repr(self)

    log_repr = deprecated_function_alias(11295, _log_repr)

    def _poly_repr(FiniteField_givaroElement self):
        r"""
        Return representation of this finite field element as a polynomial
        in the generator.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+2)._poly_repr()
            'b + 2'
        """
        return self._cache._element_poly_repr(self)

    poly_repr = deprecated_function_alias(11295, _poly_repr)

    def polynomial(FiniteField_givaroElement self, name=None):
        """
        Return self viewed as a polynomial over
        ``self.parent().prime_subfield()``.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: f = (b^2+1).polynomial(); f
            b + 4
            sage: type(f)
            <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>
            sage: parent(f)
            Univariate Polynomial Ring in b over Finite Field of size 5
        """
        cdef Cache_givaro cache = self._cache
        K = self.parent()
        quo = cache.log_to_int(self.element)
        b   = int(cache.characteristic())
        ret = []
        for i in range(K.degree()):
            coeff = quo%b
            ret.append(coeff)
            quo = quo/b
        if not name is None and K.variable_name() != name:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            return PolynomialRing(K.prime_subfield(), name)(ret)
        else:
            return K.polynomial_ring()(ret)

    def _magma_init_(self, magma):
        """
        Return a string representation of self that MAGMA can
        understand.

        EXAMPLE::

            sage: k.<a> = GF(3^5)

        String rep of parent::

            sage: k._magma_init_(magma)        # optional - magma
            'SageCreateWithNames(ext<GF(3)|_sage_[...]![GF(3)!1,GF(3)!2,GF(3)!0,GF(3)!0,GF(3)!0,GF(3)!1]>,["a"])'

        Magma repr of element::

            sage: a._magma_init_(magma)        # optional - magma
             '_sage_[...]!(_sage_[...])'

        Because of caching the string representation of an element must
        not change::

            sage: a._magma_init_(magma) == a._magma_init_(magma)   # optional - magma
            True

        We test a conversion back and forth::

            sage: k.<a> = GF(3^6)
            sage: b = magma(a^5 + 2*a^2 + 1)             # optional - magma

        Note that small fields print using a log representation in Magma
        (unlike Sage)::

            sage: b                                      # optional - magma
            a^436
            sage: b.sage()                               # optional - magma
            a^5 + 2*a^2 + 1
        """
        R = magma(self.parent())
        a = R.gen(1).name()
        return '%s!(%s)'%(R.name(), self._cache._element_poly_repr(self, a))

    def multiplicative_order(FiniteField_givaroElement self):
        """
        Return the multiplicative order of this field element.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.multiplicative_order()
            24
            sage: (b^6).multiplicative_order()
            4
        """
        # TODO -- I'm sure this can be made vastly faster
        # using how elements are represented as a power of the generator ??

        # code copy'n'pasted from element_ext_pari.py
        import sage.arith.all

        if self._multiplicative_order is not None:
            return self._multiplicative_order
        else:
            if self.is_zero():
                raise ArithmeticError("Multiplicative order of 0 not defined.")
            n = (self._cache).order_c() - 1
            order = Integer(1)
            for p, e in sage.arith.all.factor(n):
                # Determine the power of p that divides the order.
                a = self**(n/(p**e))
                while a != 1:
                    order = order * p
                    a = a**p

            self._multiplicative_order = order
            return order

    def __copy__(self):
        """
        Return a copy of this element.  Actually just returns ``self``, since
        finite field elements are immutable.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: c = copy(b); c
            b
            sage: c is b
            True
            sage: copy(5r) is 5r
            True
        """
        return self

    def _gap_init_(FiniteField_givaroElement self):
        """
        Return a string that evaluates to the GAP representation of
        this element.

        A ``NotImplementedError`` is raised if ``self.parent().modulus()`` is
        not a Conway polynomial, as the isomorphism of finite fields is
        not implemented yet.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: (4*b+3)._gap_init_()
            'Z(25)^3'
            sage: S(gap('Z(25)^3'))
            4*b + 3
        """
        #copied from element_ext_pari.py
        cdef Cache_givaro cache = self._cache
        if self == 0:
            return '0*Z(%s)'%cache.order_c()
        F = self.parent()
        if F.degree() == 1:
            # Find the root of unity used by Gap.  See _gap_init_ in sage.rings.finite_rings.integer_mod
            from sage.interfaces.all import gap        # here to reduce dependencies
            from sage.rings.finite_rings.integer_mod import mod
            g = int(gap.eval('Int(Z(%s))'%cache.order_c()))
            n = self.log(mod(g, cache.order_c()))
            return 'Z(%s)^%s'%(cache.order_c(), n)
        if not F.is_conway():
            raise NotImplementedError("conversion of (Givaro) finite field element to GAP not implemented except for fields defined by Conway polynomials.")
        if cache.order_c() > 65536:
            raise TypeError("order (=%s) must be at most 65536." % F.order_c())
        g = F.multiplicative_generator()
        n = self.log(g)
        return 'Z(%s)^%s'%(cache.order_c(), n)

    def __hash__(FiniteField_givaroElement self):
        """
        Return the hash of this finite field element.  We hash the parent
        and the underlying integer representation of this element.

        EXAMPLES::

            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: hash(a)
            5
        """
        return hash(self.integer_representation())

    def _vector_(FiniteField_givaroElement self, reverse=False):
        """
        Return a vector in ``self.parent().vector_space()`` matching
        ``self``. The most significant bit is to the right.

        INPUT:

        - ``reverse`` -- reverse the order of the bits from little endian to
          big endian.

        EXAMPLES::

            sage: k.<a> = GF(2^4)
            sage: e = a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 1, 0)
            sage: k(v)
            a^2 + 1

            sage: k.<a> = GF(3^4)
            sage: e = 2*a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 2, 0)
            sage: k(v)
            2*a^2 + 1

        You can also compute the vector in the other order::

            sage: e._vector_(reverse=True)
            (0, 2, 0, 1)
        """
        #vector(foo) might pass in ZZ
        if isinstance(reverse, Parent):
            raise TypeError("Base field is fixed to prime subfield.")
        cdef Cache_givaro cache = self._cache
        k = self.parent()

        quo = cache.log_to_int(self.element)
        b   = int(k.characteristic())

        ret = []
        for i in range(k.degree()):
            coeff = quo%b
            ret.append(coeff)
            quo = quo/b
        if reverse:
            ret = list(reversed(ret))
        return k.vector_space()(ret)

    def __reduce__(FiniteField_givaroElement self):
        """
        Used for supporting pickling of finite field elements.

        EXAMPLES::

            sage: k = GF(2**8, 'a')
            sage: e = k.random_element()
            sage: TestSuite(e).run() # indirect doctest
        """
        return unpickle_FiniteField_givaroElement,(self.parent(),self.element)

def unpickle_FiniteField_givaroElement(parent, int x):
    """
    TESTS::

        sage: k = GF(3**4, 'a')
        sage: e = k.random_element()
        sage: TestSuite(e).run() # indirect doctest
    """
    return make_FiniteField_givaroElement(parent._cache, x)

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.rings.finite_field_givaro', 'unpickle_FiniteField_givaroElement', unpickle_FiniteField_givaroElement)

cdef inline FiniteField_givaroElement make_FiniteField_givaroElement(Cache_givaro cache, int x):
    cdef FiniteField_givaroElement y

    if cache._has_array:
        return <FiniteField_givaroElement>cache._array[x]
    else:
        y = FiniteField_givaroElement.__new__(FiniteField_givaroElement)
        y._parent = <Parent> cache.parent
        y._cache = cache
        y.element = x
        return y

