r"""
Givaro Field Elements

Sage includes the Givaro finite field library, for highly optimized
arithmetic in finite fields.

.. NOTE::

    The arithmetic is performed by the Givaro C++ library which uses Zech
    logs internally to represent finite field elements. This
    implementation is the default finite extension field implementation in
    Sage for the cardinality less than `2^{16}`, as it is vastly faster than
    the PARI implementation which uses polynomials to represent finite field
    elements. Some functionality in this class however is implemented
    using the PARI implementation.

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
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.misc.randstate cimport randstate, current_randstate
from sage.rings.finite_rings.finite_field_base cimport FiniteField
from sage.rings.ring cimport Ring
from sage.rings.finite_rings.element_ext_pari import FiniteField_ext_pariElement
from sage.structure.sage_object cimport SageObject
import operator
import sage.rings.arith
import constructor as finite_field
import finite_field_ext_pari

import sage.interfaces.gap
from sage.libs.pari.all import pari
from sage.libs.pari.gen import gen

from sage.structure.parent  cimport Parent
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens

cdef object is_IntegerMod
cdef object Integer
cdef object Rational
cdef object is_Polynomial
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
           is_Polynomial, \
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

    import sage.rings.polynomial.polynomial_element
    is_Polynomial = sage.rings.polynomial.polynomial_element.is_Polynomial

    import sage.databases.conway
    ConwayPolynomials = sage.databases.conway.ConwayPolynomials

    import sage.rings.finite_rings.constructor
    conway_polynomial = sage.rings.finite_rings.conway_polynomials.conway_polynomial

    import sage.rings.polynomial.multi_polynomial_element
    MPolynomial = sage.rings.polynomial.multi_polynomial_element.MPolynomial

    import sage.rings.polynomial.polynomial_element
    Polynomial = sage.rings.polynomial.polynomial_element.Polynomial

    import sage.modules.free_module_element
    FreeModuleElement = sage.modules.free_module_element.FreeModuleElement

cdef class Cache_givaro(SageObject):
    def __init__(self, parent, p, k, modulus=None, repr="poly", cache=False):
        """
        Finite Field.

        These are implemented using Zech logs and the
        cardinality must be less than `2^{16}`. By default conway polynomials
        are used as minimal polynomial.

        INPUT:

        - ``q`` -- `p^n` (must be prime power)

        - ``name`` -- variable used for poly_repr (default: ``'a'``)

        - ``modulus`` -- you may provide a polynomial to use for reduction or
          one of the following strings:

          - ``'conway'`` -- force the use of a Conway polynomial, will
            raise a ``RuntimeError`` if none is found in the database
          - ``'random'`` -- use a random irreducible polynomial
          - ``'default'`` -- a Conway polynomial is used if found. Otherwise
            a random polynomial is used

          Furthermore, for binary fields we allow two more options:

          - ``'minimal_weight'`` -- use a minimal weight polynomial, should
            result in faster arithmetic;
          - ``'first_lexicographic'`` -- use the first irreducible polynomial
            in lexicographic order.

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

        For binary fields, you may ask for a  minimal weight polynomial::

            sage: k = GF(2**10, 'a', modulus='minimal_weight')
            sage: k.modulus()
            x^10 + x^3 + 1

        Three different representations are possible::

            sage: sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro(9,repr='poly').gen()
            a
            sage: sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro(9,repr='int').gen()
            3
            sage: sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro(9,repr='log').gen()
            1
        """
        # we are calling late_import here because this constructor is
        # called at least once before any arithmetic is performed.
        late_import()

        cdef intvec cPoly
        cdef GF2X_c ntl_m, ntl_tmp
        cdef GF2_c c

        self.parent = <Parent?> parent

        if repr=='poly':
            self.repr = 0
        elif repr=='log':
            self.repr = 1
        elif repr=='int':
            self.repr = 2
        else:
            raise RuntimeError

        if isinstance(modulus,str) and p == 2:
            if modulus == "minimal_weight":
                GF2X_BuildSparseIrred(ntl_m, k)
            elif modulus == "first_lexicographic":
                GF2X_BuildIrred(ntl_m, k)
            elif modulus == "random":
                current_randstate().set_seed_ntl(False)
                GF2X_BuildSparseIrred(ntl_tmp, k)
                GF2X_BuildRandomIrred(ntl_m, ntl_tmp)
            else:
                raise ValueError, "Cannot understand modulus"

            modulus = []
            for i in range(k+1):
                c = GF2X_coeff(ntl_m, i)
                if not GF2_IsZero(c):
                    modulus.append(1)
                else:
                    modulus.append(0)

        if is_Polynomial(modulus):
            modulus = modulus.list()

        if PY_TYPE_CHECK(modulus, list) or PY_TYPE_CHECK(modulus, tuple):
            for i in modulus:
                cPoly.push_back(int( i % p ))
            sig_on()
            self.objectptr = gfq_factorypkp(p, k, cPoly)
        elif modulus == "random":
            sig_on()
            self.objectptr = gfq_factorypk(p,k)
        else:
            raise ValueError, "Cannot understand modulus"

        self._zero_element = make_FiniteField_givaroElement(self,self.objectptr.zero)
        self._one_element = make_FiniteField_givaroElement(self,self.objectptr.one)
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

            sage: k = GF(2**8, 'a')
            sage: e = k.vector_space().gen(1); e
            (0, 1, 0, 0, 0, 0, 0, 0)
            sage: k(e) #indirect doctest
            a

        For more examples, see
        ``finite_field_givaro.FiniteField_givaro._element_constructor_``
        """
        cdef int res
        cdef int g
        cdef int x
        cdef int e_int

        cdef FiniteField_givaroElement to_add
        ########

        if PY_TYPE_CHECK(e, FiniteField_givaroElement):
            if e.parent() is self.parent:
                return e
            if e.parent() == self.parent:
                return make_FiniteField_givaroElement(self,(<FiniteField_givaroElement>e).element)
            if e.parent() is self.parent.prime_subfield() or e.parent() == self.parent.prime_subfield():
                res = self.int_to_log(int(e))
            else:
                raise TypeError, "unable to coerce from a finite field other than the prime subfield"

        elif PY_TYPE_CHECK(e, int) or \
             PY_TYPE_CHECK(e, Integer) or \
             PY_TYPE_CHECK(e, long) or is_IntegerMod(e):
            try:
                e_int = e
                if e != e_int:       # overflow in Pyrex is often not detected correctly... but this is bullet proof.
                                     # sometimes it is detected correctly, so we do have to use exceptions though.
                                     # todo -- be more eloquent here!!
                    raise OverflowError
                res = self.objectptr.initi(res,e_int)
            except OverflowError:
                e = e % self.characteristic()
                res = self.objectptr.initi(res,int(e))

        elif e is None:
            e_int = 0
            res = self.objectptr.initi(res,e_int)

        elif PY_TYPE_CHECK(e, float):
            res = self.objectptr.initd(res,e)

        elif PY_TYPE_CHECK(e, str):
            return self.parent(eval(e.replace("^","**"),self.parent.gens_dict()))

        elif PY_TYPE_CHECK(e, FreeModuleElement):
            if self.parent.vector_space() != e.parent():
                raise TypeError, "e.parent must match self.vector_space"
            ret = self._zero_element
            for i in range(len(e)):
                e_entry = e[i] % self.characteristic()
                res = self.objectptr.initi(res, int(e_entry))
                to_add = make_FiniteField_givaroElement(self, res)
                ret = ret + to_add*self.parent.gen()**i
            return ret

        elif PY_TYPE_CHECK(e, MPolynomial):
            if e.is_constant():
                return self.parent(e.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif PY_TYPE_CHECK(e, Polynomial):
            if e.is_constant():
                return self.parent(e.constant_coefficient())
            else:
                return e.change_ring(self.parent)(self.parent.gen())

        elif PY_TYPE_CHECK(e, Rational):
            num = e.numer()
            den = e.denom()
            return self.parent(num)/self.parent(den)

        elif PY_TYPE_CHECK(e, gen):
            pass # handle this in next if clause

        elif PY_TYPE_CHECK(e,FiniteField_ext_pariElement):
            # reduce FiniteField_ext_pariElements to pari
            e = e._pari_()

        elif sage.interfaces.gap.is_GapElement(e):
            from sage.interfaces.gap import gfq_gap_to_sage
            return gfq_gap_to_sage(e, self.parent)

        elif isinstance(e, list):
            if len(e) > self.exponent():
                # could reduce here...
                raise ValueError, "list is too long"
            ret = self._zero_element
            for i in range(len(e)):
                e_entry = e[i] % self.characteristic()
                res = self.objectptr.initi(res, int(e_entry))
                to_add = make_FiniteField_givaroElement(self, res)
                ret = ret + to_add*self.parent.gen()**i
            return ret

        else:
            raise TypeError, "unable to coerce"

        if PY_TYPE_CHECK(e, gen):
            e = e.lift().lift()
            try:
                res = self.int_to_log(e[0])
            except TypeError:
                res = self.int_to_log(e)

            g = self.objectptr.sage_generator()
            x = self.objectptr.one

            for i from 0 < i <= e.poldegree():
                x = self.objectptr.mul(x,x,g)
                res = self.objectptr.axpyin( res, self.int_to_log(e[i]) , x)

        return make_FiniteField_givaroElement(self,res)

    cpdef FiniteField_givaroElement gen(self):
        """
        Returns a generator of the field.

        EXAMPLES::

            sage: K.<a> = GF(625)
            sage: K._cache.gen()
            a
        """
        return make_FiniteField_givaroElement(self, self.objectptr.sage_generator())

    def log_to_int(self, int n):
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
        cdef int ret

        if n<0:
            raise ArithmeticError, "Cannot serve negative exponent %d"%n
        elif n>=self.order_c():
            raise IndexError, "n=%d must be < self.order()"%n
        sig_on()
        ret = int(self.objectptr.convert(ret, n))
        sig_off()
        return ret

    def int_to_log(self, int n):
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
        ret =  int(self.objectptr.initi(r,n))
        sig_off()
        return ret

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
        cdef int a = k.sage_generator()
        cdef int at = k.one
        cdef unsigned int ch = k.characteristic()
        cdef int _n, t, i

        if n<0 or n>k.cardinality():
            raise TypeError, "n must be between 0 and self.order()"

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

            sage: k = sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro(3^4,'a', repr='int')
            sage: a = k.gen()
            sage: k._cache._element_repr(a^20)
            '74'

            sage: k = sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro(3^4,'a', repr='log')
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

        sage: for x in GF(2^2,'a'): print x
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
            sage: i.next()
            0
            sage: i.next()
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
            raise ValueError, "must be a perfect square."

    cpdef ModuleElement _add_(self, ModuleElement right):
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

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        Add two elements inplace.

        EXAMPLES::

            sage: k.<b> = GF(9**2)
            sage: b^10 + 2*b # indirect doctest
            2*b^3 + 2*b^2 + 2*b + 1
        """
        cdef int r
        self.element = self._cache.objectptr.add(r, self.element ,
                                                         (<FiniteField_givaroElement>right).element )
        return self


    cpdef RingElement _mul_(self, RingElement right):
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


    cpdef RingElement _imul_(self, RingElement right):
        """
        Multiply two elements inplace.

        EXAMPLES::

            sage: k.<c> = GF(7**4)
            sage: 3*c # indirect doctest
            3*c
            sage: c*c
            c^2
        """
        cdef int r
        self.element = self._cache.objectptr.mul(r, self.element,
                                                         (<FiniteField_givaroElement>right).element)
        return self

    cpdef RingElement _div_(self, RingElement right):
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
            raise ZeroDivisionError, 'division by zero in finite field.'
        r = self._cache.objectptr.div(r, self.element,
                                              (<FiniteField_givaroElement>right).element)
        return make_FiniteField_givaroElement(self._cache,r)

    cpdef RingElement _idiv_(self, RingElement right):
        """
        Divide two elements inplace

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
            raise ZeroDivisionError, 'division by zero in finite field.'
        self.element = self._cache.objectptr.div(r, self.element,
                                                         (<FiniteField_givaroElement>right).element)
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
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

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        Subtract two elements inplace.

        EXAMPLES::

            sage: k.<a> = GF(3**4)
            sage: k(3) - k(1) # indirect doctest
            2
            sage: 2*a - a^2
            2*a^2 + 2*a
        """
        cdef int r
        self.element = self._cache.objectptr.sub(r, self.element,
                                                         (<FiniteField_givaroElement>right).element)
        return self

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
                raise ValueError, "exponent must be an integer"
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

    def __richcmp__(left, right, int op):
        """
        EXAMPLES::

            sage: k.<a> = GF(9); k
            Finite Field in a of size 3^2
            sage: a == k('a') # indirect doctest
            True
            sage: a == a + 1
            False
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
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
        prime subfield, the integer returned is equal to ``self`` and not
        to ``log_repr``.

        Elements of this field are represented as ints in as follows:
        for `e \in \GF{p}[x]` with `e = a_0 + a_1x + a_2x^2 + \cdots`, `e` is
        represented as: `n= a_0 + a_1  p + a_2  p^2 + \cdots`.

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
        raise TypeError, "not in prime subfield"

    def log_to_int(FiniteField_givaroElement self):
        r"""
        Returns the int representation of ``self``, as a Sage integer.   Use
        ``int(self)`` to directly get a Python int.

        Elements of this field are represented as ints in as follows:
        for `e \in \GF{p}[x]` with `e = a_0 + a_1x + a_2x^2 + \cdots`, `e` is
        represented as: `n = a_0 + a_1  p + a_2  p^2 + \cdots`.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: k(4).log_to_int()
            4
            sage: b.log_to_int()
            5
            sage: type(b.log_to_int())
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(self._cache.log_to_int(self.element))

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

    def int_repr(FiniteField_givaroElement self):
        r"""
        Return the string representation of ``self`` as an int (as returned
        by :meth:`log_to_int`).

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+1).int_repr()
            '6'
        """
        return self._cache._element_int_repr(self)

    def log_repr(FiniteField_givaroElement self):
        r"""
        Return the log representation of ``self`` as a string.  See the
        documentation of the ``_element_log_repr`` function of the
        parent field.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+2).log_repr()
            '15'
        """
        return self._cache._element_log_repr(self)

    def poly_repr(FiniteField_givaroElement self):
        r"""
        Return representation of this finite field element as a polynomial
        in the generator.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: (b+2).poly_repr()
            'b + 2'
        """
        return self._cache._element_poly_repr(self)

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

    def _finite_field_ext_pari_element(FiniteField_givaroElement self, k=None):
        """
        Return an element of ``k`` supposed to match this element.

        .. WARNING::

            No checks if ``k == self.parent()`` are performed.

        INPUT:

        - ``k`` -- (optional) :class:`FiniteField_ext_pari`

        OUTPUT:

        ``k.gen()^(self.log_repr())``

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.charpoly('x')
            x^2 + 4*x + 2
            sage: P = S._finite_field_ext_pari_(); type(P)
            <class 'sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari_with_category'>
            sage: c = b._finite_field_ext_pari_element(P); c
            b
            sage: type(c)
            <class 'sage.rings.finite_rings.element_ext_pari.FiniteField_ext_pariElement_with_category'>
            sage: c.charpoly('x')
            x^2 + 4*x + 2

        The PARI field is automatically determined if it is not given::

            sage: d = b._finite_field_ext_pari_element(); d
            b
            sage: type(d)
            <class 'sage.rings.finite_rings.element_ext_pari.FiniteField_ext_pariElement_with_category'>
        """
        if k is None:
            k = self.parent()._finite_field_ext_pari_()
        elif not isinstance(k, finite_field_ext_pari.FiniteField_ext_pari):
            raise TypeError, "k must be a pari finite field."

        variable = k.gen()._pari_()

        quo = self.integer_representation()
        b   = int(self._cache.characteristic())

        ret = k._pari_one() - k._pari_one()    # TODO -- weird
        i = 0
        while quo!=0:
            coeff = quo%b
            if coeff != 0:
                ret = coeff * variable ** i + ret
            quo = quo/b
            i = i+1
        return k(ret)

    def _pari_(FiniteField_givaroElement self, var=None):
        r"""
        Return PARI representation of this finite field element.

        INPUT:

        - ``var`` -- (default: ``None``) optional variable string

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: a._pari_()
            Mod(Mod(1, 5)*a, Mod(1, 5)*a^3 + Mod(3, 5)*a + Mod(3, 5))

            sage: a._pari_('b')
            Mod(Mod(1, 5)*b, Mod(1, 5)*b^3 + Mod(3, 5)*b + Mod(3, 5))

            sage: t = 3*a^2 + 2*a + 4
            sage: t_string = t._pari_init_('y')
            sage: t_string
            'Mod(Mod(3, 5)*y^2 + Mod(2, 5)*y + Mod(4, 5), Mod(1, 5)*y^3 + Mod(3, 5)*y + Mod(3, 5))'
            sage: type(t_string)
            <type 'str'>
            sage: t_element = t._pari_('b')
            sage: t_element
            Mod(Mod(3, 5)*b^2 + Mod(2, 5)*b + Mod(4, 5), Mod(1, 5)*b^3 + Mod(3, 5)*b + Mod(3, 5))
            sage: t_element.parent()
            Interface to the PARI C library
        """
        return pari(self._pari_init_(var))

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
        import sage.rings.arith

        if self._multiplicative_order is not None:
            return self._multiplicative_order
        else:
            if self.is_zero():
                raise ArithmeticError("Multiplicative order of 0 not defined.")
            n = (self._cache).order_c() - 1
            order = 1
            for p, e in sage.rings.arith.factor(n):
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
            raise NotImplementedError, "conversion of (Givaro) finite field element to GAP not implemented except for fields defined by Conway polynomials."
        if cache.order_c() > 65536:
            raise TypeError, "order (=%s) must be at most 65536."%F.order_c()
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
        return hash(self.log_to_int())

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
        if PY_TYPE_CHECK(reverse, Parent):
            raise TypeError, "Base field is fixed to prime subfield."
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
        y = PY_NEW(FiniteField_givaroElement)
        y._parent = <Parent> cache.parent
        y._cache = cache
        y.element = x
        return y

