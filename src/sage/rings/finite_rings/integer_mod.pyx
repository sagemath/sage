r"""
Elements of `\ZZ/n\ZZ`

An element of the integers modulo `n`.

There are three types of integer_mod classes, depending on the
size of the modulus.


-  ``IntegerMod_int`` stores its value in a
   ``int_fast32_t`` (typically an ``int``);
   this is used if the modulus is less than
   `\sqrt{2^{31}-1}`.

-  ``IntegerMod_int64`` stores its value in a
   ``int_fast64_t`` (typically a ``long
   long``); this is used if the modulus is less than
   `2^{31}-1`.

-  ``IntegerMod_gmp`` stores its value in a
   ``mpz_t``; this can be used for an arbitrarily large
   modulus.


All extend ``IntegerMod_abstract``.

For efficiency reasons, it stores the modulus (in all three forms,
if possible) in a common (cdef) class
``NativeIntStruct`` rather than in the parent.

AUTHORS:

-  Robert Bradshaw: most of the work

-  Didier Deshommes: bit shifting

-  William Stein: editing and polishing; new arith architecture

-  Robert Bradshaw: implement native is_square and square_root

-  William Stein: sqrt

-  Maarten Derickx: moved the valuation code from the global
   valuation function to here


TESTS::

    sage: R = Integers(101^3)
    sage: a = R(824362); b = R(205942)
    sage: a * b
    851127

    sage: type(IntegerModRing(2^31-1).an_element())
    <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
    sage: type(IntegerModRing(2^31).an_element())
    <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
"""

#################################################################################
#       Copyright (C) 2006 Robert Bradshaw <robertwb@math.washington.edu>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/stdsage.pxi"
from cpython.int cimport *
from cpython.list cimport *
from cpython.ref cimport *

cdef extern from "math.h":
    double log(double)
    int ceil(double)

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pyintlong(mpz_t src)

import operator

cdef bint use_32bit_type(int_fast64_t modulus):
    return modulus <= INTEGER_MOD_INT32_LIMIT

## import arith
import sage.rings.rational as rational
from sage.libs.pari.all import pari, PariError
import sage.rings.integer_ring as integer_ring

import sage.rings.commutative_ring_element as commutative_ring_element
import sage.interfaces.all

import sage.rings.integer
import sage.rings.integer_ring
cimport sage.rings.integer
from sage.rings.integer cimport Integer

import sage.structure.element
cimport sage.structure.element
coerce_binop = sage.structure.element.coerce_binop
from sage.structure.element cimport RingElement, ModuleElement, Element
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map

from sage.structure.sage_object import register_unpickle_override

#from sage.structure.parent cimport Parent

cdef Integer one_Z = Integer(1)

def Mod(n, m, parent=None):
    """
    Return the equivalence class of `n` modulo `m` as
    an element of `\ZZ/m\ZZ`.

    EXAMPLES::

        sage: x = Mod(12345678, 32098203845329048)
        sage: x
        12345678
        sage: x^100
        1017322209155072

    You can also use the lowercase version::

        sage: mod(12,5)
        2

    Illustrates that trac #5971 is fixed. Consider `n` modulo `m` when
    `m = 0`. Then `\ZZ/0\ZZ` is isomorphic to `\ZZ` so `n` modulo `0` is
    is equivalent to `n` for any integer value of `n`::

        sage: Mod(10, 0)
        10
        sage: a = randint(-100, 100)
        sage: Mod(a, 0) == a
        True
    """
    # when m is zero, then ZZ/0ZZ is isomorphic to ZZ
    if m == 0:
        return n

    # m is non-zero, so return n mod m
    cdef IntegerMod_abstract x
    import integer_mod_ring
    x = IntegerMod(integer_mod_ring.IntegerModRing(m), n)
    if parent is None:
        return x
    x._parent = parent
    return x


mod = Mod

register_unpickle_override('sage.rings.integer_mod', 'Mod', Mod)
register_unpickle_override('sage.rings.integer_mod', 'mod', mod)

def IntegerMod(parent, value):
    """
    Create an integer modulo `n` with the given parent.

    This is mainly for internal use.
    """
    cdef NativeIntStruct modulus
    cdef Py_ssize_t res
    modulus = parent._pyx_order
    if modulus.table is not None:
        if PY_TYPE_CHECK(value, sage.rings.integer.Integer) or PY_TYPE_CHECK(value, int) or PY_TYPE_CHECK(value, long):
            res = value % modulus.int64
            if res < 0:
                res = res + modulus.int64
            a = modulus.lookup(res)
            if (<Element>a)._parent is not parent:
               (<Element>a)._parent = parent
#                print (<Element>a)._parent, " is not ", parent
            return a
    if modulus.int32 != -1:
        return IntegerMod_int(parent, value)
    elif modulus.int64 != -1:
        return IntegerMod_int64(parent, value)
    else:
        return IntegerMod_gmp(parent, value)

def is_IntegerMod(x):
    """
    Return ``True`` if and only if x is an integer modulo
    `n`.

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import is_IntegerMod
        sage: is_IntegerMod(5)
        False
        sage: is_IntegerMod(Mod(5,10))
        True
    """
    return PY_TYPE_CHECK(x, IntegerMod_abstract)

def makeNativeIntStruct(sage.rings.integer.Integer z):
    """
    Function to convert a Sage Integer into class NativeIntStruct.

    .. note::

       This function is only used for the unpickle override below.
    """
    return NativeIntStruct(z)

register_unpickle_override('sage.rings.integer_mod', 'makeNativeIntStruct', makeNativeIntStruct)

cdef class NativeIntStruct:
    """
    We store the various forms of the modulus here rather than in the
    parent for efficiency reasons.

    We may also store a cached table of all elements of a given ring in
    this class.
    """
    def __init__(NativeIntStruct self, sage.rings.integer.Integer z):
        self.int64 = -1
        self.int32 = -1
        self.table = None # NULL
        self.sageInteger = z
        if mpz_cmp_si(z.value, INTEGER_MOD_INT64_LIMIT) <= 0:
            self.int64 = mpz_get_si(z.value)
            if use_32bit_type(self.int64):
                self.int32 = self.int64

    def __reduce__(NativeIntStruct self):
        return sage.rings.finite_rings.integer_mod.makeNativeIntStruct, (self.sageInteger, )

    def precompute_table(NativeIntStruct self, parent, inverses=True):
        """
        Function to compute and cache all elements of this class.

        If ``inverses == True``, also computes and caches the inverses
        of the invertible elements.

        EXAMPLES:

        This is used by the :class:`sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic` constructor::

            sage: from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
            sage: R = IntegerModRing_generic(39, cache=False)
            sage: R(5)^-1
            8
            sage: R(5)^-1 is R(8)
            False
            sage: R = IntegerModRing_generic(39, cache=True)  # indirect doctest
            sage: R(5)^-1 is R(8)
            True

        Check that the inverse of 0 modulo 1 works, see :trac:`13639`::

            sage: R = IntegerModRing_generic(1, cache=True)  # indirect doctest
            sage: R(0)^-1 is R(0)
            True
        """
        self.table = PyList_New(self.int64)
        cdef Py_ssize_t i
        if self.int32 != -1:
            for i from 0 <= i < self.int32:
                z = IntegerMod_int(parent, i)
                Py_INCREF(z); PyList_SET_ITEM(self.table, i, z)
        else:
            for i from 0 <= i < self.int64:
                z = IntegerMod_int64(parent, i)
                Py_INCREF(z); PyList_SET_ITEM(self.table, i, z)

        if inverses:
            if self.int64 == 1:
                # Special case for integers modulo 1
                self.inverses = self.table
            else:
                tmp = [None] * self.int64
                for i from 1 <= i < self.int64:
                    try:
                        tmp[i] = ~self.table[i]
                    except ZeroDivisionError:
                        pass
                self.inverses = tmp

    def _get_table(self):
        return self.table

    cdef lookup(NativeIntStruct self, Py_ssize_t value):
        return <object>PyList_GET_ITEM(self.table, value)


cdef class IntegerMod_abstract(FiniteRingElement):

    def __init__(self, parent):
        """
        EXAMPLES::

            sage: a = Mod(10,30^10); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        self._parent = parent
        self.__modulus = parent._pyx_order


    cdef _new_c_from_long(self, long value):
        cdef IntegerMod_abstract x
        x = <IntegerMod_abstract>PY_NEW(<object>PY_TYPE(self))
        if PY_TYPE_CHECK(x, IntegerMod_gmp):
            mpz_init((<IntegerMod_gmp>x).value) # should be done by the new method
        x._parent = self._parent
        x.__modulus = self.__modulus
        x.set_from_long(value)
        return x

    cdef void set_from_mpz(self, mpz_t value):
        raise NotImplementedError, "Must be defined in child class."

    cdef void set_from_long(self, long value):
        raise NotImplementedError, "Must be defined in child class."

    def __abs__(self):
        """
        Raise an error message, since ``abs(x)`` makes no sense
        when ``x`` is an integer modulo `n`.

        EXAMPLES::

            sage: abs(Mod(2,3))
            Traceback (most recent call last):
            ...
            ArithmeticError: absolute valued not defined on integers modulo n.
        """
        raise ArithmeticError, "absolute valued not defined on integers modulo n."

    def __reduce__(IntegerMod_abstract self):
        """
        EXAMPLES::

            sage: a = Mod(4,5); a
            4
            sage: loads(a.dumps()) == a
            True
            sage: a = Mod(-1,5^30)^25;
            sage: loads(a.dumps()) == a
            True
        """
        return sage.rings.finite_rings.integer_mod.mod, (self.lift(), self.modulus(), self.parent())

    def is_nilpotent(self):
        r"""
        Return ``True`` if ``self`` is nilpotent,
        i.e., some power of ``self`` is zero.

        EXAMPLES::

            sage: a = Integers(90384098234^3)
            sage: factor(a.order())
            2^3 * 191^3 * 236607587^3
            sage: b = a(2*191)
            sage: b.is_nilpotent()
            False
            sage: b = a(2*191*236607587)
            sage: b.is_nilpotent()
            True

        ALGORITHM: Let `m \geq  \log_2(n)`, where `n` is
        the modulus. Then `x \in \ZZ/n\ZZ` is
        nilpotent if and only if `x^m = 0`.

        PROOF: This is clear if you reduce to the prime power case, which
        you can do via the Chinese Remainder Theorem.

        We could alternatively factor `n` and check to see if the
        prime divisors of `n` all divide `x`. This is
        asymptotically slower :-).
        """
        if self.is_zero():
            return True
        m = self.__modulus.sageInteger.exact_log(2) + 1
        return (self**m).is_zero()

    #################################################################
    # Interfaces
    #################################################################
    def _pari_init_(self):
        return 'Mod(%s,%s)'%(str(self), self.__modulus.sageInteger)

    def _pari_(self):
        return self.lift()._pari_().Mod(self.__modulus.sageInteger)

    def _gap_init_(self):
        r"""
        Return string representation of corresponding GAP object.

        EXAMPLES::

            sage: a = Mod(2,19)
            sage: gap(a)
            Z(19)
            sage: gap(Mod(3, next_prime(10000)))
            Z(10007)^6190
            sage: gap(Mod(3, next_prime(100000)))
            ZmodpZObj( 3, 100003 )
            sage: gap(Mod(4, 48))
            ZmodnZObj( 4, 48 )
        """
        return '%s*One(ZmodnZ(%s))' % (self, self.__modulus.sageInteger)

    def _magma_init_(self, magma):
        """
        Coercion to Magma.

        EXAMPLES::

            sage: a = Integers(15)(4)
            sage: b = magma(a)                # optional - magma
            sage: b.Type()                    # optional - magma
            RngIntResElt
            sage: b^2                         # optional - magma
            1
        """
        return '%s!%s'%(self.parent()._magma_init_(magma), self)

    def _axiom_init_(self):
        """
        Return a string representation of the corresponding to
        (Pan)Axiom object.

        EXAMPLES::

            sage: a = Integers(15)(4)
            sage: a._axiom_init_()
            '4 :: IntegerMod(15)'

            sage: aa = axiom(a); aa #optional - axiom
            4
            sage: aa.type()         #optional - axiom
            IntegerMod 15

            sage: aa = fricas(a); aa #optional - fricas
            4
            sage: aa.type()          #optional - fricas
            IntegerMod(15)

        """
        return '%s :: %s'%(self, self.parent()._axiom_init_())

    _fricas_init_ = _axiom_init_

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: K = GF(7)
            sage: sage_input(K(5), verify=True)
            # Verified
            GF(7)(5)
            sage: sage_input(K(5) * polygen(K), verify=True)
            # Verified
            R.<x> = GF(7)[]
            5*x
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: K(5)._sage_input_(SageInputBuilder(), False)
            {call: {call: {atomic:GF}({atomic:7})}({atomic:5})}
            sage: K(5)._sage_input_(SageInputBuilder(), True)
            {atomic:5}
        """
        v = sib.int(self.lift())
        if coerced:
            return v
        else:
            return sib(self.parent())(v)

    def log(self, b=None):
        r"""
        Return an integer `x` such that `b^x = a`, where
        `a` is ``self``.

        INPUT:


        -  ``self`` - unit modulo `n`

        -  ``b`` - a unit modulo `n`. If ``b`` is not given,
           ``R.multiplicative_generator()`` is used, where
           ``R`` is the parent of ``self``.


        OUTPUT: Integer `x` such that `b^x = a`, if this exists; a ValueError otherwise.

        .. note::

           If the modulus is prime and b is a generator, this calls Pari's ``znlog``
           function, which is rather fast. If not, it falls back on the generic
           discrete log implementation in :meth:`sage.groups.generic.discrete_log`.

        EXAMPLES::

            sage: r = Integers(125)
            sage: b = r.multiplicative_generator()^3
            sage: a = b^17
            sage: a.log(b)
            17
            sage: a.log()
            51

        A bigger example::

            sage: FF = FiniteField(2^32+61)
            sage: c = FF(4294967356)
            sage: x = FF(2)
            sage: a = c.log(x)
            sage: a
            2147483678
            sage: x^a
            4294967356

        Things that can go wrong. E.g., if the base is not a generator for
        the multiplicative group, or not even a unit.

        ::

            sage: Mod(3, 7).log(Mod(2, 7))
            Traceback (most recent call last):
            ...
            ValueError: No discrete log of 3 found to base 2
            sage: a = Mod(16, 100); b = Mod(4,100)
            sage: a.log(b)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        We check that #9205 is fixed::

            sage: Mod(5,9).log(Mod(2, 9))
            5

        We test against a bug (side effect on PARI) fixed in #9438::

            sage: R.<a, b> = QQ[]
            sage: pari(b)
            b
            sage: GF(7)(5).log()
            5
            sage: pari(b)
            b

        AUTHORS:

        - David Joyner and William Stein (2005-11)

        - William Stein (2007-01-27): update to use PARI as requested
          by David Kohel.

        - Simon King (2010-07-07): fix a side effect on PARI
        """
        if b is None:
            b = self._parent.multiplicative_generator()
        else:
            b = self._parent(b)

        if self.modulus().is_prime() and b.multiplicative_order() == b.parent().unit_group_order():

            # use PARI

            cmd = 'if(znorder(Mod(%s,%s))!=eulerphi(%s),-1,znlog(%s,Mod(%s,%s)))'%(b, self.__modulus.sageInteger,
                                                      self.__modulus.sageInteger,
                                             self, b, self.__modulus.sageInteger)
            try:
                n = Integer(pari(cmd))
                return n
            except PariError, msg:
                raise ValueError, "%s\nPARI failed to compute discrete log (perhaps base is not a generator or is too large)"%msg

        else: # fall back on slower native implementation

            from sage.groups.generic import discrete_log
            return discrete_log(self, b)

    def generalised_log(self):
        r"""
        Return integers `n_i` such that

        ..math::

            \prod_i x_i^{n_i} = \text{self},

        where `x_1, \dots, x_d` are the generators of the unit group
        returned by ``self.parent().unit_gens()``. See also :meth:`log`.

        EXAMPLES::

            sage: m = Mod(3, 1568)
            sage: v = m.generalised_log(); v
            [1, 3, 1]
            sage: prod([Zmod(1568).unit_gens()[i] ** v[i] for i in [0..2]])
            3

        """
        if not self.is_unit():
            raise ZeroDivisionError
        N = self.modulus()
        h = []
        for (p, c) in N.factor():
            if p != 2 or (p == 2 and c == 2):
                h.append((self % p**c).log())
            elif c > 2:
                m = self % p**c
                if m % 4 == 1:
                    h.append(0)
                else:
                    h.append(1)
                    m *= -1
                h.append(m.log(5))
        return h

    def modulus(IntegerMod_abstract self):
        """
        EXAMPLES::

            sage: Mod(3,17).modulus()
            17
        """
        return self.__modulus.sageInteger

    def charpoly(self, var='x'):
        """
        Returns the characteristic polynomial of this element.

        EXAMPLES::

            sage: k = GF(3)
            sage: a = k.gen()
            sage: a.charpoly('x')
            x + 2
            sage: a + 2
            0

        AUTHORS:

        - Craig Citro
        """
        R = self.parent()[var]
        return R([-self,1])

    def minpoly(self, var='x'):
        """
        Returns the minimal polynomial of this element.

        EXAMPLES:
            sage: GF(241, 'a')(1).minpoly()
            x + 240
        """
        return self.charpoly(var)

    def minimal_polynomial(self, var='x'):
        """
        Returns the minimal polynomial of this element.

        EXAMPLES:
            sage: GF(241, 'a')(1).minimal_polynomial(var = 'z')
            z + 240
        """
        return self.minpoly(var)

    def polynomial(self, var='x'):
        """
        Returns a constant polynomial representing this value.

        EXAMPLES::

            sage: k = GF(7)
            sage: a = k.gen(); a
            1
            sage: a.polynomial()
            1
            sage: type(a.polynomial())
            <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>
        """
        R = self.parent()[var]
        return R(self)

    def norm(self):
        """
        Returns the norm of this element, which is itself. (This is here
        for compatibility with higher order finite fields.)

        EXAMPLES::

            sage: k = GF(691)
            sage: a = k(389)
            sage: a.norm()
            389

        AUTHORS:

        - Craig Citro
        """
        return self

    def trace(self):
        """
        Returns the trace of this element, which is itself. (This is here
        for compatibility with higher order finite fields.)

        EXAMPLES::

            sage: k = GF(691)
            sage: a = k(389)
            sage: a.trace()
            389

        AUTHORS:

        - Craig Citro
        """
        return self

    def centerlift(self):
        r"""
        Lift ``self`` to an integer `i` such that `n/2 < i <= n/2`
        (where `n` denotes the modulus).

        EXAMPLES::

            sage: Mod(0,5).centerlift()
            0
            sage: Mod(1,5).centerlift()
            1
            sage: Mod(2,5).centerlift()
            2
            sage: Mod(3,5).centerlift()
            -2
            sage: Mod(4,5).centerlift()
            -1
            sage: Mod(50,100).centerlift()
            50
            sage: Mod(51,100).centerlift()
            -49
            sage: Mod(-1,3^100).centerlift()
            -1
        """
        n = self.modulus()
        x = self.lift()
        if 2*x <= n:
            return x
        else:
            return x - n

    cpdef bint is_one(self):
        raise NotImplementedError

    cpdef bint is_unit(self):
        raise NotImplementedError

    def is_square(self):
        r"""
        EXAMPLES::

            sage: Mod(3,17).is_square()
            False
            sage: Mod(9,17).is_square()
            True
            sage: Mod(9,17*19^2).is_square()
            True
            sage: Mod(-1,17^30).is_square()
            True
            sage: Mod(1/9, next_prime(2^40)).is_square()
            True
            sage: Mod(1/25, next_prime(2^90)).is_square()
            True

        TESTS::

            sage: Mod(1/25, 2^8).is_square()
            True
            sage: Mod(1/25, 2^40).is_square()
            True

        ALGORITHM: Calculate the Jacobi symbol
        `(\mathtt{self}/p)` at each prime `p`
        dividing `n`. It must be 1 or 0 for each prime, and if it
        is 0 mod `p`, where `p^k || n`, then
        `ord_p(\mathtt{self})` must be even or greater than
        `k`.

        The case `p = 2` is handled separately.

        AUTHORS:

        - Robert Bradshaw
        """
        return self.is_square_c()

    cdef bint is_square_c(self) except -2:
        if self.is_zero() or self.is_one():
            return 1
        moduli = self.parent().factored_order()
        cdef int val, e
        lift = self.lift()
        if len(moduli) == 1:
            p, e = moduli[0]
            if e == 1:
                return lift.jacobi(p) != -1
            elif p == 2:
                return self._pari_().issquare() # TODO: implement directly
            elif self % p == 0:
                val = lift.valuation(p)
                return val >= e or (val % 2 == 0 and (lift // p**val).jacobi(p) != -1)
            else:
                return lift.jacobi(p) != -1
        else:
            for p, e in moduli:
                if p == 2:
                    if e > 1 and not self._pari_().issquare(): # TODO: implement directly
                        return 0
                elif e > 1 and lift % p == 0:
                    val = lift.valuation(p)
                    if val < e and (val % 2 == 1 or (lift // p**val).jacobi(p) == -1):
                        return 0
                elif lift.jacobi(p) == -1:
                    return 0
            return 1

    def sqrt(self, extend=True, all=False):
        r"""
        Returns square root or square roots of ``self`` modulo
        `n`.

        INPUT:


        -  ``extend`` - bool (default: ``True``);
           if ``True``, return a square root in an extension ring,
           if necessary. Otherwise, raise a ``ValueError`` if the
           square root is not in the base ring.

        -  ``all`` - bool (default: ``False``); if
           ``True``, return {all} square roots of self, instead of
           just one.


        ALGORITHM: Calculates the square roots mod `p` for each of
        the primes `p` dividing the order of the ring, then lifts
        them `p`-adically and uses the CRT to find a square root
        mod `n`.

        See also ``square_root_mod_prime_power`` and
        ``square_root_mod_prime`` (in this module) for more
        algorithmic details.

        EXAMPLES::

            sage: mod(-1, 17).sqrt()
            4
            sage: mod(5, 389).sqrt()
            86
            sage: mod(7, 18).sqrt()
            5
            sage: a = mod(14, 5^60).sqrt()
            sage: a*a
            14
            sage: mod(15, 389).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square
            sage: Mod(1/9, next_prime(2^40)).sqrt()^(-2)
            9
            sage: Mod(1/25, next_prime(2^90)).sqrt()^(-2)
            25

        ::

            sage: a = Mod(3,5); a
            3
            sage: x = Mod(-1, 360)
            sage: x.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square
            sage: y = x.sqrt(); y
            sqrt359
            sage: y.parent()
            Univariate Quotient Polynomial Ring in sqrt359 over Ring of integers modulo 360 with modulus x^2 + 1
            sage: y^2
            359

        We compute all square roots in several cases::

            sage: R = Integers(5*2^3*3^2); R
            Ring of integers modulo 360
            sage: R(40).sqrt(all=True)
            [20, 160, 200, 340]
            sage: [x for x in R if x^2 == 40]  # Brute force verification
            [20, 160, 200, 340]
            sage: R(1).sqrt(all=True)
            [1, 19, 71, 89, 91, 109, 161, 179, 181, 199, 251, 269, 271, 289, 341, 359]
            sage: R(0).sqrt(all=True)
            [0, 60, 120, 180, 240, 300]

        ::

            sage: R = Integers(5*13^3*37); R
            Ring of integers modulo 406445
            sage: v = R(-1).sqrt(all=True); v
            [78853, 111808, 160142, 193097, 213348, 246303, 294637, 327592]
            sage: [x^2 for x in v]
            [406444, 406444, 406444, 406444, 406444, 406444, 406444, 406444]
            sage: v = R(169).sqrt(all=True); min(v), -max(v), len(v)
            (13, 13, 104)
            sage: all([x^2==169 for x in v])
            True

        ::

            sage: t = FiniteField(next_prime(2^100))(4)
            sage: t.sqrt(extend = False, all = True)
            [2, 1267650600228229401496703205651]
            sage: t = FiniteField(next_prime(2^100))(2)
            sage: t.sqrt(extend = False, all = True)
            []

        Modulo a power of 2::

            sage: R = Integers(2^7); R
            Ring of integers modulo 128
            sage: a = R(17)
            sage: a.sqrt()
            23
            sage: a.sqrt(all=True)
            [23, 41, 87, 105]
            sage: [x for x in R if x^2==17]
            [23, 41, 87, 105]
        """
        if self.is_one():
            if all:
                return list(self.parent().square_roots_of_one())
            else:
                return self

        if not self.is_square_c():
            if extend:
                y = 'sqrt%s'%self
                R = self.parent()['x']
                modulus = R.gen()**2 - R(self)
                if self._parent.is_field():
                    import constructor
                    Q = constructor.FiniteField(self.__modulus.sageInteger**2, y, modulus)
                else:
                    R = self.parent()['x']
                    Q = R.quotient(modulus, names=(y,))
                z = Q.gen()
                if all:
                    # TODO
                    raise NotImplementedError
                return z
            if all:
                return []
            raise ValueError, "self must be a square"

        F = self._parent.factored_order()
        cdef long e, exp, val
        if len(F) == 1:
            p, e = F[0]

            if all and e > 1 and not self.is_unit():
                if self.is_zero():
                    # All multiples of p^ciel(e/2) vanish
                    return [self._parent(x) for x in xrange(0, self.__modulus.sageInteger, p**((e+1)/2))]
                else:
                    z = self.lift()
                    val = z.valuation(p)/2  # square => valuation is even
                    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
                    # Find the unit part (mod the ring with appropriate precision)
                    u = IntegerModRing(p**(e-val))(z // p**(2*val))
                    # will add multiples of p^exp
                    exp = e - val
                    if p == 2:
                        exp -= 1  # note the factor of 2 below
                    if 2*exp < e:
                        exp = (e+1)/2
                    # For all a^2 = u and all integers b
                    #   (a*p^val + b*p^exp) ^ 2
                    #   = u*p^(2*val) + 2*a*b*p^(val+exp) + b^2*p^(2*exp)
                    #   = u*p^(2*val)  mod p^e
                    # whenever min(val+exp, 2*exp) > e
                    p_val = p**val
                    p_exp = p**exp
                    w = [self._parent(a.lift() * p_val + b)
                            for a in u.sqrt(all=True)
                            for b in xrange(0, self.__modulus.sageInteger, p_exp)]
                    if p == 2:
                        w = list(set(w))
                    w.sort()
                    return w

            if e > 1:
                x = square_root_mod_prime_power(mod(self, p**e), p, e)
            else:
                x = square_root_mod_prime(self, p)
            x = x._balanced_abs()

            if not all:
                return x

            v = list(set([x*a for a in self._parent.square_roots_of_one()]))
            v.sort()
            return v

        else:
            if not all:
                # Use CRT to combine together a square root modulo each prime power
                sqrts = [square_root_mod_prime(mod(self, p), p) for p, e in F if e == 1] + \
                        [square_root_mod_prime_power(mod(self, p**e), p, e) for p, e in F if e != 1]

                x = sqrts.pop()
                for y in sqrts:
                    x = x.crt(y)
                return x._balanced_abs()
            else:
                # Use CRT to combine together all square roots modulo each prime power
                vmod = []
                moduli = []
                P = self.parent()
                from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
                for p, e in F:
                    k = p**e
                    R = IntegerModRing(p**e)
                    w = [P(x) for x in R(self).sqrt(all=True)]
                    vmod.append(w)
                    moduli.append(k)
                # Now combine in all possible ways using the CRT
                from sage.rings.arith import CRT_basis
                basis = CRT_basis(moduli)
                from sage.misc.mrange import cartesian_product_iterator
                v = []
                for x in cartesian_product_iterator(vmod):
                    # x is a specific choice of roots modulo each prime power divisor
                    a = sum([basis[i]*x[i] for i in range(len(x))])
                    v.append(a)
                v.sort()
                return v

    square_root = sqrt

    def nth_root(self, n, extend = False, all = False, algorithm = None, cunningham = False):
        r"""
        Returns an `n`\th root of ``self``.

        INPUT:

        - ``n`` - integer `\geq 1`

        - ``extend`` - bool (default: True); if True, return an nth
          root in an extension ring, if necessary. Otherwise, raise a
          ValueError if the root is not in the base ring.  Warning:
          this option is not implemented!

        - ``all`` - bool (default: ``False``); if ``True``, return all `n`\th
          roots of ``self``, instead of just one.

        - ``algorithm`` - string (default: None); The algorithm for the prime modulus case.
          CRT and p-adic log techniques are used to reduce to this case.
          'Johnston' is the only currently supported option.

        - ``cunningham`` - bool (default: ``False``); In some cases,
          factorization of ``n`` is computed. If cunningham is set to ``True``,
          the factorization of ``n`` is computed using trial division for all
          primes in the so called Cunningham table. Refer to
          sage.rings.factorint.factor_cunningham for more information. You need
          to install an optional package to use this method, this can be done
          with the following command line ``sage -i cunningham_tables``

        OUTPUT:

        If self has an `n`\th root, returns one (if ``all`` is ``False``) or a
        list of all of them (if ``all`` is ``True``).  Otherwise, raises a
        ``ValueError`` (if ``extend`` is ``False``) or a ``NotImplementedError`` (if
        ``extend`` is ``True``).

        .. warning::
           The 'extend' option is not implemented (yet).

        NOTES:

        - If `n = 0`:

          - if ``all=True``:

            - if ``self=1``: all nonzero elements of the parent are returned in
              a list.  Note that this could be very expensive for large
              parents.

            - otherwise: an empty list is returned

          - if ``all=False``:

            - if ``self=1``: ``self`` is returned

            - otherwise; a ``ValueError`` is raised

        - If `n < 0`:

          - if self is invertible, the `(-n)`\th root of the inverse of self is returned

          - otherwise a ``ValueError`` is raised or empty list returned.

        EXAMPLES::


            sage: K = GF(31)
            sage: a = K(22)
            sage: K(22).nth_root(7)
            13
            sage: K(25).nth_root(5)
            5
            sage: K(23).nth_root(3)
            29
            sage: mod(225,2^5*3^2).nth_root(4, all=True)
            [225, 129, 33, 63, 255, 159, 9, 201, 105, 279, 183, 87, 81, 273, 177, 207, 111, 15, 153, 57, 249, 135, 39, 231]
            sage: mod(275,2^5*7^4).nth_root(7, all=True)
            [58235, 25307, 69211, 36283, 3355, 47259, 14331]
            sage: mod(1,8).nth_root(2,all=True)
            [1, 7, 5, 3]
            sage: mod(4,8).nth_root(2,all=True)
            [2, 6]
            sage: mod(1,16).nth_root(4,all=True)
            [1, 15, 13, 3, 9, 7, 5, 11]
            sage: (mod(22,31)^200).nth_root(200)
            5
            sage: mod(3,6).nth_root(0,all=True)
            []
            sage: mod(3,6).nth_root(0)
            Traceback (most recent call last):
            ...
            ValueError
            sage: mod(1,6).nth_root(0,all=True)
            [1, 2, 3, 4, 5]

        TESTS::

            sage: for p in [1009,2003,10007,100003]:
            ...       K = GF(p)
            ...       for r in (p-1).divisors():
            ...           if r == 1: continue
            ...           x = K.random_element()
            ...           y = x^r
            ...           if y.nth_root(r)**r != y: raise RuntimeError
            ...           if (y^41).nth_root(41*r)**(41*r) != y^41: raise RuntimeError
            ...           if (y^307).nth_root(307*r)**(307*r) != y^307: raise RuntimeError

            sage: for t in xrange(200):
            ...       n = randint(1,2^63)
            ...       K = Integers(n)
            ...       b = K.random_element()
            ...       e = randint(-2^62, 2^63)
            ...       try:
            ...           a = b.nth_root(e)
            ...           if a^e != b:
            ...               print n, b, e, a
            ...               raise NotImplementedError
            ...       except ValueError:
            ...           pass

        We check that #13172 is resolved::

            sage: mod(-1, 4489).nth_root(2, all=True)
            []

        Check that the code path cunningham might be used::

            sage: a = Mod(9,11)
            sage: a.nth_root(2, False, True, 'Johnston', cunningham = True) # optional - cunningham
            [3, 8]

        ALGORITHMS:

        - The default for prime modulus is currently an algorithm described in the following paper:

        Johnston, Anna M. A generalized qth root algorithm. Proceedings of the tenth annual ACM-SIAM symposium on Discrete algorithms. Baltimore, 1999: pp 929-930.

        AUTHORS:

        - David Roe (2010-2-13)
        """
        if extend:
            raise NotImplementedError
        K = self.parent()
        n = Integer(n)
        if n == 0:
            if self == 1:
                if all: return [K(a) for a in range(1,K.order())]
                else: return self
            else:
                if all: return []
                else: raise ValueError
        F = K.factored_order()
        if len(F) == 0:
            if all:
                return [self]
            else:
                return self
        if len(F) != 1:
            if all:
                # we should probably do a first pass to see if there are any solutions so that we don't get giant intermediate lists and waste time...
                L = []
                for p, k in F:
                    L.append(mod(self, p**k).nth_root(n, all=True, algorithm=algorithm))
                ans = L[0]
                for i in range(1, len(L)):
                    ans = [a.crt(b) for a in ans for b in L[i]]
            else:
                ans = mod(0,1)
                for p, k in F:
                    ans = ans.crt(mod(self, p**k).nth_root(n, algorithm=algorithm))
            return ans
        p, k = F[0]
        if self.is_zero():
            if n < 0:
                if all: return []
                else: raise ValueError
            if all:
                if k == 1:
                    return [self]
                else:
                    minval = max(1, (k/n).ceil())
                    return [K(a*p**minval) for a in range(p**(k-minval))]
            else:
                return self
        if n < 0:
            try:
                self = ~self
            except ZeroDivisionError:
                if all: return []
                else: raise ValueError
            n = -n
        if p == 2 and k == 1:
            if all: return [self]
            else: return self
        if k > 1:
            pval, upart = self.lift().val_unit(p)
            if not n.divides(pval):
                if all:
                    return []
                else:
                    raise ValueError, "no nth root"
            if pval > 0:
                if all:
                    return [K(a.lift()*p**(pval // n) + p**(k - (pval - pval//n)) * b) for a in mod(upart, p**(k-pval)).nth_root(n, all=True, algorithm=algorithm) for b in range(p**(pval - pval//n))]
                else:
                    return K(p**(pval // n) * mod(upart, p**(k-pval)).nth_root(n, algorithm=algorithm).lift())
            from sage.rings.padics.all import ZpFM
            R = ZpFM(p,k)
            self_orig = self
            if p == 2:
                sign = [1]
                if self % 4 == 3:
                    if n % 2 == 0:
                        if all: return []
                        else: raise ValueError, "no nth root"
                    else:
                        sign = [-1]
                        self = -self
                elif n % 2 == 0:
                    if k > 2 and self % 8 == 5:
                        if all: return []
                        else: raise ValueError, "no nth root"
                    sign = [1, -1]
                if k == 2:
                    if all: return [K(s) for s in sign[:2]]
                    else: return K(sign[0])
                if all: modp = [mod(self,8)]
                else: modp = mod(self,8)
            else:
                sign = [1]
                modp = self % p
                self = self / K(R.teichmuller(modp))
                modp = modp.nth_root(n, all=all, algorithm=algorithm)
            # now self is congruent to 1 mod 4 or 1 mod p (for odd p), so the power series for p-adic log converges.
            # Hensel lifting is probably better, but this is easier at the moment.
            plog = R(self).log()
            nval = n.valuation(p)
            if nval >= plog.valuation() + (-1 if p == 2 else 0):
                if self == 1:
                    if all:
                        return [s*K(p*k+m.lift()) for k in range(p**(k-(2 if p==2 else 1))) for m in modp for s in sign]
                    else: return self_orig
                else:
                    if all: return []
                    else: raise ValueError, "no nth root"
            if all:
                ans = [plog // n + p**(k - nval) * i for i in range(p**nval)]
                ans = [s*K(R.teichmuller(m) * a.exp()) for a in ans for m in modp for s in sign]
                return ans
            else:
                return sign[0] * K(R.teichmuller(modp) * (plog // n).exp())
        return self._nth_root_common(n, all, algorithm, cunningham)

    def _nth_root_naive(self, n):
        """
        Computes all nth roots using brute force, for doc-testing.

        TESTS::

            sage: for n in range(2,100): # long time
            ...       K=Integers(n)
            ...       elist = range(1,min(2*n+2,100))
            ...       for e in random_sublist(elist, 5/len(elist)):
            ...           for a in random_sublist(range(1,n), min((n+2)//2,10)/(n-1)):
            ...               b = K(a)
            ...               try:
            ...                   L = b.nth_root(e, all=True)
            ...                   if len(L) > 0:
            ...                       c = b.nth_root(e)
            ...               except StandardError:
            ...                   L = [-1]
            ...               M = b._nth_root_naive(e)
            ...               if sorted(L) != M:
            ...                   print "mod(%s, %s).nth_root(%s,all=True), mod(%s, %s)._nth_root_naive(%s)"%(a,n,e,a,n,e)
            ...                   raise ValueError
            ...               if len(L) > 0 and (c not in L):
            ...                   print "mod(%s, %s).nth_root(%s), mod(%s, %s).nth_root(%s,all=True)"%(a,n,e,a,n,e)
            ...                   raise ValueError
        """
        L = []
        for a in self.parent():
            if a**n == self:
                L.append(a)
        return L

    def _balanced_abs(self):
        """
        This function returns `x` or `-x`, whichever has a
        positive representative in `-n/2 < x \leq n/2`.

        This is used so that the same square root is always returned,
        despite the possibly probabalistic nature of the underlying
        algorithm.
        """
        if self.lift() > self.__modulus.sageInteger >> 1:
            return -self
        else:
            return self


    def rational_reconstruction(self):
        """
        EXAMPLES::

            sage: R = IntegerModRing(97)
            sage: a = R(2) / R(3)
            sage: a
            33
            sage: a.rational_reconstruction()
            2/3
        """
        return self.lift().rational_reconstruction(self.modulus())

    def crt(IntegerMod_abstract self, IntegerMod_abstract other):
        r"""
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to
        ``self`` and to ``other``. The modulus of
        ``other`` must be coprime to the modulus of
        ``self``.

        EXAMPLES::

            sage: a = mod(3,5)
            sage: b = mod(2,7)
            sage: a.crt(b)
            23

        ::

            sage: a = mod(37,10^8)
            sage: b = mod(9,3^8)
            sage: a.crt(b)
            125900000037

        ::

            sage: b = mod(0,1)
            sage: a.crt(b) == a
            True
            sage: a.crt(b).modulus()
            100000000

        TESTS::

            sage: mod(0,1).crt(mod(4,2^127))
            4
            sage: mod(4,2^127).crt(mod(0,1))
            4
            sage: mod(4,2^30).crt(mod(0,1))
            4
            sage: mod(0,1).crt(mod(4,2^30))
            4
            sage: mod(0,1).crt(mod(4,2^15))
            4
            sage: mod(4,2^15).crt(mod(0,1))
            4

        AUTHORS:

        - Robert Bradshaw
        """
        cdef int_fast64_t new_modulus
        if not PY_TYPE_CHECK(self, IntegerMod_gmp) and not PY_TYPE_CHECK(other, IntegerMod_gmp):

            if other.__modulus.int64 == 1: return self
            new_modulus = self.__modulus.int64 * other.__modulus.int64
            if new_modulus < INTEGER_MOD_INT32_LIMIT:
                return self.__crt(other)

            elif new_modulus < INTEGER_MOD_INT64_LIMIT:
                if not PY_TYPE_CHECK(self, IntegerMod_int64):
                    self = IntegerMod_int64(self._parent, self.lift())
                if not PY_TYPE_CHECK(other, IntegerMod_int64):
                    other = IntegerMod_int64(other._parent, other.lift())
                return self.__crt(other)

        if not PY_TYPE_CHECK(self, IntegerMod_gmp):
            if self.__modulus.int64 == 1: return other
            self = IntegerMod_gmp(self._parent, self.lift())

        if not PY_TYPE_CHECK(other, IntegerMod_gmp):
            if other.__modulus.int64 == 1: return self
            other = IntegerMod_gmp(other._parent, other.lift())

        return self.__crt(other)


    def additive_order(self):
        r"""
        Returns the additive order of self.

        This is the same as ``self.order()``.

        EXAMPLES::

            sage: Integers(20)(2).additive_order()
            10
            sage: Integers(20)(7).additive_order()
            20
            sage: Integers(90308402384902)(2).additive_order()
            45154201192451
        """
        n = self.__modulus.sageInteger
        return sage.rings.integer.Integer(n.__floordiv__(self.lift().gcd(n)))

    def is_primitive_root(self):
        """
        Determines whether this element generates the group of units modulo n.

        This is only possible if the group of units is cyclic, which occurs if
        n is 2, 4, a power of an odd prime or twice a power of an odd prime.

        EXAMPLES::

            sage: mod(1,2).is_primitive_root()
            True
            sage: mod(3,4).is_primitive_root()
            True
            sage: mod(2,7).is_primitive_root()
            False
            sage: mod(3,98).is_primitive_root()
            True
            sage: mod(11,1009^2).is_primitive_root()
            True

        TESTS::

            sage: for p in prime_range(3,12):
            ...     for k in range(1,4):
            ...         for even in [1,2]:
            ...             n = even*p^k
            ...             phin = euler_phi(n)
            ...             for _ in range(6):
            ...                 a = Zmod(n).random_element()
            ...                 if not a.is_unit(): continue
            ...                 if a.is_primitive_root().__xor__(a.multiplicative_order()==phin):
            ...                     print "mod(%s,%s) incorrect"%(a,n)
        """
        cdef Integer p1, q = Integer(2)
        m = self.modulus()
        if m == 2:
            return self == 1
        if m == 4:
            return self == 3
        pow2, odd = m.val_unit(2)
        if pow2 > 1:
            return False
        if pow2 == 1:
            if self % 2 == 0:
                return False
            self = self % odd
        p, k = odd.perfect_power()
        if not p.is_prime():
            return False
        if k > 1:
            if self**((p-1)*p**(k-2)) == 1:
                return False
            # self**(p**(k-1)*(p-1)//q) = 1 for some q
            # iff mod(self,p)**((p-1)//q) = 1 for some q
            self = self % p
        # Now self is modulo a prime and need the factorization of p-1.
        p1 = p - 1
        while mpz_cmpabs_ui(p1.value, 1):
            q = p1.trial_division(bound=1000, start=mpz_get_ui(q.value))
            if q == p1:
                break
            if self**((p-1)//q) == 1:
                return False
            mpz_remove(p1.value, p1.value, q.value)
        if q.is_prime():
            return self**((p-1)//q) != 1
        # No small factors remain: we need to do some real work.
        for qq, e in q.factor():
            if self**((p-1)//qq) == 1:
                return False
        return True

    def multiplicative_order(self):
        """
        Returns the multiplicative order of self.

        EXAMPLES::

            sage: Mod(-1,5).multiplicative_order()
            2
            sage: Mod(1,5).multiplicative_order()
            1
            sage: Mod(0,5).multiplicative_order()
            Traceback (most recent call last):
            ...
            ArithmeticError: multiplicative order of 0 not defined since it is not a unit modulo 5
        """
        try:
            return sage.rings.integer.Integer(self._pari_().order())  # pari's "order" is by default multiplicative
        except PariError:
            raise ArithmeticError, "multiplicative order of %s not defined since it is not a unit modulo %s"%(
                self, self.__modulus.sageInteger)

    def valuation(self, p):
        """
        The largest power r such that m is in the ideal generated by p^r or infinity if there is not a largest such power.
        However it is an error to take the valuation with respect to a unit.

        .. NOTE::

            This is not a valuation in the mathematical sense. As shown with the examples below.

        EXAMPLES:

        This example shows that the (a*b).valuation(n) is not always the same as a.valuation(n) + b.valuation(n)

        ::

            sage: R=ZZ.quo(9)
            sage: a=R(3)
            sage: b=R(6)
            sage: a.valuation(3)
            1
            sage: a.valuation(3) + b.valuation(3)
            2
            sage: (a*b).valuation(3)
            +Infinity

        The valuation with respect to a unit is an error

        ::

            sage: a.valuation(4)
            Traceback (most recent call last):
            ...
            ValueError: Valuation with respect to a unit is not defined.

        TESTS::

            sage: R=ZZ.quo(12)
            sage: a=R(2)
            sage: b=R(4)
            sage: a.valuation(2)
            1
            sage: b.valuation(2)
            +Infinity
            sage: ZZ.quo(1024)(16).valuation(4)
            2

        """
        p=self.__modulus.sageInteger.gcd(p)
        if p==1:
            raise ValueError("Valuation with respect to a unit is not defined.")
        r = 0
        power = p
        while not (self % power): # self % power == 0
            r += 1
            power *= p
            if not power.divides(self.__modulus.sageInteger):
                from sage.rings.all import infinity
                return infinity
        return r

    def __floordiv__(self, other):
        """
        Exact division for prime moduli, for compatibility with other fields.

        EXAMPLES:
        sage: GF(7)(3) // GF(7)(5)
        2
        """
        # needs to be rewritten for coercion
        if other.parent() is not self.parent():
            other = self.parent().coerce(other)
        if self.parent().is_field():
            return self / other
        else:
            raise TypeError, "Floor division not defined for non-prime modulus"

    def _repr_(self):
        return str(self.lift())

    def _latex_(self):
        return str(self)

    def _integer_(self, ZZ=None):
        return self.lift()

    def _rational_(self):
        return rational.Rational(self.lift())




######################################################################
#      class IntegerMod_gmp
######################################################################


cdef class IntegerMod_gmp(IntegerMod_abstract):
    """
    Elements of `\ZZ/n\ZZ` for n not small enough
    to be operated on in word size.

    AUTHORS:

    - Robert Bradshaw (2006-08-24)
    """

    def __init__(IntegerMod_gmp self, parent, value, empty=False):
        """
        EXAMPLES::

            sage: a = mod(5,14^20)
            sage: type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: loads(dumps(a)) == a
            True
        """
        mpz_init(self.value)
        IntegerMod_abstract.__init__(self, parent)
        if empty:
            return
        cdef sage.rings.integer.Integer z
        if PY_TYPE_CHECK(value, sage.rings.integer.Integer):
            z = value
        elif PY_TYPE_CHECK(value, rational.Rational):
            z = value % self.__modulus.sageInteger
        elif PY_TYPE_CHECK(value, int):
            self.set_from_long(value)
            return
        else:
            z = sage.rings.integer_ring.Z(value)
        self.set_from_mpz(z.value)

    cdef IntegerMod_gmp _new_c(self):
        cdef IntegerMod_gmp x
        x = PY_NEW(IntegerMod_gmp)
        mpz_init(x.value)
        x.__modulus = self.__modulus
        x._parent = self._parent
        return x

    def __dealloc__(self):
        mpz_clear(self.value)

    cdef void set_from_mpz(self, mpz_t value):
        cdef sage.rings.integer.Integer modulus
        modulus = self.__modulus.sageInteger
        if mpz_sgn(value) == -1 or mpz_cmp(value, modulus.value) >= 0:
            mpz_mod(self.value, value, modulus.value)
        else:
            mpz_set(self.value, value)

    cdef void set_from_long(self, long value):
        cdef sage.rings.integer.Integer modulus
        mpz_set_si(self.value, value)
        if value < 0 or mpz_cmp_si(self.__modulus.sageInteger.value, value) >= 0:
            mpz_mod(self.value, self.value, self.__modulus.sageInteger.value)

    cdef mpz_t* get_value(IntegerMod_gmp self):
        return &self.value

    def __lshift__(IntegerMod_gmp self, k):
        r"""
        Performs a left shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(19, 10^10)
            sage: e << 102
            9443608576
        """
        return self.shift(long(k))

    def __rshift__(IntegerMod_gmp self, k):
        r"""
        Performs a right shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(19, 10^10)
            sage: e >> 1
            9
        """
        return self.shift(-long(k))

    cdef shift(IntegerMod_gmp self, long k):
        r"""
        Performs a bit-shift specified by ``k`` on ``self``.

        Suppose that ``self`` represents an integer `x` modulo `n`.  If `k` is
        `k = 0`, returns `x`.  If `k > 0`, shifts `x` to the left, that is,
        multiplies `x` by `2^k` and then returns the representative in the
        range `[0,n)`.  If `k < 0`, shifts `x` to the right, that is, returns
        the integral part of `x` divided by `2^k`.

        Note that, in any case, ``self`` remains unchanged.

        INPUT:

        - ``k`` - Integer of type ``long``

        OUTPUT

        - Result of type ``IntegerMod_gmp``

        EXAMPLES::

            sage: e = Mod(19, 10^10)
            sage: e << 102
            9443608576
            sage: e >> 1
            9
            sage: e >> 4
            1
        """
        cdef IntegerMod_gmp x
        if k == 0:
            return self
        else:
            x = self._new_c()
            if k > 0:
                mpz_mul_2exp(x.value, self.value, k)
                mpz_fdiv_r(x.value, x.value, self.__modulus.sageInteger.value)
            else:
                mpz_fdiv_q_2exp(x.value, self.value, -k)
            return x

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLES::

            sage: mod(5,13^20) == mod(5,13^20)
            True
            sage: mod(5,13^20) == mod(-5,13^20)
            False
            sage: mod(5,13^20) == mod(-5,13)
            False
        """
        cdef int i
        i = mpz_cmp((<IntegerMod_gmp>left).value, (<IntegerMod_gmp>right).value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)


    cpdef bint is_one(IntegerMod_gmp self):
        """
        Returns ``True`` if this is `1`, otherwise
        ``False``.

        EXAMPLES::

            sage: mod(1,5^23).is_one()
            True
            sage: mod(0,5^23).is_one()
            False
        """
        return mpz_cmp_si(self.value, 1) == 0

    def __nonzero__(IntegerMod_gmp self):
        """
        Returns ``True`` if this is not `0`, otherwise
        ``False``.

        EXAMPLES::

            sage: mod(13,5^23).is_zero()
            False
            sage: (mod(25,5^23)^23).is_zero()
            True
        """
        return mpz_cmp_si(self.value, 0) != 0

    cpdef bint is_unit(self):
        """
        Return True iff this element is a unit.

        EXAMPLES::

            sage: mod(13, 5^23).is_unit()
            True
            sage: mod(25, 5^23).is_unit()
            False
        """
        return self.lift().gcd(self.modulus()) == 1

    def __crt(IntegerMod_gmp self, IntegerMod_gmp other):
        cdef IntegerMod_gmp lift, x
        cdef sage.rings.integer.Integer modulus, other_modulus

        modulus = self.__modulus.sageInteger
        other_modulus = other.__modulus.sageInteger
        import integer_mod_ring
        lift = IntegerMod_gmp(integer_mod_ring.IntegerModRing(modulus*other_modulus), None, empty=True)
        try:
            if mpz_cmp(self.value, other.value) > 0:
                x = (other - IntegerMod_gmp(other._parent, self.lift())) / IntegerMod_gmp(other._parent, modulus)
                mpz_mul(lift.value, x.value, modulus.value)
                mpz_add(lift.value, lift.value, self.value)
            else:
                x = (self - IntegerMod_gmp(self._parent, other.lift())) / IntegerMod_gmp(self._parent, other_modulus)
                mpz_mul(lift.value, x.value, other_modulus.value)
                mpz_add(lift.value, lift.value, other.value)
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"


    def __copy__(IntegerMod_gmp self):
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_set(x.value, self.value)
        return x

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R(7) + R(8)
            15
        """
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_add(x.value, self.value, (<IntegerMod_gmp>right).value)
        if mpz_cmp(x.value, self.__modulus.sageInteger.value)  >= 0:
            mpz_sub(x.value, x.value, self.__modulus.sageInteger.value)
        return x;

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R(7) + R(8)
            15
        """
        mpz_add(self.value, self.value, (<IntegerMod_gmp>right).value)
        if mpz_cmp(self.value, self.__modulus.sageInteger.value)  >= 0:
            mpz_sub(self.value, self.value, self.__modulus.sageInteger.value)
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R(7) - R(8)
            9999999999
        """
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_sub(x.value, self.value, (<IntegerMod_gmp>right).value)
        if mpz_sgn(x.value) == -1:
            mpz_add(x.value, x.value, self.__modulus.sageInteger.value)
        return x;

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R(7) - R(8)
            9999999999
        """
        mpz_sub(self.value, self.value, (<IntegerMod_gmp>right).value)
        if mpz_sgn(self.value) == -1:
            mpz_add(self.value, self.value, self.__modulus.sageInteger.value)
        return self

    cpdef ModuleElement _neg_(self):
        """
        EXAMPLES::

            sage: -mod(5,10^10)
            9999999995
            sage: -mod(0,10^10)
            0
        """
        if mpz_cmp_si(self.value, 0) == 0:
            return self
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_sub(x.value, self.__modulus.sageInteger.value, self.value)
        return x

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^11)
            sage: R(700000) * R(800000)
            60000000000
        """
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_mul(x.value, self.value,  (<IntegerMod_gmp>right).value)
        mpz_fdiv_r(x.value, x.value, self.__modulus.sageInteger.value)
        return x

    cpdef RingElement _imul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^11)
            sage: R(700000) * R(800000)
            60000000000
        """
        mpz_mul(self.value, self.value,  (<IntegerMod_gmp>right).value)
        mpz_fdiv_r(self.value, self.value, self.__modulus.sageInteger.value)
        return self

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^11)
            sage: R(3) / R(7)
            71428571429
        """
        return self._mul_(~right)

    def __int__(self):
        return int(self.lift())

    def __index__(self):
        """
        Needed so integers modulo `n` can be used as list indices.

        EXAMPLES::

            sage: v = [1,2,3,4,5]
            sage: v[Mod(3,10^20)]
            4
        """
        return int(self.lift())

    def __long__(self):
        return long(self.lift())

    def __mod__(self, right):
        if self.modulus() % right != 0:
            raise ZeroDivisionError, "reduction modulo right not defined."
        import integer_mod_ring
        return IntegerMod(integer_mod_ring.IntegerModRing(right), self)

    def __pow__(IntegerMod_gmp self, exp, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES:
            sage: R = Integers(10^10)
            sage: R(2)^1000
            5668069376
            sage: p = next_prime(11^10)
            sage: R = Integers(p)
            sage: R(9876)^(p-1)
            1
            sage: mod(3, 10^100)^-2
            8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889
            sage: mod(2, 10^100)^-2
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        TESTS:

        We define ``0^0`` to be unity, :trac:`13894`::

            sage: p = next_prime(11^10)
            sage: R = Integers(p)
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: type(R(0)^0) == type(R(0))
            True

        When the modulus is ``1``, the only element in the ring is
        ``0`` (and it is equivalent to ``1``), so we return that
        instead::

            sage: from sage.rings.finite_rings.integer_mod \
            ...       import IntegerMod_gmp
            sage: zero = IntegerMod_gmp(Integers(1),0)
            sage: type(zero)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: zero^0
            0

        """
        cdef IntegerMod_gmp x = self._new_c()
        sig_on()
        try:
            mpz_pow_helper(x.value, self.value, exp, self.__modulus.sageInteger.value)
            return x
        finally:
            sig_off()

    def __invert__(IntegerMod_gmp self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES::

            sage: a = mod(3,10^100); type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: ~a
            6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667
            sage: ~mod(2,10^100)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
        """
        if self.is_zero():
            raise ZeroDivisionError, "Inverse does not exist."

        cdef IntegerMod_gmp x
        x = self._new_c()
        if mpz_invert(x.value, self.value, self.__modulus.sageInteger.value):
            return x
        else:
            raise ZeroDivisionError, "Inverse does not exist."

    def lift(IntegerMod_gmp self):
        """
        Lift an integer modulo `n` to the integers.

        EXAMPLES::

            sage: a = Mod(8943, 2^70); type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: lift(a)
            8943
            sage: a.lift()
            8943
        """
        cdef sage.rings.integer.Integer z
        z = sage.rings.integer.Integer()
        z.set_from_mpz(self.value)
        return z

    def __float__(self):
        return float(self.lift())

    def __hash__(self):
        """
        EXAMPLES::

            sage: a = Mod(8943, 2^100)
            sage: hash(a)
            8943
        """
#        return mpz_pythonhash(self.value)
        return hash(self.lift())

    @coerce_binop
    def gcd(self, IntegerMod_gmp other):
        """
        Greatest common divisor

        Returns the "smallest" generator in `\ZZ / N\ZZ` of the ideal
        generated by ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of the same ring as this one.

        EXAMPLES::

            sage: mod(2^3*3^2*5, 3^3*2^2*17^8).gcd(mod(2^4*3*17, 3^3*2^2*17^8))
            12
            sage: mod(0,17^8).gcd(mod(0,17^8))
            0
        """
        cdef IntegerMod_gmp ans = self._new_c()
        sig_on()
        mpz_gcd(ans.value, self.value, self.__modulus.sageInteger.value)
        mpz_gcd(ans.value, ans.value, other.value)
        sig_off()
        if mpz_cmp(ans.value, self.__modulus.sageInteger.value) == 0:
            # self = other = 0
            mpz_set_ui(ans.value, 0)
        return ans

######################################################################
#      class IntegerMod_int
######################################################################


cdef class IntegerMod_int(IntegerMod_abstract):
    """
    Elements of `\ZZ/n\ZZ` for n small enough to
    be operated on in 32 bits

    AUTHORS:

    - Robert Bradshaw (2006-08-24)
    """

    def __init__(self, parent, value, empty=False):
        """
        EXAMPLES::

            sage: a = Mod(10,30); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        IntegerMod_abstract.__init__(self, parent)
        if empty:
            return
        if self.__modulus.int32 == 1:
            self.ivalue = 0
            return
        cdef long x
        if PY_TYPE_CHECK(value, int):
            x = value
            self.ivalue = x % self.__modulus.int32
            if self.ivalue < 0:
                self.ivalue = self.ivalue + self.__modulus.int32
            return
        elif PY_TYPE_CHECK(value, IntegerMod_int):
            self.ivalue = (<IntegerMod_int>value).ivalue % self.__modulus.int32
            return
        cdef sage.rings.integer.Integer z
        if PY_TYPE_CHECK(value, sage.rings.integer.Integer):
            z = value
        elif PY_TYPE_CHECK(value, rational.Rational):
            z = value % self.__modulus.sageInteger
        else:
            z = sage.rings.integer_ring.Z(value)
        self.set_from_mpz(z.value)

    def _make_new_with_parent_c(self, parent): #ParentWithBase parent):
        cdef IntegerMod_int x = PY_NEW(IntegerMod_int)
        x._parent = parent
        x.__modulus = parent._pyx_order
        x.ivalue = self.ivalue
        return x

    cdef IntegerMod_int _new_c(self, int_fast32_t value):
        if self.__modulus.table is not None:
            return self.__modulus.lookup(value)
        cdef IntegerMod_int x = PY_NEW(IntegerMod_int)
        x._parent = self._parent
        x.__modulus = self.__modulus
        x.ivalue = value
        return x

    cdef void set_from_mpz(self, mpz_t value):
        if mpz_sgn(value) == -1 or mpz_cmp_si(value, self.__modulus.int32) >= 0:
            self.ivalue = mpz_fdiv_ui(value, self.__modulus.int32)
        else:
            self.ivalue = mpz_get_si(value)

    cdef void set_from_long(self, long value):
        self.ivalue = value % self.__modulus.int32

    cdef void set_from_int(IntegerMod_int self, int_fast32_t ivalue):
        if ivalue < 0:
            self.ivalue = self.__modulus.int32 + (ivalue % self.__modulus.int32)
        elif ivalue >= self.__modulus.int32:
            self.ivalue = ivalue % self.__modulus.int32
        else:
            self.ivalue = ivalue

    cdef int_fast32_t get_int_value(IntegerMod_int self):
        return self.ivalue



    cdef int _cmp_c_impl(self, Element right) except -2:
        """
        EXAMPLES::

            sage: mod(5,13) == mod(-8,13)
            True
            sage: mod(5,13) == mod(8,13)
            False
            sage: mod(5,13) == mod(5,24)
            False
            sage: mod(0, 13) == 0
            True
            sage: mod(0, 13) == int(0)
            True
        """
        if self.ivalue == (<IntegerMod_int>right).ivalue:
            return 0
        elif self.ivalue < (<IntegerMod_int>right).ivalue:
            return -1
        else:
            return 1

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)


    cpdef bint is_one(IntegerMod_int self):
        """
        Returns ``True`` if this is `1`, otherwise
        ``False``.

        EXAMPLES::

            sage: mod(6,5).is_one()
            True
            sage: mod(0,5).is_one()
            False
        """
        return self.ivalue == 1

    def __nonzero__(IntegerMod_int self):
        """
        Returns ``True`` if this is not `0`, otherwise
        ``False``.

        EXAMPLES::

            sage: mod(13,5).is_zero()
            False
            sage: mod(25,5).is_zero()
            True
        """
        return self.ivalue != 0

    cpdef bint is_unit(IntegerMod_int self):
        """
        Return True iff this element is a unit

        EXAMPLES::

            sage: a=Mod(23,100)
            sage: a.is_unit()
            True
            sage: a=Mod(24,100)
            sage: a.is_unit()
            False
        """
        return gcd_int(self.ivalue, self.__modulus.int32) == 1

    def __crt(IntegerMod_int self, IntegerMod_int other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to self and
        to other. The modulus of other must be coprime to the modulus of
        self.

        EXAMPLES::

            sage: a = mod(3,5)
            sage: b = mod(2,7)
            sage: a.crt(b)
            23

        AUTHORS:

        - Robert Bradshaw
        """
        cdef IntegerMod_int lift
        cdef int_fast32_t x

        import integer_mod_ring
        lift = IntegerMod_int(integer_mod_ring.IntegerModRing(self.__modulus.int32 * other.__modulus.int32), None, empty=True)

        try:
            x = (other.ivalue - self.ivalue % other.__modulus.int32) * mod_inverse_int(self.__modulus.int32, other.__modulus.int32)
            lift.set_from_int( x * self.__modulus.int32 + self.ivalue )
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"


    def __copy__(IntegerMod_int self):
        cdef IntegerMod_int x = PY_NEW(IntegerMod_int)
        x._parent = self._parent
        x.__modulus = self.__modulus
        x.ivalue = self.ivalue
        return x

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) + R(8)
            5
        """
        cdef int_fast32_t x
        x = self.ivalue + (<IntegerMod_int>right).ivalue
        if x >= self.__modulus.int32:
            x = x - self.__modulus.int32
        return self._new_c(x)

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) + R(8)
            5
        """
        cdef int_fast32_t x
        x = self.ivalue + (<IntegerMod_int>right).ivalue
        if x >= self.__modulus.int32:
            x = x - self.__modulus.int32
        self.ivalue = x
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef int_fast32_t x
        x = self.ivalue - (<IntegerMod_int>right).ivalue
        if x < 0:
            x = x + self.__modulus.int32
        return self._new_c(x)

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef int_fast32_t x
        x = self.ivalue - (<IntegerMod_int>right).ivalue
        if x < 0:
            x = x + self.__modulus.int32
        self.ivalue = x
        return self

    cpdef ModuleElement _neg_(self):
        """
        EXAMPLES::

            sage: -mod(7,10)
            3
            sage: -mod(0,10)
            0
        """
        if self.ivalue == 0:
            return self
        return self._new_c(self.__modulus.int32 - self.ivalue)

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        return self._new_c((self.ivalue * (<IntegerMod_int>right).ivalue) % self.__modulus.int32)

    cpdef RingElement _imul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        self.ivalue = (self.ivalue * (<IntegerMod_int>right).ivalue) % self.__modulus.int32
        return self

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(2)/3
            4
        """
        if self.__modulus.inverses is not None:
            right_inverse = self.__modulus.inverses[(<IntegerMod_int>right).ivalue]
            if right_inverse is None:
                raise ZeroDivisionError, "Inverse does not exist."
            else:
                return self._new_c((self.ivalue * (<IntegerMod_int>right_inverse).ivalue) % self.__modulus.int32)

        cdef int_fast32_t x
        x = self.ivalue * mod_inverse_int((<IntegerMod_int>right).ivalue, self.__modulus.int32)
        return self._new_c(x% self.__modulus.int32)

    def __int__(IntegerMod_int self):
        return self.ivalue

    def __index__(self):
        """
        Needed so integers modulo `n` can be used as list indices.

        EXAMPLES::

            sage: v = [1,2,3,4,5]
            sage: v[Mod(10,7)]
            4
        """
        return self.ivalue

    def __long__(IntegerMod_int self):
        return self.ivalue

    def __mod__(IntegerMod_int self, right):
        right = int(right)
        if self.__modulus.int32 % right != 0:
            raise ZeroDivisionError, "reduction modulo right not defined."
        import integer_mod_ring
        return integer_mod_ring.IntegerModRing(right)(self)

    def __lshift__(IntegerMod_int self, k):
        r"""
        Performs a left shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(5, 2^10 - 1)
            sage: e << 5
            160
            sage: e * 2^5
            160
        """
        return self.shift(int(k))

    def __rshift__(IntegerMod_int self, k):
        r"""
        Performs a right shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(5, 2^10 - 1)
            sage: e << 5
            160
            sage: e * 2^5
            160
        """
        return self.shift(-int(k))

    cdef shift(IntegerMod_int self, int k):
        """
        Performs a bit-shift specified by ``k`` on ``self``.

        Suppose that ``self`` represents an integer `x` modulo `n`.  If `k` is
        `k = 0`, returns `x`.  If `k > 0`, shifts `x` to the left, that is,
        multiplies `x` by `2^k` and then returns the representative in the
        range `[0,n)`.  If `k < 0`, shifts `x` to the right, that is, returns
        the integral part of `x` divided by `2^k`.

        Note that, in any case, ``self`` remains unchanged.

        INPUT:

        - ``k`` - Integer of type ``int``

        OUTPUT:

        - Result of type ``IntegerMod_int``

        WARNING:

        For positive ``k``, if ``x << k`` overflows as a 32-bit integer, the
        result is meaningless.

        EXAMPLES::

            sage: e = Mod(5, 2^10 - 1)
            sage: e << 5
            160
            sage: e * 2^5
            160
            sage: e = Mod(8, 2^5 - 1)
            sage: e >> 3
            1
            sage: int(e)/int(2^3)
            1
        """
        if k == 0:
            return self
        elif k > 0:
            return self._new_c((self.ivalue << k) % self.__modulus.int32)
        else:
            return self._new_c(self.ivalue >> (-k))

    def __pow__(IntegerMod_int self, exp, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(2)^10
            4
            sage: R = Integers(389)
            sage: R(7)^388
            1

            sage: mod(3, 100)^-1
            67
            sage: mod(3, 100)^-100000000
            1

            sage: mod(2, 100)^-1
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
            sage: mod(2, 100)^-100000000
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        TESTS:

        We define ``0^0`` to be unity, :trac:`13894`::

            sage: R = Integers(100)
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: type(R(0)^0) == type(R(0))
            True

        When the modulus is ``1``, the only element in the ring is
        ``0`` (and it is equivalent to ``1``), so we return that
        instead::

            sage: R = Integers(1)
            sage: R(0)^0
            0

        """
        cdef long long_exp
        cdef int_fast32_t res
        cdef mpz_t res_mpz
        if PyInt_CheckExact(exp) and -100000 < PyInt_AS_LONG(exp) < 100000:
            long_exp = PyInt_AS_LONG(exp)
        elif PY_TYPE_CHECK_EXACT(exp, Integer) and mpz_cmpabs_ui((<Integer>exp).value, 100000) == -1:
            long_exp = mpz_get_si((<Integer>exp).value)
        else:
            sig_on()
            try:
                mpz_init(res_mpz)
                base = self.lift()
                mpz_pow_helper(res_mpz, (<Integer>base).value, exp, self.__modulus.sageInteger.value)
                return self._new_c(mpz_get_ui(res_mpz))
            finally:
                mpz_clear(res_mpz)
                sig_off()

        if long_exp == 0 and self.ivalue == 0:
            # Return 0 if the modulus is 1, otherwise return 1.
            return self._new_c(self.__modulus.int32 != 1)
        cdef bint invert = False
        if long_exp < 0:
            invert = True
            long_exp = -long_exp
        res = mod_pow_int(self.ivalue, long_exp, self.__modulus.int32)
        if invert:
            return ~self._new_c(res)
        else:
            return self._new_c(res)

    def __invert__(IntegerMod_int self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES::

            sage: ~mod(7,100)
            43
            sage: Mod(0,1)^-1
            0
        """
        if self.__modulus.inverses is not None:
            x = self.__modulus.inverses[self.ivalue]
            if x is None:
                raise ZeroDivisionError, "Inverse does not exist."
            else:
                return x
        else:
            return self._new_c(mod_inverse_int(self.ivalue, self.__modulus.int32))

    def lift(IntegerMod_int self):
        """
        Lift an integer modulo `n` to the integers.

        EXAMPLES::

            sage: a = Mod(8943, 2^10); type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: lift(a)
            751
            sage: a.lift()
            751
        """
        cdef sage.rings.integer.Integer z
        z = sage.rings.integer.Integer()
        mpz_set_si(z.value, self.ivalue)
        return z

    def __float__(IntegerMod_int self):
        return <double>self.ivalue

    def __hash__(self):
        """
        EXAMPLES::

            sage: a = Mod(89, 2^10)
            sage: hash(a)
            89
        """
        return hash(self.ivalue)

    cdef bint is_square_c(self) except -2:
        if self.ivalue <= 1:
            return 1
        moduli = self._parent.factored_order()
        cdef int val, e
        cdef int_fast32_t p
        if len(moduli) == 1:
            sage_p, e = moduli[0]
            p = sage_p
            if e == 1:
                return jacobi_int(self.ivalue, p) != -1
            elif p == 2:
                return self._pari_().issquare() # TODO: implement directly
            elif self.ivalue % p == 0:
                val = self.lift().valuation(sage_p)
                return val >= e or (val % 2 == 0 and jacobi_int(self.ivalue / int(sage_p**val), p) != -1)
            else:
                return jacobi_int(self.ivalue, p) != -1
        else:
            for sage_p, e in moduli:
                p = sage_p
                if p == 2:
                    if e > 1 and not self._pari_().issquare(): # TODO: implement directly
                        return 0
                elif e > 1 and self.ivalue % p == 0:
                    val = self.lift().valuation(sage_p)
                    if val < e and (val % 2 == 1 or jacobi_int(self.ivalue / int(sage_p**val), p) == -1):
                        return 0
                elif jacobi_int(self.ivalue, p) == -1:
                    return 0
            return 1

    def sqrt(self, extend=True, all=False):
        r"""
        Returns square root or square roots of ``self`` modulo
        `n`.

        INPUT:


        -  ``extend`` - bool (default: ``True``);
           if ``True``, return a square root in an extension ring,
           if necessary. Otherwise, raise a ``ValueError`` if the
           square root is not in the base ring.

        -  ``all`` - bool (default: ``False``); if
           ``True``, return {all} square roots of self, instead of
           just one.


        ALGORITHM: Calculates the square roots mod `p` for each of
        the primes `p` dividing the order of the ring, then lifts
        them `p`-adically and uses the CRT to find a square root
        mod `n`.

        See also ``square_root_mod_prime_power`` and
        ``square_root_mod_prime`` (in this module) for more
        algorithmic details.

        EXAMPLES::

            sage: mod(-1, 17).sqrt()
            4
            sage: mod(5, 389).sqrt()
            86
            sage: mod(7, 18).sqrt()
            5
            sage: a = mod(14, 5^60).sqrt()
            sage: a*a
            14
            sage: mod(15, 389).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square
            sage: Mod(1/9, next_prime(2^40)).sqrt()^(-2)
            9
            sage: Mod(1/25, next_prime(2^90)).sqrt()^(-2)
            25

        ::

            sage: a = Mod(3,5); a
            3
            sage: x = Mod(-1, 360)
            sage: x.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square
            sage: y = x.sqrt(); y
            sqrt359
            sage: y.parent()
            Univariate Quotient Polynomial Ring in sqrt359 over Ring of integers modulo 360 with modulus x^2 + 1
            sage: y^2
            359

        We compute all square roots in several cases::

            sage: R = Integers(5*2^3*3^2); R
            Ring of integers modulo 360
            sage: R(40).sqrt(all=True)
            [20, 160, 200, 340]
            sage: [x for x in R if x^2 == 40]  # Brute force verification
            [20, 160, 200, 340]
            sage: R(1).sqrt(all=True)
            [1, 19, 71, 89, 91, 109, 161, 179, 181, 199, 251, 269, 271, 289, 341, 359]
            sage: R(0).sqrt(all=True)
            [0, 60, 120, 180, 240, 300]
            sage: GF(107)(0).sqrt(all=True)
            [0]

        ::

            sage: R = Integers(5*13^3*37); R
            Ring of integers modulo 406445
            sage: v = R(-1).sqrt(all=True); v
            [78853, 111808, 160142, 193097, 213348, 246303, 294637, 327592]
            sage: [x^2 for x in v]
            [406444, 406444, 406444, 406444, 406444, 406444, 406444, 406444]
            sage: v = R(169).sqrt(all=True); min(v), -max(v), len(v)
            (13, 13, 104)
            sage: all([x^2==169 for x in v])
            True

        Modulo a power of 2::

            sage: R = Integers(2^7); R
            Ring of integers modulo 128
            sage: a = R(17)
            sage: a.sqrt()
            23
            sage: a.sqrt(all=True)
            [23, 41, 87, 105]
            sage: [x for x in R if x^2==17]
            [23, 41, 87, 105]
        """
        cdef int_fast32_t i, n = self.__modulus.int32
        if n > 100:
            moduli = self._parent.factored_order()
        # Unless the modulus is tiny, test to see if we're in the really
        # easy case of n prime, n = 3 mod 4.
        if n > 100 and n % 4 == 3 and len(moduli) == 1 and moduli[0][1] == 1:
            if jacobi_int(self.ivalue, self.__modulus.int32) == 1:
                # it's a non-zero square, sqrt(a) = a^(p+1)/4
                i = mod_pow_int(self.ivalue, (self.__modulus.int32+1)/4, n)
                if i > n/2:
                    i = n-i
                if all:
                    return [self._new_c(i), self._new_c(n-i)]
                else:
                    return self._new_c(i)
            elif self.ivalue == 0:
                return [self] if all else self
            elif not extend:
                raise ValueError, "self must be a square"
        # Now we use a heuristic to guess whether or not it will
        # be faster to just brute-force search for squares in a c loop...
        # TODO: more tuning?
        elif n <= 100 or n / (1 << len(moduli)) < 5000:
            if all:
                return [self._new_c(i) for i from 0 <= i < n if (i*i) % n == self.ivalue]
            else:
                for i from 0 <= i <= n/2:
                    if (i*i) % n == self.ivalue:
                        return self._new_c(i)
                if not extend:
                    if all:
                        return []
                    raise ValueError, "self must be a square"
        # Either it failed but extend was True, or the generic algorithm is better
        return IntegerMod_abstract.sqrt(self, extend=extend, all=all)


    def _balanced_abs(self):
        """
        This function returns `x` or `-x`, whichever has a
        positive representative in `-n/2 < x \leq n/2`.
        """
        if self.ivalue > self.__modulus.int32 / 2:
            return -self
        else:
            return self

    @coerce_binop
    def gcd(self, IntegerMod_int other):
        """
        Greatest common divisor

        Returns the "smallest" generator in `\ZZ / N\ZZ` of the ideal
        generated by ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of the same ring as this one.

        EXAMPLES::

            sage: R = Zmod(60); S = Zmod(72)
            sage: a = R(40).gcd(S(30)); a
            2
            sage: a.parent()
            Ring of integers modulo 12
            sage: b = R(17).gcd(60); b
            1
            sage: b.parent()
            Ring of integers modulo 60

            sage: mod(72*5, 3^3*2^2*17^2).gcd(mod(48*17, 3^3*2^2*17^2))
            12
            sage: mod(0,1).gcd(mod(0,1))
            0
        """
        cdef int_fast32_t g = gcd_int(self.ivalue, self.__modulus.int32)
        g = gcd_int(g, other.ivalue)
        if g == self.__modulus.int32: # self = other = 0
            g = 0
        return self._new_c(g)

### End of class


cdef int_fast32_t gcd_int(int_fast32_t a, int_fast32_t b):
    """
    Returns the gcd of a and b

    For use with IntegerMod_int

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int_fast32_t tmp
    if a < b:
        tmp = b
        b = a
        a = tmp
    while b:
        tmp = b
        b = a % b
        a = tmp
    return a


cdef int_fast32_t mod_inverse_int(int_fast32_t x, int_fast32_t n) except 0:
    """
    Returns y such that xy=1 mod n

    For use in IntegerMod_int

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int_fast32_t tmp, a, b, last_t, t, next_t, q
    if n == 1:
        return 0
    a = n
    b = x
    t = 0
    next_t = 1
    while b:
        # a = s * n + t * x
        if b == 1:
            next_t = next_t % n
            if next_t < 0:
                next_t = next_t + n
            return next_t
        q = a / b
        tmp = b
        b = a % b
        a = tmp
        last_t = t
        t = next_t
        next_t = last_t - q * t
    raise ZeroDivisionError, "Inverse does not exist."


cdef int_fast32_t mod_pow_int(int_fast32_t base, int_fast32_t exp, int_fast32_t n):
    """
    Returns base^exp mod n

    For use in IntegerMod_int

    EXAMPLES::

        sage: z = Mod(2, 256)
        sage: z^8
        0

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int_fast32_t prod, pow2
    if exp <= 5:
        if exp == 0: return 1
        if exp == 1: return base
        prod = base * base % n
        if exp == 2: return prod
        if exp == 3: return (prod * base) % n
        if exp == 4: return (prod * prod) % n

    pow2 = base
    if exp % 2: prod = base
    else: prod = 1
    exp = exp >> 1
    while(exp != 0):
        pow2 = pow2 * pow2
        if pow2 >= INTEGER_MOD_INT32_LIMIT: pow2 = pow2 % n
        if exp % 2:
            prod = prod * pow2
            if prod >= INTEGER_MOD_INT32_LIMIT: prod = prod % n
        exp = exp >> 1

    if prod >= n:
        prod = prod % n
    return prod


cdef int jacobi_int(int_fast32_t a, int_fast32_t m) except -2:
    """
    Calculates the jacobi symbol (a/n)

    For use in IntegerMod_int

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int s, jacobi = 1
    cdef int_fast32_t b

    a = a % m

    while 1:
        if a == 0:
            return 0 # gcd was nontrivial
        elif a == 1:
            return jacobi
        s = 0
        while (1 << s) & a == 0:
            s += 1
        b = a >> s
        # Now a = 2^s * b

        # factor out (2/m)^s term
        if s % 2 == 1 and (m % 8 == 3 or m % 8 == 5):
            jacobi = -jacobi

        if b == 1:
            return jacobi

        # quadratic reciprocity
        if b % 4 == 3 and m % 4 == 3:
            jacobi = -jacobi
        a = m % b
        m = b

######################################################################
#      class IntegerMod_int64
######################################################################

cdef class IntegerMod_int64(IntegerMod_abstract):
    """
    Elements of `\ZZ/n\ZZ` for n small enough to
    be operated on in 64 bits

    AUTHORS:

    - Robert Bradshaw (2006-09-14)
    """

    def __init__(self, parent, value, empty=False):
        """
        EXAMPLES::

            sage: a = Mod(10,3^10); a
            10
            sage: type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
            sage: loads(a.dumps()) == a
            True
            sage: Mod(5, 2^31)
            5
        """
        IntegerMod_abstract.__init__(self, parent)
        if empty:
            return
        cdef int_fast64_t x
        if PY_TYPE_CHECK(value, int):
            x = value
            self.ivalue = x % self.__modulus.int64
            if self.ivalue < 0:
                self.ivalue = self.ivalue + self.__modulus.int64
            return
        cdef sage.rings.integer.Integer z
        if PY_TYPE_CHECK(value, sage.rings.integer.Integer):
            z = value
        elif PY_TYPE_CHECK(value, rational.Rational):
            z = value % self.__modulus.sageInteger
        else:
            z = sage.rings.integer_ring.Z(value)
        self.set_from_mpz(z.value)

    cdef IntegerMod_int64 _new_c(self, int_fast64_t value):
        cdef IntegerMod_int64 x
        x = PY_NEW(IntegerMod_int64)
        x.__modulus = self.__modulus
        x._parent = self._parent
        x.ivalue = value
        return x

    cdef void set_from_mpz(self, mpz_t value):
        if mpz_sgn(value) == -1 or mpz_cmp_si(value, self.__modulus.int64) >= 0:
            self.ivalue = mpz_fdiv_ui(value, self.__modulus.int64)
        else:
            self.ivalue = mpz_get_si(value)

    cdef void set_from_long(self, long value):
        self.ivalue = value % self.__modulus.int64

    cdef void set_from_int(IntegerMod_int64 self, int_fast64_t ivalue):
        if ivalue < 0:
            self.ivalue = self.__modulus.int64 + (ivalue % self.__modulus.int64) # Is ivalue % self.__modulus.int64 actually negative?
        elif ivalue >= self.__modulus.int64:
            self.ivalue = ivalue % self.__modulus.int64
        else:
            self.ivalue = ivalue

    cdef int_fast64_t get_int_value(IntegerMod_int64 self):
        return self.ivalue


    cdef int _cmp_c_impl(self, Element right) except -2:
        """
        EXAMPLES::

            sage: mod(5,13^5) == mod(13^5+5,13^5)
            True
            sage: mod(5,13^5) == mod(8,13^5)
            False
            sage: mod(5,13^5) == mod(5,13)
            True
            sage: mod(0, 13^5) == 0
            True
            sage: mod(0, 13^5) == int(0)
            True
        """
        if self.ivalue == (<IntegerMod_int64>right).ivalue: return 0
        elif self.ivalue < (<IntegerMod_int64>right).ivalue: return -1
        else: return 1

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)


    cpdef bint is_one(IntegerMod_int64 self):
        """
        Returns ``True`` if this is `1`, otherwise
        ``False``.

        EXAMPLES::

            sage: (mod(-1,5^10)^2).is_one()
            True
            sage: mod(0,5^10).is_one()
            False
        """
        return self.ivalue == 1

    def __nonzero__(IntegerMod_int64 self):
        """
        Returns ``True`` if this is not `0`, otherwise
        ``False``.

        EXAMPLES::

            sage: mod(13,5^10).is_zero()
            False
            sage: mod(5^12,5^10).is_zero()
            True
        """
        return self.ivalue != 0

    cpdef bint is_unit(IntegerMod_int64 self):
        """
        Return True iff this element is a unit.

        EXAMPLES::

            sage: mod(13, 5^10).is_unit()
            True
            sage: mod(25, 5^10).is_unit()
            False
        """
        return gcd_int64(self.ivalue, self.__modulus.int64) == 1

    def __crt(IntegerMod_int64 self, IntegerMod_int64 other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to self and
        to other. The modulus of other must be coprime to the modulus of
        self.

        EXAMPLES::

            sage: a = mod(3,5^10)
            sage: b = mod(2,7)
            sage: a.crt(b)
            29296878
            sage: type(a.crt(b)) == type(b.crt(a)) and type(a.crt(b)) == type(mod(1, 7 * 5^10))
            True

        ::

            sage: a = mod(3,10^10)
            sage: b = mod(2,9)
            sage: a.crt(b)
            80000000003
            sage: type(a.crt(b)) == type(b.crt(a)) and type(a.crt(b)) == type(mod(1, 9 * 10^10))
            True

        AUTHORS:

        - Robert Bradshaw
        """
        cdef IntegerMod_int64 lift
        cdef int_fast64_t x

        import integer_mod_ring
        lift = IntegerMod_int64(integer_mod_ring.IntegerModRing(self.__modulus.int64 * other.__modulus.int64), None, empty=True)

        try:
            x = (other.ivalue - self.ivalue % other.__modulus.int64) * mod_inverse_int64(self.__modulus.int64, other.__modulus.int64)
            lift.set_from_int( x * self.__modulus.int64 + self.ivalue )
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"

    def __copy__(IntegerMod_int64 self):
        return self._new_c(self.ivalue)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(7) + R(8)
            15
        """
        cdef int_fast64_t x
        x = self.ivalue + (<IntegerMod_int64>right).ivalue
        if x >= self.__modulus.int64:
            x = x - self.__modulus.int64
        return self._new_c(x)

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(7) + R(8)
            15
        """
        cdef int_fast64_t x
        x = self.ivalue + (<IntegerMod_int64>right).ivalue
        if x >= self.__modulus.int64:
            x = x - self.__modulus.int64
        self.ivalue = x
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(7) - R(8)
            99999
        """
        cdef int_fast64_t x
        x = self.ivalue - (<IntegerMod_int64>right).ivalue
        if x < 0:
            x = x + self.__modulus.int64
        return self._new_c(x)

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(7) - R(8)
            99999
        """
        cdef int_fast64_t x
        x = self.ivalue - (<IntegerMod_int64>right).ivalue
        if x < 0:
            x = x + self.__modulus.int64
        self.ivalue = x
        return self

    cpdef ModuleElement _neg_(self):
        """
        EXAMPLES::

            sage: -mod(7,10^5)
            99993
            sage: -mod(0,10^6)
            0
        """
        if self.ivalue == 0:
            return self
        return self._new_c(self.__modulus.int64 - self.ivalue)

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(700) * R(800)
            60000
        """
        return self._new_c((self.ivalue * (<IntegerMod_int64>right).ivalue) % self.__modulus.int64)


    cpdef RingElement _imul_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(700) * R(800)
            60000
        """
        self.ivalue = (self.ivalue * (<IntegerMod_int64>right).ivalue) % self.__modulus.int64
        return self

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(2)/3
            33334
        """
        return self._new_c((self.ivalue * mod_inverse_int64((<IntegerMod_int64>right).ivalue,
                                   self.__modulus.int64) ) % self.__modulus.int64)

    def __int__(IntegerMod_int64 self):
        return self.ivalue

    def __index__(self):
        """
        Needed so integers modulo `n` can be used as list indices.

        EXAMPLES::

            sage: v = [1,2,3,4,5]
            sage: v[Mod(3, 2^20)]
            4
        """
        return self.ivalue

    def __long__(IntegerMod_int64 self):
        return self.ivalue

    def __mod__(IntegerMod_int64 self, right):
        right = int(right)
        if self.__modulus.int64 % right != 0:
            raise ZeroDivisionError, "reduction modulo right not defined."
        import integer_mod_ring
        return integer_mod_ring.IntegerModRing(right)(self)

    def __lshift__(IntegerMod_int64 self, k):
        r"""
        Performs a left shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(5, 2^31 - 1)
            sage: e << 32
            10
            sage: e * 2^32
            10
        """
        return self.shift(int(k))

    def __rshift__(IntegerMod_int64 self, k):
        r"""
        Performs a right shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(5, 2^31 - 1)
            sage: e >> 1
            2
        """
        return self.shift(-int(k))

    cdef shift(IntegerMod_int64 self, int k):
        """
        Performs a bit-shift specified by ``k`` on ``self``.

        Suppose that ``self`` represents an integer `x` modulo `n`.  If `k` is
        `k = 0`, returns `x`.  If `k > 0`, shifts `x` to the left, that is,
        multiplies `x` by `2^k` and then returns the representative in the
        range `[0,n)`.  If `k < 0`, shifts `x` to the right, that is, returns
        the integral part of `x` divided by `2^k`.

        Note that, in any case, ``self`` remains unchanged.

        INPUT:

        - ``k`` - Integer of type ``int``

        OUTPUT:

        - Result of type ``IntegerMod_int64``

        WARNING:

        For positive ``k``, if ``x << k`` overflows as a 64-bit integer, the
        result is meaningless.

        EXAMPLES::

            sage: e = Mod(5, 2^31 - 1)
            sage: e << 32
            10
            sage: e * 2^32
            10
            sage: e = Mod(5, 2^31 - 1)
            sage: e >> 1
            2
        """
        if k == 0:
            return self
        elif k > 0:
            return self._new_c((self.ivalue << k) % self.__modulus.int64)
        else:
            return self._new_c(self.ivalue >> (-k))

    def __pow__(IntegerMod_int64 self, exp, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(2)^10
            4
            sage: p = next_prime(10^5)
            sage: R = Integers(p)
            sage: R(1234)^(p-1)
            1
            sage: R = Integers(17^5)
            sage: R(17)^5
            0

            sage: R(2)^-1 * 2
            1
            sage: R(2)^-1000000 * 2^1000000
            1
            sage: R(17)^-1
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
            sage: R(17)^-100000000
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        TESTS::

            sage: type(R(0))
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>

        We define ``0^0`` to be unity, :trac:`13894`::

            sage: p = next_prime(10^5)
            sage: R = Integers(p)
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: type(R(0)^0) == type(R(0))
            True

        When the modulus is ``1``, the only element in the ring is
        ``0`` (and it is equivalent to ``1``), so we return that
        instead::

            sage: from sage.rings.finite_rings.integer_mod \
            ...       import IntegerMod_int64
            sage: zero = IntegerMod_int64(Integers(1),0)
            sage: type(zero)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
            sage: zero^0
            0

        """
        cdef long long_exp
        cdef int_fast64_t res
        cdef mpz_t res_mpz
        if PyInt_CheckExact(exp) and -100000 < PyInt_AS_LONG(exp) < 100000:
            long_exp = PyInt_AS_LONG(exp)
        elif PY_TYPE_CHECK_EXACT(exp, Integer) and mpz_cmpabs_ui((<Integer>exp).value, 100000) == -1:
            long_exp = mpz_get_si((<Integer>exp).value)
        else:
            sig_on()
            try:
                mpz_init(res_mpz)
                base = self.lift()
                mpz_pow_helper(res_mpz, (<Integer>base).value, exp, self.__modulus.sageInteger.value)
                if mpz_fits_ulong_p(res_mpz):
                    res = mpz_get_ui(res_mpz)
                else:
                    res = mpz_get_pyintlong(res_mpz)
                return self._new_c(res)
            finally:
                mpz_clear(res_mpz)
                sig_off()

        if long_exp == 0 and self.ivalue == 0:
            # Return 0 if the modulus is 1, otherwise return 1.
            return self._new_c(self.__modulus.int64 != 1)
        cdef bint invert = False
        if long_exp < 0:
            invert = True
            long_exp = -long_exp
        res = mod_pow_int64(self.ivalue, long_exp, self.__modulus.int64)
        if invert:
            return self._new_c(mod_inverse_int64(res, self.__modulus.int64))
        else:
            return self._new_c(res)

    def __invert__(IntegerMod_int64 self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES::

            sage: a = mod(7,2^40); type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: ~a
            471219269047
            sage: a
            7
        """
        return self._new_c(mod_inverse_int64(self.ivalue, self.__modulus.int64))

    def lift(IntegerMod_int64 self):
        """
        Lift an integer modulo `n` to the integers.

        EXAMPLES::

            sage: a = Mod(8943, 2^25); type(a)
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
            sage: lift(a)
            8943
            sage: a.lift()
            8943
        """
        cdef sage.rings.integer.Integer z
        z = sage.rings.integer.Integer()
        mpz_set_si(z.value, self.ivalue)
        return z

    def __float__(IntegerMod_int64 self):
        """
        Coerce self to a float.

        EXAMPLES::

            sage: a = Mod(8943, 2^35)
            sage: float(a)
            8943.0
        """
        return <double>self.ivalue

    def __hash__(self):
        """
        Compute hash of self.

        EXAMPLES::

            sage: a = Mod(8943, 2^35)
            sage: hash(a)
            8943
        """

        return hash(self.ivalue)

    def _balanced_abs(self):
        """
        This function returns `x` or `-x`, whichever has a
        positive representative in `-n/2 < x \leq n/2`.
        """
        if self.ivalue > self.__modulus.int64 / 2:
            return -self
        else:
            return self

    @coerce_binop
    def gcd(self, IntegerMod_int64 other):
        """
        Greatest common divisor

        Returns the "smallest" generator in `\ZZ / N\ZZ` of the ideal
        generated by ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of the same ring as this one.

        EXAMPLES::

            sage: mod(2^3*3^2*5, 3^3*2^2*17^5).gcd(mod(2^4*3*17, 3^3*2^2*17^5))
            12
            sage: mod(0,17^5).gcd(mod(0,17^5))
            0
        """
        cdef int_fast64_t g = gcd_int64(self.ivalue, self.__modulus.int64)
        g = gcd_int64(g, other.ivalue)
        if g == self.__modulus.int64: # self = other = 0
            g = 0
        return self._new_c(g)

### Helper functions

cdef mpz_pow_helper(mpz_t res, mpz_t base, object exp, mpz_t modulus):
    cdef bint invert = False
    cdef long long_exp

    if PyInt_CheckExact(exp):
        long_exp = PyInt_AS_LONG(exp)
        if long_exp < 0:
            long_exp = -long_exp
            invert = True
        mpz_powm_ui(res, base, long_exp, modulus)
    else:
        if not PY_TYPE_CHECK_EXACT(exp, Integer):
            exp = Integer(exp)
        if mpz_sgn((<Integer>exp).value) < 0:
            exp = -exp
            invert = True
        mpz_powm(res, base, (<Integer>exp).value, modulus)
    if invert:
        if not mpz_invert(res, res, modulus):
            raise ZeroDivisionError, "Inverse does not exist."

cdef int_fast64_t gcd_int64(int_fast64_t a, int_fast64_t b):
    """
    Returns the gcd of a and b

    For use with IntegerMod_int64

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int_fast64_t tmp
    if a < b:
        tmp = b
        b = a
        a = tmp
    while b:
        tmp = b
        b = a % b
        a = tmp
    return a


cdef int_fast64_t mod_inverse_int64(int_fast64_t x, int_fast64_t n) except 0:
    """
    Returns y such that xy=1 mod n

    For use in IntegerMod_int64

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int_fast64_t tmp, a, b, last_t, t, next_t, q
    a = n
    b = x
    t = 0
    next_t = 1
    while b:
        # a = s * n + t * x
        if b == 1:
            next_t = next_t % n
            if next_t < 0:
                next_t = next_t + n
            return next_t
        q = a / b
        tmp = b
        b = a % b
        a = tmp
        last_t = t
        t = next_t
        next_t = last_t - q * t
    raise ZeroDivisionError, "Inverse does not exist."


cdef int_fast64_t mod_pow_int64(int_fast64_t base, int_fast64_t exp, int_fast64_t n):
    """
    Returns base^exp mod n

    For use in IntegerMod_int64

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int_fast64_t prod, pow2
    if exp <= 5:
        if exp == 0: return 1
        if exp == 1: return base
        prod = base * base % n
        if exp == 2: return prod
        if exp == 3: return (prod * base) % n
        if exp == 4: return (prod * prod) % n

    pow2 = base
    if exp % 2: prod = base
    else: prod = 1
    exp = exp >> 1
    while(exp != 0):
        pow2 = pow2 * pow2
        if pow2 >= INTEGER_MOD_INT64_LIMIT: pow2 = pow2 % n
        if exp % 2:
            prod = prod * pow2
            if prod >= INTEGER_MOD_INT64_LIMIT: prod = prod % n
        exp = exp >> 1

    if prod >= n:
        prod = prod % n
    return prod


cdef int jacobi_int64(int_fast64_t a, int_fast64_t m) except -2:
    """
    Calculates the jacobi symbol (a/n)

    For use in IntegerMod_int64

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int s, jacobi = 1
    cdef int_fast64_t b

    a = a % m

    while 1:
        if a == 0:
            return 0 # gcd was nontrivial
        elif a == 1:
            return jacobi
        s = 0
        while (1 << s) & a == 0:
            s += 1
        b = a >> s
        # Now a = 2^s * b

        # factor out (2/m)^s term
        if s % 2 == 1 and (m % 8 == 3 or m % 8 == 5):
            jacobi = -jacobi

        if b == 1:
            return jacobi

        # quadratic reciprocity
        if b % 4 == 3 and m % 4 == 3:
            jacobi = -jacobi
        a = m % b
        m = b


########################
# Square root functions
########################

def square_root_mod_prime_power(IntegerMod_abstract a, p, e):
    r"""
    Calculates the square root of `a`, where `a` is an
    integer mod `p^e`.

    ALGORITHM: Perform `p`-adically by stripping off even
    powers of `p` to get a unit and lifting
    `\sqrt{unit} \bmod p` via Newton's method.

    AUTHORS:

    - Robert Bradshaw

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import square_root_mod_prime_power
        sage: a=Mod(17,2^20)
        sage: b=square_root_mod_prime_power(a,2,20)
        sage: b^2 == a
        True

    ::

        sage: a=Mod(72,97^10)
        sage: b=square_root_mod_prime_power(a,97,10)
        sage: b^2 == a
        True
    """
    if a.is_zero() or a.is_one():
        return a

    if p == 2:
        if e == 1:
            return a
        # TODO: implement something that isn't totally idiotic.
        for x in a.parent():
            if x**2 == a:
                return x

    # strip off even powers of p
    cdef int i, val = a.lift().valuation(p)
    if val % 2 == 1:
        raise ValueError, "self must be a square."
    if val > 0:
        unit = a._parent(a.lift() // p**val)
    else:
        unit = a

    # find square root of unit mod p
    x = unit.parent()(square_root_mod_prime(mod(unit, p), p))

    # lift p-adically using Newton iteration
    # this is done to higher precision than necessary except at the last step
    one_half = ~(a._new_c_from_long(2))
    for i from 0 <= i <  ceil(log(e)/log(2)) - val/2:
        x = (x+unit/x) * one_half

    # multiply in powers of p (if any)
    if val > 0:
        x *= p**(val // 2)
    return x

cpdef square_root_mod_prime(IntegerMod_abstract a, p=None):
    r"""
    Calculates the square root of `a`, where `a` is an
    integer mod `p`; if `a` is not a perfect square,
    this returns an (incorrect) answer without checking.

    ALGORITHM: Several cases based on residue class of
    `p \bmod 16`.


    -  `p \bmod 2 = 0`: `p = 2` so
       `\sqrt{a} = a`.

    -  `p \bmod 4 = 3`: `\sqrt{a} = a^{(p+1)/4}`.

    -  `p \bmod 8 = 5`: `\sqrt{a} = \zeta i a` where
       `\zeta = (2a)^{(p-5)/8}`, `i=\sqrt{-1}`.

    -  `p \bmod 16 = 9`: Similar, work in a bi-quadratic
       extension of `\GF{p}` for small `p`, Tonelli
       and Shanks for large `p`.

    -  `p \bmod 16 = 1`: Tonelli and Shanks.


    REFERENCES:

    - Siguna Muller.  'On the Computation of Square Roots in Finite
      Fields' Designs, Codes and Cryptography, Volume 31, Issue 3
      (March 2004)

    - A. Oliver L. Atkin. 'Probabilistic primality testing' (Chapter
      30, Section 4) In Ph. Flajolet and P. Zimmermann, editors,
      Algorithms Seminar, 1991-1992. INRIA Research Report 1779, 1992,
      http://www.inria.fr/rrrt/rr-1779.html. Summary by F. Morain.
      http://citeseer.ist.psu.edu/atkin92probabilistic.html

    - H. Postl. 'Fast evaluation of Dickson Polynomials' Contrib. to
      General Algebra, Vol. 6 (1988) pp. 223-225

    AUTHORS:

    - Robert Bradshaw

    TESTS: Every case appears in the first hundred primes.

    ::

        sage: from sage.rings.finite_rings.integer_mod import square_root_mod_prime   # sqrt() uses brute force for small p
        sage: all([square_root_mod_prime(a*a)^2 == a*a
        ...        for p in prime_range(100)
        ...        for a in Integers(p)])
        True
    """
    if not a or a.is_one():
        return a

    if p is None:
        p = a._parent.order()
    if p < PyInt_GetMax():
        p = int(p)

    cdef int p_mod_16 = p % 16
    cdef double bits = log(float(p))/log(2)
    cdef long r, m

    cdef Integer resZ


    if p_mod_16 % 2 == 0:  # p == 2
        return a

    elif p_mod_16 % 4 == 3:
        return a ** ((p+1)//4)

    elif p_mod_16 % 8 == 5:
        two_a = a+a
        zeta = two_a ** ((p-5)//8)
        i = zeta**2 * two_a # = two_a ** ((p-1)//4)
        return zeta*a*(i-1)

    elif p_mod_16 == 9 and bits < 500:
        two_a = a+a
        s = two_a ** ((p-1)//4)
        if s.is_one():
            d = a._parent.quadratic_nonresidue()
            d2 = d*d
            z = (two_a * d2) ** ((p-9)//16)
            i = two_a * d2 * z*z
            return z*d*a*(i-1)
        else:
            z = two_a ** ((p-9)//16)
            i = two_a * z*z
            return z*a*(i-1)

    else:
        one = a._new_c_from_long(1)
        r, q = (p-one_Z).val_unit(2)
        v = a._parent.quadratic_nonresidue()**q

        x = a ** ((q-1)//2)
        b = a*x*x # a ^ q
        res = a*x # a ^ ((q-1)/2)

        while b != one:
            m = 1
            bpow = b*b
            while bpow != one:
                bpow *= bpow
                m += 1
            g = v**(one_Z << (r-m-1)) # v^(2^(r-m-1))
            res *= g
            b *= g*g
        return res


def fast_lucas(mm, IntegerMod_abstract P):
    """
    Return `V_k(P, 1)` where `V_k` is the Lucas
    function defined by the recursive relation

    `V_k(P, Q) = PV_{k-1}(P, Q) -  QV_{k-2}(P, Q)`

    with `V_0 = 2, V_1(P_Q) = P`.

    REFERENCES:

    - H. Postl. 'Fast evaluation of Dickson Polynomials' Contrib. to
      General Algebra, Vol. 6 (1988) pp. 223-225

    AUTHORS:

    - Robert Bradshaw

    TESTS::

        sage: from sage.rings.finite_rings.integer_mod import fast_lucas, slow_lucas
        sage: all([fast_lucas(k, a) == slow_lucas(k, a)
        ...        for a in Integers(23)
        ...        for k in range(13)])
        True
    """
    if mm == 0:
        return 2
    elif mm == 1:
        return P

    cdef sage.rings.integer.Integer m
    m = <sage.rings.integer.Integer>mm if PY_TYPE_CHECK(mm, sage.rings.integer.Integer) else sage.rings.integer.Integer(mm)
    two = P._new_c_from_long(2)
    d1 = P
    d2 = P*P - two

    sig_on()
    cdef int j
    for j from mpz_sizeinbase(m.value, 2)-1 > j > 0:
        if mpz_tstbit(m.value, j):
            d1 = d1*d2 - P
            d2 = d2*d2 - two
        else:
            d2 = d1*d2 - P
            d1 = d1*d1 - two
    sig_off()
    if mpz_odd_p(m.value):
        return d1*d2 - P
    else:
        return d1*d1 - two

def slow_lucas(k, P, Q=1):
    """
    Lucas function defined using the standard definition, for
    consistency testing.
    """
    if k == 0:
        return 2
    elif k == 1:
        return P
    else:
        return P*slow_lucas(k-1, P, Q) - Q*slow_lucas(k-2, P, Q)


############# Homomorphisms ###############

cdef class IntegerMod_hom(Morphism):
    cdef IntegerMod_abstract zero
    cdef NativeIntStruct modulus
    def __init__(self, parent):
        Morphism.__init__(self, parent)
        # we need to use element constructor so that we can register both coercions and conversions using these morphisms.
        self.zero = self._codomain._element_constructor_(0)
        self.modulus = self._codomain._pyx_order
    cpdef Element _call_(self, x):
        return IntegerMod(self.codomain(), x)

cdef class IntegerMod_to_IntegerMod(IntegerMod_hom):
    """
    Very fast IntegerMod to IntegerMod homomorphism.

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import IntegerMod_to_IntegerMod
        sage: Rs = [Integers(3**k) for k in range(1,30,5)]
        sage: [type(R(0)) for R in Rs]
        [<type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>]
        sage: fs = [IntegerMod_to_IntegerMod(S, R) for R in Rs for S in Rs if S is not R and S.order() > R.order()]
        sage: all([f(-1) == f.codomain()(-1) for f in fs])
        True
        sage: [f(-1) for f in fs]
        [2, 2, 2, 2, 2, 728, 728, 728, 728, 177146, 177146, 177146, 43046720, 43046720, 10460353202]
    """
    def __init__(self, R, S):
        if not S.order().divides(R.order()):
            raise TypeError, "No natural coercion from %s to %s" % (R, S)
        import sage.categories.homset
        IntegerMod_hom.__init__(self, sage.categories.homset.Hom(R, S))

    cpdef Element _call_(self, x):
        cdef IntegerMod_abstract a
        if PY_TYPE_CHECK(x, IntegerMod_int):
            return (<IntegerMod_int>self.zero)._new_c((<IntegerMod_int>x).ivalue % self.modulus.int32)
        elif PY_TYPE_CHECK(x, IntegerMod_int64):
            return self.zero._new_c_from_long((<IntegerMod_int64>x).ivalue  % self.modulus.int64)
        else: # PY_TYPE_CHECK(x, IntegerMod_gmp)
            a = self.zero._new_c_from_long(0)
            a.set_from_mpz((<IntegerMod_gmp>x).value)
            return a

    def _repr_type(self):
        return "Natural"

cdef class Integer_to_IntegerMod(IntegerMod_hom):
    r"""
    Fast `\ZZ \rightarrow \ZZ/n\ZZ`
    morphism.

    EXAMPLES:

    We make sure it works for every type.

    ::

        sage: from sage.rings.finite_rings.integer_mod import Integer_to_IntegerMod
        sage: Rs = [Integers(10), Integers(10^5), Integers(10^10)]
        sage: [type(R(0)) for R in Rs]
        [<type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>]
        sage: fs = [Integer_to_IntegerMod(R) for R in Rs]
        sage: [f(-1) for f in fs]
        [9, 99999, 9999999999]
    """
    def __init__(self, R):
        import sage.categories.homset
        IntegerMod_hom.__init__(self, sage.categories.homset.Hom(integer_ring.ZZ, R))

    cpdef Element _call_(self, x):
        cdef IntegerMod_abstract a
        cdef Py_ssize_t res
        if self.modulus.table is not None:
            res = x % self.modulus.int64
            if res < 0:
                res += self.modulus.int64
            a = self.modulus.lookup(res)
            if a._parent is not self._codomain:
               a._parent = self._codomain
#                print (<Element>a)._parent, " is not ", parent
            return a
        else:
            a = self.zero._new_c_from_long(0)
            a.set_from_mpz((<Integer>x).value)
            return a

    def _repr_type(self):
        return "Natural"

    def section(self):
        return IntegerMod_to_Integer(self._codomain)

cdef class IntegerMod_to_Integer(Map):
    def __init__(self, R):
        import sage.categories.homset
        from sage.categories.all import Sets
        Morphism.__init__(self, sage.categories.homset.Hom(R, integer_ring.ZZ, Sets()))

    cpdef Element _call_(self, x):
        cdef Integer ans = PY_NEW(Integer)
        if PY_TYPE_CHECK(x, IntegerMod_gmp):
            mpz_set(ans.value, (<IntegerMod_gmp>x).value)
        elif PY_TYPE_CHECK(x, IntegerMod_int):
            mpz_set_si(ans.value, (<IntegerMod_int>x).ivalue)
        elif PY_TYPE_CHECK(x, IntegerMod_int64):
            mpz_set_si(ans.value, (<IntegerMod_int64>x).ivalue)
        return ans

    def _repr_type(self):
        return "Lifting"

cdef class Int_to_IntegerMod(IntegerMod_hom):
    """
    EXAMPLES:

    We make sure it works for every type.

    ::

        sage: from sage.rings.finite_rings.integer_mod import Int_to_IntegerMod
        sage: Rs = [Integers(2**k) for k in range(1,50,10)]
        sage: [type(R(0)) for R in Rs]
        [<type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>, <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>]
        sage: fs = [Int_to_IntegerMod(R) for R in Rs]
        sage: [f(-1) for f in fs]
        [1, 2047, 2097151, 2147483647, 2199023255551]
    """
    def __init__(self, R):
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        IntegerMod_hom.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), R))

    cpdef Element _call_(self, x):
        cdef IntegerMod_abstract a
        cdef long res = PyInt_AS_LONG(x)
        if PY_TYPE_CHECK(self.zero, IntegerMod_gmp):
            if 0 <= res < INTEGER_MOD_INT64_LIMIT:
                return self.zero._new_c_from_long(res)
            else:
                return IntegerMod_gmp(self.zero._parent, x)
        else:
            res %= self.modulus.int64
            if res < 0:
                res += self.modulus.int64
            if self.modulus.table is not None:
                a = self.modulus.lookup(res)
                if a._parent is not self._codomain:
                   a._parent = self._codomain
    #                print (<Element>a)._parent, " is not ", parent
                return a
            else:
                return self.zero._new_c_from_long(res)

    def _repr_type(self):
        return "Native"
