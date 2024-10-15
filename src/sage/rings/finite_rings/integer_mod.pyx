r"""
Elements of `\ZZ/n\ZZ`

An element of the integers modulo `n`.

There are three types of integer_mod classes, depending on the
size of the modulus.


- ``IntegerMod_int`` stores its value in a
  ``int_fast32_t`` (typically an ``int``);
  this is used if the modulus is less than
  `\sqrt{2^{31}-1}`.

- ``IntegerMod_int64`` stores its value in a
  ``int_fast64_t`` (typically a ``long
  long``); this is used if the modulus is less than
  `2^{31}-1`. In many places, we assume that the values and the modulus
  actually fit inside an ``unsigned long``.

- ``IntegerMod_gmp`` stores its value in a
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
    <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
    sage: type(IntegerModRing(2^31).an_element())
    <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
"""
# ****************************************************************************
#       Copyright (C) 2006 Robert Bradshaw <robertwb@math.washington.edu>
#                     2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off, sig_check

from cpython.long cimport *
from cpython.list cimport *
from cpython.ref cimport *

from libc.math cimport log2, ceil

from sage.libs.gmp.all cimport *

cdef bint use_32bit_type(int_fast64_t modulus) noexcept:
    return modulus <= INTEGER_MOD_INT32_LIMIT

from sage.arith.long cimport (
    integer_check_long, integer_check_long_py, is_small_python_int)

import sage.rings.rational as rational

try:
    from sage.libs.pari.all import pari, PariError
except ImportError:
    class PariError(Exception):
        pass

import sage.rings.integer_ring as integer_ring
import sage.rings.rational_field

import sage.rings.integer
cimport sage.rings.integer
from sage.rings.integer cimport Integer

from sage.structure.coerce cimport py_scalar_to_element
from sage.structure.richcmp cimport rich_to_bool_sgn, rich_to_bool
import sage.structure.element
cimport sage.structure.element
coerce_binop = sage.structure.element.coerce_binop
from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map

from sage.misc.persist import register_unpickle_override

from sage.structure.parent cimport Parent

from sage.arith.misc import CRT as crt
from sage.arith.functions import lcm


cdef Integer one_Z = Integer(1)


def Mod(n, m, parent=None):
    r"""
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

    Illustrates that :issue:`5971` is fixed. Consider `n` modulo `m` when
    `m = 0`. Then `\ZZ/0\ZZ` is isomorphic to `\ZZ` so `n` modulo `0`
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

    # m is nonzero, so return n mod m
    if parent is None:
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
        parent = IntegerModRing(m)
    return IntegerMod(parent, n)


mod = Mod

register_unpickle_override('sage.rings.integer_mod', 'Mod', Mod)
register_unpickle_override('sage.rings.integer_mod', 'mod', mod)


def IntegerMod(parent, value):
    """
    Create an integer modulo `n` with the given parent.

    This is mainly for internal use.

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import IntegerMod
        sage: R = IntegerModRing(100)
        sage: type(R._pyx_order.table)
        <class 'list'>
        sage: IntegerMod(R, 42)
        42
        sage: IntegerMod(R, 142)
        42
        sage: IntegerMod(R, 10^100 + 42)
        42
        sage: IntegerMod(R, -9158)
        42
    """
    cdef NativeIntStruct modulus = parent._pyx_order

    cdef long val = 0
    cdef int err

    if modulus.table is not None:
        # Try to return an element from the precomputed table
        integer_check_long(value, &val, &err)
        if not err:
            val = (<int_fast64_t>val) % modulus.int64
            if val < 0:
                val += modulus.int64
            a = <Element>modulus.table[val]
            assert a._parent is parent
            return a
    t = modulus.element_class()
    return t(parent, value)


def is_IntegerMod(x):
    """
    Return ``True`` if and only if x is an integer modulo
    `n`.

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import is_IntegerMod
        sage: is_IntegerMod(5)
        doctest:warning...
        DeprecationWarning: The function is_IntegerMod is deprecated;
        use 'isinstance(..., IntegerMod_abstract)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        False
        sage: is_IntegerMod(Mod(5,10))
        True
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_IntegerMod is deprecated; "
                       "use 'isinstance(..., IntegerMod_abstract)' instead.")
    return isinstance(x, IntegerMod_abstract)


cdef inline inverse_or_None(x):
    try:
        return ~x
    except ArithmeticError:
        return None


cdef class NativeIntStruct:
    """
    We store the various forms of the modulus here rather than in the
    parent for efficiency reasons.

    We may also store a cached table of all elements of a given ring in
    this class.
    """
    def __cinit__(self):
        self.int32 = -1
        self.int64 = -1

    def __init__(self, m):
        self.sageInteger = Integer(m)
        cdef mpz_srcptr z = self.sageInteger.value
        if mpz_cmp_si(z, INTEGER_MOD_INT64_LIMIT) <= 0:
            self.int64 = mpz_get_si(z)
            if use_32bit_type(self.int64):
                self.int32 = self.int64

    def __repr__(self):
        return f"{type(self).__name__}({self.sageInteger})"

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.rings.finite_rings.integer_mod import NativeIntStruct
            sage: M = NativeIntStruct(12345); M
            NativeIntStruct(12345)
            sage: loads(dumps(M))
            NativeIntStruct(12345)
        """
        return type(self), (self.sageInteger, )

    def precompute_table(self, parent):
        """
        Function to compute and cache all elements of this class.

        If ``inverses == True``, also computes and caches the inverses
        of the invertible elements.

        EXAMPLES::

            sage: from sage.rings.finite_rings.integer_mod import NativeIntStruct
            sage: R = IntegerModRing(10)
            sage: M = NativeIntStruct(R.order())
            sage: M.precompute_table(R)
            sage: M.table
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: M.inverses
            [None, 1, None, 7, None, None, None, 3, None, 9]

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

        Check that the inverse of 0 modulo 1 works, see :issue:`13639`::

            sage: R = IntegerModRing_generic(1, cache=True)  # indirect doctest
            sage: R(0)^-1 is R(0)
            True

        TESTS::

            sage: R = IntegerModRing(10^50)
            sage: M = NativeIntStruct(R.order())
            sage: M.precompute_table(R)
            Traceback (most recent call last):
            ...
            OverflowError: precompute_table() is only supported for small moduli
        """
        cdef Py_ssize_t i, m = self.int64

        # Verify that the modulus m fits in a Py_ssize_t
        if m <= 0 or (<int_fast64_t>m != self.int64):
            raise OverflowError("precompute_table() is only supported for small moduli")

        t = self.element_class()
        self.table = [t(parent, i) for i in range(m)]

        if m == 1:
            # Special case for integers modulo 1
            self.inverses = self.table
        else:
            self.inverses = [inverse_or_None(x) for x in self.table]


# For unpickling
makeNativeIntStruct = NativeIntStruct
register_unpickle_override('sage.rings.integer_mod', 'makeNativeIntStruct', NativeIntStruct)


cdef class IntegerMod_abstract(FiniteRingElement):

    def __init__(self, parent, value=None):
        """
        EXAMPLES::

            sage: a = Mod(10, 30^10); a
            10
            sage: type(a)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: loads(a.dumps()) == a
            True

        TESTS::

            sage: TestSuite(Zmod(1)).run()
            sage: TestSuite(Zmod(2)).run()
            sage: TestSuite(Zmod(3)).run()
            sage: TestSuite(Zmod(4)).run()
            sage: TestSuite(Zmod(5)).run()
            sage: TestSuite(Zmod(6)).run()
            sage: TestSuite(Zmod(2^10 * 3^5)).run()
            sage: TestSuite(Zmod(2^30 * 3^50 * 5^20)).run()

            sage: GF(29)(SR(1/3))                                                       # needs sage.rings.finite_rings sage.symbolic
            10
            sage: Integers(30)(QQ['x'](1/7))
            13
            sage: Integers(30)(SR(1/4))                                                 # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(4, 30) does not exist
        """
        self._parent = parent
        self._modulus = parent._pyx_order

        if value is None:
            return

        cdef long longval = 0
        cdef int err = 0
        cdef Integer z

        if isinstance(value, Integer):
            z = <Integer>value
        elif isinstance(value, rational.Rational):
            z = value % self._modulus.sageInteger
        elif integer_check_long_py(value, &longval, &err) and not err:
            self.set_from_long(longval)
            return
        else:
            try:
                z = integer_ring.Z(value)
            except (TypeError, ValueError):
                from sage.structure.element import Expression
                if isinstance(value, Expression):
                    value = value.pyobject()
                else:
                    value = py_scalar_to_element(value)
                if isinstance(value, Element) and value.parent().is_exact():
                    value = sage.rings.rational_field.QQ(value)
                    z = value % self._modulus.sageInteger
                else:
                    raise
        self.set_from_mpz(z.value)

    cdef IntegerMod_abstract _new_c_fast(self, unsigned long value):
        cdef type t = type(self)
        x = <IntegerMod_abstract>t.__new__(t)
        x._parent = self._parent
        x._modulus = self._modulus
        x.set_from_ulong_fast(value)
        return x

    cdef _new_c_from_long(self, long value):
        cdef type t = type(self)
        cdef IntegerMod_abstract x = <IntegerMod_abstract>t.__new__(t)
        x._parent = self._parent
        x._modulus = self._modulus
        x.set_from_long(value)
        return x

    cdef void set_from_mpz(self, mpz_t value) noexcept:
        raise NotImplementedError("must be defined in child class")

    cdef void set_from_long(self, long value) noexcept:
        raise NotImplementedError("must be defined in child class")

    cdef void set_from_ulong_fast(self, unsigned long value) noexcept:
        """
        Set ``self`` to the value in ``value`` where ``value`` is
        assumed to be less than the modulus
        """
        raise NotImplementedError("must be defined in child class")

    def __abs__(self):
        """
        Raise an error message, since ``abs(x)`` makes no sense
        when ``x`` is an integer modulo `n`.

        EXAMPLES::

            sage: abs(Mod(2,3))
            Traceback (most recent call last):
            ...
            ArithmeticError: absolute value not defined on integers modulo n.
        """
        raise ArithmeticError("absolute value not defined on integers modulo n.")

    def __reduce__(IntegerMod_abstract self):
        """
        EXAMPLES::

            sage: a = Mod(4,5); a
            4
            sage: loads(a.dumps()) == a
            True
            sage: a = Mod(-1,5^30)^25
            sage: loads(a.dumps()) == a
            True
        """
        return sage.rings.finite_rings.integer_mod.mod, (self.lift(), self.modulus(), self.parent())

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the map that sends the
        generators of the parent to ``im_gens``.

        EXAMPLES::

            sage: a = Mod(7, 10)
            sage: R = ZZ.quotient(5)
            sage: a._im_gens_(R, (R(1),))
            2
        """
        # The generators are irrelevant (Zmod(n) is its own base), so we ignore base_map
        return codomain.coerce(self)

    def __mod__(self, modulus):
        """
        Coerce this element to the ring `Z/(modulus)`.

        If the new ``modulus`` does not divide the current modulus,
        an :exc:`ArithmeticError` is raised.

        EXAMPLES::

            sage: a = Mod(14, 35)
            sage: a % 5
            4
            sage: parent(a % 5)
            Ring of integers modulo 5
            sage: a % 350
            Traceback (most recent call last):
            ...
            ArithmeticError: reduction modulo 350 not defined
            sage: a % 35
            14
            sage: int(1) % a
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for %: 'int' and 'sage.rings.finite_rings.integer_mod.IntegerMod_int'
        """
        if not isinstance(self, IntegerMod_abstract):
            # something % Mod(x,y) makes no sense
            return NotImplemented
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
        R = IntegerModRing(modulus)
        if (<Element>self)._parent._IntegerModRing_generic__order % R.order():
            raise ArithmeticError(f"reduction modulo {modulus!r} not defined")
        return R(self)

    def is_nilpotent(self):
        r"""
        Return ``True`` if ``self`` is nilpotent,
        i.e., some power of ``self`` is zero.

        EXAMPLES::

            sage: a = Integers(90384098234^3)
            sage: factor(a.order())                                                     # needs sage.libs.pari
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
        m = self._modulus.sageInteger.exact_log(2) + 1
        return (self**m).is_zero()

    #################################################################
    # Interfaces
    #################################################################
    def _pari_init_(self):
        return 'Mod(%s,%s)'%(str(self), self._modulus.sageInteger)

    def __pari__(self):
        return self.lift().__pari__().Mod(self._modulus.sageInteger)

    def _gap_init_(self):
        r"""
        Return string representation of corresponding GAP object.

        EXAMPLES::

            sage: # needs sage.libs.gap
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
        return '%s*One(ZmodnZ(%s))' % (self, self._modulus.sageInteger)

    def _magma_init_(self, magma):
        """
        Coercion to Magma.

        EXAMPLES::

            sage: # optional - magma
            sage: a = Integers(15)(4)
            sage: b = magma(a)
            sage: b.Type()
            RngIntResElt
            sage: b^2
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

            sage: aa = axiom(a); aa             # optional - axiom
            4
            sage: aa.type()                     # optional - axiom
            IntegerMod 15

            sage: aa = fricas(a); aa            # optional - fricas
            4
            sage: aa.typeOf()                   # optional - fricas
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

    def log(self, b=None, order=None, check=False):
        r"""
        Compute the discrete logarithm of this element to base `b`,
        that is,
        return an integer `x` such that `b^x = a`, where
        `a` is ``self``.

        INPUT:

        - ``self`` -- unit modulo `n`

        - ``b`` -- a unit modulo `n`. If ``b`` is not given,
          ``R.multiplicative_generator()`` is used, where
          ``R`` is the parent of ``self``.

        - ``order`` -- integer (unused), the order of ``b``.
          This argument is normally unused, only there for
          coherence of apis with finite field elements.

        - ``check`` -- boolean (default: ``False``); if set,
          test whether the given ``order`` is correct

        OUTPUT:

        Integer `x` such that `b^x = a`, if this exists; a :exc:`ValueError`
        otherwise.

        .. NOTE::

           The algorithm first factors the modulus, then invokes Pari's :pari:`znlog`
           function for each odd prime power in the factorization of the modulus.
           This method can be quite slow for large moduli.

        EXAMPLES::

            sage: # needs sage.libs.pari sage.modules
            sage: r = Integers(125)
            sage: b = r.multiplicative_generator()^3
            sage: a = b^17
            sage: a.log(b)
            17
            sage: a.log()
            51

        A bigger example::

            sage: # needs sage.rings.finite_rings
            sage: FF = FiniteField(2^32 + 61)
            sage: c = FF(4294967356)
            sage: x = FF(2)
            sage: a = c.log(x)
            sage: a
            2147483678
            sage: x^a
            4294967356

        An example with a highly composite modulus::

            sage: m = 2^99 * 77^7 * 123456789 * 13712923537615486607^2
            sage: (Mod(5,m)^5735816763073854953388147237921).log(5)                     # needs sage.libs.pari
            5735816763073854953388147237921

        Errors are generated if the logarithm doesn't exist
        or the inputs are not units::

            sage: Mod(3, 7).log(Mod(2, 7))                                              # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: no logarithm of 3 found to base 2 modulo 7
            sage: a = Mod(16, 100); b = Mod(4, 100)
            sage: a.log(b)
            Traceback (most recent call last):
            ...
            ValueError: logarithm of 16 is not defined since it is not a unit modulo 100

        TESTS:

        We check that :issue:`9205` is fixed::

            sage: Mod(5, 9).log(Mod(2, 9))                                              # needs sage.libs.pari
            5

        We test against a bug (side effect on PARI) fixed in :issue:`9438`::

            sage: # needs sage.libs.pari
            sage: R.<a, b> = QQ[]
            sage: pari(b)
            b
            sage: GF(7)(5).log()
            5
            sage: pari(b)
            b

        We test that :issue:`23927` is fixed::

            sage: x = mod(48475563673907791151, 10^20 + 763)^2
            sage: e = 25248843418589594761
            sage: (x^e).log(x) == e                                                     # needs sage.libs.pari
            True

        Examples like this took extremely long before :issue:`32375`::

            sage: (Mod(5, 123337052926643**4) ^ (10^50-1)).log(5)                       # needs sage.libs.pari
            99999999999999999999999999999999999999999999999999

        We check that non-existence of solutions is detected:

        No local solutions::

            sage: Mod(1111, 1234567).log(1111**3)                                       # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: no logarithm of 1111 found to base 961261 modulo 1234567 (no solution modulo 9721)

        Incompatible local solutions::

            sage: Mod(230, 323).log(173)                                                # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: no logarithm of 230 found to base 173 modulo 323 (incompatible local solutions)

        We test that :issue:`12419` is fixed::

            sage: R.<x,y> = GF(2)[]
            sage: R(1).factor()
            1

        An example for ``check=True``::

            sage: F = GF(127, impl='modn')
            sage: t = F.primitive_element()                                             # needs sage.libs.pari
            sage: t.log(t, 57, check=True)                                              # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: base does not have the provided order

        AUTHORS:

        - David Joyner and William Stein (2005-11)

        - William Stein (2007-01-27): update to use PARI as requested
          by David Kohel.

        - Simon King (2010-07-07): fix a side effect on PARI

        - Lorenz Panny (2021): speedups for composite moduli
        """
        if not self.is_unit():
            raise ValueError(f"logarithm of {self} is not defined since it is not a unit modulo {self.modulus()}")

        if b is None:
            b = self._parent.multiplicative_generator()
        else:
            b = self._parent(b)
            if not b.is_unit():
                raise ValueError(f"logarithm with base {b} is not defined since it is not a unit modulo {b.modulus()}")

        if check:
            from sage.groups.generic import has_order
            if not has_order(b, order, '*'):
                raise ValueError('base does not have the provided order')

        cdef Integer n = Integer()
        cdef Integer m = one_Z
        cdef Integer q, na, nb

        for p, e in self.modulus().factor():
            q = p**e
            a_red = Mod(self.lift(), q)
            b_red = Mod(b.lift(), q)

            na = a_red.multiplicative_order()
            nb = b_red.multiplicative_order()
            if not na.divides(nb):  # cannot be a power
                raise ValueError(f"no logarithm of {self} found to base {b} modulo {self.modulus()}"
                              + (f" (no solution modulo {q})" if q != self.modulus() else ""))

            if p == 2 and e >= 3:   # (ZZ/2^e)* is not cyclic; must not give unsolvable DLPs to Pari
                try:
                    from sage.groups.generic import discrete_log
                    v = discrete_log(a_red, b_red, nb)
                except ValueError:
                    raise ValueError(f"no logarithm of {self} found to base {b} modulo {self.modulus()}"
                                  + (f" (no solution modulo {q})" if q != self.modulus() else ""))
            else:
                try:
                    v = pari(a_red).znlog(pari(b_red)).sage()
                except PariError as msg:
                    raise RuntimeError(f"{msg}\nPARI failed to compute discrete log modulo {q} (perhaps base is not a generator or is too large)")
                assert v != []  # if this happens, we've made a mistake above (or there is a Pari bug)

            try:
                n = crt(n, v, m, nb)
            except ValueError:
                raise ValueError(f"no logarithm of {self} found to base {b} modulo {self.modulus()} (incompatible local solutions)")
            m = lcm(m, nb)

#        assert b**n == self
        return n

    def generalised_log(self):
        r"""
        Return integers `[n_1, \ldots, n_d]` such that.

        .. MATH::

            \prod_{i=1}^d x_i^{n_i} = \text{self},

        where `x_1, \dots, x_d` are the generators of the unit group
        returned by ``self.parent().unit_gens()``.

        EXAMPLES::


            sage: m = Mod(3, 1568)
            sage: v = m.generalised_log(); v                                            # needs sage.libs.pari sage.modules
            [1, 3, 1]
            sage: prod([Zmod(1568).unit_gens()[i] ** v[i] for i in [0..2]])             # needs sage.libs.pari sage.modules
            3

        .. SEEALSO::

            The method :meth:`log`.

        .. warning::

            The output is given relative to the set of generators
            obtained by passing ``algorithm='sage'`` to the method
            :meth:`~sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic.unit_gens`
            of the parent (which is the default).  Specifying
            ``algorithm='pari'`` usually yields a different set of
            generators that is incompatible with this method.
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
        return self._modulus.sageInteger

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this element.

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
        Return the minimal polynomial of this element.

        EXAMPLES::

            sage: GF(241, 'a')(1).minpoly()
            x + 240
        """
        return self.charpoly(var)

    def minimal_polynomial(self, var='x'):
        """
        Return the minimal polynomial of this element.

        EXAMPLES::

            sage: GF(241, 'a')(1).minimal_polynomial(var = 'z')
            z + 240
        """
        return self.minpoly(var)

    def polynomial(self, var='x'):
        """
        Return a constant polynomial representing this value.

        EXAMPLES::

            sage: k = GF(7)
            sage: a = k.gen(); a
            1
            sage: a.polynomial()
            1
            sage: type(a.polynomial())                                                  # needs sage.libs.flint
            <class 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>
        """
        R = self.parent()[var]
        return R(self)

    def norm(self):
        """
        Return the norm of this element, which is itself. (This is here
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
        Return the trace of this element, which is itself. (This is here
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

    def lift_centered(self):
        r"""
        Lift ``self`` to a centered congruent integer.

        OUTPUT:

        The unique integer `i` such that `-n/2 < i \leq n/2` and `i = self \mod n`
        (where `n` denotes the modulus).

        EXAMPLES::

            sage: Mod(0,5).lift_centered()
            0
            sage: Mod(1,5).lift_centered()
            1
            sage: Mod(2,5).lift_centered()
            2
            sage: Mod(3,5).lift_centered()
            -2
            sage: Mod(4,5).lift_centered()
            -1
            sage: Mod(50,100).lift_centered()
            50
            sage: Mod(51,100).lift_centered()
            -49
            sage: Mod(-1,3^100).lift_centered()
            -1
        """
        n = self.modulus()
        x = self.lift()
        if 2*x <= n:
            return x
        else:
            return x - n

    cpdef bint is_one(self) noexcept:
        raise NotImplementedError

    cpdef bint is_unit(self) noexcept:
        raise NotImplementedError

    @coerce_binop
    def divides(self, other):
        r"""
        Test whether ``self`` divides ``other``.

        EXAMPLES::

            sage: R = Zmod(6)
            sage: R(2).divides(R(4))
            True
            sage: R(4).divides(R(2))
            True
            sage: R(2).divides(R(3))
            False
        """
        if not other:
            return True
        elif not self:
            return False
        mod = self.modulus()
        sl = self.lift().gcd(mod)
        if sl.is_one():
            return True
        return sl.divides(other.lift().gcd(mod))

    def is_square(self):
        r"""
        EXAMPLES::

            sage: Mod(3, 17).is_square()
            False

            sage: # needs sage.libs.pari
            sage: Mod(9, 17).is_square()
            True
            sage: Mod(9, 17*19^2).is_square()
            True
            sage: Mod(-1, 17^30).is_square()
            True
            sage: Mod(1/9, next_prime(2^40)).is_square()
            True
            sage: Mod(1/25, next_prime(2^90)).is_square()
            True

        TESTS::

            sage: Mod(1/25, 2^8).is_square()                                            # needs sage.libs.pari
            True
            sage: Mod(1/25, 2^40).is_square()                                           # needs sage.libs.pari
            True

            sage: for p,q,r in cartesian_product_iterator([[3,5],[11,13],[17,19]]):  # long time, needs sage.libs.pari
            ....:     for ep,eq,er in cartesian_product_iterator([[0,1,2,3],[0,1,2,3],[0,1,2,3]]):
            ....:         for e2 in [0, 1, 2, 3, 4]:
            ....:             n = p^ep * q^eq * r^er * 2^e2
            ....:             for _ in range(2):
            ....:                 a = Zmod(n).random_element()
            ....:                 if a.is_square().__xor__(a.__pari__().issquare()):
            ....:                     print(a, n)

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
        cdef int l2, m2
        if self.is_zero() or self.is_one():
            return 1
        # We first try to rule out self being a square without
        # factoring the modulus.
        lift = self.lift()
        m2, modd = self.modulus().val_unit(2)
        if m2 == 2:
            if lift & 2 == 2:  # lift = 2 or 3 (mod 4)
                return 0
        elif m2 > 2:
            l2, lodd = lift.val_unit(2)
            if l2 < m2 and (l2 % 2 == 1 or lodd % (1 << min(3, m2 - l2)) != 1):
                return 0
        # self is a square modulo 2^m2.  We compute the Jacobi symbol
        # modulo modd.  If this is -1, then self is not a square.
        if lift.jacobi(modd) == -1:
            return 0
        # We need to factor the modulus.  We do it here instead of
        # letting PARI do it, so that we can cache the factorisation.
        return lift.__pari__().Zn_issquare(self._parent.factored_order())

    def sqrt(self, extend=True, all=False):
        r"""
        Return square root or square roots of ``self`` modulo `n`.

        INPUT:

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the square root is not in the base ring.

        - ``all`` -- boolean (default: ``False``); if ``True``, return {all}
          square roots of self, instead of just one

        ALGORITHM: Calculates the square roots mod `p` for each of
        the primes `p` dividing the order of the ring, then lifts
        them `p`-adically and uses the CRT to find a square root
        mod `n`.

        See also :meth:`square_root_mod_prime_power` and
        :meth:`square_root_mod_prime` for more algorithmic details.

        EXAMPLES::

            sage: mod(-1, 17).sqrt()
            4
            sage: mod(5, 389).sqrt()
            86
            sage: mod(7, 18).sqrt()
            5

            sage: # needs sage.libs.pari
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

            sage: a = Mod(3, 5); a
            3
            sage: x = Mod(-1, 360)
            sage: x.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square
            sage: y = x.sqrt(); y
            sqrt359
            sage: y.parent()
            Univariate Quotient Polynomial Ring in sqrt359 over
             Ring of integers modulo 360 with modulus x^2 + 1
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

            sage: # needs sage.libs.pari
            sage: R = Integers(5*13^3*37); R
            Ring of integers modulo 406445
            sage: v = R(-1).sqrt(all=True); v
            [78853, 111808, 160142, 193097, 213348, 246303, 294637, 327592]
            sage: [x^2 for x in v]
            [406444, 406444, 406444, 406444, 406444, 406444, 406444, 406444]
            sage: v = R(169).sqrt(all=True); min(v), -max(v), len(v)
            (13, 13, 104)
            sage: all(x^2 == 169 for x in v)
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: t = FiniteField(next_prime(2^100))(4)
            sage: t.sqrt(extend=False, all=True)
            [2, 1267650600228229401496703205651]
            sage: t = FiniteField(next_prime(2^100))(2)
            sage: t.sqrt(extend=False, all=True)
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
                y = 'sqrt%s' % self
                R = self.parent()['x']
                modulus = R.gen()**2 - R(self)
                if self._parent.is_field():
                    from sage.rings.finite_rings.finite_field_constructor import FiniteField
                    Q = FiniteField(self._modulus.sageInteger**2, y, modulus)
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
            raise ValueError("self must be a square")

        F = self._parent.factored_order()
        cdef long e, exp, val
        if len(F) == 1:
            p, e = F[0]

            if all and e > 1 and not self.is_unit():
                if self.is_zero():
                    # All multiples of p^ciel(e/2) vanish
                    return [self._parent(x) for x in range(0, self._modulus.sageInteger, p**((e+1)/2))]
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
                            for b in range(0, self._modulus.sageInteger, p_exp)]
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
                from sage.arith.misc import CRT_basis
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
        Return an `n`-th root of ``self``.

        INPUT:

        - ``n`` -- integer `\geq 1`

        - ``extend`` -- boolean (default: ``True``); if ``True``, return an
          `n`-th root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the root is not in the base ring.  Warning:
          this option is not implemented!

        - ``all`` -- boolean (default: ``False``); if ``True``, return all
          `n`-th roots of ``self``, instead of just one

        - ``algorithm`` -- string (default: ``None``); the algorithm for the
          prime modulus case. CRT and `p`-adic log techniques are used to reduce
          to this case. ``'Johnston'`` is the only currently supported option.

        - ``cunningham`` -- boolean (default: ``False``); in some cases,
          factorization of `n` is computed. If cunningham is set to ``True``,
          the factorization of `n` is computed using trial division for all
          primes in the so called Cunningham table. Refer to
          ``sage.rings.factorint.factor_cunningham`` for more information. You
          need to install an optional package to use this method, this can be
          done with the following command line: ``sage -i cunningham_tables``.

        OUTPUT:

        If ``self`` has an `n`-th root, returns one (if ``all`` is ``False``) or a
        list of all of them (if ``all`` is ``True``).  Otherwise, raises a
        :exc:`ValueError` (if ``extend`` is ``False``) or a
        :exc:`NotImplementedError` (if ``extend`` is ``True``).

        .. warning::

           The 'extend' option is not implemented (yet).

        NOTE:

        - If `n = 0`:

          - if ``all=True``:

            - if ``self=1``: all nonzero elements of the parent are returned in
              a list.  Note that this could be very expensive for large
              parents.

            - otherwise: an empty list is returned

          - if ``all=False``:

            - if ``self=1``: ``self`` is returned

            - otherwise; a :exc:`ValueError` is raised

        - If `n < 0`:

          - if ``self`` is invertible, the `(-n)`\th root of the inverse of ``self`` is returned

          - otherwise a :exc:`ValueError` is raised or empty list returned.

        EXAMPLES::


            sage: K = GF(31)
            sage: a = K(22)
            sage: K(22).nth_root(7)
            13
            sage: K(25).nth_root(5)
            5
            sage: K(23).nth_root(3)
            29

            sage: # needs sage.rings.padics
            sage: mod(225, 2^5*3^2).nth_root(4, all=True)
            [225, 129, 33, 63, 255, 159, 9, 201, 105, 279, 183, 87, 81,
             273, 177, 207, 111, 15, 153, 57, 249, 135, 39, 231]
            sage: mod(275, 2^5*7^4).nth_root(7, all=True)
            [58235, 25307, 69211, 36283, 3355, 47259, 14331]
            sage: mod(1,8).nth_root(2, all=True)
            [1, 7, 5, 3]
            sage: mod(4,8).nth_root(2, all=True)
            [2, 6]
            sage: mod(1,16).nth_root(4, all=True)
            [1, 15, 13, 3, 9, 7, 5, 11]

            sage: (mod(22,31)^200).nth_root(200)
            5
            sage: mod(3,6).nth_root(0, all=True)
            []
            sage: mod(3,6).nth_root(0)
            Traceback (most recent call last):
            ...
            ValueError
            sage: mod(1,6).nth_root(0, all=True)
            [1, 2, 3, 4, 5]

        TESTS::

            sage: for p in [1009,2003,10007,100003]:                                    # needs sage.rings.finite_rings
            ....:     K = GF(p)
            ....:     for r in (p-1).divisors():
            ....:         if r == 1: continue
            ....:         x = K.random_element()
            ....:         y = x^r
            ....:         if y.nth_root(r)**r != y: raise RuntimeError
            ....:         if (y^41).nth_root(41*r)**(41*r) != y^41: raise RuntimeError
            ....:         if (y^307).nth_root(307*r)**(307*r) != y^307: raise RuntimeError

            sage: for t in range(200):                                                  # needs sage.libs.pari sage.rings.padics
            ....:     n = randint(1,2^63)
            ....:     K = Integers(n)
            ....:     b = K.random_element()
            ....:     e = randint(-2^62, 2^63)
            ....:     try:
            ....:         a = b.nth_root(e)
            ....:         if a^e != b:
            ....:             print(n, b, e, a)
            ....:             raise NotImplementedError
            ....:     except ValueError:
            ....:         pass

        We check that :issue:`13172` is resolved::

            sage: mod(-1, 4489).nth_root(2, all=True)                                   # needs sage.rings.padics
            []

        We check that :issue:`32084` is fixed::

            sage: mod(24, 25).nth_root(50)^50                                           # needs sage.rings.padics
            24

        Check that the code path cunningham might be used::

            sage: a = Mod(9,11)
            sage: a.nth_root(2, False, True, 'Johnston', cunningham=True)   # optional - cunningham_tables
            [3, 8]

        ALGORITHM:

        The default for prime modulus is currently an algorithm
        described in [Joh1999]_.

        AUTHORS:

        - David Roe (2010-02-13)
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
                    raise ValueError("no nth root")
            if pval > 0:
                if all:
                    return [K(a.lift()*p**(pval // n) + p**(k - (pval - pval//n)) * b) for a in mod(upart, p**(k-pval)).nth_root(n, all=True, algorithm=algorithm) for b in range(p**(pval - pval//n))]
                else:
                    return K(p**(pval // n) * mod(upart, p**(k-pval)).nth_root(n, algorithm=algorithm).lift())
            from sage.rings.padics.factory import ZpFM
            R = ZpFM(p,k)
            if p == 2:
                sign = [1]
                if self % 4 == 3:
                    if n % 2 == 0:
                        if all: return []
                        else: raise ValueError("no nth root")
                    else:
                        sign = [-1]
                        self = -self
                elif n % 2 == 0:
                    if k > 2 and self % 8 == 5:
                        if all: return []
                        else: raise ValueError("no nth root")
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
                        return [s*K(p*a+m.lift()) for a in range(p**(k-(2 if p==2 else 1))) for m in modp for s in sign]
                    else:
                        return K(modp.lift())
                else:
                    if all: return []
                    else: raise ValueError("no nth root")
            if all:
                ans = [plog // n + p**(k - nval) * i for i in range(p**nval)]
                ans = [s*K(R.teichmuller(m) * a.exp()) for a in ans for m in modp for s in sign]
                return ans
            else:
                return sign[0] * K(R.teichmuller(modp) * (plog // n).exp())
        return self._nth_root_common(n, all, algorithm, cunningham)

    def _nth_root_naive(self, n):
        """
        Compute all `n`-th roots using brute force, for doc-testing.

        TESTS::

            sage: for n in range(2,100):  # long time
            ....:     K = Integers(n)
            ....:     elist = list(range(1,min(2*n+2,100)))
            ....:     for e in random_sublist(elist, 5/len(elist)):
            ....:         for a in random_sublist(range(1,n), min((n+2)//2,10)/(n-1)):
            ....:             b = K(a)
            ....:             try:
            ....:                 L = b.nth_root(e, all=True)
            ....:                 if L:
            ....:                     c = b.nth_root(e)
            ....:             except Exception:
            ....:                 L = [-1]
            ....:             M = b._nth_root_naive(e)
            ....:             if sorted(L) != M:
            ....:                 print("mod(%s, %s).nth_root(%s,all=True), mod(%s, %s)._nth_root_naive(%s)" % (a,n,e,a,n,e))
            ....:             if L and (c not in L):
            ....:                 print("mod(%s, %s).nth_root(%s), mod(%s, %s).nth_root(%s,all=True)" % (a,n,e,a,n,e))
        """
        return [a for a in self.parent() if a**n == self]

    def _balanced_abs(self):
        r"""
        This function returns `x` or `-x`, whichever has a
        positive representative in `-n/2 < x \leq n/2`.

        This is used so that the same square root is always returned,
        despite the possibly probabilistic nature of the underlying
        algorithm.
        """
        if self.lift() > self._modulus.sageInteger >> 1:
            return -self
        return self

    def rational_reconstruction(self):
        """
        Use rational reconstruction to try to find a lift of this element to
        the rational numbers.

        EXAMPLES::

            sage: R = IntegerModRing(97)
            sage: a = R(2) / R(3)
            sage: a
            33
            sage: a.rational_reconstruction()
            2/3

        This method is also inherited by prime finite fields elements::

            sage: k = GF(97)
            sage: a = k(RationalField()('2/3'))
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

            sage: a = mod(3, 5)
            sage: b = mod(2, 7)
            sage: a.crt(b)
            23

        ::

            sage: a = mod(37, 10^8)
            sage: b = mod(9, 3^8)
            sage: a.crt(b)
            125900000037

        ::

            sage: b = mod(0, 1)
            sage: a.crt(b) == a
            True
            sage: a.crt(b).modulus()
            100000000

        TESTS::

            sage: mod(0, 1).crt(mod(4, 2^127))
            4
            sage: mod(4, 2^127).crt(mod(0, 1))
            4
            sage: mod(4, 2^30).crt(mod(0, 1))
            4
            sage: mod(0, 1).crt(mod(4, 2^30))
            4
            sage: mod(0, 1).crt(mod(4, 2^15))
            4
            sage: mod(4, 2^15).crt(mod(0, 1))
            4

        AUTHORS:

        - Robert Bradshaw
        """
        cdef int_fast64_t new_modulus
        if not isinstance(self, IntegerMod_gmp) and not isinstance(other, IntegerMod_gmp):

            if other._modulus.int64 == 1: return self
            new_modulus = self._modulus.int64 * other._modulus.int64
            if new_modulus < INTEGER_MOD_INT32_LIMIT:
                return self._crt(other)

            elif new_modulus < INTEGER_MOD_INT64_LIMIT:
                if not isinstance(self, IntegerMod_int64):
                    self = IntegerMod_int64(self._parent, self.lift())
                if not isinstance(other, IntegerMod_int64):
                    other = IntegerMod_int64(other._parent, other.lift())
                return self._crt(other)

        if not isinstance(self, IntegerMod_gmp):
            if self._modulus.int64 == 1: return other
            self = IntegerMod_gmp(self._parent, self.lift())

        if not isinstance(other, IntegerMod_gmp):
            if other._modulus.int64 == 1: return self
            other = IntegerMod_gmp(other._parent, other.lift())

        return self._crt(other)

    def additive_order(self):
        r"""
        Return the additive order of ``self``.

        This is the same as ``self.order()``.

        EXAMPLES::

            sage: Integers(20)(2).additive_order()
            10
            sage: Integers(20)(7).additive_order()
            20
            sage: Integers(90308402384902)(2).additive_order()
            45154201192451
        """
        n = self._modulus.sageInteger
        return sage.rings.integer.Integer(n // self.lift().gcd(n))

    def is_primitive_root(self):
        """
        Determines whether this element generates the group of units modulo n.

        This is only possible if the group of units is cyclic, which occurs if
        n is 2, 4, a power of an odd prime or twice a power of an odd prime.

        EXAMPLES::

            sage: mod(1, 2).is_primitive_root()
            True
            sage: mod(3, 4).is_primitive_root()
            True
            sage: mod(2, 7).is_primitive_root()
            False
            sage: mod(3, 98).is_primitive_root()                                        # needs sage.libs.pari
            True
            sage: mod(11, 1009^2).is_primitive_root()                                   # needs sage.libs.pari
            True

        TESTS::

            sage: for p in prime_range(3,12):                                           # needs sage.libs.pari
            ....:     for k in range(1,4):
            ....:         for even in [1,2]:
            ....:             n = even*p^k
            ....:             phin = euler_phi(n)
            ....:             for _ in range(6):
            ....:                 a = Zmod(n).random_element()
            ....:                 if not a.is_unit(): continue
            ....:                 if a.is_primitive_root().__xor__(a.multiplicative_order()==phin):
            ....:                     print("mod(%s,%s) incorrect" % (a, n))

        `0` is not a primitive root mod `n` (:issue:`23624`) except for `n=0`::

            sage: mod(0, 17).is_primitive_root()
            False
            sage: all(not mod(0, n).is_primitive_root() for n in srange(2, 20))         # needs sage.libs.pari
            True
            sage: mod(0, 1).is_primitive_root()
            True

            sage: all(not mod(p^j, p^k).is_primitive_root()                             # needs sage.libs.pari
            ....:     for p in prime_range(3, 12)
            ....:     for k in srange(1, 4)
            ....:     for j in srange(0, k))
            True
        """
        cdef Integer p1, q = Integer(2)
        m = self.modulus()
        if m == 1:
            return True
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
        if self == 0:
            return False
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
        Return the multiplicative order of ``self``.

        EXAMPLES::

            sage: Mod(-1, 5).multiplicative_order()                                     # needs sage.libs.pari
            2
            sage: Mod(1, 5).multiplicative_order()                                      # needs sage.libs.pari
            1
            sage: Mod(0, 5).multiplicative_order()                                      # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ArithmeticError: multiplicative order of 0 not defined
            since it is not a unit modulo 5
        """
        try:
            return sage.rings.integer.Integer(self.__pari__().znorder())
        except PariError:
            raise ArithmeticError("multiplicative order of %s not defined since it is not a unit modulo %s"%(
                self, self._modulus.sageInteger))

    def valuation(self, p):
        """
        The largest power `r` such that `m` is in the ideal generated by `p^r` or infinity if there is not a largest such power.
        However it is an error to take the valuation with respect to a unit.

        .. NOTE::

            This is not a valuation in the mathematical sense. As shown with the examples below.

        EXAMPLES:

        This example shows that ``(a*b).valuation(n)`` is not always the same as ``a.valuation(n) + b.valuation(n)``

        ::

            sage: R = ZZ.quo(9)
            sage: a = R(3)
            sage: b = R(6)
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

            sage: R = ZZ.quo(12)
            sage: a = R(2)
            sage: b = R(4)
            sage: a.valuation(2)
            1
            sage: b.valuation(2)
            +Infinity
            sage: ZZ.quo(1024)(16).valuation(4)
            2
        """
        p=self._modulus.sageInteger.gcd(p)
        if p==1:
            raise ValueError("Valuation with respect to a unit is not defined.")
        r = 0
        power = p
        while not (self % power): # self % power == 0
            r += 1
            power *= p
            if not power.divides(self._modulus.sageInteger):
                from sage.rings.infinity import infinity
                return infinity
        return r

    cpdef _floordiv_(self, right):
        """
        Exact division for prime moduli, for compatibility with other fields.

        EXAMPLES::

            sage: GF(7)(3) // 5
            2
        """
        return self._mul_(~right)

    def _repr_(self):
        return str(self.lift())

    def _latex_(self):
        return str(self)

    def _integer_(self, ZZ=None):
        return self.lift()

    def _rational_(self):
        return rational.Rational(self.lift())

    def _vector_(self):
        """
        Return ``self`` as a vector of its parent viewed as a one-dimensional
        vector space.

        This is to support prime finite fields, which are implemented as
        `IntegerMod` ring.

        EXAMPLES::

            sage: F.<a> = GF(13)
            sage: V = F.vector_space(map=False)                                         # needs sage.modules
            sage: V(a)                                                                  # needs sage.modules
            (1)
        """
        return self.parent().vector_space(map=False)([self])


######################################################################
#      class IntegerMod_gmp
######################################################################


cdef class IntegerMod_gmp(IntegerMod_abstract):
    r"""
    Elements of `\ZZ/n\ZZ` for n not small enough
    to be operated on in word size.

    AUTHORS:

    - Robert Bradshaw (2006-08-24)
    """

    def __cinit__(self):
        mpz_init(self.value)

    cdef IntegerMod_gmp _new_c(self):
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp.__new__(IntegerMod_gmp)
        x._modulus = self._modulus
        x._parent = self._parent
        return x

    def __dealloc__(self):
        mpz_clear(self.value)

    cdef void set_from_mpz(self, mpz_t value) noexcept:
        cdef sage.rings.integer.Integer modulus
        modulus = self._modulus.sageInteger
        mpz_mod(self.value, value, modulus.value)

    cdef void set_from_long(self, long value) noexcept:
        r"""
        EXAMPLES::

            sage: p = next_prime(2^32)                                                  # needs sage.libs.pari
            sage: GF(p)(int(p + 1))                                                     # needs sage.libs.pari sage.rings.finite_rings
            1
        """
        mpz_set_si(self.value, value)
        mpz_mod(self.value, self.value, self._modulus.sageInteger.value)

    cdef void set_from_ulong_fast(self, unsigned long value) noexcept:
        mpz_set_ui(self.value, value)

    def __lshift__(IntegerMod_gmp self, k):
        r"""
        Perform a left shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(19, 10^10)
            sage: e << 102
            9443608576
        """
        return self.shift(long(k))

    def __rshift__(IntegerMod_gmp self, k):
        r"""
        Perform a right shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(19, 10^10)
            sage: e >> 1
            9
        """
        return self.shift(-long(k))

    cdef shift(IntegerMod_gmp self, long k):
        r"""
        Perform a bit-shift specified by ``k`` on ``self``.

        Suppose that ``self`` represents an integer `x` modulo `n`.  If `k` is
        `k = 0`, returns `x`.  If `k > 0`, shifts `x` to the left, that is,
        multiplies `x` by `2^k` and then returns the representative in the
        range `[0,n)`.  If `k < 0`, shifts `x` to the right, that is, returns
        the integral part of `x` divided by `2^k`.

        Note that, in any case, ``self`` remains unchanged.

        INPUT:

        - ``k`` -- integer of type ``long``

        OUTPUT: result of type ``IntegerMod_gmp``

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
                mpz_fdiv_r(x.value, x.value, self._modulus.sageInteger.value)
            else:
                mpz_fdiv_q_2exp(x.value, self.value, -k)
            return x

    cpdef _richcmp_(left, right, int op):
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
        return rich_to_bool_sgn(op, i)

    cpdef bint is_one(IntegerMod_gmp self) noexcept:
        """
        Return ``True`` if this is `1`, otherwise ``False``.

        EXAMPLES::

            sage: mod(1,5^23).is_one()
            True
            sage: mod(0,5^23).is_one()
            False
        """
        return mpz_cmp_si(self.value, 1) == 0

    def __bool__(IntegerMod_gmp self):
        """
        Return ``True`` if this is not `0`, otherwise ``False``.

        EXAMPLES::

            sage: mod(13,5^23).is_zero()
            False
            sage: (mod(25,5^23)^23).is_zero()
            True
        """
        return mpz_cmp_si(self.value, 0) != 0

    cpdef bint is_unit(self) noexcept:
        """
        Return ``True`` iff this element is a unit.

        EXAMPLES::

            sage: mod(13, 5^23).is_unit()
            True
            sage: mod(25, 5^23).is_unit()
            False
        """
        return self.lift().gcd(self.modulus()) == 1

    def _crt(IntegerMod_gmp self, IntegerMod_gmp other):
        cdef IntegerMod_gmp lift, x
        cdef sage.rings.integer.Integer modulus, other_modulus

        modulus = self._modulus.sageInteger
        other_modulus = other._modulus.sageInteger
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
        lift = IntegerMod_gmp(IntegerModRing(modulus*other_modulus))
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
            raise ZeroDivisionError("moduli must be coprime")

    def __copy__(IntegerMod_gmp self):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R7 = R(7)
            sage: copy(R7) is R7
            True
        """
        # immutable
        return self

    def __deepcopy__(IntegerMod_gmp self, memo):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R7 = R(7)
            sage: deepcopy(R7) is R7
            True
        """
        # immutable
        return self

    cpdef _add_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R(7) + R(8)
            15
        """
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_add(x.value, self.value, (<IntegerMod_gmp>right).value)
        if mpz_cmp(x.value, self._modulus.sageInteger.value)  >= 0:
            mpz_sub(x.value, x.value, self._modulus.sageInteger.value)
        return x

    cpdef _sub_(self, right):
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
            mpz_add(x.value, x.value, self._modulus.sageInteger.value)
        return x

    cpdef _neg_(self):
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
        mpz_sub(x.value, self._modulus.sageInteger.value, self.value)
        return x

    cpdef _mul_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10^11)
            sage: R(700000) * R(800000)
            60000000000
        """
        cdef IntegerMod_gmp x
        x = self._new_c()
        mpz_mul(x.value, self.value,  (<IntegerMod_gmp>right).value)
        mpz_fdiv_r(x.value, x.value, self._modulus.sageInteger.value)
        return x

    cpdef _div_(self, right):
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

    def __pow__(IntegerMod_gmp self, exp, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES::

            sage: R = Integers(10^10)
            sage: R(2)^1000
            5668069376
            sage: p = next_prime(11^10)                                                 # needs sage.libs.pari
            sage: R = Integers(p)                                                       # needs sage.libs.pari
            sage: R(9876)^(p-1)                                                         # needs sage.libs.pari
            1
            sage: mod(3, 10^100)^-2
            8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889
            sage: mod(2, 10^100)^-2
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        TESTS:

        We define ``0^0`` to be unity, :issue:`13894`::

            sage: p = next_prime(11^10)                                                 # needs sage.libs.pari
            sage: R = Integers(p)                                                       # needs sage.libs.pari
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: type(R(0)^0) == type(R(0))
            True

        When the modulus is ``1``, the only element in the ring is
        ``0`` (and it is equivalent to ``1``), so we return that
        instead::

            sage: from sage.rings.finite_rings.integer_mod \
            ....:     import IntegerMod_gmp
            sage: zero = IntegerMod_gmp(Integers(1),0)
            sage: type(zero)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: zero^0
            0
        """
        cdef IntegerMod_gmp x = self._new_c()
        sig_on()
        try:
            mpz_pow_helper(x.value, self.value, exp, self._modulus.sageInteger.value)
        finally:
            sig_off()
        return x

    def __invert__(IntegerMod_gmp self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: a = mod(3,10^100); type(a)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: ~a
            6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667
            sage: ~mod(2,10^100)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(2, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000) does not exist
        """
        if self.is_zero():
            raise ZeroDivisionError(f"inverse of Mod(0, {self._modulus.sageInteger}) does not exist")

        cdef IntegerMod_gmp x
        x = self._new_c()
        if not mpz_invert(x.value, self.value, self._modulus.sageInteger.value):
            raise ZeroDivisionError(f"inverse of Mod({self}, {self._modulus.sageInteger}) does not exist")
        return x

    def lift(IntegerMod_gmp self):
        """
        Lift an integer modulo `n` to the integers.

        EXAMPLES::

            sage: a = Mod(8943, 2^70); type(a)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
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
        return mpz_pythonhash(self.value)

    @coerce_binop
    def gcd(self, IntegerMod_gmp other):
        r"""
        Greatest common divisor.

        Returns the "smallest" generator in `\ZZ / N\ZZ` of the ideal
        generated by ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of the same ring as this one

        EXAMPLES::

            sage: mod(2^3*3^2*5, 3^3*2^2*17^8).gcd(mod(2^4*3*17, 3^3*2^2*17^8))
            12
            sage: mod(0,17^8).gcd(mod(0,17^8))
            0
        """
        cdef IntegerMod_gmp ans = self._new_c()
        sig_on()
        mpz_gcd(ans.value, self.value, self._modulus.sageInteger.value)
        mpz_gcd(ans.value, ans.value, other.value)
        sig_off()
        if mpz_cmp(ans.value, self._modulus.sageInteger.value) == 0:
            # self = other = 0
            mpz_set_ui(ans.value, 0)
        return ans

######################################################################
#      class IntegerMod_int
######################################################################


cdef class IntegerMod_int(IntegerMod_abstract):
    r"""
    Elements of `\ZZ/n\ZZ` for n small enough to
    be operated on in 32 bits

    AUTHORS:

    - Robert Bradshaw (2006-08-24)

    EXAMPLES::

        sage: a = Mod(10,30); a
        10
        sage: loads(a.dumps()) == a
        True
    """

    cdef IntegerMod_int _new_c(self, int_fast32_t value):
        if self._modulus.table is not None:
            return self._modulus.table[value]
        cdef IntegerMod_int x = IntegerMod_int.__new__(IntegerMod_int)
        x._parent = self._parent
        x._modulus = self._modulus
        x.ivalue = value
        return x

    cdef void set_from_mpz(self, mpz_t value) noexcept:
        self.ivalue = mpz_fdiv_ui(value, self._modulus.int32)

    cdef void set_from_long(self, long value) noexcept:
        self.ivalue = value % self._modulus.int32
        if self.ivalue < 0:
            self.ivalue += self._modulus.int32

    cdef void set_from_ulong_fast(self, unsigned long value) noexcept:
        self.ivalue = value

    cdef void set_from_int(IntegerMod_int self, int_fast32_t ivalue) noexcept:
        if ivalue < 0:
            self.ivalue = self._modulus.int32 + (ivalue % self._modulus.int32)
        elif ivalue >= self._modulus.int32:
            self.ivalue = ivalue % self._modulus.int32
        else:
            self.ivalue = ivalue

    cdef int_fast32_t get_int_value(IntegerMod_int self) noexcept:
        return self.ivalue

    cpdef _richcmp_(self, right, int op):
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
            return rich_to_bool(op, 0)
        elif self.ivalue < (<IntegerMod_int>right).ivalue:
            return rich_to_bool(op, -1)
        else:
            return rich_to_bool(op, 1)

    cpdef bint is_one(IntegerMod_int self) noexcept:
        """
        Return ``True`` if this is `1`, otherwise ``False``.

        EXAMPLES::

            sage: mod(6,5).is_one()
            True
            sage: mod(0,5).is_one()
            False
            sage: mod(1, 1).is_one()
            True
            sage: Zmod(1).one().is_one()
            True
        """
        return self.ivalue == 1 or self._modulus.int32 == 1

    def __bool__(IntegerMod_int self):
        """
        Return ``True`` if this is not `0`, otherwise ``False``.

        EXAMPLES::

            sage: mod(13,5).is_zero()
            False
            sage: mod(25,5).is_zero()
            True
        """
        return self.ivalue != 0

    cpdef bint is_unit(IntegerMod_int self) noexcept:
        """
        Return ``True`` iff this element is a unit

        EXAMPLES::

            sage: a=Mod(23,100)
            sage: a.is_unit()
            True
            sage: a=Mod(24,100)
            sage: a.is_unit()
            False
        """
        return gcd_int(self.ivalue, self._modulus.int32) == 1

    def _crt(IntegerMod_int self, IntegerMod_int other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to ``self`` and
        to ``other``. The modulus of ``other`` must be coprime to the modulus
        of ``self``.

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

        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
        lift = IntegerMod_int(IntegerModRing(self._modulus.int32 * other._modulus.int32))

        try:
            x = (other.ivalue - self.ivalue % other._modulus.int32) * mod_inverse_int(self._modulus.int32, other._modulus.int32)
            lift.set_from_int( x * self._modulus.int32 + self.ivalue )
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError("moduli must be coprime")

    def __copy__(IntegerMod_int self):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R7 = R(7)
            sage: copy(R7) is R7
            True
        """
        # immutable
        return self

    def __deepcopy__(IntegerMod_int self, memo):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R7 = R(7)
            sage: deepcopy(R7) is R7
            True
        """
        # immutable
        return self

    cpdef _add_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) + R(8)
            5
        """
        cdef int_fast32_t x
        x = self.ivalue + (<IntegerMod_int>right).ivalue
        if x >= self._modulus.int32:
            x = x - self._modulus.int32
        return self._new_c(x)

    cpdef _sub_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef int_fast32_t x
        x = self.ivalue - (<IntegerMod_int>right).ivalue
        if x < 0:
            x = x + self._modulus.int32
        return self._new_c(x)

    cpdef _neg_(self):
        """
        EXAMPLES::

            sage: -mod(7,10)
            3
            sage: -mod(0,10)
            0
        """
        if self.ivalue == 0:
            return self
        return self._new_c(self._modulus.int32 - self.ivalue)

    cpdef _mul_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        return self._new_c((self.ivalue * (<IntegerMod_int>right).ivalue) % self._modulus.int32)

    cpdef _div_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(2)/3
            4
        """
        if self._modulus.inverses is not None:
            right_inverse = self._modulus.inverses[(<IntegerMod_int>right).ivalue]
            if right_inverse is None:
                raise ZeroDivisionError(f"inverse of Mod({right}, {self._modulus.sageInteger}) does not exist")
            else:
                return self._new_c((self.ivalue * (<IntegerMod_int>right_inverse).ivalue) % self._modulus.int32)

        cdef int_fast32_t x
        x = self.ivalue * mod_inverse_int((<IntegerMod_int>right).ivalue, self._modulus.int32)
        return self._new_c(x% self._modulus.int32)

    def __int__(IntegerMod_int self):
        """
        TESTS::

            sage: e = Mod(8, 31)
            sage: int(e)
            8
        """
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

    def __lshift__(IntegerMod_int self, k):
        r"""
        Perform a left shift by ``k`` bits.

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
        Perform a right shift by ``k`` bits.

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
        Perform a bit-shift specified by ``k`` on ``self``.

        Suppose that ``self`` represents an integer `x` modulo `n`.  If `k` is
        `k = 0`, returns `x`.  If `k > 0`, shifts `x` to the left, that is,
        multiplies `x` by `2^k` and then returns the representative in the
        range `[0,n)`.  If `k < 0`, shifts `x` to the right, that is, returns
        the integral part of `x` divided by `2^k`.

        Note that, in any case, ``self`` remains unchanged.

        INPUT:

        - ``k`` -- integer of type ``int``

        OUTPUT: result of type ``IntegerMod_int``

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
        """
        if k == 0:
            return self
        elif k > 0:
            return self._new_c((self.ivalue << k) % self._modulus.int32)
        else:
            return self._new_c(self.ivalue >> (-k))

    def __pow__(IntegerMod_int self, exp, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES::

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
            ZeroDivisionError: inverse of Mod(2, 100) does not exist
            sage: mod(2, 100)^-100000000
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        TESTS:

        We define ``0^0`` to be unity, :issue:`13894`::

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
        if type(exp) is int and -100000 < PyLong_AsLong(exp) < 100000:
            long_exp = PyLong_AsLong(exp)
        elif type(exp) is Integer and mpz_cmpabs_ui((<Integer>exp).value, 100000) == -1:
            long_exp = mpz_get_si((<Integer>exp).value)
        else:
            base = self.lift()
            sig_on()
            try:
                mpz_init(res_mpz)
                mpz_pow_helper(res_mpz, (<Integer>base).value, exp, self._modulus.sageInteger.value)
                res = mpz_get_ui(res_mpz)
                mpz_clear(res_mpz)
            finally:
                sig_off()
            return self._new_c(res)

        if long_exp == 0 and self.ivalue == 0:
            # Return 0 if the modulus is 1, otherwise return 1.
            return self._new_c(self._modulus.int32 != 1)
        cdef bint invert = False
        if long_exp < 0:
            invert = True
            long_exp = -long_exp
        res = mod_pow_int(self.ivalue, long_exp, self._modulus.int32)
        if invert:
            return ~self._new_c(res)
        else:
            return self._new_c(res)

    def __invert__(IntegerMod_int self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: ~mod(7,100)
            43
            sage: Mod(0,1)^-1
            0
        """
        if self._modulus.inverses is not None:
            x = self._modulus.inverses[self.ivalue]
            if x is None:
                raise ZeroDivisionError(f"inverse of Mod({self}, {self._modulus.sageInteger}) does not exist")
            else:
                return x
        else:
            return self._new_c(mod_inverse_int(self.ivalue, self._modulus.int32))

    def lift(IntegerMod_int self):
        """
        Lift an integer modulo `n` to the integers.

        EXAMPLES::

            sage: a = Mod(8943, 2^10); type(a)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
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
        cdef int_fast32_t l2, lodd, m2, modd
        if self.ivalue <= 1:
            return 1
        # We first try to rule out self being a square without
        # factoring the modulus.
        lift = self.lift()
        m2, modd = self.modulus().val_unit(2)
        if m2 == 2:
            if self.ivalue & 2 == 2:  # self.ivalue = 2 or 3 (mod 4)
                return 0
        elif m2 > 2:
            l2, lodd = lift.val_unit(2)
            if l2 < m2 and (l2 % 2 == 1 or lodd % (1 << min(3, m2 - l2)) != 1):
                return 0
        # self is a square modulo 2^m2.  We compute the Jacobi symbol
        # modulo modd.  If this is -1, then self is not a square.
        if jacobi_int(self.ivalue, modd) == -1:
            return 0
        # We need to factor the modulus.  We do it here instead of
        # letting PARI do it, so that we can cache the factorisation.
        return lift.__pari__().Zn_issquare(self._parent.factored_order())

    def sqrt(self, extend=True, all=False):
        r"""
        Return square root or square roots of ``self`` modulo `n`.

        INPUT:

        - ``extend`` -- boolean (default: ``True``);
          if ``True``, return a square root in an extension ring,
          if necessary. Otherwise, raise a :exc:`ValueError` if the
          square root is not in the base ring.

        - ``all`` -- boolean (default: ``False``); if
          ``True``, return {all} square roots of self, instead of
          just one.

        ALGORITHM: Calculates the square roots mod `p` for each of
        the primes `p` dividing the order of the ring, then lifts
        them `p`-adically and uses the CRT to find a square root
        mod `n`.

        See also :meth:`square_root_mod_prime_power` and
        :meth:`square_root_mod_prime` for more algorithmic details.

        EXAMPLES::

            sage: mod(-1, 17).sqrt()
            4
            sage: mod(5, 389).sqrt()
            86
            sage: mod(7, 18).sqrt()
            5

            sage: # needs sage.libs.pari
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
            Univariate Quotient Polynomial Ring in sqrt359
             over Ring of integers modulo 360 with modulus x^2 + 1
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

            sage: # needs sage.libs.pari
            sage: R = Integers(5*13^3*37); R
            Ring of integers modulo 406445
            sage: v = R(-1).sqrt(all=True); v
            [78853, 111808, 160142, 193097, 213348, 246303, 294637, 327592]
            sage: [x^2 for x in v]
            [406444, 406444, 406444, 406444, 406444, 406444, 406444, 406444]
            sage: v = R(169).sqrt(all=True); min(v), -max(v), len(v)
            (13, 13, 104)
            sage: all(x^2 == 169 for x in v)
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

        TESTS:

        Check for :issue:`30797`::

            sage: GF(103)(-1).sqrt(extend=False, all=True)
            []
        """
        cdef int_fast32_t i, n = self._modulus.int32
        if n > 100:
            moduli = self._parent.factored_order()
        # Unless the modulus is tiny, test to see if we're in the really
        # easy case of n prime, n = 3 mod 4.
        if n > 100 and n % 4 == 3 and len(moduli) == 1 and moduli[0][1] == 1:
            if jacobi_int(self.ivalue, self._modulus.int32) == 1:
                # it's a nonzero square, sqrt(a) = a^(p+1)/4
                i = mod_pow_int(self.ivalue, (self._modulus.int32+1)/4, n)
                if i > n / 2:
                    i = n - i
                if all:
                    return [self._new_c(i), self._new_c(n-i)]
                else:
                    return self._new_c(i)
            elif self.ivalue == 0:
                return [self] if all else self
            elif not extend:
                if all:
                    return []
                raise ValueError("self must be a square")
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
                    raise ValueError("self must be a square")
        # Either it failed but extend was True, or the generic algorithm is better
        return IntegerMod_abstract.sqrt(self, extend=extend, all=all)

    def _balanced_abs(self):
        r"""
        This function returns `x` or `-x`, whichever has a
        positive representative in `-n/2 < x \leq n/2`.
        """
        if self.ivalue > self._modulus.int32 / 2:
            return -self
        return self

    @coerce_binop
    def gcd(self, IntegerMod_int other):
        r"""
        Greatest common divisor.

        Returns the "smallest" generator in `\ZZ / N\ZZ` of the ideal
        generated by ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of the same ring as this one

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
        cdef int_fast32_t g = gcd_int(self.ivalue, self._modulus.int32)
        g = gcd_int(g, other.ivalue)
        if g == self._modulus.int32: # self = other = 0
            g = 0
        return self._new_c(g)

### End of class


cdef int_fast32_t gcd_int(int_fast32_t a, int_fast32_t b) noexcept:
    """
    Return the gcd of ``a`` and ``b``.

    For use with ``IntegerMod_int``.

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
    Return y such that xy=1 mod n.

    For use in ``IntegerMod_int``.

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
    raise ZeroDivisionError(f"inverse of Mod({x}, {n}) does not exist")


cdef int_fast32_t mod_pow_int(int_fast32_t base, int_fast32_t exp, int_fast32_t n) noexcept:
    """
    Return base^exp mod n.

    For use in ``IntegerMod_int``.

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
    Calculate the jacobi symbol (a/n).

    For use in ``IntegerMod_int``.

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int s, jacobi = 1
    cdef int_fast32_t b

    a = a % m

    while True:
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
    r"""
    Elements of `\ZZ/n\ZZ` for n small enough to
    be operated on in 64 bits

    EXAMPLES::

        sage: a = Mod(10,3^10); a
        10
        sage: type(a)
        <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
        sage: loads(a.dumps()) == a
        True
        sage: Mod(5, 2^31)
        5

    AUTHORS:

    - Robert Bradshaw (2006-09-14)
    """

    cdef IntegerMod_int64 _new_c(self, int_fast64_t value):
        cdef IntegerMod_int64 x
        x = IntegerMod_int64.__new__(IntegerMod_int64)
        x._modulus = self._modulus
        x._parent = self._parent
        x.ivalue = value
        return x

    cdef void set_from_mpz(self, mpz_t value) noexcept:
        self.ivalue = mpz_fdiv_ui(value, self._modulus.int64)

    cdef void set_from_long(self, long value) noexcept:
        self.ivalue = value % self._modulus.int64
        if self.ivalue < 0:
            self.ivalue += self._modulus.int64

    cdef void set_from_ulong_fast(self, unsigned long value) noexcept:
        self.ivalue = value

    cdef void set_from_int(IntegerMod_int64 self, int_fast64_t ivalue) noexcept:
        if ivalue < 0:
            self.ivalue = self._modulus.int64 + (ivalue % self._modulus.int64) # Is ivalue % self._modulus.int64 actually negative?
        elif ivalue >= self._modulus.int64:
            self.ivalue = ivalue % self._modulus.int64
        else:
            self.ivalue = ivalue

    cdef int_fast64_t get_int_value(IntegerMod_int64 self) noexcept:
        return self.ivalue

    cpdef _richcmp_(self, right, int op):
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
        if self.ivalue == (<IntegerMod_int64>right).ivalue:
            return rich_to_bool(op, 0)
        elif self.ivalue < (<IntegerMod_int64>right).ivalue:
            return rich_to_bool(op, -1)
        else:
            return rich_to_bool(op, 1)

    cpdef bint is_one(IntegerMod_int64 self) noexcept:
        """
        Return ``True`` if this is `1`, otherwise ``False``.

        EXAMPLES::

            sage: (mod(-1,5^10)^2).is_one()
            True
            sage: mod(0,5^10).is_one()
            False
        """
        return self.ivalue == 1

    def __bool__(IntegerMod_int64 self):
        """
        Return ``True`` if this is not `0`, otherwise ``False``.

        EXAMPLES::

            sage: mod(13,5^10).is_zero()
            False
            sage: mod(5^12,5^10).is_zero()
            True
        """
        return self.ivalue != 0

    cpdef bint is_unit(IntegerMod_int64 self) noexcept:
        """
        Return ``True`` iff this element is a unit.

        EXAMPLES::

            sage: mod(13, 5^10).is_unit()
            True
            sage: mod(25, 5^10).is_unit()
            False
        """
        return gcd_int64(self.ivalue, self._modulus.int64) == 1

    def _crt(IntegerMod_int64 self, IntegerMod_int64 other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to ``self`` and
        to ``other``. The modulus of ``other`` must be coprime to the modulus
        of ``self``.

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

        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
        lift = IntegerMod_int64(IntegerModRing(self._modulus.int64 * other._modulus.int64))

        try:
            x = (other.ivalue - self.ivalue % other._modulus.int64) * mod_inverse_int64(self._modulus.int64, other._modulus.int64)
            lift.set_from_int( x * self._modulus.int64 + self.ivalue )
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError("moduli must be coprime")

    def __copy__(IntegerMod_int64 self):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R7 = R(7)
            sage: copy(R7) is R7
            True
        """
        # immutable
        return self

    def __deepcopy__(IntegerMod_int64 self, memo):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R7 = R(7)
            sage: deepcopy(R7) is R7
            True
        """
        # immutable
        return self

    cpdef _add_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(7) + R(8)
            15
        """
        cdef int_fast64_t x
        x = self.ivalue + (<IntegerMod_int64>right).ivalue
        if x >= self._modulus.int64:
            x = x - self._modulus.int64
        return self._new_c(x)

    cpdef _sub_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(7) - R(8)
            99999
        """
        cdef int_fast64_t x
        x = self.ivalue - (<IntegerMod_int64>right).ivalue
        if x < 0:
            x = x + self._modulus.int64
        return self._new_c(x)

    cpdef _neg_(self):
        """
        EXAMPLES::

            sage: -mod(7,10^5)
            99993
            sage: -mod(0,10^6)
            0
        """
        if self.ivalue == 0:
            return self
        return self._new_c(self._modulus.int64 - self.ivalue)

    cpdef _mul_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(700) * R(800)
            60000
        """
        return self._new_c((self.ivalue * (<IntegerMod_int64>right).ivalue) % self._modulus.int64)

    cpdef _div_(self, right):
        """
        EXAMPLES::

            sage: R = Integers(10^5)
            sage: R(2)/3
            33334
        """
        return self._new_c((self.ivalue * mod_inverse_int64((<IntegerMod_int64>right).ivalue,
                                   self._modulus.int64) ) % self._modulus.int64)

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

    def __lshift__(IntegerMod_int64 self, k):
        r"""
        Perform a left shift by ``k`` bits.

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
        Perform a right shift by ``k`` bits.

        For details, see :meth:`shift`.

        EXAMPLES::

            sage: e = Mod(5, 2^31 - 1)
            sage: e >> 1
            2
        """
        return self.shift(-int(k))

    cdef shift(IntegerMod_int64 self, int k):
        """
        Perform a bit-shift specified by ``k`` on ``self``.

        Suppose that ``self`` represents an integer `x` modulo `n`.  If `k` is
        `k = 0`, returns `x`.  If `k > 0`, shifts `x` to the left, that is,
        multiplies `x` by `2^k` and then returns the representative in the
        range `[0,n)`.  If `k < 0`, shifts `x` to the right, that is, returns
        the integral part of `x` divided by `2^k`.

        Note that, in any case, ``self`` remains unchanged.

        INPUT:

        - ``k`` -- integer of type ``int``

        OUTPUT: result of type ``IntegerMod_int64``

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
            return self._new_c((self.ivalue << k) % self._modulus.int64)
        else:
            return self._new_c(self.ivalue >> (-k))

    def __pow__(IntegerMod_int64 self, exp, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: R(2)^10
            4
            sage: p = next_prime(10^5)                                                  # needs sage.libs.pari
            sage: R = Integers(p)                                                       # needs sage.libs.pari
            sage: R(1234)^(p - 1)                                                       # needs sage.libs.pari
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
            ZeroDivisionError: inverse of Mod(17, 1419857) does not exist
            sage: R(17)^-100000000
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        TESTS::

            sage: type(R(0))
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>

        We define ``0^0`` to be unity, :issue:`13894`::

            sage: p = next_prime(10^5)                                                  # needs sage.libs.pari
            sage: R = Integers(p)                                                       # needs sage.libs.pari
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: type(R(0)^0) == type(R(0))
            True

        When the modulus is ``1``, the only element in the ring is
        ``0`` (and it is equivalent to ``1``), so we return that
        instead::

            sage: from sage.rings.finite_rings.integer_mod \
            ....:     import IntegerMod_int64
            sage: zero = IntegerMod_int64(Integers(1),0)
            sage: type(zero)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
            sage: zero^0
            0
        """
        cdef long long_exp
        cdef int_fast64_t res
        cdef mpz_t res_mpz
        if type(exp) is int and -100000 < PyLong_AsLong(exp) < 100000:
            long_exp = PyLong_AsLong(exp)
        elif type(exp) is Integer and mpz_cmpabs_ui((<Integer>exp).value, 100000) == -1:
            long_exp = mpz_get_si((<Integer>exp).value)
        else:
            base = self.lift()
            sig_on()
            try:
                mpz_init(res_mpz)
                mpz_pow_helper(res_mpz, (<Integer>base).value, exp, self._modulus.sageInteger.value)
                res = mpz_get_ui(res_mpz)
                mpz_clear(res_mpz)
            finally:
                sig_off()
            return self._new_c(res)

        if long_exp == 0 and self.ivalue == 0:
            # Return 0 if the modulus is 1, otherwise return 1.
            return self._new_c(self._modulus.int64 != 1)
        cdef bint invert = False
        if long_exp < 0:
            invert = True
            long_exp = -long_exp
        res = mod_pow_int64(self.ivalue, long_exp, self._modulus.int64)
        if invert:
            return self._new_c(mod_inverse_int64(res, self._modulus.int64))
        else:
            return self._new_c(res)

    def __invert__(IntegerMod_int64 self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: a = mod(7,2^40); type(a)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
            sage: ~a
            471219269047
            sage: a
            7
        """
        return self._new_c(mod_inverse_int64(self.ivalue, self._modulus.int64))

    def lift(IntegerMod_int64 self):
        """
        Lift an integer modulo `n` to the integers.

        EXAMPLES::

            sage: a = Mod(8943, 2^25); type(a)
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
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
        Coerce ``self`` to a float.

        EXAMPLES::

            sage: a = Mod(8943, 2^35)
            sage: float(a)
            8943.0
        """
        return <double>self.ivalue

    def __hash__(self):
        """
        Compute hash of ``self``.

        EXAMPLES::

            sage: a = Mod(8943, 2^35)
            sage: hash(a)
            8943
        """
        return hash(self.ivalue)

    def _balanced_abs(self):
        r"""
        This function returns `x` or `-x`, whichever has a
        positive representative in `-n/2 < x \leq n/2`.
        """
        if self.ivalue > self._modulus.int64 / 2:
            return -self
        return self

    @coerce_binop
    def gcd(self, IntegerMod_int64 other):
        r"""
        Greatest common divisor.

        Returns the "smallest" generator in `\ZZ / N\ZZ` of the ideal
        generated by ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of the same ring as this one

        EXAMPLES::

            sage: mod(2^3*3^2*5, 3^3*2^2*17^5).gcd(mod(2^4*3*17, 3^3*2^2*17^5))
            12
            sage: mod(0,17^5).gcd(mod(0,17^5))
            0
        """
        cdef int_fast64_t g = gcd_int64(self.ivalue, self._modulus.int64)
        g = gcd_int64(g, other.ivalue)
        if g == self._modulus.int64: # self = other = 0
            g = 0
        return self._new_c(g)


### Helper functions

cdef int mpz_pow_helper(mpz_t res, mpz_t base, object exp, mpz_t modulus) except -1:
    cdef bint invert = False
    cdef long long_exp
    if is_small_python_int(exp):
        long_exp = exp
        if long_exp < 0:
            long_exp = -long_exp
            invert = True
        mpz_powm_ui(res, base, long_exp, modulus)
    else:
        if type(exp) is not Integer:
            exp = Integer(exp)
        if mpz_sgn((<Integer>exp).value) < 0:
            exp = -exp
            invert = True
        mpz_powm(res, base, (<Integer>exp).value, modulus)
    if invert:
        if not mpz_invert(res, res, modulus):
            raise ZeroDivisionError("Inverse does not exist.")

cdef int_fast64_t gcd_int64(int_fast64_t a, int_fast64_t b) noexcept:
    """
    Return the gcd of ``a`` and ``b``.

    For use with IntegerMod_int64.

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
    Return y such that xy=1 mod n.

    For use in ``IntegerMod_int64``.

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
    raise ZeroDivisionError(f"inverse of Mod({x}, {n}) does not exist")


cdef int_fast64_t mod_pow_int64(int_fast64_t base, int_fast64_t exp, int_fast64_t n) noexcept:
    """
    Return base^exp mod n.

    For use in ``IntegerMod_int64``.

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
    Calculate the jacobi symbol (a/n).

    For use in ``IntegerMod_int64``.

    AUTHORS:

    - Robert Bradshaw
    """
    cdef int s, jacobi = 1
    cdef int_fast64_t b

    a = a % m

    while True:
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
    Calculate the square root of `a`, where `a` is an
    integer mod `p^e`.

    ALGORITHM: Compute `p`-adically by stripping off even powers of `p`
    to get a unit and lifting `\sqrt{unit} \bmod p` via Newton's method
    whenever `p` is odd and by a variant of Hensel lifting for `p = 2`.

    AUTHORS:

    - Robert Bradshaw
    - Lorenz Panny (2022): polynomial-time algorithm for `p = 2`

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import square_root_mod_prime_power
        sage: a = Mod(17,2^20)
        sage: b = square_root_mod_prime_power(a,2,20)
        sage: b^2 == a
        True

    ::

        sage: a = Mod(72, 97^10)
        sage: b = square_root_mod_prime_power(a, 97, 10)                                # needs sage.libs.pari
        sage: b^2 == a                                                                  # needs sage.libs.pari
        True
        sage: mod(100, 5^7).sqrt()^2                                                    # needs sage.libs.pari
        100

    TESTS:

    A big example for the binary case (:issue:`33961`)::

        sage: y = Mod(-7, 2^777)
        sage: hex(y.sqrt()^2 - y)                                                       # needs sage.libs.pari
        '0x0'

    Testing with random squares in random rings::

        sage: p = random_prime(999)
        sage: e = randrange(1, 999)
        sage: x = Zmod(p^e).random_element()
        sage: (x^2).sqrt()^2 == x^2                                                     # needs sage.libs.pari
        True
    """
    if a.is_zero() or a.is_one():
        return a

    # strip off even powers of p
    cdef int i, val = a.lift().valuation(p)
    if val % 2 == 1:
        raise ValueError("self must be a square")
    if val > 0:
        unit = a._parent(a.lift() // p**val)
    else:
        unit = a

    cdef int n

    if p == 2:
        # squares in Z/2^e are of the form 4^n*(1+8*m)
        if unit.lift() % 8 != 1:
            raise ValueError("self must be a square")

        u = unit.lift()
        x = next(i for i in range(1,8,2) if i*i & 31 == u & 31)
        t = (x*x - u) >> 5
        for i in range(4, e-1 - val//2):
            if t & 1:
                x |= one_Z << i
                t += x - (one_Z << i-1)
            t >>= 1
#            assert t << i+2 == x*x - u
        x = a.parent()(x)

    else:
        # find square root of unit mod p
        x = unit.parent()(square_root_mod_prime(mod(unit, p), p))

        # lift p-adically using Newton iteration
        # this is done to higher precision than necessary except at the last step
        one_half = ~(a._new_c_from_long(2))
        # need at least (e - val//2) p-adic digits of precision, which doubles
        # at each step
        n = <int>ceil(log2(e - val//2))
        for i in range(n):
            x = (x + unit/x) * one_half

    # multiply in powers of p (if any)
    if val > 0:
        x *= p**(val//2)
    return x


cpdef square_root_mod_prime(IntegerMod_abstract a, p=None):
    r"""
    Calculate the square root of `a`, where `a` is an
    integer mod `p`; if `a` is not a perfect square,
    this returns an (incorrect) answer without checking.

    ALGORITHM: Several cases based on residue class of
    `p \bmod 16`.


    - `p \bmod 2 = 0`: `p = 2` so `\sqrt{a} = a`.

    - `p \bmod 4 = 3`: `\sqrt{a} = a^{(p+1)/4}`.

    - `p \bmod 8 = 5`: `\sqrt{a} = \zeta i a` where `\zeta = (2a)^{(p-5)/8}`,
      `i=\sqrt{-1}`.

    - `p \bmod 16 = 9`: Similar, work in a bi-quadratic extension of `\GF{p}`
      for small `p`, Tonelli and Shanks for large `p`.

    - `p \bmod 16 = 1`: Tonelli and Shanks.


    REFERENCES:

    - [Mul2004]_

    - [Atk1992]_

    - [Pos1988]_

    AUTHORS:

    - Robert Bradshaw

    TESTS:

    Every case appears in the first hundred primes.

    ::

        sage: from sage.rings.finite_rings.integer_mod import square_root_mod_prime   # sqrt() uses brute force for small p
        sage: all(square_root_mod_prime(a*a)^2 == a*a                                   # needs sage.libs.pari
        ....:     for p in prime_range(100)
        ....:     for a in Integers(p))
        True
    """
    if not a or a.is_one():
        return a

    if p is None:
        p = a._parent.order()
    p = Integer(p)

    cdef int p_mod_16 = p % 16
    cdef double bits = log2(float(p))
    cdef long r, m

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


def lucas_q1(mm, IntegerMod_abstract P):
    """
    Return `V_k(P, 1)` where `V_k` is the Lucas
    function defined by the recursive relation.

    `V_k(P, Q) = PV_{k-1}(P, Q) -  QV_{k-2}(P, Q)`

    with `V_0 = 2, V_1(P_Q) = P`.

    REFERENCES:

    - [Pos1988]_

    AUTHORS:

    - Robert Bradshaw

    TESTS::

        sage: from sage.rings.finite_rings.integer_mod import lucas_q1
        sage: all(lucas_q1(k, a) == BinaryRecurrenceSequence(a, -1, 2, a)(k)            # needs sage.combinat sage.modules
        ....:     for a in Integers(23)
        ....:     for k in range(13))
        True
    """
    if mm == 0:
        return 2
    elif mm == 1:
        return P

    cdef sage.rings.integer.Integer m
    m = <sage.rings.integer.Integer>mm if isinstance(mm, sage.rings.integer.Integer) else sage.rings.integer.Integer(mm)
    two = P._new_c_from_long(2)
    d1 = P
    d2 = P*P - two

    cdef int j
    for j from mpz_sizeinbase(m.value, 2)-1 > j > 0:
        sig_check()
        if mpz_tstbit(m.value, j):
            d1 = d1*d2 - P
            d2 = d2*d2 - two
        else:
            d2 = d1*d2 - P
            d1 = d1*d1 - two
    if mpz_odd_p(m.value):
        return d1*d2 - P
    else:
        return d1*d1 - two


def lucas(k, P, Q=1, n=None):
    r"""
    Return `[V_k(P, Q) \mod n, Q^{\lfloor k/2 \rfloor} \mod n]` where `V_k`
    is the Lucas function defined by the recursive relation

    .. MATH::

        V_k(P, Q) = P V_{k-1}(P, Q) -  Q V_{k-2}(P, Q)

    with `V_0 = 2, V_1 = P`.

    INPUT:

    - ``k`` -- integer; index to compute

    - ``P``, ``Q`` -- integers or modular integers; initial values

    - ``n`` -- integer (optional); modulus to use if ``P`` is not a modular
      integer

    REFERENCES:

    - [IEEEP1363]_

    AUTHORS:

    - Somindu Chaya Ramanna, Shashank Singh and Srinivas Vivek Venkatesh
      (2011-09-15, ECC2011 summer school)

    - Robert Bradshaw

    TESTS::

        sage: from sage.rings.finite_rings.integer_mod import lucas
        sage: p = randint(0,100000)
        sage: q = randint(0,100000)
        sage: n = randint(1,100)
        sage: all(lucas(k, p, q, n)[0] == Mod(lucas_number2(k, p, q), n)                # needs sage.combinat sage.libs.gap
        ....:     for k in Integers(20))
        True
        sage: from sage.rings.finite_rings.integer_mod import lucas
        sage: p = randint(0,100000)
        sage: q = randint(0,100000)
        sage: n = randint(1,100)
        sage: k = randint(0,100)
        sage: lucas(k, p, q, n) == [Mod(lucas_number2(k, p, q), n),                     # needs sage.combinat sage.libs.gap
        ....:                       Mod(q^(int(k/2)), n)]
        True

    EXAMPLES::

        sage: [lucas(k,4,5,11)[0] for k in range(30)]
        [2, 4, 6, 4, 8, 1, 8, 5, 2, 5, 10, 4, 10, 9, 8, 9, 7, 5, 7, 3, 10, 3, 6, 9, 6, 1, 7, 1, 2, 3]

        sage: lucas(20,4,5,11)
        [10, 1]
    """
    cdef IntegerMod_abstract p,q

    if n is None and not isinstance(P, IntegerMod_abstract):
        raise ValueError

    if n is None:
        n = P.modulus()

    if not isinstance(P, IntegerMod_abstract):
        p = Mod(P,n)
    else:
        p = P

    if not isinstance(Q, IntegerMod_abstract):
        q = Mod(Q,n)
    else:
        q = Q

    if k == 0:
        return [2, 1]
    elif k == 1:
        return [p, 1]

    cdef sage.rings.integer.Integer m
    m = <sage.rings.integer.Integer>k if isinstance(k, sage.rings.integer.Integer) else sage.rings.integer.Integer(k)
    two = p._new_c_from_long(2)

    v0 = p._new_c_from_long(2)
    v1 = p
    q0 = p._new_c_from_long(1)
    q1 = p._new_c_from_long(1)

    cdef int j
    for j from mpz_sizeinbase(m.value, 2)-1 >= j >= 0:
        sig_check()
        q0 = q0*q1
        if mpz_tstbit(m.value, j):
            q1 = q0*Q
            v0 = v0*v1 - p*q0
            v1 = v1*v1 - two*q1
        else:
            q1 = q0
            v1 = v0*v1 - p*q0
            v0 = v0*v0 - two*q0
    return [v0,q0]


############# Homomorphisms ###############

cdef class IntegerMod_hom(Morphism):
    cdef IntegerMod_abstract zero
    cdef NativeIntStruct modulus

    def __init__(self, parent):
        Morphism.__init__(self, parent)
        # we need to use element constructor so that we can register both coercions and conversions using these morphisms.
        cdef Parent C = self._codomain
        self.zero = C._element_constructor_(0)
        self.modulus = C._pyx_order

    cdef dict _extra_slots(self):
        """
        Helper for pickling and copying.

        EXAMPLES::

            sage: R5 = IntegerModRing(5)
            sage: R15 = IntegerModRing(15)
            sage: phi = R5.coerce_map_from(R15); phi
            Natural morphism:
              From: Ring of integers modulo 15
              To:   Ring of integers modulo 5

        This method helps to implement copying::

            sage: psi = copy(phi); psi
            Natural morphism:
              From: Ring of integers modulo 15
              To:   Ring of integers modulo 5
            sage: psi(R15(7))
            2
        """
        slots = Morphism._extra_slots(self)
        slots['zero'] = self.zero
        slots['modulus'] = self.modulus
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        EXAMPLES::

            sage: R5 = IntegerModRing(5)
            sage: R15 = IntegerModRing(15)
            sage: phi = R5.coerce_map_from(R15); phi
            Natural morphism:
              From: Ring of integers modulo 15
              To:   Ring of integers modulo 5

        This method helps to implement copying.
        ::

            sage: psi = copy(phi); psi
            Natural morphism:
              From: Ring of integers modulo 15
              To:   Ring of integers modulo 5
            sage: psi(R15(7))
            2
        """
        Morphism._update_slots(self, _slots)
        self.zero = _slots['zero']
        self.modulus = _slots['modulus']

    cpdef Element _call_(self, x):
        return IntegerMod(self._codomain, x)

cdef class IntegerMod_to_IntegerMod(IntegerMod_hom):
    """
    Very fast IntegerMod to IntegerMod homomorphism.

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod import IntegerMod_to_IntegerMod
        sage: Rs = [Integers(3**k) for k in range(1,30,5)]
        sage: [type(R(0)) for R in Rs]
        [<class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>]
        sage: fs = [IntegerMod_to_IntegerMod(S, R)
        ....:       for R in Rs for S in Rs if S is not R and S.order() > R.order()]
        sage: all(f(-1) == f.codomain()(-1) for f in fs)
        True
        sage: [f(-1) for f in fs]
        [2, 2, 2, 2, 2, 728, 728, 728, 728, 177146, 177146, 177146, 43046720, 43046720, 10460353202]
    """
    def __init__(self, R, S):
        if not S.order().divides(R.order()):
            raise TypeError("No natural coercion from %s to %s" % (R, S))
        import sage.categories.homset
        IntegerMod_hom.__init__(self, sage.categories.homset.Hom(R, S))

    cpdef Element _call_(self, x):
        cdef IntegerMod_abstract a
        zero = <IntegerMod_abstract>self.zero
        cdef unsigned long value
        if isinstance(x, IntegerMod_int):
            value = (<IntegerMod_int>x).ivalue
            value %= <unsigned long>self.modulus.int32
            return zero._new_c_fast(value)
        elif isinstance(x, IntegerMod_int64):
            value = (<IntegerMod_int64>x).ivalue
            value %= <unsigned long>self.modulus.int64
            return zero._new_c_fast(value)
        a = zero._new_c_fast(0)
        a.set_from_mpz((<IntegerMod_gmp?>x).value)
        return a

    def _repr_type(self):
        return "Natural"

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: Zmod(4).hom(Zmod(2)).is_surjective()
            True
        """
        return True

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: Zmod(4).hom(Zmod(2)).is_injective()
            False
        """
        return self.domain().order() == self.codomain().order()

cdef class Integer_to_IntegerMod(IntegerMod_hom):
    r"""
    Fast `\ZZ \rightarrow \ZZ/n\ZZ` morphism.

    EXAMPLES:

    We make sure it works for every type.

    ::

        sage: from sage.rings.finite_rings.integer_mod import Integer_to_IntegerMod
        sage: Rs = [Integers(10), Integers(10^5), Integers(10^10)]
        sage: [type(R(0)) for R in Rs]
        [<class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>]
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
            a = self.modulus.table[res]
#            if a._parent is not self._codomain:
            a._parent = self._codomain
            return a
        else:
            a = self.zero._new_c_from_long(0)
            a.set_from_mpz((<Integer>x).value)
            return a

    def _repr_type(self):
        return "Natural"

    def section(self):
        return IntegerMod_to_Integer(self._codomain)

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: ZZ.hom(Zmod(2)).is_surjective()
            True
        """
        return True

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: ZZ.hom(Zmod(2)).is_injective()
            False
        """
        return False


cdef class IntegerMod_to_Integer(Map):
    """
    Map to lift elements to :class:`~sage.rings.integer.Integer`.

    EXAMPLES::

        sage: ZZ.convert_map_from(GF(2))
        Lifting map:
          From: Finite Field of size 2
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        TESTS:

        Lifting maps are morphisms in the category of sets (see
        :issue:`15618`)::

            sage: ZZ.convert_map_from(GF(2)).parent()
            Set of Morphisms from Finite Field of size 2 to Integer Ring in Category of sets
        """
        import sage.categories.homset
        from sage.categories.sets_cat import Sets
        Morphism.__init__(self, sage.categories.homset.Hom(R, integer_ring.ZZ, Sets()))

    cpdef Element _call_(self, x):
        cdef Integer ans = Integer.__new__(Integer)
        if isinstance(x, IntegerMod_gmp):
            mpz_set(ans.value, (<IntegerMod_gmp>x).value)
        elif isinstance(x, IntegerMod_int):
            mpz_set_ui(ans.value, (<IntegerMod_int>x).ivalue)
        elif isinstance(x, IntegerMod_int64):
            mpz_set_ui(ans.value, (<IntegerMod_int64>x).ivalue)
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
        [<class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>,
         <class 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>]
        sage: fs = [Int_to_IntegerMod(R) for R in Rs]
        sage: [f(-1) for f in fs]
        [1, 2047, 2097151, 2147483647, 2199023255551]
    """
    def __init__(self, R):
        import sage.categories.homset
        from sage.sets.pythonclass import Set_PythonType
        IntegerMod_hom.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), R))

    cpdef Element _call_(self, x):
        cdef IntegerMod_abstract a
        zero = <IntegerMod_abstract>self.zero

        cdef long res
        cdef int err

        if not integer_check_long_py(x, &res, &err):
            raise TypeError(f"{x} is not an integer")

        if not err:
            return zero._new_c_from_long(res)

        cdef Integer z = Integer(x)
        a = zero._new_c_fast(0)
        a.set_from_mpz(z.value)
        return a

    def _repr_type(self):
        return "Native"
