"""
Sage class for PARI's GEN type

See the ``PariInstance`` class for documentation and examples.

AUTHORS:

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-03-06): added newtonpoly

- Justin Walker: contributed some of the function definitions

- Gonzalo Tornaria: improvements to conversions; much better error
  handling.

- Robert Bradshaw, Jeroen Demeyer, William Stein (2010-08-15):
  Upgrade to PARI 2.4.3 (:trac:`9343`)

- Jeroen Demeyer (2011-11-12): rewrite various conversion routines
  (:trac:`11611`, :trac:`11854`, :trac:`11952`)

- Peter Bruin (2013-11-17): move PariInstance to a separate file
  (:trac:`15185`)

- Jeroen Demeyer (2014-02-09): upgrade to PARI 2.7 (:trac:`15767`)

- Martin von Gagern (2014-12-17): Added some Galois functions (:trac:`17519`)

- Jeroen Demeyer (2015-01-12): upgrade to PARI 2.8 (:trac:`16997`)

- Jeroen Demeyer (2015-03-17): automatically generate methods from
  ``pari.desc`` (:trac:`17631` and :trac:`17860`)

- Kiran Kedlaya (2016-03-23): implement infinity type

TESTS:

Before :trac:`15654`, this used to take a very long time.
Now it takes much less than a second::

    sage: pari.allocatemem(200000)
    PARI stack size set to 200000 bytes, maximum size set to ...
    sage: x = polygen(ZpFM(3,10))
    sage: pol = ((x-1)^50 + x)
    sage: pari(pol).poldisc()
    2*3 + 3^4 + 2*3^6 + 3^7 + 2*3^8 + 2*3^9 + O(3^10)
"""

#*****************************************************************************
#       Copyright (C) 2006,2010 William Stein <wstein@gmail.com>
#       Copyright (C) ???? Justin Walker
#       Copyright (C) ???? Gonzalo Tornaria
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2010-2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

import types
cimport cython

from cpython.string cimport PyString_AsString
from cpython.int cimport PyInt_Check
from cpython.long cimport PyLong_Check
from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.complex cimport PyComplex_RealAsDouble, PyComplex_ImagAsDouble

include "cysignals/memory.pxi"
include "cysignals/signals.pxi"

from .paridecl cimport *
from .paripriv cimport *
from .convert cimport integer_to_gen, gen_to_integer
from .pari_instance cimport (PariInstance, pari_instance,
        prec_bits_to_words, prec_words_to_bits, default_bitprec)
from .pari_instance cimport pari_instance as P  # Shorthand

from sage.structure.element cimport ModuleElement, RingElement, Element
from sage.misc.randstate cimport randstate, current_randstate
from sage.structure.sage_object cimport rich_to_bool
from sage.misc.superseded import deprecation, deprecated_function_alias
from sage.libs.pari.closure cimport objtoclosure
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.infinity import Infinity


include 'auto_gen.pxi'

@cython.final
cdef class gen(gen_auto):
    """
    Cython extension class that models the PARI GEN type.
    """
    def __init__(self):
        raise RuntimeError("PARI objects cannot be instantiated directly; use pari(x) to convert x to PARI")

    def __dealloc__(self):
        sig_free(<void*>self.b)

    def __repr__(self):
        """
        Display representation of a gen.

        OUTPUT: a Python string

        EXAMPLES::

            sage: pari('vector(5,i,i)')
            [1, 2, 3, 4, 5]
            sage: pari('[1,2;3,4]')
            [1, 2; 3, 4]
            sage: pari('Str(hello)')
            "hello"
        """
        cdef char *c
        sig_on()
        # Use sig_block(), which is needed because GENtostr() uses
        # malloc(), which is dangerous inside sig_on()
        sig_block()
        c = GENtostr(self.g)
        sig_unblock()
        sig_off()

        s = str(c)
        pari_free(c)
        return s

    def __str__(self):
        """
        Convert this gen to a string.

        Except for PARI strings, we have ``str(x) == repr(x)``.
        For strings (type ``t_STR``), the returned string is not quoted.

        OUTPUT: a Python string

        EXAMPLES::

            sage: print(pari('vector(5,i,i)'))
            [1, 2, 3, 4, 5]
            sage: print(pari('[1,2;3,4]'))
            [1, 2; 3, 4]
            sage: print(pari('Str(hello)'))
            hello
        """
        # Use __repr__ except for strings
        if typ(self.g) == t_STR:
            return GSTR(self.g)
        return repr(self)

    def __hash__(self):
        """
        Return the hash of self, computed using PARI's hash_GEN().

        TESTS::

            sage: type(pari('1 + 2.0*I').__hash__())
            <type 'int'>
        """
        cdef long h
        sig_on()
        h = hash_GEN(self.g)
        sig_off()
        return h

    def list(self):
        """
        Convert self to a list of PARI gens.

        EXAMPLES:

        A PARI vector becomes a Sage list::

            sage: L = pari("vector(10,i,i^2)").list()
            sage: L
            [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
            sage: type(L)
            <type 'list'>
            sage: type(L[0])
            <type 'sage.libs.pari.gen.gen'>

        For polynomials, list() behaves as for ordinary Sage polynomials::

            sage: pol = pari("x^3 + 5/3*x"); pol.list()
            [0, 5/3, 0, 1]

        For power series or Laurent series, we get all coefficients starting
        from the lowest degree term.  This includes trailing zeros::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: s = x^2 + O(x^8)
            sage: s.list()
            [1]
            sage: pari(s).list()
            [1, 0, 0, 0, 0, 0]
            sage: s = x^-2 + O(x^0)
            sage: s.list()
            [1]
            sage: pari(s).list()
            [1, 0]

        For matrices, we get a list of columns::

            sage: M = matrix(ZZ,3,2,[1,4,2,5,3,6]); M
            [1 4]
            [2 5]
            [3 6]
            sage: pari(M).list()
            [[1, 2, 3]~, [4, 5, 6]~]

        For "scalar" types, we get a 1-element list containing ``self``::

            sage: pari("42").list()
            [42]
        """
        if typ(self.g) == t_POL:
            return list(self.Vecrev())
        return list(self.Vec())

    def __reduce__(self):
        """
        EXAMPLES::

            sage: f = pari('x^3 - 3')
            sage: loads(dumps(f)) == f
            True
            sage: f = pari('"hello world"')
            sage: loads(dumps(f)) == f
            True
        """
        s = repr(self)
        return (objtogen, (s,))

    cpdef _add_(self, right):
        sig_on()
        return P.new_gen(gadd(self.g, (<gen>right).g))

    cpdef _sub_(self, right):
        sig_on()
        return P.new_gen(gsub(self.g, (<gen> right).g))

    cpdef _mul_(self, right):
        sig_on()
        return P.new_gen(gmul(self.g, (<gen>right).g))

    cpdef _div_(self, right):
        sig_on()
        return P.new_gen(gdiv(self.g, (<gen>right).g))

    def _add_one(gen self):
        """
        Return self + 1.

        OUTPUT: pari gen

        EXAMPLES::

            sage: n = pari(5)
            sage: n._add_one()
            6
            sage: n = pari('x^3')
            sage: n._add_one()
            x^3 + 1
        """
        sig_on()
        return P.new_gen(gaddsg(1, self.g))

    def __mod__(self, other):
        """
        Return ``self`` modulo ``other``.

        EXAMPLES::

            sage: pari(15) % pari(6)
            3
            sage: pari("x^3+x^2+x+1") % pari("x^2")
            x + 1
            sage: pari(-2) % int(3)
            1
            sage: int(-2) % pari(3)
            1
        """
        cdef gen selfgen = objtogen(self)
        cdef gen othergen = objtogen(other)
        sig_on()
        return P.new_gen(gmod(selfgen.g, othergen.g))

    def __pow__(self, n, m):
        """
        Return ``self`` to the power ``n`` (if ``m`` is ``None``) or
        ``Mod(self, m)^n`` if ``m`` is not ``None``.

        EXAMPLES::

            sage: pari(5) ^ pari(3)
            125
            sage: pari("x-1") ^ 3
            x^3 - 3*x^2 + 3*x - 1
            sage: pow(pari(5), pari(28), int(29))
            Mod(1, 29)
            sage: int(2) ^ pari(-5)
            1/32
            sage: pari(2) ^ int(-5)
            1/32
        """
        cdef gen t0 = objtogen(self)
        cdef gen t1 = objtogen(n)
        if m is not None:
            t0 = t0.Mod(m)
        sig_on()
        return P.new_gen(gpow(t0.g, t1.g, prec_bits_to_words(0)))

    def __neg__(gen self):
        sig_on()
        return P.new_gen(gneg(self.g))

    def __rshift__(self, long n):
        """
        Divide ``self`` by `2^n` (truncating or not, depending on the
        input type).

        EXAMPLES::

            sage: pari(25) >> 3
            3
            sage: pari(25/2) >> 2
            25/8
            sage: pari("x") >> 3
            1/8*x
            sage: pari(1.0) >> 100
            7.88860905221012 E-31
            sage: int(33) >> pari(2)
            8
        """
        cdef gen t0 = objtogen(self)
        sig_on()
        return P.new_gen(gshift(t0.g, -n))

    def __lshift__(self, long n):
        """
        Multiply ``self`` by `2^n`.

        EXAMPLES::

            sage: pari(25) << 3
            200
            sage: pari(25/32) << 2
            25/8
            sage: pari("x") << 3
            8*x
            sage: pari(1.0) << 100
            1.26765060022823 E30
            sage: int(33) << pari(2)
            132
        """
        cdef gen t0 = objtogen(self)
        sig_on()
        return P.new_gen(gshift(t0.g, n))

    def __invert__(gen self):
        sig_on()
        return P.new_gen(ginv(self.g))

    def getattr(gen self, attr):
        """
        Return the PARI attribute with the given name.

        EXAMPLES::

            sage: K = pari("nfinit(x^2 - x - 1)")
            sage: K.getattr("pol")
            x^2 - x - 1
            sage: K.getattr("disc")
            5

            sage: K.getattr("reg")
            Traceback (most recent call last):
            ...
            PariError: _.reg: incorrect type in reg (t_VEC)
            sage: K.getattr("zzz")
            Traceback (most recent call last):
            ...
            PariError: not a function in function call
        """
        cdef str s = "_." + attr
        cdef char *t = PyString_AsString(s)
        sig_on()
        return P.new_gen(closure_callgen1(strtofunction(t), self.g))

    def mod(self):
        """
        Given an INTMOD or POLMOD ``Mod(a,m)``, return the modulus `m`.

        EXAMPLES::

            sage: pari(4).Mod(5).mod()
            5
            sage: pari("Mod(x, x*y)").mod()
            y*x
            sage: pari("[Mod(4,5)]").mod()
            Traceback (most recent call last):
            ...
            TypeError: Not an INTMOD or POLMOD in mod()
        """
        if typ(self.g) != t_INTMOD and typ(self.g) != t_POLMOD:
            raise TypeError("Not an INTMOD or POLMOD in mod()")
        sig_on()
        # The hardcoded 1 below refers to the position in the internal
        # representation of a INTMOD or POLDMOD where the modulus is
        # stored.
        return P.new_gen(gel(self.g, 1))

    def nf_get_pol(self):
        """
        Returns the defining polynomial of this number field.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: pari(K).nf_get_pol()
            y^4 - 4*y^2 + 1
            sage: bnr = pari("K = bnfinit(x^4 - 4*x^2 + 1); bnrinit(K, 2*x)")
            sage: bnr.nf_get_pol()
            x^4 - 4*x^2 + 1

        For relative number fields, this returns the relative
        polynomial. However, beware that ``pari(L)`` returns an absolute
        number field::

            sage: L.<b> = K.extension(x^2 - 5)
            sage: pari(L).nf_get_pol()        # Absolute
            y^8 - 28*y^6 + 208*y^4 - 408*y^2 + 36
            sage: L.pari_rnf().nf_get_pol()   # Relative
            x^2 - 5

        TESTS::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: K.pari_nf().nf_get_pol()
            y^4 - 4*y^2 + 1
            sage: K.pari_bnf().nf_get_pol()
            y^4 - 4*y^2 + 1

        An error is raised for invalid input::

            sage: pari("[0]").nf_get_pol()
            Traceback (most recent call last):
            ...
            PariError: incorrect type in pol (t_VEC)

        """
        sig_on()
        return P.new_gen(member_pol(self.g))

    def nf_get_diff(self):
        """
        Returns the different of this number field as a PARI ideal.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: pari(K).nf_get_diff()
            [12, 0, 0, 0; 0, 12, 8, 0; 0, 0, 4, 0; 0, 0, 0, 4]
        """
        sig_on()
        return P.new_gen(member_diff(self.g))

    def nf_get_sign(self):
        """
        Returns a Python list ``[r1, r2]``, where ``r1`` and ``r2`` are
        Python ints representing the number of real embeddings and pairs
        of complex embeddings of this number field, respectively.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: s = K.pari_nf().nf_get_sign(); s
            [4, 0]
            sage: type(s); type(s[0])
            <type 'list'>
            <type 'int'>
            sage: CyclotomicField(15).pari_nf().nf_get_sign()
            [0, 4]
        """
        cdef long r1
        cdef long r2
        cdef GEN sign
        sig_on()
        sign = member_sign(self.g)
        r1 = itos(gel(sign, 1))
        r2 = itos(gel(sign, 2))
        sig_off()
        return [r1, r2]

    def nf_get_zk(self):
        """
        Returns a vector with a `\ZZ`-basis for the ring of integers of
        this number field. The first element is always `1`.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: pari(K).nf_get_zk()
            [1, y, y^3 - 4*y, y^2 - 2]
        """
        sig_on()
        return P.new_gen(member_zk(self.g))

    def bnf_get_no(self):
        """
        Returns the class number of ``self``, a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-65)
            sage: K.pari_bnf().bnf_get_no()
            8
        """
        sig_on()
        return P.new_gen(bnf_get_no(self.g))

    def bnf_get_cyc(self):
        """
        Returns the structure of the class group of this number field as
        a vector of SNF invariants.

        NOTE: ``self`` must be a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-65)
            sage: K.pari_bnf().bnf_get_cyc()
            [4, 2]
        """
        sig_on()
        return P.new_gen(bnf_get_cyc(self.g))

    def bnf_get_gen(self):
        """
        Returns a vector of generators of the class group of this
        number field.

        NOTE: ``self`` must be a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-65)
            sage: G = K.pari_bnf().bnf_get_gen(); G
            [[3, 2; 0, 1], [2, 1; 0, 1]]
            sage: map(lambda J: K.ideal(J), G)
            [Fractional ideal (3, a + 2), Fractional ideal (2, a + 1)]
        """
        sig_on()
        return P.new_gen(bnf_get_gen(self.g))

    def bnf_get_reg(self):
        """
        Returns the regulator of this number field.

        NOTE: ``self`` must be a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: K.pari_bnf().bnf_get_reg()
            2.66089858019037...
        """
        sig_on()
        return P.new_gen(bnf_get_reg(self.g))

    def bnfunit(self):
        sig_on()
        return P.new_gen(bnf_get_fu(self.g))

    def pr_get_p(self):
        """
        Returns the prime of `\ZZ` lying below this prime ideal.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: F = pari(K).idealfactor(K.ideal(5)); F
            [[5, [-2, 1]~, 1, 1, [2, -1; 1, 2]], 1; [5, [2, 1]~, 1, 1, [-2, -1; 1, -2]], 1]
            sage: F[0,0].pr_get_p()
            5
        """
        sig_on()
        return P.new_gen(pr_get_p(self.g))

    def pr_get_e(self):
        """
        Returns the ramification index (over `\QQ`) of this prime ideal.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: pari(K).idealfactor(K.ideal(2))[0,0].pr_get_e()
            2
            sage: pari(K).idealfactor(K.ideal(3))[0,0].pr_get_e()
            1
            sage: pari(K).idealfactor(K.ideal(5))[0,0].pr_get_e()
            1
        """
        cdef long e
        sig_on()
        e = pr_get_e(self.g)
        sig_off()
        return e

    def pr_get_f(self):
        """
        Returns the residue class degree (over `\QQ`) of this prime ideal.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: pari(K).idealfactor(K.ideal(2))[0,0].pr_get_f()
            1
            sage: pari(K).idealfactor(K.ideal(3))[0,0].pr_get_f()
            2
            sage: pari(K).idealfactor(K.ideal(5))[0,0].pr_get_f()
            1
        """
        cdef long f
        sig_on()
        f = pr_get_f(self.g)
        sig_off()
        return f

    def pr_get_gen(self):
        """
        Returns the second generator of this PARI prime ideal, where the
        first generator is ``self.pr_get_p()``.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: g = pari(K).idealfactor(K.ideal(2))[0,0].pr_get_gen(); g; K(g)
            [1, 1]~
            i + 1
            sage: g = pari(K).idealfactor(K.ideal(3))[0,0].pr_get_gen(); g; K(g)
            [3, 0]~
            3
            sage: g = pari(K).idealfactor(K.ideal(5))[0,0].pr_get_gen(); g; K(g)
            [-2, 1]~
            i - 2
        """
        sig_on()
        return P.new_gen(pr_get_gen(self.g))

    def bid_get_cyc(self):
        """
        Returns the structure of the group `(O_K/I)^*`, where `I` is the
        ideal represented by ``self``.

        NOTE: ``self`` must be a "big ideal" (``bid``) as returned by
        ``idealstar`` for example.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: J = pari(K).idealstar(K.ideal(4*i + 2))
            sage: J.bid_get_cyc()
            [4, 2]
        """
        sig_on()
        return P.new_gen(bid_get_cyc(self.g))

    def bid_get_gen(self):
        """
        Returns a vector of generators of the group `(O_K/I)^*`, where
        `I` is the ideal represented by ``self``.

        NOTE: ``self`` must be a "big ideal" (``bid``) with generators,
        as returned by ``idealstar`` with ``flag`` = 2.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: J = pari(K).idealstar(K.ideal(4*i + 2), 2)
            sage: J.bid_get_gen()
            [7, [-2, -1]~]

        We get an exception if we do not supply ``flag = 2`` to
        ``idealstar``::

            sage: J = pari(K).idealstar(K.ideal(3))
            sage: J.bid_get_gen()
            Traceback (most recent call last):
            ...
            PariError: missing bid generators. Use idealstar(,,2)
        """
        sig_on()
        return P.new_gen(bid_get_gen(self.g))

    def __getitem__(gen self, n):
        """
        Return the nth entry of self. The indexing is 0-based, like in
        Python. Note that this is *different* than the default behavior
        of the PARI/GP interpreter.

        EXAMPLES::

            sage: p = pari('1 + 2*x + 3*x^2')
            sage: p[0]
            1
            sage: p[2]
            3
            sage: p[100]
            0
            sage: p[-1]
            0
            sage: q = pari('x^2 + 3*x^3 + O(x^6)')
            sage: q[3]
            3
            sage: q[5]
            0
            sage: q[6]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: m = pari('[1,2;3,4]')
            sage: m[0]
            [1, 3]~
            sage: m[1,0]
            3
            sage: l = pari('List([1,2,3])')
            sage: l[1]
            2
            sage: s = pari('"hello, world!"')
            sage: s[0]
            'h'
            sage: s[4]
            'o'
            sage: s[12]
            '!'
            sage: s[13]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: v = pari('[1,2,3]')
            sage: v[0]
            1
            sage: c = pari('Col([1,2,3])')
            sage: c[1]
            2
            sage: sv = pari('Vecsmall([1,2,3])')
            sage: sv[2]
            3
            sage: type(sv[2])
            <type 'int'>
            sage: tuple(pari(3/5))
            (3, 5)
            sage: tuple(pari('1 + 5*I'))
            (1, 5)
            sage: tuple(pari('Qfb(1, 2, 3)'))
            (1, 2, 3)
            sage: pari(57)[0]
            Traceback (most recent call last):
            ...
            TypeError: PARI object of type 't_INT' cannot be indexed
            sage: m = pari("[[1,2;3,4],5]") ; m[0][1,0]
            3
            sage: v = pari(xrange(20))
            sage: v[2:5]
            [2, 3, 4]
            sage: v[:]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[10:]
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[:5]
            [0, 1, 2, 3, 4]
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[-1]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: v[:-3]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            sage: v[5:]
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: pari([])[::]
            []
        """
        cdef int pari_type

        pari_type = typ(self.g)

        if isinstance(n, tuple):
            if pari_type != t_MAT:
                raise TypeError("self must be of pari type t_MAT")
            if len(n) != 2:
                raise IndexError("index must be an integer or a 2-tuple (i,j)")
            i = int(n[0])
            j = int(n[1])

            if i < 0 or i >= glength(<GEN>(self.g[1])):
                raise IndexError("row index out of range")
            if j < 0 or j >= glength(self.g):
                raise IndexError("column index out of range")

            ind = (i,j)

            if self.refers_to is not None and ind in self.refers_to:
                return self.refers_to[ind]
            else:
                ## In this case, we're being asked to return
                ## a GEN that has no gen pointing to it, so
                ## we need to create such a gen, add it to
                ## self.refers_to, and return it.
                val = P.new_ref(gmael(self.g, j+1, i+1), self)
                if self.refers_to is None:
                    self.refers_to = {ind: val}
                else:
                    self.refers_to[ind] = val
                return val

        elif isinstance(n, slice):
            l = glength(self.g)
            start,stop,step = n.indices(l)
            inds = xrange(start,stop,step)
            k = len(inds)
            # fast exit
            if k==0:
                return P.vector(0)
            # fast call, beware pari is one based
            if pari_type == t_VEC:
                if step==1:
                    return self.vecextract('"'+str(start+1)+".."+str(stop)+'"')
                if step==-1:
                    return self.vecextract('"'+str(start+1)+".."+str(stop+2)+'"')
            # slow call
            v = P.vector(k)
            for i,j in enumerate(inds):
                v[i] = self[j]
            return v

        ## there are no "out of bounds" problems
        ## for a polynomial or power series, so these go before
        ## bounds testing
        if pari_type == t_POL:
            return self.polcoeff(n)

        elif pari_type == t_SER:
            bound = valp(self.g) + lg(self.g) - 2
            if n >= bound:
                raise IndexError("index out of range")
            return self.polcoeff(n)

        elif pari_type in (t_INT, t_REAL, t_PADIC, t_QUAD, t_FFELT, t_INTMOD, t_POLMOD):
            # these are definitely scalar!
            raise TypeError("PARI object of type %r cannot be indexed" % self.type())

        elif n < 0 or n >= glength(self.g):
            raise IndexError("index out of range")

        elif pari_type == t_VEC or pari_type == t_MAT:
            #t_VEC    : row vector        [ code ] [  x_1  ] ... [  x_k  ]
            #t_MAT    : matrix            [ code ] [ col_1 ] ... [ col_k ]
            if self.refers_to is not None and n in self.refers_to:
                return self.refers_to[n]
            else:
                ## In this case, we're being asked to return
                ## a GEN that has no gen pointing to it, so
                ## we need to create such a gen, add it to
                ## self.refers_to, and return it.
                val = P.new_ref(gel(self.g, n+1), self)
                if self.refers_to is None:
                    self.refers_to = {n: val}
                else:
                    self.refers_to[n] = val
                return val

        elif pari_type == t_VECSMALL:
            #t_VECSMALL: vec. small ints  [ code ] [ x_1 ] ... [ x_k ]
            return self.g[n+1]

        elif pari_type == t_STR:
            #t_STR    : string            [ code ] [ man_1 ] ... [ man_k ]
            return chr( (<char *>(self.g+1))[n] )

        elif pari_type == t_LIST:
            return self.component(n+1)

        #elif pari_type in (t_FRAC, t_RFRAC):
            # generic code gives us:
            #   [0] = numerator
            #   [1] = denominator

        #elif pari_type == t_COMPLEX:
            # generic code gives us
            #   [0] = real part
            #   [1] = imag part

        #elif type(self.g) in (t_QFR, t_QFI):
            # generic code works ok

        else:
            ## generic code, which currently handles cases
            ## as mentioned above
            return P.new_ref(gel(self.g,n+1), self)

    def __setitem__(gen self, n, y):
        r"""
        Set the nth entry to a reference to y.


            -  The indexing is 0-based, like everywhere else in Python, but
               *unlike* in PARI/GP.

            -  Assignment sets the nth entry to a reference to y, assuming y is
               an object of type gen. This is the same as in Python, but
               *different* than what happens in the gp interpreter, where
               assignment makes a copy of y.

            -  Because setting creates references it is *possible* to make
               circular references, unlike in GP. Do *not* do this (see the
               example below). If you need circular references, work at the Python
               level (where they work well), not the PARI object level.



        EXAMPLES::

            sage: v = pari(range(10))
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[0] = 10
            sage: w = pari([5,8,-20])
            sage: v
            [10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[1] = w
            sage: v
            [10, [5, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w[0] = -30
            sage: v
            [10, [-30, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: t = v[1]; t[1] = 10 ; v
            [10, [-30, 10, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[1][0] = 54321 ; v
            [10, [54321, 10, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w
            [54321, 10, -20]
            sage: v = pari([[[[0,1],2],3],4]) ; v[0][0][0][1] = 12 ; v
            [[[[0, 12], 2], 3], 4]
            sage: m = pari(matrix(2,2,range(4))) ; l = pari([5,6]) ; n = pari(matrix(2,2,[7,8,9,0])) ; m[1,0] = l ; l[1] = n ; m[1,0][1][1,1] = 1111 ; m
            [0, 1; [5, [7, 8; 9, 1111]], 3]
            sage: m = pari("[[1,2;3,4],5,6]") ; m[0][1,1] = 11 ; m
            [[1, 2; 3, 11], 5, 6]

        Finally, we create a circular reference::

            sage: v = pari([0])
            sage: w = pari([v])
            sage: v
            [0]
            sage: w
            [[0]]
            sage: v[0] = w

        Now there is a circular reference. Accessing v[0] will crash Sage.

        ::

            sage: s=pari.vector(2,[0,0])
            sage: s[:1]
            [0]
            sage: s[:1]=[1]
            sage: s
            [1, 0]
            sage: type(s[0])
            <type 'sage.libs.pari.gen.gen'>
            sage: s = pari(range(20)) ; s
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: s[0:10:2] = range(50,55) ; s
            [50, 1, 51, 3, 52, 5, 53, 7, 54, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: s[10:20:3] = range(100,150) ; s
            [50, 1, 51, 3, 52, 5, 53, 7, 54, 9, 100, 11, 12, 101, 14, 15, 102, 17, 18, 103]

        TESTS::

            sage: v = pari(xrange(10)) ; v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[:] = [20..29]
            sage: v
            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
            sage: type(v[0])
            <type 'sage.libs.pari.gen.gen'>
        """
        cdef int i, j
        cdef gen x = objtogen(y)
        cdef long l
        cdef Py_ssize_t ii, jj, step

        sig_on()
        try:
            if isinstance(n, tuple):
                if typ(self.g) != t_MAT:
                    raise TypeError("cannot index PARI type %s by tuple" % typ(self.g))

                if len(n) != 2:
                    raise ValueError("matrix index must be of the form [row, column]")

                i = int(n[0])
                j = int(n[1])
                ind = (i,j)

                if i < 0 or i >= glength(<GEN>(self.g[1])):
                    raise IndexError("row i(=%s) must be between 0 and %s" % (i, self.nrows()-1))
                if j < 0 or j >= glength(self.g):
                    raise IndexError("column j(=%s) must be between 0 and %s" % (j, self.ncols()-1))
                if self.refers_to is None:
                    self.refers_to = {ind: x}
                else:
                    self.refers_to[ind] = x

                (<GEN>(self.g)[j+1])[i+1] = <long>(x.g)
                return

            elif isinstance(n, slice):
                l = glength(self.g)
                inds = xrange(*n.indices(l))
                k = len(inds)
                if k > len(y):
                    raise ValueError("attempt to assign sequence of size %s to slice of size %s" % (len(y), k))

                # actually set the values
                for i,j in enumerate(inds):
                    self[j] = y[i]
                return

            i = int(n)

            if i < 0 or i >= glength(self.g):
                raise IndexError("index (%s) must be between 0 and %s" % (i, glength(self.g)-1))

            # so python memory manager will work correctly
            # and not free x if PARI part of self is the
            # only thing pointing to it.
            if self.refers_to is None:
                self.refers_to = {i: x}
            else:
                self.refers_to[i] = x

            ## correct indexing for t_POLs
            if typ(self.g) == t_POL:
                i = i + 1

            ## actually set the value
            (self.g)[i+1] = <long>(x.g)
            return
        finally:
            sig_off()

    def __len__(gen self):
        return glength(self.g)

    cpdef _richcmp_(left, right, int op):
        """
        Compare ``left`` and ``right`` using ``op``.

        EXAMPLES::

            sage: a = pari(5)
            sage: b = 10
            sage: a < b
            True
            sage: a <= b
            True
            sage: a <= 5
            True
            sage: a > b
            False
            sage: a >= b
            False
            sage: a >= pari(10)
            False
            sage: a == 5
            True
            sage: a is 5
            False

            sage: pari(2.5) > None
            True
            sage: pari(3) == pari(3)
            True
            sage: pari('x^2 + 1') == pari('I-1')
            False
            sage: pari(I) == pari(I)
            True

        This does not define a total order.  An error is raised when
        applying inequality operators to non-ordered types::

            sage: pari("Mod(1,3)") <= pari("Mod(2,3)")
            Traceback (most recent call last):
            ...
            PariError: forbidden comparison t_INTMOD , t_INTMOD
            sage: pari("[0]") <= pari("0")
            Traceback (most recent call last):
            ...
            PariError: forbidden comparison t_VEC (1 elts) , t_INT

        TESTS:

        Check that :trac:`16127` has been fixed::

            sage: pari(1/2) < pari(1/3)
            False
            sage: pari(1) < pari(1/2)
            False

            sage: pari('O(x)') == 0
            True
            sage: pari('O(2)') == 0
            True
        """
        cdef bint r
        cdef GEN x = (<gen>left).g
        cdef GEN y = (<gen>right).g
        sig_on()
        if op == 2:    # ==
            r = (gequal(x, y) != 0)
        elif op == 3:  # !=
            r = (gequal(x, y) == 0)
        else:
            r = rich_to_bool(op, gcmp(x, y))
        sig_off()
        return r

    cpdef int _cmp_(left, right) except -2:
        """
        Compare ``left`` and ``right``.

        This uses PARI's ``cmp_universal()`` routine, which defines
        a total ordering on the set of all PARI objects (up to the
        indistinguishability relation given by ``gidentical()``).

        .. WARNING::

            This comparison is only mathematically meaningful when
            comparing 2 integers. In particular, when comparing
            rationals or reals, this does not correspond to the natural
            ordering.

        EXAMPLES::

            sage: cmp(pari(5), 5)
            0
            sage: cmp(pari(5), 10)
            -1
            sage: cmp(pari(2.5), None)
            1
            sage: cmp(pari(3), pari(3))
            0
            sage: cmp(pari('x^2 + 1'), pari('I-1'))
            1
            sage: cmp(pari(I), pari(I))
            0

        Beware when comparing rationals or reals::

            sage: cmp(pari(2/3), pari(2/5))
            -1
            sage: two = RealField(256)(2)._pari_()
            sage: cmp(two, pari(1.0))
            1
            sage: cmp(two, pari(2.0))
            1
            sage: cmp(two, pari(3.0))
            1

        Since :trac:`17026`, different elements with the same string
        representation can be distinguished by ``cmp()``::

            sage: a = pari(0); a
            0
            sage: b = pari("0*ffgen(ffinit(29, 10))"); b
            0
            sage: cmp(a, b)
            -1

            sage: x = pari("x"); x
            x
            sage: y = pari("ffgen(ffinit(3, 5))"); y
            x
            sage: cmp(x, y)
            1

        """
        cdef int r
        sig_on()
        r = cmp_universal(left.g, (<gen>right).g)
        sig_off()
        return r

    def __copy__(gen self):
        sig_on()
        return P.new_gen(gcopy(self.g))

    def list_str(gen self):
        """
        Return str that might correctly evaluate to a Python-list.

        TESTS::

            sage: pari.primes(5).list_str()
            doctest:...: DeprecationWarning: the method list_str() is deprecated
            See http://trac.sagemath.org/20219 for details.
            [2, 3, 5, 7, 11]
        """
        deprecation(20219, "the method list_str() is deprecated")

        s = str(self)
        if s[:4] == "Mat(":
            s = "[" + s[4:-1] + "]"
        s = s.replace("~","")
        if s.find(";") != -1:
            s = s.replace(";","], [")
            s = "[" + s + "]"
            return eval(s)
        else:
            return eval(s)

    def __hex__(gen self):
        """
        Return the hexadecimal digits of self in lower case.

        EXAMPLES::

            sage: print(hex(pari(0)))
            0
            sage: print(hex(pari(15)))
            f
            sage: print(hex(pari(16)))
            10
            sage: print(hex(pari(16938402384092843092843098243)))
            36bb1e3929d1a8fe2802f083
            sage: print(hex(long(16938402384092843092843098243)))
            0x36bb1e3929d1a8fe2802f083L
            sage: print(hex(pari(-16938402384092843092843098243)))
            -36bb1e3929d1a8fe2802f083
        """
        cdef GEN x
        cdef long lx
        cdef long *xp
        cdef long w
        cdef char *s
        cdef char *sp
        cdef char *hexdigits
        hexdigits = "0123456789abcdef"
        cdef int i, j
        cdef int size
        x = self.g
        if typ(x) != t_INT:
            raise TypeError("gen must be of PARI type t_INT")
        if not signe(x):
            return "0"
        lx = lgefint(x)-2  # number of words
        size = lx*2*sizeof(long)
        s = <char *>sig_malloc(size+2) # 1 char for sign, 1 char for '\0'
        sp = s + size+1
        sp[0] = 0
        xp = int_LSW(x)
        for i from 0 <= i < lx:
            w = xp[0]
            for j from 0 <= j < 2*sizeof(long):
                sp = sp-1
                sp[0] = hexdigits[w & 15]
                w = w>>4
            xp = int_nextW(xp)
        # remove leading zeros!
        while sp[0] == c'0':
            sp = sp+1
        if signe(x) < 0:
            sp = sp-1
            sp[0] = c'-'
        k = <object>sp
        sig_free(s)
        return k

    def __int__(gen self):
        """
        Convert ``self`` to a Python integer.

        If the number is too large to fit into a Pyhon ``int``, a
        Python ``long`` is returned instead.

        EXAMPLES::

            sage: int(pari(0))
            0
            sage: int(pari(10))
            10
            sage: int(pari(-10))
            -10
            sage: int(pari(123456789012345678901234567890))
            123456789012345678901234567890L
            sage: int(pari(-123456789012345678901234567890))
            -123456789012345678901234567890L
            sage: int(pari(2^31-1))
            2147483647
            sage: int(pari(-2^31))
            -2147483648
            sage: int(pari("Pol(10)"))
            10
            sage: int(pari("Mod(2, 7)"))
            2
            sage: int(pari(RealField(63)(2^63-1)))
            9223372036854775807L  # 32-bit
            9223372036854775807   # 64-bit
            sage: int(pari(RealField(63)(2^63+2)))
            9223372036854775810L
        """
        return gen_to_integer(self)

    def python_list_small(gen self):
        """
        Return a Python list of the PARI gens. This object must be of type
        t_VECSMALL, and the resulting list contains python 'int's.

        EXAMPLES::

            sage: v=pari([1,2,3,10,102,10]).Vecsmall()
            sage: w = v.python_list_small()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type 'int'>
        """
        cdef long n
        if typ(self.g) != t_VECSMALL:
            raise TypeError("Object (=%s) must be of type t_VECSMALL." % self)
        return [self.g[n+1] for n in range(glength(self.g))]

    def python_list(gen self):
        """
        Return a Python list of the PARI gens. This object must be of type
        t_VEC or t_COL.

        INPUT: None

        OUTPUT:

        -  ``list`` - Python list whose elements are the
           elements of the input gen.


        EXAMPLES::

            sage: v = pari([1,2,3,10,102,10])
            sage: w = v.python_list()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type 'sage.libs.pari.gen.gen'>
            sage: pari("[1,2,3]").python_list()
            [1, 2, 3]

            sage: pari("[1,2,3]~").python_list()
            [1, 2, 3]
        """
        cdef long n
        cdef gen t

        if typ(self.g) != t_VEC and typ(self.g) != t_COL:
            raise TypeError("Object (=%s) must be of type t_VEC or t_COL." % self)
        return [self[n] for n in range(glength(self.g))]

    def python(self, locals=None):
        """
        Return the closest Python/Sage equivalent of the given PARI object.

        INPUT:

        - `z` -- PARI ``gen``

        - `locals` -- optional dictionary used in fallback cases that
          involve :func:`sage_eval`

        .. NOTE::

            If ``self`` is a real (type ``t_REAL``), then the result
            will be a RealField element of the equivalent precision;
            if ``self`` is a complex (type ``t_COMPLEX``), then the
            result will be a ComplexField element of precision the
            maximal precision of the real and imaginary parts.

        EXAMPLES::

            sage: pari('389/17').python()
            389/17
            sage: f = pari('(2/3)*x^3 + x - 5/7 + y'); f
            2/3*x^3 + x + (y - 5/7)
            sage: var('x,y')
            (x, y)
            sage: f.python({'x':x, 'y':y})
            2/3*x^3 + x + y - 5/7

        You can also use :meth:`.sage`, which is an alias::

            sage: f.sage({'x':x, 'y':y})
            2/3*x^3 + x + y - 5/7

        Converting a real number::

            sage: pari.set_real_precision(70)
            15
            sage: a = pari('1.234').python(); a
            1.234000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a.parent()
            Real Field with 256 bits of precision
            sage: pari.set_real_precision(15)
            70
            sage: a = pari('1.234').python(); a
            1.23400000000000000
            sage: a.parent()
            Real Field with 64 bits of precision

        For complex numbers, the parent depends on the PARI type::

            sage: a = pari('(3+I)').python(); a
            i + 3
            sage: a.parent()
            Number Field in i with defining polynomial x^2 + 1

            sage: a = pari('2^31-1').python(); a
            2147483647
            sage: a.parent()
            Integer Ring

            sage: a = pari('12/34').python(); a
            6/17
            sage: a.parent()
            Rational Field

            sage: a = pari('(3+I)/2').python(); a
            1/2*i + 3/2
            sage: a.parent()
            Number Field in i with defining polynomial x^2 + 1

            sage: z = pari(CC(1.0+2.0*I)); z
            1.00000000000000 + 2.00000000000000*I
            sage: a = z.python(); a
            1.00000000000000000 + 2.00000000000000000*I
            sage: a.parent()
            Complex Field with 64 bits of precision

            sage: I = sqrt(-1)
            sage: a = pari(1.0 + 2.0*I).python(); a
            1.00000000000000000 + 2.00000000000000000*I
            sage: a.parent()
            Complex Field with 64 bits of precision

        Vectors and matrices::

            sage: a = pari('[1,2,3,4]')
            sage: a
            [1, 2, 3, 4]
            sage: a.type()
            't_VEC'
            sage: b = a.python(); b
            [1, 2, 3, 4]
            sage: type(b)
            <type 'list'>

            sage: a = pari('[1,2;3,4]')
            sage: a.type()
            't_MAT'
            sage: b = a.python(); b
            [1 2]
            [3 4]
            sage: b.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

            sage: a = pari('Vecsmall([1,2,3,4])')
            sage: a.type()
            't_VECSMALL'
            sage: a.python()
            [1, 2, 3, 4]

        We use the locals dictionary::

            sage: f = pari('(2/3)*x^3 + x - 5/7 + y')
            sage: x,y=var('x,y')
            sage: from sage.libs.pari.gen import gentoobj
            sage: gentoobj(f, {'x':x, 'y':y})
            2/3*x^3 + x + y - 5/7
            sage: gentoobj(f)
            Traceback (most recent call last):
            ...
            NameError: name 'x' is not defined

        Conversion of p-adics::

            sage: K = Qp(11,5)
            sage: x = K(11^-10 + 5*11^-7 + 11^-6); x
            11^-10 + 5*11^-7 + 11^-6 + O(11^-5)
            sage: y = pari(x); y
            11^-10 + 5*11^-7 + 11^-6 + O(11^-5)
            sage: y.sage()
            11^-10 + 5*11^-7 + 11^-6 + O(11^-5)
            sage: pari(K(11^-5)).sage()
            11^-5 + O(11^0)

        Conversion of infinities::

            sage: pari('oo').sage()
            +Infinity
            sage: pari('-oo').sage()
            -Infinity
        """
        return gentoobj(self, locals)

    sage = _sage_ = _eval_ = python

    def __long__(gen self):
        """
        Convert ``self`` to a Python ``long``.

        EXAMPLES::

            sage: long(pari(0))
            0L
            sage: long(pari(10))
            10L
            sage: long(pari(-10))
            -10L
            sage: long(pari(123456789012345678901234567890))
            123456789012345678901234567890L
            sage: long(pari(-123456789012345678901234567890))
            -123456789012345678901234567890L
            sage: long(pari(2^31-1))
            2147483647L
            sage: long(pari(-2^31))
            -2147483648L
            sage: long(pari("Pol(10)"))
            10L
            sage: long(pari("Mod(2, 7)"))
            2L
        """
        x = gen_to_integer(self)
        if isinstance(x, long):
            return x
        else:
            return long(x)

    def __float__(gen self):
        """
        Return Python float.
        """
        cdef double d
        sig_on()
        d = gtodouble(self.g)
        sig_off()
        return d

    def __complex__(self):
        r"""
        Return ``self`` as a Python ``complex``
        value.

        EXAMPLES::

            sage: g = pari(-1.0)^(1/5); g
            0.809016994374947 + 0.587785252292473*I
            sage: g.__complex__()
            (0.8090169943749475+0.5877852522924731j)
            sage: complex(g)
            (0.8090169943749475+0.5877852522924731j)

        ::

            sage: g = pari(Integers(5)(3)); g
            Mod(3, 5)
            sage: complex(g)
            Traceback (most recent call last):
            ...
            PariError: incorrect type in greal/gimag (t_INTMOD)
        """
        cdef double re, im
        sig_on()
        re = gtodouble(greal(self.g))
        im = gtodouble(gimag(self.g))
        sig_off()
        return complex(re, im)

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: pari('1').__nonzero__()
            True
            sage: pari('x').__nonzero__()
            True
            sage: bool(pari(0))
            False
            sage: a = pari('Mod(0,3)')
            sage: a.__nonzero__()
            False
        """
        return not gequal0(self.g)

    def gequal(gen a, b):
        r"""
        Check whether `a` and `b` are equal using PARI's ``gequal``.

        EXAMPLES::

            sage: a = pari(1); b = pari(1.0); c = pari('"some_string"')
            sage: a.gequal(a)
            True
            sage: b.gequal(b)
            True
            sage: c.gequal(c)
            True
            sage: a.gequal(b)
            True
            sage: a.gequal(c)
            False

        WARNING: this relation is not transitive::

            sage: a = pari('[0]'); b = pari(0); c = pari('[0,0]')
            sage: a.gequal(b)
            True
            sage: b.gequal(c)
            True
            sage: a.gequal(c)
            False
        """
        cdef gen t0 = objtogen(b)
        sig_on()
        cdef int ret = gequal(a.g, t0.g)
        sig_off()
        return ret != 0

    def gequal0(gen a):
        r"""
        Check whether `a` is equal to zero.

        EXAMPLES::

            sage: pari(0).gequal0()
            True
            sage: pari(1).gequal0()
            False
            sage: pari(1e-100).gequal0()
            False
            sage: pari("0.0 + 0.0*I").gequal0()
            True
            sage: pari(GF(3^20,'t')(0)).gequal0()
            True
        """
        sig_on()
        cdef int ret = gequal0(a.g)
        sig_off()
        return ret != 0

    def gequal_long(gen a, long b):
        r"""
        Check whether `a` is equal to the ``long int`` `b` using PARI's ``gequalsg``.

        EXAMPLES::

            sage: a = pari(1); b = pari(2.0); c = pari('3*matid(3)')
            sage: a.gequal_long(1)
            True
            sage: a.gequal_long(-1)
            False
            sage: a.gequal_long(0)
            False
            sage: b.gequal_long(2)
            True
            sage: b.gequal_long(-2)
            False
            sage: c.gequal_long(3)
            True
            sage: c.gequal_long(-3)
            False
        """
        sig_on()
        cdef int ret = gequalsg(b, a.g)
        sig_off()
        return ret != 0

    def isprime(gen self, long flag=0):
        """
        isprime(x, flag=0): Returns True if x is a PROVEN prime number, and
        False otherwise.

        INPUT:


        -  ``flag`` - int 0 (default): use a combination of
           algorithms. 1: certify primality using the Pocklington-Lehmer Test.
           2: certify primality using the APRCL test.


        OUTPUT:


        -  ``bool`` - True or False


        EXAMPLES::

            sage: pari(9).isprime()
            False
            sage: pari(17).isprime()
            True
            sage: n = pari(561)    # smallest Carmichael number
            sage: n.isprime()      # not just a pseudo-primality test!
            False
            sage: n.isprime(1)
            False
            sage: n.isprime(2)
            False
            sage: n = pari(2^31-1)
            sage: n.isprime(1)
            (True, [2, 3, 1; 3, 5, 1; 7, 3, 1; 11, 3, 1; 31, 2, 1; 151, 3, 1; 331, 3, 1])
        """
        cdef GEN x
        sig_on()
        x = gisprime(self.g, flag)
        if typ(x) != t_INT:
            # case flag=1 with prime input: x is the certificate
            return True, P.new_gen(x)
        else:
            sig_off()
            return signe(x) != 0

    def ispseudoprime(gen self, long flag=0):
        """
        ispseudoprime(x, flag=0): Returns True if x is a pseudo-prime
        number, and False otherwise.

        INPUT:


        -  ``flag`` - int 0 (default): checks whether x is a
           Baillie-Pomerance-Selfridge-Wagstaff pseudo prime (strong
           Rabin-Miller pseudo prime for base 2, followed by strong Lucas test
           for the sequence (P,-1), P smallest positive integer such that
           `P^2 - 4` is not a square mod x). 0: checks whether x is a
           strong Miller-Rabin pseudo prime for flag randomly chosen bases
           (with end-matching to catch square roots of -1).


        OUTPUT:


        -  ``bool`` - True or False, or when flag=1, either False or a tuple
           (True, cert) where ``cert`` is a primality certificate.


        EXAMPLES::

            sage: pari(9).ispseudoprime()
            False
            sage: pari(17).ispseudoprime()
            True
            sage: n = pari(561)     # smallest Carmichael number
            sage: n.ispseudoprime(2)
            False
        """
        sig_on()
        cdef long t = ispseudoprime(self.g, flag)
        sig_off()
        return t != 0

    def ispower(gen self, k=None):
        r"""
        Determine whether or not self is a perfect k-th power. If k is not
        specified, find the largest k so that self is a k-th power.

        INPUT:


        -  ``k`` - int (optional)


        OUTPUT:


        -  ``power`` - int, what power it is

        -  ``g`` - what it is a power of


        EXAMPLES::

            sage: pari(9).ispower()
            (2, 3)
            sage: pari(17).ispower()
            (1, 17)
            sage: pari(17).ispower(2)
            (False, None)
            sage: pari(17).ispower(1)
            (1, 17)
            sage: pari(2).ispower()
            (1, 2)
        """
        cdef int n
        cdef GEN x
        cdef gen t0

        if k is None:
            sig_on()
            n = gisanypower(self.g, &x)
            if n == 0:
                sig_off()
                return 1, self
            else:
                return n, P.new_gen(x)
        else:
            t0 = objtogen(k)
            sig_on()
            n = ispower(self.g, t0.g, &x)
            if n == 0:
                sig_off()
                return False, None
            else:
                return k, P.new_gen(x)

    def isprimepower(gen self):
        r"""
        Check whether ``self`` is a prime power (with an exponent >= 1).

        INPUT:

        - ``self`` - A PARI integer

        OUTPUT:

        A tuple ``(k, p)`` where `k` is a Python integer and `p` a PARI
        integer.

        - If the input was a prime power, `p` is the prime and `k` the
          power.
        - Otherwise, `k = 0` and `p` is ``self``.

        .. SEEALSO::

            If you don't need a proof that `p` is prime, you can use
            :meth:`ispseudoprimepower` instead.

        EXAMPLES::

            sage: pari(9).isprimepower()
            (2, 3)
            sage: pari(17).isprimepower()
            (1, 17)
            sage: pari(18).isprimepower()
            (0, 18)
            sage: pari(3^12345).isprimepower()
            (12345, 3)
        """
        cdef GEN x
        cdef long n

        sig_on()
        n = isprimepower(self.g, &x)
        if n == 0:
            sig_off()
            return 0, self
        else:
            return n, P.new_gen(x)

    def ispseudoprimepower(gen self):
        r"""
        Check whether ``self`` is the power (with an exponent >= 1) of
        a pseudo-prime.

        INPUT:

        - ``self`` - A PARI integer

        OUTPUT:

        A tuple ``(k, p)`` where `k` is a Python integer and `p` a PARI
        integer.

        - If the input was a pseudoprime power, `p` is the pseudoprime
          and `k` the power.
        - Otherwise, `k = 0` and `p` is ``self``.

        EXAMPLES::

            sage: pari(3^12345).ispseudoprimepower()
            (12345, 3)
            sage: p = pari(2^1500 + 1465)         # next_prime(2^1500)
            sage: (p^11).ispseudoprimepower()[0]  # very fast
            11
        """
        cdef GEN x
        cdef long n

        sig_on()
        n = ispseudoprimepower(self.g, &x)
        if n == 0:
            sig_off()
            return 0, self
        else:
            return n, P.new_gen(x)

    def vecmax(x):
        """
        Return the maximum of the elements of the vector/matrix `x`.

        EXAMPLES::

            sage: pari([1, -5/3, 8.0]).vecmax()
            8.00000000000000
        """
        sig_on()
        return P.new_gen(vecmax(x.g))

    def vecmin(x):
        """
        Return the minimum of the elements of the vector/matrix `x`.

        EXAMPLES::

            sage: pari([1, -5/3, 8.0]).vecmin()
            -5/3
        """
        sig_on()
        return P.new_gen(vecmin(x.g))

    def Col(gen x, long n = 0):
        """
        Transform the object `x` into a column vector with minimal size `|n|`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the column vector of minimal length `|n|`. If `n > 0`,
          append zeros; if `n < 0`, prepend zeros.

        OUTPUT:

        A PARI column vector (type ``t_COL``)

        EXAMPLES::

            sage: pari(1.5).Col()
            [1.50000000000000]~
            sage: pari([1,2,3,4]).Col()
            [1, 2, 3, 4]~
            sage: pari('[1,2; 3,4]').Col()
            [[1, 2], [3, 4]]~
            sage: pari('"Sage"').Col()
            ["S", "a", "g", "e"]~
            sage: pari('x + 3*x^3').Col()
            [3, 0, 1, 0]~
            sage: pari('x + 3*x^3 + O(x^5)').Col()
            [1, 0, 3, 0]~

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Col(2)
            [1, 2, 3, 4]~
            sage: pari([1,2,3,4]).Col(-2)
            [1, 2, 3, 4]~
            sage: pari([1,2,3,4]).Col(6)
            [1, 2, 3, 4, 0, 0]~
            sage: pari([1,2,3,4]).Col(-6)
            [0, 0, 1, 2, 3, 4]~

        See also :meth:`Vec` (create a row vector) for more examples
        and :meth:`Colrev` (create a column in reversed order).
        """
        sig_on()
        return P.new_gen(_Vec_append(gtocol(x.g), gen_0, n))

    def Colrev(gen x, long n = 0):
        """
        Transform the object `x` into a column vector with minimal size `|n|`.
        The order of the resulting vector is reversed compared to :meth:`Col`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          prepend zeros; if `n < 0`, append zeros.

        OUTPUT:

        A PARI column vector (type ``t_COL``)

        EXAMPLES::

            sage: pari(1.5).Colrev()
            [1.50000000000000]~
            sage: pari([1,2,3,4]).Colrev()
            [4, 3, 2, 1]~
            sage: pari('[1,2; 3,4]').Colrev()
            [[3, 4], [1, 2]]~
            sage: pari('x + 3*x^3').Colrev()
            [0, 1, 0, 3]~

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Colrev(2)
            [4, 3, 2, 1]~
            sage: pari([1,2,3,4]).Colrev(-2)
            [4, 3, 2, 1]~
            sage: pari([1,2,3,4]).Colrev(6)
            [0, 0, 4, 3, 2, 1]~
            sage: pari([1,2,3,4]).Colrev(-6)
            [4, 3, 2, 1, 0, 0]~
        """
        sig_on()
        # Create a non-reversed column vector
        cdef GEN v = _Vec_append(gtocol(x.g), gen_0, n)
        # Reverse it in-place
        cdef GEN L = v + 1
        cdef GEN R = v + (lg(v)-1)
        cdef long t
        while (L < R):
            t = L[0]
            L[0] = R[0]
            R[0] = t
            L += 1
            R -= 1
        return P.new_gen(v)

    def Ser(gen f, v=-1, long precision=-1):
        """
        Return a power series or Laurent series in the variable `v`
        constructed from the object `f`.

        INPUT:

        - ``f`` -- PARI gen

        - ``v`` -- PARI variable (default: `x`)

        - ``precision`` -- the desired relative precision (default:
          the value returned by ``pari.get_series_precision()``).
          This is the absolute precision minus the `v`-adic valuation.

        OUTPUT:

        - PARI object of type ``t_SER``

        The series is constructed from `f` in the following way:

        - If `f` is a scalar, a constant power series is returned.

        - If `f` is a polynomial, it is converted into a power series
          in the obvious way.

        - If `f` is a rational function, it will be expanded in a
          Laurent series around `v = 0`.

        - If `f` is a vector, its coefficients become the coefficients
          of the power series, starting from the constant term.  This
          is the convention used by the function ``Polrev()``, and the
          reverse of that used by ``Pol()``.

        .. warning::

           This function will not transform objects containing
           variables of higher priority than `v`.

        EXAMPLES::

            sage: pari(2).Ser()
            2 + O(x^16)
            sage: pari(Mod(0, 7)).Ser()
            Mod(0, 7)*x^15 + O(x^16)

            sage: x = pari([1, 2, 3, 4, 5])
            sage: x.Ser()
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + O(x^16)
            sage: f = x.Ser('v'); print(f)
            1 + 2*v + 3*v^2 + 4*v^3 + 5*v^4 + O(v^16)
            sage: pari(1)/f
            1 - 2*v + v^2 + 6*v^5 - 17*v^6 + 16*v^7 - 5*v^8 + 36*v^10 - 132*v^11 + 181*v^12 - 110*v^13 + 25*v^14 + 216*v^15 + O(v^16)

            sage: pari('x^5').Ser(precision=20)
            x^5 + O(x^25)
            sage: pari('1/x').Ser(precision=1)
            x^-1 + O(x^0)

        """
        if precision < 0:
            precision = P.get_series_precision()
        sig_on()
        cdef long vn = P.get_var(v)
        if typ(f.g) == t_VEC:
            # The precision flag is ignored for vectors, so we first
            # convert the vector to a polynomial.
            return P.new_gen(gtoser(gtopolyrev(f.g, vn), vn, precision))
        else:
            return P.new_gen(gtoser(f.g, vn, precision))

    def Str(self):
        """
        Str(self): Return the print representation of self as a PARI
        object.

        INPUT:


        -  ``self`` - gen


        OUTPUT:


        -  ``gen`` - a PARI gen of type t_STR, i.e., a PARI
           string


        EXAMPLES::

            sage: pari([1,2,['abc',1]]).Str()
            "[1, 2, [abc, 1]]"
            sage: pari([1,1, 1.54]).Str()
            "[1, 1, 1.54000000000000]"
            sage: pari(1).Str()       # 1 is automatically converted to string rep
            "1"
            sage: x = pari('x')       # PARI variable "x"
            sage: x.Str()             # is converted to string rep.
            "x"
            sage: x.Str().type()
            't_STR'
        """
        cdef char* c
        sig_on()
        # Use sig_block(), which is needed because GENtostr() uses
        # malloc(), which is dangerous inside sig_on()
        sig_block()
        c = GENtostr(self.g)
        sig_unblock()
        v = P.new_gen(strtoGENstr(c))
        pari_free(c)
        return v

    def Strexpand(gen x):
        """
        Concatenate the entries of the vector `x` into a single string,
        then perform tilde expansion and environment variable expansion
        similar to shells.

        INPUT:

        - ``x`` -- PARI gen. Either a vector or an element which is then
          treated like `[x]`.

        OUTPUT:

        - PARI string (type ``t_STR``)

        EXAMPLES::

            sage: pari('"~/subdir"').Strexpand()     # random
            "/home/johndoe/subdir"
            sage: pari('"$SAGE_LOCAL"').Strexpand()  # random
            "/usr/local/sage/local"

        TESTS::

            sage: a = pari('"$HOME"')
            sage: a.Strexpand() != a
            True
        """
        if typ(x.g) != t_VEC:
            x = P.vector(1, [x])
        sig_on()
        return P.new_gen(Strexpand(x.g))

    def Strtex(gen x):
        r"""
        Strtex(x): Translates the vector x of PARI gens to TeX format and
        returns the resulting concatenated strings as a PARI t_STR.

        INPUT:

        - ``x`` -- PARI gen. Either a vector or an element which is then
          treated like `[x]`.

        OUTPUT:

        - PARI string (type ``t_STR``)

        EXAMPLES::

            sage: v=pari('x^2')
            sage: v.Strtex()
            "x^2"
            sage: v=pari(['1/x^2','x'])
            sage: v.Strtex()
            "\\frac{1}{x^2}x"
            sage: v=pari(['1 + 1/x + 1/(y+1)','x-1'])
            sage: v.Strtex()
            "\\frac{ \\left(y\n + 2\\right) \\*x\n + \\left(y\n + 1\\right) }{ \\left(y\n + 1\\right) \\*x}x\n - 1"
        """
        if typ(x.g) != t_VEC:
            x = P.vector(1, [x])
        sig_on()
        return P.new_gen(Strtex(x.g))

    printtex = deprecated_function_alias(20219, Strtex)

    def Vec(gen x, long n = 0):
        """
        Transform the object `x` into a vector with minimal size `|n|`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          append zeros; if `n < 0`, prepend zeros.

        OUTPUT:

        A PARI vector (type ``t_VEC``)

        EXAMPLES::

            sage: pari(1).Vec()
            [1]
            sage: pari('x^3').Vec()
            [1, 0, 0, 0]
            sage: pari('x^3 + 3*x - 2').Vec()
            [1, 0, 3, -2]
            sage: pari([1,2,3]).Vec()
            [1, 2, 3]
            sage: pari('[1, 2; 3, 4]').Vec()
            [[1, 3]~, [2, 4]~]
            sage: pari('"Sage"').Vec()
            ["S", "a", "g", "e"]
            sage: pari('2*x^2 + 3*x^3 + O(x^5)').Vec()
            [2, 3, 0]
            sage: pari('2*x^-2 + 3*x^3 + O(x^5)').Vec()
            [2, 0, 0, 0, 0, 3, 0]

        Note the different term ordering for polynomials and series::

            sage: pari('1 + x + 3*x^3 + O(x^5)').Vec()
            [1, 1, 0, 3, 0]
            sage: pari('1 + x + 3*x^3').Vec()
            [3, 0, 1, 1]

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Vec(2)
            [1, 2, 3, 4]
            sage: pari([1,2,3,4]).Vec(-2)
            [1, 2, 3, 4]
            sage: pari([1,2,3,4]).Vec(6)
            [1, 2, 3, 4, 0, 0]
            sage: pari([1,2,3,4]).Vec(-6)
            [0, 0, 1, 2, 3, 4]

        See also :meth:`Col` (create a column vector) and :meth:`Vecrev`
        (create a vector in reversed order).
        """
        sig_on()
        return P.new_gen(_Vec_append(gtovec(x.g), gen_0, n))

    def Vecrev(gen x, long n = 0):
        """
        Transform the object `x` into a vector with minimal size `|n|`.
        The order of the resulting vector is reversed compared to :meth:`Vec`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          prepend zeros; if `n < 0`, append zeros.

        OUTPUT:

        A PARI vector (type ``t_VEC``)

        EXAMPLES::

            sage: pari(1).Vecrev()
            [1]
            sage: pari('x^3').Vecrev()
            [0, 0, 0, 1]
            sage: pari('x^3 + 3*x - 2').Vecrev()
            [-2, 3, 0, 1]
            sage: pari([1, 2, 3]).Vecrev()
            [3, 2, 1]
            sage: pari('Col([1, 2, 3])').Vecrev()
            [3, 2, 1]
            sage: pari('[1, 2; 3, 4]').Vecrev()
            [[2, 4]~, [1, 3]~]
            sage: pari('"Sage"').Vecrev()
            ["e", "g", "a", "S"]

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Vecrev(2)
            [4, 3, 2, 1]
            sage: pari([1,2,3,4]).Vecrev(-2)
            [4, 3, 2, 1]
            sage: pari([1,2,3,4]).Vecrev(6)
            [0, 0, 4, 3, 2, 1]
            sage: pari([1,2,3,4]).Vecrev(-6)
            [4, 3, 2, 1, 0, 0]
        """
        sig_on()
        return P.new_gen(_Vec_append(gtovecrev(x.g), gen_0, -n))

    def Vecsmall(gen x, long n = 0):
        """
        Transform the object `x` into a ``t_VECSMALL`` with minimal size `|n|`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          append zeros; if `n < 0`, prepend zeros.

        OUTPUT:

        A PARI vector of small integers (type ``t_VECSMALL``)

        EXAMPLES::

            sage: pari([1,2,3]).Vecsmall()
            Vecsmall([1, 2, 3])
            sage: pari('"Sage"').Vecsmall()
            Vecsmall([83, 97, 103, 101])
            sage: pari(1234).Vecsmall()
            Vecsmall([1234])
            sage: pari('x^2 + 2*x + 3').Vecsmall()
            Vecsmall([1, 2, 3])

        We demonstate the `n` argument::

            sage: pari([1,2,3]).Vecsmall(2)
            Vecsmall([1, 2, 3])
            sage: pari([1,2,3]).Vecsmall(-2)
            Vecsmall([1, 2, 3])
            sage: pari([1,2,3]).Vecsmall(6)
            Vecsmall([1, 2, 3, 0, 0, 0])
            sage: pari([1,2,3]).Vecsmall(-6)
            Vecsmall([0, 0, 0, 1, 2, 3])
        """
        sig_on()
        return P.new_gen(_Vec_append(gtovecsmall(x.g), <GEN>0, n))


    def bittest(gen x, long n):
        """
        bittest(x, long n): Returns bit number n (coefficient of
        `2^n` in binary) of the integer x. Negative numbers behave
        as if modulo a big power of 2.

        INPUT:


        -  ``x`` - gen (pari integer)


        OUTPUT:


        -  ``bool`` - a Python bool


        EXAMPLES::

            sage: x = pari(6)
            sage: x.bittest(0)
            False
            sage: x.bittest(1)
            True
            sage: x.bittest(2)
            True
            sage: x.bittest(3)
            False
            sage: pari(-3).bittest(0)
            True
            sage: pari(-3).bittest(1)
            False
            sage: [pari(-3).bittest(n) for n in range(10)]
            [True, False, True, True, True, True, True, True, True, True]
        """
        sig_on()
        cdef long b = bittest(x.g, n)
        sig_off()
        return b != 0

    lift_centered = gen_auto.centerlift

    def padicprime(gen x):
        """
        The uniformizer of the p-adic ring this element lies in, as a t_INT.

        INPUT:

        - ``x`` - gen, of type t_PADIC

        OUTPUT:

        - ``p`` - gen, of type t_INT

        EXAMPLES::

            sage: K = Qp(11,5)
            sage: x = K(11^-10 + 5*11^-7 + 11^-6)
            sage: y = pari(x)
            sage: y.padicprime()
            11
            sage: y.padicprime().type()
            't_INT'
        """
        sig_on()
        return P.new_gen(gel(x.g, 2))

    def precision(gen x, long n=-1):
        """
        Change the precision of `x` to be `n`, where `n` is an integer.
        If `n` is omitted, output the real precision of `x`.

        INPUT:

        -  ``x`` - gen

        -  ``n`` - (optional) int

        OUTPUT: gen
        """
        if n <= 0:
            return precision(x.g)
        sig_on()
        return P.new_gen(precision0(x.g, n))

    def round(gen x, estimate=False):
        """
        round(x,estimate=False): If x is a real number, returns x rounded
        to the nearest integer (rounding up). If the optional argument
        estimate is True, also returns the binary exponent e of the
        difference between the original and the rounded value (the
        "fractional part") (this is the integer ceiling of log_2(error)).

        When x is a general PARI object, this function returns the result
        of rounding every coefficient at every level of PARI object. Note
        that this is different than what the truncate function does (see
        the example below).

        One use of round is to get exact results after a long approximate
        computation, when theory tells you that the coefficients must be
        integers.

        INPUT:


        -  ``x`` - gen

        -  ``estimate`` - (optional) bool, False by default


        OUTPUT:

        - if estimate is False, return a single gen.

        - if estimate is True, return rounded version of x and error
          estimate in bits, both as gens.

        EXAMPLES::

            sage: pari('1.5').round()
            2
            sage: pari('1.5').round(True)
            (2, -1)
            sage: pari('1.5 + 2.1*I').round()
            2 + 2*I
            sage: pari('1.0001').round(True)
            (1, -14)
            sage: pari('(2.4*x^2 - 1.7)/x').round()
            (2*x^2 - 2)/x
            sage: pari('(2.4*x^2 - 1.7)/x').truncate()
            2.40000000000000*x
        """
        cdef int n
        cdef long e
        cdef gen y
        sig_on()
        if not estimate:
            return P.new_gen(ground(x.g))
        y = P.new_gen(grndtoi(x.g, &e))
        return y, e

    def sizeword(gen x):
        """
        Return the total number of machine words occupied by the
        complete tree of the object x.  A machine word is 32 or
        64 bits, depending on the computer.

        INPUT:

        -  ``x`` - gen

        OUTPUT: int (a Python int)

        EXAMPLES::

            sage: pari('0').sizeword()
            2
            sage: pari('1').sizeword()
            3
            sage: pari('1000000').sizeword()
            3
            sage: pari('10^100').sizeword()
            13      # 32-bit
            8       # 64-bit
            sage: pari(RDF(1.0)).sizeword()
            4       # 32-bit
            3       # 64-bit
            sage: pari('x').sizeword()
            9
            sage: pari('x^20').sizeword()
            66
            sage: pari('[x, I]').sizeword()
            20
        """
        return gsizeword(x.g)

    def sizebyte(gen x):
        """
        Return the total number of bytes occupied by the complete tree
        of the object x. Note that this number depends on whether the
        computer is 32-bit or 64-bit.

        INPUT:

        -  ``x`` - gen

        OUTPUT: int (a Python int)

        EXAMPLE::

            sage: pari('1').sizebyte()
            12           # 32-bit
            24           # 64-bit
        """
        return gsizebyte(x.g)

    def truncate(gen x, estimate=False):
        """
        truncate(x,estimate=False): Return the truncation of x. If estimate
        is True, also return the number of error bits.

        When x is in the real numbers, this means that the part after the
        decimal point is chopped away, e is the binary exponent of the
        difference between the original and truncated value (the
        "fractional part"). If x is a rational function, the result is the
        integer part (Euclidean quotient of numerator by denominator) and
        if requested the error estimate is 0.

        When truncate is applied to a power series (in X), it transforms it
        into a polynomial or a rational function with denominator a power
        of X, by chopping away the `O(X^k)`. Similarly, when
        applied to a p-adic number, it transforms it into an integer or a
        rational number by chopping away the `O(p^k)`.

        INPUT:


        -  ``x`` - gen

        -  ``estimate`` - (optional) bool, which is False by
           default


        OUTPUT:

        - if estimate is False, return a single gen.

        - if estimate is True, return rounded version of x and error
          estimate in bits, both as gens.

        EXAMPLES::

            sage: pari('(x^2+1)/x').round()
            (x^2 + 1)/x
            sage: pari('(x^2+1)/x').truncate()
            x
            sage: pari('1.043').truncate()
            1
            sage: pari('1.043').truncate(True)
            (1, -5)
            sage: pari('1.6').truncate()
            1
            sage: pari('1.6').round()
            2
            sage: pari('1/3 + 2 + 3^2 + O(3^3)').truncate()
            34/3
            sage: pari('sin(x+O(x^10))').truncate()
            1/362880*x^9 - 1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x
            sage: pari('sin(x+O(x^10))').round()   # each coefficient has abs < 1
            x + O(x^10)
        """
        cdef long e
        cdef gen y
        sig_on()
        if not estimate:
            return P.new_gen(gtrunc(x.g))
        y = P.new_gen(gcvtoi(x.g, &e))
        return y, e

    def _valp(gen x):
        """
        Return the valuation of x where x is a p-adic number (t_PADIC)
        or a Laurent series (t_SER).  If x is a different type, this
        will give a bogus number.

        EXAMPLES::

            sage: pari('1/x^2 + O(x^10)')._valp()
            -2
            sage: pari('O(x^10)')._valp()
            10
            sage: pari('(1145234796 + O(3^10))/771966234')._valp()
            -2
            sage: pari('O(2^10)')._valp()
            10
            sage: pari('x')._valp()   # random
            -35184372088832
        """
        # This is a simple macro, so we don't need sig_on()
        return valp(x.g)

    def bernfrac(x):
        r"""
        The Bernoulli number `B_x`, where `B_0 = 1`,
        `B_1 = -1/2`, `B_2 = 1/6,\ldots,` expressed as a
        rational number. The argument `x` should be of type
        integer.

        EXAMPLES::

            sage: pari(18).bernfrac()
            43867/798
            sage: [pari(n).bernfrac() for n in range(10)]
            [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0]
        """
        return P.bernfrac(x)

    def bernreal(x, unsigned long precision=0):
        r"""
        The Bernoulli number `B_x`, as for the function bernfrac,
        but `B_x` is returned as a real number (with the current
        precision).

        EXAMPLES::

            sage: pari(18).bernreal()
            54.9711779448622
            sage: pari(18).bernreal(precision=192).sage()
            54.9711779448621553884711779448621553884711779448621553885
        """
        return P.bernreal(x, precision)

    def besselk(gen nu, x, flag=None, unsigned long precision=0):
        """
        nu.besselk(x): K-Bessel function (modified Bessel function
        of the second kind) of index nu, which can be complex, and argument
        x.

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        INPUT:


        -  ``nu`` - a complex number

        -  ``x`` - real number (positive or negative)

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(2+i).besselk(3)
            0.0455907718407551 + 0.0289192946582081*I

        ::

            sage: pari(2+i).besselk(-3)
            -4.34870874986752 - 5.38744882697109*I

        ::

            sage: pari(2+i).besselk(300)
            3.74224603319728 E-132 + 2.49071062641525 E-134*I
            sage: pari(2+i).besselk(300, flag=1)
            doctest:...: DeprecationWarning: The flag argument to besselk() is deprecated and not used anymore
            See http://trac.sagemath.org/20219 for details.
            3.74224603319728 E-132 + 2.49071062641525 E-134*I
        """
        if flag is not None:
            deprecation(20219, 'The flag argument to besselk() is deprecated and not used anymore')
        cdef gen t0 = objtogen(x)
        sig_on()
        return P.new_gen(kbessel(nu.g, t0.g, prec_bits_to_words(precision)))

    def eint1(gen x, long n=0, unsigned long precision=0):
        r"""
        x.eint1(n): exponential integral E1(x):

        .. math::

                         \int_{x}^{\infty} \frac{e^{-t}}{t} dt


        If n is present, output the vector [eint1(x), eint1(2\*x), ...,
        eint1(n\*x)]. This is faster than repeatedly calling eint1(i\*x).

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        REFERENCE:

        - See page 262, Prop 5.6.12, of Cohen's book "A Course in
          Computational Algebraic Number Theory".

        EXAMPLES:
        """
        sig_on()
        if n <= 0:
            return P.new_gen(eint1(x.g, prec_bits_to_words(precision)))
        else:
            return P.new_gen(veceint1(x.g, stoi(n), prec_bits_to_words(precision)))

    log_gamma = gen_auto.lngamma

    def polylog(gen x, long m, long flag=0, unsigned long precision=0):
        """
        x.polylog(m,flag=0): m-th polylogarithm of x. flag is optional, and
        can be 0: default, 1: D_m -modified m-th polylog of x, 2:
        D_m-modified m-th polylog of x, 3: P_m-modified m-th polylog of
        x.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        TODO: Add more explanation, copied from the PARI manual.

        EXAMPLES::

            sage: pari(10).polylog(3)
            5.64181141475134 - 8.32820207698027*I
            sage: pari(10).polylog(3,0)
            5.64181141475134 - 8.32820207698027*I
            sage: pari(10).polylog(3,1)
            0.523778453502411
            sage: pari(10).polylog(3,2)
            -0.400459056163451
        """
        sig_on()
        return P.new_gen(polylog0(m, x.g, flag, prec_bits_to_words(precision)))

    def sqrtn(gen x, n, unsigned long precision=0):
        r"""
        x.sqrtn(n): return the principal branch of the n-th root of x,
        i.e., the one such that
        `\arg(\sqrt(x)) \in ]-\pi/n, \pi/n]`. Also returns a second
        argument which is a suitable root of unity allowing one to recover
        all the other roots. If it was not possible to find such a number,
        then this second return value is 0. If the argument is present and
        no square root exists, return 0 instead of raising an error.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        .. note::

           intmods (modulo a prime) and `p`-adic numbers are
           allowed as arguments.

        INPUT:


        -  ``x`` - gen

        -  ``n`` - integer


        OUTPUT:


        -  ``gen`` - principal n-th root of x

        -  ``gen`` - root of unity z that gives the other
           roots


        EXAMPLES::

            sage: s, z = pari(2).sqrtn(5)
            sage: z
            0.309016994374947 + 0.951056516295154*I
            sage: s
            1.14869835499704
            sage: s^5
            2.00000000000000
            sage: z^5
            1.00000000000000 - 2.710505431 E-20*I       # 32-bit
            1.00000000000000 - 2.71050543121376 E-20*I  # 64-bit
            sage: (s*z)^5
            2.00000000000000 + 0.E-19*I
        """
        cdef GEN zetan
        cdef gen t0 = objtogen(n)
        sig_on()
        ans = P.new_gen_noclear(gsqrtn(x.g, t0.g, &zetan, prec_bits_to_words(precision)))
        return ans, P.new_gen(zetan)

    phi = deprecated_function_alias(20219, gen_auto.eulerphi)

    def ffprimroot(self):
        r"""
        Return a primitive root of the multiplicative group of the
        definition field of the given finite field element.

        INPUT:

        - ``self`` -- a PARI finite field element (``FFELT``)

        OUTPUT:

        - A generator of the multiplicative group of the finite field
          generated by ``self``.

        EXAMPLES::

            sage: x = polygen(GF(3))
            sage: k.<a> = GF(9, modulus=x^2+1)
            sage: b = pari(a).ffprimroot()
            sage: b  # random
            a + 1
            sage: b.fforder()
            8
        """
        sig_on()
        return P.new_gen(ffprimroot(self.g, NULL))

    def fibonacci(self):
        r"""
        Return the Fibonacci number of index x.

        EXAMPLES::

            sage: pari(18).fibonacci()
            2584
            sage: [pari(n).fibonacci() for n in range(10)]
            [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
        """
        return P.fibonacci(self)

    def issquare(gen x, find_root=False):
        """
        issquare(x,n): ``True`` if x is a square, ``False`` if not. If
        ``find_root`` is given, also returns the exact square root.
        """
        cdef GEN G
        cdef long t
        cdef gen g
        sig_on()
        if find_root:
            t = itos(gissquareall(x.g, &G))
            if t:
                return True, P.new_gen(G)
            else:
                P.clear_stack()
                return False, None
        else:
            t = itos(gissquare(x.g))
            sig_off()
            return t != 0

    def issquarefree(gen self):
        """
        EXAMPLES::

            sage: pari(10).issquarefree()
            True
            sage: pari(20).issquarefree()
            False
        """
        sig_on()
        cdef long t = issquarefree(self.g)
        sig_off()
        return t != 0

    def sumdiv(gen n):
        """
        Return the sum of the divisors of `n`.

        EXAMPLES::

            sage: pari(10).sumdiv()
            18
        """
        sig_on()
        return P.new_gen(sumdiv(n.g))

    def sumdivk(gen n, long k):
        """
        Return the sum of the k-th powers of the divisors of n.

        EXAMPLES::

            sage: pari(10).sumdivk(2)
            130
        """
        sig_on()
        return P.new_gen(sumdivk(n.g, k))

    def Zn_issquare(gen self, n):
        """
        Return ``True`` if ``self`` is a square modulo `n`, ``False``
        if not.

        INPUT:

        - ``self`` -- integer

        - ``n`` -- integer or factorisation matrix

        EXAMPLES::

            sage: pari(3).Zn_issquare(4)
            False
            sage: pari(4).Zn_issquare(30.factor())
            True

        """
        cdef gen t0 = objtogen(n)
        sig_on()
        cdef long t = Zn_issquare(self.g, t0.g)
        sig_off()
        return t != 0

    def Zn_sqrt(gen self, n):
        """
        Return a square root of ``self`` modulo `n`, if such a square
        root exists; otherwise, raise a ``ValueError``.

        INPUT:

        - ``self`` -- integer

        - ``n`` -- integer or factorisation matrix

        EXAMPLES::

            sage: pari(3).Zn_sqrt(4)
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a square modulo 4
            sage: pari(4).Zn_sqrt(30.factor())
            22

        """
        cdef gen t0 = objtogen(n)
        cdef GEN s
        sig_on()
        s = Zn_sqrt(self.g, t0.g)
        if s == NULL:
            sig_off()
            raise ValueError("%s is not a square modulo %s" % (self, n))
        return P.new_gen(s)

    def ellan(self, long n, python_ints=False):
        """
        Return the first `n` Fourier coefficients of the modular
        form attached to this elliptic curve. See ellak for more details.

        INPUT:


        -  ``n`` - a long integer

        -  ``python_ints`` - bool (default is False); if True,
           return a list of Python ints instead of a PARI gen wrapper.


        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellan(3)
            [1, -2, -1]
            sage: e.ellan(20)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2]
            sage: e.ellan(-1)
            []
            sage: v = e.ellan(10, python_ints=True); v
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2]
            sage: type(v)
            <type 'list'>
            sage: type(v[0])
            <type 'int'>
        """
        sig_on()
        cdef GEN g = anell(self.g, n)
        if python_ints:
            v = [gtolong(gel(g, i+1)) for i in range(glength(g))]
            P.clear_stack()
            return v
        else:
            return P.new_gen(g)

    def ellaplist(self, long n, python_ints=False):
        r"""
        e.ellaplist(n): Returns a PARI list of all the prime-indexed
        coefficients `a_p` (up to n) of the `L`-function
        of the elliptic curve `e`, i.e. the Fourier coefficients of
        the newform attached to `e`.

        INPUT:

        - ``self`` -- an elliptic curve

        - ``n`` -- a long integer

        - ``python_ints`` -- bool (default is False); if True,
          return a list of Python ints instead of a PARI gen wrapper.

        .. WARNING::

            The curve e must be a medium or long vector of the type given by
            ellinit. For this function to work for every n and not just those
            prime to the conductor, e must be a minimal Weierstrass equation.
            If this is not the case, use the function ellminimalmodel first
            before using ellaplist (or you will get INCORRECT RESULTS!)

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: v = e.ellaplist(10); v
            [-2, -1, 1, -2]
            sage: type(v)
            <type 'sage.libs.pari.gen.gen'>
            sage: v.type()
            't_VEC'
            sage: e.ellan(10)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2]
            sage: v = e.ellaplist(10, python_ints=True); v
            [-2, -1, 1, -2]
            sage: type(v)
            <type 'list'>
            sage: type(v[0])
            <type 'int'>

        TESTS::

            sage: v = e.ellaplist(1)
            sage: v, type(v)
            ([], <type 'sage.libs.pari.gen.gen'>)
            sage: v = e.ellaplist(1, python_ints=True)
            sage: v, type(v)
            ([], <type 'list'>)
        """
        if python_ints:
            return [int(x) for x in self.ellaplist(n)]

        if n < 2:
            sig_on()
            return P.new_gen(zerovec(0))

        # 1. Make a table of primes up to n.
        P.init_primes(n+1)
        cdef gen t0 = objtogen(n)
        sig_on()
        cdef GEN g = primes(gtolong(primepi(t0.g)))

        # 2. Replace each prime in the table by ellap of it.
        cdef long i
        for i from 0 <= i < glength(g):
            set_gel(g, i + 1, ellap(self.g, gel(g, i + 1)))
        return P.new_gen(g)

    def ellisoncurve(self, x):
        """
        e.ellisoncurve(x): return True if the point x is on the elliptic
        curve e, False otherwise.

        If the point or the curve have inexact coefficients, an attempt is
        made to take this into account.

        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: e.ellisoncurve([1,0])
            True
            sage: e.ellisoncurve([1,1])
            False
            sage: e.ellisoncurve([1,0.00000000000000001])
            False
            sage: e.ellisoncurve([1,0.000000000000000001])
            True
            sage: e.ellisoncurve([0])
            True
        """
        cdef gen t0 = objtogen(x)
        sig_on()
        cdef int t = oncurve(self.g, t0.g)
        sig_off()
        return t != 0

    def ellminimalmodel(self):
        """
        ellminimalmodel(e): return the standard minimal integral model of
        the rational elliptic curve e and the corresponding change of
        variables. INPUT:


        -  ``e`` - gen (that defines an elliptic curve)


        OUTPUT:


        -  ``gen`` - minimal model

        -  ``gen`` - change of coordinates


        EXAMPLES::

            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: F, ch = e.ellminimalmodel()
            sage: F[:5]
            [1, -1, 0, 4, 3]
            sage: ch
            [1, -1, 0, -1]
            sage: e.ellchangecurve(ch)[:5]
            [1, -1, 0, 4, 3]
        """
        cdef GEN x, y
        cdef gen model, change
        cdef pari_sp t
        sig_on()
        x = ellminimalmodel(self.g, &y)
        change = P.new_gen_noclear(y)
        model = P.new_gen(x)
        return model, change

    def elltors(self, flag=None):
        """
        Return information about the torsion subgroup of the given
        elliptic curve.

        INPUT:

        -  ``e`` - elliptic curve over `\QQ`

        OUTPUT:


        -  ``gen`` - the order of the torsion subgroup, a.k.a.
           the number of points of finite order

        -  ``gen`` - vector giving the structure of the torsion
           subgroup as a product of cyclic groups, sorted in non-increasing
           order

        -  ``gen`` - vector giving points on e generating these
           cyclic groups


        EXAMPLES::

            sage: e = pari([1,0,1,-19,26]).ellinit()
            sage: e.elltors()
            [12, [6, 2], [[1, 2], [3, -2]]]
        """
        if flag is not None:
            deprecation(20219, 'The flag argument to elltors() is deprecated and not used anymore')
        sig_on()
        return P.new_gen(elltors(self.g))

    def omega(self, unsigned long precision=0):
        """
        Return the basis for the period lattice of this elliptic curve.

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.omega()
            [1.26920930427955, 0.634604652139777 - 1.45881661693850*I]
        """
        sig_on()
        return P.new_gen(ellR_omega(self.g, prec_bits_to_words(precision)))

    def disc(self):
        """
        Return the discriminant of this object.

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.disc()
            -161051
            sage: _.factor()
            [-1, 1; 11, 5]
        """
        sig_on()
        return P.new_gen(member_disc(self.g))

    def j(self):
        """
        Return the j-invariant of this object.

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.j()
            -122023936/161051
            sage: _.factor()
            [-1, 1; 2, 12; 11, -5; 31, 3]
        """
        sig_on()
        return P.new_gen(member_j(self.g))

    def _eltabstorel(self, x):
        """
        Return the relative number field element corresponding to `x`.

        The result is a ``t_POLMOD`` with ``t_POLMOD`` coefficients.

        .. WARNING::

            This is a low-level version of :meth:`rnfeltabstorel` that
            only needs the output of :meth:`_nf_rnfeq`, not a full
            PARI ``rnf`` structure.  This method may raise errors or
            return undefined results if called with invalid arguments.

        TESTS::

            sage: K = pari('y^2 + 1').nfinit()
            sage: rnfeq = K._nf_rnfeq(x^2 + 2)
            sage: f_abs = rnfeq[0]; f_abs
            x^4 + 6*x^2 + 1
            sage: x_rel = rnfeq._eltabstorel(x); x_rel
            Mod(x + Mod(-y, y^2 + 1), x^2 + 2)
            sage: f_abs(x_rel)
            Mod(0, x^2 + 2)

        """
        cdef gen t0 = objtogen(x)
        sig_on()
        return P.new_gen(eltabstorel(self.g, t0.g))

    def _eltabstorel_lift(self, x):
        """
        Return the relative number field element corresponding to `x`.

        The result is a ``t_POL`` with ``t_POLMOD`` coefficients.

        .. WARNING::

            This is a low-level version of :meth:`rnfeltabstorel` that
            only needs the output of :meth:`_nf_rnfeq`, not a full
            PARI ``rnf`` structure.  This method may raise errors or
            return undefined results if called with invalid arguments.

        TESTS::

            sage: K = pari('y^2 + 1').nfinit()
            sage: rnfeq = K._nf_rnfeq(x^2 + 2)
            sage: rnfeq._eltabstorel_lift(x)
            x + Mod(-y, y^2 + 1)

        """
        cdef gen t0 = objtogen(x)
        sig_on()
        return P.new_gen(eltabstorel_lift(self.g, t0.g))

    def _eltreltoabs(self, x):
        """
        Return the absolute number field element corresponding to `x`.

        The result is a ``t_POL``.

        .. WARNING::

            This is a low-level version of :meth:`rnfeltreltoabs` that
            only needs the output of :meth:`_nf_rnfeq`, not a full
            PARI ``rnf`` structure.  This method may raise errors or
            return undefined results if called with invalid arguments.

        TESTS::

            sage: K = pari('y^2 + 1').nfinit()
            sage: rnfeq = K._nf_rnfeq(x^2 + 2)
            sage: rnfeq._eltreltoabs(x)
            1/2*x^3 + 7/2*x
            sage: rnfeq._eltreltoabs('y')
            1/2*x^3 + 5/2*x

        """
        cdef gen t0 = objtogen(x)
        sig_on()
        return P.new_gen(eltreltoabs(self.g, t0.g))

    def galoissubfields(self, long flag=0, v=-1):
        """
        List all subfields of the Galois group ``self``.

        This wraps the `galoissubfields`_ function from PARI.

        This method is essentially the same as applying
        :meth:`galoisfixedfield` to each group returned by
        :meth:`galoissubgroups`.

        INPUT:

        - ``self`` -- A Galois group as generated by :meth:`galoisinit`.

        - ``flag`` -- Has the same meaning as in :meth:`galoisfixedfield`.

        - ``v`` -- Has the same meaning as in :meth:`galoisfixedfield`.

        OUTPUT:

        A vector of all subfields of this group.  Each entry is as
        described in the :meth:`galoisfixedfield` method.

        EXAMPLES::

            sage: G = pari(x^6 + 108).galoisinit()
            sage: G.galoissubfields(flag=1)
            [x, x^2 + 972, x^3 + 54, x^3 + 864, x^3 - 54, x^6 + 108]
            sage: G = pari(x^4 + 1).galoisinit()
            sage: G.galoissubfields(flag=2, v='z')[3]
            [x^2 + 2, Mod(x^3 + x, x^4 + 1), [x^2 - z*x - 1, x^2 + z*x - 1]]

        .. _galoissubfields: http://pari.math.u-bordeaux.fr/dochtml/html.stable/Functions_related_to_general_number_fields.html#galoissubfields
        """
        sig_on()
        return P.new_gen(galoissubfields(self.g, flag, P.get_var(v)))

    idealintersection = deprecated_function_alias(20219, gen_auto.idealintersect)

    def nfeltval(self, x, p):
        """
        Return the valuation of the number field element `x` at the prime `p`.

        EXAMPLES::

            sage: nf = pari('x^2 + 1').nfinit()
            sage: p = nf.idealprimedec(5)[0]
            sage: nf.nfeltval('50 - 25*x', p)
            3
        """
        cdef gen t0 = objtogen(x)
        cdef gen t1 = objtogen(p)
        sig_on()
        v = nfval(self.g, t0.g, t1.g)
        sig_off()
        return v

    elementval = deprecated_function_alias(20219, nfeltval)

    def nfbasis(self, long flag=0, fa=None):
        """
        Integral basis of the field `\QQ[a]`, where ``a`` is a root of
        the polynomial x.

        INPUT:

        - ``flag``: if set to 1 and ``fa`` is not given: assume that no
          square of a prime > 500000 divides the discriminant of ``x``.

        - ``fa``: If present, encodes a subset of primes at which to
          check for maximality. This must be one of the three following
          things:

            - an integer: check all primes up to ``fa`` using trial
              division.

            - a vector: a list of primes to check.

            - a matrix: a partial factorization of the discriminant
              of ``x``.

        .. NOTE::

            In earlier versions of Sage, other bits in ``flag`` were
            defined but these are now simply ignored.

        EXAMPLES::

            sage: pari('x^3 - 17').nfbasis()
            [1, x, 1/3*x^2 - 1/3*x + 1/3]

        We test ``flag`` = 1, noting it gives a wrong result when the
        discriminant (-4 * `p`^2 * `q` in the example below) has a big square
        factor::

            sage: p = next_prime(10^10); q = next_prime(p)
            sage: x = polygen(QQ); f = x^2 + p^2*q
            sage: pari(f).nfbasis(1)   # Wrong result
            [1, x]
            sage: pari(f).nfbasis()    # Correct result
            [1, 1/10000000019*x]
            sage: pari(f).nfbasis(fa=10^6)   # Check primes up to 10^6: wrong result
            [1, x]
            sage: pari(f).nfbasis(fa="[2,2; %s,2]"%p)    # Correct result and faster
            [1, 1/10000000019*x]
            sage: pari(f).nfbasis(fa=[2,p])              # Equivalent with the above
            [1, 1/10000000019*x]
        """
        cdef gen t0
        cdef GEN g0
        if fa is not None:
            t0 = objtogen(fa)
            g0 = t0.g
        elif flag:
            g0 = utoi(500000)
        else:
            g0 = NULL
        sig_on()
        return P.new_gen(nfbasis(self.g, NULL, g0))

    def nfbasis_d(self, long flag=0, fa=None):
        """
        Like :meth:`nfbasis`, but return a tuple ``(B, D)`` where `B`
        is the integral basis and `D` the discriminant.

        EXAMPLES::

            sage: F = NumberField(x^3-2,'alpha')
            sage: F._pari_()[0].nfbasis_d()
            ([1, y, y^2], -108)

        ::

            sage: G = NumberField(x^5-11,'beta')
            sage: G._pari_()[0].nfbasis_d()
            ([1, y, y^2, y^3, y^4], 45753125)

        ::

            sage: pari([-2,0,0,1]).Polrev().nfbasis_d()
            ([1, x, x^2], -108)
        """
        cdef gen t0
        cdef GEN g0
        cdef GEN disc
        if fa is not None:
            t0 = objtogen(fa)
            g0 = t0.g
        elif flag & 1:
            g0 = utoi(500000)
        else:
            g0 = NULL
        sig_on()
        B = P.new_gen_noclear(nfbasis(self.g, &disc, g0))
        D = P.new_gen(disc)
        return B, D

    def nfbasistoalg_lift(nf, x):
        r"""
        Transforms the column vector ``x`` on the integral basis into a
        polynomial representing the algebraic number.

        INPUT:

         - ``nf`` -- a number field
         - ``x`` -- a column of rational numbers of length equal to the
           degree of ``nf`` or a single rational number

        OUTPUT:

         - ``nf.nfbasistoalg(x).lift()``

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 - 17)
            sage: Kpari = K.pari_nf()
            sage: Kpari.getattr('zk')
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
            sage: Kpari.nfbasistoalg_lift(42)
            42
            sage: Kpari.nfbasistoalg_lift("[3/2, -5, 0]~")
            -5/3*y^2 + 5/3*y - 1/6
            sage: Kpari.getattr('zk') * pari("[3/2, -5, 0]~")
            -5/3*y^2 + 5/3*y - 1/6
        """
        cdef gen t0 = objtogen(x)
        sig_on()
        return P.new_gen(gel(basistoalg(nf.g, t0.g), 2))

    def nfgenerator(self):
        f = self[0]
        x = f.variable()
        return x.Mod(f)

    def _nf_rnfeq(self, relpol):
        """
        Return data for converting number field elements between
        absolute and relative representation.

        .. NOTE::

            The output of this method is suitable for the methods
            :meth:`_eltabstorel`, :meth:`_eltabstorel_lift` and
            :meth:`_eltreltoabs`.

        TESTS::

            sage: K = pari('y^2 + 1').nfinit()
            sage: K._nf_rnfeq(x^2 + 2)
            [x^4 + 6*x^2 + 1, 1/2*x^3 + 5/2*x, -1, y^2 + 1, x^2 + 2]

        """
        cdef gen t0 = objtogen(relpol)
        sig_on()
        return P.new_gen(nf_rnfeq(self.g, t0.g))

    def _nf_nfzk(self, rnfeq):
        """
        Return data for constructing relative number field elements
        from elements of the base field.

        INPUT:

        - ``rnfeq`` -- relative number field data as returned by
          :meth:`_nf_rnfeq`

        .. NOTE::

            The output of this method is suitable for the method
            :meth:`_nfeltup`.

        TESTS::

            sage: nf = pari('nfinit(y^2 - 2)')
            sage: nf._nf_nfzk(nf._nf_rnfeq('x^2 - 3'))
            ([2, -x^3 + 9*x], 1/2)

        """
        cdef GEN zknf, czknf
        cdef gen t0 = objtogen(rnfeq)
        cdef gen zk, czk
        sig_on()
        nf_nfzk(self.g, t0.g, &zknf, &czknf)
        zk = P.new_gen_noclear(zknf)
        czk = P.new_gen(czknf)
        return zk, czk

    def _nfeltup(self, x, zk, czk):
        """
        Construct a relative number field element from an element of
        the base field.

        INPUT:

        - ``x`` -- element of the base field

        - ``zk``, ``czk`` -- relative number field data as returned by
          :meth:`_nf_nfzk`

        .. WARNING::

            This is a low-level version of :meth:`rnfeltup` that only
            needs the output of :meth:`_nf_nfzk`, not a full PARI
            ``rnf`` structure.  This method may raise errors or return
            undefined results if called with invalid arguments.

        TESTS::

            sage: nf = pari('nfinit(y^2 - 2)')
            sage: zk, czk = nf._nf_nfzk(nf._nf_rnfeq('x^2 - 3'))
            sage: nf._nfeltup('y', zk, czk)
            -1/2*x^3 + 9/2*x

        """
        cdef gen t0 = objtogen(x)
        cdef gen t1 = objtogen(zk)
        cdef gen t2 = objtogen(czk)
        sig_on()
        return P.new_gen(nfeltup(self.g, t0.g, t1.g, t2.g))

    reverse = deprecated_function_alias(20219, gen_auto.polrecip)

    def eval(self, *args, **kwds):
        """
        Evaluate ``self`` with the given arguments.

        This is currently implemented in 3 cases:

        - univariate polynomials, rational functions, power series and
          Laurent series (using a single unnamed argument or keyword
          arguments),
        - any PARI object supporting the PARI function ``substvec``
          (in particular, multivariate polynomials) using keyword
          arguments,
        - objects of type ``t_CLOSURE`` (functions in GP bytecode form)
          using unnamed arguments.

        In no case is mixing unnamed and keyword arguments allowed.

        EXAMPLES::

            sage: f = pari('x^2 + 1')
            sage: f.type()
            't_POL'
            sage: f.eval(I)
            0
            sage: f.eval(x=2)
            5
            sage: (1/f).eval(x=1)
            1/2

        The notation ``f(x)`` is an alternative for ``f.eval(x)``::

            sage: f(3) == f.eval(3)
            True

        Evaluating power series::

            sage: f = pari('1 + x + x^3 + O(x^7)')
            sage: f(2*pari('y')^2)
            1 + 2*y^2 + 8*y^6 + O(y^14)

        Substituting zero is sometimes possible, and trying to do so
        in illegal cases can raise various errors::

            sage: pari('1 + O(x^3)').eval(0)
            1
            sage: pari('1/x').eval(0)
            Traceback (most recent call last):
            ...
            PariError: impossible inverse in gdiv: 0
            sage: pari('1/x + O(x^2)').eval(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: substituting 0 in Laurent series with negative valuation
            sage: pari('1/x + O(x^2)').eval(pari('O(x^3)'))
            Traceback (most recent call last):
            ...
            PariError: impossible inverse in gdiv: O(x^3)
            sage: pari('O(x^0)').eval(0)
            Traceback (most recent call last):
            ...
            PariError: domain error in polcoeff: t_SER = O(x^0)

        Evaluating multivariate polynomials::

            sage: f = pari('y^2 + x^3')
            sage: f(1)    # Dangerous, depends on PARI variable ordering
            y^2 + 1
            sage: f(x=1)  # Safe
            y^2 + 1
            sage: f(y=1)
            x^3 + 1
            sage: f(1, 2)
            Traceback (most recent call last):
            ...
            TypeError: evaluating PARI t_POL takes exactly 1 argument (2 given)
            sage: f(y='x', x='2*y')
            x^2 + 8*y^3
            sage: f()
            x^3 + y^2

        It's not an error to substitute variables which do not appear::

            sage: f.eval(z=37)
            x^3 + y^2
            sage: pari(42).eval(t=0)
            42

        We can define and evaluate closures as follows::

            sage: T = pari('n -> n + 2')
            sage: T.type()
            't_CLOSURE'
            sage: T.eval(3)
            5

            sage: T = pari('() -> 42')
            sage: T()
            42

            sage: pr = pari('s -> print(s)')
            sage: pr.eval('"hello world"')
            hello world

            sage: f = pari('myfunc(x,y) = x*y')
            sage: f.eval(5, 6)
            30

        Default arguments work, missing arguments are treated as zero
        (like in GP)::

            sage: f = pari("(x, y, z=1.0) -> [x, y, z]")
            sage: f(1, 2, 3)
            [1, 2, 3]
            sage: f(1, 2)
            [1, 2, 1.00000000000000]
            sage: f(1)
            [1, 0, 1.00000000000000]
            sage: f()
            [0, 0, 1.00000000000000]

        Variadic closures are supported as well (:trac:`18623`)::

            sage: f = pari("(v[..])->length(v)")
            sage: f('a', f)
            2
            sage: g = pari("(x,y,z[..])->[x,y,z]")
            sage: g(), g(1), g(1,2), g(1,2,3), g(1,2,3,4)
            ([0, 0, []], [1, 0, []], [1, 2, []], [1, 2, [3]], [1, 2, [3, 4]])

        Using keyword arguments, we can substitute in more complicated
        objects, for example a number field::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: nf = K._pari_()
            sage: nf
            [y^2 + 1, [0, 1], -4, 1, [Mat([1, 0.E-38 + 1.00000000000000*I]), [1, 1.00000000000000; 1, -1.00000000000000], [1, 1; 1, -1], [2, 0; 0, -2], [2, 0; 0, 2], [1, 0; 0, -1], [1, [0, -1; 1, 0]], []], [0.E-38 + 1.00000000000000*I], [1, y], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, 0]]
            sage: nf(y='x')
            [x^2 + 1, [0, 1], -4, 1, [Mat([1, 0.E-38 + 1.00000000000000*I]), [1, 1.00000000000000; 1, -1.00000000000000], [1, 1; 1, -1], [2, 0; 0, -2], [2, 0; 0, 2], [1, 0; 0, -1], [1, [0, -1; 1, 0]], []], [0.E-38 + 1.00000000000000*I], [1, x], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, 0]]
        """
        cdef long t = typ(self.g)
        cdef gen t0
        cdef GEN result
        cdef long arity
        cdef long nargs = len(args)
        cdef long nkwds = len(kwds)

        # Closure must be evaluated using *args
        if t == t_CLOSURE:
            if nkwds > 0:
                raise TypeError("cannot evaluate a PARI closure using keyword arguments")
            if closure_is_variadic(self.g):
                arity = closure_arity(self.g) - 1
                args = list(args[:arity]) + [0]*(arity-nargs) + [args[arity:]]
            t0 = objtogen(args)
            sig_on()
            result = closure_callgenvec(self.g, t0.g)
            if result == gnil:
                P.clear_stack()
                return None
            return P.new_gen(result)

        # Evaluate univariate polynomials, rational functions and
        # series using *args
        if nargs > 0:
            if nkwds > 0:
                raise TypeError("mixing unnamed and keyword arguments not allowed when evaluating a PARI object")
            if not (t == t_POL or t == t_RFRAC or t == t_SER):
                raise TypeError("cannot evaluate PARI %s using unnamed arguments" % self.type())
            if nargs != 1:
                raise TypeError("evaluating PARI %s takes exactly 1 argument (%d given)"
                                % (self.type(), nargs))

            t0 = objtogen(args[0])
            sig_on()
            if t == t_POL or t == t_RFRAC:
                return P.new_gen(poleval(self.g, t0.g))
            else:  # t == t_SER
                if isexactzero(t0.g):
                    # Work around the fact that PARI currently doesn't
                    # support substituting exact 0 in a power series.
                    # We don't try to imitate this when using keyword
                    # arguments, and hope this will be fixed in a
                    # future PARI version.
                    if valp(self.g) < 0:
                        sig_off()
                        raise ZeroDivisionError('substituting 0 in Laurent series with negative valuation')
                    elif valp(self.g) == 0:
                        return P.new_gen(polcoeff0(self.g, 0, -1))
                    else:
                        sig_off()
                        return P.PARI_ZERO
                return P.new_gen(gsubst(self.g, varn(self.g), t0.g))

        # Call substvec() using **kwds
        vstr = kwds.keys()            # Variables as Python strings
        t0 = objtogen(kwds.values())  # Replacements

        sig_on()
        cdef GEN v = cgetg(nkwds+1, t_VEC)  # Variables as PARI polynomials
        cdef long i
        for i in range(nkwds):
            set_gel(v, i+1, pol_x(P.get_var(vstr[i])))
        return P.new_gen(gsubstvec(self.g, v, t0.g))

    def __call__(self, *args, **kwds):
        """
        Evaluate ``self`` with the given arguments.

        This has the same effect as :meth:`eval`.

        EXAMPLES::

            sage: R.<x> = GF(3)[]
            sage: f = (x^2 + x + 1)._pari_()
            sage: f.type()
            't_POL'
            sage: f(2)
            Mod(1, 3)

        TESTS::

            sage: T = pari('n -> 1/n')
            sage: T.type()
            't_CLOSURE'
            sage: T(0)
            Traceback (most recent call last):
            ...
            PariError: _/_: impossible inverse in gdiv: 0
            sage: pari('() -> 42')(1,2,3)
            Traceback (most recent call last):
            ...
            PariError: too many parameters in user-defined function call
            sage: pari('n -> n')(n=2)
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate a PARI closure using keyword arguments
            sage: pari('x + y')(4, y=1)
            Traceback (most recent call last):
            ...
            TypeError: mixing unnamed and keyword arguments not allowed when evaluating a PARI object
            sage: pari("12345")(4)
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate PARI t_INT using unnamed arguments
        """
        return self.eval(*args, **kwds)

    def factorpadic(self, p, long r=20):
        """
        p-adic factorization of the polynomial ``pol`` to precision ``r``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: pol = (x^2 - 1)^2
            sage: pari(pol).factorpadic(5)
            [(1 + O(5^20))*x + (1 + O(5^20)), 2; (1 + O(5^20))*x + (4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20)), 2]
            sage: pari(pol).factorpadic(5,3)
            [(1 + O(5^3))*x + (1 + O(5^3)), 2; (1 + O(5^3))*x + (4 + 4*5 + 4*5^2 + O(5^3)), 2]
        """
        cdef gen t0 = objtogen(p)
        sig_on()
        return P.new_gen(factorpadic(self.g, t0.g, r))

    def poldegree(self, var=-1):
        """
        Return the degree of this polynomial.
        """
        sig_on()
        n = poldegree(self.g, P.get_var(var))
        sig_off()
        return n

    def polisirreducible(self):
        """
        f.polisirreducible(): Returns True if f is an irreducible
        non-constant polynomial, or False if f is reducible or constant.
        """
        sig_on()
        t = isirreducible(self.g)
        P.clear_stack()
        return t != 0

    def polroots(self, unsigned long precision=0):
        """
        Complex roots of the given polynomial using Schonhage's method,
        as modified by Gourdon.
        """
        sig_on()
        return P.new_gen(cleanroots(self.g, prec_bits_to_words(precision)))

    def rnfisnorm(self, T, long flag=0):
        cdef gen t0 = objtogen(T)
        sig_on()
        return P.new_gen(rnfisnorm(t0.g, self.g, flag))

    def ncols(self):
        """
        Return the number of columns of self.

        EXAMPLES::

            sage: pari('matrix(19,8)').ncols()
            8
        """
        cdef long n
        sig_on()
        n = glength(self.g)
        sig_off()
        return n

    def nrows(self):
        """
        Return the number of rows of self.

        EXAMPLES::

            sage: pari('matrix(19,8)').nrows()
            19
        """
        cdef long n
        sig_on()
        # if this matrix has no columns
        # then it has no rows.
        if self.ncols() == 0:
            sig_off()
            return 0
        n = glength(<GEN>(self.g[1]))
        sig_off()
        return n

    def mattranspose(self):
        """
        Transpose of the matrix self.

        EXAMPLES::

            sage: pari('[1,2,3; 4,5,6; 7,8,9]').mattranspose()
            [1, 4, 7; 2, 5, 8; 3, 6, 9]

        Unlike PARI, this always returns a matrix::

            sage: pari('[1,2,3]').mattranspose()
            [1; 2; 3]
            sage: pari('[1,2,3]~').mattranspose()
            Mat([1, 2, 3])
        """
        sig_on()
        return P.new_gen(gtrans(self.g)).Mat()

    def lllgram(self):
        return self.qflllgram(0)

    def lllgramint(self):
        return self.qflllgram(1)

    def qfrep(self, B, long flag=0):
        """
        Vector of (half) the number of vectors of norms from 1 to `B`
        for the integral and definite quadratic form ``self``.
        Binary digits of flag mean 1: count vectors of even norm from
        1 to `2B`, 2: return a ``t_VECSMALL`` instead of a ``t_VEC``
        (which is faster).

        EXAMPLES::

            sage: M = pari("[5,1,1;1,3,1;1,1,1]")
            sage: M.qfrep(20)
            [1, 1, 2, 2, 2, 4, 4, 3, 3, 4, 2, 4, 6, 0, 4, 6, 4, 5, 6, 4]
            sage: M.qfrep(20, flag=1)
            [1, 2, 4, 3, 4, 4, 0, 6, 5, 4, 12, 4, 4, 8, 0, 3, 8, 6, 12, 12]
            sage: M.qfrep(20, flag=2)
            Vecsmall([1, 1, 2, 2, 2, 4, 4, 3, 3, 4, 2, 4, 6, 0, 4, 6, 4, 5, 6, 4])
        """
        # PARI 2.7 always returns a t_VECSMALL, but for backwards
        # compatibility, we keep returning a t_VEC (unless flag & 2)
        cdef gen t0 = objtogen(B)
        cdef GEN r
        sig_on()
        r = qfrep0(self.g, t0.g, flag & 1)
        if (flag & 2) == 0:
            r = vecsmall_to_vec(r)
        return P.new_gen(r)

    def matkerint(self, long flag=0):
        """
        Return the integer kernel of a matrix.

        This is the LLL-reduced Z-basis of the kernel of the matrix x with
        integral entries.

        EXAMPLES::

            sage: pari('[2,1;2,1]').matker()
            [-1/2; 1]
            sage: pari('[2,1;2,1]').matkerint()
            [1; -2]
            sage: pari('[2,1;2,1]').matkerint(1)
            doctest:...: DeprecationWarning: The flag argument to matkerint() is deprecated by PARI
            See http://trac.sagemath.org/18203 for details.
            [1; -2]
        """
        if flag:
            # Keep this deprecation warning as long as PARI supports
            # this deprecated flag
            deprecation(18203, "The flag argument to matkerint() is deprecated by PARI")
        sig_on()
        return P.new_gen(matkerint0(self.g, flag))

    def factor(self, long limit=-1, proof=None):
        """
        Return the factorization of x.

        INPUT:

        -  ``limit`` -- (default: -1) is optional and can be set
           whenever x is of (possibly recursive) rational type. If limit is
           set, return partial factorization, using primes up to limit.

        - ``proof`` -- optional flag. If ``False`` (not the default),
          returned factors larger than `2^{64}` may only be pseudoprimes.
          If ``True``, always check primality. If not given, use the
          global PARI default ``factor_proven`` which is ``True`` by
          default in Sage.

        EXAMPLES::

            sage: pari('x^10-1').factor()
            [x - 1, 1; x + 1, 1; x^4 - x^3 + x^2 - x + 1, 1; x^4 + x^3 + x^2 + x + 1, 1]
            sage: pari(2^100-1).factor()
            [3, 1; 5, 3; 11, 1; 31, 1; 41, 1; 101, 1; 251, 1; 601, 1; 1801, 1; 4051, 1; 8101, 1; 268501, 1]
            sage: pari(2^100-1).factor(proof=True)
            [3, 1; 5, 3; 11, 1; 31, 1; 41, 1; 101, 1; 251, 1; 601, 1; 1801, 1; 4051, 1; 8101, 1; 268501, 1]
            sage: pari(2^100-1).factor(proof=False)
            [3, 1; 5, 3; 11, 1; 31, 1; 41, 1; 101, 1; 251, 1; 601, 1; 1801, 1; 4051, 1; 8101, 1; 268501, 1]

        We illustrate setting a limit::

            sage: pari(next_prime(10^50)*next_prime(10^60)*next_prime(10^4)).factor(10^5)
            [10007, 1; 100000000000000000000000000000000000000000000000151000000000700000000000000000000000000000000000000000000001057, 1]

        Setting a limit is invalid when factoring polynomials::

            sage: pari('x^11 + 1').factor(limit=17)
            Traceback (most recent call last):
            ...
            PariError: incorrect type in boundfact (t_POL)

        PARI doesn't have an algorithm for factoring multivariate
        polynomials::

            sage: pari('x^3 - y^3').factor()
            Traceback (most recent call last):
            ...
            PariError: sorry, factor for general polynomials is not yet implemented

        TESTS::

            sage: pari(2^1000+1).factor(limit=0)
            doctest:...: DeprecationWarning: factor(..., lim=0) is deprecated, use an explicit limit instead
            See http://trac.sagemath.org/20205 for details.
            [257, 1; 1601, 1; 25601, 1; 76001, 1; 133842787352016..., 1]
        """
        cdef GEN g
        if limit == 0:
            deprecation(20205, "factor(..., lim=0) is deprecated, use an explicit limit instead")
            limit = maxprime()
        global factor_proven
        cdef int saved_factor_proven = factor_proven
        try:
            if proof is not None:
                factor_proven = 1 if proof else 0
            sig_on()
            if limit >= 0:
                g = boundfact(self.g, limit)
            else:
                g = factor(self.g)
            return P.new_gen(g)
        finally:
            factor_proven = saved_factor_proven

    multiplicative_order = gen_auto.znorder
    order = deprecated_function_alias(20219, multiplicative_order)

    def __abs__(self):
        return self.abs()

    def nextprime(gen self, bint add_one=0):
        """
        nextprime(x): smallest pseudoprime greater than or equal to `x`.
        If ``add_one`` is non-zero, return the smallest pseudoprime
        strictly greater than `x`.

        EXAMPLES::

            sage: pari(1).nextprime()
            2
            sage: pari(2).nextprime()
            2
            sage: pari(2).nextprime(add_one = 1)
            3
            sage: pari(2^100).nextprime()
            1267650600228229401496703205653
        """
        sig_on()
        if add_one:
            return P.new_gen(nextprime(gaddsg(1, self.g)))
        return P.new_gen(nextprime(self.g))

    def change_variable_name(self, var):
        """
        In ``self``, which must be a ``t_POL`` or ``t_SER``, set the
        variable to ``var``.  If the variable of ``self`` is already
        ``var``, then return ``self``.

        .. WARNING::

            You should be careful with variable priorities when
            applying this on a polynomial or series of which the
            coefficients have polynomial components.  To be safe, only
            use this function on polynomials with integer or rational
            coefficients.  For a safer alternative, use :meth:`subst`.

        EXAMPLES::

            sage: f = pari('x^3 + 17*x + 3')
            sage: f.change_variable_name("y")
            y^3 + 17*y + 3
            sage: f = pari('1 + 2*y + O(y^10)')
            sage: f.change_variable_name("q")
            1 + 2*q + O(q^10)
            sage: f.change_variable_name("y") is f
            True

        In PARI, ``I`` refers to the square root of -1, so it cannot be
        used as variable name.  Note the difference with :meth:`subst`::

            sage: f = pari('x^2 + 1')
            sage: f.change_variable_name("I")
            Traceback (most recent call last):
            ...
            PariError: I already exists with incompatible valence
            sage: f.subst("x", "I")
            0
        """
        cdef long n = P.get_var(var)
        if varn(self.g) == n:
            return self
        if typ(self.g) != t_POL and typ(self.g) != t_SER:
            raise TypeError("set_variable() only works for polynomials or power series")
        # Copy self and then change the variable in place
        cdef gen newg = P.new_gen_noclear(self.g)
        setvarn(newg.g, n)
        return newg

    def nf_subst(self, z):
        """
        Given a PARI number field ``self``, return the same PARI
        number field but in the variable ``z``.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K = NumberField(x^2 + 5, 'a')

        We can substitute in a PARI ``nf`` structure::

            sage: Kpari = K.pari_nf()
            sage: Kpari.nf_get_pol()
            y^2 + 5
            sage: Lpari = Kpari.nf_subst('a')
            sage: Lpari.nf_get_pol()
            a^2 + 5

        We can also substitute in a PARI ``bnf`` structure::

            sage: Kpari = K.pari_bnf()
            sage: Kpari.nf_get_pol()
            y^2 + 5
            sage: Kpari.bnf_get_cyc()  # Structure of class group
            [2]
            sage: Lpari = Kpari.nf_subst('a')
            sage: Lpari.nf_get_pol()
            a^2 + 5
            sage: Lpari.bnf_get_cyc()  # We still have a bnf after substituting
            [2]
        """
        cdef gen t0 = objtogen(z)
        sig_on()
        return P.new_gen(gsubst(self.g, gvar(self.g), t0.g))

    def type(gen self):
        """
        Return the PARI type of self as a string.

        .. note::

           In Cython, it is much faster to simply use typ(self.g) for
           checking PARI types.

        EXAMPLES::

            sage: pari(7).type()
            't_INT'
            sage: pari('x').type()
            't_POL'
            sage: pari('oo').type()
            't_INFINITY'
        """
        # The following original code leaks memory:
        #        return str(type_name(typ(self.g)))
        #
        # This code is the usual workaround:
        #        cdef char* s= <char*>type_name(typ(self.g))
        #        t=str(s)
        #        free(s)
        #        return(t)
        # However, it causes segfaults with t_INTs on some
        # machines, and errors about freeing non-aligned
        # pointers on others. So we settle for the following
        # fast but ugly code. Note that should the list of
        # valid PARI types ever be updated, this code would
        # need to be updated accordingly.
        #
        cdef long t = typ(self.g)

        if   t == t_INT:      return 't_INT'
        elif t == t_REAL:     return 't_REAL'
        elif t == t_INTMOD:   return 't_INTMOD'
        elif t == t_FRAC:     return 't_FRAC'
        elif t == t_FFELT:    return 't_FFELT'
        elif t == t_COMPLEX:  return 't_COMPLEX'
        elif t == t_PADIC:    return 't_PADIC'
        elif t == t_QUAD:     return 't_QUAD'
        elif t == t_POLMOD:   return 't_POLMOD'
        elif t == t_POL:      return 't_POL'
        elif t == t_SER:      return 't_SER'
        elif t == t_RFRAC:    return 't_RFRAC'
        elif t == t_QFR:      return 't_QFR'
        elif t == t_QFI:      return 't_QFI'
        elif t == t_VEC:      return 't_VEC'
        elif t == t_COL:      return 't_COL'
        elif t == t_MAT:      return 't_MAT'
        elif t == t_LIST:     return 't_LIST'
        elif t == t_STR:      return 't_STR'
        elif t == t_VECSMALL: return 't_VECSMALL'
        elif t == t_CLOSURE:  return 't_CLOSURE'
        elif t == t_INFINITY: return 't_INFINITY'
        else:
            raise TypeError("Unknown PARI type: %s" % t)

    def polinterpolate(self, ya, x):
        """
        self.polinterpolate(ya,x,e): polynomial interpolation at x
        according to data vectors self, ya (i.e. return P such that
        P(self[i]) = ya[i] for all i). Also return an error estimate on the
        returned value.
        """
        cdef gen t0 = objtogen(ya)
        cdef gen t1 = objtogen(x)
        cdef GEN dy, g
        sig_on()
        g = polint(self.g, t0.g, t1.g, &dy)
        dif = P.new_gen_noclear(dy)
        return P.new_gen(g), dif

    def ellwp(gen self, z='z', long n=20, long flag=0, unsigned long precision=0):
        """
        Return the value or the series expansion of the Weierstrass
        `P`-function at `z` on the lattice `self` (or the lattice
        defined by the elliptic curve `self`).

        INPUT:

        -  ``self`` -- an elliptic curve created using ``ellinit`` or a
           list ``[om1, om2]`` representing generators for a lattice.

        -  ``z`` -- (default: 'z') a complex number or a variable name
           (as string or PARI variable).

        -  ``n`` -- (default: 20) if 'z' is a variable, compute the
           series expansion up to at least `O(z^n)`.

        -  ``flag`` -- (default = 0): If ``flag`` is 0, compute only
           `P(z)`.  If ``flag`` is 1, compute `[P(z), P'(z)]`.

        OUTPUT:

        - `P(z)` (if ``flag`` is 0) or `[P(z), P'(z)]` (if ``flag`` is 1).
           numbers

        EXAMPLES:

        We first define the elliptic curve X_0(11)::

            sage: E = pari([0,-1,1,-10,-20]).ellinit()

        Compute P(1)::

            sage: E.ellwp(1)
            13.9658695257485

        Compute P(1+i), where i = sqrt(-1)::

            sage: C.<i> = ComplexField()
            sage: E.ellwp(pari(1+i))
            -1.11510682565555 + 2.33419052307470*I
            sage: E.ellwp(1+i)
            -1.11510682565555 + 2.33419052307470*I

        The series expansion, to the default `O(z^20)` precision::

            sage: E.ellwp()
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)

        Compute the series for wp to lower precision::

            sage: E.ellwp(n=4)
            z^-2 + 31/15*z^2 + O(z^4)

        Next we use the version where the input is generators for a
        lattice::

            sage: pari([1.2692, 0.63 + 1.45*i]).ellwp(1)
            13.9656146936689 + 0.000644829272810...*I

        With flag=1, compute the pair P(z) and P'(z)::

            sage: E.ellwp(1, flag=1)
            [13.9658695257485, 50.5619300880073]
        """
        cdef gen t0 = objtogen(z)
        cdef GEN g0 = t0.g

        # Emulate toser_i() but with given precision
        sig_on()
        if typ(g0) == t_POL:
            g0 = RgX_to_ser(g0, n+4)
        elif typ(g0) == t_RFRAC:
            g0 = rfrac_to_ser(g0, n+4)
        return P.new_gen(ellwp0(self.g, g0, flag, prec_bits_to_words(precision)))

    def debug(gen self, long depth = -1):
        r"""
        Show the internal structure of self (like the ``\x`` command in gp).

        EXAMPLE::

            sage: pari('[1/2, 1.0*I]').debug()  # random addresses
            [&=0000000004c5f010] VEC(lg=3):2200000000000003 0000000004c5eff8 0000000004c5efb0
              1st component = [&=0000000004c5eff8] FRAC(lg=3):0800000000000003 0000000004c5efe0 0000000004c5efc8
                num = [&=0000000004c5efe0] INT(lg=3):0200000000000003 (+,lgefint=3):4000000000000003 0000000000000001
                den = [&=0000000004c5efc8] INT(lg=3):0200000000000003 (+,lgefint=3):4000000000000003 0000000000000002
              2nd component = [&=0000000004c5efb0] COMPLEX(lg=3):0c00000000000003 00007fae8a2eb840 0000000004c5ef90
                real = gen_0
                imag = [&=0000000004c5ef90] REAL(lg=4):0400000000000004 (+,expo=0):6000000000000000 8000000000000000 0000000000000000
        """
        sig_on()
        dbgGEN(self.g, depth)
        sig_off()
        return

    ####################################################################
    # Functions deprecated by upstream PARI
    #
    # NOTE: these should remain in Sage as long as PARI supports them,
    # do not just delete these methods!
    ####################################################################

    def bezout(x, y):
        deprecation(18203, "bezout() is deprecated in PARI, use gcdext() instead (note that the output is in a different order!)")
        u, v, g = x.gcdext(y)
        return g, u, v

    def sizedigit(x):
        """
        sizedigit(x): Return a quick estimate for the maximal number of
        decimal digits before the decimal point of any component of x.

        INPUT:

        -  ``x`` - gen

        OUTPUT: Python integer

        EXAMPLES::

            sage: x = pari('10^100')
            sage: x.Str().length()
            101
            sage: x.sizedigit()
            doctest:...: DeprecationWarning: sizedigit() is deprecated in PARI
            See http://trac.sagemath.org/18203 for details.
            101

        Note that digits after the decimal point are ignored::

            sage: x = pari('1.234')
            sage: x
            1.23400000000000
            sage: x.sizedigit()
            1

        The estimate can be one too big::

            sage: pari('7234.1').sizedigit()
            4
            sage: pari('9234.1').sizedigit()
            5
        """
        deprecation(18203, "sizedigit() is deprecated in PARI")
        return sizedigit(x.g)

    def bernvec(x):
        r"""
        Creates a vector containing, as rational numbers, the Bernoulli
        numbers `B_0, B_2,\ldots, B_{2x}`. This routine is
        obsolete. Use bernfrac instead each time you need a Bernoulli
        number in exact form.

        Note: this routine is implemented using repeated independent calls
        to bernfrac, which is faster than the standard recursion in exact
        arithmetic.

        EXAMPLES::

            sage: pari(8).bernvec()
            doctest:...: DeprecationWarning: bernvec() is deprecated, use repeated calls to bernfrac() instead
            See http://trac.sagemath.org/15767 for details.
            [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
            sage: [pari(2*n).bernfrac() for n in range(9)]
            [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
        """
        deprecation(15767, 'bernvec() is deprecated, use repeated calls to bernfrac() instead')
        sig_on()
        return P.new_gen(bernvec(x))

    bezoutres = deprecated_function_alias(18203, gen_auto.polresultantext)

    ellbil = deprecated_function_alias(18203, gen_auto.ellheight)

    ellpow = deprecated_function_alias(18203, gen_auto.ellmul)

    def rnfpolred(*args, **kwds):
        deprecation(18203, "rnfpolred() is deprecated in PARI, port your code to use rnfpolredbest() instead")
        return gen_auto.rnfpolred(*args, **kwds)

    def rnfpolredabs(*args, **kwds):
        deprecation(18203, "rnfpolredabs() is deprecated in PARI, port your code to use rnfpolredbest() instead")
        return gen_auto.rnfpolredabs(*args, **kwds)


cpdef gen objtogen(s):
    """
    Convert any Sage/Python object to a PARI gen.

    For Sage types, this uses the `_pari_()` method on the object.
    Basic Python types like ``int`` are converted directly. For other
    types, the string representation is used.

    EXAMPLES::

        sage: pari([2,3,5])
        [2, 3, 5]
        sage: pari(Matrix(2,2,range(4)))
        [0, 1; 2, 3]
        sage: pari(x^2-3)
        x^2 - 3

    ::

        sage: a = pari(1); a, a.type()
        (1, 't_INT')
        sage: a = pari(1/2); a, a.type()
        (1/2, 't_FRAC')
        sage: a = pari(1/2); a, a.type()
        (1/2, 't_FRAC')

    Conversion from reals uses the real's own precision::

        sage: a = pari(1.2); a, a.type(), a.precision()
        (1.20000000000000, 't_REAL', 4) # 32-bit
        (1.20000000000000, 't_REAL', 3) # 64-bit

    Conversion from strings uses the current PARI real precision.
    By default, this is 64 bits::

        sage: a = pari('1.2'); a, a.type(), a.precision()
        (1.20000000000000, 't_REAL', 4)  # 32-bit
        (1.20000000000000, 't_REAL', 3)  # 64-bit

    But we can change this precision::

        sage: pari.set_real_precision(35)  # precision in decimal digits
        15
        sage: a = pari('1.2'); a, a.type(), a.precision()
        (1.2000000000000000000000000000000000, 't_REAL', 6)  # 32-bit
        (1.2000000000000000000000000000000000, 't_REAL', 4)  # 64-bit

    Set the precision to 15 digits for the remaining tests::

        sage: pari.set_real_precision(15)
        35

    Conversion from matrices and vectors is supported::

        sage: a = pari(matrix(2,3,[1,2,3,4,5,6])); a, a.type()
        ([1, 2, 3; 4, 5, 6], 't_MAT')
        sage: v = vector([1.2, 3.4, 5.6])
        sage: pari(v)
        [1.20000000000000, 3.40000000000000, 5.60000000000000]

    Some more exotic examples::

        sage: K.<a> = NumberField(x^3 - 2)
        sage: pari(K)
        [y^3 - 2, [1, 1], -108, 1, [[1, 1.25992104989487, 1.58740105196820; 1, -0.629960524947437 + 1.09112363597172*I, -0.793700525984100 - 1.37472963699860*I], [1, 1.25992104989487, 1.58740105196820; 1, 0.461163111024285, -2.16843016298270; 1, -1.72108416091916, 0.581029111014503], [1, 1, 2; 1, 0, -2; 1, -2, 1], [3, 0, 0; 0, 0, 6; 0, 6, 0], [6, 0, 0; 0, 6, 0; 0, 0, 3], [2, 0, 0; 0, 0, 1; 0, 1, 0], [2, [0, 0, 2; 1, 0, 0; 0, 1, 0]], []], [1.25992104989487, -0.629960524947437 + 1.09112363597172*I], [1, y, y^2], [1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 0, 0, 0, 2, 0, 2, 0; 0, 1, 0, 1, 0, 0, 0, 0, 2; 0, 0, 1, 0, 1, 0, 1, 0, 0]]

        sage: E = EllipticCurve('37a1')
        sage: pari(E)
        [0, 0, 1, -1, 0, 0, -2, 1, -1, 48, -216, 37, 110592/37, Vecsmall([1]), [Vecsmall([64, 1])], [0, 0, 0, 0, 0, 0, 0, 0]]

    Conversion from basic Python types::

        sage: pari(int(-5))
        -5
        sage: pari(long(2**150))
        1427247692705959881058285969449495136382746624
        sage: pari(float(pi))
        3.14159265358979
        sage: pari(complex(exp(pi*I/4)))
        0.707106781186548 + 0.707106781186548*I
        sage: pari(False)
        0
        sage: pari(True)
        1

    Some commands are just executed without returning a value::

        sage: pari("dummy = 0; kill(dummy)")
        sage: type(pari("dummy = 0; kill(dummy)"))
        <type 'NoneType'>

    TESTS::

        sage: pari(None)
        Traceback (most recent call last):
        ...
        ValueError: Cannot convert None to pari
    """
    cdef GEN g
    cdef Py_ssize_t length, i
    cdef gen v

    if isinstance(s, gen):
        return s
    try:
        return s._pari_()
    except AttributeError:
        pass

    # Check basic Python types. Start with strings, which are a very
    # common case.
    if isinstance(s, str):
        sig_on()
        g = gp_read_str(PyString_AsString(s))
        if g == gnil:
            P.clear_stack()
            return None
        return P.new_gen(g)
    # This generates slightly more efficient code than
    # isinstance(s, (int, long))
    if PyInt_Check(s) | PyLong_Check(s):
        return integer_to_gen(s)
    if isinstance(s, bool):
        return P.PARI_ONE if s else P.PARI_ZERO
    if isinstance(s, float):
        sig_on()
        return P.new_gen(dbltor(PyFloat_AS_DOUBLE(s)))
    if isinstance(s, complex):
        sig_on()
        g = cgetg(3, t_COMPLEX)
        set_gel(g, 1, dbltor(PyComplex_RealAsDouble(s)))
        set_gel(g, 2, dbltor(PyComplex_ImagAsDouble(s)))
        return P.new_gen(g)

    if isinstance(s, (list, types.XRangeType,
                        tuple, types.GeneratorType)):
        length = len(s)
        v = P._empty_vector(length)
        for i from 0 <= i < length:
            v[i] = objtogen(s[i])
        return v

    if callable(s):
        return objtoclosure(s)

    if s is None:
        raise ValueError("Cannot convert None to pari")

    # Simply use the string representation
    return objtogen(str(s))


cpdef gentoobj(gen z, locals={}):
    """
    Convert a PARI gen to a Sage/Python object.

    See the ``python`` method of :class:`gen` for documentation and
    examples.
    """
    cdef GEN g = z.g
    cdef long t = typ(g)
    cdef long tx, ty
    cdef gen real, imag
    cdef Py_ssize_t i, j, nr, nc

    if t == t_INT:
         return Integer(z)
    elif t == t_FRAC:
         return Rational(z)
    elif t == t_REAL:
        from sage.rings.all import RealField
        prec = prec_words_to_bits(z.precision())
        return RealField(prec)(z)
    elif t == t_COMPLEX:
        real = z.real()
        imag = z.imag()
        tx = typ(real.g)
        ty = typ(imag.g)
        if tx in [t_INTMOD, t_PADIC] or ty in [t_INTMOD, t_PADIC]:
            raise NotImplementedError("No conversion to python available for t_COMPLEX with t_INTMOD or t_PADIC components")
        if tx == t_REAL or ty == t_REAL:
            xprec = real.precision()  # will be 0 if exact
            yprec = imag.precision()  # will be 0 if exact
            if xprec == 0:
                prec = prec_words_to_bits(yprec)
            elif yprec == 0:
                prec = prec_words_to_bits(xprec)
            else:
                prec = max(prec_words_to_bits(xprec), prec_words_to_bits(yprec))

            from sage.rings.all import RealField, ComplexField
            R = RealField(prec)
            C = ComplexField(prec)
            return C(R(real), R(imag))
        else:
            from sage.rings.all import QuadraticField
            K = QuadraticField(-1, 'i')
            return K([gentoobj(real), gentoobj(imag)])
    elif t == t_VEC or t == t_COL:
        return [gentoobj(x, locals) for x in z.python_list()]
    elif t == t_VECSMALL:
        return z.python_list_small()
    elif t == t_MAT:
        nc = lg(g)-1
        nr = 0 if nc == 0 else lg(gel(g,1))-1
        L = [gentoobj(z[i,j], locals) for i in range(nr) for j in range(nc)]
        from sage.matrix.constructor import matrix
        return matrix(nr, nc, L)
    elif t == t_PADIC:
        from sage.rings.padics.factory import Qp
        p = z.padicprime()
        K = Qp(Integer(p), precp(g))
        return K(z.lift())
    elif t == t_INFINITY:
        if inf_get_sign(g) >= 0:
            return Infinity
        else:
            return -Infinity
    
    # Fallback (e.g. polynomials): use string representation
    from sage.misc.sage_eval import sage_eval
    return sage_eval(str(z), locals=locals)


cdef GEN _Vec_append(GEN v, GEN a, long n):
    """
    This implements appending zeros (or another constant GEN ``a``) to
    the result of :meth:`Vec` and similar functions.

    This is a shallow function, copying ``a`` and entries of ``v`` to
    the result.  The result is simply stored on the PARI stack.

    INPUT:

    - ``v`` -- GEN of type ``t_VEC`` or ``t_COL``

    - ``a`` -- GEN which will be used for the added entries.
      Normally, this would be ``gen_0``.

    - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
      append zeros; if `n < 0`, prepend zeros.

    OUTPUT:

    A GEN of the same type as ``v``.
    """
    cdef long lenv = lg(v)-1
    cdef GEN w
    cdef long i
    # Do we need to extend the vector with zeros?
    if n > lenv:
        w = cgetg(n+1, typ(v))
        for i from 1 <= i <= lenv:
            set_gel(w, i, gel(v, i))
        for i from 1 <= i <= n-lenv:
            set_gel(w, i+lenv, a)
        return w
    elif n < -lenv:
        n = -n  # Make n positive
        w = cgetg(n+1, typ(v))
        for i from 1 <= i <= lenv:
            set_gel(w, i+(n-lenv), gel(v, i))
        for i from 1 <= i <= n-lenv:
            set_gel(w, i, a)
        return w
    else:
        return v
