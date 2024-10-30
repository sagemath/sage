# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
PowComputer

A class for computing and caching powers of the same integer.

This class is designed to be used as a field of `p`-adic rings and
fields.  Since elements of `p`-adic rings and fields need to use powers
of p over and over, this class precomputes and stores powers of p.
There is no reason that the base has to be prime however.

EXAMPLES::

    sage: X = PowComputer(3, 4, 10)
    sage: X(3)
    27
    sage: X(10) == 3^10
    True

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import weakref
from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.rings.infinity import infinity
from sage.libs.gmp.mpz cimport *
from sage.structure.richcmp cimport richcmp_not_equal, richcmp
from cpython.object cimport Py_EQ, Py_NE

from sage.ext.stdsage cimport PY_NEW

cdef long maxpreccap = (1L << (sizeof(long) * 8 - 2)) - 1

cdef class PowComputer_class(SageObject):
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Memory allocation.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC.pow_Integer_Integer(2)
            9
        """
        sig_on()
        mpz_init(self.temp_m)
        sig_off()
        self._allocated = 1

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Initialize ``self``.

        INPUT:

        - ``prime`` -- the prime that is the base of the exponentials
          stored in this ``pow_computer``

        - ``cache_limit`` -- how high to cache powers of prime

        - ``prec_cap`` -- data stored for `p`-adic elements using this
          ``pow_computer`` (so they have C-level access to fields
          common to all elements of the same parent)

        - ``ram_prec_cap`` -- prec_cap * e

        - ``in_field`` -- same idea as prec_cap

        - ``poly`` -- same idea as prec_cap

        - ``shift_seed`` -- same idea as prec_cap

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC.pow_Integer_Integer(2)
            9
        """
        self.prime = prime
        self.p2 = prime // 2
        self.in_field = in_field
        self.cache_limit = cache_limit
        self.prec_cap = prec_cap
        self.ram_prec_cap = ram_prec_cap

    def __richcmp__(self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: P = PowComputer(3, 4, 9)
            sage: P == 7
            False
            sage: Q = PowComputer(3, 6, 9)
            sage: P == Q
            False
            sage: Q = PowComputer(3, 4, 9)
            sage: P == Q
            True
            sage: P is Q
            True
        """
        if not isinstance(other, PowComputer_class):
            if op in [Py_EQ, Py_NE]:
                return (op == Py_NE)
            return NotImplemented

        cdef PowComputer_class s = self
        cdef PowComputer_class o = other

        lx = s.prime
        rx = o.prime
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = s.prec_cap
        rx = o.prec_cap
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = s.cache_limit
        rx = o.cache_limit
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        return richcmp(s.in_field, o.in_field, op)

    cdef Integer pow_Integer(self, long n):
        """
        Return ``self.prime^n``.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC.pow_Integer_Integer(2)  # indirect doctest
            9
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.pow_mpz_t_tmp(n))
        return ans

    def pow_Integer_Integer(self, n):
        """
        Test the ``pow_Integer`` function.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC.pow_Integer_Integer(4)
            81
            sage: PC.pow_Integer_Integer(6)
            729
            sage: PC.pow_Integer_Integer(0)
            1
            sage: PC.pow_Integer_Integer(10)
            59049

            sage: # needs sage.libs.ntl
            sage: PC = PowComputer_ext_maker(3, 5, 10, 20, False, ntl.ZZ_pX([-3,0,1], 3^10), 'big','e',ntl.ZZ_pX([1],3^10))
            sage: PC.pow_Integer_Integer(4)
            81
            sage: PC.pow_Integer_Integer(6)
            729
            sage: PC.pow_Integer_Integer(0)
            1
            sage: PC.pow_Integer_Integer(10)
            59049
        """
        cdef Integer _n = Integer(n)
        if _n < 0:
            if mpz_fits_ulong_p((<Integer>-_n).value) == 0:
                raise ValueError("result too big")
            return ~self.pow_Integer(mpz_get_ui((<Integer>-_n).value))
        else:
            if mpz_fits_ulong_p(_n.value) == 0:
                raise ValueError("result too big")
            return self.pow_Integer(mpz_get_ui(_n.value))

    cdef mpz_srcptr pow_mpz_t_tmp(self, long n) except NULL:
        """
        Provides fast access to an ``mpz_srcptr`` pointing to self.prime^n.

        The location pointed to depends on the underlying
        representation.  In no circumstances should you mpz_clear the
        result.  The value pointed to may be an internal temporary
        variable for the class.  In particular, you should not try to
        refer to the results of two pow_mpz_t_tmp calls at the same
        time, because the second call may overwrite the memory pointed
        to by the first.

        See pow_mpz_t_tmp_demo for an example of this phenomenon.
        """
        # READ THE DOCSTRING
        raise NotImplementedError

    def _pow_mpz_t_tmp_demo(self, m, n):
        """
        This function demonstrates a danger in using pow_mpz_t_tmp.

        EXAMPLES::

            sage: PC = PowComputer(5, 5, 10)

        When you call pow_mpz_t_tmp with an input that is not stored
        (ie n > self.cache_limit and n != self.prec_cap),
        it stores the result in self.temp_m and returns a pointer
        to that mpz_t.  So if you try to use the results of two
        calls at once, things will break. ::

            sage: PC._pow_mpz_t_tmp_demo(6, 8)  # 244140625 on some architectures and 152587890625 on others: random
            244140625
            sage: 5^6*5^8
            6103515625
            sage: 5^6*5^6
            244140625

        Note that this does not occur if you try a stored value,
        because the result of one of the calls points to that
        stored value. ::

            sage: PC._pow_mpz_t_tmp_demo(6, 10)
            152587890625
            sage: 5^6*5^10
            152587890625
        """
        m = Integer(m)
        n = Integer(n)
        if m < 0 or n < 0:
            raise ValueError("m, n must be nonnegative")
        cdef Integer ans = PY_NEW(Integer)
        mpz_mul(ans.value, self.pow_mpz_t_tmp(mpz_get_ui((<Integer>m).value)), self.pow_mpz_t_tmp(mpz_get_ui((<Integer>n).value)))
        return ans

    def _pow_mpz_t_tmp_test(self, n):
        """
        Test the ``pow_mpz_t_tmp`` function.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC._pow_mpz_t_tmp_test(4)
            81
            sage: PC._pow_mpz_t_tmp_test(6)
            729
            sage: PC._pow_mpz_t_tmp_test(0)
            1
            sage: PC._pow_mpz_t_tmp_test(10)
            59049

            sage: # needs sage.libs.ntl
            sage: PC = PowComputer_ext_maker(3, 5, 10, 20, False, ntl.ZZ_pX([-3,0,1], 3^10), 'big','e',ntl.ZZ_pX([1],3^10))
            sage: PC._pow_mpz_t_tmp_test(4)
            81
            sage: PC._pow_mpz_t_tmp_test(6)
            729
            sage: PC._pow_mpz_t_tmp_test(0)
            1
            sage: PC._pow_mpz_t_tmp_test(10)
            59049
        """
        cdef Integer _n = Integer(n)
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.pow_mpz_t_tmp(mpz_get_si(_n.value)))
        return ans

    cdef mpz_srcptr pow_mpz_t_top(self) noexcept:
        """
        Return a pointer to ``self.prime^self.prec_cap`` as an ``mpz_srcptr``.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC._pow_mpz_t_top_test()  # indirect doctest
            59049
        """
        raise NotImplementedError

    def _pow_mpz_t_top_test(self):
        """
        Test the ``pow_mpz_t_top`` function.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC._pow_mpz_t_top_test()
            59049

            sage: # needs sage.libs.ntl
            sage: PC = PowComputer_ext_maker(3, 5, 10, 20, False, ntl.ZZ_pX([-3,0,1], 3^10), 'big','e',ntl.ZZ_pX([1],3^10))
            sage: PC._pow_mpz_t_top_test()
            59049
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.pow_mpz_t_top())
        return ans

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10); PC
            PowComputer for 3
        """
        return "PowComputer for %s" % (self.prime)

    def _prime(self):
        """
        Return the base that the ``PowComputer`` is exponentiating.

        EXAMPLES::

            sage: P = PowComputer(6, 10, 15)
            sage: P._prime()
            6
        """
        return self.prime

    def _in_field(self):
        """
        Return whether or not ``self`` is attached to a field.

        EXAMPLES::

            sage: P = PowComputer(3, 5, 10)
            sage: P._in_field()
            False
        """
        return self.in_field

    def _cache_limit(self):
        """
        Return the limit to which powers of prime are computed.

        EXAMPLES::

            sage: P = PowComputer(3, 5, 10)
            sage: P._cache_limit()
            5
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.cache_limit)
        return ans

    def _prec_cap(self):
        """
        Return ``prec_cap``, a single value that for which
        ``self._prime()^prec_cap`` is stored.

        EXAMPLES::

            sage: P = PowComputer(3, 5, 10)
            sage: P._prec_cap()
            10
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.prec_cap)
        return ans

    def _top_power(self):
        """
        Return ``self._prime()^self._prec_cap()``.

        EXAMPLES::

            sage: P = PowComputer(3, 4, 6)
            sage: P._top_power()
            729
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.pow_mpz_t_top())
        return ans

    def __call__(self, n):
        """
        Return ``self.prime^n``.

        EXAMPLES::

            sage: P = PowComputer(3, 4, 6)
            sage: P(3)
            27
            sage: P(6)
            729
            sage: P(5)
            243
            sage: P(7)
            2187
            sage: P(0)
            1
            sage: P(-2)
            1/9
        """
        cdef Integer _n
        if n is infinity:
            return Integer(0)
        if not isinstance(n, Integer):
            _n = Integer(n)
        else:
            _n = <Integer>n
        if mpz_fits_slong_p(_n.value) == 0:
            raise ValueError("n too big")
        if _n < 0:
            return ~self.pow_Integer(-mpz_get_si(_n.value))
        else:
            return self.pow_Integer(mpz_get_ui(_n.value))


cdef class PowComputer_base(PowComputer_class):
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Allocate a ``PowComputer_base``.

        EXAMPLES::

            sage: PC = PowComputer(5, 7, 10)
            sage: PC(3)
            125
        """
        cdef Py_ssize_t i

        sig_on()
        try:
            self.small_powers = <mpz_t *>sig_malloc(sizeof(mpz_t) * (cache_limit + 1))
            if self.small_powers == NULL:
                raise MemoryError("out of memory allocating power storing")
            try:
                mpz_init(self.top_power)
                try:
                    mpz_init(self.powhelper_oneunit)
                    try:
                        mpz_init(self.powhelper_teichdiff)
                        try:
                            mpz_init(self.shift_rem)
                            try:
                                mpz_init(self.aliasing)
                                try:
                                    for i in range(cache_limit + 1):
                                        try:
                                            mpz_init(self.small_powers[i])
                                        except BaseException:
                                            while i:
                                                i-=1
                                                mpz_clear(self.small_powers[i])
                                            raise
                                except BaseException:
                                    mpz_clear(self.aliasing)
                                    raise
                            except BaseException:
                                mpz_clear(self.shift_rem)
                                raise
                        except BaseException:
                            mpz_clear(self.powhelper_teichdiff)
                            raise
                    except BaseException:
                        mpz_clear(self.powhelper_oneunit)
                        raise
                except BaseException:
                    mpz_clear(self.top_power)
                    raise
            except BaseException:
                sig_free(self.small_powers)
                raise
        finally:
            sig_off()

        self._allocated = 2

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Initialization.

        TESTS::

            sage: PC = PowComputer(5, 7, 10)
            sage: PC(3)
            125
        """
        PowComputer_class.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)

        cdef Py_ssize_t i

        mpz_set_ui(self.small_powers[0], 1)
        if cache_limit > 0:
            mpz_set(self.small_powers[1], prime.value)
        for i in range(2, cache_limit + 1):
            mpz_mul(self.small_powers[i], self.small_powers[i - 1], prime.value)
        sig_on()
        mpz_pow_ui(self.top_power, prime.value, prec_cap)
        sig_off()
        self.deg = 1
        self.e = 1
        self.f = 1
        self.ram_prec_cap = prec_cap

    def __dealloc__(self):
        """
        Deletion.

        EXAMPLES::

            sage: P = PowComputer(5, 7, 10)
            sage: del P
            sage: PowComputer(5, 7, 10)
            PowComputer for 5
        """
        cdef Py_ssize_t i

        if self._allocated >= 2:
            for i in range(self.cache_limit + 1):
                mpz_clear(self.small_powers[i])
            mpz_clear(self.top_power)
            mpz_clear(self.powhelper_oneunit)
            mpz_clear(self.powhelper_teichdiff)
            mpz_clear(self.shift_rem)
            mpz_clear(self.aliasing)
            mpz_clear(self.temp_m)
            sig_free(self.small_powers)

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: P = PowComputer(5, 7, 10)
            sage: R = loads(dumps(P))
            sage: P == R
            True
        """
        return PowComputer, (self.prime, self.cache_limit, self.prec_cap, self.in_field)

    cdef mpz_srcptr pow_mpz_t_top(self) noexcept:
        """
        Return a pointer to ``self.prime^self.prec_cap`` as an ``mpz_srcptr``.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC._pow_mpz_t_top_test()  # indirect doctest
            59049
        """
        return self.top_power

    cdef mpz_srcptr pow_mpz_t_tmp(self, long n) except NULL:
        """
        Compute ``self.prime^n``.

        EXAMPLES::

            sage: PC = PowComputer(3, 5, 10)
            sage: PC._pow_mpz_t_tmp_test(4)
            81
            sage: PC._pow_mpz_t_tmp_test(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be nonnegative
        """
        if n < 0:
            raise ValueError("n must be nonnegative")
        if n <= self.cache_limit:
            return self.small_powers[n]
        if n == self.prec_cap:
            return self.top_power
        # n may exceed self.prec_cap. Very large values can, however, lead to
        # out-of-memory situations in the following computation. This
        # sig_on()/sig_off() prevents sage from crashing in such cases.
        # It does not have a significant impact on performance. For small
        # values of n the powers are taken from self.small_powers, for large
        # values, the computation dominates the cost of the sig_on()/sig_off().
        sig_on()
        mpz_pow_ui(self.temp_m, self.prime.value, n)
        sig_off()
        return self.temp_m

pow_comp_cache = {}
cdef PowComputer_base PowComputer_c(Integer m, Integer cache_limit, Integer prec_cap, in_field, prec_type=None):
    """
    Return a ``PowComputer``.

    EXAMPLES::

        sage: PC = PowComputer(3, 5, 10)  # indirect doctest
        sage: PC(4)
        81
    """
    if cache_limit < 0:
        raise ValueError("cache_limit must be nonnegative")
    if prec_cap < 0:
        raise ValueError("prec_cap must be nonnegative")
    if mpz_cmp_si((<Integer>prec_cap).value, maxpreccap) >= 0:
        raise ValueError("cannot create p-adic parents with precision cap larger than (1 << (sizeof(long)*8 - 2))")

    key = (m, cache_limit, prec_cap, in_field, prec_type)
    if key in pow_comp_cache:
        PC = pow_comp_cache[key]()
        if PC is not None:
            return PC
    if prec_type == 'capped-rel':
        from sage.rings.padics.padic_capped_relative_element import PowComputer_ as PC_class
    elif prec_type == 'capped-abs':
        from sage.rings.padics.padic_capped_absolute_element import PowComputer_ as PC_class
    elif prec_type == 'fixed-mod':
        from sage.rings.padics.padic_fixed_mod_element import PowComputer_ as PC_class
    elif prec_type == 'floating-point':
        from sage.rings.padics.padic_floating_point_element import PowComputer_ as PC_class
    else:
        PC_class = PowComputer_base
    PC = PC_class(m, mpz_get_ui(cache_limit.value), mpz_get_ui(prec_cap.value), mpz_get_ui(prec_cap.value), in_field)
    pow_comp_cache[key] = weakref.ref(PC)
    return PC


# To speed up the creation of PowComputers with the same m, we might eventually want to copy over data from an existing PowComputer.

def PowComputer(m, cache_limit, prec_cap, in_field = False, prec_type=None):
    r"""
    Return a ``PowComputer`` that caches the values `1, m, m^2, \ldots, m^{C}`,
    where `C` is ``cache_limit``.

    Once you create a ``PowComputer``, merely call it to get values out.

    You can input any integer, even if it's outside of the precomputed range.

    INPUT:

    - ``m`` -- integer; the base that you want to exponentiate
    - ``cache_limit`` -- positive integer; that you want to cache powers up to

    EXAMPLES::

        sage: PC = PowComputer(3, 5, 10)
        sage: PC
        PowComputer for 3
        sage: PC(4)
        81
        sage: PC(6)
        729
        sage: PC(-1)
        1/3
    """
    if not isinstance(m, Integer):
        m = Integer(m)
    if not isinstance(cache_limit, Integer):
        cache_limit = Integer(cache_limit)
    if not isinstance(prec_cap, Integer):
        prec_cap = Integer(prec_cap)
    return PowComputer_c(m, cache_limit, prec_cap, in_field, prec_type)
