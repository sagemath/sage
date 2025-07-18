# sage.doctest: optional - primecountpy
"""
Utility classes for multi-modular algorithms
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_allocarray, check_reallocarray, sig_free

from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer, smallInteger
from sage.arith.misc import random_prime
from types import GeneratorType
from sage.ext.stdsage cimport PY_NEW
from cpython.object cimport PyObject_RichCompare

# should I have mod_int versions of these functions?
# c_inverse_mod_longlong modular inverse used exactly once in _refresh_precomputations
from sage.rings.fast_arith cimport arith_llong
cdef arith_llong ai
ai = arith_llong()

# This is the maximum modulus for the code in this module, i.e. the
# largest prime that can be used as modulus is
# previous_prime(MAX_MODULUS)
MAX_MODULUS = MOD_INT_MAX


cdef class MultiModularBasis_base():
    r"""
    This class stores a list of machine-sized prime numbers,
    and can do reduction and Chinese Remainder Theorem lifting
    modulo these primes.

    Lifting implemented via Garner's algorithm, which has the advantage
    that all reductions are word-sized. For each `i`, precompute
    `\prod_j=1^{i-1} m_j` and `\prod_j=1^{i-1} m_j^{-1} (mod m_i)`.

    This class can be initialized in two ways, either with a list of prime
    moduli or an upper bound for the product of the prime moduli. The prime
    moduli are generated automatically in the second case.

    EXAMPLES::

        sage: from sage.arith.multi_modular import MultiModularBasis_base
        sage: mm = MultiModularBasis_base([3, 5, 7]); mm
        MultiModularBasis with moduli [3, 5, 7]

        sage: height = 52348798724
        sage: mm = MultiModularBasis_base(height); mm
        MultiModularBasis with moduli [...]
        sage: mm.prod() >= 2*height
        True

    TESTS::

        sage: mm = MultiModularBasis_base((3,5,7)); mm
        MultiModularBasis with moduli [3, 5, 7]
        sage: mm = MultiModularBasis_base(primes(10,20)); mm
        MultiModularBasis with moduli [11, 13, 17, 19]

    There is no overflow if the modulus is below ``MAX_MODULUS``::

        sage: from sage.arith.multi_modular import MAX_MODULUS
        sage: p0 = previous_prime(MAX_MODULUS)
        sage: p1 = previous_prime(p0)
        sage: MultiModularBasis_base([p0, p1]).crt([p0-1, p1-1])
        -1

    If we add another bit to the prime length then there is an
    overflow, as expected::

        sage: p0 = previous_prime(2*MAX_MODULUS)
        sage: p1 = previous_prime(p0)
        sage: MultiModularBasis_base([p0, p1]).crt([p0-1, p1-1])
        Traceback (most recent call last):
        ...
        OverflowError: given modulus 6074000981 is larger than 3037000498
    """
    def __cinit__(self):
        r"""
        Allocate the space for the moduli and precomputation lists
        and initialize the first element of that list.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([1099511627791])
            Traceback (most recent call last):
            ...
            OverflowError: given modulus 1099511627791 is larger than 3037000498
        """
        mpz_init(self.product)
        mpz_init(self.half_product)

    cdef _realloc_to_new_count(self, new_count):
        self.moduli = <mod_int*>check_reallocarray(self.moduli, new_count, sizeof(mod_int))
        self.partial_products = <mpz_t*>check_reallocarray(self.partial_products, new_count, sizeof(mpz_t))
        self.C = <mod_int*>check_reallocarray(self.C, new_count, sizeof(mod_int))

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base(1099511627791); mm
            MultiModularBasis with moduli [...]
            sage: del mm
        """
        sig_free(self.moduli)
        for i in range(self.n):
            mpz_clear(self.partial_products[i])
        sig_free(self.partial_products)
        sig_free(self.C)
        mpz_clear(self.product)
        mpz_clear(self.half_product)

    def __init__(self, val, unsigned long l_bound=2**10, unsigned long u_bound=2**15):
        r"""
        Initialize a multi-modular basis and perform precomputations.

        INPUT:

        - ``val`` -- as integer
                        determines how many primes are computed
                        (their product will be at least 2*val)
                    as list, tuple or generator
                        a list of prime moduli to start with
        - ``l_bound`` -- integer (default: 2^10); lower bound for the random primes
        - ``u_bound`` -- integer (default: 2^15); upper bound for the random primes

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([1009, 10007]); mm
            MultiModularBasis with moduli [1009, 10007]
            sage: mm.prod()
            10097063

            sage: height = 10097063
            sage: mm = MultiModularBasis_base(height); mm
            MultiModularBasis with moduli [...]

            sage: mm.prod()//height >= 2
            True

            sage: mm = MultiModularBasis_base([1000000000000000000000000000057])
            Traceback (most recent call last):
            ...
            OverflowError: given modulus 1000000000000000000000000000057 is larger than 3037000498

            sage: mm = MultiModularBasis_base(0); mm
            MultiModularBasis with moduli [...]

            sage: mm = MultiModularBasis_base([6, 10])
            Traceback (most recent call last):
            ...
            ArithmeticError: The inverse of 6 modulo 10 is not defined.
        """
        if l_bound < 2:
            raise ValueError(f"minimum value for lower bound is 2, given: {l_bound}")
        if u_bound > MAX_MODULUS:
            raise ValueError(f"upper bound cannot be greater than {MAX_MODULUS}, given: {u_bound}")

        self._l_bound = l_bound
        self._u_bound = u_bound

        from primecountpy.primecount import prime_pi
        self._num_primes = prime_pi(self._u_bound) - prime_pi(self._l_bound-1)

        if isinstance(val, (list, tuple, GeneratorType)):
            self.extend_with_primes(val, check=True)
        else:
            self._extend_moduli_to_height(val)

    cdef mod_int _new_random_prime(self, set known_primes) except 1:
        """
        Choose a new random prime for inclusion in the list of moduli,
        or raise a :exc:`RuntimeError` if there are no more primes.

        INPUT:

        - ``known_primes`` -- Python set of already known primes in
          the allowed interval; we will not return a prime in
          known_primes
        """
        cdef mod_int p
        while True:
            if len(known_primes) >= self._num_primes:
                raise RuntimeError("there are not enough primes in the interval [%s, %s] to complete this multimodular computation" % (self._l_bound, self._u_bound))
            p = random_prime(self._u_bound, lbound=self._l_bound)
            if p not in known_primes:
                return p

    def extend_with_primes(self, plist, partial_products=None, check=True):
        """
        Extend the stored list of moduli with the given primes in ``plist``.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([1009, 10007]); mm
            MultiModularBasis with moduli [1009, 10007]
            sage: mm.extend_with_primes([10037, 10039])
            4
            sage: mm
            MultiModularBasis with moduli [1009, 10007, 10037, 10039]
        """
        if isinstance(plist, GeneratorType):
            plist = list(plist)
        elif not isinstance(plist, (tuple, list)):
            raise TypeError("plist should be a list, tuple or a generator, got: %s" % (str(type(plist))))

        cdef Py_ssize_t len_plist = len(plist)

        if len_plist == 0:
            return self.n
        if check:
            for p in plist:
                if p > MAX_MODULUS:
                    raise OverflowError(f"given modulus {p} is larger than {MAX_MODULUS}")
        self._realloc_to_new_count(self.n + len_plist)

        cdef Py_ssize_t i
        for i in range(len_plist):
            self.moduli[self.n+i] = plist[i]
            mpz_init(self.partial_products[self.n + i])
            if partial_products:
                mpz_set(self.partial_products[self.n + i],
                        (<Integer>partial_products[i]).value)

        cdef int old_count = self.n
        self.n += len_plist
        if not partial_products:
            self._refresh_products(old_count)
        else:
            self._refresh_prod()
        self._refresh_precomputations(old_count)
        return self.n

    def __richcmp__(self, other, int op):
        """
        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007])
            sage: nn = MultiModularBasis_base([10007])
            sage: mm == nn
            True
        """
        if not isinstance(other, MultiModularBasis_base):
            return NotImplemented
        left = self.__getstate__()
        right = other.__getstate__()
        return PyObject_RichCompare(left, right, op)

    def __setstate__(self, state):
        """
        Initialize a new :class:`MultiModularBasis_base` object from a
        state stored in a pickle.

        TESTS::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007, 10009])
            sage: mm == loads(dumps(mm))
            True

            sage: mm = MultiModularBasis_base([])
            sage: mm.__setstate__(([10007, 10009], 2^10, 2^15))

            sage: mm
            MultiModularBasis with moduli [10007, 10009]
        """
        nmoduli, lbound, ubound = state
        self._realloc_to_new_count(len(nmoduli))
        self._l_bound = lbound
        self._u_bound = ubound
        self.extend_with_primes(nmoduli, check=False)

    def __getstate__(self):
        """
        Return a tuple describing the state of this object for pickling.

        TESTS::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007, 10009])
            sage: mm.__getstate__()
            ([10007, 10009], 1024, 32768)
        """
        return (self.list(), self._l_bound, self._u_bound)

    def _extend_moduli_to_height(self, height):
        """
        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base(0); mm
            MultiModularBasis with moduli [...]
            sage: p = mm[0]

            sage: mm._extend_moduli_to_height(70000)
            sage: mm
            MultiModularBasis with moduli [...]
            sage: p == mm[0]
            True
            sage: mm.prod() >= 2*70000
            True

            sage: mm = MultiModularBasis_base([46307]); mm
            MultiModularBasis with moduli [46307]

            sage: mm._extend_moduli_to_height(10^30); mm
            MultiModularBasis with moduli [...]
            sage: mm.prod() >= 2*10^30
            True

        TESTS:

        Verify that :issue:`11358` is fixed::

            sage: set_random_seed(0); m = sage.arith.multi_modular.MultiModularBasis_base(0)
            sage: m._extend_moduli_to_height(prod(prime_range(50)))
            sage: m = sage.arith.multi_modular.MultiModularBasis_base([],2,100)
            sage: m._extend_moduli_to_height(prod(prime_range(90)))
            sage: m._extend_moduli_to_height(prod(prime_range(150)))
            Traceback (most recent call last):
            ...
            RuntimeError: there are not enough primes in the interval [2, 100] to complete this multimodular computation

        Another check (which fails horribly before :issue:`11358` is fixed)::

            sage: set_random_seed(0); m = sage.arith.multi_modular.MultiModularBasis_base(0); m._extend_moduli_to_height(10**10000)
            sage: len(set(m)) == len(m)
            True
            sage: len(m)
            2440
        """
        cdef Integer h = Integer(height)
        if h < self._l_bound:
            h = Integer(self._l_bound)
        self._extend_moduli_to_height_c(h.value)

    cdef int _extend_moduli_to_height_c(self, mpz_t height0) except -1:
        r"""
        Expand the list of primes and perform precomputations.

        INPUT:

        - ``height`` -- determines how many primes are computed
                       (their product must be at least 2*height)
        """
        # real height we use is twice the given, set height to 2*height0
        cdef mpz_t height
        mpz_init(height)
        mpz_mul_2exp(height, height0, 1)
        # check if we already have enough prime moduli
        if self.n > 0 and mpz_cmp(height, self.partial_products[self.n-1]) <= 0:
            mpz_clear(height)
            return self.n

        # find new prime moduli
        new_moduli = []
        new_partial_products = []
        cdef Integer M  # keeps current height
        cdef mod_int p  # keeps current prime moduli

        if self.n == 0:
            M = smallInteger(1)
        else:
            M = PY_NEW(Integer)
            mpz_set(M.value, self.partial_products[self.n-1])

        known_primes = set(self)
        while mpz_cmp(height, M.value) > 0:
            p = self._new_random_prime(known_primes)
            new_moduli.append(p)
            known_primes.add(p)
            M *= p
            new_partial_products.append(M)
        mpz_clear(height)
        return self.extend_with_primes(new_moduli, new_partial_products,
                                       check=False)

    def _extend_moduli_to_count(self, int count):
        r"""
        Expand the list of primes and perform precomputations.

        INPUT:

        - ``count`` -- the minimum number of moduli in the resulting list

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([46307]); mm
            MultiModularBasis with moduli [46307]
            sage: mm._extend_moduli_to_count(3)
            3
            sage: mm
            MultiModularBasis with moduli [...]
            sage: len(mm)
            3
        """
        if count <= self.n:
            return self.n
        new_moduli = []

        cdef int i
        cdef mod_int p
        known_primes = set(self)
        for i in range(self.n, count):
            p = self._new_random_prime(known_primes)
            known_primes.add(p)
            new_moduli.append(p)

        return self.extend_with_primes(new_moduli, check=False)

    def _extend_moduli(self, count):
        """
        Expand the list of prime moduli with `count` new random primes.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([46307]); mm
            MultiModularBasis with moduli [46307]
            sage: mm._extend_moduli(2); mm  # random
            MultiModularBasis with moduli [46307, 31051, 16981]
            sage: mm[0]
            46307
            sage: all(p.is_prime() for p in mm)
            True
        """
        self._extend_moduli_to_count(self.n + count)

    cdef void _refresh_products(self, int start) noexcept:
        r"""
        Compute and store `\prod_j=1^{i-1} m_j` for i > start.
        """
        cdef mpz_t z
        mpz_init(z)
        if start == 0:
            mpz_set_si(self.partial_products[0], self.moduli[0])
            start += 1
        for i in range(start, self.n):
            mpz_set_si(z, self.moduli[i])
            mpz_mul(self.partial_products[i], self.partial_products[i-1], z)
        mpz_clear(z)
        self._refresh_prod()

    cdef void _refresh_prod(self) noexcept:
        # record the product and half product for balancing the lifts.
        mpz_set(self.product, self.partial_products[self.n-1])
        mpz_fdiv_q_ui(self.half_product, self.product, 2)

    cdef void _refresh_precomputations(self, int start) except *:
        r"""
        Compute and store `\prod_j=1^{i-1} m_j^{-1} (mod m_i)` for i >= start.
        """
        if start == 0:
            start = 1  # first one is trivial, never used
            self.C[0] = 1
        for i in range(start, self.n):
            self.C[i] = ai.c_inverse_mod_longlong(mpz_fdiv_ui(self.partial_products[i-1], self.moduli[i]), self.moduli[i])

    cdef int min_moduli_count(self, mpz_t height) except -1:
        r"""
        Compute the minimum number of primes needed to uniquely determine
        an integer mod height.
        """
        self._extend_moduli_to_height_c(height)

        cdef int count
        count = self.n * mpz_sizeinbase(height, 2) / mpz_sizeinbase(self.partial_products[self.n-1], 2)  # an estimate
        count = max(min(count, self.n), 1)
        while count > 1 and mpz_cmp(height, self.partial_products[count-1]) < 0:
            count -= 1
        while mpz_cmp(height, self.partial_products[count-1]) > 0:
            count += 1

        return count

    cdef mod_int last_prime(self) noexcept:
        return self.moduli[self.n-1]

    cdef int mpz_reduce_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1:
        r"""
        Perform reduction mod `m_i` for offset <= i < len.

        `b[i] = z mod m_{i+offset}` for 0 <= i < len

        INPUT:

        - ``z`` -- the integer being reduced
        - ``b`` -- array to hold the reductions mod each `m_i`;
          it *must* be allocated and have length at least len
        - ``offset`` -- first prime in list to reduce against
        - ``len`` -- number of primes in list to reduce against
        """
        cdef int i
        cdef mod_int* m
        m = self.moduli + offset
        for i in range(len):
            b[i] = mpz_fdiv_ui(z, m[i])
        return 0

    cdef int mpz_reduce_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1:
        r"""
        Perform reduction mod `m_i` for offset <= i < len.

        `b[i][j] = z[j] mod m_{i+offset}` for 0 <= i < len

        INPUT:

        - ``z`` -- an array of integers being reduced
        - ``b`` -- array to hold the reductions mod each m_i.
          It *must* be fully allocated and each have length at least len.
        - ``vn`` -- length of z and each b[i]
        - ``offset`` -- first prime in list to reduce against
        - ``len`` -- number of primes in list to reduce against
        """
        cdef int i, j
        cdef mod_int* m
        m = self.moduli + offset
        for i in range(len):
            mi = m[i]
            for j in range(vn):
                b[i][j] = mpz_fdiv_ui(z[j], mi)
        return 0

    cdef int mpz_crt_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1:
        r"""
        Calculate lift mod `\prod_{i=0}^{offset+len-1} m_i`.

        z = b[i] mod `m_{i+offset}` for 0 <= i < len

        In the case that offset > 0,
        z remains unchanged mod `\prod_{i=0}^{offset-1} m_i`

        INPUT:

        - ``z`` -- a placeholder for the constructed integer; *must* be
          initialized if and only if ``offset > 0``
        - ``b`` -- array holding the reductions mod each m_i; *must* have
          length at least ``len``
        - ``offset`` -- first prime in list to reduce against
        - ``len`` -- number of primes in list to reduce against
        """
        cdef int i, s
        cdef mpz_t u
        cdef mod_int* m
        m = self.moduli + offset
        mpz_init(u)
        if offset == 0:
            s = 1
            mpz_init_set_si(z, b[0])
            if b[0] == 0:
                while s < len and b[s] == 0:  # fast forward to first nonzero
                    s += 1
        else:
            s = 0
        for i in range(s, len):
            mpz_set_si(u, ((b[i] + m[i] - mpz_fdiv_ui(z, m[i])) * self.C[i]) % m[i])
            mpz_mul(u, u, self.partial_products[i-1])
            mpz_add(z, z, u)

        # normalize to be between -prod/2 and prod/2.
        if mpz_cmp(z, self.half_product) > 0:
            mpz_sub(z, z, self.product)
        mpz_clear(u)
        return 0

    cdef int mpz_crt_vec_tail(self, mpz_t* z, mod_int** b, int vc, int offset, int len) except -1:
        r"""
        Calculate lift mod `\prod_{i=0}^{offset+len-1} m_i`.

        `z[j] = b[i][j] mod m_{i+offset}` for 0 <= i < len

        In the case that offset > 0,
        z[j] remains unchanged mod `\prod_{i=0}^{offset-1} m_i`

        INPUT:

        - ``z`` -- a placeholder for the constructed integers; *must* be
          allocated and have length at least ``vc``. z[j] *must* be initialized
          if and only if ``offset > 0``.
        - ``b`` -- array holding the reductions mod each m_i. *must* have
          length at least ``len``.
        - ``vn`` -- length of z and each b[i]
        - ``offset`` -- first prime in list to reduce against
        - ``len`` -- number of primes in list to reduce against
        """
        cdef int i, j
        cdef mpz_t u
        cdef mod_int* m

        m = self.moduli + offset
        mpz_init(u)
        if offset == 0:
            s = 1
        else:
            s = 0

        for j in range(vc):
            i = s
            if offset == 0:
                mpz_set_si(z[j], b[0][j])
                if b[0][j] == 0:
                    while i < len and b[i][j] == 0:  # fast forward to first nonzero
                        i += 1
            while i < len:
                mpz_set_si(u, ((b[i][j] + m[i] - mpz_fdiv_ui(z[j], m[i])) * self.C[i]) % m[i])  # u = ((b_i - z) * C_i) % m_i
                mpz_mul(u, u, self.partial_products[i-1])
                mpz_add(z[j], z[j], u)
                i += 1

            # normalize to be between -prod/2 and prod/2.
            if mpz_cmp(z[j], self.half_product) > 0:
                mpz_sub(z[j], z[j], self.product)

        cdef Integer zz
        zz = PY_NEW(Integer)
        mpz_set(zz.value, self.half_product)

        mpz_clear(u)
        return 0

    def crt(self, b):
        r"""
        Calculate lift mod `\prod_{i=0}^{len(b)-1} m_i`.

        In the case that offset > 0,
        z[j] remains unchanged mod `\prod_{i=0}^{offset-1} m_i`

        INPUT:

        - ``b`` -- list of length at most self.n

        OUTPUT:

            Integer z where `z = b[i] mod m_i` for 0 <= i < len(b)

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007, 10009, 10037, 10039, 17351])
            sage: res = mm.crt([3,5,7,9]); res
            8474803647063985
            sage: res % 10007
            3
            sage: res % 10009
            5
            sage: res % 10037
            7
            sage: res % 10039
            9
        """
        cdef int i, n
        n = len(b)
        if n > self.n:
            raise IndexError("beyond bound for multi-modular prime list")
        cdef mod_int* bs
        bs = <mod_int*>check_allocarray(n, sizeof(mod_int))
        for i in range(n):
            bs[i] = b[i]
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_crt_tail(z.value, bs, 0, n)
        sig_free(bs)
        return z

    def precomputation_list(self):
        r"""
        Return a list of the precomputed coefficients
        `\prod_j=1^{i-1} m_j^{-1} (mod m_i)`
        where `m_i` are the prime moduli.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([46307, 10007]); mm
            MultiModularBasis with moduli [46307, 10007]
            sage: mm.precomputation_list()
            [1, 4013]
        """
        return [Integer(self.C[i]) for i in range(self.n)]

    def partial_product(self, n):
        """
        Return a list containing precomputed partial products.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([46307, 10007]); mm
            MultiModularBasis with moduli [46307, 10007]
            sage: mm.partial_product(0)
            46307
            sage: mm.partial_product(1)
            463394149

        TESTS::

            sage: mm.partial_product(2)
            Traceback (most recent call last):
            ...
            IndexError: beyond bound for multi-modular prime list
            sage: mm.partial_product(-2)
            Traceback (most recent call last):
            ...
            IndexError: negative index not valid
        """
        if n >= self.n:
            raise IndexError("beyond bound for multi-modular prime list")
        if n < 0:
            raise IndexError("negative index not valid")
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_set(z.value, self.partial_products[n])
        return z

    def prod(self):
        """
        Return the product of the prime moduli.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([46307]); mm
            MultiModularBasis with moduli [46307]
            sage: mm.prod()
            46307
            sage: mm = MultiModularBasis_base([46307, 10007]); mm
            MultiModularBasis with moduli [46307, 10007]
            sage: mm.prod()
            463394149

        TESTS::

            sage: mm = MultiModularBasis_base([]); mm
            MultiModularBasis with moduli []
            sage: len(mm)
            0
            sage: mm.prod()
            1
        """
        if self.n == 0:
            return 1
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_set(z.value, self.partial_products[self.n-1])
        return z

    def list(self):
        """
        Return a list with the prime moduli.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([46307, 10007])
            sage: mm.list()
            [46307, 10007]
        """
        return [Integer(self.moduli[i]) for i in range(self.n)]

    def __len__(self):
        """
        Return the number of moduli stored.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007])
            sage: len(mm)
            1
            sage: mm._extend_moduli_to_count(2)
            2
            sage: len(mm)
            2
        """
        return self.n

    def __iter__(self):
        """
        Return an iterator over the prime moduli.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007, 10009])
            sage: t = iter(mm); t
            <list...iterator object at ...>
            sage: list(mm.__iter__())
            [10007, 10009]
        """
        return iter(self.list())

    def __getitem__(self, ix):
        """
        Return the moduli stored at index `ix` as a Python long.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: mm = MultiModularBasis_base([10007, 10009])
            sage: mm[1]
            10009
            sage: mm[-1]
            Traceback (most recent call last):
            ...
            IndexError: index out of range

            sage: mm[:1]
            MultiModularBasis with moduli [10007]
        """
        if isinstance(ix, slice):
            return self.__class__(self.list()[ix], l_bound = self._l_bound,
                                  u_bound = self._u_bound)

        cdef Py_ssize_t i = ix
        if i != ix:
            raise TypeError("index must be an integer")
        if i < 0 or i >= self.n:
            raise IndexError("index out of range")
        return self.moduli[i]

    def __repr__(self):
        """
        Return a string representation of this object.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MultiModularBasis_base
            sage: MultiModularBasis_base([10007])
            MultiModularBasis with moduli [10007]
        """
        return "MultiModularBasis with moduli " + str(self.list())


cdef class MultiModularBasis(MultiModularBasis_base):
    """
    Class used for storing a MultiModular bases of a fixed length.
    """
    cdef int mpz_reduce(self, mpz_t z, mod_int* b) except -1:
        r"""
        Perform reduction mod `m_i` for each modulus `m_i`.

        `b[i] = z mod m_i` for 0 <= i < len(self)

        INPUT:

        - ``z`` -- the integer being reduced
        - ``b`` -- array to hold the reductions mod each `m_i`; it *must* be
          allocated and have length at least len
        """
        self.mpz_reduce_tail(z, b, 0, self.n)

    cdef int mpz_reduce_vec(self, mpz_t* z, mod_int** b, int vn) except -1:
        r"""
        Perform reduction mod `m_i` for each modulus `m_i`.

        `b[i][j] = z[j] mod m_i` for 0 <= i < len(self)

        INPUT:

        - ``z`` -- an array of integers being reduced
        - ``b`` -- array to hold the reductions mod each `m_i`; it *must* be
          fully allocated and each have length at least len
        - ``vn`` -- length of ``z`` and each ``b[i]``
        """
        self.mpz_reduce_vec_tail(z, b, vn, 0, self.n)

    cdef int mpz_crt(self, mpz_t z, mod_int* b) except -1:
        r"""
        Calculate lift mod `\prod m_i`.

        `z = b[i] mod m_{i+offset}` for 0 <= i < len(self)

        INPUT:

        - ``z`` -- a placeholder for the constructed integer; *must not* be
          initialized
        - ``b`` -- array holding the reductions mod each `m_i`; it *must* have
          length at least ``len(self)``
        """
        self.mpz_crt_tail(z, b, 0, self.n)

    cdef int mpz_crt_vec(self, mpz_t* z, mod_int** b, int vn) except -1:
        r"""
        Calculate lift mod `\prod m_i`.

        `z[j] = b[i][j] mod m_i` for 0 <= i < len(self)

        INPUT:

        - ``z`` -- a placeholder for the constructed integers; *must* be
          allocated and have length at least vn,but each z[j] *must not* be
          initialized
        - ``b`` -- array holding the reductions mod each `m_i`; it *must* have
          length at least ``len(self)``
        - ``vn`` -- length of ``z`` and each ``b[i]``
        """
        self.mpz_crt_vec_tail(z, b, vn, 0, self.n)


cdef class MutableMultiModularBasis(MultiModularBasis):
    """
    Class used for performing multi-modular methods,
    with the possibility of removing bad primes.
    """
    cpdef mod_int next_prime(self) except -1:
        """
        Pick a new random prime between the bounds given during the
        initialization of this object, update the precomputed data,
        and return the new prime modulus.

        EXAMPLES::

            sage: from sage.arith.multi_modular import MutableMultiModularBasis
            sage: mm = MutableMultiModularBasis([10007])
            sage: p = mm.next_prime()
            sage: 1024 < p < 32768
            True
            sage: p != 10007
            True
            sage: mm.list() == [10007, p]
            True
        """
        self._extend_moduli(1)
        return self.moduli[self.n-1]

    cpdef mod_int replace_prime(self, int ix) except -1:
        """
        Replace the prime moduli at the given index with a different one,
        update the precomputed data accordingly, and return the new prime.

        INPUT:

        - ``ix`` -- index into list of moduli

        OUTPUT: the new prime modulus

        EXAMPLES::

            sage: from sage.arith.multi_modular import MutableMultiModularBasis
            sage: mm = MutableMultiModularBasis([10007, 10009, 10037, 10039])
            sage: mm
            MultiModularBasis with moduli [10007, 10009, 10037, 10039]
            sage: prev_prod = mm.prod(); prev_prod
            10092272478850909
            sage: mm.precomputation_list()
            [1, 5004, 6536, 6060]
            sage: mm.partial_product(2)
            1005306552331
            sage: p = mm.replace_prime(1)
            sage: mm.list() == [10007, p, 10037, 10039]
            True
            sage: mm.prod()*10009 == prev_prod*p
            True
            sage: precomputed = mm.precomputation_list()
            sage: precomputed == [prod(Integers(mm[i])(1 / mm[j])
            ....:                      for j in range(i))
            ....:                 for i in range(4)]
            True
            sage: mm.partial_product(2) == prod(mm.list()[:3])
            True
        """
        cdef mod_int new_p

        if ix < 0 or ix >= self.n:
            raise IndexError("index out of range")

        new_p = self._new_random_prime(set(self))
        self.moduli[ix] = new_p

        self._refresh_products(ix)
        self._refresh_precomputations(ix)
        return new_p
