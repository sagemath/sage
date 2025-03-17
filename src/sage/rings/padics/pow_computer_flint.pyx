# distutils: libraries = gmp NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.libs.gmp.mpz cimport mpz_init, mpz_clear, mpz_pow_ui
from sage.libs.flint.padic cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.nmod_vec cimport *
from sage.libs.flint.fmpz_vec cimport *
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_one, fmpz_mul, fmpz_set, fmpz_get_mpz, fmpz_clear, fmpz_pow_ui, fmpz_set_mpz, fmpz_fdiv_q_2exp

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport richcmp_not_equal
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint


cdef class PowComputer_flint(PowComputer_class):
    """
    A PowComputer for use in `p`-adics implemented via FLINT.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint
        sage: PowComputer_flint(5, 20, 20, 20, False)
        FLINT PowComputer for 5
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Memory initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint
            sage: type(PowComputer_flint(5, 20, 20, 20, False))
            <class 'sage.rings.padics.pow_computer_flint.PowComputer_flint'>
        """
        # fmpz_init does not allocate memory
        fmpz_init(self.fprime)
        fmpz_init(self.half_prime)
        fmpz_set_mpz(self.fprime, prime.value)
        fmpz_init(self._fpow_variable)
        fmpz_fdiv_q_2exp(self.half_prime, self.fprime, 1)
        fmpz_init(self.tfmpz)

        sig_on()
        try:
            mpz_init(self.top_power)
            padic_ctx_init(self.ctx, self.fprime, 0, prec_cap, PADIC_SERIES)
        finally:
            sig_off()

        self._allocated = 4

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f, 'capped-rel')  # indirect doctest
            sage: TestSuite(A).run()
        """
        PowComputer_class.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)

        mpz_pow_ui(self.top_power, prime.value, prec_cap)

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint
            sage: A = PowComputer_flint(5, 20, 20, 20, False)
            sage: del A
        """
        if self._allocated >= 4:
            fmpz_clear(self.fprime)
            fmpz_clear(self.half_prime)
            fmpz_clear(self._fpow_variable)
            fmpz_clear(self.tfmpz)
            mpz_clear(self.top_power)
            padic_ctx_clear(self.ctx)

    def __reduce__(self):
        """
        Pickling.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f, 'capped-rel')  # indirect doctest
            sage: A._test_pickling()  # indirect doctest
        """
        return PowComputer_flint_maker, (self.prime, self.cache_limit, self.prec_cap, self.ram_prec_cap, self.in_field, self.polynomial(), self._prec_type)

    def _repr_(self):
        """
        String representation of this powcomputer.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint
            sage: A = PowComputer_flint(5, 20, 20, 20, False); A
            FLINT PowComputer for 5
        """
        return "FLINT PowComputer for %s" % self.prime

    cdef fmpz_t* pow_fmpz_t_tmp(self, unsigned long n) except NULL:
        """
        Return a pointer to a FLINT ``fmpz_t`` holding `p^n`.

        Analogous to
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`
        but with FLINT ``fmpz_t`` rather than GMP ``mpz_t``.  The same
        important warnings apply.
        """
        cdef padic_ctx_struct ctx = (<padic_ctx_struct*>self.ctx)[0]
        if ctx.min <= n and n < ctx.max:
            self._fpow_array[0] = (ctx.pow + (n - ctx.min))[0]
            return &self._fpow_array
        else:
            fmpz_pow_ui(self._fpow_variable, self.fprime, n)
            return &self._fpow_variable

    cdef mpz_srcptr pow_mpz_t_tmp(self, long n) except NULL:
        """
        Return a pointer to an ``mpz_t`` holding `p^n`.

        See
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`
        for important warnings.
        """
        fmpz_get_mpz(self.temp_m, self.pow_fmpz_t_tmp(n)[0])
        return self.temp_m

    cdef mpz_srcptr pow_mpz_t_top(self) noexcept:
        """
        Return a pointer to an ``mpz_t`` holding `p^N`, where `N` is
        the precision cap.
        """
        return self.top_power

    cdef unsigned long capdiv(self, unsigned long n) noexcept:
        """
        Return ceil(n / e).
        """
        if self.e == 1: return n
        if n == 0: return 0
        return (n-1) / self.e + 1

    def polynomial(self, n=None, var='x'):
        """
        Return ``None``.

        For consistency with subclasses.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint
            sage: A = PowComputer_flint(5, 20, 20, 20, False, None)
            sage: A.polynomial() is None
            True
        """
        return None

cdef class PowComputer_flint_1step(PowComputer_flint):
    """
    A PowComputer for a `p`-adic extension defined by a single polynomial.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_1step
        sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
        sage: A = PowComputer_flint_1step(5, 20, 20, 20, False, f); A
        FLINT PowComputer for 5 with polynomial x^3 - 8*x - 2
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, _poly):
        """
        Memory initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_1step
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_1step(5, 20, 20, 20, False, f)
            sage: type(A)
            <class 'sage.rings.padics.pow_computer_flint.PowComputer_flint_1step'>
        """
        cdef Polynomial_integer_dense_flint poly = _poly
        cdef long length = fmpz_poly_length(poly._poly)

        cdef Py_ssize_t i

        # fmpz_init does not allocate memory
        fmpz_init(self.q)

        sig_on()
        try:
            self._moduli = <fmpz_poly_t*>sig_malloc(sizeof(fmpz_poly_t) * (cache_limit + 2))
            if self._moduli == NULL:
                raise MemoryError
            try:
                fmpz_poly_init2(self.modulus, length)
                try:
                    fmpz_poly_init2(self.powhelper_oneunit, length)
                    try:
                        fmpz_poly_init2(self.powhelper_teichdiff, length)
                        try:
                            for i in range(1,cache_limit+2):
                                try:
                                    fmpz_poly_init2(self._moduli[i], length)
                                except BaseException:
                                    i-=1
                                    while i:
                                        fmpz_poly_clear(self._moduli[i])
                                        i-=1
                                    raise
                        except BaseException:
                            fmpz_poly_clear(self.powhelper_teichdiff)
                            raise
                    except BaseException:
                        fmpz_poly_clear(self.powhelper_oneunit)
                        raise
                except BaseException:
                    fmpz_poly_clear(self.modulus)
                    raise
            except BaseException:
                sig_free(self._moduli)
                raise
        finally:
            sig_off()

        self._allocated = 8

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, _poly, shift_seed=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f, 'capped-rel')
            sage: TestSuite(A).run()
        """
        PowComputer_flint.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, _poly, shift_seed)

        cdef Polynomial_integer_dense_flint poly = _poly
        cdef long length = fmpz_poly_length(poly._poly)
        self.deg = length - 1

        fmpz_poly_set(self.modulus, poly._poly)

        cdef Py_ssize_t i
        cdef fmpz* coeffs = self.modulus.coeffs
        fmpz_one(self.tfmpz)
        for i in range(1,cache_limit+1):
            fmpz_mul(self.tfmpz, self.tfmpz, self.fprime)
            _fmpz_vec_scalar_mod_fmpz((<fmpz_poly_struct*>self._moduli[i])[0].coeffs, coeffs, length, self.tfmpz)
            _fmpz_poly_set_length(self._moduli[i], length)

        _fmpz_poly_set_length(self._moduli[cache_limit+1], length)

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_1step
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_1step(5, 20, 20, 20, False, f)
            sage: del A
        """
        cdef Py_ssize_t i

        if self._allocated >= 8:
            fmpz_clear(self.q)
            fmpz_poly_clear(self.modulus)
            fmpz_poly_clear(self.powhelper_oneunit)
            fmpz_poly_clear(self.powhelper_teichdiff)
            for i in range(1, self.cache_limit + 1):
                fmpz_poly_clear(self._moduli[i])
            sig_free(self._moduli)

    def _repr_(self):
        """
        String representation of this powcomputer.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_1step
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_1step(5, 20, 20, 20, False, f); A
            FLINT PowComputer for 5 with polynomial x^3 - 8*x - 2
        """
        return "FLINT PowComputer for %s with polynomial %s" % (self.prime, self.polynomial())

    def __richcmp__(self, other, int op):
        """
        Comparison.

        Lexicographic on class, prime, precision cap, cache_limit and polynomial.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_1step
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2; g = x^3 - (8 + 5^22)*x - 2
            sage: A = PowComputer_flint_1step(5, 20, 20, 20, False, f)
            sage: B = PowComputer_flint_1step(5, 20, 20, 20, False, g)
            sage: A == B
            False
        """
        if not isinstance(other, PowComputer_flint_1step):
            if op in [Py_EQ, Py_NE]:
                return (op == Py_NE)
            return NotImplemented

        cdef PowComputer_flint_1step s = self
        cdef PowComputer_flint_1step o = other

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

        lx = s.in_field
        rx = o.in_field
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        if fmpz_poly_equal(s.modulus, o.modulus):
            return (op == Py_EQ)
        if op != Py_NE:
            return NotImplemented
        return False

    cdef fmpz_poly_t* get_modulus(self, unsigned long k) noexcept:
        """
        Return the defining polynomial reduced modulo `p^k`.

        The same warnings apply as for
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`.
        """
        cdef long c
        if k <= self.cache_limit:
            return &(self._moduli[k])
        else:
            c = self.cache_limit+1
            _fmpz_vec_scalar_mod_fmpz((<fmpz_poly_struct*>self._moduli[c])[0].coeffs,
                                      (<fmpz_poly_struct*>self.modulus)[0].coeffs,
                                      self.deg + 1,
                                      self.pow_fmpz_t_tmp(k)[0])
            return &(self._moduli[c])

    cdef fmpz_poly_t* get_modulus_capdiv(self, unsigned long k) noexcept:
        """
        Return the defining polynomial reduced modulo `p^a`, where
        `a` is the ceiling of `k/e`.

        The same warnings apply as for
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`.
        """
        return self.get_modulus(self.capdiv(k))

    def polynomial(self, _n=None, var='x'):
        """
        Return the polynomial attached to this ``PowComputer``, possibly
        reduced modulo a power of `p`.

        INPUT:

        - ``_n`` -- (default: ``None``) an integer, the power of `p`
          modulo which to reduce

        - ``var`` -- (default: ``'x'``) the variable for the returned polynomial

        .. NOTE::

            From Cython you should use :meth:`get_modulus` instead for
            speed.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_1step
            sage: R.<y> = ZZ[]; f = y^3 - 8*y - 2
            sage: A = PowComputer_flint_1step(5, 20, 20, 20, False, f)
            sage: A.polynomial()
            x^3 - 8*x - 2
            sage: A.polynomial(var='y')
            y^3 - 8*y - 2
            sage: A.polynomial(2)
            x^3 + 17*x + 23
        """
        R = ZZ[var]
        x = R.gen()
        cdef Polynomial_integer_dense_flint ans = (<Polynomial_integer_dense_flint?>x)._new()
        if _n is None:
            fmpz_poly_set(ans._poly, self.modulus)
        else:
            fmpz_poly_set(ans._poly, self.get_modulus(_n)[0])
        return ans

    cdef _new_fmpz_poly(self, fmpz_poly_t value, var='x'):
        """
        Return a polynomial with the value stored in ``value`` and
        variable name ``var``.
        """
        R = ZZ[var]
        x = R.gen()
        cdef Polynomial_integer_dense_flint ans = (<Polynomial_integer_dense_flint?>x)._new()
        fmpz_poly_set(ans._poly, value)
        return ans

cdef class PowComputer_flint_unram(PowComputer_flint_1step):
    """
    A PowComputer for a `p`-adic extension defined by an unramified polynomial.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_unram
        sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
        sage: A = PowComputer_flint_unram(5, 20, 20, 20, False, f); A
        FLINT PowComputer for 5 with polynomial x^3 - 8*x - 2
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, _poly, shift_seed=None):
        """
        Memory initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_unram
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_unram(5, 20, 20, 20, False, f)
            sage: type(A)
            <class 'sage.rings.padics.pow_computer_flint.PowComputer_flint_unram'>
        """
        # fmpz_init does not allocate memory
        fmpz_init(self.fmpz_ccmp)
        fmpz_init(self.fmpz_cval)
        fmpz_init(self.fmpz_cinv)
        fmpz_init(self.fmpz_cinv2)
        fmpz_init(self.fmpz_cexp)
        fmpz_init(self.fmpz_ctm)
        fmpz_init(self.fmpz_cconv)

        # While the following allocations have the potential to leak
        # small amounts of memory when interrupted or when one of the
        # init methods raises a MemoryError, the only leak-free
        # solution we could devise used 11-nested try blocks.  We
        # choose readable code in this case.
        sig_on()
        fmpz_poly_init(self.poly_cconv)
        fmpz_poly_init(self.poly_ctm)
        fmpz_poly_init(self.poly_ccmp)
        fmpz_poly_init(self.poly_cinv)
        fmpz_poly_init(self.poly_cisunit)
        fmpz_poly_init(self.poly_cinv2)
        fmpz_poly_init(self.poly_flint_rep)
        fmpz_poly_init(self.poly_matmod)
        fmpz_poly_init(self.shift_rem)
        fmpz_poly_init(self.aliasing)
        mpz_init(self.mpz_cpow)
        mpz_init(self.mpz_ctm)
        mpz_init(self.mpz_cconv)
        mpz_init(self.mpz_matmod)
        sig_off()

        self._allocated = 16

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_unram
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_unram(5, 20, 20, 20, False, f)
            sage: del A
        """
        if self._allocated >= 16:
            fmpz_clear(self.fmpz_ccmp)
            fmpz_clear(self.fmpz_cval)
            fmpz_clear(self.fmpz_cinv)
            fmpz_clear(self.fmpz_cinv2)
            fmpz_clear(self.fmpz_cexp)
            fmpz_clear(self.fmpz_ctm)
            fmpz_clear(self.fmpz_cconv)
            mpz_clear(self.mpz_cconv)
            mpz_clear(self.mpz_ctm)
            mpz_clear(self.mpz_cpow)
            mpz_clear(self.mpz_matmod)
            fmpz_poly_clear(self.poly_cconv)
            fmpz_poly_clear(self.poly_ctm)
            fmpz_poly_clear(self.poly_ccmp)
            fmpz_poly_clear(self.poly_cinv)
            fmpz_poly_clear(self.poly_cisunit)
            fmpz_poly_clear(self.poly_cinv2)
            fmpz_poly_clear(self.poly_flint_rep)
            fmpz_poly_clear(self.poly_matmod)
            fmpz_poly_clear(self.shift_rem)
            fmpz_poly_clear(self.aliasing)

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f, 'capped-rel')
            sage: TestSuite(A).run()
        """
        PowComputer_flint_1step.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

        self.e = 1
        self.f = fmpz_poly_degree(self.modulus)
        fmpz_pow_ui(self.q, self.fprime, self.f)

cdef class PowComputer_flint_eis(PowComputer_flint_1step):
    """
    A PowComputer for a `p`-adic extension defined by an Eisenstein polynomial.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_eis
        sage: R.<x> = ZZ[]; f = x^3 - 25*x + 5
        sage: A = PowComputer_flint_eis(5, 20, 20, 60, False, f); A
        FLINT PowComputer for 5 with polynomial x^3 - 25*x + 5
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_eis
            sage: R.<x> = ZZ[]; f = x^3 - 25*x + 5
            sage: A = PowComputer_flint_eis(5, 20, 20, 60, False, f)
            sage: type(A)
            <class 'sage.rings.padics.pow_computer_flint.PowComputer_flint_eis'>
        """
        PowComputer_flint_1step.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

        self.e = fmpz_poly_degree(self.modulus)
        self.f = 1
        fmpz_set(self.q, self.fprime)


def PowComputer_flint_maker(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, prec_type):
    """
    Return an appropriate FLINT PowComputer for the given input.

    INPUT:

    - ``prime`` -- an integral prime

    - ``cache_limit`` -- nonnegative integer, controlling the
      caching.  Powers of ``prime``, reductions of ``poly`` modulo
      different powers of ``prime`` and inverses of the leading
      coefficient modulo different powers of ``prime`` are cached.
      Additional data is cached for ramified extensions.

    - ``prec_cap`` -- the power of `p` modulo which elements of
      largest precision are defined

    - ``ram_prec_cap`` -- approximately ``e*prec_cap``, where ``e`` is
      the ramification degree of the extension.  For a ramified
      extension this is what Sage calls the precision cap of the ring.
      In fact, it's possible to have rings with precision cap not a
      multiple of `e`, in which case the actual relationship between
      ``ram_prec_cap`` and ``prec_cap`` is that
      ``prec_cap = ceil(n/e)``

    - ``in_field`` -- (boolean) whether the associated ring is
      actually a field

    - ``poly`` -- the polynomial defining the extension

    - ``prec_type`` -- one of ``'capped-rel'``, ``'capped-abs'`` or
      ``'fixed-mod'``; the precision type of the ring

    .. NOTE::

        Because of the way templates work, this function imports the
        class of its return value from the appropriate element files.
        This means that the returned PowComputer will have the
        appropriate compile-time-type for Cython.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
        sage: R.<x> = ZZ[]
        sage: A = PowComputer_flint_maker(3, 20, 20, 20, False, x^3 + 2*x + 1, 'capped-rel'); type(A)
        <class 'sage.rings.padics.qadic_flint_CR.PowComputer_'>
        sage: A = PowComputer_flint_maker(3, 20, 20, 20, False, x^3 + 2*x + 1, 'capped-abs'); type(A)
        <class 'sage.rings.padics.qadic_flint_CA.PowComputer_'>
        sage: A = PowComputer_flint_maker(3, 20, 20, 20, False, x^3 + 2*x + 1, 'fixed-mod'); type(A)
        <class 'sage.rings.padics.qadic_flint_FM.PowComputer_'>
    """
    if prec_type == 'capped-rel':
        from sage.rings.padics.qadic_flint_CR import PowComputer_
    elif prec_type == 'capped-abs':
        from sage.rings.padics.qadic_flint_CA import PowComputer_
    elif prec_type == 'fixed-mod':
        from sage.rings.padics.qadic_flint_FM import PowComputer_
    elif prec_type == 'floating-point':
        from sage.rings.padics.qadic_flint_FP import PowComputer_
    else:
        raise ValueError("unknown prec_type `%s`" % prec_type)
    return PowComputer_(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)
