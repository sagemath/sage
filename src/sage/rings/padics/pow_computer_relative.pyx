# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
r"""
A ``PowComputer`` for relative extensions

This module provides helper classes for the various kinds of relative `p`-adic
extensions. You should never have to access these directly, unless you are
working on linkages or other low-level `p`-adics code within the Sage library.

AUTHORS:

- David Roe, Julian Rüth (2017-06-11): initial version
"""
# ****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer cimport Integer
from sage.misc.cachefunc import cached_method


cdef class PowComputer_relative(PowComputer_class):
    r"""
    Base class for a ``PowComputer`` for use in `p`-adics implemented by Sage
    Polynomials.

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
        sage: R.<a> = ZqFM(25)
        sage: S.<x> = R[]
        sage: f = x^3 - 5*x - 5*a
        sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
        sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
        sage: PC = PowComputer_relative_maker(3, 20, 20, 60, False, f, shift_seed, 'fixed-mod')

    TESTS::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative
        sage: isinstance(PC, PowComputer_relative)                                      # needs sage.libs.flint
        True
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly, shift_seed):
        r"""
        TESTS::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]
            sage: f = x^3 - 5*x - 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(3, 20, 20, 60, False, f, shift_seed, 'fixed-mod')
        """
        self._allocated = 4

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly, shift_seed):
        r"""
        TESTS::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')  # indirect doctest
            sage: TestSuite(PC).run()
        """
        PowComputer_class.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)
        self.e = poly.degree() * poly.base_ring().absolute_e()
        self.f = poly.base_ring().absolute_f()

        self.modulus = poly

        self.tmp_cconv_out = poly.parent()()
        self.tmp_ccoeffs = poly.parent()()
        self.tmp_ccmp_a = poly.parent()()
        self.tmp_ccmp_b = poly.parent()()
        self.shift_rem = poly.parent()()
        self.aliasing = poly.parent()()

        self.base_ring = poly.base_ring()
        self.poly_ring = poly.parent()
        self._shift_seed = shift_seed

    def __dealloc__(self):
        r"""
        TESTS::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]
            sage: f = x^3 - 5*x - 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')
            sage: del PC
        """

    def __reduce__(self):
        r"""
        Return a picklable representation of this ``PowComputer``.

        TESTS::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]
            sage: f = x^3 - 5*x - 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')  # indirect doctest
            sage: loads(dumps(PC)) == PC
            True
        """
        return PowComputer_relative_maker, (self.prime, self.cache_limit, self.prec_cap, self.ram_prec_cap, self.in_field, self.polynomial(), self._shift_seed, self._prec_type)

    def _repr_(self):
        r"""
        Return a string representation of this ``PowComputer``.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25,print_pos=False,show_prec=False)
            sage: S.<x> = R[]
            sage: f = x^3 + 5*x + 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')  # indirect doctest
            Relative PowComputer for modulus x^3 + 5*x + a*5
        """
        return "Relative PowComputer for modulus %s" % (self.modulus,)

    cdef unsigned long capdiv(self, unsigned long n) noexcept:
        r"""
        Return `\lceil n/e \rceil`.
        """
        if self.e == 1:
            return n
        if n == 0:
            return 0
        return (n - 1)/self.e + 1

    def polynomial(self, n=None, var='x'):
        r"""
        Return the modulus of the `p`-adic extension that is handled by this
        ``PowComputer``.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')  # indirect doctest
            sage: PC.polynomial() is f
            True
        """
        return self.modulus


cdef class PowComputer_relative_eis(PowComputer_relative):
    r"""
    A ``PowComputer`` for a relative extension defined by an Eisenstein polynomial.

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_eis, PowComputer_relative_maker
        sage: R.<a> = ZqFM(25)
        sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
        sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
        sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
        sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')

    TESTS::

        sage: isinstance(PC, PowComputer_relative_eis)                                  # needs sage.libs.flint
        True
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly, shift_seed):
        r"""
        TESTS::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, shift_seed, 'fixed-mod')
            sage: TestSuite(PC).run()
        """
        PowComputer_relative.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)
        self._inv_shift_seed = self.invert(shift_seed, self.ram_prec_cap)

    cpdef Polynomial_generic_dense invert(self, Polynomial_generic_dense a, long prec):
        r"""
        Return the inverse of ``a``.

        INPUT:

        - ``a`` -- a `p`-adic element, represented as a reduced
          Polynomial in ``poly_ring``

        - ``prec`` -- a ``long``, the required precision

        OUTPUT: a polynomial ``b`` such that ``a*b`` is one modulo `π^\mathrm{prec}`

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25,3)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: W.<w> = R.ext(f)
            sage: g = 1 + 2*w; ginv = ~g
            sage: ginv
            1 + 3*w + 4*w^2 + 2*w^3 + (3*a + 3)*w^4 + ... + (3*a + 2)*w^8
            sage: RFP = R.change(field=False, show_prec=False, type='floating-point')
            sage: shift_seed = (-f[:3] // 5).change_ring(RFP)
            sage: PC = PowComputer_relative_maker(5, 3, 3, 9, False, f, shift_seed, 'fixed-mod')
            sage: g = 1 + 2*x
            sage: ginv = PC.invert(g, 5); ginv
            (4 + (3*a + 1)*5 + (2*a + 2)*5^2)*x^2 + (3 + (a + 1)*5 + (3*a + 2)*5^2)*x + 1 + 2*a*5 + 2*5^2
        """
        k = self.base_ring.residue_field()
        a0 = k(a[0])
        if a0.is_zero():
            raise ValueError("element has no inverse")
        K = self.base_ring.change(field=True)
        Qpmodulus = self.modulus.change_ring(K)
        Qpa = a.change_ring(K)
        R = Qpa.parent()
        inv = R([~a0])
        curprec = 1
        while curprec < prec:
            # Newton iteration
            inv = 2*inv - inv**2 * Qpa
            curprec *= 2
            inv = inv % Qpmodulus
        return inv.change_ring(self.base_ring)

    @cached_method
    def px_pow(self, r):
        r"""
        Return `p/π^r` where π is the uniformizer and `p` is the uniformizer of
        the base ring (not necessarily an integer.)

        INPUT:

        - ``r`` -- integer with 0 <= r < e

        OUTPUT: a reduced polynomial in π

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Zq(25, prec=3)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: W.<w> = R.ext(f)
            sage: elt = W.prime_pow.px_pow(2); elt
            ((4*a + 4) + (4*a + 1)*5 + (4*a + 2)*5^2)*x^2 + ((2*a + 3) + (2*a + 4)*5 + (2*a + 4)*5^2)*x + (a + 1)*5 + 3*5^2 + 2*5^3
        """
        if r < 0:
            raise ValueError("r must be nonnegative")
        elif r == 0:
            return self.poly_ring(self.base_ring.uniformizer())
        elif r >= self.e:
            raise NotImplementedError
        else:
            return (self._inv_shift_seed << (self.e-r)) % self.modulus

    @cached_method
    def pxe_pow(self, r):
        r"""
        Return the ``r``-th power of the unit `p/π^e` where `e` is the relative
        ramification index and `p` is the uniformizer of the base ring (not necessarily an integer.)

        INPUT:

        - ``r`` -- nonnegative integer

        OUTPUT: a reduced polynomial in π

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Zq(25, prec=3)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: W.<w> = R.ext(f)
            sage: elt = W.prime_pow.pxe_pow(2); elt
            ((4*a + 2) + (a + 4)*5 + 2*a*5^2)*x^2 + ((a + 2) + (a + 2)*5 + (2*a + 4)*5^2)*x + (a + 1) + (3*a + 2)*5 + (2*a + 2)*5^2
        """
        if r < 0:
            raise ValueError("r must be nonnegative")
        elif r == 0:
            return self.poly_ring.one()
        elif r == 1:
            return self._inv_shift_seed
        elif r%2:
            return (self.pxe_pow(r-1) * self.pxe_pow(1)) % self.modulus
        else:
            return (self.pxe_pow(r//2)*self.pxe_pow(r//2)) % self.modulus

    @cached_method
    def uniformizer_pow(self, r):
        r"""
        Return the ``r``-th power of the uniformizer.

        INPUT:

        - ``r`` -- nonnegative integer

        OUTPUT: a reduced polynomial in π

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Zq(25, prec=3)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: W.<w> = R.ext(f)
            sage: W.prime_pow.uniformizer_pow(2)
            x^2
        """
        if r < 0:
            raise ValueError("r must be nonnegative")
        elif r == 0:
            return self.poly_ring.one()
        elif r < self.e:
            return self.poly_ring.one() << r
        elif r % 2:
            return (self.uniformizer_pow(r-1) << 1) % self.modulus
        else:
            return (self.uniformizer_pow(r//2) * self.uniformizer_pow(r//2)) % self.modulus


def PowComputer_relative_maker(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed, prec_type):
    r"""
    Create a ``PowComputer``.

    INPUT:

    - ``prime`` -- a uniformizer in the base ring

    - ``cache_limit`` -- nonnegative integer, controlling the caching. The
      ``PowComputer`` caches frequently used things like powers of ``prime``
      This parameter, e.g., controls up to which power these are cached.

    - ``prec_cap`` -- the power of ``prime`` modulo which elements of largest
      precision are defined

    - ``ram_prec_cap`` -- approximately ``e*prec_cap``, where ``e`` is
      the relative ramification degree of the extension.  For a ramified
      extension this is what Sage calls the precision cap of the ring.  In
      fact, it's possible to have rings with precision cap not a multiple of
      `e`, in which case the actual relationship between ``ram_prec_cap`` and
      ``prec_cap`` is that ``prec_cap = ceil(n/e)``

    - ``in_field`` -- boolean; whether the associated ring is actually a
      field

    - ``poly`` -- the polynomial defining the extension

    - ``prec_type`` -- one of ``'capped-rel'``, ``'capped-abs'`` or
      ``'fixed-mod'``, ``'floating-point'``; the precision type of the ring

    .. NOTE::

        Because of the way templates work, this function imports the class of
        its return value from the appropriate element files.  This means that
        the returned PowComputer will have the appropriate compile-time-type
        for Cython.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
        sage: R.<a> = ZqFM(25, prec=2)
        sage: S.<x> = R[]
        sage: f = x^3 - 5*x - 5*a
        sage: W.<w> = R.extension(f)
        sage: PC = W.prime_pow  # indirect doctest
        sage: PC
        Relative PowComputer for modulus x^3 + (4*5 + 4*5^2)*x + 4*a*5 + 4*a*5^2
    """
    PC = PowComputer_relative_eis(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)
    # We have to set this here because the signature of __cinit__ in PowComputer_base doesn't allow for prec_type to be passed.
    PC._prec_type = prec_type
    return PC
