# distutils: libraries = NTL_LIBRARIES gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
r"""
Dense univariate polynomials over  `\ZZ/n\ZZ`, implemented using NTL

This implementation is generally slower than the FLINT implementation in
:mod:`~sage.rings.polynomial.polynomial_zmod_flint`, so we use FLINT by
default when the modulus is small enough; but NTL does not require that `n` be
``int``-sized, so we use it as default when `n` is too large for FLINT.

Note that the classes :class:`Polynomial_dense_modn_ntl_zz` and
:class:`Polynomial_dense_modn_ntl_ZZ` are different; the former is limited to
moduli less than a certain bound, while the latter supports arbitrarily large
moduli.

AUTHORS:

- Robert Bradshaw: Split off from polynomial_element_generic.py (2007-09)
- Robert Bradshaw: Major rewrite to use NTL directly (2007-09)
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.rings.polynomial.polynomial_element cimport Polynomial, _dict_to_list

from sage.libs.pari.all import pari, pari_gen

from sage.rings.integer cimport smallInteger

from sage.libs.ntl.all import ZZX, ZZ_pX
from sage.rings.integer_ring import ZZ

from sage.rings.fraction_field_element import FractionFieldElement

from sage.rings.infinity import infinity

from sage.structure.element import canonical_coercion, coerce_binop, have_same_parent

from sage.libs.ntl.types cimport NTL_SP_BOUND
from sage.libs.ntl.ZZ_p cimport *
from sage.libs.ntl.lzz_p cimport *
from sage.libs.ntl.lzz_pX cimport *
from sage.libs.ntl.ZZ_pX cimport *


def make_element(parent, args):
    return parent(*args)


zz_p_max = NTL_SP_BOUND

cdef class Polynomial_dense_mod_n(Polynomial):
    """
    A dense polynomial over the integers modulo n, where n is composite, with
    the underlying arithmetic done using NTL.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(16), implementation='NTL')
        sage: f = x^3 - x + 17
        sage: f^2
        x^6 + 14*x^4 + 2*x^3 + x^2 + 14*x + 1

        sage: loads(f.dumps()) == f
        True

        sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
        sage: p = 3*x
        sage: q = 7*x
        sage: p + q
        10*x
        sage: R.<x> = PolynomialRing(Integers(8), implementation='NTL')
        sage: parent(p)
        Univariate Polynomial Ring in x over Ring of integers modulo 100 (using NTL)
        sage: p + q
        10*x
        sage: R({10:-1})
        7*x^10

    TESTS::

        sage: f = Integers(5*2^100)['x'].random_element()
        sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n
        sage: isinstance(f, Polynomial_dense_mod_n)
        True
    """
    def __init__(self, parent, x=None, check=True,
                 is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        cdef Polynomial_dense_mod_n numer, denom

        if construct:
            if isinstance(x, ZZ_pX):
                self._poly = x
                return
            self._poly = ZZ_pX(x, parent.modulus())
            return

        self._poly = ZZ_pX([], parent.modulus())

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self._poly = (<Polynomial_dense_modn_ntl_zz>x)._poly.__copy__()
                return
            else:
                R = parent.base_ring()
                x = [ZZ(R(a)) for a in x.list()]
                check = False

        elif isinstance(x, dict):
            R = parent.base_ring()
            x = _dict_to_list(x, R(0))

        elif isinstance(x, ZZX):
            self._poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.list()]
            check = False

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_dense_mod_n):
            if x.denominator().is_unit():
                numer = x.numerator()
                denom = x.denominator().inverse_of_unit()
                x = numer._poly * denom._poly
                check = False
            else:
                raise TypeError("Denominator not a unit.")

        elif not isinstance(x, list) and not isinstance(x, tuple):
            x = [x]   # constant polynomials

        if check:
            R = parent.base_ring()
            x = [ZZ(R(a)) for a in x]

        self._poly = ZZ_pX(x, parent.modulus())

    def __reduce__(self):
        return make_element, (self.parent(), (self.list(), False, self.is_gen()))

    def int_list(self):
        return eval(str(self._poly).replace(' ',','))

    def __pari__(self, variable=None):
        """
        EXAMPLES::

            sage: t = PolynomialRing(IntegerModRing(17), "t", implementation='NTL').gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            Mod(1, 17)*t^3 + Mod(3, 17)*t
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.int_list()).Polrev(variable) * \
               pari(1).Mod(self.parent().base_ring().order())

    def ntl_ZZ_pX(self):
        r"""
        Return underlying NTL representation of this polynomial.  Additional
        ''bonus'' functionality is available through this function.

        .. warning::

            You must call ``ntl.set_modulus(ntl.ZZ(n))`` before doing
            arithmetic with this object!
        """
        return self._poly

    cdef get_unsafe(self, Py_ssize_t n):
        """
        Return the `n`-th coefficient of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n
            sage: f = Polynomial_dense_mod_n(R, [5,10,13,1,4]); f
            4*x^4 + x^3 + 13*x^2 + 10*x + 5
            sage: f[2]
            13
            sage: f[:3]
            13*x^2 + 10*x + 5
        """
        return self._parent._base((<ntl_ZZ_pX> self._poly)[n]._integer_())

    def _unsafe_mutate(self, n, value):
        n = int(n)
        if n < 0:
            raise IndexError("n must be >= 0")
        self._poly[n] = int(value)

    def _pow(self, n):
        n = int(n)
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self._poly**n, construct=True)

    cpdef _add_(self, right):
        return self.parent()(self._poly + (<Polynomial_dense_mod_n>right)._poly, construct=True)

    cpdef _mul_(self, right):
        """
        EXAMPLES::

            sage: x = PolynomialRing(Integers(100), 'x', implementation='NTL').0
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 + 90*x^2 + 32*x + 68
        """
        return self.parent()(self._poly * (<Polynomial_dense_mod_n>right)._poly, construct=True)

    cpdef _lmul_(self, Element c):
        try:
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self._poly, construct=True)
        except RuntimeError as msg: # should this really be a TypeError
            raise TypeError(msg)

    @coerce_binop
    def quo_rem(self, right):
        """
        Return a tuple ``(quotient, remainder)`` where ``self = quotient*other +
        remainder``.
        """
        v = self._poly.quo_rem((<Polynomial_dense_mod_n>right)._poly)
        P = self.parent()
        return (P(v[0], construct=True), P(v[1], construct=True) )

    def shift(self, n):
        r"""
        Return this polynomial multiplied by the power `x^n`. If `n` is negative,
        terms below `x^n` will be discarded. Does not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(12345678901234567890), implementation='NTL')
            sage: p = x^2 + 2*x + 4
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(-5)
             0
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True

        AUTHOR:

        - David Harvey (2006-08-06)
        """
        if n == 0 or self.degree() < 0:
            return self
        return self.parent()(self._poly.left_shift(n),
                             construct=True)

    cpdef _sub_(self, right):
        return self.parent()(self._poly - (<Polynomial_dense_mod_n>right)._poly, construct=True)

    def __floordiv__(self, right):
        q, _ = self.quo_rem(right)
        return q

    def degree(self, gen=None):
        """
        Return the degree of this polynomial.

        The zero polynomial has degree -1.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x^3 + 3*x - 17).degree()
            3
            sage: R.zero().degree()
            -1

        TESTS:

        Check output type (see :issue:`25182`)::

            sage: R.<x> = PolynomialRing(Integers(3), implementation='NTL')
            sage: isinstance(x.degree(), Integer)
            True
        """
        return smallInteger(max(self._poly.degree(), -1))

    cpdef list list(self, bint copy=True):
        """
        Return a new copy of the list of the underlying
        elements of ``self``.

        EXAMPLES::

            sage: _.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [83, 3, 0, 1]
        """
        return [self.get_unsafe(n) for n from 0 <= n <= self.degree()]

    def ntl_set_directly(self, v):
        r"""
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers modulo n are stored internally using
        NTL's ``ZZ_pX`` class.  Use this function to set the value of this
        polynomial using the NTL constructor, which is potentially *very* fast.
        The input v is either a vector of ints or a string of the form ``[ n1
        n2 n3 ... ]`` where the ni are integers and there are no commas between
        them. The optimal input format is the string format, since that's what
        NTL uses by default.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n as poly_modn_dense
            sage: poly_modn_dense(R, ([1,-2,3]))
            3*x^2 + 98*x + 1
            sage: f = poly_modn_dense(R, 0)
            sage: f.ntl_set_directly([1,-2,3])
            sage: f
            3*x^2 + 98*x + 1
            sage: f.ntl_set_directly('[1 -2 3 4]')
            sage: f
            4*x^3 + 3*x^2 + 98*x + 1
        """
        if self.is_gen():
            raise TypeError("Cannot change the value of the generator.")
        self._poly = ZZ_pX(v, self.parent().modulus())

    # Polynomial_singular_repr stuff, copied due to lack of multiple inheritance
    def _singular_(self, singular=None, force=False):
        self.parent()._singular_(singular, force=force).set_ring()  # this is expensive
        if self.__singular is not None:
            try:
                self.__singular._check_valid()
                if self.__singular.parent() is singular:
                    return self.__singular
            except (AttributeError, ValueError):
                pass
        return self._singular_init_(singular)

    def _singular_init_(self, singular=None, force=False):
        self.parent()._singular_(singular, force=force).set_ring()  # this is expensive
        self.__singular = singular(str(self))
        return self.__singular

    @coerce_binop
    def minpoly_mod(self, other):
        r"""
        Compute the minimal polynomial of this polynomial modulo another
        polynomial in the same ring.

        ALGORITHM:

        NTL's ``MinPolyMod()``, which uses Shoup's algorithm [Sho1999]_.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(101), implementation='NTL')
            sage: f = x^17 + x^2 - 1
            sage: (x^2).minpoly_mod(f)
            x^17 + 100*x^2 + 2*x + 100

        TESTS:

        Random testing::

            sage: p = random_prime(2^99)
            sage: R.<x> = PolynomialRing(GF(p), implementation='NTL')
            sage: d = randrange(1,50)
            sage: f = R.random_element(d)
            sage: g = R.random_element((-1,5*d))
            sage: poly = g.minpoly_mod(f)
            sage: poly(R.quotient(f)(g))
            0
        """
        mod = other.ntl_ZZ_pX()
        elt = self.ntl_ZZ_pX() % mod
        res = elt.minpoly_mod(mod)
        return self.parent()(res, construct=True)

    def small_roots(self, *args, **kwds):
        r"""
        See :func:`sage.rings.polynomial.polynomial_modn_dense_ntl.small_roots`
        for the documentation of this function.

        EXAMPLES::

            sage: N = 10001
            sage: K = Zmod(10001)
            sage: P.<x> = PolynomialRing(K, implementation='NTL')
            sage: f = x^3 + 10*x^2 + 5000*x - 222
            sage: f.small_roots()
            [4]
        """
        return small_roots(self, *args, **kwds)

    def compose_mod(self, other, modulus):
        r"""
        Compute `f(g) \pmod h`.

        To be precise about the order fo compostion, given ``self``, ``other``
        and ``modulus`` as `f(x)`, `g(x)` and `h(x)` compute `f(g(x)) \bmod h(x)`.

        INPUT:

        - ``other`` -- a polynomial `g(x)`
        - ``modulus`` -- a polynomial `h(x)`

        EXAMPLES::

            sage: R.<x> = GF(2**127 - 1)[]
            sage: f = R.random_element()
            sage: g = R.random_element()
            sage: g.compose_mod(g, f) == g(g) % f
            True

            sage: R.<x> = GF(163)[]
            sage: f = R([i for i in range(100)])
            sage: g = R([i**2 for i in range(100)])
            sage: h = 1 + x + x**5
            sage: f.compose_mod(g, h)
            82*x^4 + 56*x^3 + 45*x^2 + 60*x + 127
            sage: f.compose_mod(g, h) == f(g) % h
            True

        AUTHORS:

        - Giacomo Pope (2024-08) initial implementation
        """
        elt = self.ntl_ZZ_pX()
        mod = modulus.ntl_ZZ_pX()
        other = other.ntl_ZZ_pX()
        res = elt.compose_mod(other, mod)
        return self.parent()(res, construct=True)

    # compose_mod is the natural name from the NTL bindings, but polynomial_gf2x
    # has modular_composition as the method name so here we allow both
    modular_composition = compose_mod


def small_roots(self, X=None, beta=1.0, epsilon=None, algorithm="pari", **kwds):
    r"""
    Let `N` be the characteristic of the base ring this polynomial
    is defined over: ``N = self.base_ring().characteristic()``.
    This method returns small roots of this polynomial modulo some
    factor `b` of `N` with the constraint that `b >= N^\beta`.
    Small in this context means that if `x` is a root of `f`
    modulo `b` then `|x| < X`. This `X` is either provided by the
    user or the maximum `X` is chosen such that this algorithm
    terminates in polynomial time. If `X` is chosen automatically
    it is `X = ceil(1/2 N^{\beta^2/\delta - \epsilon})`.
    The algorithm may also return some roots which are larger than `X`.
    'This algorithm' in this context means Coppersmith's algorithm for finding
    small roots using the LLL algorithm. The implementation of this algorithm
    follows Alexander May's PhD thesis referenced below.

    INPUT:

    - ``X`` -- an absolute bound for the root (default: see above)
    - ``beta`` -- compute a root mod `b` where `b` is a factor of `N` and `b
      \ge N^\beta` (default: 1.0, so `b = N`.)
    - ``epsilon`` -- the parameter `\epsilon` described above. (default: `\beta/8`)
    - ``algorithm`` -- ``"sage"`` or ``"pari"`` (default)
    - ``**kwds`` -- passed through to method :meth:`Matrix_integer_dense.LLL() <sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL>`

    EXAMPLES:

    First consider a small example::

        sage: N = 10001
        sage: K = Zmod(10001)
        sage: P.<x> = PolynomialRing(K, implementation='NTL')
        sage: f = x^3 + 10*x^2 + 5000*x - 222

    This polynomial has no roots without modular reduction (i.e. over `\ZZ`)::

        sage: f.change_ring(ZZ).roots()
        []

    To compute its roots we need to factor the modulus `N` and use the Chinese
    remainder theorem::

        sage: p, q = N.prime_divisors()
        sage: f.change_ring(GF(p)).roots()
        [(4, 1)]
        sage: f.change_ring(GF(q)).roots()
        [(4, 1)]

        sage: crt(4, 4, p, q)
        4

    This root is quite small compared to `N`, so we can attempt to
    recover it without factoring `N` using Coppersmith's small root
    method::

        sage: f.small_roots()
        [4]

    An application of this method is to consider RSA. We are using 512-bit RSA
    with public exponent `e=3` to encrypt a 56-bit DES key. Because it would be
    easy to attack this setting if no padding was used we pad the key `K` with
    1s to get a large number::

        sage: Nbits, Kbits = 512, 56
        sage: e = 3

    We choose two primes of size 256-bit each::

        sage: p = 2^256 + 2^8 + 2^5 + 2^3 + 1
        sage: q = 2^256 + 2^8 + 2^5 + 2^3 + 2^2 + 1
        sage: N = p*q
        sage: ZmodN = Zmod( N )

    We choose a random key::

        sage: K = ZZ.random_element(0, 2^Kbits)

    and pad it with `512-56=456` 1s::

        sage: Kdigits = K.digits(2)
        sage: M = [0]*Kbits + [1]*(Nbits-Kbits)
        sage: for i in range(len(Kdigits)): M[i] = Kdigits[i]

        sage: M = ZZ(M, 2)

    Now we encrypt the resulting message::

        sage: C = ZmodN(M)^e

    To recover `K` we consider the following polynomial modulo `N`::

        sage: P.<x> = PolynomialRing(ZmodN, implementation='NTL')
        sage: f = (2^Nbits - 2^Kbits + x)^e - C

    and recover its small roots::

        sage: Kbar = f.small_roots()[0]
        sage: K == Kbar
        True

    The same algorithm can be used to factor `N = pq` if partial
    knowledge about `q` is available. This example is from the
    Magma handbook:

    First, we set up `p`, `q` and `N`::

        sage: set_random_seed(1337)
        sage: length = 512
        sage: hidden = 110
        sage: p = next_prime(2^int(round(length/2)))
        sage: q = next_prime(round(pi.n()*p))                                           # needs sage.symbolic
        sage: N = p*q                                                                   # needs sage.symbolic

    Now we disturb the low 110 bits of `q` and try to recover `q` from it::

        sage: qbar = q + ZZ.random_element(0, 2^hidden - 1)                             # needs sage.symbolic
        sage: F.<x> = PolynomialRing(Zmod(N), implementation='NTL')                     # needs sage.symbolic
        sage: f = x - qbar                                                              # needs sage.symbolic

    We know that the error is `\le 2^{\text{hidden}}-1` and that the modulus
    we are looking for is `\ge \sqrt{N}`::

        sage: from sage.misc.verbose import set_verbose
        sage: set_verbose(2)
        sage: d = f.small_roots(X=2^hidden-1, beta=0.5)[0]  # time random               # needs sage.symbolic
        verbose 2 (<module>) m = 4
        verbose 2 (<module>) t = 4
        verbose 2 (<module>) X = 1298074214633706907132624082305023
        verbose 1 (<module>) LLL of 8x8 matrix (algorithm fpLLL:wrapper)
        verbose 1 (<module>) LLL finished (time = 0.006998)
        sage: q == qbar - d                                                             # needs sage.symbolic
        True

    In general, using ``algorithm="pari"`` will return better results and be quicker. For instance,
    it is able to recover `q` even if ``hidden`` is set to `120`::

        sage: # needs sage.symbolic
        sage: hidden = 120
        sage: qbar = q + ZZ.random_element(0, 2^hidden - 1)
        sage: f = x - qbar
        sage: set_verbose(0)
        sage: f.small_roots(X=2^hidden-1, beta=0.5, algorithm="sage")  # time random
        []
        sage: f.small_roots(X=2^hidden-1, beta=0.5, algorithm="pari")  # time random
        [1203913112977791332288506012179577388]
        sage: qbar - q == _[0]
        True

    .. TODO::

        Implement improved lattice constructions for small_roots.

    REFERENCES:

    Don Coppersmith. *Finding a small root of a univariate modular equation.*
    In Advances in Cryptology, EuroCrypt 1996, volume 1070 of Lecture
    Notes in Computer Science, p. 155--165. Springer, 1996.
    http://cr.yp.to/bib/2001/coppersmith.pdf

    Alexander May. *New RSA Vulnerabilities Using Lattice Reduction Methods.*
    PhD thesis, University of Paderborn, 2003.
    http://www.cs.uni-paderborn.de/uploads/tx_sibibtex/bp.pdf
    """
    from sage.misc.verbose import verbose
    from sage.matrix.constructor import Matrix
    from sage.rings.real_mpfr import RR

    N = self.parent().characteristic()

    if not self.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    beta = RR(beta)
    if beta <= 0.0 or beta > 1.0:
        raise ValueError("0.0 < beta <= 1.0 not satisfied.")

    f = self.change_ring(ZZ)

    P,(x,) = f.parent().objgens()

    delta = f.degree()

    if epsilon is None:
        epsilon = beta/8
    verbose("epsilon = %f"%epsilon, level=2)

    m = max(beta**2/(delta * epsilon), 7*beta/delta).ceil()
    verbose("m = %d"%m, level=2)

    t = int( ( delta*m*(1/beta -1) ).floor() )
    verbose("t = %d"%t, level=2)

    if X is None:
        X = (0.5 * N**(beta**2/delta - epsilon)).ceil()
    verbose("X = %s"%X, level=2)

    Nbeta = N**beta
    ZmodN = self.base_ring()

    if algorithm == "pari":
        roots = set(map(ZmodN, pari.zncoppersmith(f, N, X=X, B=Nbeta.floor())))

    elif algorithm == "sage":
        # we could do this much faster, but this is a cheap step
        # compared to LLL
        g  = [x**j * N**(m-i) * f**i for i in range(m) for j in range(delta) ]
        g.extend([x**i * f**m for i in range(t)]) # h

        B = Matrix(ZZ, len(g), delta*m + max(delta,t) )
        for i in range(B.nrows()):
            for j in range( g[i].degree()+1 ):
                B[i,j] = g[i][j]*X**j

        B =  B.LLL(**kwds)

        f = sum([ZZ(B[0,i]//X**i)*x**i for i in range(B.ncols())])
        R = f.roots()
        roots = set([ZmodN(r) for r,m in R if abs(r) <= X])

    else:
        raise ValueError('algorithm must be "pari" or "sage"')

    return [root for root in roots if N.gcd(ZZ(self(root))) >= Nbeta]


cdef class Polynomial_dense_modn_ntl_zz(Polynomial_dense_mod_n):
    r"""
    Polynomial on `\ZZ/n\ZZ` implemented via NTL.

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _lmul_
    .. automethod:: _rmul_
    .. automethod:: _mul_
    .. automethod:: _mul_trunc_
    """
    def __init__(self, parent, v=None, check=True, is_gen=False, construct=False):
        r"""
        EXAMPLES::

            sage: R = Integers(5**21)
            sage: S.<x> = PolynomialRing(R, implementation='NTL')
            sage: S(1/4)
            357627868652344
        """
        if isinstance(v, Polynomial):
            if (<Element>v)._parent == parent:
                Polynomial.__init__(self, parent, is_gen=is_gen)
                self.x = (<Polynomial_dense_modn_ntl_zz>v).x
                self.c = (<Polynomial_dense_modn_ntl_zz>v).c
                return

        Polynomial_dense_mod_n.__init__(self, parent, v, check=check, is_gen=is_gen, construct=construct)
        v = [a for a in self._poly.list()]
        self._poly = None # this will eventually go away
        cdef ntl_zz_pX ntl = ntl_zz_pX(v, parent.modulus()) # let it handle the hard work
        self.x = ntl.x
        self.c = ntl.c

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()

    def ntl_set_directly(self, v):
        # TODO: Get rid of this
        Polynomial_dense_mod_n.ntl_set_directly(self, v)
        # verbatim from __init__
        v = [int(a) for a in self._poly.list()]
        self._poly = None # this will eventually go away
        cdef ntl_zz_pX ntl = ntl_zz_pX(v, self._parent.modulus()) # let it handle the hard work
        self.x = ntl.x
        self.c = ntl.c

    cdef Polynomial_dense_modn_ntl_zz _new(self):
        cdef Polynomial_dense_modn_ntl_zz y = <Polynomial_dense_modn_ntl_zz>Polynomial_dense_modn_ntl_zz.__new__(Polynomial_dense_modn_ntl_zz)
        y.c = self.c
        y._parent = self._parent
        return y

    def int_list(self):
        """
        Return the coefficients of ``self`` as efficiently as possible as a
        list of python ints.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n as poly_modn_dense
            sage: f = poly_modn_dense(R,[5,0,0,1])
            sage: f.int_list()
            [5, 0, 0, 1]
            sage: [type(a) for a in f.int_list()]
            [<... 'int'>, <... 'int'>, <... 'int'>, <... 'int'>]
        """
        cdef long i
        return [ zz_p_rep(zz_pX_GetCoeff(self.x, i)) for i from 0 <= i <= zz_pX_deg(self.x) ]

    cdef get_unsafe(self, Py_ssize_t n):
        """
        Return the `n`-th coefficient of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_modn_ntl_zz
            sage: f = Polynomial_dense_modn_ntl_zz(R,[2, 1])^7
            sage: f[3]
            60
            sage: f[:6]
            84*x^5 + 80*x^4 + 60*x^3 + 72*x^2 + 48*x + 28
            sage: f[:50] == f
            True
        """
        return self._parent._base(zz_p_rep(zz_pX_GetCoeff(self.x, n)))

    def _unsafe_mutate(self, n, value):
        self.c.restore_c()
        zz_pX_SetCoeff_long(self.x, n, value)

    cpdef _add_(self, _right):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x + 5) + (x^2 - 6)
            x^2 + x + 99
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 1000000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_add(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef _sub_(self, _right):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x + 5) - (x^2 - 6)
            99*x^2 + x + 11
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 1000000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_sub(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef _mul_(self, _right):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x + 5) * (x^2 - 1)
            x^3 + 5*x^2 + 99*x + 95
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 10000
        if do_sig: sig_on()
        self.c.restore_c()
        if self is right:
            zz_pX_sqr(r.x, self.x)
        else:
            zz_pX_mul(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n):
        r"""
        Return the product of ``self`` and ``right`` truncated to the
        given length `n`

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation="NTL")
            sage: f = x - 2
            sage: g = x^2 - 8*x + 16
            sage: f*g
            x^3 + 90*x^2 + 32*x + 68
            sage: f._mul_trunc_(g, 42)
            x^3 + 90*x^2 + 32*x + 68
            sage: f._mul_trunc_(g, 3)
            90*x^2 + 32*x + 68
            sage: f._mul_trunc_(g, 2)
            32*x + 68
            sage: f._mul_trunc_(g, 1)
            68
            sage: f._mul_trunc_(g, 0)
            0
            sage: f = x^2 - 8*x + 16
            sage: f._mul_trunc_(f, 2)
            44*x + 56
        """
        cdef Polynomial_dense_modn_ntl_zz op2 = <Polynomial_dense_modn_ntl_zz> right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(op2.x) > 10000
        if do_sig: sig_on()
        self.c.restore_c()
        if self is op2:
            zz_pX_SqrTrunc(r.x, self.x, n)
        else:
            zz_pX_MulTrunc(r.x, self.x, op2.x, n)
        if do_sig: sig_off()
        return r

    cpdef _rmul_(self, Element c):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x+5) * 3
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) > 100000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_rmul(r.x, self.x, c)
        if do_sig: sig_off()
        return r

    cpdef _lmul_(self, Element c):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: 3 * (x+5)
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) > 100000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_lmul(r.x, c, self.x)
        if do_sig: sig_off()
        return r

    def __pow__(Polynomial_dense_modn_ntl_zz self, ee, modulus):
        r"""
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x-1)^5
            x^5 + 95*x^4 + 10*x^3 + 90*x^2 + 5*x + 99

        We define ``0^0`` to be unity, :issue:`13895`::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: type(R(0)^0) == type(R(0))
            True

        Negative powers work (over prime fields) but use the generic
        implementation of fraction fields::

            sage: R.<x> = PolynomialRing(Integers(101), implementation='NTL')
            sage: f = (x-1)^(-5)
            sage: type(f)
            <class 'sage.rings.fraction_field_element.FractionFieldElement_1poly_field'>
            sage: (f + 2).numerator()
            2*x^5 + 91*x^4 + 20*x^3 + 81*x^2 + 10*x + 100

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: (x-1)^(-5)
            Traceback (most recent call last):
            ...
            TypeError: ...
        """
        cdef bint recip = 0, do_sig
        cdef long e = ee
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1
            e = -e
        if self == 0 and e == 0:
            return self.parent(1)

        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef zz_pX_Modulus_c mod

        self.c.restore_c()

        if modulus is None:
            if zz_pX_IsX(self.x):
                zz_pX_LeftShift(r.x, self.x, e-1)
            else:
                do_sig = zz_pX_deg(self.x) *e > 1000
                if do_sig: sig_on()
                zz_pX_power(r.x, self.x, e)
                if do_sig: sig_off()
        else:
            if not isinstance(modulus, Polynomial_dense_modn_ntl_zz):
                modulus = self.parent().coerce(modulus)
            zz_pX_Modulus_build(mod, (<Polynomial_dense_modn_ntl_zz>modulus).x)

            do_sig = zz_pX_deg(self.x) * e * self.c.p_bits > 1e5
            if do_sig: sig_on()
            zz_pX_PowerMod_long_pre(r.x, self.x, e, mod)
            if do_sig: sig_off()

        if recip:
            return ~r
        else:
            return r

    @coerce_binop
    def quo_rem(self, right):
        r"""
        Return `q` and `r`, with the degree of `r` less than the degree of ``right``,
        such that `q \cdot` ``right`` `{}+ r =` ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(125), implementation='NTL')
            sage: f = x^5+1; g = (x+1)^2
            sage: q, r = f.quo_rem(g)
            sage: q
            x^3 + 123*x^2 + 3*x + 121
            sage: r
            5*x + 5
            sage: q*g + r
            x^5 + 1
        """
        cdef Polynomial_dense_modn_ntl_zz q = self._new()
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        sig_on()
        self.c.restore_c()
        zz_pX_divrem(q.x, r.x, self.x, denom.x)
        sig_off()
        return q, r

    cpdef _floordiv_(self, right):
        """
        Return the whole part of ``self``/``right``, without remainder.

        For q = n // d, we have deg(n - q*d) < deg(d)

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(25), implementation='NTL')
            sage: f = x^7 + 1; g = x^2 - 1
            sage: q = f // g; q
            x^5 + x^3 + x
            sage: f - q*g
            x + 1
        """
        cdef Polynomial_dense_modn_ntl_zz numer = <Polynomial_dense_modn_ntl_zz>self
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef Polynomial_dense_modn_ntl_zz q = numer._new()
        sig_on()
        numer.c.restore_c()
        zz_pX_div(q.x, numer.x, denom.x)
        sig_off()
        return q

    cpdef _mod_(self, right):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(81), implementation='NTL')
            sage: f = x^7 + x + 1; g = x^3
            sage: r = f % g; r
            x + 1
            sage: g * x^4 + r
            x^7 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_zz numer = <Polynomial_dense_modn_ntl_zz>self
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef Polynomial_dense_modn_ntl_zz r = numer._new()
        sig_on()
        numer.c.restore_c()
        zz_pX_mod(r.x, numer.x, denom.x)
        sig_off()
        return r

    def shift(self, n):
        """
        Shift ``self`` to left by `n`, which is multiplication by `x^n`,
        truncating if `n` is negative.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^7 + x + 1
            sage: f.shift(1)
            x^8 + x^2 + x
            sage: f.shift(-1)
            x^6 + 1
            sage: f.shift(10).shift(-10) == f
            True

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True
        """
        return self << n

    def __lshift__(Polynomial_dense_modn_ntl_zz self, long n):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^5 + 2*x + 1
            sage: f << 3
            x^8 + 2*x^4 + x^3
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(Polynomial_dense_modn_ntl_zz self, long n):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^5 + 2*x + 1
            sage: f >> 3
            x^2
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_RightShift(r.x, self.x, n)
        return r

    def _derivative(self, var=None):
        r"""
        Return the formal derivative of ``self`` with respect to ``var``.

        ``var`` must be either the generator of the polynomial ring to which
        this polynomial belongs, or ``None`` (either way the behaviour is the
        same).

        .. SEEALSO:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^4 - x - 1
            sage: f._derivative()
            4*x^3 + 76
            sage: f._derivative(None)
            4*x^3 + 76

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

        TESTS::

            sage: y = var("y")                                                          # needs sage.symbolic
            sage: f._derivative(y)                                                      # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var != self._parent.gen():
            raise ValueError('cannot differentiate with respect to {}'.format(var))

        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_diff(r.x, self.x)
        return r

    def reverse(self, degree=None):
        """
        Return the reverse of the input polynomial thought as a polynomial of
        degree ``degree``.

        If `f` is a degree-`d` polynomial, its reverse is `x^d f(1/x)`.

        INPUT:

        - ``degree`` -- ``None`` or integer; if specified, truncate or zero
          pad the list of coefficients to this degree before reversing it

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^4 - x - 1
            sage: f.reverse()
            76*x^4 + 76*x^3 + 1
            sage: f.reverse(2)
            76*x^2 + 76*x
            sage: f.reverse(5)
            76*x^5 + 76*x^4 + x
            sage: g = x^3 - x
            sage: g.reverse()
            76*x^2 + 1

        TESTS:

        We check that this implementation is compatible with the generic one::

            sage: all(f.reverse(d) == Polynomial.reverse(f, d)
            ....:     for d in [None, 0, 1, 2, 3, 4, 5])
            True
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()

        if degree is None:
            zz_pX_reverse(r.x, self.x)
        else:
            zz_pX_reverse_hi(r.x, self.x, degree)
        return r

    def is_gen(self):
        return zz_pX_IsX(self.x)

    def __bool__(self):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^4 - x - 1
            sage: not f
            False
            sage: not (x-x)
            True
        """
        return not zz_pX_IsZero(self.x)

    def valuation(self):
        """
        Return the valuation of ``self``, that is, the power of the
        lowest nonzero monomial of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10), implementation='NTL')
            sage: x.valuation()
            1
            sage: f = x-3; f.valuation()
            0
            sage: f = x^99; f.valuation()
            99
            sage: f = x-x; f.valuation()
            +Infinity
        """
        cdef long n
        for n from 0 <= n <= zz_pX_deg(self.x):
            if zz_p_rep(zz_pX_GetCoeff(self.x, n)):
                return n
        return infinity

    def degree(self):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = x^4 - x - 1
            sage: f.degree()
            4
            sage: f = 77*x + 1
            sage: f.degree()
            0
        """
        return zz_pX_deg(self.x)

    cpdef Polynomial truncate(self, long n):
        """
        Return this polynomial mod `x^n`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(77), implementation='NTL')
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_trunc(r.x, self.x, n)
        return r

    def __call__(self, *args, **kwds):
        """
        Evaluate ``self`` at ``x``. If ``x`` is a single argument coercible into
        the base ring of ``self``, this is done directly in NTL, otherwise
        the generic Polynomial call code is used.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: f = x^3 + 7
            sage: f(5)
            32
            sage: f(5r)
            32
            sage: f(mod(5, 1000))
            32
            sage: f(x)
            x^3 + 7
            sage: S.<y> = PolynomialRing(Integers(5), implementation='NTL')
            sage: f(y)
            y^3 + 2

        TESTS::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: f = x^3 + 7
            sage: f(1).parent() == R.base_ring()
            True
            sage: f(int(1)).parent() == R.base_ring()
            True
            sage: f(x + 1).parent() == f.parent()
            True

            sage: R.<x> = PolynomialRing(Zmod(12), 'x', implementation='NTL')
            sage: u = Zmod(4)(3)
            sage: x(u).parent() == u.parent()
            True
        """
        if len(args) != 1 or len(kwds) != 0:
            return Polynomial.__call__(self, *args, **kwds)
        arg = args[0]
        cdef ntl_zz_p fx = ntl_zz_p(0, self.c), x = None
        if isinstance(arg, int):
            x = ntl_zz_p(arg, self.c)
        elif isinstance(arg, Integer):
            x = ntl_zz_p(arg, self.c)
        elif isinstance(arg, Element):
            if <void *>self._parent._base == <void *>(<Element>arg)._parent: # c++ pointer hack
                x = ntl_zz_p(arg, self.c)
            else:
                map = self._parent._base.coerce_map_from((<Element>arg)._parent)
                if map is not None:
                    x = ntl_zz_p(map(arg), self.c)
        if x is None:
            return Polynomial.__call__(self, *args, **kwds)
        else:
            zz_pX_eval(fx.x, self.x, x.x)
            return self._parent._base(int(fx))


cdef class Polynomial_dense_modn_ntl_ZZ(Polynomial_dense_mod_n):

    def __init__(self, parent, v=None, check=True, is_gen=False, construct=False):
        if isinstance(v, Polynomial):
            if (<Element>v)._parent == parent:
                Polynomial.__init__(self, parent, is_gen=is_gen)
                self.x = (<Polynomial_dense_modn_ntl_ZZ>v).x
                self.c = (<Polynomial_dense_modn_ntl_ZZ>v).c
                return

        Polynomial_dense_mod_n.__init__(self, parent, v, check=check, is_gen=is_gen, construct=construct)
        cdef ntl_ZZ_pX ntl = self._poly
        self._poly = None # this will eventually go away
        self.x = ntl.x
        self.c = ntl.c

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()

    cdef Polynomial_dense_modn_ntl_ZZ _new(self):
        cdef Polynomial_dense_modn_ntl_ZZ y = <Polynomial_dense_modn_ntl_ZZ>Polynomial_dense_modn_ntl_ZZ.__new__(Polynomial_dense_modn_ntl_ZZ)
        y.c = self.c
        y._parent = self._parent
        return y

    cpdef list list(self, bint copy=True):
        return [self.get_unsafe(n) for n from 0 <= n <= self.degree()]

    cdef get_unsafe(self, Py_ssize_t n):
        """
        Return the `n`-th coefficient of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_modn_ntl_ZZ
            sage: f = Polynomial_dense_modn_ntl_ZZ(R,[2,1])^7
            sage: f[3]
            560
            sage: f[:6]
            84*x^5 + 280*x^4 + 560*x^3 + 672*x^2 + 448*x + 128
            sage: f[:50] == f
            True
        """
        self.c.restore_c()
        # TODO, make this faster
        cdef ntl_ZZ_p ntl = ntl_ZZ_p(0, self.c)
        ntl.x = ZZ_pX_coeff(self.x, n)
        return self._parent._base(ntl._integer_())

    def _unsafe_mutate(self, n, value):
        self.c.restore_c()
        cdef Integer a
        if isinstance(value, Integer):
            a = <Integer>value
        else:
            a = ZZ(value)
        cdef ntl_ZZ_p val = ntl_ZZ_p(a, self.c)
        ZZ_pX_SetCoeff(self.x, n, val.x)

    cpdef _add_(self, _right):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: (x + 5) + (x^2 - 6)
            x^2 + x + 999999999999999999999999999999
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e7
        if do_sig: sig_on()
        self.c.restore_c()
        ZZ_pX_add(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef _sub_(self, _right):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: (x + 5) - (x^2 - 6)
            999999999999999999999999999999*x^2 + x + 11
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e7
        if do_sig: sig_on()
        self.c.restore_c()
        ZZ_pX_sub(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef _mul_(self, _right):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: (x+5) * (x^2 - 1)
            x^3 + 5*x^2 + 999999999999999999999999999999*x + 999999999999999999999999999995
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e5
        if do_sig: sig_on()
        self.c.restore_c()
        if self is right:
            ZZ_pX_sqr(r.x, self.x)
        else:
            ZZ_pX_mul(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n):
        """
        Return the product of ``self`` and ``right`` truncated to the
        given length `n`, only return terms of degree less than `n`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation="NTL")
            sage: f = x - 2
            sage: g = x^2 - 8*x + 16
            sage: f*g
            x^3 + 999999999999999999999999999990*x^2 + 32*x + 999999999999999999999999999968
            sage: f._mul_trunc_(g, 42)
            x^3 + 999999999999999999999999999990*x^2 + 32*x + 999999999999999999999999999968
            sage: f._mul_trunc_(g, 3)
            999999999999999999999999999990*x^2 + 32*x + 999999999999999999999999999968
            sage: f._mul_trunc_(g, 2)
            32*x + 999999999999999999999999999968
            sage: f._mul_trunc_(g, 1)
            999999999999999999999999999968
            sage: f._mul_trunc_(g, 0)
            0
            sage: f = x^2 - 8*x + 16
            sage: f._mul_trunc_(f, 2)
            999999999999999999999999999744*x + 256
        """
        cdef Polynomial_dense_modn_ntl_ZZ op2 = <Polynomial_dense_modn_ntl_ZZ> right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(op2.x)) * self.c.p_bits > 1e5
        if do_sig: sig_on()
        self.c.restore_c()
        if self is op2:
            ZZ_pX_SqrTrunc(r.x, self.x, n)
        else:
            ZZ_pX_MulTrunc(r.x, self.x, op2.x, n)
        if do_sig: sig_off()
        return r

    cpdef _rmul_(self, Element c):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: (x+5) * 3
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = ZZ_pX_deg(self.x) * self.c.p_bits > 1e7
        if do_sig: sig_on()
        self.c.restore_c()
        cdef ntl_ZZ_p value = ntl_ZZ_p(c, self.c)
        ZZ_pX_rmul(r.x, self.x, value.x)
        if do_sig: sig_off()
        return r

    cpdef _lmul_(self, Element c):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: 3 * (x+5)
            3*x + 15
        """
        return self._rmul_(c)

    def __pow__(Polynomial_dense_modn_ntl_ZZ self, ee, modulus):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: (x+1)^5
            x^5 + 5*x^4 + 10*x^3 + 10*x^2 + 5*x + 1

        We define ``0^0`` to be unity, :issue:`13895`::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: type(R(0)^0) == type(R(0))
            True
        """
        cdef bint recip = 0, do_sig
        cdef long e = ee
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1 # delay because powering frac field elements is slow
            e = -e
        if self == 0 and e == 0:
            return self.parent(1)
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef ZZ_pX_Modulus_c mod

        self.c.restore_c()

        if modulus is None:
            if ZZ_pX_IsX(self.x):
                ZZ_pX_LeftShift(r.x, self.x, e - 1)
            else:
                do_sig = ZZ_pX_deg(self.x) * e * self.c.p_bits > 1e5
                if do_sig: sig_on()
                ZZ_pX_power(r.x, self.x, e)
                if do_sig: sig_off()
        else:
            if not isinstance(modulus, Polynomial_dense_modn_ntl_ZZ):
                modulus = self.parent().coerce(modulus)
            ZZ_pX_Modulus_build(mod, (<Polynomial_dense_modn_ntl_ZZ>modulus).x)

            do_sig = ZZ_pX_deg(self.x) * e * self.c.p_bits > 1e5
            if do_sig: sig_on()
            ZZ_pX_PowerMod_long_pre(r.x, self.x, e, mod)
            if do_sig: sig_off()
        if recip:
            return ~r
        else:
            return r

    @coerce_binop
    def quo_rem(self, right):
        r"""
        Return `q` and `r`, with the degree of `r` less than the degree of ``right``,
        such that `q \cdot` ``right`` `+ r =` ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: f = x^5+1; g = (x+1)^2
            sage: q, r = f.quo_rem(g)
            sage: q
            x^3 + 999999999999999999999999999998*x^2 + 3*x + 999999999999999999999999999996
            sage: r
            5*x + 5
            sage: q*g + r
            x^5 + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ q = self._new()
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        sig_on()
        self.c.restore_c()
        ZZ_pX_DivRem(q.x, r.x, self.x, denom.x)
        sig_off()
        return q, r

    cpdef _floordiv_(self, right):
        """
        Return the whole part of ``self`` / ``right``, without remainder.

        For q = n // d, we have deg(n - q*d) < deg(d)

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: f = x^7 + 1; g = x^2 - 1
            sage: q = f // g; q
            x^5 + x^3 + x
            sage: f - q*g
            x + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ numer = <Polynomial_dense_modn_ntl_ZZ>self
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef Polynomial_dense_modn_ntl_ZZ q = numer._new()
        sig_on()
        numer.c.restore_c()
        ZZ_pX_div(q.x, numer.x, denom.x)
        sig_off()
        return q

    cpdef _mod_(self, right):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(9^30), implementation='NTL')
            sage: f = x^7 + x + 1; g = x^3 - 1
            sage: r = f % g; r
            2*x + 1
            sage: g * (x^4 + x) + r
            x^7 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ numer = <Polynomial_dense_modn_ntl_ZZ>self
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef Polynomial_dense_modn_ntl_ZZ r = numer._new()
        sig_on()
        numer.c.restore_c()
        ZZ_pX_rem(r.x, numer.x, denom.x)
        sig_off()
        return r

    def shift(self, n):
        """
        Shift ``self`` to left by `n`, which is multiplication by `x^n`,
        truncating if `n` is negative.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(12^30), implementation='NTL')
            sage: f = x^7 + x + 1
            sage: f.shift(1)
            x^8 + x^2 + x
            sage: f.shift(-1)
            x^6 + 1
            sage: f.shift(10).shift(-10) == f
            True

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True
        """
        return self << n

    def __lshift__(Polynomial_dense_modn_ntl_ZZ self, long n):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(14^30), implementation='NTL')
            sage: f = x^5 + 2*x + 1
            sage: f << 3
            x^8 + 2*x^4 + x^3
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(Polynomial_dense_modn_ntl_ZZ self, long n):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(15^30), implementation='NTL')
            sage: f = x^5 + 2*x + 1
            sage: f >> 3
            x^2
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_RightShift(r.x, self.x, n)
        return r

    def _derivative(self, var=None):
        r"""
        Return the formal derivative of ``self`` with respect to ``var``.

        ``var`` must be either the generator of the polynomial ring to which
        this polynomial belongs, or None (either way the behaviour is the
        same).

        .. SEEALSO:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(12^29), implementation='NTL')
            sage: f = x^4 + x + 5
            sage: f._derivative()
            4*x^3 + 1
            sage: f._derivative(None)
            4*x^3 + 1

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

        TESTS::

            sage: y = var("y")                                                          # needs sage.symbolic
            sage: f._derivative(y)                                                      # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var != self._parent.gen():
            raise ValueError("cannot differentiate with respect to {}".format(var))

        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_diff(r.x, self.x)
        return r

    def reverse(self, degree=None):
        """
        Return the reverse of the input polynomial thought as a polynomial of
        degree ``degree``.

        If `f` is a degree-`d` polynomial, its reverse is `x^d f(1/x)`.

        INPUT:

        - ``degree`` -- ``None`` or integer; if specified, truncate or zero
          pad the list of coefficients to this degree before reversing it

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(12^29), implementation='NTL')
            sage: f = x^4 + 2*x + 5
            sage: f.reverse()
            5*x^4 + 2*x^3 + 1
            sage: f = x^3 + x
            sage: f.reverse()
            x^2 + 1
            sage: f.reverse(1)
            1
            sage: f.reverse(5)
            x^4 + x^2

        TESTS:

        We check that this implementation is compatible with the generic one::

            sage: all(f.reverse(d) == Polynomial.reverse(f, d)
            ....:     for d in [None, 0, 1, 2, 3, 4, 5])
            True
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()

        if degree is None:
            ZZ_pX_reverse(r.x, self.x)
        else:
            ZZ_pX_reverse_hi(r.x, self.x, degree)
        return r

    def is_gen(self):
        return ZZ_pX_IsX(self.x)

    def valuation(self):
        """
        Return the valuation of ``self``, that is, the power of the
        lowest nonzero monomial of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^50), implementation='NTL')
            sage: x.valuation()
            1
            sage: f = x - 3; f.valuation()
            0
            sage: f = x^99; f.valuation()
            99
            sage: f = x - x; f.valuation()
            +Infinity
        """
        cdef long n
        cdef ZZ_p_c coeff
        for n from 0 <= n <= ZZ_pX_deg(self.x):
            coeff = ZZ_pX_coeff(self.x, n)
            if not ZZ_p_IsZero(coeff):
                return n
        return infinity

    def __bool__(self):
        """
        TESTS::

            sage: R.<x> = PolynomialRing(Integers(12^29), implementation='NTL')
            sage: f = x^4 + 1
            sage: not f
            False
            sage: not (x-x)
            True
        """
        return not ZZ_pX_IsZero(self.x)

    def degree(self):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(14^34), implementation='NTL')
            sage: f = x^4 - x - 1
            sage: f.degree()
            4
            sage: f = 14^43*x + 1
            sage: f.degree()
            0
        """
        return ZZ_pX_deg(self.x)

    cpdef Polynomial truncate(self, long n):
        """
        Return this polynomial mod `x^n`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(15^30), implementation='NTL')
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_trunc(r.x, self.x, n)
        return r

    def __call__(self, *args, **kwds):
        """
        Evaluate ``self`` at ``x``. If ``x`` is a single argument coercible into
        the base ring of ``self``, this is done directly in NTL, otherwise
        the generic Polynomial call code is used.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: f = x^3 + 7
            sage: f(5)
            132
            sage: f(5r)
            132
            sage: f(mod(5, 10^50))
            132
            sage: f(x)
            x^3 + 7
            sage: S.<y> = PolynomialRing(Integers(5), implementation='NTL')
            sage: f(y)
            y^3 + 2

        TESTS::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: f = x^3 + 7
            sage: f(1).parent() == R.base_ring()
            True
            sage: f(int(1)).parent() == R.base_ring()
            True
            sage: f(x + 1).parent() == f.parent()
            True

            sage: R.<x> = PolynomialRing(Zmod(10^30), 'x', implementation='NTL')
            sage: u = Zmod(10^29)(3)
            sage: x(u).parent() == u.parent()
            True
            sage: v = Zmod(10)(3)
            sage: x(v).parent() == v.parent()
            True
        """
        if len(args) != 1 or len(kwds) != 0:
            return Polynomial.__call__(self, *args, **kwds)
        arg = args[0]
        cdef ntl_ZZ_p fx = ntl_ZZ_p(0, self.c), x = None
        if isinstance(arg, int) or isinstance(arg, Integer):
            x = ntl_ZZ_p(arg, self.c)
        elif isinstance(arg, Element):
            if <void *>self._parent._base == <void *>(<Element>arg)._parent: # c++ pointer hack
                x = ntl_ZZ_p(arg, self.c)
            else:
                map = self._parent._base.coerce_map_from((<Element>arg)._parent)
                if map is not None:
                    x = ntl_ZZ_p(map(arg), self.c)
        if x is None:
            return Polynomial.__call__(self, *args, **kwds)
        else:
            ZZ_pX_eval(fx.x, self.x, x.x)
            return self._parent._base(fx._integer_())


cdef class Polynomial_dense_mod_p(Polynomial_dense_mod_n):
    """
    A dense polynomial over the integers modulo `p`, where `p` is prime.
    """

    def __pow__(self, n, modulus):
        """
        Exponentiation of ``self``.

        If ``modulus`` is not ``None``, the exponentiation is performed
        modulo the polynomial ``modulus``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(101), implementation='NTL')
            sage: (x-1)^5
            x^5 + 96*x^4 + 10*x^3 + 91*x^2 + 5*x + 100
            sage: pow(x-1, 15, x^3+x+1)
            55*x^2 + 6*x + 46

        TESTS:

        Negative powers work but use the generic
        implementation of fraction fields::

            sage: R.<x> = PolynomialRing(Integers(101), implementation='NTL')
            sage: f = (x-1)^(-5)
            sage: type(f)
            <class 'sage.rings.fraction_field_element.FractionFieldElement_1poly_field'>
            sage: (f + 2).numerator()
            2*x^5 + 91*x^4 + 20*x^3 + 81*x^2 + 10*x + 100

        We define ``0^0`` to be unity, :issue:`13895`::

            sage: R.<x> = PolynomialRing(Integers(101), implementation='NTL')
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: R.<x> = PolynomialRing(Integers(101), implementation='NTL')
            sage: type(R(0)^0) == type(R(0))
            True

        The modulus can have smaller degree than ``self``::

            sage: R.<x> = PolynomialRing(GF(101), implementation="NTL")
            sage: pow(x^4 + 1, 100, x^2 + x + 1)
            100*x + 100

        Canonical coercion should apply::

            sage: R.<x> = PolynomialRing(GF(101), implementation="FLINT")
            sage: x_ZZ = ZZ["x"].gen()
            sage: pow(x+1, 100, 2)
            0
            sage: pow(x+1, 100, x_ZZ^2 + x_ZZ + 1)
            100*x + 100
            sage: pow(x+1, int(100), x_ZZ^2 + x_ZZ + 1)
            100*x + 100
            sage: xx = polygen(GF(97))
            sage: _ = pow(x + 1, 3, xx^3 + xx + 1)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: ...
        """
        n = Integer(n)
        parent = (<Element>self)._parent

        if modulus is not None:
            # Similar to coerce_binop
            if not have_same_parent(self, modulus):
                a, m = canonical_coercion(self, modulus)
                if a is not self:
                    return pow(a, n, m)
                modulus = m
            self = self % modulus
            return parent(pow(self.ntl_ZZ_pX(), n, modulus.ntl_ZZ_pX()), construct=True)
        else:
            if n < 0:
                return ~(self**(-n))
            elif self.degree() <= 0:
                return parent(self[0]**n)
            else:
                return parent(self.ntl_ZZ_pX()**n, construct=True)

    @coerce_binop
    def gcd(self, right):
        """
        Return the greatest common divisor of this polynomial and ``other``, as
        a monic polynomial.

        INPUT:

        - ``other`` -- a polynomial defined over the same ring as ``self``

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(3), implementation="NTL")
            sage: f, g = x + 2, x^2 - 1
            sage: f.gcd(g)
            x + 2
        """
        g = self.ntl_ZZ_pX().gcd(right.ntl_ZZ_pX())
        return self.parent()(g, construct=True)

    @coerce_binop
    def xgcd(self, other):
        r"""
        Compute the extended gcd of this element and ``other``.

        INPUT:

        - ``other`` -- an element in the same polynomial ring

        OUTPUT:

        A tuple ``(r,s,t)`` of elements in the polynomial ring such
        that ``r = s*self + t*other``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(3), implementation='NTL')
            sage: x.xgcd(x)
            (x, 0, 1)
            sage: (x^2 - 1).xgcd(x - 1)
            (x + 2, 0, 1)
            sage: R.zero().xgcd(R.one())
            (1, 0, 1)
            sage: (x^3 - 1).xgcd((x - 1)^2)
            (x^2 + x + 1, 0, 1)
            sage: ((x - 1)*(x + 1)).xgcd(x*(x - 1))
            (x + 2, 1, 2)
        """
        r, s, t = self.ntl_ZZ_pX().xgcd(other.ntl_ZZ_pX())
        return self.parent()(r, construct=True), self.parent()(s, construct=True), self.parent()(t, construct=True)

    @coerce_binop
    def resultant(self, other):
        """
        Return the resultant of ``self`` and ``other``, which must lie in the same
        polynomial ring.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(19), implementation='NTL')
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            11
            sage: r.parent() is GF(19)
            True
        """
        other = self.parent().coerce(other)
        return self.base_ring()(str(self.ntl_ZZ_pX().resultant(other.ntl_ZZ_pX())))

    def discriminant(self):
        """
        EXAMPLES::

            sage: _.<x> = PolynomialRing(GF(19), implementation='NTL')
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            12
        """
        return self.base_ring()(str(self.ntl_ZZ_pX().discriminant()))
