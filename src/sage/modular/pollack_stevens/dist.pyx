# distutils: libraries = gmp
# distutils: extra_compile_args = -D_XPG6
"""
`p`-adic distributions spaces

This module implements `p`-adic distributions, a `p`-adic Banach
space dual to locally analytic functions on a disc.

EXAMPLES::

    sage: D = OverconvergentDistributions(5, 7, 15)
    sage: v = D([7,14,21,28,35]); v
    (7 + O(7^5), 2*7 + O(7^4), 3*7 + O(7^3), 4*7 + O(7^2), O(7))

REFERENCES:

- [PS2011]_
"""

# ****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import operator

from sage.arith.misc import binomial, bernoulli
from sage.categories.fields import Fields
from sage.matrix.constructor import matrix
from sage.matrix.matrix cimport Matrix
from sage.misc.verbose import verbose
from sage.modular.pollack_stevens.sigma0 import Sigma0
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.infinity import Infinity
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.element cimport Element
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool

#from sage.libs.flint.ulong_extras cimport *

cdef long overflow = 1 << (4 * sizeof(long) - 1)
cdef long underflow = -overflow
cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1


def get_dist_classes(p, prec_cap, base, symk, implementation):
    r"""
    Determine the element and action classes to be used for given inputs.

    INPUT:

    - ``p`` -- prime

    - ``prec_cap`` -- the `p`-adic precision cap

    - ``base`` -- the base ring

    - ``symk`` -- an element of Symk

    - ``implementation`` -- string; if not ``None``, override the
      automatic choice of implementation. May be 'long' or 'vector',
      otherwise raise a :exc:`NotImplementedError`.

    OUTPUT:

    - Either a Dist_vector and WeightKAction_vector, or a Dist_vector_long
      and WeightKAction_vector_long

    EXAMPLES::

        sage: D = OverconvergentDistributions(2, 3, 5); D # indirect doctest
        Space of 3-adic distributions with k=2 action and precision cap 5
    """
    if implementation is not None:
        if implementation == 'long':
            raise NotImplementedError('The optimized implementation -using longs- has been disabled and may return wrong results.')
        elif implementation == 'vector':
            return Dist_vector, WeightKAction_vector
        else:
            raise NotImplementedError('The implementation "%s" does not exist yet' % (implementation))

    return Dist_vector, WeightKAction_vector


cdef class Dist(ModuleElement):
    r"""
        The main `p`-adic distribution class, implemented as per the paper [PS2011]__.
    """
    def moment(self, n):
        r"""
        Return the `n`-th moment.

        INPUT:

        - ``n`` -- integer or slice, to be passed on to moments

        OUTPUT:

        - the `n`-th moment, or a list of moments in the case that `n`
          is a slice.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 7, 10)
            sage: v = D([7,14,21,28,35])
            sage: v.moment(3)
            4*7 + O(7^2)
            sage: v.moment(0)
            7 + O(7^5)
        """
        return self.parent().prime() ** (self.ordp) * self._unscaled_moment(n)

    def moments(self):
        r"""
        Return the vector of moments.

        OUTPUT: the vector of moments

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 5, 10, base = Qp(5))
            sage: v = D([1,7,4,2,-1])
            sage: v = 1/5^3 * v
            sage: v
            5^-3 * (1 + O(5^5), 2 + 5 + O(5^4), 4 + O(5^3), 2 + O(5^2), 4 + O(5))
            sage: v.moments()
            (5^-3 + O(5^2), 2*5^-3 + 5^-2 + O(5), 4*5^-3 + O(5^0), 2*5^-3 + O(5^-1), 4*5^-3 + O(5^-2))
        """
        return self.parent().prime() ** (self.ordp) * self._moments

    cpdef normalize(self, include_zeroth_moment=True):
        r"""
        Normalize so that the precision of the `i`-th moment is `n-i`,
        where `n` is the number of moments.

        OUTPUT: normalized entries of the distribution

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15); D
            Space of 7-adic distributions with k=5 action and precision cap 15
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.normalize()
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
        """
        raise NotImplementedError

    cdef long _relprec(self) noexcept:
        raise NotImplementedError

    cdef _unscaled_moment(self, long i):
        raise NotImplementedError

    cpdef long _ord_p(self) noexcept:
        r"""
        Return power of `p` by which the moments are shifted.

        .. NOTE::

            This is not necessarily the same as the valuation,
            since the moments could all be divisible by `p`.

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([7,14,21,28,35]); v
            (7 + O(7^5), 2*7 + O(7^4), 3*7 + O(7^3), 4*7 + O(7^2), O(7))
            sage: v._ord_p()
            0
        """
        return self.ordp

    def scale(self, left):
        r"""
        Scale the moments of the distribution by ``left``.

        INPUT:

        - ``left`` -- scalar

        OUTPUT: scales the moments by ``left``

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.scale(2)
            (2 + O(7^5), 4 + O(7^4), 6 + O(7^3), 1 + 7 + O(7^2), 3 + O(7))
        """
        # if isinstance(self, Dist_long) and isinstance(left, (Integer, pAdicCappedRelativeElement, pAdicCappedAbsoluteElement, pAdicFixedModElement)):
        #     return self._lmul_(left)
        R = left.parent()
        base = self.parent().base_ring()
        if base is R:
            return self._lmul_(left)
        elif base.has_coerce_map_from(R):
            return self._lmul_(base(left))
        else:
            from sage.categories.pushout import pushout
            new_base = pushout(base, R)
            V = self.parent().change_ring(new_base)
            scalar = new_base(left)
            return V([scalar * new_base(self.moment(i)) for i in range(self.precision_absolute())])

    def is_zero(self, p=None, M=None):
        r"""
        Return ``True`` if the `i`-th moment is zero for all `i` (case ``M`` is
        ``None``) or zero modulo `p^{M-i}` for all `i` (when ``M`` is not
        ``None``).

        Note that some moments are not known to precision ``M``, in which
        case they are only checked to be equal to zero modulo the
        precision to which they are defined.

        INPUT:

        - ``p`` -- prime

        - ``M`` -- precision

        OUTPUT: boolean

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.is_zero()
            False
            sage: v = D(5*[0])
            sage: v.is_zero()
            True

        ::

            sage: D = Symk(0)
            sage: v = D([0])
            sage: v.is_zero(5,3)
            True
        """
        n = self.precision_relative()
        aprec = self.precision_absolute()
        if M is None:
            M = n
        # elif M > aprec: # DEBUG
        #     return False
        elif M < aprec:
            n -= (aprec - M)
            M -= self.ordp
        if p is None:
            p = self.parent().prime()
        cdef bint usearg = True
        if n == 0:
            return True
        else:
            try:
                z = self._unscaled_moment(0).is_zero(M)
            except TypeError:
                z = self._unscaled_moment(0).is_zero()
                use_arg = False
            if not z:
                return False
            for a in range(1, n):
                if usearg:
                    try:
                        z = self._unscaled_moment(a).is_zero(M - a)
                    except TypeError:
                        z = self._unscaled_moment(a).is_zero()
                        use_arg = False
                else:
                    z = self._unscaled_moment(a).is_zero()
                if not z:
                    return False
            return True

    def find_scalar(self, _other, p, M=None, check=True):
        r"""
        Return an ``alpha`` with ``other = self * alpha``, or raises
        a :exc:`ValueError`.

        It will also raise a :exc:`ValueError` if this distribution is zero.

        INPUT:

        - ``other`` -- another distribution

        - ``p`` -- an integral prime (only used if the parent is not a Symk)

        - ``M`` -- (default: ``None``) an integer, the relative precision
          to which the scalar must be determined

        - ``check`` -- boolean (default: ``True``); whether to validate
          that ``other`` is actually a multiple of this element

        OUTPUT: a scalar ``alpha`` with ``other = self * alpha``

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5])
            sage: w = D([3,6,9,12,15])
            sage: v.find_scalar(w,p=7)
            3 + O(7^5)
            sage: v.find_scalar(w,p=7,M=4)
            3 + O(7^4)

            sage: u = D([1,4,9,16,25])
            sage: v.find_scalar(u,p=7)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple
        """
        cdef Dist other = _other
        i = 0
        n = self.precision_relative()
        other_pr = other.precision_relative()
        if n == 0:
            raise ValueError("self is zero")
        verbose("n = %s" % n, level=2)
        verbose("moment 0", level=2)
        a = self._unscaled_moment(i)
        verbose("a = %s" % a, level = 2)
        padic = isinstance(a.parent(), pAdicGeneric)
        if self.parent().is_symk():
            while a == 0:
                if other._unscaled_moment(i) != 0:
                    raise ValueError("not a scalar multiple")
                i += 1
                verbose("moment %s" % i, level = 2)
                try:
                    a = self._unscaled_moment(i)
                except IndexError:
                    raise ValueError("self is zero")
            alpha = other._unscaled_moment(i) / a
            if check:
                i += 1
                while i < n:
                    verbose("comparing moment %s" % i, level = 2)
                    if alpha * self._unscaled_moment(i) != other._unscaled_moment(i):
                        raise ValueError("not a scalar multiple")
                    i += 1
        else:
            p = self.parent().prime()
            v = a.valuation(p)
            while v >= n - i:
                i += 1
                verbose("p moment %s" % i, level = 2)
                try:
                    a = self._unscaled_moment(i)
                except IndexError:
                    raise ValueError("self is zero")
                v = a.valuation(p)
            relprec = n - i - v
#            verbose("p=%s, n-i=%s\nself.moment=%s, other.moment=%s" % (p, n-i, a, other._unscaled_moment(i)),level=2)
## RP: This code was crashing because other may have too few moments -- so I added this bound with other's relative precision
            if padic:
                if i < other_pr:
                    alpha = (other._unscaled_moment(i) / a).add_bigoh(n - i)
                else:
                    alpha = (0 * a).add_bigoh(other_pr - i)
            else:
                if i < other_pr:
                    alpha = (other._unscaled_moment(i) / a) % p ** (n - i)
                else:
                    alpha = 0
            verbose("alpha = %s" % alpha, level = 2)
## RP: This code was crashing because other may have too few moments -- so I added this bound with other's relative precision
            while i < other_pr - 1:
                i += 1
                verbose("comparing p moment %s" % i, level = 2)
                a = self._unscaled_moment(i)
                if check:
                    # verbose("self.moment=%s, other.moment=%s" % (a, other._unscaled_moment(i)))
                    if (padic and other._unscaled_moment(i) != alpha * a) or \
                       (not padic and other._unscaled_moment(i) % p ** (n - i) != alpha * a % p ** (n - i)):
                        raise ValueError("not a scalar multiple")
                v = a.valuation(p)
                if n - i - v > relprec:
                    verbose("Reseting alpha: relprec=%s, n-i=%s, v=%s" % (relprec, n - i, v), level = 2)
                    relprec = n - i - v
                    if padic:
                        alpha = (other._unscaled_moment(i) / a).add_bigoh(n - i)
                    else:
                        alpha = (other._unscaled_moment(i) / a) % p ** (n - i)
                    verbose("alpha=%s" % alpha, level = 2)
            if relprec < M:
                raise ValueError("result not determined to high enough precision")
        alpha = alpha * self.parent().prime() ** (other.ordp - self.ordp)
        verbose("alpha=%s" % alpha, level = 2)
        try:
            alpha = self.parent().base_ring()(alpha)
            if M is not None:
                alpha = alpha.add_bigoh(M)
        except (ValueError, AttributeError):
            pass
        return alpha

    def find_scalar_from_zeroth_moment(self, _other, p, M=None, check=True):
        r"""
        Return an ``alpha`` with ``other = self * alpha`` using only
        the zeroth moment, or raises a :exc:`ValueError`.

        It will also raise a :exc:`ValueError` if the zeroth moment of the
        distribution is zero.

        INPUT:

        - ``other`` -- another distribution

        - ``p`` -- an integral prime (only used if the parent is not a Symk)

        - ``M`` -- (default: ``None``) an integer, the relative precision
          to which the scalar must be determined

        - ``check`` -- boolean (default: ``True``); whether to validate
          that ``other`` is actually a multiple of this element

        OUTPUT: a scalar ``alpha`` with ``other = self * alpha``

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5])
            sage: w = D([3,6,9,12,15])
            sage: v.find_scalar_from_zeroth_moment(w,p=7)
            3 + O(7^5)
            sage: v.find_scalar_from_zeroth_moment(w,p=7,M=4)
            3 + O(7^4)

            sage: u = D([1,4,9,16,25])
            sage: v.find_scalar_from_zeroth_moment(u,p=7)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple
        """
        cdef Dist other = _other
        n = self.precision_relative()
        if n == 0:
            raise ValueError("zeroth moment is zero")
        verbose("n = %s" % n, level = 2)
        a = self.moment(0)
        if a.is_zero():
            raise ValueError("zeroth moment is zero")
        alpha = other.moment(0) / a
        if check:
            for i in range(1, n):
                verbose("comparing moment %s" % i, level = 2)
                if alpha * self.moment(i) != other.moment(i):
                    raise ValueError("not a scalar multiple")
        alpha = self.parent().base_ring()(alpha)
        if M is not None:
            try:
                absprec = alpha.precision_absolute()
                if absprec < M:
                    raise ValueError("result not determined to high "
                                     "enough precision")
                verbose("alpha=%s" % (alpha), level = 2)
                alpha = alpha.add_bigoh(M)
            except AttributeError:
                pass
        return alpha

    cpdef _richcmp_(_left, _right, int op):
        r"""
        Comparison.

        EXAMPLES:

        Equality of two distributions::

            sage: D = OverconvergentDistributions(0, 5, 10)
            sage: D([1, 2]) == D([1])
            True
            sage: D([1]) == D([1, 2])
            True
            sage: v = D([1+O(5^3),2+O(5^2),3+O(5)])
            sage: w = D([1+O(5^2),2+O(5)])
            sage: v == w
            True
            sage: D = Symk(0,Qp(5,5))
            sage: v = 5 * D([4*5^-1+3+O(5^2)])
            sage: w = D([4+3*5+O(5^2)])
            sage: v == w
            True
        """
        cdef Dist left = _left
        cdef Dist right = _right
        left.normalize()
        right.normalize()
        cdef long rprec = min(left._relprec(), right._relprec())
        cdef long i
        p = left.parent().prime()
        if left.ordp > right.ordp:
            shift = p ** (left.ordp - right.ordp)
            for i in range(rprec):
                lx = shift * left._unscaled_moment(i)
                rx = right._unscaled_moment(i)
                if lx != rx:
                    return richcmp_not_equal(lx, rx, op)
        elif left.ordp < right.ordp:
            shift = p ** (right.ordp - left.ordp)
            for i in range(rprec):
                lx = left._unscaled_moment(i)
                rx = shift * right._unscaled_moment(i)
                if lx != rx:
                    return richcmp_not_equal(lx, rx, op)
        else:
            for i in range(rprec):
                lx = left.moment(i)
                rx = right.moment(i)
                if lx != rx:
                    return richcmp_not_equal(lx, rx, op)
        return rich_to_bool(op, 0)

    def diagonal_valuation(self, p=None):
        """
        Return the largest `m` so that this distribution lies in `Fil^m`.

        INPUT:

        - ``p`` -- (default: ``None``) a positive integral prime

        OUTPUT:

        The largest integer `m` so that `p^m` divides the `0`-th
        moment, `p^{m-1}` divides the first moment, etc.

        EXAMPLES::

            sage: D = OverconvergentDistributions(8, 7, 15)
            sage: v = D([7^(5-i) for i in range(1,5)])
            sage: v
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v.diagonal_valuation(7)
            4
        """
        if p is None:
            p = self.parent()._p
        n = self.precision_relative()
        return self.ordp + min([n] + [a + self._unscaled_moment(a).valuation(p) for a in range(n)])

    def valuation(self, p=None):
        """
        Return the minimum valuation of any moment.

        INPUT:

        - ``p`` -- (default: ``None``) a positive integral prime

        OUTPUT: integer

        .. WARNING::

            Since only finitely many moments are computed, this valuation may
            be larger than the actual valuation of this distribution.
            Moreover, this valuation may be smaller than the actual
            valuation if all entries are zero to the known precision.

        EXAMPLES::

            sage: D = OverconvergentDistributions(8, 7, 15)
            sage: v = D([7^(5-i) for i in range(1,5)])
            sage: v
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v.valuation(7)
            4
        """
        if p is None:
            p = self.parent()._p
        n = self.precision_relative()
        if self.parent().is_symk():
            return self.ordp + min([self._unscaled_moment(a).valuation(p) for a in range(n)])
        else:
            return self.ordp + min([n] + [self._unscaled_moment(a).valuation(p) for a in range(n) if not self._unscaled_moment(a).is_zero()])

    def specialize(self, new_base_ring=None):
        """
        Return the image of this overconvergent distribution under
        the canonical projection from distributions of weight `k` to
        `Sym^k`.

        INPUT:

        - ``new_base_ring`` -- (default: ``None``) a ring giving the
          desired base ring of the result

        OUTPUT:

        An element of `Sym^k(K)`, where `K` is the specified base ring.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 13)
            sage: d = D([0,2,4,6,8,10,12])
            sage: d.specialize()
            (O(13^7), 2 + O(13^6), 4 + O(13^5), 6 + O(13^4), 8 + O(13^3))
        """
        # self.normalize() # This method should not change self
        k = self.parent()._k
        if k < 0:
            raise ValueError("negative weight")
        if self.precision_absolute() < k + 1:
            raise ValueError("not enough moments")
        V = self.parent().specialize(new_base_ring)
        new_base_ring = V.base_ring()
        if self.precision_relative() == 0:
            return V.zero()
        return V([new_base_ring.coerce(self.moment(j)) for j in range(k + 1)])

    def lift(self, p=None, M=None, new_base_ring=None):
        r"""
        Lift a distribution or element of `Sym^k` to an overconvergent distribution.

        INPUT:

        - ``p`` -- (default: ``None``) a positive integral prime.  If ``None``
          then ``p`` must be available in the parent

        - ``M`` -- (default: ``None``) a positive integer giving the
          desired number of moments. If ``None``, returns a distribution having one
          more moment than this one.

        - ``new_base_ring`` -- (default: ``None``) a ring giving the desired base
          ring of the result. If ``None``, a base ring is chosen automatically.

        OUTPUT:

        - An overconvergent distribution with `M` moments whose image
          under the specialization map is this element.

        EXAMPLES::

            sage: V = Symk(0)
            sage: x = V(1/4)
            sage: y = x.lift(17, 5)
            sage: y
            (13 + 12*17 + 12*17^2 + 12*17^3 + 12*17^4 + O(17^5), O(17^4), O(17^3), O(17^2), O(17))
            sage: y.specialize()._moments == x._moments
            True
        """
        V = self.parent().lift(p, M, new_base_ring)
        k = V._k
        p = V.prime()
        M = V.precision_cap()
        R = V.base_ring()
        moments = [R(self.moment(j)) for j in range(k + 1)]
        zero = R(0)
        moments.extend([zero] * (M - k - 1))
        mu = V(moments)
        #val = mu.valuation()
        #if val < 0:
        #    # This seems unnatural
        #    print("scaling by ", p, "^", -val, " to keep things integral")
        #    mu *= p**(-val)
        return mu

    def _is_malformed(self):
        r"""
        Check that the precision of ``self`` is sensible.

        EXAMPLES::

            sage: D = sage.modular.pollack_stevens.distributions.Symk(2, base=Qp(5))
            sage: v = D([1, 2, 3])
            sage: v._is_malformed()
            False
            sage: v = D([1 + O(5), 2, 3])
            sage: v._is_malformed()
            True
        """
        n = self.precision_absolute()
        for i in range(n):
            if self.moment(i).precision_absolute() < n - i:
                return True
        return False

    def act_right(self, gamma):
        r"""
        The image of this element under the right action by a
        `2 \times 2` matrix.

        INPUT:

        - ``gamma`` -- the matrix by which to act

        OUTPUT:

        - ``self | gamma``

        .. NOTE::

            You may also just use multiplication ``self * gamma``.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4, 7, 10)
            sage: v = D([98,49,21,28,35])
            sage: M = matrix([[1,0], [7,1]])
            sage: v.act_right(M)
            (2*7^2 + 7^3 + 5*7^4 + O(7^5), 3*7^2 + 6*7^3 + O(7^4), 3*7 + 7^2 + O(7^3), 4*7 + O(7^2), O(7))
        """
        return self.parent()._act(self, gamma)

cdef class Dist_vector(Dist):
    r"""
    A distribution is stored as a vector whose `j`-th entry is the `j`-th moment of the distribution.

    The `j`-th entry is stored modulo `p^{N-j}` where `N` is the total number of moments.
    (This is the accuracy that is maintained after acting by `\Gamma_0(p)`.)

    INPUT:

    - ``moments`` -- the list of moments.  If ``check == False`` it
      must be a vector in the appropriate approximation module.

    - ``parent`` -- a :class:`distributions.OverconvergentDistributions_class` or
      :class:`distributions.Symk_class` instance

    - ``ordp`` -- integer;  this *must* be zero in the case of Symk
      of an exact ring

    - ``check`` -- boolean (default: ``True``); whether to validate input

    EXAMPLES::

        sage: D = OverconvergentDistributions(3,5,6) # indirect doctest
        sage: v = D([1,1,1])
    """
    def __init__(self, moments, parent, ordp=0, check=True, normalize=True):
        """
        Initialization.

        TESTS::

            sage: Symk(4)(0)
            (0, 0, 0, 0, 0)
        """
        # if not hasattr(parent,'Element'):
        #     parent, moments = moments, parent

        Dist.__init__(self, parent)
        if check:
            # case 1: input is a distribution already
            if isinstance(moments, Dist):
                ordp = moments._ord_p()
                moments = moments._moments.change_ring(parent.base_ring())
            # case 2: input is a vector, or something with a len
            elif hasattr(moments, '__len__'):
                M = len(moments)
                moments = parent.approx_module(M)(moments)
            # case 3: input is zero
            elif moments == 0:
                moments = parent.approx_module(parent.precision_cap())(moments)
            # case 4: everything else
            else:
                moments = parent.approx_module(1)([moments])
            # TODO: This is not quite right if the input is an inexact zero.
            if ordp != 0 and parent.prime() == 0:
                raise ValueError("cannot specify a valuation shift for an exact ring")

        self._moments = moments
        self.ordp = ordp
        if normalize:
            self.normalize()

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: D = sage.modular.pollack_stevens.distributions.Symk(2)
            sage: x = D([2,3,4])
            sage: x.__reduce__()
            (<class 'sage.modular.pollack_stevens.dist.Dist_vector'>, ((2, 3, 4), Sym^2 Q^2, 0, False))
        """
        return (self.__class__, (self._moments, self.parent(), self.ordp, False))

    cdef Dist_vector _new_c(self):
        r"""
        Create an empty distribution.

        Note that you MUST fill in the ordp attribute on the resulting distribution.

        OUTPUT:

        - A distribution with no moments.  The moments are then filled
          in by the calling function.

        EXAMPLES::

            sage: D = OverconvergentDistributions(3,5,4) # indirect doctest
            sage: v = D([1,1,1])
        """
        cdef Dist_vector ans = Dist_vector.__new__(Dist_vector)
        ans._parent = self._parent
        return ans

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: repr(v)
            '(1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))'
        """
        valstr = ""
        if self.ordp == 1:
            valstr = "%s * " % (self.parent().prime())
        elif self.ordp != 0:
            valstr = "%s^%s * " % (self.parent().prime(), self.ordp)
        if len(self._moments) == 1:
            return valstr + repr(self._moments[0])
        else:
            return valstr + repr(self._moments)

    def _rational_(self):
        """
        Convert to a rational number.

        EXAMPLES::

            sage: D = Symk(0); d = D(4/3); d
            4/3
            sage: QQ(d)
            4/3

        We get a :exc:`TypeError` if there is more than 1 moment::

            sage: D = Symk(1); d = D([1,2]); d
            (1, 2)
            sage: QQ(d)
            Traceback (most recent call last):
            ...
            TypeError: k must be 0
        """
        if len(self._moments) == 1:
            return QQ(self.moment(0))
        raise TypeError("k must be 0")

    cdef long _relprec(self) noexcept:
        """
        Return the number of moments.

        EXAMPLES::

            sage: D = Symk(4)
            sage: d = D([1,2,3,4,5]); e = D([2,3,4,5,6])
            sage: d == e # indirect doctest
            False
        """
        return len(self._moments)

    cdef _unscaled_moment(self, long n):
        r"""
        Return the `n`-th moment, unscaled by the overall power of `p`
        stored in ``self.ordp``.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4,3,5)
            sage: d = D([3,3,3,3,3])
            sage: d.moment(2) # indirect doctest
            3 + O(3^3)
        """
        return self._moments[n]

    cdef Dist_vector _addsub(self, Dist_vector right, bint negate):
        r"""
        Common code for the sum and the difference of two distributions.

        EXAMPLES::

            sage: D = Symk(2)
            sage: u = D([1,2,3]); v = D([4,5,6])
            sage: u + v # indirect doctest
            (5, 7, 9)
            sage: u - v # indirect doctest
            (-3, -3, -3)
        """
        cdef Dist_vector ans = self._new_c()
        cdef long aprec = min(self.ordp + len(self._moments), right.ordp + len(right._moments))
        ans.ordp = min(self.ordp, right.ordp)
        cdef long rprec = aprec - ans.ordp
        # In the case of symk, rprec will always be k
        V = ans.parent().approx_module(rprec)
        R = V.base_ring()
        smoments = self._moments
        rmoments = right._moments
        # We truncate if the moments are too long; extend by zero if too short
        if smoments.parent() is not V:
            vec = smoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(smoments)) if rprec > len(smoments) else [])
            smoments = V(vec)
        if rmoments.parent() is not V:
            vec = rmoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(rmoments)) if rprec > len(rmoments) else [])
            rmoments = V(vec)
        # We multiply by the relative power of p
        if self.ordp > right.ordp:
            smoments *= self.parent().prime() ** (self.ordp - right.ordp)
        elif self.ordp < right.ordp:
            rmoments *= self.parent().prime() ** (right.ordp - self.ordp)
        if negate:
            rmoments = -rmoments
        ans._moments = smoments + rmoments
        return ans

    cpdef _add_(self, _right):
        r"""
        Sum of two distributions.

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); w = D([3,6,9,12,15])
            sage: v+w
            (4 + O(7^5), 1 + 7 + O(7^4), 5 + 7 + O(7^3), 2 + 2*7 + O(7^2), 6 + O(7))
        """
        return self._addsub(<Dist_vector>_right, False)

    cpdef _sub_(self, _right):
        r"""
        Difference of two distributions.

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); w = D([1,1,1,8,8])
            sage: v-w
            (O(7^5), 1 + O(7^4), 2 + O(7^3), 3 + 6*7 + O(7^2), 4 + O(7))
        """
        return self._addsub(<Dist_vector>_right, True)

    cpdef _lmul_(self, Element right):
        r"""
        Scalar product of a distribution with a ring element that coerces into the base ring.

        EXAMPLES::

            sage: D = OverconvergentDistributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: 3*v; 7*v
            (3 + O(7^5), 6 + O(7^4), 2 + 7 + O(7^3), 5 + 7 + O(7^2), 1 + O(7))
            7 * (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v*3; v*7
            (3 + O(7^5), 6 + O(7^4), 2 + 7 + O(7^3), 5 + 7 + O(7^2), 1 + O(7))
            7 * (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
        """
        cdef Dist_vector ans = self._new_c()
        p = self.parent().prime()
        if p == 0:
            ans._moments = self._moments * right
            ans.ordp = self.ordp
        elif right.valuation(p) == Infinity:
            ans._moments = self.parent().approx_module(0)([])
            ans.ordp += self.precision_relative()
        else:
            try:
                v, u = right.val_unit(p)
            except TypeError:  # bug in p-adics: they should accept p here
                v, u = right.val_unit()
            ans._moments = self._moments * u
            ans.ordp = self.ordp + v
            # if the relative precision of u is less than that of self, ans may not be normalized.
        return ans

    def precision_relative(self):
        r"""
        Return the relative precision of this distribution.

        The precision is just the number of moments stored, which is
        also `k+1` in the case of `Sym^k(R)`.  For overconvergent
        distributions, the precision is the integer `m` so that the
        sequence of moments is known modulo `Fil^m`.

        OUTPUT: integer giving the number of moments

        EXAMPLES::

            sage: D = OverconvergentDistributions(2, 11, 15)
            sage: v = D([1,1,10,9,6,15])
            sage: v.precision_relative()
            6
            sage: v = v.reduce_precision(4); v.precision_relative()
            4
            sage: D = Symk(10)
            sage: v = D.random_element()
            sage: v.precision_relative()
            11
        """
        return Integer(len(self._moments))

    def precision_absolute(self):
        r"""
        Return the absolute precision of this distribution.

        The absolute precision is the sum of the relative precision
        (number of moments) and the valuation.

        EXAMPLES::

            sage: D = OverconvergentDistributions(3, 7, base = Qp(7))
            sage: v = D([3,1,10,0])
            sage: v.precision_absolute()
            4
            sage: v *= 7
            sage: v.precision_absolute()
            5
            sage: v = 1/7^10 * v
            sage: v.precision_absolute()
            -5
        """
        return Integer(len(self._moments) + self.ordp)

    cpdef normalize(self, include_zeroth_moment=True):
        r"""
        Normalize by reducing modulo `Fil^N`, where `N` is the number of moments.

        If the parent is Symk, then normalize has no effect.  If the
        parent is a space of distributions, then normalize reduces the
        `i`-th moment modulo `p^{N-i}`.

        OUTPUT: this distribution, after normalizing

        .. WARNING::

            This function modifies the distribution in place as well as returning it.

        EXAMPLES::

            sage: D = OverconvergentDistributions(3,7,10)
            sage: v = D([1,2,3,4,5]) ; v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: w = v.reduce_precision(3) ; w
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3))
            sage: w.normalize()
            (1 + O(7^3), 2 + O(7^2), 3 + O(7))
            sage: w
            (1 + O(7^3), 2 + O(7^2), 3 + O(7))
            sage: v.reduce_precision(3).normalize(include_zeroth_moment=False)
            (1 + O(7^5), 2 + O(7^2), 3 + O(7))
        """
        if not self.parent().is_symk() and self._moments != 0:  # non-classical
            if not self._moments:
                return self
            V = self._moments.parent()
            R = V.base_ring()
            n = self.precision_relative()
            p = self.parent()._p
            shift = self.ordp
            if include_zeroth_moment:
                if isinstance(R, pAdicGeneric):
                    self._moments = V([self._moments[i].add_bigoh(n -shift - i) for i in range(n)])
                else:
                    self._moments = V([self._moments[i] % (p ** (n -shift - i)) for i in range(n)])
            else:
                if isinstance(R, pAdicGeneric):
                    self._moments = V([self._moments[0]] + [self._moments[i].add_bigoh(n -shift - i) for i in range(1, n)])  # Don't normalize the zeroth moment
                else:
                    self._moments = V([self._moments[0]] + [self._moments[i] % (p ** (n -shift- i)) for i in range(1, n)])  # Don't normalize the zeroth moment
        return self

    def reduce_precision(self, M):
        r"""
        Only hold on to `M` moments.

        INPUT:

        - ``M`` -- positive integer less than the precision of this
          distribution

        OUTPUT:

        - a new distribution with `M` moments equal to the first `M`
          moments of this distribution.

        EXAMPLES::

            sage: D = OverconvergentDistributions(3,7,10)
            sage: v = D([3,4,5])
            sage: v
            (3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.reduce_precision(2)
            (3 + O(7^3), 4 + O(7^2))
        """
        assert M <= self.precision_relative(), "not enough moments"

        cdef Dist_vector ans = self._new_c()
        ans._moments = self._moments[:M]
        ans.ordp = self.ordp
        return ans

    def solve_difference_equation(self):
        r"""
        Solve the difference equation. `self = v | \Delta`, where `\Delta = [1, 1; 0, 1] - 1`.

        See Theorem 4.5 and Lemma 4.4 of [PS2011]_.

        OUTPUT:

        - a distribution `v` so that `self = v | Delta` , assuming ``self.moment(0) == 0``.
          Otherwise solves the difference equation for ``self - (self.moment(0),0,...,0)``.

        EXAMPLES::

            sage: D = OverconvergentDistributions(5,7,15)
            sage: v = D(([0,2,3,4,5]))
            sage: g = D._act.actor()(Matrix(ZZ,2,2,[1,1,0,1]))
            sage: w = v.solve_difference_equation()
            sage: v - (w*g - w)
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v = D(([7,2,3,4,5]))
            sage: w = v.solve_difference_equation()
            sage: v - (w*g - w)
            (7 + O(7^4), O(7^3), O(7^2), O(7))
        """
        # assert self._moments[0][0]==0, "not total measure zero"
        # print("result accurate modulo p^",self.moment(0).valuation(self.p) )
        #v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
        M = self.precision_relative()
        R = self.parent().base_ring()
        K = R.fraction_field()
        V = self._moments.parent()
        v = [K(0) for i in range(M)]
        bern = [bernoulli(i) for i in range(0, M, 2)]
        minhalf = ~K(-2)
        for m in range(1, M):
            scalar = K(self.moment(m) / m)
            # bernoulli(1) = -1/2; the only nonzero odd Bernoulli number
            v[m] += m * minhalf * scalar
            for j in range(m - 1, M, 2):
                v[j] += binomial(j, m - 1) * bern[(j - m + 1) // 2] * scalar
        p = self.parent().prime()
        cdef Dist_vector ans
        if p == 0:
            if R in Fields():
                ans = self._new_c()
                ans.ordp = 0
                ans._moments = V(v)
            else:
                newparent = self.parent().change_ring(K)
                ans = newparent(v)
        else:
            ans = self._new_c()
            try:
                ans.ordp = min(a.valuation(p) for a in v)
            except TypeError:
                ans.ordp = 0
            if ans.ordp < 0:
                scalar = K(p) ** (-ans.ordp)
                ans._moments = V([R(a * scalar) for a in v])
            elif ans.ordp > 0:
                scalar = K(p) ** ans.ordp
                ans._moments = V([R(a // scalar) for a in v])
            else:
                ans._moments = V([R(a) for a in v])
            v = ans._moments
            N = len(ans._moments)
            prec_loss = max([N - j - v[j].precision_absolute()
                             for j in range(N)])
            #            print("precision loss = ", prec_loss)
            if prec_loss > 0:
                ans._moments = ans._moments[:(N - prec_loss)]
        return ans


cdef class WeightKAction(Action):
    r"""
    Encode the action of the monoid `\Sigma_0(N)` on the space of distributions.

    INPUT:

    - ``Dk`` -- a space of distributions
    - ``character`` -- data specifying a Dirichlet character to apply to
      the top right corner, and a power of the determinant by which to scale.
      See the documentation of
      :class:`sage.modular.pollack_stevens.distributions.OverconvergentDistributions_factory`
      for more details.
    - ``adjuster`` -- a callable object that turns matrices into 4-tuples
    - ``on_left`` -- whether this action should be on the left
    - ``dettwist`` -- a power of the determinant to twist by
    - ``padic`` -- if ``True``, define an action of `p`-adic matrices (not just integer ones)

    EXAMPLES::

        sage: D = OverconvergentDistributions(4,5,10,base = Qp(5,20)); D
        Space of 5-adic distributions with k=4 action and precision cap 10
        sage: D._act
        Right action by Monoid Sigma0(5) with coefficients in 5-adic Field with capped relative precision 20 on Space of 5-adic distributions with k=4 action and precision cap 10
    """
    def __init__(self, Dk, character, adjuster, on_left, dettwist, padic=False):
        r"""
        Initialization.

        TESTS::

            sage: D = OverconvergentDistributions(4,5,10,base = Qp(5,20)); D # indirect doctest
            Space of 5-adic distributions with k=4 action and precision cap 10
            sage: D = Symk(10) # indirect doctest
        """
        self._k = Dk._k
#        if self._k < 0: raise ValueError("k must not be negative")
        self._adjuster = adjuster
        self._character = character
        self._dettwist = dettwist
        self._p = Dk._p
        self._symk = Dk.is_symk()
        self._actmat = {}
        self._maxprecs = {}
        if character is None:
            self._Np = ZZ(1)  # all of M2Z acts
        else:
            self._Np = character.modulus()
        if not self._symk:
            self._Np = self._Np.lcm(self._p)

        if padic:
            self._Sigma0 = Sigma0(self._Np, base_ring=Dk.base_ring(), adjuster=self._adjuster)
        else:
            self._Sigma0 = Sigma0(self._Np, base_ring=ZZ, adjuster=self._adjuster)
        Action.__init__(self, self._Sigma0, Dk, on_left, operator.mul)

    def clear_cache(self):
        r"""
        Clear the cached matrices which define the action of `U_p`
        (these depend on the desired precision) and the
        dictionary that stores the maximum precisions computed so far.

        EXAMPLES::

            sage: D = OverconvergentDistributions(4,5,4)
            sage: D([1,2,5,3]) * D._act.actor()(Matrix(ZZ,2,2,[1,1,0,1]))
            (1 + O(5^4), 3 + O(5^3), 2*5 + O(5^2), O(5))
            sage: D._act.clear_cache()
        """
        self._actmat = {}
        self._maxprecs = {}

    cpdef acting_matrix(self, g, M):
        r"""
        The matrix defining the action of ``g`` at precision ``M``.

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrix.matrix_generic_dense.Matrix_generic_dense`

        - ``M`` -- positive integer giving the precision at which
          ``g`` should act

        OUTPUT:

        - An `M \times M` matrix so that the action of `g` on a
          distribution with `M` moments is given by a vector-matrix
          multiplication.

        .. NOTE::

            This function caches its results.  To clear the cache use
            :meth:`clear_cache`.

        EXAMPLES::

            sage: D = Symk(3)
            sage: v = D([5,2,7,1])
            sage: g = Matrix(ZZ,2,2,[1,3,0,1])
            sage: v * D._act.actor()(g) # indirect doctest
            (5, 17, 64, 253)
        """
        g = g.matrix()
        if g not in self._maxprecs:
            A = self._compute_acting_matrix(g, M)
            self._actmat[g] = {M: A}
            self._maxprecs[g] = M
            return A
        else:
            mats = self._actmat[g]
            if M in mats:
                return mats[M]
            maxprec = self._maxprecs[g]
            if M < maxprec:
                A = mats[maxprec][:M, :M]  # submatrix; might want to reduce precisions
                mats[M] = A
                return A
            if M < 30:  # This should not be hard-coded
                maxprec = max([M, 2 * maxprec])  # This may be wasting memory
            else:
                maxprec = M
            self._maxprecs[g] = maxprec
            mats[maxprec] = self._compute_acting_matrix(g, maxprec)  # could lift from current maxprec
            if M == maxprec:
                return mats[maxprec]
            A = mats[maxprec][:M, :M]  # submatrix; might want to reduce precisions
            mats[M] = A
            return A

    cpdef _compute_acting_matrix(self, g, M):
        r"""
        Compute the matrix defining the action of ``g`` at precision ``M``.

        INPUT:

        - ``g`` -- a `2 \times 2` instance of
          :class:`sage.matrices.matrix_integer_dense.Matrix_integer_dense`

        - ``M`` -- positive integer giving the precision at which
          ``g`` should act

        OUTPUT:

        - ``G`` -- an `M \times M` matrix. If `v `is the vector of moments of a
          distribution `\mu`, then `v*G` is the vector of moments of `\mu|[a,b;c,d]`

        EXAMPLES::

            sage: D = Symk(3)
            sage: v = D([5,2,7,1])
            sage: g = Matrix(ZZ,2,2,[-2,1,-1,0])
            sage: v * D._act.actor()(g) # indirect doctest
            (-107, 35, -12, 5)
        """
        raise NotImplementedError


cdef class WeightKAction_vector(WeightKAction):
    cpdef _compute_acting_matrix(self, g, M):
        r"""
        Compute the matrix defining the action of ``g`` at precision ``M``.

        INPUT:

        - ``g`` -- a `2 \times 2` instance of
          :class:`sage.matrix.matrix_generic_dense.Matrix_generic_dense`

        - ``M`` -- positive integer giving the precision at which
          ``g`` should act

        OUTPUT:

        - ``G`` -- an `M \times M` matrix. If `v` is the vector of moments of a
          distribution `\mu`, then `v*G` is the vector of moments of `\mu|[a,b;c,d]`

        EXAMPLES::

            sage: D = Symk(3)
            sage: v = D([5,2,7,1])
            sage: g = Matrix(ZZ,2,2,[-2,1,-1,0])
            sage: v * D._act.actor()(g) # indirect doctest
            (-107, 35, -12, 5)
        """
        #tim = verbose("Starting")
        a, b, c, d = self._adjuster(g)
        # if g.parent().base_ring().is_exact():
        #     self._check_mat(a, b, c, d)
        k = self._k
        if g.parent().base_ring() is ZZ:
            if self._symk:
                base_ring = QQ
            else:
                base_ring = Zmod(self._p ** M)
        else:
            base_ring = self.underlying_set().base_ring()
        cdef Matrix B = matrix(base_ring, M, M)
        if M == 0:
            return B.change_ring(self.codomain().base_ring())
        R = PowerSeriesRing(base_ring, 'y', default_prec=M)
        y = R.gen()
        #tim = verbose("Checked, made R",tim)
        # special case for small precision, large weight
        scale = (b + d * y) / (a + c * y)
        t = (a + c * y) ** k  # will already have precision M
        cdef long row, col
        #tim = verbose("Made matrix",tim)
        for col in range(M):
            for row in range(M):
                B.set_unsafe(row, col, t[row])
            t *= scale
        #verbose("Finished loop",tim)
        # the changering here is annoying, but otherwise we have to
        # change ring each time we multiply
        B = B.change_ring(self.codomain().base_ring())
        if self._character is not None:
            B *= self._character(a)
        if self._dettwist is not None:
            B *= (a * d - b * c) ** (self._dettwist)
        return B

    cpdef _act_(self, g, _v):
        r"""
        The right action of ``g`` on a distribution.

        INPUT:

        - ``_v`` -- a :class:`Dist_vector` instance, the distribution
          on which to act

        - ``g`` -- a `2 \times 2` instance of
          :class:`sage.matrix.matrix_integer_dense.Matrix_integer_dense`

        OUTPUT: the distribution ``_v * g``

        EXAMPLES::

            sage: D = sage.modular.pollack_stevens.distributions.Symk(2)
            sage: v = D([2,3,4])
            sage: g = Matrix(ZZ,2,2,[3,-1,1,0])
            sage: v * D._act.actor()(g) # indirect doctest
            (40, -9, 2)
        """
        # if g is a matrix it needs to be immutable
        # hashing on arithmetic_subgroup_elements is by str
        if g == 1:
            return _v
        cdef Dist_vector v = <Dist_vector?>_v
        cdef Dist_vector ans = v._new_c()

        try:
            g.set_immutable()
        except AttributeError:
            pass
        v_moments = v._moments
        ans._moments = v_moments * self.acting_matrix(g, len(v_moments))
        ans.ordp = v.ordp
        return ans
