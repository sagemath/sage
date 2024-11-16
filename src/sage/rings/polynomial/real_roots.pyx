# sage.doctest: needs numpy sage.libs.linbox
"""
Isolate Real Roots of Real Polynomials

AUTHOR:

- Carl Witty (2007-09-19): initial version

This is an implementation of real root isolation.  That is, given a
polynomial with exact real coefficients, we compute isolating intervals
for the real roots of the polynomial. (Polynomials with
integer, rational, or algebraic real coefficients are supported.)

We convert the polynomials into the Bernstein basis, and then use
de Casteljau's algorithm and Descartes' rule of signs on the Bernstein
basis polynomial (using interval arithmetic) to locate the roots. The
algorithm is similar to that in "A Descartes Algorithm for Polynomials
with Bit-Stream Coefficients", by Eigenwillig, Kettner, Krandick, Mehlhorn,
Schmitt, and Wolpert, but has three crucial optimizations over the
algorithm in that paper:

- Precision reduction: at certain points in the computation, we discard the
  low-order bits of the coefficients, widening the intervals.

- Degree reduction: at certain points in the computation, we find lower-degree
  polynomials that are approximately equal to our high-degree polynomial over
  the region of interest.

- When the intervals are too wide to continue (either because of a too-low
  initial precision, or because of precision or degree reduction), and we need
  to restart with higher precision, we recall which regions have already been
  proven not to have any roots and do not examine them again.

The best description of the algorithms used (other than this source
code itself) is in the slides for my Sage Days 4 talk, currently available
from https://wiki.sagemath.org/days4schedule .
"""

################################################################################
#       Copyright (C) 2007 Carl Witty <cwitty@newtonlabs.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

# TODO:
# These things would almost certainly improve the speed:
# * Use Anatole Ruslanov's register tiling versions of de Casteljau's
# algorithms when doing de Casteljau splitting at 1/2 in the integer case.
# * Use a register tiling version of de Casteljau's algorithm for the
# floating-point case...if you have n free FP registers (after allocating
# registers to hold r and 1-r) you should be able to do about n rows at
# once.
# * Use vectorized FP instructions for de Casteljau's algorithm.
# Using SSE2 to do 2 FP operations at once should be twice as fast.
# * Faster degree elevation
# Currently, when bernstein_expand expands from a degree-d1 to a
# degree-d2 polynomial, it does O(d2^2 - d1^2) operations.  I
# have thought of (but not implemented) an approach that takes
# O(d1*d2) operations, which should be much faster when d1 is 30
# and d2 is 1000.
# Basically, the idea is to compute the d1 derivatives of the polynomial
# at x=0.  This can be done exactly, in O(d1^2) operations.  Then
# lift the d1'th derivative (a constant polynomial) to a polynomial of
# formal degree d2-d1 (trivial; the coefficients are just d2-d1+1 copies
# of the derivative).  Then, compute the integral of that polynomial
# d1 times, using the computed derivatives to give the constant offset.
# (Computing the integral in Bernstein form involves divisions; these can
# be done ahead of time, to the derivatives at x=0.)
# The only tricky bit is tracking the error bounds and making sure you
# use appropriate precisions at various parts of the computation.

# These things are interesting ideas that may or may not help:
# * Partial degree reduction
# Build a new subclass of interval_bernstein_polynomial, which
# represents a polynomial as the sum of several
# interval_bernstein_polynomial_integer (of different degrees).
# Suppose you have a degree-1000 polynomial with 2000-bit coefficients.
# When you try degree reduction, you get a degree-30 polynomial with
# 2000-bit coefficients; but you get 300 bits of error.
# Instead of just using the degree-30 polynomial and accepting the
# loss of 300 bits of precision, you could represent the reduced polynomial
# as the sum of a degree-30 polynomial with 2000-bit coefficients and
# a degree-1000 polynomial with 300-bit coefficients.  The combined
# polynomial should still be much cheaper to work with than the original.
# * Faster degree reduction
# Perhaps somebody clever can figure out a way to do degree reduction
# that doesn't involve matrix inversion, so that degree reduction can
# target degrees larger than 30.
# * Better heuristics
# ** Better degree reduction heuristics
# When should we degree reduce?  What degrees should be targeted?
# ** Better precision reduction heuristics
# ** Better split point selection
# Currently we attempt to split at 1/2, and then at a random offset.
# We do not look at the current coefficients at all to try to select
# a good split point.  Would it be worthwhile to do so?  (For instance,
# it's a standard result that the polynomial in the region is bounded
# by the convex hull of the coefficients.)
# * Better initial transformation into Bernstein form
# If a polynomial has roots which are very close to 0, or very large,
# it may be better to rescale the roots (using a transformation like
# p -> p(1000*x)).  This can be done as part of the initial transformation
# into Bernstein form.
# If a polynomial has some roots which are very close to 0, and also
# some roots which are very large, it might be better to isolate the roots
# in two groups.
# * More kinds of interval_bernstein_polynomial
# Currently we have arbitrary-precision GMP coefficients or double-precision
# floats.  Would choices between these extremes help?  (double-double?
# quad-double?  fixed-precision multi-word integer arithmetic?)
# * Currently, some of our competitors are faster than this algorithm
# for polynomials that are "easy" in some way.  (This algorithm
# seems to be vastly faster than everything else for "hard" polynomials.)
# Maybe there is a way to start out using another algorithm
# and then switch over to this algorithm once the polynomial
# is determined to be hard.
# * standard optimization
# There's probably some more to be gained by just running the profiler
# and optimizing individual functions.  (To use more Cython features, for
# instance.)

# Extra features:
# * Specify a minimal island width: once islands are this narrow,
# stop even if roots are not isolated.
# * Do some sort of inexact version, for inexact input polynomials.
# (For example, given a polynomial with (non-point) interval coefficients,
# give a set of roots such that each root is guaranteed to be a root
# of some polynomial in the set represented by the interval polynomial.
# This may be vastly faster than the exact calculations carried out
# by this algorithm!  Is it enough faster to be faster than, say,
# Pari's floating-point algorithms?)

from copy import copy
import time

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import infinity
from sage.rings.qqbar import AA
from sage.rings.real_double import RDF
from sage.rings.real_mpfi import RealIntervalField, RIF
from sage.rings.real_mpfr import RR, RealField
from sage.arith.misc import binomial, factorial
from sage.misc.randstate import randstate
from sage.modules.free_module_element import free_module_element as vector
from sage.modules.free_module import FreeModule
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.functional import numerator, denominator
from sage.misc.misc_c import prod

from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.modules.vector_real_double_dense cimport Vector_real_double_dense
from sage.rings.integer cimport Integer
from sage.rings.real_mpfr cimport RealNumber

cimport numpy

from libc.math cimport fabs, sqrt, ldexp, frexp

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.mpfr cimport *

cdef class interval_bernstein_polynomial:
    r"""
    An interval_bernstein_polynomial is an approximation to an exact
    polynomial.  This approximation is in the form of a Bernstein
    polynomial (a polynomial given as coefficients over a Bernstein
    basis) with interval coefficients.

    The Bernstein basis of degree n over the region [a .. b] is the
    set of polynomials

    .. MATH::

      \binom{n}{k} (x-a)^k (b-x)^{n-k} / (b-a)^n

    for `0 \le k \le n`.

    A degree-n interval Bernstein polynomial P with its region [a .. b] can
    represent an exact polynomial p in two different ways: it can
    "contain" the polynomial or it can "bound" the polynomial.

    We say that P contains p if, when p is represented as a degree-n
    Bernstein polynomial over [a .. b], its coefficients are contained
    in the corresponding interval coefficients of P.  For instance,
    [0.9 .. 1.1]*x^2 (which is a degree-2 interval Bernstein polynomial
    over [0 .. 1]) contains x^2.

    We say that P bounds p if, for all a <= x <= b, there exists a
    polynomial p' contained in P such that p(x) == p'(x).  For instance,
    [0 .. 1]*x is a degree-1 interval Bernstein polynomial which bounds
    x^2 over [0 .. 1].

    If P contains p, then P bounds p; but the converse is not necessarily
    true.  In particular, if n < m, it is possible for a degree-n interval
    Bernstein polynomial to bound a degree-m polynomial; but it cannot
    contain the polynomial.

    In the case where P bounds p, we maintain extra information, the
    "slope error".  We say that P (over [a .. b]) bounds p with a
    slope error of E (where E is an interval) if there is a polynomial
    p' contained in P such that the derivative of (p - p') is bounded
    by E in the range [a .. b].  If P bounds p with a slope error of 0
    then P contains p.

    (Note that "contains" and "bounds" are not standard terminology;
    I just made them up.)

    Interval Bernstein polynomials are useful in finding real roots
    because of the following properties:

    - Given an exact real polynomial p, we can compute an interval Bernstein
      polynomial over an arbitrary region containing p.

    - Given an interval Bernstein polynomial P over [a .. c], where a < b < c,
      we can compute interval Bernstein polynomials P1 over [a .. b] and P2
      over [b .. c], where P1 and P2 contain (or bound) all polynomials that P
      contains (or bounds).

    - Given a degree-n interval Bernstein polynomial P over [a .. b], and m <
      n, we can compute a degree-m interval Bernstein polynomial P' over [a ..
      b] that bounds all polynomials that P bounds.

    - It is sometimes possible to prove that no polynomial bounded by P over [a
      .. b] has any roots in [a .. b].  (Roughly, this is possible when no
      polynomial contained by P has any complex roots near the line segment [a
      .. b], where "near" is defined relative to the length b-a.)

    - It is sometimes possible to prove that every polynomial bounded by P over
      [a .. b] with slope error E has exactly one root in [a .. b].  (Roughly,
      this is possible when every polynomial contained by P over [a .. b] has
      exactly one root in [a .. b], there are no other complex roots near the
      line segment [a .. b], and every polynomial contained in P has a
      derivative which is bounded away from zero over [a .. b] by an amount
      which is large relative to E.)

    - Starting from a sufficiently precise interval Bernstein polynomial, it is
      always possible to split it into polynomials which provably have 0 or 1
      roots (as long as your original polynomial has no multiple real roots).

    So a rough outline of a family of algorithms would be:

    - Given a polynomial p, compute a region [a .. b] in which any real roots
      must lie.
    - Compute an interval Bernstein polynomial P containing p over [a .. b].
    - Keep splitting P until you have isolated all the roots.  Optionally,
      reduce the degree or the precision of the interval Bernstein polynomials
      at intermediate stages (to reduce computation time).  If this seems not
      to be working, go back and try again with higher precision.

    Obviously, there are many details to be worked out to turn this
    into a full algorithm, like:

    - What initial precision is selected for computing P?
    - How do you decide when to reduce the degree of intermediate polynomials?
    - How do you decide when to reduce the precision of intermediate
      polynomials?
    - How do you decide where to split the interval Bernstein polynomial
      regions?
    - How do you decide when to give up and start over with higher precision?

    Each set of answers to these questions gives a different algorithm
    (potentially with very different performance characteristics), but all of
    them can use this ``interval_bernstein_polynomial`` class as their basic
    building block.

    To save computation time, all coefficients in an
    ``interval_bernstein_polynomial`` share the same interval width.
    (There is one exception: when creating an ``interval_bernstein_polynomial``,
    the first and last coefficients can be marked as "known positive"
    or "known negative".  This has some of the same effect as having
    a (potentially) smaller interval width for these two coefficients,
    although it does not affect de Casteljau splitting.)
    To allow for widely varying coefficient magnitudes, all
    coefficients in an interval_bernstein_polynomial are scaled
    by `2^n` (where `n` may be positive, negative, or zero).

    There are two representations for interval_bernstein_polynomials,
    integer and floating-point. These are the two subclasses of
    this class; ``interval_bernstein_polynomial`` itself is an abstract
    class.

    ``interval_bernstein_polynomial`` and its subclasses are not expected
    to be used outside this file.
    """

    def variations(self):
        """
        Consider a polynomial (written in either the normal power basis
        or the Bernstein basis).  Take its list of coefficients, omitting
        zeroes.  Count the number of positions in the list where the
        sign of one coefficient is opposite the sign of the next coefficient.

        This count is the number of sign variations of the polynomial.
        According to Descartes' rule of signs, the number of real
        roots of the polynomial (counted with multiplicity) in a
        certain interval is always less than or equal to the number of
        sign variations, and the difference is always even.  (If the
        polynomial is written in the power basis, the region is the
        positive reals; if the polynomial is written in the Bernstein
        basis over a particular region, then we count roots in that region.)

        In particular, a polynomial with no sign variations has no real
        roots in the region, and a polynomial with one sign variation
        has one real root in the region.

        In an interval Bernstein polynomial, we do not necessarily
        know the signs of the coefficients (if some of the coefficient
        intervals contain zero), so the polynomials contained by
        this interval polynomial may not all have the same number
        of sign variations.  However, we can compute a range of
        possible numbers of sign variations.

        This function returns the range, as a 2-tuple of integers.
        """
        return (self.min_variations, self.max_variations)

    cdef void update_variations(self, interval_bernstein_polynomial bp1, interval_bernstein_polynomial bp2) noexcept:
        """
        Update the max_variations of ``bp1`` and ``bp2`` (which are assumed to
        be the result of splitting this polynomial).

        If we knew the number of variations of ``self``, ``bp1``, and ``bp2``
        exactly, we would have
          ``self.variations == bp1.variations + bp2.variations + 2*n``
        for some nonnegative integer n.  Thus, we can use our information
        on min and max variations on ``self`` and ``bp1`` (or ``bp2``) to
        refine the range on ``bp2`` (or ``bp1``).
        """
        if self.max_variations - bp1.min_variations < bp2.max_variations:
            bp2.max_variations = self.max_variations - bp1.min_variations
        if self.max_variations - bp2.min_variations < bp1.max_variations:
            bp1.max_variations = self.max_variations - bp2.min_variations

    def try_split(self, context ctx, logging_note):
        """
        Try doing a de Casteljau split of this polynomial at 1/2, resulting
        in polynomials p1 and p2.  If we see that the sign of this polynomial
        is determined at 1/2, then return (p1, p2, 1/2); otherwise,
        return None.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([50, 20, -90, -70, 200], error=5)
            sage: bp1, bp2, _ = bp.try_split(mk_context(), None)
            sage: bp1
            <IBP: (50, 35, 0, -29, -31) + [0 .. 6) over [0 .. 1/2]>
            sage: bp2
            <IBP: (-31, -33, -8, 65, 200) + [0 .. 6) over [1/2 .. 1]>
            sage: bp = mk_ibpf([0.5, 0.2, -0.9, -0.7, 0.99], neg_err=-0.1, pos_err=0.01)
            sage: bp1, bp2, _ = bp.try_split(mk_context(), None)
            sage: bp1
            <IBP: (0.5, 0.35, 0.0, -0.2875, -0.369375) + [-0.10000000000000023 .. 0.010000000000000226] over [0 .. 1/2]>
            sage: bp2
            <IBP: (-0.369375, -0.45125, -0.3275, 0.14500000000000002, 0.99) + [-0.10000000000000023 .. 0.010000000000000226] over [1/2 .. 1]>
        """
        (p1, p2, ok) = self.de_casteljau(ctx, QQ_1_2)
        ctx.dc_log_append(("half" + self._type_code(), self.scale_log2, self.bitsize, ok, logging_note))
        if ok:
            return (p1, p2, QQ_1_2)
        return None

    def try_rand_split(self, context ctx, logging_note):
        """
        Compute a random split point r (using the random number generator
        embedded in ctx).  We require 1/4 <= r < 3/4 (to ensure that
        recursive algorithms make progress).

        Then, try doing a de Casteljau split of this polynomial at r, resulting
        in polynomials p1 and p2.  If we see that the sign of this polynomial
        is determined at r, then return (p1, p2, r); otherwise,
        return None.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([50, 20, -90, -70, 200], error=5)
            sage: bp1, bp2, _ = bp.try_rand_split(mk_context(), None)
            sage: bp1
            <IBP: (50, 29, -27, -56, -11) + [0 .. 6) over [0 .. 43/64]>
            sage: bp2
            <IBP: (-11, 10, 49, 111, 200) + [0 .. 6) over [43/64 .. 1]>
            sage: bp1, bp2, _ = bp.try_rand_split(mk_context(seed=42), None)
            sage: bp1
            <IBP: (50, 32, -11, -41, -29) + [0 .. 6) over [0 .. 583/1024]>
            sage: bp2
            <IBP: (-29, -20, 13, 83, 200) + [0 .. 6) over [583/1024 .. 1]>
            sage: bp = mk_ibpf([0.5, 0.2, -0.9, -0.7, 0.99], neg_err=-0.1, pos_err=0.01)
            sage: bp1, bp2, _ = bp.try_rand_split(mk_context(), None)
            sage: bp1  # rel tol
            <IBP: (0.5, 0.2984375, -0.2642578125, -0.5511661529541015, -0.3145806974172592) + [-0.10000000000000069 .. 0.010000000000000677] over [0 .. 43/64]>
            sage: bp2  # rel tol
            <IBP: (-0.3145806974172592, -0.19903896331787108, 0.04135986328125002, 0.43546875, 0.99) + [-0.10000000000000069 .. 0.010000000000000677] over [43/64 .. 1]>
        """

        # We want a split point which is a dyadic rational (denominator
        # is a power of 2), because that speeds up de_casteljau.  We want
        # a small denominator, because that helps us find simpler isolating
        # intervals.  However, if we require that our denominator be too
        # small, then that might keep the algorithm from terminating;
        # if the current polynomial has roots at (4/16, 5/16, ..., 12/16),
        # and we pick our split points from that same set, then we will
        # never find a split point with a determined sign.  We avoid this
        # problem by making sure we have more possible split points
        # to choose from than our polynomial has roots.

        # A different algorithm here might be more efficient.

        div = 1024
        while self.degree() >= (div // 4):
            div = div * 2
        qdiv = div // 4  # divides evenly since div = 1024*2^k
        rand = Integer(ctx.random.randrange(qdiv, 3*qdiv)) / div
        (p1, p2, ok) = self.de_casteljau(ctx, rand)
        ctx.dc_log_append(("div" + self._type_code(), self.scale_log2, self.bitsize, rand, ok, logging_note))
        if ok:
            return (p1, p2, rand)
        return None

    cdef int degree(self) noexcept:
        raise NotImplementedError()

    def region(self):
        return (self.lower, self.upper)

    def region_width(self):
        return (self.upper - self.lower)

dr_cache = {}

cdef class interval_bernstein_polynomial_integer(interval_bernstein_polynomial):
    """
    This is the subclass of :class:`interval_bernstein_polynomial` where
    polynomial coefficients are represented using integers.

    In this integer representation, each coefficient is represented by
    a GMP arbitrary-precision integer A, and a (shared) interval width
    E (which is a machine integer).  These represent the coefficients
    A*2^n <= c < (A+E)*2^n.

    (Note that :func:`mk_ibpi is a simple helper` function for creating
    elements of :class:`interval_bernstein_polynomial_integer` in doctests.)

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: bp = mk_ibpi([1, 2, 3], error=5); print(bp)
        degree 2 IBP with 2-bit coefficients
        sage: bp
        <IBP: (1, 2, 3) + [0 .. 5)>
        sage: bp.variations()
        (0, 0)
        sage: bp = mk_ibpi([-3, -1, 1, -1, -3, -1], lower=1, upper=5/4, usign=1,        # needs sage.symbolic
        ....:              error=2, scale_log2=-3, level=2, slope_err=RIF(pi)); print(bp)
        degree 5 IBP with 2-bit coefficients
        sage: bp                                                                        # needs sage.symbolic
        <IBP: ((-3, -1, 1, -1, -3, -1) + [0 .. 2)) * 2^-3 over [1 .. 5/4]; usign 1;
              level 2; slope_err 3.141592653589794?>
        sage: bp.variations()                                                           # needs sage.symbolic
        (3, 3)
    """

    def __init__(self, Vector_integer_dense coeffs, Rational lower, Rational upper, int lsign, int usign, int error, int scale_log2, int level, RealIntervalFieldElement slope_err):
        """
        Initialize an interval_bernstein_polynomial_integer.

        INPUT:

        - ``coeffs`` -- a coefficient vector for a polynomial in Bernstein form
        - ``lower`` -- the lower bound of the region over which the Bernstein basis is defined
        - ``upper`` -- the upper bound of the region over which the Bernstein basis is defined
        - ``lsign`` -- the sign of the polynomial at lower, if known
        - ``usign`` -- the sign of the polynomial at upper, if known
        - ``error`` -- the maximum error in the Bernstein coefficients
        - ``scale_log2`` -- the log2 of the scaling factor for the Bernstein coefficients
        - ``level`` -- the number of times we have performed degree reduction to get this polynomial
        - ``slope_err`` -- the maximum extra error in the derivative of this polynomial from degree reduction

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = interval_bernstein_polynomial_integer(vector(ZZ, [50, -30, -10]), -3/7, 4/7, 0, -1, 17, 3, 2, RIF(10^-30))
            sage: print(bp)
            degree 2 IBP with 6-bit coefficients
            sage: bp
            <IBP: ((50, -30, -10) + [0 .. 17)) * 2^3 over [-3/7 .. 4/7]; usign -1; level 2; slope_err 1.0000000000000000?e-30>
        """
        assert len(coeffs) > 0
        self.coeffs = coeffs
        self.lower = lower
        self.upper = upper
        self.lsign = lsign
        if self.lsign == 0:
            if mpz_sgn(coeffs._entries[0]) > 0:
                self.lsign = 1
            if mpz_cmp_si(coeffs._entries[0], -error) <= 0:
                self.lsign = -1
        self.usign = usign
        cdef int n = len(coeffs)
        if self.usign == 0:
            if mpz_sgn(coeffs._entries[n-1]) > 0:
                self.usign = 1
            if mpz_cmp_si(coeffs._entries[n-1], -error) <= 0:
                self.usign = -1
        self.error = error
        self.scale_log2 = scale_log2
        self.level = level
        self.slope_err = slope_err
        self._set_bitsize()
        self._count_variations()
        self.lft = None

    def __repr__(self):
        """
        Reveal all the internals of this interval Bernstein polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([-11, 22, -33], upper=1/9, error=20, lsign=1)
            sage: repr(bp)
            '<IBP: (-11, 22, -33) + [0 .. 20) over [0 .. 1/9]; lsign 1>'
        """
        base = "%s + [0 .. %s)" % (self.coeffs, self.error)
        if self.scale_log2 != 0:
            base = "(%s) * 2^%d" % (base, self.scale_log2)
        s = "<IBP: %s" % base
        if self.lower != 0 or self.upper != 1:
            s += " over [%s .. %s]" % (self.lower, self.upper)
        if (mpz_sgn(self.coeffs._entries[0]) <= 0 and mpz_cmp_si(self.coeffs._entries[0], -self.error) > 0) and self.lsign != 0:
            s += "; lsign %d" % self.lsign
        cdef int n = len(self.coeffs)
        if (mpz_sgn(self.coeffs._entries[n-1]) <= 0 and mpz_cmp_si(self.coeffs._entries[n-1], -self.error) > 0) and self.usign != 0:
            s += "; usign %d" % self.usign
        if self.level != 0:
            s += "; level %d" % self.level
        if not (self.slope_err == 0):
            s += "; slope_err %s" % self.slope_err
        return s + ">"

    def __str__(self):
        """
        Return a short summary of this interval Bernstein polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([-11, 22, -33], upper=1/9, error=20, lsign=1)
            sage: print(bp)
            degree 2 IBP with 6-bit coefficients
            sage: str(bp)
            'degree 2 IBP with 6-bit coefficients'
        """
        return "degree %d IBP with %d-bit coefficients" % (len(self.coeffs) - 1, self.bitsize)

    def _type_code(self):
        """
        Classifies this as either an integer representation ('i') or a
        floating-point representation ('f').
        """
        return 'i'

    cdef void _set_bitsize(self) noexcept:
        """
        A private function that computes the maximum coefficient size
        of this Bernstein polynomial (in bits).

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: print(mk_ibpi([2^12345]))
            degree 0 IBP with 12346-bit coefficients
            sage: print(mk_ibpi([2^12345 - 1]))
            degree 0 IBP with 12345-bit coefficients
        """
        self.bitsize = max_bitsize_intvec(self.coeffs)

    cdef void _count_variations(self) noexcept:
        """
        A private function that counts the number of sign variations in
        this Bernstein polynomial.  Since the coefficients are represented
        with intervals, not exactly, we cannot necessarily compute the exact
        number of sign variations; instead, we compute lower and upper
        bounds on this number.

        TESTS::

            sage: from sage.rings.polynomial.real_roots import *
            sage: mk_ibpi([-1, -2, -3], error=1).variations()
            (0, 0)
            sage: mk_ibpi([-1, -2, -3], error=2).variations()
            (0, 2)
            sage: mk_ibpi([-1, -2, -3], error=2, lsign=-1).variations()
            (0, 2)
            sage: mk_ibpi([-1, -2, -3], error=2, lsign=0).variations()
            (0, 2)
            sage: mk_ibpi([-1, -2, -3], error=2, lsign=1).variations()
            (1, 1)
            sage: mk_ibpi([-1, -2, -3], error=3).variations()
            (0, 2)
            sage: mk_ibpi([-1, -2, -3], error=3, lsign=-1).variations()
            (0, 2)
            sage: mk_ibpi([-1, -2, -3], error=3, lsign=1).variations()
            (1, 1)
            sage: mk_ibpi([-1, -2, -3], error=4).variations()
            (0, 2)
            sage: mk_ibpi([-1, -2, -3], error=4, lsign=-1, usign=1).variations()
            (1, 1)
            sage: mk_ibpi([-1, -2, -3], error=4, lsign=1, usign=-1).variations()
            (1, 1)
        """
        cdef Vector_integer_dense c = self.coeffs

        cdef int count_maybe_pos, count_maybe_neg
        cdef int sign
        cdef int count_definite = 0

        cdef int n = len(c)

        cdef int new_count_maybe_pos, new_count_maybe_neg

        cdef int lower_sgn, upper_sgn

        cdef int i

        if self.lsign > 0:
            count_maybe_pos = 0
            count_maybe_neg = -1
            sign = 1
        elif self.lsign < 0:
            count_maybe_pos = -1
            count_maybe_neg = 0
            sign = -1
        else:
            count_maybe_pos = 0
            count_maybe_neg = 0
            sign = 0

        for i in range(1, n):
            lower_sgn = mpz_sgn(c._entries[i])
            upper_sgn = mpz_cmp_si(c._entries[i], -self.error)
            new_count_maybe_pos = count_maybe_pos
            new_count_maybe_neg = count_maybe_neg
            if lower_sgn > 0:
                if sign < 0:
                    count_definite = count_definite + 1
                sign = 1
                new_count_maybe_neg = -1
            if upper_sgn <= 0:
                if sign > 0:
                    count_definite = count_definite + 1
                sign = -1
                new_count_maybe_pos = -1
            if upper_sgn >= 0 and count_maybe_neg + 1 > new_count_maybe_pos:
                new_count_maybe_pos = count_maybe_neg + 1
            if lower_sgn < 0 and count_maybe_pos + 1 > new_count_maybe_neg:
                new_count_maybe_neg = count_maybe_pos + 1

            count_maybe_pos = new_count_maybe_pos
            count_maybe_neg = new_count_maybe_neg

        if self.usign > 0 and sign < 0:
            count_definite = count_definite + 1
        if self.usign < 0 and sign > 0:
            count_definite = count_definite + 1
        self.min_variations = count_definite

        if self.usign > 0:
            self.max_variations = count_maybe_pos
        elif self.usign < 0:
            self.max_variations = count_maybe_neg
        else:
            self.max_variations = max(count_maybe_pos, count_maybe_neg)

    cdef int degree(self) noexcept:
        """
        Return the (formal) degree of this polynomial.
        """
        return len(self.coeffs) - 1

    def de_casteljau(self, context ctx, mid, msign=0):
        """
        Uses de Casteljau's algorithm to compute the representation
        of this polynomial in a Bernstein basis over new regions.

        INPUT:

        - ``mid`` -- where to split the Bernstein basis region; 0 < mid < 1
        - ``msign`` -- default 0 (unknown); the sign of this polynomial at mid

        OUTPUT:

        - ``bp1``, ``bp2`` -- the new interval Bernstein polynomials
        - ``ok`` -- boolean; ``True`` if the sign of the original polynomial at
          mid is known

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([50, 20, -90, -70, 200], error=5)
            sage: ctx = mk_context()
            sage: bp1, bp2, ok = bp.de_casteljau(ctx, 1/2)
            sage: bp1
            <IBP: (50, 35, 0, -29, -31) + [0 .. 6) over [0 .. 1/2]>
            sage: bp2
            <IBP: (-31, -33, -8, 65, 200) + [0 .. 6) over [1/2 .. 1]>
            sage: bp1, bp2, ok = bp.de_casteljau(ctx, 2/3)
            sage: bp1
            <IBP: (50, 30, -26, -55, -13) + [0 .. 6) over [0 .. 2/3]>
            sage: bp2
            <IBP: (-13, 8, 47, 110, 200) + [0 .. 6) over [2/3 .. 1]>
            sage: bp1, bp2, ok = bp.de_casteljau(ctx, 7/39)
            sage: bp1
            <IBP: (50, 44, 36, 27, 17) + [0 .. 6) over [0 .. 7/39]>
            sage: bp2
            <IBP: (17, -26, -75, -22, 200) + [0 .. 6) over [7/39 .. 1]>
        """
        (c1_, c2_, err_inc) = de_casteljau_intvec(self.coeffs, self.bitsize, mid, ctx.wordsize == 32)
        cdef Vector_integer_dense c1 = c1_
        cdef Vector_integer_dense c2 = c2_

        cdef int new_err = self.error + err_inc

        cdef int sign = 0

        if mpz_sgn(c2._entries[0]) > 0:
            sign = 1
        if mpz_cmp_si(c2._entries[0], -new_err) <= 0:
            sign = -1

        if msign == 0:
            msign = sign
        elif sign != 0:
            assert msign == sign

        cdef Rational absolute_mid = self.lower + mid * (self.upper - self.lower)

        cdef interval_bernstein_polynomial_integer bp1, bp2
        bp1 = interval_bernstein_polynomial_integer(c1, self.lower, absolute_mid, self.lsign, msign, new_err, self.scale_log2, self.level, self.slope_err)
        bp2 = interval_bernstein_polynomial_integer(c2, absolute_mid, self.upper, msign, self.usign, new_err, self.scale_log2, self.level, self.slope_err)

        if self.lft is not None:
            (a, b, c, d) = self.lft
            bp1.lft = (a * mid, b, c * mid, d)
            bp2.lft = (a * (1-mid), b + a*mid, c * (1-mid), d + c*mid)

        if msign != 0:
            self.update_variations(bp1, bp2)

        return (bp1, bp2, msign != 0)

    def as_float(self):
        """
        Compute an interval_bernstein_polynomial_float which contains
        (or bounds) all the polynomials this interval polynomial
        contains (or bounds).

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([50, 20, -90, -70, 200], error=5)
            sage: print(bp.as_float())
            degree 4 IBP with floating-point coefficients
            sage: bp.as_float()
            <IBP: ((0.1953125, 0.078125, -0.3515625, -0.2734375, 0.78125) + [-1.1275702593849246e-16 .. 0.01953125000000017]) * 2^8>
        """
        (fcoeffs, neg_err, pos_err, scale_log2_delta) = intvec_to_doublevec(self.coeffs, self.error)
        cdef interval_bernstein_polynomial_float fbp = interval_bernstein_polynomial_float(fcoeffs, self.lower, self.upper, self.lsign, self.usign, neg_err, pos_err, self.scale_log2 + scale_log2_delta, self.level, self.slope_err)
        fbp.min_variations = self.min_variations
        fbp.max_variations = self.max_variations
        return fbp

    def get_msb_bit(self):
        """
        Return an approximation of the log2 of the maximum of the
        absolute values of the coefficients, as an integer.
        """
        return self.scale_log2 + self.bitsize

    def down_degree(self, context ctx, max_err, exp_err_shift):
        """
        Compute an interval_bernstein_polynomial_integer which bounds
        all the polynomials this interval polynomial bounds, but is
        of lesser degree.

        During the computation, we find an "expected error"
        expected_err, which is the error inherent in our approach
        (this depends on the degrees involved, and is proportional
        to the error of the current polynomial).

        We require that the error of the new interval polynomial
        be bounded both by max_err, and by expected_err << exp_err_shift.
        If we find such a polynomial p, then we return a pair of p and some
        debugging/logging information.  Otherwise, we return the pair
        (None, None).

        If the resulting polynomial would have error more than 2^17,
        then it is downscaled before returning.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([0, 100, 400, 903], error=2)
            sage: ctx = mk_context()
            sage: bp
            <IBP: (0, 100, 400, 903) + [0 .. 2)>
            sage: dbp, _ = bp.down_degree(ctx, 10, 32)
            sage: dbp
            <IBP: (-1, 148, 901) + [0 .. 4); level 1; slope_err 0.?e2>
        """

        # return (None, None)

        degree = self.degree()

        try:
            (next, downmat, expected_err, (bdi, den)) = dr_cache[degree]
        except KeyError:
            precompute_degree_reduction_cache(degree)
            (next, downmat, expected_err, (bdi, den)) = dr_cache[degree]

        expected_err = expected_err * self.error
        if expected_err > max_err:
            return (None, ('$', self.scale_log2 + bitsize(expected_err)))

        if (expected_err << exp_err_shift) < max_err:
            max_err = (expected_err << exp_err_shift)

        v0 = dprod_imatrow_vec(bdi, self.coeffs, 0) // den
#         v0 = dprod_matrow_vec(downmat, self.coeffs, 0)
        v0_err = abs(v0 - self.coeffs[0])
        if v0_err > max_err:
            return (None, ('>', self.scale_log2 + bitsize(v0_err)))
        vn = dprod_imatrow_vec(bdi, self.coeffs, next) // den
#         vn = dprod_matrow_vec(downmat, self.coeffs, next)
        vn_err = abs(vn - self.coeffs[degree])
        if vn_err > max_err:
            return (None, ('>', self.scale_log2 + bitsize(vn_err)))
        ribp = (ZZ**(next+1))(0)
        ribp[0] = v0
        ribp[next] = vn
        for i in range(1, next):
            ribp[i] = dprod_imatrow_vec(bdi, self.coeffs, i) // den

        before = time.time()
        (eribp, err_inc) = bernstein_expand(ribp, self.degree())
        after = time.time()
        (base_max, base_min) = min_max_delta_intvec(self.coeffs, eribp)
        err_max = base_max + self.error
        err_min = base_min - err_inc
        for i in range(len(ribp)):
            ribp[i] = ribp[i] + err_min
        err = err_max - err_min

        ctx.be_log_append((after - before, self.scale_log2, max_bitsize_intvec(ribp), err < max_err, err_max, err_min, err_inc, self.degree(), next, self.variations()))
        if err < max_err:
            rng = RIF(-2*err, 2*err)
            rng = rng * self.degree()
            width = self.region_width()
            rng = rng / width
            if self.scale_log2 >= 0:
                rng = rng << self.scale_log2
            else:
                rng = rng >> (-self.scale_log2)
            # slope_err could be smaller by actually computing
            # the derivative of the error polynomial
            slope_err = rng
            lsb = self.scale_log2
            if err <= expected_err:
                indicator = '<'
            else:
                indicator = '='
            if bitsize(err) > 17:
                shift = bitsize(err) - 15
                for i in range(len(ribp)):
                    ribp[i] = ribp[i] >> shift
                max_err = max_err >> shift
                err = -((-err) >> shift)
                lsb = lsb + shift
            # warning: lsign and usign might have changed, invalidate them to 0
            ibp = interval_bernstein_polynomial_integer(ribp, self.lower, self.upper, 0, 0, err, lsb, self.level+1, self.slope_err + slope_err)

            return (ibp, (indicator, lsb + bitsize(err)))
        else:
            return (None, ('=', self.scale_log2 + bitsize(err)))

    def down_degree_iter(self, context ctx, max_scale):
        """
        Compute a degree-reduced version of this interval polynomial, by
        iterating down_degree.

        We stop when degree reduction would give a polynomial which is
        too inaccurate, meaning that either we think the current polynomial
        may have more roots in its region than the degree of the
        reduced polynomial, or that the least significant accurate bit
        in the result (on the absolute scale) would be larger than
        1 << max_scale.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([0, 100, 400, 903, 1600, 2500], error=2)
            sage: ctx = mk_context()
            sage: bp
            <IBP: (0, 100, 400, 903, 1600, 2500) + [0 .. 2)>
            sage: rbp = bp.down_degree_iter(ctx, 6)
            sage: rbp
            <IBP: (-4, 249, 2497) + [0 .. 9); level 2; slope_err 0.?e3>
        """

        cdef interval_bernstein_polynomial bp = self

        while True:
            next = degree_reduction_next_size(bp.degree())
            if next is None:
                return bp
            if bp.variations()[0] > next:
                return bp
            (rbp, err_info) = bp.down_degree(ctx, Integer(1) << (max_scale - bp.scale_log2), 32)
            ctx.dc_log_append(('down_degree', rbp is not None, err_info))
            if rbp is None:
                # global dd_count_no
                # dd_count_no += 1
                return bp
            else:
                # global dd_count_yes
                # dd_count_yes += 1
                bp = rbp

    def downscale(self, bits):
        """
        Compute an interval_bernstein_polynomial_integer which
        contains (or bounds) all the polynomials this interval
        polynomial contains (or bounds), but uses
        "bits" fewer bits.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([0, 100, 400, 903], error=2)
            sage: bp.downscale(5)
            <IBP: ((0, 3, 12, 28) + [0 .. 1)) * 2^5>
        """
        p = self.coeffs.__copy__()
        for i in range(len(p)):
            p[i] = p[i] >> bits
        return interval_bernstein_polynomial_integer(p, self.lower, self.upper, self.lsign, self.usign, -((-self.error) >> bits), self.scale_log2 + bits, self.level, self.slope_err)

    def slope_range(self):
        """
        Compute a bound on the derivative of this polynomial, over its region.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpi([0, 100, 400, 903], error=2)
            sage: bp.slope_range().str(style='brackets')
            '[294.00000000000000 .. 1515.0000000000000]'
        """
        width = self.region_width()
        (min_diff, max_diff) = min_max_diff_intvec(self.coeffs)
        rng = RIF(min_diff - self.error, max_diff + self.error)
        rng = rng * (len(self.coeffs) - 1)
        rng = rng / width
        if self.scale_log2 >= 0:
            rng = rng << self.scale_log2
        else:
            rng = rng >> (-self.scale_log2)
        return rng


def mk_ibpi(coeffs, lower=0, upper=1, lsign=0, usign=0, error=1, scale_log2=0,
            level=0, slope_err=RIF(0)):
    """
    A simple wrapper for creating interval_bernstein_polynomial_integer
    objects with coercions, defaults, etc.

    For use in doctests.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: print(mk_ibpi([50, 20, -90, -70, 200], error=5))
        degree 4 IBP with 8-bit coefficients
    """
    return interval_bernstein_polynomial_integer(vector(ZZ, coeffs), QQ(lower), QQ(upper), lsign, usign, error, scale_log2, level, slope_err)


def de_casteljau_intvec(Vector_integer_dense c, int c_bitsize, Rational x, int use_ints):
    """
    Given a polynomial in Bernstein form with integer coefficients
    over the region [0 .. 1], and a split point x, use de Casteljau's
    algorithm to give polynomials in Bernstein form over [0 .. x] and
    [x .. 1].

    This function will work for an arbitrary rational split point x, as
    long as 0 < x < 1; but it has specialized code paths that make
    some values of x faster than others.  If x == a/(a + b),
    there are special efficient cases for a==1, b==1, a+b fits in a machine
    word, a+b is a power of 2, a fits in a machine word, b fits in
    a machine word.  The most efficient case is x==1/2.

    Given split points x == a/(a + b) and y == c/(c + d), where
    min(a, b) and min(c, d) fit in the same number of machine words
    and a+b and c+d are both powers of two, then x and y should be
    equally fast split points.

    If use_ints is nonzero, then instead of checking whether numerators
    and denominators fit in machine words, we check whether they fit in
    ints (32 bits, even on 64-bit machines).  This slows things down, but
    allows for identical results across machines.

    INPUT:

    - ``c`` -- vector of coefficients of polynomial in Bernstein form
    - ``c_bitsize`` -- approximate size of coefficients in c (in bits)
    - ``x`` -- rational splitting point; 0 < x < 1

    OUTPUT:

    - ``c1`` -- coefficients of polynomial over range [0 .. x]
    - ``c2`` -- coefficients of polynomial over range [x .. 1]
    - ``err_inc`` -- amount by which error intervals widened

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: c = vector(ZZ, [1048576, 0, 0, 0, 0, 0])
        sage: de_casteljau_intvec(c, 20, 1/2, 1)
        ((1048576, 524288, 262144, 131072, 65536, 32768), (32768, 0, 0, 0, 0, 0), 1)
        sage: de_casteljau_intvec(c, 20, 1/3, 1)
        ((1048576, 699050, 466033, 310689, 207126, 138084), (138084, 0, 0, 0, 0, 0), 1)
        sage: de_casteljau_intvec(c, 20, 7/22, 1)
        ((1048576, 714938, 487457, 332357, 226607, 154505), (154505, 0, 0, 0, 0, 0), 1)
    """
    vs = c.parent()

    cdef Vector_integer_dense c1, c2, den_powers

    c1 = Vector_integer_dense(vs, 0)
    c2 = c.__copy__()

    cdef int n = len(c)

    cdef int i, j

    cdef mpz_t num, den, diff, tmp, tmp2

    cdef unsigned long num_ui, den_ui, diff_ui
    cdef int num_fits_ui, den_fits_ui, diff_fits_ui
    cdef int den_is_pow2, den_log2
    cdef int num_less_diff

    cdef int max_den_power

    # We want to compute (diff*a + num*b)/den.  This costs
    # 2*mul, 1*div, 1*add/sub.
    # We see below a way to trade a multiplication for a subtraction,
    # for a cost of 1*mul, 1*div, 2*add/sub.
    # Another possibility is to postpone divisions, for a cost of
    # 2*mul, 1/6*div, 1*add/sub.  Clearly, this is worthwhile if
    # 5/6*div + 1*add/sub > 1*mul.  This is very likely true whenever
    # the denominator is not a power of 2 (so that you have to do
    # real divisions, instead of shifts), and might be true in other cases.

    # If one of the multiplications is by 1, then the costs become:
    # straightforward: 1*mul, 1*div, 1*add/sub
    # trade mul->sub: 1*div, 2*add/sub
    # postpone division: 1*mul, 1/6*div, 1*add/sub
    # for the same result.

    # If both multiplications are by 1, then the costs are:
    # straightforward: 1*div, 1*add/sub
    # trade mul->sub: 1*div, 2*add/sub
    # postpone division: 1/6*div, 1*add/sub
    # with postpone division being the obvious winner.

    # For now, we use "postpone division" whenever the denominator fits
    # in an unsigned long.  This is not going to be drastically bad,
    # even when dividing by a power of two.  For full generality,
    # it would also be important to use "postpone division" for
    # non-power-of-two denominators that do not fit in an unsigned long;
    # however, the current calling code essentially never does that,
    # so we'll stick with simpler code here.

    mpz_init(num)
    mpz_init(den)
    mpz_init(diff)
    mpz_init(tmp)
    mpz_init(tmp2)

    mpq_get_num(num, x.value)
    mpq_get_den(den, x.value)
    mpz_sub(diff, den, num)

    if use_ints:
        den_fits_ui = mpz_fits_uint_p(den)
        num_fits_ui = mpz_fits_uint_p(num)
        diff_fits_ui = mpz_fits_uint_p(diff)
    else:
        den_fits_ui = mpz_fits_ulong_p(den)
        num_fits_ui = mpz_fits_ulong_p(num)
        diff_fits_ui = mpz_fits_ulong_p(diff)

    if den_fits_ui:
        den_ui = mpz_get_ui(den)
    if num_fits_ui:
        num_ui = mpz_get_ui(num)
    if diff_fits_ui:
        diff_ui = mpz_get_ui(diff)

    mpz_sub_ui(tmp, den, 1)
    mpz_and(tmp, tmp, den)
    den_is_pow2 = (mpz_sgn(tmp) == 0)
    if den_is_pow2:
        den_log2 = mpz_sizeinbase(den, 2) - 1

    cdef int max_den_bits = c_bitsize / 2
    if max_den_bits < 100:
        max_den_bits = 100
# These settings are slower than the above on laguerre(1000), but that's
# the only experiment I've done so far... more testing is needed.
#     cdef int max_den_bits = 3 * c_bitsize / 2
#     if max_den_bits < 100: max_den_bits = 300

    cdef int cur_den_steps = 0

    cdef int ndivs = 0

    if den_fits_ui:
        den_powers = FreeModule(ZZ, len(c)+1)(0)
        mpz_set_ui(den_powers._entries[0], 1)
        max_den_power = 1
        for i in range(1, n + 1):
            mpz_mul_ui(den_powers._entries[i], den_powers._entries[i-1], den_ui)
            if mpz_sizeinbase(den_powers._entries[i], 2) < max_den_bits:
                max_den_power = i
            else:
                break

        for i in range(n):
            mpz_set(c1._entries[i], c2._entries[0])
            if den_ui == 2:
                # x == 1/2
                for j in range(n-i-1):
                    mpz_add(c2._entries[j], c2._entries[j], c2._entries[j+1])
            else:
                for j in range(n-i-1):
                    if diff_ui != 1:
                        mpz_mul_ui(c2._entries[j], c2._entries[j], diff_ui)
                    if num_ui == 1:
                        mpz_add(c2._entries[j], c2._entries[j], c2._entries[j+1])
                    else:
                        mpz_addmul_ui(c2._entries[j], c2._entries[j+1], num_ui)
            cur_den_steps = cur_den_steps + 1
            if cur_den_steps == max_den_power or i == n-1:
                if den_is_pow2:
                    for j from 1 <= j < cur_den_steps:
                        mpz_fdiv_q_2exp(c1._entries[i-cur_den_steps+j+1],
                                        c1._entries[i-cur_den_steps+j+1],
                                        j*den_log2)
                        mpz_fdiv_q_2exp(c2._entries[n-i+cur_den_steps-j-2],
                                        c2._entries[n-i+cur_den_steps-j-2],
                                        j*den_log2)
                    for j from 0 <= j < n-i-1:
                        mpz_fdiv_q_2exp(c2._entries[j], c2._entries[j],
                                        cur_den_steps*den_log2)
                else:
                    for j from 1 <= j < cur_den_steps:
                        mpz_fdiv_q(c1._entries[i-cur_den_steps+j+1],
                                   c1._entries[i-cur_den_steps+j+1],
                                   den_powers._entries[j])
                        mpz_fdiv_q(c2._entries[n-i+cur_den_steps-j-2],
                                   c2._entries[n-i+cur_den_steps-j-2],
                                   den_powers._entries[j])
                    for j from 0 <= j < n-i-1:
                        mpz_fdiv_q(c2._entries[j], c2._entries[j],
                                   den_powers._entries[cur_den_steps])
                cur_den_steps = 0
                ndivs = ndivs+1
    else:
        # We want to compute (diff*a + num*b)/den, where
        # a == c2._entries[j] and b == c2._entries[j+1].
        # The result goes in c2._entries[j].  Since den == diff + num,
        # this is equal to a + num*(b-a)/den or diff*(a-b)/den + b.
        # If num<diff, we compute the former, otherwise the latter.

        num_less_diff = (mpz_cmp(num, diff) < 0)

        ndivs = n-1

        for i in range(n):
            mpz_set(c1._entries[i], c2._entries[0])
            for j in range(n-i-1):
                if num_less_diff:
                    mpz_sub(tmp, c2._entries[j+1], c2._entries[j])
                    if num_fits_ui and num_ui == 1:
                        if den_is_pow2:
                            mpz_fdiv_q_2exp(tmp2, tmp, den_log2)
                        elif den_fits_ui:
                            mpz_fdiv_q_ui(tmp2, tmp, den_ui)
                        else:
                            mpz_fdiv_q(tmp2, tmp, den)
                        mpz_add(c2._entries[j], c2._entries[j], tmp2)
                    else:
                        if num_fits_ui:
                            mpz_mul_ui(tmp2, tmp, num_ui)
                        else:
                            mpz_mul(tmp2, tmp, num)
                        if den_is_pow2:
                            mpz_fdiv_q_2exp(tmp, tmp2, den_log2)
                        elif den_fits_ui:
                            mpz_fdiv_q_ui(tmp, tmp2, den_ui)
                        else:
                            mpz_fdiv_q(tmp, tmp2, den)
                        mpz_add(c2._entries[j], c2._entries[j], tmp)
                else:
                    mpz_sub(c2._entries[j], c2._entries[j], c2._entries[j+1])
                    if diff_fits_ui:
                        if diff_ui == 1:
                            mpz_set(tmp, c2._entries[j])
                        else:
                            mpz_mul_ui(tmp, c2._entries[j], diff_ui)
                    else:
                        mpz_mul(tmp, c2._entries[j], diff)
                    if den_is_pow2:
                        mpz_fdiv_q_2exp(c2._entries[j], tmp, den_log2)
                    elif den_fits_ui:
                        mpz_fdiv_q_ui(c2._entries[j], tmp, den_ui)
                    else:
                        mpz_fdiv_q(c2._entries[j], tmp, den)
                    mpz_add(c2._entries[j], c2._entries[j], c2._entries[j+1])

    mpz_clear(num)
    mpz_clear(den)
    mpz_clear(diff)
    mpz_clear(tmp)
    mpz_clear(tmp2)

    return (c1, c2, ndivs)

# An ULP is a "unit in the last place"; it is the (varying) unit for
# how much adjacent floating-point numbers differ from each other.
# A half-ULP is half this amount; it is the maximum rounding error
# in the basic operations (+, -, *, /) in a correctly-operating IEEE
# floating-point unit.
# (Note that by default, the x86 does not use IEEE double precision;
# instead, it uses extra precision, which can (counterintuitively)
# actually give worse double-precision results in some rare cases.
# To avoid this, we change the x86 floating-point unit into true
# double-precision mode in places where it matters; that is, in
# functions that add, subtract, multiply, or divide floating-point numbers.
# Functions whose floating-point operations are limited to negation
# and comparison do not require special treatment, since those operations
# give the same results in double-precision or extended precision.)


# This is the value of a half-ULP for numbers in the range 0.5 <= x < 1.
# This is actually slightly more than a half-ULP because of possible
# double-rounding on x86 PCs.
cdef double half_ulp = ldexp(1.0 * 65/64, -54)


def intvec_to_doublevec(Vector_integer_dense b, long err):
    r"""
    Given a vector of integers A = [a1, ..., an], and an integer
    error bound E, returns a vector of floating-point numbers
    B = [b1, ..., bn], lower and upper error bounds F1 and F2, and
    a scaling factor d, such that

    .. MATH::

       (bk + F1) * 2^d \le ak

    and

    .. MATH::

        ak + E \le (bk + F2) * 2^d

    If bj is the element of B with largest absolute value, then
    0.5 <= abs(bj) < 1.0 .

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: intvec_to_doublevec(vector(ZZ, [1, 2, 3, 4, 5]), 3)
        ((0.125, 0.25, 0.375, 0.5, 0.625), -1.1275702593849246e-16, 0.37500000000000017, 3)
    """
    vs = FreeModule(RDF, len(b))

    cdef Vector_real_double_dense db = vs(0)
    cdef numpy.ndarray[double, ndim=1] dbv = db._vector_numpy

    cdef long max_exp = -10000
    cdef long cur_exp

    cdef int i

    for i in range(len(b)):
        mpz_get_d_2exp(&cur_exp, b._entries[i])
        if cur_exp > max_exp:
            max_exp = cur_exp

    cdef int delta = -max_exp
    cdef double d
    cdef int new_exp

    for i in range(len(b)):
        d = mpz_get_d_2exp(&cur_exp, b._entries[i])
        # 0.5 <= d < 1; b._entries[i] ~= d*2^cur_exp
        new_exp = cur_exp + delta
        if new_exp >= -100:
            dbv[i] = ldexp(d, new_exp)
#            db[i] = ldexp(d, new_exp)

    # The true value can be off by an ulp because mpz_get_d_2exp truncates.
    # (If we created a version of mpz_get_d_2exp that rounded instead,
    # we could do a little better.)
    # The third half_ulp in the positive direction is to compensate for
    # possible bad rounding when adding ldexp(err, delta).

    cdef double low_err = -2*half_ulp
    cdef double high_err = 3*half_ulp + ldexp(err, delta)
    # fpu_fix_end(&cwf)
    return (db, low_err, high_err, -delta)


cdef class interval_bernstein_polynomial_float(interval_bernstein_polynomial):
    """
    This is the subclass of :class:`interval_bernstein_polynomial` where
    polynomial coefficients are represented using floating-point numbers.

    In the floating-point representation, each coefficient is represented
    as an IEEE double-precision float A, and the (shared) lower and
    upper interval widths E1 and E2.  These represent the coefficients
    (A+E1)*2^n <= c <= (A+E2)*2^n.

    Note that we always have E1 <= 0 <= E2.  Also, each floating-point
    coefficient has absolute value less than one.

    (Note that :func:`mk_ibpf` is a simple helper function for creating
    elements of :class:`interval_bernstein_polynomial_float` in doctests.)

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: bp = mk_ibpf([0.1, 0.2, 0.3], pos_err=0.5); print(bp)
        degree 2 IBP with floating-point coefficients
        sage: bp
        <IBP: (0.1, 0.2, 0.3) + [0.0 .. 0.5]>
        sage: bp.variations()
        (0, 0)
        sage: bp = mk_ibpf([-0.3, -0.1, 0.1, -0.1, -0.3, -0.1],                         # needs sage.symbolic
        ....:              lower=1, upper=5/4, usign=1, pos_err=0.2,
        ....:              scale_log2=-3, level=2, slope_err=RIF(pi)); print(bp)
        degree 5 IBP with floating-point coefficients
        sage: bp                                                                        # needs sage.symbolic
        <IBP: ((-0.3, -0.1, 0.1, -0.1, -0.3, -0.1) + [0.0 .. 0.2]) * 2^-3
              over [1 .. 5/4]; usign 1; level 2; slope_err 3.141592653589794?>
        sage: bp.variations()                                                           # needs sage.symbolic
        (3, 3)
    """

    def __init__(self, Vector_real_double_dense coeffs, Rational lower, Rational upper, int lsign, int usign, double neg_err, double pos_err, int scale_log2, int level, RealIntervalFieldElement slope_err):
        """
        Initialize an interval_bernstein_polynomial_float.

        INPUT:

        - ``coeffs`` -- a coefficient vector for a polynomial in Bernstein form
          (all coefficients must have absolute value less than one)
        - ``lower`` -- the lower bound of the region over which the Bernstein basis
          is defined
        - ``upper`` -- the upper bound of the region over which the Bernstein basis
          is defined
        - ``lsign`` -- the sign of the polynomial at lower, if known
        - ``usign`` -- the sign of the polynomial at upper, if known
        - ``neg_err`` -- the minimum error in the Bernstein coefficients
        - ``pos_err`` -- the maximum error in the Bernstein coefficients (so a
          Bernstein coefficient x really represents the range [neg_err+x ..
          pos_err+x]
        - ``scale_log2`` -- the log2 of the scaling factor for the Bernstein
          coefficients
        - ``level`` -- the number of times we have performed degree reduction to
          get this polynomial
        - ``slope_err`` -- the maximum extra error in the derivative of this
          polynomial from degree reduction

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = interval_bernstein_polynomial_float(vector(RDF, [0.50, -0.30, -0.10]), -3/7, 4/7, 0, -1, -0.02, 0.17, 3, 2, RIF(10^-30))
            sage: print(bp)
            degree 2 IBP with floating-point coefficients
            sage: bp
            <IBP: ((0.5, -0.3, -0.1) + [-0.02 .. 0.17]) * 2^3 over [-3/7 .. 4/7]; usign -1; level 2; slope_err 1.0000000000000000?e-30>
        """
        assert len(coeffs) > 0
        cdef numpy.ndarray[double, ndim=1] coeffs_data = coeffs._vector_numpy
        self.coeffs = coeffs
        self.lower = lower
        self.upper = upper
        self.lsign = lsign
        if self.lsign == 0:
            if coeffs_data[0] > -neg_err:
                self.lsign = 1
            if coeffs_data[0] < -pos_err:
                self.lsign = -1
        self.usign = usign
        cdef int n = len(coeffs)
        if self.usign == 0:
            if coeffs_data[n-1] > -neg_err:
                self.usign = 1
            if coeffs_data[n-1] < -pos_err:
                self.usign = -1
        self.neg_err = neg_err
        self.pos_err = pos_err
        self.scale_log2 = scale_log2
        self.level = level
        self.slope_err = slope_err
        self._count_variations()
        max_abs = max_abs_doublevec(coeffs)
        cdef int exp
        frexp(max_abs, &exp)
        self.bitsize = exp + 53
        self.lft = None

    def __repr__(self):
        """
        Reveal all the internals of this interval Bernstein polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpf([-0.11, 0.22, -0.33], upper=1/9, neg_err=-0.3, pos_err=0.1, lsign=1)
            sage: repr(bp)
            '<IBP: (-0.11, 0.22, -0.33) + [-0.3 .. 0.1] over [0 .. 1/9]>'
        """
        base = "%s + [%r .. %r]" % (self.coeffs, self.neg_err, self.pos_err)
        if self.scale_log2 != 0:
            base = "(%s) * 2^%d" % (base, self.scale_log2)
        s = "<IBP: %s" % base
        if self.lower != 0 or self.upper != 1:
            s += " over [%s .. %s]" % (self.lower, self.upper)
        if self.coeffs.get_unsafe(0) <= -self.neg_err and self.coeffs.get_unsafe(0) >= -self.pos_err and self.lsign != 0:
            s += "; lsign %d" % self.lsign
        cdef int n = len(self.coeffs)
        if self.coeffs.get_unsafe(n-1) <= -self.neg_err and self.coeffs.get_unsafe(n-1) >= -self.pos_err and self.usign != 0:
            s += "; usign %d" % self.usign
        if self.level != 0:
            s += "; level %d" % self.level
        if self.slope_err != 0:
            s += "; slope_err %s" % self.slope_err
        return s + ">"

    def __str__(self):
        """
        Return a short summary of this interval Bernstein polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpf([-0.11, 0.22, -0.33], upper=1/9, neg_err=-0.1, pos_err=0.2, lsign=1)
            sage: print(bp)
            degree 2 IBP with floating-point coefficients
            sage: str(bp)
            'degree 2 IBP with floating-point coefficients'
        """
        return "degree %d IBP with floating-point coefficients" % (len(self.coeffs) - 1)

    def _type_code(self):
        """
        Classifies this as either an integer representation ('i') or a
        floating-point representation ('f').
        """
        return 'f'

    cdef void _count_variations(self) noexcept:
        """
        A private function that counts the number of sign variations in
        this Bernstein polynomial.  Since the coefficients are represented
        with intervals, not exactly, we cannot necessarily compute the exact
        number of sign variations; instead, we compute lower and upper
        bounds on this number.
        """
        cdef numpy.ndarray[double, ndim=1] cd = self.coeffs._vector_numpy

        cdef int count_maybe_pos
        cdef int count_maybe_neg
        cdef int sign
        cdef int count_definite = 0

        cdef int n = len(self.coeffs)

        cdef int new_count_maybe_pos, new_count_maybe_neg

        cdef int i

        cdef double val

        if self.lsign > 0:
            count_maybe_pos = 0
            count_maybe_neg = -1
            sign = 1
        elif self.lsign < 0:
            count_maybe_pos = -1
            count_maybe_neg = 0
            sign = -1
        else:
            count_maybe_pos = 0
            count_maybe_neg = 0
            sign = 0

        for i in range(1, n):
            new_count_maybe_pos = count_maybe_pos
            new_count_maybe_neg = count_maybe_neg
            val = cd[i]
            if val > -self.neg_err:
                if sign < 0:
                    count_definite = count_definite + 1
                sign = 1
                new_count_maybe_neg = -1
            if val < -self.pos_err:
                if sign > 0:
                    count_definite = count_definite + 1
                sign = -1
                new_count_maybe_pos = -1
            if val > -self.pos_err and count_maybe_neg + 1 > new_count_maybe_pos:
                new_count_maybe_pos = count_maybe_neg + 1
            if val < -self.neg_err and count_maybe_pos + 1 > new_count_maybe_neg:
                new_count_maybe_neg = count_maybe_pos + 1

            count_maybe_pos = new_count_maybe_pos
            count_maybe_neg = new_count_maybe_neg

        if self.usign > 0 and sign < 0:
            count_definite = count_definite + 1
        if self.usign < 0 and sign > 0:
            count_definite = count_definite + 1
        self.min_variations = count_definite

        if self.usign > 0:
            self.max_variations = count_maybe_pos
        elif self.usign < 0:
            self.max_variations = count_maybe_neg
        else:
            self.max_variations = max(count_maybe_pos, count_maybe_neg)

    cdef int degree(self) noexcept:
        """
        Return the (formal) degree of this polynomial.
        """
        return len(self.coeffs) - 1

    def de_casteljau(self, context ctx, mid, msign=0):
        """
        Uses de Casteljau's algorithm to compute the representation
        of this polynomial in a Bernstein basis over new regions.

        INPUT:

        - ``mid`` -- where to split the Bernstein basis region; ``0 < mid < 1``
        - ``msign`` -- default 0 (unknown); the sign of this polynomial at ``mid``

        OUTPUT:

        - ``bp1``, ``bp2`` -- the new interval Bernstein polynomials
        - ``ok`` -- boolean; ``True`` if the sign of the original polynomial at
          ``mid`` is known

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: ctx = mk_context()
            sage: bp = mk_ibpf([0.5, 0.2, -0.9, -0.7, 0.99], neg_err=-0.1, pos_err=0.01)
            sage: bp1, bp2, ok = bp.de_casteljau(ctx, 1/2)
            sage: bp1
            <IBP: (0.5, 0.35, 0.0, -0.2875, -0.369375) + [-0.10000000000000023 .. 0.010000000000000226] over [0 .. 1/2]>
            sage: bp2
            <IBP: (-0.369375, -0.45125, -0.3275, 0.14500000000000002, 0.99) + [-0.10000000000000023 .. 0.010000000000000226] over [1/2 .. 1]>
            sage: bp1, bp2, ok = bp.de_casteljau(ctx, 2/3)
            sage: bp1 # rel tol 2e-16
            <IBP: (0.5, 0.30000000000000004, -0.2555555555555555, -0.5444444444444444, -0.32172839506172846) + [-0.10000000000000069 .. 0.010000000000000677] over [0 .. 2/3]>
            sage: bp2  # rel tol 3e-15
            <IBP: (-0.32172839506172846, -0.21037037037037046, 0.028888888888888797, 0.4266666666666666, 0.99) + [-0.10000000000000069 .. 0.010000000000000677] over [2/3 .. 1]>
            sage: bp1, bp2, ok = bp.de_casteljau(ctx, 7/39)
            sage: bp1  # rel tol
            <IBP: (0.5, 0.4461538461538461, 0.36653517422748183, 0.27328680523946786, 0.1765692706232836) + [-0.10000000000000069 .. 0.010000000000000677] over [0 .. 7/39]>
            sage: bp2  # rel tol
            <IBP: (0.1765692706232836, -0.26556803047927313, -0.7802038132807364, -0.3966666666666666, 0.99) + [-0.10000000000000069 .. 0.010000000000000677] over [7/39 .. 1]>
        """
        (c1_, c2_, err_inc) = de_casteljau_doublevec(self.coeffs, mid)
        cdef Vector_real_double_dense c1 = c1_
        cdef Vector_real_double_dense c2 = c2_
        cdef numpy.ndarray[double, ndim=1] c2data = c2._vector_numpy
        cdef int sign = 0

        if c2data[0] > -self.neg_err:
            sign = 1
        if c2data[0] < -self.pos_err:
            sign = -1

        if msign == 0:
            msign = sign
        elif sign != 0:
            assert msign == sign

        # As long as new_neg and new_pos have
        # magnitudes less than 0.5, these computations
        # are exact.  This will be the case for any sensible
        # usage of this class.
        cdef double new_neg = self.neg_err - half_ulp * err_inc
        cdef double new_pos = self.pos_err + half_ulp * err_inc

        if not (-0.5 <= new_neg <= 0 <= new_pos <= 0.5):
            # Give up on this computation...it's horribly inaccurate anyway.
            msign = 0

        cdef Rational absolute_mid = self.lower + mid * (self.upper - self.lower)

        cdef interval_bernstein_polynomial_float bp1, bp2
        bp1 = interval_bernstein_polynomial_float(c1, self.lower, absolute_mid, self.lsign, msign, new_neg, new_pos, self.scale_log2, self.level, self.slope_err)
        bp2 = interval_bernstein_polynomial_float(c2, absolute_mid, self.upper, msign, self.usign, new_neg, new_pos, self.scale_log2, self.level, self.slope_err)

        if self.lft is not None:
            (a, b, c, d) = self.lft
            bp1.lft = (a * mid, b, c * mid, d)
            bp2.lft = (a * (1-mid), b + a*mid, c * (1-mid), d + c*mid)

        if msign != 0:
            self.update_variations(bp1, bp2)

        return (bp1, bp2, msign != 0)

    def as_float(self):
        return self

    def get_msb_bit(self):
        """
        Return an approximation of the log2 of the maximum of the
        absolute values of the coefficients, as an integer.
        """
        return self.scale_log2 - 53 + self.bitsize

    def slope_range(self):
        """
        Compute a bound on the derivative of this polynomial, over its region.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bp = mk_ibpf([0.5, 0.2, -0.9, -0.7, 0.99], neg_err=-0.1, pos_err=0.01)
            sage: bp.slope_range().str(style='brackets')
            '[-4.8400000000000017 .. 7.2000000000000011]'
        """
        width = self.region_width()
        (min_diff, max_diff) = min_max_diff_doublevec(self.coeffs)
        err = self.pos_err - self.neg_err
        # 2 half_ulp's because subtracting two numbers with absolute values
        # in (-1 .. 1) can give a number in (-2 .. 2), and the subtraction
        # can have an error of up to half an ulp in that range, which
        # is 2 half ulps for (-1 .. 1).
        rng = RIF(min_diff - err - 2*half_ulp, max_diff + err + 2*half_ulp)
        rng = rng * (len(self.coeffs) - 1)
        rng = rng / width
        if self.scale_log2 >= 0:
            rng = rng << self.scale_log2
        else:
            rng = rng >> (-self.scale_log2)

        # fpu_fix_end(&cwf)

        return rng


def mk_ibpf(coeffs, lower=0, upper=1, lsign=0, usign=0, neg_err=0, pos_err=0,
            scale_log2=0, level=0, slope_err=RIF(0)):
    """
    A simple wrapper for creating :class:`interval_bernstein_polynomial_float`
    objects with coercions, defaults, etc.

    For use in doctests.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: print(mk_ibpf([0.5, 0.2, -0.9, -0.7, 0.99], pos_err=0.1, neg_err=-0.01))
        degree 4 IBP with floating-point coefficients
    """
    return interval_bernstein_polynomial_float(vector(RDF, coeffs), QQ(lower), QQ(upper), lsign, usign, neg_err, pos_err, scale_log2, level, slope_err)


cdef Rational QQ_1_2 = QQ((1, 2))
cdef Rational QQ_1_32 = QQ((1, 32))
cdef Rational QQ_31_32 = QQ((31, 32))
cdef Rational zero_QQ = QQ.zero()
cdef Rational one_QQ = QQ.one()
cdef Integer zero_ZZ = ZZ.zero()
cdef Integer one_ZZ = ZZ.one()
cdef Integer ZZ_2_31 = ZZ(2) ** 31
cdef Integer ZZ_2_32 = ZZ(2) ** 32
cdef RealIntervalFieldElement zero_RIF = RIF(0)


def de_casteljau_doublevec(Vector_real_double_dense c, Rational x):
    """
    Given a polynomial in Bernstein form with floating-point coefficients
    over the region [0 .. 1], and a split point x, use de Casteljau's
    algorithm to give polynomials in Bernstein form over [0 .. x] and
    [x .. 1].

    This function will work for an arbitrary rational split point x, as
    long as 0 < x < 1; but it has a specialized code path for x==1/2.

    INPUT:

    - ``c`` -- vector of coefficients of polynomial in Bernstein form
    - ``x`` -- rational splitting point; 0 < x < 1

    OUTPUT:

    - ``c1`` -- coefficients of polynomial over range [0 .. x]
    - ``c2`` -- coefficients of polynomial over range [x .. 1]
    - ``err_inc`` -- number of half-ulps by which error intervals widened

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: c = vector(RDF, [0.7, 0, 0, 0, 0, 0])
        sage: de_casteljau_doublevec(c, 1/2)
        ((0.7, 0.35, 0.175, 0.0875, 0.04375, 0.021875), (0.021875, 0.0, 0.0, 0.0, 0.0, 0.0), 5)
        sage: de_casteljau_doublevec(c, 1/3)  # rel tol
        ((0.7, 0.4666666666666667, 0.31111111111111117, 0.20740740740740746, 0.13827160493827165, 0.09218106995884777), (0.09218106995884777, 0.0, 0.0, 0.0, 0.0, 0.0), 15)
        sage: de_casteljau_doublevec(c, 7/22)  # rel tol
        ((0.7, 0.4772727272727273, 0.3254132231404959, 0.22187265214124724, 0.15127680827812312, 0.10314327837144759), (0.10314327837144759, 0.0, 0.0, 0.0, 0.0, 0.0), 15)
    """
    vs = c.parent()

    cdef Vector_real_double_dense c1, c2

    c1 = Vector_real_double_dense(vs, 0)
    c2 = c.__copy__()

    cdef numpy.ndarray[double, ndim=1] c1d = c1._vector_numpy
    cdef numpy.ndarray[double, ndim=1] c2d = c2._vector_numpy

    cdef double half = 0.5

    cdef int n = len(c)
    cdef int i, j

    cdef double rx, rnx

    cdef int extra_err

    if x == QQ_1_2:
        for i in range(n):
            c1d[i] = c2d[0]
            for j in range(n-i-1):
                c2d[j] = (c2d[j] + c2d[j+1]) * half

# The following code lets us avoid O(n^2) floating-point multiplications
# in favor of O(n) calls to ldexp().  In one experiment, though, it seems
# to actually slow down the code by about 10%.  (On an Intel Core 2 Duo
# in 32-bit mode, on the test chebyt2(200).)
#             cur_den_steps = cur_den_steps + 1
#             if cur_den_steps == 1020 or i == n-1:
#                 for j from 1 <= j < cur_den_steps:
#                     c1d[i-cur_den_steps+j+1] = ldexp(c1d[i-cur_den_steps+j+1], -j)
#                     c2d[n-i+cur_den_steps-j-2] = ldexp(c2d[n-i+cur_den_steps-j-2], -j)
#                 for j from 0 <= j < n-i-1:
#                     c2d[j] = ldexp(c2d[j], -cur_den_steps)

        extra_err = n-1
    else:
        rx = x
        rnx = 1-rx
        for i from 0 <= i < n:
            c1d[i] = c2d[0]
            for j from 0 <= j < n-i-1:
                c2d[j] = (c2d[j]*rnx + c2d[j+1]*rx)

        extra_err = 3*(n-1)

    # fpu_fix_end(&cwf)

    return (c1, c2, extra_err)


def max_abs_doublevec(Vector_real_double_dense c):
    """
    Given a floating-point vector, return the maximum of the
    absolute values of its elements.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: max_abs_doublevec(vector(RDF, [0.1, -0.767, 0.3, 0.693]))
        0.767
    """
    cdef numpy.ndarray[double, ndim=1] cd = c._vector_numpy

    cdef double m = 0
    cdef double a

    for i in range(len(c)):
        a = fabs(cd[i])
        if a > m:
            m = a

    return m


def wordsize_rational(a, b, wordsize):
    """
    Given rationals a and b, select a de Casteljau split point r between
    a and b.

    An attempt is made to select an efficient split point
    (according to the criteria mentioned in the documentation
    for de_casteljau_intvec), with a bias towards split points near a.

    In full detail:

    This takes as input two rationals, a and b, such that 0<=a<=1, 0<=b<=1,
    and a!=b. This returns rational r, such that a<=r<=b or b<=r<=a.
    The denominator of r is a power of 2.  Let m be min(r, 1-r),
    nm be numerator(m), and dml be log2(denominator(m)).  The return value
    r is taken from the first of the following classes to have any
    members between a and b (except that if a <= 1/8, or 7/8 <= a, then
    class 2 is preferred to class 1).

    1. dml < wordsize
    2. bitsize(nm) <= wordsize
    3. bitsize(nm) <= 2*wordsize
    4. bitsize(nm) <= 3*wordsize

    ...

    k. bitsize(nm) <= (k-1)*wordsize

    From the first class to have members between a and b, r is chosen
    as the element of the class which is closest to a.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: wordsize_rational(1/5, 1/7, 32)
        429496729/2147483648
        sage: wordsize_rational(1/7, 1/5, 32)
        306783379/2147483648
        sage: wordsize_rational(1/5, 1/7, 64)
        1844674407370955161/9223372036854775808
        sage: wordsize_rational(1/7, 1/5, 64)
        658812288346769701/4611686018427387904
        sage: wordsize_rational(1/17, 1/19, 32)
        252645135/4294967296
        sage: wordsize_rational(1/17, 1/19, 64)
        1085102592571150095/18446744073709551616
        sage: wordsize_rational(1/1234567890, 1/1234567891, 32)
        933866427/1152921504606846976
        sage: wordsize_rational(1/1234567890, 1/1234567891, 64)
        4010925763784056541/4951760157141521099596496896
    """

    # assert 0 <= a <= 1

    swap_01 = False
    if b <= a:
        a = one_QQ-a
        b = one_QQ-b
        swap_01 = True
    sub_1 = False
    if a > QQ_1_2:
        a = a-one_QQ
        b = b-one_QQ
        sub_1 = True

    cur_size = wordsize
    # fld = RealField(cur_size, rnd='RNDU')
    fld = get_realfield_rndu(cur_size)
    cdef RealNumber rf
    while True:
        rf = fld(a)
        if cur_size == wordsize:
            assert mpfr_number_p(rf.value)
            exp = mpfr_get_exp(rf.value)
            if rf <= -(fld(-b)):
                if exp <= -3:
                    break
                # rf2 = RealField(cur_size + exp - 1, rnd='RNDU')(a)
                fld2 = get_realfield_rndu(cur_size + exp - 1)
                rf2 = fld2(a)
                if rf2 <= -(fld2(-b)):
                    rf = rf2
                    break
        if rf <= -(fld(-b)):
            break
        cur_size = cur_size + wordsize
        fld = RealField(cur_size, rnd='RNDU')

    r = rf.exact_rational()
    if sub_1:
        r = r + one_QQ
    if swap_01:
        r = one_QQ - r
    # assert 0 <= r <= 1
    return r


def relative_bounds(a, b):
    """
    INPUT:

    - ``(al, ah)`` -- pair of rationals
    - ``(bl, bh)`` -- pair of rationals

    OUTPUT:

    - ``(cl, ch)`` -- pair of rationals

    Computes the linear transformation that maps (al, ah) to (0, 1);
    then applies this transformation to (bl, bh) and returns the result.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: relative_bounds((1/7, 1/4), (1/6, 1/5))
        (2/9, 8/15)
    """
    (al, ah) = a
    (bl, bh) = b
    width = ah - al
    return ((bl - al) / width, (bh - al) / width)


cdef int bitsize(Integer a) noexcept:
    """
    Compute the number of bits required to write a given integer
    (the sign is ignored).

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: bitsize_doctest(0)
        1
        sage: bitsize_doctest(1)
        1
        sage: bitsize_doctest(2)
        2
        sage: bitsize_doctest(-2)
        2
        sage: bitsize_doctest(65535)
        16
        sage: bitsize_doctest(65536)
        17
    """
    return int(mpz_sizeinbase(a.value, 2))


def bitsize_doctest(n):
    return bitsize(n)


def degree_reduction_next_size(n):
    """
    Given n (a polynomial degree), returns either a smaller integer or ``None``.
    This defines the sequence of degrees followed by our degree reduction
    implementation.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: degree_reduction_next_size(1000)
        30
        sage: degree_reduction_next_size(20)
        15
        sage: degree_reduction_next_size(3)
        2
        sage: degree_reduction_next_size(2) is None
        True
    """

    # Virtually any descending sequence would be "correct" here; no
    # testing has been done to see if another sequence would be better.

    # Constraints: must return None for n <= 2; degree reduction to
    # degrees > 30 may be very, very slow.  (Part of reducing from
    # degree n to degree k is computing the exact inverse of a k by k
    # rational matrix.  Fortunately, this matrix depends only on
    # n and k, so the inverted matrix can be cached; but still,
    # computing the exact inverse of a k by k matrix seems infeasible
    # for k much larger than 30.)

    if n <= 2:
        return None
    next = n * 3 // 4
    if next > 30:
        next = 30
    return next


def precompute_degree_reduction_cache(n):
    """
    Compute and cache the matrices used for degree reduction, starting
    from degree n.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: precompute_degree_reduction_cache(5)
        sage: dr_cache[5]
        (
           [121/126    8/63    -1/9   -2/63  11/126   -2/63]
           [   -3/7   37/42   16/21    1/21    -3/7     1/6]
           [    1/6    -3/7    1/21   16/21   37/42    -3/7]
        3, [  -2/63  11/126   -2/63    -1/9    8/63 121/126], 2,
        <BLANKLINE>
        ([121  16 -14  -4  11  -4]
        [-54 111  96   6 -54  21]
        [ 21 -54   6  96 111 -54]
        [ -4  11  -4 -14  16 121], 126)
        )
    """
    while True:
        if n in dr_cache:
            return
        next = degree_reduction_next_size(n)
        if next is None:
            dr_cache[n] = (None, None, 0)
            return

        # We can compute a degree n -> degree k reduction by looking at
        # as few as k+1 of the coefficients of the degree n polynomial.
        # Using more samples reduces the error involved in degree
        # reduction, but is slower.  More testing might reveal a better
        # way to select the number of samples here.

        samps = min(n+1, next+5)
        bd = bernstein_down(next, n, samps)

        # For each coefficient in the reduced polynomial, we see how
        # it varies as the (sampled) coefficients in the original
        # polynomial change by +/- 1.  Then we take the reduced
        # coefficient which changes the most, and call that the expected
        # error.  We can multiply this by the error in the original
        # polynomial, and be fairly certain (absolutely certain?) that
        # the error in the reduced polynomial will be no better
        # than this product.
        expected_err = max([sum([abs(x) for x in bd.row(k)]) for k in range(next+1)])

        # bdd = bd.denominator()
        # bdi = MatrixSpace(ZZ, next+1, samps, sparse=False)(bd * bdd)
        (bdi, bdd) = bd._clear_denom()

        dr_cache[n] = (next, bd, expected_err.floor(), (bdi, bdd))
        n = next


def bernstein_down(d1, d2, s):
    """
    Given polynomial degrees d1 and d2 (where d1 < d2), and a number
    of samples s, computes a matrix bd.

    If you have a Bernstein polynomial of formal degree d2, and select
    s of its coefficients (according to subsample_vec), and multiply
    the resulting vector by bd, then you get the coefficients
    of a Bernstein polynomial of formal degree d1, where this second
    polynomial is a good approximation to the first polynomial over the
    region of the Bernstein basis.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: bernstein_down(3, 8, 5)
        [ 612/245 -348/245   -37/49  338/245 -172/245]
        [-724/441   132/49  395/441 -290/147  452/441]
        [ 452/441 -290/147  395/441   132/49 -724/441]
        [-172/245  338/245   -37/49 -348/245  612/245]
    """

    # The use of the pseudoinverse means that the coefficients of the
    # reduced polynomial are selected according to a least-squares fit.
    # We would prefer to minimize the maximum error, rather than the
    # sum of the squares of the errors; but I don't know how to do that.

    # Also, this pseudoinverse is very slow to compute if d1 is large.
    # Since the coefficients of bernstein_up(...) can be represented
    # with a fairly simple formula (see the implementation of
    # bernstein_up()), it seems possible that there is some sort of
    # formula for the coefficients of the pseudoinverse that would be
    # faster than the computation here.

    return pseudoinverse(bernstein_up(d1, d2, s))


def pseudoinverse(m):
    mt = m.transpose()
    return ~(mt * m) * mt


def bernstein_up(d1, d2, s=None):
    """
    Given polynomial degrees d1 and d2, where d1 < d2, compute a matrix bu.

    If you have a Bernstein polynomial of formal degree d1, and
    multiply its coefficient vector by bu, then the result is the
    coefficient vector of the same polynomial represented as a
    Bernstein polynomial of formal degree d2.

    If s is not None, then it represents a number of samples; then the
    product only gives s of the coefficients of the new Bernstein polynomial,
    selected according to subsample_vec.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: bernstein_down(3, 7, 4)
        [  12/5     -4      3   -2/5]
        [-13/15   16/3     -4   8/15]
        [  8/15     -4   16/3 -13/15]
        [  -2/5      3     -4   12/5]
    """
    if s is None:
        s = d1 + 1
    MS = MatrixSpace(QQ, s, d1+1, sparse=False)
    m = MS()
    scale = factorial(d2)/factorial(d2-d1)
    for b in range(d1 + 1):
        scale2 = scale / binomial(d1, b)
        if (d1 - b) & 1 == 1:
            scale2 = -scale2
        scale2 = ~scale2
        for a in range(s):
            ra = ZZ(subsample_vec(a, s, d2 + 1))
            m[a, b] = prod([ra-i for i in range(b)]) * prod([ra-i for i in range(d2-d1+b+1, d2+1)]) * scale2

    return m


cdef int subsample_vec(int a, int slen, int llen) noexcept:
    """
    Given a vector of length llen, and slen < llen, we want to
    select slen of the elements of the vector, evenly spaced.

    Given an index 0 <= a < slen, this function computes the index
    in the original vector of length llen corresponding to a.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: [subsample_vec_doctest(a, 10, 100) for a in range(10)]
        [5, 15, 25, 35, 45, 54, 64, 74, 84, 94]
        sage: [subsample_vec_doctest(a, 3, 4) for a in range(3)]
        [1, 2, 3]
    """
    # round((a + 0.5) * (llen - 1) / slen)
    # round((2*a + 1) * (llen - 1) / (2 * slen)
    # floor(((2*a + 1) * (llen - 1) + slen) / (2 * slen))
    return ((2*a + 1) * (llen - 1) + slen) // (2 * slen)


def subsample_vec_doctest(a, slen, llen):
    return subsample_vec(a, slen, llen)


def maximum_root_first_lambda(p):
    r"""
    Given a polynomial with real coefficients, computes an upper bound
    on its largest real root.

    This is using the first-\lambda algorithm from
    "Implementations of a New Theorem for Computing Bounds for Positive
    Roots of Polynomials", by Akritas, Strzebo\'nski, and Vigklas.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: maximum_root_first_lambda((x-1)*(x-2)*(x-3))
        6.00000000000001
        sage: maximum_root_first_lambda((x+1)*(x+2)*(x+3))
        0.000000000000000
        sage: maximum_root_first_lambda(x^2 - 1)
        1.00000000000000
    """
    n = p.degree()
    if p[n] < 0:
        p = -p
    cl = [RIF(x) for x in p.list()]
    return cl_maximum_root_first_lambda(cl)


def cl_maximum_root_first_lambda(cl):
    r"""
    Given a polynomial represented by a list of its coefficients
    (as RealIntervalFieldElements), compute an upper bound on its
    largest real root.

    Uses the first-\lambda algorithm from "Implementations of a New Theorem
    for Computing Bounds for Positive Roots of Polynomials",
    by Akritas, Strzebo\'nski, and Vigklas.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: cl_maximum_root_first_lambda([RIF(-1), RIF(0), RIF(1)])
        1.00000000000000

    TESTS::

        sage: bnd = cl_maximum_root_first_lambda(list(map(RIF, [0, 0, 0, 14, 1])))
        sage: bnd, bnd.parent()
        (0.000000000000000,
        Real Field with 53 bits of precision and rounding RNDU)
    """
    n = len(cl) - 1
    assert cl[n] > 0
    pending_pos_coeff = cl[n]
    pending_pos_exp = n
    lastPos = True
    posCounter = 1
    negCounter = 0
    pos = []
    neg = []
    for j in range(n-1, -2, -1):
        if j < 0:
            coeff = 1
        else:
            coeff = cl[j]
        if coeff < 0:
            neg += [(coeff, j)]
            lastPos = False
            negCounter = negCounter+1
        if coeff > 0:
            if negCounter > posCounter:
                chunks = negCounter - posCounter + 1
                pos += [(pending_pos_coeff / chunks, pending_pos_exp)] * chunks
                posCounter = negCounter
            else:
                pos += [(pending_pos_coeff, pending_pos_exp)]
            pending_pos_coeff = coeff
            pending_pos_exp = j
            posCounter = posCounter+1

    if not neg:
        return RIF.upper_field().zero()

    max_ub_log = RIF('-infinity')
    for j in range(len(neg)):
        cur_ub_log = (-neg[j][0] / pos[j][0]).log() / (pos[j][1] - neg[j][1])
        max_ub_log = max_ub_log.union(cur_ub_log)

    return max_ub_log.upper().exp()


def maximum_root_local_max(p):
    r"""
    Given a polynomial with real coefficients, computes an upper bound
    on its largest real root, using the local-max algorithm from
    "Implementations of a New Theorem for Computing Bounds for Positive
    Roots of Polynomials", by Akritas, Strzebo\'nski, and Vigklas.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: maximum_root_local_max((x-1)*(x-2)*(x-3))
        12.0000000000001
        sage: maximum_root_local_max((x+1)*(x+2)*(x+3))
        0.000000000000000
        sage: maximum_root_local_max(x^2 - 1)
        1.41421356237310
    """
    n = p.degree()
    if p[n] < 0:
        p = -p
    cl = [RIF(x) for x in p.list()]
    return cl_maximum_root_local_max(cl)


def cl_maximum_root_local_max(cl):
    r"""
    Given a polynomial represented by a list of its coefficients
    (as RealIntervalFieldElements), compute an upper bound on its
    largest real root.

    Uses the local-max algorithm from "Implementations of a New Theorem
    for Computing Bounds for Positive Roots of Polynomials",
    by Akritas, Strzebo\'nski, and Vigklas.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: cl_maximum_root_local_max([RIF(-1), RIF(0), RIF(1)])
        1.41421356237310
    """
    n = len(cl) - 1
    assert cl[n] > 0
    max_pos_coeff = cl[n]
    max_pos_exp = n
    max_pos_uses = 0
    max_ub_log = RIF('-infinity')

    for j in range(n-1, -1, -1):
        if cl[j] < 0:
            max_pos_uses = max_pos_uses+1
            cur_ub_log = (-cl[j] / (max_pos_coeff >> max_pos_uses)).log() / (max_pos_exp - j)
            max_ub_log = max_ub_log.union(cur_ub_log)
        elif cl[j] > max_pos_coeff:
            max_pos_coeff = cl[j]
            max_pos_exp = j
            max_pos_uses = 0

    return max_ub_log.upper().exp()


def cl_maximum_root(cl):
    r"""
    Given a polynomial represented by a list of its coefficients
    (as RealIntervalFieldElements), compute an upper bound on its
    largest real root.

    Uses two algorithms of Akritas, Strzebo\'nski, and Vigklas, and
    picks the better result.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: cl_maximum_root([RIF(-1), RIF(0), RIF(1)])
        1.00000000000000
    """
    n = len(cl) - 1
    if cl[n] < 0:
        cl = [-x for x in cl]
    return min(cl_maximum_root_first_lambda(cl),
               cl_maximum_root_local_max(cl))


def root_bounds(p):
    r"""
    Given a polynomial with real coefficients, computes a lower and
    upper bound on its real roots.  Uses algorithms of
    Akritas, Strzebo\'nski, and Vigklas.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: root_bounds((x-1)*(x-2)*(x-3))
        (0.545454545454545, 6.00000000000001)
        sage: root_bounds(x^2)
        (0.000000000000000, 0.000000000000000)
        sage: root_bounds(x*(x+1))
        (-1.00000000000000, 0.000000000000000)
        sage: root_bounds((x+2)*(x-3))
        (-2.44948974278317, 3.46410161513776)
        sage: root_bounds(x^995 * (x^2 - 9999) - 1)
        (-99.9949998749937, 141.414284992713)
        sage: root_bounds(x^995 * (x^2 - 9999) + 1)
        (-141.414284992712, 99.9949998749938)

    If we can see that the polynomial has no real roots, return None. ::

        sage: root_bounds(x^2 + 1) is None
        True
    """
    n = p.degree()
    if p[n] < 0:
        p = -p
    cl = [RIF(x) for x in p.list()]

    cdef int zero_roots = 0
    while cl[0] == 0:
        cl.pop(0)
        zero_roots = zero_roots + 1
        n = n-1

    if n == 0:
        # not RIF.zero().endpoints() because of MPFI's convention that the
        # upper bound is -0.
        return RIF.lower_field().zero(), RIF.upper_field().zero()

    ub = cl_maximum_root(cl)

    neg_cl = copy(cl)
    for j in range(n-1, -1, -2):
        neg_cl[j] = -neg_cl[j]

    lb = -cl_maximum_root(neg_cl)

    if ub == 0 and lb == 0:
        if zero_roots == 0:
            return None
        else:
            return (lb, ub)

    if ub == 0 and zero_roots == 0:
        swap_neg_cl = copy(neg_cl)
        swap_neg_cl.reverse()
        ub = (-(~RIF(cl_maximum_root(swap_neg_cl)))).upper()

    if lb == 0 and zero_roots == 0:
        swap_cl = copy(cl)
        swap_cl.reverse()
        lb = (~RIF(cl_maximum_root(swap_cl))).lower()

    return (lb, ub)


def rational_root_bounds(p):
    """
    Given a polynomial p with real coefficients, computes rationals
    a and b, such that for every real root r of p, a < r < b.
    We try to find rationals which bound the roots somewhat tightly,
    yet are simple (have small numerators and denominators).

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: rational_root_bounds((x-1)*(x-2)*(x-3))
        (0, 7)
        sage: rational_root_bounds(x^2)
        (-1/2, 1/2)
        sage: rational_root_bounds(x*(x+1))
        (-3/2, 1/2)
        sage: rational_root_bounds((x+2)*(x-3))
        (-3, 6)
        sage: rational_root_bounds(x^995 * (x^2 - 9999) - 1)
        (-100, 1000/7)
        sage: rational_root_bounds(x^995 * (x^2 - 9999) + 1)
        (-142, 213/2)

    If we can see that the polynomial has no real roots, return None.
        sage: rational_root_bounds(x^2 + 7) is None
        True
    """

    # There are two stages to the root isolation process.  First,
    # we convert the polynomial to the Bernstein basis given by
    # the root bounds, exactly; then we isolate the roots from
    # that Bernstein basis, using interval arithmetic.

    # We want the root bounds to be as small as possible, because that
    # speeds up the second phase; but we also want them to be as
    # simple as possible, because that speeds up the first phase.

    # As a heuristic, given a polynomial of degree d with floating-point
    # root bounds lb and ub, we compute simple rational root bounds
    # rlb and rub such that lb - (ub - lb)/sqrt(d) <= rlb <= lb,
    # and ub <= rub <= ub + (ub - lb)/sqrt(d).  Very little testing
    # was done to come up with this heuristic, and it can probably
    # be improved.

    if p.degree() == 0:
        return (QQ(-1), QQ(1))

    sqrtd = sqrt(RR(p.degree()))
    bounds = root_bounds(p)
    if bounds is None:
        return None
    (lb, ub) = bounds

    if lb == ub:
        # The code below would give overly precise bounds in this case,
        # giving very non-simple isolating intervals in the result.
        # We don't need tight root bounds to quickly find the roots
        # of a linear polynomial, so go for simple root bounds.

        b = QQ(lb)

        return (b - QQ_1_2, b + QQ_1_2)

    # XXX A gross hack.  Sometimes, root_bounds() gives a lower or
    # upper bound which is a root.  We want bounds which are
    # not roots, so we just spread out the bounds a tiny bit.
    # (It might be more efficient to check the bounds from root_bounds()
    # and, if they are roots, use that information to factor out a linear
    # polynomial.  However, as far as I can tell, root_bounds() only
    # gives roots as bounds on toy examples; so this is not too inefficient.)
    lb = RR(lb).nextbelow()
    ub = RR(ub).nextabove()

    # Given the current implementation of to_bernstein, we want rlb
    # and (rub/rlb - 1) to be simple rationals -- we don't really care
    # about the simplicity of rub by itself.  (Or else we want rlb to
    # be zero, and rub to be a simple rational.)

    ilb = RIF(lb - (ub - lb)/sqrtd, lb)
    rlb = ilb.simplest_rational()

    if rlb == 0:
        iub = RIF(ub, ub + (ub - lb)/sqrtd)
        rub = iub.simplest_rational()
        return (rlb, rub)

    # We want to compute an interval for the upper bound,
    # iub = [ub .. ub + (ub - lb)/sqrtd],
    # and then compute (iub/ilb - 1).simplest_rational().
    # However, we want to compute this
    # interval with inward instead of outward rounding.  Instead,
    # we compute the lower and upper bounds of the interval separately.

    iub_l = RIF(ub)
    iub_h = RIF(ub + (ub - lb)/sqrtd)

    iub_l2 = iub_l/rlb - 1
    iub_h2 = iub_h/rlb - 1

    if iub_l2 < iub_h2:
        rub = (RIF(iub_l2.upper(), iub_h2.lower()).simplest_rational() + 1)*rlb
    else:
        rub = (RIF(iub_h2.upper(), iub_l2.lower()).simplest_rational() + 1)*rlb

    return (rlb, rub)


class PrecisionError(ValueError):
    pass


class bernstein_polynomial_factory:
    """
    An abstract base class for bernstein_polynomial factories.  That
    is, elements of subclasses represent Bernstein polynomials
    (exactly), and are responsible for creating
    interval_bernstein_polynomial_integer approximations at arbitrary
    precision.

    Supports four methods, coeffs_bitsize(), bernstein_polynomial(),
    lsign(), and usign().  The coeffs_bitsize() method gives an
    integer approximation to the log2 of the max of the absolute
    values of the Bernstein coefficients.  The
    bernstein_polynomial(scale_log2) method gives an approximation
    where the maximum coefficient has approximately coeffs_bitsize() -
    scale_log2 bits.  The lsign() and usign() methods give the (exact)
    sign of the first and last coefficient, respectively.
    """

    def _sign(self, v):
        if v > 0:
            return 1
        if v < 0:
            return -1
        return 0

    def lsign(self):
        """
        Return the sign of the first coefficient of this
        Bernstein polynomial.
        """
        return self._sign(self.coeffs[0])

    def usign(self):
        """
        Return the sign of the last coefficient of this
        Bernstein polynomial.
        """
        return self._sign(self.coeffs[-1])


class bernstein_polynomial_factory_intlist(bernstein_polynomial_factory):
    """
    This class holds an exact Bernstein polynomial (represented
    as a list of integer coefficients), and returns arbitrarily-precise
    interval approximations of this polynomial on demand.
    """

    def __init__(self, coeffs):
        """
        Initialize a ``bernstein_polynomial_factory_intlist``,
        given a list of integer coefficients.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bernstein_polynomial_factory_intlist([1, 2, 3])
            degree 2 Bernstein factory with 2-bit integer coefficients
        """
        self.coeffs = coeffs

    def __repr__(self):
        """
        Return a short summary of this Bernstein polynomial factory.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bernstein_polynomial_factory_intlist([1, 2, 3, 4000])
            degree 3 Bernstein factory with 12-bit integer coefficients
        """

        return "degree %d Bernstein factory with %d-bit integer coefficients" % (len(self.coeffs) - 1, self.coeffs_bitsize())

    def coeffs_bitsize(self):
        """
        Compute the approximate log2 of the maximum of the absolute
        values of the coefficients.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bernstein_polynomial_factory_intlist([1, 2, 3, -60000]).coeffs_bitsize()
            16
        """
        b = self.coeffs
        return max([bitsize(c) for c in b])
#        return intlist_size(self.coeffs)

    def bernstein_polynomial(self, scale_log2):
        """
        Compute an interval_bernstein_polynomial_integer that approximates
        this polynomial, using the given scale_log2.  (Smaller scale_log2
        values give more accurate approximations.)

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bpf = bernstein_polynomial_factory_intlist([10, -20, 30, -40])
            sage: print(bpf.bernstein_polynomial(0))
            degree 3 IBP with 6-bit coefficients
            sage: bpf.bernstein_polynomial(20)
            <IBP: ((0, -1, 0, -1) + [0 .. 1)) * 2^20; lsign 1>
            sage: bpf.bernstein_polynomial(0)
            <IBP: (10, -20, 30, -40) + [0 .. 1)>
            sage: bpf.bernstein_polynomial(-20)
            <IBP: ((10485760, -20971520, 31457280, -41943040) + [0 .. 1)) * 2^-20>
        """
        b = self.coeffs
        if scale_log2 < 0:
            shift = ZZ(-scale_log2)
            intv_b = [c << shift for c in b]
        else:
            shift = ZZ(scale_log2)
            intv_b = [c >> shift for c in b]

        return interval_bernstein_polynomial_integer((ZZ ** len(b))(intv_b), zero_QQ, one_QQ, self.lsign(), self.usign(), 1, scale_log2, 0, zero_RIF)
#        return bp_of_intlist(self.coeffs, scale_log2)


class bernstein_polynomial_factory_ratlist(bernstein_polynomial_factory):
    """
    This class holds an exact Bernstein polynomial (represented
    as a list of rational coefficients), and returns arbitrarily-precise
    interval approximations of this polynomial on demand.
    """

    def __init__(self, coeffs):
        """
        Initialize a ``bernstein_polynomial_factory_intlist``,
        given a list of rational coefficients.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bernstein_polynomial_factory_ratlist([1, 1/2, 1/3])
            degree 2 Bernstein factory with 0-bit rational coefficients
        """
        self.coeffs = coeffs

    def __repr__(self):
        """
        Return a short summary of this Bernstein polynomial factory.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bernstein_polynomial_factory_ratlist([1, 2, 3, 4000/17])
            degree 3 Bernstein factory with 7-bit rational coefficients
        """

        return "degree %d Bernstein factory with %d-bit rational coefficients" % (len(self.coeffs) - 1, self.coeffs_bitsize())

    def coeffs_bitsize(self):
        """
        Compute the approximate log2 of the maximum of the absolute
        values of the coefficients.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bernstein_polynomial_factory_ratlist([1, 2, 3, -60000]).coeffs_bitsize()
            15
            sage: bernstein_polynomial_factory_ratlist([65535/65536]).coeffs_bitsize()
            -1
            sage: bernstein_polynomial_factory_ratlist([65536/65535]).coeffs_bitsize()
            1
        """

        # This returns an estimate of the log base 2 of the rational in the
        # list with largest absolute value.  The estimate is an integer,
        # and within +/- 1 of the correct answer.
        r = max([bitsize(c.numerator()) - bitsize(c.denominator()) for c in self.coeffs])
        return r

    def bernstein_polynomial(self, scale_log2):
        """
        Compute an interval_bernstein_polynomial_integer that approximates
        this polynomial, using the given scale_log2.  (Smaller scale_log2
        values give more accurate approximations.)

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: bpf = bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99])
            sage: print(bpf.bernstein_polynomial(0))
            degree 3 IBP with 3-bit coefficients
            sage: bpf.bernstein_polynomial(20)
            <IBP: ((0, -1, 0, -1) + [0 .. 1)) * 2^20; lsign 1>
            sage: bpf.bernstein_polynomial(0)
            <IBP: (0, -4, 2, -2) + [0 .. 1); lsign 1>
            sage: bpf.bernstein_polynomial(-20)
            <IBP: ((349525, -3295525, 2850354, -1482835) + [0 .. 1)) * 2^-20>
        """
        b = self.coeffs

        if scale_log2 < 0:
            intv_b = [(numerator(c) << (-scale_log2)) // denominator(c) for c in b]
        else:
            intv_b = [(numerator(c) >> scale_log2) // denominator(c) for c in b]

        return interval_bernstein_polynomial_integer((ZZ ** len(b))(intv_b), zero_QQ, one_QQ, self.lsign(), self.usign(), 1, scale_log2, 0, zero_RIF)
#        return bp_of_ratlist(self.coeffs, scale_log2)


class bernstein_polynomial_factory_ar(bernstein_polynomial_factory):
    """
    This class holds an exact Bernstein polynomial (represented as a
    list of algebraic real coefficients), and returns
    arbitrarily-precise interval approximations of this polynomial on
    demand.
    """

    def __init__(self, poly, neg):
        """
        Initialize a ``bernstein_polynomial_factory_ar``,
        given a polynomial with algebraic real coefficients.
        If neg is True, then gives the Bernstein polynomial for
        the negative half-line; if neg is False, the positive.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: x = polygen(AA)
            sage: p = (x - 1) * (x - sqrt(AA(2))) * (x - 2)
            sage: bernstein_polynomial_factory_ar(p, False)
            degree 3 Bernstein factory with algebraic real coefficients
        """
        coeffs = to_bernstein_warp(poly)
        if neg:
            for i in range(1, len(coeffs), 2):
                coeffs[i] = -coeffs[i]
        sizes = []
        for c in coeffs:
            size = RIF(c).magnitude()
            if size > 0:
                sizes.append(size.log2().floor())
            else:
                sizes.append(-1000000)

        self.coeffs = coeffs
        self.sizes = sizes

    def __repr__(self):
        """
        Return a short summary of this Bernstein polynomial factory.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: x = polygen(AA)
            sage: p = (x - 1) * (x - sqrt(AA(2))) * (x - 2)
            sage: bernstein_polynomial_factory_ar(p, False)
            degree 3 Bernstein factory with algebraic real coefficients
        """

        return "degree %d Bernstein factory with algebraic real coefficients" % (len(self.coeffs) - 1)

    def coeffs_bitsize(self):
        """
        Compute the approximate log2 of the maximum of the absolute
        values of the coefficients.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: x = polygen(AA)
            sage: p = (x - 1) * (x - sqrt(AA(2))) * (x - 2)
            sage: bernstein_polynomial_factory_ar(p, False).coeffs_bitsize()
            1
        """
        return max(self.sizes)

    def bernstein_polynomial(self, scale_log2):
        """
        Compute an interval_bernstein_polynomial_integer that approximates
        this polynomial, using the given scale_log2.  (Smaller scale_log2
        values give more accurate approximations.)

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: x = polygen(AA)
            sage: p = (x - 1) * (x - sqrt(AA(2))) * (x - 2)
            sage: bpf = bernstein_polynomial_factory_ar(p, False)
            sage: print(bpf.bernstein_polynomial(0))
            degree 3 IBP with 2-bit coefficients
            sage: bpf.bernstein_polynomial(-20)
            <IBP: ((-2965821, 2181961, -1542880, 1048576) + [0 .. 1)) * 2^-20>
            sage: bpf = bernstein_polynomial_factory_ar(p, True)
            sage: bpf.bernstein_polynomial(-20)
            <IBP: ((-2965821, -2181962, -1542880, -1048576) + [0 .. 1)) * 2^-20>
            sage: p = x^2 - 1
            sage: bpf = bernstein_polynomial_factory_ar(p, False)
            sage: bpf.bernstein_polynomial(-10)
            <IBP: ((-1024, 0, 1024) + [0 .. 1)) * 2^-10>
        """
        res = (ZZ ** len(self.coeffs))(0)
        max_err = 0

        for i in range(len(self.coeffs)):
            intv_c = RealIntervalField(max(self.sizes[i] - scale_log2 + 6, 6))(self.coeffs[i])
            if scale_log2 < 0:
                intv_c = intv_c << (-scale_log2)
            else:
                intv_c = intv_c >> scale_log2

            # Compute lower and upper such that lower <= intv_c < upper
            lower = intv_c.lower().floor()
            upper = intv_c.upper().floor() + 1

            res[i] = lower
            max_err = max(max_err, upper - lower)

        return interval_bernstein_polynomial_integer(res, zero_QQ, one_QQ, self.lsign(), self.usign(), max_err, scale_log2, 0, zero_RIF)


def split_for_targets(context ctx, interval_bernstein_polynomial bp, target_list, precise=False):
    """
    Given an interval Bernstein polynomial over a particular region
    (assumed to be a (not necessarily proper) subregion of [0 .. 1]),
    and a list of targets, uses de Casteljau's method to compute
    representations of the Bernstein polynomial over each target.
    Uses degree reduction as often as possible while maintaining the
    requested precision.

    Each target is of the form (lgap, ugap, b).  Suppose lgap.region()
    is (l1, l2), and ugap.region() is (u1, u2).  Then we will compute
    an interval Bernstein polynomial over a region [l .. u], where
    l1 <= l <= l2 and u1 <= u <= u2.  (split_for_targets() is free to
    select arbitrary region endpoints within these bounds; it picks
    endpoints which make the computation easier.)  The third component
    of the target, b, is the maximum allowed scale_log2 of the result;
    this is used to decide when degree reduction is allowed.

    The pair (l1, l2) can be replaced by None, meaning [-infinity .. 0];
    or, (u1, u2) can be replaced by None, meaning [1 .. infinity].

    There is another constraint on the region endpoints selected by
    split_for_targets() for a target ((l1, l2), (u1, u2), b).
    We set a size goal g, such that (u - l) <= g * (u1 - l2).
    Normally g is 256/255, but if precise is True, then g is 65536/65535.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: bp = mk_ibpi([1000000, -2000000, 3000000, -4000000, -5000000, -6000000])
        sage: ctx = mk_context()
        sage: bps = split_for_targets(ctx, bp, [(rr_gap(1/1234567893, 1/1234567892, 1), rr_gap(1/1234567891, 1/1234567890, 1), 12), (rr_gap(1/3, 1/2, -1), rr_gap(2/3, 3/4, -1), 6)])
        sage: bps[0]
        <IBP: (999992, 999992, 999992) + [0 .. 15) over [8613397477114467984778830327/10633823966279326983230456482242756608 .. 591908168025934394813836527495938294787/730750818665451459101842416358141509827966271488]; level 2; slope_err 0.?e12>
        sage: bps[1]
        <IBP: (-1562500, -1875001, -2222223, -2592593, -2969137, -3337450) + [0 .. 4) over [1/2 .. 2863311531/4294967296]>
    """
    ntargs = len(target_list)
    if ntargs == 0:
        return []

    bounds = bp.region()

    out_of_bounds = False

    cdef rr_gap l
    cdef rr_gap r

    split_targets = []
    for l, r, _ in target_list:
        if l is None:
            split_targets += [(QQ(0), None, 0)]
        else:
            lbounds = relative_bounds(bounds, l.region())
            split_targets += [(lbounds[1], lbounds[0], l.sign)]
            if lbounds[0] > 0:
                out_of_bounds = True
        if r is None:
            split_targets += [(QQ.one(), None, 0)]
        else:
            rbounds = relative_bounds(bounds, r.region())
            split_targets += [(rbounds[0], rbounds[1], r.sign)]
            if rbounds[1] < 1:
                out_of_bounds = True

    if precise:
        goal = Integer(65535)/65536
    else:
        goal = Integer(255)/256

    if ntargs == 1 and not(out_of_bounds) and split_targets[1][0] - split_targets[0][0] >= goal:
        return [bp]

    half = Integer(1)/2
    best_target = split_targets[0][0]
    best_index = 0
    for i in range(1, len(split_targets)):
        if abs(split_targets[i][0] - half) < abs(best_target - half):
            best_target = split_targets[i][0]
            best_index = i

    split = wordsize_rational(split_targets[best_index][0], split_targets[best_index][1], ctx.wordsize)

    (p1_, p2_, ok) = bp.de_casteljau(ctx, split, msign=split_targets[best_index][2])
    assert ok

    cdef interval_bernstein_polynomial p1 = p1_
    cdef interval_bernstein_polynomial p2 = p2_

    if bp.lft is not None:
        (a, b, c, d) = bp.lft
        left = b/d
        if c+d != 0:
            right = (a+b)/(c+d)
            width = right-left
            err = (a/2 + b)/(c/2 + d) - (left+right)/2
        else:
            width = RR(infinity)
            err = RR(infinity)
        slope = (a*d - b*c)/(d*d)
        ctx.dc_log_append(('lft', a, b, c, d, RR(left), RR(width), RR(slope/width), RR(err/width/width)))
    ctx.dc_log_append(('split_for_targets', best_index, split, bp.scale_log2, bp.bitsize, p1.bitsize, p2.bitsize))

    target_list_splitpoint = (best_index + 1) // 2
    tl1 = target_list[:target_list_splitpoint]
    tl2 = target_list[target_list_splitpoint:]

    if len(tl1) > 0:
        if True:  # p1.region_width() / bp.region_width() < tiny:
            max_lsb = max([t[2] for t in tl1])
            p1 = p1.down_degree_iter(ctx, max_lsb)
        r1 = split_for_targets(ctx, p1, tl1)
    else:
        r1 = []
    if len(tl2) > 0:
        if True:  # p2.region_width() / bp.region_width() < tiny:
            max_lsb = max([t[2] for t in tl2])
            p2 = p2.down_degree_iter(ctx, max_lsb)
        r2 = split_for_targets(ctx, p2, tl2)
    else:
        r2 = []

    return r1 + r2


cdef class ocean:
    """
    Given the tools we've defined so far, there are many possible root
    isolation algorithms that differ on where to select split points,
    what precision to work at when, and when to attempt degree
    reduction.

    Here we implement one particular algorithm, which I call the
    ocean-island algorithm.  We start with an interval Bernstein
    polynomial defined over the region [0 .. 1].  This region is
    the "ocean".  Using de Casteljau's algorithm and Descartes'
    rule of signs, we divide this region into subregions which may
    contain roots, and subregions which are guaranteed not to contain
    roots.  Subregions which may contain roots are "islands"; subregions
    known not to contain roots are "gaps".

    All the real root isolation work happens in class island.  See the
    documentation of that class for more information.

    An island can be told to refine itself until it contains only a
    single root.  This may not succeed, if the island's interval
    Bernstein polynomial does not have enough precision.  The ocean
    basically loops, refining each of its islands, then increasing the
    precision of islands which did not succeed in isolating a single
    root; until all islands are done.

    Increasing the precision of unsuccessful islands is done in a
    single pass using split_for_target(); this means it is possible
    to share work among multiple islands.
    """

    def __init__(self, ctx, bpf, mapping):
        """
        Initialize an ocean from a context and a Bernstein polynomial
        factory.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            ocean with precision 120 and 1 island(s)
        """

        # Islands and gaps are maintained in a doubly-linked (circular)
        # list.  Islands point to the gaps on either side, and gaps
        # point to the islands on either side.

        # The list starts with self.lgap, which is the gap that starts
        # at 0; it ends at self.endpoint, which is a fictional island
        # after the gap that ends at 1.

        # The initial gaps are unique in being 0-width; all gaps that
        # are created during execution of the algorithm have positive
        # width.

        # Note that the constructor for islands side-effects its argument
        # gaps, so that they point to the island as their neighbor.

        self.ctx = ctx
        self.bpf = bpf
        self.mapping = mapping
        self.lgap = rr_gap(zero_QQ, zero_QQ, bpf.lsign())
        rgap = rr_gap(one_QQ, one_QQ, bpf.usign())
        self.endpoint = island(None, rgap, self.lgap)
        self.msb = bpf.coeffs_bitsize()
        self.prec = 120

        cdef island isle = island(self.approx_bp(self.msb - self.prec), self.lgap, rgap)
        if isle.bp.max_variations == 0:
            self.lgap.risland = self.endpoint
            self.endpoint.lgap = self.lgap

    def __repr__(self):
        """
        Return a short summary of this root isolator.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            ocean with precision 120 and 1 island(s)
        """

        cdef int isle_count = 0
        cdef island isle = self.lgap.risland
        while isle is not self.endpoint:
            isle_count = isle_count + 1
            isle = isle.rgap.risland

        return "ocean with precision %d and %d island(s)" % (self.prec, isle_count)

    def _islands(self):
        """
        Return a list of the islands in this ocean (for debugging purposes).

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc._islands()
            [island of width 1.00000000000000]
        """

        islands = []
        cdef island isle = self.lgap.risland
        while isle is not self.endpoint:
            islands.append(isle)
            isle = isle.rgap.risland

        return islands

    def approx_bp(self, scale_log2):
        """
        Return an approximation to our Bernstein polynomial with the
        given scale_log2.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc.approx_bp(0)
            <IBP: (0, -4, 2, -2) + [0 .. 1); lsign 1>
            sage: oc.approx_bp(-20)
            <IBP: ((349525, -3295525, 2850354, -1482835) + [0 .. 1)) * 2^-20>
        """
        return self.bpf.bernstein_polynomial(scale_log2)

    def find_roots(self):
        """
        Isolate all roots in this ocean.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc
            ocean with precision 120 and 1 island(s)
            sage: oc.find_roots()
            sage: oc
            ocean with precision 120 and 3 island(s)
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1, 0, -1111/2, 0, 11108889/14, 0, 0, 0, 0, -1]), lmap)
            sage: oc.find_roots()
            sage: oc
            ocean with precision 240 and 3 island(s)
        """
        while not self.all_done():
            self.refine_all()
            self.increase_precision()

    def roots(self):
        """
        Return the locations of all islands in this ocean.  (If run
        after find_roots(), this is the location of all roots in the ocean.)

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc.find_roots()
            sage: oc.roots()
            [(1/32, 1/16), (1/2, 5/8), (3/4, 7/8)]
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1, 0, -1111/2, 0, 11108889/14, 0, 0, 0, 0, -1]), lmap)
            sage: oc.find_roots()
            sage: oc.roots()
            [(95761241267509487747625/9671406556917033397649408, 191522482605387719863145/19342813113834066795298816), (1496269395904347376805/151115727451828646838272, 374067366568272936175/37778931862957161709568), (31/32, 63/64)]
        """
        rts = []
        cdef island isle = self.lgap.risland
        while isle is not self.endpoint:
            rts.append((isle.lgap.upper, isle.rgap.lower))
            isle = isle.rgap.risland
        return rts

    def refine_all(self):
        """
        Refine all islands which are not done (which are not known to
        contain exactly one root).

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc
            ocean with precision 120 and 1 island(s)
            sage: oc.refine_all()
            sage: oc
            ocean with precision 120 and 3 island(s)
        """
        cdef island isle = self.lgap.risland
        while isle is not self.endpoint:
            if not isle.done(self.ctx):
                isle.refine(self.ctx)
            isle = isle.rgap.risland

    def all_done(self):
        """
        Return ``True`` iff all islands are known to contain exactly one root.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc.all_done()
            False
            sage: oc.find_roots()
            sage: oc.all_done()
            True
        """
        cdef island isle = self.lgap.risland
        while isle is not self.endpoint:
            if not isle.done(self.ctx):
                return False
            if not isle.has_root():
                isle.lgap.risland = isle.rgap.risland
                isle.rgap.risland.lgap = isle.lgap
                isle.lgap.upper = isle.rgap.upper
            isle = isle.rgap.risland
        return True

    def reset_root_width(self, int isle_num, target_width):
        """
        Require that the isle_num island have a width at most target_width.

        If this is followed by a call to find_roots(), then the
        corresponding root will be refined to the specified width.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([-1, -1, 1]), lmap)
            sage: oc.find_roots()
            sage: oc.roots()
            [(1/2, 3/4)]
            sage: oc.reset_root_width(0, 1/2^200)
            sage: oc.find_roots()
            sage: oc
            ocean with precision 240 and 1 island(s)
            sage: RR(RealIntervalField(300)(oc.roots()[0]).absolute_diameter()).log2()
            -232.668979560890
        """
        cdef island isle = self.lgap.risland
        cdef int n = 0
        while isle is not self.endpoint:
            if n == isle_num:
                isle.reset_root_width(target_width)

            isle = isle.rgap.risland
            n = n+1

    def increase_precision(self):
        """
        Increase the precision of the interval Bernstein polynomial held
        by any islands which are not done.  (In normal use, calls to this
        function are separated by calls to self.refine_all().)

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: oc = ocean(mk_context(), bernstein_polynomial_factory_ratlist([1/3, -22/7, 193/71, -140/99]), lmap)
            sage: oc
            ocean with precision 120 and 1 island(s)
            sage: oc.increase_precision()
            sage: oc.increase_precision()
            sage: oc.increase_precision()
            sage: oc
            ocean with precision 960 and 1 island(s)
        """
        active_islands = []
        cdef int total_islands = 0
        cdef island isle = self.lgap.risland
        while isle is not self.endpoint:
            total_islands = total_islands + 1
            if not isle.done(self.ctx) and len(isle.ancestors) == 0:
                active_islands += [isle]
            isle = isle.rgap.risland
        if len(active_islands) > 0:
            self.prec = self.prec * 2
            p = self.approx_bp(self.msb - self.prec)
            targets = []
            for isle in active_islands:
                targets += [(isle.lgap, isle.rgap, isle.bp.scale_log2)]
                self.ctx.dc_log_append(('splitting', (isle.lgap.lower, isle.lgap.upper), (isle.rgap.lower, isle.rgap.upper), isle.bp.scale_log2))
            bps = split_for_targets(self.ctx, p, targets)
            for i in range(len(active_islands)):
                bp = bps[i]
                isle = active_islands[i]
                isle.bp = bp


cdef class island:
    """
    This implements the island portion of my ocean-island root isolation
    algorithm.  See the documentation for class ocean, for more
    information on the overall algorithm.

    Island root refinement starts with a Bernstein polynomial whose
    region is the whole island (or perhaps slightly more than the
    island in certain cases).  There are two subalgorithms; one when
    looking at a Bernstein polynomial covering a whole island (so we
    know that there are gaps on the left and right), and one when
    looking at a Bernstein polynomial covering the left segment of an
    island (so we know that there is a gap on the left, but the right
    is in the middle of an island).  An important invariant of the
    left-segment subalgorithm over the region [l .. r] is that it
    always finds a gap [r0 .. r] ending at its right endpoint.

    Ignoring degree reduction, downscaling (precision reduction), and
    failures to split, the algorithm is roughly:

    Whole island:

    1. If the island definitely has exactly one root, then return.
    2. Split the island in (approximately) half.
    3. If both halves definitely have no roots, then remove this island from
       its doubly-linked list (merging its left and right gaps) and return.
    4. If either half definitely has no roots, then discard that half
       and call the whole-island algorithm with the other half, then return.
    5. If both halves may have roots, then call the left-segment algorithm
       on the left half.
    6. We now know that there is a gap immediately to the left of the
       right half, so call the whole-island algorithm on the right half,
       then return.

    Left segment:

    1. Split the left segment in (approximately) half.
    2. If both halves definitely have no roots, then extend the left gap
       over the segment and return.
    3. If the left half definitely has no roots, then extend the left gap
       over this half and call the left-segment algorithm on the right half,
       then return.
    4. If the right half definitely has no roots, then split the island
       in two, creating a new gap.  Call the whole-island algorithm on the
       left half, then return.
    5. Both halves may have roots.  Call the left-segment algorithm on
       the left half.
    6. We now know that there is a gap immediately to the left of the
       right half, so call the left-segment algorithm on the right half,
       then return.

    Degree reduction complicates this picture only slightly.  Basically,
    we use heuristics to decide when degree reduction might be likely
    to succeed and be helpful; whenever this is the case, we attempt
    degree reduction.

    Precision reduction and split failure add more complications.
    The algorithm maintains a stack of different-precision representations
    of the interval Bernstein polynomial.  The base of the stack
    is at the highest (currently known) precision; each stack entry has
    approximately half the precision of the entry below it.  When we
    do a split, we pop off the top of the stack, split it, then push
    whichever half we're interested in back on the stack (so the
    different Bernstein polynomials may be over different regions).
    When we push a polynomial onto the stack, we may heuristically decide to
    push further lower-precision versions of the same polynomial onto the
    stack.

    In the algorithm above, whenever we say "split in (approximately) half",
    we attempt to split the top-of-stack polynomial using try_split()
    and try_rand_split().  However, these will fail if the sign of the
    polynomial at the chosen split point is unknown (if the polynomial is
    not known to high enough precision, or if the chosen split point
    actually happens to be a root of the polynomial).  If this fails, then
    we discard the top-of-stack polynomial, and try again with the
    next polynomial down (which has approximately twice the precision).
    This next polynomial may not be over the same region; if not, we
    split it using de Casteljau's algorithm to get a polynomial over
    (approximately) the same region first.

    If we run out of higher-precision polynomials (if we empty out the
    entire stack), then we give up on root refinement for this island.
    The ocean class will notice this, provide the island with a
    higher-precision polynomial, and restart root refinement.  Basically
    the only information kept in that case is the lower and upper bounds
    on the island.  Since these are updated whenever we discover a "half"
    (of an island or a segment) that definitely contains no roots, we
    never need to re-examine these gaps.  (We could keep more information.
    For example, we could keep a record of split points that succeeded
    and failed.  However, a split point that failed at lower precision
    is likely to succeed at higher precision, so it's not worth avoiding.
    It could be useful to select split points that are known to succeed,
    but starting from a new Bernstein polynomial over a slightly different
    region, hitting such split points would require de Casteljau splits
    with non-power-of-two denominators, which are much much slower.)
    """

    def __init__(self, interval_bernstein_polynomial bp, rr_gap lgap, rr_gap rgap):
        """
        Initialize an island from a Bernstein polynomial, and the gaps to
        the left and right of the island.
        """
        self.bp = bp
        self.ancestors = []
        self.target_width = None
        self.lgap = lgap
        self.rgap = rgap
        self.known_done = False

        lgap.risland = self
        rgap.lisland = self

    def __repr__(self):
        """
        Return a short summary of this island.
        """
        return "island of width %s" % RR(self.bp.region_width())

    def _info(self):
        # A python accessor for the information in this island.
        # For debugging purposes.
        return (self.bp, self.ancestors, self.target_width, self.lgap, self.rgap, self.known_done)

    def shrink_bp(self, context ctx):
        """
        If the island's Bernstein polynomial covers a region much larger
        than the island itself (in particular, if either the island's
        left gap or right gap are totally contained in the polynomial's
        region) then shrink the polynomial down to cover the island more
        tightly.
        """
        while True:
            bounds = self.bp.region()
            lbounds = relative_bounds(bounds, self.lgap.region())
            rbounds = relative_bounds(bounds, self.rgap.region())

            rtarget = wordsize_rational(rbounds[0], rbounds[1], ctx.wordsize)
            ltarget = wordsize_rational(lbounds[1], lbounds[0], ctx.wordsize)

            if lbounds[0] > zero_QQ or (ltarget >= QQ_1_32 and ltarget >= one_QQ-rtarget):
                (_, self.bp, _) = self.bp.de_casteljau(ctx, ltarget)
            elif rbounds[1] < one_QQ or rtarget <= QQ_31_32:
                (self.bp, _, _) = self.bp.de_casteljau(ctx, rtarget)
            else:
                break

            self.bp.lsign = self.lgap.sign
            self.bp.usign = self.rgap.sign

    def refine(self, context ctx):
        """
        Attempts to shrink and/or split this island into sub-island
        that each definitely contain exactly one root.
        """
        self.shrink_bp(ctx)
        try:
            self.refine_recurse(ctx, self.bp, self.ancestors, [], True)
        except PrecisionError:
            pass

    def refine_recurse(self, context ctx, interval_bernstein_polynomial bp, ancestors, history, rightmost):
        """
        This implements the root isolation algorithm described in the
        class documentation for class island.  This is the implementation
        of both the whole-island and the left-segment algorithms;
        if the flag rightmost is True, then it is the whole-island algorithm,
        otherwise the left-segment algorithm.

        The precision-reduction stack is (ancestors + [bp]); that is, the
        top-of-stack is maintained separately.
        """
        cdef interval_bernstein_polynomial p1, p2
        cdef rr_gap mgap
        cdef island lisland

        # If you examine the description of this algorithm in the
        # class island class documentation, you see that several of
        # the calls involved are tail recursions (that is, they are
        # of the form "call (this algorithm) recursively, then return").
        # We perform tail-recursion optimization by hand: such calls
        # assign new values to the method parameters, and then fall off
        # the end of the following "while True" loop, and return to
        # the top of the method here.

        # This optimization is VITAL; without it, several test polynomials
        # (in particular, polynomials with roots that are very, very
        # close together) give stack overflow errors.

        while True:
            if rightmost and self.bp_done(bp):
                if bp.max_variations == 0:
                    # No roots!  Make the island disappear.
                    self.lgap.risland = self.rgap.risland
                    self.rgap.risland.lgap = self.lgap
                    self.lgap.upper = self.rgap.upper
                else:
                    self.bp = bp
                    self.ancestors = ancestors
                return

            # This is our heuristic for deciding when to do degree
            # reduction.

            # In general, given a degree-d Bernstein polynomial,
            # if you perform a de Casteljau split at 1/2, the
            # coefficients of the resulting polynomials are about d bits
            # smaller than the coefficients of the original polynomial.

            # Conversely, if you do a split and the coefficients are
            # k bits smaller than the coefficients of the original
            # polynomial, this indicates that the original polynomial
            # may be very close to some degree-k polynomial.

            # Here, we look at the amount that the coefficients got smaller
            # over the last 3 splits.  If they got more than 100 bits
            # smaller (that is, an average of more than 33 bits per split),
            # then we expect that the original polynomial is not
            # close to any polynomial of degree 30 or less.  Since degree
            # reduction works by finding a polynomial of (currently)
            # degree 30 that is close to the original polynomial,
            # then this drastic reduction in coefficient size means that
            # degree reduction is likely to fail, so we don't bother
            # attempting it.

            # Note that this does not take into account the effects of
            # random splitting; maybe it should?  For example, if the
            # last three splits were random splits with split points near
            # 1/4, and we're in the left-hand branch of all these splits,
            # then we would expect twice as much coefficient reduction;
            # so a coefficient size drop of 100 bits would actually
            # mean that we are near a degree-17 polynomial.

            # Also, we are willing to throw away up to 70 + msb_delta
            # bits of precision during degree reduction.  I don't remember
            # how I selected this heuristic...

            if len(history) >= 3:
                old = history[-3]
                old_msb = old.get_msb_bit()
                cur_msb = bp.get_msb_bit()
                msb_delta = old_msb - cur_msb
                if msb_delta <= 100 and bp.bitsize - msb_delta >= 70:
                    bp = bp.down_degree_iter(ctx, bp.scale_log2 + 70 + msb_delta)

            history = history[-3:] + [bp]

            # Heuristically push more lower-precision polynomials
            # on the "ancestors" stack
            (ancestors, bp) = self.less_bits(ancestors, bp)

            # Currently, we try a random split only once before giving up
            # and trying for higher precision; and then we try a 1/2 split
            # at the higher precision.  Maybe we should try more random
            # splits before going to higher precision?  Maybe we should
            # try less 1/2 splits?
            rv = bp.try_split(ctx, [])
            while rv is None:
                if bp.bitsize > 30:
                    rv = bp.try_rand_split(ctx, [])
                if rv is None:
                    (ancestors, bp) = self.more_bits(ctx, ancestors, bp, rightmost)
                    if rv is None:
                        rv = bp.try_split(ctx, [])

            (p1_, p2_, _) = rv
            p1 = p1_
            p2 = p2_

            # ctx.dc_log_append(('split_results', p1.variations(), p2.variations()))

            if p1.variations()[1] == 0:
                self.lgap.upper = p2.lower
                if p2.variations()[1] == 0:
                    if rightmost:
                        # No roots!  Make the island disappear.
                        self.lgap.risland = self.rgap.risland
                        self.rgap.risland.lgap = self.lgap
                        self.lgap.upper = self.rgap.upper
                        return
                    self.lgap.upper = p2.upper
                    return
                else:
                    bp = p2
                    # return to top of function (tail recursion optimization)
            else:
                if p2.variations()[1] == 0:
                    if rightmost:
                        self.rgap.lower = p2.lower
                        bp = p1
                        # return to top of function
                    else:
                        # Split the island!
                        mgap = rr_gap(p2.lower, p2.upper, p2.lsign)
                        lisland = island(p1, self.lgap, mgap)
                        lisland.target_width = self.target_width
                        self.lgap = mgap
                        mgap.risland = self
                        if not lisland.done(ctx):
                            try:
                                lisland.refine_recurse(ctx, p1, ancestors, history, True)
                            except PrecisionError:
                                pass
                        return
                else:
                    self.refine_recurse(ctx, p1, ancestors, history, False)
                    assert self.lgap.upper == p2.lower
                    bp = p2
                    # return to top of function (tail recursion optimization)

    def less_bits(self, ancestors, interval_bernstein_polynomial bp):
        """
        Heuristically push lower-precision polynomials on
        the polynomial stack.  See the class documentation for class
        island for more information.
        """
        # The current heuristic is to always push one lower-precision
        # polynomial, unless the current polynomial is floating-point.
        # If the current polynomial has bitsize < 130, then the new
        # polynomial is floating-point, otherwise it is integer,
        # with half the precision of the original.

        if bp.bitsize < 130:
            if isinstance(bp, interval_bernstein_polynomial_float):
                return (ancestors, bp)
            else:
                return (ancestors + [bp], bp.as_float())
        else:
            return (ancestors + [bp], bp.downscale(bp.bitsize // 2))

    def more_bits(self, context ctx, ancestors, interval_bernstein_polynomial bp, rightmost):
        """
        Find a Bernstein polynomial on the "ancestors" stack with
        more precision than bp; if it is over a different region,
        then shrink its region to (approximately) match that of bp.
        (If this is rightmost -- if bp covers the whole island -- then
        we only require that the new region cover the whole island
        fairly tightly; if this is not rightmost, then the new region
        will have exactly the same right boundary as bp, although the
        left boundary may vary slightly.)
        """
        cur_msb = bp.scale_log2 + bp.bitsize
        extra_bits = bp.bitsize // 2
        if extra_bits < 30:
            extra_bits = 30

        target_lsb_h = cur_msb - 3*extra_bits
        target_lsb = cur_msb - 4*extra_bits
        target_lsb_l = cur_msb - 6*extra_bits

        if bp.bitsize < 32:
            target_lsb_h = cur_msb - 48
            target_lsb = target_lsb_h - 16
            target_lsb_l = target_lsb - 32

        cdef interval_bernstein_polynomial anc
        cdef interval_bernstein_polynomial ancestor_val

        for i in range(len(ancestors)-1, -1, -1):
            anc = ancestors[i]
            if target_lsb_h >= anc.scale_log2:
                ancestor_val = anc
                first_lsb = ancestor_val.scale_log2
                first_msb = first_lsb + ancestor_val.bitsize

                ancestors = ancestors[:i]

                if bp.region() == ancestor_val.region():
                    if bp.bitsize < 32:
                        return (ancestors + [ancestor_val], ancestor_val.as_float())
                    else:
                        return (ancestors, ancestor_val)

                new_lsb = ancestor_val.scale_log2
                if new_lsb < target_lsb_l:
                    new_lsb = target_lsb

                hv_width = ancestor_val.region_width()
                rel_bounds = relative_bounds(ancestor_val.region(), bp.region())
                rel_width = rel_bounds[1] - rel_bounds[0]

                rel_width_rr = RR(rel_width)

                ctx.dc_log_append(('pulling',
                                   first_msb,
                                   ancestor_val.level,
                                   first_lsb,
                                   ancestor_val.scale_log2,
                                   rel_width_rr,
                                   new_lsb,
                                   cur_msb, bp.scale_log2,
                                   target_lsb_h, target_lsb, target_lsb_l))

                if rightmost:
                    maybe_rgap = self.rgap
                else:
                    maybe_rgap = None
                    if rel_bounds[1] < 1:
                        (ancestor_val, _, _) = ancestor_val.de_casteljau(ctx, rel_bounds[1])
                        ctx.dc_log_append(('pull_right', rel_bounds[1]))
                        if ancestor_val.region_width() / hv_width < ~Integer(32):
                            ancestor_val = ancestor_val.down_degree_iter(ctx, target_lsb_h)

                    rel_bounds = relative_bounds(ancestor_val.region(), bp.region())
                    assert rel_bounds[1] == 1

                ancestor_val = split_for_targets(ctx, ancestor_val, [(self.lgap, maybe_rgap, target_lsb_h)])[0]
#                 if rel_lbounds[1] > 0:
#                     left_split = -exact_rational(simple_wordsize_float(-rel_lbounds[1], -rel_lbounds[0]))
#                     (_, ancestor_val, _) = ancestor_val.de_casteljau(ctx, left_split)
#                     ctx.dc_log_append(('pull_left', left_split))

                ancestor_val.lsign = bp.lsign
                ancestor_val.usign = bp.usign

                new_rel_bounds = relative_bounds(ancestor_val.region(), bp.region())
                assert new_rel_bounds[1] - new_rel_bounds[0] >= Integer(255)/256

                while ancestor_val.scale_log2 < target_lsb_l:
                    ancestors = ancestors + [ancestor_val]
                    ancestor_val = ancestor_val.downscale(ancestor_val.bitsize // 2)

                if bp.bitsize < 32:
                    return (ancestors + [ancestor_val], ancestor_val.as_float())

                return (ancestors, ancestor_val)

        self.ancestors = []
        raise PrecisionError()

    def reset_root_width(self, target_width):
        """
        Modify the criteria for this island to require that it is not "done"
        until its width is less than or equal to target_width.
        """

        width = self.bp.upper - self.bp.lower

        if target_width < width:
            self.known_done = False

        if self.target_width is None or target_width < self.target_width:
            self.target_width = target_width

    def bp_done(self, interval_bernstein_polynomial bp):
        """
        Examine the given Bernstein polynomial to see if it is known
        to have exactly one root in its region.  (In addition, we require
        that the polynomial region not include 0 or 1.  This makes things
        work if the user gives explicit bounds to real_roots(),
        where the lower or upper bound is a root of the polynomial.
        real_roots() deals with this by explicitly detecting it,
        dividing out the appropriate linear polynomial, and adding the
        root to the returned list of roots; but then if the island
        considers itself "done" with a region including 0 or 1, the returned
        root regions can overlap with each other.)
        """

        variations = bp.variations()[1]

        if variations > 1:
            return False
        if bp.lower == 0:
            return False
        if bp.upper == 1:
            return False
        if variations == 0:
            return True
        if self.target_width is not None and self.bp.upper - self.bp.lower > self.target_width:
            return False
        if bp.level == 0:
            return True
        if not (0 in bp.slope_err + bp.slope_range()):
            return True
        return False

    def done(self, context ctx):
        """
        Check to see if the island is known to contain zero roots or
        is known to contain one root.
        """

        if self.known_done:
            return True

        if self.bp_done(self.bp):
            self.known_done = True
        else:
            self.shrink_bp(ctx)
            if self.bp_done(self.bp):
                self.known_done = True

        return self.known_done

    def has_root(self):
        """
        Assuming that the island is done (has either 0 or 1 roots),
        reports whether the island has a root.
        """
        assert self.known_done

        return bool(self.bp.max_variations)


cdef class rr_gap:
    """
    A simple class representing the gaps between islands, in my
    ocean-island root isolation algorithm.  Named "rr_gap" for
    "real roots gap", because "gap" seemed too short and generic.
    """

    def __init__(self, lower, upper, sign):
        """
        Initialize an rr_gap element.
        """
        self.lower = lower
        self.upper = upper
        self.sign = sign

    def region(self):
        return (self.lower, self.upper)


class linear_map:
    """
    A simple class to map linearly between original coordinates
    (ranging from [lower .. upper]) and ocean coordinates (ranging
    from [0 .. 1]).
    """

    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper
        self.width = upper - lower

    def from_ocean(self, region):
        (l, u) = region
        return (self.lower + l*self.width, self.lower + u*self.width)

    def to_ocean(self, region):
        (l, u) = region
        return ((l - self.lower) / self.width, (u - self.lower) / self.width)


lmap = linear_map(0, 1)


class warp_map:
    """
    A class to map between original coordinates and ocean coordinates.
    If neg is False, then the original->ocean transform is
    x -> x/(x+1), and the ocean->original transform is x/(1-x);
    this maps between [0 .. infinity] and [0 .. 1].
    If neg is True, then the original->ocean transform is
    x -> -x/(1-x), and the ocean->original transform is the same thing:
    -x/(1-x).  This maps between [0 .. -infinity] and [0 .. 1].
    """

    def __init__(self, neg):
        self.neg = neg

    def from_ocean(self, region):
        (l, u) = region
        if self.neg:
            return (-u/(1-u), -l/(1-l))
        else:
            return (l/(1-l), u/(1-u))

    def to_ocean(self, region):
        (l, u) = region
        if self.neg:
            return (-u/(1-u), -l/(1-l))
        else:
            return (l/(l+1), u/(u+1))


def real_roots(p, bounds=None, seed=None, skip_squarefree=False, do_logging=False, wordsize=32, retval='rational', strategy=None, max_diameter=None):
    """
    Compute the real roots of a given polynomial with exact coefficients
    (integer, rational, and algebraic real coefficients are supported).

    This returns a list of pairs of a root and its multiplicity.

    The root itself can be returned in one of three different ways.
    If retval=='rational', then it is returned as a pair of rationals
    that define a region that includes exactly one root.  If
    retval=='interval', then it is returned as a RealIntervalFieldElement
    that includes exactly one root.  If retval=='algebraic_real', then
    it is returned as an AlgebraicReal.  In the former two cases, all
    the intervals are disjoint.

    An alternate high-level algorithm can be used by selecting
    strategy='warp'.  This affects the conversion into Bernstein
    polynomial form, but still uses the same ocean-island algorithm
    as the default algorithm.  The 'warp' algorithm performs the conversion
    into Bernstein polynomial form much more quickly, but performs
    the rest of the computation slightly slower in some benchmarks.
    The 'warp' algorithm is particularly likely to be helpful for
    low-degree polynomials.

    Part of the algorithm is randomized; the seed parameter gives a seed
    for the random number generator.  (By default, the same
    seed is used for every call, so that results are repeatable.)  The
    random seed may affect the running time, or the exact intervals returned,
    but the results are correct regardless of the seed used.

    The bounds parameter lets you find roots in some proper subinterval of
    the reals; it takes a pair of a rational lower and upper bound
    and only roots within this bound will be found.  Currently, specifying
    bounds does not work if you select strategy='warp', or if you
    use a polynomial with algebraic real coefficients.

    By default, the algorithm will do a squarefree decomposition
    to get squarefree polynomials.  The skip_squarefree parameter
    lets you skip this step.  (If this step is skipped, and the polynomial
    has a repeated real root, then the algorithm will loop forever!
    However, repeated non-real roots are not a problem.)

    For integer and rational coefficients, the squarefree
    decomposition is very fast, but it may be slow for algebraic
    reals.  (It may trigger exact computation, so it might be
    arbitrarily slow.  The only other way that this algorithm might
    trigger exact computation on algebraic real coefficients is that
    it checks the constant term of the input polynomial for equality with
    zero.)

    Part of the algorithm works (approximately) by splitting numbers into
    word-size pieces (that is, pieces that fit into a machine word).
    For portability, this defaults to always selecting pieces suitable
    for a 32-bit machine; the wordsize parameter lets you make choices
    suitable for a 64-bit machine instead.  (This affects the running
    time, and the exact intervals returned, but the results are correct
    on both 32- and 64-bit machines even if the wordsize is chosen "wrong".)

    The precision of the results can be improved (at the expense of time,
    of course) by specifying the max_diameter parameter.  If specified,
    this sets the maximum diameter() of the intervals returned.
    (Sage defines diameter() to be the relative diameter for intervals
    that do not contain 0, and the absolute diameter for intervals
    containing 0.)  This directly affects the results in rational or
    interval return mode; in algebraic_real mode, it increases the
    precision of the intervals passed to the algebraic number package,
    which may speed up some operations on that algebraic real.

    Some logging can be enabled with do_logging=True.  If logging is enabled,
    then the normal values are not returned; instead, a pair of
    the internal context object and a list of all the roots in their
    internal form is returned.

    ALGORITHM: We convert the polynomial into the Bernstein basis, and
    then use de Casteljau's algorithm and Descartes' rule of signs
    (using interval arithmetic) to locate the roots.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: real_roots(x^3 - x^2 - x - 1)
        [((7/4, 19/8), 1)]
        sage: real_roots((x-1)*(x-2)*(x-3)*(x-5)*(x-8)*(x-13)*(x-21)*(x-34))
        [((11/16, 33/32), 1), ((11/8, 33/16), 1), ((11/4, 55/16), 1), ((77/16, 165/32), 1), ((11/2, 33/4), 1), ((11, 55/4), 1), ((165/8, 341/16), 1), ((22, 44), 1)]
        sage: real_roots(x^5 * (x^2 - 9999)^2 - 1)
        [((-29274496381311/9007199254740992, 419601125186091/2251799813685248), 1), ((2126658450145849453951061654415153249597/21267647932558653966460912964485513216, 4253316902721330018853696359533061621799/42535295865117307932921825928971026432), 1), ((1063329226287740282451317352558954186101/10633823966279326983230456482242756608, 531664614358685696701445201630854654353/5316911983139663491615228241121378304), 1)]
        sage: real_roots(x^5 * (x^2 - 9999)^2 - 1, seed=42)
        [((-123196838480289/18014398509481984, 293964743458749/9007199254740992), 1), ((8307259573979551907841696381986376143/83076749736557242056487941267521536, 16614519150981033789137940378745325503/166153499473114484112975882535043072), 1), ((519203723562592617581015249797434335/5192296858534827628530496329220096, 60443268924081068060312183/604462909807314587353088), 1)]
        sage: real_roots(x^5 * (x^2 - 9999)^2 - 1, wordsize=64)
        [((-62866503803202151050003/19342813113834066795298816, 901086554512564177624143/4835703278458516698824704), 1), ((544424563237337315214990987922809050101157/5444517870735015415413993718908291383296, 1088849127096660194637118845654929064385439/10889035741470030830827987437816582766592), 1), ((272212281929661439711063928866060007142141/2722258935367507707706996859454145691648, 136106141275823501959100399337685485662633/1361129467683753853853498429727072845824), 1)]
        sage: real_roots(x)
        [((-47/256, 81/512), 1)]
        sage: real_roots(x * (x-1))
        [((-47/256, 81/512), 1), ((1/2, 1201/1024), 1)]
        sage: real_roots(x-1)
        [((209/256, 593/512), 1)]
        sage: real_roots(x*(x-1)*(x-2), bounds=(0, 2))
        [((0, 0), 1), ((81/128, 337/256), 1), ((2, 2), 1)]
        sage: real_roots(x*(x-1)*(x-2), bounds=(0, 2), retval='algebraic_real')
        [(0, 1), (1, 1), (2, 1)]
        sage: v = 2^40
        sage: real_roots((x^2-1)^2 * (x^2 - (v+1)/v))
        [((-12855504354077768210885019021174120740504020581912910106032833/12855504354071922204335696738729300820177623950262342682411008, -6427752177038884105442509510587059395588605840418680645585479/6427752177035961102167848369364650410088811975131171341205504), 1), ((-1125899906842725/1125899906842624, -562949953421275/562949953421312), 2), ((62165404551223330269422781018352603934643403586760330761772204409982940218804935733653/62165404551223330269422781018352605012557018849668464680057997111644937126566671941632, 3885337784451458141838923813647037871787041539340705594199885610069035709862106085785/3885337784451458141838923813647037813284813678104279042503624819477808570410416996352), 2), ((509258994083853105745586001837045839749063767798922046787130823804169826426726965449697819/509258994083621521567111422102344540262867098416484062659035112338595324940834176545849344, 25711008708155536421770038042348240136257704305733983563630791/25711008708143844408671393477458601640355247900524685364822016), 1)]
        sage: real_roots(x^2 - 2)
        [((-3/2, -1), 1), ((1, 3/2), 1)]
        sage: real_roots(x^2 - 2, retval='interval')
        [(-2.?, 1), (2.?, 1)]
        sage: real_roots(x^2 - 2, max_diameter=1/2^30)
        [((-22506280506048041472675379598886543645348790970912519198456805737131269246430553365310109/15914343565113172548972231940698266883214596825515126958094847260581103904401068017057792, -45012561012096082945350759197773086524448972309421182031053065599548946985601579935498343/31828687130226345097944463881396533766429193651030253916189694521162207808802136034115584), 1), ((45012561012096082945350759197773086524448972309421182031053065599548946985601579935498343/31828687130226345097944463881396533766429193651030253916189694521162207808802136034115584, 22506280506048041472675379598886543645348790970912519198456805737131269246430553365310109/15914343565113172548972231940698266883214596825515126958094847260581103904401068017057792), 1)]
        sage: real_roots(x^2 - 2, retval='interval', max_diameter=1/2^500)
        [(-1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641572735013846230912297024924836055850737212644121497099935831413222665927505592755799950501152782060571470109559971605970274534596862014728517418640889198609552?, 1), (1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641572735013846230912297024924836055850737212644121497099935831413222665927505592755799950501152782060571470109559971605970274534596862014728517418640889198609552?, 1)]
        sage: ar_rts = real_roots(x^2 - 2, retval='algebraic_real'); ar_rts
        [(-1.414213562373095?, 1), (1.414213562373095?, 1)]
        sage: ar_rts[0][0]^2 - 2 == 0
        True
        sage: v = 2^40
        sage: real_roots((x-1) * (x-(v+1)/v), retval='interval')
        [(1.000000000000?, 1), (1.000000000001?, 1)]
        sage: v = 2^60
        sage: real_roots((x-1) * (x-(v+1)/v), retval='interval')
        [(1.000000000000000000?, 1), (1.000000000000000001?, 1)]
        sage: real_roots((x-1) * (x-2), strategy='warp')
        [((499/525, 1173/875), 1), ((337/175, 849/175), 1)]
        sage: real_roots((x+3)*(x+1)*x*(x-1)*(x-2), strategy='warp')
        [((-1713/335, -689/335), 1), ((-2067/2029, -689/1359), 1), ((0, 0), 1), ((499/525, 1173/875), 1), ((337/175, 849/175), 1)]
        sage: real_roots((x+3)*(x+1)*x*(x-1)*(x-2), strategy='warp', retval='algebraic_real')
        [(-3.000000000000000?, 1), (-1.000000000000000?, 1), (0, 1), (1.000000000000000?, 1), (2.000000000000000?, 1)]
        sage: ar_rts = real_roots(x-1, retval='algebraic_real')
        sage: ar_rts[0][0] == 1
        True

    If the polynomial has no real roots, we get an empty list.

    ::

        sage: (x^2 + 1).real_root_intervals()
        []

    We can compute Conway's constant
    (see http://mathworld.wolfram.com/ConwaysConstant.html) to arbitrary
    precision. ::

        sage: p = x^71 - x^69 - 2*x^68 - x^67 + 2*x^66 + 2*x^65 + x^64 - x^63 - x^62 - x^61 - x^60 - x^59 + 2*x^58 + 5*x^57 + 3*x^56 - 2*x^55 - 10*x^54 - 3*x^53 - 2*x^52 + 6*x^51 + 6*x^50 + x^49 + 9*x^48 - 3*x^47 - 7*x^46 - 8*x^45 - 8*x^44 + 10*x^43 + 6*x^42 + 8*x^41 - 5*x^40 - 12*x^39 + 7*x^38 - 7*x^37 + 7*x^36 + x^35 - 3*x^34 + 10*x^33 + x^32 - 6*x^31 - 2*x^30 - 10*x^29 - 3*x^28 + 2*x^27 + 9*x^26 - 3*x^25 + 14*x^24 - 8*x^23 - 7*x^21 + 9*x^20 + 3*x^19 - 4*x^18 - 10*x^17 - 7*x^16 + 12*x^15 + 7*x^14 + 2*x^13 - 12*x^12 - 4*x^11 - 2*x^10 + 5*x^9 + x^7 - 7*x^6 + 7*x^5 - 4*x^4 + 12*x^3 - 6*x^2 + 3*x - 6
        sage: cc = real_roots(p, retval='algebraic_real')[2][0] # long time
        sage: RealField(180)(cc)                                # long time
        1.3035772690342963912570991121525518907307025046594049

    Now we play with algebraic real coefficients. ::

        sage: x = polygen(AA)
        sage: p = (x - 1) * (x - sqrt(AA(2))) * (x - 2)
        sage: real_roots(p)
        [((499/525, 2171/1925), 1), ((1173/875, 2521/1575), 1), ((337/175, 849/175), 1)]
        sage: ar_rts = real_roots(p, retval='algebraic_real'); ar_rts
        [(1.000000000000000?, 1), (1.414213562373095?, 1), (2.000000000000000?, 1)]
        sage: ar_rts[1][0]^2 == 2
        True
        sage: ar_rts = real_roots(x*(x-1), retval='algebraic_real')
        sage: ar_rts[0][0] == 0
        True
        sage: p2 = p * (p - 1/100); p2
        x^6 - 8.82842712474619?*x^5 + 31.97056274847714?*x^4 - 60.77955262170047?*x^3 + 63.98526763257801?*x^2 - 35.37613490585595?*x + 8.028284271247462?
        sage: real_roots(p2, retval='interval')
        [(1.00?, 1), (1.1?, 1), (1.38?, 1), (1.5?, 1), (2.00?, 1), (2.1?, 1)]
        sage: p = (x - 1) * (x - sqrt(AA(2)))^2 * (x - 2)^3 * sqrt(AA(3))
        sage: real_roots(p, retval='interval')
        [(1.000000000000000?, 1), (1.414213562373095?, 2), (2.000000000000000?, 3)]

    TESTS:

    Check that :issue:`20269` is fixed::

        sage: real_roots(polygen(AA))[0][0][0].parent()
        Rational Field

    Check that :issue:`10803` is fixed::

        sage: f = 2503841067*x^13 - 15465014877*x^12 + 37514382885*x^11 - 44333754994*x^10 + 24138665092*x^9 - 2059014842*x^8 - 3197810701*x^7 + 803983752*x^6 + 123767204*x^5 - 26596986*x^4 - 2327140*x^3 + 75923*x^2 + 7174*x + 102
        sage: len(real_roots(f,seed=1))
        13
    """
    base = p.base_ring()

    ar_input = False

    if base is not ZZ:
        ZZx = PolynomialRing(ZZ, 'x')
        if base is QQ:
            p = ZZx(p * p.denominator())
        elif base is AA:
            ar_input = True
        else:
            raise ValueError("Don't know how to isolate roots for " + str(p.parent()))

    if ar_input and bounds is not None:
        raise NotImplementedError("Cannot set your own bounds with algebraic real input")

    if ar_input:
        strategy = 'warp'

    if bounds is not None and strategy == 'warp':
        raise NotImplementedError("Cannot set your own bounds with strategy=warp")

    if seed is None:
        seed = 1

    if skip_squarefree:
        factors = [(p, 1)]
    else:
        factors = p.squarefree_decomposition()

    if max_diameter is not None:
        max_diameter = QQ(max_diameter)

    ctx = context(do_logging, seed, wordsize)

    extra_roots = []
    oceans = []

    cdef ocean oc

    for factor, exp in factors:
        if strategy=='warp':
            if factor.constant_coefficient() == 0:
                x = factor.parent().gen()
                extra_roots.append(((QQ.zero(), QQ.zero()), x, exp, None, None))
                factor = factor // x
            if ar_input:
                oc = ocean(ctx, bernstein_polynomial_factory_ar(factor, False), warp_map(False))
            else:
                b = to_bernstein_warp(factor)
                oc = ocean(ctx, bernstein_polynomial_factory_ratlist(b), warp_map(False))
            oc.find_roots()
            oceans.append((oc, factor, exp))

            if ar_input:
                oc = ocean(ctx, bernstein_polynomial_factory_ar(factor, True), warp_map(True))
            else:
                b = copy(b)
                for i in range(1, len(b), 2):
                    b[i] = -b[i]
                oc = ocean(ctx, bernstein_polynomial_factory_ratlist(b), warp_map(True))
            oc.find_roots()
            oceans.append((oc, factor, exp))
        else:
            if bounds is None:
                fac_bounds = rational_root_bounds(factor)
                if fac_bounds is None:
                    continue
                else:
                    (left, right) = fac_bounds
            else:
                (left, right) = bounds
                # Bad things happen if the bounds are roots themselves.
                # Avoid this by dividing out linear polynomials if
                # the bounds are roots.
                if factor(left) == 0:
                    x = factor.parent().gen()
                    linfac = (x * left.denominator() - left.numerator())
                    extra_roots.append(((left, left), linfac, exp, None, None))
                    factor = factor // linfac
                if factor(right) == 0:
                    x = factor.parent().gen()
                    linfac = (x * right.denominator() - right.numerator())
                    extra_roots.append(((right, right), linfac, exp, None, None))
                    factor = factor // linfac

            b, _ = to_bernstein(factor, left, right)

            oc = ocean(ctx, bernstein_polynomial_factory_ratlist(b), linear_map(left, right))
            oc.find_roots()
            oceans.append((oc, factor, exp))

    while True:
        all_roots = copy(extra_roots)

        for oc, factor, exp in oceans:
            rel_roots = oc.roots()

            cur_roots = [oc.mapping.from_ocean(r) for r in rel_roots]

            all_roots.extend([(cur_roots[j], factor, exp, oc, j) for j in range(len(cur_roots))])

        all_roots.sort()

        ok = True

        target_widths = [None] * len(all_roots)

        if max_diameter is not None:
            # Check to make sure that no intervals are too wide.

            # We use half_diameter, because if we ended up with a rational
            # interval that was exactly max_diameter, we might not
            # be able to coerce it into an interval small enough.
            half_diameter = max_diameter/2

            for i in range(len(all_roots)):
                root = all_roots[i][0]
                if (root[0] <= 0) and (root[1] >= 0):
                    cur_diam = root[1] - root[0]
                    if cur_diam > half_diameter:
                        target_widths[i] = half_diameter
                        ok = False
                else:
                    cur_diam = (root[1] - root[0]) / abs((root[0] + root[1]) / 2)
                    if cur_diam > half_diameter:
                        target_widths[i] = (root[1] - root[0]) / cur_diam * half_diameter
                        ok = False

        for i in range(len(all_roots) - 1):
            # Check to be sure that all intervals are disjoint.
            if all_roots[i][0][1] >= all_roots[i+1][0][0]:
                ok = False
                cur_width = max(all_roots[i+1][0][1] - all_roots[i+1][0][0], all_roots[i][0][1] - all_roots[i][0][0])
                target_width = cur_width/16
                target_widths[i] = target_width
                target_widths[i+1] = target_width

        for i in range(len(all_roots)):
            if target_widths[i] is not None:
                root = all_roots[i][0]
                oc = all_roots[i][3]
                target_region = (root[0], root[0] + target_widths[i])
                if target_region[0] <= 0 and target_region[1] >= 0:
                    target_region = (root[1] - target_widths[i], root[1])

                ocean_target = oc.mapping.to_ocean(target_region)
                oc.reset_root_width(all_roots[i][4], ocean_target[1] - ocean_target[0])

        if ok:
            break

        for oc, factor, exp in oceans:
            oc.find_roots()

    if do_logging:
        return ctx, all_roots

    if retval == 'rational':
        return [(r[0], r[2]) for r in all_roots]

    for i in range(1000):
        intv_bits = 53 << i
        intv_fld = RealIntervalField(intv_bits)
        intv_roots = [(intv_fld(r[0]), r[1], r[2], r[3], r[4]) for r in all_roots]
        ok = True

        if max_diameter is not None:
            for rt in intv_roots:
                if rt[0].diameter() > max_diameter:
                    ok = False

        for j in range(len(intv_roots) - 1):
            # The following line should work, but does not due to a Cython
            # bug (it calls PyObject_Cmp() for comparison operators,
            # instead of PyObject_RichCompare()).

            # if not (intv_roots[j][0] < intv_roots[j+1][0]):

            if not (intv_roots[j][0].upper() < intv_roots[j+1][0].lower()):
                ok = False

        if ok:
            break

    if retval=='interval':
        return [(r[0], r[2]) for r in intv_roots]

    if retval=='algebraic_real':
        return [(AA.polynomial_root(r[1], r[0]), r[2]) for r in intv_roots]

    raise ValueError("Illegal retval parameter " + retval)


def scale_intvec_var(Vector_integer_dense c, k):
    """
    Given a vector of integers c of length n+1, and a rational
    k == kn / kd, multiplies each element c[i] by (kd^i)*(kn^(n-i)).

    Modifies the input vector; has no return value.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: v = vector(ZZ, [1, 1, 1, 1])
        sage: scale_intvec_var(v, 3/4)
        sage: v
        (64, 48, 36, 27)
    """

    kn = numerator(k)
    kd = denominator(k)

    # XXX could use direct mpz calls for more speed
    cdef Integer factor = Integer(kd) ** (len(c) - 1)

    cdef int i
    for i from 0 <= i < len(c):
        c[i] = c[i] * factor
        factor = (factor * kn) // kd


def taylor_shift1_intvec(Vector_integer_dense c):
    """
    Given a vector of integers c of length d+1, representing the
    coefficients of a degree-d polynomial p, modify the vector
    to perform a Taylor shift by 1 (that is, p becomes p(x+1)).

    This is the straightforward algorithm, which is not asymptotically
    optimal.

    Modifies the input vector; has no return value.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: p = (x-1)*(x-2)*(x-3)
        sage: v = vector(ZZ, p.list())
        sage: p, v
        (x^3 - 6*x^2 + 11*x - 6, (-6, 11, -6, 1))
        sage: taylor_shift1_intvec(v)
        sage: p(x+1), v
        (x^3 - 3*x^2 + 2*x, (0, 2, -3, 1))
    """
    cdef int degree = len(c) - 1

    cdef int i, k
    for i from 1 <= i <= degree:
        for k from degree-i <= k < degree:
            mpz_add(c._entries[k], c._entries[k], c._entries[k+1])


def reverse_intvec(Vector_integer_dense c):
    """
    Given a vector of integers, reverse the vector (like the reverse()
    method on lists).

    Modifies the input vector; has no return value.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: v = vector(ZZ, [1, 2, 3, 4]); v
        (1, 2, 3, 4)
        sage: reverse_intvec(v)
        sage: v
        (4, 3, 2, 1)
    """
    cdef int i
    cdef int c_len = len(c)
    for i from 0 <= i < c_len // 2:
        mpz_swap(c._entries[i], c._entries[c_len - 1 - i])


realfield_rndu_cache = {}


def get_realfield_rndu(n):
    """
    A simple cache for RealField fields (with rounding set to
    round-to-positive-infinity).

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: get_realfield_rndu(20)
        Real Field with 20 bits of precision and rounding RNDU
        sage: get_realfield_rndu(53)
        Real Field with 53 bits of precision and rounding RNDU
        sage: get_realfield_rndu(20)
        Real Field with 20 bits of precision and rounding RNDU
    """
    try:
        return realfield_rndu_cache[n]
    except KeyError:
        fld = RealField(n, rnd='RNDU')
        realfield_rndu_cache[n] = fld
        return fld


cdef class context:
    """
    A simple context class, which is passed through parts of the
    real root isolation algorithm to avoid global variables.

    Holds logging information, a random number generator, and
    the target machine wordsize.
    """

    def __init__(self, do_logging, seed, wordsize):
        """
        Initialize a context class.
        """
        self.seed = seed  # saved to make context printable
        self.random = randstate().python_random(seed=seed)
        self.do_logging = do_logging
        self.wordsize = wordsize
        self.dc_log = []
        self.be_log = []

    def __repr__(self):
        """
        Return a short summary of this context.

        EXAMPLES::

            sage: from sage.rings.polynomial.real_roots import *
            sage: mk_context()
            root isolation context: seed=0
            sage: mk_context(do_logging=True, seed=37, wordsize=64)
            root isolation context: seed=37; do_logging=True; wordsize=64
        """

        s = "root isolation context: seed=%d" % self.seed
        if self.do_logging:
            s = s + "; do_logging=True"
        if self.wordsize != 32:
            s = s + "; wordsize=%d" % self.wordsize
        return s

    cdef void dc_log_append(self, x) noexcept:
        """
        Optional logging for the root isolation algorithm.
        """
        if self.do_logging:
            self.dc_log.append(x)

    cdef void be_log_append(self, x) noexcept:
        """
        Optional logging for degree reduction in the root isolation algorithm.
        """
        if self.do_logging:
            self.be_log.append(x)

    def get_dc_log(self):
        return self.dc_log

    def get_be_log(self):
        return self.be_log


def mk_context(do_logging=False, seed=0, wordsize=32):
    """
    A simple wrapper for creating context objects with coercions,
    defaults, etc.

    For use in doctests.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: mk_context(do_logging=True, seed=3, wordsize=64)
        root isolation context: seed=3; do_logging=True; wordsize=64
    """
    return context(do_logging, seed, wordsize)


def to_bernstein(p, low=0, high=1, degree=None):
    """
    Given a polynomial p with integer coefficients, and rational
    bounds low and high, compute the exact rational Bernstein
    coefficients of p over the region [low .. high].  The optional
    parameter degree can be used to give a formal degree higher than
    the actual degree.

    The return value is a pair (c, scale); c represents the same
    polynomial as p*scale.  (If you only care about the roots of
    the polynomial, then of course scale can be ignored.)

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: to_bernstein(x)
        ([0, 1], 1)
        sage: to_bernstein(x, degree=5)
        ([0, 1/5, 2/5, 3/5, 4/5, 1], 1)
        sage: to_bernstein(x^3 + x^2 - x - 1, low=-3, high=3)
        ([-16, 24, -32, 32], 1)
        sage: to_bernstein(x^3 + x^2 - x - 1, low=3, high=22/7)
        ([296352, 310464, 325206, 340605], 9261)
    """
    if degree is None:
        degree = p.degree()
    elif degree < p.degree():
        raise ValueError('Bernstein degree must be at least polynomial degree')
    vs = ZZ ** (degree + 1)
    c = vs(0)
    for i in range(p.degree() + 1):
        c[i] = p[i]
    scale = ZZ(1)
    if low == 0:
        scale_intvec_var(c, high)
        scale = denominator(high) ** degree
    else:
        scale_intvec_var(c, low)
        scale = denominator(low) ** degree
        taylor_shift1_intvec(c)
        scv = QQ(high - low) / low
        scale_intvec_var(c, scv)
        scale = scale * denominator(scv) ** degree
    reverse_intvec(c)
    taylor_shift1_intvec(c)
    reverse_intvec(c)
    return ([c[k] / binomial(degree, k) for k in range(degree + 1)], scale)


def to_bernstein_warp(p):
    """
    Given a polynomial p with rational coefficients, compute the
    exact rational Bernstein coefficients of p(x/(x+1)).

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: x = polygen(ZZ)
        sage: to_bernstein_warp(1 + x + x^2 + x^3 + x^4 + x^5)
        [1, 1/5, 1/10, 1/10, 1/5, 1]
    """
    c = p.list()
    for i in range(len(c)):
        c[i] = c[i] / binomial(len(c) - 1, i)
    return c


def bernstein_expand(Vector_integer_dense c, int d2):
    """
    Given an integer vector representing a Bernstein polynomial p, and
    a degree d2, compute the representation of p as a Bernstein
    polynomial of formal degree d2.

    This is similar to multiplying by the result of bernstein_up, but
    should be faster for large d2 (this has about the same number of
    multiplies, but in this version all the multiplies are by single
    machine words).

    This returns a pair consisting of the expanded polynomial, and the maximum
    error E.  (So if an element of the returned polynomial is a, and the
    true value of that coefficient is b, then a <= b < a + E.)

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: c = vector(ZZ, [1000, 2000, -3000])
        sage: bernstein_expand(c, 3)
        ((1000, 1666, 333, -3000), 1)
        sage: bernstein_expand(c, 4)
        ((1000, 1500, 1000, -500, -3000), 1)
        sage: bernstein_expand(c, 20)
        ((1000, 1100, 1168, 1205, 1210, 1184, 1126, 1036, 915, 763, 578, 363, 115, -164, -474, -816, -1190, -1595, -2032, -2500, -3000), 1)
    """
    cdef int d1 = len(c)-1

    vs = FreeModule(ZZ, d2+1)

    cdef Vector_integer_dense c2 = vs(0)

    cdef int i, j

    cdef int ndivides = 0

    cdef mpz_t tmp
    cdef mpz_t divisor

    mpz_init(tmp)
    mpz_init_set_ui(divisor, 1)

    # XXX do experimentation here on how to decide when to divide
    cdef int max_bits = max_bitsize_intvec(c) / 2
    if max_bits < 64:
        max_bits = 64

    for i from 0 <= i <= d1:
        mpz_set(c2._entries[i], c._entries[i])

    for i from d1 <= i < d2:
        for j from i >= j >= 0:
            mpz_addmul_ui(c2._entries[j+1], c2._entries[j], j+1)

            mpz_mul_ui(c2._entries[j], c2._entries[j], i+1-j)

        mpz_mul_ui(divisor, divisor, i+1)

        if i == d2-1 or mpz_sizeinbase(divisor, 2) > max_bits:
            for j from 0 <= j <= i+1:
                mpz_fdiv_q(c2._entries[j], c2._entries[j], divisor)
            mpz_set_ui(divisor, 1)
            ndivides = ndivides + 1

    mpz_clear(tmp)
    mpz_clear(divisor)

    return (c2, ndivides)


cdef int max_bitsize_intvec(Vector_integer_dense b) noexcept:
    """
    Given an integer vector, find the approximate log2 of the maximum
    of the absolute values of the elements.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: max_bitsize_intvec_doctest(vector(ZZ, [1, 2, 3, 1024]))
        11
    """
    cdef int max_bits = 0
    cdef int i
    cdef int size

    for i in range(len(b)):
        size = mpz_sizeinbase(b._entries[i], 2)
        if size > max_bits:
            max_bits = size

    return max_bits


def max_bitsize_intvec_doctest(b):
    return max_bitsize_intvec(b)


def dprod_imatrow_vec(Matrix_integer_dense m, Vector_integer_dense v, int k):
    """
    Compute the dot product of row k of the matrix m with the vector v
    (that is, compute one element of the product m*v).

    If v has more elements than m has columns, then elements of v are
    selected using subsample_vec.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: m = matrix(3, range(9))
        sage: dprod_imatrow_vec(m, vector(ZZ, [1, 0, 0, 0]), 1)
        0
        sage: dprod_imatrow_vec(m, vector(ZZ, [0, 1, 0, 0]), 1)
        3
        sage: dprod_imatrow_vec(m, vector(ZZ, [0, 0, 1, 0]), 1)
        4
        sage: dprod_imatrow_vec(m, vector(ZZ, [0, 0, 0, 1]), 1)
        5
        sage: dprod_imatrow_vec(m, vector(ZZ, [1, 0, 0]), 1)
        3
        sage: dprod_imatrow_vec(m, vector(ZZ, [0, 1, 0]), 1)
        4
        sage: dprod_imatrow_vec(m, vector(ZZ, [0, 0, 1]), 1)
        5
        sage: dprod_imatrow_vec(m, vector(ZZ, [1, 2, 3]), 1)
        26
    """
    assert 0 <= k < m.nrows()
    assert m.ncols() <= len(v)

    cdef Integer sum = Integer(0)

    cdef mpz_t tmp

    cdef int msize = m.ncols()
    cdef int vsize = len(v)
    cdef int ra
    cdef int a
    mpz_init(tmp)
    for a in range(msize):
        ra = subsample_vec(a, msize, vsize)
        m.get_unsafe_mpz(k, a, tmp)
        mpz_addmul(sum.value, tmp, v._entries[ra])
    mpz_clear(tmp)
    return sum


def min_max_delta_intvec(Vector_integer_dense a, Vector_integer_dense b):
    """
    Given two integer vectors a and b (of equal, nonzero length), return
    a pair of the minimum and maximum values taken on by a[i] - b[i].

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: a = vector(ZZ, [10, -30])
        sage: b = vector(ZZ, [15, -60])
        sage: min_max_delta_intvec(a, b)
        (30, -5)
    """
    assert len(a) == len(b)
    assert len(a) > 0

    cdef Integer max = Integer()
    cdef Integer min = Integer()

    cdef mpz_t tmp
    mpz_init(tmp)

    cdef int i
    for i from 0 <= i < len(a):
        mpz_sub(tmp, a._entries[i], b._entries[i])
        if i == 0 or mpz_cmp(tmp, max.value) > 0:
            mpz_set(max.value, tmp)
        if i == 0 or mpz_cmp(tmp, min.value) < 0:
            mpz_set(min.value, tmp)

    mpz_clear(tmp)

    return (max, min)


def min_max_diff_intvec(Vector_integer_dense b):
    """
    Given an integer vector b = (b0, ..., bn), compute the
    minimum and maximum values of b_{j+1} - b_j.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: min_max_diff_intvec(vector(ZZ, [1, 7, -2]))
        (-9, 6)
    """
    l = len(b)
    assert l > 1

    cdef Integer min_diff = b[1] - b[0]
    cdef Integer max_diff = Integer()

    cdef Integer diff = Integer()

    mpz_set(max_diff.value, min_diff.value)

    for i from 1 <= i < l-1:
        mpz_sub(diff.value, b._entries[i+1], b._entries[i])
        if mpz_cmp(diff.value, max_diff.value) > 0:
            mpz_set(max_diff.value, diff.value)
        if mpz_cmp(diff.value, min_diff.value) < 0:
            mpz_set(min_diff.value, diff.value)

    return (min_diff, max_diff)


def min_max_diff_doublevec(Vector_real_double_dense c):
    """
    Given a floating-point vector b = (b0, ..., bn), compute the
    minimum and maximum values of b_{j+1} - b_j.

    EXAMPLES::

        sage: from sage.rings.polynomial.real_roots import *
        sage: min_max_diff_doublevec(vector(RDF, [1, 7, -2]))
        (-9.0, 6.0)
    """
    cdef numpy.ndarray[double, ndim=1] cd = c._vector_numpy

    l = len(c)
    assert l > 1

    cdef double min_diff = cd[1] - cd[0]
    cdef double max_diff = min_diff

    cdef double diff

    for i from 1 <= i < l-1:
        diff = cd[i+1] - cd[i]
        if diff < min_diff:
            min_diff = diff
        if diff > max_diff:
            max_diff = diff

    return (min_diff, max_diff)
