r"""
`p`-adic Base Leaves

Implementations of `\ZZ_p` and `\QQ_p`

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
- William Stein: doctest updates

EXAMPLES:

`p`-adic rings and fields are examples of inexact structures, as the
reals are.  That means that elements cannot generally be stored
exactly: to do so would take an infinite amount of storage.  Instead,
we store an approximation to the elements with varying precision.

There are two types of precision for a `p`-adic element.  The first is
relative precision, which gives the number of known `p`-adic digits::

    sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
    2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
    20

The second type of precision is absolute precision, which gives the
power of `p` that this element is stored modulo::

    sage: a.precision_absolute()
    22

The number of times that `p` divides the element is called the
valuation, and can be accessed with the methods :meth:`valuation` and
:meth:`ordp`:

    sage: a.valuation()
    2

The following relationship holds:

``self.valuation() + self.precision_relative() == self.precision_absolute().``

    sage: a.valuation() + a.precision_relative() == a.precision_absolute()
    True

In the capped relative case, the relative precision of an element
is restricted to be at most a certain value, specified at the
creation of the field.  Individual elements also store their own
precision, so the effect of various arithmetic operations on
precision is tracked.  When you cast an exact element into a
capped relative field, it truncates it to the precision cap of the
field.::

    sage: R = Qp(5, 5); a = R(4006); a
    1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
    4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
    5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
    4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
    4*5 + 4*5^2 + 5^4 + O(5^5)

::

    sage: R = Zp(5, 5, 'capped-rel', 'series'); a = R(4006); a
    1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
    4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
    5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
    4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
    4*5 + 4*5^2 + 5^4 + O(5^5)

In the capped absolute type, instead of having a cap on the
relative precision of an element there is instead a cap on the
absolute precision.  Elements still store their own precisions,
and as with the capped relative case, exact elements are truncated
when cast into the ring.::

    sage: R = ZpCA(5, 5); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(4025); b
    5^2 + 2*5^3 + 5^4 + O(5^5)
    sage: a * b
    5^3 + 2*5^4 + O(5^5)
    sage: (a * b) // 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) // 5^3)
    <class 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
    sage: (a * b) / 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) / 5^3)
    <class 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>

The fixed modulus type is the leanest of the `p`-adic rings: it is
basically just a wrapper around `\ZZ / p^n \ZZ`
providing a unified interface with the rest of the `p`-adics.  This is
the type you should use if your primary interest is in speed (though
it's not all that much faster than other `p`-adic types).  It does not
track precision of elements.::

    sage: R = ZpFM(5, 5); a = R(4005); a
    5 + 2*5^3 + 5^4
    sage: a // 5
    1 + 2*5^2 + 5^3

`p`-adic rings and fields should be created using the creation
functions :func:`Zp` and :func:`Qp` as above.  This will ensure that there is
only one instance of `\ZZ_p` and `\QQ_p` of a given
type, `p`, print mode and precision.  It also saves typing very long
class names.::

    sage: Qp(17,10)
    17-adic Field with capped relative precision 10
    sage: R = Qp(7, prec = 20, print_mode = 'val-unit'); S = Qp(7, prec = 20, print_mode = 'val-unit'); R is S
    True
    sage: Qp(2)
    2-adic Field with capped relative precision 20

Once one has a `p`-adic ring or field, one can cast elements into it
in the standard way.  Integers, ints, longs, Rationals, other `p`-adic
types, pari `p`-adics and elements of `\ZZ / p^n \ZZ`
can all be cast into a `p`-adic field.::

    sage: R = Qp(5, 5, 'capped-rel','series'); a = R(16); a
    1 + 3*5 + O(5^5)
    sage: b = R(23/15); b
    5^-1 + 3 + 3*5 + 5^2 + 3*5^3 + O(5^4)
    sage: S = Zp(5, 5, 'fixed-mod','val-unit'); c = S(Mod(75,125)); c
    5^2 * 3
    sage: R(c)
    3*5^2 + O(5^5)

In the previous example, since fixed-mod elements don't keep track
of their precision, we assume that it has the full precision of
the ring.  This is why you have to cast manually here.

While you can cast explicitly as above, the chains of automatic
coercion are more restricted.  As always in Sage, the following
arrows are transitive and the diagram is commutative.::

    int -> long -> Integer -> Zp capped-rel -> Zp capped_abs -> IntegerMod
    Integer -> Zp fixed-mod -> IntegerMod
    Integer -> Zp capped-abs -> Qp capped-rel

In addition, there are arrows within each type.  For capped relative
and capped absolute rings and fields, these arrows go from lower
precision cap to higher precision cap.  This works since elements
track their own precision: choosing the parent with higher precision
cap means that precision is less likely to be truncated unnecessarily.
For fixed modulus parents, the arrow goes from higher precision cap to
lower.  The fact that elements do not track precision necessitates
this choice in order to not produce incorrect results.

TESTS::

    sage: R = Qp(5, 15, print_mode='bars', print_sep='&')
    sage: repr(R(2777))[3:]
    '0&0&0&0&0&0&0&0&0&0&4&2&1&0&2'
    sage: TestSuite(R).run()                                                            # needs sage.geometry.polyhedron

    sage: R = Zp(5, 15, print_mode='bars', print_sep='&')
    sage: repr(R(2777))[3:]
    '0&0&0&0&0&0&0&0&0&0&4&2&1&0&2'
    sage: TestSuite(R).run()                                                            # needs sage.geometry.polyhedron

    sage: R = ZpCA(5, 15, print_mode='bars', print_sep='&')
    sage: repr(R(2777))[3:]
    '0&0&0&0&0&0&0&0&0&0&4&2&1&0&2'
    sage: TestSuite(R).run()                                                            # needs sage.geometry.polyhedron
"""

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.richcmp import op_LE

from .generic_nodes import pAdicFieldBaseGeneric, \
                          pAdicCappedRelativeFieldGeneric, \
                          pAdicRingBaseGeneric, \
                          pAdicCappedRelativeRingGeneric, \
                          pAdicFixedModRingGeneric, \
                          pAdicCappedAbsoluteRingGeneric, \
                          pAdicFloatingPointRingGeneric, \
                          pAdicFloatingPointFieldGeneric, \
                          pAdicLatticeGeneric, \
                          pAdicRelaxedGeneric
from .padic_capped_relative_element import pAdicCappedRelativeElement
from .padic_capped_absolute_element import pAdicCappedAbsoluteElement
from .padic_fixed_mod_element import pAdicFixedModElement
from .padic_floating_point_element import pAdicFloatingPointElement

from sage.rings.integer_ring import ZZ


class pAdicRingCappedRelative(pAdicRingBaseGeneric, pAdicCappedRelativeRingGeneric):
    r"""
    An implementation of the `p`-adic integers with capped relative
    precision.
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime

        EXAMPLES::

            sage: R = ZpCR(next_prime(10^60))  # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingCappedRelative_with_category'>

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = ZpCR(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^10)],              # long time
            ....:                  max_runs=2^12, skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpCR(3, 1)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^3)])
            sage: R = ZpCR(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^6)],               # long time
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpCR(next_prime(10^60))
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^3)],               # long time
            ....:                  max_runs=2^5, skip='_test_log')
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])          # long time
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = Zp(17)
            sage: K(1) + 1  # indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Zp(17,10))
            True
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        #if isinstance(R, pAdicRingRelaxed) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingCappedRelative) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: Zp(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring with capped relative precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

    def _magma_init_(self, magma):
        """
        Conversion to magma.

        EXAMPLES::

            sage: # optional - magma
            sage: F = Qp(5,7,"capped-rel")
            sage: magma(F)
            ?
        """
        return f"pAdicRing({self.prime()},{self.precision_cap()})"


class pAdicRingCappedAbsolute(pAdicRingBaseGeneric, pAdicCappedAbsoluteRingGeneric):
    r"""
    An implementation of the `p`-adic integers with capped absolute precision.
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime

        EXAMPLES::

            sage: R = ZpCA(next_prime(10^60))  # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingCappedAbsolute_with_category'>

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = ZpCA(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^10)],              # long time
            ....:                  max_runs=2^12, skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpCA(3, 1)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^3)])
            sage: R = ZpCA(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^6)],               # long time
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpCA(next_prime(10^60))
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^3)],               # long time
            ....:                  max_runs=2^5, skip='_test_log')
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedAbsoluteElement, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = ZpCA(17)
            sage: K(1) + 1  # indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
            sage: K.has_coerce_map_from(ZpCA(17,10))
            True
            sage: K.has_coerce_map_from(Zp(17,40))
            True
        """
        #if isinstance(R, pAdicRingRelaxed) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingCappedRelative) and R.prime() == self.prime():
            return True
        if isinstance(R, pAdicRingCappedAbsolute) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: ZpCA(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring with capped absolute precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)


class pAdicRingFloatingPoint(pAdicRingBaseGeneric, pAdicFloatingPointRingGeneric):
    r"""
    An implementation of the `p`-adic integers with floating point
    precision.
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime

        EXAMPLES::

            sage: R = ZpFP(next_prime(10^60))  # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingFloatingPoint_with_category'>

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = ZpFP(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^10)],              # long time
            ....:                  max_runs=2^12, skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpFP(3, 1)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^3)])
            sage: R = ZpFP(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^6)],               # long time
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpFP(next_prime(10^60))
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^3)],               # long time
            ....:                  max_runs=2^5, skip='_test_log')
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFloatingPointElement, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = ZpFP(17)
            sage: K(1) + 1  # indirect doctest
            2
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Zp(17,10))
            False
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        if isinstance(R, pAdicRingFloatingPoint) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            elif R.precision_cap() == self.precision_cap() and self._printer.richcmp_modes(R._printer, op_LE):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from ``R`` to this ring.

        EXAMPLES::

            sage: ZpFP(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring with floating precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)


class pAdicRingFixedMod(pAdicRingBaseGeneric, pAdicFixedModRingGeneric):
    r"""
    An implementation of the `p`-adic integers using fixed modulus.
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime

        EXAMPLES::

            sage: R = ZpFM(next_prime(10^60))  # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingFixedMod_with_category'>

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = ZpFM(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^10)],  # long time
            ....:                  max_runs=2^12, skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = ZpFM(3, 1)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^3)])
            sage: R = ZpFM(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^6)],  # long time
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)]) # long time
            sage: R = ZpFM(next_prime(10^60))
            sage: TestSuite(R).run(skip='_test_log')
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^4)],  # long time
            ....:                  max_runs=2^6, skip='_test_log')
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])

        Fraction fields work after :issue:`23510`::

            sage: R = ZpFM(5)
            sage: K = R.fraction_field(); K
            5-adic Field with floating precision 20
            sage: K(R(90))
            3*5 + 3*5^2
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFixedModElement, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = ZpFM(17)
            sage: K(1) + 1  # indirect doctest
            2
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Zp(7))
            False
            sage: K.has_coerce_map_from(ZpFM(17,40))
            True
            sage: K.has_coerce_map_from(ZpFM(17,10))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
        """
        #if isinstance(R, pAdicRingRelaxed) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingFixedMod) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: ZpFM(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring of fixed modulus 7^20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

    def _magma_init_(self, magma):
        """
        Conversion to magma.

        EXAMPLES::

            sage: # optional - magma
            sage: F = Zp(5,7,"fixed-mod")
            sage: magma(F)
            ?
        """
        return f"pAdicQuotientRing({self.prime()},{self.precision_cap()})"


class pAdicFieldCappedRelative(pAdicFieldBaseGeneric, pAdicCappedRelativeFieldGeneric):
    r"""
    An implementation of `p`-adic fields with capped relative precision.

    EXAMPLES::

        sage: K = Qp(17, 1000000)  # indirect doctest
        sage: K = Qp(101)  # indirect doctest
    """

    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime

        EXAMPLES::

            sage: K = Qp(next_prime(10^60))  # indirect doctest
            sage: type(K)
            <class 'sage.rings.padics.padic_base_leaves.pAdicFieldCappedRelative_with_category'>

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = Qp(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^10)],  # long time
            ....:                  max_runs=2^12, skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time

            sage: R = Qp(3, 1)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^6)],  # long time, needs sage.geometry.polyhedron
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])  # long time

            sage: R = Qp(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^9)],  # long time, needs sage.geometry.polyhedron
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(3^3)])

            sage: R = Qp(next_prime(10^60))
            sage: TestSuite(R).run(skip='_test_log')                                    # needs sage.geometry.polyhedron
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^3)],  # long time, needs sage.geometry.polyhedron
            ....:                  max_runs=2^5, skip='_test_log')
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = Qp(17)
            sage: K(1) + 1  # indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            True
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Qp(17,40))
            False
            sage: K.has_coerce_map_from(Qp(17,10))
            True
            sage: K.has_coerce_map_from(Zp(17,40))
            True
        """
        #if isinstance(R, pAdicRingRelaxed) or isinstance(R, pAdicFieldRelaxed) and R.prime() == self.prime():
        #    return True
        if isinstance(R, (pAdicRingCappedRelative, pAdicRingCappedAbsolute)) and R.prime() == self.prime():
            return True
        if isinstance(R, pAdicFieldCappedRelative) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: Qp(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Field with capped relative precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

    def _magma_init_(self, magma):
        """
        Conversion to magma.

        EXAMPLES::

            sage: # optional - magma
            sage: F = Qp(5,7,"capped-rel")
            sage: magma(F)
            ?
        """
        return f"pAdicField({self.prime()},{self.precision_cap()})"

    def random_element(self, algorithm='default'):
        r"""
        Return a random element of ``self``, optionally using the ``algorithm``
        argument to decide how it generates the element. Algorithms currently
        implemented:

        - ``'default'``: Choose an integer `k` using the standard
          distribution on the integers.  Then choose an integer `a`
          uniformly in the range `0 \le a < p^N` where `N` is the
          precision cap of ``self``.  Return ``self(p^k * a, absprec =
          k + self.precision_cap())``.

        EXAMPLES::

            sage: Qp(17,6).random_element().parent() is Qp(17,6)
            True
        """
        if (algorithm == 'default'):
            k = ZZ.random_element()
            a = ZZ.random_element(self.prime()**self.precision_cap())
            return self(self.prime()**k * a, absprec=k + self.precision_cap())
        else:
            raise NotImplementedError("Don't know %s algorithm" % algorithm)


class pAdicFieldFloatingPoint(pAdicFieldBaseGeneric, pAdicFloatingPointFieldGeneric):
    r"""
    An implementation of the `p`-adic rationals with floating point
    precision.
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime

        EXAMPLES::

            sage: R = QpFP(next_prime(10^60))  # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicFieldFloatingPoint_with_category'>

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = QpFP(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^10)],  # long time
            ....:                  max_runs=2^12, skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = QpFP(3, 1)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^3)])
            sage: R = QpFP(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^6)],  # long time
            ....:                  skip='_test_metric_function')
            sage: R._test_metric_function(elements=[R.random_element() for i in range(2^3)])        # long time
            sage: R = QpFP(next_prime(10^60))
            sage: TestSuite(R).run(skip='_test_log')
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(2^3)],  # long time
            ....:                  max_runs=2^5, skip='_test_log')
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFloatingPointElement, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = QpFP(17)
            sage: K(1) + 1  # indirect doctest
            2
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            True
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Qp(17,10))
            False
            sage: K.has_coerce_map_from(ZpFP(17))
            True
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        if isinstance(R, (pAdicRingFixedMod, pAdicRingFloatingPoint, pAdicFieldFloatingPoint)) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            elif R.precision_cap() == self.precision_cap() and self._printer.richcmp_modes(R._printer, op_LE):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: QpFP(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Field with floating precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

# Lattice precision
###################


class pAdicRingLattice(pAdicLatticeGeneric, pAdicRingBaseGeneric):
    """
    An implementation of the `p`-adic integers with lattice precision.

    INPUT:

    - ``p`` -- prime

    - ``prec`` -- precision cap, given as a pair (``relative_cap``, ``absolute_cap``)

    - ``subtype`` -- either ``'cap'`` or ``'float'``

    - ``print_mode`` -- dictionary with print options

    - ``names`` -- how to print the prime

    - ``label`` -- the label of this ring

    .. SEEALSO::

        :meth:`label`

    EXAMPLES::

        sage: R = ZpLC(next_prime(10^60))  # indirect doctest
        doctest:...: FutureWarning: This class/method/function is marked as experimental.
        It, its functionality or its interface might change without a formal deprecation.
        See https://github.com/sagemath/sage/issues/23505 for details.
        sage: type(R)
        <class 'sage.rings.padics.padic_base_leaves.pAdicRingLattice_with_category'>

        sage: R = ZpLC(2, label='init')  # indirect doctest
        sage: R
        2-adic Ring with lattice-cap precision (label: init)
    """
    def __init__(self, p, prec, subtype, print_mode, names, label=None, category=None):
        """
        Initialization.

        TESTS:

            sage: R = ZpLC(7, label='init')
            sage: TestSuite(R).run(skip=['_test_teichmuller', '_test_matrix_smith'])  # long time
        """
        # We need to set the subtype first, so that
        # pAdicRingBaseGeneric.__init__ can work
        self._subtype = subtype
        if isinstance(prec, tuple):
            pAdicRingBaseGeneric.__init__(self, p, prec[1], print_mode, names, None, category=category)
        else:
            pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, None, category=category)
        pAdicLatticeGeneric.__init__(self, p, prec, print_mode, names, label, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to this ring.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: R.has_coerce_map_from(ZZ)
            True
            sage: R.has_coerce_map_from(QQ)
            False

            sage: K = R.fraction_field()
            sage: K.has_coerce_map_from(R)
            True
            sage: K.has_coerce_map_from(QQ)
            True

        Note that coerce map does not exist between ``p``-adic rings with
        lattice precision and other ``p``-adic rings. ::

            sage: S = Zp(2)
            sage: R.has_coerce_map_from(S)
            False
            sage: S.has_coerce_map_from(R)
            False

        Similarly there is no coercion maps between ``p``-adic rings with
        different labels. ::

            sage: R2 = ZpLC(2, label='coerce')
            sage: R.has_coerce_map_from(R2)
            False
            sage: R2.has_coerce_map_from(R)
            False
        """
        if isinstance(R, pAdicRingLattice) and R.precision() is self.precision():
            return True

    def random_element(self, prec=None):
        """
        Return a random element of this ring.

        INPUT:

        - ``prec`` -- integer or ``None`` (default); the
          absolute precision of the generated random element

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: R.random_element()    # random
            2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^10 + 2^11 + 2^14 + 2^15 + 2^16
             + 2^17 + 2^18 + 2^19 + 2^21 + O(2^23)

            sage: R.random_element(prec=10)    # random
            1 + 2^3 + 2^4 + 2^7 + O(2^10)
        """
        p = self.prime()
        if self._subtype == 'cap':
            if prec is None:
                prec = self._prec_cap_absolute
            x = ZZ.random_element(p**prec)
            relcap = x.valuation(p) + self._prec_cap_relative
            prec = min(relcap, prec)
            return self._element_class(self, x, prec=prec)
        else:
            if prec is None:
                cap = self._prec_cap_relative
            else:
                cap = prec
            x = ZZ.random_element(p**cap)
            v = x.valuation(p)
            if prec is None and v > 0:
                x += p**cap * ZZ.random_element(p**v)
            return self._element_class(self, x, prec=prec)


class pAdicFieldLattice(pAdicLatticeGeneric, pAdicFieldBaseGeneric):
    """
    An implementation of the `p`-adic numbers with lattice precision.

    INPUT:

    - ``p`` -- prime

    - ``prec`` -- precision cap, given as a pair (``relative_cap``, ``absolute_cap``)

    - ``subtype`` -- either ``'cap'`` or ``'float'``

    - ``print_mode`` -- dictionary with print options

    - ``names`` -- how to print the prime

    - ``label`` -- the label of this ring

    .. SEEALSO::

        :meth:`label`

    EXAMPLES::

        sage: R = QpLC(next_prime(10^60))  # indirect doctest
        doctest:...: FutureWarning: This class/method/function is marked as experimental.
        It, its functionality or its interface might change without a formal deprecation.
        See https://github.com/sagemath/sage/issues/23505 for details.
        sage: type(R)
        <class 'sage.rings.padics.padic_base_leaves.pAdicFieldLattice_with_category'>

        sage: R = QpLC(2,label='init')  # indirect doctest
        sage: R
        2-adic Field with lattice-cap precision (label: init)
    """
    def __init__(self, p, prec, subtype, print_mode, names, label=None, category=None):
        """
        Initialization.

        TESTS::

            sage: R = QpLC(7, label='init')
            sage: TestSuite(R).run(skip=['_test_teichmuller', '_test_matrix_smith'])  # long time
        """
        # We need to set the subtype first, so that
        # pAdicFieldBaseGeneric.__init__ can work
        self._subtype = subtype
        if isinstance(prec, tuple):
            pAdicFieldBaseGeneric.__init__(self, p, prec[1], print_mode, names, None, category=category)
        else:
            pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, None, category=category)
        pAdicLatticeGeneric.__init__(self, p, prec, print_mode, names, label, category=category)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to this ring.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: R.has_coerce_map_from(ZZ)
            True
            sage: R.has_coerce_map_from(QQ)
            False

            sage: K = R.fraction_field()
            sage: K.has_coerce_map_from(R)
            True
            sage: K.has_coerce_map_from(QQ)
            True

        Note that coerce map does not exist between ``p``-adic fields with
        lattice precision and other ``p``-adic rings. ::

            sage: L = Qp(2)
            sage: K.has_coerce_map_from(L)
            False
            sage: L.has_coerce_map_from(K)
            False

        Similarly there is no coercion maps between ``p``-adic rings with
        different labels.

            sage: K2 = QpLC(2, label='coerce')
            sage: K.has_coerce_map_from(K2)
            False
            sage: K2.has_coerce_map_from(K)
            False
        """
        if isinstance(R, (pAdicRingLattice, pAdicFieldLattice)) and R.precision() is self.precision():
            return True

    def random_element(self, prec=None, integral=False):
        """
        Return a random element of this ring.

        INPUT:

        - ``prec`` -- integer or ``None`` (default); the
          absolute precision of the generated random element

        - ``integral`` -- boolean (default: ``False``); if ``True``,
          return an element in the ring of integers

        EXAMPLES::

            sage: K = QpLC(2)
            sage: K.random_element()   # not tested, known bug (see :issue:`32126`)
            2^-8 + 2^-7 + 2^-6 + 2^-5 + 2^-3 + 1 + 2^2 + 2^3 + 2^5 + O(2^12)
            sage: K.random_element(integral=True)    # random
            2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^10 + 2^11 + 2^14 + 2^15 + 2^16
             + 2^17 + 2^18 + 2^19 + O(2^20)

            sage: K.random_element(prec=10)    # random
            2^(-3) + 1 + 2 + 2^4 + 2^8 + O(2^10)

        If the given precision is higher than the internal cap of the
        parent, then the cap is used::

            sage: K.precision_cap_relative()
            20
            sage: K.random_element(prec=100)    # random
            2^5 + 2^8 + 2^11 + 2^12 + 2^14 + 2^18 + 2^20 + 2^24 + O(2^25)
        """
        if integral:
            val = 0
        else:
            val = ZZ.random_element()
        if prec is None:
            prec = self._prec_cap_absolute - val
        p = self.prime()
        x = ZZ.random_element(p**prec)
        relcap = x.valuation(p) + self._prec_cap_relative
        prec = min(relcap, prec)
        return self._element_class(self, x*(p**val), prec=prec)

# Relaxed
#########


class pAdicRingRelaxed(pAdicRelaxedGeneric, pAdicRingBaseGeneric):
    r"""
    An implementation of relaxed arithmetics over `\ZZ_p`.

    INPUT:

    - ``p`` -- prime

    - ``prec`` -- default precision

    - ``print_mode`` -- dictionary with print options

    - ``names`` -- how to print the prime

    EXAMPLES::

        sage: R = ZpER(5)  # indirect doctest                                           # needs sage.libs.flint
        sage: type(R)                                                                   # needs sage.libs.flint
        <class 'sage.rings.padics.padic_base_leaves.pAdicRingRelaxed_with_category'>
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R = ZpER(7)
            sage: TestSuite(R).run(skip=['_test_log', '_test_matrix_smith'])
            sage: R = ZpER(7, secure=True)
            sage: TestSuite(R).run(skip=['_test_log', '_test_matrix_smith'])
        """
        from sage.rings.padics import padic_relaxed_element
        self._default_prec, self._halting_prec, self._secure = prec
        pAdicRingBaseGeneric.__init__(self, p, self._default_prec, print_mode, names, padic_relaxed_element.pAdicRelaxedElement, category=category)
        self._element_class_module = padic_relaxed_element
        self._element_class_prefix = "pAdicRelaxedElement_"


class pAdicFieldRelaxed(pAdicRelaxedGeneric, pAdicFieldBaseGeneric):
    r"""
    An implementation of relaxed arithmetics over `\QQ_p`.

    INPUT:

    - ``p`` -- prime

    - ``prec`` -- default precision

    - ``print_mode`` -- dictionary with print options

    - ``names`` -- how to print the prime

    EXAMPLES::

        sage: R = QpER(5)  # indirect doctest                                           # needs sage.libs.flint
        sage: type(R)                                                                   # needs sage.libs.flint
        <class 'sage.rings.padics.padic_base_leaves.pAdicFieldRelaxed_with_category'>
    """
    def __init__(self, p, prec, print_mode, names, category=None):
        """
        Initialization.

        TESTS::

            sage: # needs sage.libs.flint
            sage: K = QpER(7)
            sage: TestSuite(K).run(skip=['_test_log', '_test_matrix_smith'])
            sage: K = QpER(7, secure=True)
            sage: TestSuite(K).run(skip=['_test_log', '_test_matrix_smith'])
        """
        from sage.rings.padics import padic_relaxed_element
        self._default_prec, self._halting_prec, self._secure = prec
        pAdicFieldBaseGeneric.__init__(self, p, self._default_prec, print_mode, names, padic_relaxed_element.pAdicRelaxedElement, category=category)
        self._element_class_module = padic_relaxed_element
        self._element_class_prefix = "pAdicRelaxedElement_"
