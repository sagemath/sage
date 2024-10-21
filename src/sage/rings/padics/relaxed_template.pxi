# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.libs.flint
r"""
Template for relaxed `p`-adic rings and fields.

In order to use this template you need to write a linkage file and
gluing file.

The linkage file implements a common API that is then used in the
class FMElement defined here.
See sage/libs/linkages/padics/relaxed/API.pxi for the functions needed.

The gluing file does the following:
- ctypedef's ``cdigit``, ``cdigit_ptr``, ``celement`` and ``randgen``
  to be the appropriate types
- includes the linkage file
- includes this template
- defines classes inheriting for classes defined in this file,
  and implements any desired extra methods
See padic_relaxed_element.pxd and padic_relaxed_element.pyx for an example.

AUTHOR:

- Xavier Caruso (2021-02): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import re
from sage.misc import persist

from libc.stdlib cimport malloc, free
from sage.libs.gmp.mpz cimport mpz_sizeinbase, mpz_tstbit

from sage.structure.element import coerce_binop
from sage.structure.element cimport have_same_parent
from sage.structure.coerce cimport coercion_model
from sage.misc.prandom import randint

from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.padics.padic_relaxed_errors cimport *
from sage.rings.padics.padic_relaxed_errors import raise_error

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1


cdef cdigit tmp_digit
cdef celement tmp_poly
digit_init(tmp_digit)
element_init(tmp_poly)


cdef class RelaxedElement(pAdicGenericElement):
    r"""
    Template class for relaxed `p`-adic elements.

    EXAMPLES:

        sage: from sage.rings.padics.padic_relaxed_element import RelaxedElement
        sage: R = ZpER(5)
        sage: a = R(1)
        sage: isinstance(a, RelaxedElement)
        True
    """
    def __init__(self, parent):
        r"""
        Initialize this element.

        TESTS::

            sage: R = ZpER(5)
            sage: a = R(1/2)  # indirect doctest
            sage: TestSuite(a).run()
        """
        pAdicGenericElement.__init__(self, parent)
        self.prime_pow = <PowComputer_class>self._parent.prime_pow
        self._valuation = 0
        self._precrel = 0
        self._precbound = maxordp

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle
        this element.

        TESTS::

            sage: R = ZpER(5)
            sage: a = R(0)
            sage: loads(dumps(a)) == a
            True
        """
        raise NotImplementedError("must be implemented in subclasses")

    cpdef bint _is_base_elt(self, p) except -1:
        r"""
        Return ``True`` if this element is an element of Zp or Qp (rather than
        an extension).

        INPUT:

        - ``p`` -- a prime, which is compared with the parent of this element

        EXAMPLES::

            sage: a = ZpER(5)(3)
            sage: a._is_base_elt(5)
            True
            sage: a._is_base_elt(17)
            False
        """
        return p == self._parent.prime()

    cdef cdigit_ptr _getdigit_relative(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th significant digit of this number.

        .. NOTE::

            This function does not check that the requested digit
            has been already computed.

        INPUT:

        - ``i`` -- positive integer
        """
        pass

    cdef cdigit_ptr _getdigit_absolute(self, long i) noexcept:
        r"""
        Return a pointer on the digit in position `i` of
        this number.

        .. NOTE::

            This function do not check that the requested digit
            has been already computed.

        INPUT:

        - ``i`` -- integer
        """
        pass

    cdef void _getslice_relative(self, celement slice, long start, long length) noexcept:
        r"""
        Select a slice of the digits of this number.

        INPUT:

        - ``slice`` -- a ``celement`` to store the slice

        - ``start`` -- positive integer; the starting position of the slice
          in relative precision

        - ``length`` -- positive integer; the length of the slice

        .. NOTE::

            This methods only sets up a pointer to the requested slice
            (the slice is not copied). Hence any future modification
            of the slice ``slice`` will affect this number.
        """
        pass

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT: an error code which is a superposition of the following:

        - ``0`` -- no error
        - ``ERROR_ABANDON    = 1`` -- computation has been abandoned
        - ``ERROR_NOTDEFINED = 2`` -- a number is not defined
        - ``ERROR_PRECISION  = 4`` -- not enough precision
        - ``ERROR_OVERFLOW   = 8`` -- overflow
        - ``ERROR_NOTSQUARE  = 16`` -- not a square
        - ``ERROR_INTEGRAL   = 32`` -- not in the integer ring
        - ``ERROR_DIVISION   = 64`` -- division by zero (or something indistinguishable from zero)
        - ``ERROR_CIRCULAR   = 128`` -- circular definition
        """
        raise NotImplementedError("must be implemented in subclasses")

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`_next_c` for details).
        """
        cdef int error
        cdef long pr = min(prec, self._precbound)
        while self._precrel + self._valuation < pr:
            error = self._next_c()
            if error:
                return error
        if prec > self._precbound:
            return ERROR_PRECISION
        return 0

    cdef int _jump_relative_c(self, long prec, long halt) noexcept:
        r"""
        Compute the digits of this number until the relative precision
        ``prec``.

        INPUT:

        - ``prec`` -- integer

        - ``halt`` -- integer; the absolute precision after which the
          computation is abandoned if the first significant digit has not
          been found yet

        OUTPUT:

        An error code (see :meth:`_next_c` for details).
        """
        if self._valuation >= maxordp:
            return 0
        cdef int error = 0
        if self._valuation <= -maxordp:
            error = self._next_c()
        cdef long valhalt = min(halt, self._precbound)
        while not error and self._valuation < valhalt and self._precrel == 0:
            error = self._next_c()
        if self._valuation >= self._precbound:
            error |= ERROR_PRECISION
        elif self._valuation >= halt:
            error |= ERROR_ABANDON
        if not error:
            error = self._jump_c(self._valuation + prec)
        return error

    cdef int _init_jump(self) except -1:
        r"""
        Start the computation of the digits of this number.
        This method should be called at initialization.

        OUTPUT:

        An error code (see :meth:`_next_c` for details).
        """
        cdef int error = 0
        if self._precbound < maxordp:
            error = self._jump_c(self._precbound)
            raise_error(error, True)
        return error

    def digit(self, i):
        r"""
        Return the coefficient of `p^i` in the `p`-adic expansion
        of this number.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: R = ZpER(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.digit(0)
            0
            sage: a.digit(1)
            4
            sage: a.digit(7)
            3
            sage: a.digit(100)
            1

        As a shortcut, one can use the bracket operator::

            sage: a[100]
            1

        If the digit is not known, an error is raised::

            sage: b = a.add_bigoh(10)
            sage: b.digit(20)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision

        TESTS::

            sage: a.digit(-1)
            0
        """
        cdef int error = self._jump_c(i+1)
        raise_error(error)
        cdef cdigit_ptr coeff = self._getdigit_absolute(i)
        return digit_get_sage(coeff)

    def expansion(self, n=None, lift_mode='simple', start_val=None):
        r"""
        Return an iterator over the list of coefficients in a `p`-adic
        expansion of this element, that is the list of `a_i` so that
        this element can be expressed as

        .. MATH::

            \pi^v \cdot \sum_{i=0}^\infty a_i \pi^i

        where `v` is the valuation of this element when the parent is
        a field, and `v = 0` otherwise and the `a_i` are between `0`
        and `p - 1`.

        INPUT:

        - ``n`` -- integer or ``None`` (default: ``None``); if
          given, return the corresponding entries in the expansion

        - ``lift_mode`` -- ``'simple'``, ``'smallest'`` or
          ``'teichmuller'`` (default: ``'simple'``)

        - ``start_val`` -- start at this valuation rather than the
          default (`0` or the valuation of this element)

        OUTPUT:

        - If ``n`` is ``None``, an iterable giving the `p`-adic expansion
          of this element.

        - If ``n`` is an integer, the coefficient of `p^n` in the
          `p`-adic expansion of this element.

        EXAMPLES::

            sage: R = ZpER(7, print_mode='digits')
            sage: a = R(1/2021); a
            ...23615224635636163463

        Without any argument, this method returns an iterator over the
        digits of this number::

            sage: E = a.expansion()
            sage: E
            7-adic expansion of ...23615224635636163463
            sage: next(E)
            3
            sage: next(E)
            6
            sage: next(E)
            4

        On the contrary, passing in an integer returns the digit at the
        given position::

            sage: a.expansion(5)
            1

        Over a field, the expansion starts at the valuation of the element::

            sage: K = R.fraction_field()
            sage: b = K(20/21); b
            ...2222222222222222223.2
            sage: E = b.expansion()
            sage: next(E)
            2
            sage: next(E)
            3

            sage: c = 1/b; c
            ...564356435643564356440
            sage: E = c.expansion()
            sage: next(E)
            4
            sage: next(E)
            4

        When ``start_val`` is given, the expansion starts at this position
        instead:

            sage: E = c.expansion(start_val=0)
            sage: next(E)
            0
            sage: next(E)
            4

            sage: E = c.expansion(start_val=5)
            sage: next(E)
            3
            sage: next(E)
            4
        """
        if lift_mode == 'simple':
            mode = simple_mode
        elif lift_mode == 'smallest':
            mode = smallest_mode
        elif lift_mode == 'teichmuller':
            mode = teichmuller_mode
        else:
            raise ValueError("unknown lift mode")
        if n is None:
            if start_val is not None:
                start = start_val
            elif self.prime_pow.in_field:
                start = self.valuation()
            else:
                start = 0
            return ExpansionIter(self, mode, start, self._precbound)
        else:
            n = Integer(n)
            if n >= self._precbound:
                raise_error(ERROR_PRECISION)
            if n >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            if mode == simple_mode:
                return self.digit(n)
            else:
                E = ExpansionIter(self, mode, n, n+1)
                return next(E)

    def __getitem__(self, n):
        r"""
        Return the `n`-th digit of this relaxed `p`-adic number if `n` is an integer
        or return a bounded relaxed `p`-adic corresponding to the given slice if `n` is a slice.

        INPUT:

        - ``n`` -- integer or a slice

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...

            sage: a[3]
            3
            sage: a[3:6]
            3*7^3 + 6*7^4 + 7^5 + O(7^6)

        Unbounded slices are allowed::

            sage: a[:3]
            3 + 6*7 + 4*7^2 + O(7^3)
            sage: a[3:]
            3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...

        .. SEEALSO::

            :meth:`digit`, :meth:`slice`

        TESTS::

            sage: a[3:6:2]
            Traceback (most recent call last):
            ...
            NotImplementedError: step is not allowed
        """
        if isinstance(n, slice):
            if n.step is not None:
                raise NotImplementedError("step is not allowed")
            return self.slice(n.start, n.stop, True)
        else:
            n = Integer(n)
            if n >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            return self.digit(n)

    def slice(self, start=None, stop=None, bound=False):
        r"""
        Return a slice of this number.

        INPUT:

        - ``start`` -- integer or ``None`` (default: ``None``);
          the first position of the slice

        - ``stop`` -- integer or ``None`` (default: ``None``);
          the first position not included in the slice

        - ``bound`` -- boolean (default: ``False``); whether the
          precision on the output should be bounded or unbounded

        EXAMPLES::

            sage: K = QpER(7, prec=10)
            sage: a = K(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...

            sage: s = a.slice(3, 6)
            sage: s
            3*7^3 + 6*7^4 + 7^5 + ...

        In the above example, the precision on `b` remains unbounded::

            sage: s.precision_absolute()
            +Infinity

        Passing in ``bound=True`` changes this behaviour::

            sage: a.slice(3, 6, bound=True)
            3*7^3 + 6*7^4 + 7^5 + O(7^6)

        When ``start`` is omitted, the slice starts at the beginning of the number:

            sage: a.slice(stop=6)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + ...

            sage: b = K(20/21); b
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + 2*7^6 + 2*7^7 + 2*7^8 + ...
            sage: b.slice(stop=6)
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + ...

        As a shortcut, one can use the bracket operator.
        However, in this case, the precision is bounded::

            sage: a[3:6]
            3*7^3 + 6*7^4 + 7^5 + O(7^6)
            sage: b[:6]
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + O(7^6)

        TESTS::

            sage: x = a[:3]
            sage: x.slice(stop=3)
            3 + 6*7 + 4*7^2 + ...
            sage: x.slice(stop=4)
            3 + 6*7 + 4*7^2 + O(7^3)

        Taking slices of slices work as expected::

            sage: a1 = a.slice(5, 10); a1
            7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a2 = a1.slice(3, 8); a2
            7^5 + 6*7^6 + 3*7^7 + ...
            sage: a2 == a.slice(5, 8)
            True
        """
        if start is None:
            start = -maxordp
        elif start >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        if stop is None:
            stop = maxordp
        elif stop >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        cdef RelaxedElement x = element_class_slice(self._parent, self, start, stop, 0)
        if bound and x._precbound > stop:
            x = element_class_bound(self._parent, x, stop)
        return x

    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES:

        For unbounded elements, the number of printed terms is given by
        the default precision of the ring::

            sage: R = ZpER(5, 10)
            sage: a = R(1/2021)
            sage: a   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + ...

            sage: S = ZpER(5, 5)
            sage: b = S(1/2021)
            sage: b   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + ...

        On the contrary, bounded elements are printed until their own
        precision, regardless the default precision of the parent::

            sage: a[:15]   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + 2*5^10 + 5^11 + 5^12 + 5^13 + O(5^15)

            sage: b[:15]   # indirect doctest
            1 + 5 + 3*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + 2*5^10 + 5^11 + 5^12 + 5^13 + O(5^15)
        """
        # This code should be integrated to the p-adic printer
        if self._valuation <= -maxordp:
            return "O(%s^-Infinity)" % self._parent.prime()
        if self._is_exact_zero():
            return "0"
        if self._precbound >= maxordp:
            unbounded = True
            try:
                if self.prime_pow.in_field:
                    x = self.at_precision_relative(permissive=False)
                else:
                    x = self.at_precision_absolute(permissive=False)
            except (ValueError, PrecisionError, RecursionError):
                unbounded = False
                x = element_class_bound(self._parent, self, self._valuation + self._precrel)
        else:
            unbounded = False
            x = self
        s = pAdicGenericElement._repr_(x)
        mode = self._parent._printer.dict()['mode']
        if unbounded:
            s = re.sub(r'O\(.*\)$', '...', s)
        elif mode == "digits":
            s = re.sub(r'^\.\.\.\??', '...?', s)
        elif mode == "bars":
            s = re.sub(r'^\.\.\.\|?', '...?|', s)
        if s == "" or s == "...":
            s = "0 + ..."
        return s

    cdef bint _is_equal(self, RelaxedElement right, long prec, bint permissive) except -1:
        r"""
        C function for checking equality at a given precision.

        INPUT:

        - ``right`` -- the second element involved in the comparison

        - ``prec`` -- integer; the precision at which the equality is checked

        - ``permissive`` -- boolean; if ``True``, be silent if the precision
          on one input is less than ``prec``. Otherwise, raise an error.
        """
        cdef int error
        cdef long i
        if self._valuation <= -maxordp:
            error = self._next_c()
            raise_error(error)
        if right._valuation <= -maxordp:
            error = right._next_c()
            raise_error(error)
        for i in range(min(self._valuation, right._valuation), prec):
            error = self._jump_c(i+1)
            raise_error(error, permissive)
            if error:   # not enough precision
                return True
            error = right._jump_c(i+1)
            raise_error(error, permissive)
            if error:   # not enough precision
                return True
            if not digit_equal(self._getdigit_absolute(i), right._getdigit_absolute(i)):
                return False
        return True

    @coerce_binop
    def is_equal_at_precision(self, right, prec):
        r"""
        Compare this element with ``right`` at precision ``prec``.

        INPUT:

        - ``right`` -- a relaxed `p`-adic number

        - ``prec`` -- integer

        EXAMPLES::

            sage: R = ZpER(7, prec=10)
            sage: a = R(1/2); a
            4 + 3*7 + 3*7^2 + 3*7^3 + 3*7^4 + 3*7^5 + 3*7^6 + 3*7^7 + 3*7^8 + 3*7^9 + ...
            sage: b = R(99/2); b
            4 + 3*7 + 4*7^2 + 3*7^3 + 3*7^4 + 3*7^5 + 3*7^6 + 3*7^7 + 3*7^8 + 3*7^9 + ...

            sage: a.is_equal_at_precision(b, 1)
            True
            sage: a.is_equal_at_precision(b, 2)
            True
            sage: a.is_equal_at_precision(b, 3)
            False
        """
        return self._is_equal(right, min(prec, maxordp), False)

    @coerce_binop
    def is_equal_to(self, RelaxedElement right, prec=None, secure=None):
        r"""
        Compare this element with ``right``.

        INPUT:

        - ``right`` -- a relaxed `p`-adic number

        - ``prec`` -- integer or ``None`` (default: ``None``); if
          given, compare the two elements at this precision. Otherwise
          use the default halting precision of the parent.

        - ``secure`` -- boolean (default: ``False`` if ``prec`` is given,
          ``True`` otherwise); when the elements cannot be distinguished
          at the given precision, raise an error if ``secure`` is ``True``,
          return ``True`` otherwise.

        EXAMPLES::

            sage: R = ZpER(7)
            sage: a = R(1/2)
            sage: b = R(1/3)
            sage: c = R(1/6)

            sage: a.is_equal_to(b)
            False

        When equality indeed holds, it is not possible to conclude by
        comparing more and more accurate approximations.
        In this case, an error is raised::

            sage: a.is_equal_to(b + c)
            Traceback (most recent call last):
            ...
            PrecisionError: unable to decide equality; try to bound precision

        You can get around this behaviour by passing ``secure=False``::

            sage: a.is_equal_to(b + c, secure=False)
            True

        Another option (which is actually recommended) is to provide an explicit
        bound on the precision::

            sage: s = b + c + 7^50
            sage: a.is_equal_to(s, prec=20)
            True
            sage: a.is_equal_to(s, prec=100)
            False
        """
        cdef long halt
        if self is right:
            return True
        if self._valuation >= maxordp and right._valuation >= maxordp:
            return True
        if prec is None:
            if secure is None:
                secure = True
            prec = min(self._precbound, right._precbound)
        else:
            if secure is None:
                secure = False
            prec = Integer(prec)
        if prec < maxordp:
            return self._is_equal(right, prec, True)
        prec = min(self._valuation + self._precrel, right._valuation + right._precrel)
        halt = min(self._parent.halting_prec(), maxordp)
        eq = self._is_equal(right, max(prec, halt), True)
        if secure and eq:
            raise PrecisionError("unable to decide equality; try to bound precision")
        return eq

    def __eq__(self, other):
        r"""
        Return ``True`` of this element is equal to ``other``.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(1/2)
            sage: y = R(1/3)
            sage: z = R(1/6)

            sage: x == y + z
            True

        We illustrate the effect of the keyword ``secure``::

            sage: R.is_secure()
            False
            sage: s = y + z + 5^50
            sage: x == s
            True

            sage: S = ZpER(5, secure=True)
            sage: S(x) == S(s)
            Traceback (most recent call last):
            ...
            PrecisionError: unable to decide equality; try to bound precision

        Note that, when ``secure=False``, once more digits have been
        computed, the answer can change::

            sage: x[:100] == s
            False
            sage: x == s
            False

        .. SEEALSO::

            :meth:`is_equal_to`
        """
        if not have_same_parent(self, other):
            try:
                a, b = coercion_model.canonical_coercion(self, other)
            except TypeError:
                return False
            return a == b
        return self.is_equal_to(other, secure=self._parent.is_secure())

    def __bool__(self):
        r"""
        Return ``True`` if this element is distinguishable from zero.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(1)
            sage: bool(x)
            True

        In the next example, only `40` digits (which is the default halting
        precision) are computed so `x` is considered as indistinguishable from
        `0`::

            sage: x = R(5^41)
            sage: bool(x)
            False
        """
        cdef int error = 0
        cdef long prec, halt
        if self._precrel:
            return True
        prec = self._precbound
        if prec >= maxordp:
            halt = min(self._parent.halting_prec(), maxordp)
            prec = max(self._valuation + self._precrel, halt)
        while not error and self._precrel == 0 and self._valuation < prec:
            error = self._next_c()
        return bool(self._precrel)

    cpdef bint _is_exact_zero(self) except -1:
        r"""
        Return ``True`` if this element is an exact zero.

        EXAMPLES::

            sage: R = ZpER(5)
            sage: a = R(0); a
            0
            sage: a._is_exact_zero()
            True

            sage: b = a.add_bigoh(20)
            sage: b
            O(5^20)
            sage: b._is_exact_zero()
            False
        """
        return self._valuation >= maxordp

    cpdef bint _is_inexact_zero(self) except -1:
        r"""
        Return ``True`` if, at the current stage of computations, this
        number cannot be distinguished from zero.

        EXAMPLES::

            sage: R = ZpER(5, print_mode='digits')
            sage: a = R(20/21)

        Computations have not started yet; hence we are not able
        to distinguish `a` from zero so far::

            sage: a._is_inexact_zero()
            True

        When we print `a`, the first digits are computed, after what
        we know that `a` is not zero::

            sage: a
            ...34010434010434010440
            sage: a._is_inexact_zero()
            False
        """
        return self._precrel == 0

    def is_exact(self):
        r"""
        Return ``True`` if this element is exact, that is if its
        precision is unbounded.

        EXAMPLES::

            sage: R = ZpER(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.is_exact()
            True

            sage: b = a.add_bigoh(10)
            sage: b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + O(5^10)
            sage: b.is_exact()
            False
        """
        return self._precbound >= maxordp

    def precision_absolute(self):
        r"""
        Return the absolute precision of this element.

        This is the power of `p` modulo which this element is known.
        For unbounded elements, this methods return `+\infty`.

        EXAMPLES::

            sage: R = ZpER(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.precision_absolute()
            +Infinity

            sage: b = a.add_bigoh(10)
            sage: b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + O(5^10)
            sage: b.precision_absolute()
            10

        TESTS::

            sage: s = R.unknown()
            sage: (1/s).precision_absolute()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known
        """
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if self._precbound >= maxordp:
            return Infinity
        return Integer(self._precrel + self._valuation)

    def precision_relative(self):
        r"""
        Return the relative precision of this element.

        This is the power of `p` modulo which the unit part of this
        element is known.

        For unbounded nonzero elements, this methods return `+\infty`.

        EXAMPLES::

            sage: R = ZpER(5, prec=10)
            sage: a = R(20/21)
            sage: a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.precision_relative()
            +Infinity

            sage: b = a.add_bigoh(10)
            sage: b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + O(5^10)
            sage: b.precision_relative()
            9

        The relative precision of (exact and inexact) `0` is `0`::

            sage: x = R(0); x
            0
            sage: x.precision_relative()
            0

            sage: y = R(0, 10); y
            O(5^10)
            sage: y.precision_relative()
            0

        TESTS::

            sage: s = R.unknown()
            sage: (1/s).precision_relative()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known
        """
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if self._precbound >= maxordp and self._valuation < maxordp:
            return Infinity
        return Integer(self._precrel)

    def precision_current(self):
        r"""
        Return the internal absolute precision at which this relaxed `p`-adic
        number is known at the current stage of the computation.

        EXAMPLES::

            sage: R = ZpER(5, prec=10)
            sage: x = R(20/21)
            sage: y = R(21/22)
            sage: z = x + y

        When the elements are just defined, the computation has not started::

            sage: x.precision_current()
            0
            sage: y.precision_current()
            0
            sage: z.precision_current()
            0

        When elements are printed, the relevant digits are computed::

            sage: x
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: x.precision_current()
            10

        If we ask for more digits of `z`, the current precision of `z`
        increases accordingly::

            sage: z[:15]
            3 + 2*5 + 2*5^3 + 5^4 + 2*5^5 + 2*5^6 + 4*5^7 + 5^9 + 3*5^10 + 3*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + O(5^15)
            sage: z.precision_current()
            15

        and similarly the current precision of `x` and `y` increases because
        those digits are needed to carry out the computation::

            sage: x.precision_current()
            15
            sage: y.precision_current()
            15
        """
        if self._valuation <= -maxordp:
            return -Infinity
        if self._valuation >= maxordp:
            return Infinity
        return Integer(self._valuation + self._precrel)

    def at_precision_absolute(self, prec=None, permissive=None):
        r"""
        Return this element bounded at the given precision.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``);
          if ``None``, use the default precision of the parent

        - ``permissive`` -- boolean (default: ``False`` if ``prec``
          is given, ``True`` otherwise); if ``False``, raise an error
          if the precision of this element is not sufficient

        EXAMPLES::

            sage: R = ZpER(7, prec=10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a.at_precision_absolute(5)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)
            sage: a.at_precision_absolute(10)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + O(7^10)

            sage: b = a.add_bigoh(5)
            sage: b.at_precision_absolute(10)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision

            sage: b.at_precision_absolute(10, permissive=True)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)

        Bounding at a negative precision is not permitted over `\ZZ_p`::

            sage: a.at_precision_absolute(-1)
            Traceback (most recent call last):
            ...
            ValueError: precision must be nonnegative

        but, of course, it is over `\QQ_p`::

            sage: K = R.fraction_field()
            sage: K(a).at_precision_absolute(-1)
            O(7^-1)

        .. SEEALSO::

            :meth:`at_precision_relative`, :meth:`add_bigoh`
        """
        if prec is Infinity:
            if not permissive and self._precbound < maxordp:
                raise_error(ERROR_PRECISION)
            return self
        if prec is None:
            if permissive is None:
                permissive = True
            prec = self._parent.default_prec()
        else:
            prec = Integer(prec)
            if not self.prime_pow.in_field and prec < 0:
                raise ValueError("precision must be nonnegative")
            if prec > maxordp:
                raise OverflowError("beyond maximal precision (which is %s)" % maxordp)
        error = self._jump_c(prec)
        if permissive is None:
            permissive = False
        raise_error(error, permissive)
        return element_class_bound((<RelaxedElement>self)._parent, self, prec)

    def add_bigoh(self, absprec):
        r"""
        Return a new element with absolute precision decreased to ``absprec``.

        INPUT:

        - ``absprec`` -- integer or infinity

        EXAMPLES::

            sage: R = ZpER(7, prec=10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a.add_bigoh(5)
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)

        When ``absprec`` is negative, we return an element in the fraction
        field::

            sage: b = a.add_bigoh(-1)
            sage: b
            O(7^-1)
            sage: b.parent()
            7-adic Field handled with relaxed arithmetics

        .. SEEALSO::

            :meth:`at_precision_absolute`, :meth:`at_precision_relative`
        """
        if absprec < 0 and (not self.prime_pow.in_field):
            self = element_class_bound(self._parent.fraction_field(), self)
        return self.at_precision_absolute(absprec, True)

    def at_precision_relative(self, prec=None, halt=True, permissive=None):
        r"""
        Return this element bounded at the given precision.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``);
          if ``None``, use the default precision of the parent

        - ``halt`` -- integer or boolean (default: ``True``);
          the absolute precision after which the computation is abandoned
          if the first significant digit has not been found yet.
          If ``True``, the default halting precision of the parent is used.
          If ``False``, the computation is never abandoned.

        - ``permissive`` -- boolean (default: ``False`` if ``prec``
          is given, ``True`` otherwise); if ``False``, raise an error
          if the precision of this element is not sufficient

        EXAMPLES::

            sage: R = ZpER(5, prec=10, halt=10)
            sage: a = R(20/21); a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.at_precision_relative(5)
            4*5 + 4*5^2 + 5^4 + O(5^6)

        We illustrate the behaviour of the parameter ``halt``.
        We create a very small number whose first significant is far beyond
        the default precision::

            sage: b = R(5^20)
            sage: b
            0 + ...

        Without any help, Sage does not run the computation far enough to determine
        the valuation and an error is raised::

            sage: b.at_precision_relative(5)
            Traceback (most recent call last):
            ...
            PrecisionError: computation has been abandoned; try to increase precision

        By setting the argument ``halt``, one can force the computation to continue
        until a prescribed limit::

            sage: b.at_precision_relative(5, halt=20)   # not enough to find the valuation
            Traceback (most recent call last):
            ...
            PrecisionError: computation has been abandoned; try to increase precision

            sage: b.at_precision_relative(5, halt=21)   # now, we're okay
            5^20 + O(5^25)

        .. NOTE::

            It is also possible to pass in ``halt=False`` but it is not recommended
            because the computation can hang forever if this element is `0`.

        .. SEEALSO::

            :meth:`at_precision_absolute`, :meth:`add_bigoh`

        TESTS::

            sage: a.at_precision_relative(-1)
            Traceback (most recent call last):
            ...
            ValueError: precision must be nonnegative
        """
        if prec is Infinity:
            if not permissive and self._precbound < maxordp:
                raise_error(ERROR_PRECISION)
            return self
        default_prec = self._parent.default_prec()
        if prec is None:
            if permissive is None:
                permissive = True
            prec = default_prec
        else:
            prec = Integer(prec)
            if prec < 0:
                raise ValueError("precision must be nonnegative")
            if prec > maxordp:
                raise OverflowError("beyond maximal precision (which is %s)" % maxordp)
        if halt is True:
            halt = self._parent.halting_prec()
        elif halt is False:
            halt = maxordp
        else:
            halt = min(maxordp, halt)
        error = self._jump_relative_c(prec, halt)
        if permissive is None:
            permissive = False
        raise_error(error, permissive)
        return element_class_bound((<RelaxedElement>self)._parent, self, self._valuation + prec)

    def lift_to_precision(self, absprec=None):
        """
        Return another element of the same parent, lifting this element
        and having absolute precision at least ``absprec``.

        INPUT:

        - ``absprec`` -- integer or ``None`` (default: ``None``); the
          absolute precision of the result. If ``None``, the default
          precision of the parent is used.

        EXAMPLES::

            sage: R = ZpER(5, prec=10)
            sage: a = R(20/21, 5); a
            4*5 + 4*5^2 + 5^4 + O(5^5)

            sage: a.lift_to_precision(20)
            4*5 + 4*5^2 + 5^4 + O(5^20)

        When the precision is omitted, the default precision of the parent
        is used::

            sage: a.lift_to_precision()
            4*5 + 4*5^2 + 5^4 + O(5^10)

        When the parent is a field, the behaviour is slightly different since
        the default precision of the parent becomes the relative precision
        of the lifted element::

            sage: K = R.fraction_field()
            sage: K(a).lift_to_precision()
            4*5 + 4*5^2 + 5^4 + O(5^11)

        Note that the precision never decreases::

            sage: a.lift_to_precision(2)
            4*5 + 4*5^2 + 5^4 + O(5^5)

        In particular, unbounded element are not affected by this method::

            sage: b = R(20/21); b
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: b.lift_to_precision()
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: b.lift_to_precision(2)
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
        """
        if self._precbound >= maxordp:
            return self
        cdef long prec
        cdef long default_prec = self._parent.default_prec()
        if absprec is None:
            if self.prime_pow.in_field:
                self._jump_relative_c(1, self._precbound)
                if self._precrel == 0:
                    return self._parent.zero()
                prec = self._valuation + default_prec
            else:
                prec = default_prec
        else:
            if absprec > maxordp:
                raise OverflowError("beyond the maximal precision (which is %s)" % maxordp)
            prec = absprec
        if prec <= self._precbound:
            return self
        cdef RelaxedElement ans = element_class_slice(self._parent, self, -maxordp, self._precbound, 0)
        ans._precbound = prec
        ans._init_jump()
        return ans

    cdef long valuation_c(self, long halt=-maxordp) noexcept:
        r"""
        Return the best lower bound we have on the valuation of
        this element at the current stage of the computation.

        INPUT:

        - ``halt`` -- integer; if given, allow to increase the
          absolute precision on this element up to ``halt`` in order
          to get a better lower bound
        """
        cdef int error = 0
        while not error and self._precrel == 0 and self._valuation < halt:
            error = self._next_c()
        return self._valuation

    def valuation(self, halt=True, secure=None):
        r"""
        Return the valuation of this element.

        INPUT:

        - ``halt`` -- integer or boolean (default: ``True``);
          the absolute precision after which the computation is abandoned
          if the first significant digit has not been found yet.
          If ``True``, the default halting precision of the parent is used.
          If ``False``, the computation is never abandoned.

        - ``secure`` -- boolean (default: the value given at the creation
          of the parent); when the valuation cannot be determined for sure,
          raise an error if ``secure`` is ``True``, return the best known
          lower bound on the valuation otherwise

        EXAMPLES::

            sage: R = ZpER(5, prec=10, halt=10)
            sage: a = R(2001); a
            1 + 5^3 + 3*5^4 + ...
            sage: a.valuation()
            0

            sage: b = a - 1/a; b
            2*5^3 + 5^4 + 5^5 + 4*5^6 + 3*5^7 + 4*5^8 + 3*5^9 + 3*5^10 + 3*5^11 + 5^12 + ...
            sage: b.valuation()
            3

        The valuation of an exact zero is `+\infty`::

            sage: R(0).valuation()
            +Infinity

        The valuation of an inexact zero is its absolute precision::

            sage: R(0, 20).valuation()
            20

        We illustrate the behaviour of the parameter ``halt``.
        We create a very small number whose first significant is far beyond
        the default precision::

            sage: z = R(5^20)
            sage: z
            0 + ...

        Without any help, Sage does not run the computation far enough to determine
        the valuation and outputs only a lower bound::

            sage: z.valuation()
            10

        With ``secure=True``, an error is raised::

            sage: z.valuation(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: cannot determine the valuation; try to increase the halting precision

        By setting the argument ``halt``, one can force the computation to continue
        until a prescribed limit::

            sage: z.valuation(halt=15)   # not enough to find the correct valuation
            15
            sage: z.valuation(halt=20, secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: cannot determine the valuation; try to increase the halting precision

            sage: z.valuation(halt=21)   # now, we're okay
            20

        .. NOTE::

            It is also possible to pass in ``halt=False`` but it is not recommended
            because the computation can hang forever if this element is `0`.

        TESTS::

            sage: x = R.unknown()
            sage: (~x).valuation()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known
        """
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if halt is True:
            halt = self._parent.halting_prec()
        elif halt is False:
            halt = maxordp
        else:
            halt = min(maxordp, halt)
        cdef val = self.valuation_c(halt)
        if secure is None:
            secure = self._parent.is_secure()
        if secure and self._precbound >= maxordp and self._precrel == 0:
            raise PrecisionError("cannot determine the valuation; try to increase the halting precision")
        return Integer(val)

    def unit_part(self, halt=True):
        r"""
        Return the unit part of this element.

        INPUT:

        - ``halt`` -- integer or boolean (default: ``True``);
          the absolute precision after which the computation is abandoned
          if the first significant digit has not been found yet.
          If ``True``, the default halting precision of the parent is used.
          If ``False``, the computation is never abandoned.

        EXAMPLES::

            sage: R = ZpER(5, prec=10, halt=10)
            sage: a = R(20/21); a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.unit_part()
            4 + 4*5 + 5^3 + 4*5^5 + 3*5^6 + 4*5^7 + 5^9 + ...

            sage: b = 1/a; b
            4*5^-1 + 4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + ...
            sage: b.unit_part()
            4 + 4*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + ...

        The unit part of `0` is not defined::

            sage: R(0).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

            sage: R(0, 20).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

        See :meth:`valuation` for more details on the parameter ``halt``.
        """
        val = self.valuation(halt)
        if self._valuation >= self._precbound:
            raise ValueError("unit part of 0 not defined")
        return self >> val

    def val_unit(self, halt=True):
        r"""
        Return the valuation and the unit part of this element.

        INPUT:

        - ``halt`` -- integer or boolean (default: ``True``);
          the absolute precision after which the computation is abandoned
          if the first significant digit has not been found yet.
          If ``True``, the default halting precision of the parent is used.
          If ``False``, the computation is never abandoned.

        EXAMPLES::

            sage: R = ZpER(5, 10)
            sage: a = R(20/21); a
            4*5 + 4*5^2 + 5^4 + 4*5^6 + 3*5^7 + 4*5^8 + ...
            sage: a.val_unit()
            (1, 4 + 4*5 + 5^3 + 4*5^5 + 3*5^6 + 4*5^7 + 5^9 + ...)

            sage: b = 1/a; b
            4*5^-1 + 4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + ...
            sage: b.val_unit()
            (-1, 4 + 4*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + ...)

        If this element is indistinguishable from zero, an error is raised
        since the unit part of `0` is not defined::

            sage: R(0).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

            sage: R(0, 20).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined

        See :meth:`valuation` for more details on the parameter ``halt``.
        """
        val = self.valuation(halt)
        if self._valuation >= self._precbound:
            raise ValueError("unit part of 0 not defined")
        return val, self >> val

    def residue(self, absprec=1, field=True, check_prec=True):
        r"""
        Return the image of this element in the quotient
        `\ZZ/p^\mathrm{absprec}\ZZ`.

        INPUT:

        - ``absprec`` -- nonnegative integer (default: 1)

        - ``field`` -- boolean (default: ``True``); when ``absprec`` is ``1``,
          whether to return an element of GF(p) or Zmod(p)

        - ``check_prec`` -- boolean (default: ``True``); whether to raise an error
          if this element has insufficient precision to determine the reduction

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: a = R(1/2021); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + 7^5 + 6*7^6 + 3*7^7 + 6*7^8 + 5*7^9 + ...
            sage: a.residue()
            3
            sage: a.residue(2)
            45

        If this element has negative valuation, an error is raised::

            sage: K = R.fraction_field()
            sage: b = K(20/21)
            sage: b.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue
        """
        if absprec >= maxordp:
            raise OverflowError
        if absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of p")
        error = self._jump_c(absprec)
        raise_error(error, not check_prec)
        if self._valuation < 0:
            raise ValueError("element must have nonnegative valuation in order to compute residue")
        cdef celement digits
        cdef Integer ans
        if absprec <= self._valuation:
            ans = ZZ(0)
        else:
            self._getslice_relative(digits, 0, min(self._precrel, absprec - self._valuation))
            ans = element_get_sage(digits, self.prime_pow) * self.prime_pow(self._valuation)
        if field and absprec == 1:
            return self._parent.residue_class_field()(ans)
        else:
            return self._parent.residue_ring(absprec)(ans)

    def lift(self, absprec=None):
        r"""
        Return a rational number which is congruent to this element modulo
        `p^\mathrm{prec}`.

        INPUT:

        - ``absprec`` -- integer or ``None`` (default: ``None``); if ``None``,
          the absolute precision of this element is used

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: a = R(1/2021, 5); a
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)
            sage: a.lift()
            15676
            sage: a.lift(2)
            45

        Here is another example with an element of negative valuation::

            sage: K = R.fraction_field()
            sage: b = K(20/21, 5); b
            2*7^-1 + 3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + O(7^5)
            sage: b.lift()
            39223/7

        For unbounded elements, we must specify a precision::

            sage: c = R(1/2021)
            sage: c.lift()
            Traceback (most recent call last):
            ...
            ValueError: you must specify a precision for unbounded elements

            sage: c.lift(5)
            15676
        """
        if absprec is None:
            if self._precbound < maxordp:
                absprec = self._precbound
            else:
                raise ValueError("you must specify a precision for unbounded elements")
        else:
            absprec = Integer(absprec)
        cdef int error = self._jump_c(absprec)
        raise_error(error)
        if absprec < self._valuation:
            return Integer(0)
        cdef celement digits
        self._getslice_relative(digits, 0, absprec - self._valuation)
        ans = element_get_sage(digits, self.prime_pow)
        if self._valuation:
            ans *= self._parent.prime() ** self._valuation
        return ans

    def __rshift__(self, s):
        r"""
        Return this element divided by `\pi^s`, and truncated
        if the parent is not a field.

        EXAMPLES::

            sage: R = ZpER(997)
            sage: K = R.fraction_field()
            sage: a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + ...

        Shifting to the right divides by a power of `p`, but drops
        terms with negative valuation::

            sage: a >> 3
            124 + ...

        If the parent is a field no truncation is performed::

            sage: K(a) >> 3
            964*997^-2 + 572*997^-1 + 124 + ...

        A negative shift multiplies by that power of `p`::

            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6 + ...
        """
        cdef long start
        cdef long shift = long(s)
        if shift:
            if (<RelaxedElement>self)._parent.is_field():
                start = -maxordp
            else:
                start = shift
            return element_class_slice((<pAdicRelaxedElement>self)._parent, self, start, maxordp, shift)
        else:
            return self

    def __lshift__(self, s):
        r"""
        Return this element multiplied by `\pi^s`.

        If `s` is negative and this element does not lie in a field,
        digits may be truncated.  See ``__rshift__`` for details.

        EXAMPLES::

            sage: R = ZpER(997)
            sage: K = R.fraction_field()
            sage: a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + ...

        Shifting to the right divides by a power of `p`, but drops
        terms with negative valuation::

            sage: a << 2
            964*997^3 + 572*997^4 + 124*997^5 + ...

        A negative shift may result in a truncation when the base
        ring is not a field::

            sage: a << -3
            124 + ...

            sage: K(a) << -3
            964*997^-2 + 572*997^-1 + 124 + ...
        """
        return self.__rshift__(-s)

    cpdef _add_(self, other):
        r"""
        Return the sum of this element with ``other``.

        TESTS::

            sage: R = ZpER(7, 10)
            sage: R(1/3) + R(1/6)
            4 + 3*7 + 3*7^2 + 3*7^3 + 3*7^4 + 3*7^5 + 3*7^6 + 3*7^7 + 3*7^8 + 3*7^9 + ...

            sage: R(1/3, 5) + R(1/6, 10)
            4 + 3*7 + 3*7^2 + 3*7^3 + 3*7^4 + O(7^5)
        """
        if isinstance(self, RelaxedElement_zero):
            return other
        if isinstance(other, RelaxedElement_zero):
            return self
        return element_class_add(self._parent, self, <RelaxedElement>other)

    cpdef _sub_(self, other):
        r"""
        Return the difference of this element and ``other``.

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: R(1/3) - R(1/6)
            6 + 5*7 + 5*7^2 + 5*7^3 + 5*7^4 + 5*7^5 + 5*7^6 + 5*7^7 + 5*7^8 + 5*7^9 + ...

            sage: R(1/3, 5) - R(1/6, 10)
            6 + 5*7 + 5*7^2 + 5*7^3 + 5*7^4 + O(7^5)
        """
        if self is other:
            ans = self._parent.zero()
            if self._precbound < maxordp:
                ans = element_class_bound(self._parent, ans, self._precbound)
            return ans
        if isinstance(other, RelaxedElement_zero):
            return self
        return element_class_sub(self._parent, self, <RelaxedElement>other)

    cpdef _neg_(self):
        r"""
        Return the opposite of this element.

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: -R(1)
            6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + 6*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + 6*7^9 + ...

            sage: -R(1,5)
            6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)
        """
        if isinstance(self, RelaxedElement_zero):
            return self
        return element_class_sub(self._parent, self._parent.zero(), self)

    cpdef _mul_(self, other):
        r"""
        Return the product of this element with ``other``.

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: R(1/2) * R(2/3)
            5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + 4*7^6 + 4*7^7 + 4*7^8 + 4*7^9 + ...

            sage: R(1/2, 5) * R(2/3, 10)
            5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + O(7^5)

            sage: R(49, 5) * R(14, 4)
            2*7^3 + O(7^6)
        """
        if isinstance(self, RelaxedElement_zero) or isinstance(other, RelaxedElement_one):
            return self
        if isinstance(self, RelaxedElement_one) or isinstance(other, RelaxedElement_zero):
            return other
        return element_class_mul(self._parent, self, <RelaxedElement>other)

    cpdef _div_(self, other):
        r"""
        Return the quotient if this element by ``other``.

        .. NOTE::

            The result always lives in the fraction field, even if ``other`` is
            a unit.

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: x = R(2) / R(3)
            sage: x
            3 + 2*7 + 2*7^2 + 2*7^3 + 2*7^4 + 2*7^5 + 2*7^6 + 2*7^7 + 2*7^8 + 2*7^9 + ...
            sage: x.parent()
            7-adic Field handled with relaxed arithmetics

        TESTS::

            sage: x / R(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by zero

            sage: x / R(0, 10)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by something indistinguishable from zero

            sage: y = R.unknown()
            sage: x / y
            O(7^-Infinity)
        """
        if isinstance(other, RelaxedElement_one):
            if self.prime_pow.in_field:
                return self
            else:
                return element_class_bound(self._parent.fraction_field(), self)
        return element_class_div(self._parent.fraction_field(), self, <RelaxedElement>other, -maxordp)

    def __invert__(self):
        r"""
        Return the multiplicative inverse of this element.

        .. NOTE::

            The result always lives in the fraction field, even if this element
            is a unit.

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: x = ~R(3)
            sage: x
            5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + 4*7^6 + 4*7^7 + 4*7^8 + 4*7^9 + ...
            sage: x.parent()
            7-adic Field handled with relaxed arithmetics

        TESTS::

            sage: ~R(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by zero

            sage: ~R(0, 10)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by something indistinguishable from zero

            sage: y = R.unknown()
            sage: ~y
            O(7^-Infinity)
        """
        if isinstance(self, RelaxedElement_one):
            return self
        return element_class_div(self._parent.fraction_field(), self._parent.one(), self, -maxordp)

    def inverse_of_unit(self):
        r"""
        Return the multiplicative inverse of this element if
        it is a unit.

        EXAMPLES::

            sage: R = ZpER(3, 5)
            sage: a = R(2)
            sage: b = a.inverse_of_unit()
            sage: b
            2 + 3 + 3^2 + 3^3 + 3^4 + ...

        A :exc:`ZeroDivisionError` is raised if an element has no inverse in the
        ring::

            sage: R(3).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: denominator is not invertible

        Unlike the usual inverse of an element, the result is in the same ring
        as this element and not in its fraction field (for fields this does of
        course not make any difference)::

            sage: c = ~a; c
            2 + 3 + 3^2 + 3^3 + 3^4 + ...
            sage: a.parent()
            3-adic Ring handled with relaxed arithmetics
            sage: b.parent()
            3-adic Ring handled with relaxed arithmetics
            sage: c.parent()
            3-adic Field handled with relaxed arithmetics

        This method also works for self-referent numbers
        (see :meth:`sage.rings.padics.generic_nodes.pAdicRelaxedGeneric.unknown`)::

            sage: x = R.unknown(); x
            O(3^0)
            sage: x.inverse_of_unit()
            O(3^0)
            sage: x.set(1 + 3 * x.inverse_of_unit())
            True
            sage: x
            1 + 3 + 2*3^2 + 3^3 + 3^4 + ...

        Actually, in many cases, it is preferable to use it than an actual
        division. Indeed, compare::

            sage: y = R.unknown()
            sage: y.set(1 + 3/y)
            Traceback (most recent call last):
            ...
            RecursionError: definition looks circular
        """
        if isinstance(self, RelaxedElement_one):
            return self
        if self.prime_pow.in_field:
            return element_class_div(self._parent, self._parent.one(), self, -maxordp)
        else:
            return element_class_div(self._parent, self._parent.one(), self, 0)

    def sqrt(self):
        r"""
        Return the square root of this element.

        EXAMPLES::

            sage: R = ZpER(7, 10)
            sage: x = R(8)
            sage: x.sqrt()
            1 + 4*7 + 2*7^2 + 7^3 + 3*7^4 + 2*7^5 + 4*7^6 + 2*7^7 + 5*7^8 + ...

        When the element is not a square, an error is raised::

            sage: x = R(10)
            sage: x.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: not a square

        For bounded elements, the precision is tracked::

            sage: x = R(8, 5); x
            1 + 7 + O(7^5)
            sage: x.sqrt()
            1 + 4*7 + 2*7^2 + 7^3 + 3*7^4 + O(7^5)

        Note that, when `p = 2`, a digit of precision is lost::

            sage: S = ZpER(2)
            sage: x = S(17, 5)
            sage: x.sqrt()
            1 + 2^3 + O(2^4)

        This method also work for self-referent numbers
        (see :meth:`sage.rings.padics.generic_nodes.pAdicRelaxedGeneric.unknown`)::

            sage: x = R.unknown(); x
            O(7^0)
            sage: x.sqrt()
            O(7^0)
            sage: x.set(1 + 7*sqrt(x))
            True
            sage: x
            1 + 7 + 4*7^2 + 4*7^3 + 2*7^4 + 3*7^8 + 3*7^9 + ...

        TESTS::

            sage: for p in [ 7, 11, 1009 ]:
            ....:     R = ZpER(p)
            ....:     x = 1 + p * R.random_element()
            ....:     y = x.sqrt()
            ....:     assert(x == y^2)
        """
        return element_class_sqrt(self._parent, self)

    def _test_pickling(self, **options):
        r"""
        Check that this object can be pickled and unpickled properly.

        TESTS::

            sage: R = ZpER(7)
            sage: x = R.random_element()
            sage: x._test_pickling()

            sage: x[:20]._test_pickling()

        .. SEEALSO::

            :func:`dumps`, :func:`loads`
        """
        tester = self._tester(**options)
        from sage.misc.persist import loads, dumps
        if self._precbound >= maxordp:
            tester.assertEqual(loads(dumps(self)), self.at_precision_relative())
        else:
            tester.assertEqual(loads(dumps(self)), self)

    def _test_nonzero_equal(self, **options):
        r"""
        Test that ``.__bool__()`` behave consistently with `` == 0``.

        TESTS::

            sage: R = ZpER(5)
            sage: #R(0)._test_nonzero_equal()
            sage: #R(5^30)._test_nonzero_equal()
            sage: #R.unknown()._test_nonzero_equal()
        """
        tester = self._tester(**options)
        try:
            tester.assertEqual(self != self.parent().zero(), bool(self))
            tester.assertEqual(self == self.parent().zero(), not self)
        except PrecisionError:
            pass


cdef class RelaxedElement_abandon(RelaxedElement):
    r"""
    A special class for relaxed `p`-adic with all digits unknown.

    This class is used for setting temporary definition of
    some self-referent numbers.
    """
    def __init__(self):
        r"""
        Initialize this element.

        TESTS::

            sage: from sage.rings.padics.padic_relaxed_element import RelaxedElement_abandon
            sage: x = RelaxedElement_abandon()
            sage: x.valuation()
            Traceback (most recent call last):
            ...
            PrecisionError: no lower bound on the valuation is known

            sage: x[0]
            Traceback (most recent call last):
            ...
            PrecisionError: computation has been abandoned; try to increase precision
        """
        self._valuation = -maxordp

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this element.

        Here, we just abandon the computation.
        """
        return ERROR_ABANDON

cdef relaxedelement_abandon = RelaxedElement_abandon()


cdef class RelaxedElementWithDigits(RelaxedElement):
    r"""
    A generic class for relaxed `p`-adic elements that stores
    the sequence of its digits.
    """
    def __cinit__(self):
        r"""
        Allocate memory for storing digits.
        """
        element_init(self._digits)

    def __dealloc__(self):
        r"""
        Deallocate memory used for digits.
        """
        element_clear(self._digits)

    cdef cdigit_ptr _getdigit_relative(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in relative precision.
        """
        return element_get_digit(self._digits, i)

    cdef cdigit_ptr _getdigit_absolute(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in absolute precision.
        """
        return element_get_digit(self._digits, i - self._valuation)

    cdef void _getslice_relative(self, celement slice, long start, long length) noexcept:
        r"""
        Select a slice of the sequence of digits of this element.

        INPUT:

        - ``slice`` -- a ``celement`` to store the slice

        - ``start`` -- integer; the start position of the slice

        - ``length`` -- integer; the length of the slice

        .. NOTE::

            This function only sets up a pointer to the requested slice
            (the slice is not copied). Hence any future modification
            of the slice will modify this element as well.
        """
        element_get_slice(slice, self._digits, start, length)


# Assignment
############

# Zero

cdef class RelaxedElement_zero(RelaxedElement):
    r"""
    A class for representation a relaxed `p`-adic number which is
    exactly zero.

    TESTS::

        sage: R = ZpER(7)
        sage: a = R.zero()
        sage: TestSuite(a).run()
    """
    def __init__(self, parent):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(0)  # indirect doctest
            sage: x
            0

            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_zero'>
        """
        RelaxedElement.__init__(self, parent)
        self._valuation = maxordp

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: a = ZpER(5)(0)
            sage: type(a)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_zero'>
            sage: loads(dumps(a)) == a   # indirect doctest
            True
        """
        return self.__class__, (self._parent,)

    cdef cdigit_ptr _getdigit_relative(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in relative precision.
        """
        return digit_zero

    cdef cdigit_ptr _getdigit_absolute(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in absolute precision.
        """
        return digit_zero

    cdef void _getslice_relative(self, celement slice, long start, long length) noexcept:
        r"""
        Select a slice of the sequence of digits of this element.

        INPUT:

        - ``slice`` -- a ``celement`` to store the slice

        - ``start`` -- integer; the start position of the slice

        - ``length`` -- integer; the length of the slice
        """
        element_init(slice)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        return 0

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        return 0


# One

cdef class RelaxedElement_one(RelaxedElementWithDigits):
    r"""
    A class for representation a relaxed `p`-adic number which is
    exactly one.

    TESTS::

        sage: R = ZpER(7)
        sage: a = R.one()
        sage: TestSuite(a).run()
    """
    def __init__(self, parent):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(1)  # indirect doctest
            sage: x
            1 + ...

            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_one'>
        """
        RelaxedElement.__init__(self, parent)
        element_set_digit_ui(self._digits, 1, 0)
        self._precrel = 1

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: a = ZpER(5)(1)
            sage: type(a)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_one'>
            sage: a[:20] == loads(dumps(a))   # indirect doctest
            True
        """
        return self.__class__, (self._parent,)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        if prec > self._precrel:
            self._precrel = prec
        return 0

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        self._precrel += 1
        return 0


# Bound

cdef class RelaxedElement_bound(RelaxedElement):
    r"""
    A class for `p`-adic relaxed elements which are defined by bounding the
    precision of another `p`-adic relaxed element.

    TESTS::

        sage: R = ZpER(5)
        sage: x = R.random_element()
        sage: y = x[:20]
        sage: TestSuite(y).run()
    """
    def __init__(self, parent, RelaxedElement x, precbound=None):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adics, the element to bound

        - ``precbound`` -- integer or ``None`` (default: ``None``);
          the bound on the precision

        .. NOTE::

            The digits of ``x`` are not copied!

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(20/21)
            sage: y = x.add_bigoh(20)
            sage: type(y)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_bound'>
        """
        RelaxedElement.__init__(self, parent)
        self._x = x
        if precbound is None:
            self._precbound = x._precbound
        else:
            self._precbound = min(x._precbound, precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: a = ZpER(5)(1).add_bigoh(20)
            sage: type(a)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_bound'>
            sage: a == loads(dumps(a))   # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._x, self._precbound)

    cdef cdigit_ptr _getdigit_relative(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in relative precision.
        """
        return self._x._getdigit_relative(i)

    cdef cdigit_ptr _getdigit_absolute(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in absolute precision.
        """
        return self._x._getdigit_absolute(i)

    cdef void _getslice_relative(self, celement slice, long start, long length) noexcept:
        r"""
        Select a slice of the digits of this number.

        INPUT:

        - ``slice`` -- a ``celement`` to store the slice

        - ``start`` -- positive integer; the starting position of the slice
          in relative precision

        - ``length`` -- positive integer; the length of the slice

        .. NOTE::

            This methods only sets up a pointer to the requested slice
            (the slice is not copied). Hence any future modification
            of the slice ``slice`` will affect this number.
        """
        self._x._getslice_relative(slice, start, length)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Jump to the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef RelaxedElement x = self._x
        cdef int error
        if prec > self._precbound:
            error = ERROR_PRECISION | x._jump_c(self._precbound)
        else:
            error = x._jump_c(prec)
        self._precbound = min(self._precbound, x._precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        return error

    cdef int _next_c(self) noexcept:
        r"""
        Jump to the next digit.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef RelaxedElement x = self._x
        if self._valuation + self._precrel >= self._precbound:
            return ERROR_PRECISION
        cdef int error = x._next_c()
        self._precbound = min(self._precbound, x._precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        return error


# Value

cdef class RelaxedElement_value(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adics defined by the datum of a value in
    the exact subring.

    TESTS::

        sage: R = ZpER(5)
        sage: x = R(2)
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, value, long shift=0, precbound=None):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``value`` -- the value in the exact subring

        - ``shift`` -- integer (default: `0`); the position at which
          the given value is written

        - ``precbound`` -- integer or ``None`` (default: ``None``);
          the bound on the precision

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(2)
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_value'>

            sage: y = R(2, 10)
            sage: type(y)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_value'>
        """
        RelaxedElement.__init__(self, parent)
        element_set_digit_sage(self._digits, value, 0)
        self._value = value
        self._shift = shift
        self._valuation = -shift
        if precbound is not None and precbound is not Infinity:
            self._precbound = min(maxordp, precbound)
        self._valuebound = maxordp
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: a = ZpER(5)(2, 20)
            sage: type(a)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_value'>
            sage: a == loads(dumps(a))   # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._value, self._shift, self._precbound)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        if self._valuebound >= maxordp:
            return RelaxedElement._jump_c(self, prec)
        if (self._precbound is not None) and (prec > self._precbound):
            self._precrel = self._precbound - self._valuation
            return ERROR_PRECISION
        cdef long precrel = min(prec, maxordp) - self._valuation
        if precrel > self._precrel:
            self._precrel = precrel
        if prec >= maxordp:
            return ERROR_OVERFLOW
        return 0

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).

        .. NOTE::

            The algorithm is not optimal (quadratic in the precision),
            but it sounds okay...
        """
        if (self._precbound is not None) and (self._valuation + self._precrel >= self._precbound):
            return ERROR_PRECISION
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        if digit_is_zero(self._getdigit_relative(self._precrel)):
            self._valuebound = self._valuation + self._precrel
            if self._precrel == 0:
                self._valuation = self._precbound
            elif self._precbound < maxordp:
                self._precrel = self._precbound - self._valuation
        return 0


# Random

cdef class RelaxedElement_random(RelaxedElementWithDigits):
    r"""
    A class for random relaxed `p`-adic numbers.

    TESTS::

        sage: R = ZpER(5)
        sage: x = R.random_element()
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, valuation, precbound=None, seed=None):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``valuation`` -- integer or ``None``; the position from which
          random digits are picked;
          if ``None``, it is randomly chosen if the parent is a field and
          set to `0` otherwise

        - ``precbound`` -- integer or ``None`` (default: ``None``);
          the bound on the precision

        - ``seed`` -- integer or ``None`` (default: ``None``); the
          seed of the random generator

        .. NOTE::

            The argument ``valuation`` can be different from the real
            valuation of this number since the first randomly picked
            digit could vanish.

        TESTS::

            sage: R = ZpER(7)
            sage: x = R.random_element()
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_random'>
        """
        RelaxedElement.__init__(self, parent)
        if seed is None:
            self._seed = randint(0, 2*maxordp)
        else:
            self._seed = seed
        digit_random_init(self._generator, self._seed)
        if valuation is None:
            self._initialvaluation = ZZ.random_element()
        else:
            self._initialvaluation = valuation
        self._valuation = self._initialvaluation
        if precbound is not None:
            self._precbound = min(maxordp, precbound)
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5, print_mode='digits')
            sage: a = R.random_element()
            sage: a   # random
            ...32220241412003314311

            sage: b = loads(dumps(a))
            sage: b   # random
            ...32220241412003314311

        It is guaranteed that `a` and `b` are equal at any precision::

            sage: a[:30]   # random
            ...?343214211432220241412003314311
            sage: b[:30]   # random
            ...?343214211432220241412003314311

            sage: a == b
            True
        """
        return self.__class__, (self._parent, self._initialvaluation, self._precbound, self._seed)

    cdef int _next_c(self) noexcept:
        r"""
        Generate the next digit of this number at random.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef cdigit r
        digit_init(r)
        digit_random(r, self.prime_pow, self._generator)
        if self._precrel == 0 and digit_is_zero(r):
            self._valuation += 1
        else:
            element_set_digit(self._digits, r, self._precrel)
            self._precrel += 1
        digit_clear(r)
        return 0


# Operations
############

# Slice and shift

cdef class RelaxedElement_slice(RelaxedElement):
    r"""
    A class for relaxed `p`-adic numbers defined as slices.

    TESTS::

        sage: R = ZpER(5)
        sage: x = R(20/21)
        sage: y = x.slice(3, 6)
        sage: TestSuite(y).run()
    """
    def __init__(self, parent, RelaxedElement x, long start, long stop, long shift):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adic element, the element from which the
          slice is extracted

        - ``start`` -- integer; the position of the first digit of `x`
          in the slice

        - ``stop`` -- integer; the position of the first digit of `x`
          after the slice

        - ``shift`` -- integer such that ``self[i] = x[i+shift]``

        .. NOTE::

            The digits of ``x`` are not copied!

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(20/21)
            sage: y = x.slice(3, 6)
            sage: type(y)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_slice'>
        """
        RelaxedElement.__init__(self, parent)
        self._x = x
        self._start = start
        self._shift = shift
        self._valuation = max(x._valuation, start) - shift
        self._precrel = min(x._precrel + x._valuation, stop) - self._valuation - shift
        if x._precbound < stop:
            self._precbound = min(maxordp, x._precbound - shift)
        self._stop = stop
        if self._precrel < 0:
            self._precrel = 0
        while self._precrel > 0 and digit_is_zero(self._getdigit_relative(0)):
            self._precrel -= 1
            self._valuation += 1
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5, print_mode='digits')
            sage: x = R(20/21)
            sage: y = x.slice(3, 6)
            sage: y == loads(dumps(y))  # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._x, self._start, self._stop, self._shift)

    cdef cdigit_ptr _getdigit_relative(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in relative precision.
        """
        return self._getdigit_absolute(i + self._valuation)

    cdef cdigit_ptr _getdigit_absolute(self, long i) noexcept:
        r"""
        Return a pointer on the `i`-th digit of this number
        in absolute precision.
        """
        cdef long j = i + self._shift
        if j < self._start or j >= self._stop:
            return digit_zero
        else:
            return self._x._getdigit_absolute(j)

    cdef void _getslice_relative(self, celement slice, long start, long length) noexcept:
        r"""
        Select a slice of the sequence of digits of this element.

        INPUT:

        - ``slice`` -- a ``celement`` to store the slice

        - ``start`` -- integer; the start position of the slice

        - ``length`` -- integer; the length of the slice

        .. NOTE::

            This function only sets up a pointer to the requested slice
            (the slice is not copied). Hence any future modification
            of the slice will modify this element as well.
        """
        cdef RelaxedElement x = self._x
        cdef long s = start + self._valuation + self._shift
        cdef long start_absolute = max(self._start, s)
        cdef long stop_absolute = min(self._stop, s + length)
        x._getslice_relative(slice, start_absolute - x._valuation, stop_absolute - start_absolute)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Jump to the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef RelaxedElement x = self._x
        cdef int error = 0
        cdef long pr
        if prec <= self._valuation + self._precrel:
            return 0
        if prec > self._precbound:
            prec = self._precbound
            error = ERROR_PRECISION
        cdef int errorx = x._jump_c(min(prec + self._shift, self._stop))
        pr = max(self._valuation, x._valuation + x._precrel - self._shift)
        if self._precrel == 0:
            while self._valuation < pr and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
        if errorx:
            self._precrel = pr - self._valuation
            return errorx
        else:
            self._precrel = prec - self._valuation
            return error

    cdef int _next_c(self) noexcept:
        r"""
        Jump to the next digit.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef int error
        cdef long n = self._precrel + self._valuation
        n += self._shift
        if n <= self._stop:
            error = self._x._jump_c(n+1)
            if error:
                return error
        if self._precrel == 0 and (n > self._stop or digit_is_zero(self._getdigit_relative(0))):
            self._valuation += 1
        else:
            self._precrel += 1
        return 0


# Addition

cdef class RelaxedElement_add(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as sums.

    TESTS::

        sage: R = ZpER(11)
        sage: x = R.random_element() + R.random_element()
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, RelaxedElement x, RelaxedElement y):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adic element, the first summand

        - ``y`` -- a relaxed `p`-adic element, the second summand

        TESTS::

            sage: R = ZpER(11)
            sage: x = R.random_element() + R.random_element()
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_add'>
        """
        RelaxedElement.__init__(self, parent)
        self._x = x
        self._y = y
        cdef halt = self._parent.default_prec()
        self._valuation = min(x.valuation_c(0), y.valuation_c(0))
        self._precbound = min(x._precbound, y._precbound)
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R.random_element() + R.random_element()
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._x, self._y)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        # We reimplement _jump_c for better performance
        cdef long n = self._valuation + self._precrel
        cdef RelaxedElement x = self._x
        cdef RelaxedElement y = self._y
        cdef int error = x._jump_c(prec) | y._jump_c(prec)
        prec = min(prec, x._valuation + x._precrel, y._valuation + y._precrel)
        while n < prec:
            element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
            element_iadd_digit(self._digits, y._getdigit_absolute(n), self._precrel)
            element_reducesmall_digit(self._digits, self._precrel, self.prime_pow)
            if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
                element_shift_right(self._digits)
            else:
                self._precrel += 1
            n += 1
        return error

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef long n = self._valuation + self._precrel
        cdef RelaxedElement x = self._x
        cdef RelaxedElement y = self._y
        cdef int error = x._jump_c(n+1) | y._jump_c(n+1)
        if error:
            return error
        element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
        element_iadd_digit(self._digits, y._getdigit_absolute(n), self._precrel)
        element_reducesmall_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        return 0


# Subtraction

cdef class RelaxedElement_sub(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as differences.

    TESTS::

        sage: R = ZpER(11)
        sage: x = R.random_element() - R.random_element()
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, RelaxedElement x, RelaxedElement y):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adic element, the minuend

        - ``y`` -- a relaxed `p`-adic element, the subtrahend

        TESTS::

            sage: R = ZpER(11)
            sage: x = R.random_element() - R.random_element()
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_sub'>
        """
        RelaxedElement.__init__(self, parent)
        self._x = x
        self._y = y
        self._valuation = min(x.valuation_c(0), y.valuation_c(0))
        self._precbound = min(x._precbound, y._precbound)
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R.random_element() - R.random_element()
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._x, self._y)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        # We reimplement _jump_c for better performances
        cdef long n = self._valuation + self._precrel
        cdef RelaxedElement x = self._x
        cdef RelaxedElement y = self._y
        cdef int error = x._jump_c(prec) | y._jump_c(prec)
        prec = min(prec, x._valuation + x._precrel, y._valuation + y._precrel)
        while n < prec:
            element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
            element_isub_digit(self._digits, y._getdigit_absolute(n), self._precrel)
            element_reduceneg_digit(self._digits, self._precrel, self.prime_pow)
            if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
                element_shift_right(self._digits)
            else:
                self._precrel += 1
            n += 1
        return error

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef long n = self._valuation + self._precrel
        cdef RelaxedElement x = self._x
        cdef RelaxedElement y = self._y
        cdef int error = x._jump_c(n+1) | y._jump_c(n+1)
        if error:
            return error
        element_iadd_digit(self._digits, x._getdigit_absolute(n), self._precrel)
        element_isub_digit(self._digits, y._getdigit_absolute(n), self._precrel)
        element_reduceneg_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        return 0


# Multiplication

cdef class RelaxedElement_mul(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as products.

    ALGORITHM:

    We compute digits using relaxed arithmetic by var der Hoeven et al.,
    whose cost is quasi-linear with respect to the precision.

    The algorithm uses the entries behind the current position in the table
    ``self._digits`` to store carries.

    TESTS::

        sage: R = ZpER(11)
        sage: x = R.random_element() * R.random_element()
        sage: TestSuite(x).run()
    """
    def __cinit__(self):
        r"""
        Allocate memory for temporary variables.
        """
        digit_init(self._lastdigit_x)
        digit_init(self._lastdigit_y)

    def __dealloc__(self):
        r"""
        Deallocate memory for temporary variables.
        """
        digit_clear(self._lastdigit_x)
        digit_clear(self._lastdigit_y)

    def __init__(self, parent, RelaxedElement x, RelaxedElement y):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adic element, the first factor

        - ``y`` -- a relaxed `p`-adic element, the second factor

        TESTS::

            sage: R = ZpER(11)
            sage: x = R.random_element() * R.random_element()
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_mul'>
        """
        RelaxedElement.__init__(self, parent)
        self._x = x
        self._y = y
        cdef halt = self._parent.default_prec()
        self._valuation = min(maxordp, x.valuation_c(0) + y.valuation_c(0))
        if x._precbound < maxordp:
            y._jump_relative_c(1, y._valuation + self._parent.default_prec())
            self._precbound = min(self._precbound, y._valuation + x._precbound)
        if y._precbound < maxordp:
            x._jump_relative_c(1, x._valuation + self._parent.default_prec())
            self._precbound = min(self._precbound, x._valuation + y._precbound)
        self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R.random_element() * R.random_element()
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._x, self._y)

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        global tmp_digit, tmp_poly
        cdef RelaxedElement x = self._x
        cdef RelaxedElement y = self._y
        cdef long n = self._valuation + self._precrel

        cdef int errorx = x._jump_c(n - y._valuation + 1)
        cdef int errory = y._jump_c(n - x._valuation + 1)
        cdef int error = errorx | errory
        if self._precrel == 0:
            self._valuation = x._valuation + y._valuation
            if self._valuation > n:
                return 0
            if self._valuation < n or x._precrel == 0 or y._precrel == 0:
                return error | ERROR_PRECISION
        elif error:
            return error

        n = self._precrel
        digit_set(self._lastdigit_x, x._getdigit_relative(n))
        digit_set(self._lastdigit_y, y._getdigit_relative(n))
        digit_mul(tmp_digit, x._getdigit_relative(0), self._lastdigit_y)
        element_iadd_digit(self._digits, tmp_digit, n)
        if n:
            digit_mul(tmp_digit, self._lastdigit_x, y._getdigit_relative(0))
            element_iadd_digit(self._digits, tmp_digit, n)

        cdef long m = n + 2
        cdef long len = 1
        cdef celement slicex, slicey
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
            x._getslice_relative(slicex, len - 1, len)
            y._getslice_relative(slicey, (m-1)*len - 1, len)
            element_mul(tmp_poly, slicex, slicey)
            element_iadd_slice(self._digits, tmp_poly, n)
            if m > 2:
                x._getslice_relative(slicex, (m-1)*len - 1, len)
                y._getslice_relative(slicey, len - 1, len)
                element_mul(tmp_poly, slicex, slicey)
                element_iadd_slice(self._digits, tmp_poly, n)

        element_reduce_digit(self._digits, n, self.prime_pow)
        self._precrel += 1
        return 0

    cdef int _update_last_digit(self) noexcept:
        r"""
        Redo the computation of the last digit and update carries
        accordingly.

        This method is used for computing Teichmüller representatives.
        """
        if self._precrel == 0:
            return ERROR_UNEXPECTED
        cdef RelaxedElement x = self._x
        cdef RelaxedElement y = self._y
        cdef long n = self._precrel - 1
        cdef long m = n + 2
        cdef long len = 2
        cdef celement slice
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
        len -= 1

        digit_sub(tmp_digit, x._getdigit_relative(n), self._lastdigit_x)
        y._getslice_relative(slice, 0, len)
        element_scalarmul(tmp_poly, slice, tmp_digit)
        element_iadd_slice(self._digits, tmp_poly, n)
        if m == 2:
            len -= 1
        digit_sub(tmp_digit, y._getdigit_relative(n), self._lastdigit_y)
        x._getslice_relative(slice, 0, len)
        element_scalarmul(tmp_poly, slice, tmp_digit)
        element_iadd_slice(self._digits, tmp_poly, n)
        if m == 2:
            digit_mul(tmp_digit, tmp_digit, self._lastdigit_x)
            element_iadd_digit(self._digits, tmp_digit, 2*len)
        element_reduce_digit(self._digits, n, self.prime_pow)

        digit_set(self._lastdigit_x, x._getdigit_relative(n))
        digit_set(self._lastdigit_y, y._getdigit_relative(n))
        return 0


cdef class RelaxedElement_muldigit(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as products
    of a relaxed `p`-adic number by a digit.

    This class is not exposed to the user; it is only used
    internally for division.
    """
    def __init__(self, parent, RelaxedElement_div x, RelaxedElement y):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adic element, whose first significant
          digit is the first factor

        - ``y`` -- a relaxed `p`-adic element, the second factor
        """
        RelaxedElement.__init__(self, parent)
        self._x = <cdigit_ptr>x._inverse
        self._y = y
        self._valuation = y._valuation
        self._init_jump()

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef long n = self._valuation + self._precrel
        cdef int error = self._y._jump_c(n+1)
        if error:
            return error
        digit_mul(tmp_digit, self._x, self._y._getdigit_absolute(n))
        element_iadd_digit(self._digits, tmp_digit, self._precrel)
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        self._precrel += 1
        return 0


# Division

cdef class RelaxedElement_div(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as quotients.

    ALGORITHM:

    We compute the quotient `x = a/b` as the self-referent number defined by

    .. MATH::

        x = ac + (1 - bc) x

    where `c` is congruent to `b^{-1}` modulo the uniformizer.

    TESTS::

        sage: R = ZpER(5)
        sage: x = R(20) / R(21)
        sage: TestSuite(x).run()
    """
    def __cinit__(self):
        r"""
        Allocate memory for temporary variables.
        """
        digit_init(self._inverse)

    def __dealloc__(self):
        r"""
        Deallocate memory for temporary variables.
        """
        digit_clear(self._inverse)

    def __init__(self, parent, RelaxedElement num, RelaxedElement denom, long minval=-maxordp, precbound=None):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``num`` -- a relaxed `p`-adic element, the dividend

        - ``denom`` -- a relaxed `p`-adic element, the divisor

        - ``minval`` -- integer; the minimal valuation allowed for this element

        - ``precbound`` -- integer or ``None`` (default: ``None``);
          the bound on the precision

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(20) / R(21)
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_div'>

            sage: y = R.unknown()
            sage: 1/y
            O(5^-Infinity)
            sage: y.inverse_of_unit()
            O(5^0)
        """
        RelaxedElement.__init__(self, parent)
        if denom._valuation >= maxordp:
            raise ZeroDivisionError("cannot divide by zero")
        if denom._precbound < maxordp and denom._precrel == 0:
            raise ZeroDivisionError("cannot divide by something indistinguishable from zero")
        self._num = num
        self._denom = denom
        if denom._valuation <= -maxordp:
            self._maxprec = maxordp + 1
        else:
            self._maxprec = denom._valuation + max(1, self._parent.default_prec())
        self._valuation = minval
        cdef int error = self._bootstrap_c()
        if precbound is not None:
            self._precbound = min(maxordp, precbound)
        if num._precbound < maxordp:
            self._precbound = min(self._precbound, num._precbound - denom._valuation)
        if denom._precbound < maxordp:
            self._precbound = min(self._precbound, denom._precbound + num._valuation - 2*denom._valuation)
        raise_error(error, permissive=True)
        if not error:
            self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(20) / R(21)
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._num, self._denom, self._valuation, self._precbound)

    cdef int _bootstrap_c(self) noexcept:
        r"""
        Bootstrap the computation of the digits of this element, that is:

        - find the valuation
        - compute the first digit, and
        - set up the recursive definition of the next digits.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef int error
        cdef RelaxedElement num = self._num
        cdef RelaxedElement denom = self._denom

        while denom._valuation < self._maxprec and denom._precrel == 0:
            error = denom._next_c()
            if error:
                if error & ERROR_PRECISION:
                    error |= ERROR_DIVISION
                return error
        if self._maxprec < maxordp and denom._valuation > -maxordp:
            self._maxprec = denom._valuation + max(1, self._parent.default_prec())
        if denom._precrel == 0:
            return ERROR_ABANDON

        cdef long valuation = num._valuation - denom._valuation
        if valuation < self._valuation:
            error = num._jump_c(self._valuation + denom._valuation)
            if error:
                return error
            valuation = num._valuation - denom._valuation
            if valuation < self._valuation:
                return ERROR_DIVISION
        self._valuation = valuation
        digit_inv(self._inverse, denom._getdigit_relative(0), self.prime_pow)
        self._definition = relaxedelement_abandon
        cdef parent = self._parent
        cdef RelaxedElement a = element_class_muldigit(parent, self, num)
        cdef RelaxedElement b = element_class_muldigit(parent, self, denom)
        cdef RelaxedElement c = element_class_slice(parent, b, denom._valuation + 1, maxordp, 0)
        cdef RelaxedElement d = element_class_mul(parent, c, self)
        self._definition = element_class_sub(parent, a, d)
        return 0

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef RelaxedElement definition = self._definition
        if definition is None:
            return self._bootstrap_c()
        cdef long val = self._valuation + self._denom._valuation
        cdef int error = definition._jump_c(val + self._precrel + 1)
        if error:
            return error
        if definition._valuation > val:
            self._valuation = min(self._precbound, definition._valuation - self._denom._valuation)
            if definition._precbound < maxordp:
                self._precbound = min(self._precbound, definition._precbound - self._denom._valuation)
        else:
            digit = definition._getdigit_relative(self._precrel)
            element_set_digit(self._digits, digit, self._precrel)
            self._precrel += 1
        return 0


# Square root

cdef class RelaxedElement_sqrt(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as square roots.

    ALGORITHM:

    When `p \neq 2`, we compute `y = \sqrt{x}` as the self-referent number
    defined by

    .. MATH::

        y = \frac{x - (y-a)^2 + a^2}{2a}

    where `a^2` is congruent to `x` modulo the uniformizer.

    When `p = 2`, we use a variant of this construction.

    TESTS::

        sage: R = ZpER(5)
        sage: x = R(6).sqrt()
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, RelaxedElement x):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``x`` -- a relaxed `p`-adic element

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(6).sqrt()
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_sqrt'>
        """
        RelaxedElement.__init__(self, parent)
        self._x = x
        if x._valuation <= -maxordp:
            self._valuation = -maxordp
        else:
            self._valuation = x._valuation >> 1
        cdef int error = self._bootstrap_c()
        if error & ERROR_NOTSQUARE:
            raise ValueError("not a square")
        if not error:
            if x._precbound < maxordp:
                self._precbound = x._precbound - x._valuation / 2
                if self._parent.prime() == 2:
                    self._precbound -= 1
                self._precbound = min(maxordp, self._precbound)
            self._init_jump()

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(5)
            sage: x = R(6).sqrt()
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        return self.__class__, (self._parent, self._x)

    cdef int _bootstrap_c(self) noexcept:
        r"""
        Bootstrap the computation of the digits of this element, that is:

        - find the valuation
        - compute the first digit, and
        - set up the recursive definition of the next digits.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).

        .. NOTE::

            This code does not work for nontrivial extensions of `\QQ_2`.
        """
        cdef RelaxedElement x = self._x
        cdef long maxprec
        if x._valuation <= -maxordp:
            return ERROR_ABANDON
        if self._valuation <= -maxordp:
            maxprec = x._valuation
        else:
            maxprec = (self._valuation + 1) << 1
        while x._valuation < maxprec and x._precrel == 0:
            error = x._next_c()
            if error:
                return error
        if x._valuation & 1:
            self._valuation = (x._valuation + 1) >> 1
        else:
            self._valuation = x._valuation >> 1
        if x._precrel == 0:
            return 0
        if x._valuation & 1 != 0:
            return ERROR_NOTSQUARE

        cdef parent = self._parent
        cdef long val = self._valuation
        cdef cdigit digit
        cdef Integer zd, p = self.prime_pow.prime
        cdef RelaxedElement u, y, c, d

        if p == 2:
            element_set_digit_ui(self._digits, 1, 0)
            self._precrel = 1
            if x._precrel == 1:
                error = x._next_c()
                if error:
                    return error
            if not digit_equal_ui(x._getdigit_relative(1), 0):
                return ERROR_NOTSQUARE
            if x._precrel == 2:
                error = x._next_c()
                if error:
                    return error
            if not digit_equal_ui(x._getdigit_relative(2), 0):
                return ERROR_NOTSQUARE
            zd = Integer(1)
            self._definition = relaxedelement_abandon
            u = element_class_slice(parent, self, val + 2, maxordp, val)
            y = element_class_slice(parent, x, -maxordp, maxordp, val + 1)
            c = element_class_value(parent, zd, shift=-val+1)
            d = element_class_slice(parent, u*u, -maxordp, maxordp, -val + 1)
            self._definition = y + c - d
        else:
            digit_init(digit)
            if digit_sqrt(digit, x._getdigit_relative(0), self.prime_pow):
                digit_clear(digit)
                return ERROR_NOTSQUARE
            element_set_digit(self._digits, digit, 0)
            self._precrel = 1
            zd = digit_get_sage(digit)
            self._definition = relaxedelement_abandon
            u = element_class_slice(parent, self, val + 1, maxordp, val)
            y = element_class_slice(parent, x, -maxordp, maxordp, 2*val)
            c = element_class_value(parent, zd*zd)
            d = element_class_value(parent, 2*zd, shift=val)
            self._definition = (y + c - u*u) / d
        return 0

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef RelaxedElement definition = self._definition
        if definition is None:
            return self._bootstrap_c()
        cdef long n = self._valuation + self._precrel
        cdef int error = definition._jump_c(n+1)
        if error:
            return error
        element_set_digit(self._digits, definition._getdigit_relative(self._precrel), self._precrel)
        self._precrel += 1
        return 0


# Teichmüller lifts

cdef class RelaxedElement_teichmuller(RelaxedElementWithDigits):
    r"""
    A class for relaxed `p`-adic numbers defined as teichmüller representatives.

    ALGORITHM:

    We compute `x = [a]` as the unique self-referent number with last
    digit `a` and `x = x^p`.
    Note that `x^p` is known with one more digit than `x` itself.

    TESTS::

        sage: R = ZpER(7)
        sage: x = R.teichmuller(2)
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, xbar):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``xbar`` -- an element in the exact subring, which is congruent
          to this Teichmüller modulo this uniformizer

        TESTS::

            sage: R = ZpER(7)
            sage: x = R.teichmuller(2)
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_teichmuller'>
        """
        RelaxedElement.__init__(self, parent)
        cdef cdigit digit
        digit_init(digit)
        digit_set_sage(digit, xbar)
        digit_mod(digit, digit, self.prime_pow)
        if digit_equal_ui(digit, 0):
            digit_clear(digit)
            self._trivial = True
            self._valuation = maxordp
        else:
            element_set_digit(self._digits, digit, 0)
            self._trivial = digit_equal_ui(digit, 1)
        self._precrel += 1

        cdef RelaxedElement xn
        cdef Integer p
        cdef int size, i
        if not self._trivial:
            xn = self
            p = self.prime_pow.prime
            size = mpz_sizeinbase(p.value, 2)
            i = size - 2
            self._xns = []
            while i >= 0:
                xn = element_class_mul(parent, xn, xn)
                self._xns.append(xn)
                if mpz_tstbit(p.value, i):
                    xn = element_class_mul(parent, xn, self)
                    self._xns.append(xn)
                i -= 1
            self._xp = xn
        self._ready = True

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(7)
            sage: x = R.teichmuller(2)
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        xbar = digit_get_sage(element_get_digit(self._digits, 0))
        return self.__class__, (self._parent, xbar)

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Compute the digits of this number until the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        if not self._ready:
            return ERROR_ABANDON
        if self._trivial:
            if self._valuation == 0 and self._precrel < prec:
                self._precrel = prec
            return 0
        return RelaxedElement._jump_c(<RelaxedElement>self, prec)

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        if self._trivial:
            if self._valuation:
                self._precrel += 1
            return 0
        cdef int error
        cdef RelaxedElement xp = self._xp
        cdef RelaxedElement_mul xn
        self._precrel += 1
        xp._jump_c(self._precrel)
        element_set_digit(self._digits, xp._getdigit_relative(self._precrel - 1), self._precrel - 1)
        for xn in self._xns:
            error = xn._update_last_digit()
            if error:
                return error | ERROR_UNEXPECTED
        return 0


# Self-referent definitions
###########################

cdef class RelaxedElement_unknown(RelaxedElementWithDigits):
    r"""
    A class for self-referent relaxed `p`-adic numbers.

    TESTS::

        sage: R = ZpER(7)
        sage: x = R.unknown()
        sage: TestSuite(x).run()

        sage: x.set(1 + 7*x^2)
        True
        sage: TestSuite(x).run()
    """
    def __init__(self, parent, long valuation, digits=None):
        r"""
        Initialize this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``valuation`` -- integer; a lower bound on the valuation of
          this number

        - ``digits`` -- list or ``None`` (default: ``None``); the first
          significant digits of this number

        TESTS::

            sage: R = ZpER(7)
            sage: x = R.unknown()
            sage: type(x)
            <class 'sage.rings.padics.padic_relaxed_element.pAdicRelaxedElement_unknown'>
        """
        RelaxedElement.__init__(self, parent)
        if valuation >= maxordp:
            raise OverflowError("valuation is too large (maximum is %s)" % maxordp)
        self._valuation = valuation
        self._definition = None
        self._next = maxordp
        cdef cdigit digit
        if digits is not None:
            digits = [Integer(d) for d in digits]
            for d in digits:
                digit_set_sage(digit, d)
                element_iadd_digit(self._digits, digit, self._precrel)
                element_reduce_digit(self._digits, self._precrel, self.prime_pow)
                if self._precrel == 0 and digit_is_zero(element_get_digit(self._digits, 0)):
                    self._valuation += 1
                    element_shift_right(self._digits)
                else:
                    self._precrel += 1
        self._initialvaluation = self._valuation
        self._initialprecrel = self._precrel

    def __reduce__(self):
        r"""
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: R = ZpER(7)
            sage: x = R.unknown()
            sage: x.set(1 + 7*x^2)
            True
            sage: x == loads(dumps(x))  # indirect doctest
            True
        """
        digits = []
        for i in range(self._initialprecrel):
            digits.append(digit_get_sage(element_get_digit(self._digits, i)))
        definition = None
        if id(self) not in persist.already_pickled:
            persist.already_pickled[id(self)] = True
            definition = self._definition
        return unpickle_unknown, (id(self), self.__class__, self._parent, self._initialvaluation, digits, definition)

    cpdef set(self, RelaxedElement definition):
        r"""
        Set the recursive definition of this self-referent number.

        INPUT:

        - ``definition`` -- a relaxed `p`-adic number, to which this
          number is equal

        OUTPUT:

        A boolean indicating if the definition is coherent with the
        already known digits of this number.

        EXAMPLES::

            sage: R = ZpER(5, 10)
            sage: x = R.unknown()
            sage: x.set(1 + 5*x)
            True
            sage: x
            1 + 5 + 5^2 + 5^3 + 5^4 + 5^5 + 5^6 + 5^7 + 5^8 + 5^9 + ...

        The previous construction works because the relation we gave defines
        the `n`-th digit of `x` in terms of its digits at precision strictly
        less than `n` (this is due to the multiplication by `5`).

        On the contrary, the following does not work::

            sage: y = R.unknown()
            sage: y.set(1 + 3*y)
            True
            sage: y
            O(5^0)
            sage: y[:20]
            Traceback (most recent call last):
            ...
            RecursionError: definition looks circular

        In the next example, we give explicit values for the first digits
        and then a recursive definition for the next digits. However, the
        recursive definition does not hold for the first digits; that is the
        reason why the call to :meth:`set` returns ``False``::

            sage: z = R.unknown(digits=[2])
            sage: z
            2 + O(5)
            sage: z.set(1 + 5*z)
            False
            sage: z
            2 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + ...

        .. SEEALSO::

            :meth:`sage.rings.padics.generic_nodes.pAdicRelaxedGeneric.unknown`
        """
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = definition
        self._precbound = max(self._valuation + self._precrel, definition._precbound)
        eq = self._is_equal(definition, self._valuation + self._precrel, True)
        self._init_jump()
        return eq

    cdef int _next_c(self) noexcept:
        r"""
        Compute the next digit of this number.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        cdef RelaxedElement definition = self._definition
        cdef cdigit_ptr digit
        cdef long n = self._valuation + self._precrel
        cdef long diffval
        cdef int error

        if definition is None:
            return ERROR_NOTDEFINED
        if n >= self._next:
            return ERROR_CIRCULAR

        cdef long svenext = self._next
        self._next = n
        error = definition._jump_c(n+1)
        if not error:
            digit = definition._getdigit_absolute(n)
            if self._precrel == 0 and digit_is_zero(digit):
                self._valuation += 1
            else:
                element_set_digit(self._digits, digit, self._precrel)
                self._precrel += 1
        self._next = svenext
        return error


def unpickle_unknown(uid, cls, parent, valuation, digits, definition):
    r"""
    Unpickle a self-referent relaxed `p`-adic number.

    TESTS:

    Cross definitions involving several self-referent numbers are
    handled correctly::

        sage: R = ZpER(7)
        sage: x = R.unknown()
        sage: y = R.unknown()
        sage: x.set(1 + 2*y + 7*x*y)
        True
        sage: y.set(3 + 14*x^2)
        True

        sage: x == loads(dumps(x))  # indirect doctest
        True
        sage: y == loads(dumps(y))  # indirect doctest
        True
    """
    if uid in persist.already_unpickled:
        elt = persist.already_unpickled[uid]
    else:
        elt = cls(parent, valuation, digits)
        persist.already_unpickled[uid] = elt
    if definition is not None:
        elt.set(definition)
    return elt


# Expansion
###########

cdef class RelaxedElement_zeroone(RelaxedElementWithDigits):
    r"""
    A special class for `p`-adic relaxed elements with only
    `0` and `1` as digits.

    This class is used for computing expansion in Teichmuller mode.
    It is not supposed to be instantiated in other situations.
    """
    def __init__(self, parent, long valuation):
        r"""
        Instantiate this element.

        INPUT:

        - ``parent`` -- the parent of this element

        - ``valuation`` -- the valuation of this number
        """
        RelaxedElement.__init__(self, parent)
        self._valuation = valuation

    cdef void _setdigit_to_zero(self) noexcept:
        r"""
        Append `0` to the list of digits of this element.
        """
        self._precrel += 1

    cdef void _setdigit_to_one(self) noexcept:
        r"""
        Append `1` to the list of digits of this element.
        """
        element_set_digit_ui(self._digits, 1, self._precrel)
        self._precrel += 1

    cdef int _jump_c(self, long prec) noexcept:
        r"""
        Jump to the absolute precision ``prec``.

        INPUT:

        - ``prec`` -- integer

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        if prec > self._valuation + self._precrel:
            return ERROR_NOTDEFINED
        return 0

    cdef int _next_c(self) noexcept:
        r"""
        Jump to the next digit.

        OUTPUT:

        An error code (see :meth:`RelaxedElement._next_c` for details).
        """
        return ERROR_NOTDEFINED


cdef class ExpansionIter():
    """
    An iterator over a `p`-adic expansion.

    This class should not be instantiated directly, but instead using
    :meth:`RelaxedElement.expansion`.
    """
    def __cinit__(self):
        r"""
        Allocate memory for temporary variables.
        """
        digit_init(self.carry)

    def __dealloc__(self):
        r"""
        Deallocate memory for temporary variables.
        """
        digit_clear(self.carry)

    def __init__(self, RelaxedElement elt, expansion_mode mode, long start, long stop):
        r"""
        Initialize this iterator.

        INPUT:

        - ``elt`` -- a relaxed `p`-adic number

        - ``mode`` -- either ``simple_mode``, ``smallest_mode`` or ``teichmuller_mode``

        - ``start`` -- integer; the position where the expansion starts

        - ``stop`` -- integer; the position where the expansion stops

        TESTS::

            sage: E = ZpER(5,4)(373).expansion()
            sage: I = iter(E)   # indirect doctest
            sage: type(I)
            <class 'sage.rings.padics.padic_relaxed_element.ExpansionIter'>
        """
        self.elt = elt
        self.mode = mode
        self.start = start
        self.stop = stop
        self.current = min(start, elt._valuation)
        digit_init(self.digit)
        # Compute first digits if needed
        if self.mode == simple_mode:
            self.current = self.start
        elif self.mode == smallest_mode:
            while self.current < self.start:
                self._next_smallest()
        elif self.mode == teichmuller_mode:
            self.tail = elt
            self.coefficients = {}
            while self.current < self.start:
                self._next_teichmuller()

    def __repr__(self):
        r"""
        Return a string representation of this iterator.

        EXAMPLES::

            sage: R = ZpER(7, 5)
            sage: x = R(1/2021)
            sage: x.expansion()   # indirect doctest
            7-adic expansion of 3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + ...

            sage: x.expansion(lift_mode='smallest')   # indirect doctest
            7-adic expansion of 3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + ... (balanced)

            sage: x.expansion(lift_mode='teichmuller')   # indirect doctest
            7-adic expansion of 3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + ... (teichmuller)
        """
        s = "%s-adic expansion of %s" % (self.elt._parent.prime(), self.elt)
        if self.mode == smallest_mode:
            s += " (balanced)"
        elif self.mode == teichmuller_mode:
            s += " (teichmuller)"
        return s

    def __len__(self):
        r"""
        Return the length of this expansion.

        EXAMPLES::

            sage: R = ZpER(7)
            sage: x = R(1/2021, 5)
            sage: x
            3 + 6*7 + 4*7^2 + 3*7^3 + 6*7^4 + O(7^5)
            sage: E = x.expansion()
            sage: len(E)
            5

        For unbounded elements, the expansion is infinite and this method
        raises an error::

            sage: y = R(1/2021)
            sage: E = y.expansion()
            sage: len(E)
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite sequence
        """
        if self.stop >= maxordp:
            raise NotImplementedError("infinite sequence")
        return Integer(self.stop - self.start)

    def __iter__(self):
        r"""
        Return itself (as any iterator is supposed to do).

        TESTS::

            sage: E = ZpER(5)(373).expansion()
            sage: I = iter(E)
            sage: I is iter(I)
            True
        """
        return self

    cdef _next_simple(self):
        r"""
        Return the next digit of this expansion (simple mode).
        """
        cdef RelaxedElement elt = self.elt
        elt._jump_c(self.current + 1)
        digit_set(self.digit, elt._getdigit_absolute(self.current))
        self.current += 1
        return digit_get_sage(self.digit)

    cdef _next_smallest(self):
        r"""
        Return the next digit of this expansion (smallest mode).
        """
        cdef RelaxedElement elt = self.elt
        elt._jump_c(self.current + 1)
        digit_add(self.digit, elt._getdigit_absolute(self.current), self.carry)
        digit_smallest(self.digit, self.carry, self.digit, elt.prime_pow)
        self.current += 1
        return digit_get_sage(self.digit)

    cdef _next_teichmuller(self):
        r"""
        Return the next digit of this expansion (Teichmüller mode).
        """
        cdef RelaxedElement teichmuller, tail = self.tail
        cdef RelaxedElement_zeroone coeff
        cdef digit
        tail._jump_c(self.current + 1)
        digit_set(self.digit, tail._getdigit_absolute(self.current))
        digit = digit_get_sage(self.digit)
        if digit != 0 and digit != 1 and digit not in self.coefficients:
            parent = tail._parent
            self.coefficients[digit] = coeff = RelaxedElement_zeroone(parent, self.current)
            teichmuller = element_class_teichmuller(parent, digit)
            self.tail = tail - coeff * element_class_slice(parent, teichmuller, 1, maxordp, 0)
        for d, coeff in self.coefficients.items():
            if d == digit:
                coeff._setdigit_to_one()
            else:
                coeff._setdigit_to_zero()
        self.current += 1
        return digit

    def __next__(self):
        r"""
        Return the next digit of this expansion.

        EXAMPLES::

            sage: R = ZpER(11, 10)
            sage: x = R(20/21); x
            2 + 2*11 + 4*11^2 + 8*11^3 + 5*11^4 + 11^6 + 2*11^7 + 4*11^8 + 8*11^9 + ...
            sage: E = x.expansion()
            sage: next(E)
            2
            sage: next(E)
            2

        TESTS::

            sage: def check_expansion(x, mode):
            ....:     R = x.parent()
            ....:     E = x.expansion(lift_mode=mode)
            ....:     y = 0
            ....:     for i in range(len(E)):
            ....:         digit = next(E)
            ....:         if mode == 'teichmuller':
            ....:             y += R.teichmuller(digit) << i
            ....:         else:
            ....:             y += R(digit) << i
            ....:     assert(x == y)

            sage: for p in primes(100):
            ....:     x = ZpER(p).random_element()[:20]
            ....:     for mode in [ 'simple', 'smallest', 'teichmuller' ]:
            ....:         check_expansion(x, mode)
        """
        if self.current >= self.stop:
            raise StopIteration
        if self.mode == simple_mode:
            return self._next_simple()
        elif self.mode == smallest_mode:
            return self._next_smallest()
        elif self.mode == teichmuller_mode:
            return self._next_teichmuller()
