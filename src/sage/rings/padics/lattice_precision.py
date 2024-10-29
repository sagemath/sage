r"""
Lattice precision for the parents ``ZpLC``/``QpLC`` and ``ZpLF``/``QpLF``

AUTHOR:

- Xavier Caruso (2018-02): initial version

TESTS::

    sage: R = ZpLC(2)
    doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
    See https://github.com/sagemath/sage/issues/23505 for details.
    sage: prec = R.precision()
    sage: prec
    Precision lattice on 0 objects

    sage: S = ZpLF(2)
    sage: prec = S.precision()
    sage: prec
    Precision module on 0 objects
"""

# ****************************************************************************
#       Copyright (C) 2018 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections import defaultdict

from sage.misc.timing import walltime

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

from sage.rings.padics.precision_error import PrecisionError


# The default minimal size after which re-echelonization is not performed,
# i.e., when a variable is not referenced anymore and could be deleted but its
# corresponding column is further than this threshold from the right end of the
# matrix representing the precision lattice, then the column is not removed
# from the matrix because the re-echelonization would be too costly.
DEFAULT_THRESHOLD_DELETION = 50

# The number of additional digits used for internal computations
STARTING_ADDITIONAL_PREC = 5


class pRational:
    r"""
    This class implements rational numbers viewed as elements of ``Qp``.
    In particular, it provides additional methods which are specific to
    `p`-adics (as `p`-adic valuation).

    Only for internal use.

    INPUT:

    - ``p`` -- a prime number

    - ``x`` -- a rational number

    - ``exponent`` -- integer (default: 0)

    - ``valuation`` -- integer or ``None`` (default: ``None``);
      the `p`-adic valuation of this element

    If not ``None``, this method trusts the given value to the
    attribute ``valuation``.

    TESTS::

        sage: from sage.rings.padics.lattice_precision import pRational
        sage: x = pRational(2, 5); x
        5
        sage: y = pRational(2, 5/3, 2); y
        2^2 * 5/3

        sage: x + y
        35/3
        sage: x - y
        -5/3
        sage: x * y
        2^2 * 25/3
        sage: x / y
        2^-2 * 3

        sage: x.valuation()
        0
        sage: y.valuation()
        2

        sage: z = pRational(2, 1024, valuation=4)
        sage: z
        1024
        sage: z.valuation()
        4
    """
    def __init__(self, p, x, exponent=0, valuation=None):
        r"""
        Construct the element ``x * p^exponent``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: pRational(2, 5)
            5
            sage: pRational(2, 5/3, 2)
            2^2 * 5/3
        """
        self.p = p
        if x in ZZ:
            self.x = ZZ(x)
        else:
            self.x = x
        self.exponent = exponent
        self._valuation = valuation

    def __repr__(self):
        r"""
        Return a string representation of this element.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: pRational(2, 5, 2)  # indirect doctest
            2^2 * 5
        """
        if self.exponent == 0:
            return repr(self.x)
        else:
            return "%s^%s * %s" % (self.p, self.exponent, self.x)

    def reduce(self, prec):
        r"""
        Return this element reduced modulo ``p^prec``.

        INPUT:

        - ``prec`` -- integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 1234567); x
            1234567
            sage: x.reduce(12)
            1671

            sage: x = pRational(2, 1234/567); x
            1234/567
            sage: x.reduce(12)
            190
        """
        if prec is Infinity:
            return self
        x = self.x
        exp = self.exponent
        if x.parent() is ZZ:
            if prec > exp:
                x = x % (self.p ** (prec-exp))
            else:
                x = 0
        elif x.parent() is QQ:
            num = x.numerator()
            denom = x.denominator()
            valdenom = denom.valuation(self.p)
            denom //= self.p ** valdenom
            exp -= valdenom
            if prec > exp:
                modulo = self.p ** (prec - exp)
                # probably we should use Newton iteration instead
                # (but it is actually slower for now - Python implementation)
                _, inv, _ = denom.xgcd(modulo)
                x = (num*inv) % modulo
            else:
                x = 0
        if self.x == 0:
            val = Infinity
        else:
            val = self._valuation
        return self.__class__(self.p, x, exp, valuation=val)

    def reduce_relative(self, prec):
        r"""
        Return this element reduced modulo ``p^n`` where ``n = prec + val(x)``.

        INPUT:

        - ``prec`` -- nonnegative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 1234567); x
            1234567
            sage: x.reduce_relative(12)
            1671

            sage: x = pRational(2, 1234/567); x
            1234/567
            sage: x.reduce_relative(12)
            190
        """
        v = self.valuation()
        if v is Infinity:
            return self
        return self.reduce(prec+v)

    def normalize(self):
        r"""
        Normalize this element, i.e. write it as ``p^v * u`` where
        ``u`` is coprime to `p`.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.normalize(); x
            2^13 * 1929
        """
        if self.x == 0:
            self.exponent = 0
        else:
            val = self.valuation()
            exp = self.exponent
            self.x /= self.p ** (val-exp)
            if self.x in ZZ:
                self.x = ZZ(self.x)
            self.exponent = val

    def valuation(self):
        r"""
        Return the `p`-adic valuation of this element.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.valuation()
            13
        """
        if self._valuation is None:
            valx = self.x.valuation(self.p)
            self._valuation = self.exponent + valx
        return self._valuation

    def is_p_power(self):
        r"""
        Return ``True`` if this element is a power of `p`.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 1024, 2); x
            2^2 * 1024
            sage: x.is_p_power()
            True

            sage: y = pRational(2, 123456, 7); y
            2^7 * 123456
            sage: y.is_p_power()
            False
        """
        self.normalize()
        return self.x == 1

    def is_zero(self):
        r"""
        Return ``True`` if this element vanishes.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.is_zero()
            False

            sage: (x-x).is_zero()
            True
        """
        return self.x == 0

    def __add__(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x + y
            2^7 * 28635808
        """
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if sexp is Infinity:
            return other
        if oexp is Infinity:
            return self
        if self._valuation is None or other._valuation is None:
            val = None
        elif self._valuation < other._valuation:
            val = self._valuation
        elif self._valuation > other._valuation:
            val = other._valuation
        else:
            val = None
        if sexp < oexp:
            return self.__class__(p, self.x + other.x * p**(oexp-sexp), sexp, valuation=val)
        else:
            return self.__class__(p, self.x * p**(sexp-oexp) + other.x, oexp, valuation=val)

    def __sub__(self, other):
        r"""
        Return the subtraction of ``self`` by ``other``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x - y
            2^7 * -28388896
        """
        return self + (-other)

    def __neg__(self):
        r"""
        Return the opposite of this element.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: -x
            2^7 * -123456
        """
        return self.__class__(self.p, -self.x, self.exponent, valuation=self._valuation)

    def __mul__(self, other):
        r"""
        Return the product of ``self`` and ``other``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x * y
            2^19 * 110000654016
        """
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = self._valuation + other._valuation
        return self.__class__(self.p, self.x * other.x, self.exponent + other.exponent, valuation=val)

    def __truediv__(self, other):
        r"""
        Return the quotient of ``self`` by ``other``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x / y
            2^-5 * 123456/891011
        """
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = self._valuation - other._valuation
        return self.__class__(self.p, self.x / other.x, self.exponent - other.exponent, valuation=val)

    def _quo_rem(self, other):
        """
        Quotient with remainder.

        Returns a pair `q`, `r` where `r` has the `p`-adic expansion of this element,
        truncated at the valuation of other.

        EXAMPLES::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: a = pRational(2, 123456, 3)
            sage: b = pRational(2, 654321, 2)
            sage: q,r = a._quo_rem(b); q, r
            (2^7 * 643/218107, 0)
            sage: q*b+r - a
            0
            sage: q,r = b._quo_rem(a); q, r
            (5111/1929, 2^2 * 113)
            sage: q*a+r - b
            2^2 * 0
        """
        other.normalize()
        ox = other.x
        if ox == 0:
            raise ZeroDivisionError
        self.normalize()
        oval = other.exponent
        sx = self.x
        sval = self.exponent
        diff = sval - oval
        if sx == 0:
            return (self.__class__(self.p, 0, 0, valuation=Infinity),
                    self.__class__(self.p, 0, 0, valuation=Infinity))
        elif sval >= oval:
            return (self.__class__(self.p, sx / ox, diff, valuation=diff),
                    self.__class__(self.p, 0, 0, valuation=Infinity))
        else:
            pd = self.p**(-diff)
            sred = sx % pd
            return (self.__class__(self.p, (sx - sred)/(pd*ox), 0),
                    self.__class__(self.p, sred, sval, valuation=sval))

    def __lshift__(self, n):
        r"""
        Return the product of this element by ``p^n``.

        INPUT:

        - ``n`` -- relative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x << 10
            2^17 * 123456
        """
        if self._valuation is None:
            val = None
        else:
            val = self._valuation + n
        return self.__class__(self.p, self.x, self.exponent + n, valuation=val)

    def __rshift__(self, n):
        r"""
        Return the quotient of this element by ``p^n``.

        INPUT:

        - ``n`` -- relative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x >> 10
            2^-3 * 123456
        """
        return self << (-n)

    def unit_part(self):
        r"""
        Return the unit part of this element, that is the part ``u``
        in the writing ``u * p^v`` with ``u`` coprime to `p`.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.unit_part()
            1929
        """
        if self.is_zero():
            raise ValueError("the unit part of zero is not defined")
        p = self.p
        val = self.valuation()
        x = self.x / (p ** (val-self.exponent))
        return self.__class__(p, x, 0, valuation=0)

    def xgcd(self, other):
        r"""
        Return the gcd of ``self`` and ``other`` together with two
        elements ``u`` and ``v`` such that ``u*self + v*other = gcd``.

        The ``gcd`` is normalized so that it is a power of `p`.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011

            sage: d, u, v = x.xgcd(y)
            sage: d
            2^7 * 32
            sage: d.normalize(); d
            2^12 * 1

            sage: u*x + v*y
            2^7 * 32
        """
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if sexp < oexp:
            a = ZZ(self.x)
            b = ZZ(other.x * (p ** (oexp-sexp)))
            exp = sexp
        else:
            a = ZZ(self.x * (p ** (sexp-oexp)))
            b = ZZ(other.x)
            exp = oexp
        d, u, v = a.xgcd(b)
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = min(self._valuation, other._valuation)
        d = self.__class__(p, d, exp, valuation=val)
        u = self.__class__(p, u)
        v = self.__class__(p, v)
        return d, u, v

    def value(self):
        r"""
        Return this element as a rational number.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.value()
            15802368
        """
        return (self.p ** self.exponent) * self.x

    def list(self, prec):
        r"""
        Return the list of the digits of this element (written in radix
        `p`) up to position ``prec``.

        The first zeros are omitted.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456); x
            123456
            sage: x.list(5)
            []
            sage: x.list(20)
            [1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]

            sage: y = pRational(2, 123/456); y
            41/152
            sage: y.list(10)
            [1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1]

            sage: z = pRational(2, 0)
            sage: z.list(10)
            []
            sage: z.list(100)
            []
        """
        if self.x not in ZZ:
            self = self.reduce(prec)
        val = self.valuation()
        if val is Infinity:
            return []
        p = self.p
        x = ZZ(self.x * p**(self.exponent - val))
        l = []
        for _ in range(val, prec):
            x, digit = x.quo_rem(p)
            l.append(digit)
        return l


class DifferentialPrecisionGeneric(SageObject):
    r"""
    A generic class for precision objects obtained by automatic
    differentiation.

    INPUT:

    - ``p`` -- a prime number

    - ``label`` -- string; the label of the parents to which the elements
      belong that are tracked by this precision module

    .. NOTE::

        This object is used internally by the parent ring. You should not
        create instances of this class on your own.

    EXAMPLES::

        sage: R = ZpLC(2, label='init')
        sage: R.precision()
        Precision lattice on 0 objects (label: init)
    """
    def __init__(self, p, label):
        r"""
        TESTS::

            sage: prec = ZpLC(2, label='init').precision()
            sage: from sage.rings.padics.lattice_precision import DifferentialPrecisionGeneric
            sage: isinstance(prec, DifferentialPrecisionGeneric)
            True
        """
        self._p = p
        self._label = label
        self._elements = [ ]
        self._matrix = { } # A dictionary whose keys are weak references to tracked elements
                           # and values corresponding columns in the matrix
                           # representing the precision lattice
        self._collected_references = [ ]
        self._marked_for_deletion = [ ]
        self._approx_zero = pRational(p, ZZ(0))
        self._threshold_deletion = DEFAULT_THRESHOLD_DELETION
        self._history_init = None
        self._history = None

    def __reduce__(self):
        r"""
        TESTS::

            sage: R = ZpLF(2)
            sage: prec = R.precision()
            sage: dumps(prec)
            Traceback (most recent call last):
            ...
            NotImplementedError: pickling/unpickling precision modules is not implemented yet
        """
        raise NotImplementedError("pickling/unpickling precision modules is not implemented yet")

    def _repr_(self):
        r"""
        Return a string representation of this precision object.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: R.precision()
            Precision lattice on ... objects

        If a label has been specified, it is included in the representation::

            sage: R = ZpLC(2, label='mylabel')
            sage: R.precision()
            Precision lattice on 0 objects (label: mylabel)
        """
        label = "" if self._label is None else " (label: %s)" % (self._label,)
        count = "1 object" if len(self._elements) == 1 else "%s objects" % len(self._elements)
        return "%s on %s%s" % (self._repr_type, count, label)

    def threshold_deletion(self, threshold=None):
        r"""
        Return (and set) the threshold for column deletion.

        When a variable dies, i.e., goes out of scope, the ambient space in
        which the precision module lives can be reduced (by projection onto the
        hyperplane defined by the dead variable).
        This reduction has a cost because it leads to re-echelonization
        of a part of the matrix that encodes the precision. The size of this
        part is roughly measured by the number of columns between the last
        column and the one corresponding to the dead variable.

        This threshold returned by this method is the maximal distance until
        which a column of a dead variable is removed and the matrix
        re-echelonized. Beyond the threshold, the column of the dead variable
        is kept in this matrix as if the variable were not destroyed.

        INPUT:

        - ``threshold`` -- nonnegative integer, ``Infinity`` or ``None``
          (default: ``None``); if not ``None`` set the threshold to the given
          value.

        .. NOTE::

            Setting the threshold to ``0`` disables the dimension reduction.

            Setting the threshold to ``Infinity`` forces the dimension reduction
            after each deletion.

        EXAMPLES::

            sage: R = ZpLC(2, label='threshold_deletion')
            sage: prec = R.precision()
            sage: prec.threshold_deletion()
            50

            sage: prec.threshold_deletion(20)
            20
            sage: prec.threshold_deletion()
            20

            sage: prec.threshold_deletion(-2)
            Traceback (most recent call last):
            ...
            ValueError: The threshold must be a nonnegative integer or Infinity
        """
        if threshold is not None:
            if threshold is Infinity or (threshold in ZZ and threshold >= 0):
                self._threshold_deletion = threshold
            else:
                raise ValueError("The threshold must be a nonnegative integer or Infinity")
        return self._threshold_deletion

    def prime(self):
        r"""
        Return the underlying prime number attached to this precision lattice.

        EXAMPLES::

            sage: R = ZpLC(2, label='mylabel')
            sage: R.precision().prime()
            2
        """
        return self._p

    def _index(self, ref):
        r"""
        Return the index of the column in the precision matrix that
        corresponds to ``ref``.

        Only for internal use.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLC(2, label='index')
            sage: prec = R.precision()
            sage: x = R(1, 10)
            sage: y = R(1, 5)

            sage: prec._index(pAdicLatticeElementWeakProxy(x))
            0
            sage: prec._index(pAdicLatticeElementWeakProxy(y))
            1

            sage: del x
            sage: prec.del_elements()
            sage: prec._index(pAdicLatticeElementWeakProxy(y))
            0
        """
        return self._elements.index(ref)

    def ambient_dimension(self):
        r"""
        Return the dimension of the vector space in which the precision
        module/lattice lives.

        EXAMPLES::

            sage: R = ZpLC(2, label='ambient_dim')
            sage: prec = R.precision()

            sage: x, y = R(1, 10), R(1, 5)
            sage: prec.ambient_dimension()
            2
            sage: prec.dimension()
            2

            sage: u = x + y
            sage: prec.ambient_dimension()
            3
            sage: prec.dimension()
            3

        In the case of ``ZpLC`` (lattice-cap precision), it is always
        equal to the dimension of the lattice.

        In the case of ``ZpLF`` (lattice-float precision), the precision
        object is not necessarily a lattice and then may have smaller
        dimension::

            sage: R = ZpLF(2, label='ambient_dim')
            sage: prec = R.precision()

            sage: x, y = R(1, 10), R(1, 5)
            sage: prec.ambient_dimension()
            2
            sage: prec.dimension()
            2

            sage: u = x + y
            sage: prec.ambient_dimension()
            3
            sage: prec.dimension()
            2
        """
        return len(self._matrix)

    @abstract_method
    def dimension(self):
        r"""
        Return the dimension of this precision module.

        EXAMPLES::

            sage: R = ZpLC(5, label='dim')
            sage: prec = R.precision()
            sage: prec.dimension()
            0

            sage: x = R(1, 10)
            sage: prec.dimension()
            1
        """
        pass

    @abstract_method
    def _new_element(self, *args, **kwargs):
        r"""
        Insert a new element in this precision module.

        TESTS::

            sage: R = ZpLC(2)
            sage: x = R.random_element()
            sage: y = R.random_element()
            sage: z = x*y  # indirect doctest
        """
        pass

    def _record_collected_element(self, ref):
        r"""
        Record that the element with weak reference ``ref``
        has been collected by the garbage collector.

        INPUT:

        - ``ref`` -- the weak reference of the collected element

        TESTS::

            sage: R = ZpLC(2, label='gc')
            sage: prec = R.precision()

            sage: x = R.random_element()
            sage: prec._collected_references
            []

            sage: del x
            sage: prec._collected_references
            [WeakProxy#...]
        """
        self._collected_references.append(ref)

    @abstract_method
    def del_elements(self, threshold=None):
        r"""
        Delete (or mark for future deletion) the columns of precision
        matrix corresponding to elements that were collected by the
        garbage collector.

        INPUT:

        - ``threshold`` -- integer or ``None`` (default: ``None``);
          a column whose distance to the right is greater than the
          threshold is not erased but marked for deletion.
          If ``None``, always erase (never mark for deletion).

        EXAMPLES::

            sage: R = ZpLC(2, label='del_elements')
            sage: prec = R.precision()

            sage: x = R(1, 10)
            sage: prec
            Precision lattice on 1 object (label: del_elements)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024]

            sage: del x
            sage: prec
            Precision lattice on 1 object (label: del_elements)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024]

            sage: prec.del_elements()
            sage: prec
            Precision lattice on 0 objects (label: del_elements)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            []
        """
        pass

    @abstract_method
    def _precision_absolute(self, x):
        r"""
        Return the absolute precision of the given element.

        INPUT:

        - ``x`` -- the element whose absolute precision is requested

        .. NOTE::

            The absolute precision is obtained by projecting the precision
            lattice onto the line of coordinate ``dx``.

            This function is not meant to be called directly. Call
            ``x.precision_absolute()`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(1, 5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)
            sage: z.precision_absolute()  # indirect doctest
            5
        """
        pass

    @abstract_method
    def precision_lattice(self, elements=None):
        r"""
        Return a lattice modeling the precision on the given set of elements
        or, if not given, on the whole set of elements tracked by the precision
        module.

        INPUT:

        - ``elements`` -- list of elements or ``None`` (default: ``None``)

        EXAMPLES::

            sage: R = ZpLC(2, label='precision_lattice')
            sage: prec = R.precision()
            sage: x = R(1, 10); y = R(1, 5)
            sage: u = x + y
            sage: v = x - y
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [         1024             0          1024          1024]
            [            0            32            32 1099511627744]
            [            0             0       2097152             0]
            [            0             0             0 1099511627776]
            sage: prec.precision_lattice([u, v])                                        # needs sage.geometry.polyhedron
            [  32 2016]
            [   0 2048]

        If the precision module does not project to a lattice,
        an error is raised. ::

            sage: R = ZpLF(2, label='precision_lattice')
            sage: prec = R.precision()
            sage: x = R(1, 10); y = R(1, 5)
            sage: u = x + y
            sage: v = x - y
            sage: prec.precision_lattice([x,y,u,v])                                     # needs sage.geometry.polyhedron
            Traceback (most recent call last):
            ...
            PrecisionError: the differential is not surjective
        """
        pass

    def diffused_digits(self, elements=None):
        r"""
        Return the number of diffused digits of precision within a
        subset of elements.

        A diffused digit of precision is a known digit which is not
        located on a single variable but only appears on a suitable
        linear combination of variables.

        The number of diffused digits of precision quantifies the
        quality of the approximation of the lattice precision by a
        jagged precision (that is a precision which is split over
        all variables).

        We refer to [CRV2018]_ for a detail exposition of the notion of
        diffused digits.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: prec = R.precision()
            sage: x = R(1, 10); y = R(1, 5)
            sage: u = x + y
            sage: v = x - y

            sage: prec.diffused_digits([x, y])                                          # needs sage.geometry.polyhedron
            0
            sage: prec.diffused_digits([u, v])                                          # needs sage.geometry.polyhedron
            6

        The elements `u` and `v` are known at absolute precision `O(2^5)`.
        However, the sum `u + v = 2x` is known at precision `O(2^11)`, that
        is with `6` more digits.
        That is where the `6` diffused digits of precision comes from.

        Here is another example with matrices::

            sage: M = matrix(R, 2, 2, [R(3, 5), R(7, 5), R(1, 5), R(11, 1)])            # needs sage.modules
            sage: N = M^10                                                              # needs sage.modules

        The next syntax provides as easy way to select an interesting
        subset of variables (the selected subset consists of the four
        entries of the matrix ``N``)::

            sage: prec.diffused_digits(N)                                               # needs sage.geometry.polyhedron sage.modules
            17

        Note that, in some cases, the number of diffused digits can be
        infinite::

            sage: R = ZpLF(2)
            sage: prec = R.precision()
            sage: x = R(1, 10)
            sage: y = x
            sage: prec.diffused_digits([x, y])                                          # needs sage.geometry.polyhedron
            +Infinity
        """
        try:
            M = self.precision_lattice(elements)
        except PrecisionError:
            return Infinity
        n = M.nrows()
        p = self._p
        return sum(M[i, i].valuation(p) - min(M[j, i].valuation(p) for j in range(i + 1)) for i in range(n))

    def tracked_elements(self, values=True, dead=True):
        r"""
        Return the list of tracked elements.

        INPUT:

        - ``values`` -- boolean (default: ``True``); if ``False``,
          the method returns a list of weak references on tracked
          elements instead

        - ``dead`` -- boolean (default: ``True``); whether dead
          elements for which the corresponding column is still not
          erased should be listed or not

        EXAMPLES::

            sage: R = ZpLC(2, label='tracked')
            sage: prec = R.precision()
            sage: x = R(1, 10); y = R(1, 5)
            sage: prec.tracked_elements()
            [1 + O(2^10), 1 + O(2^5)]
            sage: prec.tracked_elements(values=False)
            [WeakProxy#...,
             WeakProxy#...,
             WeakProxy#...]
            sage: prec.tracked_elements(values=False, dead=False)
            [WeakProxy#...,
             WeakProxy#...]

            sage: # needs sage.rings.padics
            sage: u = x + y
            sage: v = x - y
            sage: prec.tracked_elements()
            [1 + O(2^10), 1 + O(2^5), 2 + O(2^5), O(2^5)]
            sage: prec.tracked_elements(values=False)
            [WeakProxy#...,
             WeakProxy#...,
             WeakProxy#...,
             WeakProxy#...,
             WeakProxy#...]
            sage: del x; del y
            sage: prec.tracked_elements()
            [None, None, 2 + O(2^5), O(2^5), None]
            sage: prec.tracked_elements(values=False)
            [WeakProxy#...,
             WeakProxy#...,
             WeakProxy#...]
        """
        ret = [ ref for ref in self._elements if dead or ref() is not None]
        if values:
            ret = [ ref() for ref in ret ]
        return ret

    # History

    def history_enable(self):
        r"""
        Enable history.

        We refer to the documentation of the method :meth:`history` for
        a complete documentation (including examples) about history.

        TESTS::

            sage: R = ZpLC(2, label='history_en')
            sage: prec = R.precision()

            sage: print(prec.history())  # history is disabled by default
            Traceback (most recent call last):
            ...
            ValueError: History is not tracked

            sage: prec.history_enable()
            sage: print(prec.history())
             Timings
               ---

        .. SEEALSO::

            :meth:`history`, :meth:`history_disable`, :meth:`history_clear`
        """
        if self._history is None:
            self._history_init = ( len(self._elements), list(self._marked_for_deletion) )
            self._history = [ ]

    def history_disable(self):
        r"""
        Disable history.

        We refer to the documentation of the method :meth:`history` for
        a complete documentation (including examples) about history.

        TESTS::

            sage: R = ZpLC(2, label='history_dis')
            sage: prec = R.precision()

            sage: print(prec.history())  # history is disabled by default
            Traceback (most recent call last):
            ...
            ValueError: History is not tracked

            sage: prec.history_enable()
            sage: print(prec.history())
             Timings
               ---

            sage: prec.history_disable()
            sage: print(prec.history())
            Traceback (most recent call last):
            ...
            ValueError: History is not tracked

        .. SEEALSO::

            :meth:`history`, :meth:`history_enable`, :meth:`history_clear`
        """
        self._history = self._history_init = None

    def history_clear(self):
        r"""
        Clear history.

        We refer to the documentation of the method :meth:`history` for
        a complete documentation (including examples) about history.

        TESTS::

            sage: R = ZpLC(2, label='history_clear')
            sage: prec = R.precision()
            sage: prec.history_enable()

            sage: x = R(1, 10); y = R(1, 5)
            sage: x, y = x+y, x-y
            sage: print(prec.history())  # somewhat random
             Timings
            0.000213s  oooo

        When we clear history, only the last line is kept::

            sage: prec.history_clear()
            sage: print(prec.history())
             Timings   oooo
               ---     oooo

            sage: prec.del_elements()

            sage: print(prec.history())  # somewhat random
             Timings   oooo
            0.000005s  ~~oo
            0.000285s  oo

        .. SEEALSO::

            :meth:`history`, :meth:`history_enable`, :meth:`history_disable`
        """
        if self._history is None:
            raise ValueError("History is not tracked")
        self._history_init = ( len(self._elements), list(self._marked_for_deletion) )
        self._history = [ ]

    def _format_history(self, time, status, timings):
        r"""
        Return a formatted output for the history.

        This is a helper function for the method :meth:`history`.

        TESTS::

            sage: R = ZpLC(2, label='history_en')
            sage: prec = R.precision()
            sage: prec._format_history(1.23456789, ['o', 'o', 'o', 'o', 'o', 'o', '~', 'o', 'o'], true)
            '1.234568s  oooooo~oo'
            sage: prec._format_history(1.23456789, ['o', 'o', 'o', 'o', 'o', 'o', '~', 'o', 'o'], false)
            'oooooo~oo'

            sage: prec._format_history(12.3456789, ['o', 'o', 'o', 'o', 'o', 'o', '~', 'o', 'o'], true)
            '  >= 10s   oooooo~oo'
            sage: prec._format_history(10^(-10), ['o', 'o', 'o', 'o', 'o', 'o', '~', 'o', 'o'], true)
            '   ---     oooooo~oo'
            sage: prec._format_history(-1, ['o', 'o', 'o', 'o', 'o', 'o', '~', 'o', 'o'], true)
            ' Timings   oooooo~oo'
        """
        status = ''.join(status)
        if timings:
            if time < 0:
                s = " Timings "
            elif time < 0.000001:
                s = "   ---   "
            elif time >= 10:
                s = "  >= 10s "
            else:
                s = "%.6fs" % time
            return s + "  " + status
        else:
            return status

    def history(self, compact=True, separate_reduce=False, timings=True, output_type='asciiart'):
        r"""
        Show history.

        The history records creations and deletions of elements attached
        to this precision lattice, together with many timings.

        INPUT:

        - ``compact`` -- boolean (default: ``True``); if ``True``, all
          consecutive operations of the same type appear on a single row

        - ``separate_reduce`` -- boolean (default: ``False``); specify
          whether partial/full Hermite reduction should be displayed
          separately

        - ``timings`` -- boolean (default: ``True``); specify whether
          timings should be displayed

        - ``output_type`` -- only ``asciiart`` is implemented for now

        IMPORTANT NOTE:

        History is disabled by default.
        It should then be enabled (through a call to the method :meth:`history_enable`)
        before use.

        EXAMPLES::

            sage: R = ZpLC(2, label='history_en')
            sage: prec = R.precision()

        We first enable history::

            sage: prec.history_enable()

        At the beginning, the history is of course empty::

            sage: print(prec.history())
             Timings
               ---

        Now we start creating and deleting elements::

            sage: L = [ R.random_element() for _ in range(20) ]
            sage: for p in range(20):
            ....:    if is_prime(p): L[p] = None
            sage: prec.del_elements()

            sage: print(prec.history())  # somewhat random
             Timings
            0.001108s  oooooooooooooooooooo
            0.000009s  oo~~o~o~ooo~o~ooo~o~
            0.014250s  oooooooooooo

        The legend is the following:

        - the symbol ``o`` represents a tracked element,
        - the symbol ``~`` represents an element which is marked for deletion.

        On the history, we see:

        - 1st line: twenty new elements were created
          (this corresponds to the affectation of the list ``L``);
        - 2nd line: elements at prime positions were marked for deletion
          (this corresponds to the ``for`` loop);
        - 3rd line: the above elements are indeed deleted
          (this corresponds to the call of the method :meth:`del_elements`.

        Here are some variants::

            sage: print(prec.history(timings=False))
            oooooooooooooooooooo
            oo~~o~o~ooo~o~ooo~o~
            oooooooooooo

            sage: print(prec.history(separate_reduce=True))  # somewhat random
             Timings
            0.001063s  oooooooooooooooooooo
            0.000014s  oo~~o~o~ooo~o~ooo~o~
            0.000798s  oo~~o~o~ooo~ooooo
            0.000233s  oo~~o~o~ooo~orrrr
            0.000824s  oo~~o~o~oooooooo
            0.000375s  oo~~o~o~ooorrrrr
            0.001724s  oo~~o~ooooooooo
            0.001020s  oo~~o~orrrrrrrr
            0.001989s  oo~~oooooooooo
            0.001303s  oo~~orrrrrrrrr
            0.002352s  oo~oooooooooo
            0.001632s  oo~rrrrrrrrrr
            0.002265s  oooooooooooo
            0.001630s  oorrrrrrrrrr
               ---     oooooooooooo

        The symbol ``r`` represents a column of the precision matrix which is
        currently under partial Hermite reduction.

        Timings for automatic reduction do not appear because they are included
        in the timings for deletion.

        The symbol ``R`` is used to symbolize a column which is under full
        Hermite reduction. Note that full Hermite reduction are never performed
        automatically but needs to be called by hand::

            sage: prec.reduce()
            sage: print(prec.history(separate_reduce=True))  # somewhat random
             Timings
            0.001063s  oooooooooooooooooooo
            0.000014s  oo~~o~o~ooo~o~ooo~o~
            0.000798s  oo~~o~o~ooo~ooooo
            0.000233s  oo~~o~o~ooo~orrrr
            0.000824s  oo~~o~o~oooooooo
            0.000375s  oo~~o~o~ooorrrrr
            0.001724s  oo~~o~ooooooooo
            0.001020s  oo~~o~orrrrrrrr
            0.001989s  oo~~oooooooooo
            0.001303s  oo~~orrrrrrrrr
            0.002352s  oo~oooooooooo
            0.001632s  oo~rrrrrrrrrr
            0.002265s  oooooooooooo
            0.001630s  oorrrrrrrrrr
            0.001486s  RRRRRRRRRRRR
               ---     oooooooooooo

        Here is a more common example with matrices::

            sage: R = ZpLC(3)
            sage: prec = R.precision()
            sage: prec.history_enable()
            sage: M = random_matrix(R, 5)                                               # needs sage.geometry.polyhedron
            sage: d = M.determinant()                                                   # needs sage.geometry.polyhedron
            sage: print(prec.history())  # somewhat random
               ---
            0.004212s  oooooooooooooooooooooooooooooooooooo
            0.000003s  oooooooooooooooooooooooooooooooooo~~
            0.000010s  oooooooooooooooooooooooooooooooooo
            0.001560s  ooooooooooooooooooooooooooooooooooooooooo
            0.000004s  ooooooooooooooooooooooooooooo~oooo~oooo~o
            0.002168s  oooooooooooooooooooooooooooooooooooooo
            0.001787s  ooooooooooooooooooooooooooooooooooooooooo
            0.000004s  oooooooooooooooooooooooooooooooooooooo~~o
            0.000198s  ooooooooooooooooooooooooooooooooooooooo
            0.001152s  ooooooooooooooooooooooooooooooooooooooooo
            0.000005s  ooooooooooooooooooooooooooooooooo~oooo~~o
            0.000853s  oooooooooooooooooooooooooooooooooooooo
            0.000610s  ooooooooooooooooooooooooooooooooooooooo
             [...]
            0.003879s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.000006s  oooooooooooooooooooooooooooooooooooooooooooooooooooo~~~~~
            0.000036s  oooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.006737s  oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.000005s  oooooooooooooooooooooooooooooooooooooooooooooooooooo~~~~~ooooo
            0.002637s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.007118s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.000008s  oooooooooooooooooooooooooooooooooooooooooooooooooooo~~~~o~~~~oooo
            0.003504s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.005371s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.000006s  ooooooooooooooooooooooooooooooooooooooooooooooooooooo~~~o~~~ooo
            0.001858s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.003584s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.000004s  oooooooooooooooooooooooooooooooooooooooooooooooooooooo~~o~~oo
            0.000801s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.001916s  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            0.000022s  ooooooooooooooooooooooooooooo~~~~~~~~~~~~~~~~~~~~~~oooo~o~o
            0.014705s  ooooooooooooooooooooooooooooooooooo
            0.001292s  ooooooooooooooooooooooooooooooooooooo

        We observe that deleted variables appear mostly on the right.
        This is the so-called principal of temporal locality.

        .. SEEALSO::

            :meth:`history_enable`, :meth:`history_disable`, :meth:`history_clear`
        """
        if self._history is None:
            raise ValueError("History is not tracked")
        total_time = 0
        if output_type == 'asciiart':
            # Legend:
            #  o : tracked element
            #  ~ : element marked for deletion
            #  r : partial reduction
            #  R : full Hermite reduction
            (n, mark) = self._history_init
            status = n*['o']
            for index in mark:
                status[index] = '~'
            hist = [ self._format_history(-1, status, timings) ]
            oldevent = ''
            total_time = 0
            for (event, index, tme) in self._history:
                if event == 'partial reduce' or event == 'full reduce':
                    if separate_reduce:
                        if status:
                            hist.append(self._format_history(total_time, status, timings))
                        if event == 'partial reduce':
                            code = 'r'
                        else:
                            code = 'R'
                        status_red = status[:index] + (len(status) - index) * [code]
                        hist.append(self._format_history(tme, status_red, timings))
                        total_time = 0
                        oldevent = ''
                    else:
                        total_time += tme
                    continue
                if not compact or event != oldevent:
                    if status:
                        hist.append(self._format_history(total_time, status, timings))
                    total_time = 0
                    oldevent = event
                total_time += tme
                if event == 'add':
                    if index is None:
                        status.append('o')
                    else:
                        status = status[:index] + ['o'] + status[index:]
                elif event == 'mark':
                    status[index] = '~'
                elif event == 'del':
                    del status[index]
            if status or oldevent == '':
                hist.append(self._format_history(total_time, status, timings))
            return '\n'.join(hist)
        else:
            raise NotImplementedError

    def timings(self, action=None):
        r"""
        Return cumulated timings (grouped by actions) since the last
        time history has been cleared.

        INPUT:

        - ``action`` -- ``None`` (default), ``'add'``, ``'mark'``, ``'del'``,
          ``'partial reduce'`` or ``'full reduce'``; if not ``None``, return the
          cumulated timing corresponding to this action; otherwise, return
          a dictionary

        Here are the meanings of the keywords above:

        - ``'add'``: time spent in adding new columns to the precision matrix
          (corresponding to the creation of new elements)
        - ``'mark'``: time spent in marking elements for deletion
        - ``'del'``: time spent in deleting columns of the precision matrix
          and re-echelonizing the matrix
        - ``'partial reduce'``: time spent in partial Hermite reduction
        - ``'full reduce'``: time spent in full Hermite reduction.

        EXAMPLES::

            sage: R = ZpLC(2, label='timings')
            sage: prec = R.precision()
            sage: prec.history_enable()
            sage: M = random_matrix(R, 5, 5)                                            # needs sage.geometry.polyhedron
            sage: N = M^10                                                              # needs sage.geometry.polyhedron
            sage: prec.timings()    # somewhat random
            {'add': 1.0530245304107666,
             'del': 0.24358701705932617,
             'mark': 0.0013289451599121094,
             'partial reduce': 0.21604204177856445
             'full reduce': 0}

        TESTS::

            sage: prec.history_clear()
            sage: prec.timings()
            {'add': 0, 'del': 0, 'full reduce': 0, 'mark': 0, 'partial reduce': 0}
        """
        if self._history is None:
            raise ValueError("History is not tracked")
        tme_by_event = { 'add': 0, 'del': 0, 'mark': 0, 'partial reduce': 0, 'full reduce': 0 }
        for (event, _, tme) in self._history:
            tme_by_event[event] += tme
        if action is None:
            return tme_by_event
        if action in tme_by_event:
            return tme_by_event[action]
        else:
            raise ValueError("invalid event")


class PrecisionLattice(UniqueRepresentation, DifferentialPrecisionGeneric):
    r"""
    A class for handling precision lattices which are used to
    track precision in the ZpLC model.

    The precision lattice is stored as a triangular matrix whose
    rows are generators of the lattice.

    INPUT:

    - ``p`` -- a prime number

    - ``label`` -- string; the label of the parents to which the elements
      tracked by this lattice belong

    .. NOTE::

        You should not create instances of this class directly. The precision
        lattice is automatically initialized at the creation of the parent.

    EXAMPLES::

        sage: R = ZpLC(2, label='init')
        sage: R.precision()
        Precision lattice on 0 objects (label: init)
    """
    def __init__(self, p, label):
        r"""
        TESTS::

            sage: from sage.rings.padics.lattice_precision import PrecisionLattice
            sage: R = ZpLC(2)
            sage: isinstance(R.precision(), PrecisionLattice)
            True
        """
        DifferentialPrecisionGeneric.__init__(self, p, label)
        self._repr_type = "Precision lattice"
        self._capped = { }

    # We need to copy this method.
    # Indeed otherwise it is inherited from UniqueRepresentation
    def __reduce__(self):
        r"""
        TESTS::

            sage: R = ZpLC(2)
            sage: prec = R.precision()
            sage: dumps(prec)
            Traceback (most recent call last):
            ...
            NotImplementedError: pickling/unpickling precision modules is not implemented yet
        """
        raise NotImplementedError("pickling/unpickling precision modules is not implemented yet")

    def _index(self, ref):
        r"""
        Return the index of the element whose reference is ``ref``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLC(2, label='index')
            sage: prec = R.precision()
            sage: x = R(1, 10)
            sage: y = R(1, 5)

            sage: prec._index(pAdicLatticeElementWeakProxy(x))
            0
            sage: prec._index(pAdicLatticeElementWeakProxy(y))
            1

            sage: del x
            sage: prec.del_elements()
            sage: prec._index(pAdicLatticeElementWeakProxy(y))
            0
        """
        return len(self._matrix[ref]) - 1

    def dimension(self):
        r"""
        Return the dimension of this lattice.

        EXAMPLES::

            sage: R = ZpLC(5, label='dimension')
            sage: prec = R.precision()
            sage: prec.dimension()
            0

            sage: x = R(1, 10)
            sage: prec.dimension()
            1
        """
        return len(self._matrix)

    def reduce(self, index=0, partial=False):
        r"""
        Reduce the size of the entries above the diagonal of the precision matrix.

        INPUT:

        - ``index`` -- integer; the starting row for which the reduction
          is performed

        - ``partial`` -- boolean (default: ``False``); specifying whether a
          partial or a full Hermite reduction should be performed

        NOTE:

        The partial reduction has cost `O(m^2)` where `m` is the number of
        rows that need to be reduced (that is the difference between the
        total number of rows and ``index``).

        The full Hermite reduction has cost `O(m^3)`.

        .. NOTE::

            The software ensures that the precision lattice is always partially
            reduced.  Calling the function manually with the argument
            ``partial=True`` should then just do nothing.

        TESTS::

            sage: R = ZpLC(2)
            sage: x = R.random_element()
            sage: del x
            sage: R.precision().del_elements()   # indirect doctest
        """
        n = len(self._elements)
        if index >= n-1:
            return
        if partial:
            # Partial reduction
            # Cost: O(m^2) with m = n-index
            tme = walltime()
            diffval = (n-index) * [0]
            for j in range(n-1, index, -1):
                col = self._matrix[self._elements[j]]
                prec = col[j].valuation() - diffval[j-index]
                for i in range(index, j):
                    col[i] = col[i].reduce(prec)
                    col[i].normalize()
                    dval = col[i].valuation() - prec
                    diffval[i-index] = min(dval, diffval[i-index])
            # We update history
            if self._history is not None:
                self._history.append(('partial reduce', index, walltime(tme)))
        else:
            # Full Hermite reduction
            # Cost: O(m^3) with m = n-index
            tme = walltime()
            for j in range(index+1, n):
                # In what follows, we assume that col[j] is a power of p
                col = self._matrix[self._elements[j]]
                valpivot = col[j].valuation()
                for i in range(index, j):
                    reduced = col[i].reduce(valpivot)
                    scalar = (col[i] - reduced) >> valpivot
                    if scalar.is_zero():
                        continue
                    col[i] = reduced
                    col[i].normalize()
                    for j2 in range(j+1, n):
                        col2 = self._matrix[self._elements[j2]]
                        col2[i] -= scalar*col2[i]
                        col2[i].normalize()
            # We update history
            if self._history is not None:
                self._history.append(('full reduce', index, walltime(tme)))

    def _new_element(self, x, dx, bigoh, dx_mode='linear_combination', capped=False):
        r"""
        Update the lattice when a new element is created.

        This function is not meant to be called manually.
        It is automatically called by the parent when a new
        element is created.

        INPUT:

        - ``x`` -- the newly created element

        - ``dx`` -- dictionary representing the differential of ``x``

        - ``bigoh`` -- integer or ``None`` (default: ``None``); the
          bigoh to be added to the precision of ``x``; if ``None``, the
          default cap is used.

        - ``dx_mode`` -- string; either ``linear_combination`` (the default)
          or ``values``

        - ``capped`` -- boolean; whether this element has been capped
          according to the parent's cap

        If ``dx_mode`` is ``linear_combination``, the dictionary ``dx``
        encodes the expression of the differential of ``x``.
        For example, if ``x`` was defined as ``x = y*z`` then:

        .. MATH::

            dx = y dz + z dy

        and the corresponding dictionary is ``{z: y, y: z}`` (except
        that the keys are not the elements themselves but weak references
        to them).

        If ``dx_mode`` is ``values``, the dictionary ``dx`` directly
        specifies the entries that have to be stored in the precision lattice.
        This mode is only used for multiple conversion between different
        parents (see :meth:`multiple_conversion`).

        TESTS::

            sage: R = ZpLC(2)
            sage: x = R.random_element()
            sage: y = R.random_element()
            sage: z = x*y    # indirect doctest
        """
        # First we delete some elements marked for deletion
        self.del_elements(threshold=self._threshold_deletion)

        # Then we add the new element
        tme = walltime()
        p = self._p
        n = len(self._elements)
        x_ref = pAdicLatticeElementWeakProxy(x, self._record_collected_element)
        self._elements.append(x_ref)
        col = n * [self._approx_zero]
        if dx_mode == 'linear_combination':
            for elt, scalar in dx:
                ref = pAdicLatticeElementWeakProxy(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                c = self._matrix[ref]
                for i in range(len(c)):
                    col[i] += scalar * c[i]
        elif dx_mode == 'values':
            for elt, scalar in dx:
                ref = pAdicLatticeElementWeakProxy(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                i = self._index(ref)
                col[i] = scalar
        else:
            raise ValueError("dx_mode must be either 'linear_combination' or 'values'")
        for i in range(n):
            col[i] = col[i].reduce(bigoh)
        col.append(pRational(p, ZZ(1), bigoh))
        self._matrix[x_ref] = col
        self._capped[x_ref] = capped

        # We update history
        if self._history is not None:
            self._history.append(('add', None, walltime(tme)))

    def del_elements(self, threshold=None):
        r"""
        Erase columns of the lattice precision matrix corresponding to
        elements which are marked for deletion and echelonize the matrix
        in order to keep it upper triangular.

        INPUT:

        - ``threshold`` -- integer or ``None`` (default: ``None``);
          a column whose distance to the right is greater than the
          threshold is not erased

        EXAMPLES::

            sage: R = ZpLC(2, label='delelts')
            sage: prec = R.precision()

            sage: x = R(1, 10)
            sage: prec
            Precision lattice on 1 object (label: delelts)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024]

            sage: del x
            sage: prec
            Precision lattice on 1 object (label: delelts)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024]

            sage: prec.del_elements()
            sage: prec
            Precision lattice on 0 objects (label: delelts)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            []
        """
        n = len(self._elements)

        # We mark new collected elements for deletion
        # The list self._collected_references can be updated while
        # the loop runs.
        # However, we do not need to copy it because Python supports
        # iteration over a list to which elements are added.
        count = 0
        for ref in self._collected_references:
            count += 1
            tme = walltime()
            index = self._index(ref)
            self._marked_for_deletion.append(index)
            if self._history is not None:
                self._history.append(('mark', index, walltime(tme)))
        del self._collected_references[:count]

        # We erase corresponding columns and echelonize
        self._marked_for_deletion.sort(reverse=True)
        count = 0
        for index in self._marked_for_deletion:
            if threshold is not None and index < n - threshold:
                break
            n -= 1
            count += 1

            tme = walltime()
            ref = self._elements[index]
            del self._elements[index]
            del self._matrix[ref]
            capped = self._capped[ref]
            del self._capped[ref]

            # Now, we echelonize
            for i in range(index, n):
                ref = self._elements[i]
                col = self._matrix[ref]
                if col[i].valuation() < col[i+1].valuation():
                    self._capped[ref], capped = capped, capped or self._capped[ref]
                else:
                    capped = capped or self._capped[ref]

                d, u, v = col[i].xgcd(col[i+1])
                up, vp = col[i+1]/d, col[i]/d
                col[i] = d
                del col[i+1]
                for j in range(i+1, n):
                    col = self._matrix[self._elements[j]]
                    col[i], col[i+1] = u*col[i] + v*col[i+1], up*col[i] - vp*col[i+1]

            # We update history
            if self._history is not None:
                self._history.append(('del', index, walltime(tme)))

            # And we reduce a bit
            # (we do not perform a complete reduction because it is costly)
            self.reduce(index, partial=True)

        del self._marked_for_deletion[:count]

    def _lift_to_precision(self, x, prec):
        r"""
        Lift the specified element to the specified precision.

        INPUT:

        - ``x`` -- the element whose precision has to be lifted

        - ``prec`` -- the new precision

        NOTE:

        The new precision lattice is computed as the intersection
        of the current precision lattice with the subspace

        .. MATH::

            p^{prec} \Z_p dx \oplus \bigoplus_{y \neq x} \Q_p dy

        This function may change at the same time the precision of
        other elements having the same parent.

        .. NOTE::

            This function is not meant to be called directly. Use
            ``x.lift_to_precision`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(1, 5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)

            sage: prec = R.precision()
            sage: prec._lift_to_precision(z, 12)
            sage: z
            2 + O(2^12)
            sage: y
            1 + O(2^10)
        """
        ref = pAdicLatticeElementWeakProxy(x)
        col = self._matrix[ref]
        n = len(self._elements)

        rows_by_val = defaultdict(list)
        for i in range(len(col)):
            v = col[i].valuation()
            if v >= prec:
                continue
            rows_by_val[v].append(i)
        vals = sorted(rows_by_val)
        vals.append(prec)

        for t in range(len(vals)-1):
            v, w = vals[t], vals[t+1]
            rows = rows_by_val[v]
            piv = max(rows)
            for i in rows:
                if i == piv:
                    continue
                # We clear the entry on the i-th row
                scalar = (col[i]/col[piv]).reduce(prec-v)
                for j in range(piv,n):
                    col_cur = self._matrix[self._elements[j]]
                    col_cur[i] -= scalar*col_cur[piv]
            # We rescale the piv-th row
            for j in range(piv,n):
                col_cur = self._matrix[self._elements[j]]
                col_cur[piv] <<= w - v
            # Now the entry on the piv-th row has valuation w
            # We update the dictionary accordingly
            if w < prec:
                rows_by_val[w].append(piv)

        self._precision_absolute_data.clear_cache()

    @cached_method(key=lambda self, x: pAdicLatticeElementWeakProxy(x))
    def _precision_absolute_data(self, x):
        r"""
        Return absolute precision data for ``x``.

        .. NOTE::

            Helper method for :meth:`_precision_absolute` and
            :meth:`_is_precision_capped`.

        TESTS::

            sage: R = ZpLC(2)
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(1, 5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)
            sage: z.precision_absolute()  # indirect doctest
            5
        """
        ref = pAdicLatticeElementWeakProxy(x)
        col = self._matrix[ref]
        absprec = Infinity
        capped = False
        for i in range(len(col)):
            v = col[i].valuation()
            if v < absprec:
                absprec = v
                capped = self._capped[self._elements[i]]
            elif v == absprec:
                capped = capped and self._capped[self._elements[i]]
        return (absprec, capped)

    def _precision_absolute(self, x):
        r"""
        Return the absolute precision of the given element.

        INPUT:

        - ``x`` -- an element in the parent corresponding to this lattice

        .. NOTE::

            The absolute precision is obtained by projecting the precision
            lattice onto the line of coordinate ``dx``.

            This function is not meant to be called directly. Call
            ``x.precision_absolute()`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(1, 5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)
            sage: z.precision_absolute()  # indirect doctest
            5
        """
        return self._precision_absolute_data(x)[0]

    def _is_precision_capped(self, x):
        r"""
        Return whether the absolute precision of ``x`` results from a cap
        coming from the parent.

        INPUT:

        - ``x`` -- an element in the parent corresponding to this lattice

        .. NOTE::

            This function is not meant to be called directly. Call
            ``x.is_precision_capped`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: x.is_precision_capped()  # indirect doctest
            False

            sage: y = x-x; y
            O(2^40)
            sage: y.is_precision_capped()  # indirect doctest
            True
        """
        return self._precision_absolute_data(x)[1]

    def precision_lattice(self, elements=None):
        r"""
        Return a matrix representing the precision lattice on a
        subset of elements.

        INPUT:

        - ``elements`` -- list of elements or ``None`` (default: ``None``)

        - ``echelon`` -- boolean (default: ``True``); whether the result
          should be in echelon form

        EXAMPLES::

            sage: R = ZpLC(2, label='preclattice')
            sage: prec = R.precision()
            sage: x = R(1, 10); y = R(1, 5)
            sage: u = x + y
            sage: v = x - y
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [         1024             0          1024          1024]
            [            0            32            32 1099511627744]
            [            0             0       2097152             0]
            [            0             0             0 1099511627776]
            sage: prec.precision_lattice([u, v])                                        # needs sage.geometry.polyhedron
            [  32 2016]
            [   0 2048]

        Here is another example with matrices::

            sage: M = matrix(R, 2, 2, [R(3, 5), R(7, 5), R(1, 5), R(11, 1)])            # needs sage.modules
            sage: N = M^10                                                              # needs sage.modules
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron sage.modules
            23 x 23 dense matrix over Integer Ring (use the '.str()' method to see the entries)

        The next syntax provides as easy way to select an interesting
        subset of variables (the selected subset consists of the four
        entries of the matrix ``N``)::

            sage: prec.precision_lattice(N)                                             # needs sage.modules
            [  2048    512  28160 230400]
            [     0   2048  14336 258048]
            [     0      0  65536  65536]
            [     0      0      0 262144]

        We can give a list of matrices as well::

            sage: prec.precision_lattice([M, N])                                        # needs sage.modules
            [       32         0         0         0 226115584  96788480  52174848  82804736]
            [        0        32         0         0  52174848 121765888  11829248  28516352]
            [        0         0        32         0  96788480  42762240 121765888 199614464]
            [        0         0         0         2   5175296  12475904   1782272   4045824]
            [        0         0         0         0 268435456         0         0         0]
            [        0         0         0         0         0 268435456         0         0]
            [        0         0         0         0         0         0 268435456         0]
            [        0         0         0         0         0         0         0 268435456]
        """
        if elements is None:
            elements = self._elements
        else:
            elements = list_of_padics(elements)
        n = len(self._elements)
        rows = []
        val = 0
        for ref in elements:
            col = self._matrix[ref]
            row = [ x.value() for x in col ]
            valcol = min([ x.valuation() for x in col ])
            val = min(valcol, val)
            row += (n-len(row)) * [ZZ(0)]
            rows.append(row)
        from sage.matrix.constructor import matrix
        M = matrix(rows).transpose()
        if val < 0:
            M *= self._p ** (-val)
        M = M.change_ring(ZZ)
        M.echelonize()
        n = len(elements)
        M = M.submatrix(0, 0, n, n)
        if val < 0:
            M *= self._p ** val
        return M


class PrecisionModule(UniqueRepresentation, DifferentialPrecisionGeneric):
    r"""
    A class for handling precision modules which are used to
    track precision in the ZpLF model.

    The precision module (which is not necessarily a lattice)
    is stored as a matrix whose rows are generators.
    """
    def __init__(self, p, label, prec):
        r"""
        Initialize this precision module.

        INPUT:

        - ``p`` -- a prime number

        - ``label`` -- string; the label of the parents to which belong
          the elements tracked by this precision module

        NOTE:

        The precision module is automatically initialized at the
        creation of the parent.

        TESTS::

            sage: R = ZpLF(2, label='init')
            sage: R.precision()
            Precision module on 0 objects (label: init)
        """
        DifferentialPrecisionGeneric.__init__(self, p, label)
        # elements whose valuation are not less than self._zero_cap are assumed to vanish
        self._zero_cap = prec
        self._internal_prec = prec + STARTING_ADDITIONAL_PREC
        self._count = 0
        self._threshold = 1
        self._repr_type = "Precision module"

    # We need to copy this method.
    # Indeed otherwise it is inherited from UniqueRepresentation
    def __reduce__(self):
        r"""
        TESTS::

            sage: R = ZpLF(2)
            sage: prec = R.precision()
            sage: dumps(prec)
            Traceback (most recent call last):
            ...
            NotImplementedError: pickling/unpickling precision modules is not implemented yet
        """
        raise NotImplementedError("pickling/unpickling precision modules is not implemented yet")

    def internal_prec(self):
        r"""
        Return the relative precision at which computations is handled
        internally.

        It is slightly greater than the actual precision and increases
        a bit (at a logarithmic rate) when new elements are created
        and/or computed.

        EXAMPLES::

            sage: R = ZpLF(5, prec=20, label='internal_prec')
            sage: prec = R.precision()

            sage: prec.internal_prec()
            25

            sage: L = [ R.random_element() for _ in range(50) ]
            sage: prec.internal_prec()
            28
        """
        return self._internal_prec

    def dimension(self):
        r"""
        Return the dimension of this precision module.

        EXAMPLES:

        In general, the dimension increases by 1 when a new
        element with a given precision is created::

            sage: R = ZpLF(2, label='dimension')
            sage: prec = R.precision()

            sage: prec.dimension()
            0
            sage: x = R.random_element(prec=10)
            sage: prec.dimension()
            1
            sage: y = R.random_element(prec=10)
            sage: prec.dimension()
            2

        However in general it does not increase while
        doing computations::

            sage: u = x + y
            sage: v = x^2 + 3*y + x*y + y^3
            sage: prec.dimension()
            2

        Of course, it may also decrease when a sufficient
        number of variables are collected::

            sage: del x, y, u
            sage: prec.del_elements()
            sage: prec.dimension()
            1

            sage: del v
            sage: prec.del_elements()
            sage: prec.dimension()
            0
        """
        if len(self._elements) == 0:
            return 0
        return len(self._matrix[self._elements[-1]])

    def is_lattice(self):
        r"""
        Return ``True`` if this precision module is a lattice
        (i.e. has maximal dimension).

        EXAMPLES::

            sage: R = ZpLF(2, label='is_lattice')
            sage: prec = R.precision()

            sage: x = R(1, 10)
            sage: y = R(1, 5)
            sage: prec.is_lattice()
            True

            sage: u = x + y
            sage: prec.is_lattice()
            False

        .. SEEALSO::

            :meth:`dimension`
        """
        return self.dimension() == len(self._elements)

    def _new_element(self, x, dx, bigoh, dx_mode='linear_combination'):
        r"""
        Update the lattice when a new element is created.

        This function is not meant to be called manually.
        It is automatically called by the parent when a new
        element is created.

        INPUT:

        - ``x`` -- the newly created element

        - ``dx`` -- dictionary representing the differential of ``x``

        - ``bigoh`` -- integer or ``None`` (default: ``None``); the
          bigoh to be added to the precision of ``x``. If ``None``, the
          default cap is used.

        - ``dx_mode`` -- string; either ``'linear_combination'`` (the
          default) or ``'values'``

        If ``dx_mode`` is ``'linear_combination'``, the dictionary ``dx``
        encodes the expression of the differential of ``x``.  For example, if
        ``x`` was defined as ``x = y*z`` then:

        .. MATH::

            dx = y dz + z dy

        and the corresponding dictionary is ``{z: y, y: z}`` (except
        that the keys are not the elements themselves but weak references
        to them).

        If ``dx_mode`` is ``'values'``, the dictionary ``dx`` directly
        specifies the entries that have to stored in the precision module.
        This mode is only used for multiple conversion between different
        parents (see :meth:`multiple_conversion`).

        TESTS::

            sage: R = ZpLF(2)
            sage: x = R.random_element()
            sage: y = R.random_element()
            sage: z = x*y    # indirect doctest
        """
        # First we delete some elements marked for deletion
        if self._marked_for_deletion:
            self.del_elements(threshold=self._threshold_deletion)

        # We increase the internal prec
        # The heuristic behind this is the following: when computing
        # with N digits of precision, we except that about N-log_p(c)
        # of them are correct after c elementary operations.
        self._count += 1
        if self._count > self._threshold:
            self._internal_prec += 1
            self._threshold *= self._p

        tme = walltime()
        p = self._p
        n = self.dimension()
        x_ref = pAdicLatticeElementWeakProxy(x, self._record_collected_element)
        col = n * [self._approx_zero]
        if dx_mode == 'linear_combination':
            expected_vals = n * [ Infinity ]
            for elt, scalar in dx:
                ref = pAdicLatticeElementWeakProxy(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                c = self._matrix[ref]
                for i in range(len(c)):
                    summand = scalar * c[i]
                    expected_vals[i] = min(expected_vals[i], summand.valuation())
                    col[i] += summand
            for i in range(n):
                if col[i].valuation() >= expected_vals[i] + self._zero_cap:
                    col[i] = self._approx_zero
        elif dx_mode == 'values':
            for elt, scalar in dx:
                ref = pAdicLatticeElementWeakProxy(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                i = self._index(ref)
                col[i] = scalar
        else:
            raise ValueError("dx_mode must be either 'linear_combination' or 'values'")

        for i in range(n):
            col[i] = col[i].reduce_relative(self._internal_prec)
        if bigoh is not None:
            col.append(pRational(p, ZZ(1), bigoh))

        self._elements.append(x_ref)
        self._matrix[x_ref] = col

        # We update history
        if self._history is not None:
            self._history.append(('add', None, walltime(tme)))

    def del_elements(self, threshold=None):
        r"""
        Erase columns of the lattice precision matrix corresponding to
        elements which were collected by the garbage collector.
        Then reduce the matrix in order to keep it in echelon form.

        INPUT:

        - ``threshold`` -- integer or ``None`` (default: ``None``);
          a non-pivot column whose distance to the right is greater than
          the threshold is not erased but only marked for future deletion

        EXAMPLES::

            sage: R = ZpLF(2, label='delelts')
            sage: prec = R.precision()

            sage: x = R(1, 10)
            sage: prec
            Precision module on 1 object (label: delelts)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024]

            sage: del x
            sage: prec
            Precision module on 1 object (label: delelts)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024]

            sage: prec.del_elements()
            sage: prec
            Precision module on 0 objects (label: delelts)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            []
        """
        # We mark new collected elements for deletion
        # The list self._collected_references can be updated while
        # the loop runs.
        # However, we do not need to copy it because Python supports
        # iteration over a list to which elements are added.
        count = 0
        for ref in self._collected_references:
            count += 1
            tme = walltime()
            index = self._index(ref)
            if index == 0:
                length_before = 0
            else:
                length_before = len(self._matrix[self._elements[index-1]])
            length = len(self._matrix[ref])
            if length > length_before:
                self._marked_for_deletion.append(index)
                if self._history is not None:
                    self._history.append(('mark', index, walltime(tme)))
            else:
                # if the column is not a pivot, we erase it without delay
                # (btw, is it a good idea?)
                del self._elements[index]
                self._marked_for_deletion = [i if i < index else i - 1
                                             for i in self._marked_for_deletion]
                if self._history is not None:
                    self._history.append(('del', index, walltime(tme)))
        del self._collected_references[:count]

        # We erase corresponding columns and echelonize
        n = len(self._elements)
        self._marked_for_deletion.sort(reverse=True)
        count = 0
        for index in self._marked_for_deletion:
            if threshold is not None and index < n - threshold:
                break
            n -= 1
            count += 1

            tme = walltime()

            length = len(self._matrix[self._elements[index]])
            del self._matrix[self._elements[index]]
            del self._elements[index]
            start = index
            while start < n:
                i = start
                val = Infinity
                end = n
                while i < n:
                    col = self._matrix[self._elements[i]]
                    if len(col) > length:
                        end = i
                        break
                    v = col[-1].valuation()
                    if v < val:
                        val = v
                        piv = i
                    i += 1
                if val < Infinity:
                    # another pivot has been found, we place it in front
                    self._elements[start], self._elements[piv] = self._elements[piv], self._elements[start]
                    break

                # No pivot was found. We re-echelonize
                for i in range(start, end):
                    del self._matrix[self._elements[i]][-1]
                if end == n:
                    break
                # col is the column of index "end"
                # its size is (length + 1)
                d, u, v = col[length-1].xgcd(col[length])
                up, vp = col[length]/d, col[length-1]/d
                col[length-1] = d.reduce_relative(self._internal_prec)
                del col[length]
                start = end + 1
                for j in range(start, n):
                    col = self._matrix[self._elements[j]]
                    a1 = u*col[length-1]
                    a2 = v*col[length]
                    a = a1 + a2
                    b1 = up*col[length-1]
                    b2 = vp * col[length]
                    b = b1 + b2
                    if a.valuation() > min(a1.valuation(), a2.valuation()) + self._zero_cap:
                        col[length-1] = self._approx_zero
                    else:
                        col[length-1] = a.reduce_relative(self._internal_prec)
                    if b.valuation() > min(b1.valuation(), b2.valuation()) + self._zero_cap:
                        col[length] = self._approx_zero
                    else:
                        col[length] = b.reduce_relative(self._internal_prec)
                length += 1

            # We update history
            if self._history is not None:
                self._history.append(('del', index, walltime(tme)))

        del self._marked_for_deletion[:count]

    def _lift_to_precision(self, x, prec):
        r"""
        Lift the specified element to the specified precision.

        INPUT:

        - ``x`` -- the element whose precision has to be lifted

        - ``prec`` -- the new precision

        .. NOTE::

            The new precision lattice is computed as the intersection
            of the current precision lattice with the subspace.

        .. MATH::

            p^{prec} \Z_p dx \oplus \bigoplus_{y \neq x} \Q_p dy

        This function may change at the same time the precision of
        other elements having the same parent.

        .. NOTE::

            This function is not meant to be called directly. Use
            ``x.lift_to_precision`` instead.

        EXAMPLES::

            sage: R = ZpLF(2)
            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(1, 5); y
            1 + O(2^5)
            sage: u = x^2 + x*y
            sage: v = y^2 + x*y
            sage: w = u + v

            sage: prec = R.precision()
            sage: prec._lift_to_precision(w, 11)
            sage: w
            2^2 + O(2^11)
            sage: y
            1 + O(2^9)
        """
        ref = pAdicLatticeElementWeakProxy(x)
        col = self._matrix[ref]
        n = len(self._elements)

        rows_by_val = defaultdict(list)
        for i in range(len(col)):
            v = col[i].valuation()
            if v >= prec:
                continue
            rows_by_val[v].append(i)
        vals = sorted(rows_by_val)
        vals.append(prec)

        for t in range(len(vals)-1):
            v, w = vals[t], vals[t+1]
            rows = rows_by_val[v]
            piv = max(rows)
            for i in rows:
                if i == piv:
                    continue
                # We clear the entry on the i-th row
                scalar = (col[i]/col[piv]).reduce(prec-v)
                for j in range(n):
                    col_cur = self._matrix[self._elements[j]]
                    if len(col_cur) > piv:
                        col_cur[i] -= scalar*col_cur[piv]
                        col_cur[i] = col_cur[i].reduce_relative(self._internal_prec)
            # We rescale the piv-th row
            # (if w is Infinity, we delete it)
            for j in range(n):
                col_cur = self._matrix[self._elements[j]]
                if len(col_cur) > piv:
                    if w is Infinity:
                        del col_cur[piv]
                    else:
                        col_cur[piv] <<= w - v
            # Now the entry on the piv-th row has valuation w
            # We update the dictionary accordingly
            if w < prec:
                rows_by_val[w].append(piv)

        self._precision_absolute.clear_cache()

    @cached_method(key=lambda self, x: pAdicLatticeElementWeakProxy(x))
    def _precision_absolute(self, x):
        r"""
        Return the absolute precision of the given element.

        INPUT:

        - ``x`` -- the element whose absolute precision is requested

        .. NOTE::

            The absolute precision is obtained by projecting the precision
            module onto the line of coordinate ``dx``.

            This function is not meant to be called directly. Call
            ``x.precision_absolute()`` instead.

        EXAMPLES::

            sage: R = ZpLF(2)
            sage: prec = R.precision()

            sage: x = R(1, 10); x
            1 + O(2^10)
            sage: y = R(1, 5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)
            sage: z.precision_absolute()  # indirect doctest
            5

        In some cases, the absolute precision returned by this function
        may be infinite::

            sage: y = R(1)
            sage: prec._precision_absolute(y)
            +Infinity

        However calling the method :meth:`absolute_precision` of the
        element itself reintroduces a cap::

            sage: y.precision_absolute()
            20
        """
        ref = pAdicLatticeElementWeakProxy(x)
        col = self._matrix[ref]
        if len(col) == 0:
            return Infinity
        else:
            return min( [ c.valuation() for c in col ] )

    def precision_lattice(self, elements=None):
        r"""
        Return a matrix representing the precision lattice on a
        subset of elements.

        INPUT:

        - ``elements`` -- list of elements or ``None`` (default: ``None``)

        EXAMPLES::

            sage: R = ZpLF(2, label='preclattice')
            sage: prec = R.precision()
            sage: x = R(1, 10); y = R(1, 5)
            sage: prec.precision_lattice()                                              # needs sage.geometry.polyhedron
            [1024    0]
            [   0   32]

            sage: u = x + y
            sage: v = x - y
            sage: prec.precision_lattice([u, v])                                        # needs sage.geometry.polyhedron
            [  32 2016]
            [   0 2048]

        If the precision module does not project to a lattice,
        an error is raised. ::

            sage: prec.precision_lattice([x, y, u, v])                                  # needs sage.geometry.polyhedron
            Traceback (most recent call last):
            ...
            PrecisionError: the differential is not surjective

        Here is another example with matrices::

            sage: M = matrix(R, 2, 2, [R(3, 5), R(7, 5), R(1, 5), R(11, 1)])            # needs sage.modules
            sage: N = M^10                                                              # needs sage.modules

        The next syntax provides as easy way to select an interesting
        subset of variables (the selected subset consists of the four
        entries of the matrix ``N``)::

            sage: prec.precision_lattice(N)                                             # needs sage.geometry.polyhedron sage.modules
            [  2048    512  28160 230400]
            [     0   2048  14336 258048]
            [     0      0  65536  65536]
            [     0      0      0 262144]
        """
        if elements is None:
            elements = self._elements
        else:
            elements = list_of_padics(elements)
        n = len(self._elements)
        rows = [ ]
        val = 0
        for ref in elements:
            col = self._matrix[ref]
            row = [ x.value() for x in col ]
            valcol = min([ x.valuation() for x in col ])
            val = min(valcol, val)
            row += (n-len(row)) * [ZZ(0)]
            rows.append(row)
        from sage.matrix.constructor import matrix
        M = matrix(rows).transpose()
        if val < 0:
            M *= self._p ** (-val)
        M = M.change_ring(ZZ)
        M.echelonize()
        n = len(elements)
        if len(M.pivots()) < n:
            raise PrecisionError("the differential is not surjective")
        for i in range(n):
            v = M[i, i].valuation(self._p)
            M[i, i] = self._p ** v
        M.echelonize()
        M = M.submatrix(0, 0, n, n)
        if val < 0:
            M *= self._p ** val
        return M


class pAdicLatticeElementWeakProxy:
    r"""
    The implementations of :class:`DifferentialPrecisionGeneric` hold
    weak references to :class:`pAdicLatticeElement`. They are stored in
    dictionaries, e.g., a dictionary that maps an element to the corresponding
    column in the precision lattice matrix.
    However, weak references as implemented by Python are tricky to use as
    dictionary keys. Their equality depends on the equality of the element they
    point to (as long as that element is alive) and then on the equality by
    ``id``. This means that statements such as: ``ref in D == ref in D`` could
    be false if the garbage collector kicks in between the two invocations.
    To prevent very subtle and hardly reproducible bugs, we wrap weak
    references in a proxy that gives every lattice element a unique increasing
    id and uses that id for comparisons.

    EXAMPLES:

    Proxy elements exist only internally and are not usually exposed to the user::

        sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
        sage: R = ZpLF(2, label='proxy')
        sage: p = R(2)
        sage: prec = R.precision()
        sage: proxy = prec._elements[0]
        sage: isinstance(proxy, pAdicLatticeElementWeakProxy)
        True
    """
    _next_id = 0

    def __init__(self, element, callback=None):
        r"""
        TESTS::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLF(2, label='proxy')
            sage: p = R(2)
            sage: pAdicLatticeElementWeakProxy(p) == pAdicLatticeElementWeakProxy(p)
            True
            sage: pAdicLatticeElementWeakProxy(p) is pAdicLatticeElementWeakProxy(p)
            False
        """
        if not hasattr(element, '_proxy_id'):
            element._proxy_id = pAdicLatticeElementWeakProxy._next_id
            pAdicLatticeElementWeakProxy._next_id += 1
        self._id = element._proxy_id
        from weakref import ref
        proxy_callback = callback
        if callback is not None:
            proxy_callback = lambda _: callback(self)
        self._weakref = ref(element, proxy_callback)

    def __hash__(self):
        r"""
        Return a hash value for this proxy.

        EXAMPLES::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLF(2, label='proxy')
            sage: p = R(2)
            sage: hash(pAdicLatticeElementWeakProxy(p)) == hash(pAdicLatticeElementWeakProxy(p))
            True
        """
        return self._id

    def __eq__(self, other):
        r"""
        Return whether this proxy is undistinguishable from ``other``.

        EXAMPLES::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLF(2, label='proxy')
            sage: p = R(2)
            sage: q = R(2)
            sage: pAdicLatticeElementWeakProxy(p) == pAdicLatticeElementWeakProxy(p)
            True
            sage: pAdicLatticeElementWeakProxy(q) == pAdicLatticeElementWeakProxy(p)
            False
        """
        return isinstance(other, pAdicLatticeElementWeakProxy) and self._id == other._id

    def __call__(self):
        r"""
        Return the lattice element this proxy points to, or ``None`` if the
        target has already been finalized.

        EXAMPLES::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLF(2, label='proxy')
            sage: p = R(2)
            sage: pAdicLatticeElementWeakProxy(p)()
            2 + O(2^21)
        """
        return self._weakref()

    def __repr__(self):
        r"""
        Return a printable representation of this proxy.

        EXAMPLES::

            sage: from sage.rings.padics.lattice_precision import pAdicLatticeElementWeakProxy
            sage: R = ZpLF(2, label='proxy_repr')
            sage: p = R(2)
            sage: R.precision()._elements  # indirect doctest
            [WeakProxy#...]
        """
        return "WeakProxy#%s" % (self._id,)


def list_of_padics(elements):
    r"""
    Convert a list of `p`-adic composed elements (such as polynomials, matrices)
    to a list of weak references of their `p`-adic coefficients.

    This is a helper function for the method :meth:`precision_lattice`.

    TESTS::

        sage: from sage.rings.padics.lattice_precision import list_of_padics
        sage: R = ZpLC(2)
        sage: M = random_matrix(R, 2, 2)                                                # needs sage.geometry.polyhedron
        sage: list_of_padics(M)                                                         # needs sage.geometry.polyhedron
        [WeakProxy#...,
         WeakProxy#...,
         WeakProxy#...,
         WeakProxy#...]
    """
    from sage.rings.padics.padic_lattice_element import pAdicLatticeElement
    if isinstance(elements, pAdicLatticeElement):
        return [ pAdicLatticeElementWeakProxy(elements) ]
    try:
        if elements.parent().is_sparse():
            elements = elements.coefficients()
    except AttributeError:
        pass
    if not isinstance(elements, list):
        elements = list(elements)
    ans = [ ]
    for x in elements:
        ans += list_of_padics(x)
    return ans
