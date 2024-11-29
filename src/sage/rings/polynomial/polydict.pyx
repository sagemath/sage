"""
Generic data structures for multivariate polynomials

This module provides an implementation of a generic data structure
:class:`PolyDict` and the underlying arithmetic for multi-variate polynomial
rings. It uses a sparse representation of polynomials encoded as a Python
dictionary where keys are exponents and values coefficients.

``{(e1,...,er):c1,...} <-> c1*x1^e1*...*xr^er+...``,

The exponent ``(e1,...,er)`` in this representation is an instance of the class
:class:`ETuple`.

AUTHORS:

- William Stein
- David Joyner
- Martin Albrecht (ETuple)
- Joel B. Mohler (2008-03-17) -- ETuple rewrite as sparse C array
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2022 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.string cimport memcpy
from cpython.dict cimport *
cimport cython
from cpython.object cimport (Py_EQ, Py_NE, Py_LT, Py_LE, Py_GT, Py_GE)
from cysignals.memory cimport sig_malloc, sig_free

from sage.structure.richcmp cimport rich_to_bool

from functools import reduce
from pprint import pformat

from sage.misc.latex import latex


cpdef int gen_index(PolyDict x) noexcept:
    r"""
    Return the index of the variable represented by ``x`` or ``-1`` if ``x``
    is not a monomial of degree one.

    EXAMPLES::

        sage: from sage.rings.polynomial.polydict import PolyDict, gen_index
        sage: gen_index(PolyDict({(1, 0): 1}))
        0
        sage: gen_index(PolyDict({(0, 1): 1}))
        1
        sage: gen_index(PolyDict({}))
        -1
    """
    if len(x.__repn) != 1:
        return -1
    cdef ETuple e = next(iter(x.__repn))
    if e._nonzero != 1 or e._data[1] != 1:
        return -1
    if not next(iter(x.__repn.values())).is_one():
        return -1
    return e._data[0]


cpdef ETuple monomial_exponent(PolyDict p):
    r"""
    Return the unique exponent of ``p`` if it is a monomial or raise a
    :exc:`ValueError`.

    EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict, monomial_exponent
            sage: monomial_exponent(PolyDict({(2, 3): 1}))
            (2, 3)
            sage: monomial_exponent(PolyDict({(2, 3): 3}))
            Traceback (most recent call last):
            ...
            ValueError: not a monomial
            sage: monomial_exponent(PolyDict({(1, 0): 1, (0, 1): 1}))
            Traceback (most recent call last):
            ...
            ValueError: not a monomial
    """
    if len(p.__repn) != 1 or not next(iter(p.__repn.values())).is_one():
        raise ValueError('not a monomial')
    return next(iter(p.__repn))


cdef class PolyDict:
    r"""
    Data structure for multivariate polynomials.

    A PolyDict holds a dictionary all of whose keys are :class:`ETuple` and
    whose values are coefficients on which it is implicitely assumed that
    arithmetic operations can be performed.

    No arithmetic operation on :class:`PolyDict` clear zero coefficients as of
    now there is no reliable way of testing it in the most general setting, see
    :issue:`35319`. For removing zero coefficients from a :class:`PolyDict` you
    can use the method :meth:`remove_zeros` which can be parametrized by a zero
    test.
    """
    def __init__(self, pdict, zero=None, remove_zero=None, force_int_exponents=None, force_etuples=None, bint check=True):
        """
        INPUT:

        - ``pdict`` -- dictionary or list, which represents a multi-variable
          polynomial with the distribute representation (a copy is made)

        - ``zero`` -- deprecated

        - ``remove_zero`` -- deprecated

        - ``force_int_exponents`` -- deprecated

        - ``force_etuples`` -- deprecated

        - ``check`` -- if set to ``False`` then assumes that the exponents are
          all valid ``ETuple``; in that case the construction is a bit faster

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            PolyDict with representation {(1, 2): 3, (2, 1): 4, (2, 3): 2}

            sage: PolyDict({(2, 3): 0, (1, 2): 3, (2, 1): 4})
            PolyDict with representation {(1, 2): 3, (2, 1): 4, (2, 3): 0}

            sage: PolyDict({(0, 0): RIF(-1,1)})                                         # needs sage.rings.real_interval_field
            PolyDict with representation {(0, 0): 0.?}

        TESTS::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: len(f)
            3
            sage: f = PolyDict({}, zero=3, force_int_exponents=True, force_etuples=True)
            doctest:warning
            ...
            DeprecationWarning: the arguments "zero", "forced_int_exponents"
            and "forced_etuples" of PolyDict constructor are deprecated
            See https://github.com/sagemath/sage/issues/34000 for details.
            sage: f = PolyDict({}, remove_zero=False)
            doctest:warning
            ...
            DeprecationWarning: the argument "remove_zero" of PolyDict
            constructor is deprecated; call the method remove_zeros
            See https://github.com/sagemath/sage/issues/34000 for details.
        """
        if zero is not None or force_int_exponents is not None or force_etuples is not None:
            from sage.misc.superseded import deprecation
            deprecation(34000, 'the arguments "zero", "forced_int_exponents" and "forced_etuples" of PolyDict constructor are deprecated')

        self.__repn = {}
        cdef bint has_only_etuple = True
        if isinstance(pdict, (tuple, list)):
            for coeff, exp in pdict:
                if check and type(exp) is not ETuple:
                    exp = ETuple(exp)
                self.__repn[exp] = coeff
        elif isinstance(pdict, dict):
            if check:
                for k in (<dict> pdict):
                    if type(k) is not ETuple:
                        has_only_etuple = False
                        break
            if has_only_etuple:
                self.__repn = (<dict> pdict).copy()
            else:
                self.__repn = {}
                for exp, coeff in pdict.items():
                    if type(exp) is not ETuple:
                        exp = ETuple(exp)
                    self.__repn[exp] = coeff
        else:
            raise TypeError("pdict must be a dict or a list of pairs (coeff, exponent)")

        if remove_zero is not None:
            from sage.misc.superseded import deprecation
            deprecation(34000, 'the argument "remove_zero" of PolyDict constructor is deprecated; call the method remove_zeros')
            if remove_zero:
                self.remove_zeros()

    cdef PolyDict _new(self, dict pdict):
        cdef PolyDict ans = PolyDict.__new__(PolyDict)
        ans.__repn = pdict
        return ans

    cpdef remove_zeros(self, zero_test=None):
        r"""
        Remove the entries with zero coefficients.

        INPUT:

        - ``zero_test`` -- (optional) function that performs test to zero of a coefficient

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3):0})
            sage: f
            PolyDict with representation {(2, 3): 0}
            sage: f.remove_zeros()
            sage: f
            PolyDict with representation {}

        The following example shows how to remove only exact zeros from a ``PolyDict``
        containing univariate power series::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = PolyDict({(1, 1): O(t), (1, 0): R.zero()})
            sage: f.remove_zeros(lambda s: s.is_zero() and s.prec() is Infinity)
            sage: f
            PolyDict with representation {(1, 1): O(t^1)}
        """
        if not self.__repn:
            return
        # NOTE: in each of the conditional statements below, what the first
        # loop does is equivalent to
        #
        #     if all(self.__repn.values()):
        #         return
        #
        # and
        #
        #     if all(not zero_test(coeff) for coeff in self.__repn.values()):
        #         return
        #
        # However, 'all(...)' is badly handled by the Cython compiler and we
        # rather unfold it for efficiency.
        cdef bint has_zero_coefficient = False
        if zero_test is None:
            for coeff in self.__repn.values():
                if not coeff:
                    has_zero_coefficient = True
                    break
            if not has_zero_coefficient:
                return
            for k in list(self.__repn):
                if not self.__repn[k]:
                    del self.__repn[k]
        else:
            for coeff in self.__repn.values():
                if zero_test(coeff):
                    has_zero_coefficient = True
                    break
            if not has_zero_coefficient:
                return
            for k in list(self.__repn):
                if zero_test(self.__repn[k]):
                    del self.__repn[k]

    def apply_map(self, f):
        r"""
        Apply the map ``f`` on the coefficients (inplace).

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(1, 0): 1, (1, 1): -2})
            sage: f.apply_map(lambda x: x^2)
            sage: f
            PolyDict with representation {(1, 0): 1, (1, 1): 4}
        """
        for k, v in self.__repn.items():
            self.__repn[k] = f(v)

    def coerce_coefficients(self, A):
        r"""
        Coerce the coefficients in the parent ``A``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 0})
            sage: f
            PolyDict with representation {(2, 3): 0}
            sage: f.coerce_coefficients(QQ)
            doctest:warning
            ...
            DeprecationWarning: coerce_cefficients is deprecated; use apply_map instead
            See https://github.com/sagemath/sage/issues/34000 for details.
            sage: f
            PolyDict with representation {(2, 3): 0}
        """
        from sage.misc.superseded import deprecation
        deprecation(34000, 'coerce_cefficients is deprecated; use apply_map instead')
        self.apply_map(A.coerce)

    def __hash__(self):
        """
        Return the hash.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PD1 = PolyDict({(2, 3): 0, (1, 2): 3, (2, 1): 4})
            sage: PD2 = PolyDict({(2, 3): 0, (1, 2): 3, (2, 1): 4})
            sage: PD3 = PolyDict({(2, 3): 1, (1, 2): 3, (2, 1): 4})
            sage: hash(PD1) == hash(PD2)
            True
            sage: hash(PD1) == hash(PD3)
            False
        """
        return hash(frozenset(self.__repn.items()))

    def __bool__(self):
        """
        Return whether the PolyDict is empty.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PD1 = PolyDict({(2, 3): 0, (1, 2): 3, (2, 1): 4})
            sage: bool(PD1)
            True
        """
        return bool(self.__repn)

    def __len__(self):
        """
        Return the number of terms of this polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PD1 = PolyDict({(2, 3): 0, (1, 2): 3, (2, 1): 4})
            sage: len(PD1)
            3
        """
        return len(self.__repn)

    def __richcmp__(PolyDict left, PolyDict right, int op):
        """
        Implement the ``__richcmp__`` protocol for `PolyDict`s.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: p1 = PolyDict({(0,): 1})
            sage: p2 = PolyDict({(0,): 2})
            sage: p1 == p2
            False
            sage: p1 < p2
            Traceback (most recent call last):
            ...
            TypeError: unsupported comparison between PolyDict
        """
        if op == Py_EQ:
            return left.__repn == right.__repn
        elif op == Py_NE:
            return left.__repn != right.__repn

        raise TypeError('unsupported comparison between PolyDict')

    def rich_compare(self, PolyDict other, int op, sortkey=None):
        """
        Compare two `PolyDict`s using a specified term ordering ``sortkey``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: from sage.structure.richcmp import op_EQ, op_NE, op_LT
            sage: p1 = PolyDict({(0,): 1})
            sage: p2 = PolyDict({(0,): 2})
            sage: O = TermOrder()
            sage: p1.rich_compare(PolyDict({(0,): 1}), op_EQ, O.sortkey)
            True
            sage: p1.rich_compare(p2, op_EQ, O.sortkey)
            False
            sage: p1.rich_compare(p2, op_NE, O.sortkey)
            True
            sage: p1.rich_compare(p2, op_LT, O.sortkey)
            True

            sage: p3 = PolyDict({(3, 2, 4): 1, (3, 2, 5): 2})
            sage: p4 = PolyDict({(3, 2, 4): 1, (3, 2, 3): 2})
            sage: p3.rich_compare(p4, op_LT, O.sortkey)
            False

        TESTS::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: from sage.structure.richcmp import op_EQ, op_NE, op_LT
            sage: p = PolyDict({})
            sage: ans = p.rich_compare(p, op_EQ)
            doctest:warning
            ...
            DeprecationWarning: the argument "sortkey" will become mandatory in future sage versions
            See https://github.com/sagemath/sage/issues/34000 for details.
            sage: ans
            True
        """
        if sortkey is None:
            from sage.misc.superseded import deprecation
            deprecation(34000, 'the argument "sortkey" will become mandatory in future sage versions')

        if op == Py_EQ:
            return self.__repn == other.__repn
        elif op == Py_NE:
            return self.__repn != other.__repn

        if sortkey is None:
            raise TypeError("ordering of PolyDicts requires a sortkey")

        # start with biggest
        cdef list left = sorted(self.__repn, key=sortkey, reverse=True)
        cdef list right = sorted(other.__repn, key=sortkey, reverse=True)

        cdef size_t i
        for i in range(min(len(left), len(right))):
            m = left[i]
            n = right[i]
            keym = sortkey(m)
            keyn = sortkey(n)

            # first compare the leading monomials
            if keym > keyn:
                return rich_to_bool(op, 1)
            elif keym < keyn:
                return rich_to_bool(op, -1)

            # same leading monomial, compare their coefficients
            coefm = self.__repn[m]
            coefn = other.__repn[n]
            if coefm > coefn:
                return rich_to_bool(op, 1)
            elif coefm < coefn:
                return rich_to_bool(op, -1)

        return rich_to_bool(op, (len(left) > len(right)) - (len(left) < len(right)))

    def list(self):
        """
        Return a list that defines ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: sorted(f.list())
            [[2, [2, 3]], [3, [1, 2]], [4, [2, 1]]]
        """
        return [[c, list(e)] for e, c in self.__repn.items()]

    def dict(self):
        """
        Return a copy of the dict that defines ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.dict()
            {(1, 2): 3, (2, 1): 4, (2, 3): 2}
        """
        return self.__repn.copy()

    def coefficients(self):
        """
        Return the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: sorted(f.coefficients())
            [2, 3, 4]
        """
        return list(self.__repn.values())

    def exponents(self):
        """
        Return the exponents of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: sorted(f.exponents())
            [(1, 2), (2, 1), (2, 3)]
        """
        return list(self.__repn)

    def __getitem__(self, e):
        """
        Return a coefficient of the polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f[1, 2]
            3
            sage: f[(2, 1)]
            4
        """
        if type(e) is not ETuple:
            e = ETuple(e)
        return self.__repn[e]

    def get(self, ETuple e, default=None):
        r"""
        Return the coefficient of the ETuple ``e`` if present and ``default`` otherwise.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict, ETuple
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.get(ETuple([1,2]))
            3
            sage: f.get(ETuple([1,1]), 'hello')
            'hello'
        """
        return self.__repn.get(e, default)

    def __repr__(self):
        r"""
        String representation.
        """
        repn = ' '.join(pformat(self.__repn).splitlines())
        return 'PolyDict with representation %s' % repn

    def degree(self, PolyDict x=None):
        r"""
        Return the total degree or the maximum degree in the variable ``x``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.degree()
            5
            sage: f.degree(PolyDict({(1, 0): 1}))
            2
            sage: f.degree(PolyDict({(0, 1): 1}))
            3
        """
        if x is None:
            return self.total_degree()
        cdef int i = gen_index(x)
        if i < 0:
            raise ValueError('x must be a generator')
        if not self.__repn:
            return -1
        return max((<ETuple> e).get_exp(i) for e in self.__repn)

    def total_degree(self, tuple w=None):
        r"""
        Return the total degree.

        INPUT:

        - ``w`` -- (optional) a tuple of weights

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.total_degree()
            5
            sage: f.total_degree((3, 1))
            9
            sage: PolyDict({}).degree()
            -1
        """
        if not self.__repn:
            return -1
        if w is None:
            return max((<ETuple> e).unweighted_degree() for e in self.__repn)
        else:
            return max((<ETuple> e).weighted_degree(w) for e in self.__repn)

    def monomial_coefficient(self, mon):
        """
        Return the coefficient of the monomial ``mon``.

        INPUT:

        - ``mon`` -- a PolyDict with a single key

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.monomial_coefficient(PolyDict({(2,1):1}).dict())
            doctest:warning
            ...
            DeprecationWarning: PolyDict.monomial_coefficient is deprecated; use PolyDict.get instead
            See https://github.com/sagemath/sage/issues/34000 for details.
            4
        """
        from sage.misc.superseded import deprecation
        deprecation(34000, 'PolyDict.monomial_coefficient is deprecated; use PolyDict.get instead')
        K, = mon.keys()
        if K not in self.__repn:
            return 0
        return self.__repn[K]

    def polynomial_coefficient(self, degrees):
        """
        Return a polydict that defines the coefficient in the current
        polynomial viewed as a tower of polynomial extensions.

        INPUT:

        - ``degrees`` -- list of degree restrictions; list elements are ``None``
          if the variable in that position should be unrestricted

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.polynomial_coefficient([2, None])
            PolyDict with representation {(0, 1): 4, (0, 3): 2}
            sage: f = PolyDict({(0, 3): 2, (0, 2): 3, (2, 1): 4})
            sage: f.polynomial_coefficient([0, None])
            PolyDict with representation {(0, 2): 3, (0, 3): 2}
        """
        nz = []
        cdef int i
        for i in range(len(degrees)):
            if degrees[i] is not None:
                nz.append(i)
        cdef dict ans = {}
        cdef bint exactly_divides
        for S in self.__repn:
            exactly_divides = True
            for j in nz:
                if S[j] != degrees[j]:
                    exactly_divides = False
                    break
            if exactly_divides:
                t = list(S)
                for m in nz:
                    t[m] = 0
                ans[ETuple(t)] = self.__repn[S]
        return self._new(ans)

    def coefficient(self, mon):
        """
        Return a polydict that defines a polynomial in 1 less number
        of variables that gives the coefficient of mon in this
        polynomial.

        The coefficient is defined as follows.  If f is this
        polynomial, then the coefficient is the sum T/mon where the
        sum is over terms T in f that are exactly divisible by mon.
        """
        K, = mon.keys()
        nz = K.nonzero_positions()  # set([i for i in range(len(K)) if K[i] != 0])
        ans = {}
        for S in self.__repn:
            exactly_divides = True
            for j in nz:
                if S[j] != K[j]:
                    exactly_divides = False
                    break
            if exactly_divides:
                t = list(S)
                for m in nz:
                    t[m] = 0
                ans[ETuple(t)] = self.__repn[S]
        return self._new(ans)

    def is_homogeneous(self, tuple w=None):
        r"""
        Return whether this polynomial is homogeneous.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PolyDict({}).is_homogeneous()
            True
            sage: PolyDict({(1, 2): 1, (0, 3): -2}).is_homogeneous()
            True
            sage: PolyDict({(1, 0): 1, (1, 2): 3}).is_homogeneous()
            False
        """
        if not self.__repn:
            return True
        cdef size_t s
        it = iter(self.__repn)
        if w is None:
            s = (<ETuple> next(it)).unweighted_degree()
            for elt in it:
                if (<ETuple> elt).unweighted_degree() != s:
                    return False
            return True
        else:
            s = (<ETuple> next(it)).weighted_degree(w)
            for elt in it:
                if (<ETuple> elt).weighted_degree(w) != s:
                    return False
            return True

    def is_constant(self):
        """
        Return whether this polynomial is constant.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.is_constant()
            False
            sage: g = PolyDict({(0, 0): 2})
            sage: g.is_constant()
            True
            sage: h = PolyDict({})
            sage: h.is_constant()
            True
        """
        if not self.__repn:
            return True
        if len(self.__repn) > 1:
            return False
        return not any(self.__repn)

    def homogenize(self, size_t var):
        r"""
        Return the homogeneization of ``self`` by increasing the degree of the
        variable ``var``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(0, 0): 1, (2, 1): 3, (1, 1): 5})
            sage: f.homogenize(0)
            PolyDict with representation {(2, 1): 8, (3, 0): 1}
            sage: f.homogenize(1)
            PolyDict with representation {(0, 3): 1, (1, 2): 5, (2, 1): 3}

            sage: PolyDict({(0, 1): 1, (1, 1): -1}).homogenize(0)
            PolyDict with representation {(1, 1): 0}
        """
        cdef dict H = {}
        cdef int deg = self.degree()
        cdef int shift
        for e, val in self.__repn.items():
            shift = deg - (<ETuple> e).unweighted_degree()
            if shift:
                f = (<ETuple> e).eadd_p(shift, var)
            else:
                f = e
            if f in H:
                H[f] += val
            else:
                H[f] = val
        return self._new(H)

    def latex(self, vars, atomic_exponents=True,
              atomic_coefficients=True, sortkey=None):
        r"""
        Return a nice polynomial latex representation of this PolyDict, where
        the vars are substituted in.

        INPUT:

        - ``vars`` -- list
        - ``atomic_exponents`` -- boolean (default: ``True``)
        - ``atomic_coefficients`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.latex(['a', 'WW'])
            '2 a^{2} WW^{3} + 4 a^{2} WW + 3 a WW^{2}'

        TESTS:

        We check that the issue on :issue:`9478` is resolved::

            sage: R2.<a> = QQ[]
            sage: R3.<xi, x> = R2[]
            sage: print(latex(xi*x))
            \xi x

        TESTS:

        Check that :issue:`29604` is fixed::

            sage: PolyDict({(1, 0): GF(2)(1)}).latex(['x', 'y'])
            'x'
        """
        if not self:
            return "0"

        poly = ""

        sort_kwargs = {'reverse': True}
        if sortkey:
            sort_kwargs['key'] = sortkey

        E = sorted(self.__repn, **sort_kwargs)

        try:
            ring = self.__repn[E[0]].parent()
            pos_one = ring.one()
            neg_one = -pos_one
        except (AttributeError, ArithmeticError):
            # AritchmeticError occurs when self.__repn[E[0]] is a tropical
            # semiring element
            # probably self.__repn[E[0]] is not a ring element
            pos_one = 1
            neg_one = -1

        is_characteristic_2 = bool(pos_one == neg_one)

        for e in E:
            c = self.__repn[e]
            sign_switch = False
            # First determine the multinomial:
            multi = " ".join([vars[j] +
                              ("^{%s}" % e[j] if e[j] != 1 else "")
                             for j in e.nonzero_positions(sort=True)])
            # Next determine coefficient of multinomial
            if len(multi) == 0:
                multi = latex(c)
            elif c == neg_one and not is_characteristic_2:
                # handle -1 specially because it's a pain
                if len(poly) > 0:
                    sign_switch = True
                else:
                    multi = "-%s" % multi
            elif c != pos_one:
                c = latex(c)
                if (not atomic_coefficients and multi and
                        ('+' in c or '-' in c or ' ' in c)):
                    c = "\\left(%s\\right)" % c
                multi = "%s %s" % (c, multi)

            # Now add on coefficiented multinomials
            if len(poly) > 0:
                if sign_switch:
                    poly = poly + " - "
                else:
                    poly = poly + " + "
            poly = poly + multi
        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -", "- ")
        if len(poly) == 0:
            return "0"
        return poly

    def poly_repr(self, vars, atomic_exponents=True,
                  atomic_coefficients=True, sortkey=None):
        """
        Return a nice polynomial string representation of this PolyDict, where
        the vars are substituted in.

        INPUT:

        - ``vars`` -- list
        - ``atomic_exponents`` -- boolean (default: ``True``)
        - ``atomic_coefficients`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.poly_repr(['a', 'WW'])
            '2*a^2*WW^3 + 4*a^2*WW + 3*a*WW^2'

        We check to make sure that when we are in characteristic two, we
        don't put negative signs on the generators. ::

            sage: Integers(2)['x, y'].gens()
            (x, y)

        We make sure that intervals are correctly represented. ::

            sage: f = PolyDict({(2, 3): RIF(1/2,3/2), (1, 2): RIF(-1,1)})               # needs sage.rings.real_interval_field
            sage: f.poly_repr(['x', 'y'])                                               # needs sage.rings.real_interval_field
            '1.?*x^2*y^3 + 0.?*x*y^2'

        TESTS:

        Check that :issue:`29604` is fixed::

            sage: PolyDict({(1, 0): GF(4)(1)}).poly_repr(['x', 'y'])                    # needs sage.rings.finite_rings
            'x'

            sage: # needs sage.modules
            sage: P.<x,y> = LaurentPolynomialRing(GF(2), 2)
            sage: P.gens()
            (x, y)
            sage: -x - y
            x + y
        """
        poly = ""
        sort_kwargs = {'reverse': True}
        if sortkey:
            sort_kwargs['key'] = sortkey

        E = sorted(self.__repn, **sort_kwargs)

        if not E:
            return "0"
        try:
            ring = self.__repn[E[0]].parent()
            pos_one = ring.one()
            neg_one = -pos_one
        except (AttributeError, ArithmeticError):
            # AritchmeticError occurs when self.__repn[E[0]] is a tropical
            # semiring element
            # probably self.__repn[E[0]] is not a ring element
            pos_one = 1
            neg_one = -1

        is_characteristic_2 = bool(pos_one == neg_one)

        for e in E:
            c = self.__repn[e]
            sign_switch = False
            # First determine the multinomial:
            multi = ""
            for j in e.nonzero_positions(sort=True):
                if multi:
                    multi = multi + "*"
                multi = multi + vars[j]
                if e[j] != 1:
                    if atomic_exponents:
                        multi = multi + "^%s" % e[j]
                    else:
                        multi = multi + "^(%s)" % e[j]
            # Next determine coefficient of multinomial
            if not multi:
                multi = str(c)
            elif c == neg_one and not is_characteristic_2:
                # handle -1 specially because it's a pain
                if poly:
                    sign_switch = True
                else:
                    multi = "-%s" % multi
            elif not c == pos_one:
                if not atomic_coefficients:
                    c = str(c)
                    if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                        c = "(%s)" % c
                multi = "%s*%s" % (c, multi)

            # Now add on coefficiented multinomials
            if poly:
                if sign_switch:
                    poly = poly + " - "
                else:
                    poly = poly + " + "
            poly = poly + multi
        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -", "- ")
        if not poly:
            return "0"
        return poly

    def __iadd__(PolyDict self, PolyDict other):
        r"""
        Inplace addition.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: g = PolyDict({(1, 5): -3, (2, 3): -2, (1, 1): 3})
            sage: f += g
            sage: f
            PolyDict with representation {(1, 1): 3, (1, 2): 3, (1, 5): -3, (2, 1): 4, (2, 3): 0}
        """
        cdef dict D = self.__repn
        if self is other:
            for e in D:
                v = D[e]
                D[e] += v
        else:
            for e, c in other.__repn.items():
                if e in D:
                    D[e] += c
                else:
                    D[e] = c
        return self

    def __neg__(self):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: -PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            PolyDict with representation {(1, 2): -3, (2, 1): -4, (2, 3): -2}
        """
        cdef dict D = self.__repn.copy()
        for e, c in D.items():
            D[e] = -c
        return self._new(D)

    def __add__(PolyDict self, PolyDict other):
        """
        Add two PolyDict's in the same number of variables.

        EXAMPLES:

        We add two polynomials in 2 variables::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: g = PolyDict({(1, 5): -3, (2, 3): -2, (1, 1): 3})
            sage: f + g
            PolyDict with representation {(1, 1): 3, (1, 2): 3, (1, 5): -3, (2, 1): 4, (2, 3): 0}

            sage: K = GF(2)
            sage: f = PolyDict({(1, 1): K(1)})
            sage: f + f
            PolyDict with representation {(1, 1): 0}
        """
        cdef dict D = self.__repn
        cdef dict R = other.__repn
        if len(D) < len(R):
            D, R = R, D
        D = D.copy()
        for e, c in R.items():
            if e in D:
                D[e] += c
            else:
                D[e] = c
        return self._new(D)

    def __mul__(PolyDict self, PolyDict right):
        """
        Multiply two PolyDict's in the same number of variables.

        The algorithm do not test whether a product of coefficients is zero
        or whether a final coefficient is zero because there is no reliable way
        to do so in general (eg power series ring or `p`-adic rings).

        EXAMPLES:

        Multiplication of polynomials in 2 variables with rational coefficients::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: g = PolyDict({(1, 5): -3, (2, 3): -2, (1, 1): 3})
            sage: f * g
            PolyDict with representation {(2, 3): 9, (2, 7): -9, (3, 2): 12, (3, 4): 6, (3, 5): -6, (3, 6): -12, (3, 8): -6, (4, 4): -8, (4, 6): -4}

        Multiplication of polynomials in 2 variables with power series coefficients::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = PolyDict({(1, 0): 1 + O(t), (0, 1): 1 + O(t^2)})
            sage: g = PolyDict({(1, 0): 1, (0, 1): -1})
            sage: f * g
            PolyDict with representation {(0, 2): -1 + O(t^2), (1, 1): O(t^1), (2, 0): 1 + O(t)}
            sage: f = PolyDict({(1, 0): O(t), (0, 1): O(t)})
            sage: g = PolyDict({(1, 0): 1, (0, 1): O(t)})
            sage: f * g
            PolyDict with representation {(0, 2): O(t^2), (1, 1): O(t^1), (2, 0): O(t^1)}
        """
        cdef PyObject *cc
        cdef dict newpoly
        if not self.__repn or not right.__repn:
            return self._new({})
        newpoly = {}
        for e0, c0 in self.__repn.items():
            for e1, c1 in right.__repn.items():
                e = (<ETuple> e0).eadd(<ETuple> e1)
                c = c0 * c1
                cc = PyDict_GetItem(newpoly, e)
                if cc == <PyObject*> 0:
                    PyDict_SetItem(newpoly, e, c)
                else:
                    PyDict_SetItem(newpoly, e, <object> cc + c)
        return self._new(newpoly)

    def scalar_rmult(self, s):
        """
        Return the right scalar multiplication of ``self`` by ``s``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict

            sage: x, y = FreeMonoid(2, 'x, y').gens()  # a strange object to live in a polydict, but non-commutative!   # needs sage.combinat
            sage: f = PolyDict({(2, 3): x})                                             # needs sage.combinat
            sage: f.scalar_rmult(y)                                                     # needs sage.combinat
            PolyDict with representation {(2, 3): x*y}

            sage: f = PolyDict({(2,3):2, (1, 2): 3, (2, 1): 4})
            sage: f.scalar_rmult(-2)
            PolyDict with representation {(1, 2): -6, (2, 1): -8, (2, 3): -4}
            sage: f.scalar_rmult(RIF(-1,1))                                             # needs sage.rings.real_interval_field
            PolyDict with representation {(1, 2): 0.?e1, (2, 1): 0.?e1, (2, 3): 0.?e1}
        """
        cdef dict v = {}
        for e, c in self.__repn.items():
            v[e] = c * s
        return self._new(v)

    def scalar_lmult(self, s):
        """
        Return the left scalar multiplication of ``self`` by ``s``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict

            sage: x, y = FreeMonoid(2, 'x, y').gens()  # a strange object to live in a polydict, but non-commutative!   # needs sage.combinat
            sage: f = PolyDict({(2,3):x})                                               # needs sage.combinat
            sage: f.scalar_lmult(y)                                                     # needs sage.combinat
            PolyDict with representation {(2, 3): y*x}

            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.scalar_lmult(-2)
            PolyDict with representation {(1, 2): -6, (2, 1): -8, (2, 3): -4}
            sage: f.scalar_lmult(RIF(-1,1))                                             # needs sage.rings.real_interval_field
            PolyDict with representation {(1, 2): 0.?e1, (2, 1): 0.?e1, (2, 3): 0.?e1}
        """
        cdef dict v = {}
        for e, c in self.__repn.items():
            v[e] = s * c
        return self._new(v)

    def term_lmult(self, exponent, s):
        """
        Return this element multiplied by ``s`` on the left and with exponents
        shifted by ``exponent``.

        INPUT:

        - ``exponent`` -- a ETuple

        - ``s`` -- a scalar

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple, PolyDict

            sage: x, y = FreeMonoid(2, 'x, y').gens()  # a strange object to live in a polydict, but non-commutative!   # needs sage.combinat
            sage: f = PolyDict({(2, 3): x})                                             # needs sage.combinat
            sage: f.term_lmult(ETuple((1, 2)), y)                                       # needs sage.combinat
            PolyDict with representation {(3, 5): y*x}

            sage: f = PolyDict({(2,3): 2, (1,2): 3, (2,1): 4})
            sage: f.term_lmult(ETuple((1, 2)), -2)
            PolyDict with representation {(2, 4): -6, (3, 3): -8, (3, 5): -4}
        """
        cdef dict v = {}
        for e, c in self.__repn.items():
            v[(<ETuple> e).eadd(exponent)] = s*c
        return self._new(v)

    def term_rmult(self, exponent, s):
        """
        Return this element multiplied by ``s`` on the right and with exponents
        shifted by ``exponent``.

        INPUT:

        - ``exponent`` -- a ETuple

        - ``s`` -- a scalar

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple, PolyDict

            sage: x, y = FreeMonoid(2, 'x, y').gens()  # a strange object to live in a polydict, but non-commutative!   # needs sage.combinat
            sage: f = PolyDict({(2, 3): x})                                             # needs sage.combinat
            sage: f.term_rmult(ETuple((1, 2)), y)                                       # needs sage.combinat
            PolyDict with representation {(3, 5): x*y}

            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.term_rmult(ETuple((1, 2)), -2)
            PolyDict with representation {(2, 4): -6, (3, 3): -8, (3, 5): -4}
        """
        cdef dict v = {}
        for e, c in self.__repn.items():
            v[(<ETuple> e).eadd(exponent)] = c * s
        return self._new(v)

    def __sub__(PolyDict self, PolyDict other):
        """
        Subtract two PolyDict's.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(2,3):2, (1,1):-10})
            sage: f - g
            PolyDict with representation {(1, 1): 10, (1, 2): 3, (2, 1): 4, (2, 3): 0}
            sage: g - f
            PolyDict with representation {(1, 1): -10, (1, 2): -3, (2, 1): -4, (2, 3): 0}
            sage: f - f
            PolyDict with representation {(1, 2): 0, (2, 1): 0, (2, 3): 0}
        """
        cdef dict D = self.__repn.copy()
        cdef dict R = other.__repn
        for e, c in R.items():
            if e in D:
                D[e] -= c
            else:
                D[e] = -c
        return self._new(D)

    def derivative_i(self, size_t i):
        r"""
        Return the derivative of ``self`` with respect to the ``i``-th variable.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PolyDict({(1, 1): 1}).derivative_i(0)
            PolyDict with representation {(0, 1): 1}
        """
        cdef dict D = {}
        for e, c in self.__repn.items():
            D[(<ETuple> e).eadd_p(-1, i)] = c * (<ETuple> e).get_exp(i)
        return self._new(D)

    def derivative(self, PolyDict x):
        r"""
        Return the derivative of ``self`` with respect to ``x``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.derivative(PolyDict({(1, 0): 1}))
            PolyDict with representation {(0, 2): 3, (1, 1): 8, (1, 3): 4}
            sage: f.derivative(PolyDict({(0, 1): 1}))
            PolyDict with representation {(1, 1): 6, (2, 0): 4, (2, 2): 6}

            sage: PolyDict({(-1,): 1}).derivative(PolyDict({(1,): 1}))
            PolyDict with representation {(-2,): -1}
            sage: PolyDict({(-2,): 1}).derivative(PolyDict({(1,): 1}))
            PolyDict with representation {(-3,): -2}

            sage: PolyDict({}).derivative(PolyDict({(1, 1): 1}))
            Traceback (most recent call last):
            ...
            ValueError: x must be a generator
        """
        cdef int i = gen_index(x)
        if i < 0:
            raise ValueError('x must be a generator')
        return self.derivative_i(i)

    def integral_i(self, size_t i):
        r"""
        Return the derivative of ``self`` with respect to the ``i``-th variable.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PolyDict({(1, 1): 1}).integral_i(0)
            PolyDict with representation {(2, 1): 1/2}
        """
        cdef dict D = {}
        cdef int exp
        for e, c in self.__repn.items():
            exp = (<ETuple> e).get_exp(i)
            if exp == -1:
                raise ArithmeticError('integral of monomial with exponent -1')
            D[(<ETuple> e).eadd_p(1, i)] = c / (exp + 1)
        return self._new(D)

    def integral(self, PolyDict x):
        r"""
        Return the integral of ``self`` with respect to ``x``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.integral(PolyDict({(1, 0): 1}))
            PolyDict with representation {(2, 2): 3/2, (3, 1): 4/3, (3, 3): 2/3}
            sage: f.integral(PolyDict({(0, 1): 1}))
            PolyDict with representation {(1, 3): 1, (2, 2): 2, (2, 4): 1/2}

            sage: PolyDict({(-1,): 1}).integral(PolyDict({(1,): 1}))
            Traceback (most recent call last):
            ...
            ArithmeticError: integral of monomial with exponent -1
            sage: PolyDict({(-2,): 1}).integral(PolyDict({(1,): 1}))
            PolyDict with representation {(-1,): -1}
            sage: PolyDict({}).integral(PolyDict({(1, 1): 1}))
            Traceback (most recent call last):
            ...
            ValueError: x must be a generator
        """
        cdef int i = gen_index(x)
        if i < 0:
            raise ValueError('x must be a generator')
        return self.integral_i(i)

    def lcmt(self, greater_etuple):
        """
        Provides functionality of lc, lm, and lt by calling the tuple
        compare function on the provided term order T.

        INPUT:

        - ``greater_etuple`` -- a term order
        """
        try:
            return ETuple(reduce(greater_etuple, self.__repn))
        except KeyError:
            raise ArithmeticError("%s not supported", greater_etuple)

    def __reduce__(self):
        """

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: loads(dumps(f)) == f
            True
        """
        return PolyDict, (self.__repn,)

    def min_exp(self):
        """
        Return an ETuple containing the minimum exponents appearing.  If
        there are no terms at all in the PolyDict, it returns None.

        The nvars parameter is necessary because a PolyDict doesn't know it
        from the data it has (and an empty PolyDict offers no clues).

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.min_exp()
            (1, 1)
            sage: PolyDict({}).min_exp() # returns None
        """
        cdef ETuple r
        if self.__repn:
            it = iter(self.__repn)
            r = next(it)
            for e in it:
                r = r.emin(e)
            return r
        else:
            return None

    def max_exp(self):
        """
        Return an ETuple containing the maximum exponents appearing.  If
        there are no terms at all in the PolyDict, it returns None.

        The nvars parameter is necessary because a PolyDict doesn't know it
        from the data it has (and an empty PolyDict offers no clues).

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2, 3): 2, (1, 2): 3, (2, 1): 4})
            sage: f.max_exp()
            (2, 3)
            sage: PolyDict({}).max_exp() # returns None
        """
        cdef ETuple r
        if self.__repn:
            it = iter(self.__repn)
            r = next(it)
            for e in it:
                r = r.emax(e)
            return r
        else:
            return None

cdef inline bint dual_etuple_iter(ETuple self, ETuple other, size_t *ind1, size_t *ind2, size_t *index, int *exp1, int *exp2) noexcept:
    """
    This function is a crucial helper function for a number of methods of
    the ETuple class.

    This is a rather fragile function.  Perhaps some Cython guru could make
    it appear a little less stilted -- a principal difficulty is passing
    C types by reference.  In any case, the complicated features of looping
    through two ETuple _data members is all packaged up right here and
    shouldn't be allowed to spread.
    """
    if ind1[0] >= self._nonzero and ind2[0] >= other._nonzero:
        return 0
    if ind1[0] < self._nonzero and ind2[0] < other._nonzero:
        if self._data[2*ind1[0]] == other._data[2*ind2[0]]:
            exp1[0] = self._data[2*ind1[0]+1]
            exp2[0] = other._data[2*ind2[0]+1]
            index[0] = self._data[2*ind1[0]]
            ind1[0] += 1
            ind2[0] += 1
        elif self._data[2*ind1[0]] > other._data[2*ind2[0]]:
            exp1[0] = 0
            exp2[0] = other._data[2*ind2[0]+1]
            index[0] = other._data[2*ind2[0]]
            ind2[0] += 1
        else:
            exp1[0] = self._data[2*ind1[0]+1]
            exp2[0] = 0
            index[0] = self._data[2*ind1[0]]
            ind1[0] += 1
    else:
        if ind2[0] >= other._nonzero:
            exp1[0] = self._data[2*ind1[0]+1]
            exp2[0] = 0
            index[0] = self._data[2*ind1[0]]
            ind1[0] += 1
        elif ind1[0] >= self._nonzero:
            exp1[0] = 0
            exp2[0] = other._data[2*ind2[0]+1]
            index[0] = other._data[2*ind2[0]]
            ind2[0] += 1
    return 1

cdef class ETuple:
    """
    Representation of the exponents of a polydict monomial. If
    (0,0,3,0,5) is the exponent tuple of x_2^3*x_4^5 then this class
    only stores {2:3, 4:5} instead of the full tuple. This sparse
    information may be obtained by provided methods.

    The index/value data is all stored in the _data C int array member
    variable.  For the example above, the C array would contain
    2,3,4,5.  The indices are interlaced with the values.

    This data structure is very nice to work with for some functions
    implemented in this class, but tricky for others.  One reason that
    I really like the format is that it requires a single memory
    allocation for all of the values.  A hash table would require more
    allocations and presumably be slower.  I didn't benchmark this
    question (although, there is no question that this is much faster
    than the prior use of python dicts).
    """
    cdef ETuple _new(self):
        """
        Quickly create a new initialized ETuple with the
        same length as ``self``.
        """
        cdef type t = type(self)
        cdef ETuple x = <ETuple>t.__new__(t)
        x._length = self._length
        return x

    def __init__(self, data=None, length=None):
        """
        - ``ETuple()`` -> an empty ETuple
        - ``ETuple(sequence)`` -> ETuple initialized from sequence's items

        If the argument is an ETuple, the return value is the same object.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1, 1, 0])
            (1, 1, 0)
            sage: ETuple({int(1): int(2)}, int(3))
            (0, 2, 0)
            sage: ETuple([1, -1, 0])
            (1, -1, 0)

        TESTS:

        Iterators are not accepted::

            sage: ETuple(iter([2, 3, 4]))
            Traceback (most recent call last):
            ...
            TypeError: Error in ETuple((), <list... object at ...>, None)
        """
        if data is None:
            return
        cdef size_t ind
        cdef int v
        if isinstance(data, ETuple):
            self._length = (<ETuple>data)._length
            self._nonzero = (<ETuple>data)._nonzero
            self._data = <int*>sig_malloc(sizeof(int)*self._nonzero*2)
            memcpy(self._data, (<ETuple>data)._data, sizeof(int)*self._nonzero*2)
        elif isinstance(data, dict) and isinstance(length, int):
            self._length = length
            self._nonzero = len(data)
            self._data = <int*>sig_malloc(sizeof(int)*self._nonzero*2)
            nz_elts = sorted(data.items())
            ind = 0
            for index, exp in nz_elts:
                self._data[2*ind] = index
                self._data[2*ind+1] = exp
                ind += 1
        elif isinstance(data, (list, tuple)):
            self._length = len(data)
            self._nonzero = 0
            for v in data:
                if v != 0:
                    self._nonzero += 1
            ind = 0
            self._data = <int*>sig_malloc(sizeof(int)*self._nonzero*2)
            for i from 0 <= i < self._length:
                v = data[i]
                if v != 0:
                    self._data[ind] = i
                    self._data[ind+1] = v
                    ind += 2
        else:
            raise TypeError("Error in ETuple(%s, %s, %s)" % (self, data, length))

    def __cinit__(self):
        self._data = <int*>0

    def __dealloc__(self):
        if self._data != <int*>0:
            sig_free(self._data)

    def __bool__(self):
        r"""
        Return whether ``self`` is nonzero.

        TESTS::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: bool(ETuple([1]))
            True
            sage: bool(ETuple([]))
            False
            sage: bool(ETuple([0, 0, 0]))
            False
        """
        return bool(self._nonzero)

    # methods to simulate tuple

    def __add__(ETuple self, ETuple other):
        """
        ``x.__add__(n) <==> x+n``.

        Concatenate two ETuples.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1, 1, 0]) + ETuple({int(1): int(2)}, int(3))
            (1, 1, 0, 0, 2, 0)
        """
        cdef size_t index = 0
        cdef ETuple result = <ETuple>ETuple.__new__(ETuple)
        result._length = self._length+other._length
        result._nonzero = self._nonzero+other._nonzero
        result._data = <int*>sig_malloc(sizeof(int)*result._nonzero*2)
        for index from 0 <= index < self._nonzero:
            result._data[2*index] = self._data[2*index]
            result._data[2*index+1] = self._data[2*index+1]
        for index from 0 <= index < other._nonzero:
            result._data[2*(index+self._nonzero)] = other._data[2*index]+self._length  # offset the second tuple (append to end!)
            result._data[2*(index+self._nonzero)+1] = other._data[2*index+1]
        return result

    def __mul__(ETuple self, factor):
        """
        ``x.__mul__(n) <==> x*n``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1, 2, 3])*2
            (1, 2, 3, 1, 2, 3)
        """
        cdef int _factor = factor
        cdef ETuple result = <ETuple>ETuple.__new__(ETuple)
        if factor <= 0:
            result._length = 0
            result._nonzero = 0
            return result
        cdef size_t index
        cdef size_t f
        result._length = self._length * factor
        result._nonzero = self._nonzero * factor
        result._data = <int*>sig_malloc(sizeof(int)*result._nonzero*2)
        for index from 0 <= index < self._nonzero:
            for f from 0 <= f < factor:
                result._data[2*(f*self._nonzero+index)] = self._data[2*index]+f*self._length
                result._data[2*(f*self._nonzero+index)+1] = self._data[2*index+1]
        return result

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: m = ETuple([1, 2, 0, 3])
            sage: m[2]
            0
            sage: m[1]
            2
            sage: e = ETuple([1, 2, 3])
            sage: e[1:]
            (2, 3)
            sage: e[:1]
            (1,)
        """
        cdef size_t ind
        if isinstance(i, slice):
            start, stop = i.start, i.stop
            if start is None:
                start = 0
            elif start < 0:
                start = start % self._length
            elif start > self._length:
                start = self._length

            if stop is None or stop > self._length:
                stop = self._length
            elif stop < 0:
                stop = stop % self._length

            # this is not particularly fast, but I doubt many people care
            # if you do, feel free to tweak!
            d = [self[ind] for ind from start <= ind < stop]
            return ETuple(d)
        else:
            return self.get_exp(i)

    cdef size_t get_position(self, size_t i, size_t start, size_t end) noexcept:
        r"""
        Return where to insert ``i`` in the data between ``start`` and ``end``.
        """
        if end <= start:
            return start
        cdef size_t left = start
        cdef size_t right = end - 1
        cdef size_t mid
        if self._data[2 * left] >= i:
            return left
        if self._data[2 * right] < i:
            return end
        if self._data[2 * right] == i:
            return right
        while right - left > 1:
            mid = (left + right) / 2
            if self._data[2 * mid] == i:
                return mid
            if self._data[2 * mid] > i:
                right = mid
            else:
                left = mid
        return right

    cdef int get_exp(self, size_t i) noexcept:
        """
        Return the exponent for the ``i``-th variable.
        """
        cdef size_t ind = self.get_position(i, 0, self._nonzero)
        if ind != self._nonzero and self._data[2 * ind] == i:
            return self._data[2 * ind + 1]
        return 0

    def __hash__(self):
        """
        x.__hash__() <==> hash(x)
        """
        cdef int i
        cdef int result = 0
        for i in range(self._nonzero):
            result += (1000003 * result) ^ self._data[2*i]
            result += (1000003 * result) ^ self._data[2*i+1]
        result = (1000003 * result) ^ self._length
        return result

    def __len__(self):
        """
        ``x.__len__() <==> len(x)``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2, 0, 3])
            sage: len(e)
            5
        """
        return self._length

    def __contains__(self, elem):
        """
        ``x.__contains__(n) <==> n in x``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple({int(1): int(2)}, int(3))
            sage: e
            (0, 2, 0)
            sage: 1 in e
            False
            sage: 2 in e
            True
        """
        if elem==0:
            return self._length > self._nonzero

        cdef size_t ind = 0
        for ind in range(self._nonzero):
            if elem == self._data[2 * ind + 1]:
                return True
        return False

    def __richcmp__(ETuple self, ETuple other, op):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1, 1, 0])<ETuple([1, 1, 0])
            False
            sage: ETuple([1, 1, 0])<ETuple([1, 0, 0])
            False
            sage: ETuple([1, 1, 0])<ETuple([1, 2, 0])
            True
            sage: ETuple([1, 1, 0])<ETuple([1, -1, 0])
            False
            sage: ETuple([0, -2, 0])<ETuple([1, -1, 0])
            True
            sage: ETuple([1, 1, 0])>ETuple([1, 1, 0])
            False
            sage: ETuple([1, 1, 0])>ETuple([1, 0, 0])
            True
            sage: ETuple([1, 1, 0])>ETuple([1, 2, 0])
            False
            sage: ETuple([1, 1, 0])>ETuple([1, -1, 0])
            True
            sage: ETuple([0, -2, 0])>ETuple([1, -1, 0])
            False
        """
        cdef size_t ind = 0
        if op == Py_EQ:  # ==
            if self._nonzero != other._nonzero:
                return False
            for ind in range(self._nonzero):
                if self._data[2 * ind] != other._data[2 * ind]:
                    return False
                if self._data[2 * ind + 1] != other._data[2 * ind + 1]:
                    return False
            return self._length == other._length

        if op == Py_LT:  # <
            while ind < self._nonzero and ind < other._nonzero:
                if self._data[2*ind] < other._data[2*ind]:
                    return self._data[2*ind+1] < 0
                if self._data[2*ind] > other._data[2*ind]:
                    return other._data[2*ind+1] > 0
                if self._data[2*ind] == other._data[2*ind] and self._data[2*ind+1] != other._data[2*ind+1]:
                    return self._data[2*ind+1] < other._data[2*ind+1]
                ind += 1
            if ind < self._nonzero and ind == other._nonzero:
                return self._data[2*ind+1] < 0
            if ind < other._nonzero and ind == self._nonzero:
                return other._data[2*ind+1] > 0
            return self._length < other._length

        if op == Py_GT:  # >
            while ind < self._nonzero and ind < other._nonzero:
                if self._data[2*ind] < other._data[2*ind]:
                    return self._data[2*ind+1] > 0
                if self._data[2*ind] > other._data[2*ind]:
                    return other._data[2*ind+1] < 0
                if self._data[2*ind] == other._data[2*ind] and self._data[2*ind+1] != other._data[2*ind+1]:
                    return self._data[2*ind+1] > other._data[2*ind+1]
                ind += 1
            if ind < self._nonzero and ind == other._nonzero:
                return self._data[2*ind+1] > 0
            if ind < other._nonzero and ind == self._nonzero:
                return other._data[2*ind+1] < 0
            return self._length < other._length

        # the rest of these are not particularly fast

        if op == Py_LE:  # <=
            return tuple(self) <= tuple(other)

        if op == Py_NE:  # !=
            return tuple(self) != tuple(other)

        if op == Py_GE:  # >=
            return tuple(self) >= tuple(other)

    def __iter__(self):
        """
        ``x.__iter__() <==> iter(x)``.

        TESTS::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple((4, 0, 0, 2, 0))
            sage: list(e)
            [4, 0, 0, 2, 0]

        Check that :issue:`28178` is fixed::

            sage: it = iter(e)
            sage: iter(it) is it
            True
        """
        cdef size_t i
        cdef size_t ind = 0

        for i in range(self._length):
            if ind >= self._nonzero:
                yield 0
            elif self._data[2*ind] == i:
                yield self._data[2*ind + 1]
                ind += 1
            else:
                yield 0

    def __str__(self):
        return repr(self)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple((0,))
            (0,)
            sage: ETuple((1,))
            (1,)
            sage: ETuple((0, 1, 2))
            (0, 1, 2)
        """
        if self._length == 1:
            if self._nonzero:
                return '(%d,)' % self._data[1]
            else:
                return '(0,)'
        else:
            return '(' + ', '.join(map(str, self)) + ')'

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 1, 0])
            sage: e == loads(dumps(e))
            True
        """
        cdef size_t ind
        d = {self._data[2 * ind]: self._data[2 * ind + 1] for ind in range(self._nonzero)}
        return ETuple, (d, int(self._length))

    # additional methods

    cdef int _unweighted_degree(self) except *:
        r"""
        Return the sum of entries.

        EXAMPLES::

             sage: from sage.rings.polynomial.polydict import ETuple
             sage: ETuple([1, 1, 0, 2, 0]).unweighted_degree()
             4
             sage: ETuple([-1, 1]).unweighted_degree()
             0
        """
        cdef int degree = 0
        cdef size_t i
        for i in range(self._nonzero):
            degree += self._data[2 * i + 1]
        return degree

    cpdef int unweighted_degree(self) except *:
        r"""
        Return the sum of entries.

        EXAMPLES::

             sage: from sage.rings.polynomial.polydict import ETuple
             sage: ETuple([1, 1, 0, 2, 0]).unweighted_degree()
             4
             sage: ETuple([-1, 1]).unweighted_degree()
             0
        """
        return self._unweighted_degree()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef int weighted_degree(self, tuple w) except *:
        r"""
        Return the weighted sum of entries.

        INPUT:

        - ``w`` -- tuple of nonnegative integers

        EXAMPLES::

             sage: from sage.rings.polynomial.polydict import ETuple
             sage: e = ETuple([1, 1, 0, 2, 0])
             sage: e.weighted_degree((1, 2, 3, 4, 5))
             11
             sage: ETuple([-1, 1]).weighted_degree((1, 2))
             1

             sage: ETuple([1, 0]).weighted_degree((1, 2, 3))
             Traceback (most recent call last):
             ...
             ValueError: w must be of the same length as the ETuple
        """
        if len(w) != self._length:
            raise ValueError('w must be of the same length as the ETuple')

        cdef size_t i
        cdef int deg = 0
        if len(w) != self._length:
            raise ValueError
        # NOTE: cython does optimize range(a) and range(a, b) but not range(a, b, c)
        for i in range(self._nonzero):
            deg += self._data[2 * i + 1] * <int> w[self._data[2 * i]]
        return deg

    cpdef int unweighted_quotient_degree(self, ETuple other) except *:
        """
        Return the degree of ``self`` divided by its gcd with ``other``.

        It amounts to counting the nonnegative entries of
        ``self.esub(other)``.
        """
        cdef size_t ind1 = 0    # both ind1 and ind2 will be increased in double steps.
        cdef size_t ind2 = 0
        cdef int exponent
        cdef int position
        cdef size_t selfnz = 2 * self._nonzero
        cdef size_t othernz = 2 * other._nonzero

        cdef int deg = 0
        while ind1 < selfnz:
            position = self._data[ind1]
            exponent = self._data[ind1 + 1]
            while ind2 < othernz and other._data[ind2] < position:
                ind2 += 2
            if ind2 == othernz:
                while ind1 < selfnz:
                    deg += self._data[ind1 + 1]
                    ind1 += 2
                return deg
            if other._data[ind2] > position:
                # other[position] = 0
                deg += exponent
            elif other._data[ind2 + 1] < exponent:
                # There is a positive difference that we have to insert
                deg += (exponent - other._data[ind2 + 1])
            ind1 += 2
        return deg

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef int weighted_quotient_degree(self, ETuple other, tuple w) except *:
        r"""
        Return the weighted degree of ``self`` divided by its gcd with ``other``.

        INPUT:

        - ``other`` -- an :class:`~sage.rings.polynomial.polydict.ETuple`
        - ``w`` -- tuple of nonnegative integers
        """
        if len(w) != self._length:
            raise ValueError('w must be of the same length as the ETuple')

        cdef size_t ind1 = 0    # both ind1 and ind2 will be increased in double steps.
        cdef size_t ind2 = 0
        cdef size_t exponent
        cdef int position
        cdef size_t selfnz = 2 * self._nonzero
        cdef size_t othernz = 2 * other._nonzero

        cdef size_t deg = 0
        assert len(w) == self._length
        while ind1 < selfnz:
            position = self._data[ind1]
            exponent = self._data[ind1+1]
            while ind2 < othernz and other._data[ind2] < position:
                ind2 += 2
            if ind2 == othernz:
                while ind1 < selfnz:
                    deg += <size_t>self._data[ind1+1] * <size_t> w[self._data[ind1]]
                    ind1 += 2
                return deg
            if other._data[ind2] > position:
                # other[position] = 0
                deg += exponent * <size_t>w[position]
            elif other._data[ind2+1] < exponent:
                # There is a positive difference that we have to insert
                deg += <size_t> (exponent - other._data[ind2+1]) * <size_t>w[position]
            ind1 += 2
        return deg

    cpdef ETuple eadd(self, ETuple other):
        """
        Return the vector addition of ``self`` with ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 1, 1])
            sage: e.eadd(f)
            (1, 1, 3)

        Verify that :issue:`6428` has been addressed::

            sage: # needs sage.libs.singular
            sage: R.<y, z> = Frac(QQ['x'])[]
            sage: type(y)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: y^(2^32)
            Traceback (most recent call last):
            ...
            OverflowError: exponent overflow (...)   # 64-bit
            OverflowError: Python int too large to convert to C unsigned long  # 32-bit
        """
        if self._length != other._length:
            raise ArithmeticError('ETuple of different lengths')

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef int s  # sum
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sig_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self, other, &ind1, &ind2, &index, &exp1, &exp2):
            s = exp1 + exp2
            # Check for overflow and underflow
            if (exp2 > 0 and s < exp1) or (exp2 < 0 and s > exp1):
                raise OverflowError("exponent overflow (%s)" % (int(exp1)+int(exp2)))
            if s != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = s
                result._nonzero += 1
        return result

    cpdef ETuple eadd_p(self, int other, size_t pos):
        """
        Add ``other`` to ``self`` at position ``pos``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: e.eadd_p(5, 1)
            (1, 5, 2)
            sage: e = ETuple([0]*7)
            sage: e.eadd_p(5, 4)
            (0, 0, 0, 0, 5, 0, 0)

            sage: ETuple([0,1]).eadd_p(1, 0) == ETuple([1,1])
            True

            sage: e = ETuple([0, 1, 0])
            sage: e.eadd_p(0, 0).nonzero_positions()
            [1]
            sage: e.eadd_p(0, 1).nonzero_positions()
            [1]
            sage: e.eadd_p(0, 2).nonzero_positions()
            [1]

        TESTS:

        Test segmentation faults occurring as described in :issue:`34000`::

            sage: ETuple([0, 1, 1]).eadd_p(1, 0)
            (1, 1, 1)
            sage: ETuple([0, 2, 4, 3]).eadd_p(5, 0)
            (5, 2, 4, 3)
            sage: ETuple([0, 2]).eadd_p(5, 0)
            (5, 2)
            sage: e = ETuple([0, 1, 0])
            sage: e.eadd_p(0, 0).nonzero_positions()
            [1]
            sage: e.eadd_p(0, 1).nonzero_positions()
            [1]
            sage: e.eadd_p(0, 2).nonzero_positions()
            [1]
            sage: e.eadd_p(-1, 1).nonzero_positions()
            []
        """
        cdef size_t sindex = 0
        cdef size_t rindex = 0
        cdef int new_value
        if pos >= self._length:
            raise ValueError("pos must be between 0 and %s" % self._length)

        cdef ETuple result = self._new()
        result._nonzero = self._nonzero
        if not other:
            # return a copy
            result._data = <int*> sig_malloc(sizeof(int) * self._nonzero * 2)
            memcpy(result._data, self._data, sizeof(int) * self._nonzero * 2)
            return result

        result._data = <int*> sig_malloc(sizeof(int) * (self._nonzero + 1) * 2)
        while sindex < self._nonzero and self._data[2 * sindex] < pos:
            result._data[2 * sindex] = self._data[2 * sindex]
            result._data[2 * sindex + 1] = self._data[2 * sindex + 1]
            sindex += 1

        if sindex < self._nonzero and self._data[2 * sindex] == pos:
            new_value = self._data[2 * sindex + 1] + other
            if new_value:
                result._data[2 * sindex] = pos
                result._data[2 * sindex + 1] = new_value
                sindex += 1
                rindex = sindex
            else:
                result._nonzero -= 1
                rindex = sindex
                sindex += 1
        else:
            result._data[2 * sindex] = pos
            result._data[2 * sindex + 1] = other
            result._nonzero += 1
            rindex = sindex + 1

        memcpy(result._data + 2 * rindex,
               self._data + 2 * sindex,
               sizeof(int) * 2 * (self._nonzero - sindex))

        return result

    cpdef ETuple eadd_scaled(self, ETuple other, int scalar):
        """
        Vector addition of ``self`` with ``scalar * other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 1, 1])
            sage: e.eadd_scaled(f, 3)
            (1, 3, 5)
        """
        if self._length != other._length:
            raise ArithmeticError('ETuple of different lengths')

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef int s  # sum
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sig_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self, other, &ind1, &ind2, &index, &exp1, &exp2):
            exp2 *= scalar
            s = exp1 + exp2
            # Check for overflow and underflow
            if (exp2 > 0 and s < exp1) or (exp2 < 0 and s > exp1):
                raise OverflowError("exponent overflow (%s)" % (int(exp1)+int(exp2)))
            if s != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = s
                result._nonzero += 1
        return result

    cpdef ETuple esub(self, ETuple other):
        """
        Vector subtraction of ``self`` with ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 1, 1])
            sage: e.esub(f)
            (1, -1, 1)
        """
        if self._length!=other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef int d  # difference
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sig_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self, other, &ind1, &ind2, &index, &exp1, &exp2):
            # Check for overflow and underflow
            d = exp1 - exp2
            if (exp2 > 0 and d > exp1) or (exp2 < 0 and d < exp1):
                raise OverflowError("Exponent overflow (%s)" % (int(exp1)-int(exp2)))
            if d != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = d
                result._nonzero += 1
        return result

    cpdef ETuple emul(self, int factor):
        """
        Scalar Vector multiplication of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: e.emul(2)
            (2, 0, 4)
        """
        cdef size_t ind
        cdef ETuple result = <ETuple>self._new()
        if factor == 0:
            result._nonzero = 0  # all zero, no nonzero entries!
            result._data = <int*>sig_malloc(sizeof(int) * result._nonzero * 2)
        else:
            result._nonzero = self._nonzero
            result._data = <int*>sig_malloc(sizeof(int) * result._nonzero * 2)
            for ind in range(self._nonzero):
                result._data[2 * ind] = self._data[2 * ind]
                result._data[2 * ind + 1] = self._data[2 * ind + 1] * factor
        return result

    cpdef ETuple emax(self, ETuple other):
        """
        Vector of maximum of components of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 1, 1])
            sage: e.emax(f)
            (1, 1, 2)
            sage: e = ETuple((1, 2, 3, 4))
            sage: f = ETuple((4, 0, 2, 1))
            sage: f.emax(e)
            (4, 2, 3, 4)
            sage: e = ETuple((1, -2, -2, 4))
            sage: f = ETuple((4, 0, 0, 0))
            sage: f.emax(e)
            (4, 0, 0, 4)
            sage: f.emax(e).nonzero_positions()
            [0, 3]
        """
        if self._length!=other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sig_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self, other, &ind1, &ind2, &index, &exp1, &exp2):
            if exp1 >= exp2 and exp1 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp1
                result._nonzero += 1
            elif exp2 >= exp1 and exp2 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp2
                result._nonzero += 1
        return result

    cpdef ETuple emin(self, ETuple other):
        """
        Vector of minimum of components of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 1, 1])
            sage: e.emin(f)
            (0, 0, 1)
            sage: e = ETuple([1, 0, -1])
            sage: f = ETuple([0, -2, 1])
            sage: e.emin(f)
            (0, -2, -1)
        """
        if self._length != other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sig_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self, other, &ind1, &ind2, &index, &exp1, &exp2):
            if exp1 <= exp2 and exp1 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp1
                result._nonzero += 1
            elif exp2 <= exp1 and exp2 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp2
                result._nonzero += 1
        return result

    cpdef int dotprod(self, ETuple other) except *:
        """
        Return the dot product of this tuple by ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 1, 1])
            sage: e.dotprod(f)
            2
            sage: e = ETuple([1, 1, -1])
            sage: f = ETuple([0, -2, 1])
            sage: e.dotprod(f)
            -3
        """
        if self._length != other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef int result = 0
        while dual_etuple_iter(self, other, &ind1, &ind2, &index, &exp1, &exp2):
            result += exp1 * exp2
        return result

    cpdef ETuple escalar_div(self, int n):
        r"""
        Divide each exponent by ``n``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1, 0, 2]).escalar_div(2)
            (0, 0, 1)
            sage: ETuple([0, 3, 12]).escalar_div(3)
            (0, 1, 4)

            sage: ETuple([1, 5, 2]).escalar_div(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        TESTS:

        Checking that memory allocation works fine::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: t = ETuple(list(range(2048)))
            sage: for n in range(1, 9):
            ....:     t = t.escalar_div(n)
            sage: assert t.is_constant()
        """
        if not n:
            raise ZeroDivisionError
        cdef size_t i
        cdef ETuple result = self._new()
        result._data = <int*> sig_malloc(sizeof(int) * 2 * self._nonzero)
        result._nonzero = 0
        # NOTE: cython does optimize range(a) and range(a, b) but not range(a, b, c)
        for i in range(self._nonzero):
            result._data[2 * result._nonzero + 1] = self._data[2 * i + 1] / n
            if result._data[2 * result._nonzero + 1]:
                result._data[2 * result._nonzero] = self._data[2 * i]
                result._nonzero += 1
        return result

    cpdef ETuple divide_by_gcd(self, ETuple other):
        """
        Return ``self / gcd(self, other)``.

        The entries of the result are the maximum of 0 and the
        difference of the corresponding entries of ``self`` and ``other``.
        """
        if self._length != other._length:
            raise ArithmeticError('ETuple of different lengths')
        cdef size_t ind1 = 0    # both ind1 and ind2 will be increased in 2-steps.
        cdef size_t ind2 = 0
        cdef int exponent
        cdef int position
        cdef size_t selfnz = 2 * self._nonzero
        cdef size_t othernz = 2 * other._nonzero
        cdef ETuple result = <ETuple> self._new()
        result._nonzero = 0
        result._data = <int*> sig_malloc(sizeof(int)*self._nonzero*2)
        while ind1 < selfnz:
            position = self._data[ind1]
            exponent = self._data[ind1+1]
            while ind2 < othernz and other._data[ind2] < position:
                ind2 += 2
            if ind2 == othernz:
                while ind1 < selfnz:
                    result._data[2*result._nonzero] = self._data[ind1]
                    result._data[2*result._nonzero+1] = self._data[ind1+1]
                    result._nonzero += 1
                    ind1 += 2
                return result
            if other._data[ind2] > position:
                # other[position] == 0
                result._data[2*result._nonzero] = position
                result._data[2*result._nonzero+1] = exponent
                result._nonzero += 1
            elif other._data[ind2+1] < exponent:
                # There is a positive difference that we have to insert
                result._data[2*result._nonzero] = position
                result._data[2*result._nonzero+1] = exponent - other._data[ind2+1]
                result._nonzero += 1
            ind1 += 2
        return result

    cpdef ETuple divide_by_var(self, size_t pos):
        """
        Return division of ``self`` by the variable with index ``pos``.

        If ``self[pos] == 0`` then a :exc:`ArithmeticError` is raised. Otherwise,
        an :class:`~sage.rings.polynomial.polydict.ETuple` is returned that is
        zero in position ``pos`` and coincides with ``self`` in the other
        positions.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 2, 0, 1])
            sage: e.divide_by_var(0)
            (0, 2, 0, 1)
            sage: e.divide_by_var(1)
            (1, 1, 0, 1)
            sage: e.divide_by_var(3)
            (1, 2, 0, 0)
            sage: e.divide_by_var(2)
            Traceback (most recent call last):
            ...
            ArithmeticError: not divisible by this variable
        """
        cdef int exp1
        cdef ETuple result
        cdef size_t ind

        ind = self.get_position(pos, 0, self._nonzero)
        if ind == self._nonzero or self._data[2 * ind] != pos:
            raise ArithmeticError('not divisible by this variable')
        result = <ETuple> self._new()
        result._data = <int*> sig_malloc(sizeof(int) * 2 * self._nonzero)
        exp1 = self._data[2 * ind + 1]
        if exp1 > 1:
            # division doesn't change the number of nonzero positions
            result._nonzero = self._nonzero
            memcpy(result._data, self._data, sizeof(int) * 2 * self._nonzero)
            result._data[2 * ind + 1] = exp1 - 1
        else:
            # var(pos) disappears from self
            result._nonzero = self._nonzero - 1
            memcpy(result._data, self._data, sizeof(int) * 2 * ind)
            if ind + 1 < self._nonzero:
                memcpy(result._data + 2 * ind, self._data + 2 * (ind + 1), sizeof(int) * 2 * (self._nonzero - ind - 1))
        return result

    cpdef bint divides(self, ETuple other) except *:
        """
        Return whether ``self`` divides ``other``, i.e., no entry of ``self``
        exceeds that of ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1, 1, 0, 1, 0]).divides(ETuple([2, 2, 2, 2, 2]))
            True
            sage: ETuple([0, 3, 0, 1, 0]).divides(ETuple([2, 2, 2, 2, 2]))
            False
            sage: ETuple([0, 3, 0, 1, 0]).divides(ETuple([0, 3, 2, 2, 2]))
            True
            sage: ETuple([0, 0, 0, 0, 0]).divides(ETuple([2, 2, 2, 2, 2]))
            True

            sage: ETuple({104: 18, 256: 25, 314:78}, length=400r).divides(ETuple({104: 19, 105: 20, 106: 21}, length=400r))
            False
            sage: ETuple({104: 18, 256: 25, 314:78}, length=400r).divides(ETuple({104: 19, 105: 20, 106: 21, 255: 2, 256: 25, 312: 5, 314: 79, 315: 28}, length=400r))
            True
        """
        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef int pos1

        if self._length != other._length:
            raise ArithmeticError('ETuple of different length')

        while ind1 < self._nonzero:
            if self._nonzero - ind1 > other._nonzero - ind2:
                return False
            pos1 = self._data[2 * ind1]
            ind2 = other.get_position(pos1, ind2, other._nonzero)
            if ind2 == other._nonzero or other._data[2 * ind2] != pos1 or other._data[2 * ind2 + 1] < self._data[2 * ind1 + 1]:
                return False
            ind1 += 1
            ind2 += 1

        return True

    cpdef bint is_constant(self) noexcept:
        """
        Return if all exponents are zero in the tuple.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: e.is_constant()
            False
            sage: e = ETuple([0, 0])
            sage: e.is_constant()
            True
        """
        return self._nonzero == 0

    cpdef bint is_multiple_of(self, int n) except *:
        r"""
        Test whether each entry is a multiple of ``n``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple

            sage: ETuple([0, 0]).is_multiple_of(3)
            True
            sage: ETuple([0, 3, 12, 0, 6]).is_multiple_of(3)
            True
            sage: ETuple([0, 0, 2]).is_multiple_of(3)
            False
        """
        if not n:
            raise ValueError('n should not be zero')
        cdef size_t i
        for i in range(self._nonzero):
            if self._data[2 * i + 1] % n:
                return False
        return True

    cpdef list nonzero_positions(self, bint sort=False):
        """
        Return the positions of nonzero exponents in the tuple.

        INPUT:

        - ``sort`` -- boolean (default: ``False``); if ``True`` a sorted list is
          returned; if ``False`` an unsorted list is returned

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: e.nonzero_positions()
            [0, 2]
        """
        cdef size_t ind
        return [self._data[2 * ind] for ind in range(self._nonzero)]

    cpdef common_nonzero_positions(self, ETuple other, bint sort=False):
        """
        Return an optionally sorted list of nonzero positions either
        in ``self`` or other, i.e. the only positions that need to be
        considered for any vector operation.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2])
            sage: f = ETuple([0, 0, 1])
            sage: e.common_nonzero_positions(f)
            {0, 2}
            sage: e.common_nonzero_positions(f, sort=True)
            [0, 2]
        """
        # TODO:  we should probably make a fast version of this!
        res = set(self.nonzero_positions()).union(other.nonzero_positions())
        if sort:
            return sorted(res)
        else:
            return res

    cpdef list nonzero_values(self, bint sort=True):
        """
        Return the nonzero values of the tuple.

        INPUT:

        - ``sort`` -- boolean (default: ``True``); if ``True`` the values are
          sorted by their indices. Otherwise the values are returned unsorted.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([2, 0, 1])
            sage: e.nonzero_values()
            [2, 1]
            sage: f = ETuple([0, -1, 1])
            sage: f.nonzero_values(sort=True)
            [-1, 1]
        """
        cdef size_t ind
        return [self._data[2 * ind + 1] for ind in range(self._nonzero)]

    cpdef ETuple reversed(self):
        """
        Return the reversed ETuple of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 2, 3])
            sage: e.reversed()
            (3, 2, 1)
        """
        cdef size_t ind
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = self._nonzero
        result._data = <int*>sig_malloc(sizeof(int) * result._nonzero * 2)
        for ind in range(self._nonzero):
            result._data[2 * (result._nonzero - ind - 1)] = self._length - self._data[2 * ind] - 1
            result._data[2 * (result._nonzero - ind - 1) + 1] = self._data[2 * ind + 1]
        return result

    def sparse_iter(self):
        """
        Iterator over the elements of ``self`` where the elements are returned
        as ``(i, e)`` where ``i`` is the position of ``e`` in the tuple.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1, 0, 2, 0, 3])
            sage: list(e.sparse_iter())
            [(0, 1), (2, 2), (4, 3)]
        """
        cdef size_t ind
        for ind in range(self._nonzero):
            yield (self._data[2 * ind], self._data[2 * ind + 1])

    def combine_to_positives(self, ETuple other):
        """
        Given a pair of ETuples (self, other), returns a triple of
        ETuples (a, b, c) so that self = a + b, other = a + c and b and c
        have all positive entries.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([-2, 1, -5, 3, 1, 0])
            sage: f = ETuple([1, -3, -3, 4, 0, 2])
            sage: e.combine_to_positives(f)
            ((-2, -3, -5, 3, 0, 0), (0, 4, 0, 0, 1, 0), (3, 0, 2, 1, 0, 2))
        """
        m = self.emin(other)
        return m, self.esub(m), other.esub(m)


def make_PolyDict(data):
    r"""
    Ensure support for pickled data from older sage versions.
    """
    return PolyDict(data)


def make_ETuple(data, length):
    r"""
    Ensure support for pickled data from older sage versions.
    """
    return ETuple(data, length)
