# distutils: libraries = flint
# distutils: depends = flint/fmpq_mpoly.h flint/fmpq_mpoly_factor.h
r"""
Multivariate polynomials over `\QQ`, implemented using FLINT

EXAMPLES::

    sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
    sage: f = 3*x^2*y/2 - z + 1; f
    3/2*x^2*y - z + 1
    sage: f * f
    9/4*x^4*y^2 - 3*x^2*y*z + 3*x^2*y + z^2 - 2*z + 1

.. automethod:: MPolynomial_rational_flint._add_
.. automethod:: MPolynomial_rational_flint._sub_
.. automethod:: MPolynomial_rational_flint._mul_
.. automethod:: MPolynomial_rational_flint._neg_
"""

#*****************************************************************************
#       Copyright (C) 2025 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re

from libc.stdlib cimport free

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport sig_malloc, sig_free

from cpython.object cimport Py_EQ, Py_NE

from sage.libs.flint.fmpq cimport (
    fmpq_init, fmpq_clear, fmpq_set_mpq, fmpq_get_mpq, fmpq_is_zero)
from sage.libs.flint.fmpq_mpoly cimport *
from sage.libs.flint.fmpq_mpoly_factor cimport (
    fmpq_mpoly_factor_init, fmpq_mpoly_factor_clear,
    fmpq_mpoly_factor_length, fmpq_mpoly_factor_get_constant_fmpq,
    fmpq_mpoly_factor_get_base, fmpq_mpoly_factor_get_exp_si,
    fmpq_mpoly_factor)
from sage.libs.flint.types cimport (
    fmpq_t, fmpq_mpoly_t, fmpq_mpoly_ctx_t, fmpq_mpoly_struct,
    fmpq_mpoly_factor_t,
    ordering_t, ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX, ulong, slong)

from sage.cpython.string cimport str_to_bytes, char_to_str

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational cimport Rational
from sage.rings.rational_field import QQ
from sage.structure.element cimport Element
from sage.structure.factorization import Factorization
from sage.structure.richcmp cimport rich_to_bool

from sage.rings.polynomial.multi_polynomial cimport MPolynomial_flint as MPolynomial_flint_base
from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial_rational_flint cimport (
    MPolynomialRing_rational_flint,
    MPolynomial_rational_flint)
from sage.rings.polynomial.polydict cimport ETuple


cdef ordering_t _term_order_to_flint(order) except? ORD_LEX:
    name = order.name() if hasattr(order, 'name') else str(order)
    if name == 'lex':
        return ORD_LEX
    elif name == 'deglex':
        return ORD_DEGLEX
    elif name == 'degrevlex':
        return ORD_DEGREVLEX
    raise NotImplementedError(
        "term order '{}' is not supported by FLINT fmpq_mpoly".format(name))


def _unpickle_MPolynomialRing_rational_flint(n, names, order):
    r"""
    Helper for unpickling :class:`MPolynomialRing_rational_flint`.

    TESTS::

        sage: R = PolynomialRing(QQ, 'x,y', implementation='flint')
        sage: loads(dumps(R)) is R
        True
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    return PolynomialRing(QQ, n, names=names, order=order, implementation='flint')


def _unpickle_MPolynomial_rational_flint(parent, coeffs):
    r"""
    Helper for unpickling :class:`MPolynomial_rational_flint`.

    TESTS::

        sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
        sage: f = 3*x^2*y/2 - x + 1
        sage: loads(dumps(f)) == f
        True
    """
    return parent(coeffs)


cdef class MPolynomialRing_rational_flint(MPolynomialRing_base):
    r"""
    Multivariate polynomial ring over `\QQ`, implemented via FLINT.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
        sage: R
        Multivariate Polynomial Ring in x, y over Rational Field (using FLINT)
    """

    def __cinit__(self):
        # nvars = -1 signals that fmpq_mpoly_ctx_init has not yet been called.
        self._ctx[0].zctx[0].minfo[0].nvars = -1

    def __init__(self, base_ring, n, names, order='degrevlex'):
        r"""
        Construct a multivariate polynomial ring over `\QQ` using FLINT.

        INPUT:

        - ``base_ring`` -- must be `\QQ`
        - ``n`` -- number of variables (positive integer)
        - ``names`` -- variable names
        - ``order`` -- term order (``'lex'``, ``'deglex'``, or ``'degrevlex'``)

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_rational_flint import MPolynomialRing_rational_flint
            sage: R = MPolynomialRing_rational_flint(QQ, 3, ('x','y','z'), 'degrevlex')
            sage: R
            Multivariate Polynomial Ring in x, y, z over Rational Field (using FLINT)
        """
        if base_ring is not QQ:
            raise TypeError("base ring must be QQ")
        MPolynomialRing_base.__init__(self, base_ring, n, names, order)
        cdef ordering_t ord = _term_order_to_flint(self._term_order)
        fmpq_mpoly_ctx_init(self._ctx, n, ord)

    def __dealloc__(self):
        if self._ctx[0].zctx[0].minfo[0].nvars != -1:
            fmpq_mpoly_ctx_clear(self._ctx)

    Element = MPolynomial_rational_flint

    def __hash__(self):
        r"""
        Return a hash of this ring.

        EXAMPLES::

            sage: R = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: hash(R) == hash(R)
            True
        """
        from sage.structure.category_object import CategoryObject
        return CategoryObject.__hash__(self)

    def __reduce__(self):
        r"""
        Return data for pickling this ring.

        TESTS::

            sage: R = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: loads(dumps(R)) is R
            True
        """
        return (_unpickle_MPolynomialRing_rational_flint,
                (self._ngens, self.variable_names(), self._term_order))

    def _repr_(self):
        r"""
        Return a string representation of this ring.

        EXAMPLES::

            sage: PolynomialRing(QQ, 'x,y', implementation='flint')
            Multivariate Polynomial Ring in x, y over Rational Field (using FLINT)
        """
        return "Multivariate Polynomial Ring in {} over Rational Field (using FLINT)".format(
            ", ".join(self.variable_names()))

    def gen(self, int n=0):
        r"""
        Return the ``n``-th generator of this ring.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: R.gen(0)
            x
            sage: R.gen(2)
            z
        """
        if n < 0 or n >= self._ngens:
            raise ValueError("generator index out of range")
        cdef MPolynomial_rational_flint g = self._new_element()
        sig_on()
        fmpq_mpoly_gen(g._poly, n, self._ctx)
        sig_off()
        return g

    cdef MPolynomial_rational_flint _new_element(self):
        cdef MPolynomial_rational_flint f = \
            MPolynomial_rational_flint.__new__(MPolynomial_rational_flint)
        f._parent = self
        fmpq_mpoly_init(f._poly, self._ctx)
        return f

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into an element of this ring.

        Accepted inputs include integers, rationals, strings, dictionaries
        mapping exponent tuples to coefficients, and multivariate
        polynomials over a compatible ring.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: R(3)
            3
            sage: R(2/3)
            2/3
            sage: R("x^2*y - 2*x + 1/2")
            x^2*y - 2*x + 1/2
            sage: R({(2,1): 5/3, (0,0): -1})
            5/3*x^2*y - 1
            sage: R(x + y)
            x + y
        """
        cdef MPolynomial_rational_flint f
        cdef Rational c
        cdef ETuple e
        cdef fmpq_t coeff
        cdef ulong *exp
        cdef slong n = self._ngens
        cdef int i
        cdef bytes bstr
        cdef const char **cnames

        if isinstance(x, MPolynomial_rational_flint) and x.parent() is self:
            f = self._new_element()
            sig_on()
            fmpq_mpoly_set(f._poly, (<MPolynomial_rational_flint>x)._poly, self._ctx)
            sig_off()
            return f

        if isinstance(x, (int, Integer)):
            x = Rational(x)
        if isinstance(x, Rational):
            f = self._new_element()
            fmpq_init(coeff)
            fmpq_set_mpq(coeff, (<Rational>x).value)
            sig_on()
            fmpq_mpoly_set_fmpq(f._poly, coeff, self._ctx)
            sig_off()
            fmpq_clear(coeff)
            return f

        if isinstance(x, str):
            f = self._new_element()
            bnames = [str_to_bytes(v) for v in self.variable_names()]
            cnames = <const char **>sig_malloc(n * sizeof(char *))
            if cnames == NULL:
                raise MemoryError
            for i in range(n):
                cnames[i] = bnames[i]
            bstr = str_to_bytes(x)
            sig_on()
            ok = fmpq_mpoly_set_str_pretty(f._poly, bstr, cnames, self._ctx)
            sig_off()
            sig_free(cnames)
            if ok != 0:
                raise ValueError("could not parse '{}' as a polynomial".format(x))
            return f

        if isinstance(x, dict):
            f = self._new_element()
            exp = <ulong *>sig_malloc(n * sizeof(ulong))
            if exp == NULL:
                raise MemoryError
            fmpq_init(coeff)
            try:
                for key, val in x.items():
                    if not isinstance(val, Rational):
                        val = QQ(val)
                    c = <Rational>val
                    for i in range(n):
                        exp[i] = key[i]
                    fmpq_set_mpq(coeff, c.value)
                    sig_on()
                    fmpq_mpoly_set_coeff_fmpq_ui(f._poly, coeff, exp, self._ctx)
                    sig_off()
            finally:
                fmpq_clear(coeff)
                sig_free(exp)
            return f

        # try via _mpoly_dict_recursive for MPolynomial types
        from sage.rings.polynomial.multi_polynomial import MPolynomial
        if isinstance(x, MPolynomial):
            return self._element_constructor_(
                x._mpoly_dict_recursive(self.variable_names(), QQ))

        # last resort: coerce to QQ
        return self._element_constructor_(QQ(x))

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into this ring.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: R.has_coerce_map_from(ZZ)
            True
            sage: R.has_coerce_map_from(QQ)
            True
        """
        if R is QQ or R is ZZ:
            return True
        return MPolynomialRing_base._coerce_map_from_(self, R)


cdef class MPolynomial_rational_flint(MPolynomial_flint_base):
    r"""
    A multivariate polynomial over `\QQ`, implemented via FLINT.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
        sage: f = 3*x^2/2 - y + 1; f
        3/2*x^2 - y + 1
        sage: type(f)
        <class 'sage.rings.polynomial.multi_polynomial_rational_flint.MPolynomial_rational_flint'>

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _mul_
    .. automethod:: _neg_
    """

    def __cinit__(self):
        # Python zero-initialises the memory; fmpq_mpoly_clear on a struct
        # with alloc=0 is a no-op, so __dealloc__ is safe even if
        # fmpq_mpoly_init was never called.
        pass

    def __dealloc__(self):
        cdef MPolynomialRing_rational_flint R = self._parent
        if R is not None:
            fmpq_mpoly_clear(self._poly, R._ctx)

    cdef MPolynomial_rational_flint _new(self):
        return (<MPolynomialRing_rational_flint>self._parent)._new_element()

    cpdef _new_constant_poly(self, scalar, parent):
        r"""
        Return a new constant polynomial with value ``scalar`` in ``parent``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: x._new_constant_poly(3/4, R)
            3/4
        """
        cdef MPolynomialRing_rational_flint R = parent
        cdef MPolynomial_rational_flint f = R._new_element()
        cdef fmpq_t coeff
        if not isinstance(scalar, Rational):
            scalar = QQ(scalar)
        fmpq_init(coeff)
        fmpq_set_mpq(coeff, (<Rational>scalar).value)
        fmpq_mpoly_set_fmpq(f._poly, coeff, R._ctx)
        fmpq_clear(coeff)
        return f

    def __copy__(self):
        r"""
        Return a copy of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2/2 - y + 1
            sage: g = copy(f)
            sage: f == g
            True
            sage: f is g
            False
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint f = R._new_element()
        sig_on()
        fmpq_mpoly_set(f._poly, self._poly, R._ctx)
        sig_off()
        return f

    def __deepcopy__(self, memo=None):
        r"""
        Return a deep copy of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2 + y
            sage: g = deepcopy(f)
            sage: f == g
            True
            sage: f is g
            False
        """
        cpy = self.__copy__()
        if memo is not None:
            memo[id(self)] = cpy
        return cpy

    def __reduce__(self):
        r"""
        Return data for pickling this polynomial.

        TESTS::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: f = 7*x^3*y/2 - 2*y*z + 1
            sage: loads(dumps(f)) == f
            True
        """
        return (_unpickle_MPolynomial_rational_flint,
                (self._parent, self.monomial_coefficients()))

    def _repr_(self):
        r"""
        Return a string representation of this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: 3*x^2*y/2 - z + 1
            3/2*x^2*y - z + 1
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef int n = R._ngens
        cdef char *raw
        cdef const char **cnames = <const char **>sig_malloc(n * sizeof(char *))
        if cnames == NULL:
            raise MemoryError
        bnames = [str_to_bytes(v) for v in R.variable_names()]
        for i in range(n):
            cnames[i] = bnames[i]
        sig_on()
        raw = fmpq_mpoly_get_str_pretty(self._poly, cnames, R._ctx)
        sig_off()
        sig_free(cnames)
        result = char_to_str(raw)
        free(raw)
        # FLINT omits spaces around binary + and -; add them for readability
        result = re.sub(r'([a-zA-Z0-9])\s*([+-])', r'\1 \2 ', result)
        return result

    def __hash__(self):
        r"""
        Return a hash of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: hash(x) == hash(x)
            True
            sage: hash(x) == hash(y)
            False
        """
        return self._hash_c()

    def __bool__(self):
        r"""
        Return ``True`` if this polynomial is nonzero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: bool(x)
            True
            sage: bool(R(0))
            False
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return not fmpq_mpoly_is_zero(self._poly, R._ctx)

    # ===== Predicates =====

    def is_zero(self):
        r"""
        Return ``True`` if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: R(0).is_zero()
            True
            sage: x.is_zero()
            False
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return bool(fmpq_mpoly_is_zero(self._poly, R._ctx))

    def is_one(self):
        r"""
        Return ``True`` if this polynomial equals 1.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: R(1).is_one()
            True
            sage: x.is_one()
            False
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return bool(fmpq_mpoly_is_one(self._poly, R._ctx))

    def is_constant(self):
        r"""
        Return ``True`` if this polynomial is a constant.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: R(0).is_constant()
            True
            sage: R(5/3).is_constant()
            True
            sage: x.is_constant()
            False
            sage: (x + 1).is_constant()
            False

        Constant polynomials can be converted to `\QQ`::

            sage: QQ(R(5/3))
            5/3
            sage: QQ(x)
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return bool(fmpq_mpoly_is_fmpq(self._poly, R._ctx))

    def is_term(self):
        r"""
        Return ``True`` if this polynomial is a single nonzero term.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: x.is_term()
            True
            sage: (3*x^2*y/2).is_term()
            True
            sage: (x + y).is_term()
            False
            sage: R(0).is_term()
            False
            sage: R(5/3).is_term()
            True
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return fmpq_mpoly_length(self._poly, R._ctx) == 1

    def is_monomial(self):
        r"""
        Return ``True`` if this polynomial is a monomial.

        A *monomial* has coefficient 1.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: x.is_monomial()
            True
            sage: (x^2*y).is_monomial()
            True
            sage: (3*x).is_monomial()
            False
            sage: R(1).is_monomial()
            True
            sage: R(0).is_monomial()
            False
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef fmpq_t c
        if fmpq_mpoly_length(self._poly, R._ctx) != 1:
            return False
        fmpq_init(c)
        fmpq_mpoly_get_term_coeff_fmpq(c, self._poly, 0, R._ctx)
        cdef Rational rc = Rational.__new__(Rational)
        fmpq_get_mpq(rc.value, c)
        fmpq_clear(c)
        return rc.is_one()

    def is_homogeneous(self):
        r"""
        Return ``True`` if this polynomial is homogeneous.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: (x^2 + y*z).is_homogeneous()
            True
            sage: (x^2 + y).is_homogeneous()
            False
            sage: R(0).is_homogeneous()
            True
            sage: R(5/3).is_homogeneous()
            True
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong length = fmpq_mpoly_length(self._poly, R._ctx)
        cdef slong i, j
        cdef ulong *exp
        if length <= 1:
            return True
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        try:
            fmpq_mpoly_get_term_exp_ui(exp, self._poly, 0, R._ctx)
            ref_deg = 0
            for j in range(n):
                ref_deg += exp[j]
            for i in range(1, length):
                fmpq_mpoly_get_term_exp_ui(exp, self._poly, i, R._ctx)
                deg = 0
                for j in range(n):
                    deg += exp[j]
                if deg != ref_deg:
                    return False
        finally:
            sig_free(exp)
        return True

    def is_square(self):
        r"""
        Return ``True`` if this polynomial is a perfect square.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: (x^2 + 2*x*y + y^2).is_square()
            True
            sage: x.is_square()
            False
            sage: R(0).is_square()
            True
            sage: R(4).is_square()
            True
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return bool(fmpq_mpoly_is_square(self._poly, R._ctx))

    # ===== Information =====

    cpdef long number_of_terms(self) noexcept:
        r"""
        Return the number of nonzero terms.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2 + 3*x*y - y + 2
            sage: f.number_of_terms()
            4
            sage: R(0).number_of_terms()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        return fmpq_mpoly_length(self._poly, R._ctx)

    def degree(self, x=None):
        r"""
        Return the degree of this polynomial in the variable ``x``, or the
        total degree if ``x`` is ``None``.

        INPUT:

        - ``x`` -- (optional) a generator of the parent ring

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: f = x^2*y + y^3*z^4
            sage: f.degree()
            7
            sage: f.degree(x)
            2
            sage: f.degree(y)
            3
            sage: f.degree(z)
            4
            sage: R(0).degree()
            -1
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        if fmpq_mpoly_is_zero(self._poly, R._ctx):
            return -1
        if x is None:
            return self.total_degree()
        cdef slong var_idx
        names = R.variable_names()
        var_name = str(x)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(x, R))
        var_idx = names.index(var_name)
        return Integer(fmpq_mpoly_degree_si(self._poly, var_idx, R._ctx))

    def total_degree(self):
        r"""
        Return the total degree of this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: (x^2*y + y^3*z^4 + 1).total_degree()
            7
            sage: R(5).total_degree()
            0
            sage: R(0).total_degree()
            -1
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        if fmpq_mpoly_is_zero(self._poly, R._ctx):
            return -1
        return Integer(fmpq_mpoly_total_degree_si(self._poly, R._ctx))

    def degrees(self):
        r"""
        Return a tuple containing the degree of this polynomial in each
        variable.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: f = x^2*y + y^3*z^4
            sage: f.degrees()
            (2, 3, 4)
            sage: R(5).degrees()
            (0, 0, 0)
            sage: R(0).degrees()
            (0, 0, 0)
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef slong *degs = <slong *>sig_malloc(n * sizeof(slong))
        if degs == NULL:
            raise MemoryError
        try:
            fmpq_mpoly_degrees_si(degs, self._poly, R._ctx)
            result = tuple(Integer(degs[i]) if degs[i] >= 0 else Integer(0)
                           for i in range(n))
        finally:
            sig_free(degs)
        return result

    def variables(self):
        r"""
        Return a tuple of the variables (generators) actually occurring in
        this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: (x^2*z + 1).variables()
            (x, z)
            sage: R(5).variables()
            ()
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef int *used = <int *>sig_malloc(n * sizeof(int))
        if used == NULL:
            raise MemoryError
        try:
            fmpq_mpoly_used_vars(used, self._poly, R._ctx)
            result = tuple(R.gen(i) for i in range(n) if used[i])
        finally:
            sig_free(used)
        return result

    def variable(self, int i=0):
        r"""
        Return the ``i``-th variable occurring in this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: (x^2*z + 1).variable(0)
            x
            sage: (x^2*z + 1).variable(1)
            z
        """
        return self.variables()[i]

    def nvariables(self):
        r"""
        Return the number of variables actually occurring in this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: (x^2*z + 1).nvariables()
            2
            sage: R(5).nvariables()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef int *used = <int *>sig_malloc(n * sizeof(int))
        if used == NULL:
            raise MemoryError
        try:
            fmpq_mpoly_used_vars(used, self._poly, R._ctx)
            result = sum(1 for i in range(n) if used[i])
        finally:
            sig_free(used)
        return result

    # ===== Coefficients =====

    def monomial_coefficients(self):
        r"""
        Return a dictionary of ``{exponent: coefficient}`` pairs.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 2*x^3/3 - x*y + 5/7
            sage: sorted(f.monomial_coefficients().items())
            [((0, 0), 5/7), ((1, 1), -1), ((3, 0), 2/3)]
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong length = fmpq_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef fmpq_t coeff
        cdef ulong *exp
        cdef Rational c
        fmpq_init(coeff)
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            fmpq_clear(coeff)
            raise MemoryError
        result = {}
        try:
            for i in range(length):
                fmpq_mpoly_get_term_coeff_fmpq(coeff, self._poly, i, R._ctx)
                fmpq_mpoly_get_term_exp_ui(exp, self._poly, i, R._ctx)
                c = Rational.__new__(Rational)
                fmpq_get_mpq(c.value, coeff)
                result[ETuple([exp[j] for j in range(n)])] = c
        finally:
            fmpq_clear(coeff)
            sig_free(exp)
        return result

    def coefficients(self):
        r"""
        Return the list of coefficients of this polynomial, in the order
        induced by the term order.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y/2 - x + 5/7
            sage: f.coefficients()
            [3/2, -1, 5/7]
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong length = fmpq_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef fmpq_t coeff
        cdef Rational c
        fmpq_init(coeff)
        result = []
        try:
            for i in range(length):
                fmpq_mpoly_get_term_coeff_fmpq(coeff, self._poly, i, R._ctx)
                c = Rational.__new__(Rational)
                fmpq_get_mpq(c.value, coeff)
                result.append(c)
        finally:
            fmpq_clear(coeff)
        return result

    def exponents(self, as_ETuples=True):
        r"""
        Return the list of exponent tuples of the nonzero terms of this
        polynomial, in the order induced by the term order.

        INPUT:

        - ``as_ETuples`` -- (default: ``True``) if ``True``, return a list of
          :class:`ETuple`; otherwise, return a list of plain tuples

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y/2 - x + 5/7
            sage: f.exponents()
            [(2, 1), (1, 0), (0, 0)]
            sage: f.exponents(as_ETuples=False)
            [(2, 1), (1, 0), (0, 0)]
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong length = fmpq_mpoly_length(self._poly, R._ctx)
        cdef slong i, j
        cdef ulong *exp
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        result = []
        try:
            for i in range(length):
                fmpq_mpoly_get_term_exp_ui(exp, self._poly, i, R._ctx)
                t = tuple(int(exp[j]) for j in range(n))
                if as_ETuples:
                    result.append(ETuple(t))
                else:
                    result.append(t)
        finally:
            sig_free(exp)
        return result

    def monomials(self):
        r"""
        Return the list of monomials of this polynomial, in the order
        induced by the term order.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y/2 - x + 5/7
            sage: f.monomials()
            [x^2*y, x, 1]
            sage: R(0).monomials()
            []
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong length = fmpq_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef MPolynomial_rational_flint m
        result = []
        for i in range(length):
            m = R._new_element()
            sig_on()
            fmpq_mpoly_get_term_monomial(m._poly, self._poly, i, R._ctx)
            sig_off()
            result.append(m)
        return result

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: (3*x^2 - 2*y + 7/3).constant_coefficient()
            7/3
            sage: x.constant_coefficient()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef ulong *exp
        cdef fmpq_t coeff
        cdef Rational c
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        fmpq_init(coeff)
        try:
            for i in range(n):
                exp[i] = 0
            fmpq_mpoly_get_coeff_fmpq_ui(coeff, self._poly, exp, R._ctx)
            c = Rational.__new__(Rational)
            fmpq_get_mpq(c.value, coeff)
        finally:
            fmpq_clear(coeff)
            sig_free(exp)
        return c

    def monomial_coefficient(self, mon):
        r"""
        Return the coefficient of the monomial ``mon`` in this polynomial.

        INPUT:

        - ``mon`` -- a monomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y/2 - x*y + 5/7
            sage: f.monomial_coefficient(x^2*y)
            3/2
            sage: f.monomial_coefficient(x*y)
            -1
            sage: f.monomial_coefficient(R(1))
            5/7
            sage: f.monomial_coefficient(x)
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint m
        cdef fmpq_t coeff
        cdef Rational c
        if isinstance(mon, MPolynomial_rational_flint) and (<MPolynomial_rational_flint>mon)._parent is R:
            m = <MPolynomial_rational_flint>mon
        else:
            m = <MPolynomial_rational_flint>R(mon)
        if fmpq_mpoly_is_zero(m._poly, R._ctx):
            raise ValueError("mon must not be equal to 0")
        fmpq_init(coeff)
        sig_on()
        fmpq_mpoly_get_coeff_fmpq_monomial(coeff, self._poly, m._poly, R._ctx)
        sig_off()
        c = Rational.__new__(Rational)
        fmpq_get_mpq(c.value, coeff)
        fmpq_clear(coeff)
        return c

    def coefficient(self, degrees):
        r"""
        Return the coefficient of a monomial specified by ``degrees``.

        INPUT:

        - ``degrees`` -- a dictionary mapping generators to nonnegative
          integers, or a tuple/list of nonnegative integers giving the
          exponent of each variable

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y/2 - x*y + 5/7
            sage: f.coefficient({x: 2, y: 1})
            3/2
            sage: f.coefficient((2, 1))
            3/2
            sage: f.coefficient((0, 0))
            5/7
            sage: f.coefficient((1, 0))
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef ulong *exp
        cdef fmpq_t coeff
        cdef Rational c

        if isinstance(degrees, dict):
            names = R.variable_names()
            d = [0] * n
            for k, v in degrees.items():
                key = str(k)
                if key not in names:
                    raise ValueError("{} is not a variable of {}".format(k, R))
                d[names.index(key)] = int(v)
            degrees = d
        else:
            degrees = list(degrees)
            if len(degrees) != n:
                raise ValueError("degrees must have length {}".format(n))

        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        fmpq_init(coeff)
        try:
            for i in range(n):
                exp[i] = degrees[i]
            fmpq_mpoly_get_coeff_fmpq_ui(coeff, self._poly, exp, R._ctx)
            c = Rational.__new__(Rational)
            fmpq_get_mpq(c.value, coeff)
        finally:
            fmpq_clear(coeff)
            sig_free(exp)
        return c

    # ===== Leading term/coefficient/monomial =====

    def lc(self):
        r"""
        Return the leading coefficient of this polynomial.

        Returns 0 if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint', order='lex')
            sage: f = 3*x^2*y/2 - x*y^3 + 5/7
            sage: f.lc()
            3/2
            sage: R(0).lc()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef fmpq_t coeff
        cdef Rational c
        if fmpq_mpoly_is_zero(self._poly, R._ctx):
            return Rational(0)
        fmpq_init(coeff)
        fmpq_mpoly_get_term_coeff_fmpq(coeff, self._poly, 0, R._ctx)
        c = Rational.__new__(Rational)
        fmpq_get_mpq(c.value, coeff)
        fmpq_clear(coeff)
        return c

    def lm(self):
        r"""
        Return the leading monomial of this polynomial.

        Returns 0 if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint', order='lex')
            sage: f = 3*x^2*y/2 - x*y^3 + 5/7
            sage: f.lm()
            x^2*y
            sage: R(0).lm()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint m
        if fmpq_mpoly_is_zero(self._poly, R._ctx):
            return R._new_element()
        m = R._new_element()
        sig_on()
        fmpq_mpoly_get_term_monomial(m._poly, self._poly, 0, R._ctx)
        sig_off()
        return m

    def lt(self):
        r"""
        Return the leading term of this polynomial.

        Returns 0 if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint', order='lex')
            sage: f = 3*x^2*y/2 - x*y^3 + 5/7
            sage: f.lt()
            3/2*x^2*y
            sage: R(0).lt()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint t
        if fmpq_mpoly_is_zero(self._poly, R._ctx):
            return R._new_element()
        t = R._new_element()
        sig_on()
        fmpq_mpoly_get_term(t._poly, self._poly, 0, R._ctx)
        sig_off()
        return t

    def term_content(self):
        r"""
        Return the GCD of the (nonzero) terms of this polynomial.

        Over `\QQ`, the result is a monomial with coefficient `1`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: (6*x^2*y^2 - 4*x^3*y).term_content()
            x^2*y
            sage: R(0).term_content()
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        sig_on()
        fmpq_mpoly_term_content(result._poly, self._poly, R._ctx)
        sig_off()
        return result

    def monic(self):
        r"""
        Return this polynomial divided by its leading coefficient.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: (3*x^2/2 - y).monic()
            x^2 - 2/3*y

        TESTS::

            sage: R(0).monic()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot make the zero polynomial monic
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        if self.is_zero():
            raise ZeroDivisionError("cannot make the zero polynomial monic")
        cdef MPolynomial_rational_flint result = R._new_element()
        sig_on()
        fmpq_mpoly_make_monic(result._poly, self._poly, R._ctx)
        sig_off()
        return result

    # ===== Arithmetic =====

    cpdef _add_(self, other):
        r"""
        Return the sum of this polynomial and ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: x + y
            x + y
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint res = self._new()
        sig_on()
        fmpq_mpoly_add(res._poly, self._poly,
                       (<MPolynomial_rational_flint>other)._poly, R._ctx)
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        Return the difference of this polynomial and ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: x - y
            x - y
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint res = self._new()
        sig_on()
        fmpq_mpoly_sub(res._poly, self._poly,
                       (<MPolynomial_rational_flint>other)._poly, R._ctx)
        sig_off()
        return res

    cpdef _neg_(self):
        r"""
        Return the negation of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: -(x + y)
            -x - y
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint res = self._new()
        sig_on()
        fmpq_mpoly_neg(res._poly, self._poly, R._ctx)
        sig_off()
        return res

    cpdef _lmul_(self, Element scalar):
        r"""
        Return this polynomial multiplied by the scalar ``scalar``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: 3*x
            3*x
            sage: x * (2/3)
            2/3*x
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint res = self._new()
        cdef fmpq_t c
        cdef Rational s = scalar if isinstance(scalar, Rational) else QQ(scalar)
        fmpq_init(c)
        fmpq_set_mpq(c, s.value)
        sig_on()
        fmpq_mpoly_scalar_mul_fmpq(res._poly, self._poly, c, R._ctx)
        sig_off()
        fmpq_clear(c)
        return res

    cpdef _mul_(self, other):
        r"""
        Return the product of this polynomial and ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: (x + 1) * (y - 1)
            x*y - x + y - 1
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint res = self._new()
        sig_on()
        fmpq_mpoly_mul(res._poly, self._poly,
                       (<MPolynomial_rational_flint>other)._poly, R._ctx)
        sig_off()
        return res

    def __pow__(self, exp, mod):
        r"""
        Return this polynomial raised to the power ``exp``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: (x + y)^3
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
        """
        if mod is not None:
            raise NotImplementedError("modular exponentiation not supported")
        if exp < 0:
            raise ValueError("exponent must be non-negative")
        cdef MPolynomial_rational_flint self_ = self
        cdef MPolynomialRing_rational_flint R = self_._parent
        cdef MPolynomial_rational_flint res = self_._new()
        cdef int ok
        sig_on()
        ok = fmpq_mpoly_pow_ui(res._poly, self_._poly, <ulong>exp, R._ctx)
        sig_off()
        if not ok:
            raise ArithmeticError("power failed (exponent too large?)")
        return res

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare this polynomial with ``other``.

        Only equality and inequality are supported.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: x == x
            True
            sage: x == y
            False
            sage: x != y
            True
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef bint eq = fmpq_mpoly_equal(
            self._poly,
            (<MPolynomial_rational_flint>other)._poly,
            R._ctx)
        if op == Py_EQ:
            return bool(eq)
        if op == Py_NE:
            return not eq
        return NotImplemented

    # ===== Division and related =====

    def divides(self, other):
        r"""
        Return ``True`` if this polynomial divides ``other``.

        INPUT:

        - ``other`` -- a polynomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2 - y^2
            sage: g = x - y
            sage: g.divides(f)
            True
            sage: (x + 1).divides(f)
            False
            sage: R(0).divides(f)
            False
            sage: g.divides(R(0))
            True
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        if self.is_zero():
            return R(other).is_zero()
        cdef MPolynomial_rational_flint q = R._new_element()
        cdef MPolynomial_rational_flint b = R(other)
        sig_on()
        result = fmpq_mpoly_divides(q._poly, b._poly, self._poly, R._ctx)
        sig_off()
        return bool(result)

    def quo_rem(self, other):
        r"""
        Return the quotient and remainder of the division of this polynomial
        by ``other``.

        INPUT:

        - ``other`` -- a nonzero polynomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2*y + 3*x - 1
            sage: g = x*y - 1
            sage: q, r = f.quo_rem(g)
            sage: q * g + r == f
            True
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint q = R._new_element()
        cdef MPolynomial_rational_flint r = R._new_element()
        cdef MPolynomial_rational_flint b = R(other)
        if b.is_zero():
            raise ZeroDivisionError("cannot divide by zero")
        sig_on()
        fmpq_mpoly_divrem(q._poly, r._poly, self._poly, b._poly, R._ctx)
        sig_off()
        return q, r

    def gcd(self, other):
        r"""
        Return the greatest common divisor of this polynomial and ``other``.

        The result is monic (over `\QQ`).

        INPUT:

        - ``other`` -- a polynomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2 - y^2
            sage: g = x - y
            sage: f.gcd(g)
            x - y
            sage: (2*x).gcd(4*y)
            1
            sage: f.gcd(R(0))
            x^2 - y^2
            sage: R(0).gcd(R(0))
            0
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        cdef MPolynomial_rational_flint b = R(other)
        sig_on()
        ok = fmpq_mpoly_gcd(result._poly, self._poly, b._poly, R._ctx)
        sig_off()
        if not ok:
            raise ArithmeticError("GCD computation failed")
        return result

    def sqrt(self):
        r"""
        Return the square root of this polynomial if it is a perfect square,
        otherwise raise a :class:`ValueError`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = (x^2 + 2*x*y + y^2)/4
            sage: f.sqrt()
            1/2*x + 1/2*y
            sage: R(9/4).sqrt()
            3/2
            sage: x.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: x is not a perfect square
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        sig_on()
        ok = fmpq_mpoly_sqrt(result._poly, self._poly, R._ctx)
        sig_off()
        if not ok:
            raise ValueError("{} is not a perfect square".format(self))
        return result

    # ===== Calculus =====

    def _derivative(self, var):
        r"""
        Return the partial derivative of this polynomial with respect to
        ``var``.

        INPUT:

        - ``var`` -- a generator of the parent ring (or its name)

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^3*y^2/4 + 5*y^2 + 3*x + 2
            sage: f._derivative(x)
            9/4*x^2*y^2 + 3
            sage: f._derivative(y)
            3/2*x^3*y + 10*y
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        cdef slong var_idx
        if var is None:
            raise ValueError("you must specify the variable with respect to which to differentiate")
        names = R.variable_names()
        var_name = str(var)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(var, R))
        var_idx = names.index(var_name)
        sig_on()
        fmpq_mpoly_derivative(result._poly, self._poly, var_idx, R._ctx)
        sig_off()
        return result

    def derivative(self, *args):
        r"""
        Return the (iterated) partial derivative of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^3 * y^2
            sage: f.derivative(x)
            3*x^2*y^2
            sage: f.derivative(x, x)
            6*x*y^2
            sage: f.derivative(x, y)
            6*x^2*y
        """
        result = self
        for v in args:
            result = result._derivative(v)
        return result

    def integral(self, var):
        r"""
        Return the formal indefinite integral of this polynomial with respect
        to ``var``, with zero constant term.

        INPUT:

        - ``var`` -- a generator of the parent ring (or its name)

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y + y
            sage: f.integral(x)
            x^3*y + x*y
            sage: f.integral(y)
            3/2*x^2*y^2 + 1/2*y^2
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        cdef slong var_idx
        if var is None:
            raise ValueError("you must specify the variable with respect to which to integrate")
        names = R.variable_names()
        var_name = str(var)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(var, R))
        var_idx = names.index(var_name)
        sig_on()
        fmpq_mpoly_integral(result._poly, self._poly, var_idx, R._ctx)
        sig_off()
        return result

    # ===== Factorization =====

    def factor(self):
        r"""
        Return the factorization of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = (x^2 - y^2) * (x + 1)
            sage: F = f.factor()
            sage: F.unit() * prod(g^e for g, e in F) == f
            True

        TESTS::

            sage: R(0).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: factorization of 0 is not defined
            sage: R(3/4).factor()
            3/4
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        if self.is_zero():
            raise ArithmeticError("factorization of 0 is not defined")

        cdef fmpq_mpoly_factor_t fac
        cdef fmpq_t const_fmpq
        cdef fmpq_mpoly_t base_poly
        cdef MPolynomial_rational_flint base
        cdef Rational constant
        cdef slong i, length, exp
        cdef int ok

        fmpq_mpoly_factor_init(fac, R._ctx)
        fmpq_init(const_fmpq)
        fmpq_mpoly_init(base_poly, R._ctx)
        factors = []
        try:
            sig_on()
            ok = fmpq_mpoly_factor(fac, self._poly, R._ctx)
            sig_off()
            if not ok:
                raise ArithmeticError("factorization failed")
            fmpq_mpoly_factor_get_constant_fmpq(const_fmpq, fac, R._ctx)
            constant = Rational.__new__(Rational)
            fmpq_get_mpq(constant.value, const_fmpq)
            length = fmpq_mpoly_factor_length(fac, R._ctx)
            for i in range(length):
                fmpq_mpoly_factor_get_base(base_poly, fac, i, R._ctx)
                exp = fmpq_mpoly_factor_get_exp_si(fac, i, R._ctx)
                base = R._new_element()
                sig_on()
                fmpq_mpoly_set(base._poly, base_poly, R._ctx)
                sig_off()
                factors.append((base, int(exp)))
        finally:
            fmpq_mpoly_clear(base_poly, R._ctx)
            fmpq_clear(const_fmpq)
            fmpq_mpoly_factor_clear(fac, R._ctx)

        return Factorization(factors, unit=constant, sort=False)

    # ===== Resultant, discriminant =====

    def resultant(self, other, variable=None):
        r"""
        Return the resultant of this polynomial and ``other`` with respect
        to ``variable``.

        If ``variable`` is not provided, the first variable of the ring is
        used.

        INPUT:

        - ``other`` -- a polynomial in the same ring
        - ``variable`` -- (optional) a generator of the ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2 - y
            sage: g = x - 1
            sage: f.resultant(g)
            -y + 1
            sage: f.resultant(g, y)
            x - 1
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        cdef MPolynomial_rational_flint b = R(other)
        cdef slong var_idx

        if variable is None:
            var_idx = 0
        else:
            names = R.variable_names()
            var_name = str(variable)
            if var_name not in names:
                raise ValueError("{} is not a variable of {}".format(variable, R))
            var_idx = names.index(var_name)

        sig_on()
        ok = fmpq_mpoly_resultant(result._poly, self._poly, b._poly, var_idx, R._ctx)
        sig_off()
        if not ok:
            raise ArithmeticError("resultant computation failed")
        return result

    def discriminant(self, variable):
        r"""
        Return the discriminant of this polynomial with respect to
        ``variable``.

        INPUT:

        - ``variable`` -- a generator of the parent ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 'x,y', implementation='flint')
            sage: f = x^2 + y*x + 1
            sage: f.discriminant(x)
            y^2 - 4
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef MPolynomial_rational_flint result = R._new_element()
        cdef slong var_idx
        names = R.variable_names()
        var_name = str(variable)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(variable, R))
        var_idx = names.index(var_name)
        sig_on()
        ok = fmpq_mpoly_discriminant(result._poly, self._poly, var_idx, R._ctx)
        sig_off()
        if not ok:
            raise ArithmeticError("discriminant computation failed")
        return result

    # ===== Evaluation and substitution =====

    def __call__(self, *args, **kwds):
        r"""
        Evaluate this polynomial at the given values.

        Positional arguments are assigned to the variables in order. Keyword
        arguments name variables explicitly. When all assigned values are
        rationals and all variables are assigned, the result is a
        :class:`~sage.rings.rational.Rational`; otherwise the result is a
        polynomial in the same ring.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: f = x^2 + y*z + 1
            sage: f(1, 2, 3)
            8
            sage: f(1, 2, 3).parent()
            Rational Field
            sage: f(x, y, z) == f
            True
            sage: f(y, x, z)
            y^2 + x*z + 1
            sage: f(x=2)
            y*z + 5
            sage: f(x=1, y=2, z=3)
            8
            sage: f(1/2, 1/3, 1)
            19/12
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef fmpq_t res_fmpq
        cdef fmpq **vals
        cdef fmpq_t *fmpq_storage
        cdef Rational result_q

        if kwds and args:
            raise TypeError("cannot mix positional and keyword arguments")

        # allow f((a, b, c)) as a shortcut for f(a, b, c)
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = tuple(args[0])

        names = R.variable_names()
        if args:
            if len(args) != n:
                raise TypeError("expected {} positional arguments, got {}".format(n, len(args)))
            values = list(args)
        else:
            values = [None] * n
            for k, v in kwds.items():
                if k not in names:
                    raise ValueError("{} is not a variable of {}".format(k, R))
                values[names.index(k)] = v

        # fast path: all values are rationals (or integers) and all assigned
        if all(v is not None and isinstance(v, (int, Integer, Rational)) for v in values):
            fmpq_storage = <fmpq_t *>sig_malloc(n * sizeof(fmpq_t))
            vals = <fmpq **>sig_malloc(n * sizeof(fmpq *))
            if fmpq_storage == NULL or vals == NULL:
                if fmpq_storage != NULL:
                    sig_free(fmpq_storage)
                if vals != NULL:
                    sig_free(vals)
                raise MemoryError
            for i in range(n):
                fmpq_init(fmpq_storage[i])
                fmpq_set_mpq(fmpq_storage[i], (<Rational>Rational(values[i])).value)
                vals[i] = fmpq_storage[i]
            fmpq_init(res_fmpq)
            try:
                sig_on()
                ok = fmpq_mpoly_evaluate_all_fmpq(res_fmpq, self._poly, vals, R._ctx)
                sig_off()
                if not ok:
                    raise ArithmeticError("evaluation failed")
                result_q = Rational.__new__(Rational)
                fmpq_get_mpq(result_q.value, res_fmpq)
            finally:
                fmpq_clear(res_fmpq)
                for i in range(n):
                    fmpq_clear(fmpq_storage[i])
                sig_free(fmpq_storage)
                sig_free(vals)
            return result_q

        # general path: do substitution (possibly partial)
        fixed = {}
        for i in range(n):
            if values[i] is not None:
                fixed[names[i]] = values[i]
        return self.subs(fixed)

    def subs(self, fixed=None, **kw):
        r"""
        Substitute variables in this polynomial and return the result.

        Each substitution value may be a rational or a polynomial in the
        same ring. Variables not mentioned are left unchanged.

        INPUT:

        - ``fixed`` -- (optional) dict mapping generators or variable names
          to values
        - ``**kw`` -- variable names as keyword arguments

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 'x,y,z', implementation='flint')
            sage: f = x^2 + y^2 + z^2
            sage: f.subs(x=1, y=2)
            z^2 + 5
            sage: f.subs({x: y + 1})
            2*y^2 + z^2 + 2*y + 1
            sage: f.subs(x=1, y=2, z=3)
            14
            sage: f.subs({x: y, y: x})
            x^2 + y^2 + z^2
        """
        cdef MPolynomialRing_rational_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef fmpq_mpoly_struct **C
        cdef MPolynomial_rational_flint res

        sub = {}
        if fixed is not None:
            for k, v in fixed.items():
                if isinstance(k, MPolynomial_rational_flint):
                    k = str(k)
                sub[str(k)] = v
        for k, v in kw.items():
            sub[k] = v

        var_names = R.variable_names()
        polys = []
        for i in range(n):
            name = var_names[i]
            if name in sub:
                v = sub[name]
                if not isinstance(v, MPolynomial_rational_flint):
                    v = R(v)
                polys.append(<MPolynomial_rational_flint>v)
            else:
                polys.append(<MPolynomial_rational_flint>R.gen(i))

        C = <fmpq_mpoly_struct **>sig_malloc(n * sizeof(fmpq_mpoly_struct *))
        if C == NULL:
            raise MemoryError
        for i in range(n):
            C[i] = &(<MPolynomial_rational_flint>polys[i])._poly[0]
        res = R._new_element()
        sig_on()
        ok = fmpq_mpoly_compose_fmpq_mpoly(res._poly, self._poly, C, R._ctx, R._ctx)
        sig_off()
        sig_free(C)
        if not ok:
            raise ArithmeticError("composition failed")
        return res
