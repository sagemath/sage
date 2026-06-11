# distutils: libraries = flint
# distutils: depends = flint/fmpz_mpoly.h flint/fmpz_mpoly_factor.h
r"""
Multivariate polynomials over `\ZZ`, implemented using FLINT

EXAMPLES::

    sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
    sage: f = 3*x^2*y - z + 1; f
    3*x^2*y - z + 1
    sage: f * f
    9*x^4*y^2 - 6*x^2*y*z + 6*x^2*y + z^2 - 2*z + 1

.. automethod:: MPolynomial_integer_flint._add_
.. automethod:: MPolynomial_integer_flint._sub_
.. automethod:: MPolynomial_integer_flint._mul_
.. automethod:: MPolynomial_integer_flint._neg_
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

from sage.libs.flint.fmpz cimport (
    fmpz_init, fmpz_clear, fmpz_set, fmpz_set_mpz, fmpz_get_mpz,
    fmpz_zero, fmpz_one, fmpz_is_zero, fmpz_is_one,
    fmpz_gcd, fmpz_abs, fmpz_sgn, fmpz_cmp_si)
from sage.libs.flint.fmpz_mpoly cimport *
from sage.libs.flint.fmpz_mpoly_factor cimport (
    fmpz_mpoly_factor_init, fmpz_mpoly_factor_clear,
    fmpz_mpoly_factor_length, fmpz_mpoly_factor_get_constant_fmpz,
    fmpz_mpoly_factor_get_base, fmpz_mpoly_factor_get_exp_si,
    fmpz_mpoly_factor)
from sage.libs.flint.types cimport (
    fmpz_t, fmpz_mpoly_t, fmpz_mpoly_ctx_t, fmpz_mpoly_struct,
    fmpz_mpoly_factor_t,
    ordering_t, ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX, ulong, slong)

from sage.cpython.string cimport str_to_bytes, char_to_str

from sage.rings.integer cimport Integer, _Integer_from_mpz
from sage.rings.integer_ring import ZZ
from sage.structure.element cimport Element
from sage.structure.factorization import Factorization
from sage.structure.richcmp cimport rich_to_bool

from sage.rings.polynomial.multi_polynomial cimport MPolynomial_flint as MPolynomial_flint_base
from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial_integer_flint cimport (
    MPolynomialRing_integer_flint,
    MPolynomial_integer_flint)
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
        "term order '{}' is not supported by FLINT fmpz_mpoly".format(name))


def _unpickle_MPolynomialRing_integer_flint(n, names, order):
    r"""
    Helper for unpickling :class:`MPolynomialRing_integer_flint`.

    TESTS::

        sage: R = PolynomialRing(ZZ, 'x,y', implementation='flint')
        sage: loads(dumps(R)) is R
        True
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    return PolynomialRing(ZZ, n, names=names, order=order, implementation='flint')


def _unpickle_MPolynomial_integer_flint(parent, coeffs):
    r"""
    Helper for unpickling :class:`MPolynomial_integer_flint`.

    TESTS::

        sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
        sage: f = 3*x^2*y - x + 1
        sage: loads(dumps(f)) == f
        True
    """
    return parent(coeffs)


cdef class MPolynomialRing_integer_flint(MPolynomialRing_base):
    r"""
    Multivariate polynomial ring over `\ZZ`, implemented via FLINT.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
        sage: R
        Multivariate Polynomial Ring in x, y over Integer Ring (using FLINT)
    """

    def __cinit__(self):
        # nvars = -1 signals that fmpz_mpoly_ctx_init has not yet been called.
        self._ctx[0].minfo[0].nvars = -1

    def __init__(self, base_ring, n, names, order='degrevlex'):
        r"""
        Construct a multivariate polynomial ring over `\ZZ` using FLINT.

        INPUT:

        - ``base_ring`` -- must be `\ZZ`
        - ``n`` -- number of variables (positive integer)
        - ``names`` -- variable names
        - ``order`` -- term order (``'lex'``, ``'deglex'``, or ``'degrevlex'``)

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_flint import MPolynomialRing_integer_flint
            sage: R = MPolynomialRing_integer_flint(ZZ, 3, ('x','y','z'), 'degrevlex')
            sage: R
            Multivariate Polynomial Ring in x, y, z over Integer Ring (using FLINT)
        """
        if base_ring is not ZZ:
            raise TypeError("base ring must be ZZ")
        MPolynomialRing_base.__init__(self, base_ring, n, names, order)
        cdef ordering_t ord = _term_order_to_flint(self._term_order)
        fmpz_mpoly_ctx_init(self._ctx, n, ord)

    def __dealloc__(self):
        if self._ctx[0].minfo[0].nvars != -1:
            fmpz_mpoly_ctx_clear(self._ctx)

    Element = MPolynomial_integer_flint

    def __hash__(self):
        r"""
        Return a hash of this ring.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: hash(R) == hash(R)
            True
        """
        from sage.structure.category_object import CategoryObject
        return CategoryObject.__hash__(self)

    def __reduce__(self):
        r"""
        Return data for pickling this ring.

        TESTS::

            sage: R = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: loads(dumps(R)) is R
            True
        """
        return (_unpickle_MPolynomialRing_integer_flint,
                (self._ngens, self.variable_names(), self._term_order))

    def _repr_(self):
        r"""
        Return a string representation of this ring.

        EXAMPLES::

            sage: PolynomialRing(ZZ, 'x,y', implementation='flint')
            Multivariate Polynomial Ring in x, y over Integer Ring (using FLINT)
        """
        return "Multivariate Polynomial Ring in {} over Integer Ring (using FLINT)".format(
            ", ".join(self.variable_names()))

    def gen(self, int n=0):
        r"""
        Return the ``n``-th generator of this ring.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: R.gen(0)
            x
            sage: R.gen(2)
            z
        """
        if n < 0 or n >= self._ngens:
            raise ValueError("generator index out of range")
        cdef MPolynomial_integer_flint g = self._new_element()
        sig_on()
        fmpz_mpoly_gen(g._poly, n, self._ctx)
        sig_off()
        return g

    cdef MPolynomial_integer_flint _new_element(self):
        cdef MPolynomial_integer_flint f = \
            MPolynomial_integer_flint.__new__(MPolynomial_integer_flint)
        f._parent = self
        fmpz_mpoly_init(f._poly, self._ctx)
        return f

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into an element of this ring.

        Accepted inputs include integers, strings, dictionaries mapping exponent
        tuples to coefficients, and multivariate polynomials over a compatible
        ring.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R(3)
            3
            sage: R("x^2*y - 2*x + 1")
            x^2*y - 2*x + 1
            sage: R({(2,1): 5, (0,0): -1})
            5*x^2*y - 1
            sage: R(x + y)
            x + y
        """
        cdef MPolynomial_integer_flint f
        cdef Integer c
        cdef ETuple e
        cdef fmpz_t coeff
        cdef ulong *exp
        cdef slong n = self._ngens
        cdef int i
        cdef bytes bstr
        cdef const char **cnames

        if isinstance(x, MPolynomial_integer_flint) and x.parent() is self:
            f = self._new_element()
            sig_on()
            fmpz_mpoly_set(f._poly, (<MPolynomial_integer_flint>x)._poly, self._ctx)
            sig_off()
            return f

        if isinstance(x, int):
            x = Integer(x)
        if isinstance(x, Integer):
            f = self._new_element()
            fmpz_init(coeff)
            fmpz_set_mpz(coeff, (<Integer>x).value)
            sig_on()
            fmpz_mpoly_set_fmpz(f._poly, coeff, self._ctx)
            sig_off()
            fmpz_clear(coeff)
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
            ok = fmpz_mpoly_set_str_pretty(f._poly, bstr, cnames, self._ctx)
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
            fmpz_init(coeff)
            try:
                for key, val in x.items():
                    if not isinstance(val, Integer):
                        val = ZZ(val)
                    c = <Integer>val
                    for i in range(n):
                        exp[i] = key[i]
                    fmpz_set_mpz(coeff, c.value)
                    sig_on()
                    fmpz_mpoly_set_coeff_fmpz_ui(f._poly, coeff, exp, self._ctx)
                    sig_off()
            finally:
                fmpz_clear(coeff)
                sig_free(exp)
            return f

        # try via _mpoly_dict_recursive for MPolynomial types
        from sage.rings.polynomial.multi_polynomial import MPolynomial
        if isinstance(x, MPolynomial):
            return self._element_constructor_(
                x._mpoly_dict_recursive(self.variable_names(), ZZ))

        # last resort: coerce to ZZ
        return self._element_constructor_(ZZ(x))

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into this ring.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R.has_coerce_map_from(ZZ)
            True

        Conversion from the libsingular backend and vice versa::

            sage: R_f = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R_s = PolynomialRing(ZZ, 'x,y', implementation='singular')
            sage: x_f, y_f = R_f.gens()
            sage: x_s, y_s = R_s.gens()
            sage: p_f = 2*x_f^2 - x_f*y_f + 3
            sage: p_s = 2*x_s^2 - x_s*y_s + 3
            sage: R_f(p_s)
            2*x^2 - x*y + 3
            sage: R_s(p_f)
            2*x^2 - x*y + 3
            sage: R_f(p_s) == p_f
            True
            sage: R_s(p_f) == p_s
            True
        """
        if R is ZZ:
            return True
        return MPolynomialRing_base._coerce_map_from_(self, R)


cdef class MPolynomial_integer_flint(MPolynomial_flint_base):
    r"""
    A multivariate polynomial over `\ZZ`, implemented via FLINT.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
        sage: f = 3*x^2 - y + 1; f
        3*x^2 - y + 1
        sage: type(f)
        <class 'sage.rings.polynomial.multi_polynomial_integer_flint.MPolynomial_integer_flint'>

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _mul_
    .. automethod:: _neg_
    """

    def __cinit__(self):
        # Python zero-initialises the memory; fmpz_mpoly_clear on a struct with
        # alloc=0 is a no-op, so __dealloc__ is safe even if fmpz_mpoly_init
        # was never called (e.g. if __new__ is used directly without _parent).
        pass

    def __dealloc__(self):
        cdef MPolynomialRing_integer_flint R = self._parent
        if R is not None:
            fmpz_mpoly_clear(self._poly, R._ctx)

    cdef MPolynomial_integer_flint _new(self):
        return (<MPolynomialRing_integer_flint>self._parent)._new_element()

    cpdef _new_constant_poly(self, scalar, parent):
        r"""
        Return a new constant polynomial with value ``scalar`` in ``parent``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: x._new_constant_poly(3, R)
            3
        """
        cdef MPolynomialRing_integer_flint R = parent
        cdef MPolynomial_integer_flint f = R._new_element()
        cdef fmpz_t coeff
        if not isinstance(scalar, Integer):
            scalar = ZZ(scalar)
        fmpz_init(coeff)
        fmpz_set_mpz(coeff, (<Integer>scalar).value)
        fmpz_mpoly_set_fmpz(f._poly, coeff, R._ctx)
        fmpz_clear(coeff)
        return f

    def __copy__(self):
        r"""
        Return a copy of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^2 - y + 1
            sage: g = copy(f)
            sage: f == g
            True
            sage: f is g
            False
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint f = R._new_element()
        sig_on()
        fmpz_mpoly_set(f._poly, self._poly, R._ctx)
        sig_off()
        return f

    def __deepcopy__(self, memo=None):
        r"""
        Return a deep copy of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
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

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: f = 7*x^3*y - 2*y*z + 1
            sage: loads(dumps(f)) == f
            True
        """
        return (_unpickle_MPolynomial_integer_flint,
                (self._parent, self.monomial_coefficients()))

    def _repr_(self):
        r"""
        Return a string representation of this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: 3*x^2*y - z + 1
            3*x^2*y - z + 1
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef int n = R._ngens
        cdef char *raw
        cdef const char **cnames = <const char **>sig_malloc(n * sizeof(char *))
        if cnames == NULL:
            raise MemoryError
        # keep byte strings alive for the duration of the call
        bnames = [str_to_bytes(v) for v in R.variable_names()]
        for i in range(n):
            cnames[i] = bnames[i]
        sig_on()
        raw = fmpz_mpoly_get_str_pretty(self._poly, cnames, R._ctx)
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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: bool(x)
            True
            sage: bool(R(0))
            False
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return not fmpz_mpoly_is_zero(self._poly, R._ctx)

    # ===== Predicates =====

    def is_zero(self):
        r"""
        Return ``True`` if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R(0).is_zero()
            True
            sage: (x - x).is_zero()
            True
            sage: x.is_zero()
            False
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return bool(fmpz_mpoly_is_zero(self._poly, R._ctx))

    def is_one(self):
        r"""
        Return ``True`` if this polynomial equals 1.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R(1).is_one()
            True
            sage: x.is_one()
            False
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return bool(fmpz_mpoly_is_one(self._poly, R._ctx))

    def is_constant(self):
        r"""
        Return ``True`` if this polynomial is a constant (an element of `\ZZ`).

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R(0).is_constant()
            True
            sage: R(5).is_constant()
            True
            sage: x.is_constant()
            False
            sage: (x + 1).is_constant()
            False

        Constant polynomials can be converted to `\ZZ`::

            sage: ZZ(R(5))
            5
            sage: ZZ(R(0))
            0
            sage: int(R(7))
            7
            sage: ZZ(x)
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return bool(fmpz_mpoly_is_fmpz(self._poly, R._ctx))

    def is_term(self):
        r"""
        Return ``True`` if this polynomial is a single (nonzero) term.

        A *term* is a product of a nonzero coefficient and a monomial in the
        variables.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: x.is_term()
            True
            sage: (3*x^2*y).is_term()
            True
            sage: (x + y).is_term()
            False
            sage: R(0).is_term()
            False
            sage: R(5).is_term()
            True
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return fmpz_mpoly_length(self._poly, R._ctx) == 1

    def is_monomial(self):
        r"""
        Return ``True`` if this polynomial is a monomial.

        A *monomial* is a product of variables with coefficient 1.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
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
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef fmpz_t c
        if fmpz_mpoly_length(self._poly, R._ctx) != 1:
            return False
        fmpz_init(c)
        fmpz_mpoly_get_term_coeff_fmpz(c, self._poly, 0, R._ctx)
        result = fmpz_is_one(c)
        fmpz_clear(c)
        return bool(result)

    def is_homogeneous(self):
        r"""
        Return ``True`` if this polynomial is homogeneous.

        A polynomial is homogeneous if all its terms have the same total degree.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: (x^2 + y*z).is_homogeneous()
            True
            sage: (x^2 + y).is_homogeneous()
            False
            sage: R(0).is_homogeneous()
            True
            sage: R(5).is_homogeneous()
            True
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong length = fmpz_mpoly_length(self._poly, R._ctx)
        cdef slong i, j
        cdef ulong *exp
        if length <= 1:
            return True
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        try:
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, 0, R._ctx)
            ref_deg = 0
            for j in range(n):
                ref_deg += exp[j]
            for i in range(1, length):
                fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, R._ctx)
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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (x^2 + 2*x*y + y^2).is_square()
            True
            sage: x.is_square()
            False
            sage: R(0).is_square()
            True
            sage: R(4).is_square()
            True
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return bool(fmpz_mpoly_is_square(self._poly, R._ctx))

    # ===== Information =====

    cpdef long number_of_terms(self) noexcept:
        r"""
        Return the number of nonzero terms.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2 + 3*x*y - y + 2
            sage: f.number_of_terms()
            4
            sage: R(0).number_of_terms()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        return fmpz_mpoly_length(self._poly, R._ctx)

    def degree(self, x=None):
        r"""
        Return the degree of this polynomial in the variable ``x``, or the
        total degree if ``x`` is ``None``.

        INPUT:

        - ``x`` -- (optional) a generator of the parent ring

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
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
        cdef MPolynomialRing_integer_flint R = self._parent
        if fmpz_mpoly_is_zero(self._poly, R._ctx):
            return -1
        if x is None:
            return self.total_degree()
        cdef slong var_idx
        names = R.variable_names()
        var_name = str(x)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(x, R))
        var_idx = names.index(var_name)
        return Integer(fmpz_mpoly_degree_si(self._poly, var_idx, R._ctx))

    def total_degree(self):
        r"""
        Return the total degree of this polynomial.

        The total degree is the maximum, over all nonzero terms, of the sum of
        the exponents.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: (x^2*y + y^3*z^4 + 1).total_degree()
            7
            sage: R(5).total_degree()
            0
            sage: R(0).total_degree()
            -1
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        if fmpz_mpoly_is_zero(self._poly, R._ctx):
            return -1
        return Integer(fmpz_mpoly_total_degree_si(self._poly, R._ctx))

    def degrees(self):
        r"""
        Return a tuple containing the degree of this polynomial in each
        variable.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: f = x^2*y + y^3*z^4
            sage: f.degrees()
            (2, 3, 4)
            sage: R(5).degrees()
            (0, 0, 0)
            sage: R(0).degrees()
            (0, 0, 0)
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef slong *degs = <slong *>sig_malloc(n * sizeof(slong))
        if degs == NULL:
            raise MemoryError
        try:
            fmpz_mpoly_degrees_si(degs, self._poly, R._ctx)
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

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: (x^2*z + 1).variables()
            (x, z)
            sage: R(5).variables()
            ()
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef int *used = <int *>sig_malloc(n * sizeof(int))
        if used == NULL:
            raise MemoryError
        try:
            fmpz_mpoly_used_vars(used, self._poly, R._ctx)
            result = tuple(R.gen(i) for i in range(n) if used[i])
        finally:
            sig_free(used)
        return result

    def variable(self, int i=0):
        r"""
        Return the ``i``-th variable occurring in this polynomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
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

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: (x^2*z + 1).nvariables()
            2
            sage: R(5).nvariables()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef int *used = <int *>sig_malloc(n * sizeof(int))
        if used == NULL:
            raise MemoryError
        try:
            fmpz_mpoly_used_vars(used, self._poly, R._ctx)
            result = sum(1 for i in range(n) if used[i])
        finally:
            sig_free(used)
        return result

    # ===== Coefficients =====

    def monomial_coefficients(self):
        r"""
        Return a dictionary of ``{exponent: coefficient}`` pairs.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 2*x^3 - x*y + 5
            sage: sorted(f.monomial_coefficients().items())
            [((0, 0), 5), ((1, 1), -1), ((3, 0), 2)]
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong length = fmpz_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef fmpz_t coeff
        cdef ulong *exp
        cdef Integer c
        fmpz_init(coeff)
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            fmpz_clear(coeff)
            raise MemoryError
        result = {}
        try:
            for i in range(length):
                fmpz_mpoly_get_term_coeff_fmpz(coeff, self._poly, i, R._ctx)
                fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, R._ctx)
                c = Integer.__new__(Integer)
                fmpz_get_mpz(c.value, coeff)
                result[ETuple([exp[j] for j in range(n)])] = c
        finally:
            fmpz_clear(coeff)
            sig_free(exp)
        return result

    def coefficients(self):
        r"""
        Return the list of coefficients of this polynomial, in the order
        induced by the term order.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y - x + 5
            sage: f.coefficients()
            [3, -1, 5]
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong length = fmpz_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef fmpz_t coeff
        cdef Integer c
        fmpz_init(coeff)
        result = []
        try:
            for i in range(length):
                fmpz_mpoly_get_term_coeff_fmpz(coeff, self._poly, i, R._ctx)
                c = Integer.__new__(Integer)
                fmpz_get_mpz(c.value, coeff)
                result.append(c)
        finally:
            fmpz_clear(coeff)
        return result

    def exponents(self, as_ETuples=True):
        r"""
        Return the list of exponent tuples of the nonzero terms of this
        polynomial, in the order induced by the term order.

        INPUT:

        - ``as_ETuples`` -- (default: ``True``) if ``True``, return a list of
          :class:`ETuple`; otherwise, return a list of plain tuples

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y - x + 5
            sage: f.exponents()
            [(2, 1), (1, 0), (0, 0)]
            sage: f.exponents(as_ETuples=False)
            [(2, 1), (1, 0), (0, 0)]
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong length = fmpz_mpoly_length(self._poly, R._ctx)
        cdef slong i, j
        cdef ulong *exp
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        result = []
        try:
            for i in range(length):
                fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, R._ctx)
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
        Return the list of monomials of this polynomial, in the order induced
        by the term order.

        Monomials are returned with coefficient 1.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y - x + 5
            sage: f.monomials()
            [x^2*y, x, 1]
            sage: R(0).monomials()
            []
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong length = fmpz_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef MPolynomial_integer_flint m
        result = []
        for i in range(length):
            m = R._new_element()
            sig_on()
            fmpz_mpoly_get_term_monomial(m._poly, self._poly, i, R._ctx)
            sig_off()
            result.append(m)
        return result

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (3*x^2 - 2*y + 7).constant_coefficient()
            7
            sage: x.constant_coefficient()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef ulong *exp
        cdef fmpz_t coeff
        cdef Integer c
        exp = <ulong *>sig_malloc(n * sizeof(ulong))
        if exp == NULL:
            raise MemoryError
        fmpz_init(coeff)
        try:
            for i in range(n):
                exp[i] = 0
            fmpz_mpoly_get_coeff_fmpz_ui(coeff, self._poly, exp, R._ctx)
            c = Integer.__new__(Integer)
            fmpz_get_mpz(c.value, coeff)
        finally:
            fmpz_clear(coeff)
            sig_free(exp)
        return c

    def monomial_coefficient(self, mon):
        r"""
        Return the coefficient of the monomial ``mon`` in this polynomial.

        INPUT:

        - ``mon`` -- a monomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y - x*y + 5
            sage: f.monomial_coefficient(x^2*y)
            3
            sage: f.monomial_coefficient(x*y)
            -1
            sage: f.monomial_coefficient(R(1))
            5
            sage: f.monomial_coefficient(x)
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint m
        cdef fmpz_t coeff
        cdef Integer c
        if isinstance(mon, MPolynomial_integer_flint) and (<MPolynomial_integer_flint>mon)._parent is R:
            m = <MPolynomial_integer_flint>mon
        else:
            m = <MPolynomial_integer_flint>R(mon)
        if fmpz_mpoly_is_zero(m._poly, R._ctx):
            raise ValueError("mon must not be equal to 0")
        fmpz_init(coeff)
        sig_on()
        fmpz_mpoly_get_coeff_fmpz_monomial(coeff, self._poly, m._poly, R._ctx)
        sig_off()
        c = Integer.__new__(Integer)
        fmpz_get_mpz(c.value, coeff)
        fmpz_clear(coeff)
        return c

    def coefficient(self, degrees):
        r"""
        Return the coefficient of a monomial specified by ``degrees``.

        INPUT:

        - ``degrees`` -- a dictionary mapping generators to nonnegative
          integers, or a tuple/list of nonnegative integers giving the
          exponent of each variable

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^2*y - x*y + 5
            sage: f.coefficient({x: 2, y: 1})
            3
            sage: f.coefficient((2, 1))
            3
            sage: f.coefficient((0, 0))
            5
            sage: f.coefficient((1, 0))
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef ulong *exp
        cdef fmpz_t coeff
        cdef Integer c

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
        fmpz_init(coeff)
        try:
            for i in range(n):
                exp[i] = degrees[i]
            fmpz_mpoly_get_coeff_fmpz_ui(coeff, self._poly, exp, R._ctx)
            c = Integer.__new__(Integer)
            fmpz_get_mpz(c.value, coeff)
        finally:
            fmpz_clear(coeff)
            sig_free(exp)
        return c

    # ===== Leading term/coefficient/monomial =====

    def lc(self):
        r"""
        Return the leading coefficient of this polynomial, with respect to
        the term order of the parent ring.

        Returns 0 if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint', order='lex')
            sage: f = 3*x^2*y - x*y^3 + 5
            sage: f.lc()
            3
            sage: R(0).lc()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef fmpz_t coeff
        cdef Integer c
        if fmpz_mpoly_is_zero(self._poly, R._ctx):
            return Integer(0)
        fmpz_init(coeff)
        fmpz_mpoly_get_term_coeff_fmpz(coeff, self._poly, 0, R._ctx)
        c = Integer.__new__(Integer)
        fmpz_get_mpz(c.value, coeff)
        fmpz_clear(coeff)
        return c

    def lm(self):
        r"""
        Return the leading monomial of this polynomial, with respect to the
        term order of the parent ring.

        The leading monomial is the leading term divided by its coefficient.
        Returns 0 if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint', order='lex')
            sage: f = 3*x^2*y - x*y^3 + 5
            sage: f.lm()
            x^2*y
            sage: R(0).lm()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint m
        if fmpz_mpoly_is_zero(self._poly, R._ctx):
            return R._new_element()
        m = R._new_element()
        sig_on()
        fmpz_mpoly_get_term_monomial(m._poly, self._poly, 0, R._ctx)
        sig_off()
        return m

    def lt(self):
        r"""
        Return the leading term of this polynomial, with respect to the
        term order of the parent ring.

        Returns 0 if this polynomial is zero.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint', order='lex')
            sage: f = 3*x^2*y - x*y^3 + 5
            sage: f.lt()
            3*x^2*y
            sage: R(0).lt()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint t
        if fmpz_mpoly_is_zero(self._poly, R._ctx):
            return R._new_element()
        t = R._new_element()
        sig_on()
        fmpz_mpoly_get_term(t._poly, self._poly, 0, R._ctx)
        sig_off()
        return t

    # ===== Content =====

    def content(self):
        r"""
        Return the content of this polynomial, that is, the positive gcd of
        its coefficients.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (6*x^2 - 4*x*y + 10).content()
            2
            sage: (x + y).content()
            1
            sage: R(0).content()
            0
            sage: R(-5).content()
            5
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong length = fmpz_mpoly_length(self._poly, R._ctx)
        cdef slong i
        cdef fmpz_t g, c
        cdef Integer result
        if length == 0:
            return Integer(0)
        fmpz_init(g)
        fmpz_init(c)
        try:
            fmpz_mpoly_get_term_coeff_fmpz(g, self._poly, 0, R._ctx)
            fmpz_abs(g, g)
            for i in range(1, length):
                fmpz_mpoly_get_term_coeff_fmpz(c, self._poly, i, R._ctx)
                fmpz_gcd(g, g, c)
                if fmpz_is_one(g):
                    break
            result = Integer.__new__(Integer)
            fmpz_get_mpz(result.value, g)
        finally:
            fmpz_clear(g)
            fmpz_clear(c)
        return result

    def primitive_part(self):
        r"""
        Return the primitive part of this polynomial.

        The primitive part is the polynomial divided by its content. By
        convention, the result has a positive leading coefficient.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (6*x^2 - 4*x*y + 10).primitive_part()
            3*x^2 - 2*x*y + 5
            sage: (-2*x - 4).primitive_part()
            x + 2
            sage: R(0).primitive_part()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        sig_on()
        fmpz_mpoly_primitive_part(result._poly, self._poly, R._ctx)
        sig_off()
        return result

    def term_content(self):
        r"""
        Return the GCD of the (nonzero) terms of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (6*x^2*y^2 - 4*x^3*y).term_content()
            2*x^2*y
            sage: (x + y).term_content()
            1
            sage: R(0).term_content()
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        sig_on()
        fmpz_mpoly_term_content(result._poly, self._poly, R._ctx)
        sig_off()
        return result

    # ===== Arithmetic =====

    cpdef _add_(self, other):
        r"""
        Return the sum of this polynomial and ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: x + y
            x + y
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint res = self._new()
        sig_on()
        fmpz_mpoly_add(res._poly, self._poly,
                       (<MPolynomial_integer_flint>other)._poly, R._ctx)
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        Return the difference of this polynomial and ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: x - y
            x - y
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint res = self._new()
        sig_on()
        fmpz_mpoly_sub(res._poly, self._poly,
                       (<MPolynomial_integer_flint>other)._poly, R._ctx)
        sig_off()
        return res

    cpdef _neg_(self):
        r"""
        Return the negation of this polynomial.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: -(x + y)
            -x - y
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint res = self._new()
        sig_on()
        fmpz_mpoly_neg(res._poly, self._poly, R._ctx)
        sig_off()
        return res

    cpdef _lmul_(self, Element scalar):
        r"""
        Return this polynomial multiplied by the scalar ``scalar``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: 3 * x
            3*x
            sage: x * 3
            3*x
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint res = self._new()
        cdef fmpz_t c
        cdef Integer s = scalar if isinstance(scalar, Integer) else ZZ(scalar)
        fmpz_init(c)
        fmpz_set_mpz(c, s.value)
        sig_on()
        fmpz_mpoly_scalar_mul_fmpz(res._poly, self._poly, c, R._ctx)
        sig_off()
        fmpz_clear(c)
        return res

    cpdef _mul_(self, other):
        r"""
        Return the product of this polynomial and ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (x + 1) * (y - 1)
            x*y - x + y - 1
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint res = self._new()
        sig_on()
        fmpz_mpoly_mul(res._poly, self._poly,
                       (<MPolynomial_integer_flint>other)._poly, R._ctx)
        sig_off()
        return res

    def __pow__(self, exp, mod):
        r"""
        Return this polynomial raised to the power ``exp``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: (x + y)^3
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
        """
        if mod is not None:
            raise NotImplementedError("modular exponentiation not supported")
        if exp < 0:
            raise ValueError("exponent must be non-negative")
        cdef MPolynomial_integer_flint self_ = self
        cdef MPolynomialRing_integer_flint R = self_._parent
        cdef MPolynomial_integer_flint res = self_._new()
        cdef int ok
        sig_on()
        ok = fmpz_mpoly_pow_ui(res._poly, self_._poly, <ulong>exp, R._ctx)
        sig_off()
        if not ok:
            raise ArithmeticError("power failed (exponent too large?)")
        return res

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare this polynomial with ``other``.

        Only equality and inequality are supported.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: x == x
            True
            sage: x == y
            False
            sage: x != y
            True
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef bint eq = fmpz_mpoly_equal(
            self._poly,
            (<MPolynomial_integer_flint>other)._poly,
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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
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
        cdef MPolynomialRing_integer_flint R = self._parent
        if self.is_zero():
            return R(other).is_zero()
        cdef MPolynomial_integer_flint q = R._new_element()
        cdef MPolynomial_integer_flint b = R(other)
        sig_on()
        result = fmpz_mpoly_divides(q._poly, b._poly, self._poly, R._ctx)
        sig_off()
        return bool(result)

    def quo_rem(self, other):
        r"""
        Return the quotient and remainder of the division of this polynomial
        by ``other``.

        Division is performed with respect to the term order of the ring.
        The result satisfies ``self == q * other + r`` where no monomial of
        ``r`` is divisible by the leading monomial of ``other``.

        INPUT:

        - ``other`` -- a nonzero polynomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2*y + 3*x - 1
            sage: g = x*y - 1
            sage: q, r = f.quo_rem(g)
            sage: q
            x
            sage: r
            4*x - 1
            sage: q * g + r == f
            True
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint q = R._new_element()
        cdef MPolynomial_integer_flint r = R._new_element()
        cdef MPolynomial_integer_flint b = R(other)
        if b.is_zero():
            raise ZeroDivisionError("cannot divide by zero")
        sig_on()
        fmpz_mpoly_divrem(q._poly, r._poly, self._poly, b._poly, R._ctx)
        sig_off()
        return q, r

    def gcd(self, other):
        r"""
        Return the greatest common divisor of this polynomial and ``other``.

        Over `\ZZ`, the result is the primitive GCD with positive leading
        coefficient.

        INPUT:

        - ``other`` -- a polynomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2 - y^2
            sage: g = x - y
            sage: f.gcd(g)
            x - y
            sage: (2*x).gcd(4*y)
            2
            sage: f.gcd(R(0))
            x^2 - y^2
            sage: R(0).gcd(R(0))
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        cdef MPolynomial_integer_flint b = R(other)
        sig_on()
        ok = fmpz_mpoly_gcd(result._poly, self._poly, b._poly, R._ctx)
        sig_off()
        if not ok:
            raise ArithmeticError("GCD computation failed")
        return result

    def lcm(self, other):
        r"""
        Return the least common multiple of this polynomial and ``other``.

        Defined as ``self * other / gcd(self, other)``. The result has
        positive leading coefficient (the primitive LCM over `\ZZ`).

        INPUT:

        - ``other`` -- a polynomial in the same ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2 - y^2
            sage: g = x - y
            sage: f.lcm(g)
            x^2 - y^2
            sage: (2*x).lcm(3*y)
            6*x*y
            sage: R(0).lcm(x)
            0
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint b = R(other)
        if self.is_zero() or b.is_zero():
            return R._new_element()
        g = self.gcd(b)
        q, r = (self * b).quo_rem(g)
        return q

    def sqrt(self):
        r"""
        Return the square root of this polynomial if it is a perfect square,
        otherwise raise a :class:`ValueError`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2 + 2*x*y + y^2
            sage: f.sqrt()
            x + y
            sage: R(4).sqrt()
            2
            sage: x.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: x is not a perfect square
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        sig_on()
        ok = fmpz_mpoly_sqrt(result._poly, self._poly, R._ctx)
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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f._derivative(x)
            9*x^2*y^2 + 3
            sage: f._derivative(y)
            6*x^3*y + 10*y
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        cdef slong var_idx
        if var is None:
            raise ValueError("you must specify the variable with respect to which to differentiate")
        names = R.variable_names()
        var_name = str(var)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(var, R))
        var_idx = names.index(var_name)
        sig_on()
        fmpz_mpoly_derivative(result._poly, self._poly, var_idx, R._ctx)
        sig_off()
        return result

    def derivative(self, *args):
        r"""
        Return the (iterated) partial derivative of this polynomial.

        Repeated arguments give higher-order partial derivatives.

        INPUT:

        - ``*args`` -- a sequence of variables; if a variable is repeated, the
          derivative is iterated

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
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

    # ===== Factorization =====

    def factor(self, proof=None):
        r"""
        Return the factorization of this polynomial.

        The factorization is returned as a :class:`Factorization` object,
        with each irreducible factor paired with its multiplicity, and a
        unit ``\pm 1`` collecting the integer content sign.

        INPUT:

        - ``proof`` -- ignored

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = (x^2 - y^2) * (x + 1)
            sage: F = f.factor()
            sage: F.unit() * prod(g^e for g, e in F) == f
            True
            sage: g = 6 * (x - y)^2
            sage: g.factor()
            2 * 3 * (x - y)^2

        TESTS::

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: R(0).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: factorization of 0 is not defined
            sage: R(1).factor()
            1
            sage: R(-12).factor()
            (-1) * 2^2 * 3
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        if self.is_zero():
            raise ArithmeticError("factorization of 0 is not defined")

        cdef fmpz_mpoly_factor_t fac
        cdef fmpz_t const_fmpz
        cdef fmpz_mpoly_t base_poly
        cdef MPolynomial_integer_flint base
        cdef Integer constant
        cdef slong i, length, exp
        cdef int ok

        fmpz_mpoly_factor_init(fac, R._ctx)
        fmpz_init(const_fmpz)
        fmpz_mpoly_init(base_poly, R._ctx)
        polynomial_factors = []
        try:
            sig_on()
            ok = fmpz_mpoly_factor(fac, self._poly, R._ctx)
            sig_off()
            if not ok:
                raise ArithmeticError("factorization failed")
            fmpz_mpoly_factor_get_constant_fmpz(const_fmpz, fac, R._ctx)
            constant = Integer.__new__(Integer)
            fmpz_get_mpz(constant.value, const_fmpz)
            length = fmpz_mpoly_factor_length(fac, R._ctx)
            for i in range(length):
                fmpz_mpoly_factor_get_base(base_poly, fac, i, R._ctx)
                exp = fmpz_mpoly_factor_get_exp_si(fac, i, R._ctx)
                base = R._new_element()
                sig_on()
                fmpz_mpoly_set(base._poly, base_poly, R._ctx)
                sig_off()
                polynomial_factors.append((base, int(exp)))
        finally:
            fmpz_mpoly_clear(base_poly, R._ctx)
            fmpz_clear(const_fmpz)
            fmpz_mpoly_factor_clear(fac, R._ctx)

        # split integer constant into sign and prime factorization
        unit = Integer(1)
        if constant < 0:
            unit = -unit
            constant = -constant
        factors = []
        if constant != 1:
            for p, e in constant.factor():
                factors.append((R(p), int(e)))
        factors.extend(polynomial_factors)

        return Factorization(factors, unit=unit, sort=False)

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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2 - y
            sage: g = x - 1
            sage: f.resultant(g)
            -y + 1
            sage: f.resultant(g, y)
            x - 1
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        cdef MPolynomial_integer_flint b = R(other)
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
        ok = fmpz_mpoly_resultant(result._poly, self._poly, b._poly, var_idx, R._ctx)
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

            sage: R.<x,y> = PolynomialRing(ZZ, 'x,y', implementation='flint')
            sage: f = x^2 + y*x + 1
            sage: f.discriminant(x)
            y^2 - 4
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef MPolynomial_integer_flint result = R._new_element()
        cdef slong var_idx
        names = R.variable_names()
        var_name = str(variable)
        if var_name not in names:
            raise ValueError("{} is not a variable of {}".format(variable, R))
        var_idx = names.index(var_name)
        sig_on()
        ok = fmpz_mpoly_discriminant(result._poly, self._poly, var_idx, R._ctx)
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
        integers and all variables are assigned, the result is an
        :class:`~sage.rings.integer.Integer`; otherwise the result is a
        polynomial in the same ring.

        INPUT:

        - ``*args`` -- values for the variables, in order
        - ``**kwds`` -- variable-name keyword arguments

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
            sage: f = x^2 + y*z + 1
            sage: f(1, 2, 3)
            8
            sage: f(1, 2, 3).parent()
            Integer Ring
            sage: f(x, y, z) == f
            True
            sage: f(y, x, z)
            y^2 + x*z + 1
            sage: f(x=2)
            y*z + 5
            sage: f(x=1, y=2, z=3)
            8
        """
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef fmpz_t res_fmpz
        cdef fmpz **vals
        cdef fmpz_t *fmpz_storage
        cdef MPolynomial_integer_flint p
        cdef Integer result_int

        if kwds and args:
            raise TypeError("cannot mix positional and keyword arguments")

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

        # fast path: all values are integers and all assigned
        if all(v is not None and isinstance(v, (int, Integer)) for v in values):
            fmpz_storage = <fmpz_t *>sig_malloc(n * sizeof(fmpz_t))
            vals = <fmpz **>sig_malloc(n * sizeof(fmpz *))
            if fmpz_storage == NULL or vals == NULL:
                if fmpz_storage != NULL:
                    sig_free(fmpz_storage)
                if vals != NULL:
                    sig_free(vals)
                raise MemoryError
            for i in range(n):
                fmpz_init(fmpz_storage[i])
                fmpz_set_mpz(fmpz_storage[i], (<Integer>Integer(values[i])).value)
                vals[i] = fmpz_storage[i]
            fmpz_init(res_fmpz)
            try:
                sig_on()
                ok = fmpz_mpoly_evaluate_all_fmpz(res_fmpz, self._poly, vals, R._ctx)
                sig_off()
                if not ok:
                    raise ArithmeticError("evaluation failed")
                result_int = Integer.__new__(Integer)
                fmpz_get_mpz(result_int.value, res_fmpz)
            finally:
                fmpz_clear(res_fmpz)
                for i in range(n):
                    fmpz_clear(fmpz_storage[i])
                sig_free(fmpz_storage)
                sig_free(vals)
            return result_int

        # general path: do substitution (possibly partial)
        fixed = {}
        for i in range(n):
            if values[i] is not None:
                fixed[names[i]] = values[i]
        return self.subs(fixed)

    def subs(self, fixed=None, **kw):
        r"""
        Substitute variables in this polynomial and return the result.

        Each substitution value may be an integer or a polynomial in the same
        ring. Variables not mentioned are left unchanged. All substitutions
        are applied simultaneously (not sequentially).

        INPUT:

        - ``fixed`` -- (optional) dict mapping generators or variable names to
          values
        - ``**kw`` -- variable names as keyword arguments

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 'x,y,z', implementation='flint')
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
        cdef MPolynomialRing_integer_flint R = self._parent
        cdef slong n = R._ngens
        cdef slong i
        cdef fmpz_mpoly_struct **C
        cdef MPolynomial_integer_flint res

        # build var_name → value mapping
        sub = {}
        if fixed is not None:
            for k, v in fixed.items():
                if isinstance(k, MPolynomial_integer_flint):
                    k = str(k)
                sub[str(k)] = v
        for k, v in kw.items():
            sub[k] = v

        # build list of substituting polynomials: one per variable
        var_names = R.variable_names()
        polys = []
        for i in range(n):
            name = var_names[i]
            if name in sub:
                v = sub[name]
                if not isinstance(v, MPolynomial_integer_flint):
                    v = R(v)
                polys.append(<MPolynomial_integer_flint>v)
            else:
                polys.append(<MPolynomial_integer_flint>R.gen(i))

        # call fmpz_mpoly_compose_fmpz_mpoly
        C = <fmpz_mpoly_struct **>sig_malloc(n * sizeof(fmpz_mpoly_struct *))
        if C == NULL:
            raise MemoryError
        for i in range(n):
            C[i] = &(<MPolynomial_integer_flint>polys[i])._poly[0]
        res = R._new_element()
        sig_on()
        ok = fmpz_mpoly_compose_fmpz_mpoly(res._poly, self._poly, C, R._ctx, R._ctx)
        sig_off()
        sig_free(C)
        if not ok:
            raise ArithmeticError("composition failed")
        return res
