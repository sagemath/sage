"""
Tate algebra element

A class for series in Tate algebras.

AUTHOR:

- Xavier Caruso, Thibaut Verron (2018-09)
"""

# ***************************************************************************
#    Copyright (C) 2018 Xavier Caruso <xavier.caruso@normalesup.org>
#                       Thibaut Verron <thibaut.verron@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport rich_to_bool_sgn
from sage.structure.element import coerce_binop

from sage.structure.element cimport Element
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.rings.polynomial.polydict cimport PolyDict
from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError


def _pushout_family(elements, initial=ZZ):
    """
    Return a parent in which ``initial`` and all elements
    in ``elements`` coerce.

    INPUT:

    - ``elements`` -- list of elements

    - ``initial`` -- a parent

    EXAMPLES::

        sage: from sage.rings.tate_algebra_element import _pushout_family

        sage: R = Zp(2)
        sage: A.<x,y> = TateAlgebra(R)
        sage: S.<a> = Zq(4)

        sage: _pushout_family([a, x, 3])
        Tate Algebra in x (val >= 0), y (val >= 0)
         over 2-adic Unramified Extension Field in a defined by x^2 + x + 1
    """
    from sage.structure.coerce_exceptions import CoercionException
    from sage.categories.pushout import pushout
    A = initial
    for elt in elements:
        if not isinstance(elt, Element):
            raise TypeError("cannot coerce all the elements to the same parent")
        try:
            A = pushout(A, elt.parent())
        except CoercionException:
            raise TypeError("cannot coerce all the elements to the same parent")
    return A


cdef class TateAlgebraTerm(MonoidElement):
    r"""
    A class for Tate algebra terms.

    A term in `K\{X_1,\dots,X_n\}` is the product of a coefficient in `K` and a
    monomial in the variables `X_1,\dots,X_n`.

    Those terms form a partially ordered monoid, with term multiplication and the
    term order of the parent Tate algebra.

    INPUT:

    - ``coeff`` -- an element in the base field

    - ``exponent`` -- tuple of length ``n``

    EXAMPLES::

        sage: R = Zp(2, print_mode='digits', prec=10)
        sage: A.<x,y> = TateAlgebra(R)
        sage: T = A.monoid_of_terms(); T
        Monoid of terms in x (val >= 0), y (val >= 0)
         over 2-adic Field with capped relative precision 10

        sage: T(2*x*y)
        ...00000000010*x*y
        sage: T(0)
        Traceback (most recent call last):
        ...
        TypeError: a term cannot be zero
    """
    def __init__(self, parent, coeff, exponent=None):
        """
        Initialize a Tate algebra term.

        INPUT:

        - ``coeff`` -- an element in the base field

        - ``exponent`` -- tuple

        TESTS::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()

            sage: t = T(x)
            sage: TestSuite(t).run()
        """
        MonoidElement.__init__(self, parent)
        field = parent.base_ring().fraction_field()
        if isinstance(coeff, TateAlgebraElement):
            coeff = coeff.leading_term()
        if isinstance(coeff, TateAlgebraTerm):
            if coeff.parent().variable_names() != self._parent.variable_names():
                raise ValueError("the variable names do not match")
            self._coeff = field((<TateAlgebraTerm>coeff)._coeff)
            self._exponent = (<TateAlgebraTerm>coeff)._exponent
        else:
            self._coeff = field(coeff)
            if self._coeff.is_zero():
                raise TypeError("a term cannot be zero")
            self._exponent = ETuple([0] * parent.ngens())
        if exponent is not None:
            if not isinstance(exponent, ETuple):
                exponent = ETuple(exponent)
            self._exponent = self._exponent.eadd(exponent)
        if len(self._exponent) != parent.ngens():
            raise ValueError("the length of the exponent does not match the number of variables")
        for i in self._exponent.nonzero_positions():
            if self._exponent[i] < 0:
                raise ValueError("only nonnegative exponents are allowed")
        if not parent.base_ring().is_field() and self.valuation() < 0:
            raise ValueError("this term is not in the ring of integers")

    def __hash__(self):
        """
        Return a hash of this term.

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R);
            sage: T = A.monoid_of_terms()
            sage: t = T(x^2)

            sage: hash(t) == hash((t.coefficient(), t.exponent()))
            True
        """
        return hash((self._coeff, self._exponent))

    cdef TateAlgebraTerm _new_c(self):
        r"""
        Fast creation of a Tate algebra term.

        TESTS::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R);
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: t = T(x*y); t
            ...0000000001*x*y
            sage: 2*t  # indirect doctest
            ...00000000010*x*y
        """
        cdef TateAlgebraTerm ans = TateAlgebraTerm.__new__(TateAlgebraTerm)
        ans._parent = self._parent
        return ans

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: A.<x,y> = TateAlgebra(Zp(2))
            sage: t = x.leading_term()
            sage: loads(dumps(t)) == t  # indirect doctest
            True
        """
        return TateAlgebraTerm, (self.parent(), self._coeff, self._exponent)

    def __bool__(self):
        r"""
        Return ``True`` if this term is nonzero, ``False`` otherwise.

        TESTS::

            sage: A.<x,y> = TateAlgebra(Zp(2))
            sage: t = x.leading_term()
            sage: bool(t)
            True
        """
        return bool(self._coeff)

    def _repr_(self):
        r"""
        Return a string representation of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T(2*x*y)  # indirect doctest
            ...00000000010*x*y
        """
        parent = self._parent
        if self._coeff._is_atomic() or (-self._coeff)._is_atomic():
            s = repr(self._coeff)
            if s == "1": s = ""
        else:
            s = "(%s)" % self._coeff
        for i in range(parent._ngens):
            if self._exponent[i] == 1:
                s += "*%s" % parent._names[i]
            elif self._exponent[i] > 1:
                s += "*%s^%s" % (parent._names[i], self._exponent[i])
        if s[0] == "*":
            return s[1:]
        else:
            return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T(2*x*y)
            ...00000000010*x*y
            sage: T(2*x*y)._latex_()
            '...00000000010xy'
        """
        parent = self._parent
        s = ""
        if self._coeff._is_atomic() or (-self._coeff)._is_atomic():
            s = self._coeff._latex_()
            if s == "1": s = ""
        else:
            s = "\\left(%s\\right)" % self._coeff._latex_()
        for i in range(parent._ngens):
            if self._exponent[i] == 1:
                s += "%s" % parent._latex_names[i]
            elif self._exponent[i] > 1:
                s += "%s^{%s}" % (parent._latex_names[i], self._exponent[i])
        if s[0] == "*":
            return s[1:]
        else:
            return s

    def coefficient(self):
        r"""
        Return the coefficient of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2,prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(2*x*y); t
            (2 + O(2^11))*x*y
            sage: t.coefficient()
            2 + O(2^11)
        """
        return self._coeff

    def exponent(self):
        r"""
        Return the exponents of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0)
             over 2-adic Field with capped relative precision 10
            sage: t = T(2,(1,1))
            sage: t.exponent()
            (1, 1)
        """
        return self._exponent

    cpdef _mul_(self, other):
        r"""
        Return the product of this Tate algebra term with ``other``.

        INPUT:

        - ``other`` -- a Tate algebra term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(2*x*y); s
            ...00000000010*x*y
            sage: t = T(3*x^2*y); t
            ...0000000011*x^2*y
            sage: s*t  # indirect doctest
            ...00000000110*x^3*y^2
        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.eadd((<TateAlgebraTerm>other)._exponent)
        ans._coeff = self._coeff * (<TateAlgebraTerm>other)._coeff
        return ans

    #def _div_(self, other):
    #    """
    #    Division of Tate series.
    #
    #    It is currently not implemented because fraction fields
    #    of Tate algebras are not implemented.
    #
    #    EXAMPLES::
    #
    #        sage: A.<x,y> = TateAlgebra(Zp(2))
    #        sage! x / y
    #        Traceback (most recent call last):
    #        ...
    #        NotImplementedError: fraction fields of Tate algebras are not implemented; try inverse_of_unit()
    #
    #        sage: ~x  # indirect doctest
    #        Traceback (most recent call last):
    #        ...
    #        NotImplementedError: fraction fields of Tate algebras are not implemented; try inverse_of_unit()
    #
    #    """
    #    raise NotImplementedError("fraction fields of Tate algebras are not implemented; try inverse_of_unit()")

    cdef long _cmp_c(self, TateAlgebraTerm other) except? 300:
        r"""
        Compare the Tate algebra term with ``other``.

        INPUT:

        - ``other`` -- a term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(x^2*y^2)
            sage: t = T(x*y^3)
            sage: s > t  # indirect doctest
            True
        """
        cdef long c = other._valuation_c() - self._valuation_c()
        if not c:
            skey = self._parent._sortkey
            ks = skey(self._exponent)
            ko = skey(other._exponent)
            c = (ks > ko) - (ks < ko)
        return c

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare the Tate algebra term with ``other`` according to
        the rich comparison operator ``op``.

        The term `a*X^A` is smaller or equal to `b*X^B` if its valuation is
        greater than the valuation of `b*X^B`, or when equality occurs, when
        `A` is smaller or equal to `B` for the Tate algebra monomial order.

        INPUT:

        - ``other`` -- a Tate algebra term

        - ``op`` -- the comparison operator

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T.term_order()
            Degree reverse lexicographic term order
            sage: s = T(x^2*y^2); s
            ...0000000001*x^2*y^2
            sage: t = T(x*y^3); t
            ...0000000001*x*y^3
            sage: s < t  # indirect doctest
            False
            sage: s > t  # indirect doctest
            True
            sage: s <= t  # indirect doctest
            False
            sage: s >= t  # indirect doctest
            True

        Elements with the same valuation and monomial are equivalent
        for the preorder::

            sage: ss = T(3*x^2*y^2); ss
            ...0000000011*x^2*y^2
            sage: s < ss  # indirect doctest
            False
            sage: s > ss  # indirect doctest
            False
            sage: s <= ss  # indirect doctest
            True
            sage: s >= ss  # indirect doctest
            True
            sage: s == ss
            False
        """
        if op == Py_EQ:
            return ((<TateAlgebraTerm>self)._coeff == (<TateAlgebraTerm>other)._coeff
                and (<TateAlgebraTerm>self)._exponent == (<TateAlgebraTerm>other)._exponent)
        if op == Py_NE:
            return ((<TateAlgebraTerm>self)._coeff != (<TateAlgebraTerm>other)._coeff
                 or (<TateAlgebraTerm>self)._exponent != (<TateAlgebraTerm>other)._exponent)
        c = (<TateAlgebraTerm>self)._cmp_c(<TateAlgebraTerm>other)
        return rich_to_bool_sgn(op, c)

    cpdef TateAlgebraTerm monomial(self):
        r"""
        Return this term divided by its coefficient.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(3*x^2*y^2); s
            ...0000000011*x^2*y^2
            sage: s.monomial()
            ...0000000001*x^2*y^2
        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._coeff = self._parent._field(1)
        ans._exponent = self._exponent
        return ans

    cpdef TateAlgebraTerm monic(self):
        r"""
        Return this term normalized so that it has valuation 0
        and its coefficient is a power of the uniformizer.

        EXAMPLES:

        When the log radii of convergence are all zero, the
        coefficient of the returned term is `1`. In this case,
        this method does the same thing as :meth:`monomial`::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(3*x^2*y^2); s
            ...0000000011*x^2*y^2
            sage: s.monic()
            ...0000000001*x^2*y^2
            sage: s.monomial()
            ...0000000001*x^2*y^2

        However, when log radii do not vanish, behaviors might
        be different::

            sage: A.<x,y> = TateAlgebra(R, log_radii=1)
            sage: T = A.monoid_of_terms()
            sage: s = T(3*x^2*y^2); s
            ...0000000011*x^2*y^2
            sage: s.monic()
            ...00000000010000*x^2*y^2
            sage: s.monomial()
            ...0000000001*x^2*y^2

        We compare the valuations::

            sage: s.monic().valuation()
            0
            sage: s.monomial().valuation()
            -4
        """
        cdef TateAlgebraTerm ans = self._new_c()
        cdef long v = self._exponent.dotprod(self._parent._log_radii)
        ans._coeff = self._parent._field.uniformizer_pow(v)
        ans._exponent = self._exponent
        return ans

    def valuation(self):
        r"""
        Return the valuation of this term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            ...000000000100*x^2*y^2
            sage: t.valuation()
            2

        In case of nonzero log radii, the valuations of the variables
        contribute::

            sage: A.<x,y> = TateAlgebra(R, log_radii=1)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            ...000000000100*x^2*y^2
            sage: t.valuation()
            -2
        """
        return ZZ(self._valuation_c())

    cdef long _valuation_c(self) noexcept:
        r"""
        Return the valuation of this term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            ...000000000100*x^2*y^2
            sage: t.valuation()  # indirect doctest
            2
        """
        return (<pAdicGenericElement>self._coeff).valuation_c() - <long>self._exponent.dotprod(self._parent._log_radii)

    cdef Element _call_c(self, list arg):
        """
        Return this term evaluated at ``args``.

        INPUT:

        - ``args`` -- elements

        TESTS::

            sage: R = Zp(2)
            sage: A.<x,y> = TateAlgebra(R)
            sage: M = A.monoid_of_terms()
            sage: t = M(x*y)
            sage: t(0, 1)
            0
        """
        cdef Element ans = self._coeff
        cdef ETuple exponent = self._exponent
        cdef size_t ind
        for ind in range(exponent._nonzero):
            ans *= arg[exponent._data[2*ind]] ** exponent._data[2*ind+1]
        return ans

    def __call__(self, *args):
        """
        Return this term evaluated at ``args``.

        INPUT:

        - ``args`` -- elements

        EXAMPLES::

            sage: R = Zp(2)
            sage: A.<x,y> = TateAlgebra(R)
            sage: M = A.monoid_of_terms()
            sage: t = M(x*y)
            sage: t(1, 2)
            2 + O(2^21)

        An error is raised if we ask for the evaluation at one
        point which is outside the domain of convergence::

            sage: t(1/2, 1)
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence

        TESTS::

            sage: t(1/2, 1, 0)
            Traceback (most recent call last):
            ...
            TypeError: wrong number of arguments

            sage: t(1/2, GF(3)(2))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce all the elements to the same parent
        """
        parent = self._parent
        if len(args) != parent._ngens:
            raise TypeError("wrong number of arguments")
        A = _pushout_family(args, parent._field)
        args = [ A(arg) for arg in args ]
        ratio = A.absolute_e() // parent._base.absolute_e()
        for i in range(parent._ngens):
            if args[i].valuation() < -ratio * parent._log_radii[i]:
                raise ValueError("not in the domain of convergence")
        res = self._call_c(args)
        if parent._integral:
            try:
                res = res.parent().integer_ring()(res)
            except AttributeError:
                pass
        return res

    @coerce_binop
    def is_coprime_with(self, other):
        r"""
        Return ``True`` if this term is coprime with ``other``.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            ...000000000100*x^2*y^2
            sage: s = T(y^3); s
            ...0000000001*y^3
            sage: s.is_coprime_with(t)
            False
            sage: t.is_coprime_with(s)
            False

            sage: tt = T(3*x^2); tt
            ...0000000011*x^2
            sage: s.is_coprime_with(tt)
            True
            sage: tt.is_coprime_with(s)
            True

        When working over a rational Tate algebra, only the
        monomial part of terms are compared::

            sage: t = T(2*x^2); t
            ...00000000010*x^2
            sage: s = T(4*y^3); s
            ...000000000100*y^3
            sage: s.is_coprime_with(t)
            True

        But coefficients play a role when we are working over
        the ring of integers of the Tate Algebra::

            sage: Ao = A.integer_ring()
            sage: To = Ao.monoid_of_terms()
            sage: To(s).is_coprime_with(To(t))
            False
        """
        for i in range(self._parent.ngens()):
            if self._exponent[i] > 0 and other.exponent()[i] > 0:
                return False
        if self._parent.base_ring().is_field():
            return True
        else:
            return self.valuation() == 0 or other.valuation() == 0

    @coerce_binop
    def gcd(self, other):
        r"""
        Return the greatest common divisor of this term and ``other``.

        The result is normalized so that:

        - its valuation is equal to the smallest valuation of
          this term and ``other``

        - its coefficient is a power of the uniformizer.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            ...0000000001000*x^2*y^2
            sage: t = T(4*x*y^3); t
            ...000000000100*x*y^3
            sage: s.gcd(t)
            ...000000000100*x*y^2
        """
        return self._gcd_c(other)

    cdef TateAlgebraTerm _gcd_c(self, TateAlgebraTerm other):
        r"""
        Return the greatest common divisor of this term and ``other``.

        The result is normalized so that:

        - its valuation is equal to the smallest valuation of
          this term and ``other``

        - its coefficient is a power of the uniformizer.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            ...0000000001000*x^2*y^2
            sage: t = T(4*x*y^3); t
            ...000000000100*x*y^3
            sage: s.gcd(t)  # indirect doctest
            ...000000000100*x*y^2

        ::

            sage: A.<x,y> = TateAlgebra(R, log_radii=1)
            sage: T = A.monoid_of_terms()
            sage: T(x^5).gcd(T(y^5))
            ...00000.00001
        """
        cdef TateAlgebraTerm ans = self._new_c()
        cdef long val
        ans._exponent = self._exponent.emin(other._exponent)
        val = min(self._valuation_c(), other._valuation_c()) + ans._exponent.dotprod(self._parent._log_radii)
        ans._coeff = self._parent._field.uniformizer_pow(val)
        return ans

    @coerce_binop
    def lcm(self, other):
        r"""
        Return the least common multiple of two Tate terms.

        The result is normalized so that `\gcd(a,b) \lcm(a,b) = ab`.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

        In a Tate algebra over a field:

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            ...0000000001000*x^2*y^2
            sage: t = T(4*x*y^3); t
            ...000000000100*x*y^3
            sage: s.lcm(t)
            ...0000000001000*x^2*y^3
        """
        return self._lcm_c(other)

    cdef TateAlgebraTerm _lcm_c(self, TateAlgebraTerm other):
        r"""
        Return the least common multiple of two Tate terms.

        The result is normalized so that `\gcd(a,b) \lcm(a,b) = ab`.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

        In a Tate algebra over a field:

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            ...0000000001000*x^2*y^2
            sage: t = T(12*x*y^3); t
            ...000000001100*x*y^3
            sage: s.lcm(t)  # indirect doctest
            ...0000000011000*x^2*y^3

        TESTS::

            sage: s.gcd(t) * s.lcm(t) == s * t
            True
        """
        cdef TateAlgebraTerm ans = self._new_c()
        cdef long val
        ans._exponent = self._exponent.emax(other._exponent)
        val = max(self._valuation_c(), other._valuation_c()) + ans._exponent.dotprod(self._parent._log_radii)
        ans._coeff = (self._coeff.unit_part() * other._coeff.unit_part()) << val
        return ans

    @coerce_binop
    def is_divisible_by(self, other, integral=False):
        r"""
        Return ``True`` if this term is divisible by ``other``.

        INPUT:

        - ``other`` -- a Tate term

        - ``integral`` -- (default: ``False``) if ``True``, test
          for divisibility in the ring of integers of the Tate algebra

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(4*x^2*y^2); s
            ...000000000100*x^2*y^2
            sage: t = T(4*x*y^3); t
            ...000000000100*x*y^3
            sage: s.is_divisible_by(t)
            False

            sage: t = T(4*x*y^2); t
            ...000000000100*x*y^2
            sage: s.is_divisible_by(t)
            True

            sage: t = T(16); t
            ...00000000010000
            sage: s.is_divisible_by(t)
            True
            sage: s.is_divisible_by(t, integral=True)
            False

        If you are working over the ring of integers of the Tate algebra,
        divisibility is always checked in the ring of integers (even if
        ``integral`` is set to ``False``)::

            sage: Ao = A.integer_ring()
            sage: To = Ao.monoid_of_terms()
            sage: so = To(s)
            sage: to = To(t)
            sage: so.is_divisible_by(to)
            False
            sage: so.is_divisible_by(to, integral=False)
            False

        Be careful that coercion between the Tate algebra and its ring of
        integers can be done silently::

            sage: s.is_divisible_by(to)
            True
        """
        return (<TateAlgebraTerm?>other)._divides_c(self, integral)

    @coerce_binop
    def divides(self, other, integral=False):
        r"""
        Return ``True`` if this term divides ``other``.

        INPUT:

        - ``other`` -- a Tate term

        - ``integral`` -- (default: ``False``) if ``True``, test for
          divisibility in the ring of integers of the Tate algebra

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(4*x^2*y^2); s
            ...000000000100*x^2*y^2
            sage: t = T(4*x*y^3); t
            ...000000000100*x*y^3
            sage: t.divides(s)
            False

            sage: t = T(4*x*y^2); t
            ...000000000100*x*y^2
            sage: t.divides(s)
            True

            sage: t = T(16); t
            ...00000000010000
            sage: t.divides(s)
            True
            sage: t.divides(s, integral=True)
            False

        If you are working over the ring of integers of the Tate algebra,
        divisibility is always checked in the ring of integers (even if
        ``integral`` is set to ``False``)::

            sage: Ao = A.integer_ring()
            sage: To = Ao.monoid_of_terms()
            sage: so = To(s)
            sage: to = To(t)
            sage: to.divides(so)
            False
            sage: to.divides(so, integral=False)
            False

        Be careful that coercion between the Tate algebra and its ring of
        integers can be done silently::

            sage: to.divides(s)
            True
        """
        return self._divides_c(other, integral)

    cdef bint _divides_c(self, TateAlgebraTerm other, bint integral) noexcept:
        r"""
        Return ``True`` if this term divides ``other``.

        INPUT:

        - ``other`` -- a Tate term

        - ``integral`` -- boolean (default: ``False``); if ``True``, test for
          divisibility in the ring of integers of the Tate algebra

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(4*x^2*y^2); s
            ...000000000100*x^2*y^2
            sage: t = T(4*x*y^3); t
            ...000000000100*x*y^3
            sage: t.divides(s)  # indirect doctest
            False
        """
        parent = self._parent
        if (integral or not parent.base_ring().is_field()) and self.valuation() > other.valuation():
            return False
        for i in range(parent._ngens):
            if self._exponent[i] > other._exponent[i]:
                return False
        return True

    cpdef _floordiv_(self, other):
        r"""
        Return the result of the exact division of this term by ``other``.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(2*x^2*y^3); t
            ...00000000010*x^2*y^3
            sage: s = T(6*x*y^2); s
            ...00000000110*x*y^2
            sage: t // s
            ...1010101011*x*y

        If the Tate terms are not divisible, an error is raised::

            sage: s // t
            Traceback (most recent call last):
            ...
            ValueError: the division is not exact
        """
        if not self.is_divisible_by(other):
            raise ValueError("the division is not exact")
        return (<TateAlgebraTerm>self)._floordiv_c(<TateAlgebraTerm>other)

    cdef TateAlgebraTerm _floordiv_c(self, TateAlgebraTerm other):
        r"""
        Return the result of the exact division of this term by ``other``.

        INPUT:

        - ``other`` -- a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(2*x^2*y^3); t
            ...00000000010*x^2*y^3
            sage: s = T(6*x*y^2); s
            ...00000000110*x*y^2
            sage: t // s
            ...1010101011*x*y

        If the Tate terms are not divisible, an error is raised::

            sage: s // t
            Traceback (most recent call last):
            ...
            ValueError: the division is not exact
        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self.exponent().esub(other.exponent())
        ans._coeff = self._coeff / other._coeff
        return ans


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    r"""
    A class for Tate series, elements of Tate algebras.

    EXAMPLES::

        sage: R = Zp(2, prec=10, print_mode='digits')
        sage: A.<x,y> = TateAlgebra(R)
        sage: A(2*x+1)
        ...0000000001 + ...00000000010*x
        sage: A(2*x+1, prec=5)
        ...00001 + ...00010*x + O(2^5 * <x, y>)
        sage: A(2*x+1, prec=20)
        ...0000000001 + ...00000000010*x + O(2^20 * <x, y>)
    """
    def __init__(self, parent, x, prec=None, reduce=True):
        r"""
        Initialize a Tate algebra series.

        TESTS::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: TestSuite(x).run()
        """
        cdef TateAlgebraElement xc
        CommutativeAlgebraElement.__init__(self, parent)
        self._prec = Infinity  # TODO: replace infinity by a big int
        if isinstance(x, TateAlgebraElement):
            xc = <TateAlgebraElement>x
            xparent = x.parent()
            if xparent is parent:
                self._poly = PolyDict(xc._poly.__repn, None)
                self._prec = xc._prec
            elif xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("cannot restrict to a bigger domain")
                self._poly = PolyDict({ e: parent._field(v) for (e,v) in xc._poly.__repn.items() }, None)
                if xc._prec is not Infinity:
                    self._prec = (xc._prec * ratio).ceil()
            else:
                raise TypeError("variable names do not match")
        elif isinstance(x, TateAlgebraTerm):
            xparent = x.parent()
            if xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("cannot restrict to a bigger domain")
                self._poly = PolyDict({(<TateAlgebraTerm>x)._exponent: parent._field((<TateAlgebraTerm>x)._coeff)}, None)
        else:
            try:
                poly = parent._polynomial_ring(x)
                self._poly = PolyDict(poly.monomial_coefficients(), None)
            except TypeError:
                # last chance: we first try to convert to the rational Tate series
                if parent._integral:
                    xc = parent._rational_ring(x)
                    self._poly = xc._poly
                    self._prec = xc._prec
                else:
                    raise
        if prec is not None:
            self._prec = min(self._prec, prec)
        self._normalize()
        self._terms = self._terms_nonzero = None
        if not parent.base_ring().is_field() and self.valuation() < 0:
            raise ValueError("this series is not in the ring of integers")

    cdef TateAlgebraElement _new_c(self):
        """
        Fast creation of a new Tate series.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: x + y  # indirect doctest
            ...0000000001*x + ...0000000001*y
        """
        cdef TateAlgebraElement ans = TateAlgebraElement.__new__(TateAlgebraElement)
        ans._parent = self._parent
        ans._is_normalized = False
        ans._terms = ans._terms_nonzero = None
        return ans

    cdef _normalize(self):
        """
        Normalize this series.

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: A(78612, prec=3)  # indirect doctest
            ...100 + O(2^3 * <x, y>)
        """
        self._is_normalized = True
        if self._prec is Infinity:
            return
        cdef int v
        cdef pAdicGenericElement coeff
        for (e, c) in list(self._poly.__repn.items()):
            v = (<ETuple>self._parent._log_radii).dotprod(<ETuple>e)
            coeff = self._poly.__repn[e]
            if coeff.precision_absolute() > self._prec - v:
                coeff = coeff.add_bigoh(self._prec - v)
            if coeff.valuation() >= self._prec - v:
                del self._poly.__repn[e]
            else:
                self._poly.__repn[e] = coeff

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: A.<x,y> = TateAlgebra(Zp(2))
            sage: loads(dumps(x)) == x  # indirect doctest
            True
        """
        return TateAlgebraElement, (self.parent(), self._poly, self._prec)

    def _repr_(self):
        r"""
        Return a string representation of this series.

        The terms are ordered with decreasing term order
        (increasing valuation, then the monomial order of the parent algebra).

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: x + 2*x^2 + x^3
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2

            sage: A(x + 2*x^2 + x^3, prec=5)
            ...00001*x^3 + ...00001*x + ...00010*x^2 + O(2^5 * <x, y>)
        """
        base = self._parent.base_ring()
        nvars = self._parent.ngens()
        vars = self._parent.variable_names()
        s = ""
        for t in self._terms_c():
            if t.valuation() >= self._prec:
                continue
            st = repr(t)
            if s == "":
                s += st
            elif st[0] == "-":
                s += " - " + st[1:]
            else:
                s += " + " + st
        if self._prec is not Infinity:
            if s != "":
                s += " + "
            su = self._parent._uniformizer_repr
            lr = self._parent.log_radii()
            sv = [ ]
            for i in range(len(vars)):
                if lr[i] == 0:
                    sv.append(vars[i])
                elif lr[i] == -1:
                    sv.append("%s*%s" % (su, vars[i]))
                elif lr[i] == 1:
                    sv.append("%s/%s" % (vars[i], su))
                elif lr[i] < 0:
                    sv.append("%s^%s*%s" % (su, -lr[i], vars[i]))
                else:
                    sv.append("%s/%s^%s" % (vars[i], su, lr[i]))
            sv = ", ".join(sv)
            if self._prec == 0:
                s += "O(<%s>)" % sv
            elif self._prec == 1:
                s += "O(%s * <%s>)" % (self._parent._uniformizer_repr, sv)
            else:
                s += "O(%s^%s * <%s>)" % (self._parent._uniformizer_repr, self._prec, sv)
        if s == "":
            return "0"
        return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this series.

        The terms are ordered with decreasing term order
        (increasing valuation of the coefficients, then
        the monomial order of the parent algebra).

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: f._latex_()
            '...0000000001x^{3} + ...0000000001x + ...00000000010x^{2}'
        """
        base = self._parent.base_ring()
        nvars = self._parent.ngens()
        vars = self._parent.variable_names()
        s = ""
        for t in self.terms():
            if t.valuation() >= self._prec:
                continue
            st = t._latex_()
            if s == "":
                s += st
            elif st[0] == "-":
                s += " - " + st[1:]
            else:
                s += " + " + st
        if self._prec is not Infinity:
            if s != "":
                s += " + "
            sv = ",".join(vars)
            if self._prec == 0:
                s += "O\\left(%s\\right)" % self._parent.integer_ring()._latex_()
            elif self._prec == 1:
                s += "O\\left(%s %s\\right)" % (self._parent._uniformizer_latex, self._parent.integer_ring()._latex_())
            else:
                s += "O\\left(%s^{%s} %s\\right)" % (self._parent._uniformizer_latex, self._prec, self._parent.integer_ring()._latex_())
        return s

    cpdef _add_(self, other):
        r"""
        Return the sum of this series and ``other``.

        The precision of the output is adjusted to the minimum of both precisions.

        INPUT:

        - ``other`` -- a Tate series

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: g = A(x + 2*y^2, prec=5); g
            ...00001*x + ...00010*y^2 + O(2^5 * <x, y>)
            sage: h = f + g; h  # indirect doctest
            ...00001*x^3 + ...00010*x^2 + ...00010*y^2 + ...00010*x + O(2^5 * <x, y>)

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            5
        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly + (<TateAlgebraElement>other)._poly
        ans._prec = min(self._prec, (<TateAlgebraElement>other)._prec)
        ans._normalize()
        return ans

    cpdef _neg_(self):
        r"""
        Return the opposite of this series.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: -f  # indirect doctest
            ...1111111111*x^3 + ...1111111111*x + ...11111111110*x^2
        """
        cdef TateAlgebraElement ans = self._new_c()
        cdef Element s = self._parent.base_ring()(-1)
        ans._poly = self._poly.scalar_lmult(s)
        ans._prec = self._prec
        return ans

    cpdef _sub_(self, other):
        r"""
        Return the difference of this series and ``other``.

        The precision of the output is adjusted to the minimum of both precisions.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: g = A(x + 2*y^2, prec=5); g
            ...00001*x + ...00010*y^2 + O(2^5 * <x, y>)
            sage: h = f - g; h # indirect doctest
            ...00001*x^3 + ...00010*x^2 + ...11110*y^2 + O(2^5 * <x, y>)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            5
        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly - (<TateAlgebraElement>other)._poly
        ans._prec = min(self._prec, (<TateAlgebraElement>other)._prec)
        ans._normalize()
        return ans

    cpdef _mul_(self, other):
        r"""
        Return the product of this series with ``other``.

        The precision is adjusted to match the best precision on the output.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x + 4*x^2 + 2*x^3; f
            ...00000000010*x^3 + ...00000000010*x + ...000000000100*x^2
            sage: g = A(x + 2*x^2, prec=5); g
            ...00001*x + ...00010*x^2 + O(2^5 * <x, y>)
            sage: h = f * g; h # indirect doctest
            ...001010*x^4 + ...000010*x^2 + ...000100*x^5 + ...001000*x^3 + O(2^6 * <x, y>)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            6
        """
        cdef TateAlgebraElement ans = self._new_c()
        a = self._prec + (<TateAlgebraElement>other).valuation()
        b = self.valuation() + (<TateAlgebraElement>other)._prec
        ans._poly = self._poly * (<TateAlgebraElement>other)._poly
        ans._prec = min(a, b)
        ans._normalize()
        return ans

    cpdef _lmul_(self, Element right):
        r"""
        Return the product of this series by ``right``.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: 2*f # indirect doctest
            ...00000000010*x^3 + ...00000000010*x + ...000000000100*x^2

            sage: 6*f # indirect doctest
            ...00000000110*x^3 + ...00000000110*x + ...000000001100*x^2
        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly.scalar_lmult(right)
        ans._prec = self._prec + (<pAdicGenericElement>self._parent._base(right)).valuation_c()
        return ans

    def inverse_of_unit(self, prec=None):
        r"""
        Return the inverse of this series if it is invertible.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``);
          the precision at which the result is computed, if ``None``,
          the result is truncated according to the cap of the parent

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = A(1); f
            ...0000000001
            sage: f.inverse_of_unit()
            ...0000000001 + O(2^10 * <x, y>)

            sage: f = 2*x + 1; f
            ...0000000001 + ...00000000010*x
            sage: f.inverse_of_unit()
            ...0000000001 + ...1111111110*x + ...0000000100*x^2 + ...1111111000*x^3
             + ...0000010000*x^4 + ...1111100000*x^5 + ...0001000000*x^6
             + ...1110000000*x^7 + ...0100000000*x^8 + ...1000000000*x^9 + O(2^10 * <x, y>)

            sage: f.inverse_of_unit(prec=4)
            ...0001 + ...1110*x + ...0100*x^2 + ...1000*x^3 + O(2^4 * <x, y>)

        If the series is not invertible, an error is raised::

            sage: f = 1 + x; f
            ...0000000001*x + ...0000000001
            sage: f.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ValueError: this series in not invertible
        """
        cdef TateAlgebraTerm t
        cdef long v, curprec
        cdef pAdicGenericElement c
        cdef TateAlgebraElement inv, x
        if not self.is_unit():
            raise ValueError("this series in not invertible")
        t = self.leading_term()
        c = t.coefficient()
        parent = self._parent
        v, c = c.val_unit()
        x = self >> v

        rprec = self.precision_relative()
        if prec is not None and prec + v < rprec:
            rprec = prec + v
        if rprec is Infinity:
            rprec = self._parent.precision_cap()

        x = x.add_bigoh(rprec)
        inv = (<TateAlgebraElement>self)._new_c()
        inv._poly = PolyDict({ ETuple({}, parent._ngens): parent._field(~c.residue()).lift_to_precision() }, None)
        inv._prec = rprec
        curprec = 1
        while curprec < rprec:
            curprec *= 2
            inv = 2*inv - x*inv*inv
        return inv >> v

    def is_unit(self):
        r"""
        Return ``True`` if this series is invertible.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x + 1; f
            ...0000000001 + ...00000000010*x
            sage: f.is_unit()
            True

            sage: f = 1 + x; f
            ...0000000001*x + ...0000000001
            sage: f.is_unit()
            False

        Note that invertibility is tested in the parent of this series::

            sage: f = 4*x + 2
            sage: f.is_unit()
            True

            sage: Ao = A.integer_ring()
            sage: Ao(f).is_unit()
            False
        """
        if self.is_zero():
            return False
        t = self.leading_term()
        if max(t.exponent()) != 0:
            return False
        base = self.base_ring()
        if not base.is_field() and t.valuation() > 0:
            return False
        return True

    def __pow__(self, exponent, modulus):
        r"""
        Return this element raised to the power ``exponent``.

        INPUT:

        - ``exponent`` -- either an integer, a rational number or a
          Tate series

        - ``modulus`` -- discarded

        EXAMPLES::

            sage: R = Zp(3, prec=4, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: (x + y)^3
            ...0001*x^3 + ...0001*y^3 + ...0010*x^2*y + ...0010*x*y^2

        This function can be used to compute the inverse of a Tate series::

            sage: f = 1 + 6*x^2 + 9*y^2
            sage: f^(-1)
            ...0001 + ...2210*x^2 + ...1100*x^4 + ...2200*y^2 + ...1000*x^6
             + ...1000*x^2*y^2 + O(3^4 * <x, y>)

        or a square root (or more generally an n-th root)::

            sage: g = f^(1/2); g
            ...0001 + ...0010*x^2 + ...1100*x^4 + ...1200*y^2 + ...2000*x^6
             + ...1000*x^2*y^2 + O(3^4 * <x, y>)
            sage: g^2 == f
            True

        When the exponent is not an integer, `f^e` is computed as `\exp(e \log(f))`.
        This computation fails if `f` is outside the domain of the logarithm or if
        `e \log(f)` is outside the domain of convergence of the exponential::

            sage: f^(1/3)
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence

        The exponent can be a series as well::

            sage: g = f^x; g
            ...0001 + ...0020*x^3 + ...1100*x^9 + ...2200*x^7 + ... + O(3^4 * <x, y>)

            sage: x0 = R.random_element()
            sage: y0 = R.random_element()
            sage: g(x0, y0) == f(x0, y0)^x0
            True

        TESTS::

            sage: f^(x + y) == f^x * f^y
            True
            sage: f^(x*y) == (f^x)^y
            True
            sage: f^(x*y) == (f^y)^x
            True
        """
        cdef TateAlgebraElement temp
        cdef long p, v, e
        if exponent in ZZ:
            if exponent == 0:
                return (<TateAlgebraElement>self)._parent.one()
            elif exponent > 0:
                # We handle precision
                temp = (<TateAlgebraElement>self)._new_c()
                temp._poly = (<TateAlgebraElement>self)._poly
                if (<TateAlgebraElement>self)._prec is Infinity:
                    temp._prec = Infinity
                else:
                    p = (<TateAlgebraElement>self)._parent.prime()
                    e = (<TateAlgebraElement>self)._parent.absolute_e()
                    v = (<TateAlgebraElement>self).valuation()
                    temp._prec = p*v + (<TateAlgebraElement>self)._prec
                    q = ZZ(exponent)
                    while True:
                        q, r = q.quo_rem(p)
                        if r != 0: break
                        temp._prec = min(p * temp._prec, temp._prec + e)
                # and perform the exponentiation
                return CommutativeAlgebraElement.__pow__(temp, exponent, None)
            else:
                return self.inverse_of_unit() ** (-exponent)
        elif exponent in QQ:
            exponent = QQ(exponent)
            return self.nth_root(exponent.denominator()) ** exponent.numerator()
        else:
            return ((<TateAlgebraElement>self).log() * exponent).exp()

    def square_root(self, prec=None):
        r"""
        Return the square root of this series.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``);
          the precision at which the result is computed, if ``None``,
          the result is truncated according to the cap of the parent

        EXAMPLES::

            sage: R = Zp(3, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 6*x^2 + 9*y^2
            sage: g = f.sqrt(); g
            ...0000000001 + ...0000000010*x^2 + ...1111111100*x^4 + ...1111111200*y^2
             + ...1111112000*x^6 + ...1111111000*x^2*y^2 + ... + O(3^10 * <x, y>)

            sage: f.square_root(prec=4)
            ...0001 + ...0010*x^2 + ...1100*x^4 + ...1200*y^2 + ...2000*x^6
             + ...1000*x^2*y^2 + O(3^4 * <x, y>)

            sage: g^2 == f
            True

        It's possible that `f` has a trivial square root (which is analytic on
        the correct domain) but that it takes its values outside the domain of
        convergence of the square root function.
        In this case, an error is raised::

            sage: f = x^2
            sage: f.square_root()
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence
        """
        return self.nth_root(2, prec)

    def sqrt(self, prec=None):
        r"""
        Return the square root of this series.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``);
          the precision at which the result is computed, if ``None``,
          the result is truncated according to the cap of the parent

        EXAMPLES::

            sage: R = Zp(3, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 6*x^2 + 9*y^2
            sage: g = f.sqrt(); g
            ...0000000001 + ...0000000010*x^2 + ...1111111100*x^4 + ...1111111200*y^2
             + ...1111112000*x^6 + ...1111111000*x^2*y^2 + ... + O(3^10 * <x, y>)

            sage: f.sqrt(prec=4)
            ...0001 + ...0010*x^2 + ...1100*x^4 + ...1200*y^2 + ...2000*x^6
             + ...1000*x^2*y^2 + O(3^4 * <x, y>)

            sage: g^2 == f
            True

        It's possible that `f` has a trivial square root (which is analytic on
        the correct domain) but that it takes its values outside the domain of
        convergence of the square root function.
        In this case, an error is raised::

            sage: f = x^2
            sage: f.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence
        """
        return self.nth_root(2, prec)

    def nth_root(self, n=2, prec=None):
        r"""
        Return the `n`-th root of this series.

        INPUT:

        - ``n`` -- integer (default: `2`)

        - ``prec`` -- integer or ``None`` (default: ``None``);
          the precision at which the result is computed, if ``None``,
          the result is truncated according to the cap of the parent

        .. NOTE::

            The `n`-th root is computed as `\exp(\frac 1 n \log(f))`.

        EXAMPLES::

            sage: R = Zp(3, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 9*x^2 + 9*y^2
            sage: g = f.nth_root(3, prec=3); g
            ...001 + ...010*x^2 + ...010*y^2 + ...200*x^6 + ...200*y^6 + ...200*x^4
             + ...100*x^2*y^2 + ...200*y^4 + O(3^3 * <x, y>)
            sage: g^3 == f
            True

            sage: for n in range(2, 9):
            ....:     if f.nth_root(n)^n != f: raise RuntimeError

        It's possible that `f` has a trivial ``n``-th root (which is analytic on
        the correct domain) but that `\exp(\frac 1 n \log(f))` does not converge.
        In this case, an error is raised::

            sage: f = x^3
            sage: f.nth_root(3)
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence
        """
        if n not in ZZ or n == 0:
            raise ValueError("n must be a nonzero integer")
        n = ZZ(n)

        cdef TateAlgebraElement a, root
        cdef long v, e, ep, curprec
        cdef pAdicGenericElement scalar

        parent = self._parent
        p = parent.prime()
        v = n.valuation(p)
        e = parent.absolute_e()
        ep = e // (p - 1)
        if v == 0:
            curprec = 1
        else:
            curprec = v*e + ep + 1
        if (self - 1).valuation() < curprec:
            raise ValueError("not in the domain of convergence")

        aprec = self.precision_absolute()
        if prec is not None and prec < aprec:
            aprec = prec
        if aprec is Infinity:
            aprec = parent.precision_cap()

        aprec += v*e
        root = parent(1).add_bigoh(v*e+1)
        if n > 0:
            a = self.inverse_of_unit(aprec)
        else:
            n = -n
            a = self.add_bigoh(aprec)
        scalar = ~(parent._field(n))
        while curprec < aprec:
            if v == 0:
                curprec *= 2
            else:
                curprec -= v*e
                curprec = min(2*curprec + v*e, p*curprec + (v-1)*e)
            root = root.lift_to_precision(min(aprec, curprec))
            root += scalar * root * (1 - a * root**n)
        return root

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare this series with ``other`` according to
        the rich comparison operator ``op``.

        INPUT:

        - ``other`` -- a Tate algebra element

        - ``op`` -- the comparison operator

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: A.term_order()
            Degree reverse lexicographic term order

        For terms, we first compare the valuation and ties are broken
        using the term ordering on `A`::

            sage: 2*x^2*y^2 > x*y^3
            False
            sage: x^2*y^2 > x*y^3
            True

        For general series, leading terms are first compared. In case
        of tie, second leading terms are compared, and so on::

            sage: x^2*y^2 > x*y^3 + y^4
            True
            sage: x^2*y^2 > x^2*y^2 + x*y^3
            False

        TESTS::

            sage: f = x^4 + 4*x*y + 1; f
            ...0000000001*x^4 + ...0000000001 + ...000000000100*x*y
            sage: g = f + 2
            sage: f == g
            False
            sage: f == g - 2
            True
        """
        diff = self - other
        c = None
        if diff.is_zero():
            c = 0
        if op == Py_EQ:
            return c is not None
        if op == Py_NE:
            return c is None
        if c is None:
            ts = self.terms()
            to = other.terms()
            for s, o in zip(ts, to):
                c = (<TateAlgebraTerm>s)._cmp_c(<TateAlgebraTerm>o)
                if c: break
            else:
                c = len(ts) - len(to)
        return rich_to_bool_sgn(op, c)

    def __call__(self, *args):
        """
        Return this term evaluated at ``args``.

        INPUT:

        - ``args`` -- elements

        EXAMPLES::

            sage: R = Zp(2)
            sage: A.<u,v> = TateAlgebra(R, log_radii=[0,-1])
            sage: A
            Tate Algebra in u (val >= 0), v (val >= 1) over 2-adic Field with capped relative precision 20

            sage: f = u^2 + v^2
            sage: f(1, 0)
            1 + O(2^20)

        An error is raised if we ask for the evaluation at one
        point which is outside the domain of convergence::

            sage: f(1, 1)
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence

        Evaluation at points in extensions is allowed::

            sage: S.<a> = Zq(2^3)
            sage: f(a, 2)
            a^2 + 2^2 + O(2^20)

            sage: x = polygen(ZZ, 'x')
            sage: T.<pi> = S.extension(x^2 - 2)
            sage: f(pi, 2)
            pi^2 + pi^4 + O(pi^42)

            sage: f(pi, pi)
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence

        This method can also be used to compose Tate series::

            sage: f(u + v, 2*u)
            (1 + 2^2 + O(2^20))*u^2 + (2 + O(2^20))*u*v + (1 + O(2^20))*v^2

        or for partial evaluation::

            sage: f(pi, v)
            (pi^2 + O(pi^42)) + (1 + O(pi^40))*v^2
        """
        cdef TateAlgebraTerm t
        parent = self._parent
        if len(args) != parent._ngens:
            raise TypeError("wrong number of arguments")
        A = _pushout_family(args, parent._field)
        args = [ A(arg) for arg in args ]
        ratio = A.absolute_e() // parent._base.absolute_e()
        for i in range(parent._ngens):
            if args[i].valuation() < -ratio * parent._log_radii[i]:
                raise ValueError("not in the domain of convergence")
        res = A(0, ratio * self._prec)
        for t in self._terms_c():
            if t._valuation_c() >= res.precision_absolute():
                break
            res += t._call_c(args)
        if parent._integral:
            try:
                res = res.parent().integer_ring()(res)
            except AttributeError:
                pass
        return res

    cdef TateAlgebraElement _term_mul_c(self, TateAlgebraTerm term):
        r"""
        Return the product of this series by the term ``term``.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: t = A.monoid_of_terms()(3*x^2); t
            ...0000000011*x^2
            sage: f = x^4 + 4*x*y + 1; f
            ...0000000001*x^4 + ...0000000001 + ...000000000100*x*y
            sage: t*f  # indirect doctest
            ...0000000011*x^6 + ...0000000011*x^2 + ...000000001100*x^3*y
        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly.term_lmult(term._exponent, term._coeff)
        ans._prec = self._prec + term._valuation_c()
        return ans

    cdef TateAlgebraElement _positive_lshift_c(self, n):
        r"""
        Return the product of this series by the ``n``-th power
        of the uniformizer.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: f << 2  # indirect doctest
            ...000000000100*x^3 + ...000000000100*x + ...0000000001000*x^2
        """
        cdef dict coeffs = { }
        cdef ETuple e
        cdef Element c
        cdef TateAlgebraElement ans = self._new_c()
        for (e,c) in self._poly.__repn.items():
            coeffs[e] = c << n
        ans._poly = PolyDict(coeffs, None)
        ans._prec = self._prec + n
        return ans

    cdef TateAlgebraElement _lshift_c(self, n):
        r"""
        Return the product of this series by the ``n``-th power
        of the uniformizer.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: f << 2  # indirect doctest
            ...000000000100*x^3 + ...000000000100*x + ...0000000001000*x^2
            sage: Ao = A.integer_ring()
            sage: g = Ao(f).add_bigoh(5); g
            ...00001*x^3 + ...00001*x + ...00010*x^2 + O(2^5 * <x, y>)
            sage: g << 2
            ...0000100*x^3 + ...0000100*x + ...0001000*x^2 + O(2^7 * <x, y>)
        """
        cdef dict coeffs = { }
        cdef ETuple e
        cdef Element c
        cdef TateAlgebraElement ans = self._new_c()
        parent = self._parent
        base = parent.base_ring()
        if base.is_field():
            for (e,c) in self._poly.__repn.items():
                coeffs[e] = c << n
            ans._prec = self._prec + n
        else:
            field = base.fraction_field()
            ngens = parent.ngens()
            for (e,c) in self._poly.__repn.items():
                minval = ZZ(e.dotprod(<ETuple>parent._log_radii)).ceil()
                coeffs[e] = field(base(c) >> (minval-n)) << minval
            ans._prec = max(ZZ(0), self._prec + n)
        ans._poly = PolyDict(coeffs, None)
        return ans

    def __lshift__(self, n):
        r"""
        Return the product of this series by the ``n``-th power
        of the uniformizer.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: f << 2  # indirect doctest
            ...000000000100*x^3 + ...000000000100*x + ...0000000001000*x^2
            sage: f << -1  # indirect doctest
            ...000000000.1*x^3 + ...000000000.1*x + ...0000000001*x^2

        If we're shifting by a negative number of digits over the ring of
        integers of a Tate algebra, the result is truncated -- that is, the
        output is the result of the integer division of the Tate series by
        `\pi^{-n}` where `\pi` is a uniformizer.

            sage: Ao = A.integer_ring()
            sage: Ao(f) << -1
            ...0000000001*x^2 + ...000000000*x^3 + ...000000000*x
        """
        return (<TateAlgebraElement>self)._lshift_c(n)

    def __rshift__(self, n):
        r"""
        Return the quotient in the division of this series by
        the ``n``-th power of the uniformizer.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: f << 2
            ...000000000100*x^3 + ...000000000100*x + ...0000000001000*x^2
            sage: f << -1  # indirect doctest
            ...000000000.1*x^3 + ...000000000.1*x + ...0000000001*x^2

        If we're working over the ring of integers of a Tate algebra, the
        result is truncated -- that is, the output is the result of the integer
        division of the Tate series by `\pi^n` where `\pi` is a uniformizer.

            sage: Ao = A.integer_ring()
            sage: Ao(f) << -1
            ...0000000001*x^2 + ...000000000*x^3 + ...000000000*x
        """
        return (<TateAlgebraElement>self)._lshift_c(-n)

    def __bool__(self):
        r"""
        Return ``True`` if this term is nonzero, ``False`` otherwise.

        TESTS::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + y^3
            sage: bool(f)
            True
            sage: bool(f-f)
            False
        """
        cdef list terms = self._terms_c(include_zero=False)
        return bool(terms) and (<TateAlgebraTerm>terms[0])._valuation_c() < self._prec

    def is_zero(self, prec=None):
        r"""
        Return ``True`` if this series is indistinguishable from zero.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``),
          the precision at which the series should be compared to zero

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            ...0000000001*x^3 + ...0000000001*x + ...00000000010*x^2
            sage: f.is_zero()
            False

            sage: g = f << 4; g
            ...00000000010000*x^3 + ...00000000010000*x + ...000000000100000*x^2
            sage: g.is_zero()
            False
            sage: g.is_zero(5)
            False
            sage: g.is_zero(4)
            True
        """
        cdef list terms = self._terms_c(include_zero=False)
        if prec is None:
            prec = self._prec
        if not terms:
            return True
        return (<TateAlgebraTerm>terms[0])._valuation_c() >= prec

    def restriction(self, log_radii):
        r"""
        Return the restriction of this series to a smaller domain.

        INPUT:

        - ``log_radii`` -- integer or a tuple; the log-radii of
          convergence of the smaller domain (see :class:`TateAlgebra`
          for more details)

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x + y^2; f
            ...0000000001*y^2 + ...00000000010*x

            sage: g = f.restriction(-1); g
            ...0000000001*y^2 + ...00000000010*x
            sage: g.parent()
            Tate Algebra in x (val >= 1), y (val >= 1)
             over 2-adic Field with capped relative precision 10

        Note that restricting may change the order of the terms::

            sage: f.restriction([-1,-2])
            ...00000000010*x + ...0000000001*y^2
        """
        parent = self._parent
        from sage.rings.tate_algebra import TateAlgebra
        ring = TateAlgebra(self.base_ring(), names=parent.variable_names(), log_radii=log_radii,
                           prec=parent.precision_cap(), order=parent.term_order())
        return ring(self)

    def terms(self):
        r"""
        Return a list of the terms of this series sorted in descending order.

        .. NOTE::

            The order on the terms is defined as follows: first we
            compare the valuation and, in case of equality, we compare
            the monomials with respect to the order given at the
            creation of the parent.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.terms()
            [...0000000001*x, ...00000000010*x^2]
        """
        if not self._is_normalized:
            self._normalize()
            self._terms = None
        return self._terms_c()

    cdef list _terms_c(self, bint include_zero=True):
        r"""
        Return a list of the terms of this series sorted in descending order.

        INPUT:

        - ``include_zero`` -- boolean (default: ``True``); if ``True``,
          include terms which are indistinguishable from zero

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.terms()  # indirect doctest
            [...0000000001*x, ...00000000010*x^2]
        """
        cdef pAdicGenericElement c
        cdef ETuple e
        cdef TateAlgebraTerm oneterm = self._parent._oneterm
        cdef TateAlgebraTerm term
        if self._terms is None:
            self._terms = []
            for (e,c) in self._poly.__repn.items():
                term = oneterm._new_c()
                term._coeff = c
                term._exponent = e
                if term._valuation_c() < self._prec:
                    self._terms.append(term)
            self._terms.sort(reverse=True)
            self._terms_nonzero = [ term for term in self._terms if not term.coefficient().is_zero() ]
        if include_zero:
            return self._terms
        else:
            return self._terms_nonzero

    def monomials(self):
        r"""
        Return a list of the monomials of this series.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.monomials()  # indirect doctest
            [...0000000001*x, ...0000000001*x^2]
        """
        return [ t.monomial() for t in self.terms() ]

    def monomial_coefficients(self):
        """
        Return a dictionary whose keys are the exponents and whose values
        are the corresponding coefficients of this series.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.monomial_coefficients()
            {(1, 0): ...0000000001, (2, 0): ...00000000010}

        ``dict`` is an alias::

            sage: f.dict()
            {(1, 0): ...0000000001, (2, 0): ...00000000010}
        """
        self._normalize()
        return dict(self._poly.__repn)

    dict = monomial_coefficients

    def coefficient(self, exponent):
        r"""
        Return the coefficient corresponding to the given exponent.

        INPUT:

        - ``exponent`` -- tuple of integers

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='terse')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + 53*x*y + y^3

            sage: f.coefficient((2,0))   # coeff in x^2
            2 + O(2^11)
            sage: f.coefficient((1,1))   # coeff in x*y
            53 + O(2^10)
            sage: f.coefficient((3,0))   # coeff in x^3
            0

            sage: g = f.add_bigoh(5)
            sage: g.coefficient((2,0))   # coeff in x^2
            2 + O(2^5)
            sage: g.coefficient((1,1))   # coeff in x*y
            21 + O(2^5)
            sage: g.coefficient((3,0))   # coeff in x^3
            O(2^5)
        """
        if not self._is_normalized:
            self._normalize()
        try:
            e = ETuple(exponent)
        except TypeError:
            raise IndexError("%s is not a correct exponent" % exponent)
        if len(e) != self.parent().ngens():
            raise IndexError("lengths do not match")
        if e in self._poly.__repn:
            return self._poly.__repn[e]
        else:
            return self.base_ring()(0, self.precision_absolute())

    def __getitem__(self, exponent):
        r"""
        Return the coefficient corresponding to the given exponent.

        INPUT:

        - ``exponent`` -- tuple of integers

        TESTS::

            sage: R = Zp(2, prec=10, print_mode='terse')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + 53*x*y + y^3

            sage: f['hello']
            Traceback (most recent call last):
            ...
            IndexError: hello is not a correct exponent

            sage: f[1,2,3]
            Traceback (most recent call last):
            ...
            IndexError: lengths do not match
        """
        return self.coefficient(exponent)

    def coefficients(self):
        r"""
        Return the list of coefficients of this series.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2
            sage: f.coefficients()
            [...0000000001, ...00000000010]
        """
        return [ t.coefficient() for t in self.terms() ]

    def add_bigoh(self, n):
        r"""
        Return this series truncated at precision ``n``.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 32*x + 64*x^2; f
            ...000000000100000*x + ...0000000001000000*x^2
            sage: f.add_bigoh(5)
            O(2^5 * <x, y>)

            sage: g = f.add_bigoh(6); g
            ...100000*x + O(2^6 * <x, y>)
            sage: g.precision_absolute()
            6
        """
        return self._parent(self, prec=n)

    def lift_to_precision(self, prec=None):
        """
        Return a lift of this series at precision ``prec``.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``); if
          ``None``, the cap of the parent is used if it is higher than
          the current precision

        EXAMPLES::

            sage: R = Zp(2, prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = R(1,4)*x*y + R(1,5)*x + R(1,8)*y
            sage: f
            (1 + O(2^4))*x*y + (1 + O(2^5))*x + (1 + O(2^8))*y

        This method lifts the precision of the coefficients::

            sage: f.lift_to_precision()
            (1 + O(2^10))*x*y + (1 + O(2^10))*x + (1 + O(2^10))*y

        and also acts on the global ``O(.)`` of the series::

            sage: g = f.add_bigoh(7)
            sage: g
            (1 + O(2^4))*x*y + (1 + O(2^5))*x + (1 + O(2^7))*y + O(2^7 * <x, y>)
            sage: g.lift_to_precision()
            (1 + O(2^10))*x*y + (1 + O(2^10))*x + (1 + O(2^10))*y + O(2^10 * <x, y>)

            sage: g.lift_to_precision(9)
            (1 + O(2^9))*x*y + (1 + O(2^9))*x + (1 + O(2^9))*y + O(2^9 * <x, y>)

        In the next example, the precision on the coefficient is only lifted
        to ``O(2^10)`` because it is limited by the cap of the underlying
        `p`-adic ring::

            sage: g.lift_to_precision(20)
            (1 + O(2^10))*x*y + (1 + O(2^10))*x + (1 + O(2^10))*y + O(2^20 * <x, y>)
        """
        cdef TateAlgebraElement ans = self._new_c()
        # Hmm, shouldn't we add a keyword argument to lift_to_precision()
        # to specify that we don't want it to raise an error

        def lift_without_error(elt):
            try:
                return elt.lift_to_precision(prec)
            except PrecisionError:
                return elt.lift_to_precision()
        ans._poly = PolyDict({e: lift_without_error(c)
                              for e, c in self._poly.__repn.items()},
                             None)
        if prec is None:
            prec = self._parent.precision_cap()
        ans._prec = max(self._prec, prec)
        return ans

    def precision_absolute(self):
        r"""
        Return the maximal precision at which a term of this series is known.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2; f
            ...0000000001*x + ...00000000010*x^2
            sage: f.precision_absolute()
            +Infinity

            sage: g = f.add_bigoh(5); g
            ...00001*x + ...00010*x^2 + O(2^5 * <x, y>)
            sage: g.precision_absolute()
            5

        The absolute precision may be higher than the precision of some
        individual coefficients::

            sage: g = f.add_bigoh(20); g
            ...0000000001*x + ...00000000010*x^2 + O(2^20 * <x, y>)
            sage: g.precision_absolute()
            20
        """
        return self._prec

    cpdef valuation(self):
        r"""
        Return the valuation of this series.

        .. NOTE::

            The valuation of a series `f` is defined as the minimal
            valuation of `f(x)` for `x` varying in the domain of convergence
            (specified in the parent).

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + 4*x*y + 1; f
            ...0000000001*x^4 + ...0000000001 + ...000000000100*x*y
            sage: f.valuation()
            0

            sage: g = 2*f; g
            ...00000000010*x^4 + ...00000000010 + ...0000000001000*x*y
            sage: g.valuation()
            1

        When the radius of convergence is not 1, the variables themselves
        have a nontrivial valuation::

            sage: A.<x,y> = TateAlgebra(R, log_radii=(1,2))
            sage: x.valuation()
            -1
            sage: y.valuation()
            -2

            sage: f = x^4 + 4*x*y + 1
            sage: f.valuation()
            -4
        """
        cdef TateAlgebraTerm t
        cdef list terms = self._terms_c()
        if terms:
            return min(terms[0].valuation(), self._prec)
        else:
            return self._prec

    def precision_relative(self):
        """
        Return the relative precision of this series.

        The relative precision is defined as the difference
        between the absolute precision and the valuation.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R.fraction_field())
            sage: f = x^4 + 4*x*y + 1; f
            ...0000000001*x^4 + ...0000000001 + ...000000000100*x*y
            sage: f.precision_relative()
            +Infinity

            sage: g = f.add_bigoh(5)
            sage: g.precision_relative()
            5
            sage: g.precision_absolute()
            5
            sage: g.valuation()
            0

            sage: h = g + 1/2 ; h
            ...00001.1 + ...00001*x^4 + ...00100*x*y + O(2^5 * <x, y>)
            sage: h.precision_relative()
            6
            sage: h.precision_absolute()
            5
            sage: h.valuation()
            -1
        """
        return self._prec - self.valuation()

    def log(self, prec=None):
        r"""
        Return the logarithm of this series.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``); the
          absolute precision at which the result is computed, if ``None``
          the cap of the Tate algebra is used

        EXAMPLES::

            sage: R = Zp(3, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 3*x + 9*y^2
            sage: f.log()
            ...0000000010*x + ...0000000100*x^3 + ...1111111100*x^2 + ...0000000100*y^2
             + ...2222222000*x*y^2 + ... + O(3^10 * <x, y>)

            sage: f.log(prec=4)
            ...0010*x + ...0100*x^3 + ...1100*x^2 + ...0100*y^2 + ...2000*x*y^2
             + O(3^4 * <x, y>)

        If the precision on the input is not enough to determine the
        result at precision ``prec``, a result with smaller precision
        is returned::

            sage: g = f.add_bigoh(4); g
            ...0001 + ...0010*x + ...0100*y^2 + O(3^4 * <x, y>)
            sage: g.log()
            ...0010*x + ...0100*x^3 + ...1100*x^2 + ...0100*y^2 + ...2000*x*y^2
             + O(3^4 * <x, y>)
            sage: g.log(prec=10)
            ...0010*x + ...0100*x^3 + ...1100*x^2 + ...0100*y^2 + ...2000*x*y^2
             + O(3^4 * <x, y>)

        When the input value is outside the domain of convergence, an
        error is raised::

            sage: f = 1 + x
            sage: f.log()
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence

        However `\log(1+x)` converges on a smaller disk::

            sage: f.restriction(-1).log()
            ...0000000001*x + ...000000000.1*x^3 + ...111111111*x^2 + ...
             + O(3^10 * <3*x, 3*y>)

        TESTS::

            sage: f = 1 + 3 * A.random_element(integral=True)
            sage: logf = f.log()

            sage: x0 = 3 * R.random_element()
            sage: y0 = 3 * R.random_element()
            sage: f(x0, y0).log() == logf(x0, y0)
            True

            sage: logf.exp() == f
            True
        """
        # This code is mostly copied from sage.rings.padics.padic_generic_element
        # (should we find a way to share it?)
        R = self.parent()
        p = R.base_ring().prime()
        e = R.absolute_e()
        x = 1 - self
        alpha = x.valuation()
        if alpha <= 0:
            raise ValueError("not in the domain of convergence")

        aprec = self.precision_absolute()
        if prec is not None and prec < aprec:
            aprec = prec
        if aprec is Infinity:
            aprec = R.precision_cap()

        mina = 0
        if e != 1:
            lamb = aprec - alpha
            if lamb > 0 and lamb*(p-1) <= e:
                # This is the precision region where the absolute
                # precision of the answer might be less than the
                # absolute precision of the input

                # kink is the number of times we multiply the relative
                # precision by p before starting to add e instead.
                kink = (e // (lamb * (p-1))).exact_log(p) + 1

                # deriv0 is within 1 of the n yielding the minimal
                # absolute precision
                tmp = (e / (aprec * p.log(prec=53))).floor()
                if tmp > 0:
                    deriv0 = tmp.exact_log(p)
                else:
                    deriv0 = 0

                # These are the absolute precisions of x^(p^n) at potential minimum points
                L = [(aprec * p**n - n * e, n) for n in [0, kink, deriv0, deriv0+1]]
                L.sort()
                aprec = L[0][0]
                mina = L[0][1]

        total = R.zero()
        if mina == 0 and alpha*p - e > aprec:
            # The value of x^p/p is not needed in that case
            x2p_p = R(0)
        else:
            x2p_p = x**p / p

        # To get result right to precision aprec, we need all terms for which
        # the valuation of x^n/n is strictly smaller than aprec.
        # If we rewrite n=u*p^a with u a p-adic unit, then these are the terms
        # for which u<(aprec+a*v(p))/(v(x)*p^a).
        # Two sum over these terms, we run two nested loops, the outer one
        # iterates over the possible values for a, the inner one iterates over
        # the possible values for u.
        upper_u = (aprec/alpha).floor()
        if mina > 0 or upper_u > 0:
            a=0
            p2a=1       # p^a
            x2pa = x    # x^(p^a)

            # In the unramified case, we can stop summing terms as soon as
            # there are no u for a given a to sum over. In the ramified case,
            # it can happen that for some initial a there are no such u but
            # later in the series there are such u again. mina can be set to
            # take care of this by summing at least to a=mina-1
            while True:
                # we compute the sum for the possible values for u using Horner's method
                inner_sum = R.zero()
                for u in range(upper_u,0,-1):
                    # We want u to be a p-adic unit
                    if u % p==0:
                        new_term = R.zero()
                    else:
                        new_term = ~R.base_ring()(u)

                    # This hack is to deal with rings that don't lift to fields
                    if u > 1 or x2p_p.is_zero():
                        inner_sum = (inner_sum+new_term)*x2pa
                    else:
                        inner_sum = (inner_sum+new_term)*(x2p_p**a)*(x**(p2a-a*p))

                total -= inner_sum

                # Now increase a and check if a new iteration of the loop is needed
                a += 1
                p2a *= p
                upper_u = ((aprec + a*e)/(alpha * p2a)).floor()
                if a >= mina and upper_u <= 0: break

                # We perform this last operation after the test
                # because it is costly and may raise OverflowError
                x2pa = x2pa**p

        return total.add_bigoh(aprec)

    def exp(self, prec=None):
        r"""
        Return the exponential of this series.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``); the
          absolute precision at which the result is computed, if ``None``
          the cap of the Tate algebra is used

        EXAMPLES::

            sage: R = Zp(3, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 3*x^2 + 9*y
            sage: f.exp()
            ...0000000001 + ...0000000010*x^2 + ...1111111200*x^6 + ...1111111200*x^4
             + ...0000000100*y + ... + O(3^10 * <x, y>)

            sage: f.exp(prec=3)
            ...001 + ...010*x^2 + ...200*x^6 + ...200*x^4 + ...100*y + O(3^3 * <x, y>)

        If the precision on the input is not enough to determine the
        result at precision ``prec``, a result with smaller precision
        is returned::

            sage: g = f.add_bigoh(3); g
            ...010*x^2 + ...100*y + O(3^3 * <x, y>)
            sage: g.exp()
            ...001 + ...010*x^2 + ...200*x^6 + ...200*x^4 + ...100*y + O(3^3 * <x, y>)
            sage: g.exp(prec=10)
            ...001 + ...010*x^2 + ...200*x^6 + ...200*x^4 + ...100*y + O(3^3 * <x, y>)

        When the input value is outside the domain of convergence, an
        error is raised::

            sage: f = x
            sage: f.exp()
            Traceback (most recent call last):
            ...
            ValueError: not in the domain of convergence

        However `\exp(x)` converges on a smaller disk::

            sage: f.restriction(-1).exp()
            ...0000000001 + ...0000000001*x + ...111111111.2*x^3 + ...111111112*x^2
             + ... + O(3^10 * <3*x, 3*y>)

        TESTS::

            sage: f = 3 * A.random_element(integral=True)
            sage: expf = f.exp()

            sage: x0 = 3 * R.random_element()
            sage: y0 = 3 * R.random_element()
            sage: f(x0, y0).exp() == expf(x0, y0)
            True

            sage: expf.log() == f  # long time
            True
        """
        # This code is mostly copied from sage.rings.padics.padic_generic_element
        # (should we find a way to share it?)
        A = self.parent()
        R = A.base_ring()
        p = R.prime()
        e = R.absolute_e()
        x_val = self.valuation()

        if (p-1) * x_val <= e:
            raise ValueError("not in the domain of convergence")

        aprec = self.precision_absolute()
        if aprec is Infinity:
            aprec = R.precision_cap()
        if prec is not None and prec < aprec:
            aprec = prec

        # the valuation of n! is bounded by e*n/(p-1), therefore the valuation
        # of self^n/n! is bigger or equal to n*x_val - e*n/(p-1). So, we only
        # have to sum terms for which n does not exceed N
        N = (aprec // (x_val - e/(p-1))).floor()

        # We evaluate the exponential series:
        # We compute the value of x^N + N*x^(N-1) + ... + (N-1)!*x + N!
        # by Horner's method. Then, we divide by N!.
        series = A.one()
        nfactorial = R.one()
        x = self.add_bigoh(aprec)
        for n in range(N, 0, -1):
            series *= x
            nfactorial *= n
            series += nfactorial
        return series / nfactorial

    def leading_term(self, secure=False):
        r"""
        Return the leading term of this series.

        .. NOTE::

            The order on the terms is defined as follows: first we
            compare the valuation and, in case of equality, we compare
            the monomials with respect to the order given at the
            creation of the parent.

        INPUT:

        - ``secure`` -- boolean (default: ``False``); if ``True``,
          raises an error if the leading term cannot be determined
          due to the existence of terms which are indistinguishable
          from zero. If ``False``, discard silently these terms.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + x*y + 1; f
            ...0000000001*x^4 + ...0000000001*x*y + ...0000000001
            sage: f.leading_term()
            ...0000000001*x^4

            sage: g = f + x^4; g
            ...0000000001*x*y + ...0000000001 + ...0000000010*x^4
            sage: g.leading_monomial()
            ...0000000001*x*y

        Observe that the leading term may change after restriction::

            sage: f.restriction(-1).leading_term()
            ...0000000001

        TESTS::

            sage: f = 1 + 64*x
            sage: f -= R(1, 5)
            sage: f
            ...00000 + ...0000000001000000*x

            sage: f.leading_term()
            ...0000000001000000*x
            sage: f.leading_term(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to determine the leading term

            sage: g = A(0, 10); g
            O(2^10 * <x, y>)
            sage: g.leading_term()
            Traceback (most recent call last):
            ...
            ValueError: zero has no leading term
            sage: g.leading_term(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to determine the leading term

        .. SEEALSO::

            :meth:`leading_coefficient`, :meth:`leading_monomial`
        """
        cdef list terms
        cdef TateAlgebraTerm term
        if secure:
            terms = self._terms_c()
        else:
            terms = self._terms_c(include_zero=False)
        if terms:
            term = terms[0]
            if term.coefficient().is_zero():
                raise PrecisionError("not enough precision to determine the leading term")
            return term
        if secure and self._prec is not Infinity:
            raise PrecisionError("not enough precision to determine the leading term")
        else:
            raise ValueError("zero has no leading term")

    def leading_coefficient(self, secure=False):
        """
        Return the leading coefficient of this series.

        .. NOTE::

            The leading coefficient is the coefficient of the leading term.

        INPUT:

        - ``secure`` -- boolean (default: ``False``); if ``True``,
          raises an error if the leading term cannot be determined
          due to the existence of terms which are indistinguishable
          from zero. If ``False``, discard silently these terms.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='terse')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + 3*x*y + 1; f
            (1 + O(2^10))*x^4 + (3 + O(2^10))*x*y + (1 + O(2^10))
            sage: f.leading_coefficient()
            1 + O(2^10)

            sage: g = f + x^4; g
            (3 + O(2^10))*x*y + (1 + O(2^10)) + (2 + O(2^10))*x^4
            sage: g.leading_coefficient()
            3 + O(2^10)

        .. SEEALSO::

            :meth:`leading_term`, :meth:`leading_monomial`
        """
        return self.leading_term(secure=secure).coefficient()

    def leading_monomial(self, secure=False):
        """
        Return the leading coefficient of this series.

        .. NOTE::

            The leading monomial is the monomial of the leading term.

        INPUT:

        - ``secure`` -- boolean (default: ``False``); if ``True``,
          raises an error if the leading term cannot be determined
          due to the existence of terms which are indistinguishable
          from zero. If ``False``, discard silently these terms.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + x*y + 1; f
            ...0000000001*x^4 + ...0000000001*x*y + ...0000000001
            sage: f.leading_monomial()
            ...0000000001*x^4

            sage: g = f + x^4; g
            ...0000000001*x*y + ...0000000001 + ...0000000010*x^4
            sage: g.leading_monomial()
            ...0000000001*x*y

        .. SEEALSO::

            :meth:`leading_term`, :meth:`leading_coefficient`
        """
        return self.leading_term(secure=secure).monomial()

    cpdef TateAlgebraElement monic(self):
        r"""
        Return this series normalized so that it has valuation 0
        and its leading coefficient is a power of the uniformizer.

        EXAMPLES:

        When the log radii of convergence are all zero, the
        leading coefficient of the returned series is `1`::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R, order='lex')
            sage: f = 3*x^2*y^2 + 4*x^2*y^3 + y^5; f
            ...0000000011*x^2*y^2 + ...0000000001*y^5 + ...000000000100*x^2*y^3
            sage: f.monic()
            ...0000000001*x^2*y^2 + ...1010101011*y^5 + ...101010101100*x^2*y^3

        However, when log radii do not vanish, behaviors might
        be different::

            sage: g = f.restriction(-1); g
            ...0000000011*x^2*y^2 + ...0000000001*y^5 + ...000000000100*x^2*y^3
            sage: g.monic()
            ...000000.0001*x^2*y^2 + ...101010.1011*y^5 + ...10101010.11*x^2*y^3
            sage: g.monic().valuation()
            0

        TESTS::

            sage: A(0).monic()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        cdef TateAlgebraElement ans = self._new_c()
        cdef TateAlgebraTerm t
        cdef long shi
        cdef pAdicGenericElement u
        cdef list terms = self._terms_c()
        if not terms or terms[0].coefficient().is_zero():
            raise ZeroDivisionError("rational division by zero")
        t = terms[0]
        shi, u = t._coeff.val_unit()
        shi -= t._exponent.dotprod(self._parent._log_radii)
        ans._poly = self._poly.scalar_lmult((~u) >> shi)
        ans._prec = self._prec - shi
        return ans

    def is_monic(self):
        """
        Return ``True`` if this series is monic, in the sense
        that it has valuation 0 and its leading coefficient is
        a power of the uniformizer.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x*y + 2; f
            ...0000000001*x*y + ...00000000010
            sage: f.is_monic()
            True

            sage: g = f.restriction(-1); g
            ...00000000010 + ...0000000001*x*y
            sage: g.is_monic()
            False
        """
        if self.valuation() != 0:
            return False
        c = self.leading_coefficient()
        return c != 0 and c.unit_part() == 1

    def weierstrass_degree(self):
        r"""
        Return the Weierstrass degree of this Tate series.

        .. NOTE::

            The Weierstrass degree is the total degree of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x*y + 2*x^4 + y; f
            ...0000000001*x*y + ...0000000001*y + ...00000000010*x^4
            sage: f.weierstrass_degree()
            2
        """
        v = self.valuation()
        return self.residue(v+1).degree()

    def degree(self):
        r"""
        Return the Weierstrass degree of this Tate series.

        .. NOTE::

            The Weierstrass degree is the total degree of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x*y + 2*x^4 + y; f
            ...0000000001*x*y + ...0000000001*y + ...00000000010*x^4
            sage: f.degree()
            2
        """
        return self.weierstrass_degree()

    def weierstrass_degrees(self):
        r"""
        Return the Weierstrass degrees of this Tate series.

        .. NOTE::

            The Weierstrass degrees are the partial degrees of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^2 + y^2 + 2*x^3 + y; f
            ...0000000001*x^2 + ...0000000001*y^2 + ...0000000001*y + ...00000000010*x^3
            sage: f.weierstrass_degrees()
            (2, 2)
        """
        v = self.valuation()
        return self.residue(v+1).degrees()

    def degrees(self):
        r"""
        Return the Weierstrass degrees of this series.

        .. NOTE::

            The Weierstrass degrees are the partial degrees of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^2 + y^2 + 2*x^3 + y; f
            ...0000000001*x^2 + ...0000000001*y^2 + ...0000000001*y + ...00000000010*x^3
            sage: f.degrees()
            (2, 2)
        """
        return self.weierstrass_degrees()

    def residue(self, n=None):
        r"""
        Return this series modulo the ``n``-th power of the uniformizer.

        Note that by definition of Tate series, the output is a polynomial.

        INPUT:

        - ``n`` -- integer (default: `1`)

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^2 + y^2 + 6*x^3 + 3*y
            sage: f.residue()
            x^2 + y^2 + y
            sage: f.residue().parent()
            Multivariate Polynomial Ring in x, y over Finite Field of size 2

            sage: f.residue(2)
            2*x^3 + x^2 + y^2 - y
            sage: f.residue(2).parent()
            Multivariate Polynomial Ring in x, y over Ring of integers modulo 4

        The residue can only be computed for series with nonnegative valuation.

            sage: g = f >> 2; g
            ...00000000.01*x^2 + ...00000000.01*y^2 + ...00000000.11*y + ...000000001.1*x^3
            sage: g.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue

        The residue is not implemented for series with convergence radius different from 1.

            sage: A.<x,y> = TateAlgebra(R, log_radii=(2,-1))
            sage: f = x^2 + y^2 + 6*x^3 + 3*y
            sage: f.residue()
            Traceback (most recent call last):
            ...
            NotImplementedError: residues are only implemented for radius 1
        """
        for r in self._parent.log_radii():
            if r != 0:
                raise NotImplementedError("residues are only implemented for radius 1")
        if n is None:
            Rn = self.base_ring().residue_field()
        else:
            try:
                Rn = self.base_ring().residue_ring(n)
            except (AttributeError, NotImplementedError):
                Rn = self.base_ring().change(field=False, type='fixed-mod', prec=n)
        poly = self._parent._polynomial_ring(self._poly)
        return poly.change_ring(Rn)

    cdef _quo_rem_c(self, list divisors, bint quo, bint rem, bint integral):
        r"""
        Perform the division of this series by ``divisors``.

        INPUT:

        - ``divisors`` -- the list of divisors

        - ``quo`` -- boolean; whether we should compute the quotients

        - ``rem`` -- boolean; whether we should compute the remainder

        TESTS::

            sage: R = Zp(2, 5, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: q, r = f.quo_rem(g)  # indirect doctest
            sage: q
            ...00011*y
            sage: r
            ...00001 + ...00010*x*y + ...00100*x*y^2 + O(2^5 * <x, y>)
        """
        cdef dict coeffs = { }
        cdef TateAlgebraElement f
        cdef TateAlgebraTerm lt
        cdef list ltds = [ (<TateAlgebraElement>d)._terms_c()[0] for d in divisors ]
        cdef list quos = [ ]
        cdef list terms = self._terms_c()
        cdef int index = 0
        cdef bint in_rem

        if quo:
            for d in divisors:
                f = self._new_c()
                f._poly = PolyDict({}, None)
                f._prec = d.precision_relative()
                quos.append(f)

        f = self._new_c()
        f._poly = PolyDict(self._poly.__repn, None)
        f._prec = self._prec
        while len(terms) > index:
            lt = terms[index]
            if lt._valuation_c() >= f._prec:
                break
            in_rem = True
            if not lt._coeff.is_zero():  # is it a good idea?
                for i in range(len(divisors)):
                    if (<TateAlgebraTerm>ltds[i])._divides_c(lt, integral=integral):
                        factor = lt._floordiv_c(<TateAlgebraTerm>ltds[i])
                        f = f - (<TateAlgebraElement>divisors[i])._term_mul_c(factor)
                        terms = f._terms_c()
                        index = 0
                        if quo:
                            quos[i] = (<TateAlgebraElement>quos[i]) + factor
                        in_rem = False
                        break
            if in_rem:
                if rem:
                    if lt._exponent in coeffs:
                        coeffs[lt._exponent] += lt._coeff
                    else:
                        coeffs[lt._exponent] = lt._coeff
                del f._poly.__repn[lt._exponent]
                index += 1
        if rem:
            f._poly += PolyDict(coeffs, None)
            f._terms = None
        return quos, f

    cdef _quo_rem_check(self, divisors, bint quo, bint rem):
        """
        Perform the division of this series by ``divisors``.

        INPUT:

        - ``divisors`` -- the list of divisors

        - ``quo`` -- boolean; whether we should compute the quotients

        - ``rem`` -- boolean; whether we should compute the remainder

        TESTS::

            sage: R = Zp(2, 10)
            sage: A.<u,v> = TateAlgebra(R)
            sage: Ao = A.integer_ring()

            sage: f = Ao(u^2 + 2*v^2)
            sage: f.quo_rem(u)  # indirect doctest
            ((1 + O(2^10))*u, (2 + O(2^10))*v^2 + O(2^10 * <u, v>))

        We check that coercion works::

            sage: f % 2  # indirect doctest
            (1 + O(2^10))*u^2 + O(2^10 * <u, v>)
            sage: f % [2,u]  # indirect doctest
            O(2^10 * <u, v>)

            sage: x = polygen(ZZ, 'x')
            sage: S.<pi> = R.extension(x^2 - 2)
            sage: f % (pi*u)  # indirect doctest
            (pi^2 + O(pi^20))*v^2 + O(pi^20 * <u, v>)
            sage: (pi*f) // (pi*u) == f // u
            True

       ::

            sage: s = Zp(3)(1)
            sage: f % s
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce all the elements to the same parent

            sage: f % ZZ
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce all the elements to the same parent
        """
        parent = self.parent()
        if self.precision_absolute() is Infinity:
            self = self.add_bigoh(self.valuation() + parent.precision_cap())
        if isinstance(divisors, (list, tuple)):
            onedivisor = False
        else:
            divisors = [divisors]
            onedivisor = True
        A = _pushout_family(divisors, self._parent)
        f = A(self)
        divisors = [A(d) for d in divisors]
        q, r = (<TateAlgebraElement>self)._quo_rem_c(divisors, quo, rem, False)
        if quo and onedivisor:
            q = q[0]
        if quo and rem:
            return q, r
        if quo:
            return q
        if rem:
            return r

    def quo_rem(self, divisors):
        """
        Return the quotient(s) and the remainder of the division of
        this series by ``divisors``.

        INPUT:

        - ``divisors`` -- a series, or a list of series

        NOTE:

        The condition on the remainder is that it has

        - no term which is greater than the leading term of the
          numerator and

        - no term which is divisible by the leading term of one
          divisor.

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: q, r = f.quo_rem(g)
            sage: q
            ...00011*y
            sage: r
            ...00001 + ...00010*x*y + ...00100*x*y^2 + O(2^5 * <x, y>)
            sage: f == g*q + r
            True

        We can also divide by a family of divisors::

            sage: g0 = x^2
            sage: g1 = x*y + 2*x
            sage: q, r = f.quo_rem([g0, g1])
            sage: q
            [...00011*y, ...11010 + ...00100*y]
            sage: r
            ...00001 + ...01100*x + O(2^5 * <x, y>)
            sage: f == g0*q[0] + g1*q[1] + r
            True
        """
        return (<TateAlgebraElement>self)._quo_rem_check(divisors, True, True)

    def __mod__(self, divisors):
        """
        Return the remainder of the division of this series by
        ``divisors``.

        INPUT:

        - ``divisors`` -- a series, or a list of series

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: f % g
            ...00001 + ...00010*x*y + ...00100*x*y^2 + O(2^5 * <x, y>)

        We can also divide by a family of divisors::

            sage: g0 = x^2
            sage: g1 = x*y + 2*x
            sage: f % [g0, g1]
            ...00001 + ...01100*x + O(2^5 * <x, y>)
        """
        return (<TateAlgebraElement>self)._quo_rem_check(divisors, False, True)

    def __floordiv__(self, divisors):
        """
        Return the quotient(s) of the division of this series by
        ``divisors``.

        INPUT:

        - ``divisors`` -- a series, or a list of series

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: f // g
            ...00011*y

        We can also divide by a family of divisors::

            sage: g0 = x^2
            sage: g1 = x*y + 2*x
            sage: f // [g0, g1]
            [...00011*y, ...11010 + ...00100*y]
        """
        return (<TateAlgebraElement>self)._quo_rem_check(divisors, True, False)

    def reduce(self, I):
        """
        Return a canonical representative of this series in the
        quotient of the Tate algebra (in which this series lives)
        by the ideal ``I``.

        EXAMPLES::

            sage: R = Zp(3, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 3*x^2 + 5*x*y^2
            sage: g = 5*x^2*y + 3
            sage: I = A.ideal([f, g])

            sage: f.reduce(I)
            O(3^9 * <x, y>)
            sage: h = (x^2 + 2*y)*f + (x^2*y^3 + 3*x*y^2 + 7)*g + 1
            sage: h.reduce(I)
            ...000000001 + O(3^9 * <x, y>)

        TESTS::

            sage: s = I.random_element(integral=True)
            sage: s.reduce(I).precision_absolute() >= 9
            True

            sage: h = A.random_element()
            sage: (h + s).reduce(I) == h.reduce(I)
            True
        """
        return self % I.groebner_basis()

    @coerce_binop
    def Spoly(self, other):
        """
        Return the S-polynomial of this series and ``other``.

        INPUT:

        - ``other`` -- a Tate series

        NOTE:

        If `f` and `g` are two Tate series with leading term
        `t_f` and `t_g` respectively, the S-polynomial of `f`
        and `g` is defined by

        .. MATH::

            S(f,g) = \frac{\text{lcm}(t_f,t_g)}{t_f}} f - \frac{\text{lcm}(t_f,t_g)}{t_g}} g

        By construction the terms in `\text{lcm}(t_f,t_g)` cancel,
        so that the leading term of `S(f,g)` is strictly smaller
        than `\text{lcm}(t_f,t_g)`.

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^3*y + 2*x*y + 4*x^2
            sage: g = 2*x*y^2 + 2*x
            sage: h = f.Spoly(g); h
            ...111110*x^3 + ...0000100*x*y^2 + ...00001000*x^2*y

            sage: h == 2*y*f - x^2*g
            True

        TESTS::

            sage: f.Spoly(0)
            Traceback (most recent call last):
            ...
            ValueError: the S-polynomial of zero is not defined
        """
        try:
            return self._Spoly_c(other)
        except IndexError:
            raise ValueError("the S-polynomial of zero is not defined")

    cdef TateAlgebraElement _Spoly_c(self, TateAlgebraElement other):
        """
        Return the S-polynomial of this series and ``other``.

        INPUT:

        - ``other`` -- a Tate series

        TESTS:

        We check that the S-polynomial of two monomials vanishes::

            sage: R = Zp(3, 5, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^3*y^2
            sage: g = x^2*y^3
            sage: f.Spoly(g)
            0
        """
        cdef TateAlgebraTerm st = self._terms_c()[0]
        cdef TateAlgebraTerm ot = other._terms_c()[0]
        cdef TateAlgebraTerm t = st._lcm_c(ot)
        cdef TateAlgebraElement ans = self._term_mul_c(t._floordiv_c(st)) - other._term_mul_c(t._floordiv_c(ot))
        ans._poly.__repn.pop(t._exponent, None)
        return ans
