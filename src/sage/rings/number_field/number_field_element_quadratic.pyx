"""
Optimized Quadratic Number Field Elements

This file defines a Cython class ``NumberFieldElement_quadratic`` to speed up
computations in quadratic extensions of `\QQ`.

AUTHORS:

- Robert Bradshaw (2007-09): Initial version
- David Harvey (2007-10): fix up a few bugs, polish around the edges
- David Loeffler (2009-05): add more documentation and tests
- Vincent Delecroix (2012-07): comparisons for quadratic number fields
  (:trac:`13213`), abs, floor and ceil functions (:trac:`13256`)

.. TODO::

    The ``_new()`` method should be overridden in this class to copy the ``D``
    and ``standard_embedding`` attributes
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/signals.pxi"
include "sage/ext/stdsage.pxi"
include "sage/libs/ntl/decl.pxi"
from cpython.object cimport Py_EQ, Py_NE, Py_LE, Py_GE, Py_LT, Py_GT

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX

from sage.structure.parent_base cimport ParentWithBase
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.structure.sage_object cimport rich_to_bool_sgn

from sage.rings.rational cimport Rational
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.categories.morphism cimport Morphism
from sage.rings.number_field.number_field_element import _inverse_mod_generic

import number_field

from sage.libs.gmp.pylong cimport mpz_pythonhash


def __make_NumberFieldElement_quadratic0(parent, a, b, denom):
    """
    Used in unpickling elements of number fields.

    TEST::

        sage: K.<a> = NumberField(x^2-x+13)
        sage: loads(dumps(a)) == a # indirect doctest
        True
    """
    return NumberFieldElement_quadratic(parent, (a, b, denom))

def __make_NumberFieldElement_quadratic1(parent, cls, a, b, denom):
    """
    Used in unpickling elements of number fields.

    TEST::

        sage: K.<a> = NumberField(x^2-x+13)
        sage: loads(dumps(a)) == a # indirect doctest
        True

    We test that :trac:`6462` is fixed::

        sage: L = QuadraticField(-11,'a'); OL = L.maximal_order(); w = OL.0
        sage: loads(dumps(w)) == w # indirect doctest
        True
    """
    return cls(parent, (a, b, denom))


cdef class NumberFieldElement_quadratic(NumberFieldElement_absolute):
    r"""
    A NumberFieldElement_quadratic object gives an efficient representation of
    an element of a quadratic extension of `\QQ`.

    Elements are represented internally as triples `(a, b, c)` of integers,
    where `{\rm gcd}(a, b, c) = 1` and `c > 0`, representing the element `(a +
    b \sqrt{D}) / c`. Note that if the discriminant `D` is `1 \bmod 4`,
    integral elements do not necessarily have `c = 1`.

    TESTS::

        sage: from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic

    We set up some fields::

        sage: K.<a> = NumberField(x^2+23)
        sage: a.parts()
        (0, 1)
        sage: F.<b> = NumberField(x^2-x+7)
        sage: b.parts()
        (1/2, 3/2)

    We construct elements of these fields in various ways - firstly, from
    polynomials::

        sage: NumberFieldElement_quadratic(K, x-1)
        a - 1
        sage: NumberFieldElement_quadratic(F, x-1)
        b - 1

    From triples of Integers::

        sage: NumberFieldElement_quadratic(K, (1,2,3))
        2/3*a + 1/3
        sage: NumberFieldElement_quadratic(F, (1,2,3))
        4/9*b + 1/9
        sage: NumberFieldElement_quadratic(F, (1,2,3)).parts()
        (1/3, 2/3)

    From pairs of Rationals::

        sage: NumberFieldElement_quadratic(K, (1/2,1/3))
        1/3*a + 1/2
        sage: NumberFieldElement_quadratic(F, (1/2,1/3))
        2/9*b + 7/18
        sage: NumberFieldElement_quadratic(F, (1/2,1/3)).parts()
        (1/2, 1/3)

    Direct from Rationals::

        sage: NumberFieldElement_quadratic(K, 2/3)
        2/3
        sage: NumberFieldElement_quadratic(F, 2/3)
        2/3

    This checks a bug when converting from lists::

        sage: w = CyclotomicField(3)([1/2,1])
        sage: w == w.__invert__().__invert__()
        True
    """

    def __init__(self, parent, f):
        """
        Standard initialisation function.

        EXAMPLE::

            sage: F.<a> = QuadraticField(-7)
            sage: c = a + 7
            sage: type(c) # indirect doctest
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
        """
        self.D = parent._D
        cdef Integer a, b, denom
        cdef Rational ad, bd

        cdef NumberFieldElement_quadratic gen

        if isinstance(f, NumberFieldElement_quadratic):
            self._parent = parent   # NOTE: We do *not* call NumberFieldElement_absolute.__init__, for speed reasons.
            mpz_set(self.a, (<NumberFieldElement_quadratic>f).a)
            mpz_set(self.b, (<NumberFieldElement_quadratic>f).b)
            mpz_set(self.denom, (<NumberFieldElement_quadratic>f).denom)

        elif type(f) is tuple and len(f) == 2:
            NumberFieldElement_absolute.__init__(self, parent, None)
            ad, bd = f
            mpz_lcm(self.denom, mpq_denref(ad.value), mpq_denref(bd.value))
            mpz_divexact(self.a, self.denom, mpq_denref(ad.value))
            mpz_mul(self.a, self.a, mpq_numref(ad.value))
            mpz_divexact(self.b, self.denom, mpq_denref(bd.value))
            mpz_mul(self.b, self.b, mpq_numref(bd.value))

        elif type(f) is tuple and len(f) == 3:
            NumberFieldElement_absolute.__init__(self, parent, None)
            a, b, denom = f
            mpz_set(self.a, a.value)
            mpz_set(self.b, b.value)
            mpz_set(self.denom, denom.value)
            self._reduce_c_()

        else:
            NumberFieldElement_absolute.__init__(self, parent, f)
            # poly is in gen (which may not be sqrt(d))
            self._ntl_coeff_as_mpz(self.a, 0)
            self._ntl_coeff_as_mpz(self.b, 1)
            if mpz_cmp_ui(self.a, 0) or mpz_cmp_ui(self.b, 0):
                gen = parent.gen()  # should this be cached?
                self._ntl_denom_as_mpz(self.denom)
                if mpz_cmp_ui(self.b, 0):
                    mpz_mul(self.a, self.a, gen.denom)
                    mpz_addmul(self.a, self.b, gen.a)
                    mpz_mul(self.b, self.b, gen.b)
                    mpz_mul(self.denom, self.denom, gen.denom)
            else:
                mpz_set_ui(self.denom, 1)
            self._reduce_c_()

        # set the attribute standard embedding which is used in the methods
        # __cmp__, sign, real, imag, floor, ceil, ...
        self.standard_embedding = parent._standard_embedding

    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement_quadratic with the
        same parent as self.

        EXAMPLE::

            sage: F.<b> = CyclotomicField(3)
            sage: b + b # indirect doctest
            2*b
        """
        cdef type t = type(self)
        cdef NumberFieldElement_quadratic x = <NumberFieldElement_quadratic>t.__new__(t)
        x._parent = self._parent
        x.standard_embedding = self.standard_embedding
        x.D = self.D
        return x

    cdef number_field(self):
        r"""
        Return the number field to which this element belongs. Since this is a
        Cython cdef method, it is not directly accessible by the user, but the
        function "_number_field" calls this one.

        EXAMPLE::

            sage: F.<b> = QuadraticField(-7)
            sage: b._number_field() # indirect doctest
            Number Field in b with defining polynomial x^2 + 7
        """
        return self._parent

    def _maxima_init_(self, I=None):
        """
        EXAMPLES::

            sage: K.<a> = QuadraticField(-1)
            sage: f = 1 + a
            sage: f._maxima_init_()
            '1+%i*1'
        """
        a = self.parent().gen()
        if a**2 == -1:
            x0, x1 = self
            return str(x0) + "+" + "%i*" + str(x1)
        else:
            NumberFieldElement_absolute._maxima_init_(self, I)

    def __copy__(self):
        r"""
        Returns a new copy of self.

        TESTS::

            sage: K.<a> = QuadraticField(-3)
            sage: b = a + 3
            sage: c = b.__copy__()
            sage: b
            a + 3
            sage: c
            a + 3
            sage: b is c
            False
            sage: b == c
            True
        """
        cdef NumberFieldElement_quadratic x = <NumberFieldElement_quadratic>self._new()
        mpz_set(x.a, self.a)
        mpz_set(x.b, self.b)
        mpz_set(x.denom, self.denom)
        return x


    def __cinit__(self):
        r"""
        Initialisation function.

        EXAMPLE::

            sage: QuadraticField(-3, 'a').gen() # indirect doctest
            a
        """
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.denom)


    def __dealloc__(self):
        mpz_clear(self.a)
        mpz_clear(self.b)
        mpz_clear(self.denom)


    def __reduce__(self):
        """
        Used for pickling.

        TEST:

            sage: K.<a> = NumberField(x^2-13)
            sage: loads(dumps(a)) == a
            True
            sage: loads(dumps(a/3+5)) == a/3+5
            True
        """
        cdef Integer a = <Integer>PY_NEW(Integer)
        cdef Integer b = <Integer>PY_NEW(Integer)
        cdef Integer denom = <Integer>PY_NEW(Integer)
        mpz_set(a.value, self.a)
        mpz_set(b.value, self.b)
        mpz_set(denom.value, self.denom)
        return __make_NumberFieldElement_quadratic1, (self._parent, type(self), a, b, denom)

    cdef int _randomize(self, num_bound, den_bound, distribution) except -1:
        cdef Integer temp, denom1, denom2

        # in theory, we could just generate two random numerators and
        # a random denominator. however, this would mean that we were
        # extraordinarily unlikely to run into results of the form
        # 1/3 + 1/5*sqrt(D), which are often some of the best examples
        # for testing out code. since this is probably the primary use
        # of the random element code, it's worth doing slightly more
        # work to make this possible.

        # normalize denominator bound
        if den_bound is None or den_bound < 1:
            den_bound = 1

        # generate denominators
        denom1 = <Integer>(ZZ.random_element(x=1,
                                             y=den_bound+1,
                                             distribution=distribution))
        denom2 = <Integer>(ZZ.random_element(x=1,
                                             y=den_bound+1,
                                             distribution=distribution))

        # set a, b
        temp = <Integer>(ZZ.random_element(x=num_bound, distribution=distribution))
        mpz_mul(self.a, temp.value, denom2.value)
        temp = <Integer>(ZZ.random_element(x=num_bound, distribution=distribution))
        mpz_mul(self.b, temp.value, denom1.value)
        # set denom
        mpz_mul(self.denom, denom1.value, denom2.value)

        self._reduce_c_()
        return 0  # No error


    def _lift_cyclotomic_element(self, new_parent, bint check=True, int rel=0):
        """
        Creates an element of the passed field from this field.  This
        is specific to creating elements in a cyclotomic field from
        elements in another cyclotomic field, in the case that
        self.number_field()._n() divides new_parent()._n().  This
        function aims to make this common coercion extremely fast!

        More general coercion (i.e. of zeta6 into CyclotomicField(3))
        is implemented in the _coerce_from_other_cyclotomic_field
        method of a CyclotomicField.

        EXAMPLES::

            sage: C.<zeta4>=CyclotomicField(4)
            sage: CyclotomicField(20)(zeta4+1)  # The function _lift_cyclotomic_element does the heavy lifting in the background
            zeta20^5 + 1
            sage: (zeta4+1)._lift_cyclotomic_element(CyclotomicField(40))  # There is rarely a purpose to call this function directly
            zeta40^10 + 1

            sage: cf3 = CyclotomicField(3) ; z3 = cf3.0
            sage: cf6 = CyclotomicField(6) ; z6 = cf6.0
            sage: z6._lift_cyclotomic_element(cf3)
            Traceback (most recent call last):
            ...
            TypeError: The zeta_order of the new field must be a multiple of the zeta_order of the original.
            sage: cf3(z6)
            zeta3 + 1
            sage: z3._lift_cyclotomic_element(cf6)
            zeta6 - 1

        Verify embeddings are respected::

            sage: cf6c = CyclotomicField(6, embedding=CDF(exp(-pi*I/3))) ; z6c = cf6c.0
            sage: cf3(z6c)
            -zeta3
            sage: cf6c(z3)
            -zeta6

        AUTHOR:

        - Joel B. Mohler (original version)

        - Craig Citro (reworked for quadratic elements)
        """
        if check:
            if not isinstance(self.number_field(), number_field.NumberField_cyclotomic) \
                   or not isinstance(new_parent, number_field.NumberField_cyclotomic):
                raise TypeError("The field and the new parent field must both be cyclotomic fields.")

        if rel == 0:
            small_order = self.number_field()._n()
            large_order = new_parent._n()

            try:
                rel = ZZ(large_order / small_order)
            except TypeError:
                raise TypeError("The zeta_order of the new field must be a multiple of the zeta_order of the original.")

        cdef NumberFieldElement_quadratic x2
        cdef int n = self._parent._n()

        if new_parent.degree() == 2:
            ## since self is a *quadratic* element, we can only get
            ## here if self.parent() and new_parent are:
            ## - CyclotomicField(3) and CyclotomicField(6)
            ## - CyclotomicField(3) and CyclotomicField(3)
            ## - CyclotomicField(6) and CyclotomicField(6)
            ## - CyclotomicField(4) and CyclotomicField(4)
            ## In all cases, conversion of elements is trivial!
            if n == <int>new_parent._n():
                conjugate = rel != 1
            else:
                # n = 3, new_n = 6
                conjugate = rel == 4
            x2 = <NumberFieldElement_quadratic>(self._new())
            x2._parent = new_parent
            mpz_set(x2.a, self.a)
            if conjugate:
                mpz_neg(x2.b, self.b)
            else:
                mpz_set(x2.b, self.b)
            mpz_set(x2.denom, self.denom)
            x2.D = self.D
            return x2

        cdef NumberFieldElement x
        cdef ZZX_c elt_num
        cdef ZZ_c elt_den, tmp_coeff
        cdef mpz_t tmp_mpz
        cdef long tmp_const

        x = <NumberFieldElement_absolute>NumberFieldElement_absolute.__new__(NumberFieldElement_absolute)

        mpz_to_ZZ(&elt_den, self.denom)

        mpz_init(tmp_mpz)

        ## set the two terms in the polynomial
        if n == 4:
            mpz_to_ZZ(&tmp_coeff, self.a)
            ZZX_SetCoeff(elt_num, 0, tmp_coeff)
            mpz_to_ZZ(&tmp_coeff, self.b)
            ZZX_SetCoeff(elt_num, 1, tmp_coeff)

        elif n == 3:
            ## num[0] = a + b
            mpz_add(tmp_mpz, tmp_mpz, self.a)
            mpz_add(tmp_mpz, tmp_mpz, self.b)
            mpz_to_ZZ(&tmp_coeff, tmp_mpz)
            ZZX_SetCoeff(elt_num, 0, tmp_coeff)

            ## num[1] = 2*b
            mpz_sub(tmp_mpz, tmp_mpz, self.a)
            tmp_const = 2
            mpz_mul_si(tmp_mpz, tmp_mpz, tmp_const)
            mpz_to_ZZ(&tmp_coeff, tmp_mpz)
            ZZX_SetCoeff(elt_num, 1, tmp_coeff)

        elif n == 6:
            ## num[0] = a - b
            mpz_add(tmp_mpz, tmp_mpz, self.a)
            mpz_sub(tmp_mpz, tmp_mpz, self.b)
            mpz_to_ZZ(&tmp_coeff, tmp_mpz)
            ZZX_SetCoeff(elt_num, 0, tmp_coeff)

            ## num[1] = 2*b
            mpz_sub(tmp_mpz, tmp_mpz, self.a)
            tmp_const = -2
            mpz_mul_si(tmp_mpz, tmp_mpz, tmp_const)
            mpz_to_ZZ(&tmp_coeff, tmp_mpz)
            ZZX_SetCoeff(elt_num, 1, tmp_coeff)

        mpz_clear(tmp_mpz)

        x._parent = <ParentWithBase>new_parent
        x.__fld_numerator, x.__fld_denominator = new_parent.polynomial_ntl()
        x.__denominator = elt_den
        cdef ZZX_c result
        cdef ZZ_c tmp
        cdef int i
        cdef ntl_ZZX _num
        cdef ntl_ZZ _den
        _num, _den = new_parent.polynomial_ntl()
        for i from 0 <= i <= ZZX_deg(elt_num):
            tmp = ZZX_coeff(elt_num, i)
            ZZX_SetCoeff(result, i*rel, tmp)
        ZZX_rem(x.__numerator, result, _num.x)
        (<NumberFieldElement_absolute>x)._reduce_c_()
        return x

    def parts(self):
        """
        This function returns a pair of rationals `a` and `b` such that self `=
        a+b\sqrt{D}`.

        This is much closer to the internal storage format of the
        elements than the polynomial representation coefficients (the output of
        ``self.list()``), unless the generator with which this number field was
        constructed was equal to `\sqrt{D}`. See the last example below.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-13)
            sage: K.discriminant()
            13
            sage: a.parts()
            (0, 1)
            sage: (a/2-4).parts()
            (-4, 1/2)
            sage: K.<a> = NumberField(x^2-7)
            sage: K.discriminant()
            28
            sage: a.parts()
            (0, 1)
            sage: K.<a> = NumberField(x^2-x+7)
            sage: a.parts()
            (1/2, 3/2)
            sage: a._coefficients()
            [0, 1]
        """
        cdef Rational ad = <Rational>Rational.__new__(Rational)
        cdef Rational bd = <Rational>Rational.__new__(Rational)
        if mpz_cmp_ui(self.a, 0) == 0:
            mpq_set_ui(ad.value, 0, 1)
        else:
            mpz_set(mpq_numref(ad.value), self.a)
            mpz_set(mpq_denref(ad.value), self.denom)
            mpq_canonicalize(ad.value)
        if mpz_cmp_ui(self.b, 0) == 0:
            mpq_set_ui(bd.value, 0, 1)
        else:
            mpz_set(mpq_numref(bd.value), self.b)
            mpz_set(mpq_denref(bd.value), self.denom)
            mpq_canonicalize(bd.value)

        return ad, bd

    cdef bint is_sqrt_disc(self):
        r"""
        Returns true if self is `\sqrt{D}`.

        EXAMPLE::

            sage: F.<b> = NumberField(x^2 - x + 7)
            sage: b.denominator() # indirect doctest
            1
        """
        return mpz_cmp_ui(self.denom, 1)==0 and mpz_cmp_ui(self.a, 0)==0 and mpz_cmp_ui(self.b, 1)==0

#########################################################
# Comparisons
#########################################################

    def sign(self):
        r"""
        Returns the sign of self (0 if zero, +1 if positive and -1 if negative).

        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2, name='sqrt2')
            sage: K(0).sign()
            0
            sage: sqrt2.sign()
            1
            sage: (sqrt2+1).sign()
            1
            sage: (sqrt2-1).sign()
            1
            sage: (sqrt2-2).sign()
            -1
            sage: (-sqrt2).sign()
            -1
            sage: (-sqrt2+1).sign()
            -1
            sage: (-sqrt2+2).sign()
            1

            sage: K.<a> = QuadraticField(2, embedding=-1.4142)
            sage: K(0).sign()
            0
            sage: a.sign()
            -1
            sage: (a+1).sign()
            -1
            sage: (a+2).sign()
            1
            sage: (a-1).sign()
            -1
            sage: (-a).sign()
            1
            sage: (-a-1).sign()
            1
            sage: (-a-2).sign()
            -1

            sage: K.<b> = NumberField(x^2 + 2*x + 7, 'b', embedding=CC(-1,-sqrt(6)))
            sage: b.sign()
            Traceback (most recent call last):
            ...
            ValueError: a complex number has no sign!
            sage: K(1).sign()
            1
            sage: K(0).sign()
            0
            sage: K(-2/3).sign()
            -1
        """
        cdef mpz_t i, j
        cdef int s = 1, test

        if mpz_sgn(self.b) == 0:
            return mpz_sgn(self.a)

        if mpz_sgn(self.D.value) == -1:
            raise ValueError("a complex number has no sign!")

        if not self.standard_embedding:
            s = -1

        if mpz_sgn(self.a) == 0:
            return s*mpz_sgn(self.b)

        if mpz_sgn(self.a) == 1:
            if mpz_sgn(self.b) == s:
                return 1

        elif mpz_sgn(self.b) == -s:
            return -1

        mpz_init_set(i,self.a)
        mpz_mul(i,i,i)
        mpz_init_set(j,self.b)
        mpz_mul(j,j,j)
        mpz_mul(j,j,self.D.value)
        test = mpz_cmp(i,j)
        mpz_clear(i)
        mpz_clear(j)
        if mpz_sgn(self.a) == 1 and mpz_sgn(self.b) == -s:
            return test
        return -test

    cpdef _richcmp_(left, _right, int op):
        r"""
        Rich comparison of elements.

        TESTS::

            sage: K.<i> = QuadraticField(-1)
            sage: sorted([5*i+1, 2, 3*i+1, 2-i])
            [3*i + 1, 5*i + 1, -i + 2, 2]

        Make some random tests to check that the order is compatible with the
        ones of the real field (RR) and complex field (CC)::

            sage: K1 = NumberField(x^2 - 2, 'a', embedding=RR(1.4))
            sage: K2 = NumberField(x^2 - 2, 'a', embedding=RR(-1.4))
            sage: for _ in xrange(500):
            ....:     for K in K1, K2:
            ....:         a = K.random_element()
            ....:         b = K.random_element()
            ....:         assert (a < b) == (RR(a) < RR(b))
            ....:         assert (a > b) == (RR(a) > RR(b))
            ....:         assert (a == b) == (RR(a) == RR(b))
            ....:         assert (a != b) == (RR(a) != RR(b))
            ....:         assert (a >= b) == (RR(a) >= RR(b))
            ....:         assert (a <= b) == (RR(a) <= RR(b))

        ::

            sage: K1 = NumberField(x^2 + 2, 'a', embedding=CC(0,1))
            sage: K2 = NumberField(x^2 + 2, 'a', embedding=CC(0,-1))
            sage: for _ in xrange(500):
            ....:     for K in K1, K2:
            ....:         a = K.random_element()
            ....:         b = K.random_element()
            ....:         assert (a < b) == (CC(a) < CC(b))
            ....:         assert (a > b) == (CC(a) > CC(b))
            ....:         assert (a == b) == (CC(a) == CC(b))
            ....:         assert (a != b) == (CC(a) != CC(b))
            ....:         assert (a >= b) == (CC(a) >= CC(b))
            ....:         assert (a <= b) == (CC(a) <= CC(b))

        The following is tested because of the implementation of
        func:`Q_to_quadratic_field_element` which was the cause of some problems
        with :trac:`13213`::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: 1/2 + sqrt2 > 0
            True
        """
        # When D > 0 and standard embedding, we compare (a + b * sqrt(D)) / d and (aa +
        # bb * sqrt(D)) / dd using the comparison of (dd*a - d * aa)^2 and (d*bb - dd*b)^2 * D
        # mpz_sgn: returns 1 if > 0, 0 if 0 and -1 if < 0
        cdef mpz_t i, j
        cdef NumberFieldElement_quadratic right = <NumberFieldElement_quadratic> _right
        cdef int test

        # inequality and equality
        if mpz_cmp(left.a, right.a) or mpz_cmp(left.b, right.b) or mpz_cmp(left.denom, right.denom):
            if op == Py_EQ:
                return False
            elif op == Py_NE:
                return True
        else: # equality
            if op == Py_EQ or op == Py_LE or op == Py_GE:
                return True
            if op == Py_NE or op == Py_LT or op == Py_GT:
                return False

        # comparisons are valid only in *real* quadratic number field
        # when no embedding is specified or in the case of complex embeddings we
        # use a lexicographic order.
        if mpz_sgn(left.D.value) == -1:
            mpz_init(i)
            mpz_init(j)
            mpz_mul(i, left.a, right.denom)
            mpz_mul(j, right.a, left.denom)
            test = mpz_cmp(i,j)
            if test:
                mpz_clear(i)
                mpz_clear(j)
                return rich_to_bool_sgn(op, test)
            mpz_mul(i, left.b, right.denom)
            mpz_mul(j, right.b, left.denom)
            test = mpz_cmp(i,j)
            if test:
                if not left.standard_embedding:
                    test = -test
                mpz_clear(i)
                mpz_clear(j)
                return rich_to_bool_sgn(op, test)
            test = mpz_cmp(left.denom, right.denom)
            mpz_clear(i)
            mpz_clear(j)
            return rich_to_bool_sgn(op, test)

        # comparison in the real case
        mpz_init(i)
        mpz_mul(i,  right.denom, left.a)
        mpz_submul(i, left.denom, right.a)

        mpz_init(j)
        mpz_mul(j, left.denom, right.b)
        mpz_submul(j, right.denom, left.b)

        if not left.standard_embedding:
            mpz_neg(j, j)

        if mpz_sgn(i) == 1:
            if mpz_sgn(j) == 1:
                mpz_mul(i, i, i)
                mpz_mul(j, j, j)
                mpz_mul(j, j, left.D.value)
                test = mpz_cmp(i, j)
            else:
                test = 1

        else:
            if mpz_sgn(j) == -1:
                mpz_mul(i, i, i)
                mpz_mul(j, j, j)
                mpz_mul(j, j, left.D.value)
                test = mpz_cmp(j, i)
            else:
                test = -1

        mpz_clear(i)
        mpz_clear(j)
        return rich_to_bool_sgn(op, test)

    cpdef int _cmp_(left, _right) except -2:
        """
        Comparisons of elements.

        When there is a real embedding defined, the comparisons uses comparison
        induced from the reals. Otherwise, comparison is a lexicographic
        comparison on coefficients.

        EXAMPLES:

        Two examples from the same number field with its two possible real
        embeddings::

            sage: K.<phi> = NumberField(x^2-x-1, 'phi', embedding=1.618)
            sage: phi > 0
            True
            sage: -phi > 0
            False
            sage: phi - 3 == 2*phi + 1
            False
            sage: fibonacci(10)*phi < fibonacci(11)
            True
            sage: RDF(fibonacci(10)*phi)
            88.99186938124421
            sage: fibonacci(11)
            89
            sage: l = [-2, phi+3, 2*phi-1, 2*phi-5, 0, -phi+2, fibonacci(20)*phi - fibonacci(21)]
            sage: l.sort()
            sage: l
            [-2, 2*phi - 5, 6765*phi - 10946, 0, -phi + 2, 2*phi - 1, phi + 3]
            sage: map(RDF, l)
            [-2.0, -1.7639320225002102, -6.610696073039435e-05, 0.0, 0.3819660112501051, 2.23606797749979, 4.618033988749895]

        ::

            sage: L.<psi> = NumberField(x^2-x-1, 'psi', embedding=-0.618)
            sage: psi < 0
            True
            sage: 2*psi + 3 == 2*psi + 3
            True
            sage: fibonacci(10)*psi < -fibonacci(9)
            False
            sage: RDF(fibonacci(10)*psi)
            -33.99186938124422
            sage: fibonacci(9)
            34
            sage: l = [-1, psi, 0, fibonacci(20)*psi + fibonacci(19), 3*psi+2]
            sage: l.sort()
            sage: l
            [-1, psi, 0, 6765*psi + 4181, 3*psi + 2]
            sage: map(RDF, l)
            [-1.0, -0.6180339887498949, 0.0, 6.610696073039435e-05, 0.1458980337503153]

        For a field with no specified embedding the comparison uses the standard
        embedding::

            sage: K.<sqrt2> = NumberField(x^2-2, 'sqrt2')
            sage: sqrt2 > 1 and sqrt2 < 2
            True

        The following examples illustrate the same behavior for a complex
        quadratic field::

            sage: K.<i> = QuadraticField(-1)
            sage: l = [-2, i-3, 2*i-2, 2*i+2, 5*i, 1-3*i, -1+i, 1]
            sage: l.sort()
            sage: l
            [i - 3, -2, 2*i - 2, i - 1, 5*i, -3*i + 1, 1, 2*i + 2]
            sage: map(CDF, l)
            [-3.0 + 1.0*I, -2.0, -2.0 + 2.0*I, -1.0 + 1.0*I, 5.0*I, 1.0 - 3.0*I, 1.0, 2.0 + 2.0*I]
            sage: map(CDF, l) == sorted(map(CDF, l))
            True
        """
        cdef NumberFieldElement_quadratic right = <NumberFieldElement_quadratic> _right
        cdef int test

        if left == right:
             return 0
        if left > right:
            return 1
        return -1

    def continued_fraction_list(self):
        r"""
        Return the preperiod and the period of the continued fraction expansion
        of ``self``.

        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: sqrt2.continued_fraction_list()
            ((1,), (2,))
            sage: (1/2+sqrt2/3).continued_fraction_list()
            ((0, 1, 33), (1, 32))

        For rational entries a pair of tuples is also returned but the second
        one is empty::

            sage: K(123/567).continued_fraction_list()
            ((0, 4, 1, 1, 1, 1, 3, 2), ())
        """
        cdef NumberFieldElement_quadratic x

        if mpz_sgn(self.b) == 0:
            return tuple(Rational(self).continued_fraction_list()),()

        if mpz_sgn(self.D.value) < 0:
            raise ValueError("the method is only available for positive discriminant")

        x = self
        orbit = []
        quots = []
        while x not in orbit:
            quots.append(x.floor())
            orbit.append(x)
            x = ~(x - quots[-1])

        i = orbit.index(x)

        return tuple(quots[:i]), tuple(quots[i:])

    def continued_fraction(self):
        r"""
        Return the (finite or ultimately periodic) continued fraction of ``self``.

        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: cf = sqrt2.continued_fraction(); cf
            [1; (2)*]
            sage: cf.n()
            1.41421356237310
            sage: sqrt2.n()
            1.41421356237310
            sage: cf.value()
            sqrt2

            sage: (sqrt2/3 + 1/4).continued_fraction()
            [0; 1, (2, 1, 1, 2, 3, 2, 1, 1, 2, 5, 1, 1, 14, 1, 1, 5)*]
        """
        t1,t2 = self.continued_fraction_list()
        from sage.rings.continued_fraction import ContinuedFraction_periodic
        return ContinuedFraction_periodic(t1,t2)

#########################################################
# Arithmetic
#########################################################

    cdef void _reduce_c_(self):
        r"""
        Reduces into canonical form.

        WARNING: this mutates self.
        """
        cdef mpz_t gcd
        # cancel out common factors
        mpz_init(gcd)
        mpz_gcd(gcd, self.a, self.denom)
        mpz_gcd(gcd, gcd, self.b)
        if mpz_cmp_si(gcd, 1): # != 0 (i.e. it is not 1)
            mpz_divexact(self.a, self.a, gcd)
            mpz_divexact(self.b, self.b, gcd)
            mpz_divexact(self.denom, self.denom, gcd)
        # make denominator positive
        if mpz_sgn(self.denom) < 0:
            mpz_neg(self.denom, self.denom)
            mpz_neg(self.a, self.a)
            mpz_neg(self.b, self.b)
        mpz_clear(gcd)


    cpdef _add_(self, other_m):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2-5)
            sage: K.discriminant()
            5
            sage: a+a # indirect doctest
            2*a
            sage: s = (a+2)/6; s
            1/6*a + 1/3
            sage: s+a
            7/6*a + 1/3
            sage: s+10
            1/6*a + 31/3
            sage: s+(2*a+5)/7
            19/42*a + 22/21
            sage: s+(1+a)/2
            2/3*a + 5/6
            sage: s+(1+a)/8
            7/24*a + 11/24
            sage: s+(a+5)/6
            1/3*a + 7/6
            sage: (a/3+2/3) + (2*a/3+1/3)
            a + 1
        """
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        cdef mpz_t gcd, tmp
        if mpz_cmp(self.denom, other.denom) == 0:
            mpz_add(res.a, self.a, other.a)
            mpz_add(res.b, self.b, other.b)
            mpz_set(res.denom, self.denom)
        else:
            mpz_init(gcd)
            mpz_gcd(gcd, self.denom, other.denom)
            if mpz_cmp_ui(gcd, 1) == 0:
                mpz_mul(res.a, self.a, other.denom)
                mpz_addmul(res.a, self.denom, other.a)
                mpz_mul(res.b, self.b, other.denom)
                mpz_addmul(res.b, self.denom, other.b)
                mpz_mul(res.denom, self.denom, other.denom)
            else:
                mpz_init(tmp)
                mpz_divexact(tmp, other.denom, gcd)
                mpz_mul(res.a, self.a, tmp)
                mpz_mul(res.b, self.b, tmp)
                mpz_divexact(tmp, self.denom, gcd)
                mpz_addmul(res.a, other.a, tmp)
                mpz_addmul(res.b, other.b, tmp)
                mpz_mul(res.denom, other.denom, tmp)
                mpz_clear(tmp)
            mpz_clear(gcd)
        res._reduce_c_()
        return res


    cpdef _sub_(self, other_m):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2-13)
            sage: b = (a-3)/10; b # indirect doctest
            1/10*a - 3/10
            sage: b-1
            1/10*a - 13/10
            sage: b-a
            -9/10*a - 3/10
            sage: b-1/2
            1/10*a - 4/5
            sage: b-a/15
            1/30*a - 3/10
        """
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        cdef mpz_t gcd, tmp
        if mpz_cmp(self.denom, other.denom) == 0:
            mpz_sub(res.a, self.a, other.a)
            mpz_sub(res.b, self.b, other.b)
            mpz_set(res.denom, self.denom)
        else:
            mpz_init(gcd)
            mpz_gcd(gcd, self.denom, other.denom)
            if mpz_cmp_ui(gcd, 1) == 0:
                mpz_mul(res.a, self.a, other.denom)
                mpz_submul(res.a, self.denom, other.a)
                mpz_mul(res.b, self.b, other.denom)
                mpz_submul(res.b, self.denom, other.b)
                mpz_mul(res.denom, self.denom, other.denom)
            else:
                mpz_init(tmp)
                mpz_divexact(tmp, other.denom, gcd)
                mpz_mul(res.a, self.a, tmp)
                mpz_mul(res.b, self.b, tmp)
                mpz_divexact(tmp, self.denom, gcd)
                mpz_submul(res.a, other.a, tmp)
                mpz_submul(res.b, other.b, tmp)
                mpz_mul(res.denom, other.denom, tmp)
                mpz_clear(tmp)
            mpz_clear(gcd)
        res._reduce_c_()
        return res


    def __neg__(self):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2+163)
            sage: -a
            -a
            sage: -(a+4)
            -a - 4
            sage: b = (a-3)/2
            sage: -b
            -1/2*a + 3/2
        """
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        mpz_neg(res.a, self.a)
        mpz_neg(res.b, self.b)
        mpz_set(res.denom, self.denom)
        return res


    cpdef _mul_(self, other_m):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+23)
            sage: a*a # indirect doctest
            -23
            sage: (a+1)*(a-1)
            -24
            sage: (a+1)*(a+2)
            3*a - 21
            sage: (a+1)/2 * (a+2)
            3/2*a - 21/2
            sage: (a+1)/2 * (a+2)/3
            1/2*a - 7/2
            sage: (2*a+4) * (3*a)/2
            6*a - 69

        Verify Karatsuba
            sage: K.<a> = NumberField(x^2-41)
            sage: (10^1000 * (a+1)) * K(2+3*a) == 10^1000 * ((a+1) * K(2+3*a))
            True
        """
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        cdef mpz_t tmp

        if mpz_size(self.a) + mpz_size(self.b) < 8: # could I use a macro instead?
            # Do it the traditional way
            mpz_mul(res.a, self.b, other.b)
            mpz_mul(res.a, res.a, self.D.value)
            mpz_addmul(res.a, self.a, other.a)

            mpz_mul(res.b, self.a, other.b)
            mpz_addmul(res.b, self.b, other.a)

        else:
            # Karatsuba
            sig_on()
            mpz_init(tmp)
            mpz_add(res.a, self.a, self.b) # using res.a as tmp
            mpz_add(tmp, other.a, other.b)
            mpz_mul(res.b, res.a, tmp) # res.b = (self.a + self.b)(other.a + other.b)

            mpz_mul(res.a, self.a, other.a)
            mpz_sub(res.b, res.b, res.a)
            mpz_mul(tmp, self.b, other.b)
            mpz_sub(res.b, res.b, tmp)
            mpz_mul(tmp, tmp, self.D.value)
            mpz_add(res.a, res.a, tmp)
            mpz_clear(tmp)
            sig_off()

        mpz_mul(res.denom, self.denom, other.denom)
        res._reduce_c_()
        return res


    cpdef _rmul_(self, RingElement _c):
        """
        EXAMPLE:
            sage: K.<a> = NumberField(x^2+43)
            sage: (1+a)*3 # indirect doctest
            3*a + 3
        """
        cdef Rational c =  <Rational>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        mpz_mul(res.a, self.a, mpq_numref(c.value))
        mpz_mul(res.b, self.b, mpq_numref(c.value))
        mpz_mul(res.denom, self.denom, mpq_denref(c.value))
        res._reduce_c_()
        return res


    cpdef _lmul_(self, RingElement _c):
        """
        EXAMPLE:
            sage: K.<a> = NumberField(x^2+43)
            sage: 5*(a-1/5) # indirect doctest
            5*a - 1
        """
        cdef Rational c =  <Rational>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        mpz_mul(res.a, self.a, mpq_numref(c.value))
        mpz_mul(res.b, self.b, mpq_numref(c.value))
        mpz_mul(res.denom, self.denom, mpq_denref(c.value))
        res._reduce_c_()
        return res


    cpdef _div_(self, other):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2-5)
            sage: 2/a # indirect doctest
            2/5*a
            sage: (a+2)/(a+1)
            1/4*a + 3/4
            sage: (a+1)*(a+2)/(a+1)
            a + 2
            sage: (a+1/3)*(5*a+2/7)/(a+1/3)
            5*a + 2/7
        """
        return self * ~other


    def __invert__(self):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2-5)
            sage: ~a
            1/5*a
            sage: ~(a+1)
            1/4*a - 1/4
            sage: (a-1)*(a+1)
            4
            sage: b = ~(5*a-3); b
            5/116*a + 3/116
            sage: b*(5*a-3)
            1
            sage: b = ~((3*a-2)/7); b
            21/41*a + 14/41
            sage: (3*a-2)/7 * b
            1

        This fixes ticket :trac:`9357`::

            sage: K.<a> = NumberField(x^2+1)
            sage: d = K(0)
            sage: ~d
            Traceback (most recent call last):
            ...
            ZeroDivisionError
            sage: K.random_element()/d
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        if mpz_cmp_ui(self.a, 0) == 0 and mpz_cmp_ui(self.b, 0) == 0:
            raise ZeroDivisionError
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        cdef mpz_t tmp, gcd
        mpz_init(tmp)
        mpz_init(gcd)

        mpz_gcd(gcd, self.a, self.b)
        if mpz_cmp_si(gcd, 1): # != 0 (i.e. it is not 1)
            # cancel out g (g(a'-b'd)) / (g^2 (a'^2-b'^2d^2))
            mpz_divexact(res.a, self.a, gcd)
            mpz_divexact(res.b, self.b, gcd)
            mpz_neg(res.b, res.b)
        else:
            mpz_set(res.a, self.a)
            mpz_neg(res.b, self.b)

        mpz_pow_ui(res.denom, res.a, 2)
        mpz_pow_ui(tmp, res.b, 2)
        mpz_mul(tmp, tmp, self.D.value)
        mpz_sub(res.denom, res.denom, tmp)
        # need to multiply the leftover g back in
        mpz_mul(res.denom, res.denom, gcd)

        mpz_mul(res.a, res.a, self.denom)
        mpz_mul(res.b, res.b, self.denom)

        mpz_clear(tmp)
        mpz_clear(gcd)

        res._reduce_c_()
        return res


    cpdef NumberFieldElement galois_conjugate(self):
        """
        Return the image of this element under action of the nontrivial
        element of the Galois group of this field.

        EXAMPLES::

            sage: K.<a> = QuadraticField(23)
            sage: a.galois_conjugate()
            -a

            sage: K.<a> = NumberField(x^2 - 5*x + 1)
            sage: a.galois_conjugate()
            -a + 5
            sage: b = 5*a + 1/3
            sage: b.galois_conjugate()
            -5*a + 76/3
            sage: b.norm() ==  b * b.galois_conjugate()
            True
            sage: b.trace() ==  b + b.galois_conjugate()
            True
        """
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        mpz_set(res.a, self.a)
        mpz_neg(res.b, self.b)
        mpz_set(res.denom, self.denom)
        return res

#################################################################################
# We must override everything that makes uses of self.__numerator/__denominator
#################################################################################

    def __hash__(self):
        """
        Return hash of this number field element.

        For elements in `\ZZ` or `\QQ` the hash coincides with the one in the
        native `\ZZ` or `\QQ`.

        EXAMPLES::

            sage: L.<a> = QuadraticField(-7)
            sage: hash(a)
            42082631
            sage: hash(L(1))
            1
            sage: hash(L(-3))
            -3
            sage: hash(L(-32/118)) == hash(-32/118)
            True
        """
        # 1. compute the hash of a/denom as if it was a rational
        # (see the corresponding code in sage/rings/rational.pyx)
        cdef Py_hash_t n = mpz_pythonhash(self.a)
        cdef Py_hash_t d = mpz_pythonhash(self.denom)
        cdef Py_hash_t h = n + (d - 1) * <Py_hash_t>(7461864723258187525)

        # 2. mix the hash together with b
        h += 42082631 * mpz_pythonhash(self.b)
        return h

    def __nonzero__(self):
        """
        Check whether this element is not zero.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+163)
            sage: not a
            False
            sage: not (a-a)
            True
        """
        return mpz_cmp_ui(self.a, 0) != 0 or mpz_cmp_ui(self.b, 0) != 0

    def _integer_(self, Z=None):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: (a+1-a)._integer_()
            1
            sage: (a+1/2-a)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 1/2 to an integer
        """
        cdef Integer res
        if mpz_cmp_ui(self.b, 0) != 0 or mpz_cmp_ui(self.denom, 1) != 0:
            raise TypeError("Unable to coerce %s to an integer" % self)
        else:
            res = PY_NEW(Integer)
            mpz_set(res.value, self.a)
            return res

    def _rational_(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: (a+1/2-a)._rational_()
            1/2
            sage: (a+1/2)._rational_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a + 1/2 to a rational
        """
        cdef Rational res
        if mpz_cmp_ui(self.b, 0)!=0:
            raise TypeError("Unable to coerce %s to a rational" % self)
        else:
            res = <Rational>Rational.__new__(Rational)
            mpz_set(mpq_numref(res.value), self.a)
            mpz_set(mpq_denref(res.value), self.denom)
            mpq_canonicalize(res.value)
            return res

    def _algebraic_(self, parent):
        r"""
        Convert this element to an algebraic number, if possible.

        EXAMPLES::

            sage: NF.<i> = QuadraticField(-1)
            sage: QQbar(1+i)
            I + 1
            sage: NF.<sqrt3> = QuadraticField(2)
            sage: AA(sqrt3)
            1.414213562373095?
        """
        import sage.rings.qqbar as qqbar
        if (parent is qqbar.QQbar
                and list(self._parent.polynomial()) == [1, 0, 1]):
            # AlgebraicNumber.__init__ does a better job than
            # NumberFieldElement._algebraic_ in this case, but
            # QQbar._element_constructor_ calls the latter first.
            return qqbar.AlgebraicNumber(self)
        else:
            return NumberFieldElement._algebraic_(self, parent)

    cpdef bint is_one(self):
        r"""
        Check whether this number field element is `1`.

        EXAMPLES::

            sage: K = QuadraticField(-2)
            sage: K(1).is_one()
            True
            sage: K(-1).is_one()
            False
            sage: K(2).is_one()
            False
            sage: K(0).is_one()
            False
            sage: K(1/2).is_one()
            False
            sage: K.gen().is_one()
            False
        """
        return mpz_cmp_ui(self.a, 1) == 0 and \
               mpz_cmp_ui(self.b, 0) == 0 and \
               mpz_cmp_ui(self.denom, 1) == 0

    cpdef bint is_rational(self):
        r"""
        Check whether this number field element is a rational number.

        .. SEEALSO:

            - :meth:`is_integer` to test if this element is an integer
            - :meth:`is_integral` to test if this element is an algebraic integer

        EXAMPLES::

            sage: K.<sqrt3> = QuadraticField(3)
            sage: sqrt3.is_rational()
            False
            sage: (sqrt3-1/2).is_rational()
            False
            sage: K(0).is_rational()
            True
            sage: K(-12).is_rational()
            True
            sage: K(1/3).is_rational()
            True
        """
        return mpz_cmp_ui(self.b, 0) == 0

    def is_integer(self):
        r"""
        Check whether this number field element is an integer.

        .. SEEALSO:

            - :meth:`is_rational` to test if this element is a rational number
            - :meth:`is_integral` to test if this element is an algebraic integer

        EXAMPLES::

            sage: K.<sqrt3> = QuadraticField(3)
            sage: sqrt3.is_integer()
            False
            sage: (sqrt3-1/2).is_integer()
            False
            sage: K(0).is_integer()
            True
            sage: K(-12).is_integer()
            True
            sage: K(1/3).is_integer()
            False
        """
        return mpz_cmp_ui(self.b, 0) == mpz_cmp_ui(self.denom, 1) == 0

    def real(self):
        r"""
        Returns the real part of self, which is either self (if self lives
        it a totally real field) or a rational number.

        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: sqrt2.real()
            sqrt2
            sage: K.<a> = QuadraticField(-3)
            sage: a.real()
            0
            sage: (a + 1/2).real()
            1/2
            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: a.real()
            -1/2
            sage: parent(a.real())
            Rational Field
            sage: K.<i> = QuadraticField(-1)
            sage: i.real()
            0
        """
        cdef Rational res
        if mpz_sgn(self.D.value) > 0:
            return self  # totally real
        else:
            res = <Rational>Rational.__new__(Rational)
            mpz_set(mpq_numref(res.value), self.a)
            mpz_set(mpq_denref(res.value), self.denom)
            mpq_canonicalize(res.value)
            return res

    def imag(self):
        r"""
        Returns the imaginary part of self.

        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: sqrt2.imag()
            0
            sage: parent(sqrt2.imag())
            Rational Field

            sage: K.<i> = QuadraticField(-1)
            sage: i.imag()
            1
            sage: parent(i.imag())
            Rational Field

            sage: K.<a> = NumberField(x^2 + x + 1, embedding=CDF.0)
            sage: a.imag()
            1/2*sqrt3
            sage: a.real()
            -1/2
            sage: SR(a)
            1/2*I*sqrt(3) - 1/2
            sage: bool(I*a.imag() + a.real() == a)
            True

        TESTS::

            sage: K.<a> = QuadraticField(-9, embedding=-CDF.0)
            sage: a.imag()
            -3
            sage: parent(a.imag())
            Rational Field

        """
        if mpz_sgn(self.D.value) > 0:
            return Rational.__new__(Rational)  # = 0
        embedding =  self._parent.coerce_embedding()
        cdef Integer negD = -self.D
        cdef NumberFieldElement_quadratic q = <NumberFieldElement_quadratic>self._new()
        mpz_set_ui(q.b, 1)
        mpz_set_ui(q.denom, 1)
        from sage.rings.complex_double import CDF
        cdef Rational res
        if mpz_cmp_ui(negD.value, 1) == 0 or negD.is_square():
            # D = -1 is the most common case we'll see here
            if embedding is None:
                raise ValueError("Embedding must be specified.")
            res = <Rational>Rational.__new__(Rational)
            if mpz_cmp_ui(negD.value, 1) == 0:
                mpz_set(mpq_numref(res.value), self.b)
            else:
                mpz_sqrt(mpq_numref(res.value), negD.value)
                mpz_mul(mpq_numref(res.value), mpq_numref(res.value), self.b)
            mpz_set(mpq_denref(res.value), self.denom)
            mpq_canonicalize(res.value)
            if not self.standard_embedding:
                mpq_neg(res.value, res.value)
            return res
        else:
            # avoid circular import
            if embedding is None:
                from number_field import NumberField
                K = NumberField(QQ['x'].gen()**2 - negD, 'sqrt%s' % negD)
            else:
                from number_field import QuadraticField
                K = QuadraticField(negD, 'sqrt%s' % negD)
            q = K(0)
            mpz_set(q.denom, self.denom)
            if self.standard_embedding:
                mpz_set(q.b, self.b)
            else:
                mpz_neg(q.b, self.b)
            return q

    cpdef list _coefficients(self):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2+41)
            sage: a._coefficients()
            [0, 1]
            sage: K.<a> = NumberField(x^2+x+41)
            sage: a._coefficients()
            [0, 1]
            sage: b = 3*a+1/5
            sage: b._coefficients()
            [1/5, 3]
        """
        # In terms of the generator...
        cdef NumberFieldElement_quadratic gen = self.number_field().gen()  # should this be cached?
        cdef Rational const = <Rational>Rational.__new__(Rational)
        cdef Rational lin = <Rational>Rational.__new__(Rational)
        ad, bd = self.parts()
        if not self:
            return []
        if not bd:
            return [ad]
        if gen.is_sqrt_disc():
            return [ad,bd]
        else:
            alpha, beta = gen.parts()
            scale = bd/beta
            return [ad - scale*alpha, scale]

    def denominator(self):
        """
        Return the denominator of self. This is the LCM of the denominators of
        the coefficients of self, and thus it may well be `> 1` even when the
        element is an algebraic integer.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+x+41)
            sage: a.denominator()
            1
            sage: b = (2*a+1)/6
            sage: b.denominator()
            6
            sage: K(1).denominator()
            1
            sage: K(1/2).denominator()
            2
            sage: K(0).denominator()
            1

            sage: K.<a> = NumberField(x^2 - 5)
            sage: b = (a + 1)/2
            sage: b.denominator()
            2
            sage: b.is_integral()
            True
        """
        # In terms of the generator...
        cdef NumberFieldElement_quadratic gen = self.number_field().gen()  # should this be cached?
        cdef Integer denom
        if gen.is_sqrt_disc():
            denom = PY_NEW(Integer)
            mpz_set(denom.value, self.denom)
            return denom
        else:
            c = self._coefficients()
            if len(c) == 2:
                const, lin = c
            elif len(c) == 1:
                const = c[0]
                lin = Rational(0)
            else:
                const = lin = Rational(0)
            return const.denominator().lcm(lin.denominator())

    def numerator(self):
        """
        Return self*self.denominator().

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+x+41)
            sage: b = (2*a+1)/6
            sage: b.denominator()
            6
            sage: b.numerator()
            2*a + 1
        """
        return self*self.denominator()


#########################################################
# Some things are so much easier to compute
#########################################################

    def trace(self):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2+x+41)
            sage: a.trace()
            -1
            sage: a.matrix()
            [  0   1]
            [-41  -1]

        The trace is additive::

            sage: K.<a> = NumberField(x^2+7)
            sage: (a+1).trace()
            2
            sage: K(3).trace()
            6
            sage: (a+4).trace()
            8
            sage: (a/3+1).trace()
            2
        """
        # trace = 2*self.a / self.denom
        cdef Rational res = <Rational>Rational.__new__(Rational)
        if mpz_odd_p(self.denom):
            mpz_mul_2exp(mpq_numref(res.value), self.a, 1)
            mpz_set(mpq_denref(res.value), self.denom)
        else:
            mpz_set(mpq_numref(res.value), self.a)
            mpz_divexact_ui(mpq_denref(res.value), self.denom, 2)
        mpq_canonicalize(res.value)
        return res


    def norm(self, K=None):
        """
        Return the norm of self. If the second argument is None, this is the
        norm down to `\QQ`. Otherwise, return the norm down to K (which had
        better be either `\QQ` or this number field).

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-x+3)
            sage: a.norm()
            3
            sage: a.matrix()
            [ 0  1]
            [-3  1]
            sage: K.<a> = NumberField(x^2+5)
            sage: (1+a).norm()
            6

        The norm is multiplicative::

            sage: K.<a> = NumberField(x^2-3)
            sage: a.norm()
            -3
            sage: K(3).norm()
            9
            sage: (3*a).norm()
            -27

        We test that the optional argument is handled sensibly::

            sage: (3*a).norm(QQ)
            -27
            sage: (3*a).norm(K)
            3*a
            sage: (3*a).norm(CyclotomicField(3))
            Traceback (most recent call last):
            ...
            ValueError: no way to embed L into parent's base ring K

        """
        cdef Rational res = <Rational>Rational.__new__(Rational)

        if K is None or K == QQ:
        # norm = (a^2 - d b^2) / self.denom^2
            mpz_pow_ui(mpq_numref(res.value), self.a, 2)
            mpz_pow_ui(mpq_denref(res.value), self.b, 2)  # use as temp
            mpz_mul(mpq_denref(res.value), mpq_denref(res.value), self.D.value)
            mpz_sub(mpq_numref(res.value), mpq_numref(res.value), mpq_denref(res.value))
            mpz_pow_ui(mpq_denref(res.value), self.denom, 2)
            mpq_canonicalize(res.value)
            return res
        else:
            return NumberFieldElement.norm(self, K)


    def is_integral(self):
        r"""
        Returns whether this element is an algebraic integer.

        TESTS::

            sage: K.<a> = QuadraticField(-1)
            sage: a.is_integral()
            True
            sage: K(1).is_integral()
            True
            sage: K(1/2).is_integral()
            False
            sage: K(a/2).is_integral()
            False
            sage: K((a+1)/2).is_integral()
            False
            sage: K(a/3).is_integral()
            False

            sage: K.<a> = QuadraticField(-3)
            sage: a.is_integral()
            True
            sage: K(1).is_integral()
            True
            sage: K(1/2).is_integral()
            False
            sage: K(a/2).is_integral()
            False
            sage: ((a+1)/2).is_integral()
            True
        """
        if mpz_cmp_ui(self.denom, 1) == 0:
            return True
        else:
            return self.norm().denom() == 1 and self.trace().denom() == 1

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this element over `\QQ`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-x+13)
            sage: a.charpoly()
            x^2 - x + 13
            sage: b = 3-a/2
            sage: f = b.charpoly(); f
            x^2 - 11/2*x + 43/4
            sage: f(b)
            0
        """
        R = QQ[var]
        return R([self.norm(), -self.trace(), 1])

    def minpoly(self, var='x'):
        r"""
        The minimal polynomial of this element over `\QQ`.

        INPUT:

        -  ``var`` -- the minimal polynomial is defined over a polynomial ring
           in a variable with this name. If not specified this defaults to
           ``x``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+13)
            sage: a.minpoly()
            x^2 + 13
            sage: a.minpoly('T')
            T^2 + 13
            sage: (a+1/2-a).minpoly()
            x - 1/2
        """
        if self.is_rational():
            R = QQ[var]
            return R([-self._rational_(), 1])
        else:
            return self.charpoly(var)

    def __abs__(self):
        """
        EXAMPLES::

            sage: K.<a> = QuadraticField(2, 'a', embedding=-1.4142)
            sage: abs(a)    # indirect test
            -a
            sage: abs(a+1)  # indirect test
            -a - 1
            sage: abs(a+2)  # indirect test
            a + 2

            sage: K.<a> = NumberField(x^2+1, embedding=CDF.gen())
            sage: abs(a+1)
            sqrt(2)
        """
        if mpz_sgn(self.D.value) == 1:
            if self.sign() >= 0:
                return self
            return -self

        # doing that way the parent is the symbolic ring (or IntegerRing if the
        # norm of self is a square). On the other hand, it is coherent with
        # sage.rings.integer.Integer.sqrt
        return (self.real()**2 + self.imag()**2).sqrt()

    def floor(self):
        r"""
        Returns the floor of x.

        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2,name='sqrt2')
            sage: sqrt2.floor()
            1
            sage: (-sqrt2).floor()
            -2
            sage: (13/197 + 3702/123*sqrt2).floor()
            42
            sage: (13/197-3702/123*sqrt2).floor()
            -43

        TESTS::

            sage: K2.<sqrt2> = QuadraticField(2)
            sage: K3.<sqrt3> = QuadraticField(3)
            sage: K5.<sqrt5> = QuadraticField(5)
            sage: for _ in xrange(100):
            ....:    a = QQ.random_element(1000,20)
            ....:    b = QQ.random_element(1000,20)
            ....:    assert floor(a+b*sqrt(2.)) == floor(a+b*sqrt2)
            ....:    assert floor(a+b*sqrt(3.)) == floor(a+b*sqrt3)
            ....:    assert floor(a+b*sqrt(5.)) == floor(a+b*sqrt5)

            sage: K = QuadraticField(-2)
            sage: l = [K(52), K(-3), K(43/12), K(-43/12)]
            sage: [x.floor() for x in l]
            [52, -3, 3, -4]
        """
        cdef mpz_t x
        cdef Integer result

        if mpz_sgn(self.b) == 0:
            mpz_init_set(x,self.a)
            mpz_fdiv_q(x,x,self.denom)
            result = PY_NEW(Integer)
            mpz_set(result.value,x)
            mpz_clear(x)
            return result

        if not mpz_sgn(self.D.value) == 1:
            raise ValueError("floor is not defined for complex quadratic number field")

        mpz_init(x)
        mpz_mul(x,self.b,self.b)
        mpz_mul(x,x,self.D.value)
        mpz_sqrt(x,x)
        if mpz_sgn(self.b) == -1:
            if self.standard_embedding:
                mpz_neg(x,x)
                mpz_sub_ui(x,x,1)
        elif not self.standard_embedding:
                mpz_neg(x,x)
                mpz_sub_ui(x,x,1)

        mpz_add(x,x,self.a)    # here x = a + floor(sqrt(b^2 D)) or a + floor(-sqrt(b^2 D))
        mpz_fdiv_q(x,x,self.denom)
        result = PY_NEW(Integer)
        mpz_set(result.value,x)
        mpz_clear(x)
        return result

    def ceil(self):
        r"""
        Returns the ceil.

        EXAMPLES::

            sage: K.<sqrt7> = QuadraticField(7, name='sqrt7')
            sage: sqrt7.ceil()
            3
            sage: (-sqrt7).ceil()
            -2
            sage: (1022/313*sqrt7 - 14/23).ceil()
            9

        TESTS::

            sage: K2.<sqrt2> = QuadraticField(2)
            sage: K3.<sqrt3> = QuadraticField(3)
            sage: K5.<sqrt5> = QuadraticField(5)
            sage: for _ in xrange(100):
            ....:    a = QQ.random_element(1000,20)
            ....:    b = QQ.random_element(1000,20)
            ....:    assert ceil(a+b*sqrt(2.)) == ceil(a+b*sqrt2)
            ....:    assert ceil(a+b*sqrt(3.)) == ceil(a+b*sqrt3)
            ....:    assert ceil(a+b*sqrt(5.)) == ceil(a+b*sqrt5)

            sage: K = QuadraticField(-2)
            sage: l = [K(52), K(-3), K(43/12), K(-43/12)]
            sage: [x.ceil() for x in l]
            [52, -3, 4, -3]
        """
        x = self.floor()
        if mpz_sgn(self.b) == 0 and mpz_cmp_ui(self.denom,1) == 0:
            return x
        return x+1


cdef class OrderElement_quadratic(NumberFieldElement_quadratic):
    """
    Element of an order in a quadratic field.

    EXAMPLES::

        sage: K.<a> = NumberField(x^2 + 1)
        sage: O2 = K.order(2*a)
        sage: w = O2.1; w
        2*a
        sage: parent(w)
        Order in Number Field in a with defining polynomial x^2 + 1
    """
    def __init__(self, order, f):
        r"""
        Standard initialisation function.

        EXAMPLE::

            sage: OK.<y> = EquationOrder(x^2 + 5)
            sage: v = OK.1 # indirect doctest
            sage: type(v)
            <type 'sage.rings.number_field.number_field_element_quadratic.OrderElement_quadratic'>
        """
        K = order.number_field()
        NumberFieldElement_quadratic.__init__(self, K, f)
        (<Element>self)._parent = order

    def norm(self):
        """
        The norm of an element of the ring of integers is an Integer.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 3)
            sage: O2 = K.order(2*a)
            sage: w = O2.gen(1); w
            2*a
            sage: w.norm()
            12
            sage: parent(w.norm())
            Integer Ring
        """
        return ZZ(NumberFieldElement_quadratic.norm(self))

    def trace(self):
        """
        The trace of an element of the ring of integers is an Integer.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: R = K.ring_of_integers()
            sage: b = R((1+a)/2)
            sage: b.trace()
            1
            sage: parent(b.trace())
            Integer Ring
        """
        return ZZ(NumberFieldElement_quadratic.trace(self))

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this element, which is over `\ZZ`
        because this element is an algebraic integer.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: R = K.ring_of_integers()
            sage: b = R((5+a)/2)
            sage: f = b.charpoly('x'); f
            x^2 - 5*x + 5
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: f(b)
            0
        """
        R = ZZ[var]
        return R([self.norm(), -self.trace(), 1])

    def minpoly(self, var='x'):
        r"""
        The minimal polynomial of this element over `\ZZ`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 163)
            sage: R = K.ring_of_integers()
            sage: f = R(a).minpoly('x'); f
            x^2 + 163
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: R(5).minpoly()
            x - 5
        """
        if self.is_rational():
            R = ZZ[var]
            return R([-self._rational_(), 1])
        else:
            return self.charpoly()

    cdef number_field(self):
        # So few functions actually use self.number_field() for quadratic elements, so
        # it is better *not* to return a cached value (since the call to _parent.number_field())
        # is expensive.
        return self._parent.number_field()

    # We must override these since the basering is now ZZ not QQ.
    cpdef _rmul_(self, RingElement _c):
        """
        EXAMPLE:
            sage: K.<a> = NumberField(x^2-27)
            sage: R = K.ring_of_integers()
            sage: aa = R.gen(1); aa
            1/3*a
            sage: 5 * aa # indirect doctest
            5/3*a
        """
        cdef Integer c = <Integer>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        mpz_mul(res.a, self.a, c.value)
        mpz_mul(res.b, self.b, c.value)
        mpz_set(res.denom, self.denom)
        res._reduce_c_()
        return res


    cpdef _lmul_(self, RingElement _c):
        """
        EXAMPLE:
            sage: K.<a> = NumberField(x^2+43)
            sage: R = K.ring_of_integers()
            sage: aa = R.gen(0); aa
            1/2*a + 1/2
            sage: aa*3 # indirect doctest
            3/2*a + 3/2
        """
        cdef Integer c = <Integer>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        mpz_mul(res.a, self.a, c.value)
        mpz_mul(res.b, self.b, c.value)
        mpz_set(res.denom, self.denom)
        res._reduce_c_()
        return res

    cpdef _div_(self, other):
        r"""
        Implement division, checking that the result has the
        right parent. It's not so crucial what the parent actually
        is, but it is crucial that the returned value really is an
        element of its supposed parent! This fixes trac \#4190.

        EXAMPLES::

            sage: K = NumberField(x^2 - 17, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: (17/a).parent() is K # indirect doctest
            True
            sage: 17/a in OK
            True
            sage: (17/(2*a)).parent() is OK
            False
            sage: (17/(2*a)) in OK
            False
            sage: (17/(2*a)).parent() is K
            True
        """
        cdef NumberFieldElement_quadratic x
        x = NumberFieldElement_quadratic._div_(self, other)
        return self._parent.number_field()(x)

    def __invert__(self):
        r"""
        Implement inversion, checking that the return value has the right parent.
        See trac \#4190.

        EXAMPLES::

            sage: K = NumberField(x^2 -x + 2, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: (~a).parent() is K
            True
            sage: (~a) in OK
            False
            sage: a**(-1) in OK
            False
        """
        return self._parent.number_field()(NumberFieldElement_quadratic.__invert__(self))

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the given ideal.

        INPUT:


        -  ``I`` - may be an ideal of self.parent(), or an
           element or list of elements of self.parent() generating a nonzero
           ideal. A ValueError is raised if I is non-integral or is zero.
           A ZeroDivisionError is raised if I + (x) != (1).


        EXAMPLES::

            sage: OE.<w> = EquationOrder(x^2 - x + 2)
            sage: w.inverse_mod(13) == 6*w - 6
            True
            sage: w*(6*w - 6) - 1
            -13
            sage: w.inverse_mod(13).parent() == OE
            True
            sage: w.inverse_mod(2*OE)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: w is not invertible modulo Fractional ideal (2)
        """
        R = self.parent()
        return R(_inverse_mod_generic(self, I))

cdef class Z_to_quadratic_field_element(Morphism):
    """
    Morphism that coerces from integers to elements of a quadratic number
    field K.
    """
    cdef NumberFieldElement_quadratic zero_element    # the zero element of K

    # TODO: implement __cmp__ properly so we can have a loads/dumps doctest

    def __init__(self, K):
        """
        ``K`` is the target quadratic field

        EXAMPLE::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(ZZ) # indirect doctest
            sage: type(phi)
            <type 'sage.rings.number_field.number_field_element_quadratic.Z_to_quadratic_field_element'>
            sage: phi == loads(dumps(phi)) # todo: comparison not implemented
            True

            sage: R.<b> = CyclotomicField(6)
            sage: psi = R.coerce_map_from(ZZ) # indirect doctest
            sage: type(psi)
            <type 'sage.rings.number_field.number_field_element_quadratic.Z_to_quadratic_field_element'>
            sage: psi == loads(dumps(psi)) # todo: comparison not implemented
            True
        """
        import sage.categories.homset
        Morphism.__init__(self, sage.categories.homset.Hom(ZZ, K))
        self.zero_element = K.zero()

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        TESTS::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(ZZ) # indirect doctest
            sage: phi(4)
            4
            sage: phi(5).parent() is K
            True

        """
        _slots['zero_element'] = self.zero_element
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        TESTS::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(ZZ) # indirect doctest
            sage: phi(4)
            4
            sage: phi(5).parent() is K
            True

        """
        Morphism._update_slots(self, _slots)
        self.zero_element = _slots['zero_element']

    cpdef Element _call_(self, x):
        r"""
        Evaluate at an integer ``x``.

        EXAMPLE::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(ZZ)
            sage: a = phi(2); a # indirect doctest
            2
            sage: a.parent() is K
            True

            sage: R.<b> = CyclotomicField(6)
            sage: psi = R.coerce_map_from(ZZ)
            sage: b = psi(-42); b # indirect doctest
            -42
            sage: b.parent() is R
            True
        """
        cdef NumberFieldElement_quadratic y
        if mpz_sgn((<Integer> x).value) == 0:
            return self.zero_element
        y = self.zero_element._new()
        mpz_set(y.a, (<Integer> x).value)

        # we need to set the denominator to 1 as it is 0 for
        # the zero element of K... (because gcd(0,0) = 0).
        mpz_set_ui(y.denom, 1)

        return y

    def _repr_type(self):
        r"""
        Return a short name for this morphism.

        EXAMPLES::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(ZZ)
            sage: phi # indirect doctest
            Natural morphism:
              From: Integer Ring
              To:   Number Field in a with defining polynomial x^2 - 3

            sage: R.<b> = CyclotomicField(6)
            sage: psi = R.coerce_map_from(ZZ)
            sage: psi # indirect doctest
            Natural morphism:
              From: Integer Ring
              To:   Cyclotomic Field of order 6 and degree 2
        """
        return "Natural"

cdef class Q_to_quadratic_field_element(Morphism):
    """
    Morphism that coerces from rationals to elements of a quadratic number
    field K.
    """
    cdef NumberFieldElement_quadratic zero_element    # the zero element of K

    # TODO: implement __cmp__ properly so we can have a loads/dumps doctest

    def __init__(self, K):
        """
        ``K`` is the target quadratic field

        EXAMPLE::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(QQ) # indirect doctest
            sage: type(phi)
            <type 'sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element'>
            sage: phi == loads(dumps(phi))  # todo: comparison not implemented
            True

            sage: R.<b> = CyclotomicField(6)
            sage: psi = R.coerce_map_from(QQ)
            sage: type(psi)
            <type 'sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element'>
            sage: psi == loads(dumps(psi))  # todo: comparison not implemented
            True
        """
        import sage.categories.homset
        Morphism.__init__(self, sage.categories.homset.Hom(QQ, K))
        self.zero_element = K.zero()

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        TESTS::

            sage: K.<a> = QuadraticField(-3) ## line 8983 ##
            sage: f = K.coerce_map_from(QQ); f # indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3
            sage: f(3/1)
            3
            sage: f(1/2).parent() is K
            True

        """
        _slots['zero_element'] = self.zero_element
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        TESTS::

            sage: K.<a> = QuadraticField(-3) ## line 8983 ##
            sage: f = K.coerce_map_from(QQ); f # indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3
            sage: f(3/1)
            3
            sage: f(1/2).parent() is K
            True

        """
        Morphism._update_slots(self, _slots)
        self.zero_element = _slots['zero_element']

    cpdef Element _call_(self, x):
        r"""
        Evaluate at a rational ``x``.

        EXAMPLE::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(QQ)
            sage: a = phi(2/3); a # indirect doctest
            2/3
            sage: a.parent() is K
            True

            sage: R.<b> = CyclotomicField(6)
            sage: psi = R.coerce_map_from(QQ)
            sage: b = psi(-23/15); b # indirect doctest
            -23/15
            sage: b.parent() is R
            True
        """
        cdef NumberFieldElement_quadratic y = self.zero_element._new()
        mpz_set(y.a, mpq_numref((<Rational>x).value))
        mpz_set(y.denom, mpq_denref((<Rational>x).value))
        return y

    def _repr_type(self):
        r"""
        Return a short name for this morphism.

        EXAMPLES::

            sage: K.<a> = QuadraticField(3)
            sage: phi = K.coerce_map_from(QQ)
            sage: phi # indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 - 3

            sage: R.<b> = CyclotomicField(6)
            sage: psi = R.coerce_map_from(QQ)
            sage: psi # indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Cyclotomic Field of order 6 and degree 2
        """
        return "Natural"
