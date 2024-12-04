r"""
Multivariate Polynomials via libSINGULAR

This module implements specialized and optimized implementations for
multivariate polynomials over many coefficient rings, via a shared library
interface to SINGULAR. In particular, the following coefficient rings are
supported by this implementation:

- the rational numbers `\QQ`,

- the ring of integers `\ZZ`,

- `\ZZ/n\ZZ` for any integer `n`,

- finite fields `\GF{p^n}` for `p` prime and `n > 0`,

- and absolute number fields `\QQ(a)`.

EXAMPLES:

We show how to construct various multivariate polynomial rings::

    sage: P.<x,y,z> = QQ[]
    sage: P
    Multivariate Polynomial Ring in x, y, z over Rational Field

    sage: f = 27/113 * x^2 + y*z + 1/2; f
    27/113*x^2 + y*z + 1/2

    sage: P.term_order()
    Degree reverse lexicographic term order

    sage: P = PolynomialRing(GF(127), 3, names='abc', order='lex'); P
    Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
    sage: a,b,c = P.gens()
    sage: f = 57 * a^2*b + 43 * c + 1; f
    57*a^2*b + 43*c + 1
    sage: P.term_order()
    Lexicographic term order

    sage: z = QQ['z'].0
    sage: K.<s> = NumberField(z^2 - 2)                                                  # needs sage.rings.number_field
    sage: P.<x,y> = PolynomialRing(K, 2)                                                # needs sage.rings.number_field
    sage: 1/2*s*x^2 + 3/4*s                                                             # needs sage.rings.number_field
    (1/2*s)*x^2 + (3/4*s)

    sage: P.<x,y,z> = ZZ[]; P
    Multivariate Polynomial Ring in x, y, z over Integer Ring

    sage: P.<x,y,z> = Zmod(2^10)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

    sage: P.<x,y,z> = Zmod(3^10)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

    sage: P.<x,y,z> = Zmod(2^100)[]; P
    Multivariate Polynomial Ring in x, y, z over
     Ring of integers modulo 1267650600228229401496703205376

    sage: P.<x,y,z> = Zmod(2521352)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352
    sage: type(P)
    <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

    sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
    Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232
    sage: type(P)
    <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>

We construct the Frobenius morphism on `\GF{5}[x,y,z]` over `\GF{5}`::

    sage: R.<x,y,z> = PolynomialRing(GF(5), 3)
    sage: frob = R.hom([x^5, y^5, z^5])
    sage: frob(x^2 + 2*y - z^4)
    -z^20 + x^10 + 2*y^5
    sage: frob((x + 2*y)^3)                                                             # needs sage.rings.finite_rings
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15
    sage: (x^5 + 2*y^5)^3                                                               # needs sage.rings.finite_rings
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15

We make a polynomial ring in one variable over a polynomial ring in
two variables::

    sage: R.<x, y> = PolynomialRing(QQ, 2)
    sage: S.<t> = PowerSeriesRing(R)
    sage: t*(x+y)
    (x + y)*t

TESTS::

    sage: P.<x,y,z> = QQ[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = GF(2^8,'a')[]                                                     # needs sage.rings.finite_rings
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = GF(127)[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True
    sage: P.<x,y,z> = GF(127)[]
    sage: loads(dumps(P)) == P
    True
    sage: loads(dumps(x)) == x
    True

    sage: Rt.<t> = PolynomialRing(QQ, implementation='singular')
    sage: p = 1+t
    sage: R.<u,v> = PolynomialRing(QQ, 2)
    sage: p(u/v)
    (u + v)/v

Check if :issue:`6160` is fixed::

    sage: x = polygen(ZZ, 'x')
    sage: K.<j> = NumberField(x - 1728)                                                 # needs sage.rings.number_field
    sage: R.<b,c> = K[]                                                                 # needs sage.rings.number_field
    sage: b - j*c                                                                       # needs sage.rings.number_field
    b - 1728*c

.. TODO::

    Implement Real, Complex coefficient rings via libSINGULAR

AUTHORS:

- Martin Albrecht (2007-01): initial implementation

- Joel Mohler (2008-01): misc improvements, polishing

- Martin Albrecht (2008-08): added `\QQ(a)` and `\ZZ` support

- Simon King (2009-04): improved coercion

- Martin Albrecht (2009-05): added `\ZZ/n\ZZ` support, refactoring

- Martin Albrecht (2009-06): refactored the code to allow better
  re-use

- Simon King (2011-03): use a faster way of conversion from the base
  ring.

- Volker Braun (2011-06): major cleanup, refcount singular rings, bugfixes.
"""

# ********************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ********************************************************************

# The Singular API is as follows:
#
#   pXXX does assume the currRing to be set
#   p_XXX does not.
#
# However, sometimes there are bugs, i.e. you call p_XXX and it
# crashes unless currRing is set.
#
# Notable exceptions:
#   * pNext and pIter don't need currRing
#   * p_Normalize apparently needs currRing

from cpython.object cimport Py_NE
from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.cpython.string cimport char_to_str, str_to_bytes

# singular types
from sage.libs.singular.decl cimport (ring, poly, ideal, intvec, number,
    currRing, n_unknown, n_Z, n_Zn, n_Znm, n_Z2m, sBucket, sBucketCreate,
    sBucketDestroy, sBucket_Merge_m, sBucketClearMerge, sBucketDeleteAndDestroy)

# singular functions
from sage.libs.singular.decl cimport (
    errorreported,
    n_Invers, n_GetChar,
    p_ISet, rChangeCurrRing, p_Copy, p_Init, p_SetCoeff, p_Setm, p_SetExp, p_Add_q,
    p_NSet, p_GetCoeff, p_Delete, p_GetExp, pNext, rRingVar, omAlloc0, omStrDup,
    omFree, p_Divide, p_SetCoeff0, n_Init, p_DivisibleBy, pLcm, p_LmDivisibleBy,
    pMDivide, p_MDivide, p_IsConstant, p_ExpVectorEqual, p_String, p_LmInit, n_Copy,
    p_IsUnit, p_IsOne, p_Head, idInit, fast_map_common_subexp, id_Delete,
    p_IsHomogeneous, p_Homogen, p_Totaldegree, singclap_pdivide, singclap_factorize,
    idLift, IDELEMS, On, Off, SW_USE_CHINREM_GCD, SW_USE_EZGCD,
    p_LmIsConstant, pTakeOutComp, singclap_gcd, pp_Mult_qq, p_GetMaxExp,
    pLength, kNF, p_Neg, p_Minus_mm_Mult_qq, p_Plus_mm_Mult_qq,
    pDiff, singclap_resultant, p_Normalize,
    prCopyR, prCopyR_NoSort)

# singular conversion routines
from sage.libs.singular.singular cimport si2sa, sa2si, overflow_check

# singular poly arith
from sage.libs.singular.polynomial cimport (
    singular_polynomial_call, singular_polynomial_cmp, singular_polynomial_add,
    singular_polynomial_sub, singular_polynomial_neg, singular_polynomial_rmul,
    singular_polynomial_mul, singular_polynomial_div_coeff, singular_polynomial_pow,
    singular_polynomial_str, singular_polynomial_latex,
    singular_polynomial_str_with_changed_varnames, singular_polynomial_deg,
    singular_polynomial_length_bounded, singular_polynomial_subst )

# singular rings
from sage.libs.singular.ring cimport singular_ring_new, singular_ring_reference, singular_ring_delete

# polynomial imports
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict, MPolynomialRing_polydict_domain
from sage.rings.polynomial.multi_polynomial_ring_base import BooleanPolynomialRing_base
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic

# base ring imports
import sage.rings.abc
import sage.structure.element

from sage.rings.rational cimport Rational
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import IntegerRing_class, ZZ
from sage.rings.integer cimport Integer
from sage.rings.number_field.number_field_base cimport NumberField

from sage.categories.number_fields import NumberFields

from sage.structure.element import coerce_binop

from sage.structure.parent cimport Parent
from sage.structure.category_object cimport CategoryObject

from sage.structure.coerce cimport coercion_model
from sage.structure.element cimport Element

from sage.structure.richcmp cimport rich_to_bool, richcmp
from sage.structure.factorization import Factorization
from sage.structure.sequence import Sequence

from sage.rings.fraction_field import FractionField
from sage.rings.real_mpfr import RealField

import sage.interfaces.abc

from sage.misc.misc_c import prod as mul
from sage.misc.sage_eval import sage_eval

import sage.rings.polynomial.polynomial_singular_interface

cimport cypari2.gen
from sage.rings.polynomial import polynomial_element

permstore=[]
cdef class MPolynomialRing_libsingular(MPolynomialRing_base):

    def __cinit__(self):
        """
        The Cython constructor.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: MPolynomialRing_libsingular(QQ, 3, ('x', 'y', 'z'), TermOrder('degrevlex', 3))
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: type(_)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

            sage: P.<x,y,z> = QQ[]; P
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: type(P)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
        """
        self._ring = NULL

    def __init__(self, base_ring, n, names, order='degrevlex'):
        """
        Construct a multivariate polynomial ring subject to the
        following conditions:

        INPUT:

        - ``base_ring`` -- base ring (must be either GF(q), ZZ, ZZ/nZZ,
                          QQ or absolute number field)

        - ``n`` -- number of variables (must be at least 1)

        - ``names`` -- names of ring variables, may be string of list/tuple

        - ``order`` -- term order (default: ``degrevlex``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]; P
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2
            sage: P.term_order()
            Degree reverse lexicographic term order

            sage: P = PolynomialRing(GF(127), 3, names='abc', order='lex'); P
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1; f
            57*a^2*b + 43*c + 1
            sage: P.term_order()
            Lexicographic term order

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)                                          # needs sage.rings.number_field
            sage: P.<x,y> = PolynomialRing(K, 2)                                        # needs sage.rings.number_field
            sage: 1/2*s*x^2 + 3/4*s                                                     # needs sage.rings.number_field
            (1/2*s)*x^2 + (3/4*s)

            sage: P.<x,y,z> = ZZ[]; P
            Multivariate Polynomial Ring in x, y, z over Integer Ring

            sage: P.<x,y,z> = Zmod(2^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

            sage: P.<x,y,z> = Zmod(3^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

            sage: P.<x,y,z> = Zmod(2^100)[]; P
            Multivariate Polynomial Ring in x, y, z over
             Ring of integers modulo 1267650600228229401496703205376

            sage: P.<x,y,z> = Zmod(2521352)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352
            sage: type(P)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

            sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
            Multivariate Polynomial Ring in x, y, z over
             Ring of integers modulo 25213521351515232
            sage: type(P)
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>

            sage: P.<x,y,z> = PolynomialRing(Integers(2^32), order='lex')
            sage: P(2^32-1)
            4294967295

        TESTS:

        Make sure that a faster coercion map from the base ring is used;
        see :issue:`9944`::

            sage: R.<x,y> = PolynomialRing(ZZ)
            sage: R.coerce_map_from(R.base_ring())
            Polynomial base injection morphism:
              From: Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Integer Ring

        Check some invalid input::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: MPolynomialRing_libsingular(Zmod(1), 1, ["x"], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials over Ring of integers modulo 1 are not supported in Singular
            sage: MPolynomialRing_libsingular(SR, 1, ["x"], "lex")                      # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials over Symbolic Ring are not supported in Singular
            sage: MPolynomialRing_libsingular(QQ, 0, [], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials in 0 variables are not supported in Singular
            sage: MPolynomialRing_libsingular(QQ, -1, [], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials in -1 variables are not supported in Singular
        """
        self._ngens = n
        self._ring = singular_ring_new(base_ring, n, names, order)
        self._zero_element = new_MP(self, NULL)
        cdef MPolynomial_libsingular one = new_MP(self, p_ISet(1, self._ring))
        self._one_element = one
        self._one_element_poly = one._poly
        MPolynomialRing_base.__init__(self, base_ring, n, names, order)
        self._has_singular = True
        #permanently store a reference to this ring until deallocation works reliably
        permstore.append(self)

    def __dealloc__(self):
        r"""
        Deallocate the ring without changing ``currRing``.

        TESTS:

        This example caused a segmentation fault with a previous version
        of this method::

            sage: import gc
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R1 = MPolynomialRing_libsingular(GF(5), 2, ('x', 'y'),
            ....:                                  TermOrder('degrevlex', 2))
            sage: R2 = MPolynomialRing_libsingular(GF(11), 2, ('x', 'y'),
            ....:                                  TermOrder('degrevlex', 2))
            sage: R3 = MPolynomialRing_libsingular(GF(13), 2, ('x', 'y'),
            ....:                                  TermOrder('degrevlex', 2))
            sage: _ = gc.collect()
            sage: foo = R1.gen(0)
            sage: del foo
            sage: del R1
            sage: _ = gc.collect()
            sage: del R2
            sage: _ = gc.collect()
            sage: del R3
            sage: _ = gc.collect()
        """
        if self._ring != NULL:  # the constructor did not raise an exception
            singular_ring_delete(self._ring)

    def __copy__(self):
        """
        Copy ``self``.

        The ring is unique and immutable, so we do not copy.

        TESTS::

            sage: # needs sage.rings.function_field
            sage: import gc
            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: from sage.libs.singular.ring import ring_refcount_dict
            sage: gc.collect()  # random output
            sage: n = len(ring_refcount_dict)
            sage: R = MPolynomialRing_libsingular(GF(547), 2, ('x', 'y'),
            ....:                                 TermOrder('degrevlex', 2))
            sage: len(ring_refcount_dict) == n + 1
            True
            sage: Q = copy(R)   # indirect doctest
            sage: p = R.gen(0)^2 + R.gen(1)^2
            sage: q = copy(p)
            sage: del R
            sage: del Q
            sage: del p
            sage: del q
            sage: gc.collect()  # random output
            sage: len(ring_refcount_dict) == n
            False
        """
        return self

    def __deepcopy__(self, memo):
        """
        Deep copy ``self``.

        The ring should be immutable, so we do not copy.

        TESTS::

            sage: R.<x,y> = GF(547)[]
            sage: R is deepcopy(R)
            True
        """
        memo[id(self)] = self
        return self

    cpdef _coerce_map_from_(self, other):
        """
        Return ``True`` if and only if there exists a coercion map from
        ``other`` to ``self``.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: type(R)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
            sage: R.has_coerce_map_from(ZZ['t'])
            False
            sage: R.coerce_map_from(ZZ['x'])
            Coercion map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Rational Field
        """
        base_ring = self.base_ring()
        if other is base_ring:
            # Because this parent class is a Cython class, the method
            # UnitalAlgebras.ParentMethods.__init_extra__(), which normally
            # registers the coercion map from the base ring, is called only
            # when inheriting from this class in Python (cf. Issue #26958).
            return self._coerce_map_from_base_ring()
        f = self._coerce_map_via([base_ring], other)
        if f is not None:
            return f

        if isinstance(other, MPolynomialRing_libsingular):
            if self is other:
                return True
            n = other.ngens()
            if(other.base_ring is base_ring and self.ngens() >= n and
               self.variable_names()[:n] == other.variable_names()):
                return True
            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
                return True
        elif isinstance(other, MPolynomialRing_polydict):
            if self == other:
                return True
            elif other.ngens() == 0:
                return True
            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
                return True
        elif isinstance(other, PolynomialRing_generic):
            if base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
                return True

    Element = MPolynomial_libsingular

    def _element_constructor_(self, element, check=True):
        """
        Construct a new element in this polynomial ring by converting
        ``element`` into ``self`` if possible.

        INPUT:

        - ``element`` -- several types are supported, see below

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]

        We can coerce elements of ``self`` to ``self``::

            sage: P.coerce(x*y + 1/2)
            x*y + 1/2

        We can coerce elements for a ring with the same algebraic properties::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P == R
            True

            sage: P is R
            False

            sage: P.coerce(x*y + 1)
            x*y + 1

        We can coerce base ring elements::

            sage: P.coerce(3/2)
            3/2

        and all kinds of integers::

            sage: P.coerce(ZZ(1))
            1

            sage: P.coerce(int(1))
            1

            sage: k.<a> = GF(2^8)                                                       # needs sage.rings.finite_rings
            sage: P.<x,y> = PolynomialRing(k,2)                                         # needs sage.rings.finite_rings
            sage: P.coerce(a)                                                           # needs sage.rings.finite_rings
            a

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)                                          # needs sage.rings.number_field
            sage: P.<x,y> = PolynomialRing(K, 2)                                        # needs sage.rings.number_field
            sage: P.coerce(1/2*s)                                                       # needs sage.rings.number_field
            (1/2*s)

        TESTS::

            sage: P.<x,y> = PolynomialRing(GF(127))
            sage: P("111111111111111111111111111111111111111111111111111111111")
            21
            sage: P.<x,y> = PolynomialRing(QQ)
            sage: P("111111111111111111111111111111111111111111111111111111111")
            111111111111111111111111111111111111111111111111111111111
            sage: P("31367566080")
            31367566080

        Check if :issue:`7582` is fixed::

            sage: R.<x,y,z> = PolynomialRing(CyclotomicField(2), 3)                     # needs sage.rings.number_field
            sage: R.coerce(1)                                                           # needs sage.rings.number_field
            1

        Check if :issue:`6160` is fixed::

            sage: x = polygen(ZZ, 'x')
            sage: K.<j> = NumberField(x - 1728)                                         # needs sage.rings.number_field
            sage: R.<b,c> = K[]                                                         # needs sage.rings.number_field
            sage: R.coerce(1)                                                           # needs sage.rings.number_field
            1

        Check if coercion from zero variable polynomial rings work
        (:issue:`7951`)::

            sage: P = PolynomialRing(QQ,0,'')
            sage: R.<x,y> = QQ[]
            sage: P(5)*x
            5*x
            sage: P = PolynomialRing(ZZ,0,'')
            sage: R.<x,y> = GF(127)[]
            sage: R.coerce(P(5))
            5

        Conversion from strings::

            sage: P.<x,y,z> = QQ[]
            sage: P('x+y + 1/4')
            x + y + 1/4

        Coercion from SINGULAR elements::

            sage: P._singular_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_().set_ring()
            sage: P(singular('x + 3/4'))                                                # needs sage.rings.function_field
            x + 3/4

        Coercion from symbolic variables::

            sage: R = QQ['x,y,z']
            sage: var('x')                                                              # needs sage.symbolic
            x
            sage: R(x)                                                                  # needs sage.symbolic
            x

        Coercion from 'similar' rings, which maps by index::

            sage: P.<x,y,z> = QQ[]
            sage: R.<a,b,c> = ZZ[]
            sage: P(a)
            x

        ::

            sage: P.<x,y> = QQ[]
            sage: R.<a,b,c> = QQ[]
            sage: R(x)
            a

        Coercion from PARI objects::

            sage: P.<x,y,z> = QQ[]
            sage: P(pari('x^2 + y'))                                                    # needs sage.libs.pari
            x^2 + y
            sage: P(pari('x*y'))                                                        # needs sage.libs.pari
            x*y

        Coercion from boolean polynomials, also by index::

            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.<x,y,z> = QQ[]
            sage: P(B.gen(0))
            x

        If everything else fails, we try to convert to the base ring::

            sage: R.<x,y,z> = GF(3)[]
            sage: R(1/2)                                                                # needs sage.rings.finite_rings
            -1

        Finally, conversions from other polynomial rings which are not
        coercions are provided. Variables are mapped as follows. Say,
        we are mapping an element from `P` to `Q` (this ring). If the
        variables of `P` are a subset of `Q`, we perform a name
        preserving conversion::

            sage: P.<y_2, y_1, z_3, z_2, z_1> = GF(3)[]
            sage: Q = GF(3)['y_4', 'y_3', 'y_2', 'y_1',
            ....:           'z_5', 'z_4', 'z_3', 'z_2', 'z_1']
            sage: Q(y_1*z_2^2*z_1)
            y_1*z_2^2*z_1

        Otherwise, if `P` has less than or equal the number of
        variables as `Q`, we perform a conversion by index::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a), f(b), f(c)                                                      # needs sage.rings.finite_rings
            (c, b, d)

        ::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d']
            sage: f = Q.convert_map_from(P)
            sage: f(a),f(b),f(c)
            (c, b, d)

        In all other cases, we fail::

            sage: P.<a,b,c,f> = GF(2)[]
            sage: Q = GF(2)['c','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a)                                                                  # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: Could not find a mapping of the passed element to this ring.

        Coerce in a polydict where a coefficient reduces to 0 but isn't 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]
            sage: S((5*x*y + x + 17*y)._mpoly_dict_recursive())
            xx + 2*yy

        Coerce in a polynomial one of whose coefficients reduces to 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S(5*x*y + x + 17*y)
            xx + 2*yy

        Some other examples that illustrate the same coercion idea::

            sage: R.<x,y> = ZZ[]
            sage: S.<xx,yy> = GF(25,'a')[]                                              # needs sage.rings.finite_rings
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

            sage: S.<xx,yy> = Integers(5)[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

        See :issue:`5292`::

            sage: R.<x> = QQ[]; S.<q,t> = QQ[]; F = FractionField(S)
            sage: x in S
            False
            sage: x in F
            False

        Check if :issue:`8228` is fixed::

            sage: P.<x,y> = Zmod(10)[]; P(0)
            0
            sage: P.<x,y> = Zmod(2^10)[]; P(0)
            0

        And :issue:`7597` is fixed if this does not segfault::

            sage: F2 = GF(2)
            sage: F.<x> = GF(2^8)                                                       # needs sage.rings.finite_rings
            sage: R4.<a,b> = PolynomialRing(F)
            sage: R.<u,v> = PolynomialRing(F2)
            sage: P = a
            sage: (P(0,0).polynomial()[0])*u                                            # needs sage.rings.finite_rings
            0
            sage: P(a,b)
            a

        Check that :issue:`15746` is fixed::

            sage: R.<x,y> = GF(7)[]
            sage: R(2^31)
            2

        Check that :issue:`17964` is fixed::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(17)
            sage: Q.<x,y> = K[]
            sage: f = (-3*a)*y + (5*a)
            sage: p = K.primes_above(5)[0]
            sage: R = K.residue_field(p)
            sage: S = R['x','y']
            sage: S(f)
            (2*abar)*y

        Check that creating element from strings works for transcendental extensions::

            sage: T.<c,d> = QQ[]
            sage: F = FractionField(T)
            sage: R.<x,y,z> = F[]
            sage: R('d*z+x^2*y')
            x^2*y + d*z
        """
        cdef poly *_p
        cdef poly *mon
        cdef poly *El_poly
        cdef ring *_ring = self._ring
        cdef ring *El_ring
        cdef MPolynomial_libsingular Element
        cdef MPolynomialRing_libsingular El_parent
        cdef int i, j
        cdef int e
        cdef list ind_map = []
        cdef sBucket *bucket

        if _ring != currRing: rChangeCurrRing(_ring)

        base_ring = self._base

        if isinstance(element, MPolynomial_libsingular):
            Element = <MPolynomial_libsingular>element
            El_parent = Element._parent
            if El_parent is self:
                return element
            El_poly = Element._poly
            El_ring = Element._parent_ring
            El_base = El_parent._base
            El_n = El_parent.ngens()
            if (base_ring is El_base and self.ngens() >= El_n
                    and self.variable_names()[:El_n] == El_parent.variable_names()):
                if self.term_order() == El_parent.term_order():
                    _p = prCopyR_NoSort(El_poly, El_ring, _ring)
                else:
                    _p = prCopyR(El_poly, El_ring, _ring)
                return new_MP(self, _p)
            variable_names_t = self.variable_names()
            if base_ring.has_coerce_map_from(El_parent._mpoly_base_ring(variable_names_t)):
                return self(element._mpoly_dict_recursive(variable_names_t, base_ring))
            else:
                variable_names_s = El_parent.variable_names()
                if set(variable_names_s).issubset(variable_names_t):
                    for v in variable_names_s:
                        ind_map.append(variable_names_t.index(v)+1)
                else:
                    ind_map = [i+1 for i in range(_ring.N)]

                if El_n <= self.ngens():
                    # Map the variables by indices
                    _p = p_ISet(0, _ring)

                    #this loop needs improvement
                    while El_poly:
                        c = si2sa(p_GetCoeff(El_poly, El_ring), El_ring, El_base)
                        if check:
                            try:
                                c = base_ring(c)
                            except TypeError:
                                p_Delete(&_p, _ring)
                                raise
                        if c:
                            mon = p_Init(_ring)
                            p_SetCoeff(mon, sa2si(c, _ring), _ring)
                            for j from 1 <= j <= El_ring.N:
                                e = p_GetExp(El_poly, j, El_ring)
                                if e:
                                    p_SetExp(mon, ind_map[j-1], e, _ring)
                            p_Setm(mon, _ring)
                            _p = p_Add_q(_p, mon, _ring)
                        El_poly = pNext(El_poly)
                    return new_MP(self, _p)

        elif isinstance(element, MPolynomial_polydict):
            if element.parent() == self:
                bucket = sBucketCreate(_ring)
                try:
                    for m, c in element.element().dict().items():
                        mon = p_Init(_ring)
                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
                        for pos in m.nonzero_positions():
                            overflow_check(m[pos], _ring)
                            p_SetExp(mon, pos+1, m[pos], _ring)
                        p_Setm(mon, _ring)
                        #we can use "_m" because we're merging a monomial and
                        #"Merge" because this monomial is different from the rest
                        sBucket_Merge_m(bucket, mon)
                    e=0
                    #we can use "Merge" because the monomials are distinct
                    sBucketClearMerge(bucket, &_p, &e)
                    sBucketDestroy(&bucket)
                except Exception:
                    sBucketDeleteAndDestroy(&bucket)
                    raise
                return new_MP(self, _p)
            elif element.parent().ngens() == 0:
                # zero variable polynomials
                _p = p_NSet(sa2si(base_ring(element[tuple()]), _ring),
                        _ring)
                return new_MP(self, _p)
            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
            else:
                variable_names_s = element.parent().variable_names()
                variable_names_t = self.variable_names()

                if set(variable_names_s).issubset(variable_names_t):
                    for v in variable_names_s:
                        ind_map.append(variable_names_t.index(v)+1)
                else:
                    ind_map = [i+1 for i in range(_ring.N)]

                if element.parent().ngens() <= self.ngens():
                    bucket = sBucketCreate(_ring)
                    try:
                        for m, c in element.element().dict().items():
                            if check:
                                c = base_ring(c)
                            if not c:
                                continue
                            mon = p_Init(_ring)
                            p_SetCoeff(mon, sa2si(c, _ring), _ring)
                            for pos in m.nonzero_positions():
                                overflow_check(m[pos], _ring)
                                p_SetExp(mon, ind_map[pos], m[pos], _ring)
                            p_Setm(mon, _ring)
                            sBucket_Merge_m(bucket, mon)
                        e=0
                        sBucketClearMerge(bucket, &_p, &e)
                        sBucketDestroy(&bucket)
                    except TypeError:
                        sBucketDeleteAndDestroy(&bucket)
                        raise
                    return new_MP(self, _p)

        elif isinstance(element, polynomial_element.Polynomial):
            if base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))

        elif isinstance(element, dict):
            if not element:
                _p = p_ISet(0, _ring)
            else:
                bucket = sBucketCreate(_ring)
                try:
                    for m, c in element.items():
                        if check:
                            c = base_ring(c)
                        if not c:
                            continue
                        mon = p_Init(_ring)
                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
                        if len(m) != self.ngens():
                            raise TypeError("tuple key must have same length as ngens")
                        for pos from 0 <= pos < len(m):
                            if m[pos]:
                                overflow_check(m[pos], _ring)
                                p_SetExp(mon, pos+1, m[pos], _ring)
                        p_Setm(mon, _ring)
                        sBucket_Merge_m(bucket, mon)
                    e=0
                    sBucketClearMerge(bucket, &_p, &e)
                    sBucketDestroy(&bucket)
                except TypeError:
                    sBucketDeleteAndDestroy(&bucket)
                    raise
            return new_MP(self, _p)

        if isinstance(element, sage.structure.element.Element) and isinstance(element.parent(), BooleanPolynomialRing_base):
            if element.constant():
                if element:
                    return self._one_element
                else:
                    return self._zero_element

            variable_names_s = set(element.parent().variable_names())
            variable_names_t = self.variable_names()

            if variable_names_s.issubset(variable_names_t):
                return eval(str(element),self.gens_dict(copy=False))

            elif element.parent().ngens() <= self.ngens():
                Q = element.parent()
                gens_map = dict(zip(Q.variable_names(),self.gens()[:Q.ngens()]))
                return eval(str(element),gens_map)

        if isinstance(element, (sage.interfaces.abc.SingularElement, cypari2.gen.Gen)):
            element = str(element)
        elif isinstance(element, sage.interfaces.abc.Macaulay2Element):
            element = element.external_string()

        if isinstance(element, str):
            # let python do the parsing
            d = self.gens_dict()
            if self.base_ring().gen() != 1:
                if hasattr(self.base_ring(), 'gens'):
                    for gen in self.base_ring().gens():
                        d[str(gen)] = gen
                else:
                    d[str(self.base_ring().gen())] = self.base_ring_gen()
            try:
                if '/' in element:
                    element = sage_eval(element,d)
                else:
                    element = element.replace("^","**")
                    element = eval(element, d, {})
            except (SyntaxError, NameError):
                raise TypeError("Could not find a mapping of the passed element to this ring.")

            # we need to do this, to make sure that we actually get an
            # element in self.
            return self.coerce(element)

        if hasattr(element,'_polynomial_'): # symbolic.expression.Expression
            return element._polynomial_(self)

        try:
            return self(str(element))
        except TypeError:
            pass

        try:
            # now try calling the base ring's __call__ methods
            element = self.base_ring()(element)
            _p = p_NSet(sa2si(element,_ring), _ring)
            return new_MP(self,_p)
        except (TypeError, ValueError):
            raise TypeError("Could not find a mapping of the passed element to this ring.")

    def _repr_(self):
        """
        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: P # indirect doctest
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        varstr = ", ".join(char_to_str(rRingVar(i,self._ring))
                           for i in range(self._ngens))
        return "Multivariate Polynomial Ring in %s over %s" % (varstr, self.base_ring())

    def ngens(self):
        """
        Return the number of variables in this multivariate polynomial ring.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: P.ngens()
            2

            sage: k.<a> = GF(2^16)                                                      # needs sage.rings.finite_rings
            sage: P = PolynomialRing(k, 1000, 'x')                                      # needs sage.rings.finite_rings
            sage: P.ngens()                                                             # needs sage.rings.finite_rings
            1000
        """
        return int(self._ngens)

    def gen(self, int n=0):
        """
        Return the ``n``-th generator of this multivariate polynomial
        ring.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.gen(), P.gen(1)
            (x, y)

            sage: P = PolynomialRing(GF(127), 1000, 'x')
            sage: P.gen(500)
            x500

            sage: P.<SAGE,SINGULAR> = QQ[] # weird names
            sage: P.gen(1)
            SINGULAR
        """
        cdef poly *_p
        cdef ring *_ring = self._ring

        if n < 0 or n >= self._ngens:
            raise ValueError("Generator not defined.")

        rChangeCurrRing(_ring)
        _p = p_ISet(1, _ring)
        p_SetExp(_p, n+1, 1, _ring)
        p_Setm(_p, _ring)

        return new_MP(self, _p)

    def ideal(self, *gens, **kwds):
        """
        Create an ideal in this polynomial ring.

        INPUT:

        - ``*gens`` -- list or tuple of generators (or several input arguments)

        - ``coerce`` -- boolean (default: ``True``); this must be a
          keyword argument. Only set it to ``False`` if you are certain
          that each generator is already in the ring.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: sage.rings.ideal.Katsura(P)                                           # needs sage.rings.function_field
            Ideal (x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y + 2*y*z - y)
             of Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: P.ideal([x + 2*y + 2*z-1, 2*x*y + 2*y*z-y, x^2 + 2*y^2 + 2*z^2-x])
            Ideal (x + 2*y + 2*z - 1, 2*x*y + 2*y*z - y, x^2 + 2*y^2 + 2*z^2 - x)
             of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        coerce = kwds.get('coerce', True)
        if len(gens) == 1:
            gens = gens[0]
        from sage.rings.ideal import Ideal_generic
        if isinstance(gens, Ideal_generic):
            if gens.ring() is self:
                return gens
            gens = gens.gens()
        if isinstance(gens, (sage.interfaces.abc.SingularElement, sage.interfaces.abc.Macaulay2Element)):
            gens = list(gens)
            coerce = True
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return MPolynomialIdeal(self, gens, coerce=False)

    def _macaulay2_(self, macaulay2=None):
        """
        Create an M2 representation of this polynomial ring if
        Macaulay2 is installed.

        INPUT:

        - ``macaulay2`` -- M2 interpreter (default: ``macaulay2_default``)

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: macaulay2(R).describe()  # optional - macaulay2
            ZZ[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 16},
                                                                     {GRevLex => {2:1}  }
                                                                     {Position => Up    }
            --------------------------------------------------------------------------------
            DegreeRank => 1]

            sage: R.<x,y> = QQ[]
            sage: macaulay2(R)  # optional - macaulay2
            QQ[x...y]

            sage: R.<x,y> = GF(17)[]
            sage: macaulay2(R)                  # optional - macaulay2                  # needs sage.rings.finite_rings
            ZZ
            --[x...y]
            17
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2
        try:
            R = self.__macaulay2
            if R is None or not (R.parent() is macaulay2):
                raise ValueError
            R._check_valid()
            return R
        except (AttributeError, ValueError):
            self.__macaulay2 = macaulay2(self._macaulay2_init_(macaulay2))
        return self.__macaulay2

    def _macaulay2_init_(self, macaulay2=None):
        """
        EXAMPLES::

            sage: PolynomialRing(QQ, 'x', 2, order='deglex')._macaulay2_init_()     # optional - macaulay2
            'sage...[symbol x0,symbol x1, MonomialSize=>16, MonomialOrder=>GLex]'
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2
        return macaulay2._macaulay2_input_ring(self.base_ring(), self.gens(),
                                               self.term_order().macaulay2_str())

    def _singular_(self, singular=None):
        """
        Create a SINGULAR (as in the computer algebra system)
        representation of this polynomial ring. The result is cached.

        INPUT:

        - ``singular`` -- SINGULAR interpreter (default: ``sage.interfaces.singular.singular``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P._singular_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_() is P._singular_()
            True

            sage: P._singular_().name() == P._singular_().name()
            True

            sage: k.<a> = GF(3^3)                                                       # needs sage.rings.finite_rings
            sage: P.<x,y,z> = PolynomialRing(k, 3)                                      # needs sage.rings.finite_rings
            sage: P._singular_()                                                        # needs sage.rings.finite_rings
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/3[a]/(a^3-a+1)
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_() is P._singular_()                                      # needs sage.rings.finite_rings
            True

            sage: P._singular_().name() == P._singular_().name()                        # needs sage.rings.finite_rings
            True

        TESTS::

            sage: P.<x> = QQ[]
            sage: P._singular_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
        """
        if singular is None:
            from sage.interfaces.singular import singular
        try:
            R = self.__singular
            if R is None or not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if not self.base_ring().is_field():
                return R
            if self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite() \
                    or (isinstance(self.base_ring(), NumberField) and self.base_ring().is_absolute()):
                R.set_ring() #sorry for that, but needed for minpoly
                if singular.eval('minpoly') != "(" + self.__minpoly + ")":
                    singular.eval("minpoly=%s" % (self.__minpoly))
                    self.__minpoly = singular.eval('minpoly')[1:-1] # store in correct format
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=None):
        """
        Create a SINGULAR (as in the computer algebra system)
        representation of this polynomial ring. The result is NOT
        cached.

        INPUT:

        - ``singular`` -- SINGULAR interpreter (default: ``sage.interfaces.singular.singular``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P._singular_init_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C
            sage: P._singular_init_() is P._singular_init_()
            False

            sage: P._singular_init_().name() == P._singular_init_().name()
            False

            sage: w = polygen(ZZ, 'w')
            sage: R.<x,y> = PolynomialRing(NumberField(w^2 + 1,'s'))                    # needs sage.rings.number_field
            sage: singular(R)                                                           # needs sage.rings.number_field
            polynomial ring, over a field, global ordering
            //   coefficients: QQ[s]/(s^2+1)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(2**8,'a'),10,'x', order='invlex')               # needs sage.rings.finite_rings
            sage: singular(R)                                                           # needs sage.rings.finite_rings
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/2[a]/(a^8+a^4+a^3+a^2+1)
            //   number of vars : 10
            //        block   1 : ordering ip
            //                  : names    x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127),2,'x', order='invlex')
            sage: singular(R)                                                           # needs sage.rings.finite_rings
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/127
            //   number of vars : 2
            //        block   1 : ordering ip
            //                  : names    x0 x1
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,2,'x', order='invlex')
            sage: singular(R)                                                           # needs sage.rings.function_field
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 2
            //        block   1 : ordering ip
            //                  : names    x0 x1
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,2,'x', order='degneglex')
            sage: singular(R)                                                           # needs sage.rings.function_field
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 2
            //        block   1 : ordering a
            //                  : names    x0 x1
            //                  : weights   1  1
            //        block   2 : ordering ls
            //                  : names    x0 x1
            //        block   3 : ordering C

            sage: R = PolynomialRing(QQ,'x')
            sage: singular(R)                                                           # needs sage.rings.function_field
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127),'x')
            sage: singular(R)                                                           # needs sage.rings.finite_rings
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = ZZ['x,y']
            sage: singular(R)                                                           # needs sage.rings.function_field
            polynomial ring, over a domain, global ordering
            //   coefficients: ZZ
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = IntegerModRing(1024)['x,y']
            sage: singular(R)                                                           # needs sage.rings.function_field
            polynomial ring, over a ring (with zero-divisors), global ordering
            //   coefficients: ZZ/(2^10)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = IntegerModRing(15)['x,y']
            sage: singular(R)                                                           # needs sage.rings.function_field
            polynomial ring, over a ring (with zero-divisors), global ordering
            //   coefficients: ZZ/...(15)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

        TESTS::

            sage: P.<x> = QQ[]
            sage: P._singular_init_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
        """
        if singular is None:
            from sage.interfaces.singular import singular

        if self.ngens() == 1:
            _vars = str(self.gen())
            if "*" in _vars: # 1.000...000*x
                _vars = _vars.split("*")[1]
        else:
            _vars = str(self.gens())

        order = self.term_order().singular_str()%dict(ngens=self.ngens())

        self.__singular, self.__minpoly = \
                sage.rings.polynomial.polynomial_singular_interface._do_singular_init_(
                        singular, self.base_ring(), self.characteristic(), _vars, order)

        return self.__singular

    # It is required in cython to redefine __hash__ when __richcmp__ is
    # overloaded. Also just writing
    #         __hash__ = CategoryObject.__hash__
    # doesn't work.
    def __hash__(self):
        """
        Return a hash for this ring, that is, a hash of the string
        representation of this polynomial ring.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: hash(P)      # somewhat random output
            967902441410893180 # 64-bit
            -1767675994        # 32-bit
        """
        return CategoryObject.__hash__(self)

    def __richcmp__(left, right, int op):
        r"""
        Multivariate polynomial rings are said to be equal if:

        - their base rings match,
        - their generator names match and
        - their term orderings match.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: R.<x,y,z> = QQ[]
            sage: P == R
            True

            sage: R.<x,y,z> = GF(127)[]
            sage: P == R
            False

            sage: R.<x,y> = QQ[]
            sage: P == R
            False

            sage: R.<x,y,z> = PolynomialRing(QQ, order='invlex')
            sage: P == R
            False

        TESTS::

            sage: R = QQ['x', 'y']; R
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: R == R
            True
            sage: R == QQ['x','z']
            False
        """
        if left is right:
            return rich_to_bool(op, 0)

        if not isinstance(right, Parent) or not isinstance(left, Parent):
            # One is not a parent -- not equal and not ordered
            return op == Py_NE

        if not isinstance(right, (MPolynomialRing_libsingular,
                                  MPolynomialRing_polydict_domain)):
            return op == Py_NE

        def cmp_key(x):
            return (x.base_ring(), [str(v) for v in x.gens()], x.term_order())
        return richcmp(cmp_key(left), cmp_key(right), op)

    def __reduce__(self):
        """
        Serialize ``self``.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P.<x,y,z> = PolynomialRing(ZZ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(127), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^8,'F'), names='abc')                          # needs sage.rings.finite_rings
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^16,'B'), names='abc')                         # needs sage.rings.finite_rings
            sage: P == loads(dumps(P))
            True
            sage: z = QQ['z'].0
            sage: P = PolynomialRing(NumberField(z^2 + 3, 'B'), names='abc')            # needs sage.rings.number_field
            sage: P == loads(dumps(P))                                                  # needs sage.rings.number_field
            True
        """
        return unpickle_MPolynomialRing_libsingular, \
            (self.base_ring(), self.variable_names(), self.term_order())

    def _temporarily_change_names(self, names, latex_names):
        """
        This is used by the variable names context manager.

        EXAMPLES::

            sage: R.<x,y> = QQ[] # indirect doctest
            sage: with localvars(R, 'z,w'):
            ....:      print(x^3 + y^3 - x*y)
            z^3 + w^3 - z*w
        """
        cdef ring *_ring = self._ring
        cdef char **_names
        cdef char **_orig_names
        cdef int i

        if len(names) != _ring.N:
            raise TypeError("len(names) doesn't equal self.ngens()")

        old = self._names, self._latex_names
        (self._names, self._latex_names) = names, latex_names

        _names = <char**>omAlloc0(sizeof(char*)*_ring.N)
        for i from 0 <= i < _ring.N:
            _name = str_to_bytes(names[i])
            _names[i] = omStrDup(_name)

        _orig_names = _ring.names
        _ring.names = _names

        for i from 0 <= i < _ring.N:
            omFree(_orig_names[i])
        omFree(_orig_names)

        return old

    # The following methods are handy for implementing Groebner
    # basis algorithms. They do only superficial type/sanity checks
    # and should be called carefully.

    def monomial_quotient(self, MPolynomial_libsingular f, MPolynomial_libsingular g, coeff=False):
        r"""
        Return ``f/g``, where both ``f`` and`` ``g`` are treated as
        monomials.

        Coefficients are ignored by default.

        INPUT:

        - ``f`` -- monomial
        - ``g`` -- monomial
        - ``coeff`` -- divide coefficients as well (default: ``False``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_quotient(3/2*x*y,x)
            y

            sage: P.monomial_quotient(3/2*x*y,x,coeff=True)
            3/2*y

        Note, that `\ZZ` behaves different if ``coeff=True``::

            sage: P.monomial_quotient(2*x,3*x)
            1

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: P.monomial_quotient(2*x,3*x,coeff=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot divide these coefficients.

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_quotient(x*y,x)
            y

            sage: P.monomial_quotient(x*y,R.gen())
            y

            sage: P.monomial_quotient(P(0),P(1))
            0

            sage: P.monomial_quotient(P(1),P(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

            sage: P.monomial_quotient(P(3/2),P(2/3), coeff=True)
            9/4

            sage: P.monomial_quotient(x,y) # Note the wrong result
            x*y^65535*z^65535      # 32-bit
            x*y^1048575*z^1048575  # 64-bit

            sage: P.monomial_quotient(x,P(1))
            x

        .. warning::

           Assumes that the head term of f is a multiple of the head
           term of g and return the multiplicant m. If this rule is
           violated, funny things may happen.
        """
        cdef poly *res
        cdef ring *r = self._ring
        cdef number *n

        if self is not f._parent:
            f = self.coerce(f)
        if self is not g._parent:
            g = self.coerce(g)

        if not f._poly:
            return self._zero_element
        if not g._poly:
            raise ZeroDivisionError

        if r is not currRing:
            rChangeCurrRing(r)
        res = pMDivide(f._poly, g._poly)
        if coeff:
            if r.cf.type == n_unknown or r.cf.cfDivBy(p_GetCoeff(f._poly, r), p_GetCoeff(g._poly, r), r.cf):
                n = r.cf.cfDiv(p_GetCoeff(f._poly, r), p_GetCoeff(g._poly, r), r.cf)
                p_SetCoeff0(res, n, r)
            else:
                raise ArithmeticError("Cannot divide these coefficients.")
        else:
            p_SetCoeff0(res, n_Init(1, r.cf), r)
        return new_MP(self, res)

    def monomial_divides(self, MPolynomial_libsingular a, MPolynomial_libsingular b):
        """
        Return ``False`` if a does not divide b and ``True``
        otherwise.

        Coefficients are ignored.

        INPUT:

        - ``a`` -- monomial

        - ``b`` -- monomial

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_divides(x*y*z, x^3*y^2*z^4)
            True
            sage: P.monomial_divides(x^3*y^2*z^4, x*y*z)
            False

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_divides(P(1), P(0))
            True
            sage: P.monomial_divides(P(1), x)
            True
        """
        cdef poly *_a
        cdef poly *_b
        cdef ring *_r
        if a._parent is not b._parent:
            b = a._parent.coerce(b)

        _a = a._poly
        _b = b._poly
        _r = a._parent_ring

        if _a == NULL:
            raise ZeroDivisionError
        if _b == NULL:
            return True

        if not p_DivisibleBy(_a, _b, _r):
            return False
        else:
            return True

    def monomial_lcm(self, MPolynomial_libsingular f, MPolynomial_libsingular g):
        """
        LCM for monomials. Coefficients are ignored.

        INPUT:

        - ``f`` -- monomial

        - ``g`` -- monomial

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_lcm(3/2*x*y,x)
            x*y

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_lcm(x*y,R.gen())
            x*y

            sage: P.monomial_lcm(P(3/2),P(2/3))
            1

            sage: P.monomial_lcm(x,P(1))
            x
        """
        cdef poly *m = p_ISet(1, self._ring)

        if self is not f._parent:
            f = self.coerce(f)
        if self is not g._parent:
            g = self.coerce(g)

        if f._poly == NULL:
            if g._poly == NULL:
                return self._zero_element
            else:
                raise ArithmeticError("Cannot compute LCM of zero and nonzero element.")
        if g._poly == NULL:
            raise ArithmeticError("Cannot compute LCM of zero and nonzero element.")

        if self._ring != currRing:
            rChangeCurrRing(self._ring)

        pLcm(f._poly, g._poly, m)
        p_Setm(m, self._ring)
        return new_MP(self,m)

    def monomial_reduce(self, MPolynomial_libsingular f, G):
        """
        Try to find a ``g`` in ``G`` where ``g.lm()`` divides
        ``f``. If found ``(flt,g)`` is returned, ``(0,0)`` otherwise,
        where ``flt`` is ``f/g.lm()``.

        It is assumed that ``G`` is iterable and contains *only*
        elements in this polynomial ring.

        Coefficients are ignored.

        INPUT:

        - ``f`` -- monomial
        - ``G`` -- list/set of mpolynomials

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, 1/2  ]
            sage: P.monomial_reduce(f,G)
            (y, 1/4*x*y + 2/7)

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, 1/2  ]

            sage: P.monomial_reduce(P(0),G)
            (0, 0)

            sage: P.monomial_reduce(f,[P(0)])
            (0, 0)
        """
        cdef poly *m = f._poly
        cdef ring *r = self._ring
        cdef poly *flt

        if not m:
            return (f, f)

        for g in G:
            if isinstance(g, MPolynomial_libsingular) and g:
                h = <MPolynomial_libsingular>g
                if h._parent is self and p_LmDivisibleBy(h._poly, m, r):
                    if r is not currRing:
                        rChangeCurrRing(r)
                    flt = pMDivide(f._poly, h._poly)
                    p_SetCoeff(flt, n_Init(1, r.cf), r)
                    return (new_MP(self, flt), h)
        return (self._zero_element, self._zero_element)

    def monomial_pairwise_prime(self, MPolynomial_libsingular g, MPolynomial_libsingular h):
        """
        Return ``True`` if ``h`` and ``g`` are pairwise prime. Both
        are treated as monomials.

        Coefficients are ignored.

        INPUT:

        - ``h`` -- monomial
        - ``g`` -- monomial

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_pairwise_prime(x^2*z^3, y^4)
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, 3/4*y^3)
            False

        TESTS::

            sage: Q.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_pairwise_prime(x^2*z^3, Q('y^4'))
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, Q(0))
            True

            sage: P.monomial_pairwise_prime(P(1/2),x)
            False
        """
        cdef int i
        cdef ring *r
        cdef poly *p
        cdef poly *q

        if h._parent is not g._parent:
            g = h._parent.coerce(g)

        r = h._parent_ring
        p = g._poly
        q = h._poly

        if p == NULL:
            if q == NULL:
                return False #GCD(0,0) = 0
            else:
                return True #GCD(x,0) = 1
        elif q == NULL:
            return True # GCD(0,x) = 1
        elif p_IsConstant(p,r) or p_IsConstant(q,r): # assuming a base field
            return False

        for i from 1 <= i <= r.N:
            if p_GetExp(p,i,r) and p_GetExp(q,i,r):
                return False
        return True

    def monomial_all_divisors(self, MPolynomial_libsingular t):
        """
        Return a list of all monomials that divide ``t``.

        Coefficients are ignored.

        INPUT:

        - ``t`` -- a monomial

        OUTPUT: list of monomials

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.monomial_all_divisors(x^2*z^3)
            [x, x^2, z, x*z, x^2*z, z^2, x*z^2, x^2*z^2, z^3, x*z^3, x^2*z^3]

        ALGORITHM: addwithcarry idea by Toon Segers
        """

        M = list()

        cdef ring *_ring = self._ring
        cdef poly *maxvector = t._poly
        cdef poly *tempvector = p_ISet(1, _ring)

        pos = 1

        while not p_ExpVectorEqual(tempvector, maxvector, _ring):
            tempvector = addwithcarry(tempvector, maxvector, pos, _ring)
            M.append(new_MP(self, p_Copy(tempvector, _ring)))
        return M


def unpickle_MPolynomialRing_libsingular(base_ring, names, term_order):
    """
    Inverse function for ``MPolynomialRing_libsingular.__reduce__``.

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(P)) is P # indirect doctest
        True
    """
    from sage.rings.polynomial.polynomial_ring_constructor import _multi_variate
    # If libsingular would be replaced by a different implementation in future
    # sage version, the unpickled ring will belong the new implementation.
    return _multi_variate(base_ring, tuple(names), None, term_order, None)


cdef class MPolynomial_libsingular(MPolynomial_libsingular_base):
    """
    A multivariate polynomial implemented using libSINGULAR.
    """
    def __init__(self, MPolynomialRing_libsingular parent):
        """
        Construct a zero element in parent.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomial_libsingular
            sage: P = PolynomialRing(GF(32003), 3, 'x')                                 # needs sage.rings.finite_rings
            sage: MPolynomial_libsingular(P)                                            # needs sage.rings.finite_rings
            0
        """
        self._poly = NULL
        self._parent = parent
        self._parent_ring = singular_ring_reference(parent._ring)

    def __dealloc__(self):
        # WARNING: the Cython class self._parent is now no longer accessible!
        if self._poly==NULL:
            # e.g. MPolynomialRing_libsingular._zero_element
            singular_ring_delete(self._parent_ring)
            return
        assert self._parent_ring != NULL # the constructor has no way to raise an exception
        p_Delete(&self._poly, self._parent_ring)
        singular_ring_delete(self._parent_ring)

    def __copy__(self):
        """
        Copy ``self``.

        OUTPUT: a copy

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(7^2)
            sage: R.<x,y> = F[]
            sage: p = a*x^2 + y + a^3; p
            a*x^2 + y + (-2*a - 3)
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
            sage: lst = [p,q]
            sage: matrix(ZZ, 2, 2, lambda i,j: bool(lst[i]==lst[j]))
            [1 1]
            [1 1]
            sage: matrix(ZZ, 2, 2, lambda i,j: bool(lst[i] is lst[j]))
            [1 0]
            [0 1]
        """
        return new_MP(self._parent, p_Copy(self._poly, self._parent_ring))

    def __deepcopy__(self, memo={}):
        """
        Deep copy ``self``.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: p = x^2 + y^2
            sage: p is deepcopy(p)
            False
            sage: p == deepcopy(p)
            True
            sage: p.parent() is deepcopy(p).parent()
            True
        """
        cpy = self.__copy__()
        memo[id(self)] = cpy
        return cpy

    cpdef MPolynomial_libsingular _new_constant_poly(self, x, MPolynomialRing_libsingular P):
        r"""
        Quickly create a new constant polynomial with value x in the parent P.

        ASSUMPTION:

        The value x must be an element of the base ring. That assumption is
        not verified.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: x._new_constant_poly(2/1,R)
            2
        """
        if not x:
            return new_MP(P, NULL)
        cdef poly *_p
        singular_polynomial_rmul(&_p, P._one_element_poly, x, P._ring)
        return new_MP(P, _p)

    def __call__(self, *x, **kwds):
        r"""
        Evaluate this multi-variate polynomial at ``x``, where ``x``
        is either the tuple of values to substitute in, or one can use
        functional notation ``f(a_0,a_1,a_2, \ldots)`` to evaluate
        ``f`` with the ``i``-th variable replaced by ``a_i``.

        INPUT:

        - ``x`` -- list of elements in ``self.parent()``
        - or ``**kwds`` -- dictionary of ``variable-name:value`` pairs

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = 3/2*x^2*y + 1/7 * y^2 + 13/27
            sage: f(0,0,0)
            13/27

            sage: f(1,1,1)
            803/378
            sage: 3/2 + 1/7 + 13/27
            803/378

            sage: f(45/2,19/3,1)
            7281167/1512

            sage: f(1,2,3).parent()
            Rational Field

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: P(0)(1,2,3)
            0
            sage: P(3/2)(1,2,3)
            3/2

            sage: R.<a,b,y> = QQ[]
            sage: f = a*y^3 + b*y - 3*a*b*y
            sage: f(a=5,b=3,y=10)
            4580
            sage: f(5,3,10)
            4580

        See :issue:`8502`::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^2 + 47)
            sage: R.<X,Y,Z> = K[]
            sage: f = X + Y + Z
            sage: a = f(t, t, t); a
            3*t
            sage: a.parent() is K
            True

            sage: R.<X,Y,Z> = QQ[]
            sage: f = X+Y+Z
            sage: a = f(2,3,4); a
            9
            sage: a.parent() is QQ
            True

        See :issue:`33373`::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(2^4)
            sage: R.<x> = PolynomialRing(k, 1)
            sage: f = R(1)
            sage: S.<y> = PolynomialRing(k, 1)
            sage: f(y).parent()
            Multivariate Polynomial Ring in y over Finite Field in a of size 2^4
        """
        cdef Element sage_res

        if len(kwds) > 0:
            f = self.subs(**kwds)
            if len(x) > 0:
                return f(*x)
            else:
                return f

        cdef int l = len(x)
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring

        if l == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
            l = len(x)

        if l != parent._ring.N:
            raise TypeError("number of arguments does not match number of variables in parent")

        res_parent = coercion_model.common_parent(parent._base, *x)
        cdef poly *res    # ownership will be transferred to us in the else block
        try:
            coerced_x = [parent.coerce(e) for e in x]
        except TypeError:
            # give up, evaluate functional
            sage_res = parent.base_ring().zero()
            for m, c in self.monomial_coefficients().items():
                sage_res += c * mul([x[i] ** m[i] for i in m.nonzero_positions()])
        else:
            singular_polynomial_call(&res, self._poly, _ring, coerced_x,
                                     MPolynomial_libsingular_get_element)

            if res == NULL:
                return res_parent(0)
            if p_LmIsConstant(res, _ring):
                sage_res = si2sa(p_GetCoeff(res, _ring), _ring, parent._base)
                p_Delete(&res, _ring)            # sage_res contains copy
            else:
                sage_res = new_MP(parent, res)   # pass on ownership of res to sage_res

        if sage_res._parent is not res_parent:
            sage_res = res_parent(sage_res)
        return sage_res

    def __hash__(self):
        """
        This hash incorporates the variable name in an effort to
        respect the obvious inclusions into multi-variable polynomial
        rings.

        The tuple algorithm is borrowed from http://effbot.org/zone/python-hash.htm.

        EXAMPLES::

            sage: R.<x>=QQ[]
            sage: S.<x,y>=QQ[]
            sage: hash(S(1/2))==hash(1/2)  # respect inclusions of the rationals
            True
            sage: hash(S.0)==hash(R.0)  # respect inclusions into mpoly rings
            True
            sage: # the point is to make for more flexible dictionary look ups
            sage: d={S.0:12}
            sage: d[R.0]
            12
        """
        return self._hash_c()

    cpdef _richcmp_(left, right, int op):
        """
        Compare left and right.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,order='degrevlex')
            sage: x == x
            True

            sage: x > y
            True
            sage: y^2 > x
            True

            sage: (2/3*x^2 + 1/2*y + 3) > (2/3*x^2 + 1/4*y + 10)
            True

        TESTS::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: x > P(0)
            True

            sage: P(0) == P(0)
            True

            sage: P(0) < P(1)
            True

            sage: x > P(1)
            True

            sage: 1/2*x < 3/4*x
            True

            sage: (x+1) > x
            True

            sage: f = 3/4*x^2*y + 1/2*x + 2/7
            sage: f > f
            False
            sage: f < f
            False
            sage: f == f
            True

            sage: P.<x,y,z> = PolynomialRing(GF(127), order='degrevlex')
            sage: (66*x^2 + 23) > (66*x^2 + 2)
            True
        """
        if left is right:
            return rich_to_bool(op, 0)
        cdef poly *p = (<MPolynomial_libsingular>left)._poly
        cdef poly *q = (<MPolynomial_libsingular>right)._poly
        cdef ring *r = (<MPolynomial_libsingular>left)._parent_ring
        return rich_to_bool(op, singular_polynomial_cmp(p, q, r))

    cpdef _add_(left, right):
        """
        Add left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x + 1/2*y + 1 #indirect doctest
            3/2*x + 1/2*y + 1
        """
        cdef ring *r = (<MPolynomial_libsingular>left)._parent_ring
        cdef poly *_p
        singular_polynomial_add(&_p, left._poly,
                                 (<MPolynomial_libsingular>right)._poly, r)
        return new_MP((<MPolynomial_libsingular>left)._parent, _p)

    cpdef _sub_(left, right):
        """
        Subtract left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x - 1/2*y - 1 #indirect doctest
            3/2*x - 1/2*y - 1
        """
        cdef ring *_ring = (<MPolynomial_libsingular>left)._parent_ring
        cdef poly *_p
        singular_polynomial_sub(&_p, left._poly,
                                (<MPolynomial_libsingular>right)._poly,
                                _ring)
        return new_MP((<MPolynomial_libsingular>left)._parent, _p)

    cpdef _lmul_(self, Element left):
        """
        Multiply ``self`` with a base ring element.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x # indirect doctest
            3/2*x

        ::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2) # indirect doctest
            9/4*x - 3/4*y - 3/2
        """

        cdef ring *_ring = (<MPolynomial_libsingular>self)._parent_ring
        if not left:
            return (<MPolynomial_libsingular>self)._parent._zero_element
        cdef poly *_p
        singular_polynomial_rmul(&_p, self._poly, left, _ring)
        return new_MP((<MPolynomial_libsingular>self)._parent, _p)

    cpdef _mul_(left, right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1) # indirect doctest
            9/4*x^2 - 1/4*y^2 - y - 1

            sage: P.<x,y> = PolynomialRing(QQ,order='lex')
            sage: (x^2^32) * x^2^32
            Traceback (most recent call last):
            ...
            OverflowError: Python int too large to convert to C unsigned long  # 32-bit
            OverflowError: exponent overflow (...)  # 64-bit
        """
        # all currently implemented rings are commutative
        cdef poly *_p
        singular_polynomial_mul(&_p, left._poly,
                                 (<MPolynomial_libsingular>right)._poly,
                                 (<MPolynomial_libsingular>left)._parent_ring)
        return new_MP((<MPolynomial_libsingular>left)._parent,_p)

    cpdef _div_(left, right_ringelement):
        r"""
        Divide ``left`` by ``right``.

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = (x + y)/3 # indirect doctest
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent and right is not constant. ::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: g = x
            sage: h = f/g; h
            (x^3 + y)/x
            sage: h.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field

        If we divide over `\ZZ` the result is the same as multiplying
        by 1/3 (i.e. base extension). ::

            sage: R.<x,y> = ZZ[]
            sage: f = (x + y)/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: f = (x + y) * 1/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        But we get a true fraction field if the denominator is not in
        the fraction field of the base ring.""

            sage: f = x/y
            sage: f.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring

        Division will fail for non-integral domains::

            sage: P.<x,y> = Zmod(1024)[]
            sage: x/P(3)
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

            sage: x/3
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

        TESTS::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        Ensure that :issue:`17638` is fixed::

            sage: R.<x,y> = PolynomialRing(QQ, order='neglex')
            sage: f = 1 + y
            sage: g = 1 + x
            sage: h = f/g
            sage: h*g == f
            True
        """
        cdef poly *p
        cdef MPolynomial_libsingular right = <MPolynomial_libsingular>right_ringelement
        cdef bint is_field = left._parent._base.is_field()
        if p_IsConstant(right._poly, right._parent_ring):
            if is_field:
                singular_polynomial_div_coeff(&p, left._poly, right._poly, right._parent_ring)
                return new_MP(left._parent, p)
            else:
                return left.change_ring(left.base_ring().fraction_field())/right_ringelement
        else:
            return (left._parent).fraction_field()(left,right_ringelement)

    def __pow__(MPolynomial_libsingular self, exp, mod):
        """
        Return ``self**(exp)``.

        The exponent must be an integer or a rational such that
        the result lies in the same polynomial ring.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f^2
            x^6 + 2*x^3*y + y^2
            sage: g = f^(-1); g
            1/(x^3 + y)
            sage: type(g)
            <class 'sage.rings.fraction_field_element.FractionFieldElement'>

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: P(2)**(-1)
            1/2

            sage: P.<u,v> = PolynomialRing(QQ, 2)
            sage: u^(1/2)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power

            sage: P.<x,y> = PolynomialRing(QQ,order='lex')
            sage: (x+y^2^32)^10
            Traceback (most recent call last):
            ....
            OverflowError: Python int too large to convert to C unsigned long  # 32-bit
            OverflowError: exponent overflow (...)  # 64-bit

        Test fractional powers (:issue:`22329`)::

            sage: P.<R, S> = ZZ[]
            sage: (R^3 + 6*R^2*S + 12*R*S^2 + 8*S^3)^(1/3)
            R + 2*S
            sage: _.parent()
            Multivariate Polynomial Ring in R, S over Integer Ring
            sage: P(4)^(1/2)
            2
            sage: _.parent()
            Multivariate Polynomial Ring in R, S over Integer Ring

            sage: (R^2 + 3)^(1/2)
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a 2nd power
            sage: P(2)^P(2)
            4
            sage: (R + 1)^P(2)
            R^2 + 2*R + 1
            sage: (R + 1)^R
            Traceback (most recent call last):
            ...
            TypeError: R is neither an integer nor a rational
            sage: 2^R
            Traceback (most recent call last):
            ...
            TypeError: R is neither an integer nor a rational

        Check that using third argument raises an error::

            sage: R.<x,y,z> = PolynomialRing(ZZ)
            sage: pow(x + y + z, 2, x)
            Traceback (most recent call last):
            ...
            NotImplementedError: pow() with a modulus is not implemented for this ring
        """
        if mod is not None:
            raise NotImplementedError(
                "pow() with a modulus is not implemented for this ring"
            )

        if type(exp) is not Integer:
            try:
                exp = Integer(exp)
            except TypeError:
                try:
                    n = Rational(exp)
                except TypeError:
                    raise TypeError("{} is neither an integer nor a rational".format(exp))
                num = n.numerator()
                den = n.denominator()

                if self.degree == 0:
                    return self.parent()(
                        self.constant_coefficient().nth_root(den) ** num)
                return self.nth_root(den) ** num

        if exp < 0:
            return self._parent._one_element / (self**(-exp))
        elif exp == 0:
            return self._parent._one_element

        cdef ring *_ring = self._parent_ring
        cdef poly *_p
        singular_polynomial_pow(&_p, self._poly, exp, _ring)
        return new_MP(self._parent, _p)

    def __neg__(self):
        """
        Return ``-self``.

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: -f
            -x^3 - y
        """
        cdef ring *_ring = self._parent_ring

        cdef poly *p
        singular_polynomial_neg(&p, self._poly, _ring)
        return new_MP(self._parent, p)

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f # indirect doctest
            x^3 + y
        """
        cdef ring *_ring = self._parent_ring
        s = singular_polynomial_str(self._poly, _ring)
        return s

    cpdef _repr_short_(self):
        """
        This is a faster but less pretty way to print polynomials. If
        available it uses the short SINGULAR notation.

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f._repr_short_()
            'x3+y'
        """
        cdef ring *_ring = self._parent_ring
        rChangeCurrRing(_ring)
        if _ring.CanShortOut:
            _ring.ShortOut = 1
            s = p_String(self._poly, _ring, _ring)
            _ring.ShortOut = 0
        else:
            s = p_String(self._poly, _ring, _ring)
        return char_to_str(s)

    def _latex_(self):
        """
        Return a polynomial LaTeX representation of this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: latex(f)  # indirect doctest
            -x^{2} y - \frac{25}{27} y^{3} - z^{2}
        """
        cdef ring *_ring = self._parent_ring
        gens = self.parent().latex_variable_names()
        base = self.parent().base()
        return singular_polynomial_latex(self._poly, _ring, base, gens)

    def _repr_with_changed_varnames(self, varnames):
        """
        Return string representing this polynomial but change the
        variable names to ``varnames``.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = - 1*x^2*y - 25/27 * y^3 - z^2
            sage: print(f._repr_with_changed_varnames(['FOO', 'BAR', 'FOOBAR']))
            -FOO^2*BAR - 25/27*BAR^3 - FOOBAR^2
        """
        return singular_polynomial_str_with_changed_varnames(self._poly,
                                                             self._parent_ring,
                                                             varnames)

    def degree(self, MPolynomial_libsingular x=None, int std_grading=False):
        """
        Return the degree of this polynomial.

        INPUT:

        - ``x`` -- (default: ``None``) a generator of the parent ring

        OUTPUT:

        If ``x`` is ``None``, return the total degree of ``self``. Note that
        this result is affected by the weighting given to the generators of the
        parent ring. Otherwise, if ``x`` is (or is coercible to) a generator of
        the parent ring, the output is the maximum degree of ``x`` in ``self``.
        This is not affected by the weighting of the generators.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = y^2 - x^9 - x
            sage: f.degree(x)
            9
            sage: f.degree(y)
            2
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(x)
            3
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(y)
            10

        When the generators have a grading (weighting) then the total degree
        respects this, but the degree for a given generator is unaffected::

            sage: T = TermOrder("wdegrevlex", (2, 3))
            sage: R.<x, y> = PolynomialRing(QQ, order=T)
            sage: f = x^2 * y + y^4
            sage: f.degree()
            12
            sage: f.degree(x)
            2
            sage: f.degree(y)
            4

        The term ordering of the parent ring determines the grading of the
        generators. ::

            sage: T = TermOrder('wdegrevlex', (1,2,3,4))
            sage: R = PolynomialRing(QQ, 'x', 12, order=T+T+T)
            sage: [x.degree() for x in R.gens()]
            [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]

        A matrix term ordering determines the grading of the generators by the
        first row of the matrix. ::

            sage: m = matrix(3, [3,2,1,1,1,0,1,0,0])
            sage: m
            [3 2 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=TermOrder(m))
            sage: x.degree(), y.degree(), z.degree()
            (3, 2, 1)
            sage: f = x^3*y + x*z^4
            sage: f.degree()
            11

        If the first row contains zero, the grading becomes the standard one. ::

            sage: m = matrix(3, [3,0,1,1,1,0,1,0,0])
            sage: m
            [3 0 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=TermOrder(m))
            sage: x.degree(), y.degree(), z.degree()
            (1, 1, 1)
            sage: f = x^3*y + x*z^4
            sage: f.degree()
            5

        To get the degree with the standard grading regardless of the term
        ordering of the parent ring, use ``std_grading=True``. ::

            sage: f.degree(std_grading=True)
            5

        TESTS::

            sage: P.<x, y> = QQ[]
            sage: P(0).degree(x)
            -1
            sage: P(1).degree(x)
            0

        The following example is inspired by :issue:`11652`::

            sage: R.<p,q,t> = ZZ[]
            sage: poly = p + q^2 + t^3
            sage: poly = poly.polynomial(t)[0]
            sage: poly
            q^2 + p

        There is no canonical coercion from ``R`` to the parent of ``poly``, so
        this doesn't work::

            sage: poly.degree(q)
            Traceback (most recent call last):
            ...
            TypeError: argument is not coercible to the parent

        Using a non-canonical coercion does work, but we require this
        to be done explicitly, since it can lead to confusing results
        if done automatically::

            sage: poly.degree(poly.parent()(q))
            2
            sage: poly.degree(poly.parent()(p))
            1
            sage: T.<x,y> = ZZ[]
            sage: poly.degree(poly.parent()(x))  # noncanonical coercions can be confusing
            1

        The argument to degree has to be a generator::

            sage: pp = poly.parent().gen(0)
            sage: poly.degree(pp)
            1
            sage: poly.degree(pp+1)
            Traceback (most recent call last):
            ...
            TypeError: argument is not a generator

        Canonical coercions are used::

            sage: S = ZZ['p,q']
            sage: poly.degree(S.0)
            1
            sage: poly.degree(S.1)
            2

        Ensure that :issue:`37603` is fixed::

            sage: R.<x, y> = QQ[]
            sage: f = x + y + 1
            sage: type(f.degree())
            <class 'sage.rings.integer.Integer'>
            sage: type(f.degree(x))
            <class 'sage.rings.integer.Integer'>
            sage: type(f.degree(y))
            <class 'sage.rings.integer.Integer'>
        """
        cdef ring *r = self._parent_ring
        cdef poly *p = self._poly
        if not x:
            if std_grading:
                return Integer(self.total_degree(std_grading=True))
            return Integer(singular_polynomial_deg(p, NULL, r))

        if not x.parent() is self.parent():
            try:
                x = self.parent().coerce(x)
            except TypeError:
                raise TypeError("argument is not coercible to the parent")
        if not x.is_gen():
            raise TypeError("argument is not a generator")

        return Integer(singular_polynomial_deg(p, x._poly, r))

    def total_degree(self, int std_grading=False):
        """
        Return the total degree of ``self``, which is the maximum degree
        of all monomials in ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 2*x*y^3*z^2
            sage: f.total_degree()
            6
            sage: f = 4*x^2*y^2*z^3
            sage: f.total_degree()
            7
            sage: f = 99*x^6*y^3*z^9
            sage: f.total_degree()
            18
            sage: f = x*y^3*z^6+3*x^2
            sage: f.total_degree()
            10
            sage: f = z^3+8*x^4*y^5*z
            sage: f.total_degree()
            10
            sage: f = z^9+10*x^4+y^8*x^2
            sage: f.total_degree()
            10

        A matrix term ordering changes the grading. To get the total degree
        using the standard grading, use ``std_grading=True``::

            sage: tord = TermOrder(matrix(3, [3,2,1,1,1,0,1,0,0]))
            sage: tord
            Matrix term order with matrix
            [3 2 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=tord)
            sage: f = x^2*y
            sage: f.total_degree()
            8
            sage: f.total_degree(std_grading=True)
            3

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: R(0).total_degree()
            -1
            sage: R(1).total_degree()
            0

        Ensure that :issue:`37603` is fixed::
            sage: R.<x,y,z> = QQ[]
            sage: f = x^4 + y + z
            sage: f.total_degree()
            4
            sage: type(f.total_degree())
            <class 'sage.rings.integer.Integer'>
        """
        cdef int i, result
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring

        if std_grading:
            result = 0
            while p:
                result = max(result, sum([p_GetExp(p,i,r) for i in range(1,r.N+1)]))
                p = pNext(p)
            return Integer(result)
        return Integer(singular_polynomial_deg(p, NULL, r))

    def degrees(self):
        """
        Return a tuple with the maximal degree of each variable in
        this polynomial.  The list of degrees is ordered by the order
        of the generators.

        EXAMPLES::

            sage: R.<y0,y1,y2> = PolynomialRing(QQ,3)
            sage: q = 3*y0*y1*y1*y2; q
            3*y0*y1^2*y2
            sage: q.degrees()
            (1, 2, 1)
            sage: (q + y0^5).degrees()
            (5, 2, 1)

        TESTS:

        Ensure that :issue:`37603` is fixed::

            sage: R.<x,y,z> = QQ[]
            sage: f = x^4 + y + z
            sage: f.degrees()
            (4, 1, 1)
            sage: type(f.degrees()[0])
            <class 'sage.rings.integer.Integer'>
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int i
        cdef list d = [0 for _ in range(r.N)]
        while p:
            for i from 0 <= i < r.N:
                d[i] = max(d[i],p_GetExp(p, i+1, r))
            p = pNext(p)
        return tuple(map(Integer, d))

    def coefficient(self, degrees):
        """
        Return the coefficient of the variables with the degrees
        specified in the python dictionary ``degrees``.
        Mathematically, this is the coefficient in the base ring
        adjoined by the variables of this ring not listed in
        ``degrees``.  However, the result has the same parent as this
        polynomial.

        This function contrasts with the function
        ``monomial_coefficient`` which returns the coefficient in the
        base ring of a monomial.

        INPUT:

        - ``degrees`` -- can be any of:
          - a dictionary of degree restrictions
          - a list of degree restrictions (with ``None`` in the unrestricted variables)
          - a monomial (very fast, but not as flexible)

        OUTPUT: element of the parent of this element

        .. NOTE::

           For coefficients of specific monomials, look at :meth:`monomial_coefficient`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f=x*y+y+5
            sage: f.coefficient({x:0,y:1})
            1
            sage: f.coefficient({x:0})
            y + 5
            sage: f=(1+y+y^2)*(1+x+x^2)
            sage: f.coefficient({x:0})
            y^2 + y + 1
            sage: f.coefficient([0,None])
            y^2 + y + 1
            sage: f.coefficient(x)
            y^2 + y + 1

        Note that exponents have all variables specified::

            sage: x.coefficient(x.exponents()[0])
            1
            sage: f.coefficient([1,0])
            1
            sage: f.coefficient({x:1,y:0})
            1

        Be aware that this may not be what you think! The physical
        appearance of the variable x is deceiving -- particularly if
        the exponent would be a variable. ::

            sage: f.coefficient(x^0) # outputs the full polynomial
            x^2*y^2 + x^2*y + x*y^2 + x^2 + x*y + y^2 + x + y + 1
            sage: R.<x,y> = GF(389)[]
            sage: f = x*y + 5
            sage: c = f.coefficient({x:0, y:0}); c
            5
            sage: parent(c)
            Multivariate Polynomial Ring in x, y over Finite Field of size 389

        AUTHOR:

        - Joel B. Mohler (2007.10.31)
        """
        cdef poly *_degrees = <poly*>0
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef poly *newp = p_ISet(0, r)
        cdef poly *newptemp
        cdef int i
        cdef int flag
        cdef int gens = self._parent.ngens()
        cdef int *exps = <int*>sig_malloc(sizeof(int)*gens)
        for i from 0<=i<gens:
            exps[i] = -1

        if isinstance(degrees, MPolynomial_libsingular) \
                and self._parent is (<MPolynomial_libsingular>degrees)._parent:
            _degrees = (<MPolynomial_libsingular>degrees)._poly
            if pLength(_degrees) != 1:
                raise TypeError("degrees must be a monomial")
            for i from 0<=i<gens:
                if p_GetExp(_degrees,i+1,r)!=0:
                    exps[i] = p_GetExp(_degrees,i+1,r)
        elif isinstance(degrees, list):
            for i from 0<=i<gens:
                if degrees[i] is None:
                    exps[i] = -1
                else:
                    exps[i] = int(degrees[i])
        elif isinstance(degrees, ETuple):
            for i in range(gens):
                exps[i] = int((<ETuple>degrees).get_exp(i))
        elif isinstance(degrees, dict):
            # Extract the ordered list of degree specifications from the dictionary
            poly_vars = self.parent().gens()
            for i from 0<=i<gens:
                try:
                    exps[i] = degrees[poly_vars[i]]
                except KeyError:
                    pass
        else:
            raise TypeError("The input degrees must be a dictionary of variables to exponents.")

        # Extract the monomials that match the specifications
        # this loop needs improvement
        while p:
            flag = 0
            for i from 0<=i<gens:
                if exps[i] != -1 and p_GetExp(p,i+1,r)!=exps[i]:
                    flag = 1
            if flag == 0:
                newptemp = p_LmInit(p,r)
                p_SetCoeff(newptemp,n_Copy(p_GetCoeff(p,r),r.cf),r)
                for i from 0<=i<gens:
                    if exps[i] != -1:
                        p_SetExp(newptemp,i+1,0,r)
                p_Setm(newptemp,r)
                newp = p_Add_q(newp,newptemp,r)
            p = pNext(p)

        sig_free(exps)
        return new_MP(self._parent, newp)

    def monomial_coefficient(self, MPolynomial_libsingular mon):
        """
        Return the coefficient in the base ring of the monomial mon in
        ``self``, where mon must have the same parent as ``self``.

        This function contrasts with the function ``coefficient``
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` -- a monomial

        OUTPUT: coefficient in base ring

        .. SEEALSO::

            For coefficients in a base ring of fewer variables,
            look at ``coefficient``.

        EXAMPLES::

            sage: P.<x,y> = QQ[]

            The parent of the return is a member of the base ring.
            sage: f = 2 * x * y
            sage: c = f.monomial_coefficient(x*y); c
            2
            sage: c.parent()
            Rational Field

            sage: f = y^2 + y^2*x - x^9 - 7*x + 5*x*y
            sage: f.monomial_coefficient(y^2)
            1
            sage: f.monomial_coefficient(x*y)
            5
            sage: f.monomial_coefficient(x^9)
            -1
            sage: f.monomial_coefficient(x^10)
            0
        """
        cdef poly *p = self._poly
        cdef poly *m = mon._poly
        cdef ring *r = self._parent_ring

        if mon._parent is not self._parent:
            raise TypeError("mon must have same parent as self.")

        while p:
            if p_ExpVectorEqual(p, m, r) == 1:
                return si2sa(p_GetCoeff(p, r), r, self._parent._base)
            p = pNext(p)

        return self._parent._base._zero_element

    def monomial_coefficients(self):
        """
        Return a dictionary representing ``self``. This dictionary is in
        the same format as the generic MPolynomial: The dictionary
        consists of ``ETuple:coefficient`` pairs.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 2*x*y^3*z^2 + 1/7*x^2 + 2/3
            sage: f.monomial_coefficients()
            {(0, 0, 0): 2/3, (1, 3, 2): 2, (2, 0, 0): 1/7}

        ``dict`` is an alias::

            sage: f.dict()
            {(0, 0, 0): 2/3, (1, 3, 2): 2, (2, 0, 0): 1/7}
        """
        cdef poly *p
        cdef ring *r = self._parent_ring
        cdef int n
        cdef int v
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        p = self._poly
        cdef dict d, pd = dict()
        while p:
            d = dict()
            for v from 1 <= v <= r.N:
                n = p_GetExp(p,v,r)
                if n != 0:
                    d[v-1] = n

            pd[ETuple(d,r.N)] = si2sa(p_GetCoeff(p, r), r, base)

            p = pNext(p)
        return pd

    dict = monomial_coefficients

    def iterator_exp_coeff(self, as_ETuples=True):
        """
        Iterate over ``self`` as pairs of ((E)Tuple, coefficient).

        INPUT:

        - ``as_ETuples`` -- boolean (default: ``True``); if ``True`` iterate
          over pairs whose first element is an ETuple, otherwise as a tuples

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: f = a*c^3 + a^2*b + 2*b^4
            sage: list(f.iterator_exp_coeff())
            [((0, 4, 0), 2), ((1, 0, 3), 1), ((2, 1, 0), 1)]
            sage: list(f.iterator_exp_coeff(as_ETuples=False))
            [((0, 4, 0), 2), ((1, 0, 3), 1), ((2, 1, 0), 1)]

            sage: R.<a,b,c> = PolynomialRing(QQ, 3, order='lex')
            sage: f = a*c^3 + a^2*b + 2*b^4
            sage: list(f.iterator_exp_coeff())
            [((2, 1, 0), 1), ((1, 0, 3), 1), ((0, 4, 0), 2)]
        """
        cdef poly *p
        cdef ring *r = self._parent_ring
        cdef int n
        cdef int v
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        p = self._poly
        cdef dict d
        while p:
            d = dict()
            for v from 1 <= v <= r.N:
                n = p_GetExp(p,v,r)
                if n != 0:
                    d[v-1] = n

            exp = ETuple(d,r.N)
            if as_ETuples:
                yield (exp, si2sa(p_GetCoeff(p, r), r, base))
            else:
                yield (tuple(exp), si2sa(p_GetCoeff(p, r), r, base))
            p = pNext(p)

    cpdef long number_of_terms(self) noexcept:
        """
        Return the number of nonzero coefficients of this polynomial.

        This is also called weight, :meth:`hamming_weight` or sparsity.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: f = x^3 - y
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x+y)^100
            sage: f.number_of_terms()
            101

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        cdef long w = 0
        p = self._poly
        while p:
            p = pNext(p)
            w += 1
        return w

    hamming_weight = number_of_terms

    cdef long _hash_c(self) except -1:
        """
        See ``self.__hash__``
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int n
        cdef int v
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        var_name_hash = [hash(vn) for vn in self._parent.variable_names()]
        cdef long c_hash
        while p:
            c_hash = hash(si2sa(p_GetCoeff(p, r), r, base))
            if c_hash != 0: # this is always going to be true, because we are sparse (correct?)
                # Hash (self[i], gen_a, exp_a, gen_b, exp_b, gen_c, exp_c, ...) as a tuple according to the algorithm.
                # I omit gen,exp pairs where the exponent is zero.
                result_mon = c_hash
                for v from 1 <= v <= r.N:
                    n = p_GetExp(p,v,r)
                    if n!=0:
                        result_mon = (1000003 * result_mon) ^ var_name_hash[v-1]
                        result_mon = (1000003 * result_mon) ^ n
                result += result_mon

            p = pNext(p)
        if result == -1:
            return -2
        return result

    def __getitem__(self, x):
        """
        Same as ``self.monomial_coefficient`` but for exponent vectors.

        INPUT:

        - ``x`` -- tuple or, in case of a single-variable MPolynomial
          ring x can also be an integer

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = -10*x^3*y + 17*x*y
            sage: f[3,1]
            -10
            sage: f[1,1]
            17
            sage: f[0,1]
            0

            sage: R.<x> = PolynomialRing(GF(7), implementation='singular'); R
            Multivariate Polynomial Ring in x over Finite Field of size 7
            sage: f = 5*x^2 + 3; f
            -2*x^2 + 3
            sage: f[2]
            5
        """
        cdef poly *m
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int i

        if isinstance(x, MPolynomial_libsingular):
            return self.monomial_coefficient(x)
        if not isinstance(x, tuple):
            try:
                x = tuple(x)
            except TypeError:
                x = (x,)

        if len(x) != self._parent.ngens():
            raise TypeError("x must have length self.ngens()")

        m = p_ISet(1,r)
        i = 1
        for e in x:
            overflow_check(e, r)
            p_SetExp(m, i, int(e), r)
            i += 1
        p_Setm(m, r)

        while p:
            if p_ExpVectorEqual(p, m, r) == 1:
                p_Delete(&m,r)
                return si2sa(p_GetCoeff(p, r), r, self._parent._base)
            p = pNext(p)

        p_Delete(&m,r)
        return self._parent._base._zero_element

    def __iter__(self):
        """
        Facilitates iterating over the monomials of self,
        returning tuples of the form ``(coeff, mon)`` for each
        nonzero monomial.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: f = 23*x^6*y^7 + x^3*y+6*x^7*z
            sage: list(f)
            [(23, x^6*y^7), (6, x^7*z), (1, x^3*y)]
            sage: list(R.zero())
            []

            sage: R.<x,y,z> = PolynomialRing(QQ, order='lex')
            sage: f = 23*x^6*y^7 + x^3*y+6*x^7*z
            sage: list(f)
            [(6, x^7*z), (23, x^6*y^7), (1, x^3*y)]
        """
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring
        if _ring != currRing: rChangeCurrRing(_ring)
        base = parent._base
        cdef poly *t
        cdef poly *p = p_Copy(self._poly, _ring)

        while p:
            t = pNext(p)
            p.next = NULL
            coeff = si2sa(p_GetCoeff(p, _ring), _ring, base)
            p_SetCoeff(p, n_Init(1,_ring.cf), _ring)
            p_Setm(p, _ring)
            yield (coeff, new_MP(parent, p))
            p = t

    def exponents(self, as_ETuples=True):
        """
        Return the exponents of the monomials appearing in this
        polynomial.

        INPUT:

        - ``as_ETuples`` -- boolean (default: ``True``); if ``True`` returns
          the result as a list of ETuples, otherwise returns a list of tuples

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: f = a^3 + b + 2*b^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
            sage: f.exponents(as_ETuples=False)
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        cdef int v
        cdef list pl, ml

        pl = list()
        ml = list(range(r.N))
        if as_ETuples:
            while p:
                for v from 1 <= v <= r.N:
                    ml[v-1] = p_GetExp(p,v,r)
                pl.append(ETuple(ml))
                p = pNext(p)
        else:
            while p:
                for v from 1 <= v <= r.N:
                    ml[v-1] = p_GetExp(p,v,r)
                pl.append(tuple(ml))
                p = pNext(p)
        return pl

    def inverse_of_unit(self):
        """
        Return the inverse of this polynomial if it is a unit.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: x.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: Element is not a unit.

            sage: R(1/2).inverse_of_unit()
            2
        """
        cdef ring *_ring = self._parent_ring
        if _ring != currRing:
            rChangeCurrRing(_ring)

        if not (p_IsUnit(self._poly,_ring)):
            raise ArithmeticError("Element is not a unit.")

        sig_on()
        cdef MPolynomial_libsingular r = new_MP(self._parent, p_NSet(n_Invers(p_GetCoeff(self._poly, _ring),_ring.cf),_ring))
        sig_off()
        return r

    def is_homogeneous(self):
        """
        Return ``True`` if this polynomial is homogeneous.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(RationalField(), 2)
            sage: (x+y).is_homogeneous()
            True
            sage: (x.parent()(0)).is_homogeneous()
            True
            sage: (x+y^2).is_homogeneous()
            False
            sage: (x^2 + y^2).is_homogeneous()
            True
            sage: (x^2 + y^2*x).is_homogeneous()
            False
            sage: (x^2*y + y^2*x).is_homogeneous()
            True
        """
        cdef ring *_ring = self._parent_ring
        if _ring != currRing:
            rChangeCurrRing(_ring)
        return bool(p_IsHomogeneous(self._poly,_ring))

    cpdef _homogenize(self, int var):
        """
        Return ``self`` if ``self`` is homogeneous.  Otherwise return
        a homogenized polynomial constructed by modifying the degree
        of the variable with index ``var``.

        INPUT:

        - ``var`` -- integer indicating which variable to use to
          homogenize (``0 <= var < parent(self).ngens()``)

        OUTPUT: a multivariate polynomial

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: f = x^2 + y + 1 + 5*x*y^10
            sage: g = f.homogenize('z'); g # indirect doctest
            5*x*y^10 + x^2*z^9 + y*z^10 + z^11
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: f._homogenize(0)
            2*x^11 + x^10*y + 5*x*y^10

        SEE: ``self.homogenize``
        """
        cdef ring *_ring = self._parent_ring
        cdef MPolynomialRing_libsingular parent = self._parent

        if self.is_homogeneous():
            return self

        if _ring != currRing:
            rChangeCurrRing(_ring)

        if var < parent._ring.N:
            return new_MP(parent, p_Homogen(self._poly, var + 1, _ring))
        raise TypeError("var must be < self.parent().ngens()")

    def is_monomial(self):
        """
        Return ``True`` if this polynomial is a monomial.  A monomial
        is defined to be a product of generators with coefficient 1.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: x.is_monomial()
            True
            sage: (2*x).is_monomial()
            False
            sage: (x*y).is_monomial()
            True
            sage: (x*y + x).is_monomial()
            False
            sage: P(2).is_monomial()
            False
            sage: P.zero().is_monomial()
            False
        """
        cdef poly *_p
        cdef ring *_ring
        cdef number *_n
        _ring = self._parent_ring

        if self._poly == NULL:
            return False

        if _ring != currRing:
            rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        ret = bool((not self._poly.next) and _ring.cf.cfIsOne(_n, _ring.cf))

        p_Delete(&_p, _ring)
        return ret

    def is_term(self):
        """
        Return ``True`` if ``self`` is a term, which we define to be a
        product of generators times some coefficient, which need
        not be 1.

        Use :meth:`is_monomial()` to require that the coefficient be 1.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: x.is_term()
            True
            sage: (2*x).is_term()
            True
            sage: (x*y).is_term()
            True
            sage: (x*y + x).is_term()
            False
            sage: P(2).is_term()
            True
            sage: P.zero().is_term()
            False
        """
        return self._poly != NULL and self._poly.next == NULL

    def subs(self, fixed=None, **kw):
        """
        Fixes some given variables in a given multivariate polynomial
        and returns the changed multivariate polynomials. The
        polynomial itself is not affected.  The variable,value pairs
        for fixing are to be provided as dictionary of the form
        ``{variable:value}``.

        This is a special case of evaluating the polynomial with some
        of the variables constants and the others the original
        variables, but should be much faster if only few variables are
        to be fixed.

        INPUT:

        - ``fixed`` -- (optional) dict with variable:value pairs
        - ``**kw`` -- names parameters

        OUTPUT: a new multivariate polynomial

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f(5,y)
            25*y^2 + y + 30
            sage: f.subs({x: 5})
            25*y^2 + y + 30
            sage: f.subs(x=5)
            25*y^2 + y + 30

            sage: P.<x,y,z> = PolynomialRing(GF(2), 3)
            sage: f = x + y + 1
            sage: f.subs({x:y+1})
            0
            sage: f.subs(x=y)
            1
            sage: f.subs(x=x)
            x + y + 1
            sage: f.subs({x: z})
            y + z + 1
            sage: f.subs(x=z + 1)
            y + z

            sage: f.subs(x=1/y)
            (y^2 + y + 1)/y
            sage: f.subs({x: 1/y})
            (y^2 + y + 1)/y

        The parameters are substituted in order and without side effects::

            sage: R.<x,y>=QQ[]
            sage: g=x+y
            sage: g.subs({x:x+1,y:x*y})
            x*y + x + 1
            sage: g.subs({x:x+1}).subs({y:x*y})
            x*y + x + 1
            sage: g.subs({y:x*y}).subs({x:x+1})
            x*y + x + y + 1

        ::

            sage: R.<x,y> = QQ[]
            sage: f = x + 2*y
            sage: f.subs(x=y,y=x)
            2*x + y

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: f = y
            sage: f.subs({y:x}).subs({x:z})
            z

        We test that we change the ring even if there is nothing to do::

            sage: P = QQ['x,y']
            sage: x = var('x')                                                          # needs sage.symbolic
            sage: parent(P.zero().subs(x=x))                                            # needs sage.symbolic
            Symbolic Ring

        We are catching overflows::

            sage: R.<x,y> = QQ[]
            sage: n=100; f = x^n
            sage: try:
            ....:     f.subs(x = x^n)
            ....:     print("no overflow")
            ....: except OverflowError:
            ....:     print("overflow")
            x^10000
            no overflow

            sage: n = 100000
            sage: try:
            ....:     f = x^n
            ....:     f.subs(x = x^n)
            ....:     print("no overflow")
            ....: except OverflowError:
            ....:     print("overflow")
            overflow

        Check that there is no more segmentation fault if the polynomial gets 0
        in the middle of a substitution (:issue:`17785`)::

            sage: R.<x,y,z> = QQ[]
            sage: for vx in [0,x,y,z]:
            ....:     for vy in [0,x,y,z]:
            ....:         for vz in [0,x,y,z]:
            ....:             d = {x:vx, y:vy, z:vz}
            ....:             ds = {'x': vx, 'y': vy, 'z': vz}
            ....:             assert x.subs(d) == x.subs(**ds) == vx
            ....:             assert y.subs(d) == y.subs(**ds) == vy
            ....:             assert z.subs(d) == z.subs(**ds) == vz
            ....:             assert (x+y).subs(d) == (x+y).subs(**ds) == vx+vy

        Check that substitution doesn't crash in transcendental extensions::

            sage: F = PolynomialRing(QQ,'c,d').fraction_field()
            sage: F.inject_variables()
            Defining c, d
            sage: R.<x,y,z> = F[]
            sage: f = R(d*z^2 + c*y*z^2)
            sage: f.subs({x:z^2,y:1})
            (c + d)*z^2
            sage: f.subs({z:x+1})
            c*x^2*y + d*x^2 + (2*c)*x*y + (2*d)*x + c*y + d
        """
        cdef int mi, i
        cdef bint change_ring = False   # indicates the need to change the parent ring
        cdef bint need_map = False

        cdef unsigned long degree = 0
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring

        if (_ring != currRing):
            rChangeCurrRing(_ring)

        cdef poly *_p = p_Copy(self._poly, _ring)
        cdef poly *_f

        cdef ideal *to_id = idInit(_ring.N, 1)
        cdef ideal *from_id
        cdef ideal *res_id
        cdef dict gd

        cdef list g = list(parent.gens())
        cdef list items = []

        if not change_ring:
            if fixed is not None:
                for m, v in fixed.items():
                    if isinstance(m, (int, Integer)):
                        mi = m + 1
                    elif isinstance(m, MPolynomial_libsingular) and m.parent() is parent:
                        for i from 0 < i <= _ring.N:
                            if p_GetExp((<MPolynomial_libsingular> m)._poly, i, _ring) != 0:
                                mi = i
                                break
                        if i > _ring.N:
                            id_Delete(&to_id, _ring)
                            p_Delete(&_p, _ring)
                            raise TypeError("key does not match")
                    else:
                        id_Delete(&to_id, _ring)
                        p_Delete(&_p, _ring)
                        raise TypeError("keys do not match self's parent")

                    items.append((g[mi - 1], v))
            if kw:
                gd = parent.gens_dict(copy=False)
                for m, v in kw.items():
                    items.append((gd[m], v))
            for m, v in items:
                if _p == NULL:
                    change_ring = True
                    break
                i = g.index(m)
                try:
                    v = parent.coerce(v)
                except TypeError: # give up, evaluate symbolically
                    id_Delete(&to_id, _ring)
                    p_Delete(&_p, _ring)

                    gg = list(parent.gens())
                    for m, v in items:
                        gg[g.index(m)] = v
                    y = parent.base_ring().zero()
                    for m, c in self.monomial_coefficients().items():
                        y += c*mul([gg[i] ** m[i] for i in m.nonzero_positions()])
                    return y
                _f = (<MPolynomial_libsingular> v)._poly
                if p_IsConstant(_f, _ring):
                    singular_polynomial_subst(&_p, i, _f, _ring)
                else:
                    need_map = True
                    degree = (<unsigned long> p_GetExp(_p, i + 1, _ring)) * (<unsigned long> p_GetMaxExp(_f, _ring))
                    if degree > _ring.bitmask:
                        id_Delete(&to_id, _ring)
                        p_Delete(&_p, _ring)
                        raise OverflowError("exponent overflow (%d)" % (degree,))
                    to_id.m[i] = p_Copy(_f, _ring)

        if not change_ring and need_map:
            for mi from 0 <= mi < _ring.N:
                if to_id.m[mi] == NULL:
                    to_id.m[mi] = p_ISet(1,_ring)
                    p_SetExp(to_id.m[mi], mi + 1, 1, _ring)
                    p_Setm(to_id.m[mi], _ring)

            from_id = idInit(1, 1)
            from_id.m[0] = _p

            rChangeCurrRing(_ring)
            res_id = fast_map_common_subexp(from_id, _ring, to_id, _ring)
            _p = res_id.m[0]

            from_id.m[0] = NULL
            res_id.m[0] = NULL

            id_Delete(&from_id, _ring)
            id_Delete(&res_id, _ring)

        id_Delete(&to_id, _ring)

        if not change_ring:
            return new_MP(parent,_p)

        # finally change the parent of the result

        res_parent = coercion_model.common_parent(parent._base, *[v for _, v in items])

        if _p == NULL:
            return res_parent(0)

        if p_LmIsConstant(_p, _ring):
            res = si2sa(p_GetCoeff(_p, _ring), _ring, parent._base)
            p_Delete(&_p, _ring)  # safe to delete; res contains copy
        else:
            res = new_MP(parent, _p)
        if parent(res) is not res_parent:
            res = res_parent(res)
        return res

    def monomials(self):
        """
        Return the list of monomials in ``self``. The returned list is
        decreasingly ordered by the term ordering of ``self.parent()``.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = x + 3/2*y*z^2 + 2/3
            sage: f.monomials()
            [y*z^2, x, 1]
            sage: f = P(3/2)
            sage: f.monomials()
            [1]

        TESTS::

            sage: P.<x,y,z> = QQ[]
            sage: f = x
            sage: f.monomials()
            [x]

        Check if :issue:`12706` is fixed::

            sage: f = P(0)
            sage: f.monomials()
            []

        Check if :issue:`7152` is fixed::

            sage: x = polygen(ZZ, 'x')
            sage: K.<rho> = NumberField(x**2 + 1)
            sage: R.<x,y> = QQ[]
            sage: p = rho*x
            sage: q = x
            sage: p.monomials()
            [x]
            sage: q.monomials()
            [x]
            sage: p.monomials()
            [x]
        """
        cdef list l = []
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *_ring = parent._ring
        if _ring != currRing:
            rChangeCurrRing(_ring)
        cdef poly *p = p_Copy(self._poly, _ring)
        cdef poly *t

        while p:
            t = pNext(p)
            p.next = NULL
            p_SetCoeff(p, n_Init(1,_ring.cf), _ring)
            p_Setm(p, _ring)
            l.append( new_MP(parent,p) )
            p = t

        return l

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate
        polynomial.

        EXAMPLES::

            sage: P.<x, y> = QQ[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        cdef poly *p = self._poly
        cdef ring *r = self._parent_ring
        if p == NULL:
            return self._parent._base._zero_element

        while p.next:
            p = pNext(p)

        if p_LmIsConstant(p, r):
            return si2sa(p_GetCoeff(p, r), r, self._parent._base)
        else:
            return self._parent._base._zero_element

    def univariate_polynomial(self, R=None):
        """
        Return a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:

        - ``R`` -- (default: ``None``) PolynomialRing

        If this polynomial is not in at most one variable, then a
        :exc:`ValueError` exception is raised.  This is checked using the
        :meth:`is_univariate()` method.  The new Polynomial is over
        the same base ring as the given ``MPolynomial`` and in the
        variable ``x`` if no ring ``R`` is provided.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial ()
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial(PolynomialRing(QQ,'z'))
            700*z^2 - 2*z + 305

        Here's an example with a constant multivariate polynomial::

            sage: g = R(1)
            sage: h = g.univariate_polynomial(); h
            1
            sage: h.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        cdef poly *p = self._poly
        cdef poly *p2 = self._poly
        cdef ring *r = self._parent_ring
        cdef long pTotDegMax

        k = self.base_ring()

        if not self.is_univariate():
            raise TypeError("polynomial must involve at most one variable")

        #construct ring if none
        if R is None:
            if self.is_constant():
                R = self.base_ring()['x']
            else:
                R = self.base_ring()[str(self.variables()[0])]

        zero = k(0)

        if r != currRing:
            rChangeCurrRing(r)

        pTotDegMax = -1
        while p2:
            pTotDegMax = max(pTotDegMax, p_Totaldegree(p2, r))
            p2 = pNext(p2)

        coefficients = [zero] * (pTotDegMax + 1)
        while p:
            pTotDeg = p_Totaldegree(p, r)
            if (pTotDeg >= len(coefficients) or pTotDeg < 0):
                raise IndexError("list index("+str(pTotDeg)+" out of range(0-"+str(len(coefficients))+")")
            coefficients[pTotDeg] = si2sa(p_GetCoeff(p, r), r, k)
            p = pNext(p)

        return R(coefficients)

    def is_univariate(self):
        """
        Return ``True`` if ``self`` is a univariate polynomial, that is if
        ``self`` contains only one variable.

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x^2 + 1
            sage: f.is_univariate()
            True
            sage: f = y*x^2 + 1
            sage: f.is_univariate()
            False
            sage: f = P(0)
            sage: f.is_univariate()
            True
        """
        return bool(len(self._variable_indices_(sort=False))<2)

    def _variable_indices_(self, sort=True):
        """
        Return the indices of all variables occurring in ``self``.  This
        index is the index as Sage uses them (starting at zero), not
        as SINGULAR uses them (starting at one).

        INPUT:

        - ``sort`` -- specifies whether the indices shall be sorted

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x*z^2 + z + 1
            sage: f._variable_indices_()
            [0, 2]
        """
        cdef poly *p
        cdef ring *r = self._parent_ring
        cdef int i
        s = set()
        p = self._poly
        while p:
            for i from 1 <= i <= r.N:
                if p_GetExp(p,i,r):
                    s.add(i-1)
            p = pNext(p)
        if sort:
            return sorted(s)
        else:
            return list(s)

    def variables(self):
        """
        Return a tuple of all variables occurring in ``self``.

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x*z^2 + z + 1
            sage: f.variables()
            (x, z)
        """
        cdef poly *p
        cdef poly *v
        cdef ring *r = self._parent_ring
        if r != currRing:
            rChangeCurrRing(r)
        cdef int i
        l = list()
        si = set()
        p = self._poly
        while p:
            for i from 1 <= i <= r.N:
                if i not in si and p_GetExp(p,i,r):
                    v = p_ISet(1,r)
                    p_SetExp(v, i, 1, r)
                    p_Setm(v, r)
                    l.append(new_MP(self._parent, v))
                    si.add(i)
            p = pNext(p)
        return tuple(sorted(l,reverse=True))

    def variable(self, i=0):
        """
        Return the `i`-th variable occurring in ``self``. The index `i` is the
        index in ``self.variables()``.

        EXAMPLES::

            sage: P.<x,y,z> = GF(2)[]
            sage: f = x*z^2 + z + 1
            sage: f.variables()
            (x, z)
            sage: f.variable(1)
            z
        """
        return self.variables()[i]

    def nvariables(self):
        """
        Return the number variables in this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127))
            sage: f = x*y + z
            sage: f.nvariables()
            3
            sage: f = x + y
            sage: f.nvariables()
            2
        """
        return len(self._variable_indices_(sort=False))

    cpdef is_constant(self):
        """
        Return ``True`` if this polynomial is constant.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127))
            sage: x.is_constant()
            False
            sage: P(1).is_constant()
            True
        """
        return bool(p_IsConstant(self._poly, self._parent_ring))

    def lm(MPolynomial_libsingular self):
        """
        Return the lead monomial of ``self`` with respect to the term
        order of ``self.parent()``. In Sage a monomial is a product of
        variables in some power without a coefficient.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(GF(7), 3, order='lex')
            sage: f = x^1*y^2 + y^3*z^4
            sage: f.lm()
            x*y^2
            sage: f = x^3*y^2*z^4 + x^3*y^2*z^1
            sage: f.lm()
            x^3*y^2*z^4

            sage: R.<x,y,z>=PolynomialRing(QQ, 3, order='deglex')
            sage: f = x^1*y^2*z^3 + x^3*y^2*z^0
            sage: f.lm()
            x*y^2*z^3
            sage: f = x^1*y^2*z^4 + x^1*y^1*z^5
            sage: f.lm()
            x*y^2*z^4

            sage: R.<x,y,z>=PolynomialRing(GF(127), 3, order='degrevlex')
            sage: f = x^1*y^5*z^2 + x^4*y^1*z^3
            sage: f.lm()
            x*y^5*z^2
            sage: f = x^4*y^7*z^1 + x^4*y^2*z^3
            sage: f.lm()
            x^4*y^7*z
        """
        cdef poly *_p
        cdef ring *_ring = self._parent_ring
        if self._poly == NULL:
            return self._parent._zero_element
        _p = p_Head(self._poly, _ring)
        p_SetCoeff(_p, n_Init(1,_ring.cf), _ring)
        p_Setm(_p,_ring)
        return new_MP(self._parent, _p)

    def lc(MPolynomial_libsingular self):
        """
        Leading coefficient of this polynomial with respect to the
        term order of ``self.parent()``.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(GF(7), 3, order='lex')
            sage: f = 3*x^1*y^2 + 2*y^3*z^4
            sage: f.lc()
            3

            sage: f = 5*x^3*y^2*z^4 + 4*x^3*y^2*z^1
            sage: f.lc()
            5
        """

        cdef poly *_p
        cdef ring *_ring = self._parent_ring
        cdef number *_n

        if self._poly == NULL:
            return self._parent._base._zero_element

        if _ring != currRing:
            rChangeCurrRing(_ring)

        _p = p_Head(self._poly, _ring)
        _n = p_GetCoeff(_p, _ring)

        ret =  si2sa(_n, _ring, self._parent._base)
        p_Delete(&_p, _ring)
        return ret

    def lt(MPolynomial_libsingular self):
        """
        Leading term of this polynomial. In Sage a term is a product
        of variables in some power and a coefficient.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(GF(7), 3, order='lex')
            sage: f = 3*x^1*y^2 + 2*y^3*z^4
            sage: f.lt()
            3*x*y^2

            sage: f = 5*x^3*y^2*z^4 + 4*x^3*y^2*z^1
            sage: f.lt()
            -2*x^3*y^2*z^4
        """
        if self._poly == NULL:
            return self._parent._zero_element
        return new_MP(self._parent, p_Head(self._poly, self._parent_ring))

    def is_zero(self):
        """
        Return ``True`` if this polynomial is zero.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: x.is_zero()
            False
            sage: (x - x).is_zero()
            True
        """
        if self._poly is NULL:
            return True
        else:
            return False

    def __bool__(self):
        """
        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: bool(x) # indirect doctest
            True
            sage: bool(x-x)
            False
        """
        if self._poly:
            return True
        else:
            return False

    cpdef _floordiv_(self, right):
        """
        Perform division with remainder and return the quotient.

        INPUT:

        - ``right`` -- something coercible to an MPolynomial_libsingular
          in ``self.parent()``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y,z> = GF(32003)[]
            sage: f = y*x^2 + x + 1
            sage: f//x
            x*y + 1
            sage: f//y
            x^2

            sage: P.<x,y> = ZZ[]
            sage: x//y
            0
            sage: (x+y)//y
            1

            sage: P.<x,y> = QQ[]
            sage: (x+y)//y
            1
            sage: (x)//y
            0

            sage: P.<x,y> = Zmod(1024)[]
            sage: (x+y)//x
            1
            sage: (x+y)//(2*x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Division of multivariate polynomials over non fields by non-monomials not implemented.

        TESTS::

            sage: P.<x,y> = ZZ[]
            sage: p = 3*(-x^8*y^2 - x*y^9 + 6*x^8*y + 17*x^2*y^6 - x^3*y^2)
            sage: q = 7*(x^2 + x*y + y^2 + 1)
            sage: p*q//q == p
            True
            sage: p*q//p == q
            True
        """
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef MPolynomial_libsingular _right = <MPolynomial_libsingular>right
        cdef ring *r = self._parent_ring
        cdef poly *quo
        cdef poly *temp
        cdef poly *p

        if _right._poly == NULL:
            raise ZeroDivisionError
        elif p_IsOne(_right._poly, r):
            return self

        if n_GetChar(r.cf) > 1<<29:
            raise NotImplementedError("Division of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

        if r.cf.type != n_unknown:
            if (singular_polynomial_length_bounded(_right._poly, 2) == 1
                    and r.cf.cfIsOne(p_GetCoeff(_right._poly, r), r.cf)):
                p = self._poly
                quo = p_ISet(0,r)
                while p:
                    if p_DivisibleBy(_right._poly, p, r):
                        temp = p_MDivide(p, _right._poly, r)
                        p_SetCoeff0(temp, n_Copy(p_GetCoeff(p, r), r.cf), r)
                        quo = p_Add_q(quo, temp, r)
                    p = pNext(p)
                return new_MP(parent, quo)
            if r.cf.type == n_Znm or r.cf.type == n_Zn or r.cf.type == n_Z2m:
                raise NotImplementedError("Division of multivariate polynomials over non fields by non-monomials not implemented.")

        count = singular_polynomial_length_bounded(self._poly, 15)

        # fast in the most common case where the division is exact; returns zero otherwise
        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
            sig_on()
        quo = p_Divide(p_Copy(self._poly, r), p_Copy(_right._poly, r), r)
        if count >= 15:
            sig_off()

        if quo == NULL:
            if r.cf.type == n_Z:
                P = parent.change_ring(QQ)
                f = (<MPolynomial_libsingular>P(self))._floordiv_(P(right))
                return parent(sum([c.floor() * m for c, m in f]))
            else:
                sig_on()
                quo = singclap_pdivide(self._poly, _right._poly, r)
                sig_off()

        return new_MP(parent, quo)

    def factor(self, proof=None):
        """
        Return the factorization of this polynomial.

        INPUT:

        - ``proof`` -- ignored

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            x^5 + 2*x^3*y^2 + x^4 + 2*x^2*y^2 + x^3 + 2*x*y^2
            sage: F = f.factor(); F
            x * (x^2 + x + 1) * (x^2 + 2*y^2)

        Next we factor the same polynomial, but over the finite field
        of order 3.::

            sage: R.<x, y> = GF(3)[]
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            x^5 - x^3*y^2 + x^4 - x^2*y^2 + x^3 - x*y^2
            sage: F = f.factor()
            sage: F  # order is somewhat random
            (-1) * x * (-x + y) * (x + y) * (x - 1)^2

        Next we factor a polynomial, but over a finite field of order 9.::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(3^2)
            sage: R.<x, y> = K[]
            sage: f = (x^3 + 2*a*y^2*x) * (x^2 + x + 1); f
            x^5 + (-a)*x^3*y^2 + x^4 + (-a)*x^2*y^2 + x^3 + (-a)*x*y^2
            sage: F = f.factor(); F
            ((-a)) * x * (x - 1)^2 * ((-a + 1)*x^2 + y^2)
            sage: f - F
            0

        Next we factor a polynomial over a number field.::

            sage: # needs sage.rings.number_field
            sage: p = polygen(ZZ, 'p')
            sage: K.<s> = NumberField(p^3 - 2)
            sage: KXY.<x,y> = K[]
            sage: factor(x^3 - 2*y^3)
            (x + (-s)*y) * (x^2 + s*x*y + (s^2)*y^2)
            sage: k = (x^3-2*y^3)^5*(x+s*y)^2*(2/3 + s^2)
            sage: k.factor()
            ((s^2 + 2/3)) * (x + s*y)^2 * (x + (-s)*y)^5 * (x^2 + s*x*y + (s^2)*y^2)^5

        This shows that issue :issue:`2780` is fixed, i.e. that the unit
        part of the factorization is set correctly::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 1)
            sage: R.<y, z> = PolynomialRing(K)
            sage: f = 2*y^2 + 2*z^2
            sage: F = f.factor(); F.unit()
            2

        Another example::

            sage: R.<x,y,z> = GF(32003)[]                                               # needs sage.rings.finite_rings
            sage: f = 9*(x-1)^2*(y+z)                                                   # needs sage.rings.finite_rings
            sage: f.factor()                                                            # needs sage.rings.finite_rings
            (9) * (y + z) * (x - 1)^2

            sage: R.<x,w,v,u> = QQ['x','w','v','u']
            sage: p = (4*v^4*u^2 - 16*v^2*u^4 + 16*u^6 - 4*v^4*u + 8*v^2*u^3 + v^4)
            sage: p.factor()
            (-2*v^2*u + 4*u^3 + v^2)^2
            sage: R.<a,b,c,d> = QQ[]
            sage: f =  (-2) * (a - d) * (-a + b) * (b - d) * (a - c) * (b - c) * (c - d)
            sage: F = f.factor(); F
            (-2) * (c - d) * (-b + c) * (b - d) * (-a + c) * (-a + b) * (a - d)
            sage: F[0][0]
            c - d
            sage: F.unit()
            -2

        Constant elements are factorized in the base rings. ::

            sage: P.<x,y> = ZZ[]
            sage: P(2^3*7).factor()
            2^3 * 7
            sage: P.<x,y> = GF(2)[]
            sage: P(1).factor()
            1

        Factorization for finite prime fields with characteristic
        `> 2^{29}` is not supported ::

            sage: q = 1073741789
            sage: T.<aa, bb> = PolynomialRing(GF(q))                                    # needs sage.rings.finite_rings
            sage: f = aa^2 + 12124343*bb*aa + 32434598*bb^2                             # needs sage.rings.finite_rings
            sage: f.factor()                                                            # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NotImplementedError: Factorization of multivariate polynomials
            over prime fields with characteristic > 2^29 is not implemented.

        Factorization over the integers is now supported, see :issue:`17840`::

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: f = 12 * (3*x*y + 4) * (5*x - 2) * (2*y + 7)^2
            sage: f.factor()
            2^2 * 3 * (2*y + 7)^2 * (5*x - 2) * (3*x*y + 4)
            sage: g = -12 * (x^2 - y^2)
            sage: g.factor()
            (-1) * 2^2 * 3 * (x - y) * (x + y)
            sage: factor(-4*x*y - 2*x + 2*y + 1)
            (-1) * (2*y + 1) * (2*x - 1)

        Factorization over non-integral domains is not supported ::

            sage: R.<x,y> = PolynomialRing(Zmod(4))
            sage: f = (2*x + 1) * (x^2 + x + 1)
            sage: f.factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: Factorization of multivariate polynomials
            over Ring of integers modulo 4 is not implemented.

        TESTS:

        This shows that :issue:`10270` is fixed::

            sage: R.<x,y,z> = GF(3)[]
            sage: f = x^2*z^2+x*y*z-y^2
            sage: f.factor()
            x^2*z^2 + x*y*z - y^2

        This checks that :issue:`11838` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(4,'a')
            sage: a = K.gens()[0]
            sage: R.<x,y> = K[]
            sage: p=x^8*y^3 + x^2*y^9 + a*x^9 + a*x*y^4
            sage: q=y^11 + (a)*y^10 + (a + 1)*x*y^3
            sage: f = p*q
            sage: f.factor()
            x * y^3 * (y^8 + a*y^7 + (a + 1)*x) * (x^7*y^3 + x*y^9 + a*x^8 + a*y^4)

        We test several examples which were known to return wrong
        results in the past (see :issue:`10902`)::

            sage: R.<x,y> = GF(2)[]
            sage: p = x^3*y^7 + x^2*y^6 + x^2*y^3
            sage: q = x^3*y^5
            sage: f = p*q
            sage: p.factor()*q.factor()
            x^5 * y^8 * (x*y^4 + y^3 + 1)
            sage: f.factor()
            x^5 * y^8 * (x*y^4 + y^3 + 1)
            sage: f.factor().expand() == f
            True

        ::

            sage: R.<x,y> = GF(2)[]
            sage: p = x^8 + y^8; q=x^2*y^4 + x
            sage: f = p*q
            sage: lf = f.factor()
            sage: f-lf
            0

        ::

            sage: R.<x,y> = GF(3)[]
            sage: p = -x*y^9 + x
            sage: q = -x^8*y^2
            sage: f = p*q
            sage: f
            x^9*y^11 - x^9*y^2
            sage: f.factor()
            y^2 * (y - 1)^9 * x^9
            sage: f - f.factor()
            0

        ::

            sage: R.<x,y> = GF(5)[]
            sage: p = x^27*y^9 + x^32*y^3 + 2*x^20*y^10 - x^4*y^24 - 2*x^17*y
            sage: q = -2*x^10*y^24 + x^9*y^24 - 2*x^3*y^30
            sage: f = p*q; f - f.factor()
            0

        ::

            sage: R.<x,y> = GF(7)[]
            sage: p = -3*x^47*y^24
            sage: q = -3*x^47*y^37 - 3*x^24*y^49 + 2*x^56*y^8 + 3*x^29*y^15 - x^2*y^33
            sage: f = p*q
            sage: f - f.factor()
            0

        The following examples used to give a Segmentation Fault, see
        :issue:`12918` and :issue:`13129`::

            sage: R.<x,y> = GF(2)[]
            sage: f = x^6 + x^5 + y^5 + y^4
            sage: f.factor()
            x^6 + x^5 + y^5 + y^4
            sage: f = x^16*y + x^10*y + x^9*y + x^6*y + x^5 + x*y + y^2
            sage: f.factor()
            x^16*y + x^10*y + x^9*y + x^6*y + x^5 + x*y + y^2

        Test :issue:`12928`::

            sage: R.<x,y> = GF(2)[]
            sage: p = x^2 + y^2 + x + 1
            sage: q = x^4 + x^2*y^2 + y^4 + x*y^2 + x^2 + y^2 + 1
            sage: factor(p*q)
            (x^2 + y^2 + x + 1) * (x^4 + x^2*y^2 + y^4 + x*y^2 + x^2 + y^2 + 1)

        Check that :issue:`13770` is fixed::

            sage: U.<y,t> = GF(2)[]
            sage: f = y*t^8 + y^5*t^2 + y*t^6 + t^7 + y^6 + y^5*t + y^2*t^4 + y^2*t^2 + y^2*t + t^3 + y^2 + t^2
            sage: l = f.factor()
            sage: l[0][0] == t^2 + y + t + 1 or l[1][0] == t^2 + y + t + 1
            True

        The following used to sometimes take a very long time or get
        stuck, see :issue:`12846`. These 100 iterations should take less
        than 1 second::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(4)
            sage: R.<x,y> = K[]
            sage: f = (a + 1)*x^145*y^84 + (a + 1)*x^205*y^17 + x^32*y^112 + x^92*y^45
            sage: for i in range(100):
            ....:     assert len(f.factor()) == 4

        Test for :issue:`20435`::

            sage: x,y = polygen(ZZ,'x,y')
            sage: p = x**2 - y**2
            sage: z = factor(p); z
            (x - y) * (x + y)
            sage: z[0][0].parent()
            Multivariate Polynomial Ring in x, y over Integer Ring

        Test for :issue:`17680`::

            sage: R.<a,r,v,n,g,f,h,o> = QQ[]
            sage: f = 248301045*a^2*r^10*n^2*o^10+570807000*a^2*r^9*n*o^9-137945025*a^2*r^8*n^2*o^8+328050000*a^2*r^8*o^8-253692000*a^2*r^7*n*o^7+30654450*a^2*r^6*n^2*o^6-109350000*a^2*r^6*o^6+42282000*a^2*r^5*n*o^5-3406050*a^2*r^4*n^2*o^4-22457088*a*r^2*v*n^2*o^6+12150000*a^2*r^4*o^4-3132000*a^2*r^3*n*o^3+189225*a^2*r^2*n^2*o^2+2495232*a*v*n^2*o^4-450000*a^2*r^2*o^2+87000*a^2*r*n*o-4205*a^2*n^2
            sage: len(factor(f))
            4

        Test for :issue:`17251`::

            sage: R.<z,a,b> = PolynomialRing(QQ)
            sage: N = -a^4*z^8 + 2*a^2*b^2*z^8 - b^4*z^8 - 16*a^3*b*z^7 + 16*a*b^3*z^7 + 28*a^4*z^6 - 56*a^2*b^2*z^6 + 28*b^4*z^6 + 112*a^3*b*z^5 - 112*a*b^3*z^5 - 70*a^4*z^4 + 140*a^2*b^2*z^4 - 70*b^4*z^4 - 112*a^3*b*z^3 + 112*a*b^3*z^3 + 28*a^4*z^2 - 56*a^2*b^2*z^2 + 28*b^4*z^2 + 16*a^3*b*z - 16*a*b^3*z - a^4 + 2*a^2*b^2 - b^4
            sage: N.factor()
            (-1) * (-a + b) * (a + b) * (-z^4*a + z^4*b - 4*z^3*a - 4*z^3*b + 6*z^2*a - 6*z^2*b + 4*z*a + 4*z*b - a + b) * (z^4*a + z^4*b - 4*z^3*a + 4*z^3*b - 6*z^2*a - 6*z^2*b + 4*z*a - 4*z*b + a + b)
        """
        cdef ring *_ring = self._parent_ring
        cdef poly *ptemp
        cdef intvec *iv
        cdef int *ivv
        cdef ideal *I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef int i

        if _ring!=currRing: rChangeCurrRing(_ring)

        if p_IsConstant(self._poly, _ring):
            return self.constant_coefficient().factor()

        if not self._parent._base.is_field():
            try:
                frac_field = self._parent._base.fraction_field()
                F = self.change_ring(frac_field).factor()
                FF = [(self._parent(f[0]), f[1]) for f in F]
                U = self._parent._base(F.unit()).factor()
                return Factorization(list(U) + FF, unit=U.unit())
            except Exception:
                raise NotImplementedError("Factorization of multivariate polynomials over %s is not implemented."%self._parent._base)

        if n_GetChar(_ring.cf) > 1<<29:
            raise NotImplementedError("Factorization of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

        # I make a temporary copy of the poly in self because singclap_factorize appears to modify it's parameter
        ptemp = p_Copy(self._poly, _ring)
        iv = NULL
        sig_on()
        if _ring != currRing:
            rChangeCurrRing(_ring)   # singclap_factorize
        I = singclap_factorize(ptemp, &iv, 0, _ring)
        sig_off()

        ivv = iv.ivGetVec()
        v = [(new_MP(parent, p_Copy(I.m[i], _ring)), ivv[i])
             for i in range(1, I.ncols)]
        v = [(f, m) for f, m in v if f != 0]  # we might have zero in there
        unit = new_MP(parent, p_Copy(I.m[0], _ring))

        F = Factorization(v, unit)
        F.sort()

        del iv
        id_Delete(&I, _ring)

        return F

    def lift(self, I):
        """
        Given an ideal ``I = (f_1,...,f_r)`` and some ``g (== self)`` in ``I``,
        find ``s_1,...,s_r`` such that ``g = s_1 f_1 + ... + s_r f_r``.

        A :exc:`ValueError` exception is raised if ``g (== self)`` does not belong to ``I``.

        EXAMPLES::

            sage: A.<x,y> = PolynomialRing(QQ,2,order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I)
            sage: M
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]
            sage: sum( map( mul , zip( M, I.gens() ) ) ) == f
            True

        Check that :issue:`13671` is fixed::

            sage: R.<x1,x2> = QQ[]
            sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
            sage: f = I.gen(0) + x2*I.gen(1)
            sage: f.lift(I)
            [1, x2]
            sage: (f+1).lift(I)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not in the ideal

        Check that we can work over the integers::

            sage: R.<x1,x2> = ZZ[]
            sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
            sage: f = I.gen(0) + x2*I.gen(1)
            sage: f.lift(I)
            [1, x2]
            sage: (f+1).lift(I)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not in the ideal
            sage: A.<x,y> = PolynomialRing(ZZ,2,order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I)
            sage: M
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]


        TESTS:

        Check that :issue:`13714` is fixed::

            sage: R.<x1,x2> = QQ[]
            sage: I = R.ideal(x2**2 + x1 - 2, x1**2 - 1)
            sage: R.one().lift(I)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not in the ideal
            sage: foo = I.complete_primary_decomposition() # indirect doctest
            sage: foo[0][0]
            Ideal (x1 + 1, x2^2 - 3) of Multivariate Polynomial Ring in x1, x2 over Rational Field
        """
        global errorreported

        cdef ideal *fI = idInit(1, 1)
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef int i = 0
        cdef int j
        cdef ring *r = self._parent_ring
        cdef ideal *res

        if isinstance(I, MPolynomialIdeal):
            I = I.gens()

        _I = idInit(len(I),1)

        for f in I:
            if not (isinstance(f,MPolynomial_libsingular)
                    and (<MPolynomial_libsingular>f)._parent is parent):
                try:
                    f = parent.coerce(f)
                except TypeError as msg:
                    id_Delete(&fI, r)
                    id_Delete(&_I, r)
                    raise TypeError(msg)

            _I.m[i] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
            i += 1

        fI.m[0] = p_Copy(self._poly, r)

        if r != currRing:
            rChangeCurrRing(r)  # idLift
        sig_on()
        res = idLift(_I, fI, NULL, 0, 0, 0)
        sig_off()
        if errorreported != 0:
            errorcode = errorreported
            errorreported = 0
            if errorcode == 1:
                raise ValueError("polynomial is not in the ideal")
            raise RuntimeError

        l = []
        for i from 0 <= i < IDELEMS(res):
            for j from 1 <= j <= IDELEMS(_I):
                l.append( new_MP(parent, pTakeOutComp(&res.m[i], 1)) )

        id_Delete(&fI, r)
        id_Delete(&_I, r)
        id_Delete(&res, r)
        return Sequence(l, check=False, immutable=True)

    def reduce(self, I):
        r"""
        Return a remainder of this polynomial modulo the
        polynomials in ``I``.

        INPUT:

        - ``I`` -- an ideal or a list/set/iterable of polynomials

        OUTPUT: a polynomial ``r``  such that:

        - ``self`` -- ``r`` is in the ideal generated by ``I``.

        - No term in ``r`` is divisible by any of the leading monomials
          of ``I``.

        The result ``r`` is canonical if:

        - ``I`` is an ideal, and Sage can compute a Groebner basis of it.

        - ``I`` is a list/set/iterable that is a (strong) Groebner basis
          for the term order of ``self``. (A strong Groebner basis is
          such that for every leading term ``t`` of the ideal generated
          by ``I``, there exists an element ``g`` of ``I`` such that the
          leading term of ``g`` divides ``t``.)

        The result ``r`` is implementation-dependent (and possibly
        order-dependent) otherwise. If ``I`` is an ideal and no Groebner
        basis can be computed, its list of generators ``I.gens()`` is
        used for the reduction.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f1 = -2 * x^2 + x^3
            sage: f2 = -2 * y + x* y
            sage: f3 = -x^2 + y^2
            sage: F = Ideal([f1,f2,f3])
            sage: g = x*y - 3*x*y^2
            sage: g.reduce(F)
            -6*y^2 + 2*y
            sage: g.reduce(F.gens())
            -6*y^2 + 2*y

        `\ZZ` is also supported. ::

            sage: P.<x,y,z> = ZZ[]
            sage: f1 = -2 * x^2 + x^3
            sage: f2 = -2 * y + x* y
            sage: f3 = -x^2 + y^2
            sage: F = Ideal([f1,f2,f3])
            sage: g = x*y - 3*x*y^2
            sage: g.reduce(F)
            -6*y^2 + 2*y
            sage: g.reduce(F.gens())
            -6*y^2 + 2*y

            sage: f = 3*x
            sage: f.reduce([2*x,y])
            x

        The reduction is not canonical when ``I`` is not a Groebner
        basis::

            sage: A.<x,y> = QQ[]
            sage: (x+y).reduce([x+y, x-y])
            2*y
            sage: (x+y).reduce([x-y, x+y])
            0
        """
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef int i = 0
        cdef ring *r = self._parent_ring
        cdef poly *res

        if r!=currRing: rChangeCurrRing(r)

        if isinstance(I, MPolynomialIdeal):
            try:
                strat = I._groebner_strategy()
                return strat.normal_form(self)
            except (TypeError, NotImplementedError):
                pass
            I = I.gens()

        _I = idInit(len(I),1)
        for f in I:
            if not (isinstance(f,MPolynomial_libsingular)
                    and (<MPolynomial_libsingular>f)._parent is parent):
                try:
                    f = parent.coerce(f)
                except TypeError as msg:
                    id_Delete(&_I,r)
                    raise TypeError(msg)

            _I.m[i] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
            i+=1

        #the second parameter would be qring!
        if r!=currRing: rChangeCurrRing(r)  # kNF
        res = kNF(_I, NULL, self._poly)
        id_Delete(&_I,r)
        return new_MP(parent,res)

    def divides(self, other):
        """
        Return ``True`` if this polynomial divides ``other``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: p = 3*x*y + 2*y*z + x*z
            sage: q = x + y + z + 1
            sage: r = p * q
            sage: p.divides(r)
            True
            sage: q.divides(p)
            False
            sage: r.divides(0)
            True
            sage: R.zero().divides(r)
            False
            sage: R.zero().divides(0)
            True
        """
        if self.is_zero():
            return other.is_zero()
        cdef ideal *_I
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *r = self._parent_ring
        cdef poly *rem

        if r != currRing:
            rChangeCurrRing(r)

        _I = idInit(1, 1)
        if not (isinstance(other,MPolynomial_libsingular)
                and (<MPolynomial_libsingular>other)._parent is parent):
            try:
                other = parent.coerce(other)
            except TypeError as msg:
                id_Delete(&_I,r)
                raise TypeError(msg)

        _I.m[0] = p_Copy(self._poly, r)

        if r != currRing:
            rChangeCurrRing(r)
        sig_on()
        rem = kNF(_I, NULL, (<MPolynomial_libsingular>other)._poly, 0, 1)
        sig_off()
        id_Delete(&_I, r)
        res = new_MP(parent, rem).is_zero()
        return res

    @coerce_binop
    def gcd(self, right, algorithm=None, **kwds):
        """
        Return the greatest common divisor of ``self`` and ``right``.

        INPUT:

        - ``right`` -- polynomial
        - ``algorithm``
          - ``'ezgcd'`` -- EZGCD algorithm
          - ``'modular'`` -- multi-modular algorithm (default)
        - ``**kwds`` -- ignored

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = (x*y*z)^6 - 1
            sage: g = (x*y*z)^4 - 1
            sage: f.gcd(g)
            x^2*y^2*z^2 - 1
            sage: GCD([x^3 - 3*x + 2, x^4 - 1, x^6 -1])
            x - 1

            sage: R.<x,y> = QQ[]
            sage: f = (x^3 + 2*y^2*x)^2
            sage: g = x^2*y^2
            sage: f.gcd(g)
            x^2

        We compute a gcd over a finite field::

            sage: # needs sage.rings.finite_rings
            sage: F.<u> = GF(31^2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3
            sage: gcd(p,q)  # yes, twice -- tests that singular ring is properly set.
            x^3 + (u + 1)*y^3 + z^3

        We compute a gcd over a number field::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: F.<u> = NumberField(x^3 - 2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3

        TESTS::

            sage: Q.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P(0).gcd(Q(0))
            0
            sage: x.gcd(1)
            1

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = R.change_ring(GF(3)).gen()
            sage: g = x + y
            sage: g.gcd(f)
            1
            sage: x.gcd(R.change_ring(GF(3)).gen())
            x

            sage: Pol.<x,y,z> = ZZ[]
            sage: p = x*y - 5*y^2 + x*z - z^2 + z
            sage: q = -3*x^2*y^7*z + 2*x*y^6*z^3 + 2*x^2*y^3*z^4 + x^2*y^5 - 7*x*y^5*z
            sage: (21^3*p^2*q).gcd(35^2*p*q^2) == -49*p*q
            True
        """
        cdef poly *_res
        cdef ring *_ring = self._parent_ring
        cdef MPolynomial_libsingular _right = <MPolynomial_libsingular>right

        if algorithm is None or algorithm == "modular":
            On(SW_USE_CHINREM_GCD)
            Off(SW_USE_EZGCD)
        elif algorithm == "ezgcd":
            Off(SW_USE_CHINREM_GCD)
            On(SW_USE_EZGCD)
        else:
            raise TypeError("algorithm %s not supported" % algorithm)

        if _right._poly == NULL:
            return self
        elif self._poly == NULL:
            return right
        elif p_IsOne(self._poly, _ring):
            return self
        elif p_IsOne(_right._poly, _ring):
            return right

        if _ring.cf.type != n_unknown:
            if _ring.cf.type == n_Znm or _ring.cf.type == n_Zn or _ring.cf.type == n_Z2m:
                raise NotImplementedError("GCD over rings not implemented.")

        if n_GetChar(_ring.cf) > 1<<29:
            raise NotImplementedError("GCD of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

        cdef int count = singular_polynomial_length_bounded(self._poly, 20) \
            + singular_polynomial_length_bounded(_right._poly,20)
        if count >= 20:
            sig_on()
        _res = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_right._poly, _ring), _ring )
        if count >= 20:
            sig_off()

        res = new_MP(self._parent, _res)
        return res

    @coerce_binop
    def lcm(self, MPolynomial_libsingular g):
        """
        Return the least common multiple of ``self`` and `g`.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: p = (x+y)*(y+z)
            sage: q = (z^4+2)*(y+z)
            sage: lcm(p,q)
            x*y*z^4 + y^2*z^4 + x*z^5 + y*z^5 + 2*x*y + 2*y^2 + 2*x*z + 2*y*z

            sage: P.<x,y,z> = ZZ[]
            sage: p = 2*(x+y)*(y+z)
            sage: q = 3*(z^4+2)*(y+z)
            sage: lcm(p,q)
            6*x*y*z^4 + 6*y^2*z^4 + 6*x*z^5 + 6*y*z^5 + 12*x*y + 12*y^2 + 12*x*z + 12*y*z

            sage: # needs sage.rings.finite_rings
            sage: r.<x,y> = PolynomialRing(GF(2**8, 'a'), 2)
            sage: a = r.base_ring().0
            sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
            sage: f.lcm(x^4)
            (a^2 + a)*x^6*y + (a^4 + a^3 + a)*x^4*y + (a^5)*x^4

            sage: # needs sage.rings.number_field
            sage: w = polygen(ZZ, 'w')
            sage: r.<x,y> = PolynomialRing(NumberField(w^4 + 1, 'a'), 2)
            sage: a = r.base_ring().0
            sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
            sage: f.lcm(x^4)
            (a^2 + a)*x^6*y + (a^3 + a - 1)*x^4*y + (-a)*x^4

        TESTS::

            sage: Pol.<x,y,z> = ZZ[]
            sage: p = -x*y + x*z + 54*x - 2
            sage: q = (5*p^2).lcm(3*p)
            sage: q * q.lc().sign() == 15*p^2
            True
            sage: lcm(2*x, 2*y)
            2*x*y
            sage: lcm(2*x, 2*x*y)
            2*x*y
        """
        cdef ring *_ring = self._parent_ring
        cdef poly *ret
        cdef poly *prod
        cdef poly *gcd
        cdef MPolynomial_libsingular _g
        if _ring!=currRing: rChangeCurrRing(_ring)

        if _ring.cf.type != n_unknown:
            if _ring.cf.type == n_Znm or _ring.cf.type == n_Zn or _ring.cf.type == n_Z2m:
                raise TypeError("LCM over non-integral domains not available.")

        if self._parent is not g._parent:
            _g = self._parent.coerce(g)
        else:
            _g = <MPolynomial_libsingular>g

        if n_GetChar(_ring.cf) > 1<<29:
            raise NotImplementedError("LCM of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

        cdef int count = singular_polynomial_length_bounded(self._poly, 20) \
            + singular_polynomial_length_bounded(_g._poly,20)
        if count >= 20:
            sig_on()
        if _ring != currRing:
            rChangeCurrRing(_ring)  # singclap_gcd
        gcd = singclap_gcd(p_Copy(self._poly, _ring), p_Copy(_g._poly, _ring), _ring )
        prod = pp_Mult_qq(self._poly, _g._poly, _ring)
        ret = p_Divide(prod, gcd, _ring)
        if count >= 20:
            sig_off()
        return new_MP(self._parent, ret)

    def is_squarefree(self):
        """
        Return ``True`` if this polynomial is square free.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: f= x^2 + 2*x*y + 1/2*z
            sage: f.is_squarefree()
            True
            sage: h = f^2
            sage: h.is_squarefree()
            False
        """
        # TODO:  Use Singular (4.x) intrinsics.  (Temporary solution from #17254.)
        return all(e == 1 for (f, e) in self.factor())

    @coerce_binop
    def quo_rem(self, MPolynomial_libsingular right):
        """
        Return quotient and remainder of ``self`` and ``right``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = y*x^2 + x + 1
            sage: f.quo_rem(x)
            (x*y + 1, 1)
            sage: f.quo_rem(y)
            (x^2, x + 1)

            sage: R.<x,y> = ZZ[]
            sage: f = 2*y*x^2 + x + 1
            sage: f.quo_rem(x)
            (2*x*y + 1, 1)
            sage: f.quo_rem(y)
            (2*x^2, x + 1)
            sage: f.quo_rem(3*x)
            (0, 2*x^2*y + x + 1)

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: R(0).quo_rem(R(1))
            (0, 0)
            sage: R(1).quo_rem(R(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        cdef poly *quo
        cdef poly *rem
        cdef MPolynomialRing_libsingular parent = self._parent
        cdef ring *r = self._parent_ring
        if r!=currRing: rChangeCurrRing(r)

        if right.is_zero():
            raise ZeroDivisionError

        if not self._parent._base.is_field():
            py_quo = self//right
            py_rem = self - right*py_quo
            return py_quo, py_rem

        if n_GetChar(r.cf) > 1<<29:
            raise NotImplementedError("Division of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

        cdef int count = singular_polynomial_length_bounded(self._poly, 15)
        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
            sig_on()
        if r!=currRing: rChangeCurrRing(r)   # singclap_pdivide
        quo = singclap_pdivide( self._poly, right._poly, r )
        rem = p_Add_q(p_Copy(self._poly, r), p_Neg(pp_Mult_qq(right._poly, quo, r), r), r)
        if count >= 15:
            sig_off()
        return new_MP(parent, quo), new_MP(parent, rem)

    def _singular_init_(self, singular=None):
        """
        Return a SINGULAR (as in the computer algebra system) string
        representation for this element.

        INPUT:

        - ``singular`` -- interpreter (default: ``sage.interfaces.singular.singular``)

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127), 3)
            sage: x._singular_init_()                                                   # needs sage.libs.singular
            'x'
            sage: (x^2+37*y+128)._singular_init_()                                      # needs sage.libs.singular
            'x2+37y+1'

        TESTS::

            sage: P(0)._singular_init_()                                                # needs sage.libs.singular
            '0'
        """
        if singular is None:
            from sage.interfaces.singular import singular
        self._parent._singular_().set_ring()
        return self._repr_short_()

    def sub_m_mul_q(self, MPolynomial_libsingular m, MPolynomial_libsingular q):
        """
        Return ``self - m*q``, where ``m`` must be a monomial and
        ``q`` a polynomial.

        INPUT:

        - ``m`` -- a monomial
        - ``q`` -- a polynomial

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: x.sub_m_mul_q(y, z)
            -y*z + x

        TESTS::

            sage: Q.<x,y,z> = PolynomialRing(QQ,3)
            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: P(0).sub_m_mul_q(P(0), P(1))
            0
            sage: x.sub_m_mul_q(Q.gen(1), Q.gen(2))
            -y*z + x
         """
        cdef ring *r = self._parent_ring

        if self._parent is not m._parent:
            m = self._parent.coerce(m)
        if self._parent is not q._parent:
            q = self._parent.coerce(q)

        if m._poly and m._poly.next:
            raise ArithmeticError("m must be a monomial.")
        elif not m._poly:
            return self

        cdef int le = p_GetMaxExp(m._poly, r)
        cdef int lr = p_GetMaxExp(q._poly, r)
        cdef int esum = le + lr

        overflow_check(esum, r)

        return new_MP(self._parent, p_Minus_mm_Mult_qq(p_Copy(self._poly, r), m._poly, q._poly, r))

    def _macaulay2_(self, macaulay2=None):
        """
        Return a Macaulay2 element corresponding to this polynomial.

        .. NOTE::

           Two identical rings are not canonically isomorphic in M2,
           so we require the user to explicitly set the ring, since
           there is no way to know if the ring has been set or not,
           and setting it twice screws everything up.

        EXAMPLES::

            sage: # optional - macaulay2
            sage: R.<x,y> = PolynomialRing(GF(7), 2)
            sage: f = (x^3 + 2*y^2*x)^7; f          # indirect doctest
            x^21 + 2*x^7*y^14
            sage: h = macaulay2(f); h
             21     7 14
            x   + 2x y
            sage: k = macaulay2(x+y); k
            x + y
            sage: k + h
             21     7 14
            x   + 2x y   + x + y
            sage: R(h)
            x^21 + 2*x^7*y^14
            sage: R(h^20) == f^20
            True

        TESTS:

        Check that constant polynomials are coerced to the polynomial ring, not
        the base ring (:issue:`28574`)::

            sage: R = QQ['x,y']
            sage: macaulay2(R('4')).ring()._operator('===', R)  # optional - macaulay2
            true
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2
        m2_parent = macaulay2(self.parent())
        macaulay2.use(m2_parent)
        return macaulay2('substitute(%s,%s)' % (repr(self), m2_parent._name))

    def add_m_mul_q(self, MPolynomial_libsingular m, MPolynomial_libsingular q):
        """
        Return ``self + m*q``, where ``m`` must be a monomial and
        ``q`` a polynomial.

        INPUT:

        - ``m`` -- a monomial
        - ``q`` -- a polynomial

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3)
            sage: x.add_m_mul_q(y,z)
            y*z + x

        TESTS::

            sage: R.<x,y,z>=PolynomialRing(QQ,3)
            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: P(0).add_m_mul_q(P(0),P(1))
            0
            sage: x.add_m_mul_q(R.gen(),R.gen(1))
            x*y + x
        """
        cdef ring *r = self._parent_ring

        if self._parent is not m._parent:
            m = self._parent.coerce(m)
        if self._parent is not q._parent:
            q = self._parent.coerce(q)

        if m._poly and m._poly.next:
            raise ArithmeticError("m must be a monomial.")
        elif not m._poly:
            return self

        cdef int le = p_GetMaxExp(m._poly, r)
        cdef int lr = p_GetMaxExp(q._poly, r)
        cdef int esum = le + lr

        overflow_check(esum, r)

        return new_MP(self._parent, p_Plus_mm_Mult_qq(p_Copy(self._poly, r), m._poly, q._poly, r))

    def __reduce__(self):
        """
        Serialize this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='degrevlex')
            sage: f = 27/113 * x^2 + y*z + 1/2
            sage: f == loads(dumps(f))
            True

            sage: P = PolynomialRing(GF(127), 3, names='abc')
            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1
            sage: f == loads(dumps(f))
            True

        TESTS:

        Verify that :issue:`9220` is fixed.

            sage: R=QQ['x']
            sage: S=QQ['x','y']
            sage: h=S.0^2
            sage: parent(h(R.0,0))
            Univariate Polynomial Ring in x over Rational Field
        """
        return unpickle_MPolynomial_libsingular, (self._parent,
                                                  self.monomial_coefficients())

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        INPUT:

        - ``codomain``
        - ``im_gens``

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: f = R.hom([y,x], R)
            sage: f(x^2 + 3*y^5)
            3*x^5 + y^2

            sage: R.<a,b,c,d> = QQ[]
            sage: S.<u> = QQ[]
            sage: h = R.hom([0,0,0,u], S) # indirect doctest
            sage: h((a+d)^3)
            u^3

        You can specify a map on the base ring::

            sage: # needs sage.rings.number_field
            sage: Zx.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: R.<x,y> = K[]
            sage: phi = R.hom([y,x], base_map=cc)
            sage: phi(x + i*y)
            (-i)*x + y
        """
        #TODO: very slow
        n = self.parent().ngens()
        if n == 0:
            return codomain.coerce(self)
        y = codomain(0)
        if base_map is None:
            # Just use conversion
            base_map = codomain
        for m, c in self.monomial_coefficients().items():
            y += base_map(c) * mul([im_gens[i] ** m[i] for i in range(n) if m[i]])
        return y

    def _derivative(self, MPolynomial_libsingular var):
        """
        Differentiates this polynomial with respect to the provided
        variable. This is completely symbolic so it is also defined
        over finite fields.

        INPUT:

        - ``variable`` -- the derivative is taken with respect to variable

        .. NOTE:: See also :meth:`derivative`

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f._derivative(x)
            9*x^2*y^2 + 3
            sage: f._derivative(y)
            6*x^3*y + 10*y

        The derivative is also defined over finite fields::

            sage: R.<x,y> = PolynomialRing(GF(2**8, 'a'), 2)                            # needs sage.rings.finite_rings
            sage: f = x^3*y^2 + y^2 + x + 2
            sage: f._derivative(x)                                                      # needs sage.rings.finite_rings
            x^2*y^2 + 1
        """
        if var is None:
            raise ValueError("you must specify which variable with respect to which to differentiate")

        cdef int i, var_i

        cdef poly *p
        if var._parent is not self._parent:
            raise TypeError("provided variable is not in same ring as self")
        cdef ring *_ring = self._parent_ring
        if _ring != currRing:
            rChangeCurrRing(_ring)

        var_i = -1
        for i from 0 <= i <= _ring.N:
            if p_GetExp(var._poly, i, _ring):
                if var_i == -1:
                    var_i = i
                else:
                    raise TypeError("provided variable is not univariate")

        if var_i == -1:
            raise TypeError("provided variable is constant")

        p = pDiff(self._poly, var_i)
        return new_MP(self._parent,p)

    def integral(self, MPolynomial_libsingular var):
        r"""
        Integrate this polynomial with respect to the provided variable.

        One requires that `\QQ` is contained in the ring.

        INPUT:

        - ``variable`` -- the integral is taken with respect to variable

        EXAMPLES::

            sage: R.<x, y> = PolynomialRing(QQ, 2)
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f.integral(x)
            3/4*x^4*y^2 + 5*x*y^2 + 3/2*x^2 + 2*x
            sage: f.integral(y)
            x^3*y^3 + 5/3*y^3 + 3*x*y + 2*y

        Check that :issue:`15896` is solved::

            sage: s = x+y
            sage: s.integral(x)+x
            1/2*x^2 + x*y + x
            sage: s.integral(x)*s
            1/2*x^3 + 3/2*x^2*y + x*y^2

        TESTS::

            sage: z, w = polygen(QQ, 'z, w')
            sage: f.integral(z)
            Traceback (most recent call last):
            ...
            TypeError: the variable is not in the same ring as self

            sage: f.integral(y**2)
            Traceback (most recent call last):
            ...
            TypeError: not a variable in the same ring as self

            sage: x,y = polygen(ZZ,'x,y')
            sage: y.integral(x)
            Traceback (most recent call last):
            ...
            TypeError: the ring must contain the rational numbers
        """
        cdef int index

        ambient_ring = var.parent()
        if ambient_ring is not self._parent:
            raise TypeError("the variable is not in the same ring as self")

        if not ambient_ring.has_coerce_map_from(QQ):
            raise TypeError("the ring must contain the rational numbers")

        gens = ambient_ring.gens()
        try:
            index = gens.index(var)
        except ValueError:
            raise TypeError("not a variable in the same ring as self")

        cdef poly *_p
        cdef poly *mon
        cdef ring *_ring = self._parent_ring
        if _ring != currRing:
            rChangeCurrRing(_ring)

        v = ETuple({index: 1}, len(gens))

        _p = p_ISet(0, _ring)
        for exp, coeff in self.monomial_coefficients().items():
            nexp = exp.eadd(v)  # new exponent
            mon = p_Init(_ring)
            p_SetCoeff(mon, sa2si(coeff / (1 + exp[index]), _ring), _ring)
            for pos in nexp.nonzero_positions():
                overflow_check(nexp[pos], _ring)
                p_SetExp(mon, pos + 1, nexp[pos], _ring)
            p_Setm(mon, _ring)
            _p = p_Add_q(_p, mon, _ring)
        return new_MP(self._parent, _p)

    def resultant(self, MPolynomial_libsingular other, variable=None):
        """
        Compute the resultant of this polynomial and the first
        argument with respect to the variable given as the second
        argument.

        If a second argument is not provide the first variable of
        the parent is chosen.

        INPUT:

        - ``other`` -- polynomial

        - ``variable`` -- optional variable (default: ``None``)

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: a = x+y
            sage: b = x^3-y^3
            sage: c = a.resultant(b); c
            -2*y^3
            sage: d = a.resultant(b,y); d
            2*x^3

        The SINGULAR example::

            sage: R.<x,y,z> = PolynomialRing(GF(32003), 3)                              # needs sage.rings.finite_rings
            sage: f = 3 * (x+2)^3 + y
            sage: g = x + y + z                                                         # needs sage.rings.finite_rings
            sage: f.resultant(g, x)                                                     # needs sage.rings.finite_rings
            3*y^3 + 9*y^2*z + 9*y*z^2 + 3*z^3 - 18*y^2 - 36*y*z - 18*z^2 + 35*y + 36*z - 24

        Resultants are also supported over the Integers::

            sage: R.<x,y,a,b,u> = PolynomialRing(ZZ, 5, order='lex')
            sage: r = (x^4*y^2 + x^2*y - y).resultant(x*y - y*a - x*b + a*b + u, x)
            sage: r
            y^6*a^4 - 4*y^5*a^4*b - 4*y^5*a^3*u + y^5*a^2 - y^5 + 6*y^4*a^4*b^2 + 12*y^4*a^3*b*u - 4*y^4*a^2*b + 6*y^4*a^2*u^2 - 2*y^4*a*u + 4*y^4*b - 4*y^3*a^4*b^3 - 12*y^3*a^3*b^2*u + 6*y^3*a^2*b^2 - 12*y^3*a^2*b*u^2 + 6*y^3*a*b*u - 4*y^3*a*u^3 - 6*y^3*b^2 + y^3*u^2 + y^2*a^4*b^4 + 4*y^2*a^3*b^3*u - 4*y^2*a^2*b^3 + 6*y^2*a^2*b^2*u^2 - 6*y^2*a*b^2*u + 4*y^2*a*b*u^3 + 4*y^2*b^3 - 2*y^2*b*u^2 + y^2*u^4 + y*a^2*b^4 + 2*y*a*b^3*u - y*b^4 + y*b^2*u^2

        TESTS::

            sage: P.<x,y> = PolynomialRing(QQ, order='degrevlex')
            sage: a = x+y
            sage: b = x^3-y^3
            sage: c = a.resultant(b); c
            -2*y^3
            sage: d = a.resultant(b,y); d
            2*x^3

            sage: P.<x,y> = PolynomialRing(ZZ,2)
            sage: f = x+y
            sage: g=y^2+x
            sage: f.resultant(g,y)
            x^2 + x
        """
        cdef ring *_ring = self._parent_ring
        cdef poly *rt

        if variable is None:
            variable = self.parent().gen(0)

        if self._parent is not other._parent:
            raise TypeError("first parameter needs to be an element of self.parent()")

        if not variable.parent() is self.parent():
            raise TypeError("second parameter needs to be an element of self.parent() or None")

        if n_GetChar(_ring.cf) > 1<<29:
            raise NotImplementedError("Resultants of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")

        if isinstance(self._parent._base, IntegerRing_class):
            ret = self.change_ring(QQ).resultant(other.change_ring(QQ),
                                                 variable.change_ring(QQ))
            return ret.change_ring(ZZ)
        elif not self._parent._base.is_field():
            raise ValueError("Resultants require base fields or integer base ring.")

        cdef int count = singular_polynomial_length_bounded(self._poly, 20) \
            + singular_polynomial_length_bounded(other._poly,20)
        if count >= 20:
            sig_on()
        if _ring != currRing: rChangeCurrRing(_ring)   # singclap_resultant
        rt =  singclap_resultant(p_Copy(self._poly, _ring),
                                 p_Copy(other._poly, _ring),
                                 p_Copy((<MPolynomial_libsingular>variable)._poly, _ring ),
                                 _ring)
        if count >= 20:
            sig_off()
        return new_MP(self._parent, rt)

    def coefficients(self):
        """
        Return the nonzero coefficients of this polynomial in a list.
        The returned list is decreasingly ordered by the term ordering
        of the parent.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]

            sage: R.<x,y,z> = PolynomialRing(QQ, order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        AUTHOR:

        - Didier Deshommes
        """
        cdef poly *p
        cdef ring *r = self._parent_ring
        if r!=currRing: rChangeCurrRing(r)
        base = self._parent._base
        p = self._poly
        coeffs = list()
        while p:
            coeffs.append(si2sa(p_GetCoeff(p, r), r, base))
            p = pNext(p)
        return coeffs

    def global_height(self, prec=None):
        """
        Return the (projective) global height of the polynomial.

        This returns the absolute logarithmic height of the coefficients
        thought of as a projective point.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default :class:`RealField` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 3*x^3 + 2*x*y^2
            sage: exp(f.global_height())                                                # needs sage.symbolic
            3.00000000000000

        ::

            sage: # needs sage.rings.number_field
            sage: K.<k> = CyclotomicField(3)
            sage: R.<x,y> = PolynomialRing(K, sparse=True)
            sage: f = k*x*y + 1
            sage: exp(f.global_height())
            1.00000000000000

        Scaling should not change the result::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 1/25*x^2 + 25/3*x*y + y^2
            sage: f.global_height()                                                     # needs sage.symbolic
            6.43775164973640
            sage: g = 100 * f
            sage: g.global_height()                                                     # needs sage.symbolic
            6.43775164973640

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<k> = NumberField(x^2 + 5)
            sage: T.<t,w> = PolynomialRing(K)
            sage: f = 1/1331 * t^2 + 5 * w + 7
            sage: f.global_height()                                                     # needs sage.symbolic
            9.13959596745043

        ::

            sage: R.<x,y> = QQ[]
            sage: f = 1/123*x*y + 12
            sage: f.global_height(prec=2)                                               # needs sage.symbolic
            8.0

        ::

            sage: R.<x,y> = QQ[]
            sage: f = 0*x*y
            sage: f.global_height()
            0.000000000000000
        """
        if prec is None:
            prec = 53

        if self.is_zero():
            return RealField(prec).zero()

        K = self.base_ring()
        if K in NumberFields() or isinstance(K, sage.rings.abc.Order) or K == ZZ:
            f = self
        else:
            raise TypeError("Must be over a Numberfield or a Numberfield Order.")

        from sage.schemes.projective.projective_space import ProjectiveSpace
        P = ProjectiveSpace(K, f.number_of_terms()-1)
        return P.point(f.coefficients()).global_height(prec=prec)

    def local_height(self, v, prec=None):
        """
        Return the maximum of the local height of the coefficients of
        this polynomial.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring

        - ``prec`` -- desired floating point precision (default:
          default :class:`RealField` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 1/1331*x^2 + 1/4000*y^2
            sage: f.local_height(1331)
            7.19368581839511

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<k> = NumberField(x^2 - 5)
            sage: T.<t,w> = K[]
            sage: I = K.ideal(3)
            sage: f = 1/3*t*w + 3
            sage: f.local_height(I)
            1.09861228866811

        ::

            sage: R.<x,y> = QQ[]
            sage: f = 1/2*x*y + 2
            sage: f.local_height(2, prec=2)
            0.75
        """
        if prec is None:
            prec = 53

        K = FractionField(self.base_ring())
        if K not in NumberFields() and not isinstance(K, sage.rings.abc.Order) and K != ZZ:
            raise TypeError("must be over a Numberfield or a Numberfield order")

        return max([K(c).local_height(v, prec=prec) for c in self.coefficients()])

    def local_height_arch(self, i, prec=None):
        """
        Return the maximum of the local height at the ``i``-th infinite place
        of the coefficients of this polynomial.

        INPUT:

        - ``i`` -- integer

        - ``prec`` -- desired floating point precision (default:
          default :class:`RealField` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = 210*x*y
            sage: f.local_height_arch(0)
            5.34710753071747

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<k> = NumberField(x^2 - 5)
            sage: T.<t,w> = K[]
            sage: f = 1/2*t*w + 3
            sage: f.local_height_arch(1, prec=52)
            1.09861228866811

        ::

            sage: R.<x,y> = QQ[]
            sage: f = 1/2*x*y + 3
            sage: f.local_height_arch(0, prec=2)
            1.0
        """
        K = FractionField(self.base_ring())
        if K not in NumberFields() and not isinstance(K, sage.rings.abc.Order) and K != ZZ:
            return TypeError("must be over a Numberfield or a Numberfield Order")

        if K == QQ:
            return max([K(c).local_height_arch(prec=prec) for c in self.coefficients()])
        return max([K(c).local_height_arch(i, prec=prec) for c in self.coefficients()])

    def gradient(self):
        """
        Return a list of partial derivatives of this polynomial,
        ordered by the variables of the parent.

        EXAMPLES::

           sage: P.<x,y,z> = PolynomialRing(QQ,3)
           sage: f= x*y + 1
           sage: f.gradient()
           [y, x, 0]
        """
        cdef ring *r = self._parent_ring
        cdef int k

        if r!=currRing: rChangeCurrRing(r)
        i = []
        for k from 0 < k <= r.N:
            i.append( new_MP(self._parent, pDiff(self._poly, k)))
        return i

    def numerator(self):
        """
        Return a numerator of ``self`` computed as ``self * self.denominator()``.

        If the base_field of ``self`` is the Rational Field then the
        numerator is a polynomial whose base_ring is the Integer Ring,
        this is done for compatibility to the univariate case.

        .. WARNING::

            This is not the numerator of the rational function defined by
            ``self``, which would always be ``self`` since ``self`` is a
            polynomial.

        EXAMPLES:

        First we compute the numerator of a polynomial with
        integer coefficients, which is of course ``self``.

        ::

            sage: R.<x, y> = ZZ[]
            sage: f = x^3 + 17*y + 1
            sage: f.numerator()
            x^3 + 17*y + 1
            sage: f == f.numerator()
            True

        Next we compute the numerator of a polynomial with rational
        coefficients.

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = (1/17)*x^19 - (2/3)*y + 1/3; f
            1/17*x^19 - 2/3*y + 1/3
            sage: f.numerator()
            3*x^19 - 34*y + 17
            sage: f == f.numerator()
            False
            sage: f.numerator().base_ring()
            Integer Ring

        We check that the computation of numerator and denominator
        is valid.

        ::

            sage: K=QQ['x,y']
            sage: f=K.random_element()
            sage: f.numerator() / f.denominator() == f
            True

        The following tests against a bug fixed in :issue:`11780`::

            sage: P.<foo,bar> = ZZ[]
            sage: Q.<foo,bar> = QQ[]
            sage: f = Q.random_element()
            sage: f.numerator().parent() is P
            True
        """
        if self.base_ring() is QQ:
            # This part is for compatibility with the univariate case,
            # where the numerator of a polynomial over RationalField
            # is a polynomial over IntegerRing
            #
            # Github issue #11780: Create the polynomial ring over
            # the integers using the (cached) polynomial ring constructor:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            integer_poly_ring = PolynomialRing(ZZ,
                                               self.parent().ngens(),
                                               self.parent().gens(),
                                               order=self.parent().term_order())
            return integer_poly_ring(self * self.denominator())
        else:
            return self * self.denominator()

    def in_subalgebra(self, J, algorithm='algebra_containment'):
        """
        Return whether this polynomial is contained in the subalgebra
        generated by ``J``

        INPUT:

        - ``J`` -- list of elements of the parent polynomial ring

        - ``algorithm`` -- can be ``'algebra_containment'`` (the default),
          ``'inSubring'``, or ``'groebner'``

          - ``'algebra_containment'``: use Singular's
            ``algebra_containment`` function,
            https://www.singular.uni-kl.de/Manual/4-2-1/sing_1247.htm#SEC1328. The
            Singular documentation suggests that this is frequently
            faster than the next option.

          - ``'inSubring'``: use Singular's ``inSubring`` function,
            https://www.singular.uni-kl.de/Manual/4-2-0/sing_1240.htm#SEC1321.

          - ``'groebner'``: use the algorithm described in Singular's
            documentation, but within Sage: if the subalgebra
            generators are `y_1`, ..., `y_m`, then create a new
            polynomial algebra with the old generators along with new
            ones: `z_1`, ..., `z_m`. Create the ideal `(z_1 - y_1,
            ..., z_m - y_m)`, and reduce the polynomial modulo this
            ideal. The polynomial is contained in the subalgebra if
            and only if the remainder involves only the new variables
            `z_i`.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: J = [x^2 + y^2, x^2 + z^2]
            sage: (y^2).in_subalgebra(J)
            False
            sage: a = (x^2 + y^2) * (x^2 + z^2)
            sage: a.in_subalgebra(J, algorithm='inSubring')
            True
            sage: (a^2).in_subalgebra(J, algorithm='groebner')
            True
            sage: (a + a^2).in_subalgebra(J)
            True
        """
        R = self.parent()
        algorithm = algorithm.lower()
        from sage.libs.singular.function import (singular_function,
                                                 get_printlevel,
                                                 set_printlevel,
                                                 lib as singular_lib)
        singular_lib('algebra.lib')
        if algorithm == "algebra_containment":
            # verbose unless printlevel is -1.
            saved_printlevel = get_printlevel()
            set_printlevel(-1)
            contains = singular_function('algebra_containment')(self, R.ideal(J)) == 1
            set_printlevel(saved_printlevel)
            return contains
        elif algorithm == "insubring":
            return singular_function('inSubring')(self, R.ideal(J))[0] == 1
        elif algorithm == "groebner":
            new_gens = [f"newgens{i}" for i in range(len(J))]
            ngens = len(new_gens)
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            S = PolynomialRing(R.base_ring(), R.gens() + tuple(new_gens),
                               order=f'degrevlex({len(R.gens())}),degrevlex({ngens})')
            new_gens = S.gens()[-ngens:]
            I = S.ideal([g - S(j) for (g,j) in zip(new_gens, J)])
            z = S(self).reduce(I)
            return set(z.variables()).issubset(new_gens)
        else:
            raise ValueError("unknown algorithm")


def unpickle_MPolynomial_libsingular(MPolynomialRing_libsingular R, d):
    """
    Deserialize an ``MPolynomial_libsingular`` object.

    INPUT:

    - ``R`` -- the base ring
    - ``d`` -- a Python dictionary as returned by :meth:`MPolynomial_libsingular.dict`

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(x)) == x # indirect doctest
        True
    """
    cdef ring *r = R._ring
    cdef poly *m
    cdef poly *p
    cdef int _i, _e
    cdef int ln
    cdef sBucket *bucket

    rChangeCurrRing(r)
    bucket = sBucketCreate(r)
    try:
        for mon, c in d.items():
            m = p_Init(r)
            for i, e in mon.sparse_iter():
                _i = i
                if _i >= r.N:
                    p_Delete(&m, r)
                    raise TypeError("variable index too big")
                _e = e
                if _e <= 0:
                    p_Delete(&m, r)
                    raise TypeError("exponent too small")
                overflow_check(_e, r)
                p_SetExp(m, _i+1, _e, r)
            p_SetCoeff(m, sa2si(c, r), r)
            p_Setm(m, r)
            sBucket_Merge_m(bucket, m)
        ln=0
        sBucketClearMerge(bucket, &p, &ln)
        sBucketDestroy(&bucket)
    except Exception:
        sBucketDeleteAndDestroy(&bucket)
        raise
    return new_MP(R, p)


cdef inline poly *addwithcarry(poly *tempvector, poly *maxvector, int pos, ring *_ring) noexcept:
    if p_GetExp(tempvector, pos, _ring) < p_GetExp(maxvector, pos, _ring):
        p_SetExp(tempvector, pos, p_GetExp(tempvector, pos, _ring)+1, _ring)
    else:
        p_SetExp(tempvector, pos, 0, _ring)
        tempvector = addwithcarry(tempvector, maxvector, pos + 1, _ring)
    p_Setm(tempvector, _ring)
    return tempvector


cdef inline MPolynomial_libsingular new_MP(MPolynomialRing_libsingular parent, poly *juice):
    """
    Construct MPolynomial_libsingular from parent and SINGULAR poly.

    INPUT:

    - ``parent`` -- a :class:`MPolynomialRing_libsingular`
      instance; the parent of the polynomial to create

    - ``juice`` -- a pointer to a Singular ``poly`` C struct

    OUTPUT: a Python object :class:`MPolynomial_libsingular`

    The ownership of ``juice`` will be transferred to the Python
    object. You must not free it yourself. Singular will modify the
    polynomial, so it is your responsibility to make a copy if the
    Singular data structure is used elsewhere.
    """
    cdef MPolynomial_libsingular p = MPolynomial_libsingular.__new__(MPolynomial_libsingular)
    p._parent = parent
    p._parent_ring = singular_ring_reference(parent._ring)
    p._poly = juice
    p_Normalize(p._poly, p._parent_ring)
    return p


cdef poly *MPolynomial_libsingular_get_element(object self) noexcept:
    return (<MPolynomial_libsingular>self)._poly
