"""
libSingular: Conversion Routines and Initialisation

AUTHOR:

- Martin Albrecht <malb@informatik.uni-bremen.de>

- Miguel Marco <mmarco@unizar.es> (2021): added transcendental extensions over Q
"""

# ****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

include "sage/libs/ntl/decl.pxi"

cdef extern from "limits.h":
    long INT_MAX
    long INT_MIN

import os
from warnings import warn

from libc.stdint cimport int64_t
from sage.libs.singular.decl cimport *

from sage.rings.polynomial.polydict import ETuple
from sage.libs.singular.function cimport new_sage_polynomial, access_singular_ring

from sage.rings.rational_field import RationalField
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.polynomial.polynomial_ring import PolynomialRing_field
from sage.rings.fraction_field import FractionField_generic

from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_rings.finite_field_ntl_gf2e import FiniteField_ntl_gf2e
from sage.libs.gmp.all cimport *

from sage.cpython.string import FS_ENCODING
from sage.cpython.string cimport str_to_bytes, char_to_str, bytes_to_str

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

ctypedef struct fraction "fractionObject":
    poly *numerator
    poly *denominator
    int complexity

_saved_options = (int(0),0,0)

cdef Rational si2sa_QQ(number *n, number **nn, ring *_ring):
    """
    Create a sage rational number from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular rational number

    - ``*n`` -- a pointer to a pointer like before

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    OUTPUT: a sage Rational

    TESTS::

        sage: P.<x,y,z> = QQ[]
        sage: P(1/3).lc()
        1/3
        sage: P(1).lc()
        1
        sage: P(0).lc()
        0
        sage: P(-1/3).lc()
        -1/3
        sage: type(P(3).lc())
        <class 'sage.rings.rational.Rational'>
    """
    cdef number *nom
    cdef number *denom
    cdef mpq_t _z

    cdef mpz_t nom_z, denom_z

    cdef Rational z

    mpq_init(_z)

    ##  Immediate integers handles carry the tag 'SR_INT', i.e. the last bit is 1.
    ##  This distinguishes immediate integers from other handles which point to
    ##  structures aligned on 4 byte boundaries and therefore have last bit zero.
    ##  (The second bit is reserved as tag to allow extensions of this scheme.)
    ##  Using immediates as pointers and dereferencing them gives address errors.
    nom = nlGetNumerator(n, _ring.cf)
    mpz_init(nom_z)

    if (SR_HDL(nom) & SR_INT): mpz_set_si(nom_z, SR_TO_INT(nom))
    else: mpz_set(nom_z,nom.z)

    mpq_set_num(_z,nom_z)
    nlDelete(&nom,_ring.cf)
    mpz_clear(nom_z)

    denom = nlGetDenom(n, _ring.cf)
    mpz_init(denom_z)

    if (SR_HDL(denom) & SR_INT): mpz_set_si(denom_z, SR_TO_INT(denom))
    else: mpz_set(denom_z,denom.z)

    mpq_set_den(_z, denom_z)
    nlDelete(&denom,_ring.cf)
    mpz_clear(denom_z)

    nn[0] = n
    z = Rational()
    z.set_from_mpq(_z)
    mpq_clear(_z)
    return z

cdef Integer si2sa_ZZ(number *n, ring *_ring):
    """
    Create a sage integer number from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular integer number

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    OUTPUT: a sage Integer


    TESTS::

        sage: P.<x,y,z> = ZZ[]
        sage: P(3).lc()
        3
        sage: P(0).lc()
        0
        sage: P(-3).lc()
        -3
        sage: P(-1234567890).lc()
        -1234567890
        sage: type(P(3).lc())
        <class 'sage.rings.integer.Integer'>
    """
    cdef Integer z
    z = Integer()
    z.set_from_mpz(<mpz_ptr>n)
    return z

cdef FFgivE si2sa_GFqGivaro(number *n, ring *_ring, Cache_givaro cache):
    """
    Create a sage element of a small finite field from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in a finite field

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``cache`` -- a Givaro number field

    OUTPUT: a sage element of ``cache``

    TESTS::

        sage: K.<a> = GF(5^3)
        sage: R.<x,y,z> = PolynomialRing(K)
        sage: K( (4*R(a)^2 + R(a))^3 )
        a^2
        sage: K(R(0))
        0
    """
    cdef poly *z
    cdef int c, e
    cdef int a
    cdef int ret
    cdef int order
    cdef ring *cfRing = _ring.cf.extRing

    if _ring.cf.cfIsZero(n,_ring.cf):
        return cache._zero_element
    elif _ring.cf.cfIsOne(n,_ring.cf):
        return cache._one_element

    z = <poly*>n

    a = cache.objectptr.indeterminate()
    ret = cache.objectptr.zero
    order = cache.objectptr.cardinality() - 1

    while z:
        c = cache.objectptr.initi(c, <int64_t>p_GetCoeff(z, cfRing))
        e = p_GetExp(z, 1, cfRing)
        if e == 0:
            ret = cache.objectptr.add(ret, c, ret)
        else:
            a = ( e * cache.objectptr.indeterminate() ) % order
            ret = cache.objectptr.axpy(ret, c, a, ret)
        z = <poly*>pNext(<poly*>z)
    return (<FFgivE>cache._zero_element)._new_c(ret)

cdef FFgf2eE si2sa_GFqNTLGF2E(number *n, ring *_ring, Cache_ntl_gf2e cache):
    """
    Create a sage element of a finite field of characteristic 2 from a
    singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in a finite field

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``cache`` -- a ntl_gf2e number field

    OUTPUT: a sage element of ``cache``


    TESTS::

        sage: K.<a> = GF(2^20)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        a^11 + a^10 + a^8 + a^7 + a^6 + a^5 + a^2 + a
        sage: type(f.lc())
        <class 'sage.rings.finite_rings.element_ntl_gf2e.FiniteField_ntl_gf2eElement'>
    """
    cdef poly *z
    cdef long c
    cdef int e
    cdef FFgf2eE a
    cdef FFgf2eE ret
    cdef ring *cfRing = _ring.cf.extRing

    if _ring.cf.cfIsZero(n,_ring.cf):
        return cache._zero_element
    elif _ring.cf.cfIsOne(n,_ring.cf):
        return cache._one_element

    z = <poly*>n
    a = cache._gen
    ret = cache._zero_element

    while z:
        c = <long>p_GetCoeff(z, cfRing)
        e = p_GetExp(z, 1, cfRing)
        ret += c * a**e
        z = <poly*>pNext(<poly*>z)
    return ret

cdef object si2sa_GFq_generic(number *n, ring *_ring, object base):
    """
    Create a sage element of a generic finite field from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in a finite field

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``base`` -- a sage finite field

    OUTPUT: a sage element of ``base``

    TESTS::

        sage: K.<a> = GF(3^16)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        a^12 + a^11 + a^9 + a^8 + a^7 + 2*a^6 + a^5
        sage: type(f.lc())
        <class 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt'>

    Try the largest characteristic which Singular supports::

        sage: p = previous_prime(2^31)
        sage: F.<a> = FiniteField(p^2)
        sage: R.<x,y> = F[]
        sage: R(-1).constant_coefficient()  # indirect doctest
        2147483646
    """
    cdef poly *z
    cdef long c
    cdef int e
    cdef object a
    cdef object ret
    cdef ring *cfRing = _ring.cf.extRing

    if _ring.cf.type in (n_Zn, n_Znm):
        return si2sa_ZZmod(n, _ring, base)

    if _ring.cf.cfIsZero(n,_ring.cf):
        return base.zero()
    elif _ring.cf.cfIsOne(n,_ring.cf):
        return base.one()

    z = <poly*>n

    a = base.gen()
    ret = base.zero()

    while z:
        c = <long>p_GetCoeff(z, cfRing)
        e = p_GetExp(z, 1, cfRing)
        if e == 0:
            ret = ret + c
        elif c != 0:
            ret = ret  + c * a**e
        z = <poly*>pNext(<poly*>z)
    return ret

cdef object si2sa_transext_QQ(number *n, ring *_ring, object base):
    """
    Create a sage element of a transcendental extension of ``QQ`` from a
    singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in a transcendental extension
        of the rationals

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``base`` -- a sage FractionField

    OUTPUT: a sage element of ``base``

    TESTS::

        sage: F = PolynomialRing(QQ,'a,b').fraction_field()
        sage: F.inject_variables()
        Defining a, b
        sage: R.<x,y> = F[]
        sage: a*x
        a*x
        sage: I = R.ideal([a*x])
        sage: I
        Ideal (a*x) of Multivariate Polynomial Ring in x, y over Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
        sage: I.groebner_basis()
        [x]
        sage: I = R.ideal([a*x+b*y^2, (b+a)/(b-a)*x^3-3*y*x])
        sage: I.groebner_basis()
        [x^3 + (3*a - 3*b)/(a + b)*x*y, y^2 + a/b*x]
        sage: R.term_order()
        Degree reverse lexicographic term order
    """

    cdef poly *numer
    cdef poly *denom
    cdef number *c
    cdef int e
    cdef fraction *frac
    cdef object snumer
    cdef object sdenom

    cdef ring *cfRing = _ring.cf.extRing

    if _ring.cf.cfIsZero(n,_ring.cf):
        return base._zero_element
    elif _ring.cf.cfIsOne(n,_ring.cf):
        return base._one_element

    snumer = base(0)
    sdenom = base(0)

    frac = <fraction*>n

    numer = frac.numerator
    denom = frac.denominator

    while numer:
        c = p_GetCoeff(numer, cfRing)
        coeff = si2sa_QQ(c, &c, cfRing)
        numer.coef = c
        for i in range(base.ngens()):
            e = p_GetExp(numer, i+1, cfRing)
            if e!= 0:
                coeff *= base.gen(i)**e
        snumer += coeff
        numer = <poly*>pNext(<poly*>numer)

    if not denom:
        sdenom = base(1)
    else:
        while denom:
            c = p_GetCoeff(denom, cfRing)
            coeff = si2sa_QQ(c, &c, cfRing)
            denom.coef = c
            for i in range(base.ngens()):
                e = p_GetExp(denom, i+1, cfRing)
                if e!= 0:
                    coeff *= base.gen(i)**e
            sdenom += coeff
            denom = <poly*>pNext(<poly*>denom)

    return snumer/sdenom

cdef object si2sa_transext_FF(number *n, ring *_ring, object base):
    """
    Create a sage element of a transcendental extension of a prime field from a
    singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in a transcendental extension
      of the rationals

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``base`` -- a sage FractionField

    OUTPUT: a sage element of ``base``

    TESTS::

        sage: F = PolynomialRing(FiniteField(7),'a,b').fraction_field()
        sage: R.<x,y,z> = F[]
        sage: n = R(5)
        sage: n
        -2
        sage: type(n)
        <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
    """

    cdef poly *numer
    cdef poly *denom
    cdef number *c
    cdef int e
    cdef fraction *frac
    cdef object snumer
    cdef object sdenom

    cdef ring *cfRing = _ring.cf.extRing

    if _ring.cf.cfIsZero(n,_ring.cf):
        return base._zero_element
    elif _ring.cf.cfIsOne(n,_ring.cf):
        return base._one_element

    snumer = base(0)
    sdenom = base(0)

    frac = <fraction*>n

    numer = frac.numerator
    denom = frac.denominator

    while numer:

        c = p_GetCoeff(numer, cfRing)
        coeff = base(cfRing.cf.cfInt(c, cfRing.cf))
        numer.coef = c
        for i in range(base.ngens()):
            e = p_GetExp(numer, i+1, cfRing)
            if e!= 0:
                coeff *= base.gen(i)**e
        snumer += coeff
        numer = <poly*>pNext(<poly*>numer)

    if not denom:
        sdenom = base(1)
    else:
        while denom:
            c = p_GetCoeff(denom, cfRing)
            coeff = base(cfRing.cf.cfInt(c, cfRing.cf))
            denom.coef = c
            for i in range(base.ngens()):
                e = p_GetExp(denom, i+1, cfRing)
                if e!= 0:
                    coeff *= base.gen(i)**e
            sdenom += coeff
            denom = <poly*>pNext(<poly*>denom)

    return snumer/sdenom

cdef object si2sa_NF(number *n, ring *_ring, object base):
    """
    Create a sage element of a number field from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in an algebraic extension of
      the rationals

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``base`` -- a sage NumberField

    OUTPUT: a sage element of ``base``


    TESTS::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 2)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        1024*a
        sage: type(f.lc())
        <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic_sqrt'>
    """
    cdef poly *z
    cdef number *c
    cdef int e
    cdef object a
    cdef object ret
    cdef ring *cfRing = _ring.cf.extRing

    if _ring.cf.cfIsZero(n,_ring.cf):
        return base._zero_element
    elif _ring.cf.cfIsOne(n,_ring.cf):
        return base._one_element

    z = <poly*>n

    a = base.gen()
    ret = base(0)

    while z:
        # p_GetCoeff returns a reference
        c = p_GetCoeff(z, cfRing)
        # si2sa_QQ might modify c
        coeff = si2sa_QQ(c, &c, cfRing)
        # so we force it back.
        z.coef = c
        #pSetCoeff0(z,c)
        #p_SetCoeff(z, c, cfRing)
        # rather than trying to let Cython and C++ automagically modify it
        #coeff = si2sa_QQ(p_GetCoeff(z, cfRing), cfRing)
        e = p_GetExp(z, 1, cfRing)
        if e == 0:
            ret = ret + coeff
        elif coeff != 0:
            ret = ret + coeff * a**e
        z = <poly*>pNext(<poly*>z)
    return base(ret)

cdef inline object si2sa_ZZmod(number *n, ring *_ring, object base):
    """
    Create a sage element of a ring of integers modulo n from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number in a ring of integers modulo n

    - ``_ ring`` -- a (pointer to) a singular ring, in whose coefficient field
      lives ``n``

    - ``base`` -- a sage IntegerModRing

    OUTPUT: a sage element of ``base``

    TESTS::

        sage: P.<x,y,z> = Integers(10)[]
        sage: P(3).lc()
        3
        sage: P(13).lc()
        3

        sage: P.<x,y,z> = Integers(16)[]
        sage: P(3).lc()
        3
        sage: P(19).lc()
        3

        sage: P.<x,y,z> = Integers(3**2)[]
        sage: P(3).lc()
        3
        sage: P(12).lc()
        3

        sage: P.<x,y,z> = Integers(2^32)[]
        sage: P(2^32-1).lc()
        4294967295

        sage: P(3).lc()
        3

        sage: P.<x,y,z> = Integers(17^20)[]
        sage: P(17^19 + 3).lc()
        239072435685151324847156

        sage: P(3)
        3
    """
    cdef Integer ret
    if _ring.cf.type == n_Z2m:
        return base(<long>n)
    elif _ring.cf.type == n_Znm or _ring.cf.type == n_Zn:
        ret = Integer()
        ret.set_from_mpz(<mpz_ptr>n)
        return base(ret)

    return base(_ring.cf.cfInt(n,_ring.cf))


cdef list singular_monomial_exponents(poly *p, ring *r):
    r"""
    Return the list of exponents of monomial ``p``.
    """
    cdef int v
    cdef list ml = [None] * r.N

    for v in range(1, r.N + 1):
        ml[v-1] = p_GetExp(p, v, r)
    return ml

cpdef list si2sa_resolution(Resolution res):
    r"""
    Pull the data from Singular resolution ``res`` to construct a Sage
    resolution.

    INPUT:

    - ``res`` -- Singular resolution

    The procedure is destructive and ``res`` is not usable afterward.

    EXAMPLES::

        sage: from sage.libs.singular.singular import si2sa_resolution
        sage: from sage.libs.singular.function import singular_function
        sage: module = singular_function("module")
        sage: mres = singular_function('mres')

        sage: S.<x,y,z,w> = PolynomialRing(QQ)
        sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
        sage: mod = module(I)
        sage: r = mres(mod, 0)
        sage: si2sa_resolution(r)
        [
                                         [ y  x]
                                         [-z -y]
        [z^2 - y*w y*z - x*w y^2 - x*z], [ w  z]
        ]
    """
    cdef ring *singular_ring
    cdef syStrategy singular_res
    cdef poly *p
    cdef poly *p_iter
    cdef poly *first
    cdef poly *previous
    cdef poly *acc
    cdef resolvente mods
    cdef ideal *mod
    cdef int i, j, k, idx, rank, nrows, ncols
    cdef bint zero_mat
    cdef list degs, matdegs

    from sage.modules.free_module import FreeModule
    from sage.matrix.constructor import matrix as _matrix

    singular_res = res._resolution[0]
    sage_ring = res.base_ring
    singular_ring = access_singular_ring(res.base_ring)

    if singular_res.minres != NULL:
        mods = singular_res.minres
    elif singular_res.fullres != NULL:
        mods = singular_res.fullres
    else:
        raise ValueError('Singular resolution is not usable')

    cdef list res_mats = []

    # length is the length of fullres. The length of minres
    # can be shorter. Hence we avoid SEGFAULT by stopping
    # at NULL pointer.
    for idx in range(singular_res.length):
        mod = <ideal *> mods[idx]
        if mod == NULL:
            break
        rank = mod.rank
        free_module = FreeModule(sage_ring, rank)

        nrows = rank
        ncols = mod.ncols # IDELEMS(mod)

        mat = _matrix(sage_ring, nrows, ncols)
        matdegs = []
        zero_mat = True
        for j in range(ncols):
            p = <poly *> mod.m[j]
            degs = []
            # code below copied and modified from to_sage_vector_destructive
            # in sage.libs.singular.function.Converter
            for i in range(1, rank + 1):
                previous = NULL
                acc = NULL
                first = NULL
                p_iter = p
                while p_iter != NULL:
                    if p_GetComp(p_iter, singular_ring) == i:
                        p_SetComp(p_iter, 0, singular_ring)
                        p_Setm(p_iter, singular_ring)
                        if acc == NULL:
                            first = p_iter
                        else:
                            acc.next = p_iter
                        acc = p_iter
                        if p_iter == p:
                            p = pNext(p_iter)
                        if previous != NULL:
                            previous.next = pNext(p_iter)
                        p_iter = pNext(p_iter)
                        acc.next = NULL
                    else:
                        previous = p_iter
                        p_iter = pNext(p_iter)

                if zero_mat:
                    zero_mat = first == NULL

                mat[i - 1, j] = new_sage_polynomial(sage_ring, first)

        # Singular sometimes leaves zero matrix in the resolution. We can stop
        # when one is seen.
        if zero_mat:
            break

        res_mats.append(mat)

    return res_mats

cpdef tuple si2sa_resolution_graded(Resolution res, tuple degrees):
    """
    Pull the data from Singular resolution ``res`` to construct a Sage
    resolution.

    INPUT:

    - ``res`` -- Singular resolution

    - ``degrees`` -- list of integers or integer vectors

    The procedure is destructive, and ``res`` is not usable afterward.

    EXAMPLES::

        sage: from sage.libs.singular.singular import si2sa_resolution_graded
        sage: from sage.libs.singular.function import singular_function
        sage: module = singular_function("module")
        sage: mres = singular_function('mres')

        sage: S.<x,y,z,w> = PolynomialRing(QQ)
        sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
        sage: mod = module(I)
        sage: r = mres(mod, 0)
        sage: res_mats, res_degs = si2sa_resolution_graded(r, (1, 1, 1, 1))
        sage: res_mats
        [
                                         [ y  x]
                                         [-z -y]
        [z^2 - y*w y*z - x*w y^2 - x*z], [ w  z]
        ]
        sage: res_degs
        [[[2], [2], [2]], [[1, 1, 1], [1, 1, 1]]]
    """
    cdef ring *singular_ring
    cdef syStrategy singular_res
    cdef poly *p
    cdef poly *p_iter
    cdef poly *first
    cdef poly *previous
    cdef poly *acc
    cdef resolvente mods
    cdef ideal *mod
    cdef int i, j, k, idx, rank, nrows, ncols
    cdef int ngens = len(degrees)
    cdef bint zero_mat
    cdef list matdegs, exps

    from sage.matrix.constructor import matrix as _matrix

    singular_res = res._resolution[0]
    sage_ring = res.base_ring
    singular_ring = access_singular_ring(res.base_ring)

    if singular_res.minres != NULL:
        mods = singular_res.minres
    elif singular_res.fullres != NULL:
        mods = singular_res.fullres
    else:
        raise ValueError('Singular resolution is not usable')

    cdef list res_mats = []
    cdef list res_degs = []

    # length is the length of fullres. The length of minres
    # can be shorter. Hence we avoid SEGFAULT by stopping
    # at NULL pointer.
    for idx in range(singular_res.length):
        mod = <ideal *> mods[idx]
        if mod == NULL:
            break
        rank = mod.rank
        free_module = sage_ring ** rank

        nrows = rank
        ncols = mod.ncols # IDELEMS(mod)

        mat = _matrix(sage_ring, nrows, ncols)
        matdegs = []
        zero_mat = True
        for j in range(ncols):
            p = <poly *> mod.m[j]
            degs = []
            # code below copied and modified from to_sage_vector_destructive
            # in sage.libs.singular.function.Converter
            for i in range(1, rank + 1):
                previous = NULL
                acc = NULL
                first = NULL
                p_iter = p
                while p_iter != NULL:
                    if p_GetComp(p_iter, singular_ring) == i:
                        p_SetComp(p_iter, 0, singular_ring)
                        p_Setm(p_iter, singular_ring)
                        if acc == NULL:
                            first = p_iter
                        else:
                            acc.next = p_iter
                        acc = p_iter
                        if p_iter == p:
                            p = pNext(p_iter)
                        if previous != NULL:
                            previous.next = pNext(p_iter)
                        p_iter = pNext(p_iter)
                        acc.next = NULL
                    else:
                        previous = p_iter
                        p_iter = pNext(p_iter)

                if zero_mat:
                    zero_mat = first == NULL

                mat[i - 1, j] = new_sage_polynomial(sage_ring, first)

                # degree of a homogeneous polynomial can be computed from the
                # first monomial
                if first != NULL:
                    exps = singular_monomial_exponents(first, singular_ring)
                    deg = 0
                    for k in range(ngens):
                        deg += exps[k] * degrees[k]
                    degs.append(deg)
                else:
                    degs.append(None)

            matdegs.append(degs)  # store degrees of the column

        # Singular sometimes leaves zero matrix in the resolution. We can stop
        # when one is seen.
        if zero_mat:
            break

        res_mats.append(mat)
        res_degs.append(matdegs)

    return (res_mats, res_degs)


cdef number *sa2si_QQ(Rational r, ring *_ring) noexcept:
    """
    Create a singular number from a sage rational.

    INPUT:

    - ``r`` -- a sage rational number

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number


    TESTS::

        sage: P.<x,y,z> = QQ[]
        sage: P(0) + 1/2 - 2/4
        0
        sage: P(1/2) + 3/5 - 3/5
        1/2
        sage: P(2/3) + 1/4 - 1/4
        2/3
        sage: P(12345678901234567890/23) + 5/2 - 5/2
        12345678901234567890/23
    """
    if _ring != currRing: rChangeCurrRing(_ring)
    return nlInit2gmp( mpq_numref(r.value), mpq_denref(r.value),_ring.cf )

cdef number *sa2si_GFqGivaro(int quo, ring *_ring) noexcept:
    """
    Create a singular number in a small finite field.

    INPUT:

    - ``quo`` -- sage integer

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    Number field elements are represented as polynomials in the number field
    generator. In this case, ``quo`` is the integer resulting from evaluating
    that polynomial in the characteristic of the field.

    TESTS::

        sage: F = FiniteField(5^2)
        sage: type(F)
        <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
        sage: R.<x,y,z> = F[]
        sage: R(0) + 1
        1
        sage: R(F.gen()) + 1
        (z2 + 1)
        sage: R(F.gen()^2) + 1
        (z2 - 1)
    """
    if _ring != currRing:
        rChangeCurrRing(_ring)
    cdef number* n1
    cdef number* n2
    cdef number* a
    cdef number* coeff
    cdef number* apow1
    cdef number* apow2
    cdef int b = _ring.cf.ch

    a = _ring.cf.cfParameter(1, _ring.cf)

    apow1 = _ring.cf.cfInit(1, _ring.cf)
    n1 = _ring.cf.cfInit(0, _ring.cf)

    while quo!=0:
        coeff = _ring.cf.cfInit(quo%b, _ring.cf)

        if not _ring.cf.cfIsZero(coeff, _ring.cf):
            apow2 = _ring.cf.cfMult(coeff, apow1, _ring.cf)
            n2 = _ring.cf.cfAdd(apow2, n1, _ring.cf)
            _ring.cf.cfDelete(&apow2, _ring.cf)
            _ring.cf.cfDelete(&n1, _ring.cf)
            n1 = n2

        apow2 = _ring.cf.cfMult(apow1, a, _ring.cf)
        _ring.cf.cfDelete(&apow1, _ring.cf)
        apow1 = apow2

        quo = quo/b
        _ring.cf.cfDelete(&coeff, _ring.cf)

    _ring.cf.cfDelete(&apow1, _ring.cf)
    _ring.cf.cfDelete(&a, _ring.cf)
    return n1

cdef number *sa2si_GFqNTLGF2E(FFgf2eE elem, ring *_ring) noexcept:
    """
    Create a singular number from a sage element of a finite field of
    characteristic 2.

    INPUT:

    - ``elem`` -- a sage element of a ntl_gf2e finite field

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    TESTS::

        sage: F = FiniteField(2^20)
        sage: type(F)
        <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>
        sage: R.<x,y,z> = F[]
        sage: R(0)+1
        1
        sage: R(F.gen()) + 1
        (z20 + 1)
        sage: R(F.gen()^21) + 1
        (z20^11 + z20^10 + z20^8 + z20^7 + z20^6 + z20^5 + z20^2 + z20 + 1)
    """
    if _ring != currRing: rChangeCurrRing(_ring)
    cdef int i
    cdef number *n1
    cdef number *n2
    cdef number *a
    cdef number *coeff
    cdef number *apow1
    cdef number *apow2
    cdef GF2X_c rep = GF2E_rep(elem.x)

    if GF2X_deg(rep) >= 1:
        n1 = _ring.cf.cfInit(0, _ring.cf)
        a = _ring.cf.cfParameter(1,_ring.cf)
        apow1 = _ring.cf.cfInit(1, _ring.cf)

        for i from 0 <= i <= GF2X_deg(rep):
            coeff = _ring.cf.cfInit(GF2_conv_to_long(GF2X_coeff(rep,i)), _ring.cf)

            if not _ring.cf.cfIsZero(coeff,_ring.cf):
                apow2 = _ring.cf.cfMult(coeff, apow1,_ring.cf)
                n2 = _ring.cf.cfAdd(apow2, n1,_ring.cf)
                _ring.cf.cfDelete(&apow2, _ring.cf)
                _ring.cf.cfDelete(&n1, _ring.cf)
                n1 = n2

            apow2 = _ring.cf.cfMult(apow1, a,_ring.cf)
            _ring.cf.cfDelete(&apow1, _ring.cf)
            apow1 = apow2

            _ring.cf.cfDelete(&coeff, _ring.cf)

        _ring.cf.cfDelete(&apow1, _ring.cf)
        _ring.cf.cfDelete(&a, _ring.cf)
    else:
        n1 = _ring.cf.cfInit(GF2_conv_to_long(GF2X_coeff(rep,0)), _ring.cf)

    return n1

cdef number *sa2si_GFq_generic(object elem, ring *_ring) noexcept:
    """
    Create a singular number from a sage element of a generic finite field.

    INPUT:

    - ``elem`` -- a sage element of a generic finite field

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    TESTS::

        sage: F = FiniteField(3^20)
        sage: type(F)
        <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>
        sage: R.<x,y,z> = F[]
        sage: R(0) + 1
        1
        sage: R(F.gen()) + 1
        (z20 + 1)
        sage: R(F.gen()^21) + 1
        (z20^14 - z20^12 - z20^11 - z20^10 - z20^9 + z20^6 + z20^5 + z20^4 - z20^2 + z20 + 1)
    """
    cdef int i
    cdef number *n1
    cdef number *n2
    cdef number *a
    cdef number *coeff
    cdef number *apow1
    cdef number *apow2

    if _ring.cf.type in (n_Zn, n_Znm):
        return sa2si_ZZmod(elem, _ring)
    elem = elem.polynomial()

    if _ring != currRing: rChangeCurrRing(_ring)
    if elem.degree() > 0:
        n1 = _ring.cf.cfInit(0, _ring.cf)
        a = _ring.cf.cfParameter(1,_ring.cf)
        apow1 = _ring.cf.cfInit(1, _ring.cf)

        for i from 0 <= i <= elem.degree():
            coeff = _ring.cf.cfInit(int(elem[i]), _ring.cf)

            if not _ring.cf.cfIsZero(coeff,_ring.cf):
                apow2 = _ring.cf.cfMult(coeff, apow1,_ring.cf)
                n2 = _ring.cf.cfAdd(apow2, n1,_ring.cf)
                _ring.cf.cfDelete(&apow2, _ring.cf)
                _ring.cf.cfDelete(&n1, _ring.cf)
                n1 = n2

            apow2 = _ring.cf.cfMult(apow1, a,_ring.cf)
            _ring.cf.cfDelete(&apow1, _ring.cf)
            apow1 = apow2

            _ring.cf.cfDelete(&coeff, _ring.cf)

        _ring.cf.cfDelete(&apow1, _ring.cf)
        _ring.cf.cfDelete(&a, _ring.cf)
    else:
        n1 = _ring.cf.cfInit(int(elem), _ring.cf)

    return n1

cdef number *sa2si_transext_QQ(object elem, ring *_ring) noexcept:
    """
    Create a singular number from a sage element of a transcendental extension
    of the rationals.

    INPUT:

    - ``elem`` -- a sage element of a FractionField of polynomials over the rationals

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    TESTS::

        sage: F = PolynomialRing(QQ,'a,b').fraction_field()
        sage: F.inject_variables()
        Defining a, b
        sage: R.<x,y> = F[]
        sage: a*x
        a*x
        sage: I = R.ideal([a*x])
        sage: I
        Ideal (a*x) of Multivariate Polynomial Ring in x, y over Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
        sage: I.groebner_basis()
        [x]
        sage: I = R.ideal([a*x+b*y^2, (b+a)/(b-a)*x^3-3*y*x])
        sage: I.groebner_basis()
        [x^3 + (3*a - 3*b)/(a + b)*x*y, y^2 + a/b*x]
        sage: R.term_order()
        Degree reverse lexicographic term order

    ::

        sage: F = PolynomialRing(QQ,'a').fraction_field()
        sage: R.<x,y> = F[]
        sage: F.inject_variables()
        Defining a
        sage: a*x
        a*x
        sage: I = R.ideal([a*x+5*y^2, (1+a)/(1-a)*x^3-3*y*x])
        sage: I
        Ideal (5*y^2 + a*x, (-a - 1)/(a - 1)*x^3 - 3*x*y) of Multivariate Polynomial Ring in x, y over Fraction Field of Univariate Polynomial Ring in a over Rational Field
        sage: I.groebner_basis()
        [x^3 + (3*a - 3)/(a + 1)*x*y, y^2 + a/5*x]

    ::

        sage: F = PolynomialRing(QQ,'a,b').fraction_field()
        sage: R.<x,y> = PolynomialRing(F)
        sage: S.<x,y> = QQ[]
        sage: f = x + y + 1
        sage: R(f)
        x + y + 1
    """
    cdef int j
    cdef number *n1
    cdef number *a
    cdef number *naCoeff
    cdef number *numerator
    cdef number *denominator
    cdef number *cfnum
    cdef number *cfden
    cdef number *aux1
    cdef number *aux2
    cdef number *power
    cdef int ngens
    cdef int ex
    cdef nMapFunc nMapFuncPtr = NULL

    if _ring != currRing:
        rChangeCurrRing(_ring)

    ngens = elem.parent().ngens()

    nMapFuncPtr =  naSetMap(_ring.cf, currRing.cf) # choose correct mapping function

    if nMapFuncPtr is NULL:
        raise RuntimeError("Failed to determine nMapFuncPtr")

    numerdic = elem.numerator().monomial_coefficients()
    denomdic = elem.denominator().monomial_coefficients()

    if numerdic and not isinstance(list(numerdic)[0], (tuple, ETuple)):
        numerdic = {(k,):b for k,b in numerdic.items()}

    if denomdic and not isinstance(list(denomdic)[0], (tuple, ETuple)):
        denomdic = {(k,):b for k,b in denomdic.items()}

    if _ring != currRing:
        rChangeCurrRing(_ring)
    n1 = _ring.cf.cfInit(0, _ring.cf)
    numerator = _ring.cf.cfInit(0, _ring.cf)
    for (exponents, coef) in numerdic.items():
        numer = coef.numerator()
        cfnum = _ring.cf.cfInitMPZ((<Integer>numer).value, _ring.cf)
        denom = coef.denominator()
        cfden = _ring.cf.cfInitMPZ((<Integer>denom).value, _ring.cf)
        naCoeff = _ring.cf.cfDiv(cfnum, cfden, _ring.cf)
        _ring.cf.cfDelete(&cfnum, _ring.cf)
        _ring.cf.cfDelete(&cfden, _ring.cf)
        for (j, ex) in enumerate(exponents):
            a = _ring.cf.cfParameter(j+1, _ring.cf)
            _ring.cf.cfPower(a, ex, &power, _ring.cf)
            aux1 = naCoeff
            naCoeff = _ring.cf.cfMult(aux1, power, _ring.cf)
            _ring.cf.cfDelete(&aux1, _ring.cf)
            _ring.cf.cfDelete(&a, _ring.cf)
            _ring.cf.cfDelete(&power, _ring.cf)
        aux2 = numerator
        numerator = _ring.cf.cfAdd(aux2, naCoeff,_ring.cf)
        _ring.cf.cfDelete(&aux2, _ring.cf)

    if elem.denominator() != 1:
        denominator = _ring.cf.cfInit(0, _ring.cf)

        for (exponents, coef) in denomdic.items():
            numer = coef.numerator()
            cfnum = _ring.cf.cfInitMPZ((<Integer>numer).value, _ring.cf)
            denom = coef.denominator()
            cfden = _ring.cf.cfInitMPZ((<Integer>denom).value, _ring.cf)
            naCoeff = _ring.cf.cfDiv(cfnum, cfden, _ring.cf)
            _ring.cf.cfDelete(&cfnum, _ring.cf)
            _ring.cf.cfDelete(&cfden, _ring.cf)
            for (j, ex) in enumerate(exponents):
                a = _ring.cf.cfParameter(j+1, _ring.cf)
                _ring.cf.cfPower(a, ex, &power, _ring.cf)
                aux1 = naCoeff
                naCoeff = _ring.cf.cfMult(aux1, power, _ring.cf)
                _ring.cf.cfDelete(&aux1, _ring.cf)
                _ring.cf.cfDelete(&a, _ring.cf)
                _ring.cf.cfDelete(&power, _ring.cf)
            aux2 = denominator
            denominator = _ring.cf.cfAdd(aux2, naCoeff,_ring.cf)
            _ring.cf.cfDelete(&aux2, _ring.cf)

    else:
        denominator = _ring.cf.cfInit(1, _ring.cf)

    n1 = _ring.cf.cfDiv(numerator, denominator, _ring.cf)

    _ring.cf.cfDelete(&numerator, _ring.cf)
    _ring.cf.cfDelete(&denominator, _ring.cf)
    _ring.cf.cfDelete(&naCoeff, _ring.cf)
    _ring.cf.cfDelete(&a, _ring.cf)

    return n1

cdef number *sa2si_transext_FF(object elem, ring *_ring) noexcept:
    """
    Create a singular number from a sage element of a transcendental extension
    of a prime field.

    INPUT:

    - ``elem`` -- a sage element of a FractionField of polynomials over the rationals

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    TESTS::

        sage: F = PolynomialRing(FiniteField(7),'a,b').fraction_field()
        sage: R.<x,y,z> = F[]
        sage: n = R(5)
        sage: n + n
        3
        sage: Integer(n)
        5
    """
    cdef int j
    cdef number *n1
    cdef number *a
    cdef number *naCoeff
    cdef number *numerator
    cdef number *denominator
    cdef number *aux1
    cdef number *aux2
    cdef int ngens
    cdef int ex
    cdef nMapFunc nMapFuncPtr = NULL

    if _ring != currRing:
        rChangeCurrRing(_ring)

    ngens = elem.parent().ngens()

    nMapFuncPtr =  naSetMap(_ring.cf, currRing.cf) # choose correct mapping function

    if nMapFuncPtr is NULL:
        raise RuntimeError("Failed to determine nMapFuncPtr")

    numerdic = elem.numerator().monomial_coefficients()
    denomdic = elem.denominator().monomial_coefficients()

    if numerdic and not isinstance(list(numerdic)[0], (tuple, ETuple)):
        numerdic = {(k,):b for k,b in numerdic.items()}

    if denomdic and not isinstance(list(denomdic)[0], (tuple, ETuple)):
        denomdic = {(k,):b for k,b in denomdic.items()}

    if _ring != currRing:
        rChangeCurrRing(_ring)
    numerator = _ring.cf.cfInit(0, _ring.cf)
    for (exponents, coef) in numerdic.items():
        naCoeff = _ring.cf.cfInit(<int>coef, _ring.cf)
        for (j, ex) in enumerate(exponents):
            a = _ring.cf.cfParameter(j+1, _ring.cf)
            for k in range(ex):
                aux1 = naCoeff
                naCoeff = _ring.cf.cfMult(aux1, a, _ring.cf)
                _ring.cf.cfDelete(&aux1, _ring.cf)
            _ring.cf.cfDelete(&a, _ring.cf)
        aux2 = numerator
        numerator = _ring.cf.cfAdd(aux2, naCoeff, _ring.cf)
        _ring.cf.cfDelete(&naCoeff, _ring.cf)
        _ring.cf.cfDelete(&aux2, _ring.cf)

    if elem.denominator() != 1:
        denominator = _ring.cf.cfInit(0, _ring.cf)

        for (exponents, coef) in denomdic.items():
            naCoeff = _ring.cf.cfInit(<int>coef, _ring.cf)
            for (j, ex) in enumerate(exponents):
                a = _ring.cf.cfParameter(j+1, _ring.cf)
                for k in range(ex):
                    aux1 = naCoeff
                    naCoeff = _ring.cf.cfMult(aux1, a, _ring.cf)
                    _ring.cf.cfDelete(&aux1, _ring.cf)
                _ring.cf.cfDelete(&a, _ring.cf)
            aux2 = denominator
            denominator = _ring.cf.cfAdd(aux2, naCoeff,_ring.cf)
            _ring.cf.cfDelete(&naCoeff, _ring.cf)
            _ring.cf.cfDelete(&aux2, _ring.cf)

    else:
        denominator = _ring.cf.cfInit(1, _ring.cf)

    n1 = _ring.cf.cfDiv(numerator, denominator, _ring.cf)

    _ring.cf.cfDelete(&numerator, _ring.cf)
    _ring.cf.cfDelete(&denominator, _ring.cf)
    _ring.cf.cfDelete(&a, _ring.cf)

    return n1

cdef number *sa2si_NF(object elem, ring *_ring) noexcept:
    """
    Create a singular number from a sage element of a number field.

    INPUT:

    - ``elem`` -- a sage element of a NumberField

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    TESTS::

        sage: x = polygen(ZZ, 'x')
        sage: F = NumberField(x^3 + x + 1, 'a')
        sage: type(F)
        <class 'sage.rings.number_field.number_field.NumberField_absolute_with_category'>
        sage: R.<x,y,z> = F[]
        sage: R(0) + 1
        1
        sage: R(1)
        1
        sage: R(F.gen()) + 1
        (a + 1)
        sage: R(F.gen()^5) + 1
        (-a^2 + a + 2)
    """
    cdef int i
    cdef number *n1
    cdef number *n2
    cdef number *a
    cdef number *nlCoeff
    cdef number *naCoeff
    cdef number *apow1
    cdef number *apow2

    cdef nMapFunc nMapFuncPtr = NULL

    nMapFuncPtr =  naSetMap(_ring.cf, currRing.cf) # choose correct mapping function

    if nMapFuncPtr is NULL:
        raise RuntimeError("Failed to determine nMapFuncPtr")

    elem = list(elem)

    if _ring != currRing:
        rChangeCurrRing(_ring)
    n1 = _ring.cf.cfInit(0, _ring.cf)
    a = _ring.cf.cfParameter(1, _ring.cf)
    apow1 = _ring.cf.cfInit(1, _ring.cf)

    cdef char *_name

    # the result of nlInit2gmp() is in a plain polynomial ring over QQ (not an extension ring!),
    # so we have to get/create one:
    #
    # todo: reuse qqr/ get an existing Singular polynomial ring over Q.
    _name = omStrDup("a")
    cdef char **_ext_names
    _ext_names = <char**>omAlloc0(sizeof(char*))
    _ext_names[0] = omStrDup(_name)
    qqr = rDefault( 0, 1, _ext_names)
    rComplete(qqr,1)
    qqr.ShortOut = 0

    nMapFuncPtr =  naSetMap(qqr.cf, _ring.cf)  # choose correct mapping function
    cdef poly *_p
    for i from 0 <= i < len(elem):
        nlCoeff = nlInit2gmp( mpq_numref((<Rational>elem[i]).value), mpq_denref((<Rational>elem[i]).value),  qqr.cf )
        naCoeff = nMapFuncPtr(nlCoeff, qqr.cf, _ring.cf)
        nlDelete(&nlCoeff, _ring.cf)

        # faster would be to assign the coefficient directly
        apow2 = _ring.cf.cfMult(naCoeff, apow1,_ring.cf)
        n2 = _ring.cf.cfAdd(apow2, n1,_ring.cf)
        _ring.cf.cfDelete(&apow2, _ring.cf)
        _ring.cf.cfDelete(&n1, _ring.cf)
        _ring.cf.cfDelete(&naCoeff, _ring.cf)
        n1 = n2

        apow2 = _ring.cf.cfMult(apow1, a,_ring.cf)
        _ring.cf.cfDelete(&apow1, _ring.cf)
        apow1 = apow2

    _ring.cf.cfDelete(&apow1, _ring.cf)
    _ring.cf.cfDelete(&a, _ring.cf)

    return n1

cdef number *sa2si_ZZ(Integer d, ring *_ring) noexcept:
    """
    Create a singular number from a sage Integer.

    INPUT:

    - ``elem`` -- a sage Integer

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    OUTPUT:

    - A (pointer to) a singular number

    TESTS::

        sage: P.<x,y,z> = ZZ[]
        sage: P(0) + 1 - 1
        0
        sage: P(1) + 1 - 1
        1
        sage: P(2) + 1 - 1
        2
        sage: P(12345678901234567890) + 2 - 2
        12345678901234567890
    """
    if _ring != currRing: rChangeCurrRing(_ring)
    cdef number *n = nrzInit(0, _ring.cf)
    mpz_set(<mpz_ptr>n, d.value)
    return <number*>n

cdef inline number *sa2si_ZZmod(IntegerMod_abstract d, ring *_ring) noexcept:
    """
    Create a singular number from a sage element of a IntegerModRing.

    INPUT:

    - ``elem`` -- a sage IntegerMod

    - ``_ ring`` -- a (pointer to) a singular ring, where the result will live

    TESTS::

        sage: P.<x,y,z> = Integers(10)[]
        sage: P(3)
        3
        sage: P(13)
        3

        sage: P.<x,y,z> = Integers(16)[]
        sage: P(3)
        3
        sage: P(19)
        3

        sage: P.<x,y,z> = Integers(3^2)[]
        sage: P(3)
        3
        sage: P(12)
        3

        sage: P.<x,y,z> = Integers(2^32)[]
        sage: P(2^32-1)
        -1

        sage: P(3)
        3

        sage: P.<x,y,z> = Integers(17^20)[]
        sage: P(17^19 + 3)
        239072435685151324847156

        sage: P(3)
        3
    """
    if _ring != currRing: rChangeCurrRing(_ring)

    cdef number *nn

    cdef int64_t _d
    cdef char *_name
    cdef char **_ext_names

    cdef nMapFunc nMapFuncPtr = NULL

    if _ring.cf.type == n_unknown:
        return n_Init(int(d),_ring.cf)

    if _ring.cf.type == n_Z2m:
        _d = long(d)
        return nr2mMapZp(<number *>_d, currRing.cf, _ring.cf)
    elif _ring.cf.type == n_Zn or _ring.cf.type == n_Znm:
        lift = d.lift()

        # if I understand nrnMapGMP/nMapFuncPtr correctly we need first
        # a source value in ZZr
        # create ZZr, a plain polynomial ring over ZZ with one variable.
        #
        # todo (later): reuse ZZr
        _name = omStrDup("a")
        _ext_names = <char**>omAlloc0(sizeof(char*))
        _ext_names[0] = omStrDup(_name)
        _cf = nInitChar(n_Z, NULL)  # integer coefficient ring
        ZZr = rDefault (_cf, 1, _ext_names)
        rComplete(ZZr, 1)
        ZZr.ShortOut = 0

        nn = nrzInit(0, ZZr.cf)
        mpz_set(<mpz_ptr>nn, (<Integer>lift).value)
        nMapFuncPtr  = nrnSetMap( ZZr.cf, _ring.cf)

        return nMapFuncPtr(nn, ZZr.cf, _ring.cf)
    else:
        raise ValueError

cdef object si2sa(number *n, ring *_ring, object base):
    r"""
    Create a sage number from a singular one.

    INPUT:

    - ``n`` -- a (pointer to) a singular number

    - ``_ring`` -- a (pointer to) the singular ring where ``n`` lives

    - ``object`` -- the sage parent where the result will live

    OUTPUT:

    An element of ``base``
    """
    if isinstance(base, FiniteField_prime_modn) and _ring.cf.type == n_Zp:
        return base(_ring.cf.cfInt(n, _ring.cf))

    elif isinstance(base, RationalField):
        return si2sa_QQ(n,&n,_ring)

    elif isinstance(base, IntegerRing_class):
        return si2sa_ZZ(n,_ring)

    elif isinstance(base, FiniteField_givaro):
        return si2sa_GFqGivaro(n, _ring, base._cache)

    elif isinstance(base, FiniteField_ntl_gf2e):
        return si2sa_GFqNTLGF2E(n, _ring, <Cache_ntl_gf2e>base._cache)

    elif isinstance(base, FiniteField):
        return si2sa_GFq_generic(n, _ring, base)

    elif isinstance(base, NumberField) and base.is_absolute():
        return si2sa_NF(n, _ring, base)

    elif isinstance(base, FractionField_generic) and isinstance(base.base(), (MPolynomialRing_libsingular, PolynomialRing_field)) and isinstance(base.base_ring(), RationalField):
        return si2sa_transext_QQ(n, _ring, base)

    elif isinstance(base, FractionField_generic) and isinstance(base.base(), (MPolynomialRing_libsingular, PolynomialRing_field)) and isinstance(base.base_ring(), FiniteField_prime_modn):
        return si2sa_transext_FF(n, _ring, base)

    elif isinstance(base, IntegerModRing_generic):
        return si2sa_ZZmod(n, _ring, base)

    else:
        raise ValueError("cannot convert from SINGULAR number")

cdef number *sa2si(Element elem, ring * _ring) noexcept:
    r"""
    Create a singular number from a sage one.

    INPUT:

    - ``elem`` -- a sage element from a parent; the parent must have a
      corresponding singular coefficient type

    - ``_ring`` -- a (pointer to) the singular ring where the result will live

    OUTPUT:

    a (pointer to) a singular number
    """
    cdef int i = 0

    if isinstance(elem._parent, FiniteField_prime_modn) and _ring.cf.type == n_Zp:
        return n_Init(int(elem),_ring.cf)

    elif isinstance(elem._parent, RationalField):
        return sa2si_QQ(elem, _ring)

    elif isinstance(elem._parent, IntegerRing_class):
        return sa2si_ZZ(elem, _ring)

    elif isinstance(elem._parent, FiniteField_givaro):
        return sa2si_GFqGivaro( (<FFgivE>elem)._cache.objectptr.convert(i, (<FFgivE>elem).element ), _ring )

    elif isinstance(elem._parent, FiniteField_ntl_gf2e):
        return sa2si_GFqNTLGF2E(elem, _ring)

    elif isinstance(elem._parent, FiniteField):
        return sa2si_GFq_generic(elem, _ring)

    elif isinstance(elem._parent, NumberField) and elem._parent.is_absolute():
        return sa2si_NF(elem, _ring)
    elif isinstance(elem._parent, IntegerModRing_generic):
        return sa2si_ZZmod(elem, _ring)
    elif isinstance(elem._parent, FractionField_generic) and isinstance(elem._parent.base(), (MPolynomialRing_libsingular, PolynomialRing_field)):
        if isinstance(elem._parent.base().base_ring(), RationalField):
            return sa2si_transext_QQ(elem, _ring)
        elif isinstance(elem._parent.base().base_ring(), FiniteField_prime_modn):
            return sa2si_transext_FF(elem, _ring)

    raise ValueError("cannot convert to SINGULAR number")

cdef object si2sa_intvec(intvec *v):
    r"""
    Create a sage tuple from a singular vector of integers.

    INPUT:

    - ``v`` -- a (pointer to) singular intvec

    OUTPUT:

    a sage tuple
    """
    cdef int r
    cdef list l = list()
    for r in range(v.length()):
        l.append(v.get(r))
    return tuple(l)

cdef object si2sa_bigintvec(bigintmat *v):
    r"""
    Create a sage tuple from a singular vector of big integers.

    INPUT:

    - ``v`` -- a (pointer to) singular bigintmat

    OUTPUT:

    a sage tuple
    """
    cdef int r
    cdef list l = list()
    for r in range(v.length()):
        n = v.get(r)
        l.append(si2sa_QQ(n, &n, currRing))
    return tuple(l)

# ==============
# Initialisation
# ==============

cdef extern from *: # hack to get at cython macro
    int unlikely(int)

from posix.dlfcn cimport dlopen, dlclose, dlerror, RTLD_LAZY, RTLD_GLOBAL

cdef int overflow_check(unsigned long e, ring *_ring) except -1:
    """
    Raise an :exc:`OverflowError` if e is > max degree per variable.

    INPUT:

    - ``e`` -- some integer representing a degree

    - ``_ring`` -- a pointer to some ring

    Whether an overflow occurs or not partially depends

    on the number of variables in the ring. See github issue
    :issue:`11856`. With Singular 4, it is by default optimized
    for at least 4 variables on 64-bit and 2 variables on 32-bit,
    which in both cases makes a maximal default exponent of
    2^16-1.

    EXAMPLES::

        sage: P.<x,y> = QQ[]
        sage: y^(2^30)
        Traceback (most recent call last):             # 32-bit
        ...                                            # 32-bit
        OverflowError: exponent overflow (1073741824)  # 32-bit
        y^1073741824  # 64-bit
        sage: y^2^32
        Traceback (most recent call last):
        ...
        OverflowError: Python int too large to convert to C unsigned long  # 32-bit
        OverflowError: exponent overflow (4294967296)  # 64-bit
    """
    if unlikely(e > _ring.bitmask):
        raise OverflowError("exponent overflow (%d)"%(e))

cdef init_libsingular():
    """
    This initializes the SINGULAR library. This is a hack to some
    extent.

    SINGULAR has a concept of compiled extension modules similar to
    Sage. For this, the compiled modules need to see the symbols from
    the main program. However, SINGULAR is a shared library in this
    context these symbols are not known globally. The work around so
    far is to load the library again and to specify ``RTLD_GLOBAL``.
    """
    global singular_options
    global singular_verbose_options
    global WerrorS_callback
    global error_messages

    cdef void *handle = NULL

    # This is a workaround for https://github.com/Singular/Singular/issues/1113
    # and can be removed once that fix makes it into release of Singular that
    # is supported by sage.
    from sage.features import FeatureNotPresentError
    from sage.features.singular import Singular
    from os.path import dirname
    try:
        singular_executable = Singular().absolute_filename()
    except FeatureNotPresentError:
        pass
    else:
        os.environ["SINGULAR_BIN_DIR"] = dirname(singular_executable)

    # reload the current module to force reload of libSingular (see #33446)
    lib = str_to_bytes(__loader__.path, FS_ENCODING, "surrogateescape")
    handle = dlopen(lib, RTLD_GLOBAL|RTLD_LAZY)
    if not handle:
        err = dlerror()
        raise RuntimeError(f"Could not reload Singular library with RTLD_GLOBAL ({err})")

    # load SINGULAR
    siInit(lib)

    if handle:
        dlclose(handle)

    # we set and save some global Singular options
    singular_options = singular_options | Sy_bit(OPT_REDSB) | Sy_bit(OPT_INTSTRATEGY) | Sy_bit(OPT_REDTAIL) | Sy_bit(OPT_REDTHROUGH)
    global _saved_options
    global _saved_verbose_options
    _saved_options = (int(singular_options), 0, 0)
    _saved_verbose_options = int(singular_verbose_options)

    #On(SW_USE_NTL)
    On(SW_USE_EZGCD)
    Off(SW_USE_NTL_SORT)

    WerrorS_callback = libsingular_error_callback

    error_messages = []

# Save/restore the PATH because libSingular clobbers it:
# https://github.com/Singular/Singular/issues/1119
saved_PATH = os.environ["PATH"]
init_libsingular()
os.environ["PATH"] = saved_PATH

cdef bint catching_error = False

cdef void libsingular_error_callback(const_char_ptr s) noexcept:
    _s = char_to_str(s)
    if catching_error:
        error_messages.append(_s)
    else:
        warn(f"error in Singular ignored: {_s}")

cdef int start_catch_error() except -1:
    """
    Helper function to convert Singular errors to Python exceptions.

    Must be used as follows::

        start_catch_error()
        ...
        s = check_error()  # nonempty tuple[str, ...] (error messages) or None
        if s:
            # at this point global variable ``error_messages`` is cleared
            raise RuntimeError(...)

    Return value is ignored, only used for exception handling.

    Note that :func:`check_error` can only be called exactly once.
    """
    global errorreported, catching_error, error_messages
    if catching_error:
        raise RuntimeError("internal error: previous start_catch_error not ended with check_error")
    catching_error = True

    if errorreported:
        warn(f"error in Singular ignored: {', '.join(error_messages)}")
        errorreported = False
        error_messages.clear()
    else:
        assert not error_messages
    return 0

cdef object check_error():
    """
    See :func:`start_catch_error`.
    """
    global errorreported, catching_error, error_messages
    if not catching_error:
        raise RuntimeError("check_error must be preceded with start_catch_error")
    catching_error = False

    if errorreported:
        result = tuple(error_messages)
        assert result
        errorreported = False
        error_messages.clear()
        return result
    else:
        assert not error_messages
        return None


def get_resource(id):
    """
    Return a Singular "resource".

    INPUT:

    - ``id`` -- a single-character string; see
      https://github.com/Singular/Singular/blob/spielwiese/resources/feResource.cc

    OUTPUT: string or ``None``

    EXAMPLES::

        sage: from sage.libs.singular.singular import get_resource
        sage: get_resource('D')            # SINGULAR_DATA_DIR
        '...'
        sage: get_resource('i')            # SINGULAR_INFO_FILE
        '.../singular...'
        sage: get_resource('7') is None    # not defined
        True
    """
    cdef char *result = feGetResource(<char>ord(id))
    if result == NULL:
        return None
    return bytes_to_str(result)
