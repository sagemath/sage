"""
libSingular: Conversion Routines and Initialisation

AUTHOR:

- Martin Albrecht <malb@informatik.uni-bremen.de>
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

from libc.stdint cimport int64_t
from sage.libs.singular.decl cimport *

from sage.rings.rational_field import RationalField
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_rings.finite_field_ntl_gf2e import FiniteField_ntl_gf2e
from sage.libs.pari.all import pari
from sage.libs.gmp.all cimport *

from sage.cpython.string import FS_ENCODING
from sage.cpython.string cimport str_to_bytes, char_to_str

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular

_saved_options = (int(0),0,0)

cdef Rational si2sa_QQ(number *n, number **nn, ring *_ring):
    """
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
        <type 'sage.rings.rational.Rational'>
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
        <type 'sage.rings.integer.Integer'>
    """
    cdef Integer z
    z = Integer()
    z.set_from_mpz(<mpz_ptr>n)
    return z

cdef FFgivE si2sa_GFqGivaro(number *n, ring *_ring, Cache_givaro cache):
    """
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
    TESTS::

        sage: K.<a> = GF(2^20)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        a^11 + a^10 + a^8 + a^7 + a^6 + a^5 + a^2 + a
        sage: type(f.lc())
        <type 'sage.rings.finite_rings.element_ntl_gf2e.FiniteField_ntl_gf2eElement'>
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
    TESTS::

        sage: K.<a> = GF(3^16)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        a^12 + a^11 + a^9 + a^8 + a^7 + 2*a^6 + a^5
        sage: type(f.lc())
        <type 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt'>

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

cdef object si2sa_NF(number *n, ring *_ring, object base):
    """
    TESTS::

        sage: K.<a> = NumberField(x^2 - 2)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        1024*a
        sage: type(f.lc())
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
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

cdef number *sa2si_QQ(Rational r, ring *_ring):
    """
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

cdef number *sa2si_GFqGivaro(int quo, ring *_ring):
    """
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

cdef number *sa2si_GFqNTLGF2E(FFgf2eE elem, ring *_ring):
    """
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
                _ring.cf.cfDelete(&n1, _ring.cf);
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

cdef number *sa2si_GFq_generic(object elem, ring *_ring):
    """
    """
    cdef int i
    cdef number *n1
    cdef number *n2
    cdef number *a
    cdef number *coeff
    cdef number *apow1
    cdef number *apow2
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
                _ring.cf.cfDelete(&n1, _ring.cf);
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

cdef number *sa2si_NF(object elem, ring *_ring):
    """
    """
    cdef int i
    cdef number *n1
    cdef number *n2
    cdef number *a
    cdef number *nlCoeff
    cdef number *naCoeff
    cdef number *apow1
    cdef number *apow2

    cdef nMapFunc nMapFuncPtr = NULL;

    nMapFuncPtr =  naSetMap(_ring.cf, currRing.cf) # choose correct mapping function

    if (nMapFuncPtr is NULL):
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
    qqr = rDefault( 0, 1, _ext_names);
    rComplete(qqr,1)
    qqr.ShortOut = 0


    nMapFuncPtr =  naSetMap( qqr.cf , _ring.cf ) # choose correct mapping function
    cdef poly *_p
    for i from 0 <= i < len(elem):
        nlCoeff = nlInit2gmp( mpq_numref((<Rational>elem[i]).value), mpq_denref((<Rational>elem[i]).value),  qqr.cf )
        naCoeff = nMapFuncPtr(nlCoeff, qqr.cf , _ring.cf )
        nlDelete(&nlCoeff, _ring.cf)

        # faster would be to assign the coefficient directly
        apow2 = _ring.cf.cfMult(naCoeff, apow1,_ring.cf)
        n2 = _ring.cf.cfAdd(apow2, n1,_ring.cf)
        _ring.cf.cfDelete(&apow2, _ring.cf)
        _ring.cf.cfDelete(&n1, _ring.cf);
        _ring.cf.cfDelete(&naCoeff, _ring.cf)
        n1 = n2

        apow2 = _ring.cf.cfMult(apow1, a,_ring.cf)
        _ring.cf.cfDelete(&apow1, _ring.cf)
        apow1 = apow2

    _ring.cf.cfDelete(&apow1, _ring.cf)
    _ring.cf.cfDelete(&a, _ring.cf)

    return n1

cdef number *sa2si_ZZ(Integer d, ring *_ring):
    """
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

cdef inline number *sa2si_ZZmod(IntegerMod_abstract d, ring *_ring):
    """
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
        4294967295

        sage: P(3)
        3

        sage: P.<x,y,z> = Integers(17^20)[]
        sage: P(17^19 + 3)
        239072435685151324847156

        sage: P(3)
        3
    """
    nr2mModul = d.parent().characteristic()
    if _ring != currRing: rChangeCurrRing(_ring)

    cdef number *nn

    cdef int64_t _d
    cdef char *_name
    cdef char **_ext_names

    cdef nMapFunc nMapFuncPtr = NULL;

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
        _cf = nInitChar( n_Z, NULL) # integer coefficient ring
        ZZr = rDefault (_cf ,1, _ext_names)
        rComplete(ZZr,1)
        ZZr.ShortOut = 0

        nn = nrzInit(0, ZZr.cf)
        mpz_set(<mpz_ptr>nn, (<Integer>lift).value)
        nMapFuncPtr  = nrnSetMap( ZZr.cf, _ring.cf)

        return nMapFuncPtr(nn, ZZr.cf, _ring.cf)
    else:
        raise ValueError

cdef object si2sa(number *n, ring *_ring, object base):
    if isinstance(base, FiniteField_prime_modn):
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

    elif isinstance(base, IntegerModRing_generic):
        if _ring.cf.type == n_unknown:
            return base(_ring.cf.cfInt(n, _ring.cf))
        return si2sa_ZZmod(n, _ring, base)

    else:
        raise ValueError("cannot convert from SINGULAR number")

cdef number *sa2si(Element elem, ring * _ring):
    cdef int i = 0
    if isinstance(elem._parent, FiniteField_prime_modn):
        return n_Init(int(elem),_ring)

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
        if _ring.cf.type == n_unknown:
            return n_Init(int(elem),_ring)
        return sa2si_ZZmod(elem, _ring)
    else:
        raise ValueError("cannot convert to SINGULAR number")


cdef object si2sa_intvec(intvec *v):
    cdef int r
    cdef list l = list()
    for r in range(v.length()):
        l.append(v.get(r))
    return tuple(l)

# ==============
# Initialisation
# ==============

cdef extern from *: # hack to get at cython macro
    int unlikely(int)

cdef extern from "dlfcn.h":
    void *dlopen(char *, long)
    char *dlerror()
    void dlclose(void *handle)

cdef extern from "dlfcn.h":
    cdef long RTLD_LAZY
    cdef long RTLD_GLOBAL

cdef int overflow_check(unsigned long e, ring *_ring) except -1:
    """
    Raise an ``OverflowError`` if e is > max degree per variable.

    INPUT:

    - ``e`` -- some integer representing a degree.

    - ``_ring`` -- a pointer to some ring.

    Whether an overflow occurs or not partially depends
    on the number of variables in the ring. See trac ticket
    :trac:`11856`. With Singular 4, it is by default optimized
    for at least 4 variables on 64-bit and 2 variables on 32-bit,
    which in both cases makes a maximal default exponent of
    2^16-1.

    EXAMPLES::

        sage: P.<x,y> = QQ[]
        sage: y^(2^16-1)
        y^65535
        sage: y^2^16
        Traceback (most recent call last):
        ...
        OverflowError: exponent overflow (65536)
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

    from sage.env import SINGULAR_SO
    if not SINGULAR_SO or not os.path.exists(SINGULAR_SO):
        raise RuntimeError(
            "libSingular not found--a working Singular install "
            "is required for Sage to work")

    lib = str_to_bytes(SINGULAR_SO, FS_ENCODING, "surrogateescape")

    handle = dlopen(lib, RTLD_GLOBAL|RTLD_LAZY)
    if not handle:
        err = dlerror()
        raise ImportError(f"cannot load Singular library from {SINGULAR_SO} ({err})")

    # load SINGULAR
    siInit(lib)

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

# call the init routine
init_libsingular()

cdef void libsingular_error_callback(const_char_ptr s):
    _s = char_to_str(s)
    error_messages.append(_s)
