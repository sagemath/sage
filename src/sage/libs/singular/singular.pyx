"""
libSingular: Conversion Routines and Initialisation

AUTHOR:

- Martin Albrecht <malb@informatik.uni-bremen.de>
"""
###############################################################################
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
###############################################################################
from __future__ import print_function

include "sage/libs/ntl/decl.pxi"

cdef extern from "limits.h":
    long INT_MAX
    long INT_MIN

import os

from sage.libs.singular.decl cimport intvec
from sage.libs.singular.decl cimport SR_HDL, SR_INT, SR_TO_INT
from sage.libs.singular.decl cimport singular_options, singular_verbose_options
from sage.libs.singular.decl cimport On, Off, SW_USE_NTL, SW_USE_NTL_GCD_0, SW_USE_EZGCD, SW_USE_NTL_SORT, SW_USE_NTL_GCD_P
from sage.libs.singular.decl cimport napoly, lnumber, Sy_bit, OPT_REDSB, OPT_INTSTRATEGY, OPT_REDTAIL, OPT_REDTHROUGH
from sage.libs.singular.decl cimport nlGetNumerator, nlGetDenom, nlDelete, nlInit2gmp
from sage.libs.singular.decl cimport naIsOne, naIsOne, naIsZero, naPar, naInit, naAdd, naMult, naDelete, naMap00
from sage.libs.singular.decl cimport napGetCoeff, napGetExpFrom, pNext
from sage.libs.singular.decl cimport nrzInit, nr2mMapZp, nrnMapGMP
from sage.libs.singular.decl cimport siInit
from sage.libs.singular.decl cimport n_Init
from sage.libs.singular.decl cimport rChangeCurrRing, currRing
from sage.libs.singular.decl cimport WerrorS_callback, const_char_ptr

from sage.rings.rational_field import RationalField
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_rings.finite_field_ntl_gf2e import FiniteField_ntl_gf2e
from sage.libs.pari.all import pari
from sage.libs.gmp.all cimport *

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular

_saved_options = (int(0),0,0)

cdef Rational si2sa_QQ(number *n, ring *_ring):
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
    ##  structures aligned on 4 byte boundaries and therefor have last bit zero.
    ##  (The second bit is reserved as tag to allow extensions of this scheme.)
    ##  Using immediates as pointers and dereferencing them gives address errors.
    nom = nlGetNumerator(n, _ring)
    mpz_init(nom_z)

    if (SR_HDL(nom) & SR_INT): mpz_set_si(nom_z, SR_TO_INT(nom))
    else: mpz_set(nom_z,nom.z)

    mpq_set_num(_z,nom_z)
    nlDelete(&nom,_ring)
    mpz_clear(nom_z)

    denom = nlGetDenom(n, _ring)
    mpz_init(denom_z)

    if (SR_HDL(denom) & SR_INT): mpz_set_si(denom_z, SR_TO_INT(denom))
    else: mpz_set(denom_z,denom.z)

    mpq_set_den(_z, denom_z)
    nlDelete(&denom,_ring)
    mpz_clear(denom_z)

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
    cdef napoly *z
    cdef int c, e
    cdef int a
    cdef int ret
    cdef int order

    if naIsZero(n):
        return cache._zero_element
    elif naIsOne(n):
        return cache._one_element
    z = (<lnumber*>n).z

    a = cache.objectptr.indeterminate()
    ret = cache.objectptr.zero
    order = cache.objectptr.cardinality() - 1

    while z:
        c = cache.objectptr.initi(c, <long>napGetCoeff(z))
        e = napGetExpFrom(z,1, _ring)
        if e == 0:
            ret = cache.objectptr.add(ret, c, ret)
        else:
            a = ( e * cache.objectptr.indeterminate() ) % order
            ret = cache.objectptr.axpy(ret, c, a, ret)
        z = <napoly*>pNext(<poly*>z)
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
    cdef napoly *z
    cdef long c
    cdef int e
    cdef FFgf2eE a
    cdef FFgf2eE ret

    if naIsZero(n):
        return cache._zero_element
    elif naIsOne(n):
        return cache._one_element
    z = (<lnumber*>n).z

    a = cache._gen
    ret = cache._zero_element

    while z:
        c = <long>napGetCoeff(z)
        e = napGetExpFrom(z,1, _ring)
        ret += c * a**e
        z = <napoly*>pNext(<poly*>z)
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
    cdef napoly *z
    cdef long c
    cdef int e
    cdef object a
    cdef object ret

    if naIsZero(n):
        return base.zero()
    elif naIsOne(n):
        return base.one()
    z = (<lnumber*>n).z

    a = base.gen()
    ret = base.zero()

    while z:
        c = <long>napGetCoeff(z)
        e = napGetExpFrom(z,1, _ring)
        if e == 0:
            ret = ret + c
        elif c != 0:
            ret = ret  + c * a**e
        z = <napoly*>pNext(<poly*>z)
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
    cdef napoly *z
    cdef number *c
    cdef int e
    cdef object a
    cdef object ret

    if naIsZero(n):
        return base._zero_element
    elif naIsOne(n):
        return base._one_element
    z = (<lnumber*>n).z

    a = base.gen()
    ret = base(0)

    while z:
        c = napGetCoeff(z)
        coeff = si2sa_QQ(c, _ring)
        e = napGetExpFrom(z,1, _ring)
        if e == 0:
            ret = ret + coeff
        elif coeff != 0:
            ret = ret + coeff * a**e
        z = <napoly*>pNext(<poly*>z)
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
    if _ring.ringtype == 1:
        return base(<long>n)
    else:
        ret = Integer()
        ret.set_from_mpz(<mpz_ptr>n)
        return base(ret)

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
    return nlInit2gmp( mpq_numref(r.value), mpq_denref(r.value) )

cdef number *sa2si_GFqGivaro(int quo, ring *_ring):
    """
    """
    if _ring != currRing: rChangeCurrRing(_ring)
    cdef number *n1
    cdef number *n2
    cdef number *a
    cdef number *coeff
    cdef number *apow1
    cdef number *apow2
    cdef int b = - _ring.ch

    a = naPar(1)

    apow1 = naInit(1, _ring)
    n1 = naInit(0, _ring)

    while quo!=0:
        coeff = naInit(quo%b, _ring)

        if not naIsZero(coeff):
            apow2 = naMult(coeff, apow1)
            n2 = naAdd(apow2, n1)
            naDelete(&apow2, _ring)
            naDelete(&n1, _ring)
            n1 = n2

        apow2 = naMult(apow1, a)
        naDelete(&apow1, _ring)
        apow1 = apow2

        quo = quo/b
        naDelete(&coeff, _ring)

    naDelete(&apow1, _ring)
    naDelete(&a, _ring)
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
        n1 = naInit(0, _ring)
        a = naPar(1)
        apow1 = naInit(1, _ring)

        for i from 0 <= i <= GF2X_deg(rep):
            coeff = naInit(GF2_conv_to_long(GF2X_coeff(rep,i)), _ring)

            if not naIsZero(coeff):
                apow2 = naMult(coeff, apow1)
                n2 = naAdd(apow2, n1)
                naDelete(&apow2, _ring)
                naDelete(&n1, _ring);
                n1 = n2

            apow2 = naMult(apow1, a)
            naDelete(&apow1, _ring)
            apow1 = apow2

            naDelete(&coeff, _ring)

        naDelete(&apow1, _ring)
        naDelete(&a, _ring)
    else:
        n1 = naInit(GF2_conv_to_long(GF2X_coeff(rep,0)), _ring)

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
        n1 = naInit(0, _ring)
        a = naPar(1)
        apow1 = naInit(1, _ring)

        for i from 0 <= i <= elem.degree():
            coeff = naInit(int(elem[i]), _ring)

            if not naIsZero(coeff):
                apow2 = naMult(coeff, apow1)
                n2 = naAdd(apow2, n1)
                naDelete(&apow2, _ring)
                naDelete(&n1, _ring);
                n1 = n2

            apow2 = naMult(apow1, a)
            naDelete(&apow1, _ring)
            apow1 = apow2

            naDelete(&coeff, _ring)

        naDelete(&apow1, _ring)
        naDelete(&a, _ring)
    else:
        n1 = naInit(int(elem), _ring)

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
    elem = list(elem)

    if _ring != currRing: rChangeCurrRing(_ring)
    n1 = naInit(0, _ring)
    a = naPar(1)
    apow1 = naInit(1, _ring)

    for i from 0 <= i < len(elem):
        nlCoeff = nlInit2gmp( mpq_numref((<Rational>elem[i]).value), mpq_denref((<Rational>elem[i]).value) )
        naCoeff = naMap00(nlCoeff)
        nlDelete(&nlCoeff, _ring)

        # faster would be to assign the coefficient directly
        apow2 = naMult(naCoeff, apow1)
        n2 = naAdd(apow2, n1)
        naDelete(&apow2, _ring)
        naDelete(&n1, _ring);
        naDelete(&naCoeff, _ring)
        n1 = n2

        apow2 = naMult(apow1, a)
        naDelete(&apow1, _ring)
        apow1 = apow2

    naDelete(&apow1, _ring)
    naDelete(&a, _ring)

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
    cdef number *n = nrzInit(0, _ring)
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
    cdef int _d
    if _ring.ringtype == 1:
        _d = long(d)
        return nr2mMapZp(<number *>_d)
    else:
        lift = d.lift()
        return nrnMapGMP(<number *>((<Integer>lift).value))

cdef object si2sa(number *n, ring *_ring, object base):
    if isinstance(base, FiniteField_prime_modn):
        return base(_ring.cf.n_Int(n, _ring))

    elif isinstance(base, RationalField):
        return si2sa_QQ(n,_ring)

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
        if _ring.ringtype == 0:
            return base(_ring.cf.n_Int(n, _ring))
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
        if _ring.ringtype == 0:
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

cdef int overflow_check(long e, ring *_ring) except -1:
    """
    Raises an ``OverflowError`` if e is > max degree per variable,
    or if it is not acceptable for Singular as exponent of the
    given ring.

    INPUT:

    - ``e`` - some integer representing a degree.
    - ``_ring`` - a pointer to some ring.

    TESTS:

    Whether an overflow occurs or not, partially depends
    on the number of variables in the ring. See :trac:`11856`::

        sage: P.<x,y,z> = QQ[]
        sage: y^2^30
        Traceback (most recent call last):
        ...
        OverflowError: Exponent overflow (1073741824).
        sage: P.<x,y> = QQ[]
        sage: y^2^30
        y^1073741824                                   # 64-bit
        Traceback (most recent call last):             # 32-bit
        ...                                            # 32-bit
        OverflowError: Exponent overflow (1073741824). # 32-bit

        sage: x^2^30*x^2^30
        Traceback (most recent call last):
        ...
        OverflowError: Exponent overflow (2147483648). # 64-bit
        OverflowError: Exponent overflow (1073741824). # 32-bit

    """
    # 2^31 (pPower takes ints)
    if unlikely(e >= _ring.bitmask or e >= 2**31):
        raise OverflowError("Exponent overflow (%d)."%(e))
    return 0

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

    for extension in ["so", "dylib", "dll"]:
        lib = os.environ['SAGE_LOCAL']+"/lib/libsingular."+extension
        if os.path.exists(lib):
            handle = dlopen(lib, RTLD_GLOBAL|RTLD_LAZY)
            if not handle:
                err = dlerror()
                if err:
                    print(err)
            break

    if handle == NULL:
        raise ImportError("cannot load libSINGULAR library")

    # load SINGULAR
    siInit(lib)

    dlclose(handle)

    # we set and save some global Singular options
    singular_options = singular_options | Sy_bit(OPT_REDSB) | Sy_bit(OPT_INTSTRATEGY) | Sy_bit(OPT_REDTAIL) | Sy_bit(OPT_REDTHROUGH)
    global _saved_options
    global _saved_verbose_options
    _saved_options = (int(singular_options), 0, 0)
    _saved_verbose_options = int(singular_verbose_options)

    On(SW_USE_NTL)
    On(SW_USE_NTL_GCD_0)
    On(SW_USE_NTL_GCD_P)
    On(SW_USE_EZGCD)
    Off(SW_USE_NTL_SORT)

    WerrorS_callback = libsingular_error_callback

    error_messages = []

# call the init routine
init_libsingular()

cdef void libsingular_error_callback(const_char_ptr s):
    _s = s
    error_messages.append(_s)
