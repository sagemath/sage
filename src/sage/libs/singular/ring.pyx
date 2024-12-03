"""
Wrapper for Singular's Rings

AUTHORS:

- Martin Albrecht (2009-07): initial implementation

- Kwankyu Lee (2010-06): added matrix term order support

- Miguel Marco (2021): added transcendental extensions over Q
"""
# ****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.cpython.string cimport str_to_bytes, bytes_to_str

from sage.libs.gmp.types cimport __mpz_struct
from sage.libs.gmp.mpz cimport mpz_init_set

from sage.libs.singular.decl cimport ring, currRing
from sage.libs.singular.decl cimport rChangeCurrRing, rComplete, rDelete, idInit
from sage.libs.singular.decl cimport omAlloc0, omStrDup, omAlloc
from sage.libs.singular.decl cimport ringorder_dp, ringorder_Dp, ringorder_lp, ringorder_ip, ringorder_ds, ringorder_Ds, ringorder_ls, ringorder_M, ringorder_c, ringorder_C, ringorder_wp, ringorder_Wp, ringorder_ws, ringorder_Ws, ringorder_a, rRingOrder_t
from sage.libs.singular.decl cimport prCopyR
from sage.libs.singular.decl cimport n_unknown, n_algExt, n_transExt, n_Z, n_Zn,  n_Znm, n_Z2m
from sage.libs.singular.decl cimport n_coeffType
from sage.libs.singular.decl cimport rDefault, GFInfo, ZnmInfo, nInitChar, AlgExtInfo, TransExtInfo


from sage.rings.integer cimport Integer
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.integer_ring import ZZ
import sage.rings.abc
from sage.rings.number_field.number_field_base cimport NumberField
from sage.rings.rational_field import RationalField
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic

from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.polynomial.polynomial_ring import PolynomialRing_field

from sage.rings.fraction_field import FractionField_generic

from cpython.object cimport Py_EQ, Py_NE

from collections import defaultdict


# mapping str --> SINGULAR representation
order_dict = {
    "dp": ringorder_dp,
    "Dp": ringorder_Dp,
    "lp": ringorder_lp,
    "ip": ringorder_ip,
    "ds": ringorder_ds,
    "Ds": ringorder_Ds,
    "ls": ringorder_ls,
    "wp": ringorder_wp,
    "Wp": ringorder_Wp,
    "ws": ringorder_ws,
    "Ws": ringorder_Ws,
    "a":  ringorder_a,
}

cdef extern from "singular/Singular/libsingular.h":
    cdef char * rSimpleOrdStr(rRingOrder_t)

if bytes_to_str(rSimpleOrdStr(ringorder_ip)) == "rp":
    # compatibility for singular 4.3.2p10 and before
    order_dict["rp"] = ringorder_ip
    # also patch term_order mappings
    from sage.rings.polynomial import term_order
    term_order.singular_name_mapping['invlex'] = 'rp'
    term_order.inv_singular_name_mapping['rp'] = 'invlex'

#############################################################################
cdef ring *singular_ring_new(base_ring, n, names, term_order) except NULL:
    r"""
    Create a new Singular ring over the ``base_ring`` in ``n``
    variables with the names ``names`` and the term order
    ``term_order``.

    INPUT:

    - ``base_ring`` -- a Sage ring

    - ``n`` -- the number of variables (> 0)

    - ``names`` -- list of names of length ``n``

    - ``term_order`` -- a term ordering

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: P
        Multivariate Polynomial Ring in x, y, z over Rational Field

        sage: P.term_order()
        Degree reverse lexicographic term order

        sage: P = PolynomialRing(GF(127),3,names='abc', order='lex')
        sage: P
        Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

        sage: P.term_order()
        Lexicographic term order

        sage: z = QQ['z'].0
        sage: K.<s> = NumberField(z^2 - 2)
        sage: P.<x,y> = PolynomialRing(K, 2)

        sage: P.<x,y,z> = ZZ[]; P
        Multivariate Polynomial Ring in x, y, z over Integer Ring

        sage: P.<x,y,z> = Zmod(2^10)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

        sage: P.<x,y,z> = Zmod(3^10)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

        sage: P.<x,y,z> = Zmod(2^100)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1267650600228229401496703205376

        sage: P.<x,y,z> = Zmod(2521352)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352

        sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232

        sage: K = PolynomialRing(QQ, 's,t').fraction_field()
        sage: P.<x,y> = K[]; P
        Multivariate Polynomial Ring in x, y over Fraction Field of Multivariate Polynomial Ring in s, t over Rational Field

        sage: F = PolynomialRing(FiniteField(7),'a,b').fraction_field()
        sage: R.<x,y,z> = F[]
        sage: R
        Multivariate Polynomial Ring in x, y, z over Fraction Field of Multivariate Polynomial Ring in a, b over Finite Field of size 7

    TESTS:

    Check that ``degneglex`` and ``degrevlex`` are the same up to reversal of
    variables (:issue:`29635`)::

        sage: R = PolynomialRing(QQ, 'x', 4, order='degrevlex')
        sage: S = PolynomialRing(QQ, tuple(reversed(R.gens())), order='degneglex')
        sage: L = [v for d in (0..4) for v in IntegerVectors(d, 4)]
        sage: sorted([R.monomial(*e) for e in L]) == sorted([S.monomial(*e) for e in L])
        True

    Check that we are using the libsingular backend instead of the pexpect one::

        sage: F = PolynomialRing(FiniteField(7),'a,b').fraction_field()
        sage: R.<x,y,z> = F[]
        sage: from sage.libs.singular.function import singular_function
        sage: sing_print = singular_function('print')
        sage: print(sing_print(R))
        polynomial ring, over a field, global ordering
        // coefficients: ZZ/7(a, b)
        // number of vars : 3
        //        block   1 : ordering dp
        //                  : names    x y z
        //        block   2 : ordering C

    ::

        sage: F = PolynomialRing(QQ, 's,t').fraction_field()
        sage: R.<x,y,z> = F[]
        sage: from sage.libs.singular.function import singular_function
        sage: print(sing_print(R))
        polynomial ring, over a field, global ordering
        // coefficients: QQ(s, t)
        // number of vars : 3
        //        block   1 : ordering dp
        //                  : names    x y z
        //        block   2 : ordering C

    Small primes::

        sage: R = PolynomialRing(GF(2), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a field, global ordering
        // coefficients: ZZ/2
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C
        sage: R = PolynomialRing(GF(3), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a field, global ordering
        // coefficients: ZZ/3
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C
        sage: R = PolynomialRing(GF(1000000007), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a field, global ordering
        // coefficients: ZZ/1000000007
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Large prime (note that the print is wrong, the field in fact doesn't have zero-divisors)::

        sage: R = PolynomialRing(GF(2^128+51), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a ring (with zero-divisors), global ordering
        // coefficients: ZZ/(bigint(340282366920938463463374607431768211507)^1)
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Finite field with large degree (note that if stack size is too small and the exponent is too large
    a stack overflow may happen inside libsingular)::

        sage: R = PolynomialRing(GF(2^160), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a field, global ordering
        // coefficients: ZZ/2[z160]/(z160^160+z160^159+z160^155+z160^154+z160^153+z160^152+z160^151+z160^149+z160^148+z160^147+z160^146+z160^145+z160^144+z160^143+z160^141+z160^139+z160^137+z160^131+z160^129+z160^128+z160^127+z160^126+z160^123+z160^122+z160^121+z160^117+z160^116+z160^115+z160^113+z160^111+z160^110+z160^108+z160^106+z160^102+z160^100+z160^99+z160^97+z160^96+z160^95+z160^94+z160^93+z160^92+z160^91+z160^87+z160^86+z160^82+z160^80+z160^79+z160^78+z160^74+z160^73+z160^72+z160^71+z160^70+z160^67+z160^66+z160^65+z160^62+z160^59+z160^58+z160^57+z160^55+z160^54+z160^53+z160^52+z160^51+z160^49+z160^47+z160^44+z160^40+z160^35+z160^32+z160^30+z160^28+z160^27+z160^26+z160^24+z160^23+z160^21+z160^20+z160^18+z160^16+z160^11+z160^10+z160^8+z160^7+1)
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Integer modulo small power of 2::

        sage: R = PolynomialRing(Zmod(2^32), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a ring (with zero-divisors), global ordering
        // coefficients: ZZ/(2^32)
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Integer modulo large power of 2::

        sage: R = PolynomialRing(Zmod(2^1000), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a ring (with zero-divisors), global ordering
        // coefficients: ZZ/(bigint(2)^1000)
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Integer modulo large power of odd prime::

        sage: R = PolynomialRing(Zmod(3^300), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a ring (with zero-divisors), global ordering
        // coefficients: ZZ/(bigint(3)^300)
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Integer modulo non-prime::

        sage: R = PolynomialRing(Zmod(15^20), ("a", "b"), implementation="singular"); print(sing_print(R))
        polynomial ring, over a ring (with zero-divisors), global ordering
        // coefficients: ZZ/bigint(332525673007965087890625)
        // number of vars : 2
        //        block   1 : ordering dp
        //                  : names    a b
        //        block   2 : ordering C

    Non-prime finite field with large characteristic (not supported, see :issue:`33319`)::

        sage: PolynomialRing(GF((2^31+11)^2), ("a", "b"), implementation="singular")
        Traceback (most recent call last):
        ...
        TypeError: characteristic must be <= 2147483647.
    """
    cdef long cexponent
    cdef GFInfo* _param
    cdef ZnmInfo _info
    cdef ring* _ring
    cdef char **_names
    cdef char **_ext_names
    cdef int i,j
    cdef int nblcks
    cdef int offset
    cdef int nvars
    cdef int characteristic
    cdef Integer ch, modbase
    cdef int ringorder_column_pos
    cdef int ringorder_column_asc

    cdef int ngens

    cdef n_coeffType ringtype = n_unknown
    cdef MPolynomialRing_libsingular k
    cdef MPolynomial_libsingular minpoly
    cdef AlgExtInfo extParam
    cdef TransExtInfo trextParam
    cdef n_coeffType _type = n_unknown

    #cdef cfInitCharProc myfunctionptr;

    _ring  = NULL

    n = int(n)
    if n < 1:
        raise NotImplementedError(f"polynomials in {n} variables are not supported in Singular")

    nvars = n
    order = TermOrder(term_order, n)

    cdef nbaseblcks = len(order.blocks())
    nblcks = nbaseblcks + order.singular_moreblocks() + 1  # one block for ringorder column
    offset = 0

    if (order._singular_ringorder_column is None or
        order._singular_ringorder_column < 0 or
        order._singular_ringorder_column >= 2*nbaseblcks+2):
        ringorder_column_pos = nbaseblcks
        ringorder_column_type = ringorder_C
    else:
        ringorder_column_pos = order._singular_ringorder_column // 2
        ringorder_column_type = (ringorder_C if order._singular_ringorder_column % 2 == 0
                                             else ringorder_c)

    _names = <char**>omAlloc0(sizeof(char*)*(len(names)))
    for i from 0 <= i < n:
        _name = str_to_bytes(names[i])
        _names[i] = omStrDup(_name)

    # from the SINGULAR source code documentation for the rInit function
    ##  characteristic --------------------------------------------------
    ##  input: 0 ch=0 : Q     parameter=NULL    ffChar=FALSE   float_len (done)
    ##         0    1 : Q(a,...)        *names         FALSE             (done)
    ##         0   -1 : R               NULL           FALSE  0
    ##         0   -1 : R               NULL           FALSE  prec. >6
    ##         0   -1 : C               *names         FALSE  prec. 0..?
    ##         p    p : Fp              NULL           FALSE             (done)
    ##         p   -p : Fp(a)           *names         FALSE             (done)
    ##         q    q : GF(q=p^n)       *names         TRUE              (todo)

    _wvhdl  = <int **>omAlloc0((nblcks + 2) * sizeof(int *))
    _order  = <rRingOrder_t *>omAlloc0((nblcks + 2) * sizeof(int))
    _block0 = <int *>omAlloc0((nblcks + 2) * sizeof(int))
    _block1 = <int *>omAlloc0((nblcks + 2) * sizeof(int))

    cdef int idx = 0
    for i from 0 <= i < nbaseblcks:
        if i == ringorder_column_pos:
            _order[idx] = ringorder_column_type
            idx += 1
        s = order[i].singular_str()
        if s[0] == 'M': # matrix order
            _order[idx] = ringorder_M
            mtx = order[i].matrix().list()
            wv = <int *>omAlloc0(len(mtx)*sizeof(int))
            for j in range(len(mtx)):
                wv[j] = int(mtx[j])
            _wvhdl[idx] = wv
        elif s[0] == 'w' or s[0] == 'W': # weighted degree orders
            _order[idx] = order_dict.get(s[:2], ringorder_dp)
            wts = order[i].weights()
            wv = <int *>omAlloc0(len(wts)*sizeof(int))
            for j in range(len(wts)):
                wv[j] = int(wts[j])
            _wvhdl[idx] = wv
        elif s[0] == '(' and order[i].name() == 'degneglex':  # "(a(1:n),ls(n))"
            _order[idx] = ringorder_a
            if not order[i]:    # may be zero for arbitrary-length orders
                nlen = n
            else:
                nlen = len(order[i])

            _wvhdl[idx] = <int *>omAlloc0(len(order[i])*sizeof(int))
            for j in range(nlen):
                _wvhdl[idx][j] = 1
            _block0[idx] = offset + 1     # same like subsequent ls block
            _block1[idx] = offset + nlen

            idx += 1                      # we need one more block here
            _order[idx] = ringorder_ls

        else: # ordinary orders
            _order[idx] = order_dict.get(s, ringorder_dp)

        _block0[idx] = offset + 1
        if not order[i]: # may be zero in some cases
            _block1[idx] = offset + n
        else:
            _block1[idx] = offset + len(order[i])
        offset = _block1[idx]
        idx += 1

    if ringorder_column_pos >= nbaseblcks:
        _order[idx] = ringorder_column_type

    if isinstance(base_ring, RationalField):
        characteristic = 0
        _ring = rDefault(characteristic, nvars, _names, nblcks,
                         _order, _block0, _block1, _wvhdl)

    elif isinstance(base_ring, FractionField_generic) and isinstance(base_ring.base(), (MPolynomialRing_libsingular, PolynomialRing_field)) and isinstance(base_ring.base().base_ring(), RationalField):
        characteristic = 1
        k = PolynomialRing(RationalField(),
                           names=base_ring.variable_names(), order='lex',
                           implementation='singular')

        ngens = len(k.gens())

        _ext_names = <char**>omAlloc0(ngens*sizeof(char*))
        for i in range(ngens):
            _name = str_to_bytes(k._names[i])
            _ext_names[i] = omStrDup(_name)

        _cfr = rDefault( 0, ngens, _ext_names )
        rComplete(_cfr, 1)

        trextParam.r =  _cfr

        _cf = nInitChar(n_transExt, <void *>&trextParam)

        if (_cf is NULL):
            raise RuntimeError("Failed to allocate _cf ring.")

        _ring = rDefault (_cf, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

    elif isinstance(base_ring, FractionField_generic) and isinstance(base_ring.base(), (MPolynomialRing_libsingular, PolynomialRing_field)) and isinstance(base_ring.base().base_ring(), FiniteField_generic):
        if not base_ring.base_ring().is_prime_field():
            raise NotImplementedError("Transcental extension are not implemented for non-prime finite fields")
        characteristic = int(base_ring.characteristic())
        k = PolynomialRing(base_ring.base_ring(),
            names=base_ring.variable_names(), order='lex', implementation='singular')

        ngens = len(k.gens())

        _ext_names = <char**>omAlloc0(ngens*sizeof(char*))
        for i in range(ngens):
            _name = str_to_bytes(k._names[i])
            _ext_names[i] = omStrDup(_name)

        _cfr = rDefault( characteristic, ngens, _ext_names )
        rComplete(_cfr, 1)

        trextParam.r =  _cfr

        _cf = nInitChar(n_transExt, <void *>&trextParam)

        if (_cf is NULL):
            raise RuntimeError("Failed to allocate _cf ring.")

        _ring = rDefault (_cf, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

    elif isinstance(base_ring, NumberField) and base_ring.is_absolute():
        characteristic = 1
        k = PolynomialRing(RationalField(),
            name=base_ring.variable_name(), order='lex', implementation='singular')

        minpoly = base_ring.polynomial()(k.gen())

        _ext_names = <char**>omAlloc0(sizeof(char*))
        extname = k.gen()
        _name = str_to_bytes(k._names[0])
        _ext_names[0] = omStrDup(_name)
        _cfr = rDefault( 0, 1, _ext_names )

        _cfr.qideal = idInit(1,1)
        rComplete(_cfr, 1)
        _cfr.qideal.m[0] = prCopyR(minpoly._poly, k._ring, _cfr)
        extParam.r =  _cfr

        # _type = nRegister(n_algExt, <cfInitCharProc> naInitChar);
        _cf = nInitChar( n_algExt,  <void *>&extParam) #

        if (_cf is NULL):
            raise RuntimeError("Failed to allocate _cf ring.")

        _ring = rDefault (_cf, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

    elif isinstance(base_ring, IntegerRing_class):
        _cf = nInitChar( n_Z, NULL) # integer coefficient ring
        _ring = rDefault (_cf, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

    elif isinstance(base_ring, sage.rings.abc.IntegerModRing):

        ch = base_ring.characteristic()
        if ch < 2:
            raise NotImplementedError(f"polynomials over {base_ring} are not supported in Singular")

        isprime = ch.is_prime()

        if isprime and ch <= 2147483647:
            assert isinstance(base_ring, FiniteField_generic)
            characteristic = base_ring.characteristic()

            # example for simpler ring creation interface without monomial orderings:
            #_ring = rDefault(characteristic, nvars, _names)

            _ring = rDefault(characteristic, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

        else:
            modbase, cexponent = ch.perfect_power()

            if modbase == 2:
                _cf = nInitChar(n_Z2m, <void *>cexponent)

            elif modbase.is_prime():
                _info.base = <__mpz_struct*>omAlloc(sizeof(__mpz_struct))
                mpz_init_set(_info.base, modbase.value)
                _info.exp = cexponent
                _cf = nInitChar( n_Znm, <void *>&_info )

            else:
                _info.base = <__mpz_struct*>omAlloc(sizeof(__mpz_struct))
                mpz_init_set(_info.base, ch.value)
                _info.exp = 1
                _cf = nInitChar( n_Zn, <void *>&_info )
            _ring = rDefault(_cf, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

    elif isinstance(base_ring, FiniteField_generic):
        assert not base_ring.is_prime_field()  # would have been handled above
        if base_ring.characteristic() > 2147483647:
            raise TypeError("characteristic must be <= 2147483647.")

        # TODO: This is lazy, it should only call Singular stuff not PolynomialRing()
        k = PolynomialRing(base_ring.prime_subfield(),
            name=base_ring.variable_name(), order='lex', implementation='singular')
        minpoly = base_ring.polynomial()(k.gen())

        _ext_names = <char**>omAlloc0(sizeof(char*))
        _name = str_to_bytes(k._names[0])
        _ext_names[0] = omStrDup(_name)
        _cfr = rDefault( <int>base_ring.characteristic(), 1, _ext_names )

        _cfr.qideal = idInit(1,1)
        rComplete(_cfr, 1)
        _cfr.qideal.m[0] = prCopyR(minpoly._poly, k._ring, _cfr)
        extParam.r =  _cfr
        _cf = nInitChar( n_algExt,  <void *>&extParam)

        if (_cf is NULL):
            raise RuntimeError("Failed to allocate _cf ring.")

        _ring = rDefault (_cf, nvars, _names, nblcks, _order, _block0, _block1, _wvhdl)

    else:
        raise NotImplementedError(f"polynomials over {base_ring} are not supported in Singular")

    if _ring is NULL:
        raise ValueError("Failed to allocate Singular ring.")

    _ring.ShortOut = 0

    rChangeCurrRing(_ring)

    wrapped_ring = wrap_ring(_ring)
    if wrapped_ring in ring_refcount_dict:
        raise ValueError('newly created ring already in dictionary??')
    ring_refcount_dict[wrapped_ring] = 1

    rComplete(_ring, 1)

    _ring.ShortOut = 0

    if order.is_local():
        assert(_ring.OrdSgn == -1)
    if order.is_global():
        assert(_ring.OrdSgn == 1)

    return _ring


#############################################################################
ring_refcount_dict = defaultdict(int)


cdef class ring_wrapper_Py():
    r"""
    Python object wrapping the ring pointer.

    This is useful to store ring pointers in Python containers.

    You must not construct instances of this class yourself, use
    :func:`wrap_ring` instead.

    EXAMPLES::

        sage: from sage.libs.singular.ring import ring_wrapper_Py
        sage: ring_wrapper_Py
        <class 'sage.libs.singular.ring.ring_wrapper_Py'>
    """

    cdef ring* _ring

    def __cinit__(self):
        """
        The Cython constructor.

        EXAMPLES::

            sage: from sage.libs.singular.ring import ring_wrapper_Py
            sage: t = ring_wrapper_Py(); t
            The ring pointer 0x0

        These are just wrappers around a pointer, so it isn't really meaningful
        to pickle them::

            sage: TestSuite(t).run(skip='_test_pickling')
        """
        self._ring = NULL

    def __hash__(self):
        """
        Return a hash value so that instances can be used as dictionary keys.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.libs.singular.ring import ring_wrapper_Py
            sage: t = ring_wrapper_Py()
            sage: t.__hash__()
            0
        """
        return <long>(self._ring)

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.libs.singular.ring import ring_wrapper_Py
            sage: t = ring_wrapper_Py()
            sage: t
            The ring pointer 0x0
            sage: t.__repr__()
            'The ring pointer 0x0'
        """
        return 'The ring pointer '+hex(self.__hash__())

    # This could be written using __eq__ but that does not work
    # due to https://github.com/cython/cython/issues/2019
    def __richcmp__(ring_wrapper_Py self, other, int op):
        """
        Equality comparison between two ``ring_wrapper_Py`` instances,
        for use when hashing.

        INPUT:

        - ``right`` -- a :class:`ring_wrapper_Py`

        OUTPUT: ``True`` if both ``ring_wrapper_Py`` wrap the same pointer

        EXAMPLES::

            sage: from sage.libs.singular.ring import (ring_wrapper_Py,
            ....:     currRing_wrapper)
            sage: t = ring_wrapper_Py()
            sage: t == t
            True
            sage: P.<x,y,z> = QQ[]
            sage: t2 = currRing_wrapper()
            sage: t3 = currRing_wrapper()
            sage: t == t2
            False
            sage: t2 == t3
            True
            sage: t2 != t3
            False
            sage: t2 == None
            False
        """
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented

        if type(other) is not ring_wrapper_Py:
            return op != Py_EQ

        r = <ring_wrapper_Py>other
        return (self._ring == r._ring) == (op == Py_EQ)


cdef wrap_ring(ring* R):
    """
    Wrap a C ring pointer into a Python object.

    INPUT:

    - ``R`` -- a singular ring (a C datastructure)

    OUTPUT:

    A Python object :class:`ring_wrapper_Py` wrapping the C pointer.
    """
    cdef ring_wrapper_Py W = ring_wrapper_Py()
    W._ring = R
    return W


cdef ring *singular_ring_reference(ring *existing_ring) except NULL:
    """
    Refcount the ring ``existing_ring``.

    INPUT:

    - ``existing_ring`` -- a Singular ring

    OUTPUT:

    The same ring with its refcount increased. If ``existing_ring``
    has not been refcounted yet, it will be after calling this function.
    If initially ``existing_ring`` was refcounted once, then after
    calling this function `n` times, you need to call :func:`singular_ring_delete`
    `n+1` times to actually deallocate the ring.

    EXAMPLES::

        sage: import gc
        sage: _ = gc.collect()
        sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
        sage: from sage.libs.singular.groebner_strategy import GroebnerStrategy
        sage: from sage.libs.singular.ring import ring_refcount_dict
        sage: n = len(ring_refcount_dict)
        sage: prev_rings = set(ring_refcount_dict)
        sage: P = MPolynomialRing_libsingular(GF(541), 2, ('x', 'y'), TermOrder('degrevlex', 2))
        sage: ring_ptr = set(ring_refcount_dict).difference(prev_rings).pop()
        sage: ring_ptr  # random output
        The ring pointer 0x7f78a646b8d0
        sage: ring_refcount_dict[ring_ptr]
        3

        sage: strat = GroebnerStrategy(Ideal([P.gen(0) + P.gen(1)]))
        sage: ring_refcount_dict[ring_ptr]
        6

        sage: del strat
        sage: _ = gc.collect()
        sage: ring_refcount_dict[ring_ptr]
        4

        sage: del P
        sage: _ = gc.collect()
        sage: ring_ptr in ring_refcount_dict
        True
    """
    if existing_ring is NULL:
        raise ValueError('singular_ring_reference(ring*) called with NULL pointer.')

    cdef object r = wrap_ring(existing_ring)
    ring_refcount_dict[r] += 1
    return existing_ring


#############################################################################
cdef void singular_ring_delete(ring *doomed) noexcept:
    """
    Carefully deallocate the ring, without changing "currRing" (since
    this method can be called at unpredictable times due to garbage
    collection).

    TESTS:

    This example caused a segmentation fault with a previous version
    of this method::

        sage: import gc
        sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
        sage: R1 = MPolynomialRing_libsingular(GF(5), 2, ('x', 'y'), TermOrder('degrevlex', 2))
        sage: R2 = MPolynomialRing_libsingular(GF(11), 2, ('x', 'y'), TermOrder('degrevlex', 2))
        sage: R3 = MPolynomialRing_libsingular(GF(13), 2, ('x', 'y'), TermOrder('degrevlex', 2))
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
    if doomed is NULL:
        # When this is called with a NULL pointer, we do nothing.
        # This is analogous to the libc function free().
        return

    if not ring_refcount_dict:  # arbitrary finalization order when we shut Sage down
        return

    cdef ring_wrapper_Py r = wrap_ring(doomed)
    ring_refcount_dict[r] -= 1
    if ring_refcount_dict[r] > 0:
        return

    del ring_refcount_dict[r]

    global currRing
    cdef ring *oldRing = currRing
    if currRing == doomed:
        rDelete(doomed)
        currRing = <ring*>NULL
    else:
        rChangeCurrRing(doomed)
        rDelete(doomed)
        rChangeCurrRing(oldRing)


#############################################################################
# helpers for debugging

cpdef poison_currRing(frame, event, arg):
    """
    Poison the ``currRing`` pointer.

    This function sets the ``currRing`` to an illegal value. By
    setting it as the python debug hook, you can poison the currRing
    before every evaluated Python command (but not within Cython
    code).

    INPUT:

    - ``frame``, ``event``, ``arg`` -- the standard arguments for the
      CPython debugger hook; they are not used

    OUTPUT:

    Returns itself, which ensures that :func:`poison_currRing` will
    stay in the debugger hook.

    EXAMPLES::

        sage: previous_trace_func = sys.gettrace()   # None if no debugger running
        sage: from sage.libs.singular.ring import poison_currRing
        sage: sys.settrace(poison_currRing)
        sage: sys.gettrace()
        <built-in function poison_currRing>
        sage: sys.settrace(previous_trace_func)  # switch it off again
    """
    global currRing
    currRing = <ring*>NULL
    return poison_currRing


cpdef print_currRing():
    """
    Print the ``currRing`` pointer.

    EXAMPLES::

        sage: from sage.libs.singular.ring import print_currRing
        sage: print_currRing()   # random output
        DEBUG: currRing == 0x7fc6fa6ec480

        sage: from sage.libs.singular.ring import poison_currRing
        sage: _ = poison_currRing(None, None, None)
        sage: print_currRing()
        DEBUG: currRing == 0x0
    """
    cdef size_t addr = <size_t>currRing
    print("DEBUG: currRing == " + str(hex(addr)))


def currRing_wrapper():
    """
    Return a wrapper for the current ring, for use in debugging
    ``ring_refcount_dict``.

    EXAMPLES::

        sage: from sage.libs.singular.ring import currRing_wrapper
        sage: currRing_wrapper()
        The ring pointer ...
    """
    return wrap_ring(currRing)
