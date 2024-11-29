"""
Direct low-level access to SINGULAR's Groebner basis engine via libSINGULAR

AUTHOR:

- Martin Albrecht (2007-08-08): initial version

EXAMPLES::

    sage: x,y,z = QQ['x,y,z'].gens()
    sage: I = ideal(x^5 + y^4 + z^3 - 1,  x^3 + y^3 + z^2 - 1)
    sage: I.groebner_basis('libsingular:std')
    [y^6 + x*y^4 + 2*y^3*z^2 + x*z^3 + z^4 - 2*y^3 - 2*z^2 - x + 1,
     x^2*y^3 - y^4 + x^2*z^2 - z^3 - x^2 + 1, x^3 + y^3 + z^2 - 1]

We compute a Groebner basis for cyclic 6, which is a standard
benchmark and test ideal::

    sage: R.<x,y,z,t,u,v> = QQ['x,y,z,t,u,v']
    sage: I = sage.rings.ideal.Cyclic(R,6)
    sage: B = I.groebner_basis('libsingular:std')
    sage: len(B)
    45

Two examples from the Mathematica documentation (done in Sage):

- We compute a Groebner basis::

        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: ideal(x^2 - 2*y^2, x*y - 3).groebner_basis('libsingular:slimgb')
        [x - 2/3*y^3, y^4 - 9/2]

- We show that three polynomials have no common root::

        sage: R.<x,y> = QQ[]
        sage: ideal(x+y, x^2 - 1, y^2 - 2*x).groebner_basis('libsingular:slimgb')
        [1]
"""

#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.singular.decl cimport tHomog, number, IDELEMS, p_Copy, rChangeCurrRing
from sage.libs.singular.decl cimport idInit, id_Delete, currRing, Sy_bit, OPT_REDSB
from sage.libs.singular.decl cimport scKBase, poly, testHomog, idSkipZeroes, id_RankFreeModule, kStd
from sage.libs.singular.decl cimport OPT_REDTAIL, singular_options, kInterRed, t_rep_gb, p_GetCoeff
from sage.libs.singular.decl cimport pp_Mult_nn, p_Delete, n_Delete
from sage.libs.singular.decl cimport rIsPluralRing
from sage.libs.singular.decl cimport n_Z, n_Zn, n_Znm, n_Z2m

from sage.rings.polynomial.multi_polynomial_libsingular cimport new_MP
from sage.rings.polynomial.plural cimport new_NCP

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular
from sage.structure.sequence import Sequence

from sage.rings.polynomial.plural cimport NCPolynomialRing_plural, NCPolynomial_plural

cdef object singular_ideal_to_sage_sequence(ideal *i, ring *r, object parent):
    """
    Convert a SINGULAR ideal to a Sage Sequence (the format Sage
    stores a Groebner basis in).

    INPUT:

    - ``i`` -- a SINGULAR ideal
    - ``r`` -- a SINGULAR ring
    - ``parent`` -- a Sage ring matching r
    """
    cdef int j
    cdef MPolynomial_libsingular p
    cdef NCPolynomial_plural p_nc

    l = []

    if rIsPluralRing(r):
        for j from 0 <= j < IDELEMS(i):
            p_nc = new_NCP(parent, p_Copy(i.m[j], r))
            l.append( p_nc )
    else:
        for j from 0 <= j < IDELEMS(i):
            p = new_MP(parent, p_Copy(i.m[j], r))
            l.append( p )

    return Sequence(l, check=False, immutable=True)

cdef ideal *sage_ideal_to_singular_ideal(I) except NULL:
    """
    convert a Sage ideal to a SINGULAR ideal

    INPUT:

    - ``I`` -- a Sage ideal in a ring of type
      :class:`~sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular`
      or a list of generators

    TESTS:


    We test conversion::

        sage: P.<x,y,z> = QQ[]
        sage: sage.libs.singular.function_factory.ff.std(Sequence([x,y,z]))
        [z, y, x]
        sage: sage.libs.singular.function_factory.ff.std(Ideal([x,y,z]))
        [z, y, x]
    """
    R = I.ring()
    try:
        gens = I.gens()
    except AttributeError:
        gens = I
    cdef ideal *result
    cdef ring *r
    cdef ideal *i
    cdef int j = 0

    if isinstance(R, MPolynomialRing_libsingular):
        r = (<MPolynomialRing_libsingular>R)._ring
    elif isinstance(R, NCPolynomialRing_plural):
        r = (<NCPolynomialRing_plural>R)._ring
    else:
        raise TypeError("Ring must be of type 'MPolynomialRing_libsingular'")

    rChangeCurrRing(r)

    i = idInit(len(gens),1)
    for j,f in enumerate(gens):
        if isinstance(f, MPolynomial_libsingular):
            i.m[j] = p_Copy((<MPolynomial_libsingular>f)._poly, r)
        elif isinstance(f, NCPolynomial_plural):
            i.m[j] = p_Copy((<NCPolynomial_plural>f)._poly, r)
        else:
            id_Delete(&i, r)
            raise TypeError("All generators must be of type MPolynomial_libsingular.")
    return i


def kbase_libsingular(I, degree=None):
    """
    SINGULAR's ``kbase()`` algorithm.

    INPUT:

    - ``I`` -- a groebner basis of an ideal

    - ``degree`` -- integer (default: ``None``); if not ``None``, return
      only the monomials of the given degree

    OUTPUT:

    Computes a vector space basis (consisting of monomials) of the quotient
    ring by the ideal, resp. of a free module by the module, in case it is
    finite dimensional and if the input is a standard basis with respect to
    the ring ordering. If the input is not a standard basis, the leading terms
    of the input are used and the result may have no meaning.

    With two arguments: computes the part of a vector space basis of the
    respective quotient with degree of the monomials equal to the second
    argument. Here, the quotient does not need to be finite dimensional.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: I = R.ideal(x^2-2*y^2, x*y-3)
        sage: I.normal_basis()  # indirect doctest
        [y^3, y^2, y, 1]
        sage: J = R.ideal(x^2-2*y^2)
        sage: [J.normal_basis(d) for d in (0..4)]  # indirect doctest
        [[1], [y, x], [y^2, x*y], [y^3, x*y^2], [y^4, x*y^3]]
    """

    global singular_options

    cdef ideal *i = sage_ideal_to_singular_ideal(I)
    cdef ring *r = currRing
    cdef ideal *q = currRing.qideal

    cdef ideal *result
    singular_options = singular_options | Sy_bit(OPT_REDSB)

    degree = -1 if degree is None else int(degree)

    sig_on()
    result = scKBase(degree, i, q)
    sig_off()

    id_Delete(&i, r)
    res = singular_ideal_to_sage_sequence(result,r,I.ring())

    id_Delete(&result, r)

    return res


def std_libsingular(I):
    """
    SINGULAR's ``std()`` algorithm.

    INPUT:

    - ``I`` -- a Sage ideal
    """
    global singular_options

    cdef ideal *i = sage_ideal_to_singular_ideal(I)
    cdef ring *r = currRing
    cdef tHomog hom = testHomog
    cdef ideal *result

    singular_options = singular_options | Sy_bit(OPT_REDSB)

    sig_on()
    result =kStd(i,NULL,hom,NULL)
    sig_off()

    idSkipZeroes(result)

    id_Delete(&i,r)

    res = singular_ideal_to_sage_sequence(result,r,I.ring())

    id_Delete(&result,r)
    return res


def slimgb_libsingular(I):
    """
    SINGULAR's ``slimgb()`` algorithm.

    INPUT:

    - ``I`` -- a Sage ideal
    """
    global singular_options

    cdef tHomog hom=testHomog
    cdef ideal *i
    cdef ring *r
    cdef ideal *result

    i = sage_ideal_to_singular_ideal(I)
    r = currRing

    if r.OrdSgn != 1:
        id_Delete(&i, r)
        raise TypeError("ordering must be global for slimgb")

    if i.rank < id_RankFreeModule(i, r):
        id_Delete(&i, r)
        raise TypeError

    singular_options = singular_options | Sy_bit(OPT_REDSB)

    sig_on()
    result = t_rep_gb(r, i, i.rank, 0)
    sig_off()

    id_Delete(&i, r)

    res = singular_ideal_to_sage_sequence(result, r, I.ring())

    id_Delete(&result, r)
    return res


def interred_libsingular(I):
    """
    SINGULAR's ``interred()`` command.

    INPUT:

    - ``I`` -- a Sage ideal

    EXAMPLES::

        sage: P.<x,y,z> = PolynomialRing(ZZ)
        sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
        sage: I.interreduced_basis()
        [y*z^2 - 81*x*y - 9*y - z, z^3 - x, x^2 - 3*y, 9*y^2 - y*z + 1]

        sage: P.<x,y,z> = PolynomialRing(QQ)
        sage: I = ideal( x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1 )
        sage: I.interreduced_basis()
        [y*z^2 - 81*x*y - 9*y - z, z^3 - x, x^2 - 3*y, y^2 - 1/9*y*z + 1/9]
    """
    global singular_options

    cdef ideal *i
    cdef ideal *result
    cdef ring *r
    cdef number *n
    cdef poly *p
    cdef int j
    cdef int bck

    try:
        if len(I.gens()) == 0:
            return Sequence([], check=False, immutable=True)
    except AttributeError:
        pass

    i = sage_ideal_to_singular_ideal(I)
    r = currRing

    bck = singular_options
    singular_options = singular_options | Sy_bit(OPT_REDTAIL)|Sy_bit(OPT_REDSB)
    sig_on()
    result = kInterRed(i,NULL)
    sig_off()
    singular_options = bck

    # divide head by coefficients
    if r.cf.type != n_Z and r.cf.type != n_Znm and r.cf.type != n_Zn and r.cf.type != n_Z2m:
        for j from 0 <= j < IDELEMS(result):
            p = result.m[j]
            if p:
                n = p_GetCoeff(p, r)
                n = r.cf.cfInvers(n, r.cf)
            result.m[j] = pp_Mult_nn(p, n, r)
            p_Delete(&p, r)
            n_Delete(&n, r.cf)

    id_Delete(&i, r)

    res = sorted(singular_ideal_to_sage_sequence(result,r,I.ring()),reverse=True)

    id_Delete(&result, r)
    return res
