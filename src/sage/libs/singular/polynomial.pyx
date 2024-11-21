"""
Wrapper for Singular's Polynomial Arithmetic

AUTHOR:

- Martin Albrecht (2009-07): refactoring
"""
# ****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

cdef extern from *: # hack to get at cython macro
    int unlikely(int)

import re
plusminus_pattern = re.compile(r"([^\(^])([\+\-])")
parenthvar_pattern = re.compile(r"\(([a-zA-Z][a-zA-Z0-9]*)\)")

from sage.cpython.string cimport bytes_to_str, str_to_bytes

from sage.libs.singular.decl cimport number, ideal
from sage.libs.singular.decl cimport currRing, rChangeCurrRing
from sage.libs.singular.decl cimport p_Copy, p_Add_q, p_Neg, pp_Mult_nn, p_GetCoeff, p_IsConstant, p_Cmp, pNext
from sage.libs.singular.decl cimport p_GetMaxExp, pp_Mult_qq, pPower, p_String, p_GetExp, p_LDeg
from sage.libs.singular.decl cimport n_Delete, idInit, fast_map_common_subexp, id_Delete
from sage.libs.singular.decl cimport omAlloc0, omStrDup, omFree
from sage.libs.singular.decl cimport p_GetComp
from sage.libs.singular.decl cimport pSubst
from sage.libs.singular.decl cimport p_Normalize
from sage.libs.singular.decl cimport ndCopyMap, maMapPoly

from sage.libs.singular.singular cimport sa2si, si2sa, overflow_check

from sage.misc.latex import latex

cdef int singular_polynomial_check(poly *p, ring *r) except -1:
    """
    Run consistency checks on ``p``.
    """
    while p:
        if p_GetCoeff(p, r) == NULL:
            raise ZeroDivisionError("NULL pointer as coefficient.")
        p = p.next
    return 0

cdef int singular_polynomial_add(poly **ret, poly *p, poly *q, ring *r) noexcept:
    """
    ``ret[0] = p+q`` where ``p`` and ``p`` in ``r``.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``q`` -- a Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: x + y # indirect doctest
        x + y

        sage: x + P(0)
        x
    """
    if r != currRing:
        rChangeCurrRing(r)
    p = p_Copy(p, r)
    q = p_Copy(q, r)
    ret[0] = p_Add_q(p, q, r)
    return 0

cdef int singular_polynomial_sub(poly **ret, poly *p, poly *q, ring *r) noexcept:
    """
    ``ret[0] = p-q`` where ``p`` and ``p`` in ``r``.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``q`` -- a Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: x - y # indirect doctest
        x - y

        sage: x + P(0)
        x
    """
    if r != currRing:
        rChangeCurrRing(r)
    p = p_Copy(p, r)
    q = p_Copy(q, r)
    ret[0] = p_Add_q(p, p_Neg(q, r), r)
    return 0

cdef int singular_polynomial_rmul(poly **ret, poly *p, RingElement n, ring *r) noexcept:
    """
    ``ret[0] = n*p`` where ``n`` is a coefficient and ``p`` in ``r``.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``n`` -- a Sage coefficient
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: 2*x # indirect doctest
        2*x

        sage: P(0)*x
        0
    """
    if r != currRing:
        rChangeCurrRing(r)
    cdef number *_n = sa2si(n, r)
    ret[0] = pp_Mult_nn(p, _n, r)
    n_Delete(&_n, r.cf)
    return 0

cdef int singular_polynomial_call(poly **ret, poly *p, ring *r, list args,
                                  poly *(*get_element)(object) noexcept) noexcept:
    """
    ``ret[0] = p(*args)`` where each entry in arg  is a polynomial and ``p`` in ``r``.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``r`` -- a Singular ring
    - ``args`` -- list/tuple of elements which can be converted to Singular
      polynomials using the ``(get_element)`` function provided
    - ``(*get_element)`` -- a function to turn a Sage element into a
      Singular element

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: x(0,0,0) # indirect doctest
        0

        sage: (3*x*z)(x,x,x)
        3*x^2

    TESTS:

    Test that there is no memory leak in evaluating polynomials. Note
    that (lib)Singular has pre-allocated buckets, so we have to run a
    lot of iterations to fill those up first::

        sage: import resource
        sage: import gc
        sage: F.<a> = GF(7^2)
        sage: R.<x,y> = F[]
        sage: p = x+2*y
        sage: def leak(N):
        ....:     before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        ....:     gc.collect()
        ....:     for i in range(N):
        ....:         _ = p(a, a)
        ....:     after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        ....:     return (after - before) * 1024   # ru_maxrss is in kilobytes

    Loop (at most 30 times) until we have 6 consecutive zeros when
    calling ``leak(10000)``. Depending on the operating system, it is
    possible to have several nonzero leak values in the beginning, but
    after a while we should get only zeros. The fact that we require 6
    zeros also means that Singular's pre-allocated buckets should not
    be sufficient if there really would be a memory leak. ::

        sage: zeros = 0
        sage: for i in range(30):  # long time
        ....:     n = leak(10000)
        ....:     print("Leaked {} bytes".format(n))
        ....:     if n == 0:
        ....:         zeros += 1
        ....:         if zeros >= 6:
        ....:             break
        ....:     else:
        ....:         zeros = 0
        Leaked...
        Leaked 0 bytes
        Leaked 0 bytes
        Leaked 0 bytes
        Leaked 0 bytes
        Leaked 0 bytes
    """
    cdef long l = len(args)
    cdef ideal *to_id = idInit(l, 1)
    cdef bint constant_args = 1
    for i from 0 <= i < l:
        to_id.m[i]= p_Copy( get_element(args[i]), r)
        if not p_IsConstant(to_id.m[i], r):
            constant_args = 0

    cdef ideal *from_id
    rChangeCurrRing(r)
    cdef ideal *res_id
    if not constant_args:
        from_id = idInit(1,1)
        from_id.m[0] = p

        res_id = fast_map_common_subexp(from_id, r, to_id, r)
        ret[0] = res_id.m[0]

        from_id.m[0] = NULL
        res_id.m[0] = NULL
        id_Delete(&from_id, r)
        id_Delete(&res_id, r)
    else:
        ret[0] = maMapPoly(p, r, to_id, r, ndCopyMap)

    # Unsure why we have to normalize here. See #16958
    p_Normalize(ret[0], r)

    id_Delete(&to_id, r)

    return 0

cdef int singular_polynomial_cmp(poly *p, poly *q, ring *r) noexcept:
    """
    Compare two Singular elements ``p`` and ``q`` in ``r``.

    INPUT:

    - ``p`` -- a Singular polynomial
    - ``q`` -- a Singular polynomial
    - ``r`` -- a Singular ring

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

    TESTS:

    Ensure comparison is transitive::

        sage: P.<x, y> = GF(5)[]
        sage: P(1)>P(4)
        False
        sage: P(4)>P(2)
        True
        sage: P(2)>P(1)
        True

    More tests for :issue:`35681`::

        sage: P.<x, y> = ZZ[]
        sage: P(0) > P(-1)
        True
    """
    cdef number *h
    cdef int ret

    if r != currRing:
        rChangeCurrRing(r)

    while p!=NULL and q!=NULL:
        ret = p_Cmp( p, q, r)
        if ret != 0:
            return ret

        # compare coeffs
        ret = (-1+r.cf.cfEqual(p_GetCoeff(p, r),p_GetCoeff(q, r),r.cf)
               +2*r.cf.cfGreater(p_GetCoeff(p, r),p_GetCoeff(q, r),r.cf)) # -1: <, 0:==, 1: >
        if ret != 0:
            return ret

        p = pNext(p)
        q = pNext(q)

    if p==NULL and q != NULL:
        return 1-2*r.cf.cfGreaterZero(p_GetCoeff(q,r), r.cf) # -1: <, 1: >
    elif p!=NULL and q==NULL:
        return -1+2*r.cf.cfGreaterZero(p_GetCoeff(p,r), r.cf)

    return 0

cdef int singular_polynomial_mul(poly** ret, poly *p, poly *q, ring *r) except -1:
    """
    ``ret[0] = p*q`` where ``p`` and ``p`` in ``r``.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``q`` -- a Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: x*y # indirect doctest
        x*y

        sage: x * P(0)
        0
    """
    if r != currRing:
        rChangeCurrRing(r)
    cdef unsigned long le = p_GetMaxExp(p, r)
    cdef unsigned long lr = p_GetMaxExp(q, r)
    cdef unsigned long esum = le + lr
    overflow_check(esum, r)
    ret[0] = pp_Mult_qq(p, q, r)
    return 0

cdef int singular_polynomial_div_coeff(poly** ret, poly *p, poly *q, ring *r) except -1:
    """
    ``ret[0] = p/n`` where ``p`` and ``q`` in ``r`` and ``q`` constant.

    The last condition is not checked.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``q`` -- a constant Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: x/2 # indirect doctest
        1/2*x

        sage: x/0
        Traceback (most recent call last):
        ...
        ZeroDivisionError: rational division by zero
    """
    if q == NULL:
        raise ZeroDivisionError
    sig_on()
    cdef number *n = p_GetCoeff(q, r)
    n = r.cf.cfInvers(n,r.cf)
    ret[0] = pp_Mult_nn(p, n, r)
    n_Delete(&n, r.cf)
    sig_off()
    return 0

cdef int singular_polynomial_pow(poly **ret, poly *p, unsigned long exp, ring *r) except -1:
    """
    ``ret[0] = p**exp`` where ``p`` in ``r`` and ``exp`` > 0.

    The last condition is not checked.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``exp`` -- integer
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: f = 3*x*y + 5/2*z
        sage: f*f == f^2 # indirect doctest
        True
        sage: f^2
        9*x^2*y^2 + 15*x*y*z + 25/4*z^2
        sage: f^0
        1
        sage: f^(2^60)
        Traceback (most recent call last):
        ...
        OverflowError: ...
    """
    cdef unsigned long v = p_GetMaxExp(p, r)
    v = v * exp
    overflow_check(v, r)

    if r != currRing:
        rChangeCurrRing(r)
    cdef int count = singular_polynomial_length_bounded(p, 15)
    if count >= 15 or exp > 15:
        sig_on()
    ret[0] = pPower( p_Copy(p,r), exp)
    if count >= 15 or exp > 15:
        sig_off()
    return 0

cdef int singular_polynomial_neg(poly **ret, poly *p, ring *r) noexcept:
    """
    ``ret[0] = -p where ``p`` in ``r``.

    The last condition is not checked.

    INPUT:

    - ``ret`` -- a pointer to a Singular polynomial to store the result in
    - ``p`` -- a Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: f = 3*x*y + 5/2*z
        sage: -f # indirect doctest
        -3*x*y - 5/2*z
        sage: -P(0)
        0
    """
    if r != currRing:
        rChangeCurrRing(r)
    ret[0] = p_Neg(p_Copy(p,r),r)
    return 0

cdef object singular_polynomial_str(poly *p, ring *r):
    """
    Return the string representation of ``p``.

    INPUT:

    - ``p`` -- a Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = ZZ[]
        sage: str(x) # indirect doctest
        'x'
        sage: str(10*x)
        '10*x'
    """
    if r != currRing:
        rChangeCurrRing(r)

    s = bytes_to_str(p_String(p, r, r))
    s = plusminus_pattern.sub("\\1 \\2 ", s)
    s = parenthvar_pattern.sub("\\1", s)
    return s

cdef object singular_polynomial_latex(poly *p, ring *r, object base, object latex_gens):
    r"""
    Return the LaTeX string representation of ``p``.

    INPUT:

    - ``p`` -- a Singular polynomial
    - ``r`` -- a Singular ring

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: latex(x) # indirect doctest
        x
        sage: latex(10*x^2 + 1/2*y)
        10 x^{2} + \frac{1}{2} y

    Demonstrate that coefficients over non-atomic represented rings are
    properly parenthesized (:issue:`11186`)::

        sage: x = var('x')
        sage: K.<z> = QQ.extension(x^2 + x + 1)
        sage: P.<v,w> = K[]
        sage: latex((z+1)*v*w - z*w^2 + z*v + z^2*w + z + 1)
        \left(z + 1\right) v w - z w^{2} + z v + \left(-z - 1\right) w + z + 1

    Demonstrate that there are no extra blanks in latex expression of multivariate
    polynomial (:issue:`12908`)::

        sage: R.<X,Y> = ZZ[]
        sage: latex(X-Y)
        X - Y
        sage: latex(X^2-X)
        X^{2} - X
        sage: latex(-Y^2-Y)
        -Y^{2} - Y
    """
    poly = ""
    cdef unsigned long e
    cdef int n = r.N, j
    cdef int atomic_repr = base._repr_option('element_is_atomic')
    while p:

        # First determine the multinomial:
        multi = ""
        for j in range(1, n+1):
            e = p_GetExp(p, j, r)
            if e > 0:
                multi += " " + latex_gens[j-1]
            if e > 1:
                multi += "^{%d}" % e
        multi = multi.lstrip().rstrip()

        # Next determine coefficient of multinomial
        c =  si2sa(p_GetCoeff(p, r), r, base)
        if not multi:
            multi = latex(c)
        elif c != 1:
            if c == -1:
                multi = "-%s" % (multi)
            else:
                sc = latex(c)
                # Add parenthesis if the coefficient consists of terms divided by +, -
                # (starting with - is not enough) and is not the constant term
                if not atomic_repr and multi and (sc.find("+") != -1 or sc[1:].find("-") != -1):
                    sc = "\\left(%s\\right)" % sc
                multi = "%s %s" % (sc, multi)

        # Now add on coefficiented multinomials
        if poly:
            poly = poly + " + "
        poly = poly + multi

        p = pNext(p)

    poly = poly.lstrip().rstrip()
    poly = poly.replace("+ -","- ")

    if not poly:
        return "0"
    return poly

cdef object singular_polynomial_str_with_changed_varnames(poly *p, ring *r, object varnames):
    cdef char **_names
    cdef char **_orig_names
    cdef int i

    if len(varnames) != r.N:
        raise TypeError("len(varnames) doesn't equal self.parent().ngens()")

    _names = <char**>omAlloc0(sizeof(char*)*r.N)
    for i from 0 <= i < r.N:
        _name = str_to_bytes(varnames[i])
        _names[i] = omStrDup(_name)

    _orig_names = r.names
    r.names = _names
    s = singular_polynomial_str(p, r)
    r.names = _orig_names

    for i from 0 <= i < r.N:
        omFree(_names[i])
    omFree(_names)
    return s

cdef long singular_polynomial_deg(poly *p, poly *x, ring *r) noexcept:
    cdef long _deg, deg
    cdef int dummy

    deg = -1
    _deg = -1
    if p == NULL:
        return -1
    if r != currRing:
        rChangeCurrRing(r)
    if x == NULL:
        return p_LDeg(p, &dummy, r)

    cdef int i = 0
    for i in range(1,r.N+1):
        if p_GetExp(x, i, r):
            break
    while p:
        _deg =  p_GetExp(p,i,r)
        if _deg > deg:
            deg = _deg
        p = pNext(p)
    return deg

cdef int singular_polynomial_length_bounded(poly *p, int bound) noexcept:
    """
    Return the number of monomials in ``p`` but stop counting at
    ``bound``.

    This is useful to estimate whether a certain calculation will take
    long or not.

    INPUT:

    - ``p`` -- a Singular polynomial
    - ``bound`` -- integer > 0
    """
    cdef int count = 0
    while p != NULL and count < bound:
        p = pNext(p)
        count += 1
    return count

cdef int singular_vector_maximal_component(poly *v, ring *r) except -1:
    """
    returns the maximal module component of the vector ``v``.
    INPUT:

    - ``v`` -- a polynomial/vector
    - ``r`` -- a ring
    """
    cdef int res=0
    while v!=NULL:
        res=max(p_GetComp(v, r), res)
        v = pNext(v)
    return res

cdef int singular_polynomial_subst(poly **p, int var_index, poly *value, ring *r) except -1:
    """
    Substitute variable ``var_index`` with ``value`` in ``p``.

    INPUT:

    - ``p`` -- a polynomial
    - ``var_index`` -- integer < ngens (zero based indexing)
    - ``value`` -- a polynomial
    - ``r`` -- a ring
    """

    if r != currRing:
        rChangeCurrRing(r)

    if p_IsConstant(value, r):
        p[0] = pSubst(p[0], var_index+1, value)
        return 0

    cdef unsigned long exp = p_GetExp(p[0], var_index+1, r) * p_GetMaxExp(value, r)

    overflow_check(exp, r)

    cdef int count = singular_polynomial_length_bounded(p[0], 15)
    if unlikely(count >= 15 or exp > 15): sig_on()
    p[0] = pSubst(p[0], var_index+1, value)
    if unlikely(count >= 15 or exp > 15): sig_off()
    return 0
