# sage_setup: distribution = sagemath-flint
r"""
Compute Hilbert series of monomial ideals

This implementation was provided at :issue:`26243` and is supposed to be a way
out when Singular fails with an int overflow, which will regularly be the case
in any example with more than 34 variables.
"""
#*****************************************************************************
#
#    Tools to compute Hilbert Poincaré series of monomial ideals
#
#    Copyright (C) 2018 Simon A. King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

import sage.interfaces.abc

from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

from cysignals.memory cimport sig_malloc
from cpython.list cimport PyList_GET_ITEM

# A boilerplate class to implement a binary tree
# with some additional information

from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.types cimport fmpz_poly_t

cdef class Node:
    """
    A node of a binary tree

    It has slots for data that allow to recursively compute
    the first Hilbert series of a monomial ideal.
    """
    cdef Node Back, Left, Right
    cdef list Id
    cdef fmpz_poly_t LMult
    cdef fmpz_poly_t RMult
    cdef fmpz_poly_t LeftFHS

    def __cinit__(self):
        fmpz_poly_init(self.LMult)
        fmpz_poly_init(self.RMult)
        fmpz_poly_init(self.LeftFHS)

    def __dealloc__(self):
        fmpz_poly_clear(self.LMult)
        fmpz_poly_clear(self.RMult)
        fmpz_poly_clear(self.LeftFHS)

cdef inline size_t median(list v) noexcept:
    """
    Specialized version of :func:`from sage.stats.basic_stats.median`.
    """
    cdef list values = sorted(v)
    cdef size_t nvals = len(values)
    if nvals % 2 == 1:
        return values[(nvals+1)//2 - 1]
    else:
        lower = values[(nvals+1)//2 - 1]
        upper = values[nvals // 2]
        return (lower + upper) // 2

###
#   cdef functions related with lists of monomials
###

cdef inline bint indivisible_in_list(ETuple m, list L, size_t i) noexcept:
    """
    Return if ``m`` divisible by any monomial in ``L[:i]``.
    """
    cdef size_t j
    for j in range(i):
        if (<ETuple>PyList_GET_ITEM(L,j)).divides(m):
            return False
    return True

cdef inline list interred(list L):
    """
    Return interreduction of a list of monomials.

    .. NOTE::

        The given list will be sorted in-place.

    INPUT:

    - ``L`` -- list of :class:`~sage.rings.polynomial.polydict.ETuple`

    OUTPUT:

    The interreduced list, where we interprete each ETuple as
    a monomial in a multivariate ring.
    """
    # First, we sort L ascendingly by total unweighted degree.
    # Afterwards, no monomial in L is divisible by a monomial
    # that appears later in L.
    if not L:
        return []
    L.sort(key=ETuple._unweighted_degree)
    cdef size_t i
    cdef ETuple m
    cdef list result = [<ETuple> PyList_GET_ITEM(L, 0)]
    for i in range(1, len(L)):
        m = <ETuple> PyList_GET_ITEM(L, i)
        if indivisible_in_list(m, L, i):
            result.append(m)
    return result

cdef list quotient(list L, ETuple m):
    """
    Return the quotient of the ideal represented by ``L`` and the
    monomial represented by ``m``.
    """
    cdef list result = L[:len(L)] # creates a copy
    cdef size_t i
    for i in range(len(L)):
        result.append((<ETuple>PyList_GET_ITEM(L,i)).divide_by_gcd(m))
    return interred(result)

cdef list quotient_by_var(list L, size_t index):
    """
    Return the quotient of the ideal represented by ``L`` and the
    variable number ``index``.
    """
    cdef ETuple m_j
    cdef list result = L[:len(L)] # creates a copy
    cdef size_t i
    for i in range(len(L)):
        if (<ETuple> PyList_GET_ITEM(L, i)).get_exp(index):
            result.append((<ETuple> PyList_GET_ITEM(L, i)).divide_by_var(index))
    return interred(result)

cdef ETuple sum_from_list(list L, size_t s, size_t l):
    """
    Compute the vector sum of the ETuples in ``L[s:s+l]`` in a balanced way.
    """
    if l == 1:
        return <ETuple>PyList_GET_ITEM(L,s)
    if l == 2:
        return (<ETuple>PyList_GET_ITEM(L,s)).eadd(<ETuple>PyList_GET_ITEM(L,s+1))
    cdef size_t l2 = l // 2
    cdef ETuple m1, m2
    m1 = sum_from_list(L, s, l2)
    m2 = sum_from_list(L, s+l2, l-l2)
    return m1.eadd(m2)

cdef bint HilbertBaseCase(Polynomial_integer_dense_flint fhs, Node D, tuple w) noexcept:
    """
    Try to compute the first Hilbert series of ``D.Id``, or return
    ``NotImplemented``.

    The third parameter ``w`` is a tuple of integers, the degree weights
    to be used.

    ``fhs`` is supposed to be zero. In some base cases, ``fhs`` will be
    changed to the value of the Hilbert series and ``True`` is returned.
    Otherwiese, ``False`` is returned.
    """
    cdef size_t i, j, exp
    cdef int e
    # First, the easiest cases:
    if not D.Id: # The zero ideal
        fmpz_poly_set_coeff_si(fhs._poly, 0, 1) # = PR(1)
        return True
    cdef ETuple m = <ETuple>PyList_GET_ITEM(D.Id, len(D.Id)-1)
    if m._nonzero == 0: # The one ideal
        return True

    # Second, another reasy case: D.Id is generated by variables.
    # D.Id is sorted ascendingly. Hence, if the last generator is a single
    # variable, then ALL are.
    cdef fmpz_poly_t poly_tmp
    if m._nonzero == 1 and m._data[1] == 1:
        fmpz_poly_init(poly_tmp)
        fmpz_poly_set_coeff_si(poly_tmp, 0, 1)
        fmpz_poly_set_coeff_si(fhs._poly, 0, 1) # = PR(1)
        if w is None:
            for i in range(len(D.Id)):
                exp = (<ETuple>PyList_GET_ITEM(D.Id, i)).unweighted_degree()
                fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                fmpz_poly_mul(fhs._poly, fhs._poly, poly_tmp)
                fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
        else:
            for i in range(len(D.Id)):
                exp = (<ETuple>PyList_GET_ITEM(D.Id, i)).weighted_degree(w)
                fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                fmpz_poly_mul(fhs._poly, fhs._poly, poly_tmp)
                fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
        fmpz_poly_clear(poly_tmp)
        return True # PR.prod([(1-t**degree(m,w)) for m in D.Id])

    # Thirdly, we test for proper powers of single variables.
    cdef bint easy = True
    for i in range(len(D.Id)):
        m = <ETuple>PyList_GET_ITEM(D.Id, i) # will be used later when we are breaking
        if m._nonzero > 1: # i.e., the generator contains more than a single var
            easy = False
            break
    if easy:
        # The ideal is generated by some powers of single variables, i.e., it splits.
        fmpz_poly_init(poly_tmp)
        fmpz_poly_set_coeff_si(poly_tmp, 0, 1)
        fmpz_poly_set_coeff_si(fhs._poly, 0, 1) # = PR(1)
        if w is None:
            for i in range(len(D.Id)):
                exp = (<ETuple>PyList_GET_ITEM(D.Id, i)).unweighted_degree()
                fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                fmpz_poly_mul(fhs._poly, fhs._poly, poly_tmp)
                fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
        else:
            for i in range(len(D.Id)):
                exp = (<ETuple>PyList_GET_ITEM(D.Id, i)).weighted_degree(w)
                fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                fmpz_poly_mul(fhs._poly, fhs._poly, poly_tmp)
                fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
        fmpz_poly_clear(poly_tmp)
        return True # PR.prod([(1-t**degree(m,w)) for m in D.Id])

    easy = True
    cdef ETuple m2
    cdef list v
    for j in range(i+1, len(D.Id)):
        if (<ETuple>PyList_GET_ITEM(D.Id,j))._nonzero > 1:
            # i.e., another generator contains more than a single var
            easy = False
            break
    cdef fmpz_poly_t FirstSummand, SecondSummand
    if easy:
        # The ideal only has a single non-simple power, in position i.
        # We have preserved that monomial, it is called m.
        fmpz_poly_init(poly_tmp)
        fmpz_poly_init(FirstSummand)
        fmpz_poly_init(SecondSummand)
        fmpz_poly_set_coeff_si(poly_tmp, 0, 1)
        fmpz_poly_set_coeff_si(FirstSummand, 0, 1)
        if w is None:
            fmpz_poly_set_coeff_si(SecondSummand, m.unweighted_degree(), -1)
            # Since the ideal is interreduced and all other monomials are
            # simple powers, we have the following formula:
            #~ return prod([(1-t**degree(m2,w)) for m2 in D.Id if m2 is not m])
            #~   - t**degree(m,w)*prod((1-t**quotient_degree(m2,m,w)) for m2 in D.Id if m is not m2)
            for j in range(len(D.Id)):
                if i != j:
                    m2 = <ETuple>PyList_GET_ITEM(D.Id,j)
                    exp = m2.unweighted_degree()
                    fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                    fmpz_poly_mul(FirstSummand, FirstSummand, poly_tmp)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
                    exp = m2.unweighted_quotient_degree(m)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                    fmpz_poly_mul(SecondSummand, SecondSummand, poly_tmp)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
        else:
            fmpz_poly_set_coeff_si(SecondSummand, m.weighted_degree(w), -1)
            # Since the ideal is interreduced and all other monomials are
            # simple powers, we have the following formula:
            #~ return prod([(1-t**degree(m2,w)) for m2 in D.Id if m2 is not m])
            #~   - t**degree(m,w)*prod((1-t**quotient_degree(m2,m,w)) for m2 in D.Id if m is not m2)
            for j in range(len(D.Id)):
                if i != j:
                    m2 = <ETuple>PyList_GET_ITEM(D.Id,j)
                    exp = m2.weighted_degree(w)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                    fmpz_poly_mul(FirstSummand, FirstSummand, poly_tmp)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
                    exp = m2.weighted_quotient_degree(m,w)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, -1)
                    fmpz_poly_mul(SecondSummand, SecondSummand, poly_tmp)
                    fmpz_poly_set_coeff_si(poly_tmp, exp, 0)
        fmpz_poly_clear(poly_tmp)
        fmpz_poly_add(fhs._poly, fhs._poly, FirstSummand)
        fmpz_poly_add(fhs._poly, fhs._poly, SecondSummand)
        fmpz_poly_clear(FirstSummand)
        fmpz_poly_clear(SecondSummand)
        return True
    # We are in a truly difficult case and give up for now...
    return False

cdef make_children(Node D, tuple w):
    """
    Create child nodes in ``D`` that allow to compute the first Hilbert
    series of ``D.Id``.

    Basically, the first Hilbert series of ``D.Id`` will be
    ``D.LMult`` times the first Hilbert series of ``D.Left.Id``,
    possibly plus ``D.RMult`` times the first Hilbert series of ``D.Right.Id``
    if ``D.Right`` is not ``None``.
    """
    cdef size_t j, m
    cdef int i, ii
    # Determine the variable that appears most often in the monomials.
    # If "most often" means "only once", then instead we choose a variable that is
    # guaranteed to appear in a composed monomial.
    # We will raise it to a reasonably high power that still guarantees that
    # many monomials will be divisible by it.
    cdef ETuple all_exponents = sum_from_list(D.Id, 0, len(D.Id))
    m = 0
    cdef list max_exponents = []
    for i in range(0, 2*all_exponents._nonzero, 2):
        j = all_exponents._data[i+1]
        if j > m:
            max_exponents = [all_exponents._data[i]]
            m = j
        elif j == m:
            max_exponents.append(all_exponents._data[i])
    cdef size_t e # will be the exponent, if the monomial used for cutting is power of a variable
    cdef ETuple cut
    cdef list Id2
    cdef size_t cut_deg
    # Cases:
    # - m==1, which means that all variables occur at most once.
    #   => we cut by a variable that appears in a decomposable generator
    # - max_exponents = [j]
    #   => cut = var(j)**e, where e is the median of all j-exponents
    # - max_exponents = [j1,...,jk]
    #   => cut = prod([var(j1),...,var(jk)]) or something of that type.
    if m == 1:
        # D.Id is sorted, which means that the last generator is decomposable
        all_exponents = <ETuple> PyList_GET_ITEM(D.Id, len(D.Id)-1)
        j = all_exponents._data[2*all_exponents._nonzero-2]
        cut = all_exponents._new()
        cut._nonzero = 1
        cut._data = <int*> sig_malloc(sizeof(int)*2)
        cut._data[0] = j
        cut._data[1] = 1
        if w is None:
            cut_deg = cut.unweighted_degree()
        else:
            cut_deg = cut.weighted_degree(w)
        # var(j) *only* appears in D.Id[-1]. Hence, D.Id+var(j) will be a split case,
        # with var(j) and D.Id[:-1]. So, we do the splitting right now.
        # Only the last generator contains var(j). Hence, D.Id/var(j) is obtained
        # from D.Id by adding the quotient of its last generator divided by var(j),
        # of course followed by interreduction.
        #~ D.LMult = 1 - t**degree(cut,w)
        fmpz_poly_set_coeff_si(D.LMult, 0, 1)
        fmpz_poly_set_coeff_si(D.LMult, cut_deg, -1)
        D.Left  = Node.__new__(Node)
        D.Left.Id = D.Id[:len(D.Id)-1]
        D.Left.Back = D
        Id2 = D.Id[:len(D.Id)-1]
        Id2.append(all_exponents.divide_by_var(j))
        D.Right = Node.__new__(Node)
        D.Right.Id = interred(Id2)
        D.Right.Back = D
        #~ D.RMult = 1 - D.LMult
        fmpz_poly_set_coeff_si(D.RMult, cut_deg, 1)
    else:
        j = max_exponents[0]
        e = median([(<ETuple> mon).get_exp(j) for mon in D.Id if (<ETuple> mon).get_exp(j)])
        cut = all_exponents._new()
        cut._nonzero = 1
        cut._data = <int*> sig_malloc(sizeof(int)*2)
        cut._data[0] = j
        cut._data[1] = e
        if w is None:
            cut_deg = cut.unweighted_degree()
        else:
            cut_deg = cut.weighted_degree(w)
        try:
            i = D.Id.index(cut)
        except ValueError:
            i = -1
        if i >= 0:
            # var(j)**e is a generator. Hence, e is the maximal exponent of
            # var(j) in D.Id, by D.Id being interreduced. But it also is the
            # truncated median, hence, there cannot be smaller exponents
            # (for otherwise the median would be strictly smaller than the
            # maximum). Conclusion: var(j) only appears in the generator var(j)**e
            #  -- we have a split case.
            Id2 = D.Id[:len(D.Id)]
            Id2.pop(i)
            #~ D.LMult = 1-t**degree(cut,w)
            fmpz_poly_set_coeff_si(D.LMult, 0, 1)
            fmpz_poly_set_coeff_si(D.LMult, cut_deg, -1)
            D.Left  = Node.__new__(Node)
            D.Left.Id = Id2
            D.Left.Back = D
            D.Right = None
        else:
            cut = all_exponents._new()
            cut._nonzero = 1
            cut._data = <int*> sig_malloc(sizeof(int)*2)
            cut._data[0] = j
            cut._data[1] = e
            if e > 1:
                # D.LMult = 1
                fmpz_poly_set_coeff_si(D.LMult, 0, 1)
                Id2 = D.Id[:len(D.Id)]
                Id2.append(cut)
                D.Left  = Node.__new__(Node)
                D.Left.Id = interred(Id2)
                D.Left.Back = D
                D.Right = Node.__new__(Node)
                D.Right.Id = quotient(D.Id,cut)
                D.Right.Back = D
            else:
                # m > 1, therefore var(j) cannot be a generator (D.Id is interreduced).
                # D.Id+var(j) will be a split case. So, we do the splitting right now.
                #~ D.LMult = 1-t**(1 if w is None else w[j])
                fmpz_poly_set_coeff_si(D.LMult, 0, 1)
                fmpz_poly_set_coeff_si(D.LMult, (1 if w is None else w[j]), -1)
                D.Left  = Node.__new__(Node)
                D.Left.Id = [(<ETuple> mon) for mon in D.Id if (<ETuple> mon).get_exp(j) == 0]
                D.Left.Back = D
                D.Right = Node.__new__(Node)
                D.Right.Id = quotient_by_var(D.Id, j)
                D.Right.Back = D
            #~ D.RMult = t**(e if w is None else e*w[j])
            fmpz_poly_set_coeff_si(D.RMult, (e if w is None else e * w[j]), 1)
     #else:
     #    It may be a good idea to form the product of some of the most
     #    frequent variables. But this isn't implemented yet. TODO?


def first_hilbert_series(I, grading=None, return_grading=False):
    """
    Return the first Hilbert series of the given monomial ideal.

    INPUT:

    - ``I`` -- a monomial ideal (possibly defined in singular)
    - ``grading`` -- (optional) a list or tuple of integers used as
      degree weights
    - ``return_grading`` -- boolean (default: ``False``); whether to return the
      grading

    OUTPUT:

    A univariate polynomial, namely the first Hilbert function of ``I``, and
    if ``return_grading==True`` also the grading used to compute the series.

    EXAMPLES::

        sage: from sage.rings.polynomial.hilbert import first_hilbert_series

        sage: # needs sage.libs.singular
        sage: R = singular.ring(0,'(x,y,z)','dp')
        sage: I = singular.ideal(['x^2','y^2','z^2'])
        sage: first_hilbert_series(I)
        -t^6 + 3*t^4 - 3*t^2 + 1
        sage: first_hilbert_series(I, return_grading=True)
        (-t^6 + 3*t^4 - 3*t^2 + 1, (1, 1, 1))
        sage: first_hilbert_series(I, grading=(1,2,3))
        -t^12 + t^10 + t^8 - t^4 - t^2 + 1

    TESTS:

    We test against some corner cases::

        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: I = 0*R
        sage: first_hilbert_series(I)
        1
        sage: first_hilbert_series(singular(I))
        1
        sage: I = 1*R
        sage: first_hilbert_series(I)
        0
        sage: first_hilbert_series(singular(I))
        0
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    PR = PolynomialRing(ZZ, 't')
    cdef Node AN
    # The "active node". If a recursive computation is needed, it will be equipped
    # with a 'Left' and a 'Right' child node, and some 'Multipliers'. Later, the first Hilbert
    # series of the left child node will be stored in 'LeftFHS', and together with
    # the first Hilbert series of the right child node and the multiplier yields
    # the first Hilbert series of 'Id'.
    cdef tuple w
    cdef Polynomial_integer_dense_flint fhs = Polynomial_integer_dense_flint.__new__(Polynomial_integer_dense_flint)
    fhs._parent = PR
    fhs._is_gen = 0
    if isinstance(I, sage.interfaces.abc.SingularElement):
        S = I._check_valid()
        # First, we need to deal with quotient rings, which also covers the case
        # of graded commutative rings that arise as cohomology rings in odd characteristic.
        # We replace everything by a commutative version of the quotient ring.
        br = S('basering')
        if S.eval('isQuotientRing(basering)')=='1':
            L = S('ringlist(basering)')
            R = S('ring(list(%s[1..3],ideal(0)))'%L.name())
            R.set_ring()
            I = S('fetch(%s,%s)+ideal(%s)'%(br.name(),I.name(),br.name()))

        I = [ETuple([int(x) for x in S.eval('string(leadexp({}[{}]))'.format(I.name(), i)).split(',')])
              for i in range(1,int(S.eval('size({})'.format(I.name())))+1)]
        br.set_ring()
        if grading is None:
            w = tuple(int(S.eval('deg(var({}))'.format(i))) for i in range(1,int(S.eval('nvars(basering)'))+1))
        else:
            w = tuple(grading)
    else:
        try:
            I = [bla.exponents()[0] for bla in I if bla]
        except TypeError:
            I = [bla.exponents()[0] for bla in I.gens() if bla]
        if grading is not None:
            w = tuple(grading)
        else:
            w = None

    AN = Node.__new__(Node)
    AN.Id = interred(I)
    AN.Back = None

    # Invariant of this function:
    # At each point, got_result will be false, or got_result is true and fhs
    # is the first Hilbert series of AN.
    cdef bint got_result = HilbertBaseCase(fhs, AN, w)
    while True:
        if not got_result:
            make_children(AN, w)
            AN = AN.Left
            got_result = HilbertBaseCase(fhs, AN, w)
        else:
            if AN.Back is None: # We are back on top, i.e., fhs is the First Hilber Series of I
                if return_grading:
                    return fhs, w
                else:
                    return fhs
            if AN is AN.Back.Left: # We store fhs and proceed to the sibling
                # ... unless there is no sibling
                if AN.Back.Right is None:
                    AN = AN.Back
                    #~ fhs *= AN.LMult
                    fmpz_poly_mul(fhs._poly, fhs._poly, AN.LMult)
                    got_result = True
                else:
                    fmpz_poly_set(AN.Back.LeftFHS, fhs._poly)
                    fmpz_poly_set_si(fhs._poly, 0)
                    AN = AN.Back.Right
                    AN.Back.Left = None
                    got_result = HilbertBaseCase(fhs, AN, w)
            else: # FHS of the left sibling is stored, of the right sibling is known.
                AN = AN.Back
                AN.Right = None
                #~ fhs = AN.LMult*AN.LeftFHS + AN.RMult*fhs
                fmpz_poly_mul(AN.LMult, AN.LMult, AN.LeftFHS)
                fmpz_poly_mul(AN.RMult, AN.RMult, fhs._poly)
                fmpz_poly_add(fhs._poly, AN.LMult, AN.RMult)
                got_result = True


def hilbert_poincare_series(I, grading=None):
    r"""
    Return the Hilbert Poincaré series of the given monomial ideal.

    INPUT:

    - ``I`` -- a monomial ideal (possibly defined in Singular)
    - ``grading`` -- (optional) a tuple of degree weights

    EXAMPLES::

        sage: # needs sage.libs.singular
        sage: from sage.rings.polynomial.hilbert import hilbert_poincare_series
        sage: R = PolynomialRing(QQ,'x',9)
        sage: I = [m.lm()
        ....:      for m in ((matrix(R, 3, R.gens())^2).list() * R).groebner_basis()] * R
        sage: hilbert_poincare_series(I)
        (t^7 - 3*t^6 + 2*t^5 + 2*t^4 - 2*t^3 + 6*t^2 + 5*t + 1)/(t^4 - 4*t^3 + 6*t^2 - 4*t + 1)
        sage: hilbert_poincare_series((R * R.gens())^2, grading=range(1,10))
        t^9 + t^8 + t^7 + t^6 + t^5 + t^4 + t^3 + t^2 + t + 1

    The following example is taken from :issue:`20145`::

        sage: # needs sage.libs.singular
        sage: n=4; m=11; P = PolynomialRing(QQ, n*m, "x"); x = P.gens(); M = Matrix(n, x)
        sage: from sage.rings.polynomial.hilbert import first_hilbert_series
        sage: I = P.ideal(M.minors(2))
        sage: J = P * [m.lm() for m in I.groebner_basis()]
        sage: hilbert_poincare_series(J).numerator()
        120*t^3 + 135*t^2 + 30*t + 1
        sage: hilbert_poincare_series(J).denominator().factor()
        (t - 1)^14
    """
    cdef Polynomial_integer_dense_flint HP
    HP, grading = first_hilbert_series(I, grading=grading, return_grading=True)
    PR = HP._parent

    # If grading was None, but the ideal lives in Singular, then grading is now
    # the degree vector of Singular's basering.
    # Otherwise, it may still be None.
    if grading is None:
        if HP.leading_coefficient() >= 0:
            return HP / ((1-PR.gen())**I.ring().ngens())
        return (-HP) / (-((1-PR.gen())**I.ring().ngens()))
    else:
        t = PR.gen()
        if HP.leading_coefficient() >= 0:
            return HP / PR.prod([(1-t**d) for d in grading])
        return (-HP) / (-PR.prod([(1-t**d) for d in grading]))
