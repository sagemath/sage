# sage.doctest: needs sage.rings.finite_rings sage.schemes
# cython: cdivision=True
r"""
Orthogonal arrays (find recursive constructions)

This module implements several functions to find recursive constructions of
:mod:`Orthogonal Arrays <sage.combinat.designs.orthogonal_arrays>`.

The main function of this module, i.e. :func:`find_recursive_construction`,
queries all implemented recursive constructions of designs implemented in
:mod:`~sage.combinat.designs.orthogonal_arrays_build_recursive` in order to
obtain an `OA(k,n)`.

:func:`find_recursive_construction` is called by the
:func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array` function.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`find_recursive_construction` | Find a recursive construction of an `OA(k,n)` (calls all others ``find_*`` functions)
    :func:`find_product_decomposition` | Find `n_1n_2=n` to obtain an `OA(k,n)` by the product construction
    :func:`find_wilson_decomposition_with_one_truncated_group` | Find `rm+u=n` to obtain an `OA(k,n)` by Wilson's construction with one truncated column.
    :func:`find_wilson_decomposition_with_two_truncated_groups` | Find `rm+r_1+r_2=n` to obtain an `OA(k,n)` by Wilson's construction with two truncated columns.
    :func:`find_construction_3_3` | Find a decomposition for construction 3.3 from [AC07]_.
    :func:`find_construction_3_4` | Find a decomposition for construction 3.4 from [AC07]_.
    :func:`find_construction_3_5` | Find a decomposition for construction 3.5 from [AC07]_.
    :func:`find_construction_3_6` | Find a decomposition for construction 3.6 from [AC07]_.
    :func:`find_q_x` | Find integers `q,x` such that the `q-x` construction yields an `OA(k,n)`.
    :func:`find_thwart_lemma_3_5` | Find the values on which Lemma 3.5 from [Thwarts]_ applies.
    :func:`find_thwart_lemma_4_1` | Find a decomposition for Lemma 4.1 from [Thwarts]_.
    :func:`find_three_factor_product` | Find `n_1n_2n_3=n` to obtain an `OA(k,n)` by the three-factor product from [DukesLing14]_
    :func:`find_brouwer_separable_design` | Find `t(q^2+q+1)+x=n` to obtain an `OA(k,n)` by Brouwer's separable design construction.
    :func:`find_brouwer_van_rees_with_one_truncated_column` | Find `rm+x_1+...+x_c=n` such that the Brouwer-van Rees constructions yields a `OA(k,n)`.

REFERENCES:

.. [AC07] Concerning eight mutually orthogonal latin squares
  Julian R. Abel, Nicholas Cavenagh
  Journal of Combinatorial Designs
  Vol. 15, n.3, pp. 255-261
  2007

Functions
---------
"""

from sage.misc.cachefunc import cached_function
from sage.combinat.designs.orthogonal_arrays import orthogonal_array
from sage.rings.integer cimport smallInteger
from sage.arith.misc import prime_powers


@cached_function
def find_recursive_construction(k, n):
    r"""
    Find a recursive construction of an `OA(k,n)` (calls all others ``find_*`` functions).

    This determines whether an `OA(k,n)` can be built through the following
    constructions:

    - :func:`~sage.combinat.designs.orthogonal_arrays.wilson_construction`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_3`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_4`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_5`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_6`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_q_x`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.thwart_lemma_3_5`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.thwart_lemma_4_1`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.three_factor_product`
    - :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.brouwer_separable_design`

    INPUT:

    - ``k``, ``n`` -- integers

    OUTPUT:

    Return a pair ``f,args`` such that ``f(*args)`` returns the requested `OA`
    if possible, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_recursive_construction
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: count = 0
        sage: for n in range(10,150):
        ....:     k = designs.orthogonal_arrays.largest_available_k(n)
        ....:     if find_recursive_construction(k,n):
        ....:         count = count + 1
        ....:         f,args = find_recursive_construction(k,n)
        ....:         OA = f(*args)
        ....:         assert is_orthogonal_array(OA,k,n,2,verbose=True)
        sage: count
        56
    """
    assert k > 3

    for find_c in [find_product_decomposition,
                   find_wilson_decomposition_with_one_truncated_group,
                   find_wilson_decomposition_with_two_truncated_groups,
                   find_construction_3_3,
                   find_construction_3_4,
                   find_construction_3_5,
                   find_construction_3_6,
                   find_q_x,
                   find_thwart_lemma_3_5,
                   find_thwart_lemma_4_1,
                   find_three_factor_product,
                   find_brouwer_separable_design,
                   find_brouwer_van_rees_with_one_truncated_column]:
        res = find_c(k,n)
        if res:
            return res
    return False


cpdef find_product_decomposition(int k,int n):
    r"""
    Find `n_1n_2=n` to obtain an `OA(k,n)` by the product construction.

    If Sage can build a `OA(k,n_1)` and a `OA(k,n_2)` such that `n=n_1\times
    n_2` then a `OA(k,n)` can be built by a product construction (which
    correspond to Wilson's construction with no truncated column). This
    function look for a pair of integers `(n_1,n_2)` with `n1 \leq n_2`, `n_1
    \times n_2 = n` and such that both an `OA(k,n_1)` and an `OA(k,n_2)` are
    available.

    INPUT:

    - ``k``, ``n`` -- integers

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` is an `OA(k,n)` or ``False`` if no
    product decomposition was found.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_product_decomposition
        sage: f,args = find_product_decomposition(6, 84)
        sage: args
        (None, 6, 7, 12, (), False)
        sage: _ = f(*args)
    """
    cdef int n1,n2
    for n1 in range(2, n):
        n2 = n/n1  # n2 is decreasing along the loop
        if n2 < n1:
            break
        if n%n1:  # we want to iterate only through divisors of n1... it seems
                  # faster to use that rather than calling the divisors function
            continue
        if is_available(k, n1) and is_available(k, n2):
            from sage.combinat.designs.orthogonal_arrays import wilson_construction
            return wilson_construction, (None,k,n1,n2,(),False)
    return False

cpdef find_wilson_decomposition_with_one_truncated_group(int k,int n):
    r"""
    Find `rm+u=n` to obtain an `OA(k,n)` by Wilson's construction with one truncated column.

    This function looks for possible integers `m,t,u` satisfying that `mt+u=n` and
    such that Sage knows how to build a `OA(k,m)`, `OA(k,m+1)`, `OA(k+1,t)` and a
    `OA(k,u)`.

    INPUT:

    - ``k``, ``n`` -- integers

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` is an `OA(k,n)` or ``False`` if no
    decomposition with one truncated block was found.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_wilson_decomposition_with_one_truncated_group
        sage: f,args = find_wilson_decomposition_with_one_truncated_group(4,38)
        sage: args
        (None, 4, 5, 7, (3,), False)
        sage: _ = f(*args)

        sage: find_wilson_decomposition_with_one_truncated_group(4,20)
        False
    """
    cdef int r,u,m
    # If there exists a TD(k+1,t) then k+1 < t+2, i.e. k <= t
    for r in range(max(1,k),n-1):
        u = n%r
        # We ensure that 1<=u, and that there can exists a TD(k,u), i.e k<u+2
        # (unless u == 1)
        if u == 0 or (u>1 and k >= u+2):
            continue

        m = n // r
        # If there exists a TD(k,m) then k<m+2
        if k >= m+2:
            break

        if (is_available(k, m) and is_available(k, m + 1) and
            is_available(k + 1, r) and is_available(k, u)):
            from sage.combinat.designs.orthogonal_arrays import wilson_construction
            return wilson_construction, (None,k,r,m,(u,),False)

    return False

cpdef find_wilson_decomposition_with_two_truncated_groups(int k,int n):
    r"""
    Find `rm+r_1+r_2=n` to obtain an `OA(k,n)` by Wilson's construction with two truncated columns.

    Look for integers `r,m,r_1,r_2` satisfying `n=rm+r_1+r_2` and `1\leq r_1,r_2<r`
    and such that the following designs exist : `OA(k+2,r)`, `OA(k,r1)`,
    `OA(k,r2)`, `OA(k,m)`, `OA(k,m+1)`, `OA(k,m+2)`.

    INPUT:

    - ``k``, ``n`` -- integers

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` is an `OA(k,n)` or ``False`` if no
    decomposition with two truncated blocks was found.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_wilson_decomposition_with_two_truncated_groups
        sage: f,args = find_wilson_decomposition_with_two_truncated_groups(5,58)
        sage: args
        (None, 5, 7, 7, (4, 5), False)
        sage: _ = f(*args)
    """
    cdef int r,m_min,m_max,m,r1_min,r1_max,r1,r2,r1_p_r2
    for r in [1] + list(range(k+1, n-2)):
        # as r*1+1+1 <= n and because we need
        # an OA(k+2,r), necessarily r=1 or r >= k+1
        if not is_available(k+2,r):
            continue
        m_min = (n - (2*r-2))/r
        m_max = (n - 2)/r
        if m_min > 1:
            m_values = list(range(max(m_min, k - 1), m_max + 1))
        else:
            m_values = [1] + list(range(k - 1, m_max + 1))
        for m in m_values:
            r1_p_r2 = n-r*m # the sum of r1+r2
                            # it is automatically >= 2 since m <= m_max
            if (r1_p_r2 > 2*r-2 or
                not is_available(k,m  ) or
                not is_available(k,m+1) or
                not is_available(k,m+2)):
                continue

            r1_min = r1_p_r2 - (r-1)
            r1_max = min(r-1, r1_p_r2)
            if r1_min > 1:
                r1_values = range(max(k - 1, r1_min), r1_max + 1)
            else:
                r1_values = [1] + list(range(k-1, r1_max + 1))
            for r1 in r1_values:
                if not is_available(k,r1):
                    continue
                r2 = r1_p_r2-r1
                if is_available(k,r2):
                    assert n == r*m+r1+r2
                    from sage.combinat.designs.orthogonal_arrays import wilson_construction
                    return wilson_construction, (None,k,r,m,(r1,r2),False)
    return False

cpdef find_construction_3_3(int k,int n):
    r"""
    Find a decomposition for construction 3.3 from [AC07]_.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_3`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_3
        sage: find_construction_3_3(11,177)[1]
        (11, 11, 16, 1)
        sage: find_construction_3_3(12,11)
    """
    cdef int mm,nn,i
    for mm in range(k-1, n//2+1):
        if not(is_available(k, mm) and is_available(k, mm + 1)):
            continue

        for nn in range(2, n//mm+1):
            i = n-nn*mm
            if i <= 0:
                continue

            if is_available(k + i, nn) and is_available(k, mm + i):
                from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_3
                return construction_3_3, (k, nn, mm, i)

cpdef find_construction_3_4(int k, int n):
    r"""
    Find a decomposition for construction 3.4 from [AC07]_.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_4`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_4
        sage: find_construction_3_4(8,196)[1]
        (8, 25, 7, 12, 9)
        sage: find_construction_3_4(9,24)
    """
    cdef int mm,nn,i,r,s
    for mm in range(k-1,n/2+1):
        if (not is_available(k,mm+0) or
            not is_available(k,mm+1) or
            not is_available(k,mm+2)):
            continue

        for nn in range(2, n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            for s in range(1,min(i,nn)):
                r = i-s
                if (is_available(k + r + 1, nn) and
                    is_available(k, s) and
                    (is_available(k, mm + r) or is_available(k, mm + r + 1))):
                    from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_4
                    return construction_3_4, (k, nn, mm, r, s)

cpdef find_construction_3_5(int k, int n):
    r"""
    Find a decomposition for construction 3.5 from [AC07]_.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_5`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_5
        sage: find_construction_3_5(8,111)[1]
        (8, 13, 6, 9, 11, 13)
        sage: find_construction_3_5(9,24)
    """
    cdef int mm,i,nn,r,s,t
    for mm in range(2, n//2+1):
        if (mm+3 >= n or
            not is_available(k,mm+1) or
            not is_available(k,mm+2) or
            not is_available(k,mm+3)):
            continue

        for nn in range(2, n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            if not is_available(k+3,nn):
                continue

            # Enumerate all  r,s,t<nn such that r+s+t=i and r<=s
            for s in range(min(i+1,nn)):
                for r in range(max(0,i-nn-s), min(s+1,i-s+1,nn)):
                    t = i - r - s
                    if ((nn-r-1)*(nn-s) < t         and
                        (r==0 or is_available(k,r)) and
                        (s==0 or is_available(k,s)) and
                        (t==0 or is_available(k,t))):
                        from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_5
                        return construction_3_5, (k,nn,mm,r,s,t)

cpdef find_construction_3_6(int k,int n):
    r"""
    Find a decomposition for construction 3.6 from [AC07]_.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_6`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_6
        sage: find_construction_3_6(8,95)[1]
        (8, 13, 7, 4)
        sage: find_construction_3_6(8,98)
    """
    cdef int mm,nn,i

    for mm in range(k-1,n/2+1):
        if (not is_available(k,mm+0) or
            not is_available(k,mm+1) or
            not is_available(k,mm+2)):
            continue

        for nn in range(2, n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            if (is_available(k+i,nn) and
                smallInteger(nn).is_prime_power()):
                from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_6
                return construction_3_6, (k,nn,mm,i)

cpdef find_q_x(int k,int n):
    r"""
    Find integers `q,x` such that the `q-x` construction yields an `OA(k,n)`.

    See the documentation of :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_q_x` to find out what
    hypotheses the integers `q,x` must satisfy.

    .. WARNING::

        For efficiency reasons, this function checks that Sage can build an
        `OA(k+1,q-x-1)` and an `OA(k+1,q-x+1)`, which is stronger than checking
        that Sage can build a `OA(k,q-x-1)-(q-x-1).OA(k,1)` and a
        `OA(k,q-x+1)-(q-x+1).OA(k,1)`. The latter would trigger a lot of
        independent set computations in
        :func:`sage.combinat.designs.orthogonal_arrays.incomplete_orthogonal_array`.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_q_x`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_q_x
        sage: find_q_x(10,9)
        False
        sage: find_q_x(9,158)[1]
        (9, 16, 6)
    """
    cdef int q,x

    # n = (q-1)*(q-x) + x + 2
    #   = q^2 - q*x - q + 2*x + 2
    for q in range(max(3,k+2),n):
        # n-q**2+q-2 = 2x-qx
        #            = x(2-q)
        x = (n-q**2+q-2)/(2-q)
        if (x < q and
            0 < x and
            n == (q-1)*(q-x)+x+2             and
            is_available(k+1,q-x-1)          and
            is_available(k+1,q-x+1)          and
            # The next is always True, because q is a prime power
            # is_available(k+1,q) and
            is_available(k, x+2 )            and
            smallInteger(q).is_prime_power()):
            from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_q_x
            return construction_q_x, (k,q,x)
    return False

cpdef find_thwart_lemma_3_5(int k,int N):
    r"""
    Find the values on which Lemma 3.5 from [Thwarts]_ applies.

    OUTPUT:

    A pair ``(f,args)`` such that ``f(*args)`` returns an `OA(k,n)` or ``False``
    if the construction is not available.

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.thwart_lemma_3_5`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_thwart_lemma_3_5
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array

        sage: f,args = find_thwart_lemma_3_5(7,66)
        sage: args
        (7, 9, 7, 1, 1, 1, 0, False)
        sage: OA = f(*args)
        sage: is_orthogonal_array(OA,7,66,2)
        True

        sage: f,args = find_thwart_lemma_3_5(6,100)
        sage: args
        (6, 8, 10, 8, 7, 5, 0, True)
        sage: OA = f(*args)
        sage: is_orthogonal_array(OA,6,100,2)
        True

    Some values from [Thwarts]_::

        sage: kn = ((10,1046), (10,1048), (10,1059), (11,1524),
        ....:       (11,2164), (12,3362), (12,3992),  (12,3994))
        sage: for k,n in kn:
        ....:     print("{} {} {}".format(k,n,find_thwart_lemma_3_5(k,n)[1]))
        10 1046 (10, 13, 79, 9, 1, 0, 9, False)
        10 1048 (10, 13, 79, 9, 1, 0, 11, False)
        10 1059 (10, 13, 80, 9, 1, 0, 9, False)
        11 1524 (11, 19, 78, 16, 13, 13, 0, True)
        11 2164 (11, 27, 78, 23, 19, 16, 0, True)
        12 3362 (12, 16, 207, 13, 13, 11, 13, True)
        12 3992 (12, 19, 207, 16, 13, 11, 19, True)
        12 3994 (12, 19, 207, 16, 13, 13, 19, True)

        sage: for k,n in kn:                                                     # not tested -- too long
        ....:     assert designs.orthogonal_array(k,n,existence=True) is True
    """
    from sage.combinat.designs.orthogonal_arrays_build_recursive import thwart_lemma_3_5
    cdef int n,m,a,b,c,d,NN,na,nb,nc

    for n in prime_powers(k+2,N-2): # There must exist a OA(k+3,n) thus n>=k+2
                                    # At least 3 columns are nonempty thus n<N-2

        # we look for (m,n,a,b,c,d) with N = mn + a + b + c (+d) and
        # 0 <= a,b,c,d <= n
        # hence we have N/n-4 <= m <= N/n

        # 1. look for m,a,b,c,d with complement=False
        # (we restrict to a >= b >= c)
        for m in range(max(k-1,(N+n-1)/n-4), N/n+1):
            if not (is_available(k,m+0) and
                    is_available(k,m+1) and
                    is_available(k,m+2)):
                continue

            NN = N - n*m
            # as a >= b >= c and d <= n we can restrict the start of the loops
            for a in range(max(0, (NN-n+2)/3), min(n, NN)+1): # (NN-n+2)/3 <==> ceil((NN-n)/3)x
                if not is_available(k,a):
                    continue
                for b in range(max(0, (NN-n-a+1)/2), min(a, n+1-a, NN-a)+1):
                    if not is_available(k,b):
                        continue
                    for c in range(max(0, NN-n-a-b), min(b, n+1-a-b, NN-a-b)+1):
                        if not is_available(k,c):
                            continue

                        d = NN - (a + b + c)  # necessarily 0 <= d <= n
                        if d == 0:
                            return thwart_lemma_3_5, (k,n,m,a,b,c,0,False)
                        elif (k+4 <= n+1 and
                              is_available(k, d ) and
                              is_available(k,m+3)):
                            return thwart_lemma_3_5, (k,n,m,a,b,c,d,False)

        # 2. look for m,a,b,c,d with complement=True
        # (we restrict to a >= b >= c)
        for m in range(max(k-2,N/n-4), (N+n-1)/n):
            if not (is_available(k,m+1) and
                    is_available(k,m+2) and
                    is_available(k,m+3)):
                continue

            NN = N - n*m
            for a in range(max(0, (NN-n+2)/3), min(n, NN)+1):
                # (NN-n+2)/3 <==> ceil((NN-n)/3)
                if not is_available(k,a):
                    continue
                na = n-a
                for b in range(max(0, (NN-n-a+1)/2), min(a, NN-a)+1):
                    nb = n-b
                    if na+nb > n+1 or not is_available(k,b):
                        continue
                    for c in range(max(0, NN-n-a-b), min(b, NN-a-b)+1):
                        nc = n-c
                        if na+nb+nc > n+1 or not is_available(k,c):
                            continue

                        d = NN - (a + b + c)  # necessarily d <= n
                        if d == 0:
                            return thwart_lemma_3_5, (k,n,m,a,b,c,0,True)
                        elif (k+4 <= n+1 and
                              is_available(k, d ) and
                              is_available(k,m+4)):
                            return thwart_lemma_3_5, (k,n,m,a,b,c,d,True)

    return False

cpdef find_thwart_lemma_4_1(int k,int n):
    r"""
    Find a decomposition for Lemma 4.1 from [Thwarts]_.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.thwart_lemma_4_1`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_thwart_lemma_4_1
        sage: find_thwart_lemma_4_1(10,408)[1]
        (10, 13, 28)
        sage: find_thwart_lemma_4_1(10,50)
        False
    """
    cdef int p,i,imax,nn,mm

    #      n  = nn*mm+4(nn-2)
    # <=> n+8 = nn(mm+4)
    #
    # nn is a prime power dividing n+8
    for p,imax in smallInteger(n+8).factor():
        nn = 1
        for i in range(1,imax+1):
            nn *= p
            mm = (n+8)/nn-4
            if (k+4 > nn+1 or
                mm <= 1 or
                nn % 3 == 2 or
                not is_available(k,nn-2) or
                not is_available(k,mm+1) or
                not is_available(k,mm+3) or
                not is_available(k,mm+4)):
                continue

            from sage.combinat.designs.orthogonal_arrays_build_recursive import thwart_lemma_4_1
            return thwart_lemma_4_1,(k,nn,mm)

    return False

cpdef find_three_factor_product(int k,int n):
    r"""
    Find `n_1n_2n_3=n` to obtain an `OA(k,n)` by the three-factor product from [DukesLing14]_.

    INPUT:

    - ``k``, ``n`` -- integers

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.three_factor_product`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_three_factor_product
        sage: find_three_factor_product(10,648)[1]
        (9, 8, 9, 9)
        sage: find_three_factor_product(10,50)
        False
    """
    cdef int n1,n2,n3

    # we want to write n=n1*n2*n3 where n1<=n2<=n3 and we can build:
    # - a OA(k-1,n1)
    # - a OA( k ,n2)
    # - a OA( k ,n3)
    for n1 in smallInteger(n).divisors()[1:-1]:
        if not is_available(k-1,n1):
            continue
        for n2 in smallInteger(n/n1).divisors():
            n3 = n/n1/n2
            if (n2<n1 or
                n3<n2 or
                not is_available(k,n2) or
                not is_available(k,n3)):
                continue
            from sage.combinat.designs.orthogonal_arrays_build_recursive import three_factor_product
            return three_factor_product,(k-1,n1,n2,n3)

    return False

cpdef find_brouwer_separable_design(int k,int n):
    r"""
    Find `t(q^2+q+1)+x=n` to obtain an `OA(k,n)` by Brouwer's separable design construction.

    INPUT:

    - ``k``, ``n`` -- integers

    The assumptions made on the parameters `t,q,x` are explained in the
    documentation of
    :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.brouwer_separable_design`.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_brouwer_separable_design
        sage: find_brouwer_separable_design(5,13)[1]
        (5, 1, 3, 0)
        sage: find_brouwer_separable_design(5,14)
        False
    """
    from sage.combinat.designs.orthogonal_arrays_build_recursive import brouwer_separable_design
    cdef int q, x, baer_subplane_size, max_t, min_t, t, e1, e2, e3, e4

    for q in prime_powers(2, n):
        baer_subplane_size = q**2 + q + 1
        if baer_subplane_size > n:
            break
        #                       x <= q^2+1
        # <=>        n-t(q^2+q+1) <= q^2+1
        # <=>             n-q^2-1 <= t(q^2+q+1)
        # <=> (n-q^2-1)/(q^2+q+1) <= t

        min_t = (n - q**2 - 1) / baer_subplane_size
        max_t = min(n / baer_subplane_size, q**2 - q + 1)

        for t in range(min_t, max_t + 1):
            x = n - t * baer_subplane_size
            e1 = int(x != q**2 - q - t)
            e2 = int(x != 1)
            e3 = int(x != q**2)
            e4 = int(x != t + q + 1)

            # i)
            if (x == 0 and
                    is_available(k, t) and is_available(k, t + q)):
                return brouwer_separable_design, (k, t, q, x)

            # ii)
            elif (x == t + q and
                  is_available(k + e3, t) and
                  is_available(k, t + q) and
                  is_available(k + 1, t + q + 1)):
                return brouwer_separable_design, (k, t, q, x)

            # iii)
            elif (x == q**2 - q + 1 - t and
                  is_available(k, x) and
                  is_available(k + e2, t + 1)
                  and is_available(k + 1, t + q)):
                return brouwer_separable_design, (k, t, q, x)

            # iv)
            elif (x == q**2 + 1 and
                  is_available(k, x) and
                  is_available(k + e4, t + 1) and
                  is_available(k + 1, t + q + 1)):
                return brouwer_separable_design, (k, t, q, x)

            # v)
            elif (0 < x < q**2 - q + 1 - t and (e1 or e2) and
                  is_available(k, x) and
                  is_available(k + e1, t) and
                  is_available(k + e2, t + 1) and
                  is_available(k + 1, t + q)):
                return brouwer_separable_design, (k, t, q, x)

            # vi)
            elif (t + q < x < q**2 + 1 and (e3 or e4) and
                  is_available(k, x) and
                  is_available(k + e3, t) and
                  is_available(k + e4, t + 1) and
                  is_available(k + 1, t + q + 1)):
                return brouwer_separable_design, (k, t, q, x)

    return False

# Associates to n the list of k,x with x>1 such that there exists an
# OA(k,n+x)-OA(k,x). Useful in find_brouwer_separable_design
from sage.combinat.designs.database import QDM as __QDM
cdef dict _QDM = __QDM
cdef dict ioa_indexed_by_n_minus_x = {}
for x in _QDM.itervalues():
    for (n, _, _, u), (k, _) in x.items():
        if u > 1:
            if n not in ioa_indexed_by_n_minus_x:
                ioa_indexed_by_n_minus_x[n] = []
            ioa_indexed_by_n_minus_x[n].append((k, u))


def int_as_sum(int value, list S, int k_max):
    r"""
    Return a tuple `(s_1, s_2, \ldots, s_k)` of less then `k_max` elements of `S` such
    that `value = s_1 + s_2 + \ldots + s_k`. If there is no such tuples then the
    function returns ``None``.

    INPUT:

    - ``value`` -- integer

    - ``S`` -- list of integers

    - ``k_max`` -- integer

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import int_as_sum
        sage: D = int_as_sum(21,[5,12],100)
        sage: for k in range(20,40):
        ....:     print("{} {}".format(k, int_as_sum(k,[5,12],100)))
        20 (5, 5, 5, 5)
        21 None
        22 (12, 5, 5)
        23 None
        24 (12, 12)
        25 (5, 5, 5, 5, 5)
        26 None
        27 (12, 5, 5, 5)
        28 None
        29 (12, 12, 5)
        30 (5, 5, 5, 5, 5, 5)
        31 None
        32 (12, 5, 5, 5, 5)
        33 None
        34 (12, 12, 5, 5)
        35 (5, 5, 5, 5, 5, 5, 5)
        36 (12, 12, 12)
        37 (12, 5, 5, 5, 5, 5)
        38 None
        39 (12, 12, 5, 5, 5)
    """
    cdef int i,j,v,vv,max_value
    cdef dict D,new_D,last_D
    last_D = D = {value:tuple()}
    max_value = max(S)

    if k_max * max_value < value:
        return None

    # The answer for a given k can be easily deduced from the answer
    # for k-1. That's how we build the list, incrementally starting
    # from k=0
    for j in range(k-1,-1,-1):
        new_D = {}
        for i in S:
            for v in last_D:
                vv = v-i
                if vv == 0:
                    return D[v] + (i,)
                if (vv > 0            and   # The new integer i is too big
                    vv <= j*max_value and   # The new integer i is too small
                    vv not in D       and   # We had it in D     already
                    vv not in new_D):       # We had it in new_D already
                    new_D[vv] = D[v]+(i,)
        if not new_D:
            break
        D.update(new_D)
        last_D = new_D

    return None


cpdef find_brouwer_van_rees_with_one_truncated_column(int k,int n):
    r"""
    Find `rm+x_1+...+x_c=n` such that the Brouwer-van Rees constructions yields a `OA(k,n)`.

    Let `n=rm+\sum_{1\leq i\leq c}` such that `c\leq r`. The
    generalization of Wilson's construction found by Brouwer and van
    Rees (with one truncated column) ensures that an `OA(k,n)` exists
    if the following designs exist: `OA(k+1,r)`, `OA(k,m)`,
    `OA(k,\sum_{1\leq i\leq c} u_i)`, `OA(k,m+x_1)-OA(k,x_1)`, ...,
    `OA(k,m+x_c)-OA(k,x_c)`.

    For more information, see the documentation of
    :func:`~sage.combinat.designs.orthogonal_arrays.wilson_construction`.

    INPUT:

    - ``k``, ``n`` -- integers

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_brouwer_van_rees_with_one_truncated_column
        sage: find_brouwer_van_rees_with_one_truncated_column(5,53)[1]
        (None, 5, 7, 7, [[(2, 1), (2, 1)]])
        sage: find_brouwer_van_rees_with_one_truncated_column(6,96)[1]
        (None, 6, 7, 13, [[(3, 1), (1, 1), (1, 1)]])
    """
    cdef list available_multipliers
    cdef int kk,uu,r,m,remainder,max_multiplier
    cdef tuple values

    # We write n=rm+remainder
    for m in range(2, n//2):
        if not is_available(k,m):
            continue

        # List of x such that a OA(k,m+x)-OA(k,x) exists
        #
        # This is the list of integers that can be used as multipliers
        # for the points of the truncated column
        available_multipliers = []
        if is_available(k,m+1):
            available_multipliers.append(1)
        for kk,uu in ioa_indexed_by_n_minus_x.get(m,[]):
            if kk>=k:
                available_multipliers.append(uu)

        # We stop if there is no multiplier, or if 1 is the only
        # multiplier (those cases are handled by other functions)
        if (not available_multipliers or
            (len(available_multipliers) == 1 and available_multipliers[0] == 1)):
            continue

        max_multiplier = max(available_multipliers)
        for r in range(2, n//m+1):
            remainder = n-r*m
            if (remainder > r*max_multiplier or
                not is_available(k+1,r) or
                not is_available(k,remainder)):
                continue

            values = int_as_sum(remainder, available_multipliers, r)
            if values is not None:
                from sage.combinat.designs.orthogonal_arrays import wilson_construction
                return (wilson_construction,
                        (None,k,r,m,[[(x,1) for x in values]]))

    return False

from sage.combinat.designs.designs_pyx cimport _OA_cache, _OA_cache_size
cdef int is_available(int k,int n) except -1:
    r"""
    Return whether Sage can build an OA(k,n)

    INPUT:

    - ``k``, ``n`` -- integers
    """
    if n >= _OA_cache_size:
        return orthogonal_array(k,n,existence=True) is True
    if k <= _OA_cache[n].max_true:
        return True
    elif k >= _OA_cache[n].min_unknown:
        return False
    else:
        return orthogonal_array(k,n,existence=True) is True
