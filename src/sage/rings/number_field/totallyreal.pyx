# sage_setup: distribution = sagemath-flint

r"""
Enumeration of primitive totally real fields

This module contains functions for enumerating all primitive
totally real number fields of given degree and small discriminant.
Here a number field is called *primitive* if it contains no proper
subfields except `\QQ`.

See also :mod:`sage.rings.number_field.totallyreal_rel`, which handles the non-primitive
case using relative extensions.

ALGORITHM:

We use Hunter's algorithm ([Coh2000]_, Section 9.3) with modifications
due to Takeuchi [Tak1999]_ and the author [Voi2008]_.

We enumerate polynomials `f(x) = x^n + a_{n-1} x^{n-1} + \dots + a_0`.
Hunter's theorem gives bounds on `a_{n-1}` and `a_{n-2}`; then given
`a_{n-1}` and `a_{n-2}`, one can recursively compute bounds on `a_{n-3},
\dots, a_0`, using the fact that the polynomial is totally real by
looking at the zeros of successive derivatives and applying
Rolle's theorem. See [Tak1999]_ for more details.

EXAMPLES:

In this first simple example, we compute the totally real quadratic
fields of discriminant `\le 50`. ::

    sage: enumerate_totallyreal_fields_prim(2,50)
    [[5, x^2 - x - 1],
     [8, x^2 - 2],
     [12, x^2 - 3],
     [13, x^2 - x - 3],
     [17, x^2 - x - 4],
     [21, x^2 - x - 5],
     [24, x^2 - 6],
     [28, x^2 - 7],
     [29, x^2 - x - 7],
     [33, x^2 - x - 8],
     [37, x^2 - x - 9],
     [40, x^2 - 10],
     [41, x^2 - x - 10],
     [44, x^2 - 11]]
    sage: [d for d in range(5,50)
    ....:    if (is_squarefree(d) and d%4 == 1) or (d%4 == 0 and is_squarefree(d/4))]
    [5, 8, 12, 13, 17, 20, 21, 24, 28, 29, 33, 37, 40, 41, 44]

Next, we compute all totally real quintic fields of discriminant `\le 10^5`::

    sage: ls = enumerate_totallyreal_fields_prim(5,10^5) ; ls
    [[14641, x^5 - x^4 - 4*x^3 + 3*x^2 + 3*x - 1],
     [24217, x^5 - 5*x^3 - x^2 + 3*x + 1],
     [36497, x^5 - 2*x^4 - 3*x^3 + 5*x^2 + x - 1],
     [38569, x^5 - 5*x^3 + 4*x - 1],
     [65657, x^5 - x^4 - 5*x^3 + 2*x^2 + 5*x + 1],
     [70601, x^5 - x^4 - 5*x^3 + 2*x^2 + 3*x - 1],
     [81509, x^5 - x^4 - 5*x^3 + 3*x^2 + 5*x - 2],
     [81589, x^5 - 6*x^3 + 8*x - 1],
     [89417, x^5 - 6*x^3 - x^2 + 8*x + 3]]
     sage: len(ls)
     9

We see that there are 9 such fields (up to isomorphism!).

See also [Mar1980]_.

AUTHORS:

- John Voight (2007-09-01): initial version; various optimization tweaks
- John Voight (2007-10-09): added DSage module; added pari functions to avoid
  recomputations; separated DSage component
- Craig Citro and John Voight (2007-11-04): additional doctests and type checking
- Craig Citro and John Voight (2008-02-10): final modifications for submission
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein and John Voight
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_calloc, sig_free

import math
import sys

from sage.libs.gmp.mpz cimport *
from sage.libs.pari.all import pari
from cypari2.gen cimport Gen as pari_gen
from sage.libs.pari.misc cimport new_t_POL_from_int_star

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.rings.number_field.totallyreal_data import tr_data, int_has_small_square_divisor
from sage.rings.number_field.totallyreal_data cimport tr_data


cpdef double odlyzko_bound_totallyreal(int n) noexcept:
    r"""
    This function returns the unconditional Odlyzko bound for the root
    discriminant of a totally real number field of degree `n`.

    .. NOTE::

        The bounds for `n > 50` are not necessarily optimal.

    INPUT:

    - ``n`` -- integer; the degree

    OUTPUT: a lower bound on the root discriminant (as a real number)

    EXAMPLES::

        sage: from sage.rings.number_field.totallyreal import odlyzko_bound_totallyreal
        sage: [odlyzko_bound_totallyreal(n) for n in range(1, 5)]
        [1.0, 2.223, 3.61, 5.067]

    AUTHORS:

    - John Voight (2007-09-03)

    .. NOTE::

        The values are calculated by Martinet [Mar1980]_.
    """

    if n <= 10:
        dB = [1.,2.223,3.610,5.067,6.523,7.941,9.301,10.596,11.823,12.985][n-1]
    elif n <= 20:
        dB = [14.0831,15.1217,16.1047,17.0359,17.9192,18.7580,19.5555,20.3148,21.0386,21.7294][n-11]
    elif n <= 30:
        dB = [22.3896,23.0212,23.6261,24.2061,24.7628,25.2976,25.2976,26.3071,26.3071,27.2440][n-21]
    elif n <= 40:
        dB = [27.2440,28.1165,28.1165,28.9315,28.9315,29.6948,29.6948,30.4117,30.4117,31.0865][n-31]
    elif n <= 50:
        dB = [31.0865,31.7232,31.7232,32.3252,32.3252,32.8954,32.8954,33.4365,33.4365,33.9508][n-41]
    else:
        dB = 33.9508
    return dB


def enumerate_totallyreal_fields_prim(n, B, a = [], verbose=0, return_seqs=False,
                                      phc=False, keep_fields=False, t_2=False,
                                      just_print=False,
                                      return_pari_objects=True):
    r"""
    Enumerate primitive totally real fields of degree
    `n>1` with discriminant `d \leq B`; optionally one can specify the
    first few coefficients, where the sequence `a` corresponds to

    ::

        a[d]*x^n + ... + a[0]*x^(n-d)

    where ``length(a) = d+1``, so in particular always ``a[d] = 1``.

    .. NOTE::

        This is guaranteed to give all primitive such fields, and
        seems in practice to give many imprimitive ones.

    INPUT:

    - ``n`` -- integer; the degree
    - ``B`` -- integer; the discriminant bound
    - ``a`` -- list (default: ``[]``); the coefficient list to begin with
    - ``verbose`` -- (integer or string, default: 0) if ``verbose == 1``
      (or ``2``), then print to the screen (really) verbosely; if verbose is
      a string, then print verbosely to the file specified by verbose.
    - ``return_seqs`` -- boolean (default: ``False``); if ``True``, then return
      the polynomials as sequences (for easier exporting to a file)
    - ``phc`` -- boolean or integer (default: ``False``)
    - ``keep_fields`` -- boolean or integer (default: ``False``); if
      ``keep_fields`` is ``True``, then keep fields up to ``B*log(B)``; if
      ``keep_fields`` is an integer, then keep fields up to that integer.
    - ``t_2`` -- boolean or integer (default: ``False``); if ``t_2 = T``, then
      keep only polynomials with ``t_2 norm >= T``
    - ``just_print`` -- boolean (default: ``False``); if ``just_print`` is not
      ``False``, instead of creating a sorted list of totally real number
      fields, we simply write each totally real field we find to the file
      whose filename is given by ``just_print``. In this case, we don't
      return anything.
    - ``return_pari_objects`` -- boolean (default: ``True``); if
      both ``return_seqs`` and ``return_pari_objects`` are ``False`` then
      it returns the elements as Sage objects. Otherwise it returns PARI
      objects.

    OUTPUT:

    the list of fields with entries ``[d,f]``, where ``d`` is the
    discriminant and ``f`` is a defining polynomial, sorted by
    discriminant.


    AUTHORS:

    - John Voight (2007-09-03)
    - Craig Citro (2008-09-19): moved to Cython for speed improvement

    TESTS::

        sage: len(enumerate_totallyreal_fields_prim(2,10**4))
        3043
        sage: len(enumerate_totallyreal_fields_prim(3,3**8))
        237
        sage: len(enumerate_totallyreal_fields_prim(5,5**7))
        6
        sage: len(enumerate_totallyreal_fields_prim(2,2**15)) # long time
        9957
        sage: len(enumerate_totallyreal_fields_prim(3,3**10)) # long time
        2720
        sage: len(enumerate_totallyreal_fields_prim(5,5**8)) # long time
        103

    Each of the outputs must be elements of Sage if ``return_pari_objects``
    is set to ``False``::

        sage: enumerate_totallyreal_fields_prim(2, 10)
        [[5, x^2 - x - 1], [8, x^2 - 2]]
        sage: type(enumerate_totallyreal_fields_prim(2, 10)[0][1])
        <class 'cypari2.gen.Gen'>
        sage: enumerate_totallyreal_fields_prim(2, 10, return_pari_objects=False)[0][0].parent()
        Integer Ring
        sage: enumerate_totallyreal_fields_prim(2, 10, return_pari_objects=False)[0][1].parent()
        Univariate Polynomial Ring in x over Rational Field
        sage: enumerate_totallyreal_fields_prim(2, 10, return_seqs=True)[1][0][1][0].parent()
        Rational Field
    """

    cdef pari_gen B_pari, d, d_poly, keepB, nf, t2val, ngt2, ng
    cdef int *f_out
    cdef int counts[4]
    cdef int i, n_int, j, skip_jp
    cdef bint found, use_t2, phc_flag, verb_int, temp_bint
    cdef Py_ssize_t k0, lenS
    cdef tr_data T
    cdef Integer dB

    if not isinstance(n, Integer):
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("cannot coerce n (= %s) to an integer" % n)
    if (n < 1):
        raise ValueError("n must be at least 1.")

    # Initialize
    n_int = int(n)
    S = set()        # set of pairs (d, f)
    lenS = 0

    # This is just to quiet valgrind down
    B_pari = pari(0)
    d = B_pari
    d_poly = B_pari
    keepB = B_pari
    nf = B_pari
    t2val = B_pari
    ngt2 = B_pari
    ng = B_pari

    dB = Integer.__new__(Integer)
    dB_odlyzko = odlyzko_bound_totallyreal(n_int)
    mpz_set_d(dB.value, dB_odlyzko)
    dB = 40000*((dB+1)**n_int)
    for i from 0 <= i < 4:
        counts[i] = 0

    B_pari = pari(B)
    f_out = <int *>check_calloc(n_int + 1, sizeof(int))
    f_out[n_int] = 1

    if keep_fields:
        if isinstance(keep_fields, bool):
            keepB = pari(int(math.floor(B*math.log(B))))
        else:
            keepB = pari(keep_fields)
    else:
        keepB = pari(0)

    if B > keepB:
        keepB = B_pari

    if t_2:
        k0 = len(a)
        if isinstance(t_2, Integer):
            t2val = pari(t_2)
        else:
            t2val = pari(a[k0-2]**2-2*a[k0-3])
        use_t2 = 1
    else:
        use_t2 = 0

    if phc:
        phc_flag = 1
    else:
        phc_flag = 0

    if just_print:
        skip_jp = 0
        jp_file = open(just_print, "w")
    else:
        skip_jp = 1

    # Trivial case
    if n == 1:
        sig_free(f_out)
        if return_seqs:
            return [[0,0,0,0],[[1,[-1,1]]]]
        elif return_pari_objects:
            return [[1,pari('x-1')]]
        else:
            Px = PolynomialRing(QQ, 'x')
            return [[ZZ(1), Px.gen()-1]]

    if verbose:
        verb_int = 1
        saveout = sys.stdout
        if isinstance(verbose, str):
            fsock = open(verbose, 'w')
            sys.stdout = fsock
        # Else, print to screen
    else:
        verb_int = 0

    T = tr_data(n_int,B,a)
    if verb_int == 2:
        T.incr(f_out,verb_int,0,phc_flag)
    else:
        T.incr(f_out,0,0,phc_flag)

    while f_out[n]:
        nf = new_t_POL_from_int_star(f_out, n_int+1, 0)
        if verb_int:
            print("==>", nf, "[")
            for j from 0 <= j < n-1:
                print("%s " % f_out[j])
            print("%s]" % f_out[n - 1])
        d_poly = nf.poldisc()
        counts[0] += 1
        if d_poly > 0 and nf.polsturm() == n:
            da = int_has_small_square_divisor(Integer(d_poly))
            if d_poly > dB or d_poly <= B*da:
                counts[1] += 1
                if nf.polisirreducible():
                    counts[2] += 1
                    zk, d = nf.nfbasis_d()

                    if d <= keepB:
                        if verb_int:
                            print("has discriminant", d, end='')

                        # Find a minimal lattice element
                        counts[3] += 1
                        ng = <pari_gen>((<pari_gen>(pari([nf,zk]))).polredabs())

                        dng = (d, ng)

                        if skip_jp:
                            # Check if K is contained in the list.
                            found = dng in S
                            if found and verb_int:
                                print("but is not new")

                            ngt2 = <pari_gen>(ng[n_int-1]**2-2*ng[n_int-2])
                            if not found:
                                temp_bint = ngt2 >= t2val
                                if ((not use_t2) or temp_bint):
                                    if verb_int:
                                        print("and is new!")
                                    S.add(dng)
                                    lenS += 1
                        else:
                            if ((not use_t2) or ngt2 >= t2val):
                                jp_file.write(str([d, ng.Vecrev()]) + "\n")

                    else:
                        if verb_int:
                            print("has discriminant", abs(d), "> B")
                else:
                    if verb_int:
                        print("is not irreducible")
            else:
                if verb_int:
                    print("has discriminant", abs(d), "with no large enough square divisor")
        else:
            if verb_int:
                if d == 0:
                    print("is not squarefree")
                else:
                    print("is not totally real")

        if verb_int == 2:
            T.incr(f_out,verb_int,0,phc_flag)
        else:
            T.incr(f_out,0,0,phc_flag)

    if not skip_jp:
        if n_int == 2 and B >= 5 and ((not use_t2) or t2val <= 5):
            jp_file.write(str([2,[-1,-1,1]]) + "\n")
            if B >= 8 and B < 32:
                jp_file.write(str([2,[-2,0,1]]) + "\n")
        elif n_int == 3 and B >= 49 and ((not use_t2) or 5 >= t2val):
            jp_file.write(str([3,[1,-2,-1,1]]) + "\n")
        jp_file.close()
        sig_free(f_out)
        return

    # Convert S to a sorted list of pairs [d, f]
    # we sort only according to d
    # because we cannot compare t_POL objects (PARI polynomials)
    S = sorted([list(s) for s in S], key=lambda x: x[0])

    # In the application of Smyth's theorem above (and easy
    # irreducibility test), we exclude finitely many possibilities
    # which we must now throw back in.
    if n_int == 2 and B >= 5 and ((not use_t2) or t2val <= 5):
        S = [[5, pari('x^2-x-1')]] + S
        lenS += 1
        if B >= 8 and B < 32:
            S.insert(1, [8, pari('x^2-2')])
            lenS += 1
    elif n_int == 3 and B >= 49 and ((not use_t2) or 5 >= t2val):
        S = [[49, pari('x^3-x^2-2*x+1')]] + S
        lenS += 1
    # The polynomials with n = 4 define imprimitive number fields.

    # Now check for isomorphic fields
    lenS = weed_fields(S, lenS)

    # Output.
    if verb_int:
        print("=" * 80)
        print("Polynomials tested:", counts[0])
        print("Polynomials with sssd poldisc:", counts[1])
        print("Irreducible polynomials:", counts[2])
        print("Polynomials with nfdisc <= B:", counts[3])
        for i from 0 <= i < lenS:
            print(S[i])
        if isinstance(verbose, str):
            fsock.close()
        sys.stdout = saveout

    sig_free(f_out)
    # Make sure to return elements that belong to Sage
    if return_seqs:
        return [[ZZ(counts[i]) for i in range(4)],
                [[ZZ(s[0]), [QQ(x) for x in s[1].polrecip().Vec()]] for s in S]]
    elif return_pari_objects:
        return S
    else:
        Px = PolynomialRing(QQ, 'x')
        return [[ZZ(s[0]), Px([QQ(x) for x in s[1].list()])]
                for s in S]


def weed_fields(S, Py_ssize_t lenS=0):
    r"""
    Function used internally by the :func:`~enumerate_totallyreal_fields_prim`
    routine. (Weeds the fields listed by ``[discriminant, polynomial]``
    for isomorphism classes.) Returns the size of the resulting list.

    EXAMPLES::

        sage: ls = [[5,pari('x^2-3*x+1')],[5,pari('x^2-5')]]
        sage: sage.rings.number_field.totallyreal.weed_fields(ls)
        1
        sage: ls
        [[5, x^2 - 3*x + 1]]
    """
    cdef Py_ssize_t i, j, n
    if lenS == 0:
        lenS = len(S)
    i = 0
    if not lenS:
        return lenS
    n = len(S[0][1])-1
    while i < lenS-1:
        j = i+1
        while j < lenS and S[i][0] == S[j][0]:
            if S[i][1].nfisisom(S[j][1]):
                # Keep the one with a smallest T_2
                T_2i = S[i][1][n-1]**2 - 2*S[i][1][n-2]
                T_2j = S[j][1][n-1]**2 - 2*S[j][1][n-2]
                if T_2i <= T_2j:
                    S.pop(j)
                    lenS -= 1
                else:
                    t = S.pop(j)
                    S.pop(i)
                    S.insert(i, t)
                    lenS -= 1
            else:
                j += 1
        i += 1

    return lenS
