"""
Isolate Complex Roots of Polynomials

AUTHOR:

- Carl Witty (2007-11-18): initial version

This is an implementation of complex root isolation.  That is, given a
polynomial with exact complex coefficients, we compute isolating
intervals for the complex roots of the polynomial.  (Polynomials with
integer, rational, Gaussian rational, or algebraic coefficients are
supported.)

We use a simple algorithm.  First, we compute a squarefree decomposition
of the input polynomial; the resulting polynomials have no multiple roots.
Then, we find the roots numerically, using NumPy (at low precision) or
Pari (at high precision).  Then, we verify the roots using interval
arithmetic.

EXAMPLES::

    sage: x = polygen(ZZ)
    sage: (x^5 - x - 1).roots(ring=CIF)
    [(1.167303978261419?, 1),
     (-0.764884433600585? - 0.352471546031727?*I, 1),
     (-0.764884433600585? + 0.352471546031727?*I, 1),
     (0.181232444469876? - 1.083954101317711?*I, 1),
     (0.181232444469876? + 1.083954101317711?*I, 1)]
"""

#*****************************************************************************
#       Copyright (C) 2007 Carl Witty
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from copy import copy

from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.qqbar import AA, QQbar
from sage.arith.misc import sort_complex_numbers_for_display
from sage.rings.polynomial.refine_root import refine_root


def interval_roots(p, rts, prec):
    """
    We are given a squarefree polynomial p, a list of estimated roots,
    and a precision.

    We attempt to verify that the estimated roots are in fact distinct
    roots of the polynomial, using interval arithmetic of precision ``prec``.
    If we succeed, we return a list of intervals bounding the roots; if we
    fail, we return ``None``.

    EXAMPLES::

        sage: x = polygen(ZZ)
        sage: p = x^3 - 1
        sage: rts = [CC.zeta(3)^i for i in range(0, 3)]
        sage: from sage.rings.polynomial.complex_roots import interval_roots
        sage: interval_roots(p, rts, 53)
        [1, -0.500000000000000? + 0.866025403784439?*I,
            -0.500000000000000? - 0.866025403784439?*I]
        sage: interval_roots(p, rts, 200)
        [1, -0.500000000000000000000000000000000000000000000000000000000000?
              + 0.866025403784438646763723170752936183471402626905190314027904?*I,
            -0.500000000000000000000000000000000000000000000000000000000000?
              - 0.866025403784438646763723170752936183471402626905190314027904?*I]
    """

    CIF = ComplexIntervalField(prec)
    CIFX = CIF['x']

    ip = CIFX(p)
    ipd = CIFX(p.derivative())

    irts = []

    for rt in rts:
        irt = refine_root(ip, ipd, CIF(rt), CIF)
        if irt is None:
            return None
        irts.append(irt)

    return irts


def intervals_disjoint(intvs):
    """
    Given a list of complex intervals, check whether they are pairwise
    disjoint.

    EXAMPLES::

        sage: from sage.rings.polynomial.complex_roots import intervals_disjoint
        sage: a = CIF(RIF(0, 3), 0)
        sage: b = CIF(0, RIF(1, 3))
        sage: c = CIF(RIF(1, 2), RIF(1, 2))
        sage: d = CIF(RIF(2, 3), RIF(2, 3))
        sage: intervals_disjoint([a,b,c,d])
        False
        sage: d2 = CIF(RIF(2, 3), RIF(2.001, 3))
        sage: intervals_disjoint([a,b,c,d2])
        True
    """

    # This may be quadratic in perverse cases, but will take only
    # n log(n) time in typical cases.

    intvs = sorted(copy(intvs))

    column = []
    prev_real = None

    def column_disjoint():
        column.sort()

        row = []
        prev_imag = None

        def row_disjoint():
            for a in range(len(row)):
                for b in range(a+1, len(row)):
                    if row[a].overlaps(row[b]):
                        return False
            return True

        for (y_imag, y) in column:
            if prev_imag is not None and y_imag > prev_imag:
                if not row_disjoint():
                    return False
                row = []
            prev_imag = y_imag
            row.append(y)
        if not row_disjoint():
            return False
        return True

    for x in intvs:
        x_real = x.real()
        if prev_real is not None and x_real > prev_real:
            if not column_disjoint():
                return False
            column = []
        prev_real = x_real
        column.append((x.imag(), x))

    if not column_disjoint():
        return False
    return True


def complex_roots(p, skip_squarefree=False, retval='interval', min_prec=0):
    """
    Compute the complex roots of a given polynomial with exact
    coefficients (integer, rational, Gaussian rational, and algebraic
    coefficients are supported).  Returns a list of pairs of a root
    and its multiplicity.

    Roots are returned as a ComplexIntervalFieldElement; each interval
    includes exactly one root, and the intervals are disjoint.

    By default, the algorithm will do a squarefree decomposition
    to get squarefree polynomials.  The skip_squarefree parameter
    lets you skip this step.  (If this step is skipped, and the polynomial
    has a repeated root, then the algorithm will loop forever!)

    You can specify retval='interval' (the default) to get roots as
    complex intervals.  The other options are retval='algebraic' to
    get elements of QQbar, or retval='algebraic_real' to get only
    the real roots, and to get them as elements of AA.

    EXAMPLES::

        sage: from sage.rings.polynomial.complex_roots import complex_roots
        sage: x = polygen(ZZ)
        sage: complex_roots(x^5 - x - 1)
        [(1.167303978261419?, 1),
         (-0.764884433600585? - 0.352471546031727?*I, 1),
         (-0.764884433600585? + 0.352471546031727?*I, 1),
         (0.181232444469876? - 1.083954101317711?*I, 1),
         (0.181232444469876? + 1.083954101317711?*I, 1)]
        sage: v = complex_roots(x^2 + 27*x + 181)

    Unfortunately due to numerical noise there can be a small imaginary part to each
    root depending on CPU, compiler, etc, and that affects the printing order. So we
    verify the real part of each root and check that the imaginary part is small in
    both cases::

        sage: v  # random
        [(-14.61803398874990?..., 1), (-12.3819660112501...? + 0.?e-27*I, 1)]
        sage: sorted((v[0][0].real(),v[1][0].real()))
        [-14.61803398874989?, -12.3819660112501...?]
        sage: v[0][0].imag().upper() < 1e25
        True
        sage: v[1][0].imag().upper() < 1e25
        True

        sage: K.<im> = QuadraticField(-1)
        sage: eps = 1/2^100
        sage: x = polygen(K)
        sage: p = (x-1)*(x-1-eps)*(x-1+eps)*(x-1-eps*im)*(x-1+eps*im)

    This polynomial actually has all-real coefficients, and is very, very
    close to (x-1)^5::

        sage: [RR(QQ(a)) for a in list(p - (x-1)^5)]
        [3.87259191484932e-121, -3.87259191484932e-121]
        sage: rts = complex_roots(p)
        sage: [ComplexIntervalField(10)(rt[0] - 1) for rt in rts]
        [-7.8887?e-31, 0, 7.8887?e-31, -7.8887?e-31*I, 7.8887?e-31*I]

    We can get roots either as intervals, or as elements of QQbar or AA.

    ::

        sage: p = (x^2 + x - 1)
        sage: p = p * p(x*im)
        sage: p
        -x^4 + (im - 1)*x^3 + im*x^2 + (-im - 1)*x + 1

    Two of the roots have a zero real component; two have a zero
    imaginary component.  These zero components will be found slightly
    inaccurately, and the exact values returned are very sensitive to
    the (non-portable) results of NumPy.  So we post-process the roots
    for printing, to get predictable doctest results.

    ::

        sage: def tiny(x):
        ....:     return x.contains_zero() and x.absolute_diameter() <  1e-14
        sage: def smash(x):
        ....:     x = CIF(x[0]) # discard multiplicity
        ....:     if tiny(x.imag()): return x.real()
        ....:     if tiny(x.real()): return CIF(0, x.imag())
        sage: rts = complex_roots(p); type(rts[0][0]), sorted(map(smash, rts))
        (<class 'sage.rings.complex_interval.ComplexIntervalFieldElement'>,
         [-1.618033988749895?, -0.618033988749895?*I,
          1.618033988749895?*I, 0.618033988749895?])
        sage: rts = complex_roots(p, retval='algebraic'); type(rts[0][0]), sorted(map(smash, rts))
        (<class 'sage.rings.qqbar.AlgebraicNumber'>,
         [-1.618033988749895?, -0.618033988749895?*I,
          1.618033988749895?*I, 0.618033988749895?])
        sage: rts = complex_roots(p, retval='algebraic_real'); type(rts[0][0]), rts
        (<class 'sage.rings.qqbar.AlgebraicReal'>,
         [(-1.618033988749895?, 1), (0.618033988749895?, 1)])

    TESTS:

    Verify that :issue:`12026` is fixed::

        sage: f = matrix(QQ, 8, lambda i, j: 1/(i + j + 1)).charpoly()
        sage: from sage.rings.polynomial.complex_roots import complex_roots
        sage: len(complex_roots(f))
        8
    """

    if skip_squarefree:
        factors = [(p, 1)]
    else:
        factors = p.squarefree_decomposition()

    prec = 53
    while True:
        CC = ComplexField(prec)
        CCX = CC['x']

        all_rts = []
        ok = True

        for (factor, exp) in factors:
            cfac = CCX(factor)
            rts = cfac.roots(multiplicities=False)
            # Make sure the number of roots we found is the degree. If
            # we don't find that many roots, it's because the
            # precision isn't big enough and though the (possibly
            # exact) polynomial "factor" is squarefree, it is not
            # squarefree as an element of CCX.
            if len(rts) < factor.degree():
                ok = False
                break
            irts = interval_roots(factor, rts, max(prec, min_prec))
            if irts is None:
                ok = False
                break
            if retval != 'interval':
                factor = QQbar.common_polynomial(factor)
            all_rts.extend((irt, factor, exp) for irt in irts)

        if ok and intervals_disjoint([rt for (rt, fac, mult) in all_rts]):
            all_rts = sort_complex_numbers_for_display(all_rts)
            if retval == 'interval':
                return [(rt, mult) for (rt, fac, mult) in all_rts]
            elif retval == 'algebraic':
                return [(QQbar.polynomial_root(fac, rt), mult) for (rt, fac, mult) in all_rts]
            elif retval == 'algebraic_real':
                rts = []
                for (rt, fac, mult) in all_rts:
                    qqbar_rt = QQbar.polynomial_root(fac, rt)
                    if qqbar_rt.imag().is_zero():
                        rts.append((AA(qqbar_rt), mult))
                return rts

        prec = prec * 2
