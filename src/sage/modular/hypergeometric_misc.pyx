"""
Some utility routines for the hypergeometric motives package that benefit
significantly from Cythonization.
"""
cpdef hgm_coeffs(long long p, int f, int prec, gamma, array.array m, int D,
                 gtable, int gtable_prec, bint use_longs):
    r"""
    Compute coefficients for the hypergeometric trace formula.

    This function is not intended for direct user access.

    TESTS::

        sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
        sage: import array
        sage: from sage.modular.hypergeometric_misc import hgm_coeffs
        sage: H = Hyp(cyclotomic=([3],[4]))
        sage: H.euler_factor(2, 7, cache_p=True)
        7*T^2 - 3*T + 1
        sage: gamma = H.gamma_array()
        sage: prec, gtable = H.gauss_table(7, 1, 2)
        sage: m = array.array('i', [0]*6)
        sage: D = 1
        sage: hgm_coeffs(7, 1, 2, gamma, m, D, gtable, prec, False)
        [7, 2*7, 6*7, 7, 6, 4*7]
    """
    from sage.rings.padics.factory import Zp

    cdef int gl, j, k, l, v, gv
    cdef long long i, q1, w, w1, w2, q2, r, r1
    cdef bint flip

    q1 = p ** f - 1
    gl = len(gamma)
    cdef array.array gamma_array1 = array.array('i', gamma.keys())
    cdef array.array gamma_array2 = array.array('i', gamma.values())
    cdef array.array r_array = array.array('i', [0]) * gl
    cdef array.array digit_count = array.array('i', [0]) * q1
    cdef array.array gtab2

    try:
        R = gtable[0].parent()
    except AttributeError:
        R = Zp(p, prec, "fixed-mod")
    # In certain cases, the reciprocals of the Gauss sums are reported
    # for efficiency.
    flip = (f == 1 and prec == 1 and gtable_prec == 1)
    ans = []
    if use_longs:
        q2 = p ** prec
        try:
            gtab2 = gtable
        except TypeError:
            gtab2 = array.array('l', [0]) * q1
            for r in range(q1):
                gtab2[r] = gtable[r].lift() % q2
    if f == 1:
        for r in range(q1):
            digit_count[r] = r
    else:
        for r in range(q1):
            r1 = r
            digit_count[r] = 0
            for i in range(f):
                digit_count[r] += r1 % p
                r1 //= p
    Rz = R.zero()
    for r in range(q1):
        # First determine whether this term is forced to be zero
        # for divisibility reasons. If so, skip the p-adic arithmetic.
        i = 0
        for k in range(gl):
            v = gamma_array1[k]
            gv = gamma_array2[k]
            r1 = v * r % q1
            r_array[k] = r1
            i += digit_count[r1] * gv
        i //= (p - 1)
        l = i + f * (D + m[0] - m[r])
        if l >= prec:
            ans.append(Rz)
            continue
        # Keep numerator and denominator separate for efficiency.
        if use_longs:
            w = 1
            w1 = 1
        else:
            u = R.one()
            u1 = R.one()
        for k in range(gl):
            gv = gamma_array2[k]
            r1 = r_array[k]
            if flip:
                gv = -gv
            if use_longs:
                w2 = gtab2[r1] # cast to long long to avoid overflow
                if gv > 0: 
                    for j in range(gv): w = w * w2 % q2
                else:
                    for j in range(-gv): w1 = w1 * w2 % q2
            else:
                if gv > 0:
                    for j in range(gv): u *= gtable[r1]
                else:
                    for j in range(-gv): u1 *= gtable[r1]
        if use_longs:
            u = R(w)
            u1 = R(w1)
        if i % 2: u = -u
        ans.append((u / u1) << l)
    return ans
