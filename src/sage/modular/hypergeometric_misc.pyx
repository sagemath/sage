"""
Some utility routines for the hypergeometric motives package that benefit
significantly from Cythonization.
"""
from cpython cimport array
from cysignals.signals cimport sig_check
from sage.rings.integer cimport Integer

cpdef hgm_coeffs(long long p, unsigned int f,
                 int prec, gamma, m, int D,
                 gtable, int gtable_prec, bint use_longs):
    r"""
    Compute coefficients for the hypergeometric trace formula.

    This function is not intended for direct user access.

    TESTS::

        sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
        sage: from sage.modular.hypergeometric_misc import hgm_coeffs
        sage: H = Hyp(cyclotomic=([3],[4]))
        sage: H.euler_factor(2, 7, cache_p=True)
        7*T^2 - 3*T + 1
        sage: gamma = H.gamma_array()
        sage: prec, gtable = H.gauss_table(7, 1, 2)
        sage: D = 1
        sage: hgm_coeffs(7, 1, 2, gamma, [0]*6, D, gtable, prec, False)
        [7, 2*7, 6*7, 7, 6, 4*7]

    Check issue from :issue:`28404`::

        sage: H = Hyp(cyclotomic=[[10,2],[1,1,1,1,1]])
        sage: u = H.euler_factor(2,79) # indirect doctest
        sage: u.reverse().is_weil_polynomial()
        True
    """
    from sage.rings.padics.factory import Zp

    cdef int gl, j, k, v, gv
    cdef long long i, l, q1, w, w1, w2, q2, r, r1
    cdef bint flip, use_longlongs

    q1 = p ** f - 1
    gl = len(gamma)
    cdef array.array gamma_array1 = array.array('i', gamma.keys())
    cdef array.array gamma_array2 = array.array('i', gamma.values())
    r_array = [0] * gl
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
    Rz = R.zero()
    if use_longs:
        q2 = p ** prec
        try:
            gtab2 = gtable
        except TypeError:
            gtab2 = array.array('l', [0]) * q1
            for r in range(q1):
                gtab2[r] = gtable[r].lift() % q2
    else:
        use_longlongs = (Integer(p) ** prec < 2 ** 63)
        if use_longlongs:
            gtab2 = array.array('q', [0]) * q1
            try:
                for r in range(q1):
                    gtab2[r] = gtable[r]
            except TypeError:
                for r in range(q1):
                    gtab2[r] = gtable[r].lift()
    if f == 1:
        for r in range(q1):
            digit_count[r] = r
    else:
        for r in range(q1):
            r1 = r
            w = 0
            for i in range(f):
                w += r1 % p
                r1 //= p
            digit_count[r] = w
    ans = [None] * q1
    for r in range(q1):
        sig_check()
        # Skip this coefficient if we already have it by symmetry.
        if ans[r] is not None:
            continue
        # Determine whether this term is forced to be zero
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
            ans[r] = Rz
        else:
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
                    w2 = gtab2[r1]  # cast to long long to avoid overflow
                    if gv > 0:
                        for j in range(gv):
                            w = w * w2 % q2
                    else:
                        for j in range(-gv):
                            w1 = w1 * w2 % q2
                elif use_longlongs:
                    w2 = gtab2[r1]
                    if gv > 0:
                        for j in range(gv):
                            u *= w2
                    else:
                        for j in range(-gv):
                            u1 *= w2
                else:
                    if gv > 0:
                        for j in range(gv):
                            u *= gtable[r1]
                    else:
                        for j in range(-gv):
                            u1 *= gtable[r1]
            if use_longs:
                u = R(w)
                u1 = R(w1)
            if i % 2:
                u = -u
            ans[r] = (u / u1) << l
        if f > 1:
            r1 = r
            for j in range(f-1):
                r1 = r1 * p % q1
                if ans[r1] is not None:
                    break
                ans[r1] = ans[r]
    if f == 1:
        return ans
    # Consolidate down to p-1 terms.
    ans2 = [Rz] * (p-1)
    for r1 in range(p-1):
        for r in range(r1, q1, p-1):
            ans2[r1] += ans[r]
    return ans2
