# distutils: libraries = ratpoints
r"""
Hyperelliptic Curve Point Finding, via ratpoints (deprecated)

This module is deprecated, use PARI instead::

    sage: pari(EllipticCurve("389a1")).ellratpoints(4)
    [[-2, 0], [-2, -1], [-1, 1], [-1, -2], [0, 0], [0, -1], [1, 0], [1, -1], [3, 5], [3, -6], [4, 8], [4, -9], [-3/4, 7/8], [-3/4, -15/8]]
    sage: pari("[x^3 + x^2 - 2*x, 1]").hyperellratpoints(4)
    [[-2, 0], [-2, -1], [-1, 1], [-1, -2], [0, 0], [0, -1], [1, 0], [1, -1], [3, 5], [3, -6], [4, 8], [4, -9], [-3/4, 7/8], [-3/4, -15/8]]
"""

from cysignals.memory cimport sig_malloc, sig_realloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.misc.superseded import deprecation
deprecation(24531, "the module sage.libs.ratpoints is deprecated; use pari.ellratpoints or pari.hyperellratpoints instead")


cdef int process(long x, long z, mpz_t y, void *info0, int *quit):
    # ratpoints calls this function when it finds a point [x : y : z]
    # info0 is the pointer passed to ratpoints originally
    # if quit[0] is set to a nonzero value, ratpoints will abort immediately
    cdef point_list *plist = <point_list *> info0
    cdef long i
    if plist.array_size == plist.num_points:
        i = plist.array_size
        plist.array_size *= 2
        plist.xes = <long *> sig_realloc(plist.xes, plist.array_size * sizeof(long))
        plist.ys = <mpz_t *> sig_realloc(plist.ys, plist.array_size * sizeof(mpz_t))
        plist.zs = <long *> sig_realloc(plist.zs, plist.array_size * sizeof(long))
        while i < plist.array_size:
            mpz_init(plist.ys[i])
            i += 1
    plist.xes[plist.num_points] = x
    mpz_set(plist.ys[plist.num_points], y)
    plist.zs[plist.num_points] = z
    plist.num_points += 1
    if plist.max_num_points > 0:
        if plist.max_num_points == plist.num_points:
            quit[0] = -1
    return 1 # weight for counting the points

def ratpoints(list coeffs, long H, verbose=False, long max=0,
              min_x_denom=None, max_x_denom=None, intervals=[]):
    """
    Access the ratpoints library to find points on the hyperelliptic curve:

    `y^2 = a_n x^n + \cdots + a_1 x + a_0.`

    INPUT:

    - ``coeffs`` -- list of integer coefficients `a_0` , `a_1`, ..., `a_n`

    - ``H`` -- the bound for the denominator and the absolute value of the
      numerator of the `x`-coordinate

    - ``verbose`` -- if ``True``, ratpoints will print comments about
      its progress

    - ``max`` -- maximum number of points to find (if 0, find all of them)

    OUTPUT:

    The points output by this program are points in (1, ceil(n/2), 1)-weighted
    projective space. If n is even, then the associated homogeneous equation is
    `y^2 = a_n x^n + \cdots + a_1 x z^{n-1} + a_0 z^n` while if n is odd, it is
    `y^2 = a_n x^n z + \cdots + a_1 x z^n + a_0 z^{n+1}`.

    EXAMPLES::

        sage: from sage.libs.ratpoints import ratpoints
        doctest:...: DeprecationWarning: the module sage.libs.ratpoints is deprecated; use pari.ellratpoints or pari.hyperellratpoints instead
        See http://trac.sagemath.org/24531 for details.
        sage: for x,y,z in ratpoints([1..6], 200):
        ....:     print(-1*y^2 + 1*z^6 + 2*x*z^5 + 3*x^2*z^4 + 4*x^3*z^3 + 5*x^4*z^2 + 6*x^5*z)
        0
        0
        0
        0
        0
        0
        0
        sage: for x,y,z in ratpoints([1..5], 200):
        ....:    print(-1*y^2 + 1*z^4 + 2*x*z^3 + 3*x^2*z^2 + 4*x^3*z + 5*x^4)
        0
        0
        0
        0
        0
        0
        0
        0

        sage: for x,y,z in ratpoints([1..200], 1000):
        ....:    print("{} {} {}".format(x,y,z))
        1 0 0
        0 1 1
        0 -1 1
        201 25353012004564588029934064107520000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 200
        201 -25353012004564588029934064107520000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 200

    The denominator of `x` can be restricted, for example to find
    integral points::

        sage: from sage.libs.ratpoints import ratpoints
        sage: coeffs = [400, -112, 0, 1]
        sage: ratpoints(coeffs, 10^6, max_x_denom=1, intervals=[[-10,0],[1000,2000]])
        [(1, 0, 0), (-8, 28, 1), (-8, -28, 1), (-7, 29, 1), (-7, -29, 1),
         (-4, 28, 1), (-4, -28, 1), (0, 20, 1), (0, -20, 1), (1368, 50596, 1),
         (1368, -50596, 1), (1624, 65444, 1), (1624, -65444, 1)]

        sage: ratpoints(coeffs, 1000, min_x_denom=100, max_x_denom=200)
        [(1, 0, 0),
        (-656, 426316, 121),
        (-656, -426316, 121),
        (452, 85052, 121),
        (452, -85052, 121),
        (988, 80036, 121),
        (988, -80036, 121),
        (-556, 773188, 169),
        (-556, -773188, 169),
        (264, 432068, 169),
        (264, -432068, 169)]

    Finding the integral points on the compact component of an elliptic curve::

        sage: E = EllipticCurve([0,1,0,-35220,-1346400])
        sage: e1, e2, e3 = E.division_polynomial(2).roots(multiplicities=False)
        sage: coeffs = [E.a6(),E.a4(),E.a2(),1]
        sage: ratpoints(coeffs, 1000, max_x_denom=1, intervals=[[e3,e2]])
        [(1, 0, 0),
        (-165, 0, 1),
        (-162, 366, 1),
        (-162, -366, 1),
        (-120, 1080, 1),
        (-120, -1080, 1),
        (-90, 1050, 1),
        (-90, -1050, 1),
        (-85, 1020, 1),
        (-85, -1020, 1),
        (-42, 246, 1),
        (-42, -246, 1),
        (-40, 0, 1)]
    """
    cdef ratpoints_args args
    cdef long i, total, verby
    cdef Integer sage_int, s_x, s_y, s_z
    cdef point_list *plist


    verby = ~0 if verbose else 0

    # Set the coefficient array:
    coeffs = [Integer(a) for a in coeffs]
    args.degree = len(coeffs)-1
    args.cof = <mpz_t *> sig_malloc((args.degree+1) * sizeof(mpz_t))

    # Create an array to hold the points found:
    plist = <point_list *> sig_malloc(sizeof(point_list))
    if max == 0:
        plist.array_size = 64
    else:
        plist.array_size = max
    plist.xes = <long *> sig_malloc(plist.array_size * sizeof(long))
    plist.ys = <mpz_t *> sig_malloc(plist.array_size * sizeof(mpz_t))
    for i from 0 <= i < plist.array_size:
        mpz_init(plist.ys[i])
    plist.zs = <long *> sig_malloc(plist.array_size * sizeof(long))
    plist.num_points = 0
    plist.max_num_points = max

    # Set the height bound:
    args.height = H

    # Set the intervals to be searched, including any specified:
    args.num_inter = len(intervals)
    args.domain = <ratpoints_interval *> sig_malloc((args.num_inter + args.degree) * sizeof(ratpoints_interval))
    for i,I in enumerate(intervals):
        args.domain[i].low = I[0]
        args.domain[i].up  = I[1]

    # Set the minimum and maximum denominators:
    if not min_x_denom:  min_x_denom = 1
    if not max_x_denom:  max_x_denom = H
    args.b_low = min_x_denom
    args.b_high = max_x_denom

    # Set the remaining arguments, whose non-default use is technical
    # (see ratpoints documentation)
    args.sp1 = RATPOINTS_DEFAULT_SP1
    args.sp2 = RATPOINTS_DEFAULT_SP2
    args.array_size = RATPOINTS_ARRAY_SIZE
    args.sturm = RATPOINTS_DEFAULT_STURM
    args.num_primes = RATPOINTS_DEFAULT_NUM_PRIMES
    args.max_forbidden = RATPOINTS_DEFAULT_MAX_FORBIDDEN
    args.flags = (RATPOINTS_VERBOSE & verby)

    for i from 0 <= i <= args.degree:
        mpz_init(args.cof[i])
        sage_int = <Integer> coeffs[i]
        mpz_set(args.cof[i], sage_int.value)

    sig_on()
    total = find_points(&args, process, <void *>plist)
    sig_off()
    if total == RATPOINTS_NON_SQUAREFREE:
        raise RuntimeError('Polynomial must be square-free')
    if total == RATPOINTS_BAD_ARGS:
        raise RuntimeError('Bad arguments to ratpoints')

    for i from 0 <= i <= args.degree:
        mpz_clear(args.cof[i])

    sig_free(args.cof)
    sig_free(args.domain)

    cdef list L = []
    for i from 0 <= i < plist.num_points:
        s_x = Integer(0)
        s_y = Integer(0)
        s_z = Integer(0)
        mpz_set_si(s_x.value, plist.xes[i])
        mpz_set(s_y.value, plist.ys[i])
        mpz_set_si(s_z.value, plist.zs[i])
        L.append((s_x,s_y,s_z))

    for i from 0 <= i < plist.array_size:
        mpz_clear(plist.ys[i])
    sig_free(plist.xes)
    sig_free(plist.ys)
    sig_free(plist.zs)
    sig_free(plist)

    return L

cdef int process_exists_only(long x, long z, mpz_t y, void *info0, int *quit):
    cdef info_struct_exists_only *info_s = <info_struct_exists_only *>info0
    cdef Integer YY
    if info_s.verbose:
        YY = Integer(0); mpz_set(YY.value, y)
        print('Found point [ %d : %d : %d ], quitting' % (x, YY, z))
    quit[0] = -1
    return 1

cdef int ratpoints_mpz_exists_only(mpz_t *coeffs, long H, int degree, bint verbose) except -1:
    """
    Direct call to ratpoints to search for existence only.

    WARNING - The coefficient array will be modified by ratpoints.
    """
    cdef ratpoints_args args
    cdef info_struct_exists_only info_s
    cdef long total, verby = ~0 if verbose else 0
    info_s.verbose = verbose
    assert degree <= RATPOINTS_MAX_DEGREE
    args.degree = degree
    args.cof = coeffs
    args.domain = <ratpoints_interval *> sig_malloc(2*args.degree * sizeof(ratpoints_interval))
    args.height = H
    args.num_inter = 0
    args.b_low = 1
    args.b_high = H
    args.sp1 = RATPOINTS_DEFAULT_SP1
    args.sp2 = RATPOINTS_DEFAULT_SP2
    args.array_size = RATPOINTS_ARRAY_SIZE
    args.sturm = RATPOINTS_DEFAULT_STURM
    args.num_primes = RATPOINTS_DEFAULT_NUM_PRIMES
    args.max_forbidden = RATPOINTS_DEFAULT_MAX_FORBIDDEN
    args.flags = (RATPOINTS_VERBOSE & verby)
    sig_on()
    total = find_points(&args, process_exists_only, <void *>(&info_s))
    sig_off()
    sig_free(args.domain)
    if total == RATPOINTS_NON_SQUAREFREE:
        raise RuntimeError('Polynomial must be square-free')
    if total == RATPOINTS_BAD_ARGS:
        raise RuntimeError('Bad arguments to ratpoints')
    return 1 if (total > 0) else 0




