# sage_setup: distribution = sagemath-categories
r"""
Enumeration of rational points on projective schemes

Naive algorithms for enumerating rational points over `\QQ` or finite fields
over for general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.
    (For instance using an affine function for a projective scheme or a finite
    field function for a scheme defined over an infinite field.)

EXAMPLES:

Projective, over `\QQ`::

    sage: from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
    sage: P.<X,Y,Z> = ProjectiveSpace(2, QQ)
    sage: C = P.subscheme([X + Y - Z])
    sage: enum_projective_rational_field(C, 3)
    [(-2 : 3 : 1), (-1 : 1 : 0), (-1 : 2 : 1), (-1/2 : 3/2 : 1),
     (0 : 1 : 1), (1/3 : 2/3 : 1), (1/2 : 1/2 : 1), (2/3 : 1/3 : 1),
     (1 : 0 : 1), (3/2 : -1/2 : 1), (2 : -1 : 1), (3 : -2 : 1)]


Projective over a finite field::

    sage: from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
    sage: E = EllipticCurve('72').change_ring(GF(19))                                   # needs sage.schemes
    sage: enum_projective_finite_field(E)                                               # needs sage.schemes
    [(0 : 1 : 0), (1 : 0 : 1), (3 : 0 : 1), (4 : 9 : 1), (4 : 10 : 1),
     (6 : 6 : 1), (6 : 13 : 1), (7 : 6 : 1), (7 : 13 : 1), (9 : 4 : 1),
     (9 : 15 : 1), (12 : 8 : 1), (12 : 11 : 1), (13 : 8 : 1), (13 : 11 : 1),
     (14 : 3 : 1), (14 : 16 : 1), (15 : 0 : 1), (16 : 9 : 1), (16 : 10 : 1),
     (17 : 7 : 1), (17 : 12 : 1), (18 : 9 : 1), (18 : 10 : 1)]


AUTHORS:

- David R. Kohel <kohel@maths.usyd.edu.au>: original version.

- John Cremona and Charlie Turner <charlotteturner@gmail.com> (06-2010):
  improvements to clarity and documentation.

- Raghukul Raman <raghukul.raman01@gmail.com> (2018): Added sieve algorithm
"""

# ****************************************************************************
#       Copyright (C) 2010 William Stein, David Kohel, John Cremona, Charlie Turner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import product

from sage.arith.misc import gcd, next_prime, previous_prime, crt
from sage.arith.srange import srange
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.mrange import xmrange
from sage.schemes.generic.scheme import Scheme
from sage.parallel.ncpus import ncpus
from sage.parallel.use_fork import p_iter_fork
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.integer_ring import ZZ
from sage.schemes.generic.scheme import is_Scheme

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.rings.real_mpfr', 'RR')


def enum_projective_rational_field(X, B):
    r"""
    Enumerates projective, rational points on scheme ``X`` of height up to
    bound ``B``.

    INPUT:

    - ``X`` -- a scheme or set of abstract rational points of a scheme

    - ``B`` -- a positive integer bound

    OUTPUT:

    A list containing the projective points of ``X`` of height up to ``B``,
    sorted.

    EXAMPLES::

        sage: P.<X,Y,Z> = ProjectiveSpace(2, QQ)
        sage: C = P.subscheme([X + Y - Z])
        sage: from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
        sage: enum_projective_rational_field(C(QQ), 6)
        [(-5 : 6 : 1), (-4 : 5 : 1), (-3 : 4 : 1), (-2 : 3 : 1),
         (-3/2 : 5/2 : 1), (-1 : 1 : 0), (-1 : 2 : 1), (-2/3 : 5/3 : 1),
         (-1/2 : 3/2 : 1), (-1/3 : 4/3 : 1), (-1/4 : 5/4 : 1),
         (-1/5 : 6/5 : 1), (0 : 1 : 1), (1/6 : 5/6 : 1), (1/5 : 4/5 : 1),
         (1/4 : 3/4 : 1), (1/3 : 2/3 : 1), (2/5 : 3/5 : 1), (1/2 : 1/2 : 1),
         (3/5 : 2/5 : 1), (2/3 : 1/3 : 1), (3/4 : 1/4 : 1), (4/5 : 1/5 : 1),
         (5/6 : 1/6 : 1), (1 : 0 : 1), (6/5 : -1/5 : 1), (5/4 : -1/4 : 1),
         (4/3 : -1/3 : 1), (3/2 : -1/2 : 1), (5/3 : -2/3 : 1), (2 : -1 : 1),
         (5/2 : -3/2 : 1), (3 : -2 : 1), (4 : -3 : 1), (5 : -4 : 1),
         (6 : -5 : 1)]
        sage: enum_projective_rational_field(C,6) == enum_projective_rational_field(C(QQ),6)
        True

    ::

        sage: P3.<W,X,Y,Z> = ProjectiveSpace(3, QQ)
        sage: enum_projective_rational_field(P3, 1)
        [(-1 : -1 : -1 : 1), (-1 : -1 : 0 : 1), (-1 : -1 : 1 : 0), (-1 : -1 : 1 : 1),
         (-1 : 0 : -1 : 1), (-1 : 0 : 0 : 1), (-1 : 0 : 1 : 0), (-1 : 0 : 1 : 1),
         (-1 : 1 : -1 : 1), (-1 : 1 : 0 : 0), (-1 : 1 : 0 : 1), (-1 : 1 : 1 : 0),
         (-1 : 1 : 1 : 1), (0 : -1 : -1 : 1), (0 : -1 : 0 : 1), (0 : -1 : 1 : 0),
         (0 : -1 : 1 : 1), (0 : 0 : -1 : 1), (0 : 0 : 0 : 1), (0 : 0 : 1 : 0),
         (0 : 0 : 1 : 1), (0 : 1 : -1 : 1), (0 : 1 : 0 : 0), (0 : 1 : 0 : 1),
         (0 : 1 : 1 : 0), (0 : 1 : 1 : 1), (1 : -1 : -1 : 1), (1 : -1 : 0 : 1),
         (1 : -1 : 1 : 0), (1 : -1 : 1 : 1), (1 : 0 : -1 : 1), (1 : 0 : 0 : 0),
         (1 : 0 : 0 : 1), (1 : 0 : 1 : 0), (1 : 0 : 1 : 1), (1 : 1 : -1 : 1),
         (1 : 1 : 0 : 0), (1 : 1 : 0 : 1), (1 : 1 : 1 : 0), (1 : 1 : 1 : 1)]

    ALGORITHM:

    We just check all possible projective points in correct dimension
    of projective space to see if they lie on ``X``.

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010)
    """
    from sage.schemes.projective.projective_space import ProjectiveSpace_ring
    if isinstance(X, Scheme):
        if not isinstance(X.ambient_space(), ProjectiveSpace_ring):
            raise TypeError("ambient space must be projective space over the rational field")
        X = X(X.base_ring())
    elif not isinstance(X.codomain().ambient_space(), ProjectiveSpace_ring):
        raise TypeError("codomain must be projective space over the rational field")

    n = X.codomain().ambient_space().ngens()
    zero = (0,) * n
    pts = []
    for c in product(*[srange(-B, B + 1) for _ in range(n)]):
        if gcd(c) == 1 and c > zero:
            try:
                pts.append(X(c))
            except TypeError:
                pass
    pts.sort()
    return pts


def enum_projective_number_field(X, **kwds):
    """
    Enumerates projective points on scheme ``X`` defined over a number field.

    Simply checks all of the points of absolute height of at most ``B``
    and adds those that are on the scheme to the list.

    This algorithm computes 2 lists: L containing elements x in `K` such that
    H_k(x) <= B, and a list L' containing elements x in `K` that, due to
    floating point issues,
    may be slightly larger then the bound. This can be controlled
    by lowering the tolerance.

    ALGORITHM:

    This is an implementation of the revised algorithm (Algorithm 4) in
    [DK2013]_. Algorithm 5 is used for imaginary quadratic fields.

    INPUT: keyword arguments:

    - ``bound`` -- a real number

    - ``tolerance`` -- a rational number in (0,1] used in Doyle-Krumm
      algorithm-4

    - ``precision`` -- the precision to use for computing the elements of
      bounded height of number fields

    OUTPUT: a sorted list containing the projective points of ``X`` of absolute
    height up to ``B``

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.projective.projective_rational_point import enum_projective_number_field
        sage: u = QQ['u'].0
        sage: K = NumberField(u^3 - 5, 'v')
        sage: P.<x,y,z> = ProjectiveSpace(K, 2)
        sage: X = P.subscheme([x - y])
        sage: enum_projective_number_field(X(K), bound=RR(5^(1/3)), prec=2^10)          # needs sage.symbolic
        [(0 : 0 : 1), (1 : 1 : 0), (-1 : -1 : 1), (1 : 1 : 1)]

    ::

        sage: # needs sage.rings.number_field
        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 3, 'v')
        sage: A.<x,y> = ProjectiveSpace(K, 1)
        sage: X = A.subscheme(x - y)
        sage: enum_projective_number_field(X, bound=2)
        [(1 : 1)]
    """
    B = kwds.pop('bound')
    tol = kwds.pop('tolerance', 1e-2)
    prec = kwds.pop('precision', 53)
    from sage.schemes.projective.projective_space import ProjectiveSpace_ring
    if isinstance(X, Scheme):
        if not isinstance(X.ambient_space(), ProjectiveSpace_ring):
            raise TypeError("ambient space must be projective space over a number field")
        X = X(X.base_ring())
    else:
        if not isinstance(X.codomain().ambient_space(), ProjectiveSpace_ring):
            raise TypeError("codomain must be projective space over a number field")

    R = X.codomain().ambient_space()

    pts = []

    for P in R.points_of_bounded_height(bound=B, tolerance=tol, precision=prec):
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()
    return pts


def enum_projective_finite_field(X):
    """
    Enumerates projective points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -- a scheme defined over a finite field or a set of abstract
      rational points of such a scheme

    OUTPUT:

    A list containing the projective points of ``X`` over the finite field,
    sorted.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
        sage: F = GF(53)
        sage: P.<X,Y,Z> = ProjectiveSpace(2, F)
        sage: len(enum_projective_finite_field(P(F)))
        2863
        sage: 53^2 + 53 + 1
        2863

    ::

        sage: # needs sage.rings.finite_rings
        sage: F = GF(9, 'a')
        sage: P.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: C = Curve(X^3 - Y^3 + Z^2*Y)                                              # needs sage.schemes
        sage: enum_projective_finite_field(C(F))                                        # needs sage.schemes
        [(0 : 0 : 1), (0 : 1 : 1), (0 : 2 : 1), (1 : 1 : 0), (a + 1 : 2*a : 1),
         (a + 1 : 2*a + 1 : 1), (a + 1 : 2*a + 2 : 1), (2*a + 2 : a : 1),
         (2*a + 2 : a + 1 : 1), (2*a + 2 : a + 2 : 1)]

    ::

        sage: F = GF(5)
        sage: P2F.<X,Y,Z> = ProjectiveSpace(2, F)
        sage: enum_projective_finite_field(P2F)
        [(0 : 0 : 1), (0 : 1 : 0), (0 : 1 : 1), (0 : 2 : 1), (0 : 3 : 1), (0 : 4 : 1),
         (1 : 0 : 0), (1 : 0 : 1), (1 : 1 : 0), (1 : 1 : 1), (1 : 2 : 1), (1 : 3 : 1),
         (1 : 4 : 1), (2 : 0 : 1), (2 : 1 : 0), (2 : 1 : 1), (2 : 2 : 1), (2 : 3 : 1),
         (2 : 4 : 1), (3 : 0 : 1), (3 : 1 : 0), (3 : 1 : 1), (3 : 2 : 1), (3 : 3 : 1),
         (3 : 4 : 1), (4 : 0 : 1), (4 : 1 : 0), (4 : 1 : 1), (4 : 2 : 1), (4 : 3 : 1),
         (4 : 4 : 1)]

    ALGORITHM:

    Checks all points in projective space to see if they lie on ``X``.

    .. WARNING::

        If ``X`` is defined over an infinite field, this code will not finish!

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010).
    """
    from sage.schemes.projective.projective_space import ProjectiveSpace_ring
    if isinstance(X, Scheme):
        if not isinstance(X.ambient_space(), ProjectiveSpace_ring):
            raise TypeError("ambient space must be projective space over a finite")
        X = X(X.base_ring())
    elif not isinstance(X.codomain().ambient_space(), ProjectiveSpace_ring):
        raise TypeError("codomain must be projective space over a finite field")

    n = X.codomain().ambient_space().ngens() - 1
    F = X.value_ring()
    pts = []
    for k in range(n + 1):
        for c in product(*[F for _ in range(k)]):
            try:
                pts.append(X(list(c) + [1] + [0] * (n - k)))
            except TypeError:
                pass
    pts.sort()
    return pts


def sieve(X, bound):
    r"""
    Return the list of all projective, rational points on scheme ``X`` of
    height up to ``bound``.

    Height of a projective point `X = (x_1, x_2,\dots, x_n)` is given by
    `H_X = \max(y_1, y_2,\dots, y_n)`, where the values `y_i`
    are the normalized coordinates such that all `y_i` are integers and
    `\gcd(y_1, y_2,\dots, y_n) = 1`.

    ALGORITHM:

    Main idea behind the algorithm is to find points modulo primes
    and then reconstruct them using chinese remainder theorem.
    We find modulo primes parallelly and then lift them and apply
    LLL in parallel.

    For the algorithm to work correctly, sufficient primes need
    to be present, these are calculated using the bound given in
    this([Hutz2015]_) paper.

    INPUT:

    - ``X`` -- a scheme with ambient space defined over projective space

    - ``bound`` -- positive integer bound

    OUTPUT:

    A list containing the projective rational points of ``X`` of height
    up to ``bound``, sorted

    EXAMPLES::

        sage: from sage.schemes.projective.projective_rational_point import sieve
        sage: P.<x,y,z,q> = ProjectiveSpace(QQ, 3)
        sage: Y = P.subscheme([x^2 - 3^2*y^2 + z*q, x + z + 4*q])
        sage: sorted(sieve(Y, 12))              # long time                             # needs sage.libs.singular
        [(-4 : -4/3 : 0 : 1), (-4 : 4/3 : 0 : 1),
         (-1 : -1/3 : 1 : 0), (-1 : 1/3 : 1 : 0)]

    ::

        sage: from sage.schemes.projective.projective_rational_point import sieve
        sage: E = EllipticCurve('37a')                                                  # needs sage.schemes
        sage: sorted(sieve(E, 14))              # long time                             # needs sage.libs.singular sage.schemes
        [(0 : 1 : 0),
         (-1 : -1 : 1),
         (-1 : 0 : 1),
         (0 : -1 : 1),
         (0 : 0 : 1),
         (1/4 : -5/8 : 1),
         (1/4 : -3/8 : 1),
         (1 : -1 : 1),
         (1 : 0 : 1),
         (2 : -3 : 1),
         (2 : 2 : 1),
         (6 : 14 : 1)]

    TESTS:

    Algorithm works even if coefficients are fraction::

        sage: from sage.schemes.projective.projective_rational_point import sieve
        sage: P.<x,y,z> = ProjectiveSpace(2, QQ)
        sage: X = P.subscheme(3*x - 3/2*y)
        sage: sieve(X, 3)                                                               # needs sage.libs.singular
        [(-1 : -2 : 1), (-1/2 : -1 : 1), (-1/3 : -2/3 : 1), (0 : 0 : 1),
         (1/3 : 2/3 : 1), (1/2 : 1 : 0), (1/2 : 1 : 1), (1 : 2 : 1)]
    """
    if bound < 1: # no projective rational point with height less than 1
        return []

    modulo_points = [] # list to store point modulo primes
    len_modulo_points = [] # stores number of points with respect to each prime
    primes_list = [] # list of good primes

    X.normalize_defining_polynomials()

    P = X.ambient_space()
    N = P.dimension()
    dim_scheme = X.dimension()

    # bound as per preposition - 4, in preperiodic points paper
    B = RR(2**(N/4+1)*bound**2*(N+1).sqrt())

    m = [0 for _ in range(N + 1)]

    def sufficient_primes(x):
        r"""
        Return a list of primes whose product is > `x`.
        """
        small_primes = [2,3]
        prod_primes = 6

        while prod_primes < x:
            p = next_prime(small_primes[-1])
            small_primes.append(p)
            prod_primes *= p
        return small_primes

    def good_primes(B):
        r"""
        Given the bound returns the prime whose product is greater than ``B``
        and which would take least amount of time to run main sieve algorithm

        Complexity of finding points modulo primes is assumed to be N^2 * P_max^{N}.
        Complexity of lifting points and LLL() function is assumed to
        be close to N^5 * (alpha^dim_scheme / P_max).
        where alpha is product of all primes, and P_max is largest prime in list.
        """

        M = {}  # stores optimal list of primes, corresponding to list size
        small_primes = sufficient_primes(B)
        max_length = len(small_primes)
        M[max_length] = small_primes
        current_count = max_length - 1

        while current_count > 1:
            current_list = [] # stores prime which are bigger than least
            updated_list = []
            best_list = []

            least = (RR(B)**(1.00/current_count)).floor()
            for i in range(current_count):
                current_list.append(next_prime(least))
                least = current_list[-1]
            # improving list of primes by taking prime less than least
            # this part of algorithm is used to centralize primes around `least`
            prod_prime = prod(current_list)
            least = current_list[0]
            while least != 2 and prod_prime > B and len(updated_list) < current_count:
                best_list = updated_list + current_list[:current_count - len(updated_list)]
                updated_list.append(previous_prime(least))
                least = updated_list[-1]

                removed_prime = current_list[current_count - len(updated_list)]
                prod_prime = (prod_prime * least) / removed_prime

            M[current_count] = sorted(best_list)
            current_count = current_count - 1

        best_size = 2
        best_time = (N**2)*M[2][-1]**(N) + (N**5 * RR(prod(M[2])**dim_scheme / M[2][-1]) )
        for i in range(2, max_length + 1):
            current_time = (N**2)*M[i][-1]**(N) + (N**5 * RR(prod(M[i])**dim_scheme / M[i][-1]) )
            if current_time < best_time:
                best_size = i
                best_time = current_time

        return M[best_size]

    def parallel_function(X, p):
        r"""
        Function used in parallel computation, computes a list of
        all rational points in modulo ring.
        """
        Xp = X.change_ring(GF(p))
        L = Xp.rational_points()

        return [list(_) for _ in L]

    def points_modulo_primes(X, primes):
        r"""
        Return a list of rational points modulo all `p` in primes,
        computed parallelly.
        """
        normalized_input = []
        for p in primes_list:
            normalized_input.append(((X, p, ), {}))
        p_iter = p_iter_fork(ncpus())

        points_pair = list(p_iter(parallel_function, normalized_input))
        points_pair.sort()
        return [pair[1] for pair in points_pair]

    def parallel_function_combination(point_p_max):
        r"""
        Function used in parallel computation, computes rational
        points lifted.
        """
        rat_points = set()
        for tupl in xmrange(len_modulo_points):
            point = []
            for k in range(N + 1):
                # lift all coordinates of given point using chinese remainder theorem
                L = [modulo_points[j][tupl[j]][k].lift() for j in range(len_primes - 1)]
                L.append(point_p_max[k].lift())
                point.append( crt(L, primes_list) )

            for i in range(N+1):
                m[i] = point[i]

            M = matrix(ZZ, N+2, N+1, m)
            A = M.LLL()
            point = list(A[1])

            # check if all coordinates of this point satisfy height bound
            bound_satisfied = True
            for coordinate in point:
                if coordinate.abs() > bound:
                    bound_satisfied = False
                    break
            if not bound_satisfied:
                continue

            try:
                pt = X(list(A[1]))
            except TypeError:
                pass
            else:
                rat_points.add(pt)

        return [list(_) for _ in rat_points]

    def lift_all_points():
        r"""
        Return list of all rational points lifted parallelly.
        """
        normalized_input = []
        points = modulo_points.pop() # remove the list of points corresponding to largest prime
        len_modulo_points.pop()

        for point in points:
            normalized_input.append(( (point, ), {}))
        p_iter = p_iter_fork(ncpus())
        points_satisfying = list(p_iter(parallel_function_combination, normalized_input))

        lifted_points = set()
        for pair in points_satisfying:
            L = pair[1]
            for point in L:
                lifted_points.add(X(tuple(point)))

        return list(lifted_points)

    # start of main algorithm

    primes_list = good_primes(B.ceil())

    modulo_points = points_modulo_primes(X, primes_list)
    len_modulo_points = [len(pt) for pt in modulo_points]
    len_primes = len(primes_list)
    prod_primes = prod(primes_list)

    # stores final result

    for i in range(N + 1):
        w = [0 for _ in range(N + 1)]
        w[i] = prod_primes
        m.extend(w)

    rat_points = lift_all_points()

    return sorted(rat_points)
