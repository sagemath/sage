r"""
Enumeration of rational points on product projective schemes

Naive algorithms for enumerating rational points over `\QQ`, number fields or
finite fields over general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.
    (For instance using an affine function for a product projective scheme
    or a finite field function for a scheme defined over an infinite field.)

EXAMPLES:

Product Projective, over `\QQ`::

    sage: PP.<x,y,z> = ProductProjectiveSpaces([1,0], QQ)
    sage: from sage.schemes.product_projective.rational_point import \
            enum_product_projective_rational_field
    sage: enum_product_projective_rational_field(PP, 3)
    [(-3 : 1 , 1), (-2 : 1 , 1), (-3/2 : 1 , 1),
     (-1 : 1 , 1), (-2/3 : 1 , 1), (-1/2 : 1 , 1),
     (-1/3 : 1 , 1), (0 : 1 , 1), (1/3 : 1 , 1),
     (1/2 : 1 , 1), (2/3 : 1 , 1), (1 : 0 , 1),
     (1 : 1 , 1), (3/2 : 1 , 1), (2 : 1 , 1),
     (3 : 1 , 1)]

Product projective over finite field::

    sage: P1.<x,y,a,b> = ProductProjectiveSpaces([1, 1], GF(7))
    sage: X = P1.subscheme([2*x + 3*y])
    sage: from sage.schemes.product_projective.rational_point import \
            enum_product_projective_finite_field
    sage: enum_product_projective_finite_field(X)
    [(2 : 1 , 0 : 1), (2 : 1 , 1 : 0), (2 : 1 , 1 : 1),
     (2 : 1 , 2 : 1), (2 : 1 , 3 : 1), (2 : 1 , 4 : 1),
     (2 : 1 , 5 : 1), (2 : 1 , 6 : 1)]

AUTHORS:

- Volker Braun and Ben Hutz (2014): initial version

- Raghukul Raman (2018): code cleanup and added support for rational fields
"""

# ****************************************************************************
# Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#                    Ben Hutz <bn4941@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.generic.scheme import Scheme
from sage.schemes.product_projective.space import ProductProjectiveSpaces_ring
from sage.misc.lazy_import import lazy_import
from sage.misc.mrange import xmrange
from sage.misc.misc_c import prod
from sage.arith.misc import next_prime, previous_prime, crt
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.parallel.ncpus import ncpus
from sage.parallel.use_fork import p_iter_fork

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.rings.real_mpfr', 'RR')


def enum_product_projective_rational_field(X, B):
    r"""
    Enumerate projective, rational points on scheme ``X`` of height up to
    bound ``B``.

    INPUT:

    - ``X`` -- a scheme or set of abstract rational points of a scheme

    - ``B`` -- positive integer bound

    OUTPUT:

    A list containing the product projective points of ``X`` of height up
    to ``B``, sorted.

    EXAMPLES::

        sage: PP.<x0,x1,x2,x3,x4> = ProductProjectiveSpaces([1, 2], QQ)
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_rational_field
        sage: enum_product_projective_rational_field(PP, 1)
        [(-1 : 1 , -1 : -1 : 1), (-1 : 1 , -1 : 0 : 1), (-1 : 1 , -1 : 1 : 0),
         (-1 : 1 , -1 : 1 : 1), (-1 : 1 , 0 : -1 : 1), (-1 : 1 , 0 : 0 : 1),
         (-1 : 1 , 0 : 1 : 0), (-1 : 1 , 0 : 1 : 1), (-1 : 1 , 1 : -1 : 1),
         (-1 : 1 , 1 : 0 : 0), (-1 : 1 , 1 : 0 : 1), (-1 : 1 , 1 : 1 : 0),
         (-1 : 1 , 1 : 1 : 1), (0 : 1 , -1 : -1 : 1), (0 : 1 , -1 : 0 : 1),
         (0 : 1 , -1 : 1 : 0), (0 : 1 , -1 : 1 : 1), (0 : 1 , 0 : -1 : 1),
         (0 : 1 , 0 : 0 : 1), (0 : 1 , 0 : 1 : 0), (0 : 1 , 0 : 1 : 1),
         (0 : 1 , 1 : -1 : 1), (0 : 1 , 1 : 0 : 0), (0 : 1 , 1 : 0 : 1),
         (0 : 1 , 1 : 1 : 0), (0 : 1 , 1 : 1 : 1), (1 : 0 , -1 : -1 : 1),
         (1 : 0 , -1 : 0 : 1), (1 : 0 , -1 : 1 : 0), (1 : 0 , -1 : 1 : 1),
         (1 : 0 , 0 : -1 : 1), (1 : 0 , 0 : 0 : 1), (1 : 0 , 0 : 1 : 0),
         (1 : 0 , 0 : 1 : 1), (1 : 0 , 1 : -1 : 1), (1 : 0 , 1 : 0 : 0),
         (1 : 0 , 1 : 0 : 1), (1 : 0 , 1 : 1 : 0), (1 : 0 , 1 : 1 : 1),
         (1 : 1 , -1 : -1 : 1), (1 : 1 , -1 : 0 : 1), (1 : 1 , -1 : 1 : 0),
         (1 : 1 , -1 : 1 : 1), (1 : 1 , 0 : -1 : 1), (1 : 1 , 0 : 0 : 1),
         (1 : 1 , 0 : 1 : 0), (1 : 1 , 0 : 1 : 1), (1 : 1 , 1 : -1 : 1),
         (1 : 1 , 1 : 0 : 0), (1 : 1 , 1 : 0 : 1), (1 : 1 , 1 : 1 : 0),
         (1 : 1 , 1 : 1 : 1)]

    ::

        sage: PP.<x,y,z,u,v> = ProductProjectiveSpaces([2, 1], QQ)
        sage: X = PP.subscheme([x^2 + x*y + y*z, u*u - v*u])
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_rational_field
        sage: enum_product_projective_rational_field(X, 4)
        [(-2 : 4 : 1 , 0 : 1), (-2 : 4 : 1 , 1 : 1), (-1 : 1 : 0 , 0 : 1),
         (-1 : 1 : 0 , 1 : 1), (-2/3 : -4/3 : 1 , 0 : 1), (-2/3 : -4/3 : 1 , 1 : 1),
         (-1/2 : -1/2 : 1 , 0 : 1), (-1/2 : -1/2 : 1 , 1 : 1),
         (0 : 0 : 1 , 0 : 1), (0 : 0 : 1 , 1 : 1), (0 : 1 : 0 , 0 : 1),
         (0 : 1 : 0 , 1 : 1), (1 : -1/2 : 1 , 0 : 1), (1 : -1/2 : 1 , 1 : 1)]
    """
    if isinstance(X, Scheme):
        if not isinstance(X.ambient_space(), ProductProjectiveSpaces_ring):
            raise TypeError("ambient space must be product of projective space over the rational field")
        X = X(X.base_ring())
    else:
        if not isinstance(X.codomain().ambient_space(), ProductProjectiveSpaces_ring):
            raise TypeError("codomain must be product of projective space over the rational field")

    R = X.codomain().ambient_space()
    m = R.num_components()
    iters = [ R[i].points_of_bounded_height(bound=B) for i in range(m) ]
    dim = [R[i].dimension_relative() + 1 for i in range(m)]

    dim_prefix = [0, dim[0]] # prefixes dim list
    for i in range(1, len(dim)):
        dim_prefix.append(dim_prefix[i] + dim[i])

    pts = []
    P = []
    for i in range(m):
        pt = next(iters[i])
        for j in range(dim[i]):
            P.append(pt[j]) # initial value of P

    try: # add the initial point
        pts.append(X(P))
    except TypeError:
        pass

    i = 0
    while i < m:
        try:
            pt = next(iters[i])
            for j in range(dim[i]):
                P[dim_prefix[i] + j] = pt[j]
            try:
                pts.append(X(P))
            except TypeError:
                pass
            i = 0
        except StopIteration:
            iters[i] = R[i].points_of_bounded_height(bound=B)
            pt = next(iters[i]) # reset
            for j in range(dim[i]):
                P[dim_prefix[i] + j] = pt[j]
            i += 1
    pts.sort()

    return pts

def enum_product_projective_number_field(X, **kwds):
    r"""
    Enumerates product projective points on scheme ``X`` defined over a number field.

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

    OUTPUT:

    A list containing the product projective points of ``X`` of
    absolute height up to ``B``, sorted.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 2, 'v')
        sage: PP.<x,y,z,w> = ProductProjectiveSpaces([1, 1], K)
        sage: X = PP.subscheme([x^2 + 2*y^2])
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_number_field
        sage: enum_product_projective_number_field(X, bound=1.5)
        [(-v : 1 , -1 : 1), (-v : 1 , -v : 1), (-v : 1 , -1/2*v : 1),
         (-v : 1 , 0 : 1), (-v : 1 , 1/2*v : 1), (-v : 1 , v : 1),
         (-v : 1 , 1 : 0), (-v : 1 , 1 : 1), (v : 1 , -1 : 1),
         (v : 1 , -v : 1), (v : 1 , -1/2*v : 1), (v : 1 , 0 : 1),
         (v : 1 , 1/2*v : 1), (v : 1 , v : 1), (v : 1 , 1 : 0),
         (v : 1 , 1 : 1)]
    """
    B = kwds.pop('bound')
    tol = kwds.pop('tolerance', 1e-2)
    prec = kwds.pop('precision', 53)

    if isinstance(X, Scheme):
        if not isinstance(X.ambient_space(), ProductProjectiveSpaces_ring):
            raise TypeError("ambient space must be product of projective space over the rational field")
        X = X(X.base_ring())
    else:
        if not isinstance(X.codomain().ambient_space(), ProductProjectiveSpaces_ring):
            raise TypeError("codomain must be product of projective space over the rational field")

    R = X.codomain().ambient_space()

    pts = []

    for P in R.points_of_bounded_height(bound=B, tolerance=tol, precision=prec):
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()
    return pts

def enum_product_projective_finite_field(X):
    r"""
    Enumerates projective points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -- a scheme defined over a finite field or a set of abstract
      rational points of such a scheme

    OUTPUT:

    A list containing the projective points of ``X`` over the finite field,
    sorted.

    EXAMPLES::

        sage: PP.<x,y,z,w> = ProductProjectiveSpaces([1, 1], GF(3))
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_finite_field
        sage: enum_product_projective_finite_field(PP)
        [(0 : 1 , 0 : 1), (0 : 1 , 1 : 0), (0 : 1 , 1 : 1),
         (0 : 1 , 2 : 1), (1 : 0 , 0 : 1), (1 : 0 , 1 : 0),
         (1 : 0 , 1 : 1), (1 : 0 , 2 : 1), (1 : 1 , 0 : 1),
         (1 : 1 , 1 : 0), (1 : 1 , 1 : 1), (1 : 1 , 2 : 1),
         (2 : 1 , 0 : 1), (2 : 1 , 1 : 0), (2 : 1 , 1 : 1),
         (2 : 1 , 2 : 1)]

    ::

        sage: PP.<x0,x1,x2,x3> = ProductProjectiveSpaces([1, 1], GF(17))
        sage: X = PP.subscheme([x0^2 + 2*x1^2])
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_finite_field
        sage: len(enum_product_projective_finite_field(X))
        36
    """
    if isinstance(X, Scheme):
        if not isinstance(X.ambient_space(), ProductProjectiveSpaces_ring):
            raise TypeError("ambient space must be product of projective space over the rational field")
        X = X(X.base_ring())
    else:
        if not isinstance(X.codomain().ambient_space(), ProductProjectiveSpaces_ring):
            raise TypeError("codomain must be product of projective space over the rational field")

    R = X.codomain().ambient_space()
    pts = []

    for P in R.rational_points():
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()

    return pts


def sieve(X, bound):
    r"""
    Return the list of all rational points on scheme
    ``X`` of height up to ``bound``.

    ALGORITHM:

    Main idea behind the algorithm is to find points modulo primes
    and then reconstruct them using chinese remainder theorem.
    We compute the points modulo primes parallelly and then lift
    them via chinese remainder theorem in parallel. The LLL reduction
    algorithm is applied to each component of the points, and finally
    the result is merged and converted to a point on the subscheme.

    For the algorithm to work correctly, sufficient primes need
    to be chosen, these are determined using the bounds dependent
    on the bound given in [Hutz2015]_.

    INPUT:

    - ``X`` -- a scheme with ambient space defined over a product of projective spaces

    - ``bound`` -- positive integer bound

    OUTPUT:

    A list containing the rational points of ``X`` of height
    up to ``bound``, sorted

    EXAMPLES::

        sage: from sage.schemes.product_projective.rational_point import sieve
        sage: PP.<x,y,z,u,v> = ProductProjectiveSpaces([2, 1], QQ)
        sage: X = PP.subscheme([x^2 + y^2 - x*z, u*u - v*u])
        sage: sieve(X, 2)                                                               # needs sage.libs.singular
        [(0 : 0 : 1 , 0 : 1), (0 : 0 : 1 , 1 : 1), (1/2 : -1/2 : 1 , 0 : 1),
         (1/2 : -1/2 : 1 , 1 : 1), (1/2 : 1/2 : 1 , 0 : 1), (1/2 : 1/2 : 1 , 1 : 1),
         (1 : 0 : 1 , 0 : 1), (1 : 0 : 1 , 1 : 1)]
    """

    if bound < 1:
        return []

    modulo_points = []
    len_modulo_points = []
    primes_list = []

    P = X.ambient_space()
    N = P.ngens()
    dim_scheme = X.dimension()

    num_comp = P.num_components()
    comp_dim_relative = [P[i].dimension_relative() + 1 for i in range(num_comp)]

    dim_prefix = [0, comp_dim_relative[0]] # prefixes dim list
    for i in range(1, len(comp_dim_relative)):
        dim_prefix.append(dim_prefix[i] + comp_dim_relative[i])

    dim_max = max(P[i].dimension() for i in range(num_comp))
    B = RR(2**(dim_max/4+1)*bound**2*(dim_max+1).sqrt())
    m = []

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
        Given the bound, returns the primes whose product is greater than ``B``
        and which would take the least amount of time to run the main sieve algorithm

        Complexity of finding points modulo primes is assumed to be N^2 * P_max^{N}.
        Complexity of lifting points and the LLL() function is assumed to
        be close to (dim_max^5) * (alpha / P_max)^dim_scheme.
        where alpha is the product of all primes, P_max is the largest prime in
        the list, dim_max is the max dimension of all components, and N is the dimension
        of the ambient space.
        """

        M = {}  # stores optimal list of primes, corresponding to list size
        small_primes = sufficient_primes(B)
        max_length = len(small_primes)
        M[max_length] = small_primes
        current_count = max_length - 1
        dim = X.ambient_space().dimension()

        while current_count > 1:
            current_list = [] # stores prime which are bigger than least
            updated_list = []
            best_list = []

            least = (RR(B)**(1.00/current_count)).floor()
            for i in range(current_count):
                current_list.append(next_prime(least))
                least = current_list[-1]
            # improving list of primes by taking primes less than least
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
        best_time = (dim**2)*M[2][-1]**(dim) + (dim_max**5 * (prod(M[2])/M[2][-1])**dim_scheme)
        for i in range(2, max_length + 1):
            current_time = (dim**2)*M[i][-1]**(dim) + (dim_max**5 * (prod(M[i])/M[i][-1])**dim_scheme)
            if current_time < best_time:
                best_size = i
                best_time = current_time

        return M[best_size]

    def parallel_function(X, p):
        r"""
        Function used in parallel computation, computes a list of
        all rational points in modulo p.
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
        modulo_points = []
        for pair in points_pair:
            modulo_points.append(pair[1])

        return modulo_points

    def parallel_function_combination(point_p_max):
        r"""
        Function used in parallel computation, computes rational
        points lifted.
        """
        rat_points = set()
        for tupl in xmrange(len_modulo_points):
            point = []
            for k in range(N):
                # lift all coordinates of given point using chinese remainder theorem
                L = [modulo_points[j][tupl[j]][k].lift() for j in range(len_primes - 1)]
                L.append(point_p_max[k].lift())
                point.append(crt(L, primes_list))

            for i in range(num_comp):
                for j in range(comp_dim_relative[i]):
                    m[i][j] = point[dim_prefix[i] + j]

            # generating matrix to compute LLL reduction for each component
            M = [matrix(ZZ, comp_dim_relative[i] + 1,
                        comp_dim_relative[i], m[i])
                 for i in range(num_comp)]
            A = [M[i].LLL() for i in range(num_comp)]
            point = []
            for i in range(num_comp):
                point.extend(A[i][1])

            # check if all coordinates of this point satisfy height bound
            bound_satisfied = True
            for coordinate in point:
                if coordinate.abs() > bound:
                    bound_satisfied = False
                    break
            if not bound_satisfied:
                continue

            try:
                # checks if this point lies on X or not
                rat_points.add(X(point))
            except (TypeError, ValueError):
                pass

        return [list(pt) for pt in rat_points]

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

    primes_list = good_primes(B.ceil())
    modulo_points = points_modulo_primes(X, primes_list)
    len_modulo_points = [len(_) for _ in modulo_points]
    len_primes = len(primes_list)
    prod_primes = prod(primes_list)

    # construct m for each component projective space
    for i in range(num_comp):
        dim = comp_dim_relative[i]
        temp = [0 for _ in range(dim)]
        for j in range(dim):
            w = [0 for _ in range(dim)]
            w[j] = prod_primes
            temp.extend(w)
        m.append(temp)

    rat_points = lift_all_points()

    return sorted(rat_points)
