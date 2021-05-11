r"""
Elements of bounded height in number fields

Sage functions to list all elements of a given number field with height less
than a specified bound.

AUTHORS:

- John Doyle (2013): initial version

- David Krumm (2013): initial version

- TJ Combs (2018): added Doyle-Krumm algorithm - 4

- Raghukul Raman (2018): added Doyle-Krumm algorithm - 4

REFERENCES:

- [DK2013]

"""
# ****************************************************************************
#       Copyright (C) 2013 John Doyle and David Krumm
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from itertools import product
from sage.rings.real_mpfr import RealField
from sage.rings.number_field.unit_group import UnitGroup
from sage.modules.free_module_element import vector
from sage.matrix.constructor import column_matrix
from sage.rings.rational_field import QQ
from sage.rings.all import RR, Infinity
from sage.geometry.polyhedron.constructor import Polyhedron


def bdd_norm_pr_gens_iq(K, norm_list):
    r"""
    Compute generators for all principal ideals in an imaginary quadratic field
    `K` whose norms are in ``norm_list``.

    The only keys for the output dictionary are integers n appearing in
    ``norm_list``.

    The function will only be called with `K` an imaginary quadratic field.

    The function will return a dictionary for other number fields, but it may be
    incorrect.

    INPUT:

    - `K` -- an imaginary quadratic number field

    - ``norm_list`` -- a list of positive integers

    OUTPUT:

    - a dictionary of number field elements, keyed by norm

    EXAMPLES:

    In `QQ(i)`, there is one principal ideal of norm 4, two principal ideals of
    norm 5, but no principal ideals of norm 7::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_gens_iq
        sage: K.<g> = NumberField(x^2 + 1)
        sage: L = range(10)
        sage: bdd_pr_ideals = bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals[4]
        [2]
        sage: bdd_pr_ideals[5]
        [-g - 2, -g + 2]
        sage: bdd_pr_ideals[7]
        []

    There are no ideals in the ring of integers with negative norm::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_gens_iq
        sage: K.<g> = NumberField(x^2 + 10)
        sage: L = range(-5,-1)
        sage: bdd_pr_ideals = bdd_norm_pr_gens_iq(K,L)
        sage: bdd_pr_ideals
        {-5: [], -4: [], -3: [], -2: []}

    Calling a key that is not in the input ``norm_list`` raises a KeyError::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_gens_iq
        sage: K.<g> = NumberField(x^2 + 20)
        sage: L = range(100)
        sage: bdd_pr_ideals = bdd_norm_pr_gens_iq(K, L)
        sage: bdd_pr_ideals[100]
        Traceback (most recent call last):
        ...
        KeyError: 100
    """
    return {n: K.elements_of_norm(n) for n in norm_list}


def bdd_height_iq(K, height_bound):
    r"""
    Compute all elements in the imaginary quadratic field `K` which have
    relative multiplicative height at most ``height_bound``.

    The function will only be called with `K` an imaginary quadratic field.

    If called with `K` not an imaginary quadratic, the function will likely
    yield incorrect output.

    ALGORITHM:

    This is an implementation of Algorithm 5 in [DK2013]_.

    INPUT:

    - `K` -- an imaginary quadratic number field

    - ``height_bound`` -- a real number

    OUTPUT:

    - an iterator of number field elements

    EXAMPLES::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + 191)
        sage: for t in bdd_height_iq(K,8):
        ....:     print(exp(2*t.global_height()))
        1.00000000000000
        1.00000000000000
        1.00000000000000
        4.00000000000000
        4.00000000000000
        4.00000000000000
        4.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000
        8.00000000000000

    There are 175 elements of height at most 10 in `QQ(\sqrt(-3))`::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + 3)
        sage: len(list(bdd_height_iq(K,10)))
        175

    The only elements of multiplicative height 1 in a number field are 0 and
    the roots of unity::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + x + 1)
        sage: list(bdd_height_iq(K,1))
        [0, a + 1, a, -1, -a - 1, -a, 1]

    A number field has no elements of multiplicative height less than 1::

        sage: from sage.rings.number_field.bdd_height import bdd_height_iq
        sage: K.<a> = NumberField(x^2 + 5)
        sage: list(bdd_height_iq(K,0.9))
        []
    """
    if height_bound < 1:
        return
    yield K(0)
    roots_of_unity = K.roots_of_unity()
    for zeta in roots_of_unity:
        yield zeta

    # Get a complete set of ideal class representatives
    class_group_reps = []
    class_group_rep_norms = []
    for c in K.class_group():
        a = c.ideal()
        class_group_reps.append(a)
        class_group_rep_norms.append(a.norm())
    class_number = len(class_group_reps)

    # Find principal ideals of bounded norm
    possible_norm_set = set([])
    for n in range(class_number):
        for m in range(1, int(height_bound + 1)):
            possible_norm_set.add(m*class_group_rep_norms[n])
    bdd_ideals = bdd_norm_pr_gens_iq(K, possible_norm_set)

    # Distribute the principal ideals
    generator_lists = []
    for n in range(class_number):
        this_ideal = class_group_reps[n]
        this_ideal_norm = class_group_rep_norms[n]
        gens = []
        for i in range(1, int(height_bound + 1)):
            for g in bdd_ideals[i*this_ideal_norm]:
                if g in this_ideal:
                    gens.append(g)
        generator_lists.append(gens)

    # Build all the output numbers
    for n in range(class_number):
        gens = generator_lists[n]
        s = len(gens)
        for i in range(s):
            for j in range(i + 1, s):
                if K.ideal(gens[i], gens[j]) == class_group_reps[n]:
                    new_number = gens[i]/gens[j]
                    for zeta in roots_of_unity:
                        yield zeta * new_number
                        yield zeta / new_number


def bdd_norm_pr_ideal_gens(K, norm_list):
    r"""
    Compute generators for all principal ideals in a number field `K` whose
    norms are in ``norm_list``.

    INPUT:

    - `K` -- a number field

    - ``norm_list`` -- a list of positive integers

    OUTPUT:

    - a dictionary of number field elements, keyed by norm

    EXAMPLES:

    There is only one principal ideal of norm 1, and it is generated by the
    element 1::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
        sage: K.<g> = QuadraticField(101)
        sage: bdd_norm_pr_ideal_gens(K, [1])
        {1: [1]}

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
        sage: K.<g> = QuadraticField(123)
        sage: bdd_norm_pr_ideal_gens(K, range(5))
        {0: [0], 1: [1], 2: [-g - 11], 3: [], 4: [2]}

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
        sage: K.<g> = NumberField(x^5 - x + 19)
        sage: b = bdd_norm_pr_ideal_gens(K, range(30))
        sage: key = ZZ(28)
        sage: b[key]
        [157*g^4 - 139*g^3 - 369*g^2 + 848*g + 158, g^4 + g^3 - g - 7]

    """
    negative_norm_units = K.elements_of_norm(-1)
    gens = {}
    if not negative_norm_units:
        for n in norm_list:
            if not n:
                gens[n] = [K.zero()]
            else:
                gens[n] = K.elements_of_norm(n) + K.elements_of_norm(-n)
    else:
        for n in norm_list:
            gens[n] = K.elements_of_norm(n)
    return gens


def integer_points_in_polytope(matrix, interval_radius):
    r"""
    Return the set of integer points in the polytope obtained by acting on a
    cube by a linear transformation.

    Given an r-by-r matrix ``matrix`` and a real number ``interval_radius``,
    this function finds all integer lattice points in the polytope obtained by
    transforming the cube [-interval_radius,interval_radius]^r via the linear
    map induced by ``matrix``.

    INPUT:

    - ``matrix`` -- a square matrix of real numbers

    - ``interval_radius`` -- a real number

    OUTPUT:

    - a list of tuples of integers

    EXAMPLES:

    Stretch the interval [-1,1] by a factor of 2 and find the integers in the
    resulting interval::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([2])
        sage: r = 1
        sage: integer_points_in_polytope(m,r)
        [(-2), (-1), (0), (1), (2)]

    Integer points inside a parallelogram::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([[1, 2],[3, 4]])
        sage: r = RealField()(1.3)
        sage: integer_points_in_polytope(m,r)
        [(-3, -7), (-2, -5), (-2, -4), (-1, -3), (-1, -2), (-1, -1), (0, -1), (0, 0), (0, 1), (1, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 7)]

    Integer points inside a parallelepiped::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([[1.2,3.7,0.2],[-5.3,-.43,3],[1.2,4.7,-2.1]])
        sage: r = 2.2
        sage: L = integer_points_in_polytope(m,r)
        sage: len(L)
        4143

    If ``interval_radius`` is 0, the output should include only the zero tuple::

        sage: from sage.rings.number_field.bdd_height import integer_points_in_polytope
        sage: m = matrix([[1,2,3,7],[4,5,6,2],[7,8,9,3],[0,3,4,5]])
        sage: integer_points_in_polytope(m,0)
        [(0, 0, 0, 0)]
    """
    T = matrix
    d = interval_radius
    r = T.nrows()

    # Find the vertices of the given box
    box_vertices = [vector(x) for x in product([-d, d], repeat=r)]

    # Transform the vertices
    T_trans = T.transpose()
    transformed_vertices = [v * T_trans for v in box_vertices]

    # Create polyhedron from transformed vertices and find integer points inside
    return list(Polyhedron(transformed_vertices, base_ring=QQ).integral_points())


def bdd_height(K, height_bound, tolerance=1e-2, precision=53):
    r"""
    Compute all elements in the number field `K` which have relative
    multiplicative height at most ``height_bound``.

    The function can only be called for number fields `K` with positive unit
    rank. An error will occur if `K` is `QQ` or an imaginary quadratic field.

    This algorithm computes 2 lists: L containing elements x in `K` such that
    H_k(x) <= B, and a list L' containing elements x in `K` that, due to
    floating point issues,
    may be slightly larger then the bound. This can be controlled
    by lowering the tolerance.

    In current implementation both lists (L,L') are merged and returned in
    form of iterator.

    ALGORITHM:

    This is an implementation of the revised algorithm (Algorithm 4) in
    [DK2013]_.

    INPUT:

    - ``height_bound`` -- real number

    - ``tolerance`` -- (default: 0.01) a rational number in (0,1]

    - ``precision`` -- (default: 53) positive integer

    OUTPUT:

    - an iterator of number field elements

    EXAMPLES:

    There are no elements of negative height::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = NumberField(x^5 - x + 7)
        sage: list(bdd_height(K,-3))
        []

    The only nonzero elements of height 1 are the roots of unity::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = QuadraticField(3)
        sage: list(bdd_height(K,1))
        [0, -1, 1]

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = QuadraticField(36865)
        sage: len(list(bdd_height(K,101))) # long time (4 s)
        131

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = NumberField(x^6 + 2)
        sage: len(list(bdd_height(K,60))) # long time (5 s)
        1899

    ::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<g> = NumberField(x^4 - x^3 - 3*x^2 + x + 1)
        sage: len(list(bdd_height(K,10)))
        99

    TESTS:

    Check that :trac:`22771` is fixed::

        sage: from sage.rings.number_field.bdd_height import bdd_height
        sage: K.<v> = NumberField(x^3 + x + 1)
        sage: len(list(bdd_height(K,3)))
        23
    """
    # global values, used in internal function
    B = height_bound
    theta = tolerance
    if B < 1:
        return
    embeddings = K.places(prec=precision)
    O_K = K.ring_of_integers()
    r1, r2 = K.signature()
    r = r1 + r2 - 1
    RF = RealField(precision)
    lambda_gens_approx = {}
    class_group_rep_norm_log_approx = []
    unit_log_dict = {}

    def rational_in(x, y):
        r"""
        Compute a rational number q, such that x<q<y using Archimedes' axiom
        """
        z = y - x
        if z == 0:
            n = 1
        else:
            n = RR(1/z).ceil() + 1
        if RR(n*y).ceil() is n*y:  # WHAT !?
            m = n*y - 1
        else:
            m = RR(n*y).floor()
        return m / n

    def delta_approximation(x, delta):
        r"""
        Compute a rational number in range (x-delta, x+delta)
        """
        return rational_in(x - delta, x + delta)

    def vector_delta_approximation(v, delta):
        r"""
        Compute a rational vector w=(w1, ..., wn)
        such that |vi-wi|<delta for all i in [1, n]
        """
        return [delta_approximation(vi, delta) for vi in v]

    def log_map(number):
        r"""
        Compute the image of an element of `K` under the logarithmic map.
        """
        x = number
        x_logs = []
        for i in range(r1):
            sigma = embeddings[i]  # real embeddings
            x_logs.append(sigma(x).abs().log())
        for i in range(r1, r + 1):
            tau = embeddings[i]  # Complex embeddings
            x_logs.append(2 * tau(x).abs().log())
        return vector(x_logs)

    def log_height_for_generators_approx(alpha, beta, Lambda):
        r"""
        Compute the rational approximation of logarithmic height function.
        Return a lambda approximation h_K(alpha/beta)
        """
        delta = Lambda / (r + 2)
        norm_log = delta_approximation(RR(O_K.ideal(alpha, beta).norm()).log(), delta)
        log_ga = vector_delta_approximation(log_map(alpha), delta)
        log_gb = vector_delta_approximation(log_map(beta), delta)
        arch_sum = sum([max(log_ga[k], log_gb[k]) for k in range(r + 1)])
        return (arch_sum - norm_log)

    def packet_height(n, pair, u):
        r"""
        Compute the height of the element of `K` encoded by a given packet.
        """
        gens = generator_lists[n]
        i = pair[0]
        j = pair[1]
        Log_gi = lambda_gens_approx[gens[i]]
        Log_gj = lambda_gens_approx[gens[j]]
        Log_u_gi = vector(Log_gi) + unit_log_dict[u]
        arch_sum = sum([max(Log_u_gi[k], Log_gj[k]) for k in range(r + 1)])
        return (arch_sum - class_group_rep_norm_log_approx[n])

    # Step 1
    # Computes ideal class representative and their rational approx norm
    t = theta / (3*B)
    delta_1 = t / (6*r+12)

    class_group_reps = []
    class_group_rep_norms = []

    for c in K.class_group():
        a = c.ideal()
        a_norm = a.norm()
        log_norm = RF(a_norm).log()
        log_norm_approx = delta_approximation(log_norm, delta_1)
        class_group_reps.append(a)
        class_group_rep_norms.append(a_norm)
        class_group_rep_norm_log_approx.append(log_norm_approx)
    class_number = len(class_group_reps)

    # Step 2
    # Find generators for principal ideals of bounded norm
    possible_norm_set = set([])
    for n in range(class_number):
        for m in range(1, (B + 1).ceil()):
            possible_norm_set.add(m * class_group_rep_norms[n])
    bdd_ideals = bdd_norm_pr_ideal_gens(K, possible_norm_set)

    # Stores it in form of an dictionary and gives lambda(g)_approx for key g
    for norm in possible_norm_set:
        gens = bdd_ideals[norm]
        for g in gens:
            lambda_g_approx = vector_delta_approximation(log_map(g), delta_1)
            lambda_gens_approx[g] = lambda_g_approx

    # Step 3
    # Find a list of all generators corresponding to each ideal a_l
    generator_lists = []
    for l in range(class_number):
        this_ideal = class_group_reps[l]
        this_ideal_norm = class_group_rep_norms[l]
        gens = []
        for i in range(1, (B + 1).ceil()):
            for g in bdd_ideals[i * this_ideal_norm]:
                if g in this_ideal:
                    gens.append(g)
        generator_lists.append(gens)

    # Step 4
    # Finds all relevant pair and their height
    gen_height_approx_dictionary = {}
    relevant_pair_lists = []

    for n in range(class_number):
        relevant_pairs = []
        gens = generator_lists[n]
        l = len(gens)
        for i in range(l):
            for j in range(i+1, l):
                if K.ideal(gens[i], gens[j]) == class_group_reps[n]:
                    relevant_pairs.append([i, j])
                    gen_height_approx_dictionary[(n, i, j)] = log_height_for_generators_approx(gens[i], gens[j], t/6)
        relevant_pair_lists.append(relevant_pairs)

    # Step 5
    b = rational_in(t/12 + RR(B).log(), t/4 + RR(B).log())
    maximum = 0
    for n in range(class_number):
        for p in relevant_pair_lists[n]:
            maximum = max(maximum, gen_height_approx_dictionary[(n, p[0], p[1])])
    d_tilde = b + t/6 + maximum

    # Step 6
    # computes fundamental units and their value under log map
    fund_units = UnitGroup(K).fundamental_units()
    fund_unit_logs = [log_map(fund_units[i]) for i in range(r)]
    S = column_matrix(fund_unit_logs).delete_rows([r])
    S_inverse = S.inverse()
    S_norm = S.norm(Infinity)
    S_inverse_norm = S_inverse.norm(Infinity)

    upper_bound = (r**2) * max(S_norm, S_inverse_norm)
    m = RR(upper_bound).ceil() + 1

    # Step 7
    # Variables needed for rational approximation
    lambda_tilde = (t/12) / (d_tilde*r*(1+m))
    delta_tilde = min(lambda_tilde/((r**2)*((m**2)+m*lambda_tilde)), 1/(r**2))
    M = d_tilde * (upper_bound+lambda_tilde*RR(r).sqrt())
    M = RR(M).ceil()
    d_tilde = RR(d_tilde)
    delta_2 = min(delta_tilde, (t/6)/(r*(r+1)*M))

    # Step 8, 9
    # Computes relevant points in polytope
    fund_unit_log_approx = [vector_delta_approximation(fund_unit_logs[i], delta_2) for i in range(r)]
    S_tilde = column_matrix(fund_unit_log_approx).delete_rows([r])
    S_tilde_inverse = S_tilde.inverse()
    U = integer_points_in_polytope(S_tilde_inverse, d_tilde)

    # Step 10
    # tilde suffixed list are used for computing second list (L_primed)
    yield K(0)
    U0 = []
    U0_tilde = []
    L0 = []
    L0_tilde = []

    # Step 11
    # Computes unit height
    unit_height_dict = {}
    U_copy = copy(U)
    inter_bound = b - (5*t)/12

    for u in U:
        u_log = sum([u[j]*vector(fund_unit_log_approx[j]) for j in range(r)])
        unit_log_dict[u] = u_log
        u_height = sum([max(u_log[k], 0) for k in range(r + 1)])
        unit_height_dict[u] = u_height
        if u_height < inter_bound:
                U0.append(u)
        if inter_bound <= u_height and u_height < b - (t/12):
            U0_tilde.append(u)
        if u_height > t/12 + d_tilde:
            U_copy.remove(u)
    U = U_copy

    relevant_tuples = set(U0 + U0_tilde)

    # Step 12
    # check for relevant packets
    for n in range(class_number):
        for pair in relevant_pair_lists[n]:
            i = pair[0]
            j = pair[1]
            u_height_bound = b + gen_height_approx_dictionary[(n, i, j)] + t/4
            for u in U:
                if unit_height_dict[u] < u_height_bound:
                    candidate_height = packet_height(n, pair, u)
                    if candidate_height <= b - 7*t/12:
                        L0.append([n, pair, u])
                        relevant_tuples.add(u)
                    elif candidate_height < b + t/4:
                        L0_tilde.append([n, pair, u])
                        relevant_tuples.add(u)

    # Step 13
    # forms a dictionary of all_unit_tuples and their value
    tuple_to_unit_dict = {}
    for u in relevant_tuples:
        unit = K.one()
        for k in range(r):
            unit *= fund_units[k]**u[k]
        tuple_to_unit_dict[u] = unit

    # Step 14
    # Build all output numbers
    roots_of_unity = K.roots_of_unity()
    for u in U0 + U0_tilde:
        for zeta in roots_of_unity:
            yield zeta * tuple_to_unit_dict[u]

    # Step 15
    for p in L0 + L0_tilde:
        gens = generator_lists[p[0]]
        i = p[1][0]
        j = p[1][1]
        u = p[2]
        c_p = tuple_to_unit_dict[u] * (gens[i] / gens[j])
        for zeta in roots_of_unity:
            yield zeta * c_p
            yield zeta / c_p
