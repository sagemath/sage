# sage_setup: distribution = sagemath-categories
r"""
Points of bounded height in projective spaces

This module defines functions to compute points of bounded height of a given
number field with height less than a specified bound in projective spaces.

Sage functions to list all elements of a given number field with height less
than a specified bound.

AUTHORS:

- Jing Guo (2022): initial version based on David Krumm's code

REFERENCES:

- [Krumm2016]
"""

import itertools

from math import floor

from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.schemes.projective.projective_space import ProjectiveSpace

lazy_import('sage.geometry.polyhedron.constructor', 'Polyhedron')
lazy_import('sage.libs.pari.all', 'pari')
lazy_import('sage.matrix.constructor', ['matrix', 'column_matrix'])
lazy_import('sage.modules.free_module_element', 'vector')
lazy_import('sage.rings.number_field.unit_group', 'UnitGroup')
lazy_import('sage.rings.real_mpfr', 'RealField')


def ZZ_points_of_bounded_height(PS, dim, bound):
    r"""
    Return an iterator of the points in ``self`` of absolute multiplicative
    height of at most ``bound`` in the rational field.

    INPUT:

    - ``PS`` -- a projective space

    - ``dim`` -- positive integer

    - ``bound`` -- positive integer

    OUTPUT: an iterator of points of bounded height

    EXAMPLES::

        sage: from sage.schemes.projective.proj_bdd_height import ZZ_points_of_bounded_height
        sage: PS = ProjectiveSpace(ZZ, 1)
        sage: sorted(list(ZZ_points_of_bounded_height(PS, 1, 1)))
        [(-1 : -1), (-1 : 0), (-1 : 1), (0 : -1)]
        sage: len(list(ZZ_points_of_bounded_height(PS, 1, 5)))
        40
        sage: sorted(list(ZZ_points_of_bounded_height(PS, 1, 2)))
        [(-2 : -1), (-2 : 1), (-1 : -2), (-1 : -1),
         (-1 : 0), (-1 : 1), (-1 : 2), (0 : -1)]
        sage: PS = ProjectiveSpace(ZZ, 2)
        sage: sorted(list(ZZ_points_of_bounded_height(PS, 2, 1)))
        [(-1 : -1 : -1), (-1 : -1 : 0), (-1 : -1 : 1), (-1 : 0 : -1),
         (-1 : 0 : 0), (-1 : 0 : 1), (-1 : 1 : -1), (-1 : 1 : 0),
         (-1 : 1 : 1), (0 : -1 : -1), (0 : -1 : 0), (0 : -1 : 1),
         (0 : 0 : -1)]

    There are no points of negative height::

        sage: from sage.schemes.projective.proj_bdd_height import ZZ_points_of_bounded_height
        sage: PS = ProjectiveSpace(ZZ, 1)
        sage: list(ZZ_points_of_bounded_height(PS, 1, -3))
        []
    """
    if bound < 1:
        return iter(set())

    points_of_bounded_height = set()

    for t in itertools.product(range(-bound, bound+1), repeat=dim+1):
        if gcd(t) == 1:
            point = PS(t)
            if point not in points_of_bounded_height:
                points_of_bounded_height.add(point)
                yield point


def QQ_points_of_bounded_height(PS, dim, bound, normalize=False):
    r"""
    Return an iterator of the points in ``self`` of absolute multiplicative
    height of at most ``bound`` in the rational field.

    INPUT:

    - ``PS`` -- a projective space

    - ``dim`` -- positive integer

    - ``bound`` -- a real number

    - ``normalize`` -- boolean (default: ``False``); whether to
      normalize the coordinates of returned points

    OUTPUT: an iterator of points of bounded height

    EXAMPLES::

        sage: from sage.schemes.projective.proj_bdd_height import QQ_points_of_bounded_height
        sage: PS = ProjectiveSpace(QQ, 1)
        sage: sorted(list(QQ_points_of_bounded_height(PS, 1, 1)))
        [(-1 : 1), (0 : 1), (1 : 0), (1 : 1)]
        sage: len(list(QQ_points_of_bounded_height(PS, 1, 5)))
        40
        sage: sorted(list(QQ_points_of_bounded_height(PS, 1, 2)))
        [(-2 : 1), (-1 : 1), (-1/2 : 1), (0 : 1),
         (1/2 : 1), (1 : 0), (1 : 1), (2 : 1)]
        sage: sorted(list(QQ_points_of_bounded_height(PS, 1, 2, normalize=True)))
        [(-2 : 1), (-1 : 1), (-1 : 2), (0 : 1),
         (1 : 0), (1 : 1), (1 : 2), (2 : 1)]

    There are no points of negative height::

        sage: from sage.schemes.projective.proj_bdd_height import QQ_points_of_bounded_height
        sage: PS = ProjectiveSpace(QQ, 1)
        sage: list(QQ_points_of_bounded_height(PS, 1, -3))
        []
    """
    if bound < 1:
        return iter(set())

    unit_tuples = list(itertools.product([-1, 1], repeat=dim))
    points_of_bounded_height = set()
    increasing_tuples = itertools.combinations_with_replacement(range(floor(bound + 1)), dim + 1)
    for t in increasing_tuples:
        if gcd(t) == 1:
            for p in itertools.permutations(t):
                for u in unit_tuples:
                    point = PS([a*b for a, b in zip(u, p)] + [p[dim]])
                    if point not in points_of_bounded_height:
                        if normalize:
                            point.scale_by(lcm([point[i].denominator() for i in range(dim + 1)]))

                        points_of_bounded_height.add(point)
                        yield point


def IQ_points_of_bounded_height(PS, K, dim, bound):
    r"""
    Return an iterator of the points in ``self`` of absolute multiplicative
    height of at most ``bound`` in the imaginary quadratic field ``K``.

    INPUT:

    - ``PS`` -- a projective space

    - ``K`` -- a number field

    - ``dim`` -- a positive interger

    - ``bound`` -- a real number

    OUTPUT: an iterator of points of bounded height

    EXAMPLES:

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.projective.proj_bdd_height import IQ_points_of_bounded_height
        sage: CF.<a> = CyclotomicField(3)
        sage: P.<x,y,z> = ProjectiveSpace(CF, 2)
        sage: len(list(IQ_points_of_bounded_height(P, CF, 2, -1)))
        0
        sage: len(list(IQ_points_of_bounded_height(P, CF, 2, 1)))
        57
    """
    if bound < 1:
        return iter([])

    unit_tuples = list(itertools.product(K.roots_of_unity(), repeat=dim))

    class_group_ideals = [c.ideal() for c in K.class_group()]
    class_group_ideal_norms = [i.norm() for i in class_group_ideals]
    class_number = len(class_group_ideals)

    possible_norm_set = set()
    for i in range(class_number):
        for k in range(1, floor(bound + 1)):
            possible_norm_set.add(k*class_group_ideal_norms[i])

    coordinate_space = {}
    coordinate_space[0] = [K(0)]
    for m in possible_norm_set:
        coordinate_space[m] = K.elements_of_norm(m)

    for i in range(class_number):
        a = class_group_ideals[i]
        a_norm = class_group_ideal_norms[i]
        a_norm_bound = bound * a_norm
        a_coordinates = []

        for m in coordinate_space:
            if m <= a_norm_bound:
                for x in coordinate_space[m]:
                    if x in a:
                        a_coordinates.append(x)

        points_in_class_a = set()
        t = len(a_coordinates) - 1
        increasing_tuples = itertools.combinations_with_replacement(range(t + 1), dim + 1)
        for index_tuple in increasing_tuples:
            point_coordinates = [a_coordinates[i] for i in index_tuple]
            if a == K.ideal(point_coordinates):
                for p in itertools.permutations(point_coordinates):
                    for u in unit_tuples:
                        point = PS([i*j for i, j in zip(u, p)] + [p[dim]])

                        if point not in points_in_class_a:
                            points_in_class_a.add(point)
                            yield point


def points_of_bounded_height(PS, K, dim, bound, prec=53):
    r"""
    Return an iterator of the points in ``K`` with dimension ``dim`` of
    absolute multiplicative height of at most ``bound``.

    ALGORITHM:

    This is an implementation of Algorithm 6 in [Krumm2016]_.

    INPUT:

    - ``PS`` -- a projective space

    - ``K`` -- a number field

    - ``dim`` -- positive integer

    - ``bound`` -- a real number

    - ``prec`` -- (default: 53) a positive integer

    OUTPUT: an iterator of points of bounded height

    EXAMPLES::

        sage: from sage.schemes.projective.proj_bdd_height import points_of_bounded_height
        sage: x = polygen(ZZ, 'x')

        sage: # needs sage.geometry.polyhedron sage.libs.pari sage.rings.number_field
        sage: K.<a> = NumberField(x^3 - 7)
        sage: P.<x,y,z> = ProjectiveSpace(K, 2)
        sage: sorted(list(points_of_bounded_height(P, K, 2, 1)))
        [(0 : 0 : 1), (0 : 1 : 0), (1 : 0 : 0), (0 : -1 : 1), (0 : 1 : 1),
         (-1 : 0 : 1), (1 : 0 : 1), (1 : 1 : 0), (-1 : 1 : 0), (-1 : -1 : 1),
         (-1 : 1 : 1), (1 : -1 : 1), (1 : 1 : 1)]

    ::

        sage: # needs sage.geometry.polyhedron sage.libs.pari sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(3*x^2 + 1)
        sage: O = K.maximal_order()
        sage: P.<z,w> = ProjectiveSpace(O, 1)
        sage: len(list(P.points_of_bounded_height(bound=2)))
        44

    ::

        sage: # needs sage.geometry.polyhedron sage.libs.pari sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(3*x^2 + 1)
        sage: O = K.maximal_order()
        sage: P.<z,w> = ProjectiveSpace(O, 1)
        sage: sorted(list(P.points_of_bounded_height(bound=1)))
        [(-1 : 1), (-3/2*a - 1/2 : 1), (3/2*a - 1/2 : 1), (0 : 1),
         (-3/2*a + 1/2 : 0), (-3/2*a + 1/2 : 1), (3/2*a + 1/2 : 1), (1 : 1)]

    ::

        sage: # needs sage.geometry.polyhedron sage.libs.pari sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<z> = NumberField(x^2 - 2)
        sage: R2.<y> = K[]
        sage: L.<w> = K.extension(y^2 - 3)
        sage: P.<a,b> = ProjectiveSpace(L, 1)
        sage: len(list(P.points_of_bounded_height(bound=2)))
        256
    """
    if bound < 1:
        return iter([])

    r1, r2 = K.signature()
    r = r1 + r2 - 1

    if K.is_relative():
        K_degree = K.relative_degree()
    else:
        K_degree = K.degree()

    roots_of_unity = K.roots_of_unity()
    unit_tuples = list(itertools.product(roots_of_unity, repeat=dim))

    log_embed = K.logarithmic_embedding()

    Reals = RealField(prec)
    logB = Reals(bound).log()

    class_group_ideals = [c.ideal() for c in K.class_group()]
    class_number = len(class_group_ideals)

    if K.is_relative():
        class_group_ideal_norms = [i.absolute_norm() for i in class_group_ideals]
    else:
        class_group_ideal_norms = [i.norm() for i in class_group_ideals]

    norm_bound = bound * max(class_group_ideal_norms)
    fundamental_units = UnitGroup(K).fundamental_units()
    fund_unit_logs = list(map(log_embed, fundamental_units))
    mat = column_matrix(fund_unit_logs)

    test_matrix = mat
    try:
        test_matrix.change_ring(QQ)
    except ValueError:
        raise ValueError('prec too low.')

    cut_fund_unit_logs = mat.delete_rows([r])
    lll_fund_units = []
    for c in pari(cut_fund_unit_logs).qflll().python():
        new_unit = 1
        for i in range(r):
            new_unit *= fundamental_units[i]**c[i]
        lll_fund_units.append(new_unit)
    fundamental_units = lll_fund_units
    fund_unit_logs = list(map(log_embed, fundamental_units))

    possible_norm_set = set()
    for i in range(class_number):
        for k in range(1, floor(bound + 1)):
            possible_norm_set.add(k*class_group_ideal_norms[i])

    principal_ideal_gens = {}
    negative_norm_units = K.elements_of_norm(-1)
    if len(negative_norm_units) == 0:
        for m in possible_norm_set:
            principal_ideal_gens[m] = K.elements_of_norm(m) + K.elements_of_norm(-m)
    else:
        for m in possible_norm_set:
            principal_ideal_gens[m] = K.elements_of_norm(m)

    pr_ideal_gen_logs = {}
    for key in principal_ideal_gens:
        for y in principal_ideal_gens[key]:
            pr_ideal_gen_logs[y] = log_embed(y)

    fund_parallelotope_vertices = []
    for coefficient_tuple in itertools.product([-1/2, 1/2], repeat=r):
        vertex = sum([coefficient_tuple[i]*fund_unit_logs[i] for i in range(r)])
        fund_parallelotope_vertices.append(vertex)

    D_numbers = []
    for v in range(r + 1):
        D_numbers.append(max([vertex[v] for vertex in fund_parallelotope_vertices]))

    A_numbers = []
    for v in range(r + 1):
        A_numbers.append(min([pr_ideal_gen_logs[y][v] for y in pr_ideal_gen_logs]))

    aux_constant = (1/K_degree) * Reals(norm_bound).log()

    L_numbers = []
    for v in range(r1):
        L_numbers.append(aux_constant + D_numbers[v] - A_numbers[v])
    for v in range(r1, r + 1):
        L_numbers.append(2*aux_constant + D_numbers[v] - A_numbers[v])
    L_numbers = vector(L_numbers).change_ring(QQ)

    T = column_matrix(fund_unit_logs).delete_rows([r]).change_ring(QQ)

    # insert_row only takes integers, see https://github.com/sagemath/sage/issues/11328
    M = ((-1)*matrix.identity(r)).insert_row(r, [Integer(1) for i in range(r)])
    M = M.transpose().insert_row(0, [Integer(0) for i in range(r + 1)]).transpose()
    M = M.change_ring(QQ)
    M.set_column(0, L_numbers)
    vertices = map(vector, Polyhedron(ieqs=list(M)).vertices())

    T_it = T.inverse().transpose()
    unit_polytope = Polyhedron([v*T_it for v in vertices])

    coordinate_space = {}
    coordinate_space[0] = [[K(0), log_embed(0)]]
    int_points = unit_polytope.integral_points()

    units_with_logs = {}
    for n in int_points:
        new_unit = 1
        for j in range(r):
            new_unit *= fundamental_units[j]**n[j]
        new_unit_log = sum([n[j]*fund_unit_logs[j] for j in range(r)])
        units_with_logs[n] = [new_unit, new_unit_log]

    for norm in principal_ideal_gens:
        coordinate_list = []
        for y in principal_ideal_gens[norm]:
            for n in int_points:
                unit, unit_log = units_with_logs[n]
                y_log = pr_ideal_gen_logs[y]
                g_log = unit_log + y_log
                bool1 = all(g_log[i] <= aux_constant + D_numbers[i] for i in range(r1))
                bool2 = all(g_log[j] <= 2 * aux_constant + D_numbers[j] for j in range(r1, r + 1))
                if bool1 and bool2:
                    g = unit * y
                    coordinate_list.append([g, g_log])
        if len(coordinate_list) > 0:
            coordinate_space[norm] = coordinate_list

    for m in range(class_number):
        a = class_group_ideals[m]
        a_norm = class_group_ideal_norms[m]
        log_a_norm = Reals(a_norm).log()
        a_const = (logB + log_a_norm)/K_degree
        a_coordinates = []

        for k in range(floor(bound + 1)):
            norm = k * a_norm
            if norm in coordinate_space:
                for pair in coordinate_space[norm]:
                    g, g_log = pair
                    if g in a:
                        bool1 = all(g_log[i] <= a_const + D_numbers[i] for i in range(r1))
                        bool2 = all(g_log[j] <= 2 * a_const + D_numbers[j] for j in range(r1, r + 1))
                        if bool1 and bool2:
                            a_coordinates.append(pair)

        t = len(a_coordinates) - 1
        points_in_class_a = set()
        increasing_tuples = itertools.combinations_with_replacement(range(t + 1), dim + 1)
        log_arch_height_bound = logB + log_a_norm
        for index_tuple in increasing_tuples:
            point_coordinates = [a_coordinates[i][0] for i in index_tuple]
            point_coordinate_logs = [a_coordinates[i][1] for i in index_tuple]
            log_arch_height = sum([max([x[i] for x in point_coordinate_logs]) for i in range(r + 1)])
            if log_arch_height <= log_arch_height_bound and a == K.ideal(point_coordinates):
                for p in itertools.permutations(point_coordinates):
                    for u in unit_tuples:
                        point = PS([i*j for i, j in zip(u, p)] + [p[dim]])

                        if point not in points_in_class_a:
                            points_in_class_a.add(point)
                            yield point
