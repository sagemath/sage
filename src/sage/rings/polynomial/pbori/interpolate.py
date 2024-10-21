# sage_setup: distribution = sagemath-brial
# Copyright (c) 2005-2007 by The PolyBoRi Team
from time import process_time as clock
from random import Random

from .PyPolyBoRi import (Polynomial, Variable, Monomial,
                         BoolePolynomialVector)
from .randompoly import gen_random_poly
from .pbori import (BooleSet, add_up_polynomials, interpolate_smallest_lex,
                    interpolate)
from .blocks import Block, declare_ring


generator = Random()


def add_up_poly_list(l, init):
    v = BoolePolynomialVector()
    for p in l:
        v.append(p)
    return add_up_polynomials(v, init)


def bench_interpolate(degree, nvariables, points):
    c = points
    h = len(points) / 2
    terms = set(c.terms())
    part1 = generator.sample(terms, h)
    part1 = add_up_poly_list(part1, Polynomial(c.ring().zero()))
    part2 = c + part1
    p = part1
    q = part2
    assert part1.set().intersect(part2).empty()
    c1 = clock()
    res2 = interpolate_smallest_lex(p, q)
    c2 = clock()
    print("finished interpolate_smallest_lex(p,q),len:", len(res2),
          "time", c2 - c1)
    c1 = clock()
    res1 = interpolate(p, q)
    c2 = clock()
    print("finished interpolate(p,q)" + len("_smallest_lex") * " " + ",len:",
          res1.set().size_double(), "time:", c2 - c1)
    return res2


def nf_lex_points(f, p):
    f = Polynomial(f)
    p = BooleSet(p)
    z = f.zeros_in(p)
    return interpolate_smallest_lex(z, p.diff(z))


def gen_random_o_z(points, points_p):
    k = generator.randrange(len(points) + 1)
    ones = generator.sample(points, k)
    vec = BoolePolynomialVector()
    for p in ones:
        vec.append(p)
    ones = add_up_polynomials(vec, Polynomial(points_p.ring().zero()))
    return interpolate_smallest_lex(points_p.set().diff(ones), ones)


def variety_lex_leading_terms(points, variables):
    assert isinstance(points, BooleSet), "Points needs to be a BooleSet"
    ring = variables.ring()
    standards = BooleSet(ring.zero())
    points_tuple = tuple(points)
    myvars_div = variables.divisors()
    if points != myvars_div:
        standards = BooleSet(ring.one())
    len_standards = len(standards)
    standards_old = standards
    while len_standards < len(points):
        standards = standards.union(gen_random_o_z(points_tuple, points))

        if standards_old != standards:
            standards = BooleSet(standards).include_divisors()
            len_standards = len(standards)
            standards_old = standards

    return BooleSet(myvars_div.diff(standards)).minimal_elements()


def lex_groebner_basis_points(points, variables):
    leads = variety_lex_leading_terms(points, variables)
    return [nf_lex_points(l, points) + l for l in leads]


def lex_groebner_basis_for_polynomial_via_variety(p):
    variables = p.vars_as_monomial()
    return lex_groebner_basis_points(p.zeros_in(variables.divisors()),
                                     variables)


if __name__ == '__main__':
    nvariables = 100
    r = declare_ring([Block("x", nvariables)])
    for number_of_points in (100, 500, 1000, 2000, 3000,
                             4000, 5000, 10000,
                             20000, 50000, 100000):
        print("----------")
        print("number_of_points:", number_of_points)
        print("generate points")
        points = gen_random_poly(r, number_of_points,
                                 nvariables,
                                 [Variable(i, r) for i in range(nvariables)])
        print("points generated")
        bench_interpolate(nvariables, nvariables, points)
        vars_mon = Monomial(r)
        for i in reversed(range(nvariables)):
            vars_mon = vars_mon * Variable(i, r)
        print(len(variety_lex_leading_terms(points, vars_mon)),
              "elements in groebner basis")
