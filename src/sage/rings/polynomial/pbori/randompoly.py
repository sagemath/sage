from random import Random
from pprint import pformat

from .PyPolyBoRi import (Monomial, Polynomial, Variable)
from .pbori import random_set, set_random_seed, ll_red_nf_redsb
from .ll import ll_encode
from .blocks import declare_ring


def gen_random_poly(ring, l, deg, vars_set, seed=123):
    """
    Generate a random polynomial with coefficients in ``ring``.

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import Ring, Variable
        sage: from sage.rings.polynomial.pbori.randompoly import gen_random_poly
        sage: r = Ring(16)
        sage: vars = [Variable(i,r) for i in range(10)]
        sage: gen_random_poly(r, 4, 10, vars)  # random
        x(0)*x(1)*x(2)*x(5)*x(8)*x(9) + x(0)*x(1)*x(4)*x(6) + x(0)*x(2)*x(3)*x(7)*x(9) + x(5)*x(8)
    """
    myrange = vars_set
    r = Random(seed)

    def helper(samples):
        if samples == 0:
            return Polynomial(ring.zero())
        if samples == 1:
            d = r.randint(0, deg)
            variables = r.sample(myrange, d)
            m = Monomial(ring)
            for v in sorted(set(variables), reverse=True):
                m = m * Variable(v, ring)
            return Polynomial(m)
        assert samples >= 2
        return helper(samples // 2) + helper(samples - samples // 2)
    p = Polynomial(ring.zero())
    while len(p) < l:
        p = Polynomial(p.set().union(helper(l - len(p)).set()))
    return p


def sparse_random_system(ring, number_of_polynomials, variables_per_polynomial,
                         degree, random_seed=None):
    r"""
    Generate a sparse random system.

    Generate a system, which is sparse in the sense, that each polynomial
    contains only a small subset of variables. In each variable that occurrs
    in a polynomial it is dense in the terms up to the given degree
    (every term occurs with probability 1/2).

    The system will be satisfiable by at least one solution.

    TESTS::

        sage: from sage.rings.polynomial.pbori import Ring, groebner_basis
        sage: r = Ring(10)
        sage: from sage.rings.polynomial.pbori.randompoly import sparse_random_system
        sage: s = sparse_random_system(r, number_of_polynomials=20, variables_per_polynomial=3, degree=2, random_seed=int(123))
        sage: [p.deg() for p in s]
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: sorted(groebner_basis(s), reverse=True)
        [x(0), x(1) + 1, x(2), x(3) + 1, x(4) + 1, x(5), x(6), x(7) + 1, x(8) + 1, x(9) + 1]
    """
    if random_seed is not None:
        set_random_seed(random_seed)
    random_generator = Random(random_seed)
    variables = [ring.variable(i) for i in range(ring.n_variables())]
    solutions = [v + random_generator.randint(0, 1) for v in variables]
    solutions = ll_encode(solutions)
    res = []
    while len(res) < number_of_polynomials:
        variables_as_monomial = Monomial(
            random_generator.sample(
                variables,
                variables_per_polynomial)
        )
        p = Polynomial(random_set(variables_as_monomial, 2 ** (
            variables_per_polynomial - 1)))
        p = sum([p.graded_part(i) for i in range(degree + 1)])
        if p.deg() == degree:
            res.append(p)
    # evaluate it to guarantee a solution
    return [p + ll_red_nf_redsb(p, solutions) for p in res]


def sparse_random_system_data_file_content(number_of_variables, **kwds):
    r"""
    TESTS::

        sage: from sage.rings.polynomial.pbori.randompoly import sparse_random_system_data_file_content
        sage: sparse_random_system_data_file_content(10, number_of_polynomials=5, variables_per_polynomial=3, degree=2, random_seed=int(123))
        "declare_ring(['x'+str(i) for in range(10)])\nideal=\\\n[...]\n\n"
    """
    dummy_dict = {}
    r = declare_ring(['x' + str(i) for i in range(number_of_variables)],
                     dummy_dict)
    polynomials = sparse_random_system(r, **kwds)
    polynomials = pformat(polynomials)
    return "declare_ring(['x'+str(i) for in range({})])\nideal=\\\n{}\n\n".format(
        number_of_variables, polynomials)
