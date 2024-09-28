# sage_setup: distribution = sagemath-brial
# sage.doctest: needs sage.rings.polynomial.pbori
from .pbori import (top_index, if_then_else,
                    mod_mon_set, BooleSet, BooleConstant)
from .PyPolyBoRi import (Polynomial, Monomial, Variable)


def all_monomials_of_degree_d_old(d, variables):
    """
    Return monomials of degree d in the given variables.

    Obsolete version ?
    """
    if d == 0:
        return BooleConstant(1)

    if not variables:
        return []
    variables = sorted(set(variables), reverse=True, key=top_index)

    m = variables[-1]
    for v in variables[:-1]:
        m = v + m
    m = m.set()
    i = 0
    res = Polynomial(variables[0].ring().one()).set()
    while i < d:
        i += 1
        res = res.cartesian_product(m).diff(res)
    return res


def all_monomials_of_degree_d(d, variables):
    """
    Return monomials of degree d in the given variables.
    """
    variables = Monomial(variables)
    variables = list(variables.variables())
    if not variables:
        assert d == 0
        return BooleConstant(1)
    ring = variables[0].ring()
    if d > len(variables):
        return Polynomial(0, ring)
    if d < 0:
        return Polynomial(1, ring)

    deg_variables = variables[-d:]
    # this ensures sorting by indices
    res = Monomial(deg_variables)

    for i in range(1, len(variables) - d + 1):
        deg_variables = variables[-d - i:-i]
        res = Polynomial(res)
        nav = res.navigation()
        navs = []
        while not nav.constant():
            navs.append(BooleSet(nav, ring))
            nav = nav.then_branch()
        acc = Polynomial(1, ring)
        for (nav, v) in reversed(zip(navs, deg_variables)):
            acc = if_then_else(v, acc, nav)
        res = acc
    return res.set()


def power_set(variables):
    """
    Return all subsets of the given variables.
    """
    if not variables:
        return BooleConstant(1)
    variables = sorted(set(variables), reverse=True, key=top_index)
    res = Polynomial(1, variables[0].ring()).set()
    for v in variables:
        res = if_then_else(v, res, res)
    return res


if __name__ == '__main__':
    from .blocks import declare_ring, Block
    r = declare_ring([Block("x", 10000)], globals())
    print(list(all_monomials_of_degree_d(0, [Variable(i) for i in range(100)])))
    print(list(all_monomials_of_degree_d(1, [Variable(i) for i in range(10)])))
    print(list(all_monomials_of_degree_d(2, [Variable(i) for i in range(4)])))
    print(list(all_monomials_of_degree_d(3, [Variable(i) for i in range(4)])))
    print(list(all_monomials_of_degree_d(4, [Variable(i) for i in range(4)])))
    print(list(all_monomials_of_degree_d(0, [])))
    print(list(all_monomials_of_degree_d(1, [])))
    print(list(power_set([Variable(i) for i in range(2)])))
    print(list(power_set([Variable(i) for i in range(4)])))
    print(list(power_set([])))
    # every monomial in the first 8 var, which is at most linear in the first 5
    print(list(mod_mon_set(
        power_set([Variable(i) for i in range(8)]),
        all_monomials_of_degree_d(2, [Variable(i) for i in range(5)]))))

    # specialized normal form computation
    print(Polynomial(
        mod_mon_set(
            (x(1) * x(2) + x(1) + 1).set(),
            all_monomials_of_degree_d(2, [Variable(i) for i in range(1000)]))))
    print(list(mod_mon_set(
        power_set([Variable(i) for i in range(50)]),
        all_monomials_of_degree_d(2, [Variable(i) for i in range(1000)]))))


def monomial_from_indices(ring, indices):
    res = Monomial(ring)
    for i in sorted(indices, reverse=True):
        res = res * ring.variable(i)
    return res
