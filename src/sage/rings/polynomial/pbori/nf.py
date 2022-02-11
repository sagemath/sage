from sage.rings.polynomial.pbori.pbori import mod_mon_set
from .pbori import (BooleSet, GroebnerStrategy, ReductionStrategy,
                    parallel_reduce, easy_linear_factors)
from .PyPolyBoRi import (Monomial, Polynomial, Variable,
                         BoolePolynomialVector)
from .easy_polynomials import (easy_linear_polynomials as
                               easy_linear_polynomials_func)
from .statistics import used_vars_set
from warnings import warn
import os


class GeneratorLimitExceeded(Exception):
    r"""
    Docstring for GeneratorLimitExceeded
    """

    def __init__(self, strat):
        self.strat = strat


def pkey(p):
    return (p[0], len(p))


mat_counter = 0


def build_and_print_matrices(v, strat):
    r"""
    Old solution using PIL, the currently used implementation is done in C++
    and plots the same matrices, as being calculated
    """
    treated = BooleSet()
    v = list(v)
    rows = 0
    polys_in_mat = []
    if not v:
        return
    while v:
        rows = rows + 1
        p = v[0]
        v = v[1:]
        for m in list(p.terms()):
            m = Monomial(m)
            if m not in BooleSet(treated):
                i = strat.select(m)
                if i >= 0:
                    p2 = strat[i]
                    p2 = p2 * (m // p2.lead())
                    v.append(p2)
        polys_in_mat.append(p)
        treated = treated.union(p.set())
    m2i = dict([(v, k) for (k, v) in enumerate(list(Polynomial(BooleSet(
        treated)).terms()))])
    polys_in_mat.sort(key=Polynomial.lead, reverse=True)
    polys_in_mat = [[m2i[t] for t in p.terms()] for p in polys_in_mat]

    global mat_counter
    mat_counter = mat_counter + 1
    from PIL import Image

    rows = len(polys_in_mat)
    cols = len(m2i)
    im = Image.new("1", (cols, rows), "white")
    for i in range(len(polys_in_mat)):
        p = polys_in_mat[i]
        for j in p:
            assert i < rows
            assert j < cols

            im.putpixel((j, i), 0)

    file_name = strat.matrix_prefix + str(mat_counter) + ".png"
    if os.path.exists(file_name):
        os.remove(file_name)
    im.save(file_name)
    del im

    print("MATRIX_SIZE:", rows, "x", cols)


def multiply_polynomials(l, ring):
    r"""

    TESTS::

        sage: from sage.rings.polynomial.pbori import *
        sage: r=Ring(1000)
        sage: x=r.variable
        sage: from sage.rings.polynomial.pbori.nf import multiply_polynomials
        sage: multiply_polynomials([x(3), x(2)+x(5)*x(6), x(0), x(0)+1], r)
        0
    """
    l = [Polynomial(p) for p in l]

    def sort_key(p):
        return p.navigation().value()
    l = sorted(l, key=sort_key)
    res = Polynomial(ring.one())
    for p in l:
        res = p * res
    return res


def build_and_print_matrices_deg_colored(v, strat):
    r"""
    old PIL solution using a different color for each degree
    """
    if not v:
        return

    treated = BooleSet()
    v = list(v)
    rows = 0
    polys_in_mat = []

    while v:
        rows = rows + 1
        p = v[0]
        v = v[1:]
        for m in list(p.terms()):
            m = Monomial(m)
            if m not in BooleSet(treated):
                i = strat.select(m)
                if i >= 0:
                    p2 = strat[i]
                    p2 = p2 * (m // p2.lead())
                    v.append(p2)
        polys_in_mat.append(p)
        treated = treated.union(p.set())
    m2i = dict([(v, k) for (k, v) in enumerate(BooleSet(treated))])
    max_deg = max([m.deg() for m in BooleSet(treated)])
    if max_deg == 0:
        max_deg = 1
    i2deg = dict([(m2i[m], m.deg()) for m in BooleSet(treated)])
    polys_in_mat = [[m2i[t] for t in p.terms()] for p in polys_in_mat]
    polys_in_mat.sort(key=pkey)
    global mat_counter
    mat_counter = mat_counter + 1
    from PIL import Image
    from PIL import ImageColor

    rows = len(polys_in_mat)
    cols = len(m2i)
    im = Image.new("RGB", (cols, rows), "white")
    for i in range(len(polys_in_mat)):
        p = polys_in_mat[i]
        for j in p:
            assert i < rows
            assert j < cols
            im.putpixel((j, i), ImageColor.getrgb("hsl(" + str(270 - (270 *
                i2deg[j]) / max_deg) + ",100%,50%)"))
    file_name = strat.matrix_prefix + str(mat_counter) + ".png"
    if os.path.exists(file_name):
        os.remove(file_name)
    im.save(file_name)
    del im

    print("MATRIX_SIZE:", rows, "x", cols)


def high_probability_polynomials_trick(p, strat):
    lead_deg = p.lead_deg()
    if lead_deg <= 4:
        return

    ring = p.ring()
    factor = multiply_polynomials(easy_linear_factors(p), ring)
    p = p / factor

    # again, do it twice, it's cheap
    lead_deg = p.lead_deg()
    if lead_deg <= 3:
        return

    if lead_deg > 9:
        return
    uv = p.vars_as_monomial()

    candidates = []

    if uv.deg() <= 4:
        return

    if not uv.deg() <= lead_deg + 1:
        return

    space = uv.divisors()

    lead = p.lead()
    for v in lead.variables():
        variable_selection = lead // v
        vars_reversed = list(reversed(variable_selection.variables()))
        # it's just a way to loop over the cartesian product
        for assignment in variable_selection.divisors():
            c_p = assignment
            for var in vars_reversed:
                if not assignment.reducible_by(var):
                    c_p = (var + 1) * c_p

            points = (c_p + 1).zeros_in(space)
            if p.zeros_in(points).empty():
                candidates.append(c_p * factor)
        # there many more combinations depending on plugged in values
    for c in candidates:
        strat.add_as_you_wish(c)


def symmGB_F2_python(G, deg_bound=1000000000000, over_deg_bound=0,
        use_faugere=False, use_noro=False, opt_lazy=True, opt_red_tail=True,
        max_growth=2.0, step_factor=1.0, implications=False, prot=False,
        full_prot=False, selection_size=1000, opt_exchange=True,
        opt_allow_recursion=False, ll=False,
        opt_linear_algebra_in_last_block=True, max_generators=None,
        red_tail_deg_growth=True, matrix_prefix='mat',
        modified_linear_algebra=True, draw_matrices=False,
        easy_linear_polynomials=True):
    if use_noro and use_faugere:
        raise ValueError('both use_noro and use_faugere specified')

    def add_to_basis(strat, p):
        if p.is_zero():
            if prot:
                print("-")
        else:
            if prot:
                if full_prot:
                    print(p)
                print("Result: ", "deg:", p.deg(), "lm: ", p.lead(), "el: ", p
                    .elength())
            if easy_linear_polynomials and p.lead_deg() > 2:
                lin = easy_linear_polynomials_func(p)
                for q in lin:
                    strat.add_generator_delayed(q)
            old_len = len(strat)
            strat.add_as_you_wish(p)
            new_len = len(strat)
            if new_len == 1 + old_len:
                high_probability_polynomials_trick(p, strat)

            if prot:
                print("#Generators:", len(strat))

    if isinstance(G, list):
        if not G:
            return []
        G = [Polynomial(g) for g in G]
        current_ring = G[0].ring()
        strat = GroebnerStrategy(current_ring)
        strat.reduction_strategy.opt_red_tail = opt_red_tail
        strat.opt_lazy = opt_lazy
        strat.opt_exchange = opt_exchange
        strat.opt_allow_recursion = opt_allow_recursion
        strat.enabled_log = prot
        strat.reduction_strategy.opt_ll = ll
        strat.opt_modified_linear_algebra = modified_linear_algebra
        strat.opt_linear_algebra_in_last_block = (
            opt_linear_algebra_in_last_block)
        strat.opt_red_by_reduced = False  # True
        strat.reduction_strategy.opt_red_tail_deg_growth = red_tail_deg_growth

        strat.opt_draw_matrices = draw_matrices
        strat.matrix_prefix = matrix_prefix

        for g in G:
            if not g.is_zero():
                strat.add_generator_delayed(g)
    else:
        strat = G

    if prot:
        print("added delayed")
    i = 0
    try:
        while strat.npairs() > 0:
            if max_generators and len(strat) > max_generators:
                raise GeneratorLimitExceeded(strat)
            i = i + 1
            if prot:
                print("Current Degree:", strat.top_sugar())
            if (strat.top_sugar() > deg_bound) and (over_deg_bound <= 0):
                return strat
            if (strat.top_sugar() > deg_bound):
                ps = strat.some_spolys_in_next_degree(over_deg_bound)
                over_deg_bound -= len(ps)
            else:
                ps = strat.some_spolys_in_next_degree(selection_size)

            if ps and ps[0].ring().has_degree_order():
                ps = [strat.reduction_strategy.cheap_reductions(p) for p in ps]
                ps = [p for p in ps if not p.is_zero()]
                if ps:
                    min_deg = min((p.deg() for p in ps))
                new_ps = []
                for p in ps:
                    if p.deg() <= min_deg:
                        new_ps.append(p)
                    else:
                        strat.add_generator_delayed(p)
                ps = new_ps

            if prot:
                print("(", strat.npairs(), ")")
            if prot:
                print("start reducing")
                print("Chain Crit. : ", strat.chain_criterions, "VC:", strat.
                    variable_chain_criterions, "EASYP", strat.
                    easy_product_criterions, "EXTP", strat.
                    extended_product_criterions)
                print(len(ps), "spolys added")

            if use_noro or use_faugere:
                v = BoolePolynomialVector()

                for p in ps:
                    if not p.is_zero():
                        v.append(p)
                if use_noro:
                    res = strat.noro_step(v)
                else:
                    res = strat.faugere_step_dense(v)

            else:
                v = BoolePolynomialVector()
                for p in ps:
                    p = Polynomial(mod_mon_set(BooleSet(p.set()),
                                        strat.reduction_strategy.monomials))
                    if not p.is_zero():
                        v.append(p)
                if len(v) > 100:
                    res = parallel_reduce(v, strat, int(step_factor * 10),
                                          max_growth)
                elif len(v) > 10:
                    res = parallel_reduce(v, strat, int(step_factor * 30),
                                          max_growth)
                else:
                    res = parallel_reduce(v, strat, int(step_factor * 100),
                                          max_growth)

            if prot:
                print("end reducing")

            if res and res[0].ring().has_degree_order():
                res_min_deg = min([p.deg() for p in res])
                new_res = []
                for p in res:
                    if p.deg() == res_min_deg:
                        new_res.append(p)
                    else:
                        strat.add_generator_delayed(p)
                res = new_res

            def sort_key(p):
                return p.lead()
            res_cp = sorted(res, key=sort_key)

            for p in res_cp:
                old_len = len(strat)
                add_to_basis(strat, p)
                if implications and old_len == len(strat) - 1:
                    strat.implications(len(strat) - 1)
                if p.is_one():
                    if prot:
                        print("GB is 1")
                    return strat
                if prot:
                    print("(", strat.npairs(), ")")

            strat.clean_top_by_chain_criterion()
        return strat
    except KeyboardInterrupt:
        raise


def GPS(G, vars_start, vars_end):
    def step(strat, trace, var, val):
        print("plugin: ", var, val)
        print("npairs", strat.npairs())
        strat = GroebnerStrategy(strat)
        print("npairs", strat.npairs())
        strat.add_generator_delayed(Polynomial(
            Monomial(Variable(var, strat.r)) + val))
        strat = symmGB_F2_python(strat, prot=True, deg_bound=2,
            over_deg_bound=10)
        if var <= vars_start:
            strat = symmGB_F2_python(strat, prot=True, opt_lazy=False,
                                     opt_red_tail=False)
        if strat.containsOne():
            pass
        else:
            if var <= vars_start:
                # bug: may contain Delayed polynomials
                print("!!!!!!! SOLUTION", trace)
                raise Exception
                # yield trace
            else:
                branch(strat, trace + [(var, val)], var - 1)

    def branch(strat, trace, var):
        while strat.variableHasValue(var):
            # remember to add value to trace
            var -= 1
        step(strat, trace, var, 0)
        step(strat, trace, var, 1)
    if G:
        strat = GroebnerStrategy(G[0].ring())
        # strat.add_generator(G[0])
        for g in G[:]:
            strat.add_generator_delayed(g)
        branch(strat, [], vars_end - 1)


def GPS_with_proof_path(G, proof_path, deg_bound, over_deg_bound):
    def step(strat, trace, proof_path, pos, val):
        print(proof_path)
        print("plugin: ", pos, val, proof_path[pos])
        print("npairs", strat.npairs())
        strat = GroebnerStrategy(strat)
        print("npairs", strat.npairs())
        print("npairs", strat.npairs())
        plug_p = Polynomial(proof_path[pos]) + val
        plug_p_lead = plug_p.lead()
        if len(plug_p) == 2 and (plug_p + plug_p_lead).deg() == 0:
            for v in plug_p_lead:
                strat.add_generator_delayed(v + 1)
        else:
            strat.add_generator_delayed(plug_p)
        print("npairs", strat.npairs())
        print("pos:", pos)
        strat = symmGB_F2_python(strat, deg_bound=deg_bound, opt_lazy=False,
            over_deg_bound=over_deg_bound, prot=True)
        print("npairs", strat.npairs())
        pos = pos + 1
        if pos >= len(proof_path):
            print("OVER")
            strat = symmGB_F2_python(strat, prot=True)
        if strat.containsOne():
            pass
        else:
            if pos >= len(proof_path):
                print("npairs", strat.npairs())
                print("minimized:")
                for p in strat.minimalize_and_tail_reduce():
                    print(p)
                # bug: may contain Delayed polynomials
                print("!!!!!!! SOLUTION", trace)
                raise Exception
                # yield trace
            else:
                branch(strat, trace + [(pos, val)], proof_path, pos)

    def branch(strat, trace, proof_path, pos):

        step(strat, trace, proof_path, pos, 0)
        step(strat, trace, proof_path, pos, 1)
    strat = GroebnerStrategy(G[0].ring())
    strat.add_generator(Polynomial(G[0]))
    for g in G[1:]:
        strat.add_generator_delayed(Polynomial(g))
    branch(strat, [], proof_path, 0)


def GPS_with_suggestions(G, deg_bound, over_deg_bound, opt_lazy=True,
                         opt_red_tail=True, initial_bb=True):
    def step(strat, trace, var, val):
        print(trace)
        plug_p = val + var
        print("plugin: ", len(trace), plug_p)
        trace = trace + [plug_p]
        strat = GroebnerStrategy(strat)

        strat.add_generator_delayed(plug_p)
        print("npairs", strat.npairs())

        strat = symmGB_F2_python(strat, deg_bound=deg_bound,
            opt_lazy=opt_lazy, over_deg_bound=over_deg_bound, prot=True)

        if not strat.containsOne():
            branch(strat, trace)

    def branch(strat, trace):
        print("branching")
        index = strat.suggestPluginVariable()

        if index < 0:
            uv = set(used_vars_set(strat))
            lv = set([next(iter(p.lead())).index() for p in strat if p.
                lead_deg() == 1])
            candidates = uv.difference(lv)
            if candidates:
                index = next(iter(candidates)).index()
        if index >= 0:
            print("chosen index:", index)
            step(strat, trace, Polynomial(Monomial(Variable(index))), 0)
            step(strat, trace, Polynomial(Monomial(Variable(index))), 1)
        else:
            print("FINAL!!!", index)
            strat = symmGB_F2_python(strat, prot=True)
            if not strat.containsOne():
                print("TRACE", trace)
                print("SOLUTION")
                for p in strat.minimalize_and_tail_reduce():
                    print(p)
                raise Exception

    def sort_crit(p):
        return (p.lead(), p.deg(), p.elength())
    if not G:
        return
    strat = GroebnerStrategy(G[0].ring())
    strat.reduction_strategy.opt_red_tail = opt_red_tail  # True
    strat.opt_exchange = False
    strat.opt_allow_recursion = False
    strat.opt_lazy = opt_lazy
    first_deg_bound = 1
    G = [Polynomial(p) for p in G]
    G.sort(key=sort_crit)
    if initial_bb:
        for g in G:
            print(g)
            if g.deg() == 1:  # (index<0):
                strat.add_as_you_wish(g)
            else:
                first_deg_bound = max(first_deg_bound, g.deg())
                strat.add_generator_delayed(g)
            print(g, len(strat))
    else:
        for g in G:
            strat.add_as_you_wish(g)
    if initial_bb:
        strat = symmGB_F2_python(strat, deg_bound=max(deg_bound,
            first_deg_bound), opt_lazy=opt_lazy, over_deg_bound=0, prot=True)
    strat.opt_lazy = opt_lazy
    print("INITIALIZED")
    branch(strat, [])


def GPS_with_non_binary_proof_path(G, proof_path, deg_bound, over_deg_bound):
    def step(strat, trace, proof_path, pos, choice):
        print(proof_path)
        print("plugin: ", pos)
        print("npairs", strat.npairs())
        strat = GroebnerStrategy(strat)
        print("npairs", strat.npairs())
        print("npairs", strat.npairs())
        for p in proof_path[pos][choice]:
            print(p)
            strat.add_generator_delayed(Polynomial(p))

        print("npairs", strat.npairs())
        print("pos:", pos)
        strat = symmGB_F2_python(strat, deg_bound=deg_bound,
            over_deg_bound=over_deg_bound, prot=True)
        print("npairs", strat.npairs())
        pos = pos + 1
        if pos >= len(proof_path):
            print("OVER")
            strat = symmGB_F2_python(strat)
        if strat.containsOne():
            pass
        else:
            if pos >= len(proof_path):
                print("npairs", strat.npairs())
                # bug: may contain Delayed polynomials
                print("!!!!!!! SOLUTION", trace)
                raise Exception
                # yield trace
            else:
                branch(strat, trace + [(pos, choice)], proof_path, pos)
                # workaround because of stack depth
                # step(strat,trace+[(var,val)],var-1, 0)
                # step(strat,trace+[(var,val)],var-1, 1)

    def branch(strat, trace, proof_path, pos):
        for i in range(len(proof_path[pos])):
            step(strat, trace, proof_path, pos, i)

    strat = GroebnerStrategy(G[0].ring())
    strat.add_generator(G[0])
    for g in G[1:]:
        strat.add_generator_delayed(g)
    branch(strat, [], proof_path, 0)


def symmGB_F2_C(G, opt_exchange=True,
        deg_bound=1000000000000, opt_lazy=False,
        over_deg_bound=0, opt_red_tail=True,
        max_growth=2.0, step_factor=1.0,
        implications=False, prot=False,
        full_prot=False, selection_size=1000,
        opt_allow_recursion=False, use_noro=False, use_faugere=False,
        ll=False, opt_linear_algebra_in_last_block=True,
        max_generators=None, red_tail_deg_growth=True,
        modified_linear_algebra=True, matrix_prefix="",
        draw_matrices=False):
    if use_noro:
        raise NotImplementedError("noro not implemented for symmgb")
    if isinstance(G, list):
        if not G:
            return []

        G = [Polynomial(g) for g in G]
        strat = GroebnerStrategy(G[0].ring())
        strat.reduction_strategy.opt_red_tail = opt_red_tail
        strat.enabled_log = prot
        strat.opt_lazy = opt_lazy
        strat.opt_exchange = opt_exchange
        strat.reduction_strategy.opt_ll = ll
        strat.opt_allow_recursion = opt_allow_recursion
        strat.opt_linear_algebra_in_last_block = (
            opt_linear_algebra_in_last_block)
        strat.enabled_log = prot
        strat.opt_modified_linear_algebra = modified_linear_algebra
        strat.matrix_prefix = matrix_prefix
        strat.opt_draw_matrices = draw_matrices
        strat.reduction_strategy.opt_red_tail_deg_growth = red_tail_deg_growth

        strat.redByReduced = False  # True

        for g in G:  # [1:]:
            if not g.is_zero():
                strat.add_generator_delayed(g)
    strat.symmGB_F2()
    return strat


def normal_form(poly, ideal, reduced=True):
    r"""
    Simple normal form computation of a polynomial  against an ideal.

    TESTS::

        sage: from sage.rings.polynomial.pbori import declare_ring, normal_form
        sage: r=declare_ring(['x','y'], globals())
        sage: normal_form(x+y, [y],reduced=True)
        x
        sage: normal_form(x+y,[x,y])
        0
    """
    ring = poly.ring()
    strat = ReductionStrategy(ring)
    strat.opt_red_tail = reduced
    ideal = [Polynomial(p) for p in ideal if p != 0]
    ideal = sorted(ideal, key=Polynomial.lead)
    last = None
    for p in ideal:
        if p.lead() != last:
            strat.add_generator(p)
        else:
            warn("%s will not used for reductions" % p)
        last = p.lead()
    return strat.nf(poly)


def _test():
    import doctest
    doctest.testmod()


if __name__ == "__main__":
    _test()
