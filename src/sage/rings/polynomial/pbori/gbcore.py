from copy import copy
from itertools import chain
from inspect import getfullargspec as getargspec

from .nf import GeneratorLimitExceeded, symmGB_F2_C, symmGB_F2_python
from .pbori import GroebnerStrategy, ll_red_nf_redsb
from .PyPolyBoRi import (Monomial, Polynomial,
                         OrderCode)
from .ll import eliminate, ll_encode
from .statistics import used_vars_set
from .heuristics import dense_system, gauss_on_linear
from .easy_polynomials import easy_linear_polynomials
from .interpolate import lex_groebner_basis_for_polynomial_via_variety
from .fglm import _fglm


def get_options_from_function(f):
    (argnames, varargs, varopts, defaults) = getargspec(f)[:4]
    return dict(zip(argnames[-len(defaults):], defaults))


def filter_oldstyle_options(**options):
    filtered = dict()
    for key in options:
        newkey = key
        for prefix in ['use_', 'opt_allow_', 'opt_']:
            newkey = newkey.replace(prefix, '')
        filtered[newkey] = options[key]
    return filtered


def filter_newstyle_options(func, **options):
    allowed = get_options_from_function(func).keys()
    filtered = dict()
    for key in options.keys():
        for prefix in ['', 'use_', 'opt_', 'opt_allow_']:
            if prefix + key in allowed:
                filtered[prefix + key] = options[key]

    return filtered


def want_interpolation_gb(G):
    if not G or G[0].ring().get_order_code() != OrderCode.lp or len(G) != 1:
        return False
    p = Polynomial(G[0])
    return not (p.lead_deg() <= 1 or p.set().n_nodes() > 1000)


def ll_is_good(I):
    lex_lead = set()
    for p in I:
        if not p.is_zero():
            m = p.lex_lead()
            if m.deg() == 1:
                lex_lead.add(next(iter(m.variables())).index())
    if len(lex_lead) >= 0.8 * len(I):
        uv = used_vars_set(I).deg()  # don't use len here, which will yield 1
        if len(lex_lead) > 0.9 * uv:
            if uv - len(lex_lead) > 16:
                return "llfirstonthefly"
            else:
                return "llfirst"
    return False


def ll_heuristic(d):
    d = copy(d)
    I = d["I"]
    if ("llfirstonthefly" not in d) and ("llfirst" not in d):
        hint = ll_is_good(I)
        if hint:
            d[hint] = True
    return d


def change_order_heuristic(d):
    d_orig = d
    d = copy(d)
    I = d["I"]
    if not I:
        return d
    switch_table = {OrderCode.lp: OrderCode.dp_asc, OrderCode.dlex: OrderCode.
        dp_asc}
    if "other_ordering_first" not in d:
        # TODO after ll situation might look much different, so heuristic is on
        # wrong place
        code = next(iter(I)).ring().get_order_code()
        if code in switch_table:
            max_non_linear = len(I) // 2
            non_linear = 0
            if code == OrderCode.lp:
                for p in I:
                    if p.lead_deg() > 1:
                        non_linear = non_linear + 1
                        if non_linear > max_non_linear:
                            break
            if (non_linear > max_non_linear) or (code != OrderCode.lp):
                other_ordering_opts = copy(d_orig)
                other_ordering_opts["switch_to"] = switch_table[code]
                d["other_ordering_first"] = other_ordering_opts
    return d


def interpolation_gb_heuristic(d):
    d = copy(d)
    I = d["I"]
    if not d.get("other_ordering_opts", False) and want_interpolation_gb(I):
        d["interpolation_gb"] = True
        d["other_ordering_first"] = False
    return d


def linear_algebra_heuristic(d):
    d = copy(d)
    I = d["I"]

    def want_la():
        if not I:
            return False
        n_used_vars = None
        bound = None
        if next(iter(I)).ring().has_degree_order():
            new_bound = 200
            n_used_vars = used_vars_set(I, bound=new_bound).deg()
            if n_used_vars < new_bound:
                return True
            bound = new_bound
        if dense_system(I):
            new_bound = 100
            if not (bound and new_bound < bound):
                n_used_vars = used_vars_set(I, bound=new_bound).deg()
                bound = new_bound
            if n_used_vars < bound:
                return True
        return False
    if not (("faugere" in d and not d["faugere"]) or
            ("noro" in d and d["noro"])):
        if ("faugere" in d and d["faugere"]) or want_la():

            d["faugere"] = True
            if "red_tail" not in d:
                d["red_tail"] = False
            if "selection_size" not in d:
                d["selection_size"] = 10000
            if "ll" not in d:
                d["ll"] = True

    return d


def trivial_heuristic(d):
    return d


class HeuristicalFunction(object):
    def __call__(self, *args, **kwds):
        complete_dict = copy(kwds)
        heuristic = True
        try:
            heuristic = complete_dict["heuristic"]
        except KeyError:
            pass
        for (k, v) in zip(self.argnames, args):
            complete_dict[k] = v
        if heuristic:
            complete_dict = self.heuristicFunction(complete_dict)
        return self.f(**complete_dict)

    def __init__(self, f, heuristic_function):
        self.argnames, self.varargs, self.varopts, self.defaults = getargspec(f)[:4]
        if hasattr(f, "options"):
            self.options = f.options
        else:
            self.options = dict(zip(self.argnames[-len(self.defaults):], self.
                defaults))
        self.heuristicFunction = heuristic_function
        self.f = f
        self.__doc__ = f.__doc__


def with_heuristic(heuristic_function):
    def make_wrapper(f):
        wrapped = HeuristicalFunction(f, heuristic_function)
        wrapped.__name__ = f.__name__
        return wrapped
    return make_wrapper


def clean_polys_pre(I):
    wrap = (Polynomial(p) for p in I)
    return (list(set(p for p in wrap if not p.is_zero())), None)


def gb_with_pre_post_option(option, pre=None,
                            post=None, if_not_option=tuple(),
                            default=False):
    def make_wrapper(f):
        def wrapper(I, **kwds):
            prot = kwds.get("prot", False)
            for o in if_not_option:
                if (o in kwds and kwds[o]) or (o not in kwds and
                        groebner_basis.options[o]):
                    option_set = False
            if "option_set" not in locals():
                if option in kwds:

                    option_set = kwds[option]
                else:
                    option_set = default
            kwds = dict(((o, kwds[o]) for o in kwds if o != option))
            state = None

            if option_set:
                if pre:
                    pre_args = getargspec(pre)[0]
                    if prot:
                        print("preprocessing for option:", option)

                    local_symbols = copy(locals())
                    (I, state) = pre(**dict([(k, v) for (k, v) in
                        local_symbols.items() if k in pre_args]))
            I = f(I, **kwds)
            if option_set:
                if post:
                    post_args = getargspec(post)[0]
                    if prot:
                        print("postprocessing for option:", option)
                    local_symbols = copy(locals())
                    I = post(**{k: v for (k, v) in local_symbols.items()
                                if k in post_args})

            return I
        wrapper.__name__ = f.__name__
        wrapper.__doc__ = f.__doc__
        if hasattr(f, "options"):
            wrapper.options = copy(f.options)
        else:

            wrapper.options = get_options_from_function(f)

        wrapper.options[option] = default
        return wrapper
    return make_wrapper


def redsb_post(I, state):
    if not I:
        return []
    else:
        return I.minimalize_and_tail_reduce()


def minsb_post(I, state):
    if not I:
        return []
    else:
        return I.minimalize()


def invert_all(I):
    return [p.map_every_x_to_x_plus_one() for p in I]


def invert_all_pre(I):
    return (invert_all(I), None)


def invert_all_post(I, state):
    return invert_all(I)


def llfirst_pre(I, prot):
    (eliminated, llnf, I) = eliminate(I, on_the_fly=False, prot=prot)
    return (I, eliminated)


def ll_constants_pre(I):
    ll_res = []

    while len([p for p in I if p.lex_lead_deg() == 1 and
    (p + p.lex_lead()).constant()]) > 0:
        I_new = []
        ll = []
        leads = set()
        for p in I:
            if p.lex_lead_deg() == 1:
                l = p.lead()
                if not (l in leads) and p.is_singleton_or_pair():
                    tail = p + l
                    if tail.deg() <= 0:
                        ll.append(p)
                        leads.add(l)
                        continue
            I_new.append(p)
        encoded = ll_encode(ll)
        reduced = []
        for p in I_new:
            p = ll_red_nf_redsb(p, encoded)
            if not p.is_zero():
                reduced.append(p)
        I = reduced
        ll_res.extend(ll)
    return (I, ll_res)


def variety_size_from_gb(I):
    """
    TESTS::

        sage: from sage.rings.polynomial.pbori import Ring, Monomial, Polynomial
        sage: from sage.rings.polynomial.pbori.gbcore import variety_size_from_gb
        sage: r = Ring(100)
        sage: x = r.variable
        sage: variety_size_from_gb([])
        1
        sage: variety_size_from_gb([Polynomial(0, r)])
        1
        sage: variety_size_from_gb([Polynomial(1, r)])
        0.0
        sage: variety_size_from_gb([x(1)])
        1.0
        sage: variety_size_from_gb([x(1), x(2)])
        1.0
        sage: variety_size_from_gb([x(1), x(2)*x(3)])
        3.0
        sage: variety_size_from_gb([x(1), x(1)*x(4), x(2)*x(3)])
        6.0
        sage: variety_size_from_gb([x(1)*x(2), x(2)*x(3)])
        5.0
        sage: mons = [Monomial([r.variable(i) for i in range(100) if i!=j])\
            for j in range(100)]
        sage: variety_size_from_gb(mons)
        1.2676506002282294e+30
    """
    I = (Polynomial(p) for p in I)
    I = [p for p in I if not p.is_zero()]
    if not I:
        return 1
#     # TODO Here's something wrong! See the example with 5 solutions.
#     # (reverting for now)
#     number_of_used_vars = used_vars_set(I).deg()
#     leads = set([p.lead() for p in I])
#     minimal_leads = BooleSet(leads).minimal_elements()
#     number_of_used_vars_minimal_leads =\
#         minimal_leads.vars().deg()
#     standard_monomials =\
#         minimal_leads.include_divisors().diff(minimal_leads)
#     return standard_monomials.size_double()*\
#         2**(number_of_used_vars-number_of_used_vars_minimal_leads)

    sm = Monomial(used_vars_set(I)).divisors()
    for p in I:
        m = p.lead()
        sm = sm.diff(sm.multiples_of(m))
    return sm.size_double()


def other_ordering_pre(I, option_set, kwds):
    """
    TESTS::

        sage: from sage.rings.polynomial.pbori.blocks import declare_ring
        sage: r = declare_ring(['x0', 'x1', 'x2', 'x3', 'x4'], globals())
        sage: id = [x1*x3 + x1 + x2*x3 + x3 + x4, x0*x3 + x0 + x1*x2 + x2 + 1,  x1*x3 + x1*x4 + x3*x4 + x4 + 1, x0*x2 + x0*x4 + x1 + x3 + x4]
        sage: from sage.rings.polynomial.pbori.gbcore import groebner_basis
        sage: groebner_basis(id)
        [1]

    """
    if not I:
        return (I, None)

    main_kwds = kwds
    options = option_set

    old_ring = next(iter(I)).ring()
    try:
        new_ring = old_ring.clone(ordering=options["switch_to"])

        kwds = {k: options[k] for k in options
                if k not in ("other_ordering_first", "switch_to", "I")}
        kwds["redsb"] = True
        I = groebner_basis([new_ring(poly) for poly in I], **kwds)
        variety_size = variety_size_from_gb(I)

        fglm_bound = options.get("fglm_bound") or groebner_basis.options["fglm_bound"]
        if variety_size < fglm_bound:
            main_kwds["convert_with_fglm_from_ring"] = new_ring
            main_kwds["convert_with_fglm_to_ring"] = old_ring
        else:
            I = [old_ring(poly) for poly in I]
    finally:
        pass

    return (I, None)


def llfirstonthefly_pre(I, prot):
    (eliminated, llnf, I) = eliminate(I, on_the_fly=True)
    return (I, eliminated)


def gauss_on_linear_pre(I, prot):
    return (gauss_on_linear(I), None)


def easy_linear_polynomials_pre(I):
    res = []
    for p in I:
        res.append(p)
        res.extend(easy_linear_polynomials(p))

    return (list(set(res)), None)


def llfirst_post(I, state, prot, kwds):
    eliminated = state
    for p in I:
        if p.is_one():
            return [p]

    if eliminated:
        I = list(chain(I, eliminated))
        # redsb just for safety, as don't know how option is set
        kwds = copy(kwds)
        kwds.update(
            dict(llfirst=False,
            llfirstonthefly=False,
            ll_constants=False,
            deg_bound=False,
            other_ordering_first=False,
            eliminate_identical_variables=False, redsb=True))
        I = groebner_basis(I, **kwds)
    return I


def ll_constants_post(I, state):
    eliminated = state
    for p in I:
        if p.is_one():
            return [p]
    if eliminated:
        I = list(chain(I, eliminated))
        # redsb just for safety, as don't know how option is set
    return I


def result_to_list_post(I, state):
    return list(I)


def fix_deg_bound_post(I, state):
    if isinstance(I, GroebnerStrategy):
        return I.all_generators()
    else:
        return I


def incremental_pre(I, prot, kwds):
    def sort_key(p):
        p = Polynomial(p)
        return (p.navigation().value(), -p.deg())
    I = sorted(I, key=sort_key)
    inc_sys = []
    kwds = copy(kwds)
    kwds['incremental'] = False

    for p in I[:-1]:
        inc_sys.append(p)
        inc_sys = groebner_basis(inc_sys, **kwds)
        if prot:
            print("incrementally calculating GB, adding generator:", p)
    inc_sys.append(I[:-1])
    return (inc_sys, None)


def eliminate_identical_variables_pre(I, prot):
    changed = True
    ll_system = []
    treated_linears = set()
    while changed:
        changed = False
        rules = dict()
        for p in I:
            t = p + p.lead()
            if p.lead_deg() == 1:
                l = p.lead()
                if l in treated_linears:
                    continue
                else:
                    treated_linears.add(l)
                if t.deg() > 0:
                    rules.setdefault(t, [])
                    leads = rules[t]
                    leads.append(l)

        def my_sort_key(l):
            return l.navigation().value()
        for (t, leads) in rules.items():
            if len(leads) > 1:
                changed = True
                leads = sorted(leads, key=my_sort_key, reverse=True)
                chosen = leads[0]
                for v in leads[1:]:
                    ll_system.append(chosen + v)
    if len(ll_system) > 0:
        ll_encoded = ll_encode(ll_system, reduce=True)
        I = set([ll_red_nf_redsb(p, ll_encoded) for p in I])
    return (I, ll_system)


@gb_with_pre_post_option("clean_arguments", pre=clean_polys_pre, default=True)
@gb_with_pre_post_option("easy_linear_polynomials",
    pre=easy_linear_polynomials_pre, default=True)
@gb_with_pre_post_option("result_to_list", post=result_to_list_post,
    default=True)
@with_heuristic(interpolation_gb_heuristic)
@gb_with_pre_post_option("invert", pre=invert_all_pre,
    post=invert_all_post, default=False)
@gb_with_pre_post_option("gauss_on_linear", pre=gauss_on_linear_pre,
    default=True)
@gb_with_pre_post_option("ll_constants", pre=ll_constants_pre,
    post=ll_constants_post, default=True)
@gb_with_pre_post_option("eliminate_identical_variables",
    pre=eliminate_identical_variables_pre, post=llfirst_post, default=True)
@with_heuristic(ll_heuristic)
@gb_with_pre_post_option("llfirst", if_not_option=["llfirstonthefly"],
    pre=llfirst_pre, post=llfirst_post, default=False)
@gb_with_pre_post_option("llfirstonthefly", pre=llfirstonthefly_pre,
    post=llfirst_post, default=False)
@gb_with_pre_post_option("incremental", pre=incremental_pre)
@with_heuristic(change_order_heuristic)
@gb_with_pre_post_option("other_ordering_first", if_not_option=[
    "interpolation_gb"], pre=other_ordering_pre, default=False)
@with_heuristic(linear_algebra_heuristic)
@gb_with_pre_post_option("fix_deg_bound", if_not_option=["interpolation_gb"],
    post=fix_deg_bound_post, default=True)
@gb_with_pre_post_option("minsb", post=minsb_post, if_not_option=["redsb",
    "deg_bound", "interpolation_gb", "convert_with_fglm_from_ring"],
    default=True)
@gb_with_pre_post_option("redsb", post=redsb_post, if_not_option=["deg_bound",
    "interpolation_gb", "convert_with_fglm_from_ring"], default=True)
def groebner_basis(I, heuristic=True, unique_ideal_generator=False,
        interpolation_gb=False, clean_and_restart_algorithm=False,
        convert_with_fglm_from_ring=None, convert_with_fglm_to_ring=None,
        fglm_bound=40000,
        modified_linear_algebra=True, preprocessor=None, deg_bound=False,
        implementation="Python", full_prot=False, prot=False,
        draw_matrices=False, preprocess_only=False, **impl_options):
    """Computes a Groebner basis of a given ideal I, w.r.t options."""

    if not I:
        return I

    if full_prot:
        prot = True
    if prot:
        print("number of passed generators:", len(I))
    if convert_with_fglm_from_ring is not None:
        from_ring = convert_with_fglm_from_ring
        to_ring = convert_with_fglm_to_ring
        return _fglm(I, from_ring, to_ring)

    if interpolation_gb:
        first = next(iter(I))
        if len(I) != 1 or first.ring().get_order_code() != OrderCode.lp:
            raise ValueError
        return lex_groebner_basis_for_polynomial_via_variety(first)
    if deg_bound is False:
        deg_bound = 100000000
    I = [Polynomial(p) for p in I if not p.is_zero()]
    if unique_ideal_generator and I:
        prod = 1
        for p in I:
            prod = (p + 1) * prod
        I = [prod + 1]

    if implementation == "Python":
        implementation = symmGB_F2_python
    else:
        implementation = symmGB_F2_C

    # custom preprocessing
    if preprocessor:
        I = preprocessor(I)

    if preprocess_only:
        for p in I:
            print(p)
        import sys
        sys.exit(0)

    def call_algorithm(I, max_generators=None):
        return implementation(I,
            deg_bound=deg_bound,
            full_prot=full_prot,
            prot=prot,
            max_generators=max_generators, draw_matrices=draw_matrices,
            **filter_newstyle_options(implementation, **impl_options))

    if clean_and_restart_algorithm:
        for max_generators in [1000, 10000, 50000, 100000, 200000, 300000,
                400000, None]:
            try:
                return call_algorithm(I, max_generators=max_generators)
            except GeneratorLimitExceeded as e:
                I = list(e.strat.all_generators())
                del e.strat
                if prot:
                    print("generator limit exceeded:", max_generators,
                        "restarting algorithm")
    else:
        return call_algorithm(I)


def build_groebner_basis_doc_string():
    additional_options_from_buchberger = filter_oldstyle_options(**
        get_options_from_function(symmGB_F2_python))
    for k in list(additional_options_from_buchberger):
        if k in groebner_basis.options:
            del additional_options_from_buchberger[k]

    groebner_basis.__doc__ = (groebner_basis.__doc__ + "\nOptions are:\n" +
        "\n".join((k + "  :  " + repr(groebner_basis.options[k]) for k in
        groebner_basis.options)) + """

Turn off heuristic by setting heuristic=False
  Additional options come from the actual buchberger implementation.
  In case of our standard Python implementation these are the following:

""" + "\n".join((k + "  :  " + repr(additional_options_from_buchberger[k])
                 for k in additional_options_from_buchberger)))


build_groebner_basis_doc_string()


def _test():
    import doctest
    doctest.testmod()


if __name__ == "__main__":
    _test()
