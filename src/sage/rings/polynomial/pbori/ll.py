from .pbori import (top_index, if_then_else,
                    substitute_variables, BooleSet,
                    ll_red_nf_redsb, ll_red_nf_noredsb,
                    ll_red_nf_noredsb_single_recursive_call)
from .PyPolyBoRi import (Polynomial, Monomial, Ring, BoolePolynomialVector)
from .statistics import used_vars_set
from .rank import rank

lead_index = top_index


def combine(reductors, p, reduce=None):
    p_nav = p.navigation()
    assert p_nav.value() < reductors.navigation().value()
    p_else = BooleSet(p_nav.else_branch(), p.ring())
    if reduce:
        p_else = reduce(p_else, reductors)
    return if_then_else(p_nav.value(), reductors, p_else)


def llredsb_Cudd_style(polys):

    reductors = Polynomial(polys[0].ring().one()).set() if polys else None

    linear_lead = sorted(polys, key=lead_index, reverse=True)
    assert len({p.lex_lead() for p in linear_lead}) == len(polys)
    assert not any(p.constant() for p in polys)
    assert len([p for p in polys if p.lex_lead_deg() == 1]) == len(polys)
    assert len({p.navigation().value() for p in polys}) == len(polys)
    for p in linear_lead:
        reductors = combine(reductors, p, reduce=ll_red_nf_redsb)
    return reductors


def ll_encode(polys, reduce=False, prot=False, reduce_by_linear=True):
    polys = [Polynomial(p) for p in polys]
    linear_lead = sorted(polys, key=lead_index, reverse=True)
    assert len({p.lex_lead() for p in linear_lead}) == len(polys)
    assert not any(p.constant() for p in polys)
    assert len([p for p in polys if p.lex_lead_deg() == 1]) == len(polys)
    assert len({p.navigation().value() for p in polys}) == len(polys)
    if (not reduce) and reduce_by_linear:
        linear_polys = [p for p in polys if p.deg() == 1]
        if linear_polys:
            linear_ll = ll_encode(linear_polys, reduce=True,
                                  reduce_by_linear=False)
            polys = [p.lex_lead() + ll_red_nf_redsb(p + p.lex_lead(),
                                                    linear_ll) for p in polys]
    reduce = ll_red_nf_redsb if reduce else None

    reductors = Polynomial(polys[0].ring().one()).set() if polys else None

    last = None
    counter = 0
    for p in linear_lead:

        if prot:
            counter = counter + 1
            progress = (counter * 100) / len(linear_lead)
            if last != progress:
                print(str(progress) + "%")
            last = progress
        reductors = combine(reductors, p, reduce=reduce)
    return reductors


def eliminate(polys, on_the_fly=False, prot=False, reduction_function=None,
              optimized=True):
    r"""
    There exists an optimized variant, which reorders the variable in a different ring.
    """
    polys = [Polynomial(p) for p in polys]
    rest = []
    linear_leads = []
    linear_leading_monomials = set()
    for p in polys:
        if p.is_zero():
            continue
        lm = p.lex_lead()
        if lm.deg() == 1:

            if lm not in linear_leading_monomials:
                linear_leading_monomials.add(lm)
                linear_leads.append(p)
            else:
                rest.append(p)
        else:
            rest.append(p)
    if not linear_leads:
        def identity(p):
            return p
        return (linear_leads, identity, rest)
    if reduction_function is None:
        if on_the_fly:
            if optimized:
                reduction_function = ll_red_nf_noredsb_single_recursive_call
            else:
                reduction_function = ll_red_nf_noredsb
        else:
            reduction_function = ll_red_nf_redsb

    if optimized:
        llnf, reduced_list = eliminate_ll_ranked(linear_leads, rest,
                                                 reduction_function=reduction_function,
                                                 reduce_ll_system=(not on_the_fly),
                                                 prot=prot)
    else:
        def llnf(p):
            return reduction_function(p, reductors)
        reduced_list = []
        reductors = ll_encode(linear_leads, reduce=(not on_the_fly), prot=prot)
        for p in rest:
            rp = reduction_function(p, reductors)
            if rp.is_one():
                reduced_list = [rp]
                break
            reduced_list.append(rp)

    return (linear_leads, llnf, reduced_list)


def construct_map_by_indices(to_ring, idx_mapping):
    v = BoolePolynomialVector((max(idx_mapping.keys()) + 1) * [to_ring.zero()])
    for (from_idx, to_idx) in idx_mapping.items():
        val = to_ring.variable(to_idx)
        v[from_idx] = val
    return v


def eliminate_ll_ranked(ll_system, to_reduce,
                        reduction_function=ll_red_nf_noredsb,
                        reduce_ll_system=False, prot=False):

    assert ll_system
    from_ring = ll_system[0].ring()

    ll_ranks = rank(ll_system)
    add_vars = set(used_vars_set(to_reduce).variables()).difference(ll_ranks.
                                                                    keys())
    for v in add_vars:
        ll_ranks[v] = -1

        # pushing variables ignored by ll to the front means,
        # that the routines will quickly eliminate them
        # and they won't give any overhead
    def sort_key(v):
        return (ll_ranks[v], v.index())
    sorted_vars = sorted(ll_ranks.keys(), key=sort_key)

    def var_index(v):
        return next(iter(Monomial(v).variables())).index()

    to_ring = Ring(len(sorted_vars))
    map_back_indices = {i: var_index(v) for i, v in enumerate(sorted_vars)}
    map_from_indices = {var_index(v): i for i, v in enumerate(sorted_vars)}

    var_names = [str(v) for v in sorted_vars]
    try:
        for (i, v) in enumerate(sorted_vars):
            assert var_names[i] == str(v), (var_names[i], v, var_index(v), i)

    finally:
        pass
    try:
        map_from_vec = construct_map_by_indices(to_ring, map_from_indices)
    finally:
        pass
    map_back_vec = construct_map_by_indices(from_ring, map_back_indices)

    def map_from(p):
        res = substitute_variables(to_ring, map_from_vec, p)
        return res

    def map_back(p):
        return substitute_variables(from_ring, map_back_vec, p)

    try:
        ll_opt_encoded = ll_encode([map_from(p) for p in ll_system],
                                   prot=False,
                                   reduce=reduce_ll_system)

        def llnf(p):
            return map_back(reduction_function(map_from(p), ll_opt_encoded))
        opt_eliminated = [llnf(p) for p in to_reduce]
    finally:
        pass
    return (llnf, opt_eliminated)


class RingMap:
    r"""
    Define a mapping between two rings by common variable names.

    TESTS::

        sage: from sage.rings.polynomial.pbori.pbori import *
        sage: from sage.rings.polynomial.pbori.blocks import declare_ring, Block
        sage: to_ring = declare_ring([Block("x", 10)], globals())
        sage: from_ring = declare_ring([Block("y", 5), Block("x", 10)], globals())
        sage: from sage.rings.polynomial.pbori.ll import RingMap
        sage: mapping = RingMap(to_ring, from_ring)
        sage: (x(1)+1).navigation().value()
        6
        sage: mapping(x(1)+1)
        x(1) + 1
        sage: mapping(x(1)+1).navigation().value()
        1
        sage: mapping.invert(mapping(x(1)+1))
        x(1) + 1
        sage: mapping.invert(mapping(x(1)+1)).navigation().value()
        6
        sage: mapping(y(1)+1)
        Traceback (most recent call last):
        ...
        RuntimeError: Operands come from different manager.
    """

    def __init__(self, to_ring, from_ring):
        r"""
        Initialize map by two given rings.

        TESTS::

            sage: from sage.rings.polynomial.pbori.pbori import *
            sage: from sage.rings.polynomial.pbori.blocks import declare_ring, Block
            sage: to_ring = declare_ring([Block("x", 10)], globals())
            sage: from_ring = declare_ring([Block("y", 5), Block("x", 10)], globals())
            sage: from sage.rings.polynomial.pbori.ll import RingMap
            sage: mapping = RingMap(to_ring, from_ring)
            sage: mapping(x(1)+1)
            x(1) + 1
        """
        def vars(ring):
            return [ring.variable(i) for i in range(ring.n_variables())]

        def indices(vars):
            return {str(v): idx for idx, v in enumerate(vars)}

        self.to_ring = to_ring
        self.from_ring = from_ring
        to_vars = vars(to_ring)
        from_vars = vars(from_ring)
        to_indices = indices(to_vars)
        from_indices = indices(from_vars)
        common = list(set(to_indices.keys()) & set(from_indices.keys()))

        to_map = list(from_vars)
        for elt in common:
            to_map[from_indices[elt]] = to_vars[to_indices[elt]]

        from_map = list(to_vars)
        for elt in common:
            from_map[to_indices[elt]] = from_vars[from_indices[elt]]

        self.to_map = BoolePolynomialVector(to_map)
        self.from_map = BoolePolynomialVector(from_map)

    def __call__(self, poly):
        r"""
        Execute the map to change rings.

        TESTS::

            sage: from sage.rings.polynomial.pbori.pbori import *
            sage: from sage.rings.polynomial.pbori.blocks import declare_ring, Block
            sage: to_ring = declare_ring([Block("x", 10)], globals())
            sage: from_ring = declare_ring([Block("y", 5), Block("x", 10)], globals())
            sage: from sage.rings.polynomial.pbori.ll import RingMap
            sage: mapping = RingMap(to_ring, from_ring)
            sage: mapping(x(1)+1)
            x(1) + 1
        """
        return substitute_variables(self.to_ring, self.to_map, poly)

    def invert(self, poly):
        r"""
        Inverted map to initial ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.pbori.pbori import *
            sage: from sage.rings.polynomial.pbori.blocks import declare_ring, Block
            sage: to_ring = declare_ring([Block("x", 10)], globals())
            sage: from_ring = declare_ring([Block("y", 5), Block("x", 10)], globals())
            sage: from sage.rings.polynomial.pbori.ll import RingMap
            sage: mapping = RingMap(to_ring, from_ring)
            sage: mapping.invert(mapping(x(1)+1))
            x(1) + 1
        """
        return substitute_variables(self.from_ring, self.from_map, poly)
