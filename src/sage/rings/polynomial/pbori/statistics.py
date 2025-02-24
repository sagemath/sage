from .pbori import top_index, BooleConstant
from .PyPolyBoRi import Monomial, Polynomial


def used_vars(l, bound=None):
    if not l:
        return BooleConstant(1)
    m = Monomial(Polynomial(next(iter(l))).vars_as_monomial())
    for p in l[1:]:
        m = m * Polynomial(p).vars_as_monomial()
        if bound and len(m) > bound:
            return m
    return m


def used_vars_set(l, bound=None):
    if not l:
        return BooleConstant(1)
    s = set()
    for p in l:
        s.update(Polynomial(p).vars_as_monomial().variables())
        if bound and len(s) > bound:
            break
    sorted_s = sorted(s, key=top_index, reverse=True)
    m = Monomial(next(iter(l)).ring())
    for v in sorted_s:
        m = v * m

    return m
