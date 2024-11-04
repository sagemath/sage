# sage_setup: distribution = sagemath-brial
# sage.doctest: needs sage.rings.polynomial.pbori
def input_signals(p):
    return list((p + p.lex_lead()).vars_as_monomial().variables())


def output_signal(p):
    return next(iter(p.lex_lead().variables()))


def rank(data):
    parents = {}
    res = {}
    for p in data:
        out = output_signal(p)
        parents.setdefault(out, [])
        for v in input_signals(p):
            parents.setdefault(v, []).append(out)

    def do_rank(v):
        if v in res:
            return res[v]
        my_res = res[v] = max((do_rank(p) + 1 for p in parents[v]), default=0)
        return my_res
    for v in parents:
        do_rank(v)
    return res
