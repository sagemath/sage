from sage.all import QQ

def definitely_not_equal(expr1, expr2, samples=3):
    """
    Heuristic check to detect that two symbolic expressions are not equal.
    """
    vars1 = set(expr1.variables())
    vars2 = set(expr2.variables())
    if vars1 != vars2:
        return False

    vars = list(vars1)

    for i in range(1, samples + 1):
        subs = {v: QQ(i) for v in vars}
        try:
            v1 = expr1.subs(subs)
            v2 = expr2.subs(subs)
        except Exception:
            continue

        if v1 != v2:
            return True

    return False

