from sage.numerical.self_audit.api import n_audit

def _attach_n_audit():
    from sage.symbolic.expression import Expression
    Expression.n_audit = n_audit

_attach_n_audit()

