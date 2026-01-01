from sage.all import SR
from sage.symbolic.heuristics import definitely_not_equal

def test_simple_not_equal():
    x = SR.var('x')
    assert definitely_not_equal(x + 1, x + 2)

def test_same_expression():
    x = SR.var('x')
    assert not definitely_not_equal(x^2, x*x)

