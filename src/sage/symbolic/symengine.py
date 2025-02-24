r"""
EXAMPLES::

    sage: # optional - symengine_py
    sage: import symengine
    sage: x, y = symengine.var("x y")
    sage: expr = (x + symengine.GoldenRatio * symengine.exp(y))**2
    sage: SR(expr)
    (golden_ratio*e^y + x)^2

    sage: # optional - symengine_py
    sage: f = symengine.Lambdify([x, y], expr)
    sage: f(3, 5)
    array(59115.86131768)
    sage: g = fast_callable(SR(expr), vars=[SR(x),SR(y)], domain=RDF)
    sage: g(3, 5)
    59115.86131767523
"""
