r"""
EXAMPLES::

    sage: import symengine                                             # optional - symengine_py
    sage: x, y = symengine.var("x y")                                  # optional - symengine_py
    sage: expr = (x + symengine.GoldenRatio * symengine.exp(y))**2     # optional - symengine_py
    sage: SR(expr)                                                     # optional - symengine_py
    (golden_ratio*e^y + x)^2

    sage: f = symengine.Lambdify([x, y], expr)                         # optional - symengine_py
    sage: f(3, 5)                                                      # optional - symengine_py
    array(59115.86131768)
    sage: g = fast_callable(SR(expr), vars=[SR(x),SR(y)], domain=RDF)  # optional - symengine_py
    sage: g(3, 5)                                                      # optional - symengine_py
    59115.86131767523
"""
