../conf_sub.py


def fast_callable(expression, vars=None, domain=None, hold=False, use_cache=False):
    """
    Convert a symbolic expression into a fast callable function.

    This function takes a symbolic expression and converts it into a more 
    efficient, fast callable function. This is useful for repeated evaluation 
    of the expression, especially in performance-critical applications.

    INPUT:
        - `expression` (required): 
          A symbolic expression that you want to convert into a fast callable function.
          
        - `vars` (optional, default: None): 
          A list or tuple of variables (as symbolic variables) with respect to which 
          the expression should be made callable. If `None`, the function assumes all 
          variables in the expression.
          
        - `domain` (optional, default: None): 
          Specifies the domain of the function. Common options include:
            - `"RR"` for real numbers.
            - `"CC"` for complex numbers.
            - `"ZZ"` for integers.
          If `None`, the domain is inferred from the expression.
          
        - `hold` (optional, default: False): 
          A boolean flag. If `True`, the evaluation of the expression is delayed 
          until the function is explicitly called. If `False`, the function is 
          prepared immediately.

        - `use_cache` (optional, default: False):
          A boolean flag. If `True`, the function will cache results to improve 
          performance for repeated calls with the same input. If `False`, caching 
          is disabled.

    OUTPUT:
        - A callable function that can be used to evaluate the expression for 
          specific values of the variables.
          
        The returned function will have the following properties:
            - It accepts arguments corresponding to the variables in `vars`.
            - It returns a numerical value corresponding to the evaluation of 
              `expression` for the provided inputs.

    EXAMPLES:
        sage: f = fast_callable(x^2 + y^2, vars=[x, y], domain="RR")
        sage: f(2, 3)
        13.0

        sage: g = fast_callable(sin(x) + cos(y), vars=[x, y], domain="RR")
        sage: g(pi/2, 0)
        1.0
    """
