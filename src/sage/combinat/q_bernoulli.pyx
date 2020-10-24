"""
`q`-Bernoulli Numbers and Polynomials
"""
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.cachefunc import cached_function


@cached_function
def q_bernoulli(m, p=None):
    r"""
    Compute Carlitz's `q`-analogue of the Bernoulli numbers.

    For every nonnegative integer `m`, the `q`-Bernoulli number
    `\beta_m` is a rational function of the indeterminate `q` whose
    value at `q=1` is the usual Bernoulli number `B_m`.

    INPUT:

    - `m` -- a nonnegative integer

    - `p` (default: ``None``) -- an optional value for `q`

    OUTPUT:

    A rational function of the indeterminate `q` (if `p` is ``None``)

    Otherwise, the rational function is evaluated at `p`.

    EXAMPLES::

        sage: from sage.combinat.q_bernoulli import q_bernoulli
        sage: q_bernoulli(0)
        1
        sage: q_bernoulli(1)
        -1/(q + 1)
        sage: q_bernoulli(2)
        q/(q^3 + 2*q^2 + 2*q + 1)
        sage: all(q_bernoulli(i)(q=1) == bernoulli(i) for i in range(12))
        True

    One can evaluate the rational function by giving a second argument::

        sage: x = PolynomialRing(GF(2),'x').gen()
        sage: q_bernoulli(5,x)
        x/(x^6 + x^5 + x + 1)

    The function does not accept negative arguments::

        sage: q_bernoulli(-1)
        Traceback (most recent call last):
        ...
        ValueError: the argument must be a nonnegative integer

    REFERENCES:

    .. [Ca1948] Leonard Carlitz, "q-Bernoulli numbers and polynomials". Duke Math J. 15, 987-1000 (1948), :doi:`10.1215/S0012-7094-48-01588-9`

    .. [Ca1954] Leonard Carlitz, "q-Bernoulli and Eulerian numbers". Trans Am Soc. 76, 332-350 (1954), :doi:`10.1090/S0002-9947-1954-0060538-2`
    """
    m = ZZ(m)
    if m < 0:
        raise ValueError("the argument must be a nonnegative integer")
    cdef int i
    q = PolynomialRing(ZZ, 'q').gen()
    result = (q - 1)**(1-m) * sum((-1)**(m-i)*m.binomial(i)*(i+1)/(q**(i+1)-1)
                                  for i in range(m + 1))
    if p is None:
        return result
    else:
        return result(q=p)


@cached_function
def q_bernoulli_polynomial(m):
    r"""
    Compute Carlitz's `q`-analogue of the Bernoulli polynomials.

    For every nonnegative integer `m`, the `q`-Bernoulli polynomial
    is a polynomial in one variable `x` with coefficients in `\QQ(q)` whose
    value at `q=1` is the usual Bernoulli polynomial `B_m(x)`.

    The original `q`-Bernoulli polynomials introduced by Carlitz were
    polynomials in `q^y` with coefficients in `\QQ(q)`. This function returns
    these polynomials but expressed in the variable `x=(q^y-1)/(q-1)`. This
    allows to let `q=1` to recover the classical Bernoulli polynomials.

    INPUT:

    - `m` -- a nonnegative integer

    OUTPUT:

    A polynomial in one variable `x`.

    EXAMPLES::

        sage: from sage.combinat.q_bernoulli import q_bernoulli_polynomial, q_bernoulli
        sage: q_bernoulli_polynomial(0)
        1
        sage: q_bernoulli_polynomial(1)
        (2/(q + 1))*x - 1/(q + 1)
        sage: x = q_bernoulli_polynomial(1).parent().gen()
        sage: all(q_bernoulli_polynomial(i)(q=1)==bernoulli_polynomial(x,i) for i in range(12))
        True
        sage: all(q_bernoulli_polynomial(i)(x=0)==q_bernoulli(i) for i in range(12))
        True

    The function does not accept negative arguments::

        sage: q_bernoulli_polynomial(-1)
        Traceback (most recent call last):
        ...
        ValueError: the argument must be a nonnegative integer

    REFERENCES: [Ca1948]_, [Ca1954]_
    """
    m = ZZ(m)
    if m < 0:
        raise ValueError("the argument must be a nonnegative integer")
    cdef int i
    R = PolynomialRing(ZZ, 'q')
    q = R.gen()
    x = PolynomialRing(R, 'x').gen()
    z = 1 + (q - 1) * x
    return sum(m.binomial(i) * x ** (m - i) * z ** i * q_bernoulli(i)
               for i in range(m + 1))
