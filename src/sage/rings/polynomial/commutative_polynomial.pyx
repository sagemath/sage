# sage_setup: distribution = sagemath-categories
cdef class CommutativePolynomial(CommutativeAlgebraElement):
    r"""
    Abstract base class for commutative polynomials in any number of variables.

    It is a common base for :class:`~sage.rings.polynomial.polynomial_element.Polynomial`,
    :class:`~sage.rings.polynomial.multi_polynomial.MPolynomial`, and
    :class:`~sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial`.

    EXAMPLES::

        sage: from sage.rings.polynomial.commutative_polynomial import CommutativePolynomial
        sage: K.<x> = PolynomialRing(QQ)
        sage: isinstance(x, CommutativePolynomial)
        True
        sage: K.<x,y> = PolynomialRing(QQ)
        sage: isinstance(x, CommutativePolynomial)
        True
        sage: X.<x,y> = InfinitePolynomialRing(ZZ, implementation='sparse')
        sage: isinstance(x[2], CommutativePolynomial)
        True
    """

    pass
