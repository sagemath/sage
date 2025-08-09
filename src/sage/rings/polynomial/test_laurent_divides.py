def test_divides_basic():
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    R = LaurentPolynomialRing('y', base_ring=Zmod(4))
    y = R.gen()
    a = 2 + y
    b = 2
    c = a * b
    assert a.divides(c), f"{a} should divide {c}"

def test_divides_zero():
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    R = LaurentPolynomialRing('y', base_ring=Zmod(4))
    y = R.gen()
    a = 2 + y
    assert a.divides(R(0))
    assert not R(0).divides(a)
    assert R(0).divides(R(0))

def test_divides_unit():
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    R = LaurentPolynomialRing('y', base_ring=Zmod(4))
    y = R.gen()
    u = R(2)
    a = 2 + y
    assert u.divides(a)

def test_divides_monomial():
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    R = LaurentPolynomialRing('y', base_ring=Zmod(4))
    y = R.gen()
    a = y
    c = y**5
    assert a.divides(c)
    assert not c.divides(a) 