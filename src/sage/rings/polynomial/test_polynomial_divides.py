def test_divides_zmod4():
    from sage.rings.finite_rings.constructor import Zmod
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    R = PolynomialRing(Zmod(4), 't')
    t = R.gen()
    a = 2*t**2 + t
    b = 2*t + 2
    c = a * b
    # a should divide c
    assert a.divides(c), f"{a} should divide {c} in Zmod(4)[t]"
    # a should not divide t
    assert not a.divides(t), f"{a} should not divide {t} in Zmod(4)[t]"
    # 0 only divides 0
    assert not R(0).divides(a)
    assert R(0).divides(R(0))
    # a divides 0
    assert a.divides(R(0))

def test_divides_zmod8():
    from sage.rings.finite_rings.constructor import Zmod
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    R = PolynomialRing(Zmod(8), 't')
    t = R.gen()
    a = 4*t**2 + t + 4
    b = 2
    c = a * b // t  # Laurent-like, but in polynomial ring
    # a should divide a * b
    assert a.divides(a * b)
    # a should not divide t
    assert not a.divides(t)
    # a should divide c if c is a polynomial
    if c.parent() is R:
        assert a.divides(c)

def test_divides_units_and_zero():
    from sage.rings.finite_rings.constructor import Zmod
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    R = PolynomialRing(Zmod(4), 't')
    t = R.gen()
    u = R(2)
    a = 2 + t
    # unit divides any polynomial
    assert u.divides(a)
    # everything divides 0
    assert a.divides(R(0))
    # 0 only divides 0
    assert not R(0).divides(a)
    assert R(0).divides(R(0)) 