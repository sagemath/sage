R1.<x> = PolynomialRing(GF(163))
R2.<x> = PolynomialRing(GF(2^6))
R3.<x> = PolynomialRing(QQ)

for R in [R1, R2, R3]:
    for _ in range(100):
        f = R.random_element(degree=(-1, 10))
        h = R.random_element(degree=(-1, 10))
        failed = False
        try:
            H = HyperellipticCurve(f, h)
        except:
            failed = True
        try:
            H = HyperellipticCurveSmoothModel(f, h)
            assert not failed, "Oh no"
        except:
            assert failed, f"Whoops!"

    print(f"Test passed for {R}")
