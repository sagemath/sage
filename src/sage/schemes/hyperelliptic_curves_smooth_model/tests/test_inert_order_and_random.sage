R.<x> = PolynomialRing(GF(13))

def random_sample(J, n=1000, fast=True):
    p = [J.random_element(fast=fast) for _ in range(n)]
    return len(set(p))

def random_curve(use_h=True, genus=2):
    d = 2*genus + 2
    for _ in range(10):
        # Find a polynomial for f
        f = R.random_element(degree=d)

        # Find a polynomial for h
        if use_h:
            if R.base_ring().characteristic() == 2:
                h = R.random_element(degree=(genus + 1))
            else:
                h = R.random_element(degree=2)
        else:
            h = 0

        # Ensure that there are two points at infinity and the curve is non-singular
        try:
            H = HyperellipticCurveSmoothModel(f, h)
            if H.genus() != genus:
                continue
            if not H.is_inert():
                continue
            return f, h, H
        except:
            continue

# Test that randomly sampling gets all elements in the group
for _ in range(1):
    f, h, H = random_curve(genus=2)
    J = H.jacobian()
    o = J.order()

    fast = random_sample(J, fast=True)
    slow = random_sample(J, fast=False)

    print(f"{f = }")
    print(f"{h = }")
    print(f"Order: {o}")
    print(f"Fast method: {fast}")
    print(f"Slow method: {slow}")
    print(f"")

# Test all points have order dividing the Jacobian order
for g in [2, 4, 6]:
    print(f"Testing arithmetic for genus: {g}")
    for _ in range(5):
        f, h, H = random_curve(genus=g)
        J = H.jacobian()
        o = J.order()

        # Test order
        assert all([(o * J.random_element()).is_zero() for _ in range(100)])

        # Test order on divisor
        for _ in range(100):
            D = J.random_element()
            order = D.order()
            assert order * D == J.zero()

        # Test inversion on non-rational
        for _ in range(100):
            D = J.random_element()
            assert (D - D).is_zero()

    print(f"Passed for genus: {g}")

