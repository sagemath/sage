R.<x> = PolynomialRing(GF(5))

def random_sample(J, n=1000, fast=True):
    p = [J.random_element(fast=fast) for _ in range(n)]
    return len(set(p))

def random_curve(genus=2):
    while True:
        # pick an f polynomial
        if randint(0, 1):
            d = 2*genus + 2
        else:
            d = 2*genus + 1
        f = R.random_element(degree=d)

        # Pick an h polynomial
        coin = randint(0, 2)
        if coin == 0:
            h = 0
        elif coin == 1:
            h = R.random_element(degree=(1, genus))
        # pick genus of H from degree of h 1/3 of the time
        else:
            f = R.random_element((1, 2*genus))
            h = R.random_element(degree = genus + 1)

        # Ensure that there are two points at infinity and the curve is non-singular
        try:
            H = HyperellipticCurveSmoothModel(f, h)
            # Cannot do arithmetic with odd genus curves which are inert
            if (genus % 2) and H.is_inert():
                continue
            assert H.genus() == genus
            return f, h, H
        except Exception as e:
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
for g in [2, 3, 4]:
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

