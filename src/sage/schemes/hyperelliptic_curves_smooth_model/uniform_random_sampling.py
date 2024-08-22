# code from https://github.com/sagemath/sage/pull/37118

# TODO:
# this is merged now... can we remove this?!

def uniform_random_polynomial(R, degree=(-1, 2), monic=False, *args, **kwds):
    r"""
    Return a random polynomial of given degree (bounds).

    INPUT:

    -  ``degree`` -- (default: ``(-1, 2)``) integer for fixing the degree or
        a tuple of minimum and maximum degrees

    -  ``monic`` -- boolean (optional); indicate whether the sampled
        polynomial should be monic

    -  ``*args, **kwds`` -- additional keyword parameters passed on to the
        ``random_element`` method for the base ring
    """
    k = R.base_ring()

    if isinstance(degree, (list, tuple)):
        if len(degree) != 2:
            raise ValueError(
                "degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)"
            )
        if degree[0] > degree[1]:
            raise ValueError("minimum degree must be less or equal than maximum degree")
        if degree[1] < -1:
            raise ValueError(f"maximum degree (={degree[1]}) must be at least -1")
    else:
        if degree < -1:
            raise ValueError(f"degree (={degree}) must be at least -1")
        degree = (degree, degree)

    if degree[0] <= -2:
        degree = (-1, degree[1])

    # If the coefficient range only contains 0, then
    # * if the degree range includes -1, return the zero polynomial,
    # * otherwise raise a value error
    if args == (0, 1):
        if degree[0] == -1:
            return R.zero()
        else:
            raise ValueError("No polynomial of degree >= 0 has all coefficients zero")

    if degree == (-1, -1):
        return R.zero()

    # If `monic` is set, zero should be ignored
    if degree[0] == -1 and monic:
        if degree[1] == -1:
            raise ValueError(
                "the maximum degree of monic polynomials needs to be at least 0"
            )
        if degree[1] == 0:
            return R.one()
        degree = (0, degree[1])

    # Pick random coefficients
    end = degree[1]
    if degree[0] == -1:
        return R([k.random_element(*args, **kwds) for _ in range(end + 1)])

    nonzero = False
    coefs = [None] * (end + 1)

    while not nonzero:
        # Pick leading coefficients, if `monic` is set it's handle here.
        if monic:
            for i in range(degree[1] - degree[0] + 1):
                coefs[end - i] = k.random_element(*args, **kwds)
                if not nonzero and not coefs[end - i].is_zero():
                    coefs[end - i] = k.one()
                    nonzero = True
        else:
            # Fast path
            for i in range(degree[1] - degree[0] + 1):
                coefs[end - i] = k.random_element(*args, **kwds)
                nonzero |= not coefs[end - i].is_zero()

    # Now we pick the remaining coefficients.
    for i in range(degree[1] - degree[0] + 1, degree[1] + 1):
        coefs[end - i] = k.random_element(*args, **kwds)

    return R(coefs)
