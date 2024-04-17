# sage_setup: distribution = sagemath-flint
"Monomials"

def _monomials(gens, R, n, i):
    """
    Given two lists ``gens`` and ``n`` of exactly the same length,
    return all monomials in the elements of ``gens`` in ``R`` where
    the ``i``-th generator in the monomial appears to degree strictly
    less than ``n[i]``.

    EXAMPLES::

        sage: monomials([x], [3])  # indirect doctest                                   # needs sage.symbolic
        [1, x, x^2]
    """
    # each power of the ith generator times all products
    # not involving the ith generator.
    if len(gens) == 1:
        b = gens[0]
        v = [R(1)]
        for _ in range(n[0]-1):
            v.append(v[-1]*b)
        return v
    else:
        z = gens[i]
        w = list(gens)
        del w[i]
        nn = list(n)
        del nn[i]
        v = monomials(w, nn)
        k = len(v)
        for _ in range(n[i]-1):
            for j in range(k):
                v.append(v[j]*z)
            z *= gens[i]
        return v


from sage.structure.sequence import Sequence


def monomials(v, n):
    """
    Given two lists ``v`` and ``n``, of exactly the same length,
    return all monomials in the elements of ``v``, where
    variable ``i`` (i.e., ``v[i]``) in the monomial appears to
    degree strictly less than ``n[i]``.

    INPUT:

    - ``v`` -- list of ring elements

    - ``n`` -- list of integers

    EXAMPLES::

        sage: monomials([x], [3])                                                       # needs sage.symbolic
        [1, x, x^2]
        sage: R.<x,y,z> = QQ[]
        sage: monomials([x,y], [5,5])
        [1, y, y^2, y^3, y^4, x, x*y, x*y^2, x*y^3, x*y^4, x^2, x^2*y, x^2*y^2, x^2*y^3, x^2*y^4, x^3, x^3*y, x^3*y^2, x^3*y^3, x^3*y^4, x^4, x^4*y, x^4*y^2, x^4*y^3, x^4*y^4]
        sage: monomials([x,y,z], [2,3,2])
        [1, z, y, y*z, y^2, y^2*z, x, x*z, x*y, x*y*z, x*y^2, x*y^2*z]
    """

    if (len(v) != len(n)):
        raise ValueError("inputs must be of the same length.")
    if len(v) == 0:
        return []
    v = Sequence(v)
    R = v.universe()
    return _monomials(v, R, n, 0)
