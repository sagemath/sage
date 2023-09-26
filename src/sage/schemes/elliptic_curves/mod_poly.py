r"""
Modular polynomials for elliptic curves

For a positive integer `\ell`, the classical modular polynomial
`\Phi_\ell \in \ZZ[X,Y]` is characterized by the property that its
zero set is exactly the set of pairs of `j`-invariants connected
by a cyclic `\ell`-isogeny.

AUTHORS:

- Lorenz Panny (2023)
"""

from sage.structure.element import parent

from sage.rings.integer_ring import ZZ

from sage.libs.pari import pari
from cypari2.handle_error import PariError

from sage.databases.db_modular_polynomials import ClassicalModularPolynomialDatabase
_db = ClassicalModularPolynomialDatabase()

_cache = dict()

def classical_modular_polynomial(l, j=None):
    r"""
    Return the classical modular polynomial `\Phi_\ell`, either as a
    "generic" bivariate polynomial over `\ZZ`, or as an "instantiated"
    modular polynomial where one variable has been replaced by the
    given `j`-invariant.

    Generic polynomials are cached up to a certain size of `\ell`,
    which significantly accelerates subsequent invocations with the
    same `\ell`. The default bound is `\ell \leq 150`, which can be
    adjusted by setting ``classical_modular_polynomial.cache_bound``
    to a different value. Beware that modular polynomials are very
    large and the amount of memory consumed by the cache will grow
    rapidly when the bound is set to a large value.

    INPUT:

    - ``l`` -- positive integer.
    - ``j`` -- either ``None`` or a ring element:
      * if ``None`` is given, the original modular polynomial
        is returned as an element of `\ZZ[X,Y]`
      * if a ring element `j \in R` is given, the evaluation
        `\Phi_\ell(j,Y)` is returned as an element of the
        univariate polynomial ring `R[Y]`

    ALGORITHMS:

    - The Kohel database
      :class:`~sage.databases.db_modular_polynomials.ClassicalModularPolynomialDatabase`
    - :pari:`polmodular`

    EXAMPLES::

        sage: classical_modular_polynomial(2)
        -X^2*Y^2 + X^3 + 1488*X^2*Y + 1488*X*Y^2 + Y^3 - 162000*X^2 + 40773375*X*Y - 162000*Y^2 + 8748000000*X + 8748000000*Y - 157464000000000
        sage: j = Mod(1728, 419)
        sage: classical_modular_polynomial(3, j)
        Y^4 + 230*Y^3 + 84*Y^2 + 118*Y + 329

    TESTS::

        sage: q = random_prime(50)^randrange(1,4)
        sage: j = GF(q).random_element()
        sage: l = random_prime(50)
        sage: Y = polygen(parent(j), 'Y')
        sage: classical_modular_polynomial(l,j) == classical_modular_polynomial(l)(j,Y)
        True
    """
    l = ZZ(l)

    if j is None:
        # We are supposed to return the generic modular polynomial. First
        # check if it is already in the cache, then check the database,
        # finally compute it using PARI.
        try:
            return _cache[l]
        except KeyError:
            pass

        try:
            Phi = ZZ['X,Y'](_db[l])
        except ValueError:
            try:
                pari_Phi = pari.polmodular(l)
            except PariError:
                raise NotImplementedError('modular polynomial is not in database and computing it on the fly is not yet implemented')
            d = {(i, j): c for i,f in enumerate(pari_Phi) for j, c in enumerate(f)}
            Phi = ZZ['X,Y'](d)

        if l <= classical_modular_polynomial.cache_bound:
            _cache[l] = Phi

        return Phi

    R = parent(j)['Y']
    Y = R.gen()

    # If the generic polynomial is in the cache or the database, evaluating
    # it directly should always be faster than recomputing it from scratch.
    if l  in _cache:
        return _cache[l](j, Y)
    try:
        Phi = _db[l]
    except ValueError:
        pass
    else:
        if l <= classical_modular_polynomial.cache_bound:
            _cache[l] = ZZ['X,Y'](Phi)
        return Phi(j, Y)

    # Now try to get the instantiated modular polynomial directly from PARI.
    # This should be slightly more efficient (in particular regarding memory
    # usage) than computing and evaluating the generic modular polynomial.
    try:
        pari_Phi = pari.polmodular(l, 0, j)
    except PariError:
        pass
    else:
        return R(pari_Phi)

    # Nothing worked. Fall back to computing the generic modular polynomial
    # and simply evaluating it.
    return classical_modular_polynomial(l)(j, Y)

classical_modular_polynomial.cache_bound = 150
