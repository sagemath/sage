"""
Routines for Conway and pseudo-Conway polynomials

AUTHORS:

- David Roe

- Jean-Pierre Flori

- Peter Bruin
"""

from sage.misc.fast_methods import WithEqualityById
from sage.misc.lazy_import import lazy_import
from sage.structure.sage_object import SageObject
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.integer import Integer

lazy_import('sage.databases.conway', 'ConwayPolynomials')


def conway_polynomial(p, n):
    """
    Return the Conway polynomial of degree `n` over ``GF(p)``.

    If the requested polynomial is not known, this function raises a
    ``RuntimeError`` exception.

    INPUT:

    - ``p`` -- prime number

    - ``n`` -- positive integer

    OUTPUT:

    - the Conway polynomial of degree `n` over the finite field
      ``GF(p)``, loaded from a table.

    .. NOTE::

       The first time this function is called a table is read from
       disk, which takes a fraction of a second. Subsequent calls do
       not require reloading the table.

    See also the ``ConwayPolynomials()`` object, which is the table of
    Conway polynomials used by this function.

    EXAMPLES::

        sage: conway_polynomial(2,5)                                                    # needs conway_polynomials
        x^5 + x^2 + 1
        sage: conway_polynomial(101,5)                                                  # needs conway_polynomials
        x^5 + 2*x + 99
        sage: conway_polynomial(97,101)                                                 # needs conway_polynomials
        Traceback (most recent call last):
        ...
        RuntimeError: requested Conway polynomial not in database.
    """
    (p, n) = (int(p), int(n))
    R = FiniteField(p)['x']
    try:
        return R(ConwayPolynomials()[p][n])
    except KeyError:
        raise RuntimeError("requested Conway polynomial not in database.")

def exists_conway_polynomial(p, n):
    """
    Check whether the Conway polynomial of degree `n` over ``GF(p)``
    is known.

    INPUT:

    - ``p`` -- prime number

    - ``n`` -- positive integer

    OUTPUT:

    - boolean: ``True`` if the Conway polynomial of degree `n` over
      ``GF(p)`` is in the database, ``False`` otherwise.

    If the Conway polynomial is in the database, it can be obtained
    using the command ``conway_polynomial(p,n)``.

    EXAMPLES::

        sage: exists_conway_polynomial(2,3)                                             # needs conway_polynomials
        True
        sage: exists_conway_polynomial(2,-1)
        False
        sage: exists_conway_polynomial(97,200)
        False
        sage: exists_conway_polynomial(6,6)
        False
    """
    try:
        return ConwayPolynomials().has_polynomial(p,n)
    except ImportError:
        return False

class PseudoConwayLattice(WithEqualityById, SageObject):
    r"""
    A pseudo-Conway lattice over a given finite prime field.

    The Conway polynomial `f_n` of degree `n` over `\Bold{F}_p` is
    defined by the following four conditions:

    - `f_n` is irreducible.

    - In the quotient field `\Bold{F}_p[x]/(f_n)`, the element
      `x\bmod f_n` generates the multiplicative group.

    - The minimal polynomial of `(x\bmod f_n)^{\frac{p^n-1}{p^m-1}}`
      equals the Conway polynomial `f_m`, for every divisor `m` of
      `n`.

    - `f_n` is lexicographically least among all such polynomials,
      under a certain ordering.

    The final condition is needed only in order to make the Conway
    polynomial unique.  We define a pseudo-Conway lattice to be any
    family of polynomials, indexed by the positive integers,
    satisfying the first three conditions.

    INPUT:

    - ``p`` -- prime number

    - ``use_database`` -- boolean.  If ``True``, use actual Conway
      polynomials whenever they are available in the database.  If
      ``False``, always compute pseudo-Conway polynomials.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
        sage: PCL = PseudoConwayLattice(2, use_database=False)
        sage: PCL.polynomial(3)
        x^3 + x + 1

    TESTS::

        sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
        sage: PCL = PseudoConwayLattice(3)
        sage: hash(PCL)  # random
        8738829832350

        sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
        sage: PseudoConwayLattice(3) == PseudoConwayLattice(3)
        False
        sage: PseudoConwayLattice(3) != PseudoConwayLattice(3)
        True
        sage: P = PseudoConwayLattice(5)
        sage: P == P
        True
        sage: P != P
        False
    """
    def __init__(self, p, use_database=True):
        """
        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
            sage: PCL = PseudoConwayLattice(3)
            sage: PCL.polynomial(3)
            x^3 + 2*x + 1

            sage: # needs sage.rings.finite_rings
            sage: PCL = PseudoConwayLattice(5, use_database=False)
            sage: PCL.polynomial(12)
            x^12 + 4*x^11 + 2*x^10 + 4*x^9 + 2*x^8 + 2*x^7 + 4*x^6 + x^5 + 2*x^4 + 2*x^2 + x + 2
            sage: PCL.polynomial(6)
            x^6 + x^5 + 4*x^4 + 3*x^3 + 3*x^2 + 2*x + 2
            sage: PCL.polynomial(11)
            x^11 + x^6 + 3*x^3 + 4*x + 3
        """
        self.p = p
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        self.ring = PolynomialRing(FiniteField(p), 'x')
        if use_database:
            try:
                C = ConwayPolynomials()
            except ImportError:
                self.nodes = {}
            else:
                self.nodes = {n: self.ring(C.polynomial(p, n))
                              for n in C.degrees(p)}
        else:
            self.nodes = {}

    def polynomial(self, n):
        r"""
        Return the pseudo-Conway polynomial of degree `n` in this
        lattice.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT:

        - a pseudo-Conway polynomial of degree `n` for the prime `p`.

        ALGORITHM:

        Uses an algorithm described in [HL1999]_, modified to find
        pseudo-Conway polynomials rather than Conway polynomials.  The
        major difference is that we stop as soon as we find a
        primitive polynomial.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
            sage: PCL = PseudoConwayLattice(2, use_database=False)
            sage: PCL.polynomial(3)
            x^3 + x + 1
            sage: PCL.polynomial(4)
            x^4 + x^3 + 1
            sage: PCL.polynomial(60)
            x^60 + x^59 + x^58 + x^55 + x^54 + x^53 + x^52 + x^51 + x^48 + x^46 + x^45 + x^42 + x^41 + x^39 + x^38 + x^37 + x^35 + x^32 + x^31 + x^30 + x^28 + x^24 + x^22 + x^21 + x^18 + x^17 + x^16 + x^15 + x^14 + x^10 + x^8 + x^7 + x^5 + x^3 + x^2 + x + 1
        """
        if n in self.nodes:
            return self.nodes[n]

        p = self.p
        n = Integer(n)

        if n == 1:
            f = self.ring.gen() - FiniteField(p).multiplicative_generator()
            self.nodes[1] = f
            return f

        # Work in an arbitrary field K of order p**n.
        K = FiniteField(p**n, names='a')

        # TODO: something like the following
        # gcds = [n.gcd(d) for d in self.nodes.keys()]
        # xi = { m: (...) for m in gcds }
        xi = {q: self.polynomial(n//q).any_root(K, -n//q, assume_squarefree=True)
              for q in n.prime_divisors()}

        # The following is needed to ensure that in the concrete instantiation
        # of the "new" extension all previous choices are compatible.
        _frobenius_shift(K, xi)

        # Construct a compatible element having order the lcm of orders
        q, x = xi.popitem()
        v = p**(n//q) - 1
        for q, xitem in xi.items():
            w = p**(n//q) - 1
            g, alpha, beta = v.xgcd(w)
            x = x**beta * xitem**alpha
            v = v.lcm(w)

        r = p**n - 1
        # Get the missing part of the order to be primitive
        g = r // v
        # Iterate through g-th roots of x until a primitive one is found
        z = x.nth_root(g)
        root = K.multiplicative_generator()**v
        while z.multiplicative_order() != r:
            z *= root
        # The following should work but tries to create a huge list
        # whose length overflows Python's ints for large parameters
        #Z = x.nth_root(g, all=True)
        #for z in Z:
        #    if z.multiplicative_order() == r:
        #         break
        f = z.minimal_polynomial()
        self.nodes[n] = f
        return f

    def check_consistency(self, n):
        """
        Check that the pseudo-Conway polynomials of degree dividing
        `n` in this lattice satisfy the required compatibility
        conditions.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
            sage: PCL = PseudoConwayLattice(2, use_database=False)
            sage: PCL.check_consistency(6)
            sage: PCL.check_consistency(60)  # long time
        """
        p = self.p
        K = FiniteField(p**n, modulus=self.polynomial(n), names='a')
        a = K.gen()
        for m in n.divisors():
            assert (a**((p**n-1)//(p**m-1))).minimal_polynomial() == self.polynomial(m)


def _find_pow_of_frobenius(p, n, x, y):
    """
    Find the power of Frobenius which yields `x` when applied to `y`.

    INPUT:

    - ``p`` -- prime number

    - ``n`` -- positive integer

    - ``x`` -- an element of a field `K` of `p^n` elements so that
      the multiplicative order of `x` is `p^n - 1`.

    - ``y`` -- an element of `K` with the same minimal polynomial as
      `x`.

    OUTPUT:

    - an element `i` of the integers modulo `n` such that `x = y^{p^i}`.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.rings.finite_rings.conway_polynomials import _find_pow_of_frobenius
        sage: K.<a> = GF(3^14)
        sage: x = K.multiplicative_generator()
        sage: y = x^27
        sage: _find_pow_of_frobenius(3, 14, x, y)
        11

    """
    from .integer_mod import mod
    for i in range(n):
        if x == y:
            break
        y = y**p
    else:
        raise RuntimeError("No appropriate power of Frobenius found")
    return mod(i, n)

def _crt_non_coprime(running, a):
    """
    Extension of the ``crt`` method of ``IntegerMod`` to the case of
    non-relatively prime modulus.

    EXAMPLES::

        sage: from sage.rings.finite_rings.conway_polynomials import _crt_non_coprime
        sage: a = _crt_non_coprime(mod(14, 18), mod(20,30)); a
        50
        sage: a.modulus()
        90
        sage: _crt_non_coprime(mod(13, 18), mod(20,30))
        Traceback (most recent call last):
        ...
        AssertionError

    """
    g = running.modulus().gcd(a.modulus())
    if g == 1:
        return running.crt(a)
    else:
        assert running % g == a % g
        running_modulus = running.modulus()
        a_modulus = a.modulus()
        for qq in g.prime_divisors():
            a_val_unit = a_modulus.val_unit(qq)
            running_val_unit = running_modulus.val_unit(qq)
            if a_val_unit[0] > running_val_unit[0]:
                running_modulus = running_val_unit[1]
            else:
                a_modulus = a_val_unit[1]
        return (running % running_modulus).crt(a % a_modulus)

def _frobenius_shift(K, generators, check_only=False):
    """
    Given a field `K` of degree `n` over ``GF(p)`` and a dictionary
    holding, for each divisor `q` of `n`, an element with minimal
    polynomial a pseudo-Conway polynomial of degree `n/q`, modify
    these generators into a compatible system.

    Such a system of generators is said to be compatible if for each
    pair of prime divisors `q_1` and `q_2` and each common divisor `m`
    of `n/q_1` and `n/q_2`, the equality

    ``generators[q1]^((p^(n/q1)-1)/(p^m-1)) == generators[q2]^((p^(n/q2)-1)/(p^m-1))``

    holds.

    INPUT:

    - ``K`` -- a finite field of degree `n` over its prime field

    - ``generators`` -- a dictionary, indexed by prime divisors `q` of
      `n`, whose entries are elements of `K` satisfying the `n/q`
      pseudo-Conway polynomial.

    - ``check_only`` -- if ``True``, just check that the given
      generators form a compatible system.

    EXAMPLES::

        sage: # needs sage.libs.ntl sage.rings.finite_rings
        sage: R.<x> = GF(2)[]
        sage: f30 = x^30 + x^28 + x^27 + x^25 + x^24 + x^20 + x^19 + x^18 + x^16 + x^15 + x^12 + x^10 + x^7 + x^2 + 1
        sage: f20 = x^20 + x^19 + x^15 + x^13 + x^12 + x^11 + x^9 + x^8 + x^7 + x^4 + x^2 + x + 1
        sage: f12 = x^12 + x^10 + x^9 + x^8 + x^4 + x^2 + 1
        sage: K.<a> = GF(2^60, modulus='first_lexicographic')
        sage: x30 = f30.any_root(K)
        sage: x20 = f20.any_root(K)
        sage: x12 = f12.any_root(K)
        sage: generators = {2: x30, 3: x20, 5: x12}
        sage: from sage.rings.finite_rings.conway_polynomials import _frobenius_shift, _find_pow_of_frobenius
        sage: _frobenius_shift(K, generators)
        sage: _find_pow_of_frobenius(2, 30, x30, generators[2])
        0
        sage: _find_pow_of_frobenius(2, 20, x20, generators[3])
        13
        sage: _find_pow_of_frobenius(2, 12, x12, generators[5])
        8

    """
    if len(generators) == 1:
        return generators
    p = K.characteristic()
    n = K.degree()
    compatible = {}
    from .integer_mod import mod
    for m in n.divisors():
        compatible[m] = {}
    for q, x in generators.items():
        for m in (n//q).divisors():
            compatible[m][q] = x**((p**(n//q)-1)//(p**m-1))
    if check_only:
        for m in n.divisors():
            try:
                q, x = compatible[m].popitem()
            except KeyError:
                break
            for qq, xx in compatible[m].items():
                assert x == xx
        return
    crt = {}
    qlist = sorted(generators.keys())
    for j in range(1, len(qlist)):
        for i in range(j):
            crt[(i, j)] = []
    for m in n.divisors():
        mqlist = sorted(compatible[m].keys())
        for k in range(1,len(mqlist)):
            j = qlist.index(mqlist[k])
            i = qlist.index(mqlist[k-1])
            crt[(i,j)].append(_find_pow_of_frobenius(p, m, compatible[m][qlist[j]], compatible[m][qlist[i]]))
    for i, j in list(crt):
        L = crt[(i,j)]
        running = mod(0, 1)
        for a in L:
            running = _crt_non_coprime(running, a)
        crt[(i,j)] = [(mod(running, qq**(running.modulus().valuation(qq))),
                       running.modulus().valuation(qq)) for qq in qlist]
        crt[(j,i)] = [(-a, level) for a, level in crt[(i,j)]]
    # Let x_j be the power of Frobenius we apply to generators[qlist[j]], for 0 < j < len(qlist)
    # We have some direct conditions on the x_j: x_j reduces to each entry in crt[(0,j)].
    # But we also have the equations x_j - x_i reduces to each entry in crt[(i,j)].
    # We solve for x_j one prime at a time.  For each prime, we have an equations of the form
    # x_j - x_i = c_ij.  The modulus of the currently known value of x_j, x_i and c_ij will all be powers
    # (possibly 0, possibly different) of the same prime.

    # We can set x_0=0 everywhere, can get an initial setting of x_j from the c_0j.
    # We go through prime by prime.
    import bisect
    frob_powers = [mod(0, 1) for _ in qlist]

    def find_leveller(qindex, level, x, xleveled, searched, i):
        searched[i] = True
        crt_possibles = []
        for j in range(1,len(qlist)):
            if i==j:
                continue
            if crt[(i,j)][qindex][1] >= level:
                if xleveled[j]:
                    return [j]
                elif j not in searched:
                    crt_possibles.append(j)
        for j in crt_possibles:
            path = find_leveller(qindex, level, x, xleveled, searched, j)
            if path is not None:
                path.append(j)
                return path
        return None

    def propagate_levelling(qindex, level, x, xleveled, i):
        for j in range(1, len(qlist)):
            if i==j:
                continue
            if not xleveled[j] and crt[(i,j)][qindex][1] >= level:
                newxj = x[i][0] + crt[(i,j)][qindex][0]
                x[j] = (newxj, min(x[i][1], crt[(i,j)][qindex][1]))
                xleveled[j] = True
                propagate_levelling(qindex, level, x, xleveled, j)

    for qindex in range(len(qlist)):
        q = qlist[qindex]
        # We include the initial 0 to match up our indexing with crt.
        x = [0] + [crt[(0,j)][qindex] for j in range(1,len(qlist))]
        # We first check that our equations are consistent and
        # determine which powers of q occur as moduli.
        levels = []
        for j in range(2, len(qlist)):
            for i in range(j):
                # we need crt[(0,j)] = crt[(0,i)] + crt[(i,j)]
                if i != 0:
                    assert x[j][0] == x[i][0] + crt[(i,j)][qindex][0]
                level = crt[(i,j)][qindex][1]
                if level > 0:
                    ins = bisect.bisect_left(levels,level)
                    if ins == len(levels):
                        levels.append(level)
                    elif levels[ins] != level:
                        levels.insert(ins, level)
        for level in levels:
            xleveled = [0] + [x[i][1] >= level for i in range(1,len(qlist))]
            while True:
                try:
                    i = xleveled.index(False, 1)
                    searched = {}
                    levelling_path = find_leveller(qindex, level, x, xleveled, searched, i)
                    if levelling_path is None:
                        # Any lift will work, since there are no constraints.
                        x[i] = (mod(x[i][0].lift(), q**level), level)
                        xleveled[i] = True
                        propagate_levelling(qindex, level, x, xleveled, i)
                    else:
                        levelling_path.append(i)
                        for m in range(1,len(path)):
                            # This point on the path may have already
                            # been leveled in a previous propagation.
                            if not xleveled[path[m]]:
                                newx = x[path[m-1]][0] + crt[(path[m-1],path[m])][qindex][0]
                                x[path[m]] = (newx, min(x[path[m-1]][1], crt[(path[m-1],path[m])][qindex][1]))
                                xleveled[path[m]] = True
                                propagate_levelling(qindex, level, x, xleveled, path[m])
                except ValueError:
                    break
        for j in range(1,len(qlist)):
            frob_powers[j] = frob_powers[j].crt(x[j][0])
    for j in range(1, len(qlist)):
        generators[qlist[j]] = generators[qlist[j]]**(p**(-frob_powers[j]).lift())
    _frobenius_shift(K, generators, check_only=True)
